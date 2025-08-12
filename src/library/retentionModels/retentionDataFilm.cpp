#include "retentionDataFilm.h"
#include "soilMath.h"

using namespace Soil::RetentionModels;

namespace Soil::RetentionModels
{

    retentionDataFilm::retentionDataFilm(IOdictionary &transportProperties, volScalarField &region, fvMesh &mesh, Time &runTime, simConfig &config)
        : retentionData(transportProperties, region, mesh, runTime, config),
          shp_h_0(IOobject("shp_h_0", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("shp_h_0", dimensionSet(0, 1, 0, 0, 0, 0, 0), h_0)),
          tmp_one(IOobject("tmp_one", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("tmp_one", dimensionSet(0, 0, 0, 0, 0, 0, 0), 1.0)),
          tmp_zero(IOobject("tmp_zero", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("tmp_zero", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0)),
          shp_film_alpha(IOobject("shp_film_alpha", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("shp_film_alpha", dimensionSet(0, 0, 0, 0, 0, 0, 0), alpha)),
          shp_film_h_a(IOobject("shp_film_h_a", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("shp_film_h_a", dimensionSet(0, 1, 0, 0, 0, 0, 0), 0.0)),
          shp_film_k_sat(IOobject("shp_film_k_sat", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("shp_film_k_sat", dimensionSet(0, 1, -1, 0, 0, 0, 0), 0.0)),
          b(IOobject("b", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("b", dimensionSet(0, 0, 0, 0, 0, 0, 0), 1))
    {
        if (no_of_regions > 0)
		{
			char buffer[1000];

            shp_film_alpha_vector.reserve(no_of_regions);
			for (int i = 0; i < no_of_regions; i++)
			{
				sprintf(buffer, "shp_film_alpha_%d", i);
				std::string key_name(buffer);
				shp_film_alpha_vector.push_back(dimensionedScalar(key_name.c_str(), dimensionSet(0, 0, 0, 0, 0, 0, 0), alpha, transportProperties));
			}

            //FIXME: How to enforce usage of the transportProperties dictionary  key shp_film_h_a_N for all regions or any region?
            shp_film_h_a_vector.reserve(no_of_regions);
            for (int i = 0; i < no_of_regions; i++) 
            {
                sprintf(buffer, "shp_film_h_a_%d", i);
                std::string key_name(buffer);
                shp_film_h_a_vector.push_back(dimensionedScalar(key_name.c_str(), dimensionSet(0, 1, 0, 0, 0, 0, 0), 0.0, transportProperties));
            }

            shp_film_k_sat_vector.reserve(no_of_regions);
            for (int i = 0; i < no_of_regions; i++) 
            {
                sprintf(buffer, "shp_film_k_sat_%d", i);
                std::string key_name(buffer);
                shp_film_k_sat_vector.push_back(dimensionedScalar(key_name.c_str(), dimensionSet(0, 1, -1, 0, 0, 0, 0), 0.0, transportProperties));
            }
            
			forAll(region, index)
			{
                shp_film_alpha[index] = shp_film_alpha_vector[static_cast<int>(region[index])].value();
                shp_film_h_a[index] = shp_film_h_a_vector[static_cast<int>(region[index])].value();
                shp_film_k_sat[index] = shp_film_k_sat_vector[static_cast<int>(region[index])].value();
			}

			// initiate boundary fields
			forAll(region.boundaryFieldRef(), ipatch)
			{
				// each face in the the patch
				forAll(region.boundaryFieldRef()[ipatch], index)
				{
                    shp_film_alpha.boundaryFieldRef()[ipatch][index] = shp_film_alpha_vector[static_cast<int>(region.boundaryField()[ipatch][index])].value();
                    shp_film_h_a.boundaryFieldRef()[ipatch][index] = shp_film_h_a_vector[static_cast<int>(region.boundaryField()[ipatch][index])].value();
                    shp_film_k_sat.boundaryFieldRef()[ipatch][index] = shp_film_k_sat_vector[static_cast<int>(region.boundaryField()[ipatch][index])].value();
				}
			}
		}
		else
		{
            shp_film_alpha = volScalarField(IOobject("shp_film_alpha", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
            shp_film_h_a = volScalarField(IOobject("shp_film_h_a", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
			shp_film_k_sat = volScalarField(IOobject("shp_film_k_sat", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
		}
    }

    volScalarField& retentionDataFilm::Kh(const volScalarField &h, volScalarField &Kh)
    {
        Kh = Klc(h) + Klf(h);
        return Kh;
    } 

    volScalarField& retentionDataFilm::Cv(const volScalarField &h, volScalarField &Cv)
    {
        Cv = 0.5* ((1 + sign(h)) * shp_cv_res + (1 - sign(h))*((shp_ret_th_s - shp_ret_th_r) * CvCapillary(h) + shp_ret_th_r * CvFilm(h)));
        return Cv;
    }

    volScalarField& retentionDataFilm::Theta(const volScalarField &h, volScalarField &Theta)
    {
        Theta = 0.5*((1 + sign(h))*shp_ret_th_s + (1 - sign(h))*((shp_ret_th_s - shp_ret_th_r) * Scap(h) + shp_ret_th_r * Sad(h)));
        return Theta;
    }

    Foam::tmp<Foam::volScalarField> retentionDataFilm::Gamma(const volScalarField &h)
    {
        errInfo();
        return *static_cast<volScalarField *>(0);
    }

    Foam::tmp<Foam::volScalarField> retentionDataFilm::Scap(const volScalarField &h)
    {
        const volScalarField gamma_h_0 = Gamma(shp_h_0);
        return 0.5*((1 + sign(h)) + (1 - sign(h))*(0.5*((1 + sign(h-shp_h_0))*((Gamma(h) - gamma_h_0)/(1 - gamma_h_0)))));
    }

    Foam::tmp<Foam::volScalarField> retentionDataFilm::Sad(const volScalarField &h)
    {        
        volScalarField x_a = shp_film_h_a;
        x_a.dimensions().reset(dimensionSet(0, 0, 0, 0, 0, 0, 0));
        x_a = Foam::log10(Foam::mag(x_a));
        volScalarField x_0 = shp_h_0;   
        x_0.dimensions().reset(dimensionSet(0, 0, 0, 0, 0, 0, 0));
        x_0 = Foam::log10(Foam::mag(x_0));

        volScalarField h_clipped = Soil::Math::clipField(h);
        volScalarField x = h_clipped;
        x.dimensions().reset(dimensionSet(0, 0, 0, 0, 0, 0, 0));
        x = Foam::log10(Foam::mag(x)); //TODO: check if this is correct can h be positive in case of film flow i.e. log(|h|)?

        return Foam::max(tmp_zero, Foam::min(tmp_one, ((1 + (x - x_a + b*Foam::log(1+Foam::exp((x_a-x)/b)))/(x_a - x_0)))));
    }

    Foam::tmp<Foam::volScalarField> retentionDataFilm::CvCapillary(const volScalarField &h)
    {
        errInfo();
        return *static_cast<volScalarField *>(0);
    }

    
    Foam::tmp<Foam::volScalarField> retentionDataFilm::ThetaCapilary(const volScalarField &h)
    {
        volScalarField ret(h);
        ret.dimensions().reset(dimensionSet(0, 0, 0, 0, 0, 0, 0));
        volScalarField Cv_tmp = ret;
        ret = Theta(h, ret);
        return tmp<volScalarField>(new volScalarField(ret));
    }

    Foam::tmp<Foam::volScalarField> retentionDataFilm::CvFilm(const volScalarField &h)
    {
        volScalarField x_a = shp_film_h_a;
        x_a.dimensions().reset(dimensionSet(0, 0, 0, 0, 0, 0, 0));
        x_a = Foam::log10(Foam::mag(x_a));
        volScalarField x_0 = shp_h_0;   
        x_0.dimensions().reset(dimensionSet(0, 0, 0, 0, 0, 0, 0));
        x_0 = Foam::log10(Foam::mag(x_0));
        
		volScalarField h_clipped = Soil::Math::clipField(h);
        volScalarField x = h_clipped;
        x.dimensions().reset(dimensionSet(0, 0, 0, 0, 0, 0, 0));
        x = Foam::log10(Foam::mag(x)); //TODO: check if this is correct can h be positive in case of film flow  i.e. log(|h|)?
        volScalarField exp_fun = Foam::exp((x_a - x)/b);
        
        return 0.5 * (1 - Foam::sign(h - shp_film_h_a)) * ((1-exp_fun/(1+exp_fun))/(h_clipped*Foam::log(10.0)*(x_a - x_0))); 
        
    }

    Foam::tmp<Foam::volScalarField> retentionDataFilm::Klc(const volScalarField &h)
    {
        errInfo();
        return *static_cast<volScalarField *>(0);
    }

    Foam::tmp<Foam::volScalarField> retentionDataFilm::Klf(const volScalarField &h)
    {
        //TODO: check if this is correct, what about if h is positive? It should be zero for saturated flow?
        volScalarField h_clipped = Soil::Math::clipField(h); //FIXME: Why this was needed SIGfpe powf64 calculation error
		return 0.5 * (1 - Foam::sign(h)) * (0.5 * shp_film_k_sat * ((1 - sign(h - shp_film_h_a)) * Foam::pow(Foam::mag(h_clipped/shp_film_h_a), shp_film_alpha) + (1 + sign(h - shp_film_h_a)))); //original formulation
    }

    void retentionDataFilm::write(void)
    {
        ::retentionData::write();
        shp_h_0.write();
        shp_film_alpha.write();
        shp_film_h_a.write();
        shp_film_k_sat.write();
        b.write();
    };

    //c*θm*ha−1.5 c=1.35×10−8 m5/2 s−1 https://doi.org/10.5194/hess-27-4579-2023   eq. A6
    tmp<volScalarField> retentionDataFilm::setFilmKSat(fvMesh &mesh, Time &runTime, simConfig &config, const volScalarField &region, const volScalarField &shp_film_h_a) {
        scalar h_scalar = -1000.0;
        volScalarField tmp_field(IOobject("tmp_field", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("tmp_field", dimensionSet(0, 1, 0, 0, 0, 0, 0), h_scalar));
        volScalarField res_field(IOobject("res_field", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("res_field", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0));
        Theta(tmp_field, res_field);
        volScalarField tmp_shp_film_k_sat = 1.35e-8 * res_field * Foam::pow(Foam::mag(shp_film_h_a), -1.5);
        tmp_shp_film_k_sat.dimensions().reset(dimensionSet(0, 1, -1, 0, 0, 0, 0));

        std::vector<scalar> test_theta(no_of_regions, 0.0);
        std::vector<scalar> test_k_sat(no_of_regions, 0.0);
        forAll(region, index)
        {
            test_theta[static_cast<int>(region[index])] = res_field[index];
            test_k_sat[static_cast<int>(region[index])] = tmp_shp_film_k_sat[index];
        }
        for(int i=0; i<no_of_regions; i++)
        {
            Info << "Region: " << i << " Theta: " << test_theta[i] << " K_sat_fil: " << test_k_sat[i] <<" h_a: "<<min(shp_film_h_a)<< nl;
        }
        return tmp<volScalarField>(new volScalarField(tmp_shp_film_k_sat));
    }; 

    void retentionDataFilm::debugShpExtra(const Foam::fvMesh &mesh, const Foam::Time &runTime, retentionDataFilm* child) {
        const bool debugShp = config.algorithm.getShpDebug();
        Info<<"Initialising the inverse retention function interpolator data"<<endl;
        //Initialize inverse retention function interpolator data    
        label n_cells = 1e6;
        if (debugShp) {
            n_cells = 1e4;
        }
        const scalar h_min = -10;
        const scalar h_max = 12;
        const scalar h_step = (h_max - h_min) / n_cells; 
    
        
        std::vector<double> h_debug;
        std::vector<double> theta_debug;
        std::vector<double> se_cap_debug;
        std::vector<double> se_ad_debug;
        std::vector<double> kh_debug;
        std::vector<double> klc_debug;
        std::vector<double> klf_debug;
        std::vector<double> cv_debug;
        std::vector<double> cv_film_debug;
        std::vector<double> cv_cap_debug;
    
        // Create a standalone volScalarField with a default value of 0
        volScalarField h_debug_field (
            IOobject(
                "h_debug_field",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimensionSet(0, 1, 0, 0, 0, 0, 0), 0.0)
        );
    
        const label n_cells_mesh = mesh.nCells();
        Info<<"n_cells_mesh: "<<n_cells_mesh<<endl;
        // Initialize the field with custom values
        label i_vector=0; 
        label no_progess=20;
        label l_last = -1;
        for (label round = 0; round <= n_cells/n_cells_mesh; ++round)
        {
            label l = ceil(round/((n_cells/n_cells_mesh)/no_progess));
            if (l > l_last) {
                l_last = l;
                cout<<".";
            }
        }
        cout<<nl;
        l_last = -1;
        for (label round = 0; round <= n_cells/n_cells_mesh; ++round)
        {
            label l = ceil(round/((n_cells/n_cells_mesh)/no_progess));
            if (l > l_last) {
                l_last = l;
                cout<<"=";
            }
            
    
            for (label celli = 0; celli < n_cells_mesh; ++celli)
            {
                h_debug_field[celli] = -pow(10, h_min + h_step * (round * n_cells_mesh + celli));
            }
            Info<<"a"<<nl;
            volScalarField theta_debug_field = retentionDataInterface::Theta(h_debug_field);
            Info<<"b"<<nl;
            volScalarField inv_ret_se_cap_field = Scap(h_debug_field);
            Info<<"d"<<nl;
            volScalarField inv_ret_se_ad_field = Sad(h_debug_field);
            Info<<"e"<<nl;

            
            volScalarField inv_ret_kh_field = retentionDataInterface::Kh(h_debug_field);
            Info<<"f"<<nl;
            volScalarField inv_ret_klc_field = child->Klc(h_debug_field);
            Info<<"g"<<nl;
            volScalarField inv_ret_klf_field = Klf(h_debug_field);
            Info<<"h"<<nl;
            volScalarField inv_ret_cv_field = retentionDataInterface::Cv(h_debug_field);
            Info<<"i"<<nl;
            volScalarField inv_ret_cv_film_field = CvFilm(h_debug_field);
            Info<<"j"<<nl;
            volScalarField inv_ret_cv_cap_field = child->CvCapillary(h_debug_field);
            Info<<"k"<<nl;
            // Info<<"h_debug_field: "<<h_debug_field<<endl;
    
            for (label celli = 0; celli < n_cells_mesh; ++celli)
            {
                h_debug.push_back(h_debug_field[celli]);
                theta_debug.push_back(theta_debug_field[celli]);
                if (debugShp) {
                    se_cap_debug.push_back(inv_ret_se_cap_field[celli]);
                    se_ad_debug.push_back(inv_ret_se_ad_field[celli]);
                    kh_debug.push_back(inv_ret_kh_field[celli]);
                    klc_debug.push_back(inv_ret_klc_field[celli]);
                    klf_debug.push_back(inv_ret_klf_field[celli]);
                    cv_debug.push_back(inv_ret_cv_field[celli]);
                    cv_film_debug.push_back(inv_ret_cv_film_field[celli]);
                    cv_cap_debug.push_back(inv_ret_cv_cap_field[celli]);
                }
                i_vector++;
                if (i_vector >= n_cells) {
                    break;
                }
            }
        }
        
        cout<<nl;
    
    
        if (debugShp) {
            int vsize = h_debug.size();
            string csv_fname="debug_shp_film.csv";
            OFstream debugFile(csv_fname);
            debugFile<<"h;abs_h;theta;se_cap;se_film;kh;klc;klf;cv;cv_cap;cv_film"<<endl;
            for (int n=0; n<vsize; n++)
            {
                debugFile<<h_debug[n]<<";"<<abs(h_debug[n])<<";"<<theta_debug[n]<<";"<<se_cap_debug[n]<<";"<<se_ad_debug[n]<<";"<<kh_debug[n]<<";"<<klc_debug[n]<<";"<<klf_debug[n]<<";"<<cv_debug[n]<<";"<<cv_cap_debug[n]<<";"<<cv_film_debug[n]<<endl;
            }
            Info<<"Debug file created: "<<csv_fname<<endl;
        }
    }

    
} // namespace Soil::RetentionModels
