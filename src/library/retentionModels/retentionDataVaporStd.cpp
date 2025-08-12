#include "retentionDataVaporStd.h"
#include "soilPhysics.h"

using namespace Soil::RetentionModels;
using namespace Soil::Physics;

namespace Soil::RetentionModels
{

    retentionDataVaporStd::retentionDataVaporStd(IOdictionary &transportProperties, volScalarField &region, fvMesh &mesh, Time &runTime, simConfig &config): 
        retentionData(transportProperties, region, mesh, runTime, config),
        retentionDataVaporNs(transportProperties, region, mesh, runTime, config)
    {
       
    }

    volScalarField &retentionDataVaporStd::Cv(const volScalarField &h, volScalarField &Cv)
    {
        Info<<"retentionDataVaporStd::Cv"<<nl;
        Cv = CvFullVapor(retentionDataVaporNs::Cv(h), h);
        return Cv;
    }

    
    Foam::tmp<Foam::volScalarField>  retentionDataVaporStd::ThetaVapor(const volScalarField &h, const volScalarField &Theta_cap) {
        const volScalarField &T = getSoilTemp();
        const volScalarField &tot_porosity = getShpRetTotalPor();
        const volScalarField rho_H = saturatedVaporDensity(T)*soilPoreRelativeHumidity(T, h)/specificWaterDensity(T);
        return tmp<volScalarField>(new volScalarField((getShpRetTotalPor() - Theta_cap)*rho_H));
    }

    Foam::tmp<Foam::volScalarField> retentionDataVaporStd::CvFullVapor(const volScalarField &Cv, const volScalarField &h) {
        const volScalarField &T = getSoilTemp();
        const volScalarField rho_H = saturatedVaporDensity(T)*soilPoreRelativeHumidity(T, h)/specificWaterDensity(T);
        const volScalarField ret = Cv*(1 - rho_H) + (getShpRetTotalPor() - Theta(h))*rho_H*Soil::Physics::Constants::g/(Soil::Physics::Constants::Rsv*T);
	    return tmp<volScalarField>(new volScalarField(ret));
    }

    Foam::tmp<Foam::volScalarField> retentionDataVaporStd::Cv_vapor_debug(const volScalarField &h) {
        volScalarField ret(h);
        ret.dimensions().reset(dimensionSet(0, -1, 0, 0, 0, 0, 0));
        ret = CvFullVapor(retentionDataVaporNs::Cv(h), h) - retentionDataVaporNs::Cv(h);
        return tmp<volScalarField>(new volScalarField(ret));
    }

    void retentionDataVaporStd::debugShpExtra(const Foam::fvMesh &mesh, const Foam::Time &runTime, retentionDataFilm* child) {
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
        std::vector<double> inv_ret_h_;
        std::vector<double> inv_ret_theta_;
        std::vector<double> se_cap_debug;
        std::vector<double> se_ad_debug;
        std::vector<double> theta_vapor_debug;
        std::vector<double> kh_debug;
        std::vector<double> klc_debug;
        std::vector<double> klf_debug;
        std::vector<double> kv_debug;
        std::vector<double> cv_debug;
        std::vector<double> cv_film_debug;
        std::vector<double> cv_cap_debug;
        std::vector<double> cv_vapor_debug;

        // Create a standalone volScalarField with a default value of 0
        volScalarField inv_ret_h_field (
            IOobject(
                "inv_ret_h_field",
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
                inv_ret_h_field[celli] = -pow(10, h_min + h_step * (round * n_cells_mesh + celli));
            }
            Info<<"a"<<nl;
            volScalarField inv_ret_theta_field = retentionDataInterface::Theta(inv_ret_h_field);
            Info<<"b"<<nl;
            volScalarField inv_ret_se_cap_field = child->Scap(inv_ret_h_field);
            Info<<"c"<<nl;
            volScalarField inv_ret_se_ad_field = child->Sad(inv_ret_h_field);
            Info<<"d"<<nl;
            volScalarField inv_ret_kh_field = Kh(inv_ret_h_field);
            Info<<"e"<<nl;
            volScalarField inv_ret_klc_field = child->Klc(inv_ret_h_field);
            Info<<"f"<<nl;
            volScalarField inv_ret_klf_field = child->Klf(inv_ret_h_field);
            Info<<"g"<<nl;
            volScalarField inv_ret_kv_field = Kv(inv_ret_h_field);
            Info<<"h"<<nl;
            volScalarField inv_ret_cv_field = retentionDataInterface::Cv(inv_ret_h_field);
            Info<<"i"<<nl;
            volScalarField inv_ret_cv_film_field = child->CvFilm(inv_ret_h_field);
            Info<<"j"<<nl;
            volScalarField inv_ret_cv_cap_field = child->CvCapillary(inv_ret_h_field);
            Info<<"k"<<nl;
            volScalarField inv_ret_cv_vapor_field = Cv_vapor_debug(inv_ret_h_field);
            Info<<"l"<<nl;
            volScalarField inv_ret_theta_vapor_field = ThetaVapor(inv_ret_h_field, child->ThetaCapilary(inv_ret_h_field));

            for (label celli = 0; celli < n_cells_mesh; ++celli)
            {
                inv_ret_h_.push_back(inv_ret_h_field[celli]);
                inv_ret_theta_.push_back(inv_ret_theta_field[celli]);
                if (debugShp) {
                    se_cap_debug.push_back(inv_ret_se_cap_field[celli]);
                    se_ad_debug.push_back(inv_ret_se_ad_field[celli]);
                    kh_debug.push_back(inv_ret_kh_field[celli]);
                    klc_debug.push_back(inv_ret_klc_field[celli]);
                    klf_debug.push_back(inv_ret_klf_field[celli]);
                    kv_debug.push_back(inv_ret_kv_field[celli]);
                    cv_debug.push_back(inv_ret_cv_field[celli]);
                    cv_film_debug.push_back(inv_ret_cv_film_field[celli]);
                    cv_cap_debug.push_back(inv_ret_cv_cap_field[celli]);
                    cv_vapor_debug.push_back(inv_ret_cv_vapor_field[celli]);
                    theta_vapor_debug.push_back(inv_ret_theta_vapor_field[celli]);
                }
                i_vector++;
                if (i_vector >= n_cells) {
                    break;
                }
            }
        }
        
        cout<<nl;


        if (debugShp) {
            int vsize = inv_ret_h_.size();
            string csv_fname="debug_shp_vapor_std.csv";
            OFstream debugFile(csv_fname);
            debugFile<<"h;abs_h;theta;se_cap;se_film;theta_vapor;kh;klc;klf;kv;cv;cv_cap;cv_film;cv_vapor"<<endl;
            for (int n=0; n<vsize; n++)
            {
                debugFile<<inv_ret_h_[n]<<";"<<abs(inv_ret_h_[n])<<";"<<inv_ret_theta_[n]<<";"<<se_cap_debug[n]<<";"<<se_ad_debug[n]<<";"<<theta_vapor_debug[n]<<";"<<kh_debug[n]<<";"<<klc_debug[n]<<";"<<klf_debug[n]<<";"<<kv_debug[n]<<";"<<cv_debug[n]<<";"<<cv_cap_debug[n]<<";"<<cv_film_debug[n]<<";"<<cv_vapor_debug[n]<<endl;
            }
            Info<<"Debug file created: "<<csv_fname<<endl;
        }
    }

    void retentionDataVaporStd::write() {
        retentionDataVaporNs::write();
    }

} // namespace Soil::RetentionModels

