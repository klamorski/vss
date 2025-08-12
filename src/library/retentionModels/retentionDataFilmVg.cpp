#include "retentionDataFilmVg.h"
#include "soilMath.h"

// uncomment to disable assert()
// #define NDEBUG
#include <cassert>
using namespace Soil::RetentionModels;

namespace Soil::RetentionModels
{

    retentionDataFilmVg::retentionDataFilmVg(IOdictionary &transportProperties, volScalarField &region, fvMesh &mesh, Time &runTime, simConfig &config)
        : retentionData(transportProperties, region, mesh, runTime, config),
          retentionDataFilm(transportProperties, region, mesh, runTime, config),
          retentionDataVg(transportProperties, region, mesh, runTime, config)
    {
        b = 0.1 + 0.2 * (1 - Foam::exp(-Foam::pow(shp_ret_th_r / (shp_ret_th_s - shp_ret_th_r),2))) / Foam::pow(shp_ret_vg_n, 2); //[Iden&Durner 2014] eq. 5

        Info << "average(shp_film_h_a).value() = " << average(shp_film_h_a).value() << endl;
        if (average(shp_film_h_a).value() == 0.0) {
            shp_film_h_a = -Foam::pow((Foam::pow(2, 1/shp_ret_vg_m) - 1), 1/shp_ret_vg_n)/shp_ret_vg_alpha; //Peters 3013, [22] - full approximation
            //shp_film_h_a = -1.0/shp_ret_vg_alpha; //Peters 3013, [22] - simplified approximation
            Info<< "shp_film_h_a estimated by model: "<< average(shp_film_h_a).value() << endl;
        } else {
            Info<< "shp_film_h_a set by transportProperties: "<< average(shp_film_h_a).value() << endl;
        }
        
        assert((Foam::max(shp_film_h_a).value() < 0) && "shp_film_h_a must be negative");
        
        if (average(shp_film_k_sat).value() == 0.0 ) {
            shp_film_k_sat = setFilmKSat(mesh, runTime, config, region, shp_film_h_a);
            Info<< "shp_film_k_sat estimated by model: "<< average(shp_film_k_sat).value() << endl;
        } else {
            Info<< "shp_film_k_sat set by transportProperties: "<< average(shp_film_k_sat).value() << endl;
        }

        //debugShpExtra(mesh, runTime, this);
    }
    
    volScalarField& retentionDataFilmVg::Kh(const volScalarField &h, volScalarField &Kh) {
        return ::retentionDataFilm::Kh(h, Kh);
    }

    volScalarField& retentionDataFilmVg::Cv(const volScalarField &h, volScalarField &Cv) {
        return ::retentionDataFilm::Cv(h, Cv);
    }

    volScalarField& retentionDataFilmVg::Theta(const volScalarField &h, volScalarField &Theta) {
        return ::retentionDataFilm::Theta(h, Theta);
    }

    Foam::tmp<Foam::volScalarField> retentionDataFilmVg::Gamma(const volScalarField &h)
    {
        return 0.5 * ((1 + Foam::sign(h)) + (1 - Foam::sign(h)) * (Foam::pow((1 + Foam::pow(Foam::mag(shp_ret_vg_alpha*h), shp_ret_vg_n)), -shp_ret_vg_m)));
    }

    Foam::tmp<Foam::volScalarField> retentionDataFilmVg::CvCapillary(const volScalarField &h)
    {
		volScalarField h_clipped = Soil::Math::clipField(h);
        return 0.5 * (1 - Foam::sign(h)) * (0.5*(1 + sign(h-shp_h_0))*shp_ret_vg_m*shp_ret_vg_n*Foam::pow(Foam::mag(shp_ret_vg_alpha*h), shp_ret_vg_n)*Foam::pow(1+Foam::pow(Foam::mag(shp_ret_vg_alpha*h), shp_ret_vg_n), -shp_ret_vg_m-1)/(Foam::mag(h_clipped)*(1-Gamma(shp_h_0))));
    }

    Foam::tmp<Foam::volScalarField> retentionDataFilmVg::Klc(const volScalarField &h)
    {
        const volScalarField base = shp_mual_k_sat*Foam::pow(1-Foam::pow((1-Foam::pow(Gamma(h),1/shp_ret_vg_m))/(1-Foam::pow(Gamma(shp_h_0),1/shp_ret_vg_m)),shp_ret_vg_m),2);
        volScalarField h_clamped = h; 
        h_clamped.clamp_min(0.99*getShpH0Value());
        return Foam::pow(Scap(h_clamped), shp_mual_se_l)*base;
   }

    void retentionDataFilmVg::write(void)
    {
        ::retentionDataFilm::write();
        ::retentionDataVg::write();
    }

} // namespace Soil::RetentionModels
