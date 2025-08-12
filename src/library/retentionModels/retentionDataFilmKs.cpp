#include "retentionDataFilmKs.h"
#include "soilMath.h"


// uncomment to disable assert()
// #define NDEBUG
#include <cassert>

using namespace Soil::RetentionModels;
using namespace Soil::Math;

namespace Soil::RetentionModels
{

    retentionDataFilmKs::retentionDataFilmKs(IOdictionary &transportProperties, volScalarField &region, fvMesh &mesh, Time &runTime, simConfig &config)
        : retentionData(transportProperties, region, mesh, runTime, config),
          retentionDataFilm(transportProperties, region, mesh, runTime, config),
          retentionDataKs(transportProperties, region, mesh, runTime, config)
    {
        b = 0.1+0.07*getShpRetKsSigma()*(1-Foam::exp(-Foam::pow(getShpRetThR()/(getShpRetThS()-getShpRetThR()),2))); //[Iden&Durner 2014] eq. 4

        if (average(shp_film_alpha).value() == 0.0) {
            shp_film_h_a = getShpRetKsHM(); //Peters 3013, [22]
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
    }
    
    volScalarField& retentionDataFilmKs::Kh(const volScalarField &h, volScalarField &Kh) {
        return ::retentionDataFilm::Kh(h, Kh);
    }

    volScalarField& retentionDataFilmKs::Cv(const volScalarField &h, volScalarField &Cv) {
        return ::retentionDataFilm::Cv(h, Cv);
    }

    volScalarField& retentionDataFilmKs::Theta(const volScalarField &h, volScalarField &Theta) {
        return ::retentionDataFilm::Theta(h, Theta);
    }

    Foam::tmp<Foam::volScalarField> retentionDataFilmKs::Gamma(const volScalarField &h)
    {
        return 0.5*Foam::erfc(Foam::log(Foam::mag(h/getShpRetKsHM()))/(Foam::sqrt(2.0)*getShpRetKsSigma()));
    }

    Foam::tmp<Foam::volScalarField> retentionDataFilmKs::CvCapillary(const volScalarField &h)
    {
		volScalarField h_clipped = Soil::Math::clipField(h);
        return 0.5*(1 - sign(h)) * (0.5*(1 + sign(h-shp_h_0))*Foam::exp(-Foam::pow((Foam::log(Foam::mag(h/getShpRetKsHM()))/Foam::sqrt(2.0)*getShpRetKsSigma()),2))/((Gamma(getShpH0())-1)*Foam::sqrt(2.0*constant::mathematical::pi)*getShpRetKsSigma()*h_clipped));
    }

    Foam::tmp<Foam::volScalarField> retentionDataFilmKs::Klc(const volScalarField &h)
    {
        volScalarField er_g_h0 = Foam::erf(erfc_inv(2*Gamma(getShpH0())+getShpRetKsSigma()/Foam::sqrt(2.0)));
        volScalarField er_g_h = Foam::erf(erfc_inv(2*Gamma(h)+getShpRetKsSigma()/Foam::sqrt(2.0)));
        const volScalarField base = shp_mual_k_sat*Foam::pow((er_g_h0-er_g_h)/(1+er_g_h0), 2);
        volScalarField h_clamped = h; 
        h_clamped.clamp_min(0.99*getShpH0Value());
        return Foam::pow(Scap(h_clamped), shp_mual_se_l)*base;
    }

    void retentionDataFilmKs::write(void)
    {
        ::retentionDataFilm::write();
        ::retentionDataKs::write();
    }

} // namespace Soil::RetentionModels
