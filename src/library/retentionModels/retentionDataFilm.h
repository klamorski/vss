#ifndef RETENTIONDATAFILM_H
#define RETENTIONDATAFILM_H

#include "retentionData.h"
#include "../external/interpolate/src/libInterpolate/Interpolate.hpp"

namespace Soil::RetentionModels
{
    /**
     * @brief Base class for liquid capillary and film flow retention models
     *
     */
    class retentionDataFilm : virtual public retentionData
    {
    public:
        retentionDataFilm(IOdictionary &transportProperties, volScalarField &region, fvMesh &mesh, Time &runTime, simConfig &config);
        volScalarField& Kh(const volScalarField &h, volScalarField &Kh);
        volScalarField& Cv(const volScalarField &h, volScalarField &Cv);
        volScalarField& Theta(const volScalarField &h, volScalarField &Theta);

        inline volScalarField &getShpH0() { return shp_h_0; }
        inline const scalar &getShpH0Value() { return h_0; }
        inline volScalarField &getShpFilmHa() { return shp_film_h_a; }
        inline volScalarField &getShpFilmKSat() { return shp_film_k_sat; }
        inline volScalarField &getB() { return b; }

        Foam::tmp<Foam::volScalarField> setFilmKSat(fvMesh &mesh, Time &runTime, simConfig &config, const volScalarField &region, const volScalarField &shp_film_h_a);

        void write(void);

    protected:
        const scalar h_0 = -6.3e4;
        const scalar alpha = -1.5;

        volScalarField shp_h_0;
        volScalarField tmp_one;
        volScalarField tmp_zero;

        volScalarField shp_film_alpha;
        std::vector<dimensionedScalar> shp_film_alpha_vector;
        volScalarField shp_film_h_a;
        std::vector<dimensionedScalar> shp_film_h_a_vector;
        volScalarField shp_film_k_sat;
        std::vector<dimensionedScalar> shp_film_k_sat_vector;

        volScalarField b;

    public:
        Foam::tmp<Foam::volScalarField> Klf(const volScalarField &h);
        void debugShpExtra(const Foam::fvMesh &mesh, const Foam::Time &runTime, retentionDataFilm* child);
        virtual Foam::tmp<Foam::volScalarField> Gamma(const volScalarField &h);
        virtual Foam::tmp<Foam::volScalarField> CvCapillary(const volScalarField &h);
        virtual Foam::tmp<Foam::volScalarField> Klc(const volScalarField &h);
        Foam::tmp<Foam::volScalarField> Scap(const volScalarField &h);
        Foam::tmp<Foam::volScalarField> Sad(const volScalarField &h);
        Foam::tmp<Foam::volScalarField> CvFilm(const volScalarField &h);
        Foam::tmp<Foam::volScalarField> ThetaCapilary(const volScalarField &h);
    };

    
} // namespace Soil::RetentionModels

#endif // RETENTIONDATAFILM_H
