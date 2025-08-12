#ifndef RETENTIONDATAFILMVG_H
#define RETENTIONDATAFILMVG_H

#include "retentionDataFilm.h"
#include "retentionDataVg.h"

namespace Soil::RetentionModels
{
    /**
     * @brief Base class for liquid capillary and film flow implementing van Genuchten retention model
     *
     */
    class retentionDataFilmVg : public retentionDataFilm, public retentionDataVg
    {
    public:
        retentionDataFilmVg(IOdictionary &transportProperties, volScalarField &region, fvMesh &mesh, Time &runTime, simConfig &config);
        volScalarField& Kh(const volScalarField &h, volScalarField &Kh);
        volScalarField& Cv(const volScalarField &h, volScalarField &Cv);
        volScalarField& Theta(const volScalarField &h, volScalarField &Theta);
    
        void write(void);

    protected:
        Foam::tmp<Foam::volScalarField> Gamma(const volScalarField &h);
        Foam::tmp<Foam::volScalarField> CvCapillary(const volScalarField &h);
        Foam::tmp<Foam::volScalarField> Klc(const volScalarField &h);
    };

} // namespace Soil::RetentionModels

#endif // RETENTIONDATAFILMVG_H
