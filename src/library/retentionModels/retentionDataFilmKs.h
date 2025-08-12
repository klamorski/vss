#ifndef RETENTIONDATAFILMKS_H
#define RETENTIONDATAFILMKS_H

#include "retentionDataFilm.h"
#include "retentionDataKs.h"

namespace Soil::RetentionModels
{
    /**
     * @brief Base class for liquid capillary and film flow implementing Kosugi retention model
     *
     */
    class retentionDataFilmKs : public retentionDataFilm, public retentionDataKs
    {
    public:
        retentionDataFilmKs(IOdictionary &transportProperties, volScalarField &region, fvMesh &mesh, Time &runTime, simConfig &config);
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

#endif // RETENTIONDATAFILMKS_H
