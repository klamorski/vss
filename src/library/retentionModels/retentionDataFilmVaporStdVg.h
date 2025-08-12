#ifndef RETENTIONDATAFILMVAPORSTDVG_H
#define RETENTIONDATAFILMVAPORSTDVG_H

#include "retentionDataFilmVg.h"
#include "retentionDataVaporStd.h"

namespace Soil::RetentionModels
{
    /**
     * @brief Base class for liquid capillary and film flow implementing van Genuchten retention model
     *
     */
    class retentionDataFilmVaporStdVg : public retentionDataFilmVg, public retentionDataVaporStd
    {
    public:
        retentionDataFilmVaporStdVg(IOdictionary &transportProperties, volScalarField &region, fvMesh &mesh, Time &runTime, simConfig &config);
        volScalarField& Cv(const volScalarField &h, volScalarField &Cv);    
        Foam::tmp<Foam::volScalarField> Cv_vapor_debug(const volScalarField &h);
        volScalarField& Theta(const volScalarField &h, volScalarField &Theta);
        volScalarField& Kh(const volScalarField &h, volScalarField &Kh); 

        void write();
    };

} // namespace Soil::RetentionModels

#endif // RETENTIONDATAFILMVAPORSTDVG_H
