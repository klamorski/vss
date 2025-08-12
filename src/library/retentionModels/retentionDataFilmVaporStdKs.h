#ifndef RETENTIONDATAFILMVAPORSTDKS_H
#define RETENTIONDATAFILMVAPORSTDKS_H

#include "retentionDataFilmKs.h"
#include "retentionDataVaporStd.h"

namespace Soil::RetentionModels
{
    /**
     * @brief Base class for liquid capillary and film flow implementing van Genuchten retention model
     *
     */
    class retentionDataFilmVaporStdKs : public retentionDataFilmKs, public retentionDataVaporStd
    {
    public:
        retentionDataFilmVaporStdKs(IOdictionary &transportProperties, volScalarField &region, fvMesh &mesh, Time &runTime, simConfig &config);
        
        volScalarField& Cv(const volScalarField &h, volScalarField &Cv);    
        Foam::tmp<Foam::volScalarField> Cv_vapor_debug(const volScalarField &h);
        volScalarField& Theta(const volScalarField &h, volScalarField &Theta);
        volScalarField& Kh(const volScalarField &h, volScalarField &Kh); 

        void write();
    };

} // namespace Soil::RetentionModels

#endif // RETENTIONDATAFILMVAPORSTDKS_H
