#ifndef RETENTIONDATAFILMVAPORNSKS_H
#define RETENTIONDATAFILMVAPORNSKS_H

#include "retentionDataFilmKs.h"
#include "retentionDataVaporNs.h"

namespace Soil::RetentionModels
{
    /**
     * @brief Base class for liquid capillary and film flow implementing van Genuchten retention model
     *
     */
    class retentionDataFilmVaporNsKs : public retentionDataFilmKs, public retentionDataVaporNs
    {
    public:
        retentionDataFilmVaporNsKs(IOdictionary &transportProperties, volScalarField &region, fvMesh &mesh, Time &runTime, simConfig &config);
        
        volScalarField& Kh(const volScalarField &h, volScalarField &Kh);        

        void write();
    };

} // namespace Soil::RetentionModels

#endif // RETENTIONDATAFILMVAPORNSKS_H
