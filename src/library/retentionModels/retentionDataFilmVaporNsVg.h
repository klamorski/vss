#ifndef RETENTIONDATAFILMVAPORNSVG_H
#define RETENTIONDATAFILMVAPORNSVG_H

#include "retentionDataFilmVg.h"
#include "retentionDataVaporNs.h"

namespace Soil::RetentionModels
{
    /**
     * @brief Base class for liquid capillary and film flow implementing van Genuchten retention model
     *
     */
    class retentionDataFilmVaporNsVg : public retentionDataFilmVg, public retentionDataVaporNs
    {
    public:
        retentionDataFilmVaporNsVg(IOdictionary &transportProperties, volScalarField &region, fvMesh &mesh, Time &runTime, simConfig &config);
        
        volScalarField& Kh(const volScalarField &h, volScalarField &Kh);        

        void write();
    };

} // namespace Soil::RetentionModels

#endif // RETENTIONDATAFILMVAPORNSVG_H
