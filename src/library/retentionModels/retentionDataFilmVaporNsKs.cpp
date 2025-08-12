#include "retentionDataFilmVaporNsKs.h"
#include "soilMath.h"

// uncomment to disable assert()
// #define NDEBUG
#include <cassert>
using namespace Soil::RetentionModels;

namespace Soil::RetentionModels
{

    retentionDataFilmVaporNsKs::retentionDataFilmVaporNsKs(IOdictionary &transportProperties, volScalarField &region, fvMesh &mesh, Time &runTime, simConfig &config)
        :retentionData(transportProperties, region, mesh, runTime, config), 
        retentionDataVaporNs(transportProperties, region, mesh, runTime, config),
        retentionDataFilmKs(transportProperties, region, mesh, runTime, config)
    {
        if (config.algorithm.getShpDebug()) {
            retentionDataVaporNs::debugShpExtra(mesh, runTime, this);
        }
    }

    volScalarField &retentionDataFilmVaporNsKs::Kh(const volScalarField &h, volScalarField &Kh)
    {   
        volScalarField Kh_tmp = Kh;
        Kh = retentionDataFilmKs::Kh(h, Kh_tmp) + getShpVaporBeta()*Kv(h);
        return Kh;
    }

    void retentionDataFilmVaporNsKs::write(void)
    {
        ::retentionDataFilmKs::write();
        ::retentionDataVaporNs::write();
    }

} // namespace Soil::RetentionModels
