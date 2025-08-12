#include "retentionDataFilmVaporNsVg.h"
#include "soilMath.h"

// uncomment to disable assert()
// #define NDEBUG
#include <cassert>
using namespace Soil::RetentionModels;

namespace Soil::RetentionModels
{

    retentionDataFilmVaporNsVg::retentionDataFilmVaporNsVg(IOdictionary &transportProperties, volScalarField &region, fvMesh &mesh, Time &runTime, simConfig &config)
        :retentionData(transportProperties, region, mesh, runTime, config), 
        retentionDataVaporNs(transportProperties, region, mesh, runTime, config),
        retentionDataFilmVg(transportProperties, region, mesh, runTime, config)
    {
        if (config.algorithm.getShpDebug()) {
            retentionDataVaporNs::debugShpExtra(mesh, runTime, this);
        }
    }

    volScalarField &retentionDataFilmVaporNsVg::Kh(const volScalarField &h, volScalarField &Kh)
    {   
        volScalarField Kh_tmp = Kh;
        Kh = retentionDataFilmVg::Kh(h, Kh_tmp) + getShpVaporBeta()*Kv(h);
        return Kh;
    }

    void retentionDataFilmVaporNsVg::write(void)
    {
        ::retentionDataFilmVg::write();
        ::retentionDataVaporNs::write();
    }

} // namespace Soil::RetentionModels
