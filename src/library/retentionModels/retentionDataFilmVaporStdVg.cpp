#include "retentionDataFilmVaporStdVg.h"
#include "soilMath.h"
#include "soilPhysics.h"

// uncomment to disable assert()
// #define NDEBUG
#include <cassert>
using namespace Soil::RetentionModels;
using namespace Soil::Physics;

namespace Soil::RetentionModels
{

    retentionDataFilmVaporStdVg::retentionDataFilmVaporStdVg(IOdictionary &transportProperties, volScalarField &region, fvMesh &mesh, Time &runTime, simConfig &config)
        :retentionData(transportProperties, region, mesh, runTime, config), 
        retentionDataVaporStd(transportProperties, region, mesh, runTime, config),
        retentionDataFilmVg(transportProperties, region, mesh, runTime, config)
    {
        if (config.algorithm.getShpDebug()) {
            retentionDataVaporStd::debugShpExtra(mesh, runTime, this);
        }
    }

    volScalarField &retentionDataFilmVaporStdVg::Cv(const volScalarField &h, volScalarField &Cv)
    {
        Info<<"retentionDataFilmVaporStdVg::Cv"<<nl;
        volScalarField Cv_tmp = Cv;
        Cv = CvFullVapor(retentionDataFilmVg::Cv(h, Cv_tmp), h);
        return Cv;
    }

    Foam::tmp<Foam::volScalarField> retentionDataFilmVaporStdVg::Cv_vapor_debug(const volScalarField &h) {
        volScalarField ret(h);
        ret.dimensions().reset(dimensionSet(0, -1, 0, 0, 0, 0, 0));
        volScalarField Cv_tmp = ret;
        ret = CvFullVapor(retentionDataFilmVg::Cv(h, Cv_tmp), h) - retentionDataFilmVg::Cv(h, Cv_tmp);
        return tmp<volScalarField>(new volScalarField(ret));
    }

    volScalarField& retentionDataFilmVaporStdVg::Theta(const volScalarField &h, volScalarField &Theta) {
        volScalarField ret = retentionDataFilmVg::Theta(h, Theta);
        Theta =  ret + ThetaVapor(h, ret);
        return Theta;
    }

    volScalarField &retentionDataFilmVaporStdVg::Kh(const volScalarField &h, volScalarField &Kh)
    {   
        volScalarField Kh_tmp = Kh;
        Kh = retentionDataFilmVg::Kh(h, Kh_tmp) + getShpVaporBeta()*Kv(h);
        return Kh;
    }

    void retentionDataFilmVaporStdVg::write(void)
    {
        ::retentionDataFilmVg::write();
        ::retentionDataVaporStd::write();
    }

} // namespace Soil::RetentionModels
