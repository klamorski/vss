#include "retentionDataFilmVaporStdKs.h"
#include "soilMath.h"

// uncomment to disable assert()
// #define NDEBUG
#include <cassert>
using namespace Soil::RetentionModels;

namespace Soil::RetentionModels
{

    retentionDataFilmVaporStdKs::retentionDataFilmVaporStdKs(IOdictionary &transportProperties, volScalarField &region, fvMesh &mesh, Time &runTime, simConfig &config)
        :retentionData(transportProperties, region, mesh, runTime, config), 
        retentionDataVaporStd(transportProperties, region, mesh, runTime, config),
        retentionDataFilmKs(transportProperties, region, mesh, runTime, config)
    {
        if (config.algorithm.getShpDebug()) {
            retentionDataVaporStd::debugShpExtra(mesh, runTime, this);
        }
    }

    volScalarField &retentionDataFilmVaporStdKs::Cv(const volScalarField &h, volScalarField &Cv)
    {
        Info<<"retentionDataFilmVaporStdKs::Cv"<<nl;
        volScalarField Cv_tmp = Cv;
        Cv = CvFullVapor(retentionDataFilmKs::Cv(h, Cv_tmp), h);
        return Cv;
    }

    volScalarField& retentionDataFilmVaporStdKs::Theta(const volScalarField &h, volScalarField &Theta) {
        volScalarField ret = retentionDataFilmKs::Theta(h, Theta);
        Theta =  ret + ThetaVapor(h, ret);
        return Theta;
    }

    Foam::tmp<Foam::volScalarField> retentionDataFilmVaporStdKs::Cv_vapor_debug(const volScalarField &h) {
        volScalarField ret(h);
        ret.dimensions().reset(dimensionSet(0, -1, 0, 0, 0, 0, 0));
        volScalarField Cv_tmp = ret;
        ret = CvFullVapor(retentionDataFilmKs::Cv(h, Cv_tmp), h) - retentionDataFilmKs::Cv(h, Cv_tmp);
        return tmp<volScalarField>(new volScalarField(ret));
    }

    volScalarField &retentionDataFilmVaporStdKs::Kh(const volScalarField &h, volScalarField &Kh)
    {   
        volScalarField Kh_tmp = Kh;
        Kh = retentionDataFilmKs::Kh(h, Kh_tmp) + getShpVaporBeta()*Kv(h);
        return Kh;
    }

    void retentionDataFilmVaporStdKs::write(void)
    {
        ::retentionDataFilmKs::write();
        ::retentionDataVaporStd::write();
    }

} // namespace Soil::RetentionModels
