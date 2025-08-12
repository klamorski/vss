#include "retentionDataInterface.h"
#include "OFstream.H"

using namespace Soil::RetentionModels;

namespace Soil::RetentionModels
{

retentionDataInterface::retentionDataInterface(fvMesh &mesh, Time &runTime, simConfig &config):
    config(config)
{
}

// Foam::tmp<Foam::volScalarField> retentionDataInterface::h(volScalarField &Theta) {
//     volScalarField ret(Theta);
//     ret.dimensions().reset(dimensionSet(0, 1, 0, 0, 0, 0, 0));
//     forAll(ret, index)
//     {
//         ret[index] = inv_ret_interpolator_(Theta[index]);
//     }

//     forAll(ret.boundaryFieldRef(), ipatch)
//     {
//         forAll(ret.boundaryFieldRef()[ipatch], index)
//         {
//             ret.boundaryFieldRef()[ipatch][index] = inv_ret_interpolator_(Theta.boundaryFieldRef()[ipatch][index]);
//         }
//     }
// 	return tmp<volScalarField>(new volScalarField(ret));
// }


Foam::tmp<Foam::volScalarField> retentionDataInterface::Kh(const volScalarField &h)
{
    volScalarField ret(h);
    ret.dimensions().reset(dimensionSet(0, 1, -1, 0, 0, 0, 0));
    ret = Kh(h, ret);
	return tmp<volScalarField>(new volScalarField(ret));
}


volScalarField& retentionDataInterface::Kh(const volScalarField &h, volScalarField &Kh) {
    errInfo();
    return *static_cast<volScalarField *>(0);
}

Foam::tmp<Foam::volScalarField> retentionDataInterface::Cv(const volScalarField &h)
{
    volScalarField ret(h);
    ret.dimensions().reset(dimensionSet(0, -1, 0, 0, 0, 0, 0));
    ret = Cv(h, ret);
	return tmp<volScalarField>(new volScalarField(ret));
}

volScalarField& retentionDataInterface::Cv(const volScalarField &h, volScalarField &Cv) {
    errInfo();
    return *static_cast<volScalarField *>(0);
}

Foam::tmp<Foam::volScalarField> retentionDataInterface::Theta(const volScalarField &h)
{
    volScalarField ret(h);
    ret.dimensions().reset(dimensionSet(0, 0, 0, 0, 0, 0, 0));
    ret = Theta(h, ret);
	return tmp<volScalarField>(new volScalarField(ret));
}

volScalarField& retentionDataInterface::Theta(const volScalarField &h, volScalarField &Theta) {
    errInfo();
    return *static_cast<volScalarField *>(0);
}

volScalarField &retentionDataInterface::getDneTau()
{
    errInfo();
    return *static_cast<volScalarField *>(0);
}

volScalarField &retentionDataInterface::getSoilTemp()
{
    errInfo();
    return *static_cast<volScalarField *>(0);
}

volScalarField &retentionDataInterface::getShpCvRes()
{
    errInfo();
    return *static_cast<volScalarField *>(0);
}

volScalarField &retentionDataInterface::getShpMualSeL()
{
    errInfo();
    return *static_cast<volScalarField *>(0);
}

volScalarField &retentionDataInterface::getShpDneTau()
{
    errInfo();
    return *static_cast<volScalarField *>(0);
}

volScalarField &retentionDataInterface::getShpMualKSat()
{
    errInfo();
    return *static_cast<volScalarField *>(0);
}

volScalarField &retentionDataInterface::getShpRetThS()
{
    errInfo();
    return *static_cast<volScalarField *>(0);
}

volScalarField &retentionDataInterface::getShpRetThR()
{
    errInfo();
    return *static_cast<volScalarField *>(0);
}

volScalarField &retentionDataInterface::getShpRetTotalPor()
{
    errInfo();
    return *static_cast<volScalarField *>(0);
}

void retentionDataInterface::errInfo(void)
{
    Info << "Missing retention model declaration in dictionary \"transportProperties\". Exiting now." << endl;
    exit(0);
}


} // namespace Soil::RetentionModels
