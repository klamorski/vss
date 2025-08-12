#ifndef SOILMATH_H
#define SOILMATH_H

#include "fvCFD.H"

namespace Soil::Math
{
    namespace Private {
        double evaluate_polynomial(const double *poly, double const &z, int count);
        double erf_inv_imp(const double &p, const double &q);
    }
    
    double erfc_inv(double z);
    double erf_inv(double z);

    scalar getFieldValueAtHeight(const volScalarField &field, const fvMesh &mesh, scalar height);
    volScalarField erfc_inv(const volScalarField &z);
    volScalarField erf_inv(const volScalarField &z);    
    volScalarField &clipField(volScalarField &field);
    tmp<volScalarField> clipField(const volScalarField &field);
    volScalarField &setZeroGradientBoundaryCondition(volScalarField &field);
    tmp<volScalarField> setZeroGradientBoundaryCondition(const volScalarField &field);
    tmp<volScalarField> f1(const volScalarField &h, dimensionedScalar a);
    tmp<volScalarField> f2(const volScalarField &h, dimensionedScalar a);

}

#endif // SOILMATH_H
