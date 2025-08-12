#include "soilPhysics.h"
#include "soilGlobal.h"
#include "dimensionedScalar.H"

namespace Soil::Physics
{
    
    volScalarField vaporSoilDiffusivity(const volScalarField &T, const volScalarField &theta, const volScalarField &totPorosity)
    {
        return vaporDiffusivity(T) * Foam::pow(Foam::max(dimensionedScalar("tmp_zero", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0), totPorosity - theta), 2.5) / totPorosity;
    }
    
    scalar vaporSoilDiffusivity(const scalar &T, const scalar &theta, const scalar &totPorosity)
    {
        return vaporDiffusivity(T) * Foam::pow(Foam::max(0.0, totPorosity - theta), 2.5) / totPorosity;
    }

    // scalar vaporDiffusivity(const scalar &T)
    // {
    //     return 2.841407e-10 * Foam::pow(T, 2);
    // }

    // volScalarField vaporDiffusivity(const volScalarField &T)
    // {
    //     // Check if the units of T are temperature (K)
    //     Soil::Global::checkUnits(T, dimensionSet(0, 0, 0, 1, 0, 0, 0), "T");
    //     volScalarField ret(2.841407e-10 * Foam::pow(T, 2));
    //     ret.dimensions().reset(dimensionSet(0, 2, -1, 0, 0, 0, 0));
    //     return ret;
    // }
    
    scalar vaporDiffusivity(const scalar &T)
    {
        return 2.12e-5 * Foam::pow(T/273.15, 2);
    }

    volScalarField vaporDiffusivity(const volScalarField &T)
    {
        // Check if the units of T are temperature (K)
        Soil::Global::checkUnits(T, dimensionSet(0, 0, 0, 1, 0, 0, 0), "T");
        volScalarField ret(2.12e-5 * Foam::pow(T/273.15, 2));
        ret.dimensions().reset(dimensionSet(0, 2, -1, 0, 0, 0, 0));
        return ret;
    }

    scalar saturatedVaporDensity(const scalar &T)
    {
        return Foam::exp(31.3716 - 6014.79 / T - 7.92495e-3 * T) / (1000.0 * T);
    }

    volScalarField saturatedVaporDensity(const volScalarField &T)
    {
        // Check if the units of T are temperature (K)
        Soil::Global::checkUnits(T, dimensionSet(0, 0, 0, 1, 0, 0, 0), "T");
        volScalarField Tnu = T;
        Tnu.dimensions().reset(dimensionSet(0, 0, 0, 0, 0, 0, 0));
        volScalarField ret(Foam::exp(31.3716 - 6014.79 / Tnu - 7.92495e-3 * Tnu) / (1000.0 * Tnu));
        ret.dimensions().reset(dimensionSet(1, -3, 0, 0, 0, 0, 0));
        return ret;
    }

    volScalarField soilPoreRelativeHumidity(const volScalarField &T, const volScalarField &h)
    {
        // Check if the units of T are temperature (K)
        Soil::Global::checkUnits(T, dimensionSet(0, 0, 0, 1, 0, 0, 0), "T");
        // Check if the units of h are pressure head (m)
        Soil::Global::checkUnits(h, dimensionSet(0, 1, 0, 0, 0, 0, 0), "h");
        return Foam::exp((h * Constants::g) / (Constants::Rsv * T)); //FIXME: limit to 1 in case of full saturation
    }
    
    scalar soilPoreRelativeHumidity(const scalar &T, const scalar &h)
    {
        return Foam::exp((h * Constants::g.value()) / (Constants::Rsv.value() * T));
    }

    volScalarField specificWaterDensity(const volScalarField &T)
    {
        // Check if the units of T are temperature (K)
        Soil::Global::checkUnits(T, dimensionSet(0, 0, 0, 1, 0, 0, 0), "T");

        const fvMesh &mesh = T.mesh();
        tmp<volScalarField> tmpField(
            new volScalarField(IOobject("specificWaterDensity", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
                            mesh,
                            dimensionedScalar("specificWaterDensity", dimensionSet(1, -3, 0, 0, 0, 0, 0), 1000)));
        return tmpField; // TODO: Implement T dependence by some equation
    }

    scalar specificWaterDensity(const scalar &T)
    {
        return 1000.0; // TODO: Implement T dependence by some equation
    }
} // namespace Soil::Physics


