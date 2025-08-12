#ifndef SOILPHYSICS_H
#define SOILPHYSICS_H

#include "fvCFD.H"

namespace Soil::Physics::Constants
{
    
    const dimensionedScalar g(dimensionSet(0, 1, -2, 0, 0, 0, 0), 9.80665);      // gravitational accelleration [m/s^2] TODO: get from dictionary
    const dimensionedScalar R(dimensionSet(1, 2, -2, -1, -1, 0, 0), 8.31446261815324); // Universal gas constant [kg⋅m2⋅s−2⋅K−1⋅mol−1]
    const dimensionedScalar Rsv(dimensionSet(0, 2, -2, -1, 0, 0, 0), 4.6152e2); // specific gas constant of the water vapour R/Mwat [m^2/s^2K]
    const dimensionedScalar Rsa(dimensionSet(0, 2, -2, -1, 0, 0, 0), 2.87052874e2); // specific gas constant of the dry air [m^2/s^2K]
    const dimensionedScalar Mwat(dimensionSet(1, 0, 0, 0, -1, 0, 0), 0.018015); // molecular weight of water [m^2/s^2K]

} // namespace Soil::Physics::Constants

namespace Soil::Physics
{

     /**
     * @brief Calculate the soil vapor diffusion coefficient
     *
     * @param T temperature [K]
     * @param theta volumetric water content [m^3/m^3]
     * @param totPorosity soil total porosity [m^3/m^3]
     *
     * @return volScalarField vapor diffusivity [m^2/s]
     *
     * Calculated according to [Iden et al., 2021] equation no. 9.
     */
    volScalarField vaporSoilDiffusivity(const volScalarField &T, const volScalarField &theta, const volScalarField &totPorosity);
    scalar vaporSoilDiffusivity(const scalar &T, const scalar &theta, const scalar &totPorosity);

    /**
     * @brief Calculate the vapor diffusivity
     *
     * @param T temperature [K]
     * @return volScalarField vapor diffusivity [m^2/s]
     *
     * Calculated according to [Saito et al., 2006] equation no. 16.
     */
    volScalarField vaporDiffusivity(const volScalarField &T);
    scalar vaporDiffusivity(const scalar &T);

    /**
     * @brief Calculate the saturated vapor density
     *
     * @param T temperature [K]
     * @return volScalarField saturated vapor density [kg/m^3]
     *
     * Calculated according to [Saito et al., 2006] equation no. 17.
     */
    volScalarField saturatedVaporDensity(const volScalarField &T);
    scalar saturatedVaporDensity(const scalar &T);

    /**
     * @brief Calculate the relative humidity
     *
     * @param T temperature [K]
     * @param h volumetric water content [m^3/m^3]
     * @return volScalarField relative humidity [-]
     *
     * Calculated according to [Saito et al., 2006] equation no. 18.
     */
    volScalarField soilPoreRelativeHumidity(const volScalarField &T, const volScalarField &h);
    scalar soilPoreRelativeHumidity(const scalar &T, const scalar &h);

    /**
     * @brief Calculate the specific water density
     *
     * @param T temperature [K]
     * @return volScalarField specific water density [kg/m^3]
     *
     */
    volScalarField specificWaterDensity(const volScalarField &T);
    scalar specificWaterDensity(const scalar &T);

    /**
     * @brief Convert temperature from Celsius to Kelvin
     *
     * @param T temperature [C]
     * @return volScalarField temperature [K]
     *
     */
    inline volScalarField& CtoK(volScalarField &C)
    {
        C = C + dimensionedScalar("K", dimensionSet(0,0,0,1,0,0,0), 273.15);
        return C;
    };

    inline volScalarField CtoK(const volScalarField &C)
    {
        return C + dimensionedScalar("K", dimensionSet(0,0,0,1,0,0,0), 273.15);
    };
    
    inline scalar CtoK(const scalar &C)
    {
        return C + 273.15;
    };

    /**
     * @brief Convert temperature from Kelvin to Celsius
     *
     * @param T temperature [K]
     * @return volScalarField temperature [C]
     *
     */
    
    inline volScalarField& KtoC(volScalarField &K)
    {
        K = K - dimensionedScalar("K", dimensionSet(0,0,0,1,0,0,0), 273.15);
        return K;
    }

    inline volScalarField KtoC(const volScalarField &K)
    {
        return K - dimensionedScalar("K", dimensionSet(0,0,0,1,0,0,0), 273.15);
    }

    inline scalar KtoC(const scalar &K)
    {
        return K - 273.15;
    }

} // namespace Soil::Physics

#endif // SOILPHYSICS_H
