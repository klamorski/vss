#include "retentionDataKs.h"
#include "mathematicalConstants.H"
#include "soilMath.h"

using namespace Soil::RetentionModels;

namespace Soil::RetentionModels
{

    retentionDataKs::retentionDataKs(IOdictionary &transportProperties, volScalarField &region, fvMesh &mesh, Time &runTime, simConfig &config)
        : retentionData(transportProperties, region, mesh, runTime, config),
          shp_ret_ks_sigma(IOobject("shp_ret_ks_sigma", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("shp_ret_ks_sigma", dimensionSet(0, 0, 0, 0, 0, 0, 0), 1)),
          shp_ret_ks_h_m(IOobject("shp_ret_ks_h_m", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("shp_ret_ks_h_m", dimensionSet(0, 1, 0, 0, 0, 0, 0), 1))
    {
        if (no_of_regions > 0)
        {
            char buffer[1000];

            shp_ret_ks_sigma_vector.reserve(no_of_regions);
            for (int i = 0; i < no_of_regions; i++)
            {
                sprintf(buffer, "shp_ret_ks_sigma_%d", i);
                std::string key_name(buffer);
                shp_ret_ks_sigma_vector.push_back(dimensionedScalar(key_name.c_str(), dimensionSet(0, 0, 0, 0, 0, 0, 0), transportProperties));
            }

            shp_ret_ks_h_m_vector.reserve(no_of_regions);
            for (int i = 0; i < no_of_regions; i++)
            {
                sprintf(buffer, "shp_ret_ks_h_m%d", i);
                std::string key_name(buffer);
                shp_ret_ks_h_m_vector.push_back(dimensionedScalar(key_name.c_str(), dimensionSet(0, 1, 0, 0, 0, 0, 0), transportProperties));
            }

            forAll(region, index)
            {
                shp_ret_ks_sigma[index] = shp_ret_ks_sigma_vector[static_cast<int>(region[index])].value();
                shp_ret_ks_h_m[index] = shp_ret_ks_h_m_vector[static_cast<int>(region[index])].value();
            }

            forAll(region.boundaryFieldRef(), ipatch)
            {
                forAll(region.boundaryFieldRef()[ipatch], index)
                {
                    shp_ret_ks_sigma.boundaryFieldRef()[ipatch][index] = shp_ret_ks_sigma_vector[static_cast<int>(region.boundaryField()[ipatch][index])].value();
                    shp_ret_ks_h_m.boundaryFieldRef()[ipatch][index] = shp_ret_ks_h_m_vector[static_cast<int>(region.boundaryField()[ipatch][index])].value();
                }
            }
        }
        else
        {
            shp_ret_ks_sigma = volScalarField(IOobject("shp_ret_ks_sigma", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
            shp_ret_ks_h_m = volScalarField(IOobject("shp_ret_ks_h_m", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
        }
    }

    Foam::volScalarField &retentionDataKs::Kh(const volScalarField &h, Foam::volScalarField &Kh)
    {
        const volScalarField se = 0.5 * Foam::erfc(Foam::log(h / shp_ret_ks_h_m) / Foam::sqrt(2.0) * shp_ret_ks_sigma);
        Kh = shp_mual_k_sat * 0.5 * ((1 + sign(h)) + (1 - sign(h)) * Foam::pow(se, shp_mual_se_l) * Foam::pow(0.5 * Foam::erfc((Foam::log(h / shp_ret_ks_h_m) / shp_ret_ks_sigma + shp_ret_ks_sigma) / Foam::sqrt(2.0)), 2));
        return Kh;
    }
    Foam::volScalarField &retentionDataKs::Cv(const volScalarField &h, Foam::volScalarField &Cv)
    {
		volScalarField h_clipped = Soil::Math::clipField(h);
        Cv = 0.5 * ((1 + sign(h)) * shp_cv_res + (1 - sign(h)) * ((shp_ret_th_s - shp_ret_th_r) * 0.5 * (-(2 * Foam::exp(-Foam::pow(Foam::log(h / shp_ret_ks_h_m) / (Foam::sqrt(2.0) * shp_ret_ks_sigma), 2))) / (Foam::sqrt(constant::mathematical::pi) * Foam::sqrt(2.0) * shp_ret_ks_sigma * h_clipped))));
        return Cv;
    }

    volScalarField& retentionDataKs::Theta(const volScalarField &h, volScalarField &Theta)
    {
        Theta = 0.5 * ((1 + sign(h)) * shp_ret_th_s + (1 - sign(h)) * (shp_ret_th_r + (shp_ret_th_s - shp_ret_th_r) * (0.5 * Foam::erfc(Foam::log(Foam::mag(h / shp_ret_ks_h_m)) / Foam::sqrt(2.0) * shp_ret_ks_sigma))));
        return Theta;
    }

    void retentionDataKs::write(void)
    {
        ::retentionData::write();
        shp_ret_ks_sigma.write();
        shp_ret_ks_h_m.write();
    };

} // namespace Soil::RetentionModels
