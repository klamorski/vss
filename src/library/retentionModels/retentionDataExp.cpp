#include "retentionDataExp.h"
#include "soilMath.h"

using namespace Soil::RetentionModels;

namespace Soil::RetentionModels
{

	retentionDataExp::retentionDataExp(IOdictionary &transportProperties, volScalarField &region, fvMesh &mesh, Time &runTime, simConfig &config)
		: retentionData(transportProperties, region, mesh, runTime, config),
		  shp_ret_exp_alpha(IOobject("shp_ret_exp_alpha", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("shp_ret_exp_alpha", dimensionSet(0, -1, 0, 0, 0, 0, 0), 1))
	{
		if (no_of_regions > 0)
		{
			char buffer[1000];

			shp_ret_exp_alpha_vector.reserve(no_of_regions);
			for (int i = 0; i < no_of_regions; i++)
			{
				sprintf(buffer, "shp_ret_exp_alpha_%d", i);
				std::string key_name(buffer);
				shp_ret_exp_alpha_vector.push_back(dimensionedScalar(key_name.c_str(), dimensionSet(0, -1, 0, 0, 0, 0, 0), transportProperties));
			}

			// initiate internal fields
			forAll(region, index)
			{
				shp_ret_exp_alpha[index] = shp_ret_exp_alpha_vector[static_cast<int>(region[index])].value();
			}

			// initiate boundary fields
			forAll(region.boundaryFieldRef(), ipatch)
			{
				// each face in the the patch
				forAll(region.boundaryFieldRef()[ipatch], index)
				{
					shp_ret_exp_alpha.boundaryFieldRef()[ipatch][index] = shp_ret_exp_alpha_vector[static_cast<int>(region.boundaryField()[ipatch][index])].value();
				}
			}
		}
		else
		{
			shp_ret_exp_alpha = volScalarField(IOobject("shp_ret_exp_alpha", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
		}
	}

	Foam::volScalarField &retentionDataExp::Kh(const volScalarField &h, Foam::volScalarField &Kh)
	{
		Kh = shp_mual_k_sat * 0.5 * ((1 + sign(h)) + (1 - sign(h)) * Foam::exp(shp_ret_exp_alpha * h));
		return Kh;
	}

	Foam::volScalarField &retentionDataExp::Cv(const volScalarField &h, Foam::volScalarField &Cv)
	{
		volScalarField h_clipped = Soil::Math::clipField(h);
		Cv = 0.5 * ((1 + sign(h)) * shp_cv_res + (1 - sign(h)) * (shp_ret_exp_alpha * (shp_ret_th_s - shp_ret_th_r) * exp(shp_ret_exp_alpha * h_clipped)));
		return Cv;
	}

	Foam::volScalarField &retentionDataExp::Theta(const volScalarField &h, Foam::volScalarField &Theta)
	{
		Theta = 0.5 * ((1 + sign(h)) * shp_ret_th_s + (1 - sign(h)) * (shp_ret_th_r + (shp_ret_th_s - shp_ret_th_r) * exp(shp_ret_exp_alpha * h)));
		return Theta;
	}

	void retentionDataExp::write(void)
	{
		::retentionData::write();
		shp_ret_exp_alpha.write();
	};
} // namespace Soil::RetentionModels
