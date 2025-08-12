#include "retentionDataHvk.h"
#include "soilMath.h"

using namespace Soil::RetentionModels;

namespace Soil::RetentionModels
{

	retentionDataHvk::retentionDataHvk(IOdictionary &transportProperties, volScalarField &region, fvMesh &mesh, Time &runTime, simConfig &config)
		: retentionData(transportProperties, region, mesh, runTime, config),
		  shp_ret_hvk_alpha(IOobject("shp_ret_hvk_alpha", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("shp_ret_hvk_alpha", dimensionSet(0, 0, 0, 0, 0, 0, 0), 1)),
		  shp_ret_hvk_beta(IOobject("shp_ret_hvk_beta", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("shp_ret_hvk_beta", dimensionSet(0, 0, 0, 0, 0, 0, 0), 1)),
		  shp_ret_hvk_a(IOobject("shp_ret_hvk_a", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("shp_ret_hvk_a", dimensionSet(0, 0, 0, 0, 0, 0, 0), 1)),
		  shp_ret_hvk_gamma(IOobject("shp_ret_hvk_gamma", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("shp_ret_hvk_gamma", dimensionSet(0, 0, 0, 0, 0, 0, 0), 1))
	{
		if (no_of_regions > 0)
		{
			char buffer[1000];

			shp_ret_hvk_alpha_vector.reserve(no_of_regions);
			for (int i = 0; i < no_of_regions; i++)
			{
				sprintf(buffer, "shp_ret_hvk_alpha_%d", i);
				std::string key_name(buffer);
				shp_ret_hvk_alpha_vector.push_back(dimensionedScalar(key_name.c_str(), dimensionSet(0, 0, 0, 0, 0, 0, 0), transportProperties));
			}

			shp_ret_hvk_beta_vector.reserve(no_of_regions);
			for (int i = 0; i < no_of_regions; i++)
			{
				sprintf(buffer, "shp_ret_hvk_beta_%d", i);
				std::string key_name(buffer);
				shp_ret_hvk_beta_vector.push_back(dimensionedScalar(key_name.c_str(), dimensionSet(0, 0, 0, 0, 0, 0, 0), transportProperties));
			}

			shp_ret_hvk_a_vector.reserve(no_of_regions);
			for (int i = 0; i < no_of_regions; i++)
			{
				sprintf(buffer, "shp_ret_hvk_a_%d", i);
				std::string key_name(buffer);
				shp_ret_hvk_a_vector.push_back(dimensionedScalar(key_name.c_str(), dimensionSet(0, 0, 0, 0, 0, 0, 0), transportProperties));
			}

			shp_ret_hvk_gamma_vector.reserve(no_of_regions);
			for (int i = 0; i < no_of_regions; i++)
			{
				sprintf(buffer, "shp_ret_hvk_gamma_%d", i);
				std::string key_name(buffer);
				shp_ret_hvk_gamma_vector.push_back(dimensionedScalar(key_name.c_str(), dimensionSet(0, 0, 0, 0, 0, 0, 0), transportProperties));
			}

			forAll(region, index)
			{
				shp_ret_hvk_alpha[index] = shp_ret_hvk_alpha_vector[static_cast<int>(region[index])].value();
				shp_ret_hvk_beta[index] = shp_ret_hvk_beta_vector[static_cast<int>(region[index])].value();
				shp_ret_hvk_a[index] = shp_ret_hvk_a_vector[static_cast<int>(region[index])].value();
				shp_ret_hvk_gamma[index] = shp_ret_hvk_gamma_vector[static_cast<int>(region[index])].value();
			}

			forAll(region.boundaryFieldRef(), ipatch)
			{
				forAll(region.boundaryFieldRef()[ipatch], index)
				{
					shp_ret_hvk_alpha.boundaryFieldRef()[ipatch][index] = shp_ret_hvk_alpha_vector[static_cast<int>(region.boundaryField()[ipatch][index])].value();
					shp_ret_hvk_beta.boundaryFieldRef()[ipatch][index] = shp_ret_hvk_beta_vector[static_cast<int>(region.boundaryField()[ipatch][index])].value();
					shp_ret_hvk_a.boundaryFieldRef()[ipatch][index] = shp_ret_hvk_a_vector[static_cast<int>(region.boundaryField()[ipatch][index])].value();
					shp_ret_hvk_gamma.boundaryFieldRef()[ipatch][index] = shp_ret_hvk_gamma_vector[static_cast<int>(region.boundaryField()[ipatch][index])].value();
				}
			}
		}
		else
		{
			shp_ret_hvk_alpha = volScalarField(IOobject("shp_ret_hvk_alpha", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
			shp_ret_hvk_beta = volScalarField(IOobject("shp_ret_hvk_beta", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
			shp_ret_hvk_a = volScalarField(IOobject("shp_ret_hvk_a", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
			shp_ret_hvk_gamma = volScalarField(IOobject("shp_ret_hvk_gamma", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
		}
	}

    volScalarField& retentionDataHvk::Kh(const volScalarField &h, volScalarField &Kh)	
	{
		Kh = shp_mual_k_sat * 0.5 * ((1 + sign(h)) + (1 - sign(h)) * shp_ret_hvk_a / (shp_ret_hvk_a + pow(mag(dimensionedScalar(dimensionSet(0, -1, 0, 0, 0, 0, 0), 100.0) * h), shp_ret_hvk_gamma)));
		return Kh;
	}

	volScalarField& retentionDataHvk::Cv(const volScalarField &h, volScalarField &Cv)
	{
		const dimensionedScalar one_on_meter = dimensionedScalar(dimensionSet(0, -1, 0, 0, 0, 0, 0), 1.0);
		const dimensionedScalar one_hundredth_on_meter = dimensionedScalar(dimensionSet(0, -1, 0, 0, 0, 0, 0), 100.0);
		
		volScalarField h_clipped = Soil::Math::clipField(h);

		Cv = 0.5 * ((1 + sign(h)) * shp_cv_res + (1 - sign(h)) * (100.0 * one_on_meter * (shp_ret_hvk_alpha * shp_ret_hvk_beta * (shp_ret_th_s - shp_ret_th_r) * pow(mag(one_hundredth_on_meter * h), shp_ret_hvk_beta - 1)) / (pow(shp_ret_hvk_alpha + pow(mag(one_hundredth_on_meter * h_clipped), shp_ret_hvk_beta), 2))));
		return Cv;
	}

	volScalarField& retentionDataHvk::Theta(const volScalarField &h, volScalarField &Theta)
	{
		Theta = 0.5 * ((1 + sign(h)) * shp_ret_th_s + (1 - sign(h)) * (shp_ret_th_r + shp_ret_hvk_alpha * (shp_ret_th_s - shp_ret_th_r) / (shp_ret_hvk_alpha + pow(mag(dimensionedScalar(dimensionSet(0, -1, 0, 0, 0, 0, 0), 100.0) * h), shp_ret_hvk_beta))));
		return Theta;
	}

	void retentionDataHvk::write(void)
	{
		::retentionData::write();
		shp_ret_hvk_alpha.write();
		shp_ret_hvk_beta.write();
		shp_ret_hvk_a.write();
		shp_ret_hvk_gamma.write();
	};
	// https://www.wolframalpha.com/input?i=derivative+of+a*c%2F%28a%2B%280.01*x%29%5Eb%29
	// derivative of a*c/(a+(0.01*x)^b)

} // namespace Soil::RetentionModels
