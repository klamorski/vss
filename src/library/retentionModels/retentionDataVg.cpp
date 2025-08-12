#include "retentionDataVg.h"
#include "soilMath.h"

using namespace Soil::RetentionModels;

namespace Soil::RetentionModels
{

	retentionDataVg::retentionDataVg(IOdictionary &transportProperties, volScalarField &region, fvMesh &mesh, Time &runTime, simConfig &config)
		: retentionData(transportProperties, region, mesh, runTime, config),
		  shp_ret_vg_alpha(IOobject("shp_ret_vg_alpha", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("shp_ret_vg_alpha", dimensionSet(0, -1, 0, 0, 0, 0, 0), 1)),
		  shp_ret_vg_n(IOobject("shp_ret_vg_n", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("shp_ret_vg_n", dimensionSet(0, 0, 0, 0, 0, 0, 0), 1)),
		  shp_ret_vg_m(IOobject("shp_ret_vg_m", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("shp_ret_vg_m", dimensionSet(0, 0, 0, 0, 0, 0, 0), 1))
	{
		if (no_of_regions > 0)
		{
			char buffer[1000];

			shp_ret_vg_alpha_vector.reserve(no_of_regions);
			for (int i = 0; i < no_of_regions; i++)
			{
				sprintf(buffer, "shp_ret_vg_alpha_%d", i);
				std::string key_name(buffer);
				shp_ret_vg_alpha_vector.push_back(dimensionedScalar(key_name.c_str(), dimensionSet(0, -1, 0, 0, 0, 0, 0), transportProperties));
			}

			shp_ret_vg_n_vector.reserve(no_of_regions);
			for (int i = 0; i < no_of_regions; i++)
			{
				sprintf(buffer, "shp_ret_vg_n_%d", i);
				std::string key_name(buffer);
				shp_ret_vg_n_vector.push_back(dimensionedScalar(key_name.c_str(), dimensionSet(0, 0, 0, 0, 0, 0, 0), transportProperties));
			}

			shp_ret_vg_m_vector.reserve(no_of_regions);
			for (int i = 0; i < no_of_regions; i++)
			{
				sprintf(buffer, "shp_ret_vg_m_%d", i);
				std::string key_name(buffer);
				shp_ret_vg_m_vector.push_back(dimensionedScalar(key_name.c_str(), dimensionSet(0, 0, 0, 0, 0, 0, 0), 1.0 - 1.0 / shp_ret_vg_n_vector[i].value(), transportProperties));
			}

			//TODO: check if all cells were initialized - otherwise regions number is more likely wrong in the dictionary
			forAll(region, index)
			{
				shp_ret_vg_alpha[index] = shp_ret_vg_alpha_vector[static_cast<int>(region[index])].value();
				shp_ret_vg_n[index] = shp_ret_vg_n_vector[static_cast<int>(region[index])].value();
				shp_ret_vg_m[index] = shp_ret_vg_m_vector[static_cast<int>(region[index])].value();
			}

			forAll(region.boundaryFieldRef(), ipatch)
			{
				forAll(region.boundaryFieldRef()[ipatch], index)
				{
					shp_ret_vg_alpha.boundaryFieldRef()[ipatch][index] = shp_ret_vg_alpha_vector[static_cast<int>(region.boundaryField()[ipatch][index])].value();
					shp_ret_vg_n.boundaryFieldRef()[ipatch][index] = shp_ret_vg_n_vector[static_cast<int>(region.boundaryField()[ipatch][index])].value();
					shp_ret_vg_m.boundaryFieldRef()[ipatch][index] = shp_ret_vg_m_vector[static_cast<int>(region.boundaryField()[ipatch][index])].value();
				}
			}
		}
		else
		{
			shp_ret_vg_alpha = volScalarField(IOobject("shp_ret_vg_alpha", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
			shp_ret_vg_n = volScalarField(IOobject("shp_ret_vg_n", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
			shp_ret_vg_m = volScalarField(IOobject("shp_ret_vg_m", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
		}
	}

	volScalarField &retentionDataVg::Kh(const volScalarField &h, volScalarField &Kh)
	{
		volScalarField h_clipped = Soil::Math::clipField(h);
		const volScalarField se = 1 / Foam::pow((1 + Foam::pow(Foam::mag(shp_ret_vg_alpha * h_clipped), shp_ret_vg_n)), shp_ret_vg_m);
		Kh = shp_mual_k_sat * 0.5 * ((1 + sign(h)) + (1 - sign(h)) * Foam::pow(se, shp_mual_se_l) * Foam::pow(1 - Foam::pow(1 - Foam::pow(se, 1 / shp_ret_vg_m), shp_ret_vg_m), 2));
		Info<<"Kh estimated, min: "<<min(Kh)<<", avg: "<<average(Kh)<<", max: "<<max(Kh)<<endl;
		return Kh;
	}

	volScalarField &retentionDataVg::Cv(const volScalarField &h, volScalarField &Cv)
	{
		volScalarField h_clipped = Soil::Math::clipField(h);
		Cv = 0.5 * ((1 + sign(h)) * shp_cv_res + (1 - sign(h)) * ((shp_ret_th_s - shp_ret_th_r) * shp_ret_vg_m * shp_ret_vg_n * pow(mag(shp_ret_vg_alpha * h), (shp_ret_vg_n)) * pow((1 + pow(mag(shp_ret_vg_alpha * h), shp_ret_vg_n)), -shp_ret_vg_m - 1)) / mag(h_clipped));
		Info << "Cv estimated, min: " << min(Cv) << ", avg: " << average(Cv) << ", max: " << max(Cv) << endl;
		return Cv;
	}

	volScalarField &retentionDataVg::Theta(const volScalarField &h, volScalarField &Theta)
	{
		Theta = 0.5 * ((1 + sign(h)) * shp_ret_th_s + (1 - sign(h)) * (shp_ret_th_r + (shp_ret_th_s - shp_ret_th_r) * pow((1 + pow(mag(shp_ret_vg_alpha * h), shp_ret_vg_n)), -shp_ret_vg_m)));
		return Theta;
	}

	void retentionDataVg::write(void)
	{
		::retentionData::write();
		shp_ret_vg_alpha.write();
		shp_ret_vg_n.write();
		shp_ret_vg_m.write();
	};

} // namespace Soil::RetentionModels
