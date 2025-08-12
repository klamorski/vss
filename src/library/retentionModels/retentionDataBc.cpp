#include "retentionDataBc.h"
#include "soilMath.h"

using namespace Soil::RetentionModels;

namespace Soil::RetentionModels
{
	retentionDataBc::retentionDataBc(IOdictionary &transportProperties, volScalarField &region, fvMesh &mesh, Time &runTime, simConfig &config)
		: retentionData(transportProperties, region, mesh, runTime, config),
		  shp_ret_bc_lambda(IOobject("shp_ret_bc_lambda", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("shp_ret_bc_lambda", dimensionSet(0, 0, 0, 0, 0, 0, 0), 1)),
		  shp_ret_bc_h_e(IOobject("shp_ret_bc_h_e", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("shp_ret_bc_h_e", dimensionSet(0, 1, 0, 0, 0, 0, 0), 1))
	{
		if (no_of_regions > 0)
		{
			char buffer[1000];

			// TODO Test and eventually report error when there will be no appriopriate definitions for region in the dictionary
			// TODO Test if field region really has no_of_regions values defined. Will not have if user forgot to prepare 0/region file than field will be initiated to default 0 for all mesh

			shp_ret_bc_lambda_vector.reserve(no_of_regions);
			for (int i = 0; i < no_of_regions; i++)
			{
				sprintf(buffer, "shp_ret_bc_lambda_%d", i);
				std::string key_name(buffer);
				shp_ret_bc_lambda_vector.push_back(dimensionedScalar(key_name.c_str(), dimensionSet(0, 0, 0, 0, 0, 0, 0), transportProperties));
				// cout<<buffer<<endl;
			}

			shp_ret_bc_h_e_vector.reserve(no_of_regions);
			for (int i = 0; i < no_of_regions; i++)
			{
				sprintf(buffer, "shp_ret_bc_h_e_%d", i);
				std::string key_name(buffer);
				shp_ret_bc_h_e_vector.push_back(dimensionedScalar(key_name.c_str(), dimensionSet(0, 1, 0, 0, 0, 0, 0), transportProperties));
				// cout<<buffer<<endl;
			}

			// initiate internal fields
			forAll(region, index)
			{
				shp_ret_bc_lambda[index] = shp_ret_bc_lambda_vector[static_cast<int>(region[index])].value();
				shp_ret_bc_h_e[index] = shp_ret_bc_h_e_vector[static_cast<int>(region[index])].value();
			}

			// initiate boundary fields
			forAll(region.boundaryFieldRef(), ipatch)
			{
				// each face in the the patch
				forAll(region.boundaryFieldRef()[ipatch], index)
				{
					shp_ret_bc_lambda.boundaryFieldRef()[ipatch][index] = shp_ret_bc_lambda_vector[static_cast<int>(region.boundaryField()[ipatch][index])].value();
					shp_ret_bc_h_e.boundaryFieldRef()[ipatch][index] = shp_ret_bc_h_e_vector[static_cast<int>(region.boundaryField()[ipatch][index])].value();
				}
			}
		}
		else
		{
			shp_ret_bc_lambda = volScalarField(IOobject("shp_ret_bc_lambda", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
			shp_ret_bc_h_e = volScalarField(IOobject("shp_ret_bc_h_e", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);
		}
	}
	
    volScalarField& retentionDataBc::Kh(const volScalarField &h, volScalarField &Kh) {
		volScalarField iahn = (1 + Foam::pow(Foam::mag(shp_ret_bc_lambda * h), shp_ret_bc_h_e));
		Kh = shp_mual_k_sat * 0.5 * ((1 + sign(h - shp_ret_bc_h_e)) + (1 - sign(h - shp_ret_bc_h_e)) * Foam::pow(mag(shp_ret_bc_h_e / h), shp_ret_bc_lambda * (shp_mual_se_l + 2 + 2 / shp_ret_bc_lambda)));
		return Kh;
	}

	Foam::volScalarField &retentionDataBc::Cv(const volScalarField &h, volScalarField &Cv)
	{
		volScalarField h_clipped = Soil::Math::clipField(h);
		Cv = 0.5 * ((1 + sign(h - shp_ret_bc_h_e)) * shp_cv_res + (1 - sign(h - shp_ret_bc_h_e)) * ((shp_ret_th_s - shp_ret_th_r) * shp_ret_bc_lambda * pow(shp_ret_bc_h_e, shp_ret_bc_lambda) * pow(mag(h_clipped), -shp_ret_bc_lambda - 1)));
		return Cv;
	}

	Foam::volScalarField &retentionDataBc::Theta(const volScalarField &h, volScalarField &Theta)
	{
		Theta = 0.5 * ((1 + sign(h - shp_ret_bc_h_e)) * shp_ret_th_s + (1 - sign(h - shp_ret_bc_h_e)) * (shp_ret_th_r + (shp_ret_th_s - shp_ret_th_r) * pow(mag(shp_ret_bc_h_e / h), shp_ret_bc_lambda)));
		return Theta;
	}

	void retentionDataBc::write(void) {
		::retentionData::write();
		shp_ret_bc_lambda.write();
		shp_ret_bc_h_e.write();
	};

} // namespace Soil::RetentionModels
