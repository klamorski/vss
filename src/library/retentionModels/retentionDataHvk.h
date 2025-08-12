#ifndef RETENTIONDATAHVK_H
#define RETENTIONDATAHVK_H

#include "fvCFD.H"
#include "retentionData.h"

#include <vector>
#include <string>
#include <stdio.h>

namespace Soil::RetentionModels
{

    class retentionDataHvk : public retentionData
    {
    public:
        retentionDataHvk(IOdictionary &transportProperties, volScalarField &region, fvMesh &mesh, Time &runTime, simConfig &config);
        volScalarField& Kh(const volScalarField &h, volScalarField &Kh);
        volScalarField& Cv(const volScalarField &h, volScalarField &Cv);
        volScalarField& Theta(const volScalarField &h, volScalarField &Theta);


        inline volScalarField &getShpRetHvkAlpha() { return shp_ret_hvk_alpha; }
        inline volScalarField &getShpRetHvkBeta() { return shp_ret_hvk_beta; }
        inline volScalarField &getShpRetHvkA() { return shp_ret_hvk_a; }
        inline volScalarField &getShpRetHvkGamma() { return shp_ret_hvk_gamma; }

        void write(void);

    protected:
        volScalarField shp_ret_hvk_alpha;
        std::vector<dimensionedScalar> shp_ret_hvk_alpha_vector;
        volScalarField shp_ret_hvk_beta;
        std::vector<dimensionedScalar> shp_ret_hvk_beta_vector;
        volScalarField shp_ret_hvk_a;
        std::vector<dimensionedScalar> shp_ret_hvk_a_vector;
        volScalarField shp_ret_hvk_gamma;
        std::vector<dimensionedScalar> shp_ret_hvk_gamma_vector;
    };

} // namespace Soil::RetentionModels

#endif // RETENTIONDATAHVK_H
