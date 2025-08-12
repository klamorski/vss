#ifndef RETENTIONDATABC_H
#define RETENTIONDATABC_H

#include "fvCFD.H"
#include <vector>
#include <string>
#include <stdio.h>
#include "retentionData.h"

namespace Soil::RetentionModels
{
    /**
     * @brief Brook-Corey capillary flow retention model
     *
     */
    class retentionDataBc : public retentionData
    {
    public:
        retentionDataBc(IOdictionary &transportProperties, volScalarField &region, fvMesh &mesh, Time &runTime, simConfig &config);
        volScalarField& Kh(const volScalarField &h, volScalarField &Kh);
        volScalarField& Cv(const volScalarField &h, volScalarField &Cv);
        volScalarField& Theta(const volScalarField &h, volScalarField &Theta);

        inline volScalarField &getShpRetBcLambda() { return shp_ret_bc_lambda; }
        inline volScalarField &getShpRetBcHE() { return shp_ret_bc_h_e; }

        void write(void);
        
    protected:
    private:
        volScalarField shp_ret_bc_lambda;
        std::vector<dimensionedScalar> shp_ret_bc_lambda_vector;
        volScalarField shp_ret_bc_h_e;
        std::vector<dimensionedScalar> shp_ret_bc_h_e_vector;
    };

} // namespace Soil::RetentionModels

#endif // RETENTIONDATABC_H
