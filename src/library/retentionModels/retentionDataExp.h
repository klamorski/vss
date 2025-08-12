#ifndef RETENTIONDATAEXP_H
#define RETENTIONDATAEXP_H

#include "fvCFD.H"
#include "retentionData.h"
#include <vector>
#include <string>
#include <stdio.h>

namespace Soil::RetentionModels
{

    class retentionDataExp : public retentionData
    {
    public:
        retentionDataExp(IOdictionary &transportProperties, volScalarField &region, fvMesh &mesh, Time &runTime, simConfig &config);
        volScalarField& Kh(const volScalarField &h, volScalarField &Kh);
        volScalarField& Cv(const volScalarField &h, volScalarField &Cv);
        volScalarField& Theta(const volScalarField &h, volScalarField &Theta);

        inline volScalarField &getShpRetExpAlpha() { return shp_ret_exp_alpha; }

        void write(void);

    protected:
    private:
        volScalarField shp_ret_exp_alpha;
        std::vector<dimensionedScalar> shp_ret_exp_alpha_vector;
    };

} // namespace Soil::RetentionModels

#endif // RETENTIONDATAEXP_H
