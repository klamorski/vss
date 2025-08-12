#ifndef RETENTIONDATAMVG_H
#define RETENTIONDATAMVG_H

#include "fvCFD.H"
#include "retentionData.h"

#include <vector>
#include <string>
#include <stdio.h>

namespace Soil::RetentionModels
{

    class retentionDataVg : virtual public retentionData
    {
    public:
        retentionDataVg(IOdictionary &transportProperties, volScalarField &region, fvMesh &mesh, Time &runTime, simConfig &config);
        volScalarField &Kh(const volScalarField &h, volScalarField &Kh);
        volScalarField &Cv(const volScalarField &h, volScalarField &Cv);
        volScalarField &Theta(const volScalarField &h, volScalarField &Theta);


        inline volScalarField &getShpRetVgAlpha() { return shp_ret_vg_alpha; };
        inline volScalarField &getShpRetVgN() { return shp_ret_vg_n; }
        inline volScalarField &getShpRetVgM() { return shp_ret_vg_m; }

        void write(void);

    protected:
        volScalarField shp_ret_vg_alpha;
        std::vector<dimensionedScalar> shp_ret_vg_alpha_vector;
        volScalarField shp_ret_vg_n;
        std::vector<dimensionedScalar> shp_ret_vg_n_vector;
        volScalarField shp_ret_vg_m;
        std::vector<dimensionedScalar> shp_ret_vg_m_vector;
    };

} // namespace Soil::RetentionModels

#endif // RETENTIONDATAMVG_H
