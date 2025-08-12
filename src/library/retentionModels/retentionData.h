#ifndef RETENTIONDATA_H
#define RETENTIONDATA_H

#include "fvCFD.H"
#include "retentionDataInterface.h"
#include "simConfig.h"

using namespace Foam;

namespace Soil::RetentionModels
{
    /**
     * @brief Base class for all retention models
     *
     */
    class retentionData : public retentionDataInterface
    {
    public:
        /**
         * @brief Construct a new retention Data object
         *
         * @param transportProperties
         * @param region
         * @param mesh
         * @param runTime
         * @param config
         */
        retentionData(IOdictionary &transportProperties, volScalarField &region, fvMesh &mesh, Time &runTime, simConfig &config);

        inline volScalarField &getDneTau() { return shp_dne_tau; }
        inline volScalarField &getSoilTemp() { return soil_temp; }
        inline volScalarField &getShpCvRes() { return shp_cv_res; }
        inline volScalarField &getShpMualSeL() { return shp_mual_se_l; }
        inline volScalarField &getShpDneTau() { return shp_dne_tau; }
        inline volScalarField &getShpMualKSat() { return shp_mual_k_sat; }
        inline volScalarField &getShpRetThS() { return shp_ret_th_s; }
        inline volScalarField &getShpRetThR() { return shp_ret_th_r; }
        inline volScalarField &getShpRetTotalPor() { return shp_ret_total_por; }

        void errInfo(void);
        inline void write(void) {
            shp_mual_se_l.write();
            shp_dne_tau.write();
            shp_mual_k_sat.write();
            shp_ret_th_s.write();
            shp_ret_th_r.write();
            shp_ret_total_por.write();
        };

    protected:
        int no_of_regions;
        volScalarField soil_temp;  /// soil temperature [K]
        volScalarField shp_cv_res; /// residual value of the soil water differential capacity (defaults to 0) [1/m]
        volScalarField shp_mual_se_l;
        std::vector<dimensionedScalar> shp_mual_se_l_vector;
        volScalarField shp_dne_tau;
        std::vector<dimensionedScalar> shp_dne_tau_vector;
        volScalarField shp_mual_k_sat;
        std::vector<dimensionedScalar> shp_mual_k_sat_vector;
        volScalarField shp_ret_th_s;
        std::vector<dimensionedScalar> shp_ret_th_s_vector;
        volScalarField shp_ret_th_r;
        std::vector<dimensionedScalar> shp_ret_th_r_vector;
        volScalarField shp_ret_total_por;
        std::vector<dimensionedScalar> shp_ret_total_por_vector;
    };

} // namespace Soil::RetentionModels

#endif // RETENTIONDATA_H
