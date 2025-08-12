#ifndef RETENTIONDATAKS_H
#define RETENTIONDATAKS_H

#include "retentionData.h"

namespace Soil::RetentionModels
{
    /**
     * @brief Base class for liquid capillary flow implementing Kosugi retention model
     *
     */
    class retentionDataKs : virtual public retentionData
    {
    public:
        retentionDataKs(IOdictionary &transportProperties, volScalarField &region, fvMesh &mesh, Time &runTime, simConfig &config);
        volScalarField& Kh(const volScalarField &h, volScalarField &Kh);
        volScalarField& Cv(const volScalarField &h, volScalarField &Cv);
        volScalarField& Theta(const volScalarField &h, volScalarField &Theta);

        inline volScalarField &getShpRetKsSigma() { return shp_ret_ks_sigma; }
        inline volScalarField &getShpRetKsHM() { return shp_ret_ks_h_m; }

        void write(void);

    protected:
        volScalarField shp_ret_ks_sigma;
        std::vector<dimensionedScalar> shp_ret_ks_sigma_vector;
        volScalarField shp_ret_ks_h_m;
        std::vector<dimensionedScalar> shp_ret_ks_h_m_vector;
    };

} // namespace Soil::RetentionModels

#endif // RETENTIONDATAKS_H
