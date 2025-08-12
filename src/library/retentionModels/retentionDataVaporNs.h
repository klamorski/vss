#ifndef RETENTIONDATAVAPORNS_H
#define RETENTIONDATAVAPORNS_H

#include "retentionDataInterface.h"
#include "retentionDataFilm.h"
#include "retentionData.h"
#include "simConfig.h"

namespace Soil::RetentionModels
{
    /**
     * @brief Base class for liquid capillary and film flow retention models
     *
     */
    class retentionDataVaporNs : virtual public retentionData
    {
    public:
        retentionDataVaporNs(IOdictionary &transportProperties, volScalarField &region, fvMesh &mesh, Time &runTime, simConfig &config);
        Foam::tmp<Foam::volScalarField> Kv(const volScalarField &h);
        volScalarField& Kv(const volScalarField &h, volScalarField &Kh);

        inline volScalarField &getShpVaporBeta() { return shp_vapor_beta; }
        
        void write(void);

    protected: 
        volScalarField shp_vapor_beta;
        std::vector<dimensionedScalar> shp_vapor_beta_vector;


    public:
        void debugShpExtra(const Foam::fvMesh &mesh, const Foam::Time &runTime, retentionDataFilm* child);
    };

} // namespace Soil::RetentionModels

#endif // RETENTIONDATAVAPORNS_H
