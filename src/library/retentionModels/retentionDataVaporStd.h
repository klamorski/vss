#ifndef RETENTIONDATAVAPORSTD_H
#define RETENTIONDATAVAPORSTD_H

#include "retentionDataVaporNs.h"
#include "retentionData.h"
#include "simConfig.h"

namespace Soil::RetentionModels
{
    /**
     * @brief Base class for liquid capillary and film flow retention models
     *
     */
    class retentionDataVaporStd : public retentionDataVaporNs
    {
    public:
        retentionDataVaporStd(IOdictionary &transportProperties, volScalarField &region, fvMesh &mesh, Time &runTime, simConfig &config);
        volScalarField& Cv(const volScalarField &h, volScalarField &Cv);
        Foam::tmp<Foam::volScalarField> Cv_vapor_debug(const volScalarField &h);
        
        void write(void);
        void debugShpExtra(const Foam::fvMesh &mesh, const Foam::Time &runTime, retentionDataFilm* child);
    protected:
        Foam::tmp<Foam::volScalarField> CvFullVapor(const volScalarField &Cv, const volScalarField &h);
        Foam::tmp<Foam::volScalarField> ThetaVapor(const volScalarField &h, const volScalarField &Theta_cap);
    };

} // namespace Soil::RetentionModels

#endif // RETENTIONDATAVAPORSTD_H
