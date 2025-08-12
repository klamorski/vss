#ifndef RETENTIONDATAINTERFACE_H
#define RETENTIONDATAINTERFACE_H

#include "fvCFD.H"
#include "../external/interpolate/src/libInterpolate/Interpolate.hpp"
#include "simConfig.h"

namespace Soil::RetentionModels
{

    /**
     * @brief Base class defining interface for all retention models
     *
     */
    class retentionDataInterface
    {
    protected:
        simConfig config;
        
    public:
        retentionDataInterface(fvMesh &mesh, Time &runTime, simConfig &config);
        
        /**
         * @brief Unsaturated hydraulic conductivity
         *
         * @param h pressure head
         * @return K hydraulic conductivity
         */
        Foam::tmp<Foam::volScalarField> Kh(const volScalarField &h);
        virtual volScalarField& Kh(const volScalarField &h, volScalarField &Kh);
        Foam::tmp<Foam::volScalarField> Cv(const volScalarField &h);
        virtual volScalarField& Cv(const volScalarField &h, volScalarField &Cv);
        Foam::tmp<Foam::volScalarField> Theta(const volScalarField &h);
        virtual volScalarField& Theta(const volScalarField &h, volScalarField &Theta);
        //Foam::tmp<Foam::volScalarField> h(volScalarField &Theta);
        virtual ~retentionDataInterface() {};


        virtual volScalarField &getDneTau();
        virtual volScalarField &getSoilTemp();
        virtual volScalarField &getShpCvRes();
        virtual volScalarField &getShpMualSeL();
        virtual volScalarField &getShpDneTau();
        virtual volScalarField &getShpMualKSat();
        virtual volScalarField &getShpRetThS();
        virtual volScalarField &getShpRetThR();
        virtual volScalarField &getShpRetTotalPor();

        void errInfo(void);
        virtual void write(void) {};

        
        virtual volScalarField &getShpFilmHa() { 
            errInfo();
            return *static_cast<volScalarField *>(0);
        }

        virtual Foam::tmp<Foam::volScalarField> Klf(const volScalarField &h) {
            errInfo();
            return *static_cast<volScalarField *>(0);
        }
    };

} // namespace Soil::RetentionModels

#endif // RETENTIONDATA_H
