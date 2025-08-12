#ifndef SIMCONFIG_H
#define SIMCONFIG_H

#include "fvCFD.H"

class Processes {
public:
    Processes(const dictionary &transportProperties);
    Processes() {};
    inline bool isWaterCapillaryFlow() { return procesesWaterCapillaryFlow; };
    inline word getWaterVaporFlow() { return procesesWaterVaporFlow; };
    inline bool isWaterFilmFlow() { return procesesWaterFilmFlow; };
    inline bool isHeatConductionFlow() { return procesesHeatConductionFlow; };
    inline bool isHeatConvectionFlow() { return procesesHeatConvectionFlow; };
    inline bool isHeatVaporFlow() { return procesesHeatVaporFlow; };
    inline word getDneModel() { return dneModel; };
    inline bool isDneModelSet() { return dneModelSet; };

private:
        bool procesesWaterCapillaryFlow;
        word procesesWaterVaporFlow;
        bool procesesWaterFilmFlow;
        bool procesesHeatConductionFlow;
        bool procesesHeatConvectionFlow;
        bool procesesHeatVaporFlow;
        word dneModel;
        bool dneModelSet;
};

class Algorithm {
public:
    enum kIntEstMethod {now, nowOldArithmeticMean, nowOldHarmonicMean, nowOldGeometricMean};
public:
    Algorithm(const dictionary &transportProperties);
    Algorithm() {};
    inline bool isCorrectBcK() { return correctBcK; };
    inline bool isCorrectBcCv() { return correctBcCv; };
    inline bool isCorrectBcTheta() { return correctBcTheta; };
    inline scalar getHPicardResidual() { return hPicardResidual; };
    inline scalar getKPicardResidual() { return kPicardResidual; };
    inline scalar getCvPicardResidual() { return cvPicardResidual; };
    inline scalar getPicardMaxIter() { return picardMaxIter; };
    inline scalar getShpCvRes() { return shpCvRes; };
    inline scalar getShpKLim() { return shpKLim; };
    inline kIntEstMethod getKInterstepEstimationMethod() { return kInterstepEstimationMethod; };
    inline bool getKInterstepEstimationDoSurfaceInterpolation() { return kInterstepEstimationDoSurfaceInterpolation; };
    inline string getMbcHeaderInfo() { return mbcHeaderInfo; };
    inline string getMbcRowInfo() { return mbcRowInfo; };
    inline string getReFormulation() { return reFormulation; };
    inline bool isReFormulationPressureHead() { return reFormulation == "pressure_head"; };
    inline bool isReFormulationMixed() { return reFormulation == "mixed"; };
    inline bool getShpDebug() { return shpDebug; };
    inline void setShpDebug(bool shpDebug) { this->shpDebug = shpDebug; };

private:
    bool correctBcK;
    bool correctBcCv;
    bool correctBcTheta;
    scalar hPicardResidual;
    scalar kPicardResidual;
    scalar cvPicardResidual;
    scalar picardMaxIter;
    scalar shpCvRes;
    scalar shpKLim;
    kIntEstMethod kInterstepEstimationMethod;
    bool kInterstepEstimationDoSurfaceInterpolation;
    word mbcHeaderInfo;
    word mbcRowInfo;
    word reFormulation;
    bool shpDebug;
};

class simConfig
{
public:
    simConfig(const dictionary &transportProperties);
    Processes processes;
    Algorithm algorithm;
};

#endif // SIMCONFIG_H
