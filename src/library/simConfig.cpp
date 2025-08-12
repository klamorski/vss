#include "simConfig.h"

simConfig::simConfig(const dictionary &transportProperties)
{
    processes = Processes(transportProperties);
    algorithm = Algorithm(transportProperties);
};

Processes::Processes(const dictionary &transportProperties) {
    procesesWaterCapillaryFlow = false;
    procesesWaterVaporFlow = word("none");
    procesesWaterFilmFlow = false;
    procesesHeatConductionFlow = false;
    procesesHeatConvectionFlow = false;
    procesesHeatVaporFlow = false;
    dneModel = word("none");
    dneModelSet = false;


    if (const dictionary *processesDictPtr = transportProperties.findDict("processes"))
    {
        Info << "Found 'processes' dictionary in transportProperties file" << nl;
        if (const dictionary *waterDictPtr = processesDictPtr->findDict("water"))
        {
            Info << "Found 'water' dictionary in transportProperties file" << nl;
            procesesWaterCapillaryFlow = waterDictPtr->getOrDefault("capillary_flow", false);
            procesesWaterVaporFlow = waterDictPtr->getOrDefault("vapor_flow",  word("none"));
            procesesWaterFilmFlow = waterDictPtr->getOrDefault("film_flow", false);
            if (!procesesWaterCapillaryFlow && (procesesWaterVaporFlow == "none" || procesesWaterFilmFlow))
            {
                FatalErrorInFunction << "Vapor and film flow can be simulated only if capillary flow is enabled" << exit(FatalError);
            }
        }
        else
        {
            FatalErrorInFunction << "No 'water' dictionary found in transportProperties file" << exit(FatalError);
        }
        if (const dictionary *heatDictPtr = processesDictPtr->findDict("heat"))
        {
            Info << "Found 'heat' dictionary in transportProperties file" << nl;
            procesesHeatConductionFlow = heatDictPtr->getOrDefault("conduction_flow", false);
            procesesHeatConvectionFlow = heatDictPtr->getOrDefault("convection_flow", false);
            procesesHeatVaporFlow = heatDictPtr->getOrDefault("vapor_flow", false);
            if (!procesesHeatConductionFlow && (procesesHeatConvectionFlow || procesesHeatVaporFlow))
            {
                FatalErrorInFunction << "Convection and vapor flow can be simulated only if conduction flow is enabled" << exit(FatalError);
            }
        }
        else
        {
            FatalErrorInFunction << "No 'heat' dictionary found in transportProperties file" << exit(FatalError);
        }
        
		dneModel = word(processesDictPtr->lookupOrDefault("dne_model", word("none")));

        if (dneModel == "none")
        {
            dneModelSet = false;
        } else {
            if (dneModel == "ross_smettem")
            {
                Info << "Ross-Smettem model is not implemented yet." << nl;
                exit(0);
                dneModelSet = true;
            } else {
                Info << "Unknown DNE model: " << dneModel << nl;
                exit(0);
            }
        }
    }
    else
    {
        FatalErrorInFunction << "No 'processes' dictionary found in transportProperties file" << exit(FatalError);
    }

}


Algorithm::Algorithm(const dictionary &transportProperties) {
    correctBcK = false;
    correctBcCv = false;
    correctBcTheta = false;
    hPicardResidual = 1e-6;
    kPicardResidual = 0.0;
    cvPicardResidual = 0.0;
    picardMaxIter = 100;
    shpCvRes = 0.0;
    kInterstepEstimationMethod = kIntEstMethod::now;
    kInterstepEstimationDoSurfaceInterpolation = false;
    mbcHeaderInfo = string("");
    mbcRowInfo = string("");
    reFormulation = word("");
    
    if (const dictionary *algorithmDictPtr = transportProperties.findDict("algorithm"))
    {
        Info << "Found 'algorithm' dictionary in transportProperties file" << nl;
        correctBcK = algorithmDictPtr->getOrDefault("k_correct_bc", false);
        correctBcCv = algorithmDictPtr->getOrDefault("cv_correct_bc", false);
        correctBcTheta = algorithmDictPtr->getOrDefault("theta_correct_bc", false);
        hPicardResidual = algorithmDictPtr->getOrDefault("h_picard_residual", 1e-6);
        kPicardResidual = algorithmDictPtr->getOrDefault("k_picard_residual", 0.0);
        cvPicardResidual = algorithmDictPtr->getOrDefault("cv_picard_residual", 0.0);
        picardMaxIter = algorithmDictPtr->getOrDefault("picard_max_iter", 100);
        shpCvRes = algorithmDictPtr->getOrDefault("shp_cv_res", 0.0);
        mbcHeaderInfo = string(algorithmDictPtr->lookupOrDefault("mbc_header_info", string("")));
        mbcRowInfo = string(algorithmDictPtr->lookupOrDefault("mbc_row_info", string("")));
        reFormulation = word(algorithmDictPtr->lookupOrDefault("re_formulation", word("pressure_head")));
        cout << "mbcHeaderInfo: " << mbcHeaderInfo << nl;
        cout << "mbcRowInfo: " << mbcRowInfo << nl;
        word kInterstepEstimationMethodTmp = word(algorithmDictPtr->lookupOrDefault("k_interstep_estimation_method", word("now")));
        if (kInterstepEstimationMethodTmp == "now")
        {
            kInterstepEstimationMethod = kIntEstMethod::now;
        } else {
            if (kInterstepEstimationMethodTmp == "now-old-arithmetic-mean")
            {
                kInterstepEstimationMethod = kIntEstMethod::nowOldArithmeticMean;
            } else {
                if (kInterstepEstimationMethodTmp == "now-old-geometrical-mean")
                {
                    kInterstepEstimationMethod = kIntEstMethod::nowOldGeometricMean;
                } else {
                    if (kInterstepEstimationMethodTmp == "now-old-harmonic-mean")
                    {
                        kInterstepEstimationMethod = kIntEstMethod::nowOldHarmonicMean;
                    } else {
                        Info << "Unknown kInterstepEstimationMethod: " << kInterstepEstimationMethodTmp << nl;
                        exit(0);
                    }
                }
            }
        }
        kInterstepEstimationDoSurfaceInterpolation = algorithmDictPtr->getOrDefault("k_interstep_estimation_do_surface_interpolation", false);
        Info<<"Algorithm parameters: hPicardResidual: "<<hPicardResidual<<", kPicardResidual: "<<kPicardResidual<<", cvPicardResidual: "<<cvPicardResidual<<", shpCvRes: "<<shpCvRes<<", shpKLim: "<<shpKLim<<nl;
    }
    else
    {
        FatalErrorInFunction << "No 'algorithm' dictionary found in transportProperties file" << exit(FatalError);
    }
}
