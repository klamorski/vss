#ifndef SOLVER_FUNCTIONS_H
#define SOLVER_FUNCTIONS_H

#include "fvCFD.H"
#include "retentionDataInterface.h"

using namespace Soil::RetentionModels;

void debugShp(const Foam::fvMesh &mesh, const Foam::Time &runTime, retentionDataInterface* shpModel);

#endif // SOLVER_FUNCTIONS_H
