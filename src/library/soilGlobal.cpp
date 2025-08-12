#include "soilGlobal.h"
#include "dimensionedScalar.H"
#include "error.H"

namespace Soil::Global
{
    void checkUnits(const dimensionedScalar &arg, const dimensionSet &expectedUnits, const std::string &argName)
    {
        if (arg.dimensions() != expectedUnits)
        {
            FatalErrorInFunction
                << "Incorrect units for argument " << argName << ". Expected " << expectedUnits << " but got " << arg.dimensions() << Foam::abort(Foam::FatalError);
        }
    }

    void checkUnits(const volScalarField &arg, const dimensionSet &expectedUnits, const std::string &argName)
    {
        if (arg.dimensions() != expectedUnits)
        {
            FatalErrorInFunction
                << "Incorrect units for argument " << argName << ". Expected " << expectedUnits << " but got " << arg.dimensions() << Foam::abort(Foam::FatalError);
        }
    }

} // namespace Soil::Global


