#include "rampedFluxFvPatchField.H"
#include "dictionary.H"
#include "fvCFD.H"
#include "gravityMeshObject.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
    //- Construct from patch and internal field
    rampedFluxFvPatchField::rampedFluxFvPatchField(
        const fvPatch &p,
        const DimensionedField<scalar, volMesh> &iF)
        : fixedFluxFvPatchField(p, iF),
          hMin_(0),
          fluxToK_(0),
          rampFunction_("steep")
    {
    }

    //- Construct from patch, internal field and dictionary
    rampedFluxFvPatchField::rampedFluxFvPatchField(
        const fvPatch &p,
        const DimensionedField<scalar, volMesh> &iF,
        const dictionary &dict)
        : fixedFluxFvPatchField(p, iF, dict),
          hMin_(dict.get<scalar>("hMin")),
          fluxToK_(dict.get<scalar>("fluxToK")),
          rampFunction_(dict.lookup("rampFunction"))
    {
    }

    //- Construct by mapping the given rampedFluxFvPatchField
    //  onto a new patch
    rampedFluxFvPatchField::rampedFluxFvPatchField(
        const rampedFluxFvPatchField &ptf,
        const fvPatch &p,
        const DimensionedField<scalar, volMesh> &iF,
        const fvPatchFieldMapper &mapper)
        : fixedFluxFvPatchField(ptf, p, iF, mapper),
          hMin_(ptf.hMin_),
          fluxToK_(ptf.fluxToK_),
          rampFunction_(ptf.rampFunction_)
    {
    }

    //- Construct as copy
    rampedFluxFvPatchField::rampedFluxFvPatchField(
        const rampedFluxFvPatchField &ptf)
        : fixedFluxFvPatchField(ptf),
          hMin_(ptf.hMin_),
          fluxToK_(ptf.fluxToK_),
          rampFunction_(ptf.rampFunction_)
    {
    }

    //- Construct as copy setting internal field reference
    rampedFluxFvPatchField::rampedFluxFvPatchField(
        const rampedFluxFvPatchField &ptf,
        const DimensionedField<scalar, volMesh> &iF)
        : fixedFluxFvPatchField(ptf, iF),
          hMin_(ptf.hMin_),
          fluxToK_(ptf.fluxToK_),
          rampFunction_(ptf.rampFunction_)
    {
    }

    // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
    //- Reverse map the given fvPatchField onto this fvPatchField
    void rampedFluxFvPatchField::rmap(
        const fvPatchScalarField &ptf,
        const labelList &addr)
    {
        fvPatchScalarField::rmap(ptf, addr);

        const rampedFluxFvPatchField &fgptf =
            refCast<const rampedFluxFvPatchField>(ptf);

        flux().rmap(fgptf.flux(), addr);

        this->initiateFluxNormVert();
    }

    tmp<scalarField> rampedFluxFvPatchField::calculateGradient(void) const
    {
        const scalarField &K = patch().lookupPatchField<volScalarField>("K");
        Info << "K patch - min(" << min(K) << "), avg(" << average(K) << ")" << endl;
        const scalarField &h = patch().lookupPatchField<volScalarField>("h");
        // Info<<"h - min("<<min(h)<<"), avg("<<average(h)<<")"<<endl;
        const scalar hMinPatch = min(h);
        Info << "hMinPatch: " << hMinPatch << endl;
        const scalar hDiff = (hMin_ - hMinPatch) / hMin_;
        Info << "hDiff: " << hDiff << endl;
        scalar rampTmp = 0;
        if (rampFunction_ == "steep")
        {
            rampTmp = (hDiff > 0) ? 1 : 0;
        }
        else
        {
            if (rampFunction_ == "linear")
            {
                rampTmp = max(min(hDiff, 1), 0);
            }
            else
            {
                FatalErrorInFunction
                    << "Unknown ramping function " << rampFunction_
                    << " for rampedFluxFvPatchField" << exit(FatalError);
            }
        }
        Info << "rampTmp: " << rampTmp << endl;

        return (gradient() - fluxToK_ * rampTmp);
    }

    void rampedFluxFvPatchField::write(Ostream &os) const
    {
        zeroFluxFvPatchField::write(os);
        flux().writeEntry("flux", os);
        os.writeEntry("hMin", hMin_);
        os.writeEntry("fluxToK", fluxToK_);
        os.writeEntry("rampFunction", rampFunction_);
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam

namespace Foam
{

    makePatchTypeField(
        fvPatchScalarField,
        rampedFluxFvPatchField);

}; // End namespace Foam

// ************************************************************************* //
