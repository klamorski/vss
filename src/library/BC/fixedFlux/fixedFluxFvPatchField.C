#include "fixedFluxFvPatchField.H"
#include "dictionary.H"
#include "fvCFD.H"
#include "gravityMeshObject.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
    //- Construct from patch and internal field
    fixedFluxFvPatchField::fixedFluxFvPatchField(
        const fvPatch &p,
        const DimensionedField<scalar, volMesh> &iF)
        : zeroFluxFvPatchField(p, iF),
          flux_(p.size(), pTraits<scalar>::zero)
    {
    }

    //- Construct from patch, internal field and dictionary
    fixedFluxFvPatchField::fixedFluxFvPatchField(
        const fvPatch &p,
        const DimensionedField<scalar, volMesh> &iF,
        const dictionary &dict)
        : zeroFluxFvPatchField(p, iF, dict),
          flux_("flux", dict, p.size())
    {
        this->initiateFluxNormVert();
        evaluate();
    }

    //- Construct by mapping the given fixedFluxFvPatchField
    //  onto a new patch
    fixedFluxFvPatchField::fixedFluxFvPatchField(
        const fixedFluxFvPatchField &ptf,
        const fvPatch &p,
        const DimensionedField<scalar, volMesh> &iF,
        const fvPatchFieldMapper &mapper)
        : zeroFluxFvPatchField(ptf, p, iF, mapper),
          flux_(ptf.flux_, mapper)
    {
        this->initiateFluxNormVert();

        if (&iF && mapper.hasUnmapped())
        {
            WarningIn(
                "fixedFluxFvPatchField::fixedFluxFvPatchField\n"
                "(\n"
                "    const fixedFluxFvPatchField&,\n"
                "    const fvPatch&,\n"
                "    const DimensionedField<scalar, volMesh>&,\n"
                "    const fvPatchFieldMapper&\n"
                ")\n")
                << "On field " << iF.name() << " patch " << p.name()
                << " patchField " << this->type()
                << " : mapper does not map all values." << nl
                << "    To avoid this warning fully specify the mapping in derived"
                << " patch fields." << endl;
        }
    }

    //- Construct as copy
    fixedFluxFvPatchField::fixedFluxFvPatchField(
        const fixedFluxFvPatchField &ptf)
        : zeroFluxFvPatchField(ptf),
          flux_(ptf.flux_)
    {
        this->initiateFluxNormVert();
    }

    //- Construct as copy setting internal field reference
    fixedFluxFvPatchField::fixedFluxFvPatchField(
        const fixedFluxFvPatchField &ptf,
        const DimensionedField<scalar, volMesh> &iF)
        : zeroFluxFvPatchField(ptf, iF),
          flux_(ptf.flux_)
    {
        this->initiateFluxNormVert();
    }

    // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
    //- Map (and resize as needed) from self given a mapping object
    void fixedFluxFvPatchField::autoMap(
        const fvPatchFieldMapper &m)
    {
        fvPatchScalarField::autoMap(m);
        flux().autoMap(m);

        this->initiateFluxNormVert();
    }

    //- Reverse map the given fvPatchField onto this fvPatchField
    void fixedFluxFvPatchField::rmap(
        const fvPatchScalarField &ptf,
        const labelList &addr)
    {
        fvPatchScalarField::rmap(ptf, addr);

        const fixedFluxFvPatchField &fgptf =
            refCast<const fixedFluxFvPatchField>(ptf);

        flux().rmap(fgptf.flux(), addr);

        this->initiateFluxNormVert();
    }

    void fixedFluxFvPatchField::evaluate(const Pstream::commsTypes)
    {
        if (!this->updated())
        {
            updateCoeffs();
        }

        scalarField::operator=(
            patchInternalField() + calculateGradient() / patch().deltaCoeffs());

        // fvPatchScalarField::evaluate();
    }

    void fixedFluxFvPatchField::updateCoeffs()
    {
        if (updated())
        {
            return;
        }

        this->initiateFluxNormVert(); // FIXME: this shall not be needed

        scalarField::operator=(
            patchInternalField() + calculateGradient() / patch().deltaCoeffs());

        setUpdated(true);
        // zeroFluxNgFvPatchField::updateCoeffs();
    }

    tmp<scalarField> fixedFluxFvPatchField::valueInternalCoeffs(
        const tmp<scalarField> &) const
    {
        return tmp<scalarField>(new scalarField(this->size(), pTraits<scalar>::one));
    }

    tmp<scalarField> fixedFluxFvPatchField::valueBoundaryCoeffs(
        const tmp<scalarField> &) const
    {
        return calculateGradient() / patch().deltaCoeffs();
    }

    tmp<scalarField> fixedFluxFvPatchField::gradientInternalCoeffs() const
    {
        return tmp<scalarField>(
            new scalarField(this->size(), pTraits<scalar>::zero));
    }

    tmp<scalarField> fixedFluxFvPatchField::gradientBoundaryCoeffs() const
    {
        return calculateGradient();
    }

    void fixedFluxFvPatchField::write(Ostream &os) const
    {
        zeroFluxFvPatchField::write(os);
        flux().writeEntry("flux", os);
    }

    tmp<scalarField> fixedFluxFvPatchField::calculateGradient(void) const
    {
        // FIXME: Why these two K scalarField approaches gives different results?!
        // const scalarField &K = patch().lookupPatchField<volScalarField>("K");

        const scalarField &K = (this->db().objectRegistry::foundObject<volScalarField>("K") == true) ? static_cast<scalarField>(this->db().objectRegistry::lookupObject<volScalarField>("K").boundaryField()[this->patch().index()]) : Field<double>(this->size(), pTraits<double>::one);
        if (!this->db().objectRegistry::foundObject<volScalarField>("K"))
        {
            Info << "Warning! Dummy K field initialized to scalar(1)." << endl;
        }

        return (gradient() - flux() / K);
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam

namespace Foam
{

    makePatchTypeField(
        fvPatchScalarField,
        fixedFluxFvPatchField);

}; // End namespace Foam

// ************************************************************************* //
