#include "zeroFluxFvPatchField.H"
#include "dictionary.H"
#include "fvCFD.H"
#include "gravityMeshObject.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

    zeroFluxFvPatchField::zeroFluxFvPatchField(
        const fvPatch &p,
        const DimensionedField<scalar, volMesh> &iF)
        : zeroGradientFvPatchScalarField(p, iF)
    {
        this->initiateFluxNormVert();
    }

    zeroFluxFvPatchField::zeroFluxFvPatchField(
        const zeroFluxFvPatchField &ptf,
        const fvPatch &p,
        const DimensionedField<scalar, volMesh> &iF,
        const fvPatchFieldMapper &mapper)
        : zeroGradientFvPatchScalarField(p, iF)
    {
        this->initiateFluxNormVert();
        patchType() = ptf.patchType();
        // Map gradient. Set unmapped values and overwrite with mapped ptf
        gradient().map(ptf.gradient(), mapper);

        // Evaluate the value field from the gradient if the internal field is valid
        if (notNull(iF))
        {
            if (iF.size())
            {
                // Note: cannot ask for nf() if zero faces

                scalarField::operator=(
                    // patchInternalField() + gradient()/patch().deltaCoeffs()
                    //  ***HGW Hack to avoid the construction of mesh.deltaCoeffs
                    //  which fails for AMI patches for some mapping operations
                    patchInternalField() + gradient() * (patch().nf() & patch().delta()));
            }
        }
        else
        {
            // Enforce mapping of values so we have a valid starting value
            this->map(ptf, mapper);
        }
    }

    zeroFluxFvPatchField::zeroFluxFvPatchField(
        const fvPatch &p,
        const DimensionedField<scalar, volMesh> &iF,
        const dictionary &dict)
        : zeroGradientFvPatchScalarField(p, iF, dict)
    {
        this->initiateFluxNormVert();
        evaluate();
    }

    zeroFluxFvPatchField::zeroFluxFvPatchField(
        const zeroFluxFvPatchField &ptf)
        : zeroGradientFvPatchScalarField(ptf)
    {
        this->initiateFluxNormVert();
    }

    zeroFluxFvPatchField::zeroFluxFvPatchField(
        const zeroFluxFvPatchField &ptf,
        const DimensionedField<scalar, volMesh> &iF)
        : zeroGradientFvPatchScalarField(ptf, iF)
    {
        this->initiateFluxNormVert();
    }

    // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

    void zeroFluxFvPatchField::autoMap(
        const fvPatchFieldMapper &m)
    {
        zeroGradientFvPatchScalarField::autoMap(m);

        this->initiateFluxNormVert();
    }

    void zeroFluxFvPatchField::rmap(
        const fvPatchScalarField &ptf,
        const labelList &addr)
    {
        zeroGradientFvPatchScalarField::rmap(ptf, addr);

        this->initiateFluxNormVert();
    }

    tmp<scalarField> zeroFluxFvPatchField::valueInternalCoeffs(
        const tmp<scalarField> &) const
    {
        return tmp<scalarField>(new scalarField(this->size(), pTraits<scalar>::one));
    }

    tmp<scalarField> zeroFluxFvPatchField::valueBoundaryCoeffs(
        const tmp<scalarField> &) const
    {
        // TODO: tmp<scalarField>( new  scalarField(gradient()/this->patch().deltaCoeffs())) - here and in other places?
        return gradient() / this->patch().deltaCoeffs();
    }

    tmp<scalarField> zeroFluxFvPatchField::gradientInternalCoeffs() const
    {
        return tmp<scalarField>(
            new scalarField(this->size(), pTraits<scalar>::zero));
    }

    tmp<scalarField> zeroFluxFvPatchField::gradientBoundaryCoeffs() const
    {
        return gradient();
    }

    void zeroFluxFvPatchField::evaluate(const Pstream::commsTypes)
    {
        //Info << "zeroFluxFvPatchField::evaluate()" << endl;
        if (!updated())
        {
            updateCoeffs();
        }
        scalarField::operator=(
            patchInternalField() + gradient() / patch().deltaCoeffs());

        // zeroGradientFvPatchScalarField::evaluate();
    }

    void zeroFluxFvPatchField::updateCoeffs()
    {
        if (updated())
        {
            return;
        }

        this->initiateFluxNormVert(); // FIXME: this shall not be needed
        scalarField::operator=(
            patchInternalField() + gradient() / patch().deltaCoeffs());
        // zeroGradientFvPatchScalarField::updateCoeffs();
    }

    void zeroFluxFvPatchField::initiateFluxNormVert(void)
    {
        //Info << "XXX BC init" << nl;
        vectorField normal = this->patch().nf();
        const uniformDimensionedVectorField &g = meshObjects::gravity::New(this->db().time());
        vector g_n = g.value();
        g_n = g_n / Foam::mag(g_n);
        gradient() = Field<scalar>(this->size(), pTraits<scalar>::one) * (normal & g_n);
        //Info << "zeroFluxFvPatchField::initiateFluxNormVert()" << nl;
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

namespace Foam
{
    makePatchTypeField(
        fvPatchScalarField,
        zeroFluxFvPatchField);
};
// ************************************************************************* //
