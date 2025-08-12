#ifndef fixedFluxFvPatchField_H
#define fixedFluxFvPatchField_H

#include "fvPatchFields.H"
#include "zeroFluxFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
                     Class fixedFluxFvPatchField Declaration
    \*---------------------------------------------------------------------------*/
    class fixedFluxFvPatchField
        : public zeroFluxFvPatchField // TODO: Use zeroGradientFvPatchField as the parent class
    {
        // Private data
        scalarField flux_; // flux initiated from patch boundary definition

    public:
        //- Runtime type information
        TypeName("fixedFlux");

        // Constructors
        //- Construct from patch and internal field
        fixedFluxFvPatchField(
            const fvPatch &,
            const DimensionedField<scalar, volMesh> &);

        //- Construct from patch, internal field and dictionary
        fixedFluxFvPatchField(
            const fvPatch &,
            const DimensionedField<scalar, volMesh> &,
            const dictionary &);

        //- Construct by mapping the given fixedFluxFvPatchField
        //  onto a new patch
        fixedFluxFvPatchField(
            const fixedFluxFvPatchField &,
            const fvPatch &,
            const DimensionedField<scalar, volMesh> &,
            const fvPatchFieldMapper &);

        //- Construct as copy
        fixedFluxFvPatchField(
            const fixedFluxFvPatchField &);

        //- Construct as copy setting internal field reference
        fixedFluxFvPatchField(
            const fixedFluxFvPatchField &,
            const DimensionedField<scalar, volMesh> &);

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>(
                new fixedFluxFvPatchField(*this));
        }

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone(
            const DimensionedField<scalar, volMesh> &iF) const
        {
            return tmp<fvPatchScalarField>(
                new fixedFluxFvPatchField(*this, iF));
        }

        // Mapping functions
        //- Map (and resize as needed) from self given a mapping object
        virtual void autoMap(
            const fvPatchFieldMapper &);

        //- Reverse map the given fvPatchField onto this fvPatchField
        virtual void rmap(
            const fvPatchScalarField &,
            const labelList &);

        // Evaluation functions
        //- Return gradient at boundary
        virtual tmp<scalarField> snGrad() const
        {
            return calculateGradient();
        }

        //- Evaluate the patch field
        virtual void evaluate(
            const Pstream::commsTypes commsType = Pstream::commsTypes::blocking);

        virtual void updateCoeffs();

        //- Return the matrix diagonal coefficients corresponding to the
        //  evaluation of the value of this patchField with given weights
        virtual tmp<scalarField> valueInternalCoeffs(
            const tmp<scalarField> &) const;

        //- Return the matrix source coefficients corresponding to the
        //  evaluation of the value of this patchField with given weights
        virtual tmp<scalarField> valueBoundaryCoeffs(
            const tmp<scalarField> &) const;

        //- Return the matrix diagonal coefficients corresponding to the
        //  evaluation of the gradient of this patchField
        virtual tmp<scalarField> gradientInternalCoeffs() const;

        //- Return the matrix source coefficients corresponding to the
        //  evaluation of the gradient of this patchField
        virtual tmp<scalarField> gradientBoundaryCoeffs() const;

        virtual tmp<scalarField> calculateGradient(void) const;

        //- Write
        virtual void write(Ostream &) const;

        virtual scalarField &flux()
        {
            return flux_;
        };

        virtual const scalarField &flux() const
        {
            return flux_;
        };
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
