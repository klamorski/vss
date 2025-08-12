#ifndef zeroFluxFvPatchField_H
#define zeroFluxFvPatchField_H
#include "fvPatchFields.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
                     Class zeroFluxFvPatchField Declaration
    \*---------------------------------------------------------------------------*/
    class zeroFluxFvPatchField
        : public zeroGradientFvPatchScalarField
    {
        // Private data
        scalarField gradient_;

    public:
        //- Runtime type information
        TypeName("zeroFlux");

        // Constructors
        //- Construct from patch and internal field
        zeroFluxFvPatchField(
            const fvPatch &,
            const DimensionedField<scalar, volMesh> &);

        //- Construct from patch, internal field and dictionary
        zeroFluxFvPatchField(
            const fvPatch &,
            const DimensionedField<scalar, volMesh> &,
            const dictionary &);

        //- Construct by mapping the given zeroFluxFvPatchField
        //  onto a new patch
        zeroFluxFvPatchField(
            const zeroFluxFvPatchField &,
            const fvPatch &,
            const DimensionedField<scalar, volMesh> &,
            const fvPatchFieldMapper &);

        //- Construct as copy
        zeroFluxFvPatchField(
            const zeroFluxFvPatchField &);

        //- Construct as copy setting internal field reference
        zeroFluxFvPatchField(
            const zeroFluxFvPatchField &,
            const DimensionedField<scalar, volMesh> &);

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>(
                new zeroFluxFvPatchField(*this));
        }

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone(
            const DimensionedField<scalar, volMesh> &iF) const
        {
            return tmp<fvPatchScalarField>(
                new zeroFluxFvPatchField(*this, iF));
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
            return gradient();
        }

        //- Evaluate the patch field
        virtual void evaluate(
            const Pstream::commsTypes commsType = Pstream::commsTypes::blocking);

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

        // Member functions
        virtual scalarField &gradient()
        {
            return gradient_;
        };

        virtual const scalarField &gradient() const
        {
            return gradient_;
        };

        virtual void updateCoeffs();

        void initiateFluxNormVert(void);
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
