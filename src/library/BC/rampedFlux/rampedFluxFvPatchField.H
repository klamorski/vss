#ifndef rampedFluxFvPatchField_H
#define rampedFluxFvPatchField_H

#include "fvPatchFields.H"
#include "fixedFluxFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
                     Class fixedFluxFvPatchField Declaration
    \*---------------------------------------------------------------------------*/
    class rampedFluxFvPatchField
        : public fixedFluxFvPatchField // TODO: Use zeroGradientFvPatchField as the parent class
    {
        // Private data
        scalar hMin_;       // ramping pressure head
        scalar fluxToK_;    // conversion factor from flux to hydraulic conductivity
        word rampFunction_; // ramping function [steep | linear]

    public:
        //- Runtime type information
        TypeName("rampedFlux");

        // Constructors
        //- Construct from patch and internal field
        rampedFluxFvPatchField(
            const fvPatch &,
            const DimensionedField<scalar, volMesh> &);

        //- Construct from patch, internal field and dictionary
        rampedFluxFvPatchField(
            const fvPatch &,
            const DimensionedField<scalar, volMesh> &,
            const dictionary &);

        //- Construct by mapping the given rampedFluxFvPatchField
        //  onto a new patch
        rampedFluxFvPatchField(
            const rampedFluxFvPatchField &,
            const fvPatch &,
            const DimensionedField<scalar, volMesh> &,
            const fvPatchFieldMapper &);

        //- Construct as copy
        rampedFluxFvPatchField(
            const rampedFluxFvPatchField &);

        //- Construct as copy setting internal field reference
        rampedFluxFvPatchField(
            const rampedFluxFvPatchField &,
            const DimensionedField<scalar, volMesh> &);

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>(
                new rampedFluxFvPatchField(*this));
        };

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone(
            const DimensionedField<scalar, volMesh> &iF) const
        {
            return tmp<fvPatchScalarField>(
                new rampedFluxFvPatchField(*this, iF));
        };

        void rmap(
            const fvPatchScalarField &ptf,
            const labelList &addr);

        tmp<scalarField> calculateGradient(void) const;

        virtual void write(Ostream &) const;
    };
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
