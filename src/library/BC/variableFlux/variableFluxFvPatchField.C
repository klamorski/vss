#include "variableFluxFvPatchField.H"
#include "dictionary.H"
#include "fvCFD.H"
#include "gravityMeshObject.H"
#include "addToRunTimeSelectionTable.H"
#include "../../external/csv/csv.h"

//https://github.com/ben-strasser/fast-cpp-csv-parser
//https://medium.com/@ryan_forrester_/reading-csv-files-in-c-how-to-guide-35030eb378ad
//https://github.com/CD3/libInterpolate
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
    //- Construct from patch and internal field
    variableFluxFvPatchField::variableFluxFvPatchField(
        const fvPatch &p,
        const DimensionedField<scalar, volMesh> &iF)
        : zeroFluxFvPatchField(p, iF)//,
         // flux_(p.size(), pTraits<scalar>::zero)
    {
    }

    //- Construct from patch, internal field and dictionary
    variableFluxFvPatchField::variableFluxFvPatchField(
        const fvPatch &p,
        const DimensionedField<scalar, volMesh> &iF,
        const dictionary &dict)
        : zeroFluxFvPatchField(p, iF, dict),
         // flux_("flux", dict, p.size()),
          flux_data_csv_fname_(dict.get<string>("variable_flux_fname"))
    {
        string csv_fname = this->db().time().path() + "/simulation/" + flux_data_csv_fname_;
        io::CSVReader<2> in(csv_fname);
        in.read_header(io::ignore_extra_column, "time", "flux");
        double time; double flux;
        while(in.read_row(time, flux)){
            //Info << "time: " << time << ", flux: " << flux << nl;
            flux_data_time_.push_back(time);
            flux_data_flux_.push_back(flux);
        }
        Info<<"Loaded "<<flux_data_time_.size()<<" flux data records from "<<flux_data_csv_fname_<<"."<<nl;
        if ((this->db().time().endTime().value() > flux_data_time_.back()) || (this->db().time().startTime().value() < flux_data_time_[0]))
        {
            Info<<"The timespan of the flux data does not agree with the timespan of the simulation. Check the "<<flux_data_csv_fname_<<" file. Exiting now!"<<nl;
            std::exit(1);
        }
        interpolator_.setData(flux_data_time_,flux_data_flux_);
        this->initiateFluxNormVert();
        evaluate();
    }

    //- Construct by mapping the given variableFluxFvPatchField
    //  onto a new patch
    variableFluxFvPatchField::variableFluxFvPatchField(
        const variableFluxFvPatchField &ptf,
        const fvPatch &p,
        const DimensionedField<scalar, volMesh> &iF,
        const fvPatchFieldMapper &mapper)
        : zeroFluxFvPatchField(ptf, p, iF, mapper)//,
         // flux_(ptf.flux_, mapper)
    {
        this->initiateFluxNormVert();

        if (&iF && mapper.hasUnmapped())
        {
            WarningIn(
                "variableFluxFvPatchField::variableFluxFvPatchField\n"
                "(\n"
                "    const variableFluxFvPatchField&,\n"
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
    variableFluxFvPatchField::variableFluxFvPatchField(
        const variableFluxFvPatchField &ptf)
        : zeroFluxFvPatchField(ptf)//,
          //flux_(ptf.flux_)
    {
        this->initiateFluxNormVert();
    }

    //- Construct as copy setting internal field reference
    variableFluxFvPatchField::variableFluxFvPatchField(
        const variableFluxFvPatchField &ptf,
        const DimensionedField<scalar, volMesh> &iF)
        : zeroFluxFvPatchField(ptf, iF)//,
          //flux_(ptf.flux_)
    {
        this->initiateFluxNormVert();
    }

    // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
    //- Map (and resize as needed) from self given a mapping object
    void variableFluxFvPatchField::autoMap(
        const fvPatchFieldMapper &m)
    {
        fvPatchScalarField::autoMap(m);
        //flux().autoMap(m);

        this->initiateFluxNormVert();
    }

    //- Reverse map the given fvPatchField onto this fvPatchField
    void variableFluxFvPatchField::rmap(
        const fvPatchScalarField &ptf,
        const labelList &addr)
    {
        fvPatchScalarField::rmap(ptf, addr);

        const variableFluxFvPatchField &fgptf =
            refCast<const variableFluxFvPatchField>(ptf);

        //flux().rmap(fgptf.flux(), addr);

        this->initiateFluxNormVert();
    }

    void variableFluxFvPatchField::evaluate(const Pstream::commsTypes)
    {
        if (!this->updated())
        {
            updateCoeffs();
        }

        scalarField::operator=(
            patchInternalField() + calculateGradient() / patch().deltaCoeffs());

        // fvPatchScalarField::evaluate();
    }

    void variableFluxFvPatchField::updateCoeffs()
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

    tmp<scalarField> variableFluxFvPatchField::valueInternalCoeffs(
        const tmp<scalarField> &) const
    {
        return tmp<scalarField>(new scalarField(this->size(), pTraits<scalar>::one));
    }

    tmp<scalarField> variableFluxFvPatchField::valueBoundaryCoeffs(
        const tmp<scalarField> &) const
    {
        return calculateGradient() / patch().deltaCoeffs();
    }

    tmp<scalarField> variableFluxFvPatchField::gradientInternalCoeffs() const
    {
        return tmp<scalarField>(
            new scalarField(this->size(), pTraits<scalar>::zero));
    }

    tmp<scalarField> variableFluxFvPatchField::gradientBoundaryCoeffs() const
    {
        return calculateGradient();
    }

    void variableFluxFvPatchField::write(Ostream &os) const
    {
        zeroFluxFvPatchField::write(os);
        //flux().writeEntry("flux", os);
        os.writeEntry("variable_flux_fname", flux_data_csv_fname_);
    }

    tmp<scalarField> variableFluxFvPatchField::calculateGradient(void) const
    {
        // FIXME: Why these two K scalarField approaches gives different results?!
        // const scalarField &K = patch().lookupPatchField<volScalarField>("K");

        const scalarField &K = (this->db().objectRegistry::foundObject<volScalarField>("K") == true) ? static_cast<scalarField>(this->db().objectRegistry::lookupObject<volScalarField>("K").boundaryField()[this->patch().index()]) : Field<double>(this->size(), pTraits<double>::one);
        if (!this->db().objectRegistry::foundObject<volScalarField>("K"))
        {
            Info << "Warning! Dummy K field initialized to scalar(1)." << endl;
        }
        double val = interpolator_(this->db().time().value());
        //tmp<scalarField> flux = scalarField(this->size(), val) / K;
        
        //Info<<"FLux data - time; "<<this->db().time().value()<<"; flux ratio; "<<val<<" K.value(); "<<K<<" gradient; "<<(val/K)<<endl;
        //Info<<val/K<<endl;
        return (gradient() - val/K);
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam

namespace Foam
{

    makePatchTypeField(
        fvPatchScalarField,
        variableFluxFvPatchField);

}; // End namespace Foam

// ************************************************************************* //
