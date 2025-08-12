#include "solverFunctions.h"

void debugShp(const Foam::fvMesh &mesh, const Foam::Time &runTime, retentionDataInterface* shpModel) {
    Info<<"Saving base SHP estimations."<<endl;
    scalar n_cells = 1e4;
    const scalar h_min = -10;
    const scalar h_max = 12;
    const scalar h_step = (h_max - h_min) / n_cells; 

    std::vector<double> h_debug;
    std::vector<double> theta_debug;
    std::vector<double> kh_debug;
    std::vector<double> cv_debug;

    // Create a standalone volScalarField with a default value of 0
    volScalarField h_debug_field (
        IOobject(
            "h_debug_field",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0, 1, 0, 0, 0, 0, 0), 0.0)
    );

    const label n_cells_mesh = mesh.nCells();
    Info<<"n_cells_mesh: "<<n_cells_mesh<<endl;
    // Initialize the field with custom values
    label i_vector=0; 
    label no_progess=20;
    label l_last = -1;
    for (label round = 0; round <= n_cells/n_cells_mesh; ++round)
    {
        label l = ceil(round/((n_cells/n_cells_mesh)/no_progess));
        if (l > l_last) {
            l_last = l;
            cout<<".";
        }
    }
    cout<<nl;
    l_last = -1;
    for (label round = 0; round <= n_cells/n_cells_mesh; ++round)
    {
        label l = ceil(round/((n_cells/n_cells_mesh)/no_progess));
        if (l > l_last) {
            l_last = l;
            cout<<"=";
        }
        

        for (label celli = 0; celli < n_cells_mesh; ++celli)
        {
            h_debug_field[celli] = -pow(10, h_min + h_step * (round * n_cells_mesh + celli));
        }
        Info<<"a"<<nl;
        volScalarField theta_debug_field = shpModel->Theta(h_debug_field);
        Info<<"c"<<nl;
        volScalarField inv_ret_kh_field = shpModel->Kh(h_debug_field);
        Info<<"d"<<nl;
        volScalarField inv_ret_cv_field = shpModel->Cv(h_debug_field);
        Info<<"e"<<nl;

        for (label celli = 0; celli < n_cells_mesh; ++celli)
        {
            h_debug.push_back(h_debug_field[celli]);
            theta_debug.push_back(theta_debug_field[celli]);
            kh_debug.push_back(inv_ret_kh_field[celli]);
            cv_debug.push_back(inv_ret_cv_field[celli]);
            i_vector++;
            if (i_vector >= n_cells) {
                break;
            }
        }
    }
    
    cout<<nl;

    int vsize = h_debug.size();
    OFstream debugFile("debug_shp.csv");
    debugFile<<"h;abs_h;theta;kh;cv"<<endl;
    for (int n=0; n<vsize; n++)
    {
        debugFile<<h_debug[n]<<";"<<abs(h_debug[n])<<";"<<theta_debug[n]<<";"<<kh_debug[n]<<";"<<cv_debug[n]<<endl;
    }
    Info<<"Debug file created: debug_shp.txt"<<endl;
}


