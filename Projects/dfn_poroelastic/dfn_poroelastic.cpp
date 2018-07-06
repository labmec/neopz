

#include <iostream>

// Geometry description
#include "TPZGmshReader.h"
#include "TPZVTKGeoMesh.h"

// Computational operators
#include "TNFRElastic.h"

// Print geometric mesh
void PrintGeometry(TPZGeoMesh * geometry);


int main()
{
    
    std::string dirname = PZSOURCEDIR;
    //  Reading mesh
    std::string mesh_file;
    mesh_file = dirname + "/Projects/dfn_poroelastic/gmsh/wellbore.msh";
    
    TPZGmshReader mesh_reader;
    REAL scale = 1.0;
    mesh_reader.SetDimensionlessL(scale);
    TPZGeoMesh * geometry = mesh_reader.GeometricGmshMesh(mesh_file);
    geometry->SetName("NFR geometry description");
    PrintGeometry(geometry);
    
    
    
    //  Construction for the elastic operator
    TNFRElastic A_e;
    
    // Auxiliary objects
    std::vector<REAL> tensor_values(6,0.0);
    std::vector<REAL> vector_values(3,0.0);
    
    // Defining boundary conditions the elastic operator
    std::vector<TNFRBoundaryDescription> boundary_data;
    
    REAL wb_p = -40.0;
    tensor_values[0] = wb_p;
    tensor_values[3] = wb_p;
    tensor_values[5] = wb_p;
    
    // Wellbore bc
    TNFRBoundaryDescription Wellbore_bc;
    Wellbore_bc.SetBCId(2);
    Wellbore_bc.SetBCType(2);
    Wellbore_bc.SetBCValues(tensor_values);
    
    REAL s_xx_0 = -50.0;
    REAL s_yy_0 = -60.0;
    REAL s_zz_0 = -10.0;
    tensor_values[0] = s_xx_0;
    tensor_values[3] = s_yy_0;
    tensor_values[5] = s_zz_0;
    
    // South bc
    TNFRBoundaryDescription S_bc;
    S_bc.SetBCId(3);
    S_bc.SetBCType(2);
    S_bc.SetBCValues(tensor_values);
    
    // East bc
    TNFRBoundaryDescription E_bc;
    E_bc.SetBCId(4);
    E_bc.SetBCType(2);
    E_bc.SetBCValues(tensor_values);
    
    // North bc
    TNFRBoundaryDescription N_bc;
    N_bc.SetBCId(5);
    N_bc.SetBCType(2);
    N_bc.SetBCValues(tensor_values);
    
    // West bc
    TNFRBoundaryDescription W_bc;
    W_bc.SetBCId(6);
    W_bc.SetBCType(2);
    W_bc.SetBCValues(tensor_values);
    
    // Fixed ux bc
    TNFRBoundaryDescription fixed_ux_bc;
    fixed_ux_bc.SetBCId(7);
    fixed_ux_bc.SetBCType(3);
    fixed_ux_bc.SetBCValues(vector_values);

    // Fixed uy bc
    TNFRBoundaryDescription fixed_uy_bc;
    fixed_uy_bc.SetBCId(8);
    fixed_uy_bc.SetBCType(4);
    fixed_uy_bc.SetBCValues(vector_values);
    
    boundary_data.push_back(Wellbore_bc);
    boundary_data.push_back(S_bc);
    boundary_data.push_back(E_bc);
    boundary_data.push_back(N_bc);
    boundary_data.push_back(W_bc);
    boundary_data.push_back(fixed_ux_bc);
    boundary_data.push_back(fixed_uy_bc);
    A_e.SetBoundaryData(boundary_data);
    A_e.SetNumberOfThreads(0);
    A_e.SetEnableBandwidthReduction(true);
    
    // Build the operator with polynomial order 1
    int p_order = 2;
    bool Ae_CoherenceQ = A_e.BuildOperator(geometry,p_order);
    

    // usage of A_e
    
    A_e.ExecuteASingleTimeStep();
    A_e.PostProcess();
    
}

void PrintGeometry(TPZGeoMesh * geometry){
    //  Print Geometrical Base Mesh
    std::ofstream argument("geometry.txt");
    geometry->Print(argument);
    std::ofstream Dummyfile("geometry.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(geometry,Dummyfile, true);
}
