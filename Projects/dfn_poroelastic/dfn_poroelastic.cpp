

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
    PrintGeometry(geometry);
    
    //  Construction for the elastic operator
    TNFRElastic A_e;
    int p_order = 1;
    A_e.BuildOperator(geometry,p_order);
    
    
}

void PrintGeometry(TPZGeoMesh * geometry){
    //  Print Geometrical Base Mesh
    std::ofstream argument("geometry.txt");
    geometry->Print(argument);
    std::ofstream Dummyfile("geometry.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(geometry,Dummyfile, true);
}
