

#include <iostream>

// Geometry description
#include "TPZGmshReader.h"
#include "TPZVTKGeoMesh.h"

void PrintGeometry(TPZGeoMesh * geometry){
        //  Print Geometrical Base Mesh
        std::ofstream argument("geometry.txt");
        geometry->Print(argument);
        std::ofstream Dummyfile("geometry.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(geometry,Dummyfile, true);
}

int main()
{
    
    std::string dirname = PZSOURCEDIR;
    //  Reading mesh
    std::string mesh_file;
    mesh_file = dirname + "/Projects/dfn_poroelastic/gmsh/reservoir.msh";
    
    TPZGmshReader mesh_reader;
    REAL scale = 1.0;
    mesh_reader.SetDimensionlessL(scale);
    TPZGeoMesh * geometry = mesh_reader.GeometricGmshMesh(mesh_file);
    PrintGeometry(geometry);
    
    
}

