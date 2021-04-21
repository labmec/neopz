
#include <iostream>
#include <string>
#include <math.h>


#include "pzgmesh.h"
#include "pzcmesh.h"
#include "TPZVTKGeoMesh.h"
#include "pzcheckmesh.h"
#include "TPZGmshReader.h"
#include "pzintel.h"
#include "pzelast3d.h"
#include "pzbndcond.h"
#include "pzgeopoint.h"
#include "tpzgeoelrefpattern.h"
#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "pzskylstrmatrix.h"
#include "TPZSSpStructMatrix.h"

// Read the drill pipe section
TPZGeoMesh * Geometriy();

// Print the drill pipe section
void PrintGeometry(TPZGeoMesh * gmesh);

// Creating the computational mesh
TPZCompMesh * ComputationalMesh(TPZGeoMesh * geometry, int p_order);

//Creating the analysis
TPZAnalysis * CreateAnalysis(TPZCompMesh * cmesh);

// Post-processing solution
void PostProcess(TPZAnalysis  * an, std::string file);

int main()
{
    // Compute the geometry
    TPZGeoMesh * gmesh;
    gmesh = Geometriy();
    PrintGeometry(gmesh);
    
    int64_t p_order = 1;
    TPZCompMesh * cmesh = ComputationalMesh(gmesh, p_order);
    TPZAnalysis * analysis = CreateAnalysis(cmesh);
    analysis->Run();
    std::cout << "Number of equations = " << cmesh->NEquations() << std::endl;
    std::string file = "drill_pipe.vtk";
    PostProcess(analysis, file);
    
	return 0;
}

TPZGeoMesh * Geometriy(){

    TPZGeoMesh * geomesh = new TPZGeoMesh;
    std::string dirname = PZSOURCEDIR;
    std::string grid;
    grid = dirname + "/Projects/DrillPipe/mesh/drill_pipe.msh";
    
    TPZGmshReader Geometry;
    REAL s = 1.0;
    geomesh = Geometry.GeometricGmshMesh(grid);
    const std::string name("Drill pipe section");
    geomesh->SetName(name);
    
    TPZManVector<REAL,3> co(3,0.0);
    co[0] = 0.0;
    co[1] = -0.0635;
    co[2] = 0.0;
    TPZGeoNode * node = geomesh->FindNode(co);

    int64_t n_el = geomesh->NElements();
    // Point
    TPZManVector <int64_t,1> TopolPoint(1);
    int64_t zero_d_element_id = n_el;
    TopolPoint[0] = node->Id();
    int64_t matid = 4;
    new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> (zero_d_element_id, TopolPoint, matid, *geomesh);
    geomesh->BuildConnectivity();
    
    return geomesh;
}

void PrintGeometry(TPZGeoMesh * gmesh){
    std::stringstream text_name;
    std::stringstream vtk_name;
    text_name   << "geometry" << ".txt";
    vtk_name    << "geometry" << ".vtk";
    std::ofstream textfile(text_name.str().c_str());
    gmesh->Print(textfile);
    std::ofstream vtkfile(vtk_name.str().c_str());
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);
}

TPZCompMesh * ComputationalMesh(TPZGeoMesh * geometry, int p_order){
    
    int64_t dim = 3;
    TPZCompMesh * cmesh = new TPZCompMesh(geometry);
    cmesh->SetName("Drill pipe deflection");
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(p_order);
    cmesh->SetAllCreateFunctionsContinuous();
    
    int64_t steel_id = 1;
    
    // Acero
    STATE Eyoung  = 203.0e9;
    STATE poisson = 0.27;
    STATE rho     = 7860.0;
    
//     Hierro de fundicion
//    STATE Eyoung  = 80.0e9;
//    STATE poisson = 0.25;
//    STATE rho     = 6920.0;
    
    STATE g       = 9.81;
    TPZManVector<STATE,3> b(3,0.0);
    b[1] = - g * rho;
    
    STATE preStressXX = 0.0;
    STATE preStressYY = 0.0;
    STATE preStressZZ = 0.0;
    
    TPZMaterial * material = new TPZElasticity3D(steel_id, Eyoung,poisson, b, preStressXX, preStressYY, preStressZZ);
    cmesh->InsertMaterialObject(material);
    
    int64_t bc_id = 2;
    int64_t fixed_uz = 8;
    TPZFMatrix<STATE> val1(3,3,0.0);
    TPZFMatrix<STATE> val2(3,1,0.0);
    TPZMaterial * bc_bottom = material->CreateBC(material, bc_id, fixed_uz, val1, val2);
    cmesh->InsertMaterialObject(bc_bottom);
    
    int64_t bc_id_point = 4;
    int64_t fixed_u = 0;
    TPZMaterial * bc_encastre = material->CreateBC(material, bc_id_point, fixed_u, val1, val2);
    cmesh->InsertMaterialObject(bc_encastre);
    
    cmesh->AutoBuild();
    
    
    
    return cmesh;
}

TPZAnalysis * CreateAnalysis(TPZCompMesh * cmesh){
    
    bool UsePardisoQ = true;
    int64_t n_threads = 8;
    TPZAnalysis * analysis = new TPZAnalysis(cmesh, true);
    if (UsePardisoQ) {
        
        TPZSSpStructMatrix matrix(cmesh);
        matrix.SetNumThreads(n_threads);
        analysis->SetStructuralMatrix(matrix);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt);
        analysis->SetSolver(step);
    }
    else{
        
        TPZSkylineStructMatrix matrix(cmesh);
        matrix.SetNumThreads(n_threads);
        TPZStepSolver<STATE> step;
        step.SetDirect(ECholesky);
        analysis->SetSolver(step);
        analysis->SetStructuralMatrix(matrix);
    }
    
    return analysis;
    
}

void PostProcess(TPZAnalysis  * an, std::string file)
{
    int64_t dim = 3;
    int64_t div = 0;
    
    TPZStack<std::string,10> scalnames, vecnames;
    vecnames.Push("Displacement");
    scalnames.Push("StressX");
    scalnames.Push("StressY");
    scalnames.Push("StressZ");
    an->DefineGraphMesh(dim,scalnames,vecnames,file);
    an->PostProcess(div,dim);
    
}
