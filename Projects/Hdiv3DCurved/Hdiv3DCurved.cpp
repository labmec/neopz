#include "pzlog.h"
#include "tpzautopointer.h"
#include "TPZRefPatternTools.h"


#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "TPZRefPatternDataBase.h"
#include "TPZRefPatternTools.h"
#include "pzgeopoint.h"
#include "TPZGeoLinear.h"
#include "TPZGeoCube.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "pzgeoelside.h"
#include "tpzgeoblend.h"
#include "tpzarc3d.h"
#include "pzgeotetrahedra.h"
#include "pzgeoelrefless.h"
#include "tpzquadraticquad.h"
#include "tpzquadraticline.h"
#include "TPZQuadSphere.h"
#include "TPZTriangleSphere.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "pzfunction.h"
#include "tpzchangeel.h"

#include "pzpoisson3d.h"
#include "mixedpoisson.h"
#include "pzbndcond.h"

#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "pzskylstrmatrix.h"
#include "TPZParFrontStructMatrix.h"

#include "pzgmesh.h"
#include "TPZVTKGeoMesh.h"

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif


struct SimulationCase {
    bool            IsHdivQ;
    int             n_h_levels;
    int             n_p_levels;
    int             int_order;
    int             n_threads;
    std::string     mesh_type;
    std::string     domain_type;
    std::string     conv_summary;
    std::string     dump_folder;
    TPZStack<int>   omega_ids;
    TPZStack<int>   gamma_ids;
};

#define Solution1




void Analytic(const TPZVec<REAL> &x, TPZVec<STATE> &u,TPZFMatrix<STATE> &gradu);
void Solution(const TPZVec<REAL> &x, TPZVec<STATE> &f);
void f(const TPZVec<REAL> &p, TPZVec<STATE> &f, TPZFMatrix<STATE> &gradf);

TPZGeoMesh * GeomtricMesh(int ndiv, SimulationCase sim_data);

TPZGeoMesh * MakeSphereFromLinearQuadrilateralFaces(int ndiv, SimulationCase sim_data);
TPZGeoMesh * MakeSphereFromQuadrilateralFaces(int ndiv, SimulationCase sim_data);

TPZManVector<STATE,3> ParametricSphere(REAL radius,REAL phi,REAL theta);
void TransformToQuadratic(TPZGeoMesh *gmesh);
void RotateGeomesh(TPZGeoMesh *gmesh, REAL CounterClockwiseAngle, int &Axis);


TPZCompMesh * ComputationalMesh(TPZGeoMesh * geometry, int p, SimulationCase sim_data);

TPZCompMesh * PrimalMesh(TPZGeoMesh * geometry, int p, SimulationCase sim_data); //  Primal approximation
TPZCompMesh * DualMesh(TPZGeoMesh * geometry, int p, SimulationCase sim_data); // Dual approximation
TPZCompMesh * uMesh(TPZGeoMesh * geometry, int p, SimulationCase sim_data); // Hdiv space
TPZCompMesh * pMesh(TPZGeoMesh * geometry, int p, SimulationCase sim_data); // L2 space


TPZAnalysis * CreateAnalysis(TPZCompMesh * cmesh, SimulationCase sim_data);
void PosProcess(TPZCompMesh* cmesh, TPZAnalysis * an, std::string file, SimulationCase sim_data);

void ComputeCases(TPZStack<SimulationCase> cases);
void ComputeConvergenceRates(SimulationCase sim_data);


int main()
{
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    TPZStack<SimulationCase> simulations;
    
    // Primal Formulation over the solid sphere
    struct SimulationCase H1Case;
    H1Case.IsHdivQ = false;
    H1Case.n_h_levels = 2;
    H1Case.n_p_levels = 1;
    H1Case.int_order  = 5;
    H1Case.n_threads  = 0;
    H1Case.mesh_type = "quadratic";
    H1Case.domain_type = "sphere";
    H1Case.conv_summary = "convergence_summary";
    H1Case.dump_folder = "H1_sphere";
    H1Case.omega_ids.Push(1);     // Domain
    H1Case.gamma_ids.Push(-1);    // Gamma_D inner surface
    H1Case.gamma_ids.Push(-2);    // Gamma_D outer surface
    simulations.Push(H1Case);

    ComputeCases(simulations);
    
    return 0;
}

void ComputeCases(TPZStack<SimulationCase> cases){
    
    int n_cases = cases.size();
    for (int i = 0; i < n_cases; i++) {
        ComputeConvergenceRates(cases[i]);
    }
}

void ComputeConvergenceRates(SimulationCase sim_data){
    
    // Creating the directory
    std::string command = "mkdir " + sim_data.dump_folder;
    system(command.c_str());
    
    int n_h_levels = sim_data.n_h_levels;
    int n_p_levels = sim_data.n_p_levels;
    
    for (int p = 1; p <= n_p_levels; p++) {
        for (int h = 0; h <= n_h_levels; h++) {
            
            // Compute the geometry
            TPZGeoMesh * gmesh = GeomtricMesh(h, sim_data);
            
#ifdef PZDEBUG
            
            std::stringstream text_name;
            std::stringstream vtk_name;
            text_name   << sim_data.dump_folder << "/" "geo" << "_" << sim_data.mesh_type << "_" << sim_data.domain_type << "_" << "p" << p << "h" <<  h << ".txt";
            vtk_name    << sim_data.dump_folder << "/" "geo" << "_" << sim_data.mesh_type << "_" << sim_data.domain_type << "_" << "p" << p << "h" <<  h << ".vtk";
            ofstream textfile(text_name.str());
            gmesh->Print(textfile);
            
            std::ofstream vtkfile(vtk_name.str());
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);
        
#endif
            
            // Compute the geometry
            TPZCompMesh * cmesh = ComputationalMesh(gmesh, p, sim_data);
            
            // Create Analysis
            TPZAnalysis * analysis = CreateAnalysis(cmesh,sim_data);
            
#ifdef USING_BOOST
            boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
#endif
            analysis->Assemble();
            
#ifdef USING_BOOST
            boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
#endif
            analysis->Solve();
            
#ifdef USING_BOOST
            boost::posix_time::ptime t3 = boost::posix_time::microsec_clock::local_time();
            STATE assemble_time = boost::numeric_cast<double>((t2-t1).total_milliseconds());
            STATE solving_time  = boost::numeric_cast<double>((t3-t2).total_milliseconds());
#endif
            
            // PostProccessing
            std::stringstream sol_vtk_name;
            sol_vtk_name    << sim_data.dump_folder << "/" "sol" << "_" << sim_data.mesh_type << "_" << sim_data.domain_type << "_" << "p" << p << "h" <<  h << ".vtk";
            std::string file(sol_vtk_name.str());
            PosProcess(cmesh, analysis, file, sim_data);
            
            // compute the error
            
            
            // current output summary
            
        }
        
        // print convergence summary
        
    }
    
}


void Analytic(const TPZVec<REAL> &p, TPZVec<STATE> &u,TPZFMatrix<STATE> &gradu){
    DebugStop();
}

void Solution(const TPZVec<REAL> &p, TPZVec<STATE> &f){

    REAL x,y,z;
    x = p[0];
    y = p[1];
    z = p[2];
    
    REAL r = sqrt(x*x+y*y+z*z);
    REAL theta = acos(z/r);
    REAL phi = atan2(y,x);
    
#ifdef Solution1
    
    f[0] = r*r;
    
#endif
    
}

void f(const TPZVec<REAL> &p, TPZVec<STATE> &f, TPZFMatrix<STATE> &gradf){
    REAL x,y,z;
    x = p[0];
    y = p[1];
    z = p[2];
    
    REAL r = sqrt(x*x+y*y+z*z);
    REAL theta = acos(z/r);
    REAL phi = atan2(y,x);
    
#ifdef Solution1
    
    f[0] = -6.0;
    
#endif
}


TPZGeoMesh * GeomtricMesh(int ndiv, SimulationCase sim_data){
    
    TPZGeoMesh * geometry = NULL;
    
    if (sim_data.domain_type == "sphere") {
        
        if (sim_data.mesh_type == "linear") {
            geometry = MakeSphereFromLinearQuadrilateralFaces(ndiv, sim_data);
            return geometry;
        }
        
        if (sim_data.mesh_type == "quadratic") {
            geometry = MakeSphereFromQuadrilateralFaces(ndiv, sim_data);
            TransformToQuadratic(geometry);
            return geometry;
        }
        
        if (sim_data.mesh_type == "blended") {
            geometry = MakeSphereFromQuadrilateralFaces(ndiv, sim_data);
            return geometry;
        }
        
        std::cout << "Error:: unable to compute geometry for the given case = " << &sim_data << std::endl;
        DebugStop();
        
    }
    
    std::cout << "Error:: unable to compute geometry for the given case = " << &sim_data << std::endl;
    DebugStop();
    return geometry;
}

TPZCompMesh * ComputationalMesh(TPZGeoMesh * geometry, int p, SimulationCase sim_data){
    
    TPZCompMesh * mesh = NULL;
    
    if (sim_data.IsHdivQ) {
        
    }
    else
    {
        mesh = PrimalMesh(geometry, p, sim_data);
        return mesh;        
    }
    
    std::cout << "Error:: unable to compute the given case = " << &sim_data << std::endl;
    DebugStop();
    return mesh;
}

TPZAnalysis * CreateAnalysis(TPZCompMesh * cmesh, SimulationCase sim_data){
    
    TPZAnalysis * analysis = new TPZAnalysis(cmesh, true);
    
    // solver settings
    bool IsFrontQ      = sim_data.IsHdivQ;
    bool IsParsidoQ    = !sim_data.IsHdivQ;
    
    if (IsFrontQ) {
        TPZParFrontStructMatrix<TPZFrontSym<STATE> > matrix(cmesh);
        matrix.SetDecomposeType(ECholesky);
        matrix.SetNumThreads(sim_data.n_threads);
        analysis->SetStructuralMatrix(matrix);
        return analysis;
    }
    else{
        TPZSkylineStructMatrix matrix(cmesh);
        matrix.SetNumThreads(sim_data.n_threads);        
        TPZStepSolver<STATE> step;
        step.SetDirect(ECholesky);
        analysis->SetSolver(step);
        analysis->SetStructuralMatrix(matrix);
        return analysis;
    }
    
    if (IsParsidoQ) {
        TPZParFrontStructMatrix<TPZFrontSym<STATE> > matrix(cmesh);
        analysis->SetStructuralMatrix(matrix);
        return analysis;        
    }
    
   return analysis;
    
}

void PosProcess(TPZCompMesh* cmesh, TPZAnalysis  * an, std::string file, SimulationCase sim_data)
{
    int dim = 3;
//    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, cmesh);
    TPZStack<std::string,10> scalnames, vecnames;
    
    if (sim_data.IsHdivQ) {
        vecnames.Push("Flux");
        vecnames.Push("ExactFlux");
        scalnames.Push("Pressure");
        scalnames.Push("ExactPressure");
        scalnames.Push("Rhs");
        scalnames.Push("Divergence");
    }
    else{
        vecnames.Push("Flux");
        scalnames.Push("Solution");
    }

    int div = 0;
    an->DefineGraphMesh(dim,scalnames,vecnames,file);
    an->PostProcess(div,dim);
    
}

TPZCompMesh * PrimalMesh(TPZGeoMesh * geometry, int p, SimulationCase sim_data){
    
    int dimension = 3;
    int dirichlet = 0;
    int nvolumes = sim_data.omega_ids.size();
    int nboundaries = sim_data.gamma_ids.size();

    
#ifdef PZDEBUG
    if (nvolumes != 1) {
    std::cout << "Error:: unable to compute the given case = " << &sim_data << std::endl;
        DebugStop();
    }
#endif
    
    TPZCompMesh *cmesh = new TPZCompMesh(geometry);
    
    TPZMaterial * volume;
    TPZMaterial * face;
    
    TPZFMatrix<STATE> val1(dimension,dimension,0.),val2(dimension,1,0.);
    for (int iv = 0; iv < nvolumes ; iv++) {
        
        volume = new TPZMatPoisson3d(sim_data.omega_ids[iv],dimension);
        TPZDummyFunction<STATE> *dum = new TPZDummyFunction<STATE>(f);
        dum->SetPolynomialOrder(sim_data.int_order);
        TPZAutoPointer<TPZFunction<STATE> > rhs = dum;
        volume->SetForcingFunction(rhs);
        cmesh->InsertMaterialObject(volume);
        
        for (int ib = 0; ib < nboundaries; ib++) {
            face = volume->CreateBC(volume,sim_data.gamma_ids[ib],dirichlet,val1,val2);
            
            TPZDummyFunction<STATE> *analytic = new TPZDummyFunction<STATE>(Solution);
            analytic->SetPolynomialOrder(sim_data.int_order);
            TPZAutoPointer<TPZFunction<STATE> > solution = analytic;
            face->SetForcingFunction(solution);
            
            cmesh->InsertMaterialObject(face);
        }
        
    }
    
    cmesh->SetDefaultOrder(p);
    cmesh->SetAllCreateFunctionsContinuous();
    
    cmesh->AutoBuild();
    cmesh->ExpandSolution();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
#ifdef PZDEBUG
    std::stringstream file_name;
    file_name   << sim_data.dump_folder << "/" << "Primal_cmesh" << ".txt";
    std::ofstream sout(file_name.str());
    cmesh->Print(sout);
#endif
    
    return cmesh;
    
}

TPZCompMesh *DualMesh(TPZGeoMesh * geometry, int p, SimulationCase sim_data)
{
//    TPZCompMesh *fluxmesh = CMeshFlux(gmesh,porder,dim);
//    TPZCompMesh *pressuremesh = CMeshPressure(gmesh, porder, dim);
//    AdjustFluxPolynomialOrders(fluxmesh, hdivplusplus);
//    SetPressureOrders(fluxmesh, pressuremesh);
//    meshvec.resize(2);
//    meshvec[0] = fluxmesh;
//    meshvec[1] = pressuremesh;
//    TPZCompMesh *mixed = CMeshMixed(gmesh, meshvec);
//    TPZCompMeshTools::GroupElements(mixed);
//    TPZCompMeshTools::CreatedCondensedElements(mixed, true);
//    return mixed;
}

TPZGeoMesh * MakeSphereFromLinearQuadrilateralFaces(int ndiv, SimulationCase sim_data){
    
#ifdef PZDEBUG
    if (sim_data.omega_ids.size() != 1 || sim_data.gamma_ids.size() != 2) {
        std::cout << "Sphere:: Please pass materials ids for volume (1) and boundary domains = { (-1 -> inner surface ),  (-2 -> outer surface) }" << std::endl;
        DebugStop();
    }
#endif
    
    TPZGeoMesh * geomesh = new TPZGeoMesh;
    geomesh->SetDimension(3);
    int nl = 2;// Let it fixed
    int basenodes = 8;
    int nodes =  basenodes * (nl);
    REAL radius_o = 1.0;
    REAL radius_i = 0.25;
    
    REAL dr = (radius_o- radius_i)/REAL(nl-1);
    geomesh->SetMaxNodeId(nodes-1);
    geomesh->NodeVec().Resize(nodes);
    TPZManVector<TPZGeoNode,4> Node(nodes);
    
    TPZManVector<long,4> TopolQuad(4);
    TPZManVector<long,8> TopolCube(8);
    TPZManVector<REAL,3> coord(3,0.);
    TPZVec<REAL> xc(3,0.);
    REAL cphi = atan(sqrt(2.0));
    
    int nodeindex = 0;
    long id = 0;
    int matid = sim_data.omega_ids[0];
    
    TPZManVector< TPZManVector<REAL,3> , 8 > points(nodes,0.);
    for (int il = 0; il < nl; il++) {
        
        if (il==0) {
            matid = sim_data.gamma_ids[0];
        }
        
        if (il==nl-1) {
            matid = sim_data.gamma_ids[1];
        }
        
        REAL radius = radius_o - REAL(il)*dr;
        points[0].Resize(3, 0.0);
        points[0][0]=radius;
        points[0][1]=M_PI-cphi;
        points[0][2]=M_PI/4.0;
        
        points[1].Resize(3, 0.0);
        points[1][0]=radius;
        points[1][1]=cphi;
        points[1][2]=M_PI/4.0;
        
        points[2].Resize(3, 0.0);
        points[2][0]=radius;
        points[2][1]=cphi;
        points[2][2]=-M_PI/4.0;
        
        points[3].Resize(3, 0.0);
        points[3][0]=radius;
        points[3][1]=M_PI-cphi;
        points[3][2]=-M_PI/4.0;
        
        points[4].Resize(3, 0.0);
        points[4][0]=radius;
        points[4][1]=M_PI-cphi;
        points[4][2]=3.0*M_PI/4.0;
        
        points[5].Resize(3, 0.0);
        points[5][0]=radius;
        points[5][1]=cphi;
        points[5][2]=3.0*M_PI/4.0;
        
        points[6].Resize(3, 0.0);
        points[6][0]=radius;
        points[6][1]=cphi;
        points[6][2]=-3.0*M_PI/4.0;
        
        points[7].Resize(3, 0.0);
        points[7][0]=radius;
        points[7][1]=M_PI-cphi;
        points[7][2]=-3.0*M_PI/4.0;
        
        
        
        for (int i = 0; i < basenodes; i++) {
            coord = ParametricSphere(points[i][0],points[i][1],points[i][2]);
            geomesh->NodeVec()[nodeindex].SetCoord(coord);
            geomesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
            nodeindex++;
        }
        
        
        TopolQuad[0] = 0+il*basenodes;
        TopolQuad[1] = 1+il*basenodes;
        TopolQuad[2] = 2+il*basenodes;
        TopolQuad[3] = 3+il*basenodes;
        new TPZGeoElRefPattern<  pzgeom::TPZGeoQuad  > (id,TopolQuad,matid,*geomesh);
        id++;
        
        TopolQuad[0] = 4+il*basenodes;
        TopolQuad[1] = 5+il*basenodes;
        TopolQuad[2] = 6+il*basenodes;
        TopolQuad[3] = 7+il*basenodes;
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (id,TopolQuad,matid,*geomesh);
        id++;
        
        TopolQuad[0] = 0+il*basenodes;
        TopolQuad[1] = 4+il*basenodes;
        TopolQuad[2] = 5+il*basenodes;
        TopolQuad[3] = 1+il*basenodes;
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (id,TopolQuad,matid,*geomesh);
        id++;
        
        TopolQuad[0] = 3+il*basenodes;
        TopolQuad[1] = 7+il*basenodes;
        TopolQuad[2] = 6+il*basenodes;
        TopolQuad[3] = 2+il*basenodes;
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad  > (id,TopolQuad,matid,*geomesh);
        id++;
        
        TopolQuad[0] = 0+il*basenodes;
        TopolQuad[1] = 4+il*basenodes;
        TopolQuad[2] = 7+il*basenodes;
        TopolQuad[3] = 3+il*basenodes;
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (id,TopolQuad,matid,*geomesh);
        id++;
        
        TopolQuad[0] = 1+il*basenodes;
        TopolQuad[1] = 5+il*basenodes;
        TopolQuad[2] = 6+il*basenodes;
        TopolQuad[3] = 2+il*basenodes;
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (id,TopolQuad,matid,*geomesh);
        id++;
    }
    
    matid = 1;
    
    for (int il = 0; il < nl - 1 ; il++) {
        //      Inserting blend elements
        TopolCube[0] = 0+il*basenodes;
        TopolCube[1] = 1+il*basenodes;
        TopolCube[2] = 2+il*basenodes;
        TopolCube[3] = 3+il*basenodes;
        TopolCube[4] = 0+(il+1)*basenodes;
        TopolCube[5] = 1+(il+1)*basenodes;
        TopolCube[6] = 2+(il+1)*basenodes;
        TopolCube[7] = 3+(il+1)*basenodes;
        new TPZGeoElRefPattern< pzgeom::TPZGeoCube > (id,TopolCube,matid,*geomesh);
        id++;
        
        TopolCube[0] = 4+il*basenodes;
        TopolCube[1] = 5+il*basenodes;
        TopolCube[2] = 6+il*basenodes;
        TopolCube[3] = 7+il*basenodes;
        TopolCube[4] = 4+(il+1)*basenodes;
        TopolCube[5] = 5+(il+1)*basenodes;
        TopolCube[6] = 6+(il+1)*basenodes;
        TopolCube[7] = 7+(il+1)*basenodes;
        new TPZGeoElRefPattern< pzgeom::TPZGeoCube > (id,TopolCube,matid,*geomesh);
        id++;
        
        TopolCube[0] = 0+il*basenodes;
        TopolCube[1] = 4+il*basenodes;
        TopolCube[2] = 5+il*basenodes;
        TopolCube[3] = 1+il*basenodes;
        TopolCube[4] = 0+(il+1)*basenodes;
        TopolCube[5] = 4+(il+1)*basenodes;
        TopolCube[6] = 5+(il+1)*basenodes;
        TopolCube[7] = 1+(il+1)*basenodes;
        new TPZGeoElRefPattern< pzgeom::TPZGeoCube > (id,TopolCube,matid,*geomesh);
        id++;
        
        TopolCube[0] = 3+il*basenodes;
        TopolCube[1] = 7+il*basenodes;
        TopolCube[2] = 6+il*basenodes;
        TopolCube[3] = 2+il*basenodes;
        TopolCube[4] = 3+(il+1)*basenodes;
        TopolCube[5] = 7+(il+1)*basenodes;
        TopolCube[6] = 6+(il+1)*basenodes;
        TopolCube[7] = 2+(il+1)*basenodes;
        new TPZGeoElRefPattern< pzgeom::TPZGeoCube > (id,TopolCube,matid,*geomesh);
        id++;
        
        TopolCube[0] = 0+il*basenodes;
        TopolCube[1] = 4+il*basenodes;
        TopolCube[2] = 7+il*basenodes;
        TopolCube[3] = 3+il*basenodes;
        TopolCube[4] = 0+(il+1)*basenodes;
        TopolCube[5] = 4+(il+1)*basenodes;
        TopolCube[6] = 7+(il+1)*basenodes;
        TopolCube[7] = 3+(il+1)*basenodes;
        new TPZGeoElRefPattern< pzgeom::TPZGeoCube > (id,TopolCube,matid,*geomesh);
        id++;
        
        TopolCube[0] = 1+il*basenodes;
        TopolCube[1] = 5+il*basenodes;
        TopolCube[2] = 6+il*basenodes;
        TopolCube[3] = 2+il*basenodes;
        TopolCube[4] = 1+(il+1)*basenodes;
        TopolCube[5] = 5+(il+1)*basenodes;
        TopolCube[6] = 6+(il+1)*basenodes;
        TopolCube[7] = 2+(il+1)*basenodes;
        new TPZGeoElRefPattern< pzgeom::TPZGeoCube > (id,TopolCube,matid,*geomesh);
        id++;
        
    }
    
    geomesh->BuildConnectivity();
    
    TPZVec<TPZGeoEl *> sons;
    const int nref = ndiv;
    for (int iref = 0; iref < nref; iref++) {
        int nel = geomesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = geomesh->ElementVec()[iel];
            if(!gel->HasSubElement())
            {
                gel->Divide(sons);
            }
        }
    }
    
    int axis = 1;
    REAL angle = -45.0;
    RotateGeomesh(geomesh, angle, axis);
    return geomesh;
}

TPZGeoMesh * MakeSphereFromQuadrilateralFaces(int ndiv, SimulationCase sim_data)
{
    
#ifdef PZDEBUG
    if (sim_data.omega_ids.size() != 1 || sim_data.gamma_ids.size() != 2) {
        std::cout << "Sphere:: Please pass materials ids for volume (1) and boundary domains = { (-1 -> inner surface ),  (-2 -> outer surface) }" << std::endl;
        DebugStop();
    }
#endif
    
    TPZGeoMesh * geomesh = new TPZGeoMesh;
    geomesh->SetDimension(3);
    int nl = 2;// Let it fixed
    int basenodes = 8;
    int nodes =  basenodes * (nl);
    REAL radius_o = 1.0;
    REAL radius_i = 0.25;
    
    REAL dr = (radius_o- radius_i)/REAL(nl-1);
    geomesh->SetMaxNodeId(nodes-1);
    geomesh->NodeVec().Resize(nodes);
    TPZManVector<TPZGeoNode,4> Node(nodes);
    
    TPZManVector<long,4> TopolQuad(4);
    TPZManVector<long,8> TopolCube(8);
    TPZManVector<REAL,3> coord(3,0.);
    TPZVec<REAL> xc(3,0.);
    REAL cphi = atan(sqrt(2.0));
    
    int nodeindex = 0;
    long id = 0;
    int matid = sim_data.omega_ids[0];
    
    TPZManVector< TPZManVector<REAL,3> , 8 > points(nodes,0.);
    for (int il = 0; il < nl; il++) {
        
        if (il==0) {
            matid = sim_data.gamma_ids[0];
        }
        
        if (il==nl-1) {
            matid = sim_data.gamma_ids[1];
        }
        
        REAL radius = radius_o - REAL(il)*dr;
        points[0].Resize(3, 0.0);
        points[0][0]=radius;
        points[0][1]=M_PI-cphi;
        points[0][2]=M_PI/4.0;
        
        points[1].Resize(3, 0.0);
        points[1][0]=radius;
        points[1][1]=cphi;
        points[1][2]=M_PI/4.0;
        
        points[2].Resize(3, 0.0);
        points[2][0]=radius;
        points[2][1]=cphi;
        points[2][2]=-M_PI/4.0;
        
        points[3].Resize(3, 0.0);
        points[3][0]=radius;
        points[3][1]=M_PI-cphi;
        points[3][2]=-M_PI/4.0;
        
        points[4].Resize(3, 0.0);
        points[4][0]=radius;
        points[4][1]=M_PI-cphi;
        points[4][2]=3.0*M_PI/4.0;
        
        points[5].Resize(3, 0.0);
        points[5][0]=radius;
        points[5][1]=cphi;
        points[5][2]=3.0*M_PI/4.0;
        
        points[6].Resize(3, 0.0);
        points[6][0]=radius;
        points[6][1]=cphi;
        points[6][2]=-3.0*M_PI/4.0;
        
        points[7].Resize(3, 0.0);
        points[7][0]=radius;
        points[7][1]=M_PI-cphi;
        points[7][2]=-3.0*M_PI/4.0;
        
        
        
        for (int i = 0; i < basenodes; i++) {
            coord = ParametricSphere(points[i][0],points[i][1],points[i][2]);
            geomesh->NodeVec()[nodeindex].SetCoord(coord);
            geomesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
            nodeindex++;
        }
        
        
        TopolQuad[0] = 0+il*basenodes;
        TopolQuad[1] = 1+il*basenodes;
        TopolQuad[2] = 2+il*basenodes;
        TopolQuad[3] = 3+il*basenodes;
        TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > * quad1 = new TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > (id,TopolQuad,matid,*geomesh);
        quad1->Geom().SetData(radius, xc);
        id++;
        
        TopolQuad[0] = 4+il*basenodes;
        TopolQuad[1] = 5+il*basenodes;
        TopolQuad[2] = 6+il*basenodes;
        TopolQuad[3] = 7+il*basenodes;
        TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > * quad2 = new TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > (id,TopolQuad,matid,*geomesh);
        quad2->Geom().SetData(radius, xc);
        id++;
        
        TopolQuad[0] = 0+il*basenodes;
        TopolQuad[1] = 4+il*basenodes;
        TopolQuad[2] = 5+il*basenodes;
        TopolQuad[3] = 1+il*basenodes;
        TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > * quad3 = new TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > (id,TopolQuad,matid,*geomesh);
        quad3->Geom().SetData(radius, xc);
        id++;
        
        TopolQuad[0] = 3+il*basenodes;
        TopolQuad[1] = 7+il*basenodes;
        TopolQuad[2] = 6+il*basenodes;
        TopolQuad[3] = 2+il*basenodes;
        TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > * quad4 = new TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > (id,TopolQuad,matid,*geomesh);
        quad4->Geom().SetData(radius, xc);
        id++;
        
        TopolQuad[0] = 0+il*basenodes;
        TopolQuad[1] = 4+il*basenodes;
        TopolQuad[2] = 7+il*basenodes;
        TopolQuad[3] = 3+il*basenodes;
        TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > * quad5 = new TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > (id,TopolQuad,matid,*geomesh);
        quad5->Geom().SetData(radius, xc);
        id++;
        
        TopolQuad[0] = 1+il*basenodes;
        TopolQuad[1] = 5+il*basenodes;
        TopolQuad[2] = 6+il*basenodes;
        TopolQuad[3] = 2+il*basenodes;
        TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > * quad6 = new TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > (id,TopolQuad,matid,*geomesh);
        quad6->Geom().SetData(radius, xc);
        id++;
    }
    
    matid = 1;
    
    for (int il = 0; il < nl - 1 ; il++) {
        //      Inserting blend elements
        TopolCube[0] = 0+il*basenodes;
        TopolCube[1] = 1+il*basenodes;
        TopolCube[2] = 2+il*basenodes;
        TopolCube[3] = 3+il*basenodes;
        TopolCube[4] = 0+(il+1)*basenodes;
        TopolCube[5] = 1+(il+1)*basenodes;
        TopolCube[6] = 2+(il+1)*basenodes;
        TopolCube[7] = 3+(il+1)*basenodes;
        geomesh->CreateGeoBlendElement(ECube, TopolCube, matid, id);
        id++;
        
        TopolCube[0] = 4+il*basenodes;
        TopolCube[1] = 5+il*basenodes;
        TopolCube[2] = 6+il*basenodes;
        TopolCube[3] = 7+il*basenodes;
        TopolCube[4] = 4+(il+1)*basenodes;
        TopolCube[5] = 5+(il+1)*basenodes;
        TopolCube[6] = 6+(il+1)*basenodes;
        TopolCube[7] = 7+(il+1)*basenodes;
        geomesh->CreateGeoBlendElement(ECube, TopolCube, matid, id);
        id++;
        
        TopolCube[0] = 0+il*basenodes;
        TopolCube[1] = 4+il*basenodes;
        TopolCube[2] = 5+il*basenodes;
        TopolCube[3] = 1+il*basenodes;
        TopolCube[4] = 0+(il+1)*basenodes;
        TopolCube[5] = 4+(il+1)*basenodes;
        TopolCube[6] = 5+(il+1)*basenodes;
        TopolCube[7] = 1+(il+1)*basenodes;
        geomesh->CreateGeoBlendElement(ECube, TopolCube, matid, id);
        id++;
        
        TopolCube[0] = 3+il*basenodes;
        TopolCube[1] = 7+il*basenodes;
        TopolCube[2] = 6+il*basenodes;
        TopolCube[3] = 2+il*basenodes;
        TopolCube[4] = 3+(il+1)*basenodes;
        TopolCube[5] = 7+(il+1)*basenodes;
        TopolCube[6] = 6+(il+1)*basenodes;
        TopolCube[7] = 2+(il+1)*basenodes;
        geomesh->CreateGeoBlendElement(ECube, TopolCube, matid, id);
        id++;
        
        TopolCube[0] = 0+il*basenodes;
        TopolCube[1] = 4+il*basenodes;
        TopolCube[2] = 7+il*basenodes;
        TopolCube[3] = 3+il*basenodes;
        TopolCube[4] = 0+(il+1)*basenodes;
        TopolCube[5] = 4+(il+1)*basenodes;
        TopolCube[6] = 7+(il+1)*basenodes;
        TopolCube[7] = 3+(il+1)*basenodes;
        geomesh->CreateGeoBlendElement(ECube, TopolCube, matid, id);
        id++;
        
        TopolCube[0] = 1+il*basenodes;
        TopolCube[1] = 5+il*basenodes;
        TopolCube[2] = 6+il*basenodes;
        TopolCube[3] = 2+il*basenodes;
        TopolCube[4] = 1+(il+1)*basenodes;
        TopolCube[5] = 5+(il+1)*basenodes;
        TopolCube[6] = 6+(il+1)*basenodes;
        TopolCube[7] = 2+(il+1)*basenodes;
        geomesh->CreateGeoBlendElement(ECube, TopolCube, matid, id);
        id++;
        
    }
    
    geomesh->BuildConnectivity();
    
    TPZVec<TPZGeoEl *> sons;
    const int nref = ndiv;
    for (int iref = 0; iref < nref; iref++) {
        int nel = geomesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = geomesh->ElementVec()[iel];
            if(!gel->HasSubElement())
            {
                gel->Divide(sons);
            }
        }
    }
    
    int axis = 1;
    REAL angle = -45.0;
    RotateGeomesh(geomesh, angle, axis);

    return geomesh;
}

void TransformToQuadratic(TPZGeoMesh *gmesh)
{
    long nel = gmesh->NElements();
    for (long el=0; el<nel; el++)
    {
        TPZGeoEl *gel = gmesh->Element(el);
        if (!gel || gel->HasSubElement()) {
            continue;
        }
        TPZGeoEl *father = gel->Father();
        int whichsubel = gel->WhichSubel();
        gel = TPZChangeEl::ChangeToQuadratic(gmesh, el);
        if (whichsubel != -1) {
            father->SetSubElement(whichsubel, gel);
        }
    }
}


TPZManVector<STATE,3> ParametricSphere(REAL radius,REAL phi,REAL theta)
{
    TPZManVector<STATE,3> xcoor(3,0.0);
    xcoor[0] = radius * cos(theta) * sin(phi);
    xcoor[1] = radius * sin(theta) * sin(phi);
    xcoor[2] = radius * cos(phi) ;
    return xcoor;
}

void RotateGeomesh(TPZGeoMesh *gmesh, REAL CounterClockwiseAngle, int &Axis)
{
    REAL theta =  (M_PI/180.0)*CounterClockwiseAngle;
    // It represents a 3D rotation around the axis -> i.
    TPZFMatrix<STATE> RotationMatrix(3,3,0.0);
    
    switch (Axis) {
        case 1:
        {
            RotationMatrix(0,0) = 1.0;
            RotationMatrix(1,1) =   +cos(theta);
            RotationMatrix(1,2) =   -sin(theta);
            RotationMatrix(2,1) =   +sin(theta);
            RotationMatrix(2,2) =   +cos(theta);
        }
            break;
        case 2:
        {
            RotationMatrix(0,0) =   +cos(theta);
            RotationMatrix(0,2) =   +sin(theta);
            RotationMatrix(1,1) = 1.0;
            RotationMatrix(2,0) =   -sin(theta);
            RotationMatrix(2,2) =   +cos(theta);
        }
            break;
        case 3:
        {
            RotationMatrix(0,0) =   +cos(theta);
            RotationMatrix(0,1) =   -sin(theta);
            RotationMatrix(1,0) =   +sin(theta);
            RotationMatrix(1,1) =   +cos(theta);
            RotationMatrix(2,2) = 1.0;
        }
            break;
        default:
        {
            RotationMatrix(0,0) =   +cos(theta);
            RotationMatrix(0,1) =   -sin(theta);
            RotationMatrix(1,0) =   +sin(theta);
            RotationMatrix(1,1) =   +cos(theta);
            RotationMatrix(2,2) = 1.0;
        }
            break;
    }
    
    TPZVec<STATE> iCoords(3,0.0);
    TPZVec<STATE> iCoordsRotated(3,0.0);
    
    //RotationMatrix.Print("Rotation = ");
    
    int NumberofGeoNodes = gmesh->NNodes();
    for (int inode = 0; inode < NumberofGeoNodes; inode++)
    {
        TPZGeoNode GeoNode = gmesh->NodeVec()[inode];
        GeoNode.GetCoordinates(iCoords);
        // Apply rotation
        iCoordsRotated[0] = RotationMatrix(0,0)*iCoords[0]+RotationMatrix(0,1)*iCoords[1]+RotationMatrix(0,2)*iCoords[2];
        iCoordsRotated[1] = RotationMatrix(1,0)*iCoords[0]+RotationMatrix(1,1)*iCoords[1]+RotationMatrix(1,2)*iCoords[2];
        iCoordsRotated[2] = RotationMatrix(2,0)*iCoords[0]+RotationMatrix(2,1)*iCoords[1]+RotationMatrix(2,2)*iCoords[2];
        GeoNode.SetCoord(iCoordsRotated);
        gmesh->NodeVec()[inode] = GeoNode;
    }
    
}