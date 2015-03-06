//
//  LaplaceInSphere.cpp
//  PZ
//
//  Created by Douglas Castro on 1/29/15.
//
//

#include "LaplaceInSphere.h"
#include "tools.h"

LaplaceInSphere::LaplaceInSphere(int ordemP, int ndiv, std::map<REAL, REAL> &fDebugMapL2, std::map<REAL, REAL> &fDebugMapHdiv)
{
    
    std::cout<< " INICIO(CASCA ESFERA) - grau  polinomio " << ordemP << " numero de divisoes " << ndiv << std::endl;
    std::cout<< " Dimensao == " << fDim << std::endl;
    
    TPZGeoMesh *gmesh = GMeshSphericalShell(ndiv);
    
    gmesh->SetDimension(fDim);
    ofstream argm("gmesh2d-Esfera.txt");
    gmesh->Print(argm);
    
    TPZCompMesh *cmesh2 = CMeshPressure(gmesh, ordemP, fDim);
    TPZCompMesh *cmesh1 = CMeshFlux(gmesh, ordemP, fDim);
    
    // Um teste para a solucao via H1, sem hdiv
    if (isH1) {
        TPZCompMesh *cmeshH1 = CMeshH1(gmesh, ordemP, fDim);
        TPZAnalysis anh1(cmeshH1, true);
        
        tools::SolveSyst(anh1, cmeshH1);
        
        stringstream refh1,grauh1;
        grauh1 << ordemP;
        refh1 << ndiv;
        string strgh1 = grauh1.str();
        string strrh1 = refh1.str();
        std::string plotnameh1("OurSolutionH1");
        std::string Grauh1("P");
        std::string Refh1("H");
        std::string VTKh1(".vtk");
        std::string plotDatah1;
        plotDatah1 = plotnameh1+Grauh1+strgh1+Refh1+strrh1+VTKh1;
        std::string plotfileh1(plotDatah1);
        
        tools::PosProcess(anh1, plotfileh1,fDim);
        
        return;
    }
    // exit
    
    ofstream arg1("cmeshflux.txt");
    cmesh1->Print(arg1);
    
    ofstream arg2("cmeshpressure.txt");
    cmesh2->Print(arg2);
    
    ofstream arg4("gmesh2.txt");
    gmesh->Print(arg4);
    
    
    //malha multifisica
    TPZVec<TPZCompMesh *> meshvec(2);
    meshvec[0] = cmesh1;
    meshvec[1] = cmesh2;
    
    TPZCompMesh * mphysics = CMeshMixed(gmesh,meshvec);
    
    //TestMesh(mphysics);
    ofstream arg5("cmeshmultiphysics.txt");
    mphysics->Print(arg5);
    
    TPZAnalysis an(mphysics, true);
    
    tools::SolveSyst(an, mphysics);
    
    stringstream ref,grau;
    grau << ordemP;
    ref << ndiv;
    string strg = grau.str();
    string strr = ref.str();
    std::string plotname("OurSolutionMetaEsfera");
    std::string Grau("P");
    std::string Ref("H");
    std::string VTK(".vtk");
    std::string plotData;
    plotData = plotname+Grau+strg+Ref+strr+VTK;
    std::string plotfile(plotData);
    
    tools::PosProcessMultphysics(meshvec,  mphysics, an, plotfile,fDim);
    
    //Calculo do erro
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
    TPZVec<REAL> erros;
    
    std::cout << "Postprocessed\n";
    
    stringstream ss;
    ss << ordemP;
    string str = ss.str();
    
    std::cout<< " grau  polinomio " << ordemP << " numero de divisoes " << ndiv << std::endl;
    std::string filename("InputDataMetaEsfera");
    std::string L2("L2.txt");
    std::string Hdiv("Hdiv.txt");
    std::string HdivData,L2Data;
    HdivData = filename+str+Hdiv;
    L2Data = filename+str+L2;
    
    ErrorHDiv(cmesh1, ordemP, ndiv, fDebugMapL2, fDebugMapHdiv);
    
    ErrorL2(cmesh2, ordemP, ndiv, fDebugMapL2, fDebugMapHdiv);
    
    tools::PrintDebugMapForMathematica(HdivData, L2Data, fDebugMapL2, fDebugMapHdiv);
    
    std::cout<< " FIM (ESFERA) - grau  polinomio " << ordemP << " numero de divisoes " << ndiv << std::endl;


}

LaplaceInSphere::~LaplaceInSphere()
{
    
}

//int LaplaceInSphere::getDimension()
//{
//    return fDim;
//}

TPZGeoMesh *LaplaceInSphere::GMeshSphericalShell(int ndiv)
{
    ///INSTANCIACAO DA MALHA GEOMETRICA
    TPZGeoMesh * geomesh = new TPZGeoMesh;
    
    /// Materiais
    long materialId = fmatId;
    long arc1 = fbc1; // -1;
    long arc2 = fbc2; // -2;
    long arc3 = fbc3; // -3;
    long arc4 = fbc4; // -4;
    
    int nnodes = 9;//quantidade de nos da malha geometrica
    geomesh->NodeVec().Resize(nnodes);
    geomesh->SetDimension(fDim);
    
    
    
    ///INICIALIZACAO DA MALHA GEOMETRICA PELA INSTANCIACAO E INICIALIZACAO DOS NOS
    TPZGeoNode node;
    TPZVec<REAL> coord(3,0.);
    const REAL r = 2.;
    TPZManVector<REAL,3> xc(3,0.);
    xc[0] = 0.;
    xc[1] = 0.;
    xc[2] = 0.;
    
    int id = 0;
    //no 0
    coord[0] = xc[0] + r;
    coord[1] = xc[1] + 0.;
    coord[2] = xc[2] + 0.;
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 1
    coord[0] = xc[0] + 0.;
    coord[1] = xc[1] + r;
    coord[2] = xc[2] + 0.;
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 2
    coord[0] = xc[0] + -r;
    coord[1] = xc[1] + 0.;
    coord[2] = xc[2] + 0.;
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 3
    coord[0] = xc[0] + 0.;
    coord[1] = xc[1] + -r;
    coord[2] = xc[2] + 0.;
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 4
    coord[0] = xc[0] + 0.;
    coord[1] = xc[1] + 0.;
    coord[2] = xc[2] + r;
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 5
    coord[0] = xc[0] + r*sqrt(2.)/2.;
    coord[1] = xc[1] + r*sqrt(2.)/2.;
    coord[2] = xc[2] + 0.;
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 6
    coord[0] = xc[0] + -r*sqrt(2.)/2.;
    coord[1] = xc[1] + r*sqrt(2.)/2.;
    coord[2] = xc[2] + 0.;
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 7
    coord[0] = xc[0] + -r*sqrt(2.)/2.;
    coord[1] = xc[1] + -r*sqrt(2.)/2.;
    coord[2] = xc[2] + 0.;
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 8
    coord[0] = xc[0] + r*sqrt(2.)/2.;
    coord[1] = xc[1] + -r*sqrt(2.)/2.;
    coord[2] = xc[2] + 0.;
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    
    int elementid = 0;
    // Using triangle to sphere special map
    TPZVec<long> topology(3);
    
    // El 0
    topology[0] = 0; 
    topology[1] = 1; 
    topology[2] = 4; 
    
    TPZGeoElRefPattern< pzgeom::TPZTriangleSphere > * SphereEighth1 =
    new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere > (elementid, topology,materialId,*geomesh);
    SphereEighth1->Geom().SetData(r,xc);
    elementid++;
    
    // El 1
    topology[0] = 1; 
    topology[1] = 2; 
    topology[2] = 4; 
    
    TPZGeoElRefPattern< pzgeom::TPZTriangleSphere > * SphereEighth2 =
    new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere > (elementid, topology,materialId,*geomesh);
    SphereEighth2->Geom().SetData(r,xc);
    elementid++;//
    
    // El 2
    topology[0] = 2; 
    topology[1] = 3; 
    topology[2] = 4; 
    
    TPZGeoElRefPattern< pzgeom::TPZTriangleSphere > * SphereEighth3 =
    new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere > (elementid, topology,materialId,*geomesh);
    SphereEighth3->Geom().SetData(r,xc);
    elementid++;//
    
    // El 3
    topology[0] = 3; 
    topology[1] = 0; 
    topology[2] = 4; 
    
    TPZGeoElRefPattern< pzgeom::TPZTriangleSphere > * SphereEighth4 =
    new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere > (elementid, topology,materialId,*geomesh);
    SphereEighth4->Geom().SetData(r,xc);
    elementid++;//
    
    // El linha
    // Definition of Arc coordenates
    topology.resize(3);
    // Create Geometrical Arc #1
    topology[0] = 0;
    topology[1] = 1;
    topology[2] = 5;
    new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,topology, arc3, *geomesh);
    elementid++;
    
    // Create Geometrical Arc #2
    topology[0] = 1;
    topology[1] = 2;
    topology[2] = 6;
    new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,topology, arc4, *geomesh);
    elementid++;
    
    // Create Geometrical Arc #3
    topology[0] = 2;
    topology[1] = 3;
    topology[2] = 7;
    new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,topology, arc1, *geomesh);
    elementid++;
    
    // Create Geometrical Arc #4
    topology[0] = 3;
    topology[1] = 0;
    topology[2] = 8;
    new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,topology, arc2, *geomesh);
    
    //CONCLUINDO A CONSTRUCAO DA MALHA GEOMETRICA
    geomesh->BuildConnectivity();
    
    
    TPZVec<TPZGeoEl *> sons;
    const int nref = 1;
    for (int iref = 0; iref < nref; iref++) {
        int nel = geomesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = geomesh->ElementVec()[iel];
            if(!gel->HasSubElement()) gel->Divide(sons);
        }
    }
    
    std::ofstream outfile("malhaCascaEsfera.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(geomesh, outfile, true);
    
    return geomesh;
    
   
}


void LaplaceInSphere::SolExata(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux){
    
    solp.resize(1);
    flux.Resize(3, 1);
    
    solp[0]=0.;
    flux(0,0)=0.0;
    flux(1,0)=0.0;
    flux(2,0)=0.0;
    
    REAL x = pt[0];
    REAL y = pt[1];
    REAL z = pt[2];
    
    
    int dim = 2; //getDimension();

    // tensor de permutacao
    TPZFNMatrix<2,REAL> TP(dim,dim,0.0);
    TPZFNMatrix<2,REAL> InvTP(dim,dim,0.0);
    
    
    // Hard coded
    for (int id = 0; id < dim; id++){
        TP(id,id) = 1.0;
        InvTP(id,id) = 1.0;
    }
    
    // sobre a casca da esfera -- conferir o raio aqui usado com o da malha geometrica

    REAL r = sqrt(x*x+y*y+z*z);
    REAL theta = atan2(sqrt(x*x+y*y),z);
    REAL phi = atan2(y,x);
    REAL cot = 1.0/tan(theta);
    
    solp[0] = ( 8.0 - 4.0*M_PI*(cot) + 8.0*theta*(cot) )/(M_PI*M_PI*r*r);
    
    flux(0,0) = (4.0*(M_PI - 2.0*theta)*(-TP(0,2)*sin(theta) + cos(theta)* (TP(0,0)*cos(phi) + TP(0,1) * sin(phi))))/(M_PI*M_PI*r);
    //4.0*(Pi - 2.0*theta)*cos(theta)*cos(phi)/(Pi*Pi*r);
    
    flux(1,0) = (4.0*(M_PI - 2.0*theta)*(-TP(1,2)*sin(theta) + cos(theta)* (TP(1,0)*cos(phi) + TP(1,1) * sin(phi))))/(M_PI*M_PI*r);
    //4.0*(Pi - 2.0*theta)*cos(theta)*sin(phi)/(Pi*Pi*r);
    
    flux(2,0) = (4.0*(M_PI - 2.0*theta)*(-TP(2,2)*sin(theta) + cos(theta)* (TP(2,0)*cos(phi) + TP(2,1) * sin(phi))))/(M_PI*M_PI*r);
    //-4.0*(Pi - 2.0*theta)*sin(theta)/(Pi*Pi*r);

    
//------------
    

    
    
}

void LaplaceInSphere::Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &ff){
  
    REAL x,y,z;

    x = pt[0];
    y = pt[1];
    z = pt[2];

    REAL r = sqrt(x*x+y*y+z*z);
    REAL theta = atan2(sqrt(x*x+y*y),z);
    REAL phi = atan2(y,x);
    REAL cot = 1.0/tan(theta);
    ff[0] = ( 8.0 - 4.0*M_PI*(cot) + 8.0*theta*(cot) )/(M_PI*M_PI*r*r);


}

void LaplaceInSphere::ForcingH1(const TPZVec<REAL> &pt, TPZVec<STATE> &ff, TPZFMatrix<STATE> &flux)
{
    flux.Resize(3, 1);
    REAL x,y,z;
    
    x = pt[0];
    y = pt[1];
    z = pt[2];
    
    REAL r = sqrt(x*x+y*y+z*z);
    REAL theta = atan2(sqrt(x*x+y*y),z);
    REAL phi = atan2(y,x);
    REAL cot = 1.0/tan(theta);
    
    int dim = 2; //getDimension();
    
    // tensor de permutacao
    TPZFNMatrix<2,REAL> TP(dim,dim,0.0);
    TPZFNMatrix<2,REAL> InvTP(dim,dim,0.0);
    
    
    // Hard coded
    for (int id = 0; id < dim; id++){
        TP(id,id) = 1.0;
        InvTP(id,id) = 1.0;
    }
    
    ff[0] = ( 8.0 - 4.0*M_PI*(cot) + 8.0*theta*(cot) )/(M_PI*M_PI*r*r);
    
    flux(0,0) = (4.0*(M_PI - 2.0*theta)*(-TP(0,2)*sin(theta) + cos(theta)* (TP(0,0)*cos(phi) + TP(0,1) * sin(phi))))/(M_PI*M_PI*r);
    //4.0*(Pi - 2.0*theta)*cos(theta)*cos(phi)/(Pi*Pi*r);
    
    flux(1,0) = (4.0*(M_PI - 2.0*theta)*(-TP(1,2)*sin(theta) + cos(theta)* (TP(1,0)*cos(phi) + TP(1,1) * sin(phi))))/(M_PI*M_PI*r);
    //4.0*(Pi - 2.0*theta)*cos(theta)*sin(phi)/(Pi*Pi*r);
    
    flux(2,0) = (4.0*(M_PI - 2.0*theta)*(-TP(2,2)*sin(theta) + cos(theta)* (TP(2,0)*cos(phi) + TP(2,1) * sin(phi))))/(M_PI*M_PI*r);
    //-4.0*(Pi - 2.0*theta)*sin(theta)/(Pi*Pi*r);
    
}

void LaplaceInSphere::ForcingBC0D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){
    
    REAL x,y,z;
    
    x = pt[0];
    y = pt[1];
    z = pt[2];
    
    REAL r = sqrt(x*x+y*y+z*z);
    REAL theta = atan2(sqrt(x*x+y*y),z);
    REAL phi = atan2(y,x);
    
    solp[0] = ((M_PI-2.0*theta)*(M_PI-2.0*theta))/(M_PI*M_PI);
    
    
    
}

void LaplaceInSphere::ForcingBC1D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){
    
    REAL x,y,z;
    
    x = pt[0];
    y = pt[1];
    z = pt[2];
    
    REAL r = sqrt(x*x+y*y+z*z);
    REAL theta = atan2(sqrt(x*x+y*y),z);
    REAL phi = atan2(y,x);
    
    solp[0] = ((M_PI-2.0*theta)*(M_PI-2.0*theta))/(M_PI*M_PI);
    
    
}

void LaplaceInSphere::ForcingBC2D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){
    
    REAL x,y,z;
    
    x = pt[0];
    y = pt[1];
    z = pt[2];
    
    REAL r = sqrt(x*x+y*y+z*z);
    REAL theta = atan2(sqrt(x*x+y*y),z);
    REAL phi = atan2(y,x);
    
    solp[0] = ((M_PI-2.0*theta)*(M_PI-2.0*theta))/(M_PI*M_PI);
    
}

void LaplaceInSphere::ForcingBC3D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){
    
    REAL x,y,z;
    
    x = pt[0];
    y = pt[1];
    z = pt[2];
    
    REAL r = sqrt(x*x+y*y+z*z);
    REAL theta = atan2(sqrt(x*x+y*y),z);
    REAL phi = atan2(y,x);
    
    solp[0] = ((M_PI-2.0*theta)*(M_PI-2.0*theta))/(M_PI*M_PI);

    
}

void LaplaceInSphere::ForcingBC4D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){

    REAL x,y,z;
    
    x = pt[0];
    y = pt[1];
    z = pt[2];
    
    REAL r = sqrt(x*x+y*y+z*z);
    REAL theta = atan2(sqrt(x*x+y*y),z);
    REAL phi = atan2(y,x);
    
    solp[0] = ((M_PI-2.0*theta)*(M_PI-2.0*theta))/(M_PI*M_PI);

    
}

void LaplaceInSphere::ForcingBC5D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){

    REAL x,y,z;

    x = pt[0];
    y = pt[1];
    z = pt[2];

    REAL r = sqrt(x*x+y*y+z*z);
    REAL theta = atan2(sqrt(x*x+y*y),z);
    REAL phi = atan2(y,x);

    solp[0] = ((M_PI-2.0*theta)*(M_PI-2.0*theta))/(M_PI*M_PI);

    
}


void LaplaceInSphere::ForcingBC0N(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
    DebugStop();
}

void LaplaceInSphere::ForcingBC1N(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
    DebugStop();
    
}

void LaplaceInSphere::ForcingBC2N(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    DebugStop();
}

void LaplaceInSphere::ForcingBC3N(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
    DebugStop();
    
}

void LaplaceInSphere::ForcingBC4N(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
    DebugStop();
}

void LaplaceInSphere::ForcingBC5N(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
    DebugStop();
}

TPZCompMesh *LaplaceInSphere::CMeshH1(TPZGeoMesh *gmesh, int pOrder, int dim)
{
    /// criar materiais
    TPZMatPoisson3d *material = new TPZMatPoisson3d(fmatId,dim);
    TPZMaterial * mat(material);
    material->NStateVariables();
    
    //    //solucao exata
    TPZAutoPointer<TPZFunction<STATE> > solexata;
    solexata = new TPZDummyFunction<STATE>(SolExata);
    material->SetForcingFunctionExact(solexata);
    
    //funcao do lado direito da equacao do problema
    //    TPZAutoPointer<TPZFunction<STATE> > forcef;
    //    forcef = new TPZDummyFunction<STATE>(Forcing);
    //    material->SetForcingFunction(forcef);
    
    //funcao do lado direito da equacao do problema
    TPZDummyFunction<STATE> *dum = new TPZDummyFunction<STATE>(ForcingH1);
    TPZAutoPointer<TPZFunction<STATE> > forcef;
    dum->SetPolynomialOrder(20);
    forcef = dum;
    material->SetForcingFunction(forcef);
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(pOrder);
    
    ///Inserir condicao de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZMaterial * BCond0;
    TPZMaterial * BCond1;
    TPZMaterial * BCond2;
    TPZMaterial * BCond3;
    TPZMaterial * BCond4;
    TPZMaterial * BCond5;
    
    if( dim == 3 ) { BCond0 = material->CreateBC(mat, fbc0,fdirichlet, val1, val2);}
    BCond1 = material->CreateBC(mat, fbc1,fdirichlet, val1, val2);
    BCond2 = material->CreateBC(mat, fbc2,fdirichlet, val1, val2);
    BCond3 = material->CreateBC(mat, fbc3,fdirichlet, val1, val2);
    BCond4 = material->CreateBC(mat, fbc4,fdirichlet, val1, val2);
    if( dim == 3 ) { BCond5 = material->CreateBC(mat, fbc5,fdirichlet, val1, val2);}
    
    cmesh->InsertMaterialObject(mat);
    
    cmesh->SetAllCreateFunctionsContinuous();
    if( dim == 3 ) { cmesh->InsertMaterialObject(BCond0); }
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    if( dim == 3 ) { cmesh->InsertMaterialObject(BCond5); }
    
    //Ajuste da estrutura de dados computacional
    cmesh->AutoBuild();
    
    return cmesh;
    
}


TPZCompMesh *LaplaceInSphere::CMeshFlux(TPZGeoMesh *gmesh, int pOrder, int dim)
{
    /// criar materiais
    //TPZMatPoisson3d *material = new TPZMatPoisson3d(matId,dim);
    TPZMatPoisson3d *material = new TPZMatPoisson3d(fmatId,dim);
    TPZMaterial * mat(material);
    material->NStateVariables();
    
    //  TPZAutoPointer<TPZFunction<STATE> > force1 = new TPZDummyFunction<STATE>(Forcing1);
    //	material->SetForcingFunction(force1);
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    
    
    ///Inserir condicao de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    TPZMaterial * BCond0;
    TPZMaterial * BCond1;
    TPZMaterial * BCond2;
    TPZMaterial * BCond3;
    TPZMaterial * BCond4;
    TPZMaterial * BCond5;
    
    if( dim == 3 ) { BCond0 = material->CreateBC(mat, fbc0,fdirichlet, val1, val2);}
    BCond1 = material->CreateBC(mat, fbc1,fdirichlet, val1, val2);
    BCond2 = material->CreateBC(mat, fbc2,fdirichlet, val1, val2);
    BCond3 = material->CreateBC(mat, fbc3,fdirichlet, val1, val2);
    BCond4 = material->CreateBC(mat, fbc4,fdirichlet, val1, val2);
    if( dim == 3 ) { BCond5 = material->CreateBC(mat, fbc5,fdirichlet, val1, val2);}
    
    cmesh->InsertMaterialObject(mat);
    
    cmesh->SetDimModel(dim);
    
    cmesh->SetAllCreateFunctionsHDiv();
    
    
    if( dim == 3 ) { cmesh->InsertMaterialObject(BCond0); }
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    if( dim == 3 ) { cmesh->InsertMaterialObject(BCond5); }
    
    cmesh->SetDefaultOrder(pOrder);
    
    
//    if (!isgeoblend) {
//        TPZLagrangeMultiplier *matskelet = new TPZLagrangeMultiplier(matskeleton, dim-1, 1);
//        TPZMaterial * mat2(matskelet);
//        cmesh->InsertMaterialObject(mat2);
//    }
    
    
    
    //Ajuste da estrutura de dados computacional
    cmesh->AutoBuild();
    
    //#ifdef LOG4CXX
    //	if(logdata->isDebugEnabled())
    //	{
    //        std::stringstream sout;
    //        sout<<"\n\n Malha Computacional_1 Fluxo\n ";
    //        cmesh->Print(sout);
    //        LOGPZ_DEBUG(logdata,sout.str())
    //	}
    //#endif
    
    return cmesh;
    
}

TPZCompMesh *LaplaceInSphere::CMeshPressure(TPZGeoMesh *gmesh, int pOrder, int dim)
{
    /// criar materiais
    TPZMatPoisson3d *material = new TPZMatPoisson3d(fmatId,dim);
    material->NStateVariables();
    
    //    TPZAutoPointer<TPZFunction<STATE> > force1 = new TPZDummyFunction<STATE>(Forcing1);
    //	material->SetForcingFunction(force1);
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    
    ///Inserir condicao de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    //    TPZMaterial * BCond0 = material->CreateBC(mat, bc0,dirichlet, val1, val2);
    //    TPZMaterial * BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
    //    TPZMaterial * BCond2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);
    //    TPZMaterial * BCond3 = material->CreateBC(mat, bc3,dirichlet, val1, val2);
    //    TPZMaterial * BCond4 = material->CreateBC(mat, bc4,dirichlet, val1, val2);
    //    TPZMaterial * BCond5 = material->CreateBC(mat, bc5,dirichlet, val1, val2);
    //
    //    cmesh->InsertMaterialObject(BCond0);
    //    cmesh->InsertMaterialObject(BCond1);
    //    cmesh->InsertMaterialObject(BCond2);
    //    cmesh->InsertMaterialObject(BCond3);
    //    cmesh->InsertMaterialObject(BCond4);
    //    cmesh->InsertMaterialObject(BCond5);
    
    cmesh->SetDefaultOrder(pOrder);
    //cmesh->SetDimModel(dim);
    
    bool h1function = true;//false em triangulo
    if(h1function){
        cmesh->SetAllCreateFunctionsContinuous();
        cmesh->ApproxSpace().CreateDisconnectedElements(true);
    }
    else{
        cmesh->SetAllCreateFunctionsDiscontinuous();
    }
    
    
    //Ajuste da estrutura de dados computacional
    cmesh->AutoBuild();
    
    
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        //newnod.SetPressure(true);
        newnod.SetLagrangeMultiplier(1);
    }
    
    if(!h1function)
    {
        
        int nel = cmesh->NElements();
        for(int i=0; i<nel; i++){
            TPZCompEl *cel = cmesh->ElementVec()[i];
            TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
            celdisc->SetConstC(1.);
            celdisc->SetCenterPoint(0, 0.);
            celdisc->SetCenterPoint(1, 0.);
            celdisc->SetCenterPoint(2, 0.);
            celdisc->SetTrueUseQsiEta();
            //celdisc->SetFalseUseQsiEta();
            
            //            TPZVec<REAL> qsi(3,0.);
            //            qsi[0] = 0.5;
            //            qsi[1] = 0.5;
            //            TPZFMatrix<REAL> phi;
            //            TPZFMatrix<REAL> dphi;
            //            celdisc->Shape(qsi, phi,dphi);
            //            phi.Print("phi = ");
            
            
            if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
            {
                if(ftriang==true) celdisc->SetTotalOrderShape();
                else celdisc->SetTensorialShape();
            }
            
        }
    }
    
#ifdef DEBUG
    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }
#endif
    
    
    
    //#ifdef LOG4CXX
    //	if(logdata->isDebugEnabled())
    //	{
    //        std::stringstream sout;
    //        sout<<"\n\n Malha Computacional_2 pressure\n ";
    //        cmesh->Print(sout);
    //        LOGPZ_DEBUG(logdata,sout.str());
    //	}
    //#endif
    return cmesh;
    
}

TPZCompMesh *LaplaceInSphere::CMeshMixed(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec)
{
    
    //Creating computational mesh for multiphysic elements
    gmesh->ResetReference();
    TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    //criando material
    int dim = gmesh->Dimension();
    
    TPZMatMixedPoisson3D *material = new TPZMatMixedPoisson3D(fmatId,dim);
    
    //incluindo os dados do problema
    TPZFNMatrix<9,REAL> PermTensor(3,3,0.);
    TPZFNMatrix<9,REAL> InvPermTensor(3,3,0.);
    
    
    // tensor de permutacao
    TPZFNMatrix<9,REAL> TP(3,3,0.0);
    TPZFNMatrix<9,REAL> InvTP(3,3,0.0);
    
    // Hard coded
    for (int id = 0; id < 3; id++){
        TP(id,id) = 1.0;
        InvTP(id,id) = 1.0;
    }
    
    PermTensor = TP;
    InvPermTensor = InvTP;
    
    material->SetPermeabilityTensor(PermTensor, InvPermTensor);
    
    //solucao exata
    TPZAutoPointer<TPZFunction<STATE> > solexata;
    
    solexata = new TPZDummyFunction<STATE>(SolExata);
    material->SetForcingFunctionExact(solexata);
    mphysics->SetDimModel(dim);
    //funcao do lado direito da equacao do problema
    TPZDummyFunction<STATE> *dum = new TPZDummyFunction<STATE>(Forcing);
    TPZAutoPointer<TPZFunction<STATE> > forcef;
    dum->SetPolynomialOrder(1);
    forcef = dum;
    material->SetForcingFunction(forcef);
    
    //inserindo o material na malha computacional
    TPZMaterial *mat(material);
    mphysics->InsertMaterialObject(mat);
    
    
    //Criando condicoes de contorno
    TPZMaterial * BCond0;
    TPZMaterial * BCond1;
    TPZMaterial * BCond2;
    TPZMaterial * BCond3;
    TPZMaterial * BCond4;
    TPZMaterial * BCond5;
    
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    if (dim==3)
    {
        val2(0,0) = 0.0;
        val2(1,0) = 0.0;
        TPZAutoPointer<TPZFunction<STATE> > FBCond0 = new TPZDummyFunction<STATE>(ForcingBC0D);
        BCond0 = material->CreateBC(mat, fbc0,fdirichlet, val1, val2);
        BCond0->SetForcingFunction(FBCond0);
    }
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZAutoPointer<TPZFunction<STATE> > FBCond1 = new TPZDummyFunction<STATE>(ForcingBC1N);
    BCond1 = material->CreateBC(mat, fbc1,fneumann, val1, val2);
    BCond1->SetForcingFunction(FBCond1);
    //    BCond1 = material->CreateBC(mat, bc1,neumann, val1, val2);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZAutoPointer<TPZFunction<STATE> > FBCond2 = new TPZDummyFunction<STATE>(ForcingBC2D);
    BCond2 = material->CreateBC(mat, fbc2,fdirichlet, val1, val2);
    BCond2->SetForcingFunction(FBCond2);
    //    TPZAutoPointer<TPZFunction<STATE> > FBCond2 = new TPZDummyFunction<STATE>(ForcingBC2N);
    //    BCond2 = material->CreateBC(mat, bc2,neumann, val1, val2);
    //    BCond2->SetForcingFunction(FBCond2);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZAutoPointer<TPZFunction<STATE> > FBCond3 = new TPZDummyFunction<STATE>(ForcingBC3N);
    BCond3 = material->CreateBC(mat, fbc3,fneumann, val1, val2);
    BCond3->SetForcingFunction(FBCond3);
    //    BCond3 = material->CreateBC(mat, bc3,neumann, val1, val2);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZAutoPointer<TPZFunction<STATE> > FBCond4 = new TPZDummyFunction<STATE>(ForcingBC4D);
    BCond4 = material->CreateBC(mat, fbc4,fdirichlet, val1, val2);
    BCond4->SetForcingFunction(FBCond4);
    //    TPZAutoPointer<TPZFunction<STATE> > FBCond4 = new TPZDummyFunction<STATE>(ForcingBC4N);
    //    BCond4 = material->CreateBC(mat, bc4,neumann, val1, val2);
    //    BCond4->SetForcingFunction(FBCond4);
    
    if (dim==3)
    {
        val2(0,0) = 0.0;
        val2(1,0) = 0.0;
        TPZAutoPointer<TPZFunction<STATE> > FBCond5 = new TPZDummyFunction<STATE>(ForcingBC5D);
        BCond5 = material->CreateBC(mat, fbc5,fdirichlet, val1, val2);
        BCond5->SetForcingFunction(FBCond5);
    }
    
    
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    if( dim == 3 ) { mphysics->InsertMaterialObject(BCond0); }
    mphysics->InsertMaterialObject(BCond1);
    mphysics->InsertMaterialObject(BCond2);
    mphysics->InsertMaterialObject(BCond3);
    mphysics->InsertMaterialObject(BCond4);
    if( dim == 3 ) { mphysics->InsertMaterialObject(BCond5); }
    
    //Fazendo auto build
    mphysics->AutoBuild();
    mphysics->AdjustBoundaryElements();
    mphysics->CleanUpUnconnectedNodes();
    
    
    TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
    
    //        TPZMaterial * skeletonEl = material->CreateBC(mat, matskeleton, 3, val1, val2);
    //        mphysics->InsertMaterialObject(skeletonEl);
    
    //        TPZLagrangeMultiplier *matskelet = new TPZLagrangeMultiplier(matskeleton, dim-1, 1);
    //        TPZMaterial * mat2(matskelet);
    //        mphysics->InsertMaterialObject(mat2);
    
    int nel = mphysics->ElementVec().NElements();
    TPZStack< TPZStack< TPZMultiphysicsElement *,7> > wrapEl;
    for(int el = 0; el < nel; el++)
    {
        TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *>(mphysics->Element(el));
        if(mfcel->Dimension()==dim) TPZBuildMultiphysicsMesh::AddWrap(mfcel, fmatId, wrapEl);//criei elementos com o mesmo matId interno, portanto nao preciso criar elemento de contorno ou outro material do tipo TPZLagrangeMultiplier
    }
    meshvec[0]->CleanUpUnconnectedNodes();
    TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
    
    //------- Create and add group elements -------
    long index, nenvel;
    nenvel = wrapEl.NElements();
    for(int ienv=0; ienv<nenvel; ienv++){
        TPZElementGroup *elgr = new TPZElementGroup(*wrapEl[ienv][0]->Mesh(),index);
        nel = wrapEl[ienv].NElements();
        for(int jel=0; jel<nel; jel++){
            elgr->AddElement(wrapEl[ienv][jel]);
        }
    }
    
    
    return mphysics;
    
}


void LaplaceInSphere::ErrorHDiv(TPZCompMesh *hdivmesh, int p, int ndiv, std::map<REAL, REAL> &fDebugMapL2, std::map<REAL, REAL> &fDebugMapHdiv)
{
    long nel = hdivmesh->NElements();
    int dim = hdivmesh->Dimension();
    TPZManVector<STATE,10> globalerrors(10,0.);
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = hdivmesh->ElementVec()[el];
        if(cel->Reference()->Dimension()!=dim) continue;
        TPZManVector<STATE,10> elerror(10,0.);
        elerror.Fill(0.);
        cel->EvaluateError(SolExata, elerror, NULL);
        int nerr = elerror.size();
        for (int i=0; i<nerr; i++) {
            globalerrors[i] += elerror[i]*elerror[i];
        }
        
    }
    //    out << "Errors associated with HDiv space - ordem polinomial = " << p << "- divisoes = " << ndiv << endl;
    //    out << "L2 Norm for flux - "<< endl; //L2 Norm for divergence - Hdiv Norm for flux " << endl;
    //    out <<  setw(16) << sqrt(globalerrors[1]) <<endl;// setw(25)  << sqrt(globalerrors[2]) << setw(21)  << sqrt(globalerrors[3]) << endl;
    //
    //    out << "L2 Norm for flux = "    << sqrt(globalerrors[1]) << endl;
    //    out << "L2 Norm for divergence = "    << sqrt(globalerrors[2])  <<endl;
    //    out << "Hdiv Norm for flux = "    << sqrt(globalerrors[3])  <<endl;
    //
    fDebugMapHdiv.insert(std::pair<REAL, REAL> (ndiv,sqrt(globalerrors[1])));
}

void LaplaceInSphere::ErrorL2(TPZCompMesh *l2mesh, int p, int ndiv, std::map<REAL, REAL> &fDebugMapL2, std::map<REAL, REAL> &fDebugMapHdiv)
{
    long nel = l2mesh->NElements();
    //int dim = l2mesh->Dimension();
    TPZManVector<STATE,10> globalerrors(10,0.);
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = l2mesh->ElementVec()[el];
        TPZManVector<STATE,10> elerror(10,0.);
        cel->EvaluateError(SolExata, elerror, NULL);
        int nerr = elerror.size();
        globalerrors.resize(nerr);
        //#ifdef LOG4CXX
        //        if (logdata->isDebugEnabled()) {
        //            std::stringstream sout;
        //            sout << "L2 Error sq of element " << el << elerror[0]*elerror[0];
        //            LOGPZ_DEBUG(logdata, sout.str())
        //        }
        //#endif
        for (int i=0; i<nerr; i++) {
            globalerrors[i] += elerror[i]*elerror[i];
        }
        
    }
    //    out << "Errors associated with L2 space - ordem polinomial = " << p << "- divisoes = " << ndiv << endl;
    //    out << "L2 Norm = "    << sqrt(globalerrors[1]) << endl;
    fDebugMapL2.insert(std::pair<REAL, REAL> (ndiv,sqrt(globalerrors[1])));
}







