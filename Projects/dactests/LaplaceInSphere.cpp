//
//  LaplaceInSphere.cpp
//  PZ
//
//  Created by Douglas Castro on 1/29/15.
//
//

#include "LaplaceInSphere.h"
#include "pzcheckgeom.h"
#include "tools.h"

//#define TROPICO
#define CompleteShpere
const int  norder = 6;

LaplaceInSphere::LaplaceInSphere()
{
    fDim = 2;
    fmatId = 1;
    fdirichlet = 0;
    fneumann = 1;
    fbc0 = -1;
    fbc1 = -2;
    fbc2 = -3;
    fbc3 = -4;
    fbc4 = -5;
    fbc5 = -6;
    fmatskeleton = -7;
    fisH1 = false;
    ftriang = false;
    isgeoblend = true;
}

LaplaceInSphere::~LaplaceInSphere()
{
    
}

void LaplaceInSphere::Run(int ordemP, int ndiv, std::map<REAL, REAL> &fDebugMapL2, std::map<REAL, REAL> &fDebugMapHdiv, std::ofstream &saidaErro, bool HdivMaisMais)
{
    std::cout<< " INICIO(CASCA ESFERA) - grau  polinomio " << ordemP << " numero de divisoes " << ndiv << std::endl;
    std::cout<< " Dimensao == " << fDim << std::endl;
    
    
#ifdef TROPICO
//    TPZVec<bool> CurvesSides(4,true);
//    TPZGeoMesh *gmesh = GMeshTropicodeCancer(ndiv, CurvesSides, false, 1);
    TPZVec<bool> CurvesSides(3,true);
    TPZGeoMesh *gmesh = GMeshCirculoPolarArtico(ndiv, CurvesSides, false, 1);
#else
    //TPZGeoMesh *gmesh = GMeshSphericalShell(ndiv);
    //TPZGeoMesh *gmesh = GMeshOitavoCasca( ndiv);
    TPZGeoMesh *gmesh = MakeSphereFromQuadrilateral(2, false, ndiv);
    //TPZGeoMesh *gmesh = GMeshSphericalRingQuarter(2, false, ndiv);
    //TPZGeoMesh *gmesh = GMeshSphericalShell(2, true, ndiv);
    //TPZGeoMesh *gmesh = GMeshSphericalShell(2,false, ndiv);
    //TPZGeoMesh *gmesh = GMeshSphericalShellBlendQ(2, false, ndiv);
#endif
    
    gmesh->SetDimension(fDim);
    {
        //        ofstream argm("gmesh2d-Esfera.txt");
        //        gmesh->Print(argm);
    }
    
    
    // Um teste para a solucao via H1, sem hdiv
    if (fisH1)
    {
        TPZCompMesh *cmeshH1 = CMeshH1(gmesh, ordemP, fDim);
        TPZAnalysis anh1(cmeshH1, true);
        REAL t1,t2;
        tools::SolveSyst(anh1, cmeshH1, t1, t2);
        
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
    
    TPZCompMesh *cmesh2 = CMeshPressure(gmesh, ordemP, fDim);
    TPZCompMesh *cmesh1 = CMeshFlux(gmesh, ordemP, fDim);
    if (HdivMaisMais) {
        // para rodar P**
        ChangeExternalOrderConnects(cmesh1);
    }
    int DofCond, DoFT;
    DoFT = cmesh1->NEquations() + cmesh2->NEquations();
    
//    TPZCheckGeom *checaerro = new TPZCheckGeom(gmesh);
//    int nerro = checaerro->PerformCheck();
//    if (nerro>0) {
//        std::cout << " Erro malha geometrica " << std::endl;
//        DebugStop();
//    }
    
    
    {
        //        ofstream arg1("cmeshflux.txt");
        //        cmesh1->Print(arg1);
        //
        //        ofstream arg2("cmeshpressure.txt");
        //        cmesh2->Print(arg2);
        //
        //        ofstream arg4("gmesh2.txt");
        //        gmesh->Print(arg4);
        
    }
    
    //malha multifisica
    TPZVec<TPZCompMesh *> meshvec(2);
    meshvec[0] = cmesh1;
    meshvec[1] = cmesh2;
    
    TPZCompMesh * mphysics = CMeshMixed(gmesh,meshvec);
    
    std::ofstream sout("FluxcmeshFinal.txt");
    meshvec[0]->Print(sout);
    
    
    std::ofstream sout2("cmeshmphysics.txt");
    mphysics->Print(sout2);
    
    DofCond = mphysics->NEquations();
    
    //TestMesh(mphysics);
    {
        //        ofstream arg5("cmeshmultiphysics.txt");
        //        mphysics->Print(arg5);
    }
    
    TPZAnalysis an(mphysics, true);
    REAL t1,t2;
    tools::SolveSyst(an, mphysics, t1, t2);
    
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
    
    std::cout << "Postprocessed - inicio calculo do erro\n";
    
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
    
    ErrorPrimalDual( cmesh2, cmesh1,  ordemP, ndiv, saidaErro, DoFT, DofCond);

    
    tools::PrintDebugMapForMathematica(HdivData, L2Data, fDebugMapL2, fDebugMapHdiv);
    
    std::cout<< " FIM (ESFERA) - grau  polinomio " << ordemP << " numero de divisoes " << ndiv << std::endl;
    

}

TPZGeoMesh *LaplaceInSphere::GMeshSphericalShell(int ndiv)
{

    bool ftriangulo = true;//false;//
    
    ///INSTANCIACAO DA MALHA GEOMETRICA
    TPZGeoMesh * geomesh = new TPZGeoMesh;
    
    /// Materiais
    int64_t materialId = fmatId;
    
    int nnodes = 9;//quantidade de nos da malha geometrica
    geomesh->NodeVec().Resize(nnodes);
    geomesh->SetDimension(2);
    
    
    
    ///INICIALIZACAO DA MALHA GEOMETRICA PELA INSTANCIACAO E INICIALIZACAO DOS NOS
    TPZGeoNode node;
    TPZVec<REAL> coord(3,0.);
    const REAL r = 1.;
    
    TPZManVector<REAL,3> xc(3,0.);
    xc[0] = 0.;
    xc[1] = 0.;
    xc[2] = 0.;
    
    int Axis = 3;
    REAL angrot = 0.0;
    

    
    int id = 0;
    //no 0
    coord[0] = xc[0] + r*sin( M_PI/2.0)*cos(0.0);
    coord[1] = xc[1] + r*sin( M_PI/2.0)*sin(0.0);
    coord[2] = xc[2] + r*cos( M_PI/2.0);
    tools::RotateNode(coord, angrot, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 1
    coord[0] = xc[0] + r*sin( M_PI/2.0)*cos(M_PI/2.0);
    coord[1] = xc[1] + r*sin( M_PI/2.0)*sin(M_PI/2.0);
    coord[2] = xc[2] + r*cos( M_PI/2.0);
    tools::RotateNode(coord, angrot, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 2
    coord[0] = xc[0] + r*sin( M_PI/2.0)*cos(M_PI);
    coord[1] = xc[1] + r*sin( M_PI/2.0)*sin(M_PI);
    coord[2] = xc[2] + r*cos( M_PI/2.0);
    tools::RotateNode(coord, angrot, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 3
    coord[0] = xc[0] + r*sin( M_PI/2.0)*cos(3.0*M_PI/2.0);
    coord[1] = xc[1] + r*sin( M_PI/2.0)*sin(3.0*M_PI/2.0);
    coord[2] = xc[2] + r*cos( M_PI/2.0);
    tools::RotateNode(coord, angrot, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 4
    coord[0] = xc[0] + r*sin( M_PI/2.0)*cos(M_PI/4.0);
    coord[1] = xc[1] + r*sin( M_PI/2.0)*sin(M_PI/4.0);
    coord[2] = xc[2] + r*cos( M_PI/2.0);
    tools::RotateNode(coord, angrot, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 5
    coord[0] = xc[0] + r*sin( M_PI/2.0)*cos(3.0*M_PI/4.0);
    coord[1] = xc[1] + r*sin( M_PI/2.0)*sin(3.0*M_PI/4.0);
    coord[2] = xc[2] + r*cos( M_PI/2.0);
    tools::RotateNode(coord, angrot, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 6
    coord[0] = xc[0] + r*sin( M_PI/2.0)*cos(5.0*M_PI/4.0);
    coord[1] = xc[1] + r*sin( M_PI/2.0)*sin(5.0*M_PI/4.0);
    coord[2] = xc[2] + r*cos( M_PI/2.0);
    tools::RotateNode(coord, angrot, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 7
    coord[0] = xc[0] + r*sin( M_PI/2.0)*cos(7.0*M_PI/4.0);
    coord[1] = xc[1] + r*sin( M_PI/2.0)*sin(7.0*M_PI/4.0);
    coord[2] = xc[2] + r*cos( M_PI/2.0);
    tools::RotateNode(coord, angrot, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 8
    coord[0] = xc[0] + r*sin( 0.0)*cos(0.0);
    coord[1] = xc[1] + r*sin( 0.0)*sin(0.0);
    coord[2] = xc[2] + r*cos( 0.0);
    tools::RotateNode(coord, angrot, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    //id++;
    
    int elementid = 0;
    // Using triangle to sphere special map
    TPZVec<int64_t> topology(4);
    
    
    if (ftriangulo)
    {
        
        topology.resize(3);
        
        // El 0
        topology[0] = 0;
        topology[1] = 1;
        topology[2] = 8;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > * SphereRingT1 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > (elementid, topology,materialId,*geomesh);
            SphereRingT1->Geom().SetData(r,xc);
            elementid++;
        }
        // El 1
        topology[0] = 1;
        topology[1] = 2;
        topology[2] = 8;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > * SphereRingT2 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > (elementid, topology,materialId,*geomesh);
            SphereRingT2->Geom().SetData(r,xc);
            elementid++;
        }
        // El 2
        topology[0] = 2;
        topology[1] = 3;
        topology[2] = 8;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > * SphereRingT1 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > (elementid, topology,materialId,*geomesh);
            SphereRingT1->Geom().SetData(r,xc);
            elementid++;
        }
        // El 3
        topology[0] = 3;
        topology[1] = 0;
        topology[2] = 8;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > * SphereRingT2 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > (elementid, topology,materialId,*geomesh);
            SphereRingT2->Geom().SetData(r,xc);
            elementid++;
        }
        
    }
    else
    {
        
        topology.resize(4);
        
        // El 0
        topology[0] = 0;
        topology[1] = 4;
        topology[2] = 1;
        topology[3] = 8;
        {
            TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > * SphereRingQ =
            new TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > (elementid, topology,materialId,*geomesh);
            SphereRingQ->Geom().SetData(r,xc);
            elementid++;
        }
        // El 1
        topology[0] = 1;
        topology[1] = 5;
        topology[2] = 2;
        topology[3] = 8;
        {
            TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > * SphereRingQ =
            new TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > (elementid, topology,materialId,*geomesh);
            SphereRingQ->Geom().SetData(r,xc);
            elementid++;
        }
        // El 2
        topology[0] = 2;
        topology[1] = 6;
        topology[2] = 3;
        topology[3] = 8;
        {
            TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > * SphereRingQ =
            new TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > (elementid, topology,materialId,*geomesh);
            SphereRingQ->Geom().SetData(r,xc);
            elementid++;
        }
        // El 3
        topology[0] = 3;
        topology[1] = 7;
        topology[2] = 0;
        topology[3] = 8;
        {
            TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > * SphereRingQ =
            new TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > (elementid, topology,materialId,*geomesh);
            SphereRingQ->Geom().SetData(r,xc);
            elementid++;
        }
        
    }
    
    
    // El linha
    // Definition of Arc coordenates
    topology.resize(2);
    // Create Geometrical Arc #0
    topology[0] = 0;
    topology[1] = 1;
    {
        //TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc =
        new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, fbc1 /** fbc3 */, *geomesh);
        //arc->Geom().Print(std::cout);
    }
    elementid++;
    
    // Create Geometrical Arc #1
    topology[0] = 1;
    topology[1] = 2;
    {
        //TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc =
        new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, fbc1 /** fbc4 */, *geomesh);
        // arc->Geom().Print(std::cout);
    }
    elementid++;
    
    // Create Geometrical Arc #2
    topology[0] = 2;
    topology[1] = 3;
    {
        //TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc =
        new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, fbc1 /** fbc1 */, *geomesh);
        //arc->Geom().Print(std::cout);
    }
    elementid++;
    
    // Create Geometrical Arc #3
    topology[0] = 3;
    topology[1] = 0;
    {
        //TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc =
        new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, fbc1 /** fbc2 */, *geomesh);
        //arc->Geom().Print(std::cout);
    }
    
    
    //CONCLUINDO A CONSTRUCAO DA MALHA GEOMETRICA
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
    
    std::ofstream outfile("malhaCascaesfera.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(geomesh, outfile, true);
    
    return geomesh;
    
}


//-----------------------------------------------------------------


TPZGeoMesh *LaplaceInSphere::GMeshSphericalRingQuarter(int dimensao, bool triang, int ndiv)
{
    
    bool ftriangulo = triang;
    
    ///INSTANCIACAO DA MALHA GEOMETRICA
    TPZGeoMesh * geomesh = new TPZGeoMesh;
    
    /// Materiais
    int64_t materialId = fmatId;
    
    int nnodes = 4;//quantidade de nos da malha geometrica
    geomesh->NodeVec().Resize(nnodes);
    geomesh->SetDimension(dimensao);
    
    
    
    ///INICIALIZACAO DA MALHA GEOMETRICA PELA INSTANCIACAO E INICIALIZACAO DOS NOS
    TPZGeoNode node;
    TPZVec<REAL> coord(3,0.);
    const REAL r = 1.;

    TPZManVector<REAL,3> xc(3,0.);
    xc[0] = 0.;
    xc[1] = 0.;
    xc[2] = 0.;
    
    REAL theta = M_PI/2.0 ;
    REAL phi = M_PI/6.0;
    
    int Axis = 3;
    REAL angulo = 0.0;
    
    int id = 0;
    //no 0
    coord[0] = xc[0] + r*sin(theta + M_PI/6.0)*cos(-phi);
    coord[1] = xc[1] + r*sin(theta + M_PI/6.0)*sin(-phi);
    coord[2] = xc[2] + r*cos(theta + M_PI/6.0);
    
    tools::RotateNode(coord, angulo, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 1
    coord[0] = xc[0] + r*sin(theta + M_PI/6.0)*cos(phi);
    coord[1] = xc[1] + r*sin(theta + M_PI/6.0)*sin(phi);
    coord[2] = xc[2] + r*cos(theta + M_PI/6.0);
    tools::RotateNode(coord, angulo, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 2
    coord[0] = xc[0] + r*sin(theta - M_PI/6.0)*cos(phi);
    coord[1] = xc[1] + r*sin(theta - M_PI/6.0)*sin(phi);
    coord[2] = xc[2] + r*cos(theta - M_PI/6.0);
    tools::RotateNode(coord, angulo, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    
    //no 3
    coord[0] = xc[0] + r*sin(theta - M_PI/6.0)*cos(-phi);
    coord[1] = xc[1] + r*sin(theta - M_PI/6.0)*sin(-phi);
    coord[2] = xc[2] + r*cos(theta - M_PI/6.0);
    tools::RotateNode(coord, angulo, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    //id++;
    
    int elementid = 0;
    
    // Using triangle to sphere special map
    TPZVec<int64_t> topology(4);
    
    if (ftriangulo)
    {
        
        
        topology.Resize(3);
        
        // El 0
        topology[0] = 0;
        topology[1] = 1;
        topology[2] = 2;
        
        TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > * SphereEighth1 =
        new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > (elementid, topology,materialId,*geomesh);
        SphereEighth1->Geom().SetData(r,xc);
        elementid++;
        
        // El 0
        topology[0] = 0;
        topology[1] = 2;
        topology[2] = 3;
        
        TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > * SphereEighth2 =
        new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > (elementid, topology,materialId,*geomesh);
        SphereEighth2->Geom().SetData(r,xc);
        //elementid++;
        
        
    }
    else
    {
        
        topology.Resize(4);
        
        // El 0
        topology[0] = 0;
        topology[1] = 1;
        topology[2] = 2;
        topology[3] = 3;
        
        TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > * SphereEighth1 =
        new TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > (elementid, topology,materialId,*geomesh);
        SphereEighth1->Geom().SetData(r,xc);
        elementid++;
        
    }
    
    // El linha
    // Definition of Arc coordenates
    topology.resize(2);
    // Create Geometrical Arc #1
    topology[0] = 0;
    topology[1] = 1;
    {
        //TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc =
        new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, fbc1, *geomesh);
        //arc->Geom().Print(std::cout);
    }
    //new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,topology, fbc3, *geomesh);
    elementid++;
    
    // Create Geometrical Arc #2
    topology[0] = 1;
    topology[1] = 2;
    {
        //TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc =
        new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, fbc2, *geomesh);
        // arc->Geom().Print(std::cout);
    }
    //new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,topology, fbc4, *geomesh);
    elementid++;
    
    // Create Geometrical Arc #3
    topology[0] = 2;
    topology[1] = 3;
    {
        //TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc =
        new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, fbc3, *geomesh);
        //arc->Geom().Print(std::cout);
    }
    //new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,topology, fbc1, *geomesh);
    elementid++;
    
    // Create Geometrical Arc #4
    topology[0] = 3;
    topology[1] = 0;
    {
        //TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc =
        new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, fbc4, *geomesh);
        //arc->Geom().Print(std::cout);
    }
    //new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,topology, fbc2, *geomesh);
    //    elementid++;
    
    
    //CONCLUINDO A CONSTRUCAO DA MALHA GEOMETRICA
    geomesh->BuildConnectivity();
    
    
    
//    const unsigned int nel = 1;//geomesh->NElements();
//    for (int i = 0 ; i < nel ; i++){
//        TPZGeoEl *gel = geomesh->Element(i);
//        if (!gel) {
//            continue;
//        }
//        const int npt = 5;
//        TPZManVector<REAL,3> qsi(3,0.), x(3,0.);
//        for (int iqsi = 0; iqsi < npt; iqsi++) {
//            for (int ieta = 0; ieta < npt; ieta++) {
//                qsi[0] = -1 + iqsi * 2./(npt-1);
//                qsi[1] = -1 + ieta * 2./(npt-1);
//                gel->X(qsi, x);
//                //calcula a funcao e ve se zero
//                cout << " qsi = " << qsi << " x = " << x << endl;
//                cout << " r = " << sqrt(x[0]*x[0]+x[1]*x[1]) << endl;
//                cout << " R = " << sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]) << endl;
//                //cout << " phi
//            }
//        }
//    }
    
    
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
    
    std::ofstream outfile("malhaCascaesferaQuarto.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(geomesh, outfile, true);
    
    return geomesh;
    
    
}



TPZGeoMesh *LaplaceInSphere::GMeshCirculoPolarArtico(int ndiv , TPZVec<bool>  &CurvesSides, bool isPlane, int plane)
{
    
    
    ///INSTANCIACAO DA MALHA GEOMETRICA
    TPZGeoMesh * geomesh = new TPZGeoMesh;
    
    /// Materiais
    int64_t materialId = fmatId;

    
    int nnodes = 6;//quantidade de nos da malha geometrica
    geomesh->NodeVec().Resize(nnodes);
    geomesh->SetDimension(2);
    
    if (CurvesSides.size() != 3) {
        DebugStop();
    }
    
    
    ///INICIALIZACAO DA MALHA GEOMETRICA PELA INSTANCIACAO E INICIALIZACAO DOS NOS
    TPZGeoNode node;
    //    REAL beta = M_PI/16.0;
    //    REAL alpha = M_PI/4.0;
    //    REAL h = M_PI/4.0;
    //    REAL k = M_PI/4.0;
    
    REAL beta = M_PI/3.0;
    REAL alpha = M_PI/4.0;
    REAL h = M_PI/2.0;
    REAL k = 0.0;
    
    int id=0;
    TPZManVector<REAL,3> coord(3,0.);
    
    if (isPlane) {
        
        switch (plane) {
            case 1:
            {
                //no 0
                coord = SphereToKartesian(1.0, h-beta, k-alpha);
                node.SetNodeId(id);
                node.SetCoord(coord);
                geomesh->NodeVec()[id] = node;
                id++;
                cout << coord << endl;
                //no 1
                coord = SphereToKartesian(1.0, h-beta, k+alpha);
                node.SetNodeId(id);
                node.SetCoord(coord);
                geomesh->NodeVec()[id] = node;
                id++;
                cout << coord << endl;
                //no 2
                coord = SphereToKartesian(1.0, h-beta, k+alpha + 2.0*(h-alpha));
                node.SetNodeId(id);
                node.SetCoord(coord);
                geomesh->NodeVec()[id] = node;
                id++;
                cout << coord << endl;
                //no 3
                coord = SphereToKartesian(1.0, h-beta, k-alpha - 2.0*(h-alpha));
                node.SetNodeId(id);
                node.SetCoord(coord);
                geomesh->NodeVec()[id] = node;
                id++;
                cout << coord << endl;
            }
                break;
                
            case 2:
            {
                //no 0
                coord = SphereToKartesian(1.0, h-beta, k-alpha);
                node.SetNodeId(id);
                node.SetCoord(coord);
                geomesh->NodeVec()[id] = node;
                id++;
                cout << coord << endl;
                //no 1
                coord = SphereToKartesian(1.0, h-beta, k+alpha);
                node.SetNodeId(id);
                node.SetCoord(coord);
                geomesh->NodeVec()[id] = node;
                id++;
                cout << coord << endl;
                //no 2
                coord = SphereToKartesian(1.0, h+beta, k+alpha);
                node.SetNodeId(id);
                node.SetCoord(coord);
                geomesh->NodeVec()[id] = node;
                id++;
                cout << coord << endl;
                //no 3
                coord = SphereToKartesian(1.0, h+beta, k-alpha);
                node.SetNodeId(id);
                node.SetCoord(coord);
                geomesh->NodeVec()[id] = node;
                id++;
                cout << coord << endl;
            }
                break;
                
            case 3:
            {
                //no 0
                coord = SphereToKartesian(1.0, h-beta, k+alpha + 2.0*(h-alpha));
                node.SetNodeId(id);
                node.SetCoord(coord);
                geomesh->NodeVec()[id] = node;
                id++;
                cout << coord << endl;
                //no 1
                coord = SphereToKartesian(1.0, h-beta, k + alpha);
                node.SetNodeId(id);
                node.SetCoord(coord);
                geomesh->NodeVec()[id] = node;
                id++;
                cout << coord << endl;
                //no 2
                coord = SphereToKartesian(1.0, h+beta, k + alpha);
                node.SetNodeId(id);
                node.SetCoord(coord);
                geomesh->NodeVec()[id] = node;
                id++;
                cout << coord << endl;
                //no 3
                coord = SphereToKartesian(1.0, h+beta, k + alpha + 2.0*(h-alpha));
                node.SetNodeId(id);
                node.SetCoord(coord);
                geomesh->NodeVec()[id] = node;
                id++;
                cout << coord << endl;
            }
                break;
                
            case 4:
            {
                // domain with equal lentgh
                //no 0
                coord[0] = -M_PI/8.0;
                coord[1] = -sqrt(3.0)/2.0;
                coord[2] = 0.0;
                node.SetNodeId(id);
                node.SetCoord(coord);
                geomesh->NodeVec()[id] = node;
                id++;
                
                //no 1
                coord[0] = M_PI/8.0;
                coord[1] = -sqrt(3.0)/2.0;
                coord[2] = 0.0;
                
                node.SetNodeId(id);
                node.SetCoord(coord);
                geomesh->NodeVec()[id] = node;
                id++;
                
                //no 2
                coord[0] = M_PI/8.0;
                coord[1] = sqrt(3.0)/2.0;
                coord[2] = 0.0;
                node.SetNodeId(id);
                node.SetCoord(coord);
                geomesh->NodeVec()[id] = node;
                id++;
                
                //no 3
                coord[0] = -M_PI/8.0;
                coord[1] = sqrt(3.0)/2.0;
                coord[2] = 0.0;
                node.SetNodeId(id);
                node.SetCoord(coord);
                geomesh->NodeVec()[id] = node;
                id++;
                
            }
                break;
                
            default:
            {
                DebugStop();
            }
                break;
        }
        
    }
    else{
        //no 0
        coord = SphereToKartesian(1.0, h - beta , k - alpha );
        node.SetNodeId(id);
        node.SetCoord(coord);
        geomesh->NodeVec()[id] = node;
        id++;
        
        //no 1
        coord = SphereToKartesian(1.0, h - beta , k + alpha );
        node.SetNodeId(id);
        node.SetCoord(coord);
        geomesh->NodeVec()[id] = node;
        id++;
        
        //no 2 // POLO N
        coord = SphereToKartesian(1.0, 0.0 , 0.0 );
        node.SetNodeId(id);
        node.SetCoord(coord);
        geomesh->NodeVec()[id] = node;
        id++;
        
        //no 3
        coord = SphereToKartesian(1.0, h - beta , k);
        node.SetNodeId(id);
        node.SetCoord(coord);
        geomesh->NodeVec()[id] = node;
        id++;
        
        //no 4
        coord = SphereToKartesian(1.0, (h - beta)/2.0 , k + alpha);
        node.SetNodeId(id);
        node.SetCoord(coord);
        geomesh->NodeVec()[id] = node;
        id++;
        
        //no 5
        coord = SphereToKartesian(1.0, (h - beta)/2.0, k - alpha);
        node.SetNodeId(id);
        node.SetCoord(coord);
        geomesh->NodeVec()[id] = node;
        id++;
        
    }
    
    int elementid = 0;
    // Using triangle to sphere special map
    TPZVec<int64_t> topology(3);
    TPZVec<int64_t> topologyLine(2);
    
    
    // Side 0
    if (CurvesSides[0]) {
        topology[0] = 0;
        topology[1] = 1;
        topology[2] = 3;
        new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,topology, fbc1, *geomesh);
        elementid++;
    }
    else{
        topologyLine[0] = 0;
        topologyLine[1] = 1;
        new TPZGeoElRefPattern < pzgeom::TPZGeoLinear > (elementid,topologyLine, fbc1, *geomesh);
        elementid++;
    }
    
    // Side 1
    if (CurvesSides[1]) {
        topology[0] = 1;
        topology[1] = 2;
        topology[2] = 4;
        new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,topology, fbc2, *geomesh);
        elementid++;
    }
    else{
        topologyLine[0] = 1;
        topologyLine[1] = 2;
        new TPZGeoElRefPattern < pzgeom::TPZGeoLinear > (elementid,topologyLine, fbc2, *geomesh);
        elementid++;
    }
    
    // Side 2
    if (CurvesSides[2]) {
        topology[0] = 2;
        topology[1] = 0;
        topology[2] = 5;
        new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,topology, fbc3, *geomesh);
        elementid++;
    }
    else{
        topologyLine[0] = 2;
        topologyLine[1] = 0;
        new TPZGeoElRefPattern < pzgeom::TPZGeoLinear > (elementid,topologyLine,  fbc3, *geomesh);
        elementid++;
    }
    
    
    
    REAL r = 1.0;
    TPZManVector<REAL,3> xc(3,0.0);
    
    
    // Create Geometrical Quad #1
    topology.Resize(3);
    topology[0] = 0;
    topology[1] = 1;
    topology[2] = 2;
    TPZGeoElRefPattern< pzgeom::TPZTriangleSphere< pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle > > > *TS = new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere< pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle > > > (elementid, topology, materialId, *geomesh);
    TS->Geom().SetData(r,xc);
    
    //CONCLUINDO A CONSTRUCAO DA MALHA GEOMETRICA
    geomesh->BuildConnectivity();
    

    
    //    const unsigned int nel = geomesh->NElements();
    //    for (int i = 0 ; i < nel ; i++){
    //        TPZGeoEl *gel = geomesh->Element(i);
    //        if (!gel) {
    //            continue;
    //        }
    //        const int npt = 3;
    //        TPZManVector<REAL,3> qsi(3,0.), x(3,0.);
    //        for (int iqsi = 0; iqsi < npt; iqsi++) {
    //            for (int ieta = 0; ieta < npt; ieta++) {
    //                qsi[0] = -1 + iqsi * 2./(npt-1);
    //                qsi[1] = -1 + ieta * 2./(npt-1);
    //                gel->X(qsi, x);
    //                //calcula a funcao e ve se zero
    //                cout << " qsi = " << qsi << " x = " << x << endl;
    //                cout << " r = " << sqrt(x[0]*x[0]+x[1]*x[1]) << endl;
    //                cout << " R = " << sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]) << endl;
    //
    //            }
    //        }
    //    }
    
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
    
    std::ofstream outfileAnel("CirculoPolarArtico.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(geomesh, outfileAnel, true);
    
    return geomesh;
    
    
}


TPZGeoMesh *LaplaceInSphere::GMeshOitavoCasca(int ndiv)
{
    
    
    ///INSTANCIACAO DA MALHA GEOMETRICA
    TPZGeoMesh * geomesh = new TPZGeoMesh;
    
    /// Materiais
    int64_t materialId = fmatId;
    
    REAL r = 1.0;
    TPZManVector<REAL,3> xc(3,0.0);
    
    int nnodes = 6;//quantidade de nos da malha geometrica
    geomesh->NodeVec().Resize(nnodes);
    geomesh->SetDimension(2);
    
    
    ///INICIALIZACAO DA MALHA GEOMETRICA PELA INSTANCIACAO E INICIALIZACAO DOS NOS
    TPZGeoNode node;
    
    
    
    int id=0;
    TPZManVector<REAL,3> coord(3,0.);
    
    
    //no 0
    coord = SphereToKartesian(xc, r, M_PI/2.0, 0.0 );
    //coord = SphereToKartesian(r, h - beta , k - alpha );
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 1
    coord = SphereToKartesian(xc, r, M_PI/2.0, M_PI/2.0 );
    //coord = SphereToKartesian(r, h - beta , k + alpha );
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 2 // POLO N
    coord = SphereToKartesian(xc, r, 0.0 , 0.0 );
    //coord = SphereToKartesian( r, 0.0 , 0.0 );
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 3
    coord = SphereToKartesian(xc, r, M_PI/2.0 , M_PI/4.0);
    //coord = SphereToKartesian(r, h - beta , k);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 4
    coord = SphereToKartesian(xc, r, M_PI/4.0, M_PI/2.0);
    //coord = SphereToKartesian( r, (h - beta)/2.0 , k + alpha);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 5
    coord = SphereToKartesian(xc, r, M_PI/4.0 , 0.0);
    //coord = SphereToKartesian( r, (h - beta)/2.0, k - alpha);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    
    int elementid = 0;
    // Using triangle to sphere special map
    TPZVec<int64_t> topology(3);
    TPZVec<int64_t> topologyLine(2);
    
    
    // Side 0
    {
        topology[0] = 0;
        topology[1] = 1;
        topology[2] = 3;
        new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,topology, fbc1, *geomesh);
        elementid++;
    }
    
    
    // Side 1
    {
        topology[0] = 1;
        topology[1] = 2;
        topology[2] = 4;
        new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,topology, fbc2, *geomesh);
        elementid++;
    }
    
    
    // Side 2
    {
        topology[0] = 2;
        topology[1] = 0;
        topology[2] = 5;
        new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,topology, fbc3, *geomesh);
        elementid++;
    }
    
    
    // Create Geometrical Quad #1
    topology.Resize(3);
    topology[0] = 0;
    topology[1] = 1;
    topology[2] = 2;
    TPZGeoElRefPattern< pzgeom::TPZTriangleSphere< pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle > > > *TS = new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere< pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle > > > (elementid, topology, materialId, *geomesh);
//    TPZGeoElRefPattern< pzgeom::TPZTriangleSphere< > > *TS = new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere< > > (elementid, topology, materialId, *geomesh);
    TS->Geom().SetData(r,xc);
    
    //CONCLUINDO A CONSTRUCAO DA MALHA GEOMETRICA
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
    
    const unsigned int nel = geomesh->NElements();
    for (int i = 0 ; i < nel ; i++){
        TPZGeoEl *gel = geomesh->Element(i);
        if (!gel) {
            continue;
        }
        
        if (gel->Dimension()!=2) {
            continue;
        }
        const int npt = 3;
        TPZManVector<REAL,3> qsi(3,0.), x(3,0.);
        for (int iqsi = 0; iqsi < npt; iqsi++) {
            for (int ieta = 0; ieta < npt; ieta++) {
                qsi[0] = 0.0 + iqsi * 1./(npt-1);
                qsi[1] = 0.0 + ieta * (1.-qsi[0])/(npt-1);
                gel->X(qsi, x);
                //calcula a funcao e ve se zero
                //                cout << " qsi = " << qsi << " x = " << x << endl;
                //                cout << " r = " << sqrt(x[0]*x[0]+x[1]*x[1]) << endl;
                REAL R = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
//                cout << " R = " << R << endl;
                if (fabs(R-r)>1e-10) {
                    std::cout << "opa " << std::endl;
                }
                
            }
        }
    }
    
    std::ofstream outfile("OitavoCasca.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(geomesh, outfile, true);
    
    return geomesh;
    
    
}


// theta (0,pi) angulo que se inicia no polo norte. phi (0,2pi) o angulo no plano xy
TPZVec<REAL> LaplaceInSphere::SphereToKartesian(REAL r, REAL theta, REAL phi)
{
    TPZVec<REAL> xyz(3,0.0);
    xyz[0] = r*cos(phi)*sin(theta);
    xyz[1] = r*sin(theta)*sin(phi);
    xyz[2] = r*cos(theta);
    return xyz;
}
// theta (0,pi) angulo que se inicia no polo norte. phi (0,2pi) o angulo no plano xy
TPZVec<REAL> LaplaceInSphere::SphereToKartesian(TPZManVector<REAL> xc, REAL r, REAL theta, REAL phi)
{
    TPZVec<REAL> xyz(3,0.0);
    xyz[0] = xc[0] + r*sin(theta)*cos(phi);
    xyz[1] = xc[1] + r*sin(theta)*sin(phi);
    xyz[2] = xc[2] + r*cos(theta);
    return xyz;
}

TPZGeoMesh *LaplaceInSphere::MakeSphereFromQuadrilateral(int dimensao, bool triang, int ndiv)
{
    TPZGeoMesh * geomesh = new TPZGeoMesh;
    
    int nodes =  8;
    REAL radius = 1.0;
    geomesh->SetMaxNodeId(nodes-1);
    geomesh->NodeVec().Resize(nodes);
    TPZManVector<TPZGeoNode,4> Node(nodes);
    
    TPZManVector<int64_t,2> TopolQuad(4);
    TPZManVector<int64_t,1> TopolPoint(1);
    TPZManVector<int64_t,2> TopolLine(2);
    TPZManVector<REAL,3> coord(3,0.);
    TPZVec<REAL> xc(3,0.);
    
    REAL cphi = atan(sqrt(2.0));
    
    int nodeindex = 0;
    
    coord = ParametricSphere(radius,M_PI-cphi,M_PI/4.0);
    geomesh->NodeVec()[nodeindex].SetCoord(coord);
    geomesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    coord = ParametricSphere(radius,cphi,M_PI/4.0);
    geomesh->NodeVec()[nodeindex].SetCoord(coord);
    geomesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    coord = ParametricSphere(radius,cphi,-M_PI/4.0);
    geomesh->NodeVec()[nodeindex].SetCoord(coord);
    geomesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    coord = ParametricSphere(radius,M_PI-cphi,-M_PI/4.0);
    geomesh->NodeVec()[nodeindex].SetCoord(coord);
    geomesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    
    coord = ParametricSphere(radius,M_PI-cphi,3.0*M_PI/4.0);
    geomesh->NodeVec()[nodeindex].SetCoord(coord);
    geomesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    coord = ParametricSphere(radius,cphi,3.0*M_PI/4.0);
    geomesh->NodeVec()[nodeindex].SetCoord(coord);
    geomesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    coord = ParametricSphere(radius,cphi,-3.0*M_PI/4.0);
    geomesh->NodeVec()[nodeindex].SetCoord(coord);
    geomesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    coord = ParametricSphere(radius,M_PI-cphi,-3.0*M_PI/4.0);
    geomesh->NodeVec()[nodeindex].SetCoord(coord);
    geomesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    
    int id = 0;
    int matid = 1;
    
    TopolQuad[0] = 0;
    TopolQuad[1] = 1;
    TopolQuad[2] = 2;
    TopolQuad[3] = 3;
    TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > * quad1 = new TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > (id,TopolQuad,matid,*geomesh);
    quad1->Geom().SetData(radius, xc);
    id++;
    
    TopolQuad[0] = 4;
    TopolQuad[1] = 5;
    TopolQuad[2] = 6;
    TopolQuad[3] = 7;
    TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > * quad2 = new TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > (id,TopolQuad,matid,*geomesh);
    quad2->Geom().SetData(radius, xc);
    id++;
    
    TopolQuad[0] = 0;
    TopolQuad[1] = 4;
    TopolQuad[2] = 5;
    TopolQuad[3] = 1;
    TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > * quad3 = new TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > (id,TopolQuad,matid,*geomesh);
    quad3->Geom().SetData(radius, xc);
    id++;
    
    TopolQuad[0] = 3;
    TopolQuad[1] = 7;
    TopolQuad[2] = 6;
    TopolQuad[3] = 2;
    TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > * quad4 = new TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > (id,TopolQuad,matid,*geomesh);
    quad4->Geom().SetData(radius, xc);
    id++;
    
    TopolQuad[0] = 0;
    TopolQuad[1] = 4;
    TopolQuad[2] = 7;
    TopolQuad[3] = 3;
    TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > * quad5 = new TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > (id,TopolQuad,matid,*geomesh);
    quad5->Geom().SetData(radius, xc);
    id++;
    
    TopolQuad[0] = 1;
    TopolQuad[1] = 5;
    TopolQuad[2] = 6;
    TopolQuad[3] = 2;
    TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > * quad6 = new TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > (id,TopolQuad,matid,*geomesh);
    quad6->Geom().SetData(radius, xc);
    id++;
    
    TopolLine[0] = 0;
    TopolLine[1] = 4;
    new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (id,TopolLine, fbc0, *geomesh);
    id++;
    
    TopolLine[0] = 0;
    TopolLine[1] = 4;
    new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (id,TopolLine, fbc1, *geomesh);
    id++;

    
//    TopolLine[0] = 0;
//    TopolLine[1] = 4;
//    TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoLinear> > * arcleft = new TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoLinear> > (id,TopolLine, fbc0, *geomesh);
//    arcleft->Geom().SetData(radius, xc);
//    id++;
//    
//    TopolLine[0] = 0;
//    TopolLine[1] = 4;
//    TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoLinear> > * arcrigth = new TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoLinear> > (id,TopolLine, fbc1, *geomesh);
//    arcrigth->Geom().SetData(radius, xc);
//    id++;
    
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
    this->RotateGeomesh(geomesh, angle, axis);
    
//    TPZCheckGeom * checkgeo = new TPZCheckGeom(geomesh);
//    int badmeshQ = checkgeo->PerformCheck();
//    
//    if (badmeshQ) {
//        DebugStop();
//    }
    
    std::ofstream outfile("NiceSphere.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(geomesh, outfile, true);
    
    return geomesh;
}

TPZManVector<REAL,3> LaplaceInSphere::ParametricSphere(REAL radius,REAL phi,REAL theta)
{
    TPZManVector<REAL,3> xcoor(3,0.0);
    xcoor[0] = radius * cos(theta) * sin(phi);
    xcoor[1] = radius * sin(theta) * sin(phi);
    xcoor[2] = radius * cos(phi) ;
    return xcoor;
}

// um anel de quadrilateros e topo que pode ser quadilatero ou triangulo
TPZGeoMesh *LaplaceInSphere::GMeshSphericalShell(int dimensao, bool triang, int ndiv)
{
    
    bool topoComTriangulos = triang;
    
    ///INSTANCIACAO DA MALHA GEOMETRICA
    TPZGeoMesh * geomesh = new TPZGeoMesh;
    
    /// Materiais
    int64_t materialId = fmatId;
    
    int nnodes = 13;//quantidade de nos da malha geometrica
    geomesh->NodeVec().Resize(nnodes);
    geomesh->SetDimension(2);
    
    
    
    ///INICIALIZACAO DA MALHA GEOMETRICA PELA INSTANCIACAO E INICIALIZACAO DOS NOS
    TPZGeoNode node;
    TPZVec<REAL> coord(3,0.);
    const REAL r = 1.;

    TPZManVector<REAL,3> xc(3,0.);
    xc[0] = 0.;
    xc[1] = 0.;
    xc[2] = 0.;
    
    int Axis = 3;
    REAL angrot = 0.0;
    
    int id = 0;
    //no 0
    coord = SphereToKartesian(xc, r, M_PI/2.0, 0.0);
//    coord[0] = xc[0] + r*sin( M_PI/2.0)*cos(0.0);
//    coord[1] = xc[1] + r*sin( M_PI/2.0)*sin(0.0);
//    coord[2] = xc[2] + r*cos( M_PI/2.0);
    //cout<<"coord = "<< coord<<endl;
    tools::RotateNode(coord, angrot, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 1
    coord = SphereToKartesian(r, M_PI/2.0, M_PI/2.0);
    //cout<<"coord = "<< coord<<endl;
    tools::RotateNode(coord, angrot, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 2
    coord = SphereToKartesian(r, M_PI/2.0, M_PI);
    //cout<<"coord = "<< coord<<endl;
    tools::RotateNode(coord, angrot, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 3
    coord = SphereToKartesian(r, M_PI/2.0, 3.0*M_PI/2.0);
    //cout<<"coord = "<< coord<<endl;
    tools::RotateNode(coord, angrot, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 4
    coord = SphereToKartesian(r, M_PI/4.0, 0.0);
    //cout<<"coord = "<< coord<<endl;
    tools::RotateNode(coord, angrot, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 5
    coord = SphereToKartesian(r, M_PI/4.0, M_PI/2.0);
    //cout<<"coord = "<< coord<<endl;
    tools::RotateNode(coord, angrot, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 6
    coord = SphereToKartesian(r, M_PI/4.0, M_PI);
    //cout<<"coord = "<< coord<<endl;
    tools::RotateNode(coord, angrot, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 7
    coord = SphereToKartesian(r, M_PI/4.0, 3.0*M_PI/2.0);
    //cout<<"coord = "<< coord<<endl;
    tools::RotateNode(coord, angrot, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 8
    coord = SphereToKartesian(r, 0.0, 0.0);
    //cout<<"coord = "<< coord<<endl;
    tools::RotateNode(coord, angrot, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    //no 9
    coord = SphereToKartesian(r, M_PI/2.0, M_PI/4.0);
    tools::RotateNode(coord, angrot, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 10
    coord = SphereToKartesian(r, M_PI/2.0, 3.0*M_PI/4.0);
    tools::RotateNode(coord, angrot, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 11
    coord = SphereToKartesian(r, M_PI/2.0, 5.0*M_PI/4.0);
    tools::RotateNode(coord, angrot, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 12
    coord = SphereToKartesian(r, M_PI/2.0, 7.0*M_PI/4.0);
    tools::RotateNode(coord, angrot, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    int elementid = 0;
    
    TPZVec<int64_t> topology(4);
    
    // El 0
    topology[0] = 0;
    topology[1] = 1;
    topology[2] = 5;
    topology[3] = 4;
    {
        TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > * SphereRingQ =
        new TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > (elementid, topology,materialId,*geomesh);
        SphereRingQ->Geom().SetData(r,xc);
        elementid++;
    }
    // El 1
    topology[0] = 1;
    topology[1] = 2;
    topology[2] = 6;
    topology[3] = 5;
    {
        TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > * SphereRingQ =
        new TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > (elementid, topology,materialId,*geomesh);
        SphereRingQ->Geom().SetData(r,xc);
        elementid++;
    }
    topology.resize(4);
    
    // El 0
    topology[0] = 2;
    topology[1] = 3;
    topology[2] = 7;
    topology[3] = 6;
    {
        TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > * SphereRingQ =
        new TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > (elementid, topology,materialId,*geomesh);
        SphereRingQ->Geom().SetData(r,xc);
        elementid++;
    }
    // El 1
    topology[0] = 3;
    topology[1] = 0;
    topology[2] = 4;
    topology[3] = 7;
    {
        TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > * SphereRingQ =
        new TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > (elementid, topology,materialId,*geomesh);
        SphereRingQ->Geom().SetData(r,xc);
        elementid++;
    }
    
    if (topoComTriangulos)
    {
        
        topology.resize(3);
        
        // El 0
        topology[0] = 4;
        topology[1] = 5;
        topology[2] = 8;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > * SphereRingT1 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > (elementid, topology,materialId,*geomesh);
            SphereRingT1->Geom().SetData(r,xc);
            elementid++;
        }
        // El 1
        topology[0] = 5;
        topology[1] = 6;
        topology[2] = 8;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > * SphereRingT2 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > (elementid, topology,materialId,*geomesh);
            SphereRingT2->Geom().SetData(r,xc);
            elementid++;
        }
        // El 2
        topology[0] = 6;
        topology[1] = 7;
        topology[2] = 8;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > * SphereRingT1 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > (elementid, topology,materialId,*geomesh);
            SphereRingT1->Geom().SetData(r,xc);
            elementid++;
        }
        // El 3
        topology[0] = 7;
        topology[1] = 4;
        topology[2] = 8;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > * SphereRingT2 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > (elementid, topology,materialId,*geomesh);
            SphereRingT2->Geom().SetData(r,xc);
            elementid++;
        }
        
    }
    else
    {
        
//        topology.resize(4);
//        
//        // El 0
//        topology[0] = 4;
//        topology[1] = 5;
//        topology[2] = 6;
//        topology[3] = 7;
//        {
//            TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > * SphereRingQ =
//            new TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > (elementid, topology,materialId,*geomesh);
//            SphereRingQ->Geom().SetData(r,xc);
//            elementid++;
//        }

        topology.resize(2);
        
        topology[0] = 4;
        topology[1] = 5;
        //topology[2] = 9;
        {
            //TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc =
            //new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,topology, fbc1 /** fbc3 */, *geomesh);
            new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, fbc1 /** fbc3 */, *geomesh);
            //arc->Geom().Print(std::cout);
        }
        elementid++;
        topology[0] = 5;
        topology[1] = 6;
        //topology[2] = 9;
        {
            //TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc =
            //new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,topology, fbc1 /** fbc3 */, *geomesh);
            new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, fbc1 /** fbc3 */, *geomesh);
            //arc->Geom().Print(std::cout);
        }
        elementid++;
        topology[0] = 6;
        topology[1] = 7;
        //topology[2] = 9;
        {
            //TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc =
            //new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,topology, fbc1 /** fbc3 */, *geomesh);
            new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, fbc1 /** fbc3 */, *geomesh);
            //arc->Geom().Print(std::cout);
        }
        elementid++;
        topology[0] = 7;
        topology[1] = 4;
        //topology[2] = 9;
        {
            //TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc =
            //new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,topology, fbc1 /** fbc3 */, *geomesh);
            new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, fbc1 /** fbc3 */, *geomesh);
            //arc->Geom().Print(std::cout);
        }
        elementid++;
    }
    
        
    // El linha
    // Definition of Arc coordenates
    topology.resize(2);
    // Create Geometrical Arc #0
    topology[0] = 0;
    topology[1] = 1;
    //topology[2] = 9;
    {
        //TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc =
        //new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,topology, fbc1 /** fbc3 */, *geomesh);
        new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, fbc1 /** fbc3 */, *geomesh);
        //arc->Geom().Print(std::cout);
    }
    elementid++;
    
    // Create Geometrical Arc #1
    topology[0] = 1;
    topology[1] = 2;
    //topology[2] = 10;
    {
        //TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc =
        //new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,topology, fbc1 /** fbc4 */, *geomesh);
        new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, fbc1 /** fbc3 */, *geomesh);
        // arc->Geom().Print(std::cout);
    }
    elementid++;
    // Create Geometrical Arc #1
    topology[0] = 2;
    topology[1] = 3;
    //topology[2] = 11;
    {
        //TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc =
        //new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,topology, fbc1 /** fbc4 */, *geomesh);
        new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, fbc1 /** fbc3 */, *geomesh);
        // arc->Geom().Print(std::cout);
    }
    elementid++;
    // Create Geometrical Arc #1
    topology[0] = 3;
    topology[1] = 0;
    {
        //TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc =
        //new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,topology, fbc1 /** fbc4 */, *geomesh);
        new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, fbc1 /** fbc3 */, *geomesh);
        // arc->Geom().Print(std::cout);
    }
    elementid++;
    
    
    
    //CONCLUINDO A CONSTRUCAO DA MALHA GEOMETRICA
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
    
    std::ofstream outfile("malhaCascaesfera.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(geomesh, outfile, true);
    
    return geomesh;
    
    
}

// anel com quadrilateros e topo de quadrilateros ou triangulos
TPZGeoMesh *LaplaceInSphere::GMeshSphericalShellBlendQ(int dimensao, bool triang, int ndiv)
{
    
    bool topoComTriangulos = triang;
    
    ///INSTANCIACAO DA MALHA GEOMETRICA
    TPZGeoMesh * geomesh = new TPZGeoMesh;
    
    /// Materiais
    int64_t materialId = fmatId;
    
    int nnodes = 17;//quantidade de nos da malha geometrica
    geomesh->NodeVec().Resize(nnodes);
    geomesh->SetDimension(2);
    
    
    
    ///INICIALIZACAO DA MALHA GEOMETRICA PELA INSTANCIACAO E INICIALIZACAO DOS NOS
    TPZGeoNode node;
    TPZVec<REAL> coord(3,0.);
    const REAL r = 1.;
    
    TPZManVector<REAL,3> xc(3,0.);
    xc[0] = 0.;
    xc[1] = 0.;
    xc[2] = 0.;
    
    REAL base = M_PI/2.0;
    REAL topo = M_PI/10.0;
    
    int Axis = 3;
    REAL angrot = 0.0;
    
    int id = 0;
    //no 0
    coord = SphereToKartesian(xc, r, base, 0.0);
    tools::RotateNode(coord, angrot, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 1
    coord = SphereToKartesian(r, base, M_PI/2.0);
    tools::RotateNode(coord, angrot, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 2
    coord = SphereToKartesian(r, base, M_PI);
    //cout<<"coord = "<< coord<<endl;
    tools::RotateNode(coord, angrot, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 3
    coord = SphereToKartesian(r, base, 3.0*M_PI/2.0);
    //cout<<"coord = "<< coord<<endl;
    tools::RotateNode(coord, angrot, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 4
    coord = SphereToKartesian(r, topo, 0.0);
    //cout<<"coord = "<< coord<<endl;
    tools::RotateNode(coord, angrot, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 5
    coord = SphereToKartesian(r, topo, M_PI/2.0);
    //cout<<"coord = "<< coord<<endl;
    tools::RotateNode(coord, angrot, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 6
    coord = SphereToKartesian(r, topo, M_PI);
    //cout<<"coord = "<< coord<<endl;
    tools::RotateNode(coord, angrot, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 7
    coord = SphereToKartesian(r, topo, 3.0*M_PI/2.0);
    //cout<<"coord = "<< coord<<endl;
    tools::RotateNode(coord, angrot, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 8
    coord = SphereToKartesian(r, 0.0, 0.0);
    //cout<<"coord = "<< coord<<endl;
    tools::RotateNode(coord, angrot, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    //no 9
    coord = SphereToKartesian(r, base, M_PI/4.0);
    tools::RotateNode(coord, angrot, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 10
    coord = SphereToKartesian(r, base, 3.0*M_PI/4.0);
    tools::RotateNode(coord, angrot, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 11
    coord = SphereToKartesian(r, base, 5.0*M_PI/4.0);
    tools::RotateNode(coord, angrot, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 12
    coord = SphereToKartesian(r, base, 7.0*M_PI/4.0);
    tools::RotateNode(coord, angrot, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 13
    coord = SphereToKartesian(r, topo, M_PI/4.0);
    tools::RotateNode(coord, angrot, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 14
    coord = SphereToKartesian(r, topo, 3.0*M_PI/4.0);
    tools::RotateNode(coord, angrot, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 15
    coord = SphereToKartesian(r, topo, 5.0*M_PI/4.0);
    tools::RotateNode(coord, angrot, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    //no 16
    coord = SphereToKartesian(r, topo, 7.0*M_PI/4.0);
    tools::RotateNode(coord, angrot, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    id++;
    
    int elementid = 0;
    
    TPZVec<int64_t> topology(4);
    
    
    
    topology.resize(3);
    
    topology[0] = 0;
    topology[1] = 1;
    topology[2] = 9;
    {
        new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,topology, fbc1 /** fbc3 */, *geomesh);
    }
    elementid++;
    
    topology[0] = 1;
    topology[1] = 2;
    topology[2] = 10;
    {
        new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,topology, fbc1 /** fbc3 */, *geomesh);
    }
    elementid++;
    
    topology[0] = 2;
    topology[1] = 3;
    topology[2] = 11;
    {
        new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,topology, fbc1 /** fbc3 */, *geomesh);
    }
    elementid++;
    
    topology[0] = 3;
    topology[1] = 0;
    topology[2] = 12;
    {
        new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,topology, fbc1 /** fbc3 */, *geomesh);
    }
    elementid++;
    
    
    int fnu = fbc1;//-10;
    
    topology[0] = 4;
    topology[1] = 5;
    topology[2] = 13;
    {
        new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,topology, fnu  /** fbc1 */, *geomesh);
    }
    elementid++;
    
    topology[0] = 5;
    topology[1] = 6;
    topology[2] = 14;
    {
        new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,topology, fnu  /** fbc1 */, *geomesh);
    }
    elementid++;
    
    topology[0] = 6;
    topology[1] = 7;
    topology[2] = 15;
    {
        new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,topology, fnu  /** fbc1 */, *geomesh);
    }
    elementid++;
    
    topology[0] = 7;
    topology[1] = 4;
    topology[2] = 16;
    {
        new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,topology, fnu /** fbc1 */, *geomesh);
    }
    elementid++;
    
    
    
    
    topology.resize(4);
    // El 0
    topology[0] = 0;
    topology[1] = 1;
    topology[2] = 5;
    topology[3] = 4;
    {
        TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> > > * SphereRingQ =
        new TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> > > (elementid, topology,materialId,*geomesh);
        SphereRingQ->Geom().SetData(r,xc);
        elementid++;
    }
    // El 1
    topology[0] = 1;
    topology[1] = 2;
    topology[2] = 6;
    topology[3] = 5;
    {
        TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> > > * SphereRingQ =
        new TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> > > (elementid, topology,materialId,*geomesh);
        SphereRingQ->Geom().SetData(r,xc);
        elementid++;
    }
    topology.resize(4);
    
    // El 0
    topology[0] = 2;
    topology[1] = 3;
    topology[2] = 7;
    topology[3] = 6;
    {
        TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> > > * SphereRingQ =
        new TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> > > (elementid, topology,materialId,*geomesh);
        SphereRingQ->Geom().SetData(r,xc);
        elementid++;
    }
    // El 1
    topology[0] = 3;
    topology[1] = 0;
    topology[2] = 4;
    topology[3] = 7;
    {
        TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> > > * SphereRingQ =
        new TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> > > (elementid, topology,materialId,*geomesh);
        SphereRingQ->Geom().SetData(r,xc);
        elementid++;
    }
    
    if (topoComTriangulos)
    {
        
//        topology.resize(3);
//        
//        // El 0
//        topology[0] = 4;
//        topology[1] = 5;
//        topology[2] = 8;
//        {
//            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere< pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle > > > * SphereRingT1 =
//            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere< pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle > > > (elementid, topology,materialId,*geomesh);
//            SphereRingT1->Geom().SetData(r,xc);
//            elementid++;
//        }
//        // El 1
//        topology[0] = 5;
//        topology[1] = 6;
//        topology[2] = 8;
//        {
//            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere< pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle > > > * SphereRingT2 =
//            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere< pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle > > > (elementid, topology,materialId,*geomesh);
//            SphereRingT2->Geom().SetData(r,xc);
//            elementid++;
//        }
//        // El 2
//        topology[0] = 6;
//        topology[1] = 7;
//        topology[2] = 8;
//        {
//            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere< pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle > > > * SphereRingT1 =
//            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere< pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle > > > (elementid, topology,materialId,*geomesh);
//            SphereRingT1->Geom().SetData(r,xc);
//            elementid++;
//        }
//        // El 3
//        topology[0] = 7;
//        topology[1] = 4;
//        topology[2] = 8;
//        {
//            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere< pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle > > > * SphereRingT2 =
//            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere< pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle > > > (elementid, topology,materialId,*geomesh);
//            SphereRingT2->Geom().SetData(r,xc);
//            elementid++;
//        }
        
    }
    else
    {
        
//        topology.resize(4);
//
//        // El 0
//        topology[0] = 4;
//        topology[1] = 5;
//        topology[2] = 6;
//        topology[3] = 7;
//        {
//            TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> > > * SphereRingQ =
//            new TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> > > (elementid, topology,materialId,*geomesh);
//            SphereRingQ->Geom().SetData(r,xc);
//            elementid++;
//        }
        

    }
    
    
    
    //CONCLUINDO A CONSTRUCAO DA MALHA GEOMETRICA
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
    
    std::ofstream outfile("malhaCascaesfera.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(geomesh, outfile, true);
    
    return geomesh;
    
    
}

TPZGeoMesh *LaplaceInSphere::GMeshSphericalShell2(int dimensao, bool triang, int ndiv)
{
    
    
    bool ftriangulo = triang;
    
    ///INSTANCIACAO DA MALHA GEOMETRICA
    TPZGeoMesh * geomesh = new TPZGeoMesh;
    
    /// Materiais
    int64_t materialId = fmatId;
    
    int nnodes = 37;//quantidade de nos da malha geometrica
    geomesh->NodeVec().Resize(nnodes);
    geomesh->SetDimension(2);
    
    
    
    ///INICIALIZACAO DA MALHA GEOMETRICA PELA INSTANCIACAO E INICIALIZACAO DOS NOS
    TPZGeoNode node;
    TPZVec<REAL> coord(3,0.);
    const REAL r = 1.;

    TPZManVector<REAL,3> xc(3,0.);
    xc[0] = 0.;
    xc[1] = 0.;
    xc[2] = 0.;
    
    int Axis = 3;
    REAL theta = 0.0;
    
    int sind = (nnodes+1)/3;
    
    TPZManVector<int> indinf(sind);
    TPZManVector<int> indmid(sind);
    TPZManVector<int> indsup(sind);
    int cont = 0;
    
    int id = 0;
    for (int j = 0; j < 3; j++) {
        for (int i = 0; i < 12; i++) {
            //no id
            coord[0] = xc[0] + r*cos(i*M_PI/6.0)*sin(M_PI/2.0 - j*M_PI/6.0);
            coord[1] = xc[1] + r*sin(i*M_PI/6.0)*sin(M_PI/2.0 - j*M_PI/6.0);
            coord[2] = xc[2] + r*cos(M_PI/2.0 - j*M_PI/6.0);
            tools::RotateNode(coord, theta, Axis);
            node.SetNodeId(id);
            node.SetCoord(coord);
            geomesh->NodeVec()[id] = node;
            if (j==0) {
                indinf[cont] = id;
                cont++;
            }
            else if (j==1)
            {
                indmid[cont] = id;
                cont++;
            }
            else
            {
                indsup[cont] = id;
                cont++;
            }
            id++;
            
        }
        cont = 0;
    }
    
    //no polo norte
    coord[0] = xc[0] + 0.;
    coord[1] = xc[1] + 0.;
    coord[2] = xc[2] + r;
    tools::RotateNode(coord, theta, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    //id++;
    
    int elementid = 0;
    // Using triangle to sphere special map
    TPZVec<int64_t> topology(4);
    
    
    if (ftriangulo)
    {
        
        topology.resize(3);
        
        // El 0
        topology[0] = 0;
        topology[1] = 1;
        topology[2] = 5;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > * SphereRingT1 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > (elementid, topology,materialId,*geomesh);
            SphereRingT1->Geom().SetData(r,xc);
            elementid++;
        }
        // El 1
        topology[0] = 0;
        topology[1] = 5;
        topology[2] = 4;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > * SphereRingT2 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > (elementid, topology,materialId,*geomesh);
            SphereRingT2->Geom().SetData(r,xc);
            elementid++;
        }
        // El 2
        topology[0] = 1;
        topology[1] = 2;
        topology[2] = 6;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > * SphereRingT1 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > (elementid, topology,materialId,*geomesh);
            SphereRingT1->Geom().SetData(r,xc);
            elementid++;
        }
        // El 3
        topology[0] = 1;
        topology[1] = 6;
        topology[2] = 5;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > * SphereRingT2 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > (elementid, topology,materialId,*geomesh);
            SphereRingT2->Geom().SetData(r,xc);
            elementid++;
        }
        // El 4
        topology[0] = 2;
        topology[1] = 3;
        topology[2] = 7;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > * SphereRingT1 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > (elementid, topology,materialId,*geomesh);
            SphereRingT1->Geom().SetData(r,xc);
            elementid++;
        }
        // El 5
        topology[0] = 2;
        topology[1] = 7;
        topology[2] = 6;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > * SphereRingT2 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > (elementid, topology,materialId,*geomesh);
            SphereRingT2->Geom().SetData(r,xc);
            elementid++;
        }
        // El 6
        topology[0] = 3;
        topology[1] = 0;
        topology[2] = 4;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > * SphereRingT1 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > (elementid, topology,materialId,*geomesh);
            SphereRingT1->Geom().SetData(r,xc);
            elementid++;
        }
        // El 7
        topology[0] = 3;
        topology[1] = 4;
        topology[2] = 7;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > * SphereRingT2 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > (elementid, topology,materialId,*geomesh);
            SphereRingT2->Geom().SetData(r,xc);
            elementid++;
        }
    }
    else
    {
        
        topology.resize(4);
        
        for (int nelinf = 0; nelinf < 12; nelinf++) {
            // El nel
            
            topology[0] = indinf[nelinf];
            topology[1] = indinf[(nelinf+1)%12];
            topology[2] = indmid[(nelinf+1)%12];
            topology[3] = indmid[nelinf];
            {
                TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > * SphereRingQ =
                new TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > (elementid, topology,materialId,*geomesh);
                SphereRingQ->Geom().SetData(r,xc);
                elementid++;
            }
        }
        
        for (int nelmid = 0; nelmid < 12; nelmid++) {
            // El nel
            topology[0] = indmid[nelmid];
            topology[1] = indmid[(nelmid+1)%12];
            topology[2] = indsup[(nelmid+1)%12];
            topology[3] = indsup[nelmid];
            {
                TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > * SphereRingQ =
                new TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > (elementid, topology,materialId,*geomesh);
                SphereRingQ->Geom().SetData(r,xc);
                elementid++;
            }
        }
        
        topology.resize(3);
        
        //        for (int nelsup = 0; nelsup < 12; nelsup++) {
        //            // El nel
        ////            int a = indsup[nelsup];
        ////            int b = indsup[(nelsup+1)%12];
        ////            int c = polo;
        //            topology[0] = indsup[nelsup];
        //            topology[1] = indsup[(nelsup+1)%12];
        //            topology[2] = polo;
        //            {
        //                TPZGeoElRefPattern< pzgeom::TPZTriangleSphere > * SphereRingQ =
        //                new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere > (elementid, topology,materialId,*geomesh);
        //                SphereRingQ->Geom().SetData(r,xc);
        //                elementid++;
        //            }
        //        }
        
        
    }
    
    // El linha
    // Definition of Arc
    topology.resize(2);
    
    for (int nelinf = 0; nelinf < 12; nelinf++) {
        // El nel
        topology[0] = indinf[nelinf];
        topology[1] = indinf[(nelinf+1)%12];
        {
            //TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc =
            new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, fbc1 /** fbc3 */, *geomesh);
            //arc->Geom().Print(std::cout);
            elementid++;
        }
    }
    
    // arcos do topo
    for (int nelsup = 0; nelsup < 12; nelsup++) {
        // El nel
        topology[0] = indsup[nelsup];
        topology[1] = indsup[(nelsup+1)%12];
        {
            //TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc =
            new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, fbc1 /** fbc3 */, *geomesh);
            //arc->Geom().Print(std::cout);
            elementid++;
        }
    }
    
    
    
    //    // Create Geometrical Arc #0
    //    topology[0] = 0;
    //    topology[1] = 1;
    //    {
    //        TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc = new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, fbc1 /** fbc3 */, *geomesh);
    //        //arc->Geom().Print(std::cout);
    //    }
    //    elementid++;
    //
    //    // Create Geometrical Arc #1
    //    topology[0] = 1;
    //    topology[1] = 2;
    //    {
    //        TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc = new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, fbc1 /** fbc4 */, *geomesh);
    //        // arc->Geom().Print(std::cout);
    //    }
    //    elementid++;
    //
    //    // Create Geometrical Arc #2
    //    topology[0] = 2;
    //    topology[1] = 3;
    //    {
    //        TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc = new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, fbc1 /** fbc1 */, *geomesh);
    //        //arc->Geom().Print(std::cout);
    //    }
    //    elementid++;
    //
    //    // Create Geometrical Arc #3
    //    topology[0] = 3;
    //    topology[1] = 0;
    //    {
    //        TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc = new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, fbc1 /** fbc2 */, *geomesh);
    //        //arc->Geom().Print(std::cout);
    //    }
    
    
    //CONCLUINDO A CONSTRUCAO DA MALHA GEOMETRICA
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
    
    std::ofstream outfile("malhaCascaesfera.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(geomesh, outfile, true);
    
    return geomesh;
    
    
}

TPZGeoMesh *LaplaceInSphere::GMeshSliceSphericalShell(int dimensao, bool triang, int ndiv)
{
    
    int nfatias = 1;
    const REAL r = 1.;
    // centro da esfera
    TPZManVector<REAL,3> xc(3,0.);
    xc[0] = 0.;
    xc[1] = 0.;
    xc[2] = 0.;
    
    // para rotacao da malha, se quiser
    int Axis = 3;
    REAL theta = 0.0;
    
    
    bool ftriangulo = triang;
    
    ///INSTANCIACAO DA MALHA GEOMETRICA
    TPZGeoMesh * geomesh = new TPZGeoMesh;
    
    /// Materiais
    int64_t materialId = fmatId;
    
    int nnodes = nfatias == 12 ? 37 :(nfatias+1)*3+1;//quantidade de nos da malha geometrica
    geomesh->NodeVec().Resize(nnodes);
    geomesh->SetDimension(2);
    
    
    
    ///INICIALIZACAO DA MALHA GEOMETRICA PELA INSTANCIACAO E INICIALIZACAO DOS NOS
    TPZGeoNode node;
    TPZVec<REAL> coord(3,0.);
    
    // quantidade de indices
    int sind = (nnodes+1)/3;
    
    TPZManVector<int> indinf(sind);
    TPZManVector<int> indmid(sind);
    TPZManVector<int> indsup(sind);
    int cont = 0;
    int npts = nfatias == 12 ? 12 : nfatias + 1;
    
    int id = 0;
    for (int j = 0; j < 3; j++) {
        for (int i = 0; i < npts; i++) {
            //no id
            coord[0] = xc[0] + r*cos(i*M_PI/6.0)*sin(M_PI/2.0 - j*M_PI/6.0);
            coord[1] = xc[1] + r*sin(i*M_PI/6.0)*sin(M_PI/2.0 - j*M_PI/6.0);
            coord[2] = xc[2] + r*cos(M_PI/2.0 - j*M_PI/6.0);
            tools::RotateNode(coord, theta, Axis);
            node.SetNodeId(id);
            node.SetCoord(coord);
            geomesh->NodeVec()[id] = node;
            if (j==0) {
                indinf[cont] = id;
                cont++;
            }
            else if (j==1)
            {
                indmid[cont] = id;
                cont++;
            }
            else
            {
                indsup[cont] = id;
                cont++;
            }
            id++;
            
        }
        cont = 0;
    }
    
    //no polo norte
    coord[0] = xc[0] + 0.;
    coord[1] = xc[1] + 0.;
    coord[2] = xc[2] + r;
    tools::RotateNode(coord, theta, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    geomesh->NodeVec()[id] = node;
    //int polo = id;
    //id++;
    
    int elementid = 0;
    // Using triangle to sphere special map
    TPZVec<int64_t> topology(4);
    
    
    if (ftriangulo)
    {
        DebugStop();
        topology.resize(3);
        
        // El 0
        topology[0] = 0;
        topology[1] = 1;
        topology[2] = 5;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > * SphereRingT1 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > (elementid, topology,materialId,*geomesh);
            SphereRingT1->Geom().SetData(r,xc);
            elementid++;
        }
        // El 1
        topology[0] = 0;
        topology[1] = 5;
        topology[2] = 4;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > * SphereRingT2 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > (elementid, topology,materialId,*geomesh);
            SphereRingT2->Geom().SetData(r,xc);
            elementid++;
        }
        // El 2
        topology[0] = 1;
        topology[1] = 2;
        topology[2] = 6;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > * SphereRingT1 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > (elementid, topology,materialId,*geomesh);
            SphereRingT1->Geom().SetData(r,xc);
            elementid++;
        }
        // El 3
        topology[0] = 1;
        topology[1] = 6;
        topology[2] = 5;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > * SphereRingT2 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > (elementid, topology,materialId,*geomesh);
            SphereRingT2->Geom().SetData(r,xc);
            elementid++;
        }
        // El 4
        topology[0] = 2;
        topology[1] = 3;
        topology[2] = 7;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > * SphereRingT1 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > (elementid, topology,materialId,*geomesh);
            SphereRingT1->Geom().SetData(r,xc);
            elementid++;
        }
        // El 5
        topology[0] = 2;
        topology[1] = 7;
        topology[2] = 6;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > * SphereRingT2 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > (elementid, topology,materialId,*geomesh);
            SphereRingT2->Geom().SetData(r,xc);
            elementid++;
        }
        // El 6
        topology[0] = 3;
        topology[1] = 0;
        topology[2] = 4;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > * SphereRingT1 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > (elementid, topology,materialId,*geomesh);
            SphereRingT1->Geom().SetData(r,xc);
            elementid++;
        }
        // El 7
        topology[0] = 3;
        topology[1] = 4;
        topology[2] = 7;
        {
            TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > * SphereRingT2 =
            new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<> > (elementid, topology,materialId,*geomesh);
            SphereRingT2->Geom().SetData(r,xc);
            elementid++;
        }
    }
    else
    {
        
        topology.resize(4);
        
        for (int nelinf = 0; nelinf < nfatias; nelinf++) {
            // El nel
//            int a = indinf[nelinf];
//            int b = indinf[(nelinf+1)%npts];
//            int c = indmid[(nelinf+1)%npts];
//            int d = indmid[nelinf];
            
            topology[0] = indinf[nelinf];
            topology[1] = indinf[(nelinf+1)%npts];
            topology[2] = indmid[(nelinf+1)%npts];
            topology[3] = indmid[nelinf];
            {
                TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > * SphereRingQ =
                new TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > (elementid, topology,materialId,*geomesh);
                SphereRingQ->Geom().SetData(r,xc);
                elementid++;
            }
        }
        
        for (int nelmid = 0; nelmid < nfatias; nelmid++) {
            // El nel
//            int a = indmid[nelmid];
//            int b = indmid[(nelmid+1)%npts];
//            int c = indsup[(nelmid+1)%npts];
//            int d = indsup[nelmid];
            topology[0] = indmid[nelmid];
            topology[1] = indmid[(nelmid+1)%npts];
            topology[2] = indsup[(nelmid+1)%npts];
            topology[3] = indsup[nelmid];
            {
                TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > * SphereRingQ =
                new TPZGeoElRefPattern< pzgeom::TPZQuadSphere<> > (elementid, topology,materialId,*geomesh);
                SphereRingQ->Geom().SetData(r,xc);
                elementid++;
            }
        }
        
        topology.resize(3);
        
        //        for (int nelsup = 0; nelsup < nfatias; nelsup++) {
        //            // El nel
        ////            int a = indsup[nelsup];
        ////            int b = indsup[(nelsup+1)%npts];
        ////            int c = polo;
        //            topology[0] = indsup[nelsup];
        //            topology[1] = indsup[(nelsup+1)%npts];
        //            topology[2] = polo;
        //            {
        //                TPZGeoElRefPattern< pzgeom::TPZTriangleSphere > * SphereRingQ =
        //                new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere > (elementid, topology,materialId,*geomesh);
        //                SphereRingQ->Geom().SetData(r,xc);
        //                elementid++;
        //            }
        //        }
        
        
    }
    
    // El linha
    // Definition of Arc
    topology.resize(2);
    
    for (int nelinf = 0; nelinf < nfatias; nelinf++) {
        // El nel
        topology[0] = indinf[nelinf];
        topology[1] = indinf[(nelinf+1)%npts];
        {
            //TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc =
            new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, fbc1 /** fbc3 */, *geomesh);
            //arc->Geom().Print(std::cout);
            elementid++;
        }
    }
    
    // arcos do topo
    for (int nelsup = 0; nelsup < nfatias; nelsup++) {
        // El nel
        topology[0] = indsup[nelsup];
        topology[1] = indsup[(nelsup+1)%npts];
        {
            //TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc =
            new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, fbc1 /** fbc3 */, *geomesh);
            //arc->Geom().Print(std::cout);
            elementid++;
        }
    }
    
    
    // condicoes lado
    if (nfatias<12) {
        topology[0] = 1;
        topology[1] = 3;
        {
            //TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc =
            new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, fbc3 /** fbc3 */, *geomesh);
            //arc->Geom().Print(std::cout);
            elementid++;
        }
        topology[0] = 3;
        topology[1] = 5;
        {
            //TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc =
            new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, fbc3 /** fbc3 */, *geomesh);
            //arc->Geom().Print(std::cout);
            elementid++;
        }
        topology[0] = 0;
        topology[1] = 2;
        {
            //TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc =
            new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, fbc3 /** fbc3 */, *geomesh);
            //arc->Geom().Print(std::cout);
            elementid++;
        }
        topology[0] = 2 ;
        topology[1] = 4;
        {
            //TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc =
            new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, fbc3 /** fbc3 */, *geomesh);
            //arc->Geom().Print(std::cout);
            elementid++;
        }
    }
    
    
    
    //    // Create Geometrical Arc #0
    //    topology[0] = 0;
    //    topology[1] = 1;
    //    {
    //        TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc = new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, fbc1 /** fbc3 */, *geomesh);
    //        //arc->Geom().Print(std::cout);
    //    }
    //    elementid++;
    //
    //    // Create Geometrical Arc #1
    //    topology[0] = 1;
    //    topology[1] = 2;
    //    {
    //        TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc = new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, fbc1 /** fbc4 */, *geomesh);
    //        // arc->Geom().Print(std::cout);
    //    }
    //    elementid++;
    //
    //    // Create Geometrical Arc #2
    //    topology[0] = 2;
    //    topology[1] = 3;
    //    {
    //        TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc = new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, fbc1 /** fbc1 */, *geomesh);
    //        //arc->Geom().Print(std::cout);
    //    }
    //    elementid++;
    //
    //    // Create Geometrical Arc #3
    //    topology[0] = 3;
    //    topology[1] = 0;
    //    {
    //        TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > *arc = new TPZGeoElRefPattern < pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear> > (elementid,topology, fbc1 /** fbc2 */, *geomesh);
    //        //arc->Geom().Print(std::cout);
    //    }
    
    
    //CONCLUINDO A CONSTRUCAO DA MALHA GEOMETRICA
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
    
    std::ofstream outfile("malhaFatiaCascaEsfera.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(geomesh, outfile, true);
    
    return geomesh;
    
    
}


//-------------------------------------------------------------------------


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
    
    
#ifdef CompleteShpere
    
    REAL r = sqrt(x*x+y*y+z*z);
    REAL theta = acos(z/r);
    REAL phi = atan2(y,x);

    REAL costheta = cos(theta);
    REAL sintheta = sin(theta);
//    REAL sin2theta = sin(2.0*theta);
    

    
    REAL sinphi = sin(phi);
    REAL cosphi = cos(phi);
    REAL sin2phi = sin(2.0*phi);
    
    REAL oneminuscosphicosphi = (1.0-cosphi*cosphi);
    
    REAL npowerofsintheta = 1.0;
    
    for (int i = 1; i < norder ; i++)
    {
        npowerofsintheta *= sintheta;
    }
    
    solp[0] = npowerofsintheta*sintheta*oneminuscosphicosphi;
    
    // Gradient computations
    REAL Thetaunitx = cosphi*costheta;
    REAL Thetaunity = costheta*sinphi;
    REAL Thetaunitz = -sintheta;

    REAL Phiunitx = -sinphi;
    REAL Phiunity = cosphi;
    REAL Phiunitz = 0.0;
    
    REAL dfdTheta   = (REAL(norder)/r)*costheta*npowerofsintheta*sinphi*sinphi;
    REAL dfdPhi     = (1.0/r)*npowerofsintheta*sin2phi;

    flux(0,0) = -1.0*(dfdTheta * Thetaunitx + dfdPhi * Phiunitx);
    flux(1,0) = -1.0*(dfdTheta * Thetaunity + dfdPhi * Phiunity);
    flux(2,0) = -1.0*(dfdTheta * Thetaunitz + dfdPhi * Phiunitz);
    

    
    return;

#else
 
    // sobre a casca da esfera -- conferir o raio aqui usado com o da malha geometrica

    REAL r = sqrt(x*x+y*y+z*z);
    REAL theta = atan2(sqrt(x*x+y*y),z);
    REAL phi =  atan2(y,x);
//    REAL cot = 1.0/tan(theta);
    REAL a = M_PI;
    
    solp[0] = (a-theta)*sin(theta)*sin(theta);
    
    flux(0,0) = -(cos(phi)*cos(theta)*(2.0*(a - theta)*cos(theta) - sin(theta))*sin(theta))/r;
    
    flux(1,0) = -(cos(theta)*sin(phi)*(2.0*(a - theta)*cos(theta) - sin(theta))*sin(theta))/r;
    
    flux(2,0) = -(sin(theta)*sin(theta)*(2.0*(-a + theta)*cos(theta) + sin(theta)))/r;

    
//------------

#ifdef TROPICO
    // Anel
    solp[0] = -2.0*log(cos(theta/2.0))-log(2.0);
    flux(0,0)= (cos(phi)*cos(theta)*tan(theta/2.0))/r;
    flux(1,0)= (cos(theta)*sin(phi)*tan(theta/2.0))/r;
    flux(2,0)= (-1.0+cos(theta))/r;
    
    theta = (atan2(sqrt(x*x+y*y),z));// acos(z/r); //
    phi = atan2(y,x);
    a=5.0*M_PI/16.0;
    solp[0] = (a-theta);
    flux(0,0)= ( cos(phi)*cos(theta) )/r;
    flux(1,0)= ( cos(theta)*sin(phi) )/r;
    flux(2,0)= -( sin(theta) )/r;
#endif
    
#endif
}

void LaplaceInSphere::Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &ff){
  
    REAL x,y,z;

    x = pt[0];
    y = pt[1];
    z = pt[2];
    
#ifdef CompleteShpere
    
    REAL r = sqrt(x*x+y*y+z*z);
    REAL theta = acos(z/r);
    REAL phi = atan2(y,x);

    REAL sintheta = sin(theta);
    REAL sinphi = sin(phi);
    REAL cosphi = cos(phi);
    REAL costheta = cos(theta);
    
    REAL npowerofsintheta = 1.0;
    for (int i = 1; i < norder - 1; i++)
    {
        npowerofsintheta *= sintheta;
    }
    
    ff[0] = -(npowerofsintheta/r*r)*(- REAL(norder) * sintheta*sintheta + cosphi*cosphi*(2.0 + REAL(norder)*sintheta*sintheta)
                                        + (-2.0 + REAL(norder)*REAL(norder)*costheta*costheta)*sinphi*sinphi);

    return;
    
#else
    REAL r = sqrt(x*x+y*y+z*z);
    REAL theta = atan2(sqrt(x*x+y*y),z);
//    REAL phi = atan2(y,x);
//    REAL cot = 1.0/tan(theta);
    REAL a = M_PI;
    ff[0] = -((2.0*(a - theta)*(1.0 + 3.0*cos(2.0*theta)) - 5.0*sin(2.0*theta))/(2.0*r*r));

#ifdef TROPICO
    // anel
    ff[0] = (1.0/(r*r))*(1.0/tan(theta));
#endif
    //ff[0] = -1.0;
    
    //ff[0] = 0.0;
    
#endif
}

void LaplaceInSphere::ForcingH1(const TPZVec<REAL> &pt, TPZVec<STATE> &ff, TPZFMatrix<STATE> &flux)
{
    
    DebugStop();
    
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
    
    DebugStop();
    
}

void LaplaceInSphere::ForcingBC1D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){
    
    REAL x = pt[0];
    REAL y = pt[1];
    REAL z = pt[2];
    
    // sobre a casca da esfera -- conferir o raio aqui usado com o da malha geometrica
    
//    REAL r = sqrt(x*x+y*y+z*z);
    REAL theta = atan2(sqrt(x*x+y*y),z);
//    REAL phi = atan2(y,x);
//    REAL cot = 1.0/tan(theta);
    REAL a = M_PI;
    
    solp[0] = (a-theta)*sin(theta)*sin(theta);
    
    // Anel
#ifdef TROPICO
    solp[0] = -2.0*log(cos(theta/2.0))-log(2.0);
    
    solp[0] = -4.0;
    
    theta = (atan2(sqrt(x*x+y*y),z));// acos(z/r); //
    a=5.0*M_PI/16.0;
    solp[0] = (a-theta);
#endif
}

void LaplaceInSphere::ForcingBC2D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){
    
    REAL x = pt[0];
    REAL y = pt[1];
    REAL z = pt[2];
    
    // sobre a casca da esfera -- conferir o raio aqui usado com o da malha geometrica
    
//    REAL r = sqrt(x*x+y*y+z*z);
    REAL theta = M_PI + atan2(sqrt(x*x+y*y),z);
//    REAL phi = atan2(y,x);
//    REAL cot = 1.0/tan(theta);
    REAL a = M_PI;
    
    solp[0] = (a-theta)*sin(theta)*sin(theta);
    
#ifdef TROPICO
    solp[0] = -2.0*log(cos(theta/2.0))-log(2.0);
    
    solp[0] = -4.0;
    
    theta = (atan2(sqrt(x*x+y*y),z));// acos(z/r); //
    a=5.0*M_PI/16.0;
    solp[0] = (a-theta);
#endif
}

void LaplaceInSphere::ForcingBC3D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){
    
    REAL x = pt[0];
    REAL y = pt[1];
    REAL z = pt[2];
    
    // sobre a casca da esfera -- conferir o raio aqui usado com o da malha geometrica
    
//    REAL r = sqrt(x*x+y*y+z*z);
    REAL theta = M_PI + atan2(sqrt(x*x+y*y),z);
//    REAL phi = atan2(y,x);
//    REAL cot = 1.0/tan(theta);
    REAL a = M_PI;
    
    solp[0] = (a-theta)*sin(theta)*sin(theta);
    
    // Anel
#ifdef TROPICO
    solp[0] = -2.0*log(cos(theta/2.0))-log(2.0);
    
    solp[0] = -4.0;
    
    theta = (atan2(sqrt(x*x+y*y),z));// acos(z/r); //
    a=5.0*M_PI/16.0;
    solp[0] = (a-theta);
#endif
}

void LaplaceInSphere::ForcingBC4D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){

    REAL x = pt[0];
    REAL y = pt[1];
    REAL z = pt[2];
    
    // sobre a casca da esfera -- conferir o raio aqui usado com o da malha geometrica
    
//    REAL r = sqrt(x*x+y*y+z*z);
    REAL theta = M_PI + atan2(sqrt(x*x+y*y),z);
//    REAL phi = atan2(y,x);
//    REAL cot = 1.0/tan(theta);
    REAL a = M_PI;
    
    solp[0] = (a-theta)*sin(theta)*sin(theta);
#ifdef TROPICO
    solp[0] = -2.0*log(cos(theta/2.0))-log(2.0);
    
    solp[0] = -4.0;
    
    theta = (atan2(sqrt(x*x+y*y),z));// acos(z/r); //
    a=5.0*M_PI/16.0;
    solp[0] = (a-theta);
#endif
}

void LaplaceInSphere::ForcingBC5D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){

    REAL x = pt[0];
    REAL y = pt[1];
    REAL z = pt[2];
    
    
    REAL r = sqrt(x*x+y*y+z*z);
    REAL theta = acos(z/r);
    REAL phi = atan2(y,x);
    
//    REAL costheta = cos(theta);
    REAL sintheta = sin(theta);
    //    REAL sin2theta = sin(2.0*theta);
    
    
    
//    REAL sinphi = sin(phi);
    REAL cosphi = cos(phi);
//    REAL sin2phi = sin(2.0*phi);
    
    REAL oneminuscosphicosphi = (1.0-cosphi*cosphi);
    
    REAL npowerofsintheta = 1.0;
    
    for (int i = 1; i < norder ; i++)
    {
        npowerofsintheta *= sintheta;
    }
    
    solp[0] = npowerofsintheta*sintheta*oneminuscosphicosphi;
    
//    // sobre a casca da esfera -- conferir o raio aqui usado com o da malha geometrica
//    
////    REAL r = sqrt(x*x+y*y+z*z);
//    REAL theta = atan2(sqrt(x*x+y*y),z);
////    REAL phi = atan2(y,x);
////    REAL cot = 1.0/tan(theta);
//    REAL a = M_PI;
//    
//    solp[0] = (a-theta)*sin(theta)*sin(theta);
//    solp[0] = 0.0;
}


void LaplaceInSphere::ForcingBC0N(const TPZVec<REAL> &pt, TPZVec<STATE> &normflux){
    
    DebugStop();
}

void LaplaceInSphere::ForcingBC1N(const TPZVec<REAL> &pt, TPZVec<STATE> &normflux){
    
    DebugStop();
    
}

void LaplaceInSphere::ForcingBC2N(const TPZVec<REAL> &pt, TPZVec<STATE> &normflux){
    
    //DebugStop();
    normflux[0] = 0.0;
}

void LaplaceInSphere::ForcingBC3N(const TPZVec<REAL> &pt, TPZVec<STATE> &normflux){
    
    //DebugStop();
    normflux[0] = -1.0;
    
}

void LaplaceInSphere::ForcingBC4N(const TPZVec<REAL> &pt, TPZVec<STATE> &normflux){
    
    //DebugStop();
    normflux[0] = 0.0;
}

void LaplaceInSphere::ForcingBC5N(const TPZVec<REAL> &pt, TPZVec<STATE> &normflux){
    
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
    solexata = new TPZDummyFunction<STATE>(SolExata, 5);
    material->SetForcingFunctionExact(solexata);
    
    //funcao do lado direito da equacao do problema
    //    TPZAutoPointer<TPZFunction<STATE> > forcef;
    //    forcef = new TPZDummyFunction<STATE>(Forcing);
    //    material->SetForcingFunction(forcef);
    
    //funcao do lado direito da equacao do problema
    TPZDummyFunction<STATE> *dum = new TPZDummyFunction<STATE>(ForcingH1, 5);
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
    //TPZMatPoisson3d *material = new TPZMatPoisson3d(fmatId,fDim);
    //TPZMatPoissonD3 *material = new TPZMatPoissonD3(fmatId,fDim);
    TPZMatMixedPoisson3D *material = new TPZMatMixedPoisson3D(fmatId,fDim);
    
    TPZMaterial * mat(material);
    material->NStateVariables();
    
    //  TPZAutoPointer<TPZFunction<STATE> > force1 = new TPZDummyFunction<STATE>(Forcing1);
    //	material->SetForcingFunction(force1);
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    
    
    ///Inserir condicao de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    TPZMaterial * BCond0;
    TPZMaterial * BCond1;
//    TPZMaterial * BCond2;
//    TPZMaterial * BCond3;
//    TPZMaterial * BCond4;
//    TPZMaterial * BCond5;
    
//    if( dim == 3 ) { BCond0 = material->CreateBC(mat, fbc0,fdirichlet, val1, val2);}
//    BCond1 = material->CreateBC(mat, fbc1,fdirichlet, val1, val2);
//    BCond2 = material->CreateBC(mat, fbc2,fdirichlet, val1, val2);
//    BCond3 = material->CreateBC(mat, fbc3,fdirichlet, val1, val2);
//    BCond4 = material->CreateBC(mat, fbc4,fdirichlet, val1, val2);
//    if( dim == 3 ) { BCond5 = material->CreateBC(mat, fbc5,fdirichlet, val1, val2);}
    

    
    cmesh->InsertMaterialObject(mat);
    cmesh->SetDimModel(dim);
    cmesh->SetAllCreateFunctionsHDiv();
    
    BCond0 = material->CreateBC(mat, fbc0,fdirichlet, val1, val2);
    BCond1 = material->CreateBC(mat, fbc1,fdirichlet, val1, val2);
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
    
    
//    if( dim == 3 ) { cmesh->InsertMaterialObject(BCond0); }
//    cmesh->InsertMaterialObject(BCond1);
//    cmesh->InsertMaterialObject(BCond2);
//    cmesh->InsertMaterialObject(BCond3);
//    cmesh->InsertMaterialObject(BCond4);
//    if( dim == 3 ) { cmesh->InsertMaterialObject(BCond5); }
    
    cmesh->SetDefaultOrder(pOrder);
    
    
//    bool isgeoblend = false;
//    if (!isgeoblend) {
        TPZLagrangeMultiplier *matskelet = new TPZLagrangeMultiplier(fmatskeleton, dim-1, 1);
        TPZMaterial * mat2(matskelet);
        cmesh->InsertMaterialObject(mat2);
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
    
    
    this->SetupDisconnectedHdivboud(fbc0,fbc1,cmesh);
    
    std::ofstream sout("Fluxcmesh.txt");
    cmesh->Print(sout);
    
    
    
    return cmesh;
    
}

void LaplaceInSphere::SetupDisconnectedHdivboud(const int left,const int rigth, TPZCompMesh * cmesh)
{
    
    int64_t ncels = cmesh->NElements();
    for (int64_t icel = 0; icel < ncels; icel++) {
        TPZCompEl * cel = cmesh->Element(icel);
        if (!cel) {
            DebugStop();
        }
        
        if (cel->Dimension() != 1) {
            continue;
        }
        
        if (cel->Reference()->MaterialId() == fbc0)
        {
            TPZConnect & connect = cel->Connect(0);
            int64_t newindex = cmesh->AllocateNewConnect(connect);
            TPZGeoEl * gel = cel->Reference();
            if (!gel) {
                DebugStop();
            }
            
            TPZGeoElSide gelside(gel,gel->NSides()-1);
            
            TPZStack<TPZCompElSide > celstack;
            gelside.EqualLevelCompElementList(celstack, 1, 0);
            
            int64_t sz = celstack.size();
            if (sz == 0) {
                DebugStop();
            }
            
            for (int64_t iside=0; iside < sz; iside++) {
                TPZCompElSide cels = celstack[iside];
                
                if(cels.Element()->Dimension() == 2)
                {
                    TPZInterpolationSpace * InterpElside = dynamic_cast<TPZInterpolationSpace *>(cels.Element());
                    int nsideconnects = InterpElside->NSideConnects(cels.Side());
                    if(nsideconnects != 1 )
                    {
                        DebugStop();
                    }
                    int localindex = InterpElside->SideConnectLocId(0, cels.Side());
                    InterpElside->SetConnectIndex(localindex, newindex);
//                    int orientside = InterpElside->GetSideOrient(cels.Side());
                    InterpElside->SetSideOrient(cels.Side(),1);
                    
                    cel->SetConnectIndex(0, newindex);
                    TPZInterpolationSpace * celbound  = dynamic_cast<TPZInterpolationSpace *>(cel);
                    celbound->SetSideOrient(2,1);
                    cmesh->ExpandSolution();
//                    cmesh->CleanUpUnconnectedNodes();
                    break;
                }
            }
        }
    }

}

TPZCompMesh *LaplaceInSphere::CMeshPressure(TPZGeoMesh *gmesh, int pOrder, int dim)
{
    /// criar materiais
    //TPZMatPoisson3d *material = new TPZMatPoisson3d(fmatId,fDim);
    //TPZMatPoissonD3 *material = new TPZMatPoissonD3(fmatId,fDim);
    TPZMatMixedPoisson3D *material = new TPZMatMixedPoisson3D(fmatId,fDim);
    
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
    //    TPZMaterial * BCond5 = material->CreateBC(mat,fbc5,10, val1, val2);
    //
    //    cmesh->InsertMaterialObject(BCond0);
    //    cmesh->InsertMaterialObject(BCond1);
    //    cmesh->InsertMaterialObject(BCond2);
    //    cmesh->InsertMaterialObject(BCond3);
    //    cmesh->InsertMaterialObject(BCond4);
    //    cmesh->InsertMaterialObject(BCond5);
    
    cmesh->SetDefaultOrder(pOrder);
    //cmesh->SetDimModel(dim);
    
    bool h1function = true;//(espaco sera sempre assim quando estiver usando elementos wrap)
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
    
//#ifdef PZDEBUG
//    int ncel = cmesh->NElements();
//    for(int i =0; i<ncel; i++){
//        TPZCompEl * compEl = cmesh->ElementVec()[i];
//        if(!compEl) continue;
//        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
//        if(facel)DebugStop();
//        
//    }
//#endif
    
    
    
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
    bool condensacaoestatica = true; //false; //
    
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
    // Hard coded
    for (int id = 0; id < 3; id++){
        PermTensor(id,id) = 1.0;
        InvPermTensor(id,id) = 1.0;
    }
    
    material->SetPermeabilityTensor(PermTensor, InvPermTensor);
    
    //solucao exata
    TPZAutoPointer<TPZFunction<STATE> > solexata;
    
    solexata = new TPZDummyFunction<STATE>(SolExata, 5);
    material->SetForcingFunctionExact(solexata);
    mphysics->SetDimModel(dim);
    //funcao do lado direito da equacao do problema
    TPZDummyFunction<STATE> *dum = new TPZDummyFunction<STATE>(Forcing, 5);
    TPZAutoPointer<TPZFunction<STATE> > forcef;
    dum->SetPolynomialOrder(20);
    forcef = dum;
    material->SetForcingFunction(forcef);
    
    //inserindo o material na malha computacional
    TPZMaterial *mat(material);
    mphysics->InsertMaterialObject(mat);
    
    
    //Criando condicoes de contorno
    TPZMaterial * BCond0;
    TPZMaterial * BCond1;
//    TPZMaterial * BCond2;
//    TPZMaterial * BCond3;
//    TPZMaterial * BCond4;
//    TPZMaterial * BCond5;
    
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZDummyFunction<STATE> * FBCond0D = new TPZDummyFunction<STATE>(ForcingBC5D, 5);
    //FBCond0D->SetPolynomialOrder(20);
    TPZAutoPointer<TPZFunction<STATE> > FBCond0 = FBCond0D;
    BCond0 = material->CreateBC(mat, fbc0,fdirichlet, val1, val2);
    BCond0->SetForcingFunction(FBCond0);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZDummyFunction<STATE> * FBCond1D = new TPZDummyFunction<STATE>(ForcingBC5D, 5);
    //FBCond0D->SetPolynomialOrder(20);
    TPZAutoPointer<TPZFunction<STATE> > FBCond1 = FBCond1D;
    BCond1 = material->CreateBC(mat, fbc1,fdirichlet, val1, val2);
    BCond1->SetForcingFunction(FBCond1);
    
    mphysics->InsertMaterialObject(BCond0);
    mphysics->InsertMaterialObject(BCond1);

    
//    val2(0,0) = 0.0;
//    val2(1,0) = 0.0;
//    TPZAutoPointer<TPZFunction<STATE> > FBCond2 = new TPZDummyFunction<STATE>(ForcingBC2D);
//    BCond2 = material->CreateBC(mat, fbc2,fdirichlet, val1, val2);
//    BCond2->SetForcingFunction(FBCond2);
////    TPZAutoPointer<TPZFunction<STATE> > FBCond2 = new TPZDummyFunction<STATE>(ForcingBC2N);
////    BCond2 = material->CreateBC(mat, fbc2,fneumann, val1, val2);
////    BCond2->SetForcingFunction(FBCond2);
//    
//    val2(0,0) = 0.0;
//    val2(1,0) = 0.0;
//    TPZAutoPointer<TPZFunction<STATE> > FBCond3 = new TPZDummyFunction<STATE>(ForcingBC3D);
//    BCond3 = material->CreateBC(mat, fbc3,fdirichlet, val1, val2);
//    BCond3->SetForcingFunction(FBCond3);
////    TPZAutoPointer<TPZFunction<STATE> > FBCond3 = new TPZDummyFunction<STATE>(ForcingBC3N);
////    BCond3 = material->CreateBC(mat, fbc3,fneumann, val1, val2);
////    BCond3->SetForcingFunction(FBCond3);
//    
//    val2(0,0) = 0.0;
//    val2(1,0) = 0.0;
//    TPZAutoPointer<TPZFunction<STATE> > FBCond4 = new TPZDummyFunction<STATE>(ForcingBC4D);
//    BCond4 = material->CreateBC(mat, fbc4,fdirichlet, val1, val2);
//    BCond4->SetForcingFunction(FBCond4);
////    TPZAutoPointer<TPZFunction<STATE> > FBCond4 = new TPZDummyFunction<STATE>(ForcingBC4N);
////    BCond4 = material->CreateBC(mat, fbc4,fneumann, val1, val2);
////    BCond4->SetForcingFunction(FBCond4);
//    
//    if (dim==3)
//    {
//        val2(0,0) = 0.0;
//        val2(1,0) = 0.0;
//        TPZAutoPointer<TPZFunction<STATE> > FBCond5 = new TPZDummyFunction<STATE>(ForcingBC5D);
//        BCond5 = material->CreateBC(mat, fbc5,10, val1, val2);
//        BCond5->SetForcingFunction(FBCond5);
//    }
    
//    if( dim == 3 ) { mphysics->InsertMaterialObject(BCond0); }
//    mphysics->InsertMaterialObject(BCond0);
//    mphysics->InsertMaterialObject(BCond1);
//    mphysics->InsertMaterialObject(BCond3);
//    mphysics->InsertMaterialObject(BCond4);
//    if( dim == 3 ) { mphysics->InsertMaterialObject(BCond5); }
    
    
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    
    //Fazendo auto build
    mphysics->AutoBuild();
    mphysics->AdjustBoundaryElements();
    mphysics->CleanUpUnconnectedNodes();
    
    if (condensacaoestatica)
    {
        //Creating multiphysic elements containing skeletal elements.
        TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
        mphysics->Reference()->ResetReference();
        mphysics->LoadReferences();
        
        int64_t nel = mphysics->ElementVec().NElements();
        
        std::map<int64_t, int64_t> bctoel, eltowrap;
        for (int64_t el=0; el<nel; el++) {
            TPZCompEl *cel = mphysics->Element(el);
            TPZGeoEl *gel = cel->Reference();
            int matid = gel->MaterialId();
            if (matid < 0) {
                TPZGeoElSide gelside(gel,gel->NSides()-1);
                TPZGeoElSide neighbour = gelside.Neighbour();
                while (neighbour != gelside) {
                    if (neighbour.Element()->Dimension() == dim && neighbour.Element()->Reference()) {
                        // got you!!
                        bctoel[el] = neighbour.Element()->Reference()->Index();
                        break;
                    }
                    neighbour = neighbour.Neighbour();
                }
                if (neighbour == gelside) {
                    DebugStop();
                }
            }
        }
        
        TPZStack< TPZStack< TPZMultiphysicsElement *,7> > wrapEl;
        for(int64_t el = 0; el < nel; el++)
        {
            TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *>(mphysics->Element(el));
            if(mfcel->Dimension()==dim) TPZBuildMultiphysicsMesh::AddWrap(mfcel, fmatId, wrapEl);//criei elementos com o mesmo matId interno, portanto nao preciso criar elemento de contorno ou outro material do tipo TPZLagrangeMultiplier
        }
        
        for (int64_t el =0; el < wrapEl.size(); el++) {
            TPZCompEl *cel = wrapEl[el][0];
            int64_t index = cel->Index();
            eltowrap[index] = el;
        }
        
        meshvec[0]->CleanUpUnconnectedNodes();
        TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
        TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
        
        std::map<int64_t, int64_t>::iterator it;
        for (it = bctoel.begin(); it != bctoel.end(); it++) {
            int64_t bcindex = it->first;
            int64_t elindex = it->second;
            if (eltowrap.find(elindex) == eltowrap.end()) {
                DebugStop();
            }
            int64_t wrapindex = eltowrap[elindex];
            TPZCompEl *bcel = mphysics->Element(bcindex);
            TPZMultiphysicsElement *bcmf = dynamic_cast<TPZMultiphysicsElement *>(bcel);
            if (!bcmf) {
                DebugStop();
            }
            wrapEl[wrapindex].Push(bcmf);
            
        }
        
        //------- Create and add group elements -------
        int64_t index, nenvel;
        nenvel = wrapEl.NElements();
        TPZStack<TPZElementGroup *> elgroups;
        for(int64_t ienv=0; ienv<nenvel; ienv++){
            TPZElementGroup *elgr = new TPZElementGroup(*wrapEl[ienv][0]->Mesh(),index);
            elgroups.Push(elgr);
            nel = wrapEl[ienv].NElements();
            for(int jel=0; jel<nel; jel++){
                elgr->AddElement(wrapEl[ienv][jel]);
            }
        }
        
        mphysics->ComputeNodElCon();
        // create condensed elements
        // increase the NumElConnected of one pressure connects in order to prevent condensation
        for (int64_t ienv=0; ienv<nenvel; ienv++) {
            TPZElementGroup *elgr = elgroups[ienv];
            int nc = elgr->NConnects();
            for (int ic=0; ic<nc; ic++) {
                TPZConnect &c = elgr->Connect(ic);
                if (c.LagrangeMultiplier() > 0) {
                    c.IncrementElConnected();
                    break;
                }
            }
            new TPZCondensedCompEl(elgr);
//            TPZCondensedCompEl *condense = new TPZCondensedCompEl(elgr);
        }
        
        mphysics->CleanUpUnconnectedNodes();
        mphysics->ExpandSolution();
    }
    else
    {
        TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
        mphysics->Reference()->ResetReference();
        mphysics->LoadReferences();
        
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

//        //------- Create and add group elements -------
//        int64_t index, nenvel;
//        nenvel = wrapEl.NElements();
//        for(int ienv=0; ienv<nenvel; ienv++){
//            TPZElementGroup *elgr = new TPZElementGroup(*wrapEl[ienv][0]->Mesh(),index);
//            nel = wrapEl[ienv].NElements();
//            for(int jel=0; jel<nel; jel++){
//                elgr->AddElement(wrapEl[ienv][jel]);
//            }
//        }
        
    }

    
    
    mphysics->ComputeNodElCon();
    mphysics->CleanUpUnconnectedNodes();
    mphysics->ExpandSolution();
    
    return mphysics;
    
}


void LaplaceInSphere::ErrorHDiv(TPZCompMesh *hdivmesh, int p, int ndiv, std::map<REAL, REAL> &fDebugMapL2, std::map<REAL, REAL> &fDebugMapHdiv)
{
    int64_t nel = hdivmesh->NElements();
    int dim = hdivmesh->Dimension();
    TPZManVector<REAL,10> globalerrors(10,0.);
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = hdivmesh->ElementVec()[el];
        if(cel->Reference()->Dimension()!=dim) continue;
        TPZManVector<REAL,10> elerror(10,0.);
        elerror.Fill(0.);
        cel->EvaluateError(SolExata, elerror, false);
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
    int64_t nel = l2mesh->NElements();
    //int dim = l2mesh->Dimension();
    TPZManVector<REAL,10> globalerrors(10,0.);
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = l2mesh->ElementVec()[el];
        TPZManVector<REAL,10> elerror(10,0.);
        cel->EvaluateError(SolExata, elerror, false);
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

void LaplaceInSphere::ErrorPrimalDual(TPZCompMesh *l2mesh, TPZCompMesh *hdivmesh,  int p, int ndiv, std::ostream &out, int DoFT, int DofCond)
{
    int64_t nel = hdivmesh->NElements();
    int dim = hdivmesh->Dimension();
    TPZManVector<REAL,10> globalerrorsDual(10,0.);
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = hdivmesh->ElementVec()[el];
        if(cel->Reference()->Dimension()!=dim) continue;
        TPZManVector<REAL,10> elerror(10,0.);
        elerror.Fill(0.);
        cel->EvaluateError(SolExata, elerror, false);
        int nerr = elerror.size();
        for (int i=0; i<nerr; i++) {
            globalerrorsDual[i] += elerror[i]*elerror[i];
        }
        
        
    }
    
    
    nel = l2mesh->NElements();
    //int dim = l2mesh->Dimension();
    TPZManVector<REAL,10> globalerrorsPrimal(10,0.);
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = l2mesh->ElementVec()[el];
        TPZManVector<REAL,10> elerror(10,0.);
        cel->EvaluateError(SolExata, elerror, false);
        int nerr = elerror.size();
        globalerrorsPrimal.resize(nerr);
        //#ifdef LOG4CXX
        //        if (logdata->isDebugEnabled()) {
        //            std::stringstream sout;
        //            sout << "L2 Error sq of element " << el << elerror[0]*elerror[0];
        //            LOGPZ_DEBUG(logdata, sout.str())
        //        }
        //#endif
        for (int i=0; i<nerr; i++) {
            globalerrorsPrimal[i] += elerror[i]*elerror[i];
        }
        
    }
    
    out << ndiv << setw(10) << DoFT << setw(20) << DofCond << setw(28) << sqrt(globalerrorsPrimal[1]) << setw(35)  << sqrt(globalerrorsDual[1])  << endl;
    
}

void LaplaceInSphere::ChangeExternalOrderConnects(TPZCompMesh *mesh){
    
    int nEl= mesh-> NElements();
    int dim = mesh->Dimension();
    
    int cordermin = -1;
    for (int iel=0; iel<nEl; iel++) {
        TPZCompEl *cel = mesh->ElementVec()[iel];
        if (!cel) continue;
        int ncon = cel->NConnects();
        int corder = 0;
        int nshape = 0;
        
        if(cel->Dimension()== dim){
            for (int icon=0; icon<ncon-1; icon++){
                TPZConnect &co  = cel->Connect(icon);
                corder = co.Order();
                nshape = co.NShape();
                if(corder!=cordermin){
                    cordermin = corder-1;
                    int64_t cindex = cel->ConnectIndex(icon);
                    co.SetOrder(cordermin,cindex);
                    co.SetNShape(nshape-1);
                    mesh->Block().Set(co.SequenceNumber(),nshape-1);
                }
            }
        }
    }
    mesh->ExpandSolution();
    mesh->CleanUpUnconnectedNodes();
}



void LaplaceInSphere::RotateGeomesh(TPZGeoMesh *gmesh, REAL CounterClockwiseAngle, int &Axis)
{
    REAL theta =  (M_PI/180.0)*CounterClockwiseAngle;
    // It represents a 3D rotation around the z axis.
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
    
    TPZVec<REAL> iCoords(3,0.0);
    TPZVec<REAL> iCoordsRotated(3,0.0);
    
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






