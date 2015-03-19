//
//  LaplaceInCircle.cpp
//  PZ
//
//  Created by Douglas Castro on 1/29/15.
//
//

#include "LaplaceInCircle.h"
#include "tools.h"

#define LINEAR

LaplaceInCircle::LaplaceInCircle(int ordemP, int ndiv, std::map<REAL, REAL> &fDebugMapL2, std::map<REAL, REAL> &fDebugMapHdiv)
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
    isH1 = false;
    ftriang = false;
    isgeoblend = true;
    
    
    std::cout<< " INICIO(CIRCULO) - grau  polinomio " << ordemP << " numero de divisoes " << ndiv << std::endl;
    std::cout<< " Dimensao == " << fDim << std::endl;
    
#ifdef LINEAR
    TPZGeoMesh *gmesh = this->GmeshCirculoPorElementosRetos(ndiv);
#else
    TPZGeoMesh *gmesh = this->GMeshCirculoGeob(ndiv);
#endif
    
    
    gmesh->SetDimension(fDim);
    {
//        ofstream argm("gmesh2d-circulo.txt");
//        gmesh->Print(argm);
    }
    
    TPZCompMesh *cmesh2 = this->CMeshPressure(gmesh, ordemP, fDim);
    TPZCompMesh *cmesh1 = this->CMeshFlux(gmesh, ordemP, fDim);
    
    // Um teste para a solucao via H1, sem hdiv
    if (isH1) {
        TPZCompMesh *cmeshH1 = this->CMeshH1(gmesh, ordemP, fDim);
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
        
        tools::PosProcess(anh1, plotfileh1, fDim);
        
        return ;
    }
    // exit
    
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
    {
//        ofstream arg5("cmeshmultiphysics.txt");
//        mphysics->Print(arg5);
    }
    
    TPZAnalysis an(mphysics, true);
    
    tools::SolveSyst(an, mphysics);
    
    stringstream ref,grau;
    grau << ordemP;
    ref << ndiv;
    string strg = grau.str();
    string strr = ref.str();
#ifdef LINEAR
    std::string plotname("OurSolutionMetaCirculoLINEAR");
#else
    std::string plotname("OurSolutionMetaCirculo");
#endif
    
    std::string Grau("P");
    std::string Ref("H");
    std::string VTK(".vtk");
    std::string plotData;
    plotData = plotname+Grau+strg+Ref+strr+VTK;
    std::string plotfile(plotData);
    
    tools::PosProcessMultphysics(meshvec,  mphysics, an, plotfile, fDim);
    
    //Calculo do erro
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
    TPZVec<REAL> erros;
    
    std::cout << "Postprocessed\n";
    
    stringstream ss;
    ss << ordemP;
    string str = ss.str();
    
    std::cout<< " grau  polinomio " << ordemP << " numero de divisoes " << ndiv << std::endl;
#ifdef LINEAR
    std::string filename("InputDataMetaCirculoLINEAR");
#else
    std::string filename("InputDataMetaCirculo");
#endif
    std::string L2("L2.txt");
    std::string Hdiv("Hdiv.txt");
    std::string HdivData,L2Data;
    HdivData = filename+str+Hdiv;
    L2Data = filename+str+L2;
    
    ErrorHDiv(cmesh1, ordemP, ndiv, fDebugMapL2, fDebugMapHdiv);
    
    ErrorL2(cmesh2, ordemP, ndiv, fDebugMapL2, fDebugMapHdiv);
    
    tools::PrintDebugMapForMathematica(HdivData, L2Data, fDebugMapL2, fDebugMapHdiv);
    
    std::cout<< " FIM (CIRCULO) - grau  polinomio " << ordemP << " numero de divisoes " << ndiv << std::endl;
    
}

LaplaceInCircle::~LaplaceInCircle()
{
    
}


// malha que representa um circulo de raio R formado por 4 quadrilateros com geoblends
TPZGeoMesh *LaplaceInCircle::GMeshCirculoGeob( int ndiv)
{
    if(fDim != 2)
    {
        DebugStop();
    }
    
    TPZGeoMesh * gmesh;// = new TPZGeoMesh;
    
    // description of Geometry and application
    // 2D Cylindrical Domain boundaries
    //int matId = 1;
    int arc1 = fbc1; // -1;
    int arc2 = fbc2; // -2;
    int arc3 = fbc3; // -3;
    int arc4 = fbc4; // -4;
    
    int nodenumber = 12;
    REAL ModelRadius = 1.5;
    REAL ModelRadiusInt2 = 0.5;//ModelRadius/2.;
    
    gmesh = new TPZGeoMesh;
    gmesh->NodeVec().Resize(nodenumber);
    
    // Setting node coordantes for Arc3D 1
    int id = 0;
    //0
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,ModelRadius );//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    id++;
    //1
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,ModelRadius);//coord Y
    id++;
    //2
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-ModelRadius );//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    id++;
    //3
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,-ModelRadius);//coord Y
    id++;
    //4
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,sqrt(2.)*ModelRadius/2.);//coord X
    gmesh->NodeVec()[id].SetCoord(1,sqrt(2.)*ModelRadius/2.);//coord Y
    id++;
    //5
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-sqrt(2.)*ModelRadius/2.);//coord X
    gmesh->NodeVec()[id].SetCoord(1,sqrt(2.)*ModelRadius/2.);//coord Y
    id++;
    //6
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-sqrt(2.)*ModelRadius/2.);//coord X
    gmesh->NodeVec()[id].SetCoord(1,-sqrt(2.)*ModelRadius/2.);//coord Y
    id++;
    //7
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,sqrt(2.)*ModelRadius/2.);//coord X
    gmesh->NodeVec()[id].SetCoord(1,-sqrt(2.)*ModelRadius/2.);//coord Y
    id++;
    //8
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,ModelRadiusInt2 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    id++;
    //9
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,ModelRadiusInt2);//coord Y
    id++;
    //10
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-ModelRadiusInt2 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    id++;
    //11
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,-ModelRadiusInt2);//coord Y
    
    
    int elementid = 0;
    TPZVec < long > nodeindex(3,0.0);
    nodeindex.resize(4);
    
    // Create Geometrical Quad #1
    nodeindex[0] = 0;
    nodeindex[1] = 1;
    nodeindex[2] = 9;
    nodeindex[3] = 8;
    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoQuad > > (elementid,nodeindex, fmatId,*gmesh);
    elementid++;
    
    // Create Geometrical Quad #2
    nodeindex[0] = 1;
    nodeindex[1] = 2;
    nodeindex[2] = 10;
    nodeindex[3] = 9;
    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoQuad > > (elementid,nodeindex, fmatId,*gmesh);
    elementid++;
    
    // Create Geometrical Quad #3
    nodeindex[0] = 2;
    nodeindex[1] = 3;
    nodeindex[2] = 11;
    nodeindex[3] = 10;
    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoQuad > > (elementid,nodeindex, fmatId,*gmesh);
    elementid++;
    
    // Create Geometrical Quad #4
    nodeindex[0] = 3;
    nodeindex[1] = 0;
    nodeindex[2] = 8;
    nodeindex[3] = 11;
    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoQuad > > (elementid,nodeindex, fmatId,*gmesh);
    elementid++;
    
    // Create Geometrical Quad #5
    nodeindex[0] = 8;
    nodeindex[1] = 9;
    nodeindex[2] = 10;
    nodeindex[3] = 11;
    new TPZGeoElRefPattern < pzgeom::TPZGeoQuad >  (elementid,nodeindex, fmatId,*gmesh);
    elementid++;
    
    // Definition of Arc coordenates
    nodeindex.resize(3);
    // Create Geometrical Arc #1
    nodeindex[0] = 0;
    nodeindex[1] = 1;
    nodeindex[2] = 4;
    new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, arc3, *gmesh);
    elementid++;
    
    // Create Geometrical Arc #2
    nodeindex[0] = 1;
    nodeindex[1] = 2;
    nodeindex[2] = 5;
    new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, arc4, *gmesh);
    elementid++;
    
    // Create Geometrical Arc #3
    nodeindex[0] = 2;
    nodeindex[1] = 3;
    nodeindex[2] = 6;
    new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, arc1, *gmesh);
    elementid++;
    
    // Create Geometrical Arc #4
    nodeindex[0] = 3;
    nodeindex[1] = 0;
    nodeindex[2] = 7;
    new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, arc2, *gmesh);
    
    
    
    gmesh->BuildConnectivity();
    
    int nref = ndiv;
    TPZVec<TPZGeoEl *> sons;
    for (int iref = 0; iref < nref; iref++) {
        int nel = gmesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = gmesh->ElementVec()[iel];
            if (gel->HasSubElement()) {
                continue;
            }
            gel->Divide(sons);
        }
    }
    
    
    std::ofstream out("CurvoDAC3.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
    
    return gmesh;
}

TPZGeoMesh *LaplaceInCircle::GMeshCirculoTriangGeob(int ndiv)
{
    if(fDim != 2)
    {
        DebugStop();
    }
    
    TPZGeoMesh * gmesh;// = new TPZGeoMesh;
    
    // description of Geometry and application
    // 2D Cylindrical Domain boundaries
    //int matId = 1;
    int arc1 = fbc1; // -1;
    int arc2 = fbc2; // -2;
    int arc3 = fbc3; // -3;
    int arc4 = fbc4; // -4;
    
    int nodenumber = 13;
    REAL ModelRadius = 1.5;
    REAL ModelRadiusInt2 = 0.5;//ModelRadius/2.;
    
    gmesh = new TPZGeoMesh;
    gmesh->NodeVec().Resize(nodenumber);
    
    // Setting node coordantes for Arc3D 1
    int id = 0;
    //0
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,ModelRadius );//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    id++;
    //1
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,ModelRadius);//coord Y
    id++;
    //2
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-ModelRadius );//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    id++;
    //3
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,-ModelRadius);//coord Y
    id++;
    //4
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,sqrt(2.)*ModelRadius/2.);//coord X
    gmesh->NodeVec()[id].SetCoord(1,sqrt(2.)*ModelRadius/2.);//coord Y
    id++;
    //5
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-sqrt(2.)*ModelRadius/2.);//coord X
    gmesh->NodeVec()[id].SetCoord(1,sqrt(2.)*ModelRadius/2.);//coord Y
    id++;
    //6
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-sqrt(2.)*ModelRadius/2.);//coord X
    gmesh->NodeVec()[id].SetCoord(1,-sqrt(2.)*ModelRadius/2.);//coord Y
    id++;
    //7
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,sqrt(2.)*ModelRadius/2.);//coord X
    gmesh->NodeVec()[id].SetCoord(1,-sqrt(2.)*ModelRadius/2.);//coord Y
    id++;
    //8
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,ModelRadiusInt2 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    id++;
    //9
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,ModelRadiusInt2);//coord Y
    id++;
    //10
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-ModelRadiusInt2 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    id++;
    //11
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,-ModelRadiusInt2);//coord Y
    id++;
    //12
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    
    int elementid = 0;
    TPZVec < long > nodeindex(3,0.0);
    nodeindex.resize(3);
    
    // Create Geometrical Quad #1
    nodeindex[0] = 0;
    nodeindex[1] = 1;
    nodeindex[2] = 9;
    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoTriangle > > (elementid,nodeindex, fmatId,*gmesh);
    elementid++;
    
    // Create Geometrical Quad #2
    nodeindex[0] = 0;
    nodeindex[1] = 9;
    nodeindex[2] = 8;
    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoTriangle > > (elementid,nodeindex, fmatId,*gmesh);
    elementid++;
    
    // Create Geometrical Quad #3
    nodeindex[0] = 1;
    nodeindex[1] = 2;
    nodeindex[2] = 10;
    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoTriangle > > (elementid,nodeindex, fmatId,*gmesh);
    elementid++;
    
    // Create Geometrical Quad #4
    nodeindex[0] = 1;
    nodeindex[1] = 10;
    nodeindex[2] = 9;
    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoTriangle > > (elementid,nodeindex, fmatId,*gmesh);
    elementid++;
    
    // Create Geometrical Quad #4
    nodeindex[0] = 2;
    nodeindex[1] = 3;
    nodeindex[2] = 11;
    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoTriangle > > (elementid,nodeindex, fmatId,*gmesh);
    elementid++;
    // Create Geometrical Quad #4
    nodeindex[0] = 2;
    nodeindex[1] = 11;
    nodeindex[2] = 10;
    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoTriangle > > (elementid,nodeindex, fmatId,*gmesh);
    elementid++;
    // Create Geometrical Quad #4
    nodeindex[0] = 3;
    nodeindex[1] = 0;
    nodeindex[2] = 8;
    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoTriangle > > (elementid,nodeindex, fmatId,*gmesh);
    elementid++;
    // Create Geometrical Quad #4
    nodeindex[0] = 3;
    nodeindex[1] = 8;
    nodeindex[2] = 11;
    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoTriangle > > (elementid,nodeindex, fmatId,*gmesh);
    elementid++;
    
    // Create Geometrical Quad #5
    nodeindex.resize(3);
    nodeindex[0] = 8;
    nodeindex[1] = 9;
    nodeindex[2] = 12;
    new TPZGeoElRefPattern < pzgeom::TPZGeoTriangle >  (elementid,nodeindex, fmatId,*gmesh);
    elementid++;
    // Create Geometrical Quad #5
    nodeindex[0] = 9;
    nodeindex[1] = 10;
    nodeindex[2] = 12;
    new TPZGeoElRefPattern < pzgeom::TPZGeoTriangle >  (elementid,nodeindex, fmatId,*gmesh);
    elementid++;
    // Create Geometrical Quad #5
    nodeindex[0] = 10;
    nodeindex[1] = 11;
    nodeindex[2] = 12;
    new TPZGeoElRefPattern < pzgeom::TPZGeoTriangle >  (elementid,nodeindex, fmatId,*gmesh);
    elementid++;
    // Create Geometrical Quad #5
    nodeindex[0] = 11;
    nodeindex[1] = 8;
    nodeindex[2] = 12;
    new TPZGeoElRefPattern < pzgeom::TPZGeoTriangle >  (elementid,nodeindex, fmatId,*gmesh);
    elementid++;
    
    // Definition of Arc coordenates
    nodeindex.resize(3);
    // Create Geometrical Arc #1
    nodeindex[0] = 0;
    nodeindex[1] = 1;
    nodeindex[2] = 4;
    new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, arc3, *gmesh);
    elementid++;
    
    // Create Geometrical Arc #2
    nodeindex[0] = 1;
    nodeindex[1] = 2;
    nodeindex[2] = 5;
    new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, arc4, *gmesh);
    elementid++;
    
    // Create Geometrical Arc #3
    nodeindex[0] = 2;
    nodeindex[1] = 3;
    nodeindex[2] = 6;
    new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, arc1, *gmesh);
    elementid++;
    
    // Create Geometrical Arc #4
    nodeindex[0] = 3;
    nodeindex[1] = 0;
    nodeindex[2] = 7;
    new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, arc2, *gmesh);
    
    
    
    gmesh->BuildConnectivity();
    
    int nref = ndiv;
    TPZVec<TPZGeoEl *> sons;
    for (int iref = 0; iref < nref; iref++) {
        int nel = gmesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = gmesh->ElementVec()[iel];
            if (gel->HasSubElement()) {
                continue;
            }
            gel->Divide(sons);
        }
    }
    
    
    std::ofstream out("CurvoCirculo.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
    
    return gmesh;
}

// malha que representa um circulo de raio R formado por 5 quadrilateros com quadraticos
TPZGeoMesh *LaplaceInCircle::GMeshCirculoQuad( int ndiv)
{
    if(fDim != 2)
    {
        DebugStop();
    }
    
    TPZGeoMesh * gmesh;// = new TPZGeoMesh;
    
    // description of Geometry and application
    // 2D Cylindrical Domain boundaries
    //int matId = 1;
    int arc1 = fbc1; // -1;
    int arc2 = fbc2; // -2;
    int arc3 = fbc3; // -3;
    int arc4 = fbc4; // -4;
    
    int nodenumber = 20;
    REAL ModelRadius = 1.5;
    REAL ModelRadiusInt2 = 0.5;//ModelRadius/2.;
    
    gmesh = new TPZGeoMesh;
    gmesh->NodeVec().Resize(nodenumber);
    
    // Setting node coordantes for Arc3D 1
    int id = 0;
    //0
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,ModelRadius );//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    id++;
    //1
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,ModelRadius);//coord Y
    id++;
    //2
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-ModelRadius );//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    id++;
    //3
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,-ModelRadius);//coord Y
    id++;
    //4
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,sqrt(2.)*ModelRadius/2.);//coord X
    gmesh->NodeVec()[id].SetCoord(1,sqrt(2.)*ModelRadius/2.);//coord Y
    id++;
    //5
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-sqrt(2.)*ModelRadius/2.);//coord X
    gmesh->NodeVec()[id].SetCoord(1,sqrt(2.)*ModelRadius/2.);//coord Y
    id++;
    //6
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-sqrt(2.)*ModelRadius/2.);//coord X
    gmesh->NodeVec()[id].SetCoord(1,-sqrt(2.)*ModelRadius/2.);//coord Y
    id++;
    //7
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,sqrt(2.)*ModelRadius/2.);//coord X
    gmesh->NodeVec()[id].SetCoord(1,-sqrt(2.)*ModelRadius/2.);//coord Y
    id++;
    //8
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,ModelRadiusInt2 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    id++;
    //9
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,ModelRadiusInt2);//coord Y
    id++;
    //10
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-ModelRadiusInt2 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    id++;
    //11
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
    gmesh->NodeVec()[id].SetCoord(1,-ModelRadiusInt2);//coord Y
    id++;
    //12
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,sqrt(2.)*ModelRadiusInt2/2.0);//coord X
    gmesh->NodeVec()[id].SetCoord(1,sqrt(2.)*ModelRadiusInt2/2.0);//coord Y
    id++;
    //13
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-sqrt(2.)*ModelRadiusInt2/2.0);//coord X
    gmesh->NodeVec()[id].SetCoord(1,sqrt(2.)*ModelRadiusInt2/2.0);//coord Y
    id++;
    //14
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-sqrt(2.)*ModelRadiusInt2/2.0);//coord X
    gmesh->NodeVec()[id].SetCoord(1,-sqrt(2.)*ModelRadiusInt2/2.0);//coord Y
    id++;
    //15
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,sqrt(2.)*ModelRadiusInt2/2.0);//coord X
    gmesh->NodeVec()[id].SetCoord(1,-sqrt(2.)*ModelRadiusInt2/2.0);//coord Y
    id++;
    //16
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,(ModelRadiusInt2+ModelRadius)/2.0);//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    id++;
    //17
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0);//coord X
    gmesh->NodeVec()[id].SetCoord(1,(ModelRadiusInt2+ModelRadius)/2.0);//coord Y
    id++;
    //18
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,-(ModelRadiusInt2+ModelRadius)/2.0);//coord X
    gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
    id++;
    //19
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,0.0);//coord X
    gmesh->NodeVec()[id].SetCoord(1,-(ModelRadiusInt2+ModelRadius)/2.0);//coord Y
    
    
    
    int elementid = 0;
    TPZVec < long > nodeindex(3,0.0);
    nodeindex.resize(8);
    
    // Create Quadratic Quad #1
    nodeindex[0] = 0;
    nodeindex[1] = 1;
    nodeindex[2] = 9;
    nodeindex[3] = 8;
    nodeindex[4] = 4;
    nodeindex[5] = 17;
    nodeindex[6] = 12;
    nodeindex[7] = 16;
    new TPZGeoElRefPattern < pzgeom::TPZQuadraticQuad>  (elementid,nodeindex, fmatId,*gmesh);
    //new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoQuad > > (elementid,nodeindex, matId,*gmesh);
    elementid++;
    
    // Create Quadratic Quad #2
    nodeindex[0] = 1;
    nodeindex[1] = 2;
    nodeindex[2] = 10;
    nodeindex[3] = 9;
    nodeindex[4] = 5;
    nodeindex[5] = 18;
    nodeindex[6] = 13;
    nodeindex[7] = 17;
    new TPZGeoElRefPattern < pzgeom::TPZQuadraticQuad>  (elementid,nodeindex, fmatId,*gmesh);
    //new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoQuad > > (elementid,nodeindex, matId,*gmesh);
    elementid++;
    
    // Create Quadratic Quad #3
    nodeindex[0] = 2;
    nodeindex[1] = 3;
    nodeindex[2] = 11;
    nodeindex[3] = 10;
    nodeindex[4] = 6;
    nodeindex[5] = 19;
    nodeindex[6] = 14;
    nodeindex[7] = 18;
    new TPZGeoElRefPattern < pzgeom::TPZQuadraticQuad>  (elementid,nodeindex, fmatId,*gmesh);
    //new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoQuad > > (elementid,nodeindex, matId,*gmesh);
    elementid++;
    
    // Create Quadratic Quad #4
    nodeindex[0] = 3;
    nodeindex[1] = 0;
    nodeindex[2] = 8;
    nodeindex[3] = 11;
    nodeindex[4] = 7;
    nodeindex[5] = 16;
    nodeindex[6] = 15;
    nodeindex[7] = 19;
    new TPZGeoElRefPattern < pzgeom::TPZQuadraticQuad>  (elementid,nodeindex, fmatId,*gmesh);
    //new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoQuad > > (elementid,nodeindex, matId,*gmesh);
    elementid++;
    
    // Create Quadratic Quad #5
    nodeindex[0] = 8;
    nodeindex[1] = 9;
    nodeindex[2] = 10;
    nodeindex[3] = 11;
    nodeindex[4] = 12;
    nodeindex[5] = 13;
    nodeindex[6] = 14;
    nodeindex[7] = 15;
    new TPZGeoElRefPattern < pzgeom::TPZQuadraticQuad>  (elementid,nodeindex, fmatId,*gmesh);
    //new TPZGeoElRefPattern < pzgeom::TPZGeoQuad >  (elementid,nodeindex, matId,*gmesh);
    elementid++;
    
    // Definition of Arc coordenates
    nodeindex.resize(3);
    
    // Create Quadratic Arc #1
    nodeindex[0] = 0;
    nodeindex[1] = 1;
    nodeindex[2] = 4;
    new TPZGeoElRefPattern < pzgeom::TPZQuadraticLine > (elementid,nodeindex, arc3, *gmesh);
    //new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, arc1, *gmesh);
    elementid++;
    
    // Create Quadratic Arc #2
    nodeindex[0] = 1;
    nodeindex[1] = 2;
    nodeindex[2] = 5;
    new TPZGeoElRefPattern < pzgeom::TPZQuadraticLine > (elementid,nodeindex, arc4, *gmesh);
    //new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, arc2, *gmesh);
    elementid++;
    
    // Create Quadratic Arc #3
    nodeindex[0] = 2;
    nodeindex[1] = 3;
    nodeindex[2] = 6;
    new TPZGeoElRefPattern < pzgeom::TPZQuadraticLine > (elementid,nodeindex, arc1, *gmesh);
    //new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, arc3, *gmesh);
    elementid++;
    
    // Create Quadratic Arc #4
    nodeindex[0] = 3;
    nodeindex[1] = 0;
    nodeindex[2] = 7;
    new TPZGeoElRefPattern < pzgeom::TPZQuadraticLine > (elementid,nodeindex, arc2, *gmesh);
    //new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, arc4, *gmesh);
    
    
    
    gmesh->BuildConnectivity();
    
    int nref = ndiv;
    TPZVec<TPZGeoEl *> sons;
    for (int iref = 0; iref < nref; iref++) {
        int nel = gmesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = gmesh->ElementVec()[iel];
            if (gel->HasSubElement()) {
                continue;
            }
            gel->Divide(sons);
        }
    }
    
    
    std::ofstream out("CurvoDAC4.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
    
    return gmesh;
}

TPZGeoMesh *LaplaceInCircle::GmeshCirculoPorElementosRetos( int ndiv)
{
    if(fDim != 2)
    {
        DebugStop();
    }
    
    TPZGeoMesh * gmesh;// = new TPZGeoMesh;
    
    // description of Geometry and application
    //int matId = 1;
    int arc1 = fbc1; // -1;
    int arc2 = fbc2; // -2;
    int arc3 = fbc3; // -3;
    int arc4 = fbc4; // -4;
    
    int numberoflevels = (ndiv+1); // exceto centro
    int nodesperlevel = 4*(ndiv+1);
    int nodenumber = nodesperlevel*(ndiv+1) + 1;
     
    REAL R = 1.0;
    REAL stepR = R/((REAL)numberoflevels);
    TPZManVector<REAL,3> xc(3,0.0);
    
    
    REAL theta = (M_PI/2.0)/((REAL)numberoflevels);
    
    // Para rotacao
    int Axis = 3;
    REAL angulo = 0.0;
    
    gmesh = new TPZGeoMesh;
    gmesh->NodeVec().Resize(nodenumber);
    
    TPZManVector<REAL,3> coord(3,0.);
    
    TPZManVector<int,220> nodeindex(nodenumber,0);
    TPZManVector<int,220> firstnodeindexlevel(numberoflevels+1,0);
    
    int contador = 0;
    
    TPZGeoNode node;
    
    // Setting node coordantes
    long id = 0;
    
    for (int nivel = 0 ; nivel < numberoflevels; nivel++)
    {
        REAL r = R  - ((REAL)nivel)*stepR;
        
        firstnodeindexlevel[nivel] = id;
        
        for (int no = 0; no < nodesperlevel; no++)
        {
            
            coord = PolarToKartesian(r, no*theta, xc);
            tools::RotateNode(coord, angulo, Axis);
            node.SetNodeId(id);
            node.SetCoord(coord);
            gmesh->NodeVec()[id] = node;
            
            nodeindex[id] = id;
            
            id++;
        }
        contador++;
    }
    
    // Centro do circulo
    coord = PolarToKartesian(0.0, 0.0, xc);
    tools::RotateNode(coord, angulo, Axis);
    node.SetNodeId(id);
    node.SetCoord(coord);
    gmesh->NodeVec()[id] = node;
    nodeindex[id] = id;
    firstnodeindexlevel[contador] = id;
    
    
    int elementid = 0;
    TPZVec < long > topology(4,0);
    contador = 0;
    // Building 2d elements
    for (int nivel = 0 ; nivel < numberoflevels-1; nivel++)
    {
        contador++;
        for (int no = 0; no < nodesperlevel; no++)
        {

            // Create Geometrical Quad #1
            long nodestenivel = firstnodeindexlevel[nivel];
            long nodoproximonivel = firstnodeindexlevel[contador];
            //int nodoproximonivelp1 = firstnodeindexlevel[contador+2];
            
            long correcao1 = no < (nodesperlevel - 1) ? (no+nodestenivel+1)%nodoproximonivel : nodestenivel;
            long correcao2 = no < (nodesperlevel - 1) ? (no+nodoproximonivel+1) : nodoproximonivel;

            long a = nodeindex[  no+nodestenivel ];
            long b = nodeindex[correcao1];//(no+nodestenivel+1)%nodoproximonivel
            long c = nodeindex[correcao2];
            long d = nodeindex[(no+nodoproximonivel)];//%nodoproximonivelp1
            
            //std::cout << a << " " << b << " " << c << " " << d << std::endl;
            
            topology[0] = a;
            topology[1] = b;
            topology[2] = c;
            topology[3] = d;
            new TPZGeoElRefPattern<  pzgeom::TPZGeoQuad  > (elementid, topology, fmatId,*gmesh);
            elementid++;

        }
    }
    
    topology.resize(3);
    
    for (int no = 0; no < nodesperlevel; no++)
    {
        
        // Create Geometrical Quad #1
        long nodesteonivel = firstnodeindexlevel[contador];
        long nodoproximonivel = firstnodeindexlevel[contador+1];
        
        long correcao = no < (nodesperlevel - 1) ? (no+nodesteonivel+1) : nodesteonivel;
        
        long a = nodeindex[no+nodesteonivel];
        long b = nodeindex[ correcao ];
        long c = nodeindex[nodoproximonivel];
        
        //std::cout << a << " " << b << " " << c << std::endl;
        
        topology[0] = a;
        topology[1] = b;
        topology[2] = c;
        new TPZGeoElRefPattern<  pzgeom::TPZGeoTriangle  > (elementid, topology, fmatId,*gmesh);
        elementid++;
        
    }
    
    
    //Creating elements 1d
    topology.resize(2);
    
    for (int no = 0; no < nodesperlevel; no++)
    {
        
        // Create Geometrical Quad #1
        long nodesteonivel = firstnodeindexlevel[1];
        
        long correcao = (no+1)%nodesteonivel;
        
        long a = nodeindex[no];
        long b = nodeindex[ correcao ];

        //std::cout << a << " " << b << std::endl;
        
        topology[0] = a;
        topology[1] = b;
        new TPZGeoElRefPattern<  pzgeom::TPZGeoLinear  > (elementid, topology, arc1, *gmesh);
        elementid++;
        
    }
    
    gmesh->BuildConnectivity();
    
    {
//               ofstream argm("gmesh2d-circulo.txt");
//                gmesh->Print(argm);
    }
    
    std::ofstream out("DiscoPorElLineares.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
    
    return gmesh;
}

TPZVec<REAL> LaplaceInCircle::PolarToKartesian(REAL r, REAL theta, TPZManVector<REAL> xc)
{
    TPZVec<REAL> xyz(3,0.0);
    xyz[0] = xc[0] + r*cos(theta);
    xyz[1] = xc[1] + r*sin(theta);
    xyz[2] = xc[2];
    return xyz;
}

void LaplaceInCircle::SolExata(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux){
    
    solp.resize(1);
    solp[0]=0.;
    
    int dim = 3; //getDimension();
    
    // tensor de permutacao
    TPZFNMatrix<2,REAL> TP(dim,dim,0.0);
    TPZFNMatrix<2,REAL> InvTP(dim,dim,0.0);
    
    
    // Hard coded
    for (int id = 0; id < dim; id++){
        TP(id,id) = 1.0;
        InvTP(id,id) = 1.0;
    }
    
    flux.Resize(dim, 1);
    

    
    double x = pt[0];
    double y = pt[1];
    
    REAL r = sqrt( x*x + y*y );
    REAL theta = atan2(y,x);
    
    solp[0] = r*r*(1.0 - r*r);
    flux(0,0)= (2.0*r - 4.0*r*r*r)*cos(theta);
    flux(1,0)=  (2.0*r - 4.0*r*r*r)*sin(theta);
    

    
//    flux(0,0)=flux(1,0)=0.;
//    double x = pt[0];
//    double y = pt[1];
//    REAL raio = sqrt( x*x + y*y );
//    if (raio < 1.0)
//    {
//        
//        
//        // para a Lap p = f no circulo
//        solp[0] = 3.0*exp(1.0/(x*x + y*y-1.0));
//        flux(0,0) = 6.0*exp(1.0/(x*x + y*y-1.0))*(x*TP(0,0)+y*TP(0,1))/( (x*x + y*y-1.0)*(x*x + y*y-1.0) );
//        flux(1,0) = 6.0*exp(1.0/(x*x + y*y-1.0))*(x*TP(1,0)+y*TP(1,1))/( (x*x + y*y-1.0)*(x*x + y*y-1.0) );
//        
//    }
//    else
//    {
//        // para a Lap p = f no circulo
//        solp[0] = 0.0;
//        flux(0,0)= 0.0;
//        flux(1,0)= 0.0;
//        
//    }


    
    


    
}

void LaplaceInCircle::Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &ff){

    int dim = 3; //getDimension();
    
    // tensor de permutacao
    TPZFNMatrix<2,REAL> TP(dim,dim,0.0);
    TPZFNMatrix<2,REAL> InvTP(dim,dim,0.0);
    
    
    // Hard coded
    for (int id = 0; id < dim; id++){
        TP(id,id) = 1.0;
        InvTP(id,id) = 1.0;
    }
    
    double x = pt[0];
    double y = pt[1];
    REAL r = sqrt( x*x + y*y );
    ff[0] = -(4.0 - 16.0*r*r);
    
   
//    double x = pt[0];
//    double y = pt[1];
//    REAL raio = sqrt( x*x + y*y );
//    if (raio < 1.0)
//    {
//        ff[0] = -6.0*exp(1.0/(x*x + y*y-1.0))*((-1.0 + 3.0*x*x*x*x + 2.0*(1.0 + x*x)*y*y - y*y*y*y)*TP(0,0) +2.0*x*y*(-1.0 + 2.0*x*x + 2.0*y*y)*(TP(0,1) + TP(1,0)) + (-1.0 - x*x*x*x +3.0*y*y*y*y + 2.0*x*x*(1.0 + y*y))*TP(1,1));
//        
//    }
//    else
//    {
//        ff[0] = 0.0;
//        
//    }
    
    
}

void LaplaceInCircle::ForcingH1(const TPZVec<REAL> &pt, TPZVec<STATE> &ff, TPZFMatrix<STATE> &flux)
{
    int dim = 2; //getDimension();
    
    // tensor de permutacao
    TPZFNMatrix<2,REAL> TP(dim,dim,0.0);
    TPZFNMatrix<2,REAL> InvTP(dim,dim,0.0);
    
    
    // Hard coded
    for (int id = 0; id < dim; id++){
        TP(id,id) = 1.0;
        InvTP(id,id) = 1.0;
    }
    
    flux.Resize(dim, 1);

    double x = pt[0];
    double y = pt[1];
    REAL raio = sqrt( x*x + y*y );
    if (raio < 1.0)
    {
        ff[0] = 6.0*exp(1.0/(x*x + y*y-1.0))*((-1.0 + 3.0*x*x*x*x + 2.0*(1.0 + x*x)*y*y - y*y*y*y)*TP(0,0) +2.0*x*y*(-1.0 + 2.0*x*x + 2.0*y*y)*(TP(0,1) + TP(1,0)) + (-1.0 - x*x*x*x +3.0*y*y*y*y + 2.0*x*x*(1.0 + y*y))*TP(1,1));
        flux(0,0) = 6.0*exp(1.0/(x*x + y*y-1.0))*(x*TP(0,0)+y*TP(0,1))/( (x*x + y*y-1.0)*(x*x + y*y-1.0) );
        flux(1,0) = 6.0*exp(1.0/(x*x + y*y-1.0))*(x*TP(1,0)+y*TP(1,1))/( (x*x + y*y-1.0)*(x*x + y*y-1.0) );
        
    }
    else
    {
        ff[0] = 0.0;
        flux(0,0)= 0.0;
        flux(1,0)= 0.0;
    }
    

    
}

void LaplaceInCircle::ForcingBC0D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){


    double x = pt[0];
    double y = pt[1];
    
    REAL r = sqrt( x*x + y*y );
    
    //solp[0] = r*r*(1.0 - r*r);
    solp[0] = 0.0;


    //    flux(0,0)=flux(1,0)=0.;
    //    double x = pt[0];
    //    double y = pt[1];
    //    REAL raio = sqrt( x*x + y*y );
    //    if (raio < 1.0)
    //    {
    //
    //
    //        // para a Lap p = f no circulo
    //        solp[0] = 3.0*exp(1.0/(x*x + y*y-1.0));
    //
    //    }
    //    else
    //    {
    //        // para a Lap p = f no circulo
    //        solp[0] = 0.0;
    //        
    //    }
    

    
}

void LaplaceInCircle::ForcingBC1D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){
    
    double x = pt[0];
    double y = pt[1];
    
    REAL r = sqrt( x*x + y*y );
    
    //solp[0] = r*r*(1.0 - r*r);
    solp[0] = 0.0;
    
    
    //    flux(0,0)=flux(1,0)=0.;
    //    double x = pt[0];
    //    double y = pt[1];
    //    REAL raio = sqrt( x*x + y*y );
    //    if (raio < 1.0)
    //    {
    //
    //
    //        // para a Lap p = f no circulo
    //        solp[0] = 3.0*exp(1.0/(x*x + y*y-1.0));
    //
    //    }
    //    else
    //    {
    //        // para a Lap p = f no circulo
    //        solp[0] = 0.0;
    //
    //    }
    
}

void LaplaceInCircle::ForcingBC2D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){

    
    
    double x = pt[0];
    double y = pt[1];
    
    REAL r = sqrt( x*x + y*y );
    
    //solp[0] = r*r*(1.0 - r*r);
    solp[0] = 0.0;
    
    
    //    flux(0,0)=flux(1,0)=0.;
    //    double x = pt[0];
    //    double y = pt[1];
    //    REAL raio = sqrt( x*x + y*y );
    //    if (raio < 1.0)
    //    {
    //
    //
    //        // para a Lap p = f no circulo
    //        solp[0] = 3.0*exp(1.0/(x*x + y*y-1.0));
    //
    //    }
    //    else
    //    {
    //        // para a Lap p = f no circulo
    //        solp[0] = 0.0;
    //
    //    }

}

void LaplaceInCircle::ForcingBC3D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){
    
    
    double x = pt[0];
    double y = pt[1];
    
    REAL r = sqrt( x*x + y*y );
    
    //solp[0] = r*r*(1.0 - r*r);
    solp[0] = 0.0;
    
    
    //    flux(0,0)=flux(1,0)=0.;
    //    double x = pt[0];
    //    double y = pt[1];
    //    REAL raio = sqrt( x*x + y*y );
    //    if (raio < 1.0)
    //    {
    //
    //
    //        // para a Lap p = f no circulo
    //        solp[0] = 3.0*exp(1.0/(x*x + y*y-1.0));
    //
    //    }
    //    else
    //    {
    //        // para a Lap p = f no circulo
    //        solp[0] = 0.0;
    //
    //    }

}

void LaplaceInCircle::ForcingBC4D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){

    
    double x = pt[0];
    double y = pt[1];
    
    REAL r = sqrt( x*x + y*y );
    
    //solp[0] = r*r*(1.0 - r*r);
    solp[0] = 0.0;
    
    
    //    flux(0,0)=flux(1,0)=0.;
    //    double x = pt[0];
    //    double y = pt[1];
    //    REAL raio = sqrt( x*x + y*y );
    //    if (raio < 1.0)
    //    {
    //
    //
    //        // para a Lap p = f no circulo
    //        solp[0] = 3.0*exp(1.0/(x*x + y*y-1.0));
    //
    //    }
    //    else
    //    {
    //        // para a Lap p = f no circulo
    //        solp[0] = 0.0;
    //
    //    }
}

void LaplaceInCircle::ForcingBC5D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){

    double x = pt[0];
    double y = pt[1];
    
    REAL r = sqrt( x*x + y*y );
    
    //solp[0] = r*r*(1.0 - r*r);
    solp[0] = 0.0;
    
    //    flux(0,0)=flux(1,0)=0.;
    //    double x = pt[0];
    //    double y = pt[1];
    //    REAL raio = sqrt( x*x + y*y );
    //    if (raio < 1.0)
    //    {
    //
    //
    //        // para a Lap p = f no circulo
    //        solp[0] = 3.0*exp(1.0/(x*x + y*y-1.0));
    //
    //    }
    //    else
    //    {
    //        // para a Lap p = f no circulo
    //        solp[0] = 0.0;
    //
    //    }
}


void LaplaceInCircle::ForcingBC0N(const TPZVec<REAL> &pt, TPZVec<STATE> &normflux){
    
    DebugStop();
}

void LaplaceInCircle::ForcingBC1N(const TPZVec<REAL> &pt, TPZVec<STATE> &normflux){
    
    DebugStop();
    
}

void LaplaceInCircle::ForcingBC2N(const TPZVec<REAL> &pt, TPZVec<STATE> &normflux){
    DebugStop();
}

void LaplaceInCircle::ForcingBC3N(const TPZVec<REAL> &pt, TPZVec<STATE> &normflux){
    
    DebugStop();
    
}

void LaplaceInCircle::ForcingBC4N(const TPZVec<REAL> &pt, TPZVec<STATE> &normflux){
    
    DebugStop();
}

void LaplaceInCircle::ForcingBC5N(const TPZVec<REAL> &pt, TPZVec<STATE> &normflux){
    
    DebugStop();
}

TPZCompMesh *LaplaceInCircle::LaplaceInCircle::CMeshH1(TPZGeoMesh *gmesh, int pOrder, int dim)
{
    /// criar materiais
    TPZMatPoisson3d *material = new TPZMatPoisson3d(fmatId,fDim);
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
    cmesh->SetDimModel(fDim);
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
    
    if( fDim == 3 ) { BCond0 = material->CreateBC(mat, fbc0,fdirichlet, val1, val2);}
    BCond1 = material->CreateBC(mat, fbc1,fdirichlet, val1, val2);
    BCond2 = material->CreateBC(mat, fbc2,fdirichlet, val1, val2);
    BCond3 = material->CreateBC(mat, fbc3,fdirichlet, val1, val2);
    BCond4 = material->CreateBC(mat, fbc4,fdirichlet, val1, val2);
    if( fDim == 3 ) { BCond5 = material->CreateBC(mat, fbc5,fdirichlet, val1, val2);}
    
    cmesh->InsertMaterialObject(mat);
    
    cmesh->SetAllCreateFunctionsContinuous();
    if( fDim == 3 ) { cmesh->InsertMaterialObject(BCond0); }
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    if( fDim == 3 ) { cmesh->InsertMaterialObject(BCond5); }
    
    //Ajuste da estrutura de dados computacional
    cmesh->AutoBuild();
    
    return cmesh;
    
}


TPZCompMesh *LaplaceInCircle::CMeshFlux(TPZGeoMesh *gmesh, int pOrder, int dim)
{
    /// criar materiais
    //TPZMatPoisson3d *material = new TPZMatPoisson3d(matId,fDim);
    TPZMatPoisson3d *material = new TPZMatPoisson3d(fmatId,fDim);
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
    
    if( fDim == 3 ) { BCond0 = material->CreateBC(mat, fbc0,fdirichlet, val1, val2);}
    BCond1 = material->CreateBC(mat, fbc1,fdirichlet, val1, val2);
    BCond2 = material->CreateBC(mat, fbc2,fdirichlet, val1, val2);
    BCond3 = material->CreateBC(mat, fbc3,fdirichlet, val1, val2);
    BCond4 = material->CreateBC(mat, fbc4,fdirichlet, val1, val2);
    if( fDim == 3 ) { BCond5 = material->CreateBC(mat, fbc5,fdirichlet, val1, val2);}
    
    cmesh->InsertMaterialObject(mat);
    
    cmesh->SetDimModel(fDim);
    
    cmesh->SetAllCreateFunctionsHDiv();
    
    
    if( fDim == 3 ) { cmesh->InsertMaterialObject(BCond0); }
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    if( fDim == 3 ) { cmesh->InsertMaterialObject(BCond5); }
    
    cmesh->SetDefaultOrder(pOrder);
    
    
    if (!isgeoblend) {
        TPZLagrangeMultiplier *matskelet = new TPZLagrangeMultiplier(fmatskeleton, fDim-1, 1);
        TPZMaterial * mat2(matskelet);
        cmesh->InsertMaterialObject(mat2);
    }
    
    
    
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

TPZCompMesh *LaplaceInCircle::CMeshPressure(TPZGeoMesh *gmesh, int pOrder, int dim)
{
    /// criar materiais
    TPZMatPoisson3d *material = new TPZMatPoisson3d(fmatId,fDim);
    material->NStateVariables();
    
    //  TPZAutoPointer<TPZFunction<STATE> > force1 = new TPZDummyFunction<STATE>(Forcing1);
    //	material->SetForcingFunction(force1);
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(fDim);
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
    //cmesh->SetDimModel(fDim);
    
    bool h1function = true;// com esqueleto precisa disso
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

TPZCompMesh *LaplaceInCircle::CMeshMixed(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec)
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
    TPZAutoPointer<TPZFunction<STATE> > FBCond1 = new TPZDummyFunction<STATE>(ForcingBC1D);
    BCond1 = material->CreateBC(mat, fbc1,fdirichlet, val1, val2);
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
    TPZAutoPointer<TPZFunction<STATE> > FBCond3 = new TPZDummyFunction<STATE>(ForcingBC3D);
    BCond3 = material->CreateBC(mat, fbc3,fdirichlet, val1, val2);
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


void LaplaceInCircle::ErrorHDiv(TPZCompMesh *hdivmesh, int p, int ndiv, std::map<REAL, REAL> &fDebugMapL2, std::map<REAL, REAL> &fDebugMapHdiv)
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

void LaplaceInCircle::ErrorL2(TPZCompMesh *l2mesh, int p, int ndiv, std::map<REAL, REAL> &fDebugMapL2, std::map<REAL, REAL> &fDebugMapHdiv)
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






