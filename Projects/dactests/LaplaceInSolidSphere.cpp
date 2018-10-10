//
//  LaplaceInSolidSphere.cpp
//  PZ
//
//  Created by Omar on 5/9/16.
//
//

#include "LaplaceInSolidSphere.h"
#include "pzcheckgeom.h"
#include "tools.h"


//#define Solution1
//#define Solution2
//#define Solution3
#define Solution4

const int  norder = 6;

LaplaceInSolidSphere::LaplaceInSolidSphere()
{
    fDim = 3;
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
    fIsNonLinearMeshQ = false;
}

LaplaceInSolidSphere::~LaplaceInSolidSphere()
{
    
}

void LaplaceInSolidSphere::Run(int ordemP, int ndiv, std::map<REAL, REAL> &fDebugMapL2, std::map<REAL, REAL> &fDebugMapHdiv, std::ofstream &saidaErro, bool HdivMaisMais)
{
    std::cout<< " (SOLID SPHERE) - grau  polinomio " << ordemP << " numero de divisoes " << ndiv << std::endl;
    std::cout<< " Dimensao == " << fDim << std::endl;
    TPZGeoMesh *gmesh;
    
    if(fIsNonLinearMeshQ){
        gmesh = MakeSphereFromQuadrilateralFaces(ndiv);
    }
    else{
        gmesh = MakeSphereFromLinearQuadrilateralFaces(ndiv);
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
    
    //malha multifisica
    TPZVec<TPZCompMesh *> meshvec(2);
    meshvec[0] = cmesh1;
    meshvec[1] = cmesh2;
    
    TPZCompMesh * mphysics = CMeshMixed(gmesh,meshvec);
    
#ifdef PZDBUG
    
    std::ofstream sout("FluxcmeshFinal.txt");
    meshvec[0]->Print(sout);
    
    std::ofstream sout2("cmeshmphysics.txt");
    mphysics->Print(sout2);
    
#endif
    
    DofCond = mphysics->NEquations();
    
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
    
    tools::PosProcessMultphysics(meshvec,  mphysics, an, plotfile, fDim);
    
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
    REAL error_primal, error_dual;
    ErrorPrimalDual(cmesh2, cmesh1, error_primal, error_dual );
    
    // Printing required information
    
    saidaErro << ndiv << setw(10) << DoFT << setw(20) << DofCond << setw(20) << t1 << setw(20) << t2 << setw(20) << (t1+t2) << setw(20) << error_primal << setw(20) << error_dual << endl;
    
    std::cout<< " FIM (ESFERA) - grau  polinomio " << ordemP << " numero de divisoes " << ndiv << std::endl;
    
    
}



// theta (0,pi) angulo que se inicia no polo norte. phi (0,2pi) o angulo no plano xy
TPZVec<REAL> LaplaceInSolidSphere::SphereToKartesian(REAL r, REAL theta, REAL phi)
{
    TPZVec<REAL> xyz(3,0.0);
    xyz[0] = r*cos(phi)*sin(theta);
    xyz[1] = r*sin(theta)*sin(phi);
    xyz[2] = r*cos(theta);
    return xyz;
}
// theta (0,pi) angulo que se inicia no polo norte. phi (0,2pi) o angulo no plano xy
TPZVec<REAL> LaplaceInSolidSphere::SphereToKartesian(TPZManVector<REAL> xc, REAL r, REAL theta, REAL phi)
{
    TPZVec<REAL> xyz(3,0.0);
    xyz[0] = xc[0] + r*sin(theta)*cos(phi);
    xyz[1] = xc[1] + r*sin(theta)*sin(phi);
    xyz[2] = xc[2] + r*cos(theta);
    return xyz;
}

TPZGeoMesh *LaplaceInSolidSphere::MakeSphereFromLinearQuadrilateralFaces(int ndiv){

    TPZGeoMesh * geomesh = new TPZGeoMesh;
    geomesh->SetDimension(fDim);
    int nl = 2;// Let it fixed
    int basenodes = 8;
    int nodes =  basenodes * (nl);
    REAL radius_o = 1.0;
    REAL radius_i = 0.25;
    
    REAL dr = (radius_o- radius_i)/REAL(nl-1);
    geomesh->SetMaxNodeId(nodes-1);
    geomesh->NodeVec().Resize(nodes);
    TPZManVector<TPZGeoNode,4> Node(nodes);
    
    TPZManVector<int64_t,4> TopolQuad(4);
    TPZManVector<int64_t,8> TopolCube(8);
    TPZManVector<REAL,3> coord(3,0.);
    TPZVec<REAL> xc(3,0.);
    REAL cphi = atan(sqrt(2.0));
    
    int nodeindex = 0;
    int64_t id = 0;
    int matid = 1;
    
    TPZManVector< TPZManVector<REAL,3> , 8 > points(nodes,0.);
    for (int il = 0; il < nl; il++) {
        
        if (il==0) {
            matid = fbc0;
        }
        
        if (il==nl-1) {
            matid = fbc1;
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
    this->RotateGeomesh(geomesh, angle, axis);
    
    ofstream argm("NiceSphere.txt");
    geomesh->Print(argm);
    
    std::ofstream outfile("NiceSphere.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(geomesh, outfile, true);
    
    return geomesh;
}

TPZGeoMesh *LaplaceInSolidSphere::MakeSphereFromQuadrilateralFaces(int ndiv)
{
    TPZGeoMesh * geomesh = new TPZGeoMesh;
    geomesh->SetDimension(fDim);
    int nl = 2;// Let it fixed
    int basenodes = 8;
    int nodes =  basenodes * (nl);
    REAL radius_o = 1.0;
    REAL radius_i = 0.25;

    REAL dr = (radius_o- radius_i)/REAL(nl-1);
    geomesh->SetMaxNodeId(nodes-1);
    geomesh->NodeVec().Resize(nodes);
    TPZManVector<TPZGeoNode,4> Node(nodes);
    
    TPZManVector<int64_t,4> TopolQuad(4);
    TPZManVector<int64_t,8> TopolCube(8);
    TPZManVector<REAL,3> coord(3,0.);
    TPZVec<REAL> xc(3,0.);
    REAL cphi = atan(sqrt(2.0));
    
    int nodeindex = 0;
    int64_t id = 0;
    int matid = 1;
    
    TPZManVector< TPZManVector<REAL,3> , 8 > points(nodes,0.);
    for (int il = 0; il < nl; il++) {
        
        if (il==0) {
            matid = fbc0;
        }
        
        if (il==nl-1) {
            matid = fbc1;
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
    this->RotateGeomesh(geomesh, angle, axis);
    
    ofstream argm("NiceSphere.txt");
    geomesh->Print(argm);
    
    std::ofstream outfile("NiceSphere.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(geomesh, outfile, true);

    return geomesh;
}

TPZManVector<REAL,3> LaplaceInSolidSphere::ParametricSphere(REAL radius,REAL phi,REAL theta)
{
    TPZManVector<REAL,3> xcoor(3,0.0);
    xcoor[0] = radius * cos(theta) * sin(phi);
    xcoor[1] = radius * sin(theta) * sin(phi);
    xcoor[2] = radius * cos(phi) ;
    return xcoor;
}



//-------------------------------------------------------------------------


void LaplaceInSolidSphere::SolExata(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux){
    
    solp.resize(1);
    flux.Resize(4, 1);
    
    solp[0]=0.;
    flux(0,0)=0.0;
    flux(1,0)=0.0;
    flux(2,0)=0.0;
    flux(3,0)=0.0;
    
    REAL x = pt[0];
    REAL y = pt[1];
    REAL z = pt[2];
    
    REAL p = 0.0;
    REAL r = sqrt(x*x+y*y+z*z);
    REAL theta = acos(z/r);
    REAL phi = atan2(y,x);
    
//    REAL costheta = cos(theta);
    REAL sintheta = sin(theta);
    
//    REAL cos2phi = cos(2.0*phi);
//    REAL cos2theta = cos(2.0*theta);
    
//    REAL sinphi = sin(phi);
    REAL cosphi = cos(phi);
//    REAL sin2phi = sin(2.0*phi);
    
    
    REAL oneminuscosphicosphi = (1.0-cosphi*cosphi);
    
    REAL npowerofsintheta = 1.0;
    for (int i = 1; i <= norder ; i++)
    {
        npowerofsintheta *= sintheta;
    }
    
    p = r*r*r*r*npowerofsintheta*oneminuscosphicosphi;
    
    // Gradient computations
    
/*    REAL Radialunitx = sintheta*cosphi;
    REAL Radialunity = sintheta*sinphi;
    REAL Radialunitz = costheta;
    
    REAL Thetaunitx = cosphi*costheta;
    REAL Thetaunity = costheta*sinphi;
    REAL Thetaunitz = -sintheta;
    
    REAL Phiunitx = -sinphi;
    REAL Phiunity = cosphi;
    REAL Phiunitz = 0.0; */
    
#ifdef Solution1
    
    REAL dfdr       = 2.0*r;
    REAL dfdTheta   = 0.0;
    REAL dfdPhi     = 0.0;
    
    solp[0] = r*r;
    
    flux(0,0) = -1.0*(dfdr * Radialunitx + dfdTheta * Thetaunitx + dfdPhi * Phiunitx);
    flux(1,0) = -1.0*(dfdr * Radialunity + dfdTheta * Thetaunity + dfdPhi * Phiunity);
    flux(2,0) = -1.0*(dfdr * Radialunitz + dfdTheta * Thetaunitz + dfdPhi * Phiunitz);
    
    flux(3,0) = -6.0;
    
#endif
    
#ifdef Solution2
    
    REAL dfdr       = 4.0*r*r*r*npowerofsintheta*sinphi*sinphi;
    REAL dfdTheta   = 6.0*r*r*r*costheta*(sintheta*sintheta*sintheta*sintheta*sintheta)*sinphi*sinphi;
    REAL dfdPhi     = r*r*r*(sintheta*sintheta*sintheta*sintheta*sintheta)*sin2phi;
    
    solp[0] = p;
    
    flux(0,0) = -1.0*(dfdr * Radialunitx + dfdTheta * Thetaunitx + dfdPhi * Phiunitx);
    flux(1,0) = -1.0*(dfdr * Radialunity + dfdTheta * Thetaunity + dfdPhi * Phiunity);
    flux(2,0) = -1.0*(dfdr * Radialunitz + dfdTheta * Thetaunitz + dfdPhi * Phiunitz);
    
    flux(3,0) = -2.0*r*r*sintheta*sintheta*sintheta*sintheta*( cos2phi + 0.5*(25.0+11.0*cos2theta)*sinphi*sinphi );

#endif
    
#ifdef Solution3
    
    REAL sinlpix = sin(flambda*M_PI*x);
    REAL sinlpiy = sin(flambda*M_PI*y);
    REAL sinlpiz = sin(flambda*M_PI*z);

    REAL coslpix = cos(flambda*M_PI*x);
    REAL coslpiy = cos(flambda*M_PI*y);
    REAL coslpiz = cos(flambda*M_PI*z);

    
    solp[0] = sinlpix*sinlpiy*sinlpiz;
    
    flux(0,0) = -1.0*(flambda*M_PI*coslpix*sinlpiy*sinlpiz);
    flux(1,0) = -1.0*(flambda*M_PI*sinlpix*coslpiy*sinlpiz);
    flux(2,0) = -1.0*(flambda*M_PI*sinlpix*sinlpiy*coslpiz);
    flux(3,0) = 3.0*M_PI*M_PI*flambda*flambda*sinlpix*sinlpiy*sinlpiz;
    
#endif
    
#ifdef Solution4
    
    solp[0] = (1. - x)*x*(1. - y)*y*(1. - z)*z;
    
    flux(0,0) = -1.0*(-2.*(-0.5 + x)*(-1. + y)*y*(-1. + z)*z);
    flux(1,0) = -1.0*(-2.*(-1. + x)*x*(-0.5 + y)*(-1. + z)*z);
    flux(2,0) = -1.0*(-2.*(-1. + x)*x*(-1. + y)*y*(-0.5 + z));
    flux(3,0) = 2.*(-1. + x)*x*(-1. + y)*y + 2.*(-1. + x)*x*(-1. + z)*z + 2.*(-1. + y)*y*(-1. + z)*z;
    
#endif
    
    return;

}

void LaplaceInSolidSphere::Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &ff){
    
//    REAL flambda = 2.0;
    REAL x,y,z;
    
    x = pt[0];
    y = pt[1];
    z = pt[2];
    
//    REAL r = sqrt(x*x+y*y+z*z);
//    REAL phi = atan2(y,x);
    
#ifdef Solution1
    
    REAL sintheta = sin(theta);
    REAL sinphi = sin(phi);
    REAL cos2phi = cos(2.0*phi);
    REAL cos2theta = cos(2.0*theta);
    
    ff[0] = -6.0;
    
#endif
    
#ifdef Solution2
    
    REAL sintheta = sin(theta);
    REAL sinphi = sin(phi);
    REAL cos2phi = cos(2.0*phi);
    REAL cos2theta = cos(2.0*theta);
    
    ff[0] = -2.0*r*r*sintheta*sintheta*sintheta*sintheta*( cos2phi + 0.5*(25.0+11.0*cos2theta)*sinphi*sinphi );
    
#endif
    
    
    
#ifdef Solution3
    
    REAL sinlpix = sin(flambda*M_PI*x);
    REAL sinlpiy = sin(flambda*M_PI*y);
    REAL sinlpiz = sin(flambda*M_PI*z);
    
    ff[0] = 3.0*M_PI*M_PI*flambda*flambda*sinlpix*sinlpiy*sinlpiz;
    
#endif
    
#ifdef Solution4
    
    ff[0] = 2.*(-1. + x)*x*(-1. + y)*y + 2.*(-1. + x)*x*(-1. + z)*z + 2.*(-1. + y)*y*(-1. + z)*z;
    
#endif
    
    
    return;
}

void LaplaceInSolidSphere::ForcingH1(const TPZVec<REAL> &pt, TPZVec<STATE> &ff, TPZFMatrix<STATE> &flux)
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

void LaplaceInSolidSphere::ForcingBC0D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){
    
    DebugStop();
    
}

void LaplaceInSolidSphere::ForcingBC1D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){
    
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

void LaplaceInSolidSphere::ForcingBC2D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){
    
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

void LaplaceInSolidSphere::ForcingBC3D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){
    
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

void LaplaceInSolidSphere::ForcingBC4D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){
    
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

void LaplaceInSolidSphere::ForcingBC5D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){
    
//    REAL flambda = 2.0;
    solp.resize(1);
    
    solp[0]=0.;

    
    REAL x = pt[0];
    REAL y = pt[1];
    REAL z = pt[2];
    
    REAL p = 0.0;
    REAL r = sqrt(x*x+y*y+z*z);
    REAL theta = acos(z/r);
//    REAL phi = atan2(y,x);
    
//    REAL costheta = cos(theta);
    REAL sintheta = sin(theta);
    
//    REAL cos2phi = cos(2.0*phi);
//    REAL cos2theta = cos(2.0*theta);
    
//    REAL sinphi = sin(phi);
//    REAL cosphi = cos(phi);
//    REAL sin2phi = sin(2.0*phi);
    
//    REAL oneminuscosphicosphi = (1.0-cosphi*cosphi);
    
    REAL npowerofsintheta = 1.0;
    
    for (int i = 1; i <= norder ; i++)
    {
        npowerofsintheta *= sintheta;
    }
    
    
#ifdef Solution1
    
    solp[0] = r*r;
    
#endif
    
#ifdef Solution2
    
    p = r*r*r*r*npowerofsintheta*oneminuscosphicosphi;
    solp[0] = p;
    
#endif

#ifdef Solution3
    
    
    REAL sinlpix = sin(flambda*M_PI*x);
    REAL sinlpiy = sin(flambda*M_PI*y);
    REAL sinlpiz = sin(flambda*M_PI*z);
    
    p = sinlpix*sinlpiy*sinlpiz;
    solp[0] = p;
    
#endif
    
#ifdef Solution4
    
    p = (1. - x)*x*(1. - y)*y*(1. - z)*z;
    solp[0] = p;
    
#endif
    
    
    
    return;
}


void LaplaceInSolidSphere::ForcingBC0N(const TPZVec<REAL> &pt, TPZVec<STATE> &normflux){
    
    DebugStop();
}

void LaplaceInSolidSphere::ForcingBC1N(const TPZVec<REAL> &pt, TPZVec<STATE> &normflux){
    
    DebugStop();
    
}

void LaplaceInSolidSphere::ForcingBC2N(const TPZVec<REAL> &pt, TPZVec<STATE> &normflux){
    
    //DebugStop();
    normflux[0] = 0.0;
}

void LaplaceInSolidSphere::ForcingBC3N(const TPZVec<REAL> &pt, TPZVec<STATE> &normflux){
    
    //DebugStop();
    normflux[0] = -1.0;
    
}

void LaplaceInSolidSphere::ForcingBC4N(const TPZVec<REAL> &pt, TPZVec<STATE> &normflux){
    
    //DebugStop();
    normflux[0] = 0.0;
}

void LaplaceInSolidSphere::ForcingBC5N(const TPZVec<REAL> &pt, TPZVec<STATE> &normflux){
    
    DebugStop();
}

TPZCompMesh *LaplaceInSolidSphere::CMeshH1(TPZGeoMesh *gmesh, int pOrder, int dim)
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


TPZCompMesh *LaplaceInSolidSphere::CMeshFlux(TPZGeoMesh *gmesh, int pOrder, int dim)
{

    TPZMatMixedPoisson3D *material = new TPZMatMixedPoisson3D(fmatId,fDim);
    
    TPZMaterial * mat(material);
    material->NStateVariables();

    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    
    
    ///Inserir condicao de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    TPZMaterial * BCond0;
    TPZMaterial * BCond1;

    cmesh->InsertMaterialObject(mat);
    cmesh->SetDimModel(dim);
    cmesh->SetAllCreateFunctionsHDiv();
    
    BCond0 = material->CreateBC(mat, fbc0,fdirichlet, val1, val2);
    BCond1 = material->CreateBC(mat, fbc1,fdirichlet, val1, val2);
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);

    
    cmesh->SetDefaultOrder(pOrder);
    
    
    TPZLagrangeMultiplier *matskelet = new TPZLagrangeMultiplier(fmatskeleton, dim-1, 1);
    TPZMaterial * mat2(matskelet);
    cmesh->InsertMaterialObject(mat2);

    
    
    
    //Ajuste da estrutura de dados computacional
    cmesh->AutoBuild();
    
//    
//    
//    this->SetupDisconnectedHdivboud(fbc0,fbc1,cmesh);
    
    std::ofstream sout("Fluxcmesh.txt");
    cmesh->Print(sout);
    
    
    
    return cmesh;
    
}

void LaplaceInSolidSphere::SetupDisconnectedHdivboud(const int left,const int rigth, TPZCompMesh * cmesh)
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

TPZCompMesh *LaplaceInSolidSphere::CMeshPressure(TPZGeoMesh *gmesh, int pOrder, int dim)
{

    TPZMatMixedPoisson3D *material = new TPZMatMixedPoisson3D(fmatId,fDim);    
    material->NStateVariables();
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    
    cmesh->SetDefaultOrder(pOrder);
    
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
            
            
//            if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
//            {
//                if(ftriang==true) celdisc->SetTotalOrderShape();
//                else celdisc->SetTensorialShape();
//            }
            
        }
    }
    
    return cmesh;
    
}

TPZCompMesh *LaplaceInSolidSphere::CMeshMixed(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec)
{

    int intorder = 15;
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
    dum->SetPolynomialOrder(intorder);
    forcef = dum;
    material->SetForcingFunction(forcef);
    
    //inserindo o material na malha computacional
    TPZMaterial *mat(material);
    mphysics->InsertMaterialObject(mat);
    
    
    //Criando condicoes de contorno
    TPZMaterial * BCond0;
    TPZMaterial * BCond1;
    
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZDummyFunction<STATE> * FBCond0D = new TPZDummyFunction<STATE>(ForcingBC5D, 5);
    FBCond0D->SetPolynomialOrder(intorder);
    TPZAutoPointer<TPZFunction<STATE> > FBCond0 = FBCond0D;
    BCond0 = material->CreateBC(mat, fbc0,fdirichlet, val1, val2);
    BCond0->SetForcingFunction(FBCond0);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZDummyFunction<STATE> * FBCond1D = new TPZDummyFunction<STATE>(ForcingBC5D, 5);
    FBCond0D->SetPolynomialOrder(intorder);
    TPZAutoPointer<TPZFunction<STATE> > FBCond1 = FBCond1D;
    BCond1 = material->CreateBC(mat, fbc1,fdirichlet, val1, val2);
    BCond1->SetForcingFunction(FBCond1);
    
    mphysics->InsertMaterialObject(BCond0);
    mphysics->InsertMaterialObject(BCond1);
    
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    
    //Fazendo auto build
    mphysics->AutoBuild();
    mphysics->AdjustBoundaryElements();
    mphysics->CleanUpUnconnectedNodes();
    
    if(material->IsUsedSecondIntegration()==false){
        // Creating multiphysic elements into mphysics computational mesh
        TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
        TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
        TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
        
        
        mphysics->ComputeNodElCon();
        // create condensed elements
        // increase the NumElConnected of one pressure connects in order to prevent condensation
        for (int64_t icel=0; icel < mphysics->NElements(); icel++) {
            TPZCompEl  * cel = mphysics->Element(icel);
            if(!cel) continue;
            int nc = cel->NConnects();
            for (int ic=0; ic<nc; ic++) {
                TPZConnect &c = cel->Connect(ic);
                if (c.LagrangeMultiplier() > 0) {
                    c.IncrementElConnected();
                    break;
                }
            }
            new TPZCondensedCompEl(cel);
        }
        
        
        mphysics->CleanUpUnconnectedNodes();
        mphysics->ExpandSolution();
    }
    else{
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
        }
        
        mphysics->CleanUpUnconnectedNodes();
        mphysics->ExpandSolution();
    }
    
    return mphysics;
}



void LaplaceInSolidSphere::ErrorPrimalDual(TPZCompMesh *l2mesh, TPZCompMesh *hdivmesh,  REAL &error_primal , REAL & error_dual)
{
    int64_t nel = hdivmesh->NElements();
    int dim = hdivmesh->Dimension();
    TPZManVector<REAL,10> globalerrorsDual(10,0.   );
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

        for (int i=0; i<nerr; i++) {
            globalerrorsPrimal[i] += elerror[i]*elerror[i];
        }
        
    }
    
    error_primal    = globalerrorsPrimal[1];
    error_dual      = globalerrorsDual[1];
    
}

void LaplaceInSolidSphere::ChangeExternalOrderConnects(TPZCompMesh *mesh){
    
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



void LaplaceInSolidSphere::RotateGeomesh(TPZGeoMesh *gmesh, REAL CounterClockwiseAngle, int &Axis)
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






