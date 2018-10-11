//
//  LaplaceInCylinder.cpp
//  PZ
//
//  Created by Douglas Castro on 1/29/15.
//
//

#include "LaplaceInCylinder.h"
#include "tools.h"

LaplaceInCylinder::LaplaceInCylinder()
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


LaplaceInCylinder::~LaplaceInCylinder()
{
    
}

void LaplaceInCylinder::Run(int ordemP, int ndiv, std::map<REAL, REAL> &fDebugMapL2, std::map<REAL, REAL> &fDebugMapHdiv, std::ofstream &saidaErro, bool HdivMaisMais)
{
    std::cout<< " INICIO(CASCA CILINDRO) - grau  polinomio " << ordemP << " numero de divisoes " << ndiv << std::endl;
    std::cout<< " Dimensao == " << fDim << std::endl;
    
    TPZGeoMesh *gmesh = GMeshCilindricalMesh(ndiv);
    //TPZGeoMesh *gmesh = GMesh(fDim, ftriang, ndiv);
    
    
    gmesh->SetDimension(fDim);
    {
        //        ofstream argm("gmesh2d-circulo.txt");
        //        gmesh->Print(argm);
    }
    
    
    // Um teste para a solucao via H1, sem hdiv
    if (fisH1) {
        
        TPZCompMesh *cmeshH1 = this->CMeshH1(gmesh, ordemP, fDim);
        //        {
        //         ofstream arg1("cmeshH1.txt");
        //         cmeshH1->Print(arg1);
        //        }
        
        int dofTotal = cmeshH1->NEquations();
        
        //condensar
        for (int64_t iel=0; iel<cmeshH1->NElements(); iel++) {
            TPZCompEl *cel = cmeshH1->Element(iel);
            if(!cel) continue;
            new TPZCondensedCompEl(cel);
//            TPZCondensedCompEl *condense = new TPZCondensedCompEl(cel);
        }
        
        cmeshH1->ExpandSolution();
        cmeshH1->CleanUpUnconnectedNodes();
        
        int dofCondensed = cmeshH1->NEquations();
        
        
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
        
        tools::PosProcess(anh1, plotfileh1, fDim);
        
        ErrorH1(cmeshH1, ordemP, ndiv, saidaErro, dofTotal, dofCondensed);
        
        return ;
    }
    else
    {
        TPZCompMesh *cmesh1 = this->CMeshFlux(gmesh, ordemP, fDim);
        
        TPZCompMesh *cmesh2 = this->CMeshPressure(gmesh, ordemP, fDim);
        
        if (HdivMaisMais) {
            // para rodar P**
            ChangeExternalOrderConnects(cmesh1);
        }
        int DofCond, DoFT;
        DoFT = cmesh1->NEquations() + cmesh2->NEquations();
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
        DofCond = mphysics->NEquations();
        
        TPZAnalysis an(mphysics, true);
        REAL t1,t2;
        tools::SolveSyst(an, mphysics, t1, t2);
        
        stringstream ref,grau;
        grau << ordemP;
        ref << ndiv;
        string strg = grau.str();
        string strr = ref.str();
        std::string plotname("OurSolutionMetaCilindro");
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
        std::string filename("InputDataMetaCilindro");
        std::string L2("L2.txt");
        std::string Hdiv("Hdiv.txt");
        std::string HdivData,L2Data;
        HdivData = filename+str+Hdiv;
        L2Data = filename+str+L2;
        
        
        ErrorHDiv(cmesh1, ordemP, ndiv, fDebugMapL2, fDebugMapHdiv );
        
        ErrorL2(cmesh2, ordemP, ndiv, fDebugMapL2, fDebugMapHdiv);
        
        ErrorPrimalDual( cmesh2, cmesh1,  ordemP, ndiv, saidaErro, DoFT, DofCond);
        
        
        tools::PrintDebugMapForMathematica(HdivData, L2Data, fDebugMapL2, fDebugMapHdiv);
        
        std::cout<< " FIM (CILIND) - grau  polinomio " << ordemP << " numero de divisoes " << ndiv << std::endl;
    }
    // exit
    
    
}


TPZGeoMesh *LaplaceInCylinder::GMeshCilindricalMesh( int ndiv)
{
    ///INSTANCIACAO DA MALHA GEOMETRICA
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    
    /// Materiais
    
    int nodenumber = 10;
    REAL r = 1.0;
    // o angulo theta varia de theta0 a theta1
    REAL theta0 = 0.0;
    REAL theta1 = M_PI;
    // z varia de z0 a z1
    REAL z0 = -M_PI/2.0;
    REAL z1 = M_PI/2.0;
    
    TPZVec<REAL> coord(3,0.0);
    int Axis = 3;
    REAL angrot = 0.0;
    
    gmesh = new TPZGeoMesh;
    gmesh->NodeVec().Resize(nodenumber);
    
    // Setting node coordinates
    int id = 0;
    //0
    coord[0] = r*cos(theta0);
    coord[1] = r*sin(theta0);
    coord[2] = z0;
    tools::RotateNode(coord, angrot, Axis);
    
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,coord[0]);//coord X
    gmesh->NodeVec()[id].SetCoord(1,coord[1]);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,coord[2]);//coord Z
    id++;
    
    //1
    coord[0] = r*cos(theta1);
    coord[1] = r*sin(theta1);
    coord[2] = z0;
    tools::RotateNode(coord, angrot, Axis);
    
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,coord[0]);//coord X
    gmesh->NodeVec()[id].SetCoord(1,coord[1]);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,coord[2]);//coord Z
    id++;
    
    //2
    coord[0] = r*cos(theta1);
    coord[1] = r*sin(theta1);
    coord[2] = z1;
    tools::RotateNode(coord, angrot, Axis);
    
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,coord[0]);//coord X
    gmesh->NodeVec()[id].SetCoord(1,coord[1]);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,coord[2]);//coord Z
    id++;
    
    //3
    coord[0] = r*cos(theta0);
    coord[1] = r*sin(theta0);
    coord[2] = z1;
    tools::RotateNode(coord, angrot, Axis);
    
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,coord[0]);//coord X
    gmesh->NodeVec()[id].SetCoord(1,coord[1]);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,coord[2]);//coord Z
    id++;
    
    //4
    coord[0] = r*cos(M_PI/4.0);
    coord[1] = r*sin(M_PI/4.0);
    coord[2] = z0;
    tools::RotateNode(coord, angrot, Axis);
    
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,coord[0]);//coord X
    gmesh->NodeVec()[id].SetCoord(1,coord[1]);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,coord[2]);//coord Z
    id++;
    
    //5
    coord[0] = r*cos(M_PI/2.0);
    coord[1] = r*sin(M_PI/2.0);
    coord[2] = z0;
    tools::RotateNode(coord, angrot, Axis);
    
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,coord[0]);//coord X
    gmesh->NodeVec()[id].SetCoord(1,coord[1]);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,coord[2]);//coord Z
    id++;

    //6
    coord[0] = r*cos(3.0*M_PI/4.0);
    coord[1] = r*sin(3.0*M_PI/4.0);
    coord[2] = z0;
    tools::RotateNode(coord, angrot, Axis);
    
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,coord[0]);//coord X
    gmesh->NodeVec()[id].SetCoord(1,coord[1]);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,coord[2]);//coord Z
    id++;
    
    //7
    coord[0] = r*cos(M_PI/4.0);
    coord[1] = r*sin(M_PI/4.0);
    coord[2] = z1;
    tools::RotateNode(coord, angrot, Axis);
    
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,coord[0]);//coord X
    gmesh->NodeVec()[id].SetCoord(1,coord[1]);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,coord[2]);//coord Z
    id++;

    //8
    coord[0] = r*cos(M_PI/2.0);
    coord[1] = r*sin(M_PI/2.0);
    coord[2] = z1;
    tools::RotateNode(coord, angrot, Axis);
    
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,coord[0]);//coord X
    gmesh->NodeVec()[id].SetCoord(1,coord[1]);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,coord[2]);//coord Z
    id++;
    
    //9
    coord[0] = r*cos(3.0*M_PI/4.0);
    coord[1] = r*sin(3.0*M_PI/4.0);
    coord[2] = z1;
    tools::RotateNode(coord, angrot, Axis);
    
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,coord[0]);//coord X
    gmesh->NodeVec()[id].SetCoord(1,coord[1]);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,coord[2]);//coord Z
    id++;
    
    int elementid = 0;
    TPZVec < int64_t > nodeindex(3,0.0);
    
    // Definition of Arc coordenates
    
    
    nodeindex.resize(2);
    
    // Create Geometrical Arc #2
    nodeindex[0] = 0;
    nodeindex[1] = 3;
    new TPZGeoElRefPattern < pzgeom::TPZGeoLinear > (elementid,nodeindex, fbc1, *gmesh);
    elementid++;
    
    nodeindex.resize(3);
    // Create Geometrical Arc #1
    nodeindex[0] = 0;
    nodeindex[1] = 5;
    nodeindex[2] = 4;
    new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, fbc2, *gmesh);
    elementid++;
    // Create Geometrical Arc #1
    nodeindex[0] = 5;
    nodeindex[1] = 1;
    nodeindex[2] = 6;
    new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, fbc2, *gmesh);
    elementid++;

    
    nodeindex.resize(2);
    // Create Geometrical Arc #4
    nodeindex[0] = 1;
    nodeindex[1] = 2;
    new TPZGeoElRefPattern < pzgeom::TPZGeoLinear > (elementid,nodeindex, fbc3, *gmesh);
    
    
    nodeindex.resize(3);
    
    // Create Geometrical Arc #3
    nodeindex[0] = 3;
    nodeindex[1] = 8;
    nodeindex[2] = 7;
    new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, fbc4, *gmesh);
    elementid++;
    // Create Geometrical Arc #1
    nodeindex[0] = 8;
    nodeindex[1] = 2;
    nodeindex[2] = 9;
    new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, fbc2, *gmesh);
    elementid++;
    
    
    nodeindex.resize(4);
    
    // Create Geometrical Quad #1
    nodeindex[0] = 0;
    nodeindex[1] = 5;
    nodeindex[2] = 8;
    nodeindex[3] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoQuad > > (elementid,nodeindex, fmatId,*gmesh);
    elementid++;
    // Create Geometrical Quad #1
    nodeindex[0] = 1;
    nodeindex[1] = 2;
    nodeindex[2] = 8;
    nodeindex[3] = 5;
    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoQuad > > (elementid,nodeindex, fmatId,*gmesh);
    elementid++;
    
    
    
    
    //CONCLUINDO A CONSTRUCAO DA MALHA GEOMETRICA
    gmesh->BuildConnectivity();
    
    TPZVec<TPZGeoEl *> sons;
    const int nref = ndiv;
    for (int iref = 0; iref < nref; iref++) {
        int nel = gmesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = gmesh->ElementVec()[iel];
            if(!gel->HasSubElement()) gel->Divide(sons);
        }
    }
    
    std::ofstream outfile("malhaCilindrica.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outfile, true);
    
    return gmesh;
    
}

TPZGeoMesh *LaplaceInCylinder::GMesh(int dim, bool ftriang, int ndiv)
{
    
    if(dim!=2)
    {
        std::cout << "dimensao errada" << std::endl;
        dim = 2;
        DebugStop();
    }
    
    int Qnodes =  4;
    
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    gmesh->SetMaxNodeId(Qnodes-1);
    gmesh->NodeVec().Resize(Qnodes);
    TPZVec<TPZGeoNode> Node(Qnodes);
    
    gmesh->SetDimension(dim);
    
    TPZVec <int64_t> TopolQuad(4);
    TPZVec <int64_t> TopolTriang(3);
    TPZVec <int64_t> TopolLine(2);
    TPZVec <int64_t> TopolPoint(1);
    
    //indice dos nos
    int64_t id = 0;
    //    REAL valx;
    //    for(int xi = 0; xi < Qnodes/2; xi++)
    //    {
    //        valx = xi*Lx;
    //        Node[id].SetNodeId(id);
    //        Node[id].SetCoord(0 ,valx );//coord X
    //        Node[id].SetCoord(1 ,0. );//coord Y
    //        gmesh->NodeVec()[id] = Node[id];
    //        id++;
    //    }
    //
    //    for(int xi = 0; xi < Qnodes/2; xi++)
    //    {
    //        valx = Lx - xi*Lx;
    //        Node[id].SetNodeId(id);
    //        Node[id].SetCoord(0 ,valx );//coord X
    //        Node[id].SetCoord(1 ,Ly);//coord Y
    //        gmesh->NodeVec()[id] = Node[id];
    //        id++;
    //    }
    //
    TPZManVector<REAL,3> coord(2,0.);
    int in = 0;
    //c0
    coord[0] = 0.0;
    coord[1] = 0.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    //c1
    coord[0] =  1.0;
    coord[1] = 0.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    //c2
    coord[0] =  0.0;
    coord[1] =  1.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    //c3
    coord[0] = 1.0;
    coord[1] =  1.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    //indice dos elementos
    id = 0;
    
    if(ftriang==true) // triangulo
    {
        TopolTriang[0] = 0;
        TopolTriang[1] = 1;
        TopolTriang[2] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (id,TopolTriang,fmatId,*gmesh);
        id++;
        
        TopolTriang[0] = 0;
        TopolTriang[1] = 3;
        TopolTriang[2] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,fmatId,*gmesh);
        id++;
        
        TopolLine[0] = 0;
        TopolLine[1] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,fbc1,*gmesh);
        id++;
        
        TopolLine[0] = 1;
        TopolLine[1] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,fbc2,*gmesh);
        id++;
        
        TopolLine[0] = 3;
        TopolLine[1] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,fbc3,*gmesh);
        id++;
        
        TopolLine[0] = 2;
        TopolLine[1] = 0;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,fbc4,*gmesh);
        id++;
    }
    else{
        
        TopolQuad[0] = 0;
        TopolQuad[1] = 1;
        TopolQuad[2] = 3;
        TopolQuad[3] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,fmatId,*gmesh);
        id++;
        
        TopolLine[0] = 0;
        TopolLine[1] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,fbc1,*gmesh);
        id++;
        
        TopolLine[0] = 1;
        TopolLine[1] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,fbc2,*gmesh);
        id++;
        
        TopolLine[0] = 3;
        TopolLine[1] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,fbc3,*gmesh);
        id++;
        
        TopolLine[0] = 2;
        TopolLine[1] = 0;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,fbc4,*gmesh);
    }
    
    gmesh->BuildConnectivity();
    
    /// gmesh para aqui
    
    TPZVec<TPZGeoEl *> sons;
    for (int iref = 0; iref < ndiv; iref++) {
        int nel = gmesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = gmesh->ElementVec()[iel];
            if (gel->HasSubElement()) {
                continue;
            }
            gel->Divide(sons);
        }
    }
    
    
    //#ifdef LOG4CXX
    //	if(logdata->isDebugEnabled())
    //	{
    //        std::stringstream sout;
    //        sout<<"\n\n Malha Geometrica Inicial\n ";
    //        gmesh->Print(sout);
    //        LOGPZ_DEBUG(logdata,sout.str())
    //	}
    //#endif
    
    std::ofstream outfile("malhaQuadm.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outfile, true);
    
    return gmesh;
}

void LaplaceInCylinder::SolExata(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux){
    
    solp.resize(1);
    solp[0]=0.;
    flux.Resize(3,1);
    
    REAL x = pt[0];
    REAL y = pt[1];
    REAL z = pt[2];
    
    REAL r = sqrt(x*x+y*y);
    REAL theta = atan2(y,x);
    
    solp[0] = sin(theta)*sin(z + M_PI/2.0);
    flux(0,0) = -((cos(theta)*cos(z)*sin(theta))/r);
    flux(1,0) = (cos(theta)*cos(theta)*cos(z))/r;
    flux(2,0) = -sin(theta)*sin(z);
    flux(0,0) *= -1.;
    flux(1,0) *= -1.;
    flux(2,0) *= -1.;
    
}

void LaplaceInCylinder::SolExataH1(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux){
    
    solp.resize(1);
    solp[0]=0.;
    
    int dim = 3;
    
    flux.Resize(dim, 1);
    
    double x = pt[0];
    double y = pt[1];
    double z = pt[2];
    
    
    REAL r = sqrt(x*x+y*y);
    REAL theta = atan2(y,x);
    
    solp[0] = sin(theta)*sin(z + M_PI/2.0);
    flux(0,0) = -((cos(theta)*cos(z)*sin(theta))/r);
    flux(1,0) = (cos(theta)*cos(theta)*cos(z))/r;
    flux(2,0) = -sin(theta)*sin(z);
    
    
    flux(0,0) *= -1.;
    flux(1,0) *= -1.;
    flux(2,0) *= -1.;
    
    
    
}

void  LaplaceInCylinder::Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &ff){
    
    REAL x = pt[0];
    REAL y = pt[1];
    REAL z = pt[2];
    
    REAL r = sqrt(x*x+y*y);
    REAL theta = atan2(y,x);
    
	ff[0] = ((1.0 + r*r)*cos(z)*sin(theta))/(r*r);
}

void  LaplaceInCylinder::ForcingH1(const TPZVec<REAL> &pt, TPZVec<STATE> &ff, TPZFMatrix<STATE> &flux)
{
    flux.Resize(3, 1);

    ff[0] = 0.0;
    flux(0,0) = 0.0;
    flux(1,0) = 0.0;
    flux(2,0) = 0.0;
    // verifique os dados acima
    DebugStop();

}

void  LaplaceInCylinder::ForcingBC0D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){

    DebugStop();
    
}

void  LaplaceInCylinder::ForcingBC1D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){

    REAL x = pt[0];
    REAL y = pt[1];
    REAL z = pt[2];
    
    //REAL r = sqrt(x*x+y*y);
    REAL theta = atan2(y,x);
    
    solp[0] = sin(theta)*sin(z + M_PI/2.0);
}

void  LaplaceInCylinder::ForcingBC2D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){
    
    REAL x = pt[0];
    REAL y = pt[1];
    REAL z = pt[2];
    
    //REAL r = sqrt(x*x+y*y);
    REAL theta = atan2(y,x);
    
    solp[0] = sin(theta)*sin(z + M_PI/2.0);
    
}

void  LaplaceInCylinder::ForcingBC3D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){

    REAL x = pt[0];
    REAL y = pt[1];
    REAL z = pt[2];
    
    //REAL r = sqrt(x*x+y*y);
    REAL theta = atan2(y,x);
    
    solp[0] = sin(theta)*sin(z + M_PI/2.0);
}

void  LaplaceInCylinder::ForcingBC4D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){

    REAL x = pt[0];
    REAL y = pt[1];
    REAL z = pt[2];
    
    //REAL r = sqrt(x*x+y*y);
    REAL theta = atan2(y,x);
    
    solp[0] = sin(theta)*sin(z + M_PI/2.0);
    
}

void  LaplaceInCylinder::ForcingBC5D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){
    DebugStop();
}


void  LaplaceInCylinder::ForcingBC0N(const TPZVec<REAL> &pt, TPZVec<STATE> &normflux){
    
    DebugStop();
}

void  LaplaceInCylinder::ForcingBC1N(const TPZVec<REAL> &pt, TPZVec<STATE> &normflux){
    
    DebugStop();
}

void  LaplaceInCylinder::ForcingBC2N(const TPZVec<REAL> &pt, TPZVec<STATE> &normflux){
    
    DebugStop();
}

void  LaplaceInCylinder::ForcingBC3N(const TPZVec<REAL> &pt, TPZVec<STATE> &normflux){
    
    DebugStop();
}

void  LaplaceInCylinder::ForcingBC4N(const TPZVec<REAL> &pt, TPZVec<STATE> &normflux){
    
    DebugStop();
}


void  LaplaceInCylinder::ForcingBC5N(const TPZVec<REAL> &pt, TPZVec<STATE> &normflux){
    
    DebugStop();
}

TPZCompMesh * LaplaceInCylinder::CMeshH1(TPZGeoMesh *gmesh, int pOrder, int dim)
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


TPZCompMesh * LaplaceInCylinder::CMeshFlux(TPZGeoMesh *gmesh, int pOrder, int dim)
{
    /// criar materiais
    //TPZMatPoisson3d *material = new TPZMatPoisson3d(matId,dim);
    //TPZMatPoisson3d *material = new TPZMatPoisson3d(fmatId,dim);
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
//        TPZLagrangeMultiplier *matskelet = new TPZLagrangeMultiplier(fmatskeleton, fDim-1, 1);
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

TPZCompMesh * LaplaceInCylinder::CMeshPressure(TPZGeoMesh *gmesh, int pOrder, int dim)
{
    /// criar materiais
    //TPZMatPoisson3d *material = new TPZMatPoisson3d(fmatId,dim);
    //TPZMatPoisson3d *material = new TPZMatPoisson3d(fmatId,dim);
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
                if(ftriang) celdisc->SetTotalOrderShape();
                else celdisc->SetTensorialShape();
            }
            
        }
    }
    
#ifdef PZDEBUG
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

TPZCompMesh *LaplaceInCylinder::CMeshMixed(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec)
{
    
    bool condensacaoestatica = true;
    
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
    
    solexata = new TPZDummyFunction<STATE>(SolExata, 5);
    material->SetForcingFunctionExact(solexata);
    mphysics->SetDimModel(dim);
    //funcao do lado direito da equacao do problema
    TPZDummyFunction<STATE> *dum = new TPZDummyFunction<STATE>(Forcing, 5);
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
        TPZAutoPointer<TPZFunction<STATE> > FBCond0 = new TPZDummyFunction<STATE>(ForcingBC0D, 5);
        BCond0 = material->CreateBC(mat, fbc0,fdirichlet, val1, val2);
        BCond0->SetForcingFunction(FBCond0);
    }
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZAutoPointer<TPZFunction<STATE> > FBCond1 = new TPZDummyFunction<STATE>(ForcingBC1D, 5);
    BCond1 = material->CreateBC(mat, fbc1,fdirichlet, val1, val2);
    BCond1->SetForcingFunction(FBCond1);
    //    BCond1 = material->CreateBC(mat, bc1,neumann, val1, val2);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZAutoPointer<TPZFunction<STATE> > FBCond2 = new TPZDummyFunction<STATE>(ForcingBC2D, 5);
    BCond2 = material->CreateBC(mat, fbc2,fdirichlet, val1, val2);
    BCond2->SetForcingFunction(FBCond2);
    //    TPZAutoPointer<TPZFunction<STATE> > FBCond2 = new TPZDummyFunction<STATE>(ForcingBC2N, 5);
    //    BCond2 = material->CreateBC(mat, bc2,neumann, val1, val2);
    //    BCond2->SetForcingFunction(FBCond2);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZAutoPointer<TPZFunction<STATE> > FBCond3 = new TPZDummyFunction<STATE>(ForcingBC3D, 5);
    BCond3 = material->CreateBC(mat, fbc3,fdirichlet, val1, val2);
    BCond3->SetForcingFunction(FBCond3);
    //    BCond3 = material->CreateBC(mat, bc3,neumann, val1, val2);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZAutoPointer<TPZFunction<STATE> > FBCond4 = new TPZDummyFunction<STATE>(ForcingBC4D, 5);
    BCond4 = material->CreateBC(mat, fbc4,fdirichlet, val1, val2);
    BCond4->SetForcingFunction(FBCond4);
    //    TPZAutoPointer<TPZFunction<STATE> > FBCond4 = new TPZDummyFunction<STATE>(ForcingBC4N, 5);
    //    BCond4 = material->CreateBC(mat, bc4,neumann, val1, val2);
    //    BCond4->SetForcingFunction(FBCond4);
    
    if (dim==3)
    {
        val2(0,0) = 0.0;
        val2(1,0) = 0.0;
        TPZAutoPointer<TPZFunction<STATE> > FBCond5 = new TPZDummyFunction<STATE>(ForcingBC5D, 5);
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
    
    
    return mphysics;
    
}



void LaplaceInCylinder::ErrorHDiv(TPZCompMesh *hdivmesh, int p, int ndiv, std::map<REAL, REAL> &fDebugMapL2, std::map<REAL, REAL> &fDebugMapHdiv)
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

void LaplaceInCylinder::ErrorL2(TPZCompMesh *l2mesh, int p, int ndiv, std::map<REAL, REAL> &fDebugMapL2, std::map<REAL, REAL> &fDebugMapHdiv)
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

void LaplaceInCylinder::ErrorH1(TPZCompMesh *l2mesh, int p, int ndiv, std::ostream &out, int DoFT, int DofCond)
{
    
    int64_t nel = l2mesh->NElements();
    int dim = l2mesh->Dimension();
    TPZManVector<REAL,10> globalerrors(10,0.);
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = l2mesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel || gel->Dimension() != dim) {
            continue;
        }
        TPZManVector<REAL,10> elerror(10,0.);
        elerror.Fill(0.);
        cel->EvaluateError(SolExataH1, elerror, false);
        
        int nerr = elerror.size();
        globalerrors.resize(nerr);
        //#ifdef LOG4CXX
        //        if (logger->isDebugEnabled()) {
        //            std::stringstream sout;
        //            sout << "L2 Error sq of element " << el << elerror[0]*elerror[0];
        //            LOGPZ_DEBUG(logger, sout.str())
        //        }
        //#endif
        for (int i=0; i<nerr; i++) {
            globalerrors[i] += elerror[i]*elerror[i];
        }
    }
    
    
    //    int64_t nel = l2mesh->NElements();
    //    //int dim = l2mesh->Dimension();
    //    TPZManVector<STATE,10> globalerrors(10,0.);
    //    for (int64_t el=0; el<nel; el++) {
    //        TPZCompEl *cel = l2mesh->ElementVec()[el];
    //        TPZManVector<STATE,10> elerror(10,0.);
    //        cel->EvaluateError(SolExata, elerror, NULL);
    //        int nerr = elerror.size();
    //        globalerrors.resize(nerr);
    //        //#ifdef LOG4CXX
    //        //        if (logdata->isDebugEnabled()) {
    //        //            std::stringstream sout;
    //        //            sout << "L2 Error sq of element " << el << elerror[0]*elerror[0];
    //        //            LOGPZ_DEBUG(logdata, sout.str())
    //        //        }
    //        //#endif
    //        for (int i=0; i<nerr; i++) {
    //            globalerrors[i] += elerror[i]*elerror[i];
    //        }
    //
    //    }
    //    out << "Errors associated with L2 space - ordem polinomial = " << p << "- divisoes = " << ndiv << endl;
    //    out << "L2 Norm = "    << sqrt(globalerrors[1]) << endl;
    out << ndiv << setw(10) << DoFT << setw(20) << DofCond << setw(28) << sqrt(globalerrors[1]) << setw(35)  << sqrt(globalerrors[2])  << endl;
    
    
    //out << "\nSemi H1 Norm = "    << sqrt(globerrors[2]) << endl;
}


void LaplaceInCylinder::ErrorPrimalDual(TPZCompMesh *l2mesh, TPZCompMesh *hdivmesh,  int p, int ndiv, std::ostream &out, int DoFT, int DofCond)
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

void LaplaceInCylinder::ChangeExternalOrderConnects(TPZCompMesh *mesh){
    
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




