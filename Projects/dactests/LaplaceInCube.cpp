//
//  LaplaceInCube.cpp
//  PZ
//
//  Created by Douglas Castro on 17/03/15.
//
//

#include "LaplaceInCube.h"
#include "tools.h"



LaplaceInCube::LaplaceInCube(int ordemP, int ndiv, std::map<REAL, REAL> &fDebugMapL2, std::map<REAL, REAL> &fDebugMapHdiv)
{
    std::cout<< " INICIO(CUBO) - grau  polinomio " << ordemP << " numero de divisoes " << ndiv << std::endl;
    std::cout<< " Dimensao == " << fDim << std::endl;
    
    //TPZGeoMesh *gmesh = this->GMeshCirculoGeob(ndiv);
    TPZGeoMesh *gmesh = this->GMeshWithPrism(ndiv);
    
    gmesh->SetDimension(fDim);
    {
        //        ofstream argm("gmesh2d-cubo.txt");
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
    std::string plotname("OurSolutionMetaCubo");
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
    std::string filename("InputDataMetaCubo");
    std::string L2("L2.txt");
    std::string Hdiv("Hdiv.txt");
    std::string HdivData,L2Data;
    HdivData = filename+str+Hdiv;
    L2Data = filename+str+L2;
    
    ErrorHDiv(cmesh1, ordemP, ndiv, fDebugMapL2, fDebugMapHdiv);
    
    ErrorL2(cmesh2, ordemP, ndiv, fDebugMapL2, fDebugMapHdiv);
    
    tools::PrintDebugMapForMathematica(HdivData, L2Data, fDebugMapL2, fDebugMapHdiv);
    
    std::cout<< " FIM (CUBO) - grau  polinomio " << ordemP << " numero de divisoes " << ndiv << std::endl;
    
}

LaplaceInCube::~LaplaceInCube()
{
    
}


TPZGeoMesh *LaplaceInCube::GMeshWithPrism( int ndiv)
{
    
    int Qnodes =  8;
    
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    gmesh->SetMaxNodeId(Qnodes-1);
    gmesh->NodeVec().Resize(Qnodes);
    TPZVec<TPZGeoNode> Node(Qnodes);
    
    gmesh->SetDimension(3);
    
    TPZVec <long> TopolPrism(6);
    TPZVec <long> TopolQuad(4);
    TPZVec <long> TopolTriang(3);
    
    //indice dos nos
    long id = 0;
    
    TPZManVector<REAL,3> coord(3,0.);
    int in = 0;
    //c0
    coord[0] = 0.0;
    coord[1] = 0.0;
    coord[2] = 0.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    //c1
    coord[0] = 1.0;
    coord[1] = 0.0;
    coord[2] = 0.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    //c2
    coord[0] = 1.0;
    coord[1] = 1.0;
    coord[2] = 0.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    //c3
    coord[0] = 0.0;
    coord[1] = 1.0;
    coord[2] = 0.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    //c4
    coord[0] = 0.0;
    coord[1] = 0.0;
    coord[2] = 1.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    //c5
    coord[0] = 1.0;
    coord[1] = 0.0;
    coord[2] = 1.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    //c6
    coord[0] = 1.0;
    coord[1] = 1.0;
    coord[2] = 1.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    //c7
    coord[0] = 0.0;
    coord[1] = 1.0;
    coord[2] = 1.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    
    
    //indice dos elementos
    id = 0;
    
    TopolTriang[0] = 0;
    TopolTriang[1] = 1;
    TopolTriang[2] = 2;
    new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (id,TopolTriang,fbc0,*gmesh);
    id++;
    
    TopolTriang[0] = 0;
    TopolTriang[1] = 2;
    TopolTriang[2] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,fbc0,*gmesh);
    id++;
    
    TopolQuad[0] = 0;
    TopolQuad[1] = 1;
    TopolQuad[2] = 5;
    TopolQuad[3] = 4;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (id, TopolQuad, fbc1, *gmesh);
    id++;
    
    TopolQuad[0] = 1;
    TopolQuad[1] = 2;
    TopolQuad[2] = 6;
    TopolQuad[3] = 5;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (id, TopolQuad, fbc2, *gmesh);
    id++;
    
    TopolQuad[0] = 3;
    TopolQuad[1] = 2;
    TopolQuad[2] = 6;
    TopolQuad[3] = 7;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (id, TopolQuad, fbc3, *gmesh);
    id++;
    
    TopolQuad[0] = 0;
    TopolQuad[1] = 3;
    TopolQuad[2] = 7;
    TopolQuad[3] = 4;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (id, TopolQuad, fbc4, *gmesh);
    id++;
    
    TopolTriang[0] = 4;
    TopolTriang[1] = 5;
    TopolTriang[2] = 6;
    new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (id,TopolTriang,fbc5,*gmesh);
    id++;
    
    TopolTriang[0] = 4;
    TopolTriang[1] = 6;
    TopolTriang[2] = 7;
    new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,fbc5,*gmesh);
    id++;
    
    
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
    
    return gmesh;
}

TPZGeoMesh *LaplaceInCube::CreateOneCuboWithTetraedrons(long nelem)
{
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    GenerateNodes(gmesh,nelem);
    
    for (long i=0; i<nelem; i++) {
        for (long j=0; j<nelem; j++) {
            for (long k=0; k<nelem; k++) {
                TPZManVector<long,8> nodes(8,0);
                nodes[0] = k*(nelem+1)*(nelem+1)+j*(nelem+1)+i;
                nodes[1] = k*(nelem+1)*(nelem+1)+j*(nelem+1)+i+1;
                nodes[2] = k*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i+1;
                nodes[3] = k*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i;
                nodes[4] = (k+1)*(nelem+1)*(nelem+1)+j*(nelem+1)+i;
                nodes[5] = (k+1)*(nelem+1)*(nelem+1)+j*(nelem+1)+i+1;
                nodes[6] = (k+1)*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i+1;
                nodes[7] = (k+1)*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i;

                for (int el=0; el<6; el++)
                {
                    TPZManVector<long,4> elnodes(4);
                    long index;
                    for (int il=0; il<4; il++) {
                        elnodes[il] = nodes[tetraedra_2[el][il]];
                    }
                    gmesh->CreateGeoElement(ETetraedro, elnodes, fmatId, index);
                }
            }
        }
    }
    
    gmesh->BuildConnectivity();
    
    // Boundary Conditions
    const int numelements = gmesh->NElements();
    //    const int bczMinus = -3, bczplus = -2, bcids = -1;
    //    const int bczMinus = -1, bczplus = -1, bcids = -1;
    
    for(int el=0; el<numelements; el++)
    {
        TPZManVector <TPZGeoNode,4> Nodefinder(4);
        TPZManVector <REAL,3> nodecoord(3);
        TPZGeoEl *tetra = gmesh->ElementVec()[el];
        TPZVec<long> ncoordVec(0);
        long sizeOfVec = 0;
        
        // na face z = 0
        for (int i = 0; i < 4; i++)
        {
            long pos = tetra->NodeIndex(i);
            Nodefinder[i] = gmesh->NodeVec()[pos];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (tools::MyDoubleComparer(nodecoord[2],0.))
            {
                sizeOfVec++;
                ncoordVec.Resize(sizeOfVec);
                ncoordVec[sizeOfVec-1] = pos;
            }
        }
        if(ncoordVec.NElements() == 3)
        {
            int lado = tetra->WhichSide(ncoordVec);
            TPZGeoElSide tetraSide(tetra, lado);
            TPZGeoElBC(tetraSide,fbc0);
        }
        
        ncoordVec.clear();
        sizeOfVec = 0;
        // na face y = 0
        for (int i = 0; i < 4; i++)
        {
            long pos = tetra->NodeIndex(i);
            Nodefinder[i] = gmesh->NodeVec()[pos];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (tools::MyDoubleComparer(nodecoord[1],0.))
            {
                sizeOfVec++;
                ncoordVec.Resize(sizeOfVec);
                ncoordVec[sizeOfVec-1] = pos;
            }
        }
        if(ncoordVec.NElements() == 3)
        {
            int lado = tetra->WhichSide(ncoordVec);
            TPZGeoElSide tetraSide(tetra, lado);
            TPZGeoElBC(tetraSide,fbc1);
        }
        
        ncoordVec.clear();
        sizeOfVec = 0;
        // na face x = 1
        for (int i = 0; i < 4; i++)
        {
            long pos = tetra->NodeIndex(i);
            Nodefinder[i] = gmesh->NodeVec()[pos];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (tools::MyDoubleComparer(nodecoord[0],1.))
            {
                sizeOfVec++;
                ncoordVec.Resize(sizeOfVec);
                ncoordVec[sizeOfVec-1] = pos;
            }
        }
        if(ncoordVec.NElements() == 3)
        {
            int lado = tetra->WhichSide(ncoordVec);
            TPZGeoElSide tetraSide(tetra, lado);
            TPZGeoElBC(tetraSide,fbc2);
        }
        
        ncoordVec.clear();
        sizeOfVec = 0;
        // na face y = 1
        for (int i = 0; i < 4; i++)
        {
            long pos = tetra->NodeIndex(i);
            Nodefinder[i] = gmesh->NodeVec()[pos];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (tools::MyDoubleComparer(nodecoord[1],1.))
            {
                sizeOfVec++;
                ncoordVec.Resize(sizeOfVec);
                ncoordVec[sizeOfVec-1] = pos;
            }
        }
        if(ncoordVec.NElements() == 3)
        {
            int lado = tetra->WhichSide(ncoordVec);
            TPZGeoElSide tetraSide(tetra, lado);
            TPZGeoElBC(tetraSide,fbc3);
        }
        
        
        ncoordVec.clear();
        sizeOfVec = 0;
        // na face x = 0
        for (int i = 0; i < 4; i++)
        {
            long pos = tetra->NodeIndex(i);
            Nodefinder[i] = gmesh->NodeVec()[pos];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (tools::MyDoubleComparer(nodecoord[0],0.))
            {
                sizeOfVec++;
                ncoordVec.Resize(sizeOfVec);
                ncoordVec[sizeOfVec-1] = pos;
            }
        }
        if(ncoordVec.NElements() == 3)
        {
            int lado = tetra->WhichSide(ncoordVec);
            TPZGeoElSide tetraSide(tetra, lado);
            TPZGeoElBC(tetraSide,fbc4);
        }
        
        ncoordVec.clear();
        sizeOfVec = 0;
        // na face z = 1
        for (int i = 0; i < 4; i++)
        {
            long pos = tetra->NodeIndex(i);
            Nodefinder[i] = gmesh->NodeVec()[pos];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (tools::MyDoubleComparer(nodecoord[2],1.))
            {
                sizeOfVec++;
                ncoordVec.Resize(sizeOfVec);
                ncoordVec[sizeOfVec-1] = pos;
            }
        }
        if(ncoordVec.NElements() == 3)
        {
            int lado = tetra->WhichSide(ncoordVec);
            TPZGeoElSide tetraSide(tetra, lado);
            TPZGeoElBC(tetraSide,fbc5);
        }
        
        
        
    }
    
    return gmesh;
}

void LaplaceInCube::GenerateNodes(TPZGeoMesh *gmesh, long nelem)
{
    gmesh->NodeVec().Resize((nelem+1)*(nelem+1)*(nelem+1));
    for (long i=0; i<=nelem; i++) {
        for (long j=0; j<=nelem; j++) {
            for (long k=0; k<=nelem; k++) {
                TPZManVector<REAL,3> x(3);
                x[0] = k*1./nelem;
                x[1] = j*1./nelem;
                x[2] = i*1./nelem;
                gmesh->NodeVec()[i*(nelem+1)*(nelem+1)+j*(nelem+1)+k].Initialize(x, *gmesh);
            }
        }
    }
}


TPZGeoMesh *LaplaceInCube::CreateOneCubo(int nref)
{
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    int nnodes = 8;
    int idf0=-1;
    int idf1=-2;
    int idf2=-3;
    int idf3=-4;
    int idf4=-5;
    int idf5=-6;
    
    gmesh->SetDimension(3);
    gmesh->NodeVec().Resize(nnodes);
    
    TPZManVector<REAL,3> coord(3,0.);
    int in = 0;
    
    //cubo [0,1]ˆ3
    //c0
    coord[0] = 0.0;
    coord[1] = 0.0;
    coord[2] = 0.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    //c1
    coord[0] = 1.0;
    coord[1] = 0.0;
    coord[2] = 0.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    //c2
    coord[0] = 1.0;
    coord[1] = 1.0;
    coord[2] = 0.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    //c3
    coord[0] = 0.0;
    coord[1] = 1.0;
    coord[2] = 0.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    
    //c4
    coord[0] = 0.0;
    coord[1] = 0.0;
    coord[2] = 1.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    //c5
    coord[0] = 1.0;
    coord[1] = 0.0;
    coord[2] = 1.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    //c6
    coord[0] = 1.0;
    coord[1] = 1.0;
    coord[2] = 1.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    //c7
    coord[0] = 0.0;
    coord[1] = 1.0;
    coord[2] = 1.0;
    gmesh->NodeVec()[in].SetCoord(coord);
    gmesh->NodeVec()[in].SetNodeId(in);
    in++;
    
    // cubo [-1,1]^3
    //    //c0
    //    coord[0] = -1.0;
    //    coord[1] = -1.0;
    //    coord[2] = -1.0;
    //    gmesh->NodeVec()[in].SetCoord(coord);
    //    gmesh->NodeVec()[in].SetNodeId(in);
    //    in++;
    //    //c1
    //    coord[0] =  1.0;
    //    coord[1] = -1.0;
    //    coord[2] = -1.0;
    //    gmesh->NodeVec()[in].SetCoord(coord);
    //    gmesh->NodeVec()[in].SetNodeId(in);
    //    in++;
    //    //c2
    //    coord[0] =  1.0;
    //    coord[1] =  1.0;
    //    coord[2] = -1.0;
    //    gmesh->NodeVec()[in].SetCoord(coord);
    //    gmesh->NodeVec()[in].SetNodeId(in);
    //    in++;
    //    //c3
    //    coord[0] = -1.0;
    //    coord[1] =  1.0;
    //    coord[2] = -1.0;
    //    gmesh->NodeVec()[in].SetCoord(coord);
    //    gmesh->NodeVec()[in].SetNodeId(in);
    //    in++;
    //    //c4
    //    coord[0] = -1.0;
    //    coord[1] = -1.0;
    //    coord[2] =  1.0;
    //    gmesh->NodeVec()[in].SetCoord(coord);
    //    gmesh->NodeVec()[in].SetNodeId(in);
    //    in++;
    //    //c5
    //    coord[0] =  1.0;
    //    coord[1] = -1.0;
    //    coord[2] =  1.0;
    //    gmesh->NodeVec()[in].SetCoord(coord);
    //    gmesh->NodeVec()[in].SetNodeId(in);
    //    in++;
    //    //c6
    //    coord[0] =  1.0;
    //    coord[1] =  1.0;
    //    coord[2] =  1.0;
    //    gmesh->NodeVec()[in].SetCoord(coord);
    //    gmesh->NodeVec()[in].SetNodeId(in);
    //    in++;
    //    //c7
    //    coord[0] = -1.0;
    //    coord[1] =  1.0;
    //    coord[2] =  1.0;
    //    gmesh->NodeVec()[in].SetCoord(coord);
    //    gmesh->NodeVec()[in].SetNodeId(in);
    //    in++;
    
    
    
    int index = 0;
    
    TPZVec<long> TopologyQuad(4);
    
    // bottom
    TopologyQuad[0] = 0;
    TopologyQuad[1] = 1;
    TopologyQuad[2] = 2;
    TopologyQuad[3] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index,TopologyQuad,idf0,*gmesh);
    index++;
    
    // Front
    TopologyQuad[0] = 0;
    TopologyQuad[1] = 1;
    TopologyQuad[2] = 5;
    TopologyQuad[3] = 4;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index,TopologyQuad,idf1,*gmesh);
    index++;
    
    // Rigth
    TopologyQuad[0] = 1;
    TopologyQuad[1] = 2;
    TopologyQuad[2] = 6;
    TopologyQuad[3] = 5;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index,TopologyQuad,idf2,*gmesh);
    index++;
    // Back
    TopologyQuad[0] = 3;
    TopologyQuad[1] = 2;
    TopologyQuad[2] = 6;
    TopologyQuad[3] = 7;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index,TopologyQuad,idf3,*gmesh);
    index++;
    
    // Left
    TopologyQuad[0] = 0;
    TopologyQuad[1] = 3;
    TopologyQuad[2] = 7;
    TopologyQuad[3] = 4;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index,TopologyQuad,idf4,*gmesh);
    index++;
    
    // Top
    TopologyQuad[0] = 4;
    TopologyQuad[1] = 5;
    TopologyQuad[2] = 6;
    TopologyQuad[3] = 7;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index,TopologyQuad,idf5,*gmesh);
    index++;
    
    TPZManVector<long,8> TopolCubo(8,0);
    TopolCubo[0] = 0;
    TopolCubo[1] = 1;
    TopolCubo[2] = 2;
    TopolCubo[3] = 3;
    TopolCubo[4] = 4;
    TopolCubo[5] = 5;
    TopolCubo[6] = 6;
    TopolCubo[7] = 7;
    
    
    new TPZGeoElRefPattern< pzgeom::TPZGeoCube> (index, TopolCubo, fmatId, *gmesh);
    index++;
    
    
    gmesh->BuildConnectivity();
    
    /// gmesh para aqui
    
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
    
    
    std::ofstream out("SingleCubeWithBcs.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
    
    return gmesh;
}



void LaplaceInCube::SolExata(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux){
    
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
    
    flux(0,0)=flux(1,0)=0.;
    double x = pt[0];
    double y = pt[1];
    REAL raio = sqrt( x*x + y*y );
    if (raio < 1.0)
    {
        
        
        // para a Lap p = f no circulo
        solp[0] = 3.0*exp(1.0/(x*x + y*y-1.0));
        flux(0,0) = 6.0*exp(1.0/(x*x + y*y-1.0))*(x*TP(0,0)+y*TP(0,1))/( (x*x + y*y-1.0)*(x*x + y*y-1.0) );
        flux(1,0) = 6.0*exp(1.0/(x*x + y*y-1.0))*(x*TP(1,0)+y*TP(1,1))/( (x*x + y*y-1.0)*(x*x + y*y-1.0) );
        
    }
    else
    {
        // para a Lap p = f no circulo
        solp[0] = 0.0;
        flux(0,0)= 0.0;
        flux(1,0)= 0.0;
        
    }
    
#ifdef LINEAR
    
    double x = pt[0];
    double y = pt[1];
    REAL r = sqrt( x*x + y*y );
    solp[0] = r*r*(1.0 - r*r);
    flux(0,0)= (2.0*r - 4.0*r*r*r)*cos(theta);
    flux(1,0)=  (2.0*r - 4.0*r*r*r)*sin(theta);
    
#endif
    
    
}

void LaplaceInCube::Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &ff){
    
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
    REAL raio = sqrt( x*x + y*y );
    if (raio < 1.0)
    {
        ff[0] = -6.0*exp(1.0/(x*x + y*y-1.0))*((-1.0 + 3.0*x*x*x*x + 2.0*(1.0 + x*x)*y*y - y*y*y*y)*TP(0,0) +2.0*x*y*(-1.0 + 2.0*x*x + 2.0*y*y)*(TP(0,1) + TP(1,0)) + (-1.0 - x*x*x*x +3.0*y*y*y*y + 2.0*x*x*(1.0 + y*y))*TP(1,1));
        
    }
    else
    {
        ff[0] = 0.0;
        
    }
    
#ifdef LINEAR
    
    double x = pt[0];
    double y = pt[1];
    REAL r = sqrt( x*x + y*y );
    ff[0] = -(4.0 - 16.0*r*r);
    
#endif
    
}

void LaplaceInCube::ForcingH1(const TPZVec<REAL> &pt, TPZVec<STATE> &ff, TPZFMatrix<STATE> &flux)
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

void LaplaceInCube::ForcingBC0D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){
    
    
    double x = pt[0];
    double y = pt[1];
    REAL raio = sqrt( x*x + y*y );
    if (raio < 1.0)
    {
        // para a Lap p = f no circulo
        solp[0] = 3.0*exp(1.0/(x*x + y*y-1.0));
        
    }
    else
    {
        solp[0] = 0.0;
    }
    
    
#ifdef LINEAR
    
    solp[0] = 0.0;
    
#endif
    
    
    
}

void LaplaceInCube::ForcingBC1D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){
    
    double x = pt[0];
    double y = pt[1];
    REAL raio = sqrt( x*x + y*y );
    if (raio < 1.0)
    {
        // para a Lap p = f no circulo
        solp[0] = 3.0*exp(1.0/(x*x + y*y-1.0));
        
    }
    else
    {
        solp[0] = 0.0;
    }
    
#ifdef LINEAR
    
    solp[0] = 0.0;
    
#endif
}

void LaplaceInCube::ForcingBC2D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){
    
    
    
    double x = pt[0];
    double y = pt[1];
    REAL raio = sqrt( x*x + y*y );
    if (raio < 1.0)
    {
        // para a Lap p = f no circulo
        solp[0] = 3.0*exp(1.0/(x*x + y*y-1.0));
        
    }
    else
    {
        solp[0] = 0.0;
    }
#ifdef LINEAR
    
    solp[0] = 0.0;
    
#endif
    
}

void LaplaceInCube::ForcingBC3D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){
    
    
    double x = pt[0];
    double y = pt[1];
    REAL raio = sqrt( x*x + y*y );
    if (raio < 1.0)
    {
        // para a Lap p = f no circulo
        solp[0] = 3.0*exp(1.0/(x*x + y*y-1.0));
        
    }
    else
    {
        solp[0] = 0.0;
    }
#ifdef LINEAR
    
    solp[0] = 0.0;
    
#endif
    
}

void LaplaceInCube::ForcingBC4D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){
    
    
    double x = pt[0];
    double y = pt[1];
    REAL raio = sqrt( x*x + y*y );
    if (raio < 1.0)
    {
        // para a Lap p = f no circulo
        solp[0] = 3.0*exp(1.0/(x*x + y*y-1.0));
        
    }
    else
    {
        solp[0] = 0.0;
    }
#ifdef LINEAR
    
    solp[0] = 0.0;
    
#endif
}

void LaplaceInCube::ForcingBC5D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){
    
    
    
    double x = pt[0];
    double y = pt[1];
    REAL raio = sqrt( x*x + y*y );
    if (raio < 1.0)
    {
        // para a Lap p = f no circulo
        solp[0] = 3.0*exp(1.0/(x*x + y*y-1.0));
        
    }
    else
    {
        solp[0] = 0.0;
    }
#ifdef LINEAR
    
    solp[0] = 0.0;
    
#endif
}


void LaplaceInCube::ForcingBC0N(const TPZVec<REAL> &pt, TPZVec<STATE> &normflux){
    
    DebugStop();
}

void LaplaceInCube::ForcingBC1N(const TPZVec<REAL> &pt, TPZVec<STATE> &normflux){
    
    DebugStop();
    
}

void LaplaceInCube::ForcingBC2N(const TPZVec<REAL> &pt, TPZVec<STATE> &normflux){
    DebugStop();
}

void LaplaceInCube::ForcingBC3N(const TPZVec<REAL> &pt, TPZVec<STATE> &normflux){
    
    DebugStop();
    
}

void LaplaceInCube::ForcingBC4N(const TPZVec<REAL> &pt, TPZVec<STATE> &normflux){
    
    DebugStop();
}

void LaplaceInCube::ForcingBC5N(const TPZVec<REAL> &pt, TPZVec<STATE> &normflux){
    
    DebugStop();
}

TPZCompMesh *LaplaceInCube::CMeshH1(TPZGeoMesh *gmesh, int pOrder, int dim)
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


TPZCompMesh *LaplaceInCube::CMeshFlux(TPZGeoMesh *gmesh, int pOrder, int dim)
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

TPZCompMesh *LaplaceInCube::CMeshPressure(TPZGeoMesh *gmesh, int pOrder, int dim)
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

TPZCompMesh *LaplaceInCube::CMeshMixed(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec)
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


void LaplaceInCube::ErrorHDiv(TPZCompMesh *hdivmesh, int p, int ndiv, std::map<REAL, REAL> &fDebugMapL2, std::map<REAL, REAL> &fDebugMapHdiv)
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

void LaplaceInCube::ErrorL2(TPZCompMesh *l2mesh, int p, int ndiv, std::map<REAL, REAL> &fDebugMapL2, std::map<REAL, REAL> &fDebugMapHdiv)
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






