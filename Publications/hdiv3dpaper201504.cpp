//
//  hdiv3dpaper201504.cpp
//  PZ
//
//  Created by Douglas Castro on 22/04/15.
//
//

#include "hdiv3dpaper201504.h"


Hdiv3dPaper201504::Hdiv3dPaper201504()
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
    
    ftetra = false;
    
    fprisma = false;
    
    tetraedra_2.Resize(6, 4);
    
    tetraedra_2(0,0) = 1;
    tetraedra_2(0,1) = 2;
    tetraedra_2(0,2) = 5;
    tetraedra_2(0,3) = 4;
    
    tetraedra_2(1,0) = 4;
    tetraedra_2(1,1) = 7;
    tetraedra_2(1,2) = 3;
    tetraedra_2(1,3) = 2;
    
    tetraedra_2(2,0) = 0;
    tetraedra_2(2,1) = 1;
    tetraedra_2(2,2) = 2;
    tetraedra_2(2,3) = 4;
    
    tetraedra_2(3,0) = 0;
    tetraedra_2(3,1) = 2;
    tetraedra_2(3,2) = 3;
    tetraedra_2(3,3) = 4;

    tetraedra_2(4,0) = 4;
    tetraedra_2(4,1) = 5;
    tetraedra_2(4,2) = 6;
    tetraedra_2(4,3) = 2;
    
    tetraedra_2(5,0) = 4;
    tetraedra_2(5,1) = 6;
    tetraedra_2(5,2) = 7;
    tetraedra_2(5,3) = 2;
    
}

Hdiv3dPaper201504::~Hdiv3dPaper201504()
{
    
}

void Hdiv3dPaper201504::Run(ApproximationSpace problem, Eltype element, TPZVec<int> POrderBeginAndEnd, TPZVec<int> ndivinterval, TPZFMatrix< REAL > &errors)
{
    if (POrderBeginAndEnd.size()!=2)
    {
        std::cout << " POrderBeginAndEnd Vector must have size == 2 and contain just two values, the first and last approximation orders. " << std::endl;
        std::cout << " Some like POrderBeginAndEnd[0] = 1, POrderBeginAndEnd[1] = 3 " << std::endl;
        return;
    }
    
    if (ndivinterval.size()!=2)
    {
        std::cout << " ndivinterval Vector must have size == 2 and contain just two values, the first and last number of uniform refinements. " << std::endl;
        std::cout << " Some like ndivinterval[0] = 0, ndivinterval[1] = 3 for simulations with 0, 1, 2 and 3 refinements " << std::endl;
        std::cout << " Some like ndivinterval[0] = 1, ndivinterval[1] = 1 for simulations with 1 refinement " << std::endl;
        return;
    }
    
    int sizeErrorMatrix = (ndivinterval[1]-ndivinterval[0])*(POrderBeginAndEnd[1] - POrderBeginAndEnd[0]) + 1;
    int errorposition = 0;
    
    errors.Resize(sizeErrorMatrix, 4);
    
    for(int p = POrderBeginAndEnd[0]; p <= POrderBeginAndEnd[1]; p++)
    {
        
        for (int ndiv = ndivinterval[0]; ndiv<=ndivinterval[1]; ndiv++)
        {
            int ordemP = p;
            
            std::cout<< " BEGIN (CUBE Domain) - Polinomial degree: " << ordemP << " rerfinement size (h): " << ndiv << std::endl;
        
            TPZGeoMesh *gmesh;
            switch (element) {
                case ECub:
                {
		    gmesh = this->CreateOneCubo(ndiv);
                    //gmesh = this->CreateOneQuadraticCube(ndiv);
                }
                    break;
                case EPrism:
                {
                    gmesh = this->GMeshWithPrism( ndiv);
                }
                    break;
		case ETetra:
		{
		    double dndiv = ndiv;
		    int nref = (int) pow(2., dndiv);
		    gmesh = this->CreateOneCuboWithTetraedrons(nref);
		}
		  break;
                default:
                    break;
            }


            
            switch (problem) {
                case EH1:
                {
                    
                    gmesh->SetDimension(fDim);
                    TPZCompMesh *cmeshH1 = this->CMeshH1(gmesh, ordemP, fDim);
                    //condensar
                    for (int64_t iel=0; iel<cmeshH1->NElements(); iel++) {
                        TPZCompEl *cel = cmeshH1->Element(iel);
                        if(!cel) continue;
                        new TPZCondensedCompEl(cel);
                    }
                    
                    cmeshH1->ExpandSolution();
                    cmeshH1->CleanUpUnconnectedNodes();
                    
                    TPZAnalysis anh1(cmeshH1, true);
                    
                    SolveSyst(anh1, cmeshH1);
                    
                    ErrorH1(cmeshH1, ordemP, ndiv, errorposition, errors);

                }
                    break;
                case EHDiv:
                {
                    DebugStop(); //Mudanca na topologia
                }
                    break;
                case EHDivStar:
                {

                    TPZCompMesh *cmesh1 = this->CMeshFlux(gmesh, ordemP, fDim);
                    
                    TPZCompMesh *cmesh2 = this->CMeshPressure(gmesh, ordemP, fDim);
                    
                    // P**
                    ChangeExternalOrderConnects(cmesh1);
                    
                    
                    //malha multifisica
                    TPZVec<TPZCompMesh *> meshvec(2);
                    meshvec[0] = cmesh1;
                    meshvec[1] = cmesh2;
                    
                    TPZCompMesh * mphysics = CMeshMixed(gmesh,meshvec);
                    
                    TPZAnalysis an(mphysics, true);
                    
                    SolveSyst(an, mphysics);
                    
                    //Calculo do erro
                    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
                    TPZVec<REAL> erros;
                    
                    std::cout << "Computing Errors\n";
                    
                    ErrorPrimalDual( cmesh2, cmesh1,  ordemP, ndiv, errorposition, errors);
                }
                    break;
                    
                case EHDivStarStar:
                {

                    ordemP  =  ordemP + 1;
                    
                    TPZCompMesh *cmesh1 = this->CMeshFlux(gmesh, ordemP, fDim);
                    
                    TPZCompMesh *cmesh2 = this->CMeshPressure(gmesh, ordemP, fDim);
                    //P**
                    ChangeExternalOrderConnects(cmesh1);
                    
                    //malha multifisica
                    TPZVec<TPZCompMesh *> meshvec(2);
                    meshvec[0] = cmesh1;
                    meshvec[1] = cmesh2;
                    
                    TPZCompMesh * mphysics = CMeshMixed(gmesh,meshvec);
                    
                    TPZAnalysis an(mphysics, true);
                    
                    SolveSyst(an, mphysics);
                    
                    //Calculo do erro
                    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
                    TPZVec<REAL> erros;
                    
                    std::cout << "Computing Errors\n";
                    
                    ErrorPrimalDual( cmesh2, cmesh1,  ordemP, ndiv, errorposition, errors);
                }
                    break;
                    
                default:
                    break;
            }
            
            std::cout<< " END  - polinomial degree: " << ordemP << " refinement size (h): " << ndiv << std::endl;
            
            errorposition++;

        }
    }
    
    std::cout<< " The End " << std::endl;
    
    
}

void Hdiv3dPaper201504::PrintErrors(ApproximationSpace problem, Eltype element, TPZVec<int> POrderBeginAndEnd, TPZVec<int> ndivinterval, TPZVec<REAL> &errors, std::ostream &output)
{
    if (POrderBeginAndEnd.size()!=2)
    {
        std::cout << " POrderBeginAndEnd Vector must have size == 2 and contain just two values, the first and last approximation orders. " << std::endl;
        std::cout << " Some like POrderBeginAndEnd[0] = 1, POrderBeginAndEnd[1] = 3 " << std::endl;
        return;
    }
    
    if (ndivinterval.size()!=2)
    {
        std::cout << " ndivinterval Vector must have size == 2 and contain just two values, the first and last number of uniform refinements. " << std::endl;
        std::cout << " Some like ndivinterval[0] = 0, ndivinterval[1] = 3 for simulations with 0, 1, 2 and 3 refinements " << std::endl;
        std::cout << " Some like ndivinterval[0] = 1, ndivinterval[1] = 1 for simulations with 1 refinement " << std::endl;
        return;
    }
    
    for(int p = POrderBeginAndEnd[0]; p <= POrderBeginAndEnd[1]; p++)
    {
                output << "\n WHEN p = " << p << " \n " << endl;
                output << "ndiv " << setw(6) << "DoFTot" << setw(20) << "DofCond" << setw(28) << "PrimalL2Error" << setw(35) << "L2DualError"  << endl;
        
        for (int ndiv = ndivinterval[0]; ndiv<=ndivinterval[1]; ndiv++)
        {
            int ordemP = p;
            
            std::cout<< " BEGIN (Cube Domain) - Polinomial degree: " << ordemP << " rerfinement size (h): " << ndiv << std::endl;
            
            TPZGeoMesh *gmesh;
            switch (element) {
                case ECub:
                {
		    gmesh = this->CreateOneCubo(ndiv);
                    //gmesh = this->CreateOneQuadraticCube(ndiv);
                }
                    break;
                case EPrism:
                {
                    gmesh = this->GMeshWithPrism( ndiv);
                }
                    break;
		case ETetra:
		{
		    double dndiv = ndiv;
		    int nref = (int) pow(2., dndiv);
		    gmesh = this->CreateOneCuboWithTetraedrons(nref);
		}
		  break;
                default:
                    break;
            }
            
            gmesh->SetDimension(fDim);
            
            switch (problem) {
                case EH1:
                {
                    
                    TPZCompMesh *cmeshH1 = this->CMeshH1(gmesh, ordemP, fDim);
                    
                    int dofTotal = cmeshH1->NEquations();
                    
                    //condensar
                    for (int64_t iel=0; iel<cmeshH1->NElements(); iel++) {
                        TPZCompEl *cel = cmeshH1->Element(iel);
                        if(!cel) continue;
                        new TPZCondensedCompEl(cel);
                    }
                    
                    cmeshH1->ExpandSolution();
                    cmeshH1->CleanUpUnconnectedNodes();
                    
                    int dofCondensed = cmeshH1->NEquations();
                    
                    TPZAnalysis anh1(cmeshH1, true);
                    
                    SolveSyst(anh1, cmeshH1);
                    
                    ErrorH1(cmeshH1, ordemP, ndiv, output, dofTotal, dofCondensed);
                    
                }
                    break;
                case EHDiv:
                {
                    DebugStop(); //Mudanca na topologia
                }
                    break;
                case EHDivStar:
                {

                    TPZCompMesh *cmesh1 = this->CMeshFlux(gmesh, ordemP, fDim);
                    
                    TPZCompMesh *cmesh2 = this->CMeshPressure(gmesh, ordemP, fDim);
                    
                    // para rodar P**
                    ChangeExternalOrderConnects(cmesh1);
                    
                    int DofCond, DoFT;
                    DoFT = cmesh1->NEquations() + cmesh2->NEquations();
                    
                    //malha multifisica
                    TPZVec<TPZCompMesh *> meshvec(2);
                    meshvec[0] = cmesh1;
                    meshvec[1] = cmesh2;
                    
                    TPZCompMesh * mphysics = CMeshMixed(gmesh,meshvec);
                    
                    DofCond = mphysics->NEquations();
                    
                    TPZAnalysis an(mphysics, true);
                    
                    SolveSyst(an, mphysics);
                    
                    //Calculo do erro
                    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
                    TPZVec<REAL> erros;
                    
                    std::cout << "Computing Errors\n";
                    
                    ErrorPrimalDual( cmesh2, cmesh1,  ordemP, ndiv, output, DoFT, DofCond);
                }
                    break;
                    
                case EHDivStarStar:
                {

                    
                    ordemP  =  ordemP + 1;
                    
                    TPZCompMesh *cmesh1 = this->CMeshFlux(gmesh, ordemP, fDim);
                    
                    TPZCompMesh *cmesh2 = this->CMeshPressure(gmesh, ordemP, fDim);
                    
                    // para rodar P**
                    ChangeExternalOrderConnects(cmesh1);
                    
                    int DofCond, DoFT;
                    DoFT = cmesh1->NEquations() + cmesh2->NEquations();
                    
                    //malha multifisica
                    TPZVec<TPZCompMesh *> meshvec(2);
                    meshvec[0] = cmesh1;
                    meshvec[1] = cmesh2;
                    
                    TPZCompMesh * mphysics = CMeshMixed(gmesh,meshvec);
                    
                    DofCond = mphysics->NEquations();
                    
                    TPZAnalysis an(mphysics, true);
                    
                    SolveSyst(an, mphysics);
                    
                    //Calculo do erro
                    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
                    TPZVec<REAL> erros;
                    
                    std::cout << "Computing Errors\n";
                    
                    ErrorPrimalDual( cmesh2, cmesh1,  ordemP, ndiv, output, DoFT, DofCond);
                }
                    break;
                    
                default:
                    break;
            }

            
        }
        
        output << "\n ----------------------------------------------------------------------------- " << endl;
        std::cout<< " END (Quadrilateral Domain) - Polinomial degree: " << p << std::endl;
        
    }
    
    std::cout<< " The End " << std::endl;
    
    
}

TPZGeoMesh *Hdiv3dPaper201504::GMeshWithPrism( int ndiv)
{
    
    int Qnodes =  8;
    
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    gmesh->SetMaxNodeId(Qnodes-1);
    gmesh->NodeVec().Resize(Qnodes);
    TPZVec<TPZGeoNode> Node(Qnodes);
    
    gmesh->SetDimension(3);
    
    TPZVec <int64_t> TopolPrism(6);
    TPZVec <int64_t> TopolQuad(4);
    TPZVec <int64_t> TopolTriang(3);
    
    //indice dos nos
    int64_t id = 0;
    
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
    
    TopolPrism[0] = 0;
    TopolPrism[1] = 1;
    TopolPrism[2] = 2;
    TopolPrism[3] = 4;
    TopolPrism[4] = 5;
    TopolPrism[5] = 6;
    new TPZGeoElRefPattern< pzgeom::TPZGeoPrism > (id, TopolPrism,fmatId,*gmesh);
    id++;
    
    TopolPrism[0] = 0;
    TopolPrism[1] = 2;
    TopolPrism[2] = 3;
    TopolPrism[3] = 4;
    TopolPrism[4] = 6;
    TopolPrism[5] = 7;
    new TPZGeoElRefPattern< pzgeom::TPZGeoPrism > (id, TopolPrism,fmatId,*gmesh);
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

TPZGeoMesh *Hdiv3dPaper201504::CreateOneCuboWithTetraedrons(int64_t nelem)
{
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    GenerateNodes(gmesh,nelem);
    
    for (int64_t i=0; i<nelem; i++) {
        for (int64_t j=0; j<nelem; j++) {
            for (int64_t k=0; k<nelem; k++) {
                TPZManVector<int64_t,8> nodes(8,0);
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
                    TPZManVector<int64_t,4> elnodes(4);
                    int64_t index;
                    for (int il=0; il<4; il++) {
                        elnodes[il] = nodes[tetraedra_2(el,il)];//nodes[tetraedra_2[el][il]];
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
        TPZVec<int64_t> ncoordVec(0);
        int64_t sizeOfVec = 0;
        
        // na face z = 0
        for (int i = 0; i < 4; i++)
        {
            int64_t pos = tetra->NodeIndex(i);
            Nodefinder[i] = gmesh->NodeVec()[pos];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (MyDoubleComparer(nodecoord[2],0.))
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
            int64_t pos = tetra->NodeIndex(i);
            Nodefinder[i] = gmesh->NodeVec()[pos];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (MyDoubleComparer(nodecoord[1],0.))
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
            int64_t pos = tetra->NodeIndex(i);
            Nodefinder[i] = gmesh->NodeVec()[pos];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (MyDoubleComparer(nodecoord[0],1.))
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
            int64_t pos = tetra->NodeIndex(i);
            Nodefinder[i] = gmesh->NodeVec()[pos];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (MyDoubleComparer(nodecoord[1],1.))
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
            int64_t pos = tetra->NodeIndex(i);
            Nodefinder[i] = gmesh->NodeVec()[pos];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (MyDoubleComparer(nodecoord[0],0.))
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
            int64_t pos = tetra->NodeIndex(i);
            Nodefinder[i] = gmesh->NodeVec()[pos];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (MyDoubleComparer(nodecoord[2],1.))
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

void Hdiv3dPaper201504::GenerateNodes(TPZGeoMesh *gmesh, int64_t nelem)
{
    gmesh->NodeVec().Resize((nelem+1)*(nelem+1)*(nelem+1));
    for (int64_t i=0; i<=nelem; i++) {
        for (int64_t j=0; j<=nelem; j++) {
            for (int64_t k=0; k<=nelem; k++) {
                TPZManVector<REAL,3> x(3);
                x[0] = k*1./nelem;
                x[1] = j*1./nelem;
                x[2] = i*1./nelem;
                gmesh->NodeVec()[i*(nelem+1)*(nelem+1)+j*(nelem+1)+k].Initialize(x, *gmesh);
            }
        }
    }
}


TPZGeoMesh *Hdiv3dPaper201504::CreateOneCubo(int nref)
{
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    int nnodes = 8;
    
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
    
    TPZVec<int64_t> TopologyQuad(4);
    
    // bottom
    TopologyQuad[0] = 0;
    TopologyQuad[1] = 1;
    TopologyQuad[2] = 2;
    TopologyQuad[3] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index,TopologyQuad,fbc0,*gmesh);
    index++;
    
    // Front
    TopologyQuad[0] = 0;
    TopologyQuad[1] = 1;
    TopologyQuad[2] = 5;
    TopologyQuad[3] = 4;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index,TopologyQuad,fbc1,*gmesh);
    index++;
    
    // Rigth
    TopologyQuad[0] = 1;
    TopologyQuad[1] = 2;
    TopologyQuad[2] = 6;
    TopologyQuad[3] = 5;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index,TopologyQuad,fbc2,*gmesh);
    index++;
    // Back
    TopologyQuad[0] = 3;
    TopologyQuad[1] = 2;
    TopologyQuad[2] = 6;
    TopologyQuad[3] = 7;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index,TopologyQuad,fbc3,*gmesh);
    index++;
    
    // Left
    TopologyQuad[0] = 0;
    TopologyQuad[1] = 3;
    TopologyQuad[2] = 7;
    TopologyQuad[3] = 4;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index,TopologyQuad,fbc4,*gmesh);
    index++;
    
    // Top
    TopologyQuad[0] = 4;
    TopologyQuad[1] = 5;
    TopologyQuad[2] = 6;
    TopologyQuad[3] = 7;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad>(index,TopologyQuad,fbc5,*gmesh);
    index++;
    
    TPZManVector<int64_t,8> TopolCubo(8,0);
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



void Hdiv3dPaper201504::SolExata(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux){
    
    solp.resize(1);
    solp[0]=0.;
    flux.Resize(3, 1);
    double x = pt[0];
    double y = pt[1];
    double z = pt[2];
    for(int d=0; d<3;d++)
    {
        flux(d,0)=0.;
    }
    
    
    solp[0] = sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);
    flux(0,0) = -M_PI*(cos(M_PI*x)*sin(M_PI*y)*sin(M_PI*z));
    flux(1,0) = -M_PI*(sin(M_PI*x)*cos(M_PI*y)*sin(M_PI*z));
    flux(2,0) = -M_PI*(sin(M_PI*x)*sin(M_PI*y)*cos(M_PI*z));
    
}

void Hdiv3dPaper201504::Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &ff){
    
    double x = pt[0];
    double y = pt[1];
    double z = pt[2];
    
    ff[0] = 3.0*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);
    
}

void Hdiv3dPaper201504::SolExataH1(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux){
    
    solp.resize(1);
    solp[0]=0.;
    flux.Resize(3, 1);
    double x = pt[0];
    double y = pt[1];
    double z = pt[2];
    for(int d=0; d<3;d++)
    {
        flux(d,0)=0.;
    }
    
    
    solp[0] = sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);
    flux(0,0) = M_PI*(cos(M_PI*x)*sin(M_PI*y)*sin(M_PI*z));
    flux(1,0) = M_PI*(sin(M_PI*x)*cos(M_PI*y)*sin(M_PI*z));
    flux(2,0) = M_PI*(sin(M_PI*x)*sin(M_PI*y)*cos(M_PI*z));
    
}

void Hdiv3dPaper201504::ForcingH1(const TPZVec<REAL> &pt, TPZVec<STATE> &ff, TPZFMatrix<STATE> &flux)
{
    
    double x = pt[0];
    double y = pt[1];
    double z = pt[2];
    
    ff[0] = -3.0*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);
    
}


void Hdiv3dPaper201504::ForcingBC0D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){
    
    solp[0] = 0.0;
    
}

void Hdiv3dPaper201504::ForcingBC1D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){
    
    solp[0] = 0.0;
    
}

void Hdiv3dPaper201504::ForcingBC2D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){

    solp[0] = 0.0;
    
}

void Hdiv3dPaper201504::ForcingBC3D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){
    
    solp[0] = 0.0;
    
}

void Hdiv3dPaper201504::ForcingBC4D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){
    
    
    solp[0] = 0.0;

}

void Hdiv3dPaper201504::ForcingBC5D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){
    
    solp[0] = 0.0;
    
}

TPZCompMesh *Hdiv3dPaper201504::CMeshH1(TPZGeoMesh *gmesh, int pOrder, int dim)
{
    /// criar materiais
    TPZMatPoisson3d *material = new TPZMatPoisson3d(fmatId,fDim);
    TPZMaterial * mat(material);
    material->NStateVariables();
    
    //    //solucao exata
    TPZAutoPointer<TPZFunction<STATE> > solexata;
    solexata = new TPZDummyFunction<STATE>(SolExataH1, 5);
    material->SetForcingFunctionExact(solexata);
    
    //funcao do lado direito da equacao do problema
    TPZDummyFunction<STATE> *dum = new TPZDummyFunction<STATE>(ForcingH1, 5);
    TPZAutoPointer<TPZFunction<STATE> > forcef;
    dum->SetPolynomialOrder(20);
    forcef = dum;
    material->SetForcingFunction(forcef);
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(fDim);
    cmesh->SetDefaultOrder(pOrder);
    
    cmesh->InsertMaterialObject(mat);
    
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
    
    BCond0 = material->CreateBC(mat, fbc0,fdirichlet, val1, val2);
    BCond1 = material->CreateBC(mat, fbc1,fdirichlet, val1, val2);
    BCond2 = material->CreateBC(mat, fbc2,fdirichlet, val1, val2);
    BCond3 = material->CreateBC(mat, fbc3,fdirichlet, val1, val2);
    BCond4 = material->CreateBC(mat, fbc4,fdirichlet, val1, val2);
    BCond5 = material->CreateBC(mat, fbc5,fdirichlet, val1, val2);
    
    cmesh->InsertMaterialObject(BCond0); 
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    cmesh->InsertMaterialObject(BCond5); 
    
    
    cmesh->SetAllCreateFunctionsContinuous();
    
    //Ajuste da estrutura de dados computacional
    cmesh->AutoBuild();
    
    return cmesh;
    
}


TPZCompMesh *Hdiv3dPaper201504::CMeshFlux(TPZGeoMesh *gmesh, int pOrder, int dim)
{
    /// criar materiais
    //TPZMatPoisson3d *material = new TPZMatPoisson3d(matId,fDim);
    //TPZMatPoissonD3 *material = new TPZMatPoissonD3(fmatId,fDim);
    TPZMatMixedPoisson3D *material = new TPZMatMixedPoisson3D(fmatId,fDim);
    
    TPZMaterial * mat(material);
    material->NStateVariables();
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    
    
    ///Inserir condicao de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    TPZMaterial * BCond0;
    TPZMaterial * BCond1;
    TPZMaterial * BCond2;
    TPZMaterial * BCond3;
    TPZMaterial * BCond4;
    TPZMaterial * BCond5;

    BCond0 = material->CreateBC(mat, fbc0,fdirichlet, val1, val2);
    BCond1 = material->CreateBC(mat, fbc1,fdirichlet, val1, val2);
    BCond2 = material->CreateBC(mat, fbc2,fdirichlet, val1, val2);
    BCond3 = material->CreateBC(mat, fbc3,fdirichlet, val1, val2);
    BCond4 = material->CreateBC(mat, fbc4,fdirichlet, val1, val2);
    BCond5 = material->CreateBC(mat, fbc5,fdirichlet, val1, val2);
    
    cmesh->InsertMaterialObject(mat);
    
    cmesh->SetDimModel(fDim);
    
    cmesh->SetAllCreateFunctionsHDiv();
    
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    cmesh->InsertMaterialObject(BCond5);

    
    cmesh->SetDefaultOrder(pOrder);
    
    
    TPZLagrangeMultiplier *matskelet = new TPZLagrangeMultiplier(fmatskeleton, fDim-1, 1);
    TPZMaterial * mat2(matskelet);
    cmesh->InsertMaterialObject(mat2);
    
    

    //Ajuste da estrutura de dados computacional
    cmesh->AutoBuild();
    
    return cmesh;
    
}

TPZCompMesh *Hdiv3dPaper201504::CMeshPressure(TPZGeoMesh *gmesh, int pOrder, int dim)
{
    /// criar materiais
    //TPZMatPoisson3d *material = new TPZMatPoisson3d(fmatId,fDim);
    //TPZMatPoissonD3 *material = new TPZMatPoissonD3(fmatId,fDim);
    TPZMatMixedPoisson3D *material = new TPZMatMixedPoisson3D(fmatId,fDim);
    material->NStateVariables();
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(fDim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    
    ///Inserir condicao de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    cmesh->SetDefaultOrder(pOrder);
    
    bool h1function = true;// com esqueleto precisa disso
    if(pOrder>0/*h1function*/){
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
            
            if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
            {
                //if(ftriang==true) celdisc->SetTotalOrderShape();
                //else celdisc->SetTensorialShape();
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
    
    return cmesh;
    
}

TPZCompMesh *Hdiv3dPaper201504::CMeshMixed(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec)
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
    
    solexata = new TPZDummyFunction<STATE>(SolExata, 5);
    material->SetForcingFunctionExact(solexata);
    mphysics->SetDimModel(dim);
    //funcao do lado direito da equacao do problema
    TPZDummyFunction<STATE> *dum = new TPZDummyFunction<STATE>(Forcing, 5);
    TPZAutoPointer<TPZFunction<STATE> > forcef;
    dum->SetPolynomialOrder(10);
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

    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZAutoPointer<TPZFunction<STATE> > FBCond0 = new TPZDummyFunction<STATE>(ForcingBC0D, 5);
    BCond0 = material->CreateBC(mat, fbc0,fdirichlet, val1, val2);
    BCond0->SetForcingFunction(FBCond0);
	
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZAutoPointer<TPZFunction<STATE> > FBCond1 = new TPZDummyFunction<STATE>(ForcingBC1D, 5);
    BCond1 = material->CreateBC(mat, fbc1,fdirichlet, val1, val2);
    BCond1->SetForcingFunction(FBCond1);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZAutoPointer<TPZFunction<STATE> > FBCond2 = new TPZDummyFunction<STATE>(ForcingBC2D, 5);
    BCond2 = material->CreateBC(mat, fbc2,fdirichlet, val1, val2);
    BCond2->SetForcingFunction(FBCond2);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZAutoPointer<TPZFunction<STATE> > FBCond3 = new TPZDummyFunction<STATE>(ForcingBC3D, 5);
    BCond3 = material->CreateBC(mat, fbc3,fdirichlet, val1, val2);
    BCond3->SetForcingFunction(FBCond3);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZAutoPointer<TPZFunction<STATE> > FBCond4 = new TPZDummyFunction<STATE>(ForcingBC4D, 5);
    BCond4 = material->CreateBC(mat, fbc4,fdirichlet, val1, val2);
    BCond4->SetForcingFunction(FBCond4);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZAutoPointer<TPZFunction<STATE> > FBCond5 = new TPZDummyFunction<STATE>(ForcingBC5D, 5);
    BCond5 = material->CreateBC(mat, fbc5,fdirichlet, val1, val2);
    BCond5->SetForcingFunction(FBCond5);
    
    
    
    mphysics->SetAllCreateFunctionsMultiphysicElem();

    mphysics->InsertMaterialObject(BCond0);
    mphysics->InsertMaterialObject(BCond1);
    mphysics->InsertMaterialObject(BCond2);
    mphysics->InsertMaterialObject(BCond3);
    mphysics->InsertMaterialObject(BCond4);
    mphysics->InsertMaterialObject(BCond5);
    
    //Fazendo auto build
    mphysics->AutoBuild();
    mphysics->AdjustBoundaryElements();
    mphysics->CleanUpUnconnectedNodes();
    
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

    return mphysics;
    
}


void Hdiv3dPaper201504::ErrorH1(TPZCompMesh *l2mesh, int p, int ndiv, int pos,  TPZFMatrix< REAL > &errors)
{
    
    int64_t nel = l2mesh->NElements();
    int dim = l2mesh->Dimension();
    TPZManVector<STATE,10> globalerrors(10,0.);
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
        cel->EvaluateError(SolExataH1, elerror, 0);
        
        int nerr = elerror.size();
        globalerrors.resize(nerr);

        for (int i=0; i<nerr; i++) {
            globalerrors[i] += elerror[i]*elerror[i];
        }
    }
    // sequence of storage
    // for each
    errors.PutVal(pos, 0, p);      // polinomial order
    errors.PutVal(pos, 1, ndiv);   // refinement
    errors.PutVal(pos, 2, sqrt(globalerrors[1]));  // pressure
    errors.PutVal(pos, 3, sqrt(globalerrors[2]));  // flux
    
}

void Hdiv3dPaper201504::ErrorH1(TPZCompMesh *l2mesh, int p, int ndiv, std::ostream &out, int DoFT, int DofCond)
{
    
    int64_t nel = l2mesh->NElements();
    int dim = l2mesh->Dimension();
    TPZManVector<STATE,10> globalerrors(10,0.);
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
        cel->EvaluateError(SolExataH1, elerror, 0);
        
        int nerr = elerror.size();
        globalerrors.resize(nerr);
        
        for (int i=0; i<nerr; i++) {
            globalerrors[i] += elerror[i]*elerror[i];
        }
    }
    
    out << ndiv << setw(10) << DoFT << setw(20) << DofCond << setw(28) << sqrt(globalerrors[1]) << setw(35)  << sqrt(globalerrors[2])  << endl;
    
}


void Hdiv3dPaper201504::ErrorPrimalDual(TPZCompMesh *l2mesh, TPZCompMesh *hdivmesh,  int p, int ndiv, int pos, TPZFMatrix< REAL > &errors)
{
    int64_t nel = hdivmesh->NElements();
    int dim = hdivmesh->Dimension();
    TPZManVector<STATE,10> globalerrorsDual(10,0.);
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = hdivmesh->ElementVec()[el];
        if(cel->Reference()->Dimension()!=dim) continue;
        TPZManVector<REAL,10> elerror(10,0.);
        elerror.Fill(0.);
        cel->EvaluateError(SolExata, elerror, 0);
        int nerr = elerror.size();
        for (int i=0; i<nerr; i++) {
            globalerrorsDual[i] += elerror[i]*elerror[i];
        }
        
        
    }
    
    
    nel = l2mesh->NElements();
    
    TPZManVector<STATE,10> globalerrorsPrimal(10,0.);
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = l2mesh->ElementVec()[el];
        TPZManVector<REAL,10> elerror(10,0.);
        cel->EvaluateError(SolExata, elerror, 0);
        int nerr = elerror.size();
        globalerrorsPrimal.resize(nerr);
        
        for (int i=0; i<nerr; i++) {
            globalerrorsPrimal[i] += elerror[i]*elerror[i];
        }
        
    }
    
    // sequence of storage
    // for each
    errors.PutVal(pos, 0, p);      // polinomial order
    errors.PutVal(pos, 1, ndiv);   // refinement
    errors.PutVal(pos, 2, sqrt(globalerrorsPrimal[1]));  // pressure
    errors.PutVal(pos, 3, sqrt(globalerrorsDual[1]));  // flux
    
}

void Hdiv3dPaper201504::ErrorPrimalDual(TPZCompMesh *l2mesh, TPZCompMesh *hdivmesh,  int p, int ndiv, std::ostream &out, int DoFT, int DofCond)
{
    int64_t nel = hdivmesh->NElements();
    int dim = hdivmesh->Dimension();
    TPZManVector<STATE,10> globalerrorsDual(10,0.);
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = hdivmesh->ElementVec()[el];
        if(cel->Reference()->Dimension()!=dim) continue;
        TPZManVector<REAL,10> elerror(10,0.);
        elerror.Fill(0.);
        cel->EvaluateError(SolExata, elerror, 0);
        int nerr = elerror.size();
        for (int i=0; i<nerr; i++) {
            globalerrorsDual[i] += elerror[i]*elerror[i];
        }
        
        
    }
    
    
    nel = l2mesh->NElements();

    TPZManVector<STATE,10> globalerrorsPrimal(10,0.);
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = l2mesh->ElementVec()[el];
        TPZManVector<REAL,10> elerror(10,0.);
        cel->EvaluateError(SolExata, elerror, 0);
        int nerr = elerror.size();
        globalerrorsPrimal.resize(nerr);
        
        for (int i=0; i<nerr; i++) {
            globalerrorsPrimal[i] += elerror[i]*elerror[i];
        }
        
    }
    
    out << ndiv << setw(10) << DoFT << setw(20) << DofCond << setw(28) << sqrt(globalerrorsPrimal[1]) << setw(35)  << sqrt(globalerrorsDual[1])  << endl;
    
}

void Hdiv3dPaper201504::ChangeExternalOrderConnects(TPZCompMesh *mesh){
    
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
                    co.SetOrder(cordermin,cel->ConnectIndex(icon));
                    co.SetNShape(nshape-1);
                    mesh->Block().Set(co.SequenceNumber(),nshape-1);
                }
            }
        }
    }
    mesh->ExpandSolution();
    mesh->CleanUpUnconnectedNodes();
}



void Hdiv3dPaper201504::SolveSyst(TPZAnalysis &an, TPZCompMesh *fCmesh)
{
    std::cout <<"Numero de equacoes "<< fCmesh->NEquations()<< std::endl;
    
    bool isdirect = true;
    bool simetrico = true;
    bool isfrontal = false;
    if (isdirect)
    {
        if (simetrico)
        {
            //TPZSkylineStructMatrix strmat(fCmesh);
            if (isfrontal) {
                TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(fCmesh);
                strmat.SetDecomposeType(ELDLt);
                
                int numthreads = 1;
                
                strmat.SetNumThreads(numthreads);
                
                an.SetStructuralMatrix(strmat);
            }
            else
            {
                //TPZBandStructMatrix full(fCmesh);
                TPZSkylineStructMatrix skylstr(fCmesh); //caso simetrico
                //    TPZSkylineNSymStructMatrix full(fCmesh);
                an.SetStructuralMatrix(skylstr);
            }
            
            
            TPZStepSolver<STATE> step;
            step.SetDirect(ELDLt); //caso simetrico
            an.SetSolver(step);
            an.Run();
        }
        else
        {
            TPZBandStructMatrix full(fCmesh);
            an.SetStructuralMatrix(full);
            TPZStepSolver<STATE> step;
            step.SetDirect(ELU);
            an.SetSolver(step);
            an.Run();
        }
        
    }
    else
    {
        TPZSkylineStructMatrix skylstr(fCmesh); //caso simetrico
        skylstr.SetNumThreads(10);
        an.SetStructuralMatrix(skylstr);
        
        TPZAutoPointer<TPZMatrix<STATE> > matbeingcopied = skylstr.Create();
        TPZAutoPointer<TPZMatrix<STATE> > matClone = matbeingcopied->Clone();
        
        TPZStepSolver<STATE> *precond = new TPZStepSolver<STATE>(matClone);
        TPZStepSolver<STATE> *Solver = new TPZStepSolver<STATE>(matbeingcopied);
        precond->SetReferenceMatrix(matbeingcopied);
        precond->SetDirect(ELDLt);
        Solver->SetGMRES(20, 20, *precond, 1.e-18, 0);
        //        Solver->SetCG(10, *precond, 1.0e-10, 0);
        an.SetSolver(*Solver);
        an.Run();
    }
    
    
    
}
