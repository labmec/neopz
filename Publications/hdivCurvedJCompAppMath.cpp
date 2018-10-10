//
//  hdivCurvedJCompAppMath.cpp
//  PZ
//
//  Created by Douglas Castro on 18/06/15.
//
//

#include "hdivCurvedJCompAppMath.h"

const int  norder = 6;

bool hdivCurvedJCompAppMath::probAtCircle = false;
bool hdivCurvedJCompAppMath::probAtCylinder = false;
bool hdivCurvedJCompAppMath::probAtSphere = false;


hdivCurvedJCompAppMath::hdivCurvedJCompAppMath()
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
    
    probAtCircle = false;
    probAtCylinder = false;
    probAtSphere = false;

    ftriang = false;
    
    isgeoblend = true;
}


hdivCurvedJCompAppMath::hdivCurvedJCompAppMath(geomDomain geodomain)
{
  switch(geodomain){
    case ECircle:
    {
      probAtCircle = true;
      probAtCylinder = false;
      probAtSphere = false;
    }
    break;
    case ECylinder:
    {
      probAtCircle = false;
      probAtCylinder = true;
      probAtSphere = false;
    }
    break;
    case ESphere:
    {
      probAtCircle = false;
      probAtCylinder = false;
      probAtSphere = true;
    }
    break;
    default:
      break;
    
  }
//   switch(geodomain){
//     case ECircle:
//     {
//       probAtCircle = 1;
//       probAtCylinder = 0;
//       probAtSphere = 0;
//     }
//     break;
//     case ECylinder:
//     {
//       probAtCircle = 0;
//       probAtCylinder = 1;
//       probAtSphere = 0;
//     }
//     break;
//     case ESphere:
//     {
//       probAtCircle = 0;
//       probAtCylinder = 0;
//       probAtSphere = 1;
//     }
//     break;
//     default:
//       break;
//     
//   }
}

hdivCurvedJCompAppMath::~hdivCurvedJCompAppMath()
{
    
}

void hdivCurvedJCompAppMath::Run(geomDomain geodomain, ApproximationSpace problem, Eltype element, TPZVec<int> POrderBeginAndEnd, TPZVec<int> ndivinterval, TPZFMatrix< REAL > &errors)
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
            
            std::cout<< " BEGIN - Polinomial degree: " << ordemP << " rerfinement size (h): " << ndiv << std::endl;
        
            TPZGeoMesh *gmesh;
	    /*
            switch (element) {
                case EQuad:
                {
                    gmesh = this->GMesh(fDim, false, ndiv);
                }
                    break;
                case ETriangle:
                {
                    setTriangTrue();
                    gmesh = this->GMesh(fDim, true, ndiv);
                }
                    break;
                default:
                    break;
            }*/
            
	    hdivCurvedJCompAppMath geodomain;
	    if(probAtCircle)
	    {
	      gmesh = this->MakeCircle(ndiv);
	    }
	    else if(probAtCylinder)
	    {
	      gmesh = this->GMeshCilindricalMesh(ndiv);;
                
	    }  
	    else //if(probAtSphere)
	    {
	      gmesh = this->MakeSphereFromQuadrilateral(2, false, ndiv);
	    }
	      
//             switch (geodomain) {
//                 case ECircle:
//                 {
//                     gmesh = this->MakeCircle(ndiv);
//                 }
//                     break;
//                 case ECylinder:
//                 {
// 		  gmesh = this->GMeshCilindricalMesh(ndiv);;
//                 }
//                     break;
// 		case ESphere:
// 		{
// 		  gmesh = this->MakeSphereFromQuadrilateral(2, false, ndiv);
// 		}
// 		    break;
//                 default:
//                     break;
//             }
            
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
            
            std::cout<< " END - polinomial degree: " << ordemP << " refinement size (h): " << ndiv << std::endl;
            
            errorposition++;

        }
    }
    
    std::cout<< " The End " << std::endl;
    
    
}

void hdivCurvedJCompAppMath::PrintErrors(geomDomain geodomain, ApproximationSpace problem, Eltype element, TPZVec<int> POrderBeginAndEnd, TPZVec<int> ndivinterval, TPZVec<REAL> &errors, std::ostream &output)
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
            
            std::cout<< " BEGIN (Quadrilateral Domain) - Polinomial degree: " << ordemP << " rerfinement size (h): " << ndiv << std::endl;
            
            TPZGeoMesh *gmesh;
	    
	    hdivCurvedJCompAppMath geodomain;

	    if(probAtCircle)
	    {
	      gmesh = this->MakeCircle(ndiv);
	    }
	    else if(probAtCylinder)
	    {
	      gmesh = this->GMeshCilindricalMesh(ndiv);;
                
	    }  
	    else //if(probAtSphere)
	    {
	      gmesh = this->MakeSphereFromQuadrilateral(2, false, ndiv);
	    }
// 	     switch (geodomain) {
//                 case ECircle:
//                 {
//                     gmesh = this->MakeCircle( ndiv);
//                 }
//                     break;
//                 case ECylinder:
//                 {
// 		  gmesh = this->GMeshCilindricalMesh(ndiv);;
//                 }
//                     break;
// 		case ESphere:
// 		{
// 		  gmesh = this->MakeSphereFromQuadrilateral(2, false, ndiv);
// 		}
// 		    break;
//                 default:
//                     break;
//             }
            
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

TPZGeoMesh *hdivCurvedJCompAppMath::MakeCircle( int ndiv)
{
  
    if( fDim!= 2)
    {
        DebugStop();
    }
    
    TPZGeoMesh * geomesh = new TPZGeoMesh;
    
    int nodes =  17 + 8;
    REAL radius = 1.0;
    REAL innerradius = radius/2.0;
    geomesh->SetMaxNodeId(nodes-1);
    geomesh->NodeVec().Resize(nodes);
    TPZManVector<TPZGeoNode,7> Node(nodes);
    
    TPZManVector<int64_t,8> TopolQQuadrilateral(8);
    TPZManVector<int64_t,8> TopolQuadrilateral(4);
    TPZManVector<int64_t,6> TopolQTriangle(6);
    TPZManVector<int64_t,2> TopolLine(2);
    TPZManVector<int64_t,3> TopolArc(3);
    TPZManVector<REAL,3> coord(3,0.);
    TPZVec<REAL> xc(3,0.);
    
    
    int nodeindex = 0;
    
    for (int inode = 0; inode < 8 ; inode++) {
        // i node
        coord = ParametricCircle(innerradius, inode * M_PI/4.0);
        geomesh->NodeVec()[nodeindex].SetCoord(coord);
        geomesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
        nodeindex++;
    }

    for (int inode = 0; inode < 8 ; inode++) {
        // i node
        coord = ParametricCircle(radius, inode * M_PI/4.0);
        geomesh->NodeVec()[nodeindex].SetCoord(coord);
        geomesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
        nodeindex++;
    }
    
    for (int inode = 0; inode < 4 ; inode++) {
        // i node
        coord = ParametricCircle(0.5*innerradius, inode * M_PI/2.0);
        geomesh->NodeVec()[nodeindex].SetCoord(coord);
        geomesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
        nodeindex++;
    }
    
    for (int inode = 0; inode < 4 ; inode++) {
        // i node
        coord = ParametricCircle(1.5*innerradius, inode * M_PI/2.0);
        geomesh->NodeVec()[nodeindex].SetCoord(coord);
        geomesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
        nodeindex++;
    }
    
    // 24 node id at circle center
    coord = ParametricCircle(0.0,0.0);
    geomesh->NodeVec()[nodeindex].SetCoord(coord);
    geomesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    int elementid = 0;
    
    // outer domain
    
    TopolQuadrilateral[0] = 0;
    TopolQuadrilateral[1] = 8;
    TopolQuadrilateral[2] = 10;
    TopolQuadrilateral[3] = 2;
    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad > > (elementid,TopolQuadrilateral, fmatId,*geomesh);
    elementid++;
    
    TopolQuadrilateral[0] = 2;
    TopolQuadrilateral[1] = 10;
    TopolQuadrilateral[2] = 12;
    TopolQuadrilateral[3] = 4;
    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad > > (elementid,TopolQuadrilateral, fmatId,*geomesh);
    elementid++;
    
    TopolQuadrilateral[0] = 4;
    TopolQuadrilateral[1] = 12;
    TopolQuadrilateral[2] = 14;
    TopolQuadrilateral[3] = 6;

    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad > > (elementid,TopolQuadrilateral, fmatId,*geomesh);
    elementid++;
    
    TopolQuadrilateral[0] = 6;
    TopolQuadrilateral[1] = 14;
    TopolQuadrilateral[2] = 8;
    TopolQuadrilateral[3] = 0;
    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad > > (elementid,TopolQuadrilateral, fmatId,*geomesh);
    elementid++;
    
    
    // inner domain
    
    
    TopolQTriangle[0] = 0;
    TopolQTriangle[1] = 2;
    TopolQTriangle[2] = 24;
    TopolQTriangle[3] = 1;
    TopolQTriangle[4] = 17;
    TopolQTriangle[5] = 16;
    new TPZGeoElRefPattern<  pzgeom::TPZQuadraticTrig  > (elementid,TopolQTriangle, fmatId,*geomesh);
    elementid++;

    TopolQTriangle[0] = 2;
    TopolQTriangle[1] = 4;
    TopolQTriangle[2] = 24;
    TopolQTriangle[3] = 3;
    TopolQTriangle[4] = 18;
    TopolQTriangle[5] = 17;
    new TPZGeoElRefPattern<  pzgeom::TPZQuadraticTrig  > (elementid,TopolQTriangle, fmatId,*geomesh);
    elementid++;

    TopolQTriangle[0] = 4;
    TopolQTriangle[1] = 6;
    TopolQTriangle[2] = 24;
    TopolQTriangle[3] = 5;
    TopolQTriangle[4] = 19;
    TopolQTriangle[5] = 18;
    new TPZGeoElRefPattern<  pzgeom::TPZQuadraticTrig  > (elementid,TopolQTriangle, fmatId,*geomesh);
    elementid++;
    
    TopolQTriangle[0] = 6;
    TopolQTriangle[1] = 0;
    TopolQTriangle[2] = 24;
    TopolQTriangle[3] = 7;
    TopolQTriangle[4] = 16;
    TopolQTriangle[5] = 19;
    new TPZGeoElRefPattern<  pzgeom::TPZQuadraticTrig  > (elementid,TopolQTriangle, fmatId,*geomesh);
    elementid++;
    
    
    
    // outer arcs bc's
    
    TopolArc[0] = 8;
    TopolArc[1] = 10;
    TopolArc[2] = 9;
    new TPZGeoElRefPattern< pzgeom::TPZArc3D > (elementid,TopolArc, fbc1,*geomesh);
    elementid++;

    TopolArc[0] = 10;
    TopolArc[1] = 12;
    TopolArc[2] = 11;
    new TPZGeoElRefPattern< pzgeom::TPZArc3D > (elementid,TopolArc, fbc1,*geomesh);
    elementid++;
    
    TopolArc[0] = 12;
    TopolArc[1] = 14;
    TopolArc[2] = 13;
    new TPZGeoElRefPattern< pzgeom::TPZArc3D > (elementid,TopolArc, fbc1,*geomesh);
    elementid++;
    
    TopolArc[0] = 14;
    TopolArc[1] = 8;
    TopolArc[2] = 15;
    new TPZGeoElRefPattern< pzgeom::TPZArc3D > (elementid,TopolArc, fbc1,*geomesh);
    elementid++;
    
    
    
    geomesh->BuildConnectivity();
    
    int nref = ndiv;
    TPZVec<TPZGeoEl *> sons;
    for (int iref = 0; iref < nref; iref++) {
        int nel = geomesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = geomesh->ElementVec()[iel];
            if (gel->HasSubElement()) {
                continue;
            }
            gel->Divide(sons);
        }
    }
    
    
    std::ofstream out("CircleMixed.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(geomesh, out, true);
    
    return geomesh;
}

TPZManVector<REAL,3> hdivCurvedJCompAppMath::ParametricCircle(REAL radius,REAL theta)
{
    TPZManVector<REAL,3> xcoor(3,0.0);
    xcoor[0] = radius * cos(theta);
    xcoor[1] = radius * sin(theta);
    xcoor[2] = 0.0 ;
    return xcoor;
}

TPZGeoMesh *hdivCurvedJCompAppMath::MakeSphereFromQuadrilateral(int dimensao, bool triang, int ndiv)
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
    
//     std::ofstream outfile("NiceSphere.vtk");
//     TPZVTKGeoMesh::PrintGMeshVTK(geomesh, outfile, true);
//     
    return geomesh;
}

void hdivCurvedJCompAppMath::RotateGeomesh(TPZGeoMesh *gmesh, REAL CounterClockwiseAngle, int &Axis)
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



TPZManVector<REAL,3> hdivCurvedJCompAppMath::ParametricSphere(REAL radius,REAL phi,REAL theta)
{
    TPZManVector<REAL,3> xcoor(3,0.0);
    xcoor[0] = radius * cos(theta) * sin(phi);
    xcoor[1] = radius * sin(theta) * sin(phi);
    xcoor[2] = radius * cos(phi) ;
    return xcoor;
}

TPZGeoMesh *hdivCurvedJCompAppMath::GMeshCilindricalMesh( int ndiv)
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
    RotateNode(coord, angrot, Axis);
    
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,coord[0]);//coord X
    gmesh->NodeVec()[id].SetCoord(1,coord[1]);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,coord[2]);//coord Z
    id++;
    
    //1
    coord[0] = r*cos(theta1);
    coord[1] = r*sin(theta1);
    coord[2] = z0;
    RotateNode(coord, angrot, Axis);
    
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,coord[0]);//coord X
    gmesh->NodeVec()[id].SetCoord(1,coord[1]);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,coord[2]);//coord Z
    id++;
    
    //2
    coord[0] = r*cos(theta1);
    coord[1] = r*sin(theta1);
    coord[2] = z1;
    RotateNode(coord, angrot, Axis);
    
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,coord[0]);//coord X
    gmesh->NodeVec()[id].SetCoord(1,coord[1]);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,coord[2]);//coord Z
    id++;
    
    //3
    coord[0] = r*cos(theta0);
    coord[1] = r*sin(theta0);
    coord[2] = z1;
    RotateNode(coord, angrot, Axis);
    
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,coord[0]);//coord X
    gmesh->NodeVec()[id].SetCoord(1,coord[1]);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,coord[2]);//coord Z
    id++;
    
    //4
    coord[0] = r*cos(M_PI/4.0);
    coord[1] = r*sin(M_PI/4.0);
    coord[2] = z0;
    RotateNode(coord, angrot, Axis);
    
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,coord[0]);//coord X
    gmesh->NodeVec()[id].SetCoord(1,coord[1]);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,coord[2]);//coord Z
    id++;
    
    //5
    coord[0] = r*cos(M_PI/2.0);
    coord[1] = r*sin(M_PI/2.0);
    coord[2] = z0;
    RotateNode(coord, angrot, Axis);
    
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,coord[0]);//coord X
    gmesh->NodeVec()[id].SetCoord(1,coord[1]);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,coord[2]);//coord Z
    id++;

    //6
    coord[0] = r*cos(3.0*M_PI/4.0);
    coord[1] = r*sin(3.0*M_PI/4.0);
    coord[2] = z0;
    RotateNode(coord, angrot, Axis);
    
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,coord[0]);//coord X
    gmesh->NodeVec()[id].SetCoord(1,coord[1]);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,coord[2]);//coord Z
    id++;
    
    //7
    coord[0] = r*cos(M_PI/4.0);
    coord[1] = r*sin(M_PI/4.0);
    coord[2] = z1;
    RotateNode(coord, angrot, Axis);
    
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,coord[0]);//coord X
    gmesh->NodeVec()[id].SetCoord(1,coord[1]);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,coord[2]);//coord Z
    id++;

    //8
    coord[0] = r*cos(M_PI/2.0);
    coord[1] = r*sin(M_PI/2.0);
    coord[2] = z1;
    RotateNode(coord, angrot, Axis);
    
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,coord[0]);//coord X
    gmesh->NodeVec()[id].SetCoord(1,coord[1]);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,coord[2]);//coord Z
    id++;
    
    //9
    coord[0] = r*cos(3.0*M_PI/4.0);
    coord[1] = r*sin(3.0*M_PI/4.0);
    coord[2] = z1;
    RotateNode(coord, angrot, Axis);
    
    gmesh->NodeVec()[id].SetNodeId(id);
    gmesh->NodeVec()[id].SetCoord(0,coord[0]);//coord X
    gmesh->NodeVec()[id].SetCoord(1,coord[1]);//coord Y
    gmesh->NodeVec()[id].SetCoord(2,coord[2]);//coord Z
    id++;
    
    int elementid = 0;
    TPZVec < int64_t > nodeindex(3,0L);
    
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
    
//     std::ofstream outfile("malhaCilindrica.vtk");
//     TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outfile, true);
    
    return gmesh;
    
}

void hdivCurvedJCompAppMath::RotateNode(TPZVec<REAL> &iCoords, REAL CounterClockwiseAngle, int &Axis)
{
    REAL theta =  (M_PI/180.0)*CounterClockwiseAngle;
    // It represents a 3D rotation around the z axis.
    TPZFMatrix<REAL> RotationMatrix(3,3,0.0);
    
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
    
    TPZVec<REAL> iCoordsRotated(3,0.0);
    // Apply rotation
    iCoordsRotated[0] = RotationMatrix(0,0)*iCoords[0]+RotationMatrix(0,1)*iCoords[1]+RotationMatrix(0,2)*iCoords[2];
    iCoordsRotated[1] = RotationMatrix(1,0)*iCoords[0]+RotationMatrix(1,1)*iCoords[1]+RotationMatrix(1,2)*iCoords[2];
    iCoordsRotated[2] = RotationMatrix(2,0)*iCoords[0]+RotationMatrix(2,1)*iCoords[1]+RotationMatrix(2,2)*iCoords[2];
    iCoords = iCoordsRotated;
}


void hdivCurvedJCompAppMath::SolExata(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux){
  
  solp.resize(1);
  solp[0]=0.;
  flux.Resize(3,1);
  
  solp[0]=0.;
  flux(0,0)=0.0;
  flux(1,0)=0.0;
  flux(2,0)=0.0;
    
  REAL x = pt[0];
  REAL y = pt[1];
  REAL z = pt[2];
  
  if(probAtCircle)
  {
    //+======================================================================================
    // Circle
    
    REAL r = sqrt(x*x+y*y+z*z);
    REAL theta = atan2(y,x);
    
    REAL costheta = cos(theta);
    REAL sintheta = sin(theta);
    REAL cos2theta = cos(2.0*theta);
    REAL sin2theta = sin(2.0*theta);
    REAL r4 = r*r*r*r;
    
    
    
    solp[0] = r4*cos2theta*sin2theta;
    
    // Gradient computations
    
    REAL Runitx = costheta;
    REAL Runity = sintheta;
    REAL Runitz = 0.0;
    
    REAL Thetaunitx = -sintheta;
    REAL Thetaunity = costheta;
    REAL Thetaunitz = 0.0;
    
    REAL Zunitx = 0.0;
    REAL Zunity = 0.0;
    REAL Zunitz = 1.0;
    
    REAL dfdR           = 4.0*r*r*r*cos2theta*sin2theta;
    REAL dfdTheta       = (1.0)*(2.0*r*r*r*(cos2theta*cos2theta-sin2theta*sin2theta));
    REAL dfdZ           = 0.0;
    
    flux(0,0) = -1.0*(dfdR * Runitx + dfdTheta * Thetaunitx + dfdZ * Zunitx);
    flux(1,0) = -1.0*(dfdR * Runity + dfdTheta * Thetaunity + dfdZ * Zunity);
    flux(2,0) = -1.0*(dfdR * Runitz + dfdTheta * Thetaunitz + dfdZ * Zunitz);
  }
  else if(probAtCylinder)
  {
    //+======================================================================================
    // Cylinder
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
  else if(probAtSphere)
  {
    //+======================================================================================
    // Complete Sphere
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
  }
  else
  {
    std::cout << " Something is wrong " << std::endl;
    // One of the previous choices should be true
    DebugStop();
  }
    
}

void hdivCurvedJCompAppMath::Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &ff){
    
  REAL x = pt[0];
  REAL y = pt[1];
  REAL z = pt[2];
  ff.resize(1);
  
    if(probAtCircle)
  {
    //+======================================================================================
    // Circle
    ff[0] = 0.0;
  }
  else if(probAtCylinder)
  {
    //+======================================================================================
    // Cylinder
    REAL r = sqrt(x*x+y*y);
    REAL theta = atan2(y,x);
    ff[0] = ((1.0 + r*r)*cos(z)*sin(theta))/(r*r);
  }
  else if(probAtSphere)
  {
    //+======================================================================================
    // Complete Sphere
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
  }
  else
  {
    std::cout << " Something is wrong " << std::endl;
    // One of the previous choices should be true
    DebugStop();
  }
    
}

void hdivCurvedJCompAppMath::SolExataH1(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux){
    
  solp.resize(1);
  solp[0]=0.;
  flux.Resize(3, 1);
    
  // not implemented
  DebugStop();
    
  if(probAtCircle)
  {
    
  }
  else if(probAtCylinder)
  {
    
  }
  else if(probAtSphere)
  {
    
  }
  else
  {
    std::cout << " Something is wrong " << std::endl;
    // One of the previous choices should be true
    DebugStop();;
  }

}

void hdivCurvedJCompAppMath::ForcingH1(const TPZVec<REAL> &pt, TPZVec<STATE> &ff, TPZFMatrix<STATE> &flux)
{
  // not implemented
  DebugStop();
  
  
  if(probAtCircle)
  {
    
  }
  else if(probAtCylinder)
  {
    
  }
  else if(probAtSphere)
  {
    
  }
  else
  {
    std::cout << " Something is wrong " << std::endl;
    // One of the previous choices should be true
    DebugStop();
  }

}


void hdivCurvedJCompAppMath::ForcingBC1D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){
  solp.resize(1);
  solp[0]=0.;
    
  REAL x = pt[0];
  REAL y = pt[1];
  REAL z = pt[2];
  
  if(probAtCircle)
  {
    //+======================================================================================
    // Circle
    
    REAL r = sqrt(x*x+y*y+z*z);
    REAL theta = atan2(y,x);
    
    REAL cos2theta = cos(2.0*theta);
    REAL sin2theta = sin(2.0*theta);
    REAL r4 = r*r*r*r;
    
    solp[0] = r4*cos2theta*sin2theta;
  }
  else if(probAtCylinder)
  {
    //+======================================================================================
    // Cylinder

    REAL theta = atan2(y,x);
    
    solp[0] = sin(theta)*sin(z + M_PI/2.0);
  }
  else if(probAtSphere)
  {
    //+======================================================================================
    // Complete Sphere
    // sobre a casca da esfera -- conferir o raio aqui usado com o da malha geometrica

    REAL theta = atan2(sqrt(x*x+y*y),z);
    REAL a = M_PI;
    
    solp[0] = (a-theta)*sin(theta)*sin(theta);
  }
  else
  {
    std::cout << " Something is wrong " << std::endl;
    // One of the previous choices should be true
    DebugStop();
  }  

}

void hdivCurvedJCompAppMath::ForcingBC2D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){

  solp.resize(1);
  solp[0]=0.;
    
  REAL x = pt[0];
  REAL y = pt[1];
  REAL z = pt[2];
  
  if(probAtCircle)
  {
    //+======================================================================================
    // Circle
    
    REAL r = sqrt(x*x+y*y+z*z);
    REAL theta = atan2(y,x);
    
    REAL cos2theta = cos(2.0*theta);
    REAL sin2theta = sin(2.0*theta);
    REAL r4 = r*r*r*r;
    
    solp[0] = r4*cos2theta*sin2theta;
  }
  else if(probAtCylinder)
  {
    //+======================================================================================
    // Cylinder

    REAL theta = atan2(y,x);
    
    solp[0] = sin(theta)*sin(z + M_PI/2.0);
  }
  else if(probAtSphere)
  {
    //+======================================================================================
    // Complete Sphere
    // sobre a casca da esfera -- conferir o raio aqui usado com o da malha geometrica
    
    REAL theta = atan2(sqrt(x*x+y*y),z);
    REAL a = M_PI;
    
    solp[0] = (a-theta)*sin(theta)*sin(theta);
  }
  else
  {
    std::cout << " Something is wrong " << std::endl;
    // One of the previous choices should be true
    DebugStop();
  }  

}


void hdivCurvedJCompAppMath::ForcingBC3D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){
  
  solp.resize(1);
  solp[0]=0.;
    
  REAL x = pt[0];
  REAL y = pt[1];
  REAL z = pt[2];
  
  if(probAtCircle)
  {
    //+======================================================================================
    // Circle
    
    REAL r = sqrt(x*x+y*y+z*z);
    REAL theta = atan2(y,x);
    
    REAL cos2theta = cos(2.0*theta);
    REAL sin2theta = sin(2.0*theta);
    REAL r4 = r*r*r*r;
    
    solp[0] = r4*cos2theta*sin2theta;
  }
  else if(probAtCylinder)
  {
    //+======================================================================================
    // Cylinder

    REAL theta = atan2(y,x);
    
    solp[0] = sin(theta)*sin(z + M_PI/2.0);
  }
  else if(probAtSphere)
  {
    //+======================================================================================
    // Complete Sphere
    // sobre a casca da esfera -- conferir o raio aqui usado com o da malha geometrica

    REAL theta = atan2(sqrt(x*x+y*y),z);
    REAL a = M_PI;
    
    solp[0] = (a-theta)*sin(theta)*sin(theta);
  }
  else
  {
    std::cout << " Something is wrong " << std::endl;
    // One of the previous choices should be true
    DebugStop();
  }  

}


void hdivCurvedJCompAppMath::ForcingBC4D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){
  
  solp.resize(1);
  solp[0]=0.;
    
  REAL x = pt[0];
  REAL y = pt[1];
  REAL z = pt[2];
  
  if(probAtCircle)
  {
    //+======================================================================================
    // Circle
    
    REAL r = sqrt(x*x+y*y+z*z);
    REAL theta = atan2(y,x);
    
    REAL cos2theta = cos(2.0*theta);
    REAL sin2theta = sin(2.0*theta);
    REAL r4 = r*r*r*r;
    
    solp[0] = r4*cos2theta*sin2theta;
  }
  else if(probAtCylinder)
  {
    //+======================================================================================
    // Cylinder
    REAL theta = atan2(y,x);
    
    solp[0] = sin(theta)*sin(z + M_PI/2.0);
  }
  else if(probAtSphere)
  {
    //+======================================================================================
    // Complete Sphere
    
    REAL theta = atan2(sqrt(x*x+y*y),z);
    REAL a = M_PI;
    
    solp[0] = (a-theta)*sin(theta)*sin(theta);
  }
  else
  {
    std::cout << " Something is wrong " << std::endl;
    // One of the previous choices should be true
    DebugStop();
  }  

}
void hdivCurvedJCompAppMath::ForcingBC5D(const TPZVec<REAL> &pt, TPZVec<STATE> &solp){

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
    
}



TPZCompMesh *hdivCurvedJCompAppMath::CMeshH1(TPZGeoMesh *gmesh, int pOrder, int dim)
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

    TPZMaterial * BCond1;
    TPZMaterial * BCond2;
    TPZMaterial * BCond3;
    TPZMaterial * BCond4;

    
    BCond1 = material->CreateBC(mat, fbc1,fdirichlet, val1, val2);
    BCond2 = material->CreateBC(mat, fbc2,fdirichlet, val1, val2);
    BCond3 = material->CreateBC(mat, fbc3,fdirichlet, val1, val2);
    BCond4 = material->CreateBC(mat, fbc4,fdirichlet, val1, val2);
    
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    
    
    cmesh->SetAllCreateFunctionsContinuous();
    
    //Ajuste da estrutura de dados computacional
    cmesh->AutoBuild();
    
    return cmesh;
    
}


TPZCompMesh *hdivCurvedJCompAppMath::CMeshFlux(TPZGeoMesh *gmesh, int pOrder, int dim)
{
  if(probAtCircle)
  {
    /// criar materiais
    TPZMatMixedPoisson3D *material = new TPZMatMixedPoisson3D(fmatId,fDim);
    
    TPZMaterial * mat(material);
    material->NStateVariables();
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    
    ///Inserir condicao de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);

    TPZMaterial * BCond1;
    TPZMaterial * BCond2;
    TPZMaterial * BCond3;
    TPZMaterial * BCond4;
    BCond1 = material->CreateBC(mat, fbc1,fdirichlet, val1, val2);
    BCond2 = material->CreateBC(mat, fbc2,fdirichlet, val1, val2);
    BCond3 = material->CreateBC(mat, fbc3,fdirichlet, val1, val2);
    BCond4 = material->CreateBC(mat, fbc4,fdirichlet, val1, val2);
    
    cmesh->InsertMaterialObject(mat);
    
    cmesh->SetDimModel(fDim);
    
    cmesh->SetAllCreateFunctionsHDiv();
    
    
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    
    cmesh->SetDefaultOrder(pOrder);
    
    TPZLagrangeMultiplier *matskelet = new TPZLagrangeMultiplier(fmatskeleton, fDim-1, 1);
    TPZMaterial * mat2(matskelet);
    cmesh->InsertMaterialObject(mat2);
    
    //Ajuste da estrutura de dados computacional
    cmesh->AutoBuild();
    
    return cmesh;
  }
  else if(probAtCylinder)
  {
    /// criar materiais
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
    
    //Ajuste da estrutura de dados computacional
    cmesh->AutoBuild();
    
    return cmesh;
  }
  else if(probAtSphere)
  {
    /// criar materiais
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
    
    this->SetupDisconnectedHdivboud(fbc0,fbc1,cmesh);

    return cmesh;
  }
  else
  {
    std::cout << " Something is wrong " << std::endl;
    // One of the previous choices should be true
    DebugStop();
  }
  return new TPZCompMesh(gmesh);
}

void hdivCurvedJCompAppMath::SetupDisconnectedHdivboud(const int left,const int rigth, TPZCompMesh * cmesh)
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

TPZCompMesh *hdivCurvedJCompAppMath::CMeshPressure(TPZGeoMesh *gmesh, int pOrder, int dim)
{
  
  if(probAtCircle)
  {
   
    /// criar materiais
    TPZMatPoisson3d *material = new TPZMatPoisson3d(fmatId,fDim);
    material->NStateVariables();
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(fDim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    
    ///Inserir condicao de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    cmesh->SetDefaultOrder(pOrder);

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
            
            if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
            {
                if(ftriang==true) celdisc->SetTotalOrderShape();
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
    
//    int nel = cmesh->NElements();
//    for(int i=0; i<nel; i++){
//        TPZCompEl *cel = cmesh->ElementVec()[i];
//
//        
//    }
    return cmesh;
  }
  else if(probAtCylinder)
  {
    /// criar materiais
    TPZMatMixedPoisson3D *material = new TPZMatMixedPoisson3D(fmatId,fDim);
    material->NStateVariables();
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    
    ///Inserir condicao de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    cmesh->SetDefaultOrder(pOrder);
    
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
    
    return cmesh;
  }
  else if(probAtSphere)
  {
    TPZMatMixedPoisson3D *material = new TPZMatMixedPoisson3D(fmatId,fDim);
    
    material->NStateVariables();

    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    
    ///Inserir condicao de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
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
           
            
            if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
            {
                if(ftriang==true) celdisc->SetTotalOrderShape();
                else celdisc->SetTensorialShape();
            }
            
        }
    }
    
    //#endif
    return cmesh;
  }
  else
  {
    std::cout << " Something is wrong " << std::endl;
    // One of the previous choices should be true
    DebugStop();
  }
    return new TPZCompMesh;
}

TPZCompMesh *hdivCurvedJCompAppMath::CMeshMixed(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec)
{
    
    if(probAtCircle)
  {
   
    //Creating computational mesh for multiphysic elements
    gmesh->ResetReference();
    TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    //criando material
    int dim = gmesh->Dimension();
    bool intface;
//    TPZMatPoissonD3 *material = new TPZMatPoissonD3(fmatId,dim); intface = true; // nesse material tem que ser true
    TPZMatMixedPoisson3D *material = new TPZMatMixedPoisson3D(fmatId,dim); intface = true; // nesse material tem que ser true
    //TPZMixedPoisson *material = new TPZMixedPoisson(matId,dim); intface = false; // nesse material tem que ser false
    
    //incluindo os dados do problema
    
    //incluindo os dados do problema
    TPZFNMatrix<2,REAL> PermTensor(dim,dim,0.);
    TPZFNMatrix<2,REAL> InvPermTensor(dim,dim,0.);
    
    // tensor de permutacao
    TPZFNMatrix<9,REAL> TP(dim,dim,0.0);
    TPZFNMatrix<9,REAL> InvTP(dim,dim,0.0);
    
    // Hard coded
    for (int id = 0; id < dim; id++){
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
    dum->SetPolynomialOrder(20);
    forcef = dum;
    material->SetForcingFunction(forcef);
    
    //inserindo o material na malha computacional
    TPZMaterial *mat(material);
    mphysics->InsertMaterialObject(mat);
    
    
    //Criando condicoes de contorno
    TPZMaterial * BCond0 = NULL;
    TPZMaterial * BCond1;
    TPZMaterial * BCond2 = NULL;
    TPZMaterial * BCond3 = NULL;
    TPZMaterial * BCond4 = NULL;
    TPZMaterial * BCond5 = NULL;
    
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);

    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZAutoPointer<TPZFunction<STATE> > FBCond1 = new TPZDummyFunction<STATE>(ForcingBC1D, 5);
    BCond1 = material->CreateBC(mat, fbc1,fdirichlet, val1, val2);
    BCond1->SetForcingFunction(FBCond1);
    //    BCond1 = material->CreateBC(mat, bc1,neumann, val1, val2);
    
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
    
    // Creating multiphysic elements into mphysics computational mesh
    TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
    TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
    
    mphysics->Reference()->ResetReference();
    mphysics->LoadReferences();
    
    // Creation of interface elements
    if (intface)
    {
        int nel = mphysics->ElementVec().NElements();
        for(int el = 0; el < nel; el++)
        {
            TPZCompEl * compEl = mphysics->ElementVec()[el];
            if(!compEl) continue;
            int index = compEl ->Index();
            if(compEl->Dimension() == mphysics->Dimension())
            {
                TPZMultiphysicsElement * InterpEl = dynamic_cast<TPZMultiphysicsElement *>(mphysics->ElementVec()[index]);
                if(!InterpEl) continue;
                InterpEl->CreateInterfaces();
            }
        }
        
    }
    
    return mphysics;
  }
  else if(probAtCylinder)
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
    TPZMaterial * BCond0 = NULL;
    TPZMaterial * BCond1;
    TPZMaterial * BCond2;
    TPZMaterial * BCond3;
    TPZMaterial * BCond4;
    TPZMaterial * BCond5 = NULL;
    
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
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
  else if(probAtSphere)
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
        }
        
        mphysics->CleanUpUnconnectedNodes();
        mphysics->ExpandSolution();
    }
    else
    {
        TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
        mphysics->Reference()->ResetReference();
        mphysics->LoadReferences();
        
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
        
    }

    
    
    mphysics->ComputeNodElCon();
    mphysics->CleanUpUnconnectedNodes();
    mphysics->ExpandSolution();
    
    return mphysics;
  }
  else
  {
    std::cout << " Something is wrong " << std::endl;
    // One of the previous choices should be true
    DebugStop();
  }
    return new TPZCompMesh;
}


void hdivCurvedJCompAppMath::ErrorH1(TPZCompMesh *l2mesh, int p, int ndiv, int pos,  TPZFMatrix< REAL > &errors)
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

void hdivCurvedJCompAppMath::ErrorH1(TPZCompMesh *l2mesh, int p, int ndiv, std::ostream &out, int DoFT, int DofCond)
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


void hdivCurvedJCompAppMath::ErrorPrimalDual(TPZCompMesh *l2mesh, TPZCompMesh *hdivmesh,  int p, int ndiv, int pos, TPZFMatrix< REAL > &errors)
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

void hdivCurvedJCompAppMath::ErrorPrimalDual(TPZCompMesh *l2mesh, TPZCompMesh *hdivmesh,  int p, int ndiv, std::ostream &out, int DoFT, int DofCond)
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

void hdivCurvedJCompAppMath::ChangeExternalOrderConnects(TPZCompMesh *mesh){
    
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



void hdivCurvedJCompAppMath::SolveSyst(TPZAnalysis &an, TPZCompMesh *fCmesh)
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






