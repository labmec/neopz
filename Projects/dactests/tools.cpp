//
//  tools.cpp
//  PZ
//
//  Created by Douglas Castro on 1/29/15.
//
//

#include "tools.h"


tools::tools()
{
    DebugStop();
}

tools::~tools()
{
    
}



void tools::SolveSyst(TPZAnalysis &an, TPZCompMesh *fCmesh)
{
    cout <<"Numero de equacoes "<< fCmesh->NEquations()<< endl;
    
	bool isdirect = true;
    bool simetrico = true;
    if (isdirect)
    {
        if (simetrico)
        {
            //TPZBandStructMatrix full(fCmesh);
            TPZSkylineStructMatrix skylstr(fCmesh); //caso simetrico
            //    TPZSkylineNSymStructMatrix full(fCmesh);
            an.SetStructuralMatrix(skylstr);
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

void tools::PosProcess(TPZAnalysis &an, std::string plotfile, int dim){
    TPZManVector<std::string,10> scalnames(1), vecnames(0);
    scalnames[0] = "Solution";
    int div = 3;
    an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    an.PostProcess(div,dim);

}

void tools::PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile, int dim)
{
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
    TPZManVector<std::string,10> scalnames(2), vecnames(2);
    vecnames[0]  = "Flux";
    vecnames[1]  = "ExactFlux";
    scalnames[0] = "Pressure";
    scalnames[1] = "ExactPressure";
    
    int div = 3;
    an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    an.PostProcess(div,dim);
    
    std::ofstream out("malha.txt");
    an.Print("nothing",out);
    
    //    mphysics->Solution().Print("Solucao");
    
}


void tools::RotateGeomesh(TPZGeoMesh *gmesh, REAL CounterClockwiseAngle, int &Axis)
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

void tools::RotateNode(TPZVec<STATE> &iCoords, REAL CounterClockwiseAngle, int &Axis)
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
    
    TPZVec<STATE> iCoordsRotated(3,0.0);
    // Apply rotation
    iCoordsRotated[0] = RotationMatrix(0,0)*iCoords[0]+RotationMatrix(0,1)*iCoords[1]+RotationMatrix(0,2)*iCoords[2];
    iCoordsRotated[1] = RotationMatrix(1,0)*iCoords[0]+RotationMatrix(1,1)*iCoords[1]+RotationMatrix(1,2)*iCoords[2];
    iCoordsRotated[2] = RotationMatrix(2,0)*iCoords[0]+RotationMatrix(2,1)*iCoords[1]+RotationMatrix(2,2)*iCoords[2];
    iCoords = iCoordsRotated;
}

void tools::PrintDebugMapForMathematica(std::string filenameHdiv, std::string filenameL2, std::map<REAL, REAL> &fDebugMapL2, std::map<REAL, REAL> &fDebugMapHdiv)
{
    std::ofstream outHdiv(filenameHdiv.c_str());
    std::ofstream outL2(filenameL2.c_str());
    
    for (std::map<REAL,REAL>::const_iterator it = fDebugMapL2.begin(); it != fDebugMapL2.end(); it++) {
        outL2 << it->first << "   " << it->second << std::endl;
    }
    outL2.close();
    
    
    for (std::map<REAL,REAL>::const_iterator it = fDebugMapHdiv.begin(); it != fDebugMapHdiv.end(); it++) {
        outHdiv <<  it->first << "   " << it->second << std::endl;
    }
    outHdiv.close();
}



//void tools::AddWrap(TPZMultiphysicsElement *mfcel, int matskeleton, TPZStack< TPZStack<TPZMultiphysicsElement *,7> > &ListGroupEl)
//{
//    TPZCompMesh *multiMesh = mfcel->Mesh();
//    TPZInterpolationSpace *hdivel = dynamic_cast<TPZInterpolationSpace *> (mfcel->Element(0));
//    TPZCompElDisc *discel = dynamic_cast<TPZCompElDisc *>(mfcel->Element(1));
//    TPZGeoEl *gel = mfcel->Reference();
//    
//    int dimMesh = mfcel->Mesh()->Dimension();
//    if (!hdivel || !discel || gel->Dimension() != dimMesh) {
//        DebugStop();
//    }
//    
//    //wrap element
//    TPZStack<TPZMultiphysicsElement *, 7> wrapEl;
//    wrapEl.push_back(mfcel);
//    
//    for (int side = 0; side < gel->NSides(); side++)
//    {
//        if (gel->SideDimension(side) != gel->Dimension()-1) {
//            continue;
//        }
//        TPZGeoEl *gelbound = gel->CreateBCGeoEl(side, matskeleton);
//        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(hdivel);
//        int loccon = intel->SideConnectLocId(0,side);
//        long index;
//        
//        TPZInterpolationSpace *bound;
//        MElementType elType = gel->Type(side);
//        switch(elType)
//        {
//            case(EOned)://line
//            {
//                bound = new TPZCompElHDivBound2<pzshape::TPZShapeLinear>(* intel->Mesh(),gelbound,index);
//                int sideorient = intel->GetSideOrient(side);
//                TPZCompElHDivBound2<pzshape::TPZShapeLinear> *hdivbound = dynamic_cast< TPZCompElHDivBound2<pzshape::TPZShapeLinear> *>(bound);
//                hdivbound->SetSideOrient(sideorient);
//                break;
//            }
//            case(ETriangle)://triangle
//            {
//                bound = new TPZCompElHDivBound2<pzshape::TPZShapeTriang>(* intel->Mesh(),gelbound,index);
//                int sideorient = intel->GetSideOrient(side);
//                TPZCompElHDivBound2<pzshape::TPZShapeTriang> *hdivbound = dynamic_cast< TPZCompElHDivBound2<pzshape::TPZShapeTriang> *>(bound);
//                hdivbound->SetSideOrient(sideorient);
//                break;
//            }
//            case(EQuadrilateral)://quadrilateral
//            {
//                bound = new TPZCompElHDivBound2<pzshape::TPZShapeQuad>(* intel->Mesh(),gelbound,index);
//                int sideorient = intel->GetSideOrient(side);
//                TPZCompElHDivBound2<pzshape::TPZShapeQuad> *hdivbound = dynamic_cast< TPZCompElHDivBound2<pzshape::TPZShapeQuad> *>(bound);
//                hdivbound->SetSideOrient(sideorient);
//                break;
//            }
//                
//            default:
//            {
//                bound=0;
//                std::cout << "ElementType not found!";
//                DebugStop();
//                break;
//            }
//        }
//        
//        long sideconnectindex = intel->ConnectIndex(loccon);
//        bound->SetConnectIndex(0, sideconnectindex);
//        //bound->Print(std::cout);
//        
//        TPZCompEl *newMFBound = multiMesh->CreateCompEl(gelbound, index);
//        TPZMultiphysicsElement *locMF = dynamic_cast<TPZMultiphysicsElement *>(newMFBound);
//        
//        locMF->AddElement(bound, 0);
//        locMF->AddElement(TPZCompElSide(discel,side), 1);
//        
//        wrapEl.push_back(locMF);
//    }
//    
//    ListGroupEl.push_back(wrapEl);
//}
