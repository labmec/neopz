//
//  TPZDPGMeshControl.cpp
//  PZ
//
//  Created by Agnaldo Farias on 26/08/14.
//
//

#include "TPZDPGMeshControl.h"
#include "pzmultiphysicscompel.h"
#include "pzbuildmultiphysicsmesh.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.dpgmeshcontrol"));
#endif


TPZDPGMeshControl::TPZDPGMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, std::set<int64_t> &coarseindices) : fGMesh(gmesh),
       fPressureCoarseMesh(gmesh), fMHMControl(gmesh,coarseindices), fPOrderCoarseInternal(-1), fFinerMatId(0),fCoarseMatId(0), fSkeletMatId(0)
{
    fPressureCoarseMesh.SetDimModel(gmesh->Dimension());
    
}

TPZDPGMeshControl::TPZDPGMeshControl(const TPZDPGMeshControl &copy){
    
      this->operator=(copy);
}

TPZDPGMeshControl & TPZDPGMeshControl::operator=(const TPZDPGMeshControl &cp){
    
    fMHMControl = cp.fMHMControl;
    fGMesh = fMHMControl.GMesh();
    fPressureCoarseMesh = cp.fPressureCoarseMesh;
    fPressureCoarseMesh.SetReference(fGMesh);
    fPOrderCoarseInternal = cp.fPOrderCoarseInternal;
    fFinerMatId = cp.fFinerMatId;
    fCoarseMatId = cp.fCoarseMatId;
    fSkeletMatId = cp.fSkeletMatId;
    return *this;
}

void TPZDPGMeshControl::BuildComputationalMesh()
{
    fGMesh->ResetReference();
    //fPressureCoarseMesh.LoadReferences();
    fPressureCoarseMesh.SetAllCreateFunctionsContinuous();
    fPressureCoarseMesh.AutoBuild();
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        sout<< "\n MESH COARSE \n";
        fPressureCoarseMesh.Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    fMHMControl.SetLagrangeAveragePressure(true);
    fMHMControl.CreateCoarseInterfaces(fSkeletMatId);
    fMHMControl.BuildComputationalMesh(false);
    
    //adding elements from coarse mesh pressure to mhm mesh
    int coarsemeshindex = 3;
    fMHMControl.GMesh()->ResetReference();
    fMHMControl.CMesh()->LoadReferences();
    int64_t nelc = fPressureCoarseMesh.NElements();
    for (int64_t el=0; el<nelc; el++) {
        TPZCompEl *cel = fPressureCoarseMesh.Element(el);
        TPZGeoEl *gel = cel->Reference();
        int nsides = gel->NSides();
        TPZGeoElSide gelside(gel,nsides-1);
        TPZStack<TPZCompElSide> higherorequal;
        gelside.HigherLevelCompElementList3(higherorequal, 1, 0);
        if (higherorequal.size() < 1) {
            gelside.EqualLevelCompElementList3(higherorequal, 1, 0);
            if (higherorequal.size() < 1) {
                TPZCompElSide celside = gelside.Reference();
                if (celside.Element()) {
                    higherorequal.Push(celside);
                }
                else
                {
                    DebugStop();
                }
            }
        }
        int64_t nhigher = higherorequal.size();
        for (int64_t ih=0; ih<nhigher; ih++) {
            TPZCompElSide higher = higherorequal[ih];
            TPZCompEl *chigh = higher.Element();
            TPZMultiphysicsElement *cmult = dynamic_cast<TPZMultiphysicsElement *>(chigh);
            if (!cmult) {
                DebugStop();
            }
            if(cmult->Dimension()==cel->Dimension()){
                cmult->AddElement(cel, coarsemeshindex);
            }
        }
    }
    
    //Add the mesh 4 to the multi-physics element with 3 meshes
    int nel = fMHMControl.CMesh()->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fMHMControl.CMesh()->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZMultiphysicsElement *multel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if (!multel) {
            continue;
        }
        int nmeshes = multel->NMeshes();
        for (int im=nmeshes; im<coarsemeshindex+1; im++) {
            multel->AddElement(0, im);
        }
    }

    
    //adding connects from coarse mesh pressure to mhm mesh
    TPZCompMesh * cmesh = &fPressureCoarseMesh;
    TPZBuildMultiphysicsMesh::AppendConnects(cmesh,fMHMControl.CMesh().operator->());
    //AppendConnect
    
    
    fMHMControl.CMesh()->ExpandSolution();
    fMHMControl.CMesh()->SaddlePermute();
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        sout<< "\n\nMULTIPHYSICS MESH DPG METHOD \n";
        fMHMControl.CMesh()->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif

}