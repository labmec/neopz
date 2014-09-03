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


TPZDPGMeshControl::TPZDPGMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, std::set<long> &coarseindices) : fGMesh(gmesh),
       fPressureCoarseMesh(gmesh), fMHMControl(gmesh,coarseindices), fPOrderCoarseInternal(-1), fFinerMatId(0),fCoarseMatId(0), fSkeletMatId(0)
{
 
    
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
    fPressureCoarseMesh.LoadReferences();
    fPressureCoarseMesh.SetAllCreateFunctionsContinuous();
    fPressureCoarseMesh.AutoBuild();
    std::cout <<"Aqui ===\n";
    fPressureCoarseMesh.Print();
    
    fMHMControl.SetLagrangeAveragePressure(true);
    fMHMControl.CreateCoarseInterfaces(fSkeletMatId);
       
    fMHMControl.BuildComputationalMesh();
    TPZVec<TPZCompMesh *> meshvec(1,&fPressureCoarseMesh);
    TPZBuildMultiphysicsMesh::AddConnects(meshvec,fMHMControl.CMesh().operator->());
    int coarsemeshindex = 3;
    
    fMHMControl.GMesh()->ResetReference();
    fMHMControl.CMesh()->LoadReferences();
    long nelc = fPressureCoarseMesh.NElements();
    for (long el=0; el<nelc; el++) {
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
        long nhigher = higherorequal.size();
        for (long ih=0; ih<nhigher; ih++) {
            TPZCompElSide higher = higherorequal[ih];
            TPZCompEl *chigh = higher.Element();
            TPZMultiphysicsElement *cmult = dynamic_cast<TPZMultiphysicsElement *>(chigh);
            if (!cmult) {
                DebugStop();
            }
            cmult->AddElement(cel, coarsemeshindex);
        }
    }
    
    fMHMControl.CMesh()->ExpandSolution();
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        sout<< "\n\nMULTIPHYSICS MESH DPG METHOD \n";
        fMHMControl.CMesh()->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif

}