/**
 * @file
 * @brief Contains the implementation of the TPZElementMatrix methods.
 */

#include "pzelmat.h"
#include "pzbasematrix.h"
#include "pzcmesh.h"

#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.tpzelmat");
#endif

using namespace std;

TPZElementMatrix::TPZElementMatrix(const TPZElementMatrix &cp) : 
    fType(cp.fType), fMesh(cp.fMesh), fConnect(cp.fConnect), 
    fConstrConnect(cp.fConstrConnect),
    fDestinationIndex(cp.fDestinationIndex), fSourceIndex(cp.fSourceIndex)
{

}

TPZElementMatrix &TPZElementMatrix::operator=(const TPZElementMatrix &cp)
{
    fType = cp.fType;
    fMesh = cp.fMesh;
    fConnect = cp.fConnect;    
    fConstrConnect = cp.fConstrConnect;
    fDestinationIndex = cp.fDestinationIndex;
    fSourceIndex = cp.fSourceIndex;
    return *this;
}

void TPZElementMatrix::ComputeDestinationIndices(){
    if (!this->HasDependency()){
        TPZBaseMatrix &mat = Matrix();
		this->fSourceIndex.Resize(mat.Rows());
		this->fDestinationIndex.Resize(mat.Rows());
		int64_t destindex = 0L;
        int64_t fullmatindex = 0L;
		const int numnod = this->NConnects();
		for(int in = 0; in < numnod; in++){
			const int64_t npindex = this->ConnectIndex(in);
			TPZConnect &np = this->fMesh->ConnectVec()[npindex];
			int64_t blocknumber = np.SequenceNumber();
			int64_t firsteq = this->fMesh->Block().Position(blocknumber);
			int ndf = this->fMesh->Block().Size(blocknumber);
            if(np.HasDependency() || np.IsCondensed()) {
                fullmatindex += ndf;
                continue;
            }//for (np)
			for(int idf=0; idf<ndf; idf++){
				this->fSourceIndex[destindex] = fullmatindex++;
				this->fDestinationIndex[destindex++] = firsteq+idf;
			}//for idf
		}//for in
        this->fSourceIndex.Resize(destindex);
        this->fDestinationIndex.Resize(destindex);		
#ifdef PZ_LOG
        if (logger.isDebugEnabled())
		{
			std::stringstream sout;
			sout<<"fSourceIndex " <<fSourceIndex<< "\nfDestinationIndex "<<fDestinationIndex<<std::endl;
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
	}//if
	else{
        TPZBaseMatrix &constrMat = ConstrMatrix();
        int64_t destindex = 0L;
        int64_t fullmatindex = 0L;
        this->fDestinationIndex.Resize(constrMat.Rows());
        this->fSourceIndex.Resize(constrMat.Rows());
        int numnod = this->fConstrConnect.NElements();
        for(int in = 0; in < numnod; in++){
            const int64_t npindex = this->fConstrConnect[in];
            TPZConnect &np = this->fMesh->ConnectVec()[npindex];
            int64_t blocknumber = np.SequenceNumber();
            int64_t firsteq = this->fMesh->Block().Position(blocknumber);
            int ndf = this->fMesh->Block().Size(blocknumber);
            if(np.HasDependency() || np.IsCondensed()) {
                fullmatindex += ndf;
                continue;
            }//for (np)
            for(int idf=0; idf<ndf; idf++) {
                this->fSourceIndex[destindex] = fullmatindex++;
                this->fDestinationIndex[destindex++] = firsteq+idf;
            }//for idf
        }//for in
        this->fSourceIndex.Resize(destindex);
        this->fDestinationIndex.Resize(destindex);
    }
}//void

bool TPZElementMatrix::HasDependency()
{
    if (fOneRestraints.size()) {
        return true;
    }
	int in, nconnects = this->NConnects();
	int64_t index;
	for(in=0; in<nconnects; in++)
    {
		index = this->ConnectIndex(in);
       // bool val =this->fMesh->ConnectVec()[index].HasDependency();
		if(this->fMesh->ConnectVec()[index].HasDependency())
        {
			return true;
		}
	}
	return false;
}

void TPZElementMatrix::BuildDependencyOrder(TPZVec<int64_t> &connectlist, TPZVec<int> &DependenceOrder, TPZCompMesh &mesh) {
    // nodelist (input) : vector which contains pointers to all nodes which
    // are in the dependency chain of the nodes of the element
    int64_t totalnodes = connectlist.NElements();
    DependenceOrder.Resize(totalnodes);
    DependenceOrder.Fill(0,0);
    // initialize the vector which contains the
    // dependency order to zero
    int CurrentOrder = 0;
    // order which is currently processed
    int64_t numnodes_processed = totalnodes;

    for (std::list<TPZOneShapeRestraint>::iterator it = fOneRestraints.begin(); it != fOneRestraints.end(); it++) {
        for (int i=1; i<4; i++)
        {
            int64_t index = it->fFaces[1].first;
            TPZConnect &dfn = mesh.ConnectVec()[index];
            dfn.SetDependenceOrder(index,mesh,1,connectlist,DependenceOrder);
        }
    }

    
    // number of nodes processed during the current cycle
    while(numnodes_processed) {
        
        numnodes_processed = 0;
        int64_t i;
        for(i=0; i<totalnodes; i++) {
            int64_t dfnindex = connectlist[i];
            TPZConnect &dfn = mesh.ConnectVec()[dfnindex];
            if(dfn.HasDependency() && DependenceOrder[i] == CurrentOrder) {
                dfn.SetDependenceOrder(dfnindex,mesh,CurrentOrder,connectlist,DependenceOrder);
                // this method will fill in the DependenceOrder vector by recursively
                // calling SetDependenceOrder for the nodes upon which dfn depends
                numnodes_processed++;
            }
        }
        // force the loop to process the order one connects
        if (fOneRestraints.size() && CurrentOrder == 0) {
            numnodes_processed++;
        }
        CurrentOrder++;
    }
}
