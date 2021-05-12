/**
 * @file
 * @brief Contains the implementation of the TPZBSpStructMatrix methods. 
 */

#include "TPZBSpStructMatrix.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzelmat.h"
#include "pzysmp.h"
#include "TPZRenumbering.h"
#include "TPZGuiInterface.h"

using namespace std;

template<class TVar, class TPar>
TPZStructMatrix * TPZBSpStructMatrix<TVar,TPar>::Clone(){
    return new TPZBSpStructMatrix(*this);
}

template<class TVar,class TPar>
TPZMatrix<TVar> * TPZBSpStructMatrix<TVar,TPar>::Create(){
    //checked
    if(this->fMesh->FatherMesh()) {
		PZError << "TPZSpStructMatrix should not be called with CreateAssemble for a substructure mesh\n";
        DebugStop();
    }
    int64_t neq = this->fEquationFilter.NActiveEquations();
    TPZFYsmpMatrix<TVar> * mat = new TPZFYsmpMatrix<TVar>(neq,neq);
	
    /**Rearange elements order*/
	//    TPZVec<int> elorder(fMesh->NEquations(),0);
	
	
    /**
     *Longhin implementation
	 */
    TPZStack<int64_t> elgraph;
    TPZVec<int64_t> elgraphindex;
	//    int nnodes = 0;
    this->fMesh->ComputeElGraph(elgraph,elgraphindex);
    /**Creates a element graph*/
    TPZRenumbering metis;
    metis.SetElementsNodes(elgraphindex.NElements() -1 ,this->fMesh->NIndependentConnects());
    metis.SetElementGraph(elgraph,elgraphindex);
	
    TPZManVector<int64_t> nodegraph;
    TPZManVector<int64_t> nodegraphindex;
    /**
     *converts an element graph structure into a node graph structure
     *those vectors have size ZERO !!!
     */
    metis.ConvertGraph(elgraph,elgraphindex,nodegraph,nodegraphindex);
    /**vector sizes*/
    int64_t i;
    int64_t nblock = nodegraphindex.NElements()-1;
    int64_t totalvar = 0;
    int64_t totaleq = 0;
    for(i=0;i<nblock;i++){
		int64_t iblsize = this->fMesh->Block().Size(i);
		int64_t iblpos = this->fMesh->Block().Position(i);
        int64_t numactive = this->fEquationFilter.NumActive(iblpos, iblpos+iblsize);
        if (!numactive) {
            continue;
        }
        if (numactive != iblsize) {
            DebugStop();
        }
		totaleq += iblsize;
		int64_t icfirst = nodegraphindex[i];
		int64_t iclast = nodegraphindex[i+1];
		int64_t j;
		//longhin
		totalvar+=iblsize*iblsize;
		for(j=icfirst;j<iclast;j++) {
			int64_t col = nodegraph[j];
			int64_t colsize = this->fMesh->Block().Size(col);
			int64_t colpos = this->fMesh->Block().Position(col);
            int64_t numactive = this->fEquationFilter.NumActive(colpos, colpos+colsize);
            if (!numactive) {
                continue;
            }
			totalvar += iblsize*colsize;
		}
    }
	
    int64_t ieq = 0;
    int64_t pos = 0;
	
    nblock=this->fMesh->NIndependentConnects();
	
    int64_t * Eq = new int64_t[totaleq+1];
    int64_t * EqCol = new int64_t[totalvar/2];
    TVar * EqValue = new TVar[totalvar/2];
    for(i=0;i<nblock;i++){
		int64_t iblsize = this->fMesh->Block().Size(i);
		int64_t iblpos = this->fMesh->Block().Position(i);
        int64_t numactive = this->fEquationFilter.NumActive(iblpos, iblpos+iblsize);
        if (!numactive) {
            continue;
        }
		int64_t ibleq;
		for(ibleq=0; ibleq<iblsize; ibleq++) {
			Eq[ieq] = pos;
			if(ieq%2) {
				ieq++;
				continue;
			}
			int64_t colsize = this->fMesh->Block().Size(i);
			int64_t colpos = this->fMesh->Block().Position(i);
			int64_t jbleq;
			for(jbleq=0; jbleq<colsize; jbleq++) {
				/**It can also be implemented using half the size of both columns and data vectors*/
				EqCol[pos] = -1;//colpos;
				EqValue[pos] = 0.;
				colpos++;
				pos++;
			}
			
			int64_t icfirst = nodegraphindex[i];
			int64_t iclast = nodegraphindex[i+1];
			int64_t j;
			for(j=icfirst;j<iclast;j++) {
				int64_t col = nodegraph[j];
				colsize = this->fMesh->Block().Size(col);
				colpos = this->fMesh->Block().Position(col);
                int64_t numactive = this->fEquationFilter.NumActive(colpos, colpos+colsize);
                if (!numactive) {
                    continue;
                }
                for(jbleq=0; jbleq<colsize; jbleq++) {
					EqCol[pos] = -1;//colpos;
					EqValue[pos] = 0.;
					colpos++;
					pos++;
				}
			}
			ieq++;
		}
    }
    Eq[ieq] = pos;
	/*    for(i=0;i<totalvar;i++){
	 if(i<totaleq+1){
	 cout << i <<  " " << Eq[i] << " "<< EqCol[i] << " " << EqValue[i] << endl;
	 }else{
	 cout << i <<  " " << " "<< EqCol[i] << " " << EqValue[i] << endl;
	 }
	 }
	 */
    mat->SetData(Eq,EqCol,EqValue);
    return mat;
}

template<class TVar, class TPar>
int TPZBSpStructMatrix<TVar,TPar>::ClassId() const{
    return Hash("TPZBSpStructMatrix") ^
        TPZBSpStructMatrix<TVar,TPar>::ClassId() << 1 ^
        TPar::ClassId() << 2;
}

#include "pzstrmatrixot.h"
#include "pzstrmatrixflowtbb.h"


template class TPZBSpStructMatrix<STATE,TPZStructMatrixOR<STATE>>;
template class TPZBSpStructMatrix<STATE,TPZStructMatrixOT<STATE>>;
template class TPZBSpStructMatrix<STATE,TPZStructMatrixTBBFlow<STATE>>;

template class TPZBSpStructMatrix<CSTATE,TPZStructMatrixOR<CSTATE>>;
template class TPZBSpStructMatrix<CSTATE,TPZStructMatrixOT<CSTATE>>;
template class TPZBSpStructMatrix<CSTATE,TPZStructMatrixTBBFlow<CSTATE>>;