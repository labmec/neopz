/**
 * @file
 * @brief Contains the implementation of the TPZSpStructMatrix methods. 
 */

#include "TPZSpStructMatrix.h"
#include "pzcmesh.h"
#include "pzysmp.h"
#include "TPZRenumbering.h"
#include "TPZGuiInterface.h"

#include "TPZTimer.h"

#include "pzlog.h"
#ifdef PZ_LOG
static TPZLogger logger("pz.StrMatrix");
#endif

using namespace std;

template<class TVar, class TPar>
TPZStructMatrix * TPZSpStructMatrix<TVar,TPar>::Clone(){
    return new TPZSpStructMatrix(*this);
}
template<class TVar, class TPar>
TPZMatrix<TVar> * TPZSpStructMatrix<TVar,TPar>::CreateAssemble(TPZBaseMatrix &rhs,
                                              TPZAutoPointer<TPZGuiInterface> guiInterface){

#ifdef PZ_LOG
    if(logger.isDebugEnabled())
    {
        LOGPZ_DEBUG(logger,"TPZSpStructMatrix::CreateAssemble starting")
    }
#endif
	TPar::InitCreateAssemble();
    int64_t neq = this->fMesh->NEquations();
    if(this->fMesh->FatherMesh()) {
		cout << "TPZSpStructMatrix should not be called with CreateAssemble for a substructure mesh\n";
		return new TPZFYsmpMatrix<TVar>(0,0);
    }
    TPZMatrix<TVar> *stiff = Create();//new TPZFYsmpMatrix(neq,neq);
    TPZFYsmpMatrix<TVar> *mat = dynamic_cast<TPZFYsmpMatrix<TVar> *> (stiff);
    rhs.Redim(neq,1);
    //stiff->Print("Stiffness TPZFYsmpMatrix :: CreateAssemble()");
    TPZTimer before("Assembly of a sparse matrix");
    before.start();
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        LOGPZ_DEBUG(logger,"TPZSpStructMatrix::CreateAssemble calling Assemble()");
    }
#endif
    this->Assemble(*stiff,rhs,guiInterface);
    before.stop();
    std::cout << __PRETTY_FUNCTION__ << " " << before << std::endl;
//    mat->ComputeDiagonal();
    //    mat->ComputeDiagonal();
    //stiff->Print("Stiffness TPZFYsmpMatrix :: CreateAssemble()");
#ifdef PZ_LOG
    if(logger.isDebugEnabled()) LOGPZ_DEBUG(logger,"TPZSpStructMatrix::CreateAssemble exiting");
#endif
    return stiff;
}
template<class TVar, class TPar>
TPZMatrix<TVar> * TPZSpStructMatrix<TVar,TPar>::Create(){
    int64_t neq = this->fEquationFilter.NActiveEquations();
	/*    if(fMesh->FatherMesh()) {
	 TPZSubCompMesh *smesh = (TPZSubCompMesh *) fMesh;
	 neq = smesh->NumInternalEquations();
	 }*/
    TPZFYsmpMatrix<TVar> * mat = new TPZFYsmpMatrix<TVar>(neq,neq);
	
    /**
     *Longhin implementation
	 */
    TPZStack<int64_t> elgraph;
    TPZVec<int64_t> elgraphindex;
    //    int nnodes = 0;
    this->fMesh->ComputeElGraph(elgraph,elgraphindex);
    /**Creates a element graph*/
    TPZRenumbering metis;
    metis.SetElementsNodes(elgraphindex.NElements() -1 ,
                           this->fMesh->NIndependentConnects());
    metis.SetElementGraph(elgraph,elgraphindex);
	
    TPZManVector<int64_t> nodegraph;
    TPZManVector<int64_t> nodegraphindex;
    /**
     *converts an element graph structure into a node graph structure
     *those vectors have size ZERO !!!
     */
    metis.ConvertGraph(elgraph,elgraphindex,nodegraph,nodegraphindex);
    
#ifdef PZ_LOG2
    if(logger.isDebugEnabled()){
        std::stringstream sout;
        sout << "Node graph \n";
        metis.TPZRenumbering::Print(nodegraph, nodegraphindex);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    /**vector sizes*/
    int64_t nblock = nodegraphindex.NElements()-1;
    // number of values in the sparse matrix
    int64_t totalvar = 0;
    // number of equations
    int64_t totaleq = 0;
    for(int64_t i=0;i<nblock;i++){
		int64_t iblsize = this->fMesh->Block().Size(i);
		int64_t iblpos = this->fMesh->Block().Position(i);
        int64_t numactive = this->fEquationFilter.NumActive(iblpos, iblpos+iblsize);
        if (!numactive) {
            continue;
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
	
#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "Number of equations " << totaleq << " number of nonzero s " << totalvar;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    int64_t ieq = 0;
    // pos is the position where we will put the column value
    int64_t pos = 0;
	
    nblock=this->fMesh->NIndependentConnects();
	
    TPZManVector<int64_t,400> Eq(totaleq+1);
    TPZVec<int64_t> EqCol(totalvar);
    TPZVec<TVar> EqValue(totalvar);
    for(int64_t i=0;i<nblock;i++){
		int64_t iblsize = this->fMesh->Block().Size(i);
		int64_t iblpos = this->fMesh->Block().Position(i);
        TPZManVector<int64_t> rowdestindices(iblsize);
        for (int64_t ij=0; ij<iblsize; ij++) {
            rowdestindices[ij] = iblpos+ij;
        }
        this->fEquationFilter.Filter(rowdestindices);

		int64_t ibleq;
        // working equation by equation
        // rowdestindices contains the equation number of each element in the block number "i"
		for(ibleq=0; ibleq<rowdestindices.size(); ibleq++) {
            int rowind = rowdestindices[ibleq];
//            if (rowind != pos) {
//                DebugStop();
//            }
			Eq[ieq] = pos;
			int64_t colsize,colpos,jbleq;
			int64_t diagonalinsert = 0;
			int64_t icfirst = nodegraphindex[i];
			int64_t iclast = nodegraphindex[i+1];
			int64_t j;
			for(j=icfirst;j<iclast;j++)
            {
                // col is the block linked to block "i"
				int64_t col = nodegraph[j];
                
                // force the diagonal block to be inserted
                // the nodegraph does not contain the pointer to itself
				if(!diagonalinsert && col > i)
				{
					diagonalinsert = 1;
					int64_t colsize = this->fMesh->Block().Size(i);
					int64_t colpos = this->fMesh->Block().Position(i);
                    TPZManVector<int64_t> destindices(colsize);
                    for (int64_t i=0; i<colsize; i++) {
                        destindices[i] = colpos+i;
                    }
                    this->fEquationFilter.Filter(destindices);
					int64_t jbleq;
					for(jbleq=0; jbleq<destindices.size(); jbleq++) {
						//             if(colpos+jbleq == ieq) continue;
						EqCol[pos] = destindices[jbleq];
						EqValue[pos] = 0.;
						//            colpos++;
                        // pos is the position within EqCol or EqVal where we will assemble
						pos++;
					}
				}
				colsize = this->fMesh->Block().Size(col);
				colpos = this->fMesh->Block().Position(col);
                // optimization statement : if all equations in the range are inactive -> continue
                if (this->fEquationFilter.NumActive(colpos, colpos+colsize) == 0) {
                    continue;
                }
                TPZManVector<int64_t> destindices(colsize);
                for (int64_t i=0; i<colsize; i++) {
                    destindices[i] = colpos+i;
                }
                this->fEquationFilter.Filter(destindices);
				for(jbleq=0; jbleq<destindices.size(); jbleq++) {
					EqCol[pos] = destindices[jbleq];
					EqValue[pos] = 0.;
					colpos++;
					pos++;
				}
			}
            // all elements are below (last block certainly)
			if(!diagonalinsert)
			{
				diagonalinsert = 1;
				int64_t colsize = this->fMesh->Block().Size(i);
				int64_t colpos = this->fMesh->Block().Position(i);
                TPZManVector<int64_t> destindices(colsize);
                for (int64_t i=0; i<colsize; i++) {
                    destindices[i] = colpos+i;
                }
                this->fEquationFilter.Filter(destindices);
				int64_t jbleq;
				for(jbleq=0; jbleq<destindices.size(); jbleq++) {
					//             if(colpos+jbleq == ieq) continue;
					EqCol[pos] = destindices[jbleq];
					EqValue[pos] = 0.;
					//            colpos++;
					pos++;
				}
			}
			ieq++;
		}
    }
    Eq[ieq] = pos;
    if(pos != totalvar)
    {
        DebugStop();
    }
    mat->SetData(Eq,EqCol,EqValue);
    return mat;
}

template<class TVar, class TPar>
int TPZSpStructMatrix<TVar,TPar>::ClassId() const{
    return Hash("TPZSpStructMatrix") ^
        TPZStructMatrixT<TVar>::ClassId() << 1 ^
        TPar::ClassId() << 2;
}

template<class TVar, class TPar>
void TPZSpStructMatrix<TVar,TPar>::Read(TPZStream& buf, void* context){
    TPZStructMatrix::Read(buf,context);
    TPar::Read(buf,context);
}

template<class TVar, class TPar>
void TPZSpStructMatrix<TVar,TPar>::Write(TPZStream& buf, int withclassid) const{
    TPZStructMatrix::Write(buf,withclassid);
    TPar::Write(buf,withclassid);
}

#include "pzstrmatrixot.h"
#include "pzstrmatrixflowtbb.h"

template class TPZSpStructMatrix<STATE,TPZStructMatrixOR<STATE>>;
template class TPZSpStructMatrix<STATE,TPZStructMatrixOT<STATE>>;
template class TPZSpStructMatrix<STATE,TPZStructMatrixTBBFlow<STATE>>;
