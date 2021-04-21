/**
 * @file
 * @brief Contains the implementation of the TPZSSpStructMatrix methods. 
 */

#include "TPZSSpStructMatrix.h"
#include "pzcmesh.h"
#include "pzsysmp.h"
#include "TPZRenumbering.h"
#include "TPZGuiInterface.h"
#include "TPZTimer.h"
#include "pzlog.h"
#ifdef PZ_LOG
static TPZLogger logger("pz.StrMatrix");
#endif

using namespace std;

template<class TVar, class TPar>
TPZStructMatrix * TPZSSpStructMatrix<TVar,TPar>::Clone(){
    return new TPZSSpStructMatrix(*this);
}
template<class TVar, class TPar>
TPZMatrix<TVar> * TPZSSpStructMatrix<TVar,TPar>::CreateAssemble(TPZFMatrix<TVar> &rhs,
                                              TPZAutoPointer<TPZGuiInterface> guiInterface){
	
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        LOGPZ_DEBUG(logger,"TPZSSpStructMatrix::CreateAssemble starting");
    }
#endif
	TPar::InitCreateAssemble();
    int64_t neq = fMesh->NEquations();
    if(fMesh->FatherMesh()) {
		cout << "TPZSSpStructMatrix should not be called with CreateAssemble for a substructure mesh\n";
		return new TPZSYsmpMatrix<TVar>(0,0);
    }
//    std::cout << "Creating\n";
    TPZMatrix<TVar> *stiff = Create();//new TPZFYsmpMatrix(neq,neq);
    TPZSYsmpMatrix<TVar> *mat = dynamic_cast<TPZSYsmpMatrix<TVar> *> (stiff);
    rhs.Redim(neq,1);
    //stiff->Print("Stiffness TPZFYsmpMatrix :: CreateAssemble()");
    TPZTimer before("Assembly of a sparse matrix");
//    std::cout << "Assembling\n";
    before.start();
#ifdef PZ_LOG
    if(logger.isDebugEnabled()) LOGPZ_DEBUG(logger,"TPZSSpStructMatrix::CreateAssemble calling Assemble()");
#endif
	Assemble(*stiff,rhs,guiInterface);
    mat->ComputeDiagonal();
    
//    std::cout << "Rhs norm " << Norm(rhs) << std::endl;
    
    before.stop();
    //std::cout << __PRETTY_FUNCTION__ << " " << before << std::endl;
    //    mat->ComputeDiagonal();
    //stiff->Print("Stiffness TPZFYsmpMatrix :: CreateAssemble()");
#ifdef PZ_LOG
    if(logger.isDebugEnabled()) LOGPZ_DEBUG(logger,"TPZSSpStructMatrix::CreateAssemble exiting");
#endif
    return stiff;
}
template<class TVar, class TPar>
TPZMatrix<TVar> * TPZSSpStructMatrix<TVar,TPar>::Create(){

    /**
     *Longhin implementation
     */
    TPZStack<int64_t> elgraph;
    TPZVec<int64_t> elgraphindex;
    //    int nnodes = 0;
    fMesh->ComputeElGraph(elgraph,elgraphindex);
    TPZMatrix<TVar> * mat = SetupMatrixData(elgraph, elgraphindex);
    return mat;
}

template<class TVar, class TPar>
TPZMatrix<TVar> * TPZSSpStructMatrix<TVar,TPar>::SetupMatrixData(TPZStack<int64_t> & elgraph, TPZVec<int64_t> &elgraphindex){
    
    int64_t neq = fEquationFilter.NActiveEquations();
    TPZSYsmpMatrix<TVar> * mat = new TPZSYsmpMatrix<TVar>(neq,neq);
    
    /**Creates a element graph*/
    TPZRenumbering metis;
    metis.SetElementsNodes(elgraphindex.NElements() -1 ,fMesh->NIndependentConnects());
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
    // number of values in the sparse matrix
    int64_t totalvar = 0;
    // number of equations
    int64_t totaleq = 0;
    for(i=0;i<nblock;i++){
        int64_t iblsize = fMesh->Block().Size(i);
        int64_t iblpos = fMesh->Block().Position(i);
        int64_t numactive = fEquationFilter.NumActive(iblpos, iblpos+iblsize);
        if (!numactive) {
            continue;
        }
        totaleq += iblsize;
        int64_t icfirst = nodegraphindex[i];
        int64_t iclast = nodegraphindex[i+1];
        int64_t j;
        //longhin
        totalvar+=(iblsize*(iblsize+1))/2;
        for(j=icfirst;j<iclast;j++) {
            int64_t col = nodegraph[j];
            if (col < i) {
                continue;
            }
            
            if (col == i) {
                DebugStop();
            }
            
            int64_t colsize = fMesh->Block().Size(col);
            int64_t colpos = fMesh->Block().Position(col);
            int64_t numactive = fEquationFilter.NumActive(colpos, colpos+colsize);
            if (!numactive) {
                continue;
            }
            totalvar += iblsize*colsize;
        }
    }
    
    int64_t ieq = 0;
    // pos is the position where we will put the column value
    int64_t pos = 0;
    
    nblock=fMesh->NIndependentConnects();
    
    TPZVec<int64_t> Eq(totaleq+1);
    TPZVec<int64_t> EqCol(totalvar);
    TPZVec<TVar> EqValue(totalvar,0.);
    for(i=0;i<nblock;i++){
        int64_t iblsize = fMesh->Block().Size(i);
        int64_t iblpos = fMesh->Block().Position(i);
        TPZManVector<int64_t> rowdestindices(iblsize);
        for (int64_t i=0; i<iblsize; i++) {
            rowdestindices[i] = iblpos+i;
        }
        fEquationFilter.Filter(rowdestindices);

        int64_t ibleq;
        // working equation by equation
        for(ibleq=0; ibleq<rowdestindices.size(); ibleq++) {
            if (rowdestindices[ibleq] != ieq) {
                DebugStop();
            }
            Eq[ieq] = pos;
            int64_t colsize,colpos,jbleq;
            int64_t diagonalinsert = 0;
            int64_t icfirst = nodegraphindex[i];
            int64_t iclast = nodegraphindex[i+1];
            int64_t j;
            for(j=icfirst;j<iclast;j++)
            {
                int64_t col = nodegraph[j];
                if (col < i) {
                    continue;
                }
                // force the diagonal block to be inserted
                // the nodegraph does not contain the pointer to itself
                if(!diagonalinsert && col > i)
                {
                    diagonalinsert = 1;
                    int64_t colsize = fMesh->Block().Size(i);
                    int64_t colpos = fMesh->Block().Position(i);
                    TPZManVector<int64_t> destindices(colsize);
                    for (int64_t i=0; i<colsize; i++) {
                        destindices[i] = colpos+i;
                    }
                    fEquationFilter.Filter(destindices);
                    int64_t jbleq;
                    for(jbleq=0; jbleq<destindices.size(); jbleq++) {
                        //             if(colpos+jbleq == ieq) continue;
                        int64_t jeq = destindices[jbleq];
                        if (jeq < ieq) {
                            continue;
                        }
                        EqCol[pos] = destindices[jbleq];
                        EqValue[pos] = 0.;
                        //            colpos++;
                        pos++;
                    }
                }
                colsize = fMesh->Block().Size(col);
                colpos = fMesh->Block().Position(col);
                if (fEquationFilter.NumActive(colpos, colpos+colsize) == 0) {
                    continue;
                }
                TPZManVector<int64_t> destindices(colsize);
                for (int64_t i=0; i<colsize; i++) {
                    destindices[i] = colpos+i;
                }
                fEquationFilter.Filter(destindices);
                for(jbleq=0; jbleq<destindices.size(); jbleq++) {
                    int64_t jeq = destindices[jbleq];
                    if (jeq < ieq) {
                        continue;
                    }
                    EqCol[pos] = jeq;
                    EqValue[pos] = 0.;
                    colpos++;
                    pos++;
                }
            }
            // all elements are below (last block certainly)
            if(!diagonalinsert)
            {
                diagonalinsert = 1;
                int64_t colsize = fMesh->Block().Size(i);
                int64_t colpos = fMesh->Block().Position(i);
                TPZManVector<int64_t> destindices(colsize);
                for (int64_t i=0; i<colsize; i++) {
                    destindices[i] = colpos+i;
                }
                fEquationFilter.Filter(destindices);
                int64_t jbleq;
                for(jbleq=0; jbleq<destindices.size(); jbleq++) {
                    //             if(colpos+jbleq == ieq) continue;
                    int64_t jeq = destindices[jbleq];
                    if (jeq < ieq) {
                        continue;
                    }
                    EqCol[pos] = jeq;
                    EqValue[pos] = 0.;
                    //            colpos++;
                    pos++;
                }
            }
            ieq++;
        }
    }

    Eq[ieq] = pos;
    mat->SetData(Eq,EqCol,EqValue);
    return mat;
}
template<class TVar, class TPar>
TPZSSpStructMatrix<TVar,TPar>::TPZSSpStructMatrix(TPZCompMesh *mesh) :
    TPZStructMatrix(mesh)
{}

template<class TVar, class TPar>
TPZSSpStructMatrix<TVar,TPar>::TPZSSpStructMatrix(TPZAutoPointer<TPZCompMesh> mesh) :
    TPZStructMatrix(mesh)
{}

template<class TVar, class TPar>
int TPZSSpStructMatrix<TVar,TPar>::ClassId() const{
    return Hash("TPZSSpStructMatrix") ^
        TPZStructMatrix::ClassId() << 1 ^
        TPar::ClassId() << 2;
}

template<class TVar, class TPar>
void TPZSSpStructMatrix<TVar,TPar>::Read(TPZStream& buf, void* context){
    TPZStructMatrix::Read(buf,context);
    TPar::Read(buf,context);
}

template<class TVar, class TPar>
void TPZSSpStructMatrix<TVar,TPar>::Write(TPZStream& buf, int withclassid) const{
    TPZStructMatrix::Write(buf,withclassid);
    TPar::Write(buf,withclassid);
}

#include "pzstrmatrixot.h"
#include "pzstrmatrixflowtbb.h"

template class TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>>;
template class TPZSSpStructMatrix<STATE,TPZStructMatrixOT<STATE>>;
template class TPZSSpStructMatrix<STATE,TPZStructMatrixTBBFlow<STATE>>;
