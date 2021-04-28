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
void TPZSSpStructMatrix<TVar,TPar>::EndCreateAssemble(TPZBaseMatrix *mat){
    auto spMat = dynamic_cast<TPZSYsmpMatrix<TVar> *>(mat);
    spMat->ComputeDiagonal();
}

template<class TVar, class TPar>
TPZMatrix<TVar> * TPZSSpStructMatrix<TVar,TPar>::Create(){
    int64_t neq = this->fMesh->NEquations();
    if(this->fMesh->FatherMesh()) {
		cout << "TPZSSpStructMatrix should not be called with CreateAssemble for a substructure mesh\n";
        DebugStop();
    }
    /**
     *Longhin implementation
     */
    TPZStack<int64_t> elgraph;
    TPZVec<int64_t> elgraphindex;
    //    int nnodes = 0;
    this->fMesh->ComputeElGraph(elgraph,elgraphindex);
    
    TPZMatrix<TVar> * mat = SetupMatrixData(elgraph, elgraphindex);
    return mat;
}

template<class TVar, class TPar>
TPZMatrix<TVar> * TPZSSpStructMatrix<TVar,TPar>::SetupMatrixData(TPZStack<int64_t> & elgraph, TPZVec<int64_t> &elgraphindex){
    
    const int64_t neq = this->fEquationFilter.NActiveEquations();
    TPZSYsmpMatrix<TVar> * mat = new TPZSYsmpMatrix<TVar>(neq,neq);
    
    /**Creates a element graph*/
    TPZRenumbering graph;
    graph.SetElementsNodes(elgraphindex.NElements() -1
                           ,this->fMesh->NIndependentConnects());
    graph.SetElementGraph(elgraph,elgraphindex);
    
    TPZManVector<int64_t> nodegraph;
    TPZManVector<int64_t> nodegraphindex;
    /**
     *converts an element graph structure into a node graph structure
     *those vectors have size ZERO !!!
     */
    graph.ConvertGraph(elgraph,elgraphindex,nodegraph,nodegraphindex);
    /**vector sizes*/
    const int64_t nblock = nodegraphindex.NElements()-1;
    // number of values in the sparse matrix
    int64_t totalvar = 0;
    // number of equations
    int64_t totaleq = 0;
    for(auto i=0;i<nblock;i++){
        const int64_t iblsize = this->fMesh->Block().Size(i);
        const int64_t iblpos = this->fMesh->Block().Position(i);
        const int64_t numactive = this->fEquationFilter.NumActive(iblpos, iblpos+iblsize);
        if (!numactive) {
            continue;
        }
        totaleq += iblsize;
        const int64_t icfirst = nodegraphindex[i];
        const int64_t iclast = nodegraphindex[i+1];
        //longhin
        totalvar+=(iblsize*(iblsize+1))/2;
        for(auto j=icfirst;j<iclast;j++) {
            const int64_t col = nodegraph[j];
            if (col < i) {
                continue;
            }
            
            if (col == i) {
                DebugStop();
            }
            
            const int64_t colsize = this->fMesh->Block().Size(col);
            const int64_t colpos = this->fMesh->Block().Position(col);
            const int64_t numactive = this->fEquationFilter.NumActive(colpos, colpos+colsize);
            if (!numactive) {
                continue;
            }
            totalvar += iblsize*colsize;
        }
    }
    
    int64_t ieq = 0;
    // pos is the position where we will put the column value
    int64_t pos = 0;
    
    TPZVec<int64_t> Eq(totaleq+1);
    TPZVec<int64_t> EqCol(totalvar);
    TPZVec<TVar> EqValue(totalvar,0.);
    //lambda for avoid repeating code
    //lambda for avoid repeating code
    auto AddColEqs =
        [this,&EqCol,&EqValue,&pos](const int colsize, const int colpos, const int ieq)
        {
            TPZManVector<int64_t> destindices(colsize);
            for (int64_t i=0; i<colsize; i++) {
                destindices[i] = colpos+i;
            }
            this->fEquationFilter.Filter(destindices);
            for(auto jbleq=0; jbleq<destindices.size(); jbleq++) {
                const int64_t jeq = destindices[jbleq];
                if (jeq < ieq) {
                    continue;
                }
                EqCol[pos] = destindices[jbleq];
                EqValue[pos] = 0.;
                pos++;
            }
        };
    
    for(auto i=0;i<nblock;i++){
        const int64_t iblsize = this->fMesh->Block().Size(i);
        const int64_t iblpos = this->fMesh->Block().Position(i);
        const int64_t numactive =
            this->fEquationFilter.NumActive(iblpos, iblpos+iblsize);
        if (!numactive) {
            continue;
        }
        TPZManVector<int64_t> rowdestindices(iblsize);
        for (int64_t i=0; i<iblsize; i++) {
            rowdestindices[i] = iblpos+i;
        }
        this->fEquationFilter.Filter(rowdestindices);
        // working equation by equation
        for(auto ibleq=0; ibleq<rowdestindices.size(); ibleq++) {
            if (rowdestindices[ibleq] != ieq) {
                DebugStop();
            }
            Eq[ieq] = pos;
            bool diagonalinsert = false;
            const int64_t icfirst = nodegraphindex[i];
            const int64_t iclast = nodegraphindex[i+1];
            for(auto j=icfirst;j<iclast;j++)
            {
                const int64_t col = nodegraph[j];
                if (col < i) {
                    continue;
                }
                // force the diagonal block to be inserted
                // the nodegraph does not contain the pointer to itself
                if(!diagonalinsert && col > i)
                {
                    diagonalinsert = true;
                    const auto colsize = this->fMesh->Block().Size(i);
                    const auto colpos = this->fMesh->Block().Position(i);
                    AddColEqs(colsize,colpos,ieq);
                }
                const auto colsize = this->fMesh->Block().Size(col);
                const auto colpos = this->fMesh->Block().Position(col);
                if (this->fEquationFilter.NumActive(colpos, colpos+colsize) == 0) {
                    continue;
                }
                AddColEqs(colsize,colpos,ieq);
            }
            // all elements are below (last block certainly)
            if(!diagonalinsert)
            {
                diagonalinsert = true;
                const auto colsize = this->fMesh->Block().Size(i);
                const auto colpos = this->fMesh->Block().Position(i);
                AddColEqs(colsize,colpos,ieq);
            }
            ieq++;
        }
    }

    Eq[ieq] = pos;
    mat->SetData(Eq,EqCol,EqValue);
    return mat;
}

template<class TVar, class TPar>
int TPZSSpStructMatrix<TVar,TPar>::ClassId() const{
    return Hash("TPZSSpStructMatrix") ^
        TPZStructMatrixT<TVar>::ClassId() << 1 ^
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

template class TPZSSpStructMatrix<CSTATE,TPZStructMatrixOR<CSTATE>>;
template class TPZSSpStructMatrix<CSTATE,TPZStructMatrixOT<CSTATE>>;
template class TPZSSpStructMatrix<CSTATE,TPZStructMatrixTBBFlow<CSTATE>>;
