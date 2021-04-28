/**
 * @file
 * @brief Contains the implementation of the TPZBlockDiagonalStructMatrix methods.
 */

#include "pzbdstrmatrix.h"
#include "pzblockdiag.h"
#include "pzsubcmesh.h"

using namespace std;

template<class TVar, class TPar>
void TPZBlockDiagonalStructMatrix<TVar,TPar>::AssembleBlockDiagonal(TPZBlockDiagonal<TVar> & block){
    
    TPZVec<int> blocksizes;
    BlockSizes(blocksizes);
    block.Initialize(blocksizes);
    int nblock = blocksizes.NElements();
    TPZAdmChunkVector<TPZCompEl *> &elementvec = this->fMesh->ElementVec();
    TPZAdmChunkVector<TPZConnect> &connectvec = this->fMesh->ConnectVec();
    TPZStack<int64_t> connectlist;
    TPZBlockDiagonal<TVar> elblock;
    int64_t numel = elementvec.NElements();
    int64_t el;
    for(el=0; el<numel; el++) {
        TPZCompEl *cel = elementvec[el];
        if(!cel) continue;
        TPZBlockDiagonal<TVar> eldiag;
        cel->CalcBlockDiagonal(connectlist,elblock);
        //    elblock.Print("Element block diagonal");
        int64_t ncon = connectlist.NElements();
        int64_t c,seqnum;
        for(c=0; c<ncon; c++) {
            TPZConnect &con = connectvec[connectlist[c]];
            seqnum = con.SequenceNumber();
            int64_t eqnum = this->fMesh->Block().Position(seqnum);
            if(seqnum <0 || seqnum >= nblock) continue;
            int bsize = blocksizes[seqnum];
            int64_t numactive = this->fEquationFilter.NumActive(eqnum, eqnum+bsize);
            if (!numactive) {
                continue;
            }
            if (numactive != bsize) {
                // Please implement me
                DebugStop();
            }
            if(con.NDof(*this->fMesh) != bsize ) {
                cout << "TPZBlockDiagonalStructMatrix::AssembleBlockDiagonal wrong data structure\n";
                continue;
            }
            TPZFMatrix<TVar> temp(bsize,bsize);
            elblock.GetBlock(c,temp);
            block.AddBlock(seqnum,temp);
        }
    }
}

template<class TVar, class TPar>
void TPZBlockDiagonalStructMatrix<TVar,TPar>::BlockSizes(TPZVec < int > & blocksizes){
    
    if(this->fMesh->FatherMesh() != 0) {
        TPZSubCompMesh *mesh = (TPZSubCompMesh *) this->fMesh;
        mesh->PermuteExternalConnects();
    }
    int nblocks = 0;
    TPZAdmChunkVector<TPZConnect> &connectvec = this->fMesh->ConnectVec();
    int64_t nc = connectvec.NElements();
    int64_t c;
    for(c=0; c<nc; c++) {
        TPZConnect &con = connectvec[c];
        if(con.HasDependency() || con.IsCondensed() || con.SequenceNumber() < 0) continue;
        nblocks++;
    }
    blocksizes.Resize(nblocks);
    int64_t bl,blsize;
    for(c=0; c<nc; c++) {
        TPZConnect &con = connectvec[c];
        if(con.HasDependency() || con.IsCondensed() || con.SequenceNumber() < 0) continue;
        bl = con.SequenceNumber();
        blsize = con.NDof(*this->fMesh);
        int64_t blpos = this->fMesh->Block().Position(bl);
        int64_t numactiv = this->fEquationFilter.NumActive(blpos, blpos+blsize);
        if (numactiv && numactiv != blsize) {
            DebugStop();
        }
        if(!numactiv)
        {
            blocksizes[bl] = 0;
        }
        else
        {
            blocksizes[bl] = blsize;
        }
    }
}

template<class TVar, class TPar>
TPZStructMatrix * TPZBlockDiagonalStructMatrix<TVar,TPar>::Clone(){
    return new TPZBlockDiagonalStructMatrix(*this);
}

template<class TVar, class TPar>
void TPZBlockDiagonalStructMatrix<TVar,TPar>::EndCreateAssemble(TPZBaseMatrix *mat){
    auto *block =
        dynamic_cast<TPZBlockDiagonal<TVar>*>(mat);
    AssembleBlockDiagonal(*block);
}

template<class TVar, class TPar>
TPZMatrix<TVar> * TPZBlockDiagonalStructMatrix<TVar,TPar>::Create(){
    TPZVec<int> blocksize;
    BlockSizes(blocksize);
    return new TPZBlockDiagonal<TVar>(blocksize);
}

template<class TVar, class TPar>
int TPZBlockDiagonalStructMatrix<TVar,TPar>::ClassId() const{
    return Hash("TPZBlockDiagonalStructMatrix") ^
        TPZStructMatrix::ClassId() << 1 ^
        TPar::ClassId() << 2;
}
template<class TVar, class TPar>
void TPZBlockDiagonalStructMatrix<TVar,TPar>::Read(TPZStream& buf, void* context){
    TPZStructMatrix::Read(buf,context);
    TPar::Read(buf,context);
}

template<class TVar, class TPar>
void TPZBlockDiagonalStructMatrix<TVar,TPar>::Write(TPZStream& buf, int withclassid) const{
    TPZStructMatrix::Write(buf,withclassid);
    TPar::Write(buf,withclassid);
}

#include "pzstrmatrixot.h"
#include "pzstrmatrixflowtbb.h"

template class TPZBlockDiagonalStructMatrix<STATE,TPZStructMatrixOR<STATE>>;
template class TPZBlockDiagonalStructMatrix<STATE,TPZStructMatrixOT<STATE>>;
template class TPZBlockDiagonalStructMatrix<STATE,TPZStructMatrixTBBFlow<STATE>>;
template class TPZBlockDiagonalStructMatrix<CSTATE,TPZStructMatrixOR<CSTATE>>;
template class TPZBlockDiagonalStructMatrix<CSTATE,TPZStructMatrixOT<CSTATE>>;
template class TPZBlockDiagonalStructMatrix<CSTATE,TPZStructMatrixTBBFlow<CSTATE>>;