/**
 * @file
 * @brief Contains the implementation of the TPZBlockDiagonalStructMatrix methods.
 */

#include "pzbdstrmatrix.h"
#include "pzblockdiag.h"
#include "pzsubcmesh.h"

using namespace std;

template<class TVar, class TPar>
TPZBlockDiagonalStructMatrix<TVar,TPar>::~TPZBlockDiagonalStructMatrix(){
    
}

template<class TVar, class TPar>
void TPZBlockDiagonalStructMatrix<TVar,TPar>::AssembleBlockDiagonal(TPZBlockDiagonal<TVar> & block){
    
    TPZVec<int> blocksizes;
    BlockSizes(blocksizes);
    block.Initialize(blocksizes);
    int nblock = blocksizes.NElements();
    TPZAdmChunkVector<TPZCompEl *> &elementvec = fMesh->ElementVec();
    TPZAdmChunkVector<TPZConnect> &connectvec = fMesh->ConnectVec();
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
            int64_t eqnum = fMesh->Block().Position(seqnum);
            if(seqnum <0 || seqnum >= nblock) continue;
            int bsize = blocksizes[seqnum];
            int64_t numactive = fEquationFilter.NumActive(eqnum, eqnum+bsize);
            if (!numactive) {
                continue;
            }
            if (numactive != bsize) {
                // Please implement me
                DebugStop();
            }
            if(con.NDof(*fMesh) != bsize ) {
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
    
    if(fMesh->FatherMesh() != 0) {
        TPZSubCompMesh *mesh = (TPZSubCompMesh *) fMesh;
        mesh->PermuteExternalConnects();
    }
    int nblocks = 0;
    TPZAdmChunkVector<TPZConnect> &connectvec = fMesh->ConnectVec();
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
        blsize = con.NDof(*fMesh);
        int64_t blpos = fMesh->Block().Position(bl);
        int64_t numactiv = fEquationFilter.NumActive(blpos, blpos+blsize);
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
TPZMatrix<TVar> * TPZBlockDiagonalStructMatrix<TVar,TPar>::CreateAssemble(TPZBaseMatrix &rhs,TPZAutoPointer<TPZGuiInterface> guiInterface){
    int64_t neq = fMesh->NEquations();
    TPZBlockDiagonal<TVar> *block = new TPZBlockDiagonal<TVar>();
    rhs.Redim(neq,1);
    Assemble(rhs,guiInterface);
    AssembleBlockDiagonal(*block);
    return block;
}
template<class TVar, class TPar>
TPZMatrix<TVar> * TPZBlockDiagonalStructMatrix<TVar,TPar>::Create(){
    TPZVec<int> blocksize;
    BlockSizes(blocksize);
    return new TPZBlockDiagonal<TVar>(blocksize);
}
template<class TVar, class TPar>
TPZBlockDiagonalStructMatrix<TVar,TPar>::TPZBlockDiagonalStructMatrix(TPZCompMesh *mesh) : 
TPZRegisterClassId(&TPZBlockDiagonalStructMatrix::ClassId), TPZStructMatrix(mesh),fBlockStructure(EVertexBased),fOverlap(0)
{
}

template<class TVar, class TPar>
TPZBlockDiagonalStructMatrix<TVar,TPar>::TPZBlockDiagonalStructMatrix(TPZAutoPointer<TPZCompMesh>mesh) : 
TPZRegisterClassId(&TPZBlockDiagonalStructMatrix::ClassId), TPZStructMatrix(mesh),fBlockStructure(EVertexBased),fOverlap(0)
{
}

template<class TVar, class TPar>
TPZBlockDiagonalStructMatrix<TVar,TPar>::TPZBlockDiagonalStructMatrix() : TPZRegisterClassId(&TPZBlockDiagonalStructMatrix::ClassId), 
TPZStructMatrix(),fBlockStructure(EVertexBased),fOverlap(0)
{
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
