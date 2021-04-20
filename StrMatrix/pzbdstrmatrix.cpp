/**
 * @file
 * @brief Contains the implementation of the TPZBlockDiagonalStructMatrix methods.
 */

#include "pzbdstrmatrix.h"
#include "pzblockdiag.h"
#include "pzsubcmesh.h"

using namespace std;

TPZBlockDiagonalStructMatrix::~TPZBlockDiagonalStructMatrix(){
    
}

void TPZBlockDiagonalStructMatrix::AssembleBlockDiagonal(TPZBlockDiagonal<STATE> & block){
    
    TPZVec<int> blocksizes;
    BlockSizes(blocksizes);
    block.Initialize(blocksizes);
    int nblock = blocksizes.NElements();
    TPZAdmChunkVector<TPZCompEl *> &elementvec = fMesh->ElementVec();
    TPZAdmChunkVector<TPZConnect> &connectvec = fMesh->ConnectVec();
    TPZStack<int64_t> connectlist;
    TPZBlockDiagonal<STATE> elblock;
    int64_t numel = elementvec.NElements();
    int64_t el;
    for(el=0; el<numel; el++) {
        TPZCompEl *cel = elementvec[el];
        if(!cel) continue;
        TPZBlockDiagonal<STATE> eldiag;
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
            TPZFMatrix<STATE> temp(bsize,bsize);
            elblock.GetBlock(c,temp);
            block.AddBlock(seqnum,temp);
        }
    }
}

void TPZBlockDiagonalStructMatrix::BlockSizes(TPZVec < int > & blocksizes){
    
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

TPZStructMatrix * TPZBlockDiagonalStructMatrix::Clone(){
    return new TPZBlockDiagonalStructMatrix(*this);
}
TPZBaseMatrix * TPZBlockDiagonalStructMatrix::CreateAssemble(TPZBaseMatrix &rhs,TPZAutoPointer<TPZGuiInterface> guiInterface){
    int64_t neq = fMesh->NEquations();
    TPZBlockDiagonal<STATE> *block = new TPZBlockDiagonal<STATE>();
    rhs.Redim(neq,1);
    Assemble(rhs,guiInterface);
    AssembleBlockDiagonal(*block);
    return block;
}
TPZMatrix<STATE> * TPZBlockDiagonalStructMatrix::Create(){
    TPZVec<int> blocksize;
    BlockSizes(blocksize);
    return new TPZBlockDiagonal<STATE>(blocksize);
}
TPZBlockDiagonalStructMatrix::TPZBlockDiagonalStructMatrix(TPZCompMesh *mesh) : 
TPZRegisterClassId(&TPZBlockDiagonalStructMatrix::ClassId), TPZStructMatrix(mesh),fBlockStructure(EVertexBased),fOverlap(0)
{
}

TPZBlockDiagonalStructMatrix::TPZBlockDiagonalStructMatrix(TPZAutoPointer<TPZCompMesh>mesh) : 
TPZRegisterClassId(&TPZBlockDiagonalStructMatrix::ClassId), TPZStructMatrix(mesh),fBlockStructure(EVertexBased),fOverlap(0)
{
}

TPZBlockDiagonalStructMatrix::TPZBlockDiagonalStructMatrix() : TPZRegisterClassId(&TPZBlockDiagonalStructMatrix::ClassId), 
TPZStructMatrix(),fBlockStructure(EVertexBased),fOverlap(0)
{
}

int TPZBlockDiagonalStructMatrix::ClassId() const{
    return Hash("TPZBlockDiagonalStructMatrix") ^ TPZStructMatrix::ClassId() << 1;
}
void TPZBlockDiagonalStructMatrix::Read(TPZStream& buf, void* context){
    TPZStructMatrix::Read(buf,context);
    TPZStructMatrixOR::Read(buf,context);
}

void TPZBlockDiagonalStructMatrix::Write(TPZStream& buf, int withclassid) const{
    TPZStructMatrix::Write(buf,withclassid);
    TPZStructMatrixOR::Write(buf,withclassid);
}