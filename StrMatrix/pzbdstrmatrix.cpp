/**
 * @file
 * @brief Contains the implementation of the TPZBlockDiagonalStructMatrix methods. 
 */

#include "pzbdstrmatrix.h"
#include "pzblockdiag.h"
#include "pzvec.h"
#include "pzadmchunk.h"
#include "pzconnect.h"
#include "pzcmesh.h"
#include "pzskylmat.h"
#include "pzsubcmesh.h"
#include "pzgmesh.h"
#include "pzsolve.h"
#include "pzstepsolver.h"

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
	TPZStack<int> connectlist;
	TPZBlockDiagonal<STATE> elblock;
	int numel = elementvec.NElements();
	int el;
	for(el=0; el<numel; el++) {
		TPZCompEl *cel = elementvec[el];
		if(!cel) continue;
		TPZBlockDiagonal<STATE> eldiag;
		cel->CalcBlockDiagonal(connectlist,elblock);
		//    elblock.Print("Element block diagonal");
		int ncon = connectlist.NElements();
		int c,seqnum;
		for(c=0; c<ncon; c++) {
			TPZConnect &con = connectvec[connectlist[c]];
			seqnum = con.SequenceNumber();
			int eqnum = fMesh->Block().Position(seqnum);
			if(seqnum <0 || seqnum >= nblock) continue;
			if(HasRange() && (eqnum <fMinEq || eqnum >= fMaxEq)) continue;
			int bsize = blocksizes[seqnum];
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
    int nc = connectvec.NElements();
    int c;
    for(c=0; c<nc; c++) {
        TPZConnect &con = connectvec[c];
        if(con.HasDependency() || con.IsCondensed() || con.SequenceNumber() < 0) continue;
        nblocks++;
    }
    blocksizes.Resize(nblocks);
    int bl,blsize;
    for(c=0; c<nc; c++) {
        TPZConnect &con = connectvec[c];
        if(con.HasDependency() || con.IsCondensed() || con.SequenceNumber() < 0) continue;
        bl = con.SequenceNumber();
        blsize = con.NDof(*fMesh);
        int blpos = fMesh->Block().Position(bl);
        if(HasRange() && (blpos < fMinEq || blpos >= fMaxEq))
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
TPZMatrix<STATE> * TPZBlockDiagonalStructMatrix::CreateAssemble(TPZFMatrix<STATE> &rhs,TPZAutoPointer<TPZGuiInterface> guiInterface){
	int neq = fMesh->NEquations();
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
TPZBlockDiagonalStructMatrix::TPZBlockDiagonalStructMatrix(TPZCompMesh *mesh) : TPZStructMatrix(mesh),fBlockStructure(EVertexBased),fOverlap(0)
{
	fMinEq = 0;
	fMaxEq = mesh->NEquations();
}
