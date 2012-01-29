/**
 * @file
 * @brief Contains the implementation of the TPZBlockDiagonalStructMatrix methods. 
 */
//$Id: pzbdstrmatrix.cpp,v 1.10 2010-04-06 17:22:04 fortiago Exp $

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

void TPZBlockDiagonalStructMatrix::AssembleBlockDiagonal(TPZBlockDiagonal & block){
	
	TPZVec<int> blocksizes;
	BlockSizes(blocksizes);
	block.Initialize(blocksizes);
	int nblock = blocksizes.NElements();
	TPZAdmChunkVector<TPZCompEl *> &elementvec = fMesh->ElementVec();
	TPZAdmChunkVector<TPZConnect> &connectvec = fMesh->ConnectVec();
	TPZStack<int> connectlist;
	TPZBlockDiagonal elblock;
	int numel = elementvec.NElements();
	int el;
	for(el=0; el<numel; el++) {
		TPZCompEl *cel = elementvec[el];
		if(!cel) continue;
		TPZBlockDiagonal eldiag;
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
			TPZFMatrix temp(bsize,bsize);
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
		//        int bl = con.SequenceNumber();
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
TPZMatrix * TPZBlockDiagonalStructMatrix::CreateAssemble(TPZFMatrix &rhs,TPZAutoPointer<TPZGuiInterface> guiInterface){
	int neq = fMesh->NEquations();
	//cout << "TPZBlockDiagonalStructMatrix::CreateAssemble will not assemble the right hand side\n";
	TPZBlockDiagonal *block = new TPZBlockDiagonal();
	rhs.Redim(neq,1);
	//  TPZStructMatrix::Assemble(rhs, *fMesh,fMinEq,fMaxEq);
	Assemble(rhs,guiInterface);
	AssembleBlockDiagonal(*block);
	//  block->Print("Block Diagonal matrix");
	return block;
}
TPZMatrix * TPZBlockDiagonalStructMatrix::Create(){
	TPZVec<int> blocksize;
	BlockSizes(blocksize);
	return new TPZBlockDiagonal(blocksize);
}
TPZBlockDiagonalStructMatrix::TPZBlockDiagonalStructMatrix(TPZCompMesh *mesh) : TPZStructMatrix(mesh),fBlockStructure(EVertexBased),fOverlap(0)
{
	fMinEq = 0;
	fMaxEq = mesh->NEquations();
}
