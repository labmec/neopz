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
	TPZStack<long> connectlist;
	TPZBlockDiagonal<STATE> elblock;
	long numel = elementvec.NElements();
	long el;
	for(el=0; el<numel; el++) {
		TPZCompEl *cel = elementvec[el];
		if(!cel) continue;
		TPZBlockDiagonal<STATE> eldiag;
		cel->CalcBlockDiagonal(connectlist,elblock);
		//    elblock.Print("Element block diagonal");
		long ncon = connectlist.NElements();
		long c,seqnum;
		for(c=0; c<ncon; c++) {
			TPZConnect &con = connectvec[connectlist[c]];
			seqnum = con.SequenceNumber();
			long eqnum = fMesh->Block().Position(seqnum);
			if(seqnum <0 || seqnum >= nblock) continue;
			int bsize = blocksizes[seqnum];
            long numactive = fEquationFilter.NumActive(eqnum, eqnum+bsize);
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
    long nc = connectvec.NElements();
    long c;
    for(c=0; c<nc; c++) {
        TPZConnect &con = connectvec[c];
        if(con.HasDependency() || con.IsCondensed() || con.SequenceNumber() < 0) continue;
        nblocks++;
    }
    blocksizes.Resize(nblocks);
    long bl,blsize;
    for(c=0; c<nc; c++) {
        TPZConnect &con = connectvec[c];
        if(con.HasDependency() || con.IsCondensed() || con.SequenceNumber() < 0) continue;
        bl = con.SequenceNumber();
        blsize = con.NDof(*fMesh);
        long blpos = fMesh->Block().Position(bl);
        long numactiv = fEquationFilter.NumActive(blpos, blpos+blsize);
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
TPZMatrix<STATE> * TPZBlockDiagonalStructMatrix::CreateAssemble(TPZFMatrix<STATE> &rhs,TPZAutoPointer<TPZGuiInterface> guiInterface){
	long neq = fMesh->NEquations();
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
}
