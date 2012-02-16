/*
 *  pzbuildmultiphysicsmesh.cpp
 *  PZ
 *
 *  Created by Agnaldo on 10/31/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "pzbuildmultiphysicsmesh.h"
#include "pzmultiphysicselement.h"


TPZBuildMultiphysicsMesh::TPZBuildMultiphysicsMesh(){
}

TPZBuildMultiphysicsMesh::~TPZBuildMultiphysicsMesh(){
}

void TPZBuildMultiphysicsMesh::AddElements(TPZVec<TPZCompMesh *> cmeshVec, TPZCompMesh *MFMesh)
{
	TPZGeoMesh *gmesh = MFMesh->Reference();
	gmesh->ResetReference();
	int nMFEl = MFMesh->NElements();
	int nmesh = cmeshVec.size();
	int imesh;
	for(imesh = 0; imesh<nmesh; imesh++)
	{
		cmeshVec[imesh]->LoadReferences();
		int iel, is;
		for(iel=0; iel<nMFEl; iel++)
		{
			TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *> (MFMesh->ElementVec()[iel]);
			if(!mfcel)
			{
				DebugStop();
			}
            int found = 0;
			TPZGeoEl *gel = mfcel->Reference();
			TPZStack<TPZCompElSide> celstack;
			TPZGeoElSide gelside(gel,gel->NSides()-1);
			if (gel->Reference()) {
				celstack.Push(gelside.Reference());
			}
			gelside.ConnectedCompElementList(celstack, 0, 0);
			while (celstack.size())
            {
                int ncel = celstack.size();
				TPZGeoElSide gelside = celstack[ncel-1].Reference();
                TPZStack<TPZCompElSide> celstack2;
                celstack2.Push(celstack[ncel-1]);
                gelside.EqualLevelCompElementList(celstack2, 0, 0);
                while (celstack2.size()) 
                {
                    int ncel2 = celstack2.size();
                    TPZGeoElSide gelside2 = celstack2[ncel2-1].Reference();
                    if(gelside2.Element()->Dimension()==gel->Dimension()) {
                        mfcel->AddElement(gelside2.Element()->Reference(), imesh);
                        found = 1;
                        celstack2.Resize(0);
                        celstack.Resize(0);                        
                    }
                }
			}
            if (!found) {
                DebugStop();
            }
		}
		gmesh->ResetReference();
	}
}

void TPZBuildMultiphysicsMesh::AddConnects(TPZVec<TPZCompMesh *> cmeshVec, TPZCompMesh *MFMesh)
{
	int nmeshes = cmeshVec.size();
	TPZVec<int> FirstConnect(nmeshes,0);
	int nconnects = 0;
	int imesh;
	for (imesh=0; imesh<nmeshes; imesh++) 
	{
		FirstConnect[imesh] = nconnects;
		nconnects += cmeshVec[imesh]->ConnectVec().NElements();
	}
	MFMesh->ConnectVec().Resize(nconnects);
	MFMesh->Block().SetNBlocks(nconnects);
	int counter = 0;
	int seqnum = 0;
	for (imesh=0; imesh<nmeshes; imesh++) 
	{
		int ic;
		int nc = cmeshVec[imesh]->ConnectVec().NElements();
		for (ic=0; ic<nc; ic++) 
		{
			TPZConnect &refcon =  cmeshVec[imesh]->ConnectVec()[ic];
			MFMesh->ConnectVec()[counter] = refcon;
			if (refcon.SequenceNumber() >= 0) {
				MFMesh->ConnectVec()[counter].SetSequenceNumber(seqnum);
                MFMesh->ConnectVec()[counter].SetNState(refcon.NState());
                MFMesh->ConnectVec()[counter].SetNShape(refcon.NShape());
				int ndof = refcon.NDof(*cmeshVec[imesh]);
				MFMesh->Block().Set(seqnum,ndof);
				seqnum++;
			}
			counter++;
		}	
		// ajustar as dependencias
		for (ic=0; ic<nc; ic++) 
		{
			TPZConnect &cn = MFMesh->ConnectVec()[FirstConnect[imesh]+ic];
			if (cn.HasDependency()) 
			{
				TPZConnect::TPZDepend *dep = cn.FirstDepend();
				while (dep) {
					dep->fDepConnectIndex = dep->fDepConnectIndex+FirstConnect[imesh];
					dep = dep->fNext;
				}
			}
		}	
	}
	MFMesh->Block().SetNBlocks(seqnum);
	MFMesh->ExpandSolution();
	int iel;
	int nelem = MFMesh->NElements();
	for (iel = 0; iel < nelem; iel++) 
	{
		TPZMultiphysicsElement *cel = dynamic_cast<TPZMultiphysicsElement *> (MFMesh->ElementVec()[iel]);
		if (!cel) {
			DebugStop();
		}
		TPZStack<int> connectindexes;
		int imesh;
		for (imesh=0; imesh < nmeshes; imesh++) {
			TPZCompEl *celref = cel->ReferredElement(imesh);
			int ncon = celref->NConnects();
			int ic;
			for (ic=0; ic<ncon; ic++) {
				connectindexes.Push(celref->ConnectIndex(ic)+FirstConnect[imesh]);
			}
		}
		cel->SetConnectIndexes(connectindexes);
	}
}

void TPZBuildMultiphysicsMesh::TransferFromMeshes(TPZVec<TPZCompMesh *> &cmeshVec, TPZCompMesh *MFMesh)
{
    int imesh;
    int nmeshes = cmeshVec.size();
    TPZManVector<int> FirstConnectIndex(nmeshes+1,0);
    for (imesh = 0; imesh < nmeshes; imesh++) {
        FirstConnectIndex[imesh+1] = FirstConnectIndex[imesh]+cmeshVec[imesh]->NConnects();
    }
    TPZBlock &blockMF = MFMesh->Block();
    for (imesh = 0; imesh < nmeshes; imesh++) {
        int ncon = cmeshVec[imesh]->NConnects();
        TPZBlock &block = cmeshVec[imesh]->Block();
        int ic;
        for (ic=0; ic<ncon; ic++) {
            TPZConnect &con = cmeshVec[imesh]->ConnectVec()[ic];
            int seqnum = con.SequenceNumber();
			if(seqnum<0) continue;       // Whether connect was deleted by previous refined process
            int blsize = block.Size(seqnum);
            TPZConnect &conMF = MFMesh->ConnectVec()[FirstConnectIndex[imesh]+ic];
            int seqnumMF = conMF.SequenceNumber();
            int idf;
            for (idf=0; idf<blsize; idf++) {
                blockMF.Put(seqnumMF, idf, 0, block.Get(seqnum, idf, 0));
            }
        }
    }
}

void TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(TPZVec<TPZCompMesh *> &cmeshVec, TPZCompMesh *MFMesh)
{
    int imesh;
    int nmeshes = cmeshVec.size();
    TPZManVector<int> FirstConnectIndex(nmeshes+1,0);
    for (imesh = 0; imesh < nmeshes; imesh++) {
        FirstConnectIndex[imesh+1] = FirstConnectIndex[imesh]+cmeshVec[imesh]->NConnects();
    }
    TPZBlock &blockMF = MFMesh->Block();
    for (imesh = 0; imesh < nmeshes; imesh++) {
        int ncon = cmeshVec[imesh]->NConnects();
        TPZBlock &block = cmeshVec[imesh]->Block();
        int ic;
        for (ic=0; ic<ncon; ic++) {
            TPZConnect &con = cmeshVec[imesh]->ConnectVec()[ic];
            int seqnum = con.SequenceNumber();
			if(seqnum<0) continue;       // Whether connect was deleted by previous refined process
            int blsize = block.Size(seqnum);
            TPZConnect &conMF = MFMesh->ConnectVec()[FirstConnectIndex[imesh]+ic];
            int seqnumMF = conMF.SequenceNumber();
            int idf;
            for (idf=0; idf<blsize; idf++) {
                block.Put(seqnum, idf, 0, blockMF.Get(seqnumMF, idf, 0));
            }
        }
    }
    
}

void TPZBuildMultiphysicsMesh::RefinUniformElemComp(TPZCompMesh  *cMesh, int ndiv)
{
	TPZVec<int > subindex;
	for (int iref = 0; iref < ndiv; iref++) {
		TPZAdmChunkVector<TPZCompEl *> elvec = cMesh->ElementVec();
		int nel = elvec.NElements(); 
		for(int el=0; el < nel; el++){
			TPZCompEl * compEl = elvec[el];
			if(!compEl) continue;
			int ind = compEl->Index();
//			TPZGeoEl *geoel = compEl->Reference();
			compEl->Divide(ind, subindex, 0);
		}
	}
}

