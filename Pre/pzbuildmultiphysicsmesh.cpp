/*
 * @file
 * @brief Implementations to mesh multiphysics
 */

#include "pzbuildmultiphysicsmesh.h"
#include "pzmultiphysicselement.h"
#include "pzmaterial.h"

#include "TPZInterfaceEl.h"


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
		int iel;
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
						break;
					}
					celstack2.Pop();
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
		///ajustar as dependencias
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
    TPZBlock<REAL> &blockMF = MFMesh->Block();
    for (imesh = 0; imesh < nmeshes; imesh++) {
		int ncon = cmeshVec[imesh]->NConnects();
		TPZBlock<REAL> &block = cmeshVec[imesh]->Block();
		int ic;
		for (ic=0; ic<ncon; ic++) {
			TPZConnect &con = cmeshVec[imesh]->ConnectVec()[ic];
			int seqnum = con.SequenceNumber();
			if(seqnum<0) continue;       /// Whether connect was deleted by previous refined process
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
    TPZBlock<REAL> &blockMF = MFMesh->Block();
    for (imesh = 0; imesh < nmeshes; imesh++) {
		int ncon = cmeshVec[imesh]->NConnects();
		TPZBlock<REAL> &block = cmeshVec[imesh]->Block();
		int ic;
		for (ic=0; ic<ncon; ic++) {
			TPZConnect &con = cmeshVec[imesh]->ConnectVec()[ic];
			int seqnum = con.SequenceNumber();
			if(seqnum<0) continue;       /// Whether connect was deleted by previous refined process
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

void BuildHybridMesh(TPZCompMesh *cmesh, std::set<int> &MaterialIDs, int LagrangeMat, int InterfaceMat)
{
	TPZAdmChunkVector<TPZGeoEl *> &elvec = cmesh->Reference()->ElementVec();
	int meshdim = cmesh->Dimension();
	
	///cria todos os elementos sem conectar-se aos vizinhos
	int i, nelem = elvec.NElements();
	int neltocreate = 0;
	int index;
	for(i=0; i<nelem; i++) {
		TPZGeoEl *gel = elvec[i];
		if(!gel) continue;
		if(!gel->HasSubElement()) {
			neltocreate++;
		}
	}
	std::set<int> matnotfound;
	int nbl = cmesh->Block().NBlocks();
	if(neltocreate > nbl) cmesh->Block().SetNBlocks(neltocreate);
	cmesh->Block().SetNBlocks(nbl);
	for(i=0; i<nelem; i++) {
		TPZGeoEl *gel = elvec[i];
		if(!gel) continue;
		if(!gel->HasSubElement()) {
			int matid = gel->MaterialId();
			TPZAutoPointer<TPZMaterial> mat = cmesh->FindMaterial(matid);
			if(!mat)
			{
				matnotfound.insert(matid);
				continue;
			}
			int printing = 0;
			if (printing) {
				gel->Print(std::cout);
			}
			
			///checking material in MaterialIDs
			std::set<int>::const_iterator found = MaterialIDs.find(matid);
			if (found == MaterialIDs.end()) continue;
			
			if(!gel->Reference() && gel->NumInterfaces() == 0)
			{
				cmesh->CreateCompEl(gel,index);
				gel->ResetReference();
			}
		}
	}
    
	cmesh->LoadReferences();
    
	///Gera elementos geometricos para os elementos menores
	for (i=0; i<nelem; ++i) {
		TPZGeoEl *gel = elvec[i];
		if (!gel || gel->Dimension() != meshdim || !gel->Reference()) {
            continue;
		}
		int matid = gel->MaterialId();
		if(MaterialIDs.find(matid) == MaterialIDs.end())
		{
			continue;
		}
		///over the dimension-1 sides
		int nsides = gel->NSides();
		int is;
		for (is=0; is<nsides; ++is) {
			int sidedim = gel->SideDimension(is);
			if (sidedim != meshdim-1) {
				continue;
			}
			///check if there is a smaller element connected to this element
			TPZStack<TPZCompElSide> celsides;
			celsides.Resize(0);
			TPZGeoElSide gelside(gel,is);
			gelside.HigherLevelCompElementList2(celsides, 0, 0);
			//int ncelsid =  celsides.NElements();
			if(celsides.NElements()) continue;
			
			///check the neighboring
			TPZCompElSide celside;
			celside = gelside.LowerLevelCompElementList2(0);
			if (celside && celside.Element()->Reference()->Dimension() != meshdim) continue;
			TPZStack<TPZGeoElSide> allneigh;
			allneigh.Resize(0);	
			gelside.AllNeighbours(allneigh);
			int nneig = allneigh.NElements();
			if(allneigh.NElements()>1) continue;
			if (nneig && allneigh[0].Element()->Dimension() != meshdim) continue;
			
			gel->CreateBCGeoEl(is, LagrangeMat);
		}
	}
	
	/// now create the lagrange elements
	cmesh->Reference()->ResetReference();
	nelem = elvec.NElements();
	for(i=0; i<nelem; i++) {
		TPZGeoEl *gel = elvec[i];
		if(!gel) continue;
		if(!gel->HasSubElement()) {
			int matid = gel->MaterialId();
			TPZAutoPointer<TPZMaterial> mat = cmesh->FindMaterial(matid);
			
			if(!mat)
			{
				matnotfound.insert(matid);
				continue;
			}
			
			///checking material in MaterialIDs
			if (matid != LagrangeMat) {
				continue;
			}
			int printing = 0;
			if (printing) {
				gel->Print(std::cout);
			}
			
			if(!gel->Reference())
			{
				cmesh->CreateCompEl(gel,index);
				gel->ResetReference();
			}
		}
	}
	
	cmesh->LoadReferences();
	
	/// now create the interface elements between the lagrange elements and other elements
	nelem = elvec.NElements();
	for (i=0; i<nelem; ++i) {
		TPZGeoEl *gel = elvec[i];
		if (!gel || gel->Dimension() != meshdim-1 || !gel->Reference()) {
            continue;
		}
		int matid = gel->MaterialId();
		if(matid != LagrangeMat){
			continue;
		}
		///over the dimension-1 sides
		int nsides = gel->NSides();
		int is;
		for (is=0; is<nsides; ++is) {
			int sidedim = gel->SideDimension(is);
			if (sidedim != meshdim-1) {
				continue;
			}
			/// check if there is a smaller element connected to this element
			TPZStack<TPZCompElSide> celsides;
			TPZGeoElSide gelside(gel,is);
			gelside.EqualLevelCompElementList(celsides, 0, 0);
			if(celsides.size() < 1){
				DebugStop();
			}
			TPZCompElSide cels = gelside.LowerLevelCompElementList2(0);
			if(cels) 
			{
				celsides.Push(cels);
			}
			int nelsides = celsides.NElements();
			if(nelsides != 2) 
			{
				DebugStop();
			} 
			for (int lp=0; lp<nelsides; ++lp) {
				TPZGeoEl *interface = gel->CreateBCGeoEl(is, InterfaceMat);
				TPZCompElSide right = celsides[lp];
				TPZCompElSide left(gel->Reference(),is);
				int index;
				new TPZInterfaceElement(*cmesh,interface,index,left,right);
			}
		}
	}
	
	cmesh->InitializeBlock();
}


void TPZBuildMultiphysicsMesh::UniformRefineCompMesh(TPZCompMesh  *cMesh, int ndiv)
{
	TPZVec<int > subindex;
	for (int iref = 0; iref < ndiv; iref++) {
		TPZAdmChunkVector<TPZCompEl *> elvec = cMesh->ElementVec();
		int nel = elvec.NElements(); 
		for(int el=0; el < nel; el++){
			TPZCompEl * compEl = elvec[el];
			if(!compEl) continue;
			int ind = compEl->Index();
			compEl->Divide(ind, subindex, 0);
		}
	}
}

void TPZBuildMultiphysicsMesh::UniformRefineCompEl(TPZCompMesh  *cMesh, int indexEl){
	
	TPZVec<int > subindex; 
	int nel = cMesh->ElementVec().NElements(); 
	for(int el=0; el < nel; el++){
		TPZCompEl * compEl = cMesh->ElementVec()[el];
		if(!compEl) continue;
		int ind = compEl->Index();
		if(ind==indexEl){
			compEl->Divide(indexEl, subindex, 1);
		}
	}	
}


