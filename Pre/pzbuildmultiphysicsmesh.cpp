 /*
 * @file
 * @brief Implementations to mesh multiphysics
 */

#include "pzbuildmultiphysicsmesh.h"
#include "pzmultiphysicselement.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzmaterial.h"

#include "TPZInterfaceEl.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.tpzbuildmultiphysicsmesh"));
#endif




#include <iostream>

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
            TPZCompEl *cel = MFMesh->ElementVec()[iel];
			TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *> (cel);
            TPZMultiphysicsInterfaceElement *mfint = dynamic_cast<TPZMultiphysicsInterfaceElement *>(cel);
			if(mfcel)
			{
                int found = 0;
                TPZGeoEl *gel = mfcel->Reference();
                TPZStack<TPZCompElSide> celstack;
                TPZGeoElSide gelside(gel,gel->NSides()-1);
                // if the geometric element has a reference, it is an obvious candidate
                if (gel->Reference()) {
                    celstack.Push(gelside.Reference());
                }
                // put all large and small elements on the stack
                gelside.ConnectedCompElementList(celstack, 0, 0);
                while (celstack.size())
                {
                    int ncel = celstack.size();
                    // the last element on the stack
                    TPZGeoElSide gelside = celstack[ncel-1].Reference();
                    TPZStack<TPZCompElSide> celstack2;
                    // put te last element on the new stack
                    celstack2.Push(celstack[ncel-1]);
                    celstack.Pop();
                    // put all equal level elements on the new stack
                    gelside.EqualLevelCompElementList(celstack2, 0, 0);
                    
                    while (celstack2.size()) 
                    {
                        int ncel2 = celstack2.size();
                        TPZGeoElSide gelside2 = celstack2[ncel2-1].Reference();
                        // put all elements in the stack - if there is one element, stop the search
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
                    mfcel->AddElement(0, imesh);
                }
            }
            else if (mfint) {
                //set up interface
            }
            else {
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
                MFMesh->ConnectVec()[counter].SetPressure(refcon.IsLagrMult());
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
		TPZMultiphysicsInterfaceElement *interfacee = dynamic_cast<TPZMultiphysicsInterfaceElement *> (MFMesh->ElementVec()[iel]);
		if (interfacee) {
			continue;
		}
		if (!cel) {
			DebugStop();
		}
		TPZStack<int> connectindexes;
		int imesh;
		for (imesh=0; imesh < nmeshes; imesh++) {
			TPZCompEl *celref = cel->ReferredElement(imesh);
            if (!celref) {
                continue;
            }
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
    TPZBlock<STATE> &blockMF = MFMesh->Block();
    for (imesh = 0; imesh < nmeshes; imesh++) {
		int ncon = cmeshVec[imesh]->NConnects();
		TPZBlock<STATE> &block = cmeshVec[imesh]->Block();
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
    TPZBlock<STATE> &blockMF = MFMesh->Block();
    for (imesh = 0; imesh < nmeshes; imesh++) {
		int ncon = cmeshVec[imesh]->NConnects();
		TPZBlock<STATE> &block = cmeshVec[imesh]->Block();
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

void TPZBuildMultiphysicsMesh::BuildHybridMesh(TPZCompMesh *cmesh, std::set<int> &MaterialIDs, std::set<int> &BCMaterialIds, int LagrangeMat, int InterfaceMat)
{
	TPZAdmChunkVector<TPZGeoEl *> &elvec = cmesh->Reference()->ElementVec();
	int meshdim = cmesh->Dimension();
    
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    cmesh->ApproxSpace().BuildMesh(*cmesh, MaterialIDs);
	
    cmesh->LoadReferences();
    
    cmesh->ApproxSpace().CreateDisconnectedElements(false);
    cmesh->ApproxSpace().BuildMesh(*cmesh, BCMaterialIds);
    
//#ifdef LOG4CXX
//    {
//        std::stringstream sout;
//        cmesh->Reference()->Print(sout);
//        LOGPZ_DEBUG(logger,sout.str())
//    }
//#endif

    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    
    int nelem = cmesh->Reference()->NElements();

	//2- Generate geometric elements (with dimension (meshdim-1)) between the previous elements.
	for (int i=0; i<nelem; ++i) {
		TPZGeoEl *gel = elvec[i];
        // skip all elements which are not volumetric
		if (!gel || gel->Dimension() != meshdim || !gel->Reference()) {
            continue;
		}
		int matid = gel->MaterialId();
		if(MaterialIDs.find(matid) == MaterialIDs.end())
		{
			continue;
		}
		// over the dimension-1 sides
		int nsides = gel->NSides();
		int is;
		for (is=0; is<nsides; ++is) {
			int sidedim = gel->SideDimension(is);
			if (sidedim != meshdim-1) {
				continue;
			}
			// check if there is a smaller element connected to this element
			TPZStack<TPZCompElSide> celsides;
			celsides.Resize(0);
			TPZGeoElSide gelside(gel,is);
			gelside.HigherLevelCompElementList2(celsides, 0, 0);
			//int ncelsid =  celsides.NElements();
            
            // we only treat elements which look at equal or larger sizes
			if(celsides.NElements()) continue;
			
			// check the neighboring
			TPZCompElSide celside;
			celside = gelside.LowerLevelCompElementList2(0);
			if (celside && celside.Element()->Reference()->Dimension() != meshdim) 
            {
                // there a larger sided boundary element --- strange lets see what happened
                DebugStop();
            }
            TPZStack<TPZCompElSide> equallevel;
            gelside.EqualLevelCompElementList(equallevel, 1, 0);
            if (equallevel.size() > 1) {
                // at this point there should only be maximum one neighbour along this type of side
                DebugStop();
            }
			TPZStack<TPZGeoElSide> allneigh;
			allneigh.Resize(0);	
			gelside.AllNeighbours(allneigh);
			//int nneig = allneigh.NElements();
            // this would mean the lagrange element was already created along this side
			if(allneigh.NElements()>1) continue;
            
            // create only lagrange interfaces between internal sides
            int neighmat = allneigh[0].Element()->MaterialId();
            if (MaterialIDs.find(neighmat) == MaterialIDs.end() ) {
                continue;
            }
			//if (nneig && allneigh[0].Element()->Dimension() != meshdim) continue;//joao, comentar essa linha permite criar elementos 1D(Lagrange) entre elemento de contorno e um elemento 2D
			
			gel->CreateBCGeoEl(is, LagrangeMat);
		}
	}
    
    cmesh->Reference()->ResetReference();
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    

    std::set<int> lagrangematids;
    lagrangematids.insert(LagrangeMat);
    cmesh->AutoBuild(lagrangematids);
    

	
	cmesh->LoadReferences();
    
#ifdef LOG4CXX
    {
        std::stringstream sout;
        cmesh->Reference()->Print(sout);
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif

    
//  std::ofstream arg2("gmeshB.txt");
//	cmesh->Reference()->Print(arg2);

	
	//4- Create the interface elements between the lagrange elements and other elements
	nelem = elvec.NElements();
	for (int i=0; i<nelem; ++i) {
		TPZGeoEl *gel = elvec[i];
		if (!gel || gel->Dimension() != meshdim-1 || !gel->Reference()) {
            continue;
		}
		int matid = gel->MaterialId();
		if(matid != LagrangeMat){
			continue;
		}
		//over the dimension-1 sides
		int nsides = gel->NSides();
        if(nsides!=3)
        {
            DebugStop();// como estamos no caso 2D todas as arestas são 1D
        }
		int is = nsides-1;
        int sidedim = gel->SideDimension(is);
        if (sidedim != meshdim-1) {
            DebugStop();
        }
        //check if there is a smaller element connected to this element
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
            TPZCompElSide right = celsides[lp];
            TPZCompElSide left(gel->Reference(),is);
            
            TPZGeoEl *gelright=right.Reference().Element();
            int matidleft=gel->MaterialId();// sempre é LagrangeMat
            int matidright =gelright->MaterialId();
            ///???? o InterfaceMaterial não esta fazendo o que preciso. Por isso nao estou usando matid !
            const int matid = cmesh->Reference()->InterfaceMaterial(matidleft, matidright );
            
            // there is no need to create a lagrange multiplier between an interior element and boundary element
            if(matid == 0) 
            {
                continue;
            }
            
            TPZGeoEl *interfaceEl = gel->CreateBCGeoEl(is, matid);
            int index;
            new TPZInterfaceElement(*cmesh,interfaceEl,index,left,right);
            
        }
	}
	
//  std::ofstream arg3("gmeshC.txt");
//	cmesh->Reference()->Print(arg3);

	cmesh->InitializeBlock();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    
    //Set the connect as a Lagrange multiplier  
    int nel = cmesh->NElements();
    for(int i=0; i<nel; i++){
        TPZCompEl *cel = cmesh->ElementVec()[i];
        if(!cel) continue;
        
        int mid = cel->Material()->Id();
        
        if(mid==LagrangeMat){
            
            int nsides = cel->Reference()->NSides();
            
            for(int i = 0; i<nsides; i++){
                TPZConnect &newnod = cel->Connect(i);
                newnod.SetPressure(true);
            }
        }
    }
}

#include "TPZCompElDisc.h"
void TPZBuildMultiphysicsMesh::UniformRefineCompMesh(TPZCompMesh *cMesh, int ndiv, bool isLagrMult)
{

    TPZAdmChunkVector<TPZCompEl *> elvec = cMesh->ElementVec();
    int nel = elvec.NElements();
    for(int el=0; el < nel; el++){
        TPZCompEl * compEl = elvec[el];
        if(!compEl) continue;
        
        if(compEl->IsInterface()){
            compEl->Reference()->ResetReference();
            delete compEl;
        }
    }
    
	TPZVec<int > subindex(0);
	for (int iref = 0; iref < ndiv; iref++) {
		TPZAdmChunkVector<TPZCompEl *> elvec = cMesh->ElementVec();
		int nel = elvec.NElements(); 
		for(int el=0; el < nel; el++){
			TPZCompEl * compEl = elvec[el];
			if(!compEl) continue;
			int ind = compEl->Index();
            if(compEl->Dimension() >0/* cMesh->Dimension()*/){
                compEl->Divide(ind, subindex, 0);
            }
            
		}
	}
    
    cMesh->AdjustBoundaryElements();
	cMesh->CleanUpUnconnectedNodes();
    
    
    //When is using one mesh with L2 space for pressure  
    if(isLagrMult==true)
    {
        int ncon = cMesh->NConnects();
        for(int i=0; i<ncon; i++)
        {
            TPZConnect &newnod = cMesh->ConnectVec()[i];
            newnod.SetPressure(true);
        }
        
        int nel = cMesh->NElements();
        for(int i=0; i<nel; i++){
            TPZCompEl *cel = cMesh->ElementVec()[i];
            if(!cel) continue;
            TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
            celdisc->SetConstC(1.);
            celdisc->SetCenterPoint(0, 0.);
            celdisc->SetCenterPoint(1, 0.);
            celdisc->SetCenterPoint(2, 0.);
            celdisc->SetTrueUseQsiEta();
        }

    }
}

void TPZBuildMultiphysicsMesh::UniformRefineCompEl(TPZCompMesh  *cMesh, int indexEl, bool isLagrMult){
	
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
    
	cMesh->AdjustBoundaryElements();
	cMesh->CleanUpUnconnectedNodes();
    
    
    //When is using one mesh with L2 space for pressure
    if(isLagrMult==true)
    {
        int ncon = cMesh->NConnects();
        for(int i=0; i<ncon; i++)
        {
            TPZConnect &newnod = cMesh->ConnectVec()[i];
            newnod.SetPressure(true);
        }
        
        int nel = cMesh->NElements();
        for(int i=0; i<nel; i++){
            TPZCompEl *cel = cMesh->ElementVec()[i];
            if(!cel) continue;
            TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
            celdisc->SetConstC(1.);
            celdisc->SetCenterPoint(0, 0.);
            celdisc->SetCenterPoint(1, 0.);
            celdisc->SetCenterPoint(2, 0.);
            celdisc->SetTrueUseQsiEta();
        }
    }
}


