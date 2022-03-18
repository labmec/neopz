 /*
 * @file
 * @brief Implementations to mesh multiphysics
 */

#include "TPZMultiphysicsCompMesh.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzmultiphysicselement.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "TPZMaterial.h"
#include "TPZLinearAnalysis.h"
#include "pzstack.h"
#include "TPZInterfaceEl.h"
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"
#include "pzsubcmesh.h"

#include "pzelchdivbound2.h"
#include "pzshapequad.h"
#include "pzshapelinear.h"
#include "pzshapetriang.h"

#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.tpzbuildmultiphysicsmesh");
#endif




#include <iostream>

TPZBuildMultiphysicsMesh::TPZBuildMultiphysicsMesh(){
}

TPZBuildMultiphysicsMesh::~TPZBuildMultiphysicsMesh(){
}


void TPZBuildMultiphysicsMesh::AddElements(TPZVec<TPZCompMesh *> &cmeshVec, TPZCompMesh *MFMesh)
{
	TPZGeoMesh *gmesh = MFMesh->Reference();
	gmesh->ResetReference();
	int64_t nMFEl = MFMesh->NElements();
	int64_t nmesh = cmeshVec.size();
	int64_t imesh;
	for(imesh = 0; imesh<nmesh; imesh++)
	{
		cmeshVec[imesh]->LoadReferences();
		int64_t iel;
		for(iel=0; iel<nMFEl; iel++)
		{
            TPZCompEl * cel = MFMesh->ElementVec()[iel];
			TPZMultiphysicsElement * mfcel = dynamic_cast<TPZMultiphysicsElement *> (cel);
            TPZMultiphysicsInterfaceElement * mfint = dynamic_cast<TPZMultiphysicsInterfaceElement *>(cel);
			if(mfcel)
			{
                int64_t found = 0;
                TPZGeoEl * gel = mfcel->Reference();
                TPZStack<TPZCompElSide> celstack;
                TPZGeoElSide gelside(gel,gel->NSides()-1);
                // if the geometric element has a reference, it is an obvious candidate
                if (gel->Reference()) {
                    mfcel->AddElement(gel->Reference(), imesh);
                    continue;
                }
                else
                {
                    TPZGeoEl *gelF = gel;
                    while(gelF->Father())
                    {
                        gelF = gelF->Father();
                        if (gelF->Reference()) {
#ifdef PZDEBUG
                            if (gelF->MaterialId() != gel->MaterialId()) {
                                DebugStop();
                            }
#endif
                            mfcel->AddElement(gelF->Reference(), imesh);
                            found = true;
                            break;
                        }
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
    for (int64_t el = 0; el < nMFEl; el++) {
        TPZCompEl *cel = MFMesh->Element(el);
        TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if (!mfcel) {
            continue;
        }
        mfcel->InitializeIntegrationRule();
    }
}

void TPZBuildMultiphysicsMesh::AddConnects(TPZVec<TPZCompMesh *> &cmeshVec, TPZCompMesh *MFMesh)
{
	int64_t nmeshes = cmeshVec.size();
    MFMesh->SetNMeshes(nmeshes);
    
	TPZVec<int64_t> FirstConnect(nmeshes,0);
	int64_t nconnects = 0;
	int64_t imesh;
	for (imesh=0; imesh<nmeshes; imesh++) 
	{
		FirstConnect[imesh] = nconnects;
		nconnects += cmeshVec[imesh]->ConnectVec().NElements();
	}
	MFMesh->ConnectVec().Resize(nconnects);
	MFMesh->Block().SetNBlocks(nconnects);
	int64_t counter = 0;
	int64_t seqnum = 0;
	for (imesh=0; imesh<nmeshes; imesh++) 
	{
		int64_t ic;
		int64_t nc = cmeshVec[imesh]->ConnectVec().NElements();
		for (ic=0; ic<nc; ic++) 
		{
			TPZConnect &refcon =  cmeshVec[imesh]->ConnectVec()[ic];
			MFMesh->ConnectVec()[counter] = refcon;
			if (refcon.SequenceNumber() >= 0) {
				MFMesh->ConnectVec()[counter].SetSequenceNumber(seqnum);
				MFMesh->ConnectVec()[counter].SetNState(refcon.NState());
				MFMesh->ConnectVec()[counter].SetNShape(refcon.NShape());
                MFMesh->ConnectVec()[counter].SetLagrangeMultiplier(refcon.LagrangeMultiplier());
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
	int64_t iel;
	int64_t nelem = MFMesh->NElements();
	for (iel = 0; iel < nelem; iel++) 
	{
        TPZCompEl *celorig = MFMesh->ElementVec()[iel];
		TPZMultiphysicsElement *cel = dynamic_cast<TPZMultiphysicsElement *> (celorig);
		TPZMultiphysicsInterfaceElement *interfacee = dynamic_cast<TPZMultiphysicsInterfaceElement *> (celorig);
		if (interfacee) {
			continue;
		}
		if (!cel) {
			continue;
		}
		TPZStack<int64_t> connectindexes;
		int64_t imesh;
		for (imesh=0; imesh < nmeshes; imesh++) {
			TPZCompEl *celref = cel->ReferredElement(imesh);
            if (!celref) {
                continue;
            }
			int64_t ncon = celref->NConnects();
			int64_t ic;
			for (ic=0; ic<ncon; ic++) {
				connectindexes.Push(celref->ConnectIndex(ic)+FirstConnect[imesh]);
			}
		}
		cel->SetConnectIndexes(connectindexes);
	}
}

void TPZBuildMultiphysicsMesh::AppendConnects(TPZCompMesh *cmesh, TPZCompMesh *MFMesh)
{
    int64_t nmeshes = MFMesh->GetNMeshes();
    if(nmeshes<1) DebugStop();
    nmeshes +=1;
    
    //adding connects from cmesh to MFMesh
    int64_t nconnects_old = MFMesh->NConnects();
    int64_t nconnects = nconnects_old + cmesh->NConnects();
    
    MFMesh->ConnectVec().Resize(nconnects);
    MFMesh->Block().SetNBlocks(nconnects);
    
    int64_t counter = nconnects_old;
    int64_t seqnum  = nconnects_old;
    
    int64_t ic;
    int64_t nc = cmesh->NConnects();
    for (ic=0; ic<nc; ic++)
    {
        TPZConnect &refcon =  cmesh->ConnectVec()[ic];
        MFMesh->ConnectVec()[counter] = refcon;
        if (refcon.SequenceNumber() >= 0)
        {
            MFMesh->ConnectVec()[counter].SetSequenceNumber(seqnum);
            MFMesh->ConnectVec()[counter].SetNState(refcon.NState());
            MFMesh->ConnectVec()[counter].SetNShape(refcon.NShape());
            MFMesh->ConnectVec()[counter].SetLagrangeMultiplier(refcon.LagrangeMultiplier());
            int ndof = refcon.NDof(*cmesh);
            MFMesh->Block().Set(seqnum,ndof);
            seqnum++;
        }
        counter++;
    }
    
    ///ajustar as dependencias
    for (ic=0; ic<nc; ic++)
    {
        TPZConnect &cn = MFMesh->ConnectVec()[nconnects_old + ic];
        if (cn.HasDependency())
        {
            TPZConnect::TPZDepend *dep = cn.FirstDepend();
            while (dep)
            {
                dep->fDepConnectIndex = dep->fDepConnectIndex + nconnects_old;
                dep = dep->fNext;
            }
        }
    }

	MFMesh->Block().SetNBlocks(seqnum);
	MFMesh->ExpandSolution();
    TPZCompEl *celorig  = NULL;
    TPZMultiphysicsElement *cel = NULL;
    TPZCompEl *celref = NULL;
    TPZMultiphysicsInterfaceElement *interface1 = NULL;
    
    int64_t iel;
    TPZVec<int64_t> FirstConnect(nmeshes,0);
    int64_t nelem = MFMesh->NElements();
    for (iel = 0; iel < nelem; iel++)
    {
        celorig = MFMesh->ElementVec()[iel];
        cel = dynamic_cast<TPZMultiphysicsElement *> (celorig);
        if (!cel) {
            continue;
        }
        
        int64_t nfirstcon= 0;
        int64_t imesh;
        int nrefel = cel->NMeshes();
        if(nrefel!=nmeshes) DebugStop();
        for (imesh=0; imesh<nrefel; imesh++)
        {
            celref = cel->ReferredElement(imesh);
            if (!celref) {
                continue;
            }
            FirstConnect[imesh] = nfirstcon;
            nfirstcon += celref->Mesh()->NConnects();
        }
        break;
    }
	
    //Set connect index to multiphysics element
	for (iel = 0; iel < nelem; iel++)
	{
        celorig = MFMesh->ElementVec()[iel];
		cel = dynamic_cast<TPZMultiphysicsElement *> (celorig);
		interface1 = dynamic_cast<TPZMultiphysicsInterfaceElement *> (celorig);
		if (interface1) {
			continue;
		}
		if (!cel) {
			continue;
		}
		TPZStack<int64_t> connectindexes;
		
        int64_t imesh;
		for (imesh=0; imesh < nmeshes; imesh++) {
			celref = cel->ReferredElement(imesh);
            if (!celref) {
                continue;
            }
			int64_t ncon = celref->NConnects();
			int64_t ic;
			for (ic=0; ic<ncon; ic++) {
				connectindexes.Push(celref->ConnectIndex(ic)+FirstConnect[imesh]);
			}
		}
		cel->SetConnectIndexes(connectindexes);
	}
}

void TPZBuildMultiphysicsMesh::TransferFromMeshes(TPZVec<TPZCompMesh *> &cmeshVec,
                                                  TPZCompMesh *MFMesh)
{
    TPZBaseMatrix &solMF = MFMesh->Solution();
    
    auto *realSol = dynamic_cast<TPZFMatrix<STATE>*>(&solMF);
    auto *cplxSol = dynamic_cast<TPZFMatrix<CSTATE>*>(&solMF);
    if(realSol) TransferFromMeshesT<STATE>(cmeshVec,MFMesh);
    else TransferFromMeshesT<CSTATE>(cmeshVec,MFMesh);
}

template<class TVar>
void TPZBuildMultiphysicsMesh::TransferFromMeshesT(TPZVec<TPZCompMesh *> &cmeshVec,
                                                   TPZCompMesh *MFMesh)
{
    
    TPZVec<atomic_index> indexes;
    ComputeAtomicIndexes(MFMesh, indexes);
    int64_t nconnect = indexes.size();
    TPZBlock &blockMF = MFMesh->Block();
    TPZFMatrix<TVar> &solMF = MFMesh->Solution();
    for(int64_t connect = 0; connect < nconnect; connect++)
    {
        TPZCompMesh *atomic_mesh = indexes[connect].first;
        if(!atomic_mesh) continue;
		TPZBlock &block = atomic_mesh->Block();
        TPZFMatrix<TVar> &sol = atomic_mesh->Solution();
        TPZConnect &con = atomic_mesh->ConnectVec()[indexes[connect].second];
        int64_t seqnum = con.SequenceNumber();
        if(seqnum<0) DebugStop();       /// Whether connect was deleted by previous refined process
        int blsize = block.Size(seqnum);
        TPZConnect &conMF = MFMesh->ConnectVec()[connect];
        int64_t seqnumMF = conMF.SequenceNumber();
        if(seqnumMF < 0) DebugStop();
        for (int idf=0; idf<blsize; idf++) {
            auto getval = sol(block.Index(seqnum, idf));
            solMF(blockMF.Index(seqnumMF, idf)) = getval;
        }
	}
    
    
    // copy the solution of the submesh to father mesh
    if(1)
    {
        TPZSubCompMesh *msub = dynamic_cast<TPZSubCompMesh*>(MFMesh);
        if(msub){
            TPZCompMesh * fathermesh = msub->FatherMesh();
            //TODOCOMPLEX
            TPZFMatrix<TVar> &fathermeshSol = fathermesh->Solution();
            TPZCompEl *compel = dynamic_cast<TPZCompEl*>(msub);
            int nconnect = compel->NConnects();
            
            for(int ic=0; ic<nconnect ; ic++){
                
                int64_t fatherconIndex = compel->ConnectIndex(ic);
                int64_t submeshIndex = msub->InternalIndex(fatherconIndex);
                if(fatherconIndex == -1) DebugStop();
                //acessing the block on father mesh
                TPZBlock &blockfather = fathermesh->Block();
                TPZFMatrix<TVar> &solfather = fathermesh->Solution();
                
                TPZConnect &confather = fathermesh->ConnectVec()[fatherconIndex];
                int64_t seqnumfather = confather.SequenceNumber();
                int nblock = blockfather.Size(seqnumfather);
                //acessing the block on submesh
                TPZBlock &blocksub = msub->Block();
                TPZFMatrix<TVar> &solsub = ((TPZCompMesh *)(msub))->Solution();
                const TPZConnect &consub = msub->ConnectVec()[submeshIndex];
                const int64_t seqnumsub = consub.SequenceNumber();
                
                if(seqnumfather < 0) DebugStop();
                for(int idf=0 ; idf<nblock; idf++){
                    const int posfather = blockfather.Position(seqnumfather);
                    auto valsub = solsub(blocksub.Index(seqnumsub, idf));
                    fathermeshSol(posfather + idf) = valsub;
                }
            }
            
        }
    }
    
    int64_t nel = MFMesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = MFMesh->Element(el);
        TPZSubCompMesh *sub = dynamic_cast<TPZSubCompMesh *>(cel);
        if(sub)
        {
            TransferFromMeshes(cmeshVec, sub);
        }
    }
}

void TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(TPZVec<TPZCompMesh *> &cmeshVec, TPZCompMesh *MFMesh)
{
    TPZBaseMatrix &solMF = MFMesh->Solution();
    
    auto *realSol = dynamic_cast<TPZFMatrix<STATE>*>(&solMF);
    auto *cplxSol = dynamic_cast<TPZFMatrix<CSTATE>*>(&solMF);
    if(realSol) TransferFromMultiPhysicsT<STATE>(cmeshVec,MFMesh);
    else TransferFromMultiPhysicsT<CSTATE>(cmeshVec,MFMesh);
}

template<class TVar>
void TPZBuildMultiphysicsMesh::TransferFromMultiPhysicsT(TPZVec<TPZCompMesh *> &cmeshVec, TPZCompMesh *MFMesh)
{
    
    TPZVec<atomic_index> indexes;
    ComputeAtomicIndexes(MFMesh, indexes);
    int64_t nconnect = indexes.size();
    TPZBlock &blockMF = MFMesh->Block();
    TPZFMatrix<TVar> &solMF = MFMesh->Solution();
    for(int64_t connect = 0; connect < nconnect; connect++)
    {
        TPZCompMesh *atomic_mesh = indexes[connect].first;
        //TODOCOMPLEX
        if(!atomic_mesh) continue;
        TPZFMatrix<TVar> &atomic_mesh_sol = atomic_mesh->Solution();
        TPZBlock &block = atomic_mesh->Block();
        int64_t atomicindexconnect = indexes[connect].second;
        TPZConnect &con = atomic_mesh->ConnectVec()[atomicindexconnect];
        int64_t seqnum = con.SequenceNumber();
        if(seqnum<0) DebugStop();       /// Whether connect was deleted by previous refined process
        int blsize = block.Size(seqnum);
        TPZConnect &conMF = MFMesh->ConnectVec()[connect];
        int64_t seqnumMF = conMF.SequenceNumber();
        if(seqnumMF < 0) DebugStop();
        for (int idf=0; idf<blsize; idf++) {
            TVar val = solMF(blockMF.Index(seqnumMF, idf));
            int64_t pos = block.Position(seqnum);
            atomic_mesh_sol(pos+idf) = val;
        }
    }
    
    int64_t nel = MFMesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = MFMesh->Element(el);
        TPZSubCompMesh *sub = dynamic_cast<TPZSubCompMesh *>(cel);
        if(sub)
        {
            TransferFromMultiPhysics(cmeshVec, sub);
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
    
//#ifdef PZ_LOG
//    {
//        std::stringstream sout;
//        cmesh->Reference()->Print(sout);
//        LOGPZ_DEBUG(logger,sout.str())
//    }
//#endif

    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    
    int64_t nelem = cmesh->Reference()->NElements();

	//2- Generate geometric elements (with dimension (meshdim-1)) between the previous elements.
	for (int64_t i=0; i<nelem; ++i) {
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
    
#ifdef PZ_LOG
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
	for (int64_t i=0; i<nelem; ++i) {
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
        int64_t nelsides = celsides.NElements();
        if(nelsides != 2) 
        {
            DebugStop();
        } 
        for (int64_t lp=0; lp<nelsides; ++lp) {
            TPZCompElSide left = celsides[lp];
            TPZCompElSide right(gel->Reference(),is);
            
            TPZGeoEl *gelright=right.Reference().Element();
            TPZGeoEl *gelleft = left.Reference().Element();
            int matidleft = gelleft->MaterialId();// sempre é LagrangeMat
            int matidright = gelright->MaterialId();
            ///???? o InterfaceMaterial não esta fazendo o que preciso. Por isso nao estou usando matid !
            const int interfacematid = cmesh->Reference()->InterfaceMaterial(matidleft, matidright );
            
            // there is no need to create a lagrange multiplier between an interior element and boundary element
            if(interfacematid == 0) 
            {
                DebugStop();
            }
            
            TPZGeoEl *interfaceEl = gel->CreateBCGeoEl(is, interfacematid);            
            new TPZInterfaceElement(*cmesh,interfaceEl,left,right);
            
        }
	}
	
//  std::ofstream arg3("gmeshC.txt");
//	cmesh->Reference()->Print(arg3);

	cmesh->InitializeBlock();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    
    //Set the connect as a Lagrange multiplier  
    int64_t nel = cmesh->NElements();
    for(int64_t i=0; i<nel; i++){
        TPZCompEl *cel = cmesh->ElementVec()[i];
        if(!cel || !cel->Material())
            continue;
        
        int mid = cel->Material()->Id();
        
        if(mid==LagrangeMat){
            
            int nsides = cel->Reference()->NSides();
            
            for(int i = 0; i<nsides; i++){
                TPZConnect &newnod = cel->Connect(i);
                newnod.SetLagrangeMultiplier(2);
            }
        }
    }
}

#include "TPZCompElDisc.h"
void TPZBuildMultiphysicsMesh::UniformRefineCompMesh(TPZCompMesh *cMesh, int ndiv, bool isLagrMult)
{

    // delete the interface elements
    TPZAdmChunkVector<TPZCompEl *> elvec = cMesh->ElementVec();
    int64_t nel = elvec.NElements();
    for(int64_t el=0; el < nel; el++){
        TPZCompEl * compEl = elvec[el];
        if(!compEl) continue;
        
        if(compEl->IsInterface()){
            compEl->Reference()->ResetReference();
            delete compEl;
        }
    }
    
    // divide all elements
	TPZVec<int64_t > subindex(0);
	for (int iref = 0; iref < ndiv; iref++) {
		TPZAdmChunkVector<TPZCompEl *> elvec = cMesh->ElementVec();
		int64_t nel = elvec.NElements(); 
		for(int64_t el=0; el < nel; el++){
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
    
    
    //When the mesh is an L2 space used as lagrange multiplier  
    if(isLagrMult==true)
    {
        int64_t ncon = cMesh->NConnects();
        for(int64_t i=0; i<ncon; i++)
        {
            TPZConnect &newnod = cMesh->ConnectVec()[i];
            newnod.SetLagrangeMultiplier(1);
        }
        
        int64_t nel = cMesh->NElements();
        for(int64_t i=0; i<nel; i++){
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

void TPZBuildMultiphysicsMesh::UniformRefineCompEl(TPZCompMesh  *cMesh, int64_t indexEl, bool isLagrMult){
	
	TPZVec<int64_t> subindex; 
	int64_t nel = cMesh->ElementVec().NElements(); 
	for(int64_t el=0; el < nel; el++){
		TPZCompEl * compEl = cMesh->ElementVec()[el];
		if(!compEl) continue;
		int64_t ind = compEl->Index();
		if(ind==indexEl){
            //-------------------------------------------
            TPZStack<TPZCompElSide> neighequal;
            for(int side = compEl->Reference()->NSides()-2; side > compEl->Reference()->NCornerNodes()-1; side--)
            {
                TPZVec<int64_t> subindexneigh;
                int64_t indneigh;
                neighequal.Resize(0);
                TPZCompElSide celside(compEl,side);
                celside.EqualLevelElementList(neighequal, 0, 0);
                
                int nneighs = neighequal.size();
                if(nneighs != 0)
                {
                    for(int i =0; i<nneighs; i++)
                    {
                        TPZInterpolationSpace * InterpEl = dynamic_cast<TPZInterpolationSpace *>(neighequal[i].Element());
                        if(!InterpEl || InterpEl->Dimension() == compEl->Dimension()) continue;
                        indneigh = InterpEl->Index();
                        InterpEl->Divide(indneigh, subindexneigh, 1);
                    }
                }
            }
            //-------------------------------------------
			compEl->Divide(indexEl, subindex, 1);
            break;
		}
	}
    
	cMesh->AdjustBoundaryElements();
	cMesh->CleanUpUnconnectedNodes();
    
    
    //When is using one mesh with L2 space for pressure
    if(isLagrMult==true)
    {
        int64_t ncon = cMesh->NConnects();
        for(int64_t i=0; i<ncon; i++)
        {
            TPZConnect &newnod = cMesh->ConnectVec()[i];
            newnod.SetLagrangeMultiplier(1);
        }
        
        int64_t nel = cMesh->NElements();
        for(int64_t i=0; i<nel; i++){
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


/**
 * @brief Show shape functions associated with connects of a multiphysics mesh
 */
void TPZBuildMultiphysicsMesh::ShowShape(TPZVec<TPZCompMesh *> &cmeshVec, TPZCompMesh *MFMesh, TPZLinearAnalysis &analysis, const std::string &filename, TPZVec<int64_t> &equationindices)
{
    
    TPZStack<std::string> scalnames,vecnames;
    scalnames.Push("State");
    analysis.DefineGraphMesh(analysis.Mesh()->Dimension(), scalnames, vecnames, filename);
    int porder = analysis.Mesh()->GetDefaultOrder();
    
    int neq = equationindices.size();
    TPZFMatrix<STATE> solkeep(analysis.Solution());
    analysis.Solution().Zero();
    //TODOCOMPLEX
    TPZFMatrix<STATE> &anSol = analysis.Solution();
    for (int ieq = 0; ieq < neq; ieq++) {
        anSol(equationindices[ieq],0) = 1.;
        analysis.LoadSolution();
        TransferFromMultiPhysics(cmeshVec, MFMesh);
        analysis.PostProcess(porder);
        analysis.Solution().Zero();
    }
    analysis.Solution() = solkeep;
    analysis.LoadSolution();

}


void TPZBuildMultiphysicsMesh::AddWrap(TPZMultiphysicsElement *mfcel, int matskeleton, TPZStack< TPZStack<TPZMultiphysicsElement *,7> > &ListGroupEl)
{
    TPZCompMesh *multiMesh = mfcel->Mesh();
    TPZInterpolationSpace *hdivel = dynamic_cast<TPZInterpolationSpace *> (mfcel->Element(0));
    TPZCompEl *cel = dynamic_cast<TPZCompEl *>(mfcel->Element(1));
    MElementType celType = cel->Type();
    
    TPZGeoEl *gel = mfcel->Reference();
    
    int dimMesh = mfcel->Mesh()->Dimension();
    if (!hdivel || !cel || gel->Dimension() != dimMesh) {
        DebugStop();
    }
    
    //wrap element
    TPZStack<TPZMultiphysicsElement *, 7> wrapEl;
    wrapEl.push_back(mfcel);
    
    for (int side = 0; side < gel->NSides(); side++)
    {
        if (gel->SideDimension(side) != gel->Dimension()-1) {
            continue;
        }
        TPZGeoEl *gelbound = gel->CreateBCGeoEl(side, matskeleton);
        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(hdivel);
        int loccon = intel->SideConnectLocId(0,side);
        
        //Impose that the skeleton element has the same polynomial order  to the element of side.
        TPZConnect &conside = intel->Connect(loccon);
        int sideorder = conside.Order();
        intel->Mesh()->SetDefaultOrder(sideorder);
                
        TPZInterpolationSpace *bound;
        MElementType elType = gel->Type(side);
        switch(elType)
        {
            case(EOned)://line
            {
                bound = new TPZCompElHDivBound2<pzshape::TPZShapeLinear>(* intel->Mesh(),gelbound);
                int sideorient = intel->GetSideOrient(side);
                TPZCompElHDivBound2<pzshape::TPZShapeLinear> *hdivbound = dynamic_cast< TPZCompElHDivBound2<pzshape::TPZShapeLinear> *>(bound);
                hdivbound->SetSideOrient(pzshape::TPZShapeLinear::NSides-1,sideorient);
                break;
            }
            case(ETriangle)://triangle
            {
                bound = new TPZCompElHDivBound2<pzshape::TPZShapeTriang>(* intel->Mesh(),gelbound);
                int sideorient = intel->GetSideOrient(side);
                TPZCompElHDivBound2<pzshape::TPZShapeTriang> *hdivbound = dynamic_cast< TPZCompElHDivBound2<pzshape::TPZShapeTriang> *>(bound);
                hdivbound->SetSideOrient(pzshape::TPZShapeTriang::NSides-1,sideorient);
                break;
            }
            case(EQuadrilateral)://quadrilateral
            {
                bound = new TPZCompElHDivBound2<pzshape::TPZShapeQuad>(* intel->Mesh(),gelbound);
                int sideorient = intel->GetSideOrient(side);
                TPZCompElHDivBound2<pzshape::TPZShapeQuad> *hdivbound = dynamic_cast< TPZCompElHDivBound2<pzshape::TPZShapeQuad> *>(bound);
                hdivbound->SetSideOrient(pzshape::TPZShapeQuad::NSides-1,sideorient);
                break;
            }
                
            default:
            {
                bound=0;
                std::cout << "ElementType not found!";
                DebugStop();
                break;
            }
        }
        
        int64_t sideconnectindex = intel->ConnectIndex(loccon);
        
        TPZConnect &co = bound->Connect(0);
        if(co.HasDependency()){
            if(bound->NConnects()!=1) DebugStop();
            //int64_t cindex_bound = bound->ConnectIndex(0);
            co.RemoveDepend();
        }
        bound->SetConnectIndex(0, sideconnectindex);
        //bound->Print(std::cout);
        
        TPZCompEl *newMFBound = multiMesh->CreateCompEl(gelbound);
        TPZMultiphysicsElement *locMF = dynamic_cast<TPZMultiphysicsElement *>(newMFBound);
        
        locMF->AddElement(bound, 0);
        //locMF->Print(std::cout);
        
        if(celType==EDiscontinuous){
            TPZCompElDisc *discel = dynamic_cast<TPZCompElDisc *>(mfcel->Element(1));
            locMF->AddElement(TPZCompElSide(discel,side), 1);
        }else{
            TPZInterpolationSpace *contcel = dynamic_cast<TPZInterpolationSpace *>(mfcel->Element(1));
            locMF->AddElement(TPZCompElSide(contcel,side), 1);
        }
        
        wrapEl.push_back(locMF);
    }
    
    ListGroupEl.push_back(wrapEl);
}

static void FillAtomic(TPZCompEl *cel, TPZVec<atomic_index> &indexes);

static void FillAtomic(TPZMultiphysicsElement *mphys, TPZVec<atomic_index> &indexes)
{
    int ncon = mphys->NConnects();
    int count = 0;
    int nmeshes = mphys->NMeshes();
    for (int imesh = 0; imesh<nmeshes; imesh++) {
        if(mphys->IsActiveApproxSpaces(imesh) == false) continue;
        TPZCompEl *cel = mphys->Element(imesh);
        if(!cel) continue;
        TPZCompMesh *atomic_mesh = cel->Mesh();
        int nc = cel->NConnects();
        for (int ic=0; ic<nc; ic++) {
            int64_t atomic_conindex = cel->ConnectIndex(ic);
            int64_t conindex = mphys->ConnectIndex(count);
            indexes[conindex] = atomic_index(atomic_mesh,atomic_conindex);
            count++;
        }
    }
    if(count != ncon) DebugStop();
}

static void FillAtomic(TPZMultiphysicsInterfaceElement *intface, TPZVec<atomic_index> &indexes)
{
    int ncon = intface->NConnects();
    TPZCompElSide leftside, rightside;
    intface->GetLeftRightElement(leftside, rightside);
    TPZManVector<int64_t> leftindices, rightindices;
    intface->GetLeftRightElementIndices(leftindices, rightindices);
    TPZCompEl *left = leftside.Element();
    TPZMultiphysicsElement *mphysleft = dynamic_cast<TPZMultiphysicsElement *>(left);
    if(!mphysleft) DebugStop();
    TPZCompEl *right = rightside.Element();
    TPZMultiphysicsElement *mphysright = dynamic_cast<TPZMultiphysicsElement *>(right);
    if(!mphysright) DebugStop();

    int count = 0;
    for(int64_t i = 0; i<leftindices.size(); i++)
    {
        int64_t imesh = leftindices[i];
        TPZCompEl *cel = mphysleft->Element(imesh);
        if(!cel) continue;
        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
        int nside_connects = intel->NSideConnects(leftside.Side());
        for(int ic=0; ic<nside_connects; ic++)
        {
            int64_t atomic_conindex = intel->SideConnectIndex(ic, leftside.Side());
            int64_t mphys_index = intface->ConnectIndex(count);
            indexes[mphys_index] = atomic_index(cel->Mesh(),atomic_conindex);
            count++;
        }
    }
    for(int64_t i = 0; i<rightindices.size(); i++)
    {
        int64_t imesh = rightindices[i];
        TPZCompEl *cel = mphysright->Element(imesh);
        if(!cel) continue;
        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
        int nside_connects = intel->NSideConnects(rightside.Side());
        for(int ic=0; ic<nside_connects; ic++)
        {
            int64_t atomic_conindex = intel->SideConnectIndex(ic, rightside.Side());
            int64_t mphys_index = intface->ConnectIndex(count);
            indexes[mphys_index] = atomic_index(cel->Mesh(),atomic_conindex);
            count++;
        }
    }
#ifdef PZDEBUG
    if(count != intface->NConnects()) DebugStop();
#endif
}

static void FillAtomic(TPZElementGroup *elgr, TPZVec<atomic_index> &indexes)
{
    TPZVec<TPZCompEl *> elvec = elgr->GetElGroup();
    for (int el=0; el<elvec.size(); el++) {
        TPZCompEl *cel = elvec[el];
        FillAtomic(cel, indexes);
    }
}

static void FillAtomic(TPZCondensedCompEl *cond, TPZVec<atomic_index> &indexes)
{
    TPZCompEl *cel = cond->ReferenceCompEl();
    FillAtomic(cel, indexes);
}

static void FillAtomic(TPZCompEl *cel, TPZVec<atomic_index> &indexes)
{
    TPZMultiphysicsElement *mphys = dynamic_cast<TPZMultiphysicsElement *>(cel);
    TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *>(cel);
    TPZCondensedCompEl *condense = dynamic_cast<TPZCondensedCompEl *>(cel);
    TPZMultiphysicsInterfaceElement *intface = dynamic_cast<TPZMultiphysicsInterfaceElement *>(cel);
    if(mphys)
    {
        FillAtomic(mphys, indexes);
    }
    if(elgr)
    {
        FillAtomic(elgr, indexes);
    }
    if(condense)
    {
        FillAtomic(condense, indexes);
    }
    if(intface)
    {
        FillAtomic(intface, indexes);
    }
}

/**
 * Compute the correspondence between the connect index in the multiphysics
 * mesh and the connect indexes in the atomic meshes
 */
void TPZBuildMultiphysicsMesh::ComputeAtomicIndexes(TPZCompMesh *mesh, TPZVec<atomic_index> &indexes)
{
    int64_t ncon = mesh->NConnects();
    int64_t nel = mesh->NElements();
    atomic_index def(0,-1);
    indexes.Resize(ncon, def);
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = mesh->Element(el);
        if(cel)
        {
            FillAtomic(cel, indexes);
        }
    }
#ifdef PZDEBUG
    {
        int notfound = 0;
        for (int64_t i=0; i<indexes.size(); i++) {
            if(mesh->ConnectVec()[i].SequenceNumber() < 0) continue;
            if(indexes[i].first == 0) notfound++;
        }
        if(notfound)
        {
            std::cout << __PRETTY_FUNCTION__ << " number of missing connects " << notfound << std::endl;
        }
    }
#endif
}
