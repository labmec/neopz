 /*
 * @file
 * @brief Implementations to mesh multiphysics
 */

#include "pzbuildmultiphysicsmesh.h"
#include "pzmultiphysicselement.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzmaterial.h"
#include "pzanalysis.h"
#include "pzstack.h"
#include "TPZInterfaceEl.h"

#include "pzelchdivbound2.h"
#include "pzshapequad.h"
#include "pzshapelinear.h"
#include "pzshapetriang.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.tpzbuildmultiphysicsmesh"));
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
	long nMFEl = MFMesh->NElements();
	long nmesh = cmeshVec.size();
	long imesh;
	for(imesh = 0; imesh<nmesh; imesh++)
	{
		cmeshVec[imesh]->LoadReferences();
		long iel;
		for(iel=0; iel<nMFEl; iel++)
		{
            TPZCompEl * cel = MFMesh->ElementVec()[iel];
			TPZMultiphysicsElement * mfcel = dynamic_cast<TPZMultiphysicsElement *> (cel);
            TPZMultiphysicsInterfaceElement * mfint = dynamic_cast<TPZMultiphysicsInterfaceElement *>(cel);
			if(mfcel)
			{
                long found = 0;
                TPZGeoEl * gel = mfcel->Reference();
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
                    long ncel = celstack.size();
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
                        long ncel2 = celstack2.size();
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

void TPZBuildMultiphysicsMesh::AddConnects(TPZVec<TPZCompMesh *> &cmeshVec, TPZCompMesh *MFMesh)
{
	long nmeshes = cmeshVec.size();
    MFMesh->SetNMeshes(nmeshes);
    
	TPZVec<long> FirstConnect(nmeshes,0);
	long nconnects = 0;
	long imesh;
	for (imesh=0; imesh<nmeshes; imesh++) 
	{
		FirstConnect[imesh] = nconnects;
		nconnects += cmeshVec[imesh]->ConnectVec().NElements();
	}
	MFMesh->ConnectVec().Resize(nconnects);
	MFMesh->Block().SetNBlocks(nconnects);
	long counter = 0;
	long seqnum = 0;
	for (imesh=0; imesh<nmeshes; imesh++) 
	{
		long ic;
		long nc = cmeshVec[imesh]->ConnectVec().NElements();
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
	long iel;
	long nelem = MFMesh->NElements();
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
		TPZStack<long> connectindexes;
		long imesh;
		for (imesh=0; imesh < nmeshes; imesh++) {
			TPZCompEl *celref = cel->ReferredElement(imesh);
            if (!celref) {
                continue;
            }
			long ncon = celref->NConnects();
			long ic;
			for (ic=0; ic<ncon; ic++) {
				connectindexes.Push(celref->ConnectIndex(ic)+FirstConnect[imesh]);
			}
		}
		cel->SetConnectIndexes(connectindexes);
	}
}

void TPZBuildMultiphysicsMesh::AppendConnects(TPZCompMesh *cmesh, TPZCompMesh *MFMesh)
{
    long nmeshes = MFMesh->GetNMeshes();
    if(nmeshes<1) DebugStop();
    nmeshes +=1;
    
    //adding connects from cmesh to MFMesh
    long nconnects_old = MFMesh->NConnects();
    long nconnects = nconnects_old + cmesh->NConnects();
    
    MFMesh->ConnectVec().Resize(nconnects);
    MFMesh->Block().SetNBlocks(nconnects);
    
    long counter = nconnects_old;
    long seqnum  = nconnects_old;
    
    long ic;
    long nc = cmesh->NConnects();
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
    
    long iel;
    TPZVec<long> FirstConnect(nmeshes,0);
    long nelem = MFMesh->NElements();
    for (iel = 0; iel < nelem; iel++)
    {
        celorig = MFMesh->ElementVec()[iel];
        cel = dynamic_cast<TPZMultiphysicsElement *> (celorig);
        if (!cel) {
            continue;
        }
        
        long nfirstcon= 0;
        long imesh;
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
		TPZStack<long> connectindexes;
		
        long imesh;
		for (imesh=0; imesh < nmeshes; imesh++) {
			celref = cel->ReferredElement(imesh);
            if (!celref) {
                continue;
            }
			long ncon = celref->NConnects();
			long ic;
			for (ic=0; ic<ncon; ic++) {
				connectindexes.Push(celref->ConnectIndex(ic)+FirstConnect[imesh]);
			}
		}
		cel->SetConnectIndexes(connectindexes);
	}
}


void TPZBuildMultiphysicsMesh::TransferFromMeshes(TPZVec<TPZCompMesh *> &cmeshVec, TPZCompMesh *MFMesh)
{
    long imesh;
    long nmeshes = cmeshVec.size();
    TPZManVector<long> FirstConnectIndex(nmeshes+1,0);
    for (imesh = 0; imesh < nmeshes; imesh++) {
		FirstConnectIndex[imesh+1] = FirstConnectIndex[imesh]+cmeshVec[imesh]->NConnects();
    }
    TPZBlock<STATE> &blockMF = MFMesh->Block();
    for (imesh = 0; imesh < nmeshes; imesh++) {
		long ncon = cmeshVec[imesh]->NConnects();
		TPZBlock<STATE> &block = cmeshVec[imesh]->Block();
		long ic;
		for (ic=0; ic<ncon; ic++) {
			TPZConnect &con = cmeshVec[imesh]->ConnectVec()[ic];
			long seqnum = con.SequenceNumber();
			if(seqnum<0) continue;       /// Whether connect was deleted by previous refined process
			int blsize = block.Size(seqnum);
			TPZConnect &conMF = MFMesh->ConnectVec()[FirstConnectIndex[imesh]+ic];
			long seqnumMF = conMF.SequenceNumber();
			int idf;
			for (idf=0; idf<blsize; idf++) {
				blockMF.Put(seqnumMF, idf, 0, block.Get(seqnum, idf, 0));
			}
		}
	}
}

void TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(TPZVec<TPZCompMesh *> &cmeshVec, TPZCompMesh *MFMesh)
{
    long imesh;
    long nmeshes = cmeshVec.size();
    TPZManVector<long> FirstConnectIndex(nmeshes+1,0);
    for (imesh = 0; imesh < nmeshes; imesh++) {
		FirstConnectIndex[imesh+1] = FirstConnectIndex[imesh]+cmeshVec[imesh]->NConnects();
    }
    TPZBlock<STATE> &blockMF = MFMesh->Block();
    for (imesh = 0; imesh < nmeshes; imesh++) {
		long ncon = cmeshVec[imesh]->NConnects();
		TPZBlock<STATE> &block = cmeshVec[imesh]->Block();
		long ic;
		for (ic=0; ic<ncon; ic++) {
			TPZConnect &con = cmeshVec[imesh]->ConnectVec()[ic];
			long seqnum = con.SequenceNumber();
			if(seqnum<0) continue;       /// Whether connect was deleted by previous refined process
			int blsize = block.Size(seqnum);
			TPZConnect &conMF = MFMesh->ConnectVec()[FirstConnectIndex[imesh]+ic];
            int nelconnected = conMF.NElConnected();
            if (nelconnected == 0) {
                continue;
            }
			long seqnumMF = conMF.SequenceNumber();
			int idf;
			for (idf=0; idf<blsize; idf++) {
				block.Put(seqnum, idf, 0, blockMF.Get(seqnumMF, idf, 0));
			}
		}
	}
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        sout << "Solutions of the referred meshes";
    }
#endif
    for (imesh=0; imesh<nmeshes; imesh++) {
        cmeshVec[imesh]->LoadSolution(cmeshVec[imesh]->Solution());
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
    
    long nelem = cmesh->Reference()->NElements();

	//2- Generate geometric elements (with dimension (meshdim-1)) between the previous elements.
	for (long i=0; i<nelem; ++i) {
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
	for (long i=0; i<nelem; ++i) {
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
        long nelsides = celsides.NElements();
        if(nelsides != 2) 
        {
            DebugStop();
        } 
        for (long lp=0; lp<nelsides; ++lp) {
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
            long index;
            new TPZInterfaceElement(*cmesh,interfaceEl,index,left,right);
            
        }
	}
	
//  std::ofstream arg3("gmeshC.txt");
//	cmesh->Reference()->Print(arg3);

	cmesh->InitializeBlock();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    
    //Set the connect as a Lagrange multiplier  
    long nel = cmesh->NElements();
    for(long i=0; i<nel; i++){
        TPZCompEl *cel = cmesh->ElementVec()[i];
        if(!cel) continue;
        
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

    TPZAdmChunkVector<TPZCompEl *> elvec = cMesh->ElementVec();
    long nel = elvec.NElements();
    for(long el=0; el < nel; el++){
        TPZCompEl * compEl = elvec[el];
        if(!compEl) continue;
        
        if(compEl->IsInterface()){
            compEl->Reference()->ResetReference();
            delete compEl;
        }
    }
    
	TPZVec<long > subindex(0);
	for (int iref = 0; iref < ndiv; iref++) {
		TPZAdmChunkVector<TPZCompEl *> elvec = cMesh->ElementVec();
		long nel = elvec.NElements(); 
		for(long el=0; el < nel; el++){
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
        long ncon = cMesh->NConnects();
        for(long i=0; i<ncon; i++)
        {
            TPZConnect &newnod = cMesh->ConnectVec()[i];
            newnod.SetLagrangeMultiplier(1);
        }
        
        long nel = cMesh->NElements();
        for(long i=0; i<nel; i++){
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

void TPZBuildMultiphysicsMesh::UniformRefineCompEl(TPZCompMesh  *cMesh, long indexEl, bool isLagrMult){
	
	TPZVec<long> subindex; 
	long nel = cMesh->ElementVec().NElements(); 
	for(long el=0; el < nel; el++){
		TPZCompEl * compEl = cMesh->ElementVec()[el];
		if(!compEl) continue;
		long ind = compEl->Index();
		if(ind==indexEl){
            //-------------------------------------------
            TPZStack<TPZCompElSide> neighequal;
            for(int side = compEl->Reference()->NSides()-2; side > compEl->Reference()->NCornerNodes()-1; side--)
            {
                TPZVec<long> subindexneigh;
                long indneigh;
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
        long ncon = cMesh->NConnects();
        for(long i=0; i<ncon; i++)
        {
            TPZConnect &newnod = cMesh->ConnectVec()[i];
            newnod.SetLagrangeMultiplier(1);
        }
        
        long nel = cMesh->NElements();
        for(long i=0; i<nel; i++){
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
void TPZBuildMultiphysicsMesh::ShowShape(TPZVec<TPZCompMesh *> &cmeshVec, TPZCompMesh *MFMesh, TPZAnalysis &analysis, const std::string &filename, TPZVec<long> &equationindices)
{
    TPZStack<std::string> scalnames,vecnames;
    scalnames.Push("State");
    analysis.DefineGraphMesh(analysis.Mesh()->Dimension(), scalnames, vecnames, filename);
    int porder = analysis.Mesh()->GetDefaultOrder();
    
    int neq = equationindices.size();
    TPZFMatrix<STATE> solkeep(analysis.Solution());
    analysis.Solution().Zero();
    for (int ieq = 0; ieq < neq; ieq++) {
        analysis.Solution()(equationindices[ieq],0) = 1.;
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
        long index;
        
        TPZInterpolationSpace *bound;
        MElementType elType = gel->Type(side);
        switch(elType)
        {
            case(EOned)://line
            {
                bound = new TPZCompElHDivBound2<pzshape::TPZShapeLinear>(* intel->Mesh(),gelbound,index);
                int sideorient = intel->GetSideOrient(side);
                TPZCompElHDivBound2<pzshape::TPZShapeLinear> *hdivbound = dynamic_cast< TPZCompElHDivBound2<pzshape::TPZShapeLinear> *>(bound);
                hdivbound->SetSideOrient(sideorient);
                break;
            }
            case(ETriangle)://triangle
            {
                bound = new TPZCompElHDivBound2<pzshape::TPZShapeTriang>(* intel->Mesh(),gelbound,index);
                int sideorient = intel->GetSideOrient(side);
                TPZCompElHDivBound2<pzshape::TPZShapeTriang> *hdivbound = dynamic_cast< TPZCompElHDivBound2<pzshape::TPZShapeTriang> *>(bound);
                hdivbound->SetSideOrient(sideorient);
                break;
            }
            case(EQuadrilateral)://quadrilateral
            {
                bound = new TPZCompElHDivBound2<pzshape::TPZShapeQuad>(* intel->Mesh(),gelbound,index);
                int sideorient = intel->GetSideOrient(side);
                TPZCompElHDivBound2<pzshape::TPZShapeQuad> *hdivbound = dynamic_cast< TPZCompElHDivBound2<pzshape::TPZShapeQuad> *>(bound);
                hdivbound->SetSideOrient(sideorient);
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
        
        long sideconnectindex = intel->ConnectIndex(loccon);
        bound->SetConnectIndex(0, sideconnectindex);
        //bound->Print(std::cout);
        
        TPZCompEl *newMFBound = multiMesh->CreateCompEl(gelbound, index);
        TPZMultiphysicsElement *locMF = dynamic_cast<TPZMultiphysicsElement *>(newMFBound);
        
        locMF->AddElement(bound, 0);
        
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

