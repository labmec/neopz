/*
 *  pzdohrstructmatrix.cpp
 *  SubStruct
 *
 *  Created by Philippe Devloo on 28/06/10.
 *  Copyright 2010 UNICAMP. All rights reserved.
 *
 */

#include "pzdohrstructmatrix.h"
#include "tpzdohrmatrix.h"
#include "tpzdohrsubstructCondense.h"
#include "tpzdohrprecond.h"
#include "tpznodesetcompute.h"
#include "pzrenumbering.h"
#include "pzmetis.h"

#include "pzskylstrmatrix.h"
#include "pzmatred.h"
#include "tpzmatredstructmatrix.h"
#include "tpzpairstructmatrix.h"

#include "pzsubcmesh.h"

#include "TPZBoostGraph.h"
#include "pzvisualmatrix.h"
#include "TPZRefPatternTools.h"

#include <sstream>
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("structmatrix.dohrstructmatrix"));
static LoggerPtr loggerasm(Logger::getLogger("structmatrix.dohrstructmatrix.asm"));
#endif

/// return the number of submeshes
static int NSubMesh(TPZAutoPointer<TPZCompMesh> compmesh);

/// return a pointer to the isub submesh
static TPZSubCompMesh *SubMesh(TPZAutoPointer<TPZCompMesh> compmesh, int isub);

// This is a lengthy process which should run on the remote processor
static void InitializeMatrices(TPZSubCompMesh *submesh, TPZAutoPointer<TPZDohrSubstructCondense> substruct, TPZAutoPointer<TPZDohrAssembly> dohrassembly, pthread_mutex_t &testthread);

TPZDohrStructMatrix::TPZDohrStructMatrix(TPZAutoPointer<TPZCompMesh> cmesh) : TPZStructMatrix(cmesh.operator->()), fDohrAssembly(0),
	fDohrPrecond(0), fMesh(cmesh)
{
	pthread_mutex_init(&fAccessElement, 0);
}

TPZDohrStructMatrix::~TPZDohrStructMatrix()
{
	pthread_mutex_destroy(&fAccessElement);
}

// this will create a DohrMatrix
TPZMatrix * TPZDohrStructMatrix::Create()
{

	fMesh->ComputeNodElCon();
	TPZAutoPointer<TPZDohrAssembly> assembly = new TPZDohrAssembly;
	fDohrAssembly = assembly;

	fMesh->InitializeBlock();
	{
		TPZVec<int> perm,iperm;
		TPZStack<int> elgraph,elgraphindex;
		
		
		int nindep = fMesh->NIndependentConnects();
		fMesh->ComputeElGraph(elgraph,elgraphindex);
		int nel = elgraphindex.NElements()-1;
		TPZBoostGraph boost(nel,nindep);
		boost.fGType = TPZBoostGraph::KMCExpensive;
		boost.SetElementGraph(elgraph, elgraphindex);
		boost.Resequence(perm, iperm);
		fMesh->Permute(perm);
	}
	int nsub = NSubMesh(fMesh);
	int isub;
	for(isub=0; isub<nsub; isub++)
	{
		TPZSubCompMesh *submesh = SubMesh(fMesh, isub);
		if(!submesh) 
		{
			continue;
		}
		TPZVec<int> perm,iperm;
		TPZStack<int> elgraph,elgraphindex;
		int nindep = submesh->NIndependentConnects();
		submesh->ComputeElGraph(elgraph,elgraphindex);
		int nel = elgraphindex.NElements()-1;
		TPZBoostGraph boost(nel,nindep);
		boost.fGType = TPZBoostGraph::KMCExpensive;
		boost.SetElementGraph(elgraph, elgraphindex);
		boost.Resequence(perm, iperm);
		submesh->Permute(perm);
#ifdef DEBUG 
		std::stringstream filename;
		filename << "SubMatrix" << submesh->Index() << ".vtk";
		TPZFMatrix fillin(50,50);
		submesh->ComputeFillIn(50, fillin);
		VisualMatrix(fillin, filename.str().c_str());
#endif
	}		
	IdentifyCornerNodes();

	
	TPZDohrMatrix<TPZDohrSubstructCondense> *dohr = new TPZDohrMatrix<TPZDohrSubstructCondense>(assembly);
	dohr->SetNumThreads(this->fNumThreads);

	int neq = fMesh->NEquations();
	dohr->Resize(neq,neq);
	// fCornerEqs was initialized during the mesh generation process
	dohr->SetNumCornerEqs(this->fCornerEqs.size());

	assembly->fFineEqs.Resize(nsub);
	assembly->fCoarseEqs.Resize(nsub);
	for(isub=0; isub<nsub; isub++)
	{
		TPZSubCompMesh *submesh = SubMesh(fMesh, isub);
		if(!submesh) 
		{
			continue;
		}
		TPZAutoPointer<TPZDohrSubstructCondense> substruct = new TPZDohrSubstructCondense();
		submesh->ComputeNodElCon();
		int neq = ((TPZCompMesh *)submesh)->NEquations();
		//    int neq = substruct->fStiffness->Rows();
		
		substruct->fNEquations = neq;
		
		
		// identify the equation numbers of the submesh
		std::map<int,int> globinv,globaleqs;
		// initialize the globaleqs data structure
		// the global eqs are ordered in the sequence the connects appear
		IdentifyEqNumbers(submesh, globaleqs ,globinv);
		int next = globaleqs.size();
		substruct->fNumExternalEquations = next;
		assembly->fFineEqs[isub].Resize(next);
		std::map<int,int>::iterator it;
		int count = 0;
		for(it=globaleqs.begin(); it!=globaleqs.end(); it++)
		{
			assembly->fFineEqs[isub][count++] = it->second; 
		}
		
		// initialize the permutations from the mesh enumeration to the external enumeration
		typedef TPZDohrSubstructCondense::ENumbering ENumbering;
		typedef std::pair<ENumbering,ENumbering> Numberingpair;
		ENumbering tsub,text,tint;
		tsub = TPZDohrSubstructCondense::Submesh;
		text = TPZDohrSubstructCondense::ExternalFirst;
		tint = TPZDohrSubstructCondense::InternalFirst;
		
		TPZVec<int> &toexternal = substruct->fPermutationsScatter[Numberingpair(tsub,text)];
		TPZVec<int> &fromexternal = substruct->fPermutationsScatter[Numberingpair(text,tsub)];
		toexternal.Resize(neq,-1);
		fromexternal.Resize(neq,-1);
		int nel = globaleqs.size();
		count = 0;
		for(it=globaleqs.begin(); it!=globaleqs.end(); it++)
		{
			toexternal[it->first] = count++;
		}
		count = nel++;
		int ieq;
		for(ieq=0; ieq<neq; ieq++)
		{
			if(toexternal[ieq] == -1) toexternal[ieq] = count++;
		}
		for(ieq=0; ieq<neq; ieq++)
		{
			fromexternal[toexternal[ieq]] = ieq;
		}
		
		ComputeInternalEquationPermutation(submesh, substruct->fPermutationsScatter[Numberingpair(tsub,tint)], substruct->fPermutationsScatter[Numberingpair(tint,tsub)]);
		//		IdentifyEqNumbers(submesh, substruct->fGlobalIndex,globinv);
		
		// initialize the fC matrix
		// associate each column of the fC matrix with a coarse index
		IdentifySubCornerEqs(globinv,substruct->fCoarseNodes,assembly->fCoarseEqs[isub]);
		//		int ncoarse = substruct->fCoarseNodes.NElements();
		
		// reorder by internal nodes
		// the fInternalEqs data structure will not be filled if the connects are made internal
		
		// this permutes the nodes of the submesh
		// This is a lengthy process which should run on the remote processor
		dohr->AddSubstruct(substruct);
	}
	return dohr;
}

// this will create a DohrMatrix and compute its matrices
TPZMatrix * TPZDohrStructMatrix::CreateAssemble(TPZFMatrix &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface)
{
	std::cout << "Computing the system of equations for each substructure\n";
	TPZMatrix *dohrgeneric = Create();
	TPZDohrMatrix<TPZDohrSubstructCondense> *dohr = dynamic_cast<TPZDohrMatrix<TPZDohrSubstructCondense> *> (dohrgeneric);
	const std::list<TPZAutoPointer<TPZDohrSubstructCondense> > &sublist = dohr->SubStructures();
	
	int nsub = NSubMesh(fMesh);
	std::list<TPZAutoPointer<TPZDohrSubstructCondense> >::const_iterator it = sublist.begin();
	
	ThreadDohrmanAssemblyList worklist;

	int isub;
	for (isub=0; isub<nsub ; isub++) {
		TPZSubCompMesh *submesh = SubMesh(fMesh, isub);
		if(!submesh) 
		{
			continue;
		}
		ThreadDohrmanAssembly *work = new ThreadDohrmanAssembly(fMesh,isub,*it,fDohrAssembly);
		if (fNumThreads > 0) {
			worklist.Append(work);
		}
		else {
			work->AssembleMatrices((worklist.fTestThreads));
			delete work;
		}
		it++;
	}
	
	
	const int numthreads = fNumThreads;
	std::vector<pthread_t> allthreads(numthreads);
	int itr;
	if(guiInterface){
		if(guiInterface->AmIKilled()){
			return 0;
		}
	}
	for(itr=0; itr<numthreads; itr++)
	{
		pthread_create(&allthreads[itr], NULL,ThreadDohrmanAssemblyList::ThreadWork, &worklist);
	}
	for(itr=0; itr<numthreads; itr++)
	{
		pthread_join(allthreads[itr],NULL);
	}
	
	
	dohr->Initialize();
	TPZDohrPrecond<TPZDohrSubstructCondense> *precond = new TPZDohrPrecond<TPZDohrSubstructCondense> (*dohr,fDohrAssembly);
	precond->Initialize();
	fDohrPrecond = precond;
	return dohrgeneric;
	
}

/// identify cornernodes
void TPZDohrStructMatrix::IdentifyCornerNodes()
{
	fCornerEqs.clear();
	TPZNodesetCompute nodeset;
	TPZStack<int> elementgraph,elementgraphindex;
	//    fCompMesh->ComputeElGraph(elementgraph,elementgraphindex);
	int nindep = fMesh->NIndependentConnects();
	//  int neq = fCMesh->NEquations();
	fMesh->ComputeElGraph(elementgraph,elementgraphindex);
	int nel = elementgraphindex.NElements()-1;
	TPZRenumbering renum(nel,nindep);
    //nodeset.Print(file,elementgraphindex,elementgraph);
	renum.ConvertGraph(elementgraph,elementgraphindex,nodeset.Nodegraph(),nodeset.Nodegraphindex());
	//   cout << "nodegraphindex " << nodeset.Nodegraphindex() << endl;
	//   cout << "nodegraph " << nodeset.Nodegraph() << endl;
	nodeset.AnalyseGraph();
#ifdef LOG4CXX
	{
		std::stringstream str;
		nodeset.Print(str);
		LOGPZ_DEBUG(logger,str.str());
	}
#endif
#ifdef DEBUG
	std::set<int> cornerseqnums;
#endif
	int nnodes = nodeset.Levels().NElements();
//	int maxlev = nodeset.MaxLevel();
	int maxlev = 1;
	int in;
	for(in=0; in<nnodes; in++)
	{
		if(nodeset.Levels()[in] >= maxlev) 
		{
#ifdef DEBUG
			cornerseqnums.insert(in);
#endif
			//      this->fCornerEqs.insert(in);
			//		int seqnum = fCMesh->ConnectVec()[in].SequenceNumber();
			int pos = fMesh->Block().Position(in);
			int size = fMesh->Block().Size(in);
			int ieq;
			for(ieq=0; ieq<size; ieq++)
			{
				this->fCornerEqs.insert(pos+ieq);
			}
		}
	}
#ifdef DEBUG
	std::set<int> connectindices;
	TPZStack<int> geonodeindices;
	int ncon = fMesh->ConnectVec().NElements();
	int ic;
	for (ic=0; ic<ncon; ic++) {
		if (cornerseqnums.find(fMesh->ConnectVec()[ic].SequenceNumber()) != cornerseqnums.end()) {
			connectindices.insert(ic);
		}
	}
	int el;
	int numcel = fMesh->NElements();
	for (el=0; el<numcel; el++) {
		TPZCompEl *cel = fMesh->ElementVec()[el];
		if(!cel) continue;
		TPZSubCompMesh *submesh = dynamic_cast<TPZSubCompMesh *> (cel);
		if(!submesh) continue;
		int elsub;
		int nelsub = submesh->NElements();
		for (elsub=0; elsub<nelsub; elsub++) {
			TPZCompEl *cel = submesh->ElementVec()[elsub];
			if (!cel) {
				continue;
			}
			int ic;
			int nc = cel->NConnects();
			for (ic=0; ic<nc ; ic++) {
				int connectindex = cel->ConnectIndex(ic);
				int fatherindex = submesh->NodeIndex(connectindex,fMesh.operator->());
				if(fatherindex != -1)
				{
					if (connectindices.find(fatherindex) != connectindices.end()) 
					{
					// good one
						TPZGeoEl *gel = cel->Reference();
						int ncornernodes = gel->NCornerNodes();
						if(ic<ncornernodes)
						{
							int nodeindex = gel->NodeIndex(ic);
							geonodeindices.Push(nodeindex);
						}
						connectindices.erase(fatherindex);
					}
				}
			}
		}
	}
	TPZAutoPointer<TPZGeoMesh> pointgmesh = new TPZGeoMesh;
	pointgmesh->NodeVec() = fMesh->Reference()->NodeVec();
	TPZManVector<int> nodeindices(1,0);
	int ngeo = geonodeindices.NElements();
	int igeo;
	for (igeo=0; igeo<ngeo; igeo++) {
		nodeindices[0] = geonodeindices[igeo];
		int index;
		pointgmesh->CreateGeoElement(EPoint,nodeindices,1,index);
	}
	pointgmesh->BuildConnectivity();
	std::ofstream arquivo("PointMesh.vtk");
	TPZRefPatternTools::PrintGMeshVTK(pointgmesh.operator->(),arquivo);
#endif
	
#ifdef LOG4CXX
	{
		std::stringstream str;
		str << "number of corner indices " << fCornerEqs.size() << std::endl;
		str << " corner connect indices ";
		std::set<int>::iterator it;
		for(it=fCornerEqs.begin(); it!=fCornerEqs.end(); it++)
		{
			str << *it << " ";
		}
		LOGPZ_DEBUG(logger,str.str());
	}
#endif
}

/// get the global equation numbers of a substructure (and their inverse)
void TPZDohrStructMatrix::IdentifyEqNumbers(TPZSubCompMesh *sub, std::map<int,int> &global, std::map<int,int> &globinv)
{
	int ncon = sub->ConnectVec().NElements();
	// ncon is the number of connects of the subcompmesh
	TPZCompMesh *super = fMesh.operator->();
	int ic;
#ifdef LOG4CXX_STOP
	std::stringstream sout;
	sout << "total submesh connects/glob/loc ";
#endif
	for(ic=0; ic<ncon; ic++)
	{
		int glob = sub->NodeIndex(ic,super);
		// continue is the connect is internal
		if(glob == -1) continue;
		int locseq = sub->ConnectVec()[ic].SequenceNumber();
		int globseq = super->ConnectVec()[glob].SequenceNumber();
		int locpos = sub->Block().Position(locseq);
		int globpos = super->Block().Position(globseq);
		int locsize = sub->Block().Size(locseq);
		//    int globsize = super->Block().Size(globseq);
		int ieq;
		for(ieq =0; ieq<locsize; ieq++)
		{
#ifdef LOG4CXX_STOP
			sout << ic << "/" << globpos+ieq << "/" << locpos+ieq << " ";
#endif
			global[locpos+ieq] = globpos+ieq;
			globinv[globpos+ieq] = locpos+ieq;
		}
	}
#ifdef LOG4CXX_STOP
	LOGPZ_DEBUG(logger,sout.str())
#endif
}

/// return the number of submeshes
int NSubMesh(TPZAutoPointer<TPZCompMesh> compmesh)
{
	int nel = compmesh->NElements();
	TPZCompEl *cel;
	int iel, count = 0;
	for(iel=0; iel<nel; iel++)
	{
		cel = compmesh->ElementVec()[iel];
		if(!cel) continue;
		TPZSubCompMesh *sub = dynamic_cast<TPZSubCompMesh *>(cel);
		if(sub) count++;
	}
	return count;
}

/// return a pointer to the isub submesh
TPZSubCompMesh *SubMesh(TPZAutoPointer<TPZCompMesh> compmesh, int isub)
{
	int nel = compmesh->NElements();
	TPZCompEl *cel;
	int iel, count = 0;
	for(iel=0; iel<nel; iel++)
	{
		cel = compmesh->ElementVec()[iel];
		if(!cel) continue;
		TPZSubCompMesh *sub = dynamic_cast<TPZSubCompMesh *>(cel);
		if(sub && isub == count) return sub;
		if(sub) count++;
	}
	return NULL;
}

// computes the permutation vectors from the subcompmesh ordening to the "internal first" ordering
// the mesh is modified during this method but is returned to its original state at the end of execution
void TPZDohrStructMatrix::ComputeInternalEquationPermutation(TPZSubCompMesh *sub,
														 TPZVec<int> &scatterpermute, TPZVec<int> &gatherpermute)
{
	// This permutation vector is with respect to the blocks of the mesh
	TPZVec<int> scatterpermuteblock;
	sub->ComputePermutationInternalFirst(scatterpermuteblock);
	TPZBlock destblock = sub->Block();
	TPZBlock &origblock = sub->Block();
	int nblocks = origblock.NBlocks();
	if(scatterpermuteblock.NElements() != origblock.NBlocks())
	{
		std::cout << __PRETTY_FUNCTION__ << " something seriously wrong!!!\n";
	}
	int ib;
	for(ib=0; ib<nblocks; ib++)
	{
		destblock.Set(scatterpermuteblock[ib],origblock.Size(ib));
	}
	destblock.Resequence();
	
	int neq = ((TPZCompMesh *)sub)->NEquations();
	scatterpermute.Resize(neq);
	gatherpermute.Resize(neq);
	scatterpermute.Fill(-1);
	gatherpermute.Fill(-1);
	int ncon = sub->ConnectVec().NElements();
#ifdef LOG4CXX_STOP
	std::stringstream sout;
	sout << "internal submesh connects/glob/loc ";
#endif
	int ic;
	for(ic=0; ic<ncon; ic++)
	{
		// skip dependent connects
		if(sub->ConnectVec()[ic].HasDependency()) continue;
		int locseq = sub->ConnectVec()[ic].SequenceNumber();
		// skip unused connects
		if(locseq < 0) continue;
		int destseq = scatterpermuteblock[locseq];
		int locpos = origblock.Position(locseq);
		int destpos = destblock.Position(destseq);
		int size = origblock.Size(locseq);
		//    int globsize = super->Block().Size(globseq);
		int ieq;
		for(ieq =0; ieq<size; ieq++)
		{
#ifdef LOG4CXX_STOP
			sout << ic << "/" << locpos+ieq << "/" << destpos+ieq << " ";
#endif
			scatterpermute[locpos+ieq] = destpos+ieq;
		}
	}
	int ieq;
	for(ieq = 0; ieq < neq; ieq++)
	{
		gatherpermute[scatterpermute[ieq]] = ieq;
	}
#ifdef LOG4CXX_STOP
	LOGPZ_DEBUG(logger,sout.str())
#endif
	
}

/// Identify the corner equations associated with a substructure
void TPZDohrStructMatrix::IdentifySubCornerEqs(std::map<int,int> &globaltolocal, TPZVec<int> &cornereqs,
										   TPZVec<int> &coarseindex)
{
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Input data for IdentifySubCornerEqs \nglobaltolocal";
		std::map<int,int>::iterator mapit;
		for(mapit = globaltolocal.begin(); mapit != globaltolocal.end(); mapit++)
		{
			sout << " [" << mapit->first << " , " << mapit->second << "] ";
		}
		sout << "\nCorner equations stored in the GenSubStructure data ";
		std::set<int>::iterator setit;
		for(setit = fCornerEqs.begin(); setit != fCornerEqs.end(); setit++)
		{
			sout << *setit << " , ";
		}
		sout << "\ncornereqs " << cornereqs;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	
	
	 cornereqs.Resize(fCornerEqs.size());
	 coarseindex.Resize(fCornerEqs.size());
	 std::set<int>::iterator it;
	 int count = 0;
	 int localcount = 0;
	 for(it = fCornerEqs.begin(); it!= fCornerEqs.end(); it++,count++)
	 {
	 if(globaltolocal.find(*it) != globaltolocal.end())
	 {
	 cornereqs[localcount] = globaltolocal[*it];
	 coarseindex[localcount] = count;
	 localcount++;
	 }
	 }
	 cornereqs.Resize(localcount);
	 coarseindex.Resize(localcount);
}


// partition the mesh in submeshes
void TPZDohrStructMatrix::SubStructure(TPZAutoPointer<TPZCompMesh> cMesh, int nsub )
{
	int nel = cMesh->NElements();
	int nnodes = cMesh->NIndependentConnects();
	TPZMetis metis(nel,nnodes);
	TPZStack<int> elgraph,elgraphindex;
	cMesh->ComputeElGraph(elgraph,elgraphindex);
	metis.SetElementGraph(elgraph, elgraphindex);
	TPZVec<int> domain_index(nel,-1);
	metis.Subdivide(nsub, domain_index);
#ifdef DEBUG 
	{
		TPZGeoMesh *gmesh = cMesh->Reference();
		int nelgeo = gmesh->NElements();
		TPZVec<int> domaincolor(nelgeo,-999);
		int cel;
		for (cel=0; cel<nel; cel++) {
			TPZCompEl *compel = cMesh->ElementVec()[cel];
			if(!compel) continue;
			TPZGeoEl *gel = compel->Reference();
			if (!gel) {
				continue;
			}
			domaincolor[gel->Index()] = domain_index[cel];
		}
		ofstream vtkfile("partition.vtk");
		TPZRefPatternTools::PrintGMeshVTK(gmesh, vtkfile, domaincolor);
	}
#endif
	int isub;
	TPZManVector<TPZSubCompMesh *> submeshes(nsub,0);
	for (isub=0; isub<nsub; isub++) {
		int index;
		submeshes[isub] = new TPZSubCompMesh(cMesh,index);
		if (index < domain_index.NElements()) {
			domain_index[index] = -1;
		}
	}
	int iel;
	for (iel=0; iel<nel; iel++) {
		int domindex = domain_index[iel];
		if (domindex >= 0) {
			TPZCompEl *cel = cMesh->ElementVec()[iel];
			if (!cel) {
				continue;
			}
			submeshes[domindex]->TransferElement(cMesh.operator->(),iel);
		}
	}
	cMesh->ComputeNodElCon();
	for (isub=0; isub<nsub; isub++) {
		submeshes[isub]->MakeAllInternal();
	}
	cMesh->ComputeNodElCon();
	cMesh->CleanUpUnconnectedNodes();
}


// This is a lengthy process which should run on the remote processor
void InitializeMatrices(TPZSubCompMesh *submesh, TPZAutoPointer<TPZDohrSubstructCondense> substruct, TPZAutoPointer<TPZDohrAssembly> dohrassembly,
						pthread_mutex_t &TestThread)
{
//	static std::set<int> subindexes;
//	int index = submesh->Index();
//	if (subindexes.find(index) != subindexes.end()) {
//		DebugStop();
//	}
//	subindexes.insert(index);
	
	typedef TPZDohrSubstructCondense::ENumbering ENumbering;
	typedef std::pair<ENumbering,ENumbering> pairnumbering;
	pairnumbering fromsub(TPZDohrSubstructCondense::Submesh,TPZDohrSubstructCondense::InternalFirst);
	TPZVec<int> &permutescatter = substruct->fPermutationsScatter[fromsub];
	
	// create a skyline matrix based on the current numbering of the mesh
	// put the stiffness matrix in a TPZMatRed object to facilitate the computation of phi and zi
	TPZSkylineStructMatrix skylstr(submesh);
	skylstr.AssembleAllEquations();
	

	TPZAutoPointer<TPZMatrix> Stiffness = skylstr.Create();

	
	TPZMatRed<TPZFMatrix> *matredbig = new TPZMatRed<TPZFMatrix>(Stiffness->Rows()+substruct->fCoarseNodes.NElements(),Stiffness->Rows());


	matredbig->SetK00(Stiffness);
	substruct->fMatRedComplete = matredbig;
	
	
	TPZVec<int> permuteconnectscatter;
	
	substruct->fNumInternalEquations = submesh->NumInternalEquations();
	
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "SubMesh Index = " << submesh->Index() << " Before permutation sequence numbers ";
		int i;
		int ncon = submesh->ConnectVec().NElements();
		for (i=0; i<ncon; i++) {
			sout << i << '|' << submesh->ConnectVec()[i].SequenceNumber() << " ";
		}
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	// change the sequencing of the connects of the mesh, putting the internal connects first
	submesh->PermuteInternalFirst(permuteconnectscatter);

//	pthread_mutex_lock(&TestThread);

#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "SubMesh Index = " << submesh->Index() << " After permutation sequence numbers ";
		int i;
		int ncon = submesh->ConnectVec().NElements();
		for (i=0; i<ncon; i++) {
			sout << i << '|' << submesh->ConnectVec()[i].SequenceNumber() << " ";
		}
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "SubMesh Index = " << submesh->Index() << "\nComputed scatter vector ";
		sout << permuteconnectscatter;
		sout << "\nStored scatter vector " << permutescatter;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	

	// create a "substructure matrix" based on the submesh using a skyline matrix structure as the internal matrix
	TPZMatRedStructMatrix<TPZSkylineStructMatrix,TPZVerySparseMatrix> redstruct(submesh);
	TPZMatRed<TPZVerySparseMatrix> *matredptr = dynamic_cast<TPZMatRed<TPZVerySparseMatrix> *>(redstruct.Create());
	TPZAutoPointer<TPZMatRed<TPZVerySparseMatrix> > matred = matredptr;
	
	// create a structural matrix which will assemble both stiffnesses simultaneously
	TPZPairStructMatrix pairstructmatrix(submesh,permutescatter);
	
	// reorder the sequence numbering of the connects to reflect the original ordering
	TPZVec<int> invpermuteconnectscatter(permuteconnectscatter.NElements());
	int iel;
	for (iel=0; iel < permuteconnectscatter.NElements(); iel++) {
		invpermuteconnectscatter[permuteconnectscatter[iel]] = iel;
	}
	TPZAutoPointer<TPZMatrix> InternalStiffness = matred->K00();
	
#ifdef DEBUG 
	std::stringstream filename;
	filename << "SubMatrixInternal" << submesh->Index() << ".vtk";
	TPZFMatrix fillin(50,50);
	submesh->ComputeFillIn(50, fillin);
	VisualMatrix(fillin, filename.str().c_str());
#endif
	
	submesh->Permute(invpermuteconnectscatter);
	
	

//	pthread_mutex_unlock(&TestThread);
	
	// compute both stiffness matrices simultaneously
	substruct->fLocalLoad.Redim(Stiffness->Rows(),1);
	pairstructmatrix.Assemble(-1, -1, Stiffness.operator->(), matredptr, substruct->fLocalLoad);
	matredbig->Simetrize();
	matredptr->Simetrize();
	
	substruct->fWeights.Resize(Stiffness->Rows());
	int i;
	for(i=0; i<substruct->fWeights.NElements(); i++)
	{
		substruct->fWeights[i] = Stiffness->GetVal(i,i);
	}
	// Desingularize the matrix without affecting the solution
	int ncoarse = substruct->fCoarseNodes.NElements(), ic;
	int neq = Stiffness->Rows();
	for(ic=0; ic<ncoarse; ic++)
	{
		int coarse = substruct->fCoarseNodes[ic];
		Stiffness->operator()(coarse,coarse) += 10.;
		matredbig->operator()(neq+ic,coarse) = 1.;
		matredbig->operator()(coarse,neq+ic) = 1.;
	}
	//substruct->fStiffness = Stiffness;
	TPZStepSolver *InvertedStiffness = new TPZStepSolver(Stiffness);
    InvertedStiffness->SetMatrix(Stiffness);
    InvertedStiffness->SetDirect(ECholesky);
	matredbig->SetSolver(InvertedStiffness);
	
	
	TPZStepSolver *InvertedInternalStiffness = new TPZStepSolver(InternalStiffness);
    InvertedInternalStiffness->SetMatrix(InternalStiffness);
    InvertedInternalStiffness->SetDirect(ECholesky);
	matred->SetSolver(InvertedInternalStiffness);
	matred->SetReduced();
	substruct->fMatRed = matred;
	substruct->Initialize();

}

void ThreadDohrmanAssembly::AssembleMatrices(pthread_mutex_t &threadtest)
{	
	ThreadDohrmanAssembly *threadData = this;
	TPZSubCompMesh *submesh = SubMesh(threadData->fMesh,threadData->fSubMeshIndex);
	InitializeMatrices(submesh, threadData->fSubstruct, threadData->fAssembly,threadtest);
#ifdef LOG4CXX
	{
		std::stringstream sout;
		/*      sout << "Submesh for element " << iel << std::endl;
		 submesh->Print(sout);*/
		sout << "Substructure for submesh " << fSubMeshIndex << std::endl;
		fSubstruct->Print(sout);
		LOGPZ_DEBUG(loggerasm,sout.str())
	}
#endif
	
}

ThreadDohrmanAssemblyList::ThreadDohrmanAssemblyList()
{
	pthread_mutex_init(&fAccessElement,NULL);
	pthread_mutex_init(&fTestThreads,NULL);

}

ThreadDohrmanAssemblyList::~ThreadDohrmanAssemblyList()
{
	pthread_mutex_destroy(&fAccessElement);
	pthread_mutex_destroy(&fTestThreads);
}

void ThreadDohrmanAssemblyList::Append(TPZAutoPointer<ThreadDohrmanAssembly> object)
{
	pthread_mutex_lock(&fAccessElement);
	fList.push_back(object);
	pthread_mutex_unlock(&fAccessElement);
}

TPZAutoPointer<ThreadDohrmanAssembly> ThreadDohrmanAssemblyList::NextObject()
{
	TPZAutoPointer<ThreadDohrmanAssembly> result;
	pthread_mutex_lock(&fAccessElement);
	if (fList.begin() != fList.end()) {
		result = *fList.begin();
		fList.pop_front();
	}
	pthread_mutex_unlock(&fAccessElement);
	return result;
}

void *ThreadDohrmanAssemblyList::ThreadWork(void *voidptr)
{
	ThreadDohrmanAssemblyList *List = (ThreadDohrmanAssemblyList *)(voidptr);
	
	TPZAutoPointer<ThreadDohrmanAssembly> runner = List->NextObject();
	
	while (runner) {
		runner->AssembleMatrices((List->fTestThreads));
		runner = List->NextObject();
	}
	return 0;
}

