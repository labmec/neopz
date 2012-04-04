/**
 * @file
 * @brief Contains the implementation of the TPZDohrStructMatrix methods. 
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
#include "pzfstrmatrix.h"

#include "pzsubcmesh.h"
#include "pzintel.h"

#include "TPZBoostGraph.h"
#include "pzsloan.h"
#include "pzvisualmatrix.h"
#include "TPZRefPatternTools.h"

//#include "TPZfTime.h"

#include <sstream>
#include <map>
#include "pzlog.h"


#include "TPZfTime.h"
#include "TPZTimeTemp.h"
#include "TPZVTKGeoMesh.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("structmatrix.dohrstructmatrix"));
static LoggerPtr loggerasm(Logger::getLogger("structmatrix.dohrstructmatrix.asm"));
#endif

// return the number of submeshes
static int NSubMesh(TPZAutoPointer<TPZCompMesh> compmesh);

// return a pointer to the isub submesh
static TPZSubCompMesh *SubMesh(TPZAutoPointer<TPZCompMesh> compmesh, int isub);

// This is a lengthy process which should run on the remote processor
/*
 static void ComputeMatrices(TPZSubCompMesh *submesh, TPZAutoPointer<TPZDohrSubstructCondense> substruct, TPZAutoPointer<TPZDohrAssembly> dohrassembly, pthread_mutex_t &testthread);
 */
static void DecomposeBig(TPZSubCompMesh *submesh, TPZAutoPointer<TPZDohrSubstructCondense> substruct, TPZAutoPointer<TPZDohrAssembly> dohrassembly, pthread_mutex_t &testthread);
static void DecomposeInternal(TPZSubCompMesh *submesh, TPZAutoPointer<TPZDohrSubstructCondense> substruct, TPZAutoPointer<TPZDohrAssembly> dohrassembly, pthread_mutex_t &testthread);

TPZDohrStructMatrix::TPZDohrStructMatrix(TPZAutoPointer<TPZCompMesh> cmesh, int numthreads_compute, int numthreads_decompose) : 
TPZStructMatrix(cmesh.operator->()), fDohrAssembly(0),
fDohrPrecond(0), fMesh(cmesh)
{
	fNumThreads = numthreads_compute;
	fNumThreadsDecompose = numthreads_decompose;
	
	pthread_mutex_init(&fAccessElement, 0);
}

TPZDohrStructMatrix::TPZDohrStructMatrix(const TPZDohrStructMatrix &copy) : TPZStructMatrix(copy), fNumThreadsDecompose(copy.fNumThreadsDecompose), fDohrAssembly(copy.fDohrAssembly),
fDohrPrecond(copy.fDohrPrecond), fMesh(copy.fMesh)
{
	pthread_mutex_init(&fAccessElement, 0);
}

TPZDohrStructMatrix::~TPZDohrStructMatrix()
{
	pthread_mutex_destroy(&fAccessElement);
}

// this will create a DohrMatrix
TPZMatrix<REAL> * TPZDohrStructMatrix::Create()
{
	
	TPZfTime timeforcopute; // init of timer for compute
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
#ifdef USING_BOOST
		TPZBoostGraph boost(nel,nindep);
		boost.fGType = TPZBoostGraph::KMCExpensive;
		boost.SetElementGraph(elgraph, elgraphindex);
		boost.Resequence(perm, iperm);
#else
		TPZSloan sloan(nel,nindep);
		sloan.SetElementGraph(elgraph, elgraphindex);
		sloan.Resequence(perm, iperm);
#endif
		fMesh->Permute(perm);
	}
	int nsub = NSubMesh(fMesh);
	int isub;
	
	for(isub=0; isub<nsub; isub++)
	{
		TPZSubCompMesh *submesh = SubMesh(fMesh, isub);
		std::cout << '.'; std::cout.flush();
		if(!submesh) 
		{
			continue;
		}
		TPZVec<int> perm,iperm;
		TPZStack<int> elgraph,elgraphindex;
		int nindep = submesh->NIndependentConnects();
		submesh->ComputeElGraph(elgraph,elgraphindex);
		int nel = elgraphindex.NElements()-1;
#ifdef USING_BOOST
		TPZBoostGraph boost(nel,nindep);
		boost.fGType = TPZBoostGraph::KMCExpensive;
		boost.SetElementGraph(elgraph, elgraphindex);
		boost.Resequence(perm, iperm);
#else
		TPZSloan sloan(nel,nindep);
		sloan.SetElementGraph(elgraph, elgraphindex);
		sloan.Resequence(perm, iperm);
#endif
		
		submesh->Permute(perm);
#ifdef DEBUG 
		std::stringstream filename;
		filename << "SubMatrix" << submesh->Index() << ".vtk";
		TPZFMatrix<REAL> fillin(50,50);
		submesh->ComputeFillIn(50, fillin);
		VisualMatrix(fillin, filename.str().c_str());
#endif
	}		
	
	tempo.ft1comput = timeforcopute.ReturnTimeDouble(); //end of time for compute
	std::cout << tempo.ft1comput << std::endl;
	
	std::cout << "Identifying corner nodes\n";
	TPZfTime timefornodes; // init of timer
	
	
	IdentifyCornerNodes();
	
	tempo.ft4identcorner = timefornodes.ReturnTimeDouble();
	std::cout << "Total for Identifying Corner Nodes: " << tempo.ft4identcorner << std::endl; // end of timer
	
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
TPZMatrix<REAL> * TPZDohrStructMatrix::CreateAssemble(TPZFMatrix<REAL> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface)
{
	std::cout << "Computing the system of equations for each substructure\n";
	TPZMatrix<REAL> *dohrgeneric = Create();
	
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
			work->fTask = ThreadDohrmanAssembly::EComputeMatrix;
			work->AssembleMatrices((worklist.fTestThreads));
			work->fTask = ThreadDohrmanAssembly::EDecomposeBig;
			work->AssembleMatrices((worklist.fTestThreads));
			work->fTask = ThreadDohrmanAssembly::EDecomposeInternal;
			work->AssembleMatrices((worklist.fTestThreads));
			delete work;
		}
		it++;
	}
	
	const int numthreads_assemble = fNumThreads;
	std::vector<pthread_t> allthreads_assemble(numthreads_assemble);
	int itr;
	if(guiInterface){
		if(guiInterface->AmIKilled()){
			return 0;
		}
	}
	
	std::cout << "ThreadDohrmanAssembly\n";
	// First pass : assembling the matrices
	ThreadDohrmanAssemblyList worklistAssemble(worklist);
	std::list<TPZAutoPointer<ThreadDohrmanAssembly> >::iterator itwork = worklistAssemble.fList.begin();
	while (itwork != worklistAssemble.fList.end()) {
		(*itwork)->fTask = ThreadDohrmanAssembly::EComputeMatrix;
		itwork++;
	}
	
	TPZfTime timerforassembly; // init of timer
	
	
	for(itr=0; itr<numthreads_assemble; itr++)
	{
		pthread_create(&allthreads_assemble[itr], NULL,ThreadDohrmanAssemblyList::ThreadWork, &worklistAssemble);
	}
	for(itr=0; itr<numthreads_assemble; itr++)
	{
		pthread_join(allthreads_assemble[itr],NULL);
	}
	
	const int numthreads_decompose = this->fNumThreadsDecompose;
	std::vector<pthread_t> allthreads_decompose(numthreads_decompose);
	
	std::cout << "ThreadDohrmanCompute\n";
	// First pass : assembling the matrices
	ThreadDohrmanAssemblyList worklistDecompose;
	itwork = worklist.fList.begin();
	while (itwork != worklist.fList.end()) {
		TPZAutoPointer<ThreadDohrmanAssembly> pt1 = new ThreadDohrmanAssembly(*itwork);
		pt1->fTask = ThreadDohrmanAssembly::EDecomposeBig;
		worklistDecompose.Append(pt1);
		TPZAutoPointer<ThreadDohrmanAssembly> pt2 = new ThreadDohrmanAssembly(*itwork);
		pt2->fTask = ThreadDohrmanAssembly::EDecomposeInternal;
		worklistDecompose.Append(pt2);
		itwork++;
	}
	
	TPZfTime timerfordecompose; // init of timer
	
	
	for(itr=0; itr<numthreads_assemble; itr++)
	{
		pthread_create(&allthreads_decompose[itr], NULL,ThreadDohrmanAssemblyList::ThreadWork, &worklistDecompose);
	}
	for(itr=0; itr<numthreads_decompose; itr++)
	{
		pthread_join(allthreads_decompose[itr],NULL);
	}
	
	tempo.ft5dohrassembly = timerfordecompose.ReturnTimeDouble(); // end of timer
	std::cout << "Time to ThreadDohrmanAssembly" << tempo.ft5dohrassembly << std::endl;
	//	int isub;
	for (it=sublist.begin(), isub=0; it != sublist.end(); it++,isub++) {
		
		// const std::list<TPZAutoPointer<TPZDohrSubstructCondense> > &sublist
		// *it represents the substructure
		TPZFMatrix<REAL> rhsloc((*it)->fNumExternalEquations,1,0.);
		(*it)->ContributeRhs(rhsloc);
		fDohrAssembly->Assemble(isub,rhsloc,rhs);
	}
	
	
	dohr->Initialize();
	TPZDohrPrecond<TPZDohrSubstructCondense> *precond = new TPZDohrPrecond<TPZDohrSubstructCondense> (*dohr,fDohrAssembly);
	precond->Initialize();
	fDohrPrecond = precond;
	return dohrgeneric;
	
}

/**
 * @brief Assemble the global right hand side
 */
void TPZDohrStructMatrix::Assemble(TPZFMatrix<REAL> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface)
{
	
    TPZDohrPrecond<TPZDohrSubstructCondense> *precond = dynamic_cast<TPZDohrPrecond<TPZDohrSubstructCondense> *>(fDohrPrecond.operator->());
	const std::list<TPZAutoPointer<TPZDohrSubstructCondense> > &sublist = precond->Global();//dohr->SubStructures();
	
	int nsub = NSubMesh(fMesh);
	std::list<TPZAutoPointer<TPZDohrSubstructCondense> >::const_iterator it = sublist.begin();
	
	
	int isub;
	for (isub=0; isub<nsub ; isub++) {
		TPZSubCompMesh *submesh = SubMesh(fMesh, isub);
		if(!submesh) 
		{
            DebugStop();
			continue;
		}
        TPZFStructMatrix fullstr(submesh);
        (*it)->fLocalLoad.Zero();
        fullstr.Assemble((*it)->fLocalLoad,guiInterface);
		it++;
	}
    for (it=sublist.begin(), isub=0; it != sublist.end(); it++,isub++) {
		
		// const std::list<TPZAutoPointer<TPZDohrSubstructCondense> > &sublist
		// *it represents the substructure
		TPZFMatrix<REAL> rhsloc((*it)->fNumExternalEquations,1,0.);
		(*it)->ContributeRhs(rhsloc);
		fDohrAssembly->Assemble(isub,rhsloc,rhs);
	}

    
}


// identify cornernodes
void TPZDohrStructMatrix::IdentifyCornerNodes()
{
	fCornerEqs.clear();
	TPZStack<int> elementgraph,elementgraphindex;
	TPZStack<int> expelementgraph,expelementgraphindex;
	std::set<int> subelindexes;
	int nelem = fMesh->NElements();
	int iel;
	for (iel=0; iel<nelem ; iel++) {
		TPZSubCompMesh *sub = dynamic_cast<TPZSubCompMesh *> (fMesh->ElementVec()[iel]);
		if (sub) {
			subelindexes.insert(iel);
		}
	}
	//    fCompMesh->ComputeElGraph(elementgraph,elementgraphindex);
	int nindep = fMesh->NIndependentConnects();
	//  int neq = fCMesh->NEquations();
	fMesh->ComputeElGraph(elementgraph,elementgraphindex);
	int nel = elementgraphindex.NElements()-1;
	// expand the element graph to include a ficticious internal node to all elements
	expelementgraphindex.Push(0);
	int nelprev = nel;
	
	int count = 0;
	for (iel=0; iel<nel; iel++) {
		int nc = elementgraphindex[iel+1]-elementgraphindex[iel];
		if (nc) {
			int index = elementgraphindex[iel];
			int ic;
			for (ic=0; ic<nc; ic++) {
				expelementgraph.Push(0);
				expelementgraph[count++] = elementgraph[index++];
			}
			expelementgraph.Push(0);
			expelementgraph[count++] = nindep;
			nindep++;
		}
		expelementgraphindex.Push(count);
	}
	int next = fExternalConnectIndexes.NElements();
	
	if(next)
		//	if(0)
	{
		TPZManVector<int> externalconnect(nindep,0);
		// add the external connects
		int iext;
		for (iext=0; iext<next; iext++) {
			int extindex = fExternalConnectIndexes[iext];
			int seqnum = fMesh->ConnectVec()[extindex].SequenceNumber();
			if (seqnum >= 0) {
				externalconnect[seqnum] = 1;
			}
		}
		nel = expelementgraphindex.NElements()-1;
		for (iel=0; iel<nel; iel++) {
			bool hasext = false;
			int firstnode = expelementgraphindex[iel];
			int lastnode = expelementgraphindex[iel+1];
			int nodeindex;
			for (nodeindex= firstnode; nodeindex < lastnode; nodeindex++) {
				int node = expelementgraph[nodeindex];
				if (externalconnect[node] ==1) {
					hasext = true;
					break;
				}
			}
			if (hasext) {
				for (nodeindex= firstnode; nodeindex < lastnode; nodeindex++) {
					int node = expelementgraph[nodeindex];
					if (externalconnect[node] ==1) {
						expelementgraph.Push(node);
					}
					expelementgraph.Push(nindep++);
				}
				expelementgraphindex.Push(expelementgraph.NElements());
			}
		}
	}
	// Put a global external element on top of everything
	//	if (next) {
	if (0) {
		count = expelementgraph.NElements();
		int iext;
		for (iext=0; iext<next; iext++) {
			int extindex = fExternalConnectIndexes[iext];
			int seqnum = fMesh->ConnectVec()[extindex].SequenceNumber();
			if (seqnum >= 0) {
				expelementgraph.Push(0);
				expelementgraph[count++] = seqnum;
			}
		}
		expelementgraphindex.Push(count);
	}
	nel = expelementgraphindex.NElements()-1;
	//	nel = elementgraphindex.NElements()-1;
	TPZRenumbering renum(nel,nindep);
	renum.SetElementGraph(expelementgraph, expelementgraphindex);
	//	renum.SetElementGraph(elementgraph, elementgraphindex);
	std::set<int> othercornereqs;
	renum.CornerEqs(3,nelprev,othercornereqs);
#ifdef LOG4CXX
	{
		std::stringstream str;
		int nelem = fMesh->NElements();
		int iel;
		int sub = 0;
		for (iel=0; iel<nelem; iel++) {
			TPZCompEl *cel = fMesh->ElementVec()[iel];
			if (!cel) {
				continue;
			}
			str << "SubCMesh " << sub << std::endl;
			int nc = cel->NConnects();
			int ic;
			for (ic=0; ic<nc; ic++) {
				TPZConnect &con = cel->Connect(ic);
				int seqnum = con.SequenceNumber();
				if (othercornereqs.find(seqnum) != othercornereqs.end()) {
					str << seqnum << " ";
				}
			}
			str << std::endl;
			sub++;
		}
		LOGPZ_DEBUG(logger,str.str());
	}
#endif
#ifdef DEBUG
	std::set<int> cornerseqnums;
#endif
	int nnodes = fMesh->Block().NBlocks();
	int in;
	for (in=0; in<nnodes; in++) {
		if (othercornereqs.find(in) != othercornereqs.end()) {
#ifdef DEBUG
			cornerseqnums.insert(in);
#endif
			int pos = fMesh->Block().Position(in);
			int size = fMesh->Block().Size(in);
			int ieq;
			for(ieq=0; ieq<size; ieq++)
			{
				this->fCornerEqs.insert(pos+ieq);
			}
			
		}
	}
	std::cout << "Number cornereqs " << fCornerEqs.size() << std::endl;
#ifdef DEBUG
	cornerseqnums = othercornereqs;
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
	TPZVTKGeoMesh::PrintGMeshVTK(pointgmesh.operator->(),arquivo,true);
#endif
	
#ifdef LOG4CXX
	{
		std::stringstream str;
		str << "number of corner equations " << fCornerEqs.size() << std::endl;
		int count = 0;
		str << " corner equations ";
		std::set<int>::iterator it;
		for(it=fCornerEqs.begin(); it!=fCornerEqs.end(); it++)
		{
			str << *it << " ";
			count++;
			if (!(count%100)) {
				str << std::endl;
			}
		}
		str << std::endl;
		
		count = 0;
		
		str << "\nnumber of corner block indices after " << othercornereqs.size() << std::endl;
		for(it=othercornereqs.begin(); it!=othercornereqs.end(); it++)
		{
			str << *it << " ";
			count++;
			if (!(count%100)) {
				str << std::endl;
			}
			
		}
		LOGPZ_DEBUG(logger,str.str());
	}
#endif
}

// get the global equation numbers of a substructure (and their inverse)
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

// return the number of submeshes
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

// return a pointer to the isub submesh
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
	TPZBlock<REAL> destblock = sub->Block();
	TPZBlock<REAL> &origblock = sub->Block();
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
        TPZConnect &con = sub->ConnectVec()[ic];
		if(con.HasDependency() || con.IsCondensed() ) continue;
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

// Identify the corner equations associated with a substructure
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
void TPZDohrStructMatrix::SubStructure(int nsub )
{
	
	int nel = fMesh->NElements();
	int meshdim = fMesh->Dimension();
	int nnodes = fMesh->NIndependentConnects();
	
	TPZMetis metis(nel,nnodes);
	TPZStack<int> elgraph,elgraphindex;
	fMesh->ComputeElGraph(elgraph,elgraphindex);
	metis.SetElementGraph(elgraph, elgraphindex);
	TPZVec<int> domain_index(nel,-1);
	metis.Subdivide(nsub, domain_index);
#ifdef DEBUG 
	{
		TPZGeoMesh *gmesh = fMesh->Reference();
		int nelgeo = gmesh->NElements();
		TPZVec<int> domaincolor(nelgeo,-999);
		int cel;
		for (cel=0; cel<nel; cel++) {
			TPZCompEl *compel = fMesh->ElementVec()[cel];
			if(!compel) continue;
			TPZGeoEl *gel = compel->Reference();
			if (!gel) {
				continue;
			}
			domaincolor[gel->Index()] = domain_index[cel];
		}
		std::ofstream vtkfile("partitionbefore.vtk");
		TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, domaincolor);
	}
#endif
	if(meshdim == 3)
	{
		int nsubnew = 0;
		while (nsubnew != nsub)
		{
			nsubnew = SeparateUnconnected(domain_index,nsub,meshdim-1);
			nsub = nsubnew;
		}
		nsub = ClusterIslands(domain_index,nsub,meshdim-1);
	}	
#ifdef DEBUG 
	{
		TPZGeoMesh *gmesh = fMesh->Reference();
		int nelgeo = gmesh->NElements();
		TPZVec<int> domaincolor(nelgeo,-999);
		int cel;
		for (cel=0; cel<nel; cel++) {
			TPZCompEl *compel = fMesh->ElementVec()[cel];
			if(!compel) continue;
			TPZGeoEl *gel = compel->Reference();
			if (!gel) {
				continue;
			}
			domaincolor[gel->Index()] = domain_index[cel];
		}
		std::ofstream vtkfile("partition.vtk");
		TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, domaincolor);
	}
#endif
	int isub;
	TPZManVector<TPZSubCompMesh *> submeshes(nsub,0);
	for (isub=0; isub<nsub; isub++) {
		int index;
		std::cout << '^'; std::cout.flush();
		submeshes[isub] = new TPZSubCompMesh(fMesh,index);
		if (index < domain_index.NElements()) {
			domain_index[index] = -1;
		}
	}
	int iel;
	for (iel=0; iel<nel; iel++) {
		int domindex = domain_index[iel];
		if (domindex >= 0) {
			TPZCompEl *cel = fMesh->ElementVec()[iel];
			if (!cel) {
				continue;
			}
			submeshes[domindex]->TransferElement(fMesh.operator->(),iel);
		}
	}
	fMesh->ComputeNodElCon();
	for (isub=0; isub<nsub; isub++) {
		submeshes[isub]->MakeAllInternal();
		std::cout << '*'; std::cout.flush();
	}
	
	fMesh->ComputeNodElCon();
	fMesh->CleanUpUnconnectedNodes();
}


// This is a lengthy process which should run on the remote processor assembling all
void AssembleMatrices(TPZSubCompMesh *submesh, TPZAutoPointer<TPZDohrSubstructCondense> substruct, TPZAutoPointer<TPZDohrAssembly> dohrassembly,
					  pthread_mutex_t &TestThread)
{
	//	static std::set<int> subindexes;
	//	int index = submesh->Index();
	//	if (subindexes.find(index) != subindexes.end()) {
	//		DebugStop();
	//	}
	//	subindexes.insert(index);
	
	
	{
		typedef TPZDohrSubstructCondense::ENumbering ENumbering;
		typedef std::pair<ENumbering,ENumbering> pairnumbering;
		pairnumbering fromsub(TPZDohrSubstructCondense::Submesh,TPZDohrSubstructCondense::InternalFirst);
		TPZVec<int> &permutescatter = substruct->fPermutationsScatter[fromsub];
		// create a skyline matrix based on the current numbering of the mesh
		// put the stiffness matrix in a TPZMatRed object to facilitate the computation of phi and zi
		TPZSkylineStructMatrix skylstr(submesh);
		skylstr.AssembleAllEquations();
		
		
		TPZAutoPointer<TPZMatrix<REAL> > Stiffness = skylstr.Create();
		
		
		TPZMatRed<REAL, TPZFMatrix<REAL> > *matredbig = new TPZMatRed<REAL,TPZFMatrix<REAL> >(Stiffness->Rows()+substruct->fCoarseNodes.NElements(),Stiffness->Rows());
		
		
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
		TPZMatRedStructMatrix<TPZSkylineStructMatrix,TPZVerySparseMatrix<REAL> > redstruct(submesh);
		TPZMatRed<REAL, TPZVerySparseMatrix<REAL> > *matredptr = dynamic_cast<TPZMatRed<REAL, TPZVerySparseMatrix<REAL> > *>(redstruct.Create());
		//TPZAutoPointer<TPZMatRed<TPZVerySparseMatrix> > matred = matredptr;
		
		// create a structural matrix which will assemble both stiffnesses simultaneously
		// permutescatter will reorder the equations to internal first
		TPZPairStructMatrix pairstructmatrix(submesh,permutescatter);
		
		// reorder the sequence numbering of the connects to reflect the original ordering
		TPZVec<int> invpermuteconnectscatter(permuteconnectscatter.NElements());
		int iel;
		for (iel=0; iel < permuteconnectscatter.NElements(); iel++) {
			invpermuteconnectscatter[permuteconnectscatter[iel]] = iel;
		}
		TPZAutoPointer<TPZMatrix<REAL> > InternalStiffness = matredptr->K00();
		
#ifdef DEBUG 
		std::stringstream filename;
		filename << "SubMatrixInternal" << submesh->Index() << ".vtk";
		TPZFMatrix<REAL> fillin(50,50);
		submesh->ComputeFillIn(50, fillin);
		VisualMatrix(fillin, filename.str().c_str());
#endif
		
		// put the equation back in the optimized ordering for all equations (original ordering)
		submesh->Permute(invpermuteconnectscatter);
		
		
		
		//	pthread_mutex_unlock(&TestThread);
		
		// compute both stiffness matrices simultaneously
		substruct->fLocalLoad.Redim(Stiffness->Rows(),1);
		pairstructmatrix.Assemble(-1, -1, Stiffness.operator->(), matredptr, substruct->fLocalLoad);
		// fLocalLoad is in the original ordering of the submesh
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
		TPZStepSolver<REAL> *InvertedStiffness = new TPZStepSolver<REAL>(Stiffness);
		InvertedStiffness->SetMatrix(Stiffness);
		InvertedStiffness->SetDirect(ECholesky);
		matredbig->SetSolver(InvertedStiffness);
		
		
		TPZStepSolver<REAL> *InvertedInternalStiffness = new TPZStepSolver<REAL>(InternalStiffness);
		InvertedInternalStiffness->SetMatrix(InternalStiffness);
		InvertedInternalStiffness->SetDirect(ECholesky);
		matredptr->SetSolver(InvertedInternalStiffness);
		matredptr->SetReduced();
		TPZMatRed<REAL,TPZFMatrix<REAL> > *matredfull = new TPZMatRed<REAL,TPZFMatrix<REAL> >(*matredptr);
		substruct->fMatRed = matredfull;
		
	}
}
void DecomposeBig(TPZSubCompMesh *submesh, TPZAutoPointer<TPZDohrSubstructCondense> substruct, TPZAutoPointer<TPZDohrAssembly> dohrassembly,
				  pthread_mutex_t &TestThread)
{
	
	{
		TPZAutoPointer<TPZMatRed<REAL,TPZFMatrix<REAL> > > matredbig = substruct->fMatRedComplete;
		
		TPZAutoPointer<TPZMatrix<REAL> > Stiffness = matredbig->K00();
		
		Stiffness->Decompose_Cholesky();
		substruct->Initialize();
	}
}
void DecomposeInternal(TPZSubCompMesh *submesh, TPZAutoPointer<TPZDohrSubstructCondense> substruct, TPZAutoPointer<TPZDohrAssembly> dohrassembly,
					   pthread_mutex_t &TestThread)
{
	
	{
		TPZAutoPointer<TPZMatRed<REAL,TPZFMatrix<REAL> > > matred = substruct->fMatRed;
		TPZAutoPointer<TPZMatrix<REAL> > InternalStiffness = matred->K00();
		InternalStiffness->Decompose_Cholesky();
	}
	
}

//EComputeMatrix, EDecomposeInternal, EDecomposeBig
void ThreadDohrmanAssembly::AssembleMatrices(pthread_mutex_t &threadtest)
{	
	ThreadDohrmanAssembly *threadData = this;
	TPZSubCompMesh *submesh = SubMesh(threadData->fMesh,threadData->fSubMeshIndex);
	switch (fTask) {
		case EComputeMatrix:
			::AssembleMatrices(submesh,threadData->fSubstruct,threadData->fAssembly,threadtest);
			break;
		case EDecomposeInternal:
			DecomposeInternal(submesh,threadData->fSubstruct,threadData->fAssembly,threadtest);
			break;
		case EDecomposeBig:
			DecomposeBig(submesh,threadData->fSubstruct,threadData->fAssembly,threadtest);
			break;
		default:
			DebugStop();
			break;
	}
#ifdef LOG4CXX
	if (fTask == EComputeMatrix)
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
		std::cout << '+'; std::cout.flush();
		runner->AssembleMatrices((List->fTestThreads));
		runner = List->NextObject();
	}
	return 0;
}

// Identify the external connects
void TPZDohrStructMatrix::IdentifyExternalConnectIndexes()
{
	// for each computational element
	std::set<int> connectindexes;
	int iel;
	int nel = fMesh->NElements();
	for (iel=0; iel<nel; iel++) {
		// if it has a neighbour along its interior, skip
		TPZCompEl *cel = fMesh->ElementVec()[iel];
		if (!cel) {
			continue;
		}
		TPZGeoEl *gel = cel->Reference();
		if (!gel) {
			continue;
		}
		int is;
		int ns = gel->NSides();
		int dim = gel->Dimension();
		TPZStack<TPZCompElSide> compneigh;
		
		// if there is a neighbour along the side of dimension dim skip
		TPZGeoElSide gelside(gel,ns-1);
		gelside.ConnectedCompElementList(compneigh,0,0);
		if (compneigh.NElements()) {
			continue;
		}
		TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
		if (!intel) {
			continue;
		}
		// loop over the sides of dimension dim-1
		for (is=0; is<ns; is++) 
		{
			// if there is a neighbour of dimension >= dim skip
			// the side connects are external
			TPZGeoElSide gelside(gel,is);
			if (gelside.Dimension() != dim-1) {
				continue;
			}
			compneigh.Resize(0);
			gelside.ConnectedCompElementList(compneigh, 0, 0);
			int ncomp = compneigh.NElements();
			int ic;
			for (ic=0; ic<ncomp; ic++) {
				TPZCompElSide celside = compneigh[ic];
				TPZGeoElSide gelside = celside.Reference();
				if (gelside.Element()->Dimension() == dim) {
					break;
				}
			}
			// if no neighbour has dimension dim
			if (ic == ncomp) {
				int nsconnect = intel->NSideConnects(is);
				int isc;
				for (isc=0; isc<nsconnect; isc++) {
					int ind = intel->SideConnectIndex(isc,is);
					connectindexes.insert(ind);
				}
			}
		}
	}
	std::set<int>::iterator it;
	fExternalConnectIndexes.Resize(connectindexes.size());
	int i = 0;
	for (it=connectindexes.begin(); it != connectindexes.end(); it++,i++) {
		fExternalConnectIndexes[i] = *it;
	}
}

// Verifies if the subdomains are connected by sides of connectdimension and separate them if not
// returns the new number of subdomains
int TPZDohrStructMatrix::SeparateUnconnected(TPZVec<int> &domain_index, int nsub, int connectdimension)
{
	
	std::map<int,int> domain_index_count;
	int iel;
	int nel = fMesh->NElements();
	for (iel=0; iel<nel; iel++) {
		TPZCompEl *cel = fMesh->ElementVec()[iel];
		if (!cel) {
			continue;
		}
		int mydomainindex = domain_index[cel->Index()];
		domain_index_count[mydomainindex]++;
	}
	std::set<int> domain_check;
	
	for (iel=0; iel<nel; iel++) {
		TPZCompEl *cel = fMesh->ElementVec()[iel];
		if (!cel) {
			continue;
		}
		int mydomainindex = domain_index[cel->Index()];
		if (domain_check.find(mydomainindex) != domain_check.end()) {
			continue;
		}
		domain_check.insert(mydomainindex);
		TPZGeoEl *gel = cel->Reference();
		if (!gel) {
			continue;
		}
		TPZStack<TPZGeoEl *> gelstack;
		gelstack.Push(gel);
		std::set<TPZCompEl *> gelcluster;
		while (gelstack.NElements()) 
		{
			TPZGeoEl *gel = gelstack.Pop();
			if (gelcluster.find(gel->Reference()) != gelcluster.end()) {
				continue;
			}
			int beforesize = gelcluster.size();
			gelcluster.insert(gel->Reference());
			int checksize = gelcluster.size();
			if (checksize == beforesize) {
				DebugStop();
			}
			
			int nsides = gel->NSides();
			int is;
			for (is=0; is<nsides; is++) {
				int sidedim = gel->SideDimension(is);
				if (sidedim != connectdimension) {
					continue;
				}
				TPZGeoElSide gelside(gel,is);
				TPZStack<TPZCompElSide> elsidevec;
				gelside.ConnectedCompElementList(elsidevec, 0, 0);
				int nneigh = elsidevec.NElements();
				int neigh;
				for (neigh = 0; neigh <nneigh; neigh++) {
					TPZCompElSide celside = elsidevec[neigh];
					TPZCompEl *celloc = celside.Element();
					if (domain_index[celloc->Index()] != mydomainindex) {
						continue;
					}
					if (gelcluster.find(celloc) == gelcluster.end()) {
						gelstack.Push(celloc->Reference());
					}
				}
			}
		}
		if (gelcluster.size() != (std::set<TPZCompEl *>::size_type)domain_index_count[mydomainindex]) {
			if (gelcluster.size() > (std::set<TPZCompEl *>::size_type)domain_index_count[mydomainindex]) {
				DebugStop();
			}
			domain_index_count[mydomainindex] -= gelcluster.size();
			domain_index_count[nsub] = gelcluster.size();
			std::set<TPZCompEl *>::iterator it;
			for (it=gelcluster.begin(); it!=gelcluster.end(); it++) {
				domain_index[(*it)->Index()]=nsub;
			}
			nsub++;
		}
	}
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Number of elements per domain ";
		std::map<int,int>::iterator it;
		int count = 0;
		for (it=domain_index_count.begin(); it != domain_index_count.end(); it++) {
			if (! (count++ %40)) {
				sout << std::endl;
			}
			sout << it->first << " " << it->second << " ";
		}
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	return nsub;
}

// Eliminate subdomains who are embedded in other subdomains
// returns the number of subdomains
int TPZDohrStructMatrix::ClusterIslands(TPZVec<int> &domain_index,int nsub,int connectdimension)
{
	int meshdim = fMesh->Dimension();
	int nel = fMesh->NElements();
	int mincount = nel/nsub/20;
	TPZVec<std::set<int> > domain_neighbours(nsub);
	std::map<int,int> domain_index_count;
	int iel;
	for (iel=0; iel<nel; iel++) {
		TPZCompEl *cel = fMesh->ElementVec()[iel];
		if (!cel) {
			continue;
		}
		int mydomainindex = domain_index[cel->Index()];
		domain_index_count[mydomainindex]++;
		TPZGeoEl *gel = cel->Reference();
		if (!gel) {
			continue;
		}
		int nsides = gel->NSides();
		int is;
		for (is=0; is<nsides; is++) {
			int sidedim = gel->SideDimension(is);
			if (sidedim != connectdimension) {
				continue;
			}
			TPZGeoElSide gelside(gel,is);
			TPZStack<TPZCompElSide> elsidevec;
			gelside.ConnectedCompElementList(elsidevec, 0, 0);
			int nneigh = elsidevec.NElements();
			int neigh;
			int nneighvalid = 0;
			for (neigh = 0; neigh <nneigh; neigh++) {
				TPZCompElSide celside = elsidevec[neigh];
				TPZCompEl *celloc = celside.Element();
				TPZGeoEl *gelloc = celloc->Reference();
				if (gelloc->Dimension() != meshdim) {
					continue;
				}
				nneighvalid++;
				int celdomain = domain_index[celloc->Index()];
				if (celdomain != mydomainindex) 
				{
					domain_neighbours[mydomainindex].insert(celdomain);
				}
			}
			if (nneighvalid == 0) 
			{
				// include a boundary as a neighbour
				domain_neighbours[mydomainindex].insert(-1);
			}
		}
	}
	int isub;
	TPZManVector<int> domain_dest(nsub,-1);
	int count = 0;
	for (isub=0; isub < nsub; isub++) 
	{
		if (domain_neighbours[isub].size() == 1 ) 
		{
			int target = *(domain_neighbours[isub].begin());
			if (domain_dest[target] == -1 && domain_dest[isub] == -1) 
			{
				domain_dest[isub] = count;
				domain_dest[target] = count;
				count++;
			}
			else if (domain_dest[target] == -1)
			{
				domain_dest[target] = domain_dest[isub];
			}
			else 
			{
				domain_dest[isub] = domain_dest[target];
			}
			
		}
		else if(domain_dest[isub] == -1 && domain_index_count[isub] < mincount)
		{
			std::map<int,int> sizeDomain;
			std::set<int>::iterator it;
			for (it = domain_neighbours[isub].begin(); it != domain_neighbours[isub].end(); it++) {
				if (*it != -1) {
					sizeDomain[domain_index_count[isub]] = *it;
				}
			}
			int domaintargetindex = sizeDomain.rbegin()->second;
			int destdomainindexcount = domain_index_count[domaintargetindex];
			int domainshrinkcount = domain_index_count[isub];
			domain_index_count[domaintargetindex] = destdomainindexcount+domainshrinkcount;
			domain_index_count[isub] = 0;
			if(domain_dest[domaintargetindex] == -1)
			{
				domain_dest[domaintargetindex] = count;
				domain_dest[isub] = count;
				count++;
			}
			else {
				domain_dest[isub] = domain_dest[domaintargetindex];
			}
			
		}
		else if (domain_dest[isub] == -1) 
		{
			domain_dest[isub] = count++;
		}
	}
#ifdef LOG4CXX
	{
		std::stringstream sout;
		int isub;
		for (isub=0; isub < nsub; isub++) {
			sout << "isub = " << isub << " number of neighbours " << domain_neighbours[isub].size() << " domains ";
			std::set<int>::iterator it;
			for (it = domain_neighbours[isub].begin(); it != domain_neighbours[isub].end(); it++) {
				sout << *it << " ";
			}
			sout << std::endl;
		}
		sout << "Destination domain " << domain_dest << std::endl;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	int domsize = domain_index.NElements();
	int d;
	for (d=0; d<domsize; d++) {
		domain_index[d] = domain_dest[domain_index[d]];
	}
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Number of elements per domain ";
		std::map<int,int>::iterator it;
		int count = 0;
		for (it=domain_index_count.begin(); it != domain_index_count.end(); it++) {
			if (! (count++ %40)) {
				sout << std::endl;
			}
			sout << it->first << " " << it->second << " ";
		}
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	
	return count;
}