/**
 * @file
 * @brief Contains the implementation of the TPZDohrStructMatrix methods.
 */

#include "pzdohrstructmatrix.h"
#include "tpzdohrmatrix.h"
#include "tpzdohrsubstructCondense.h"
#include "tpzdohrprecond.h"
#include "tpznodesetcompute.h"
#include "TPZRenumbering.h"
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
#include "tpzverysparsematrix.h"

#include <sstream>
#include <map>
#include "pzlog.h"

#include "TPZfTime.h"
#include "TPZTimeTemp.h"
#include "TPZVTKGeoMesh.h"
#include <stdlib.h>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("structmatrix.dohrstructmatrix"));
static LoggerPtr loggerasm(Logger::getLogger("structmatrix.dohrstructmatrix.asm"));
#endif

#include "pz_pthread.h"
#include "clock_timer.h"
#include "timing_analysis.h"
#include "arglib.h"
#include "run_stats_table.h"

#ifdef USING_TBB
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
using namespace tbb;
#endif

#include <rcm.h>

#ifdef USING_PAPI
#include <papi.h>

static float stiff_sum = 0;
#endif

/** @brief Return the number of submeshes */
static int64_t NSubMesh(TPZAutoPointer<TPZCompMesh> compmesh);

/** @brief return a pointer to the isub submesh */
static TPZSubCompMesh *SubMesh(TPZAutoPointer<TPZCompMesh> compmesh, int isub);

static void AssembleMatrices(TPZSubCompMesh *submesh, TPZAutoPointer<TPZDohrSubstructCondense<STATE> > substruct, TPZAutoPointer<TPZDohrAssembly<STATE> > dohrassembly,
                             pthread_mutex_t* TestThread);

static void DecomposeBig(TPZAutoPointer<TPZDohrSubstructCondense<STATE> > substruct, int numa_node);
static void DecomposeInternal(TPZAutoPointer<TPZDohrSubstructCondense<STATE> > substruct, int numa_node);

TPZDohrStructMatrix::TPZDohrStructMatrix() :
TPZStructMatrix(), fDohrAssembly(0), fDohrPrecond(0)
{
	PZ_PTHREAD_MUTEX_INIT(&fAccessElement, 0, "TPZDohrStructMatrix::TPZDohrStructMatrix()");
}

TPZDohrStructMatrix::TPZDohrStructMatrix(TPZAutoPointer<TPZCompMesh> cmesh) :
TPZStructMatrix(cmesh), fDohrAssembly(0),
fDohrPrecond(0)
{
	PZ_PTHREAD_MUTEX_INIT(&fAccessElement, 0, "TPZDohrStructMatrix::TPZDohrStructMatrix()");
}

TPZDohrStructMatrix::TPZDohrStructMatrix(const TPZDohrStructMatrix &copy) :
TPZStructMatrix(copy), fDohrAssembly(copy.fDohrAssembly), fDohrPrecond(copy.fDohrPrecond)
{
	PZ_PTHREAD_MUTEX_INIT(&fAccessElement, 0, "TPZDohrStructMatrix::TPZDohrStructMatrix(copy)");
}

TPZDohrStructMatrix::~TPZDohrStructMatrix()
{
	PZ_PTHREAD_MUTEX_DESTROY(&fAccessElement, "TPZDohrStructMatrix::~TPZDohrStructMatrix()");
}


// this will create a DohrMatrix
TPZMatrix<STATE> * TPZDohrStructMatrix::Create()
{
    
	TPZfTime timeforcopute; // init of timer for compute
	fMesh->ComputeNodElCon();
	TPZAutoPointer<TPZDohrAssembly<STATE> > assembly = new TPZDohrAssembly<STATE>;
	fDohrAssembly = assembly;
	
	fMesh->InitializeBlock();
	{
		TPZVec<int64_t> perm,iperm;
		TPZStack<int64_t> elgraph,elgraphindex;
		
		
		int nindep = fMesh->NIndependentConnects();
		fMesh->ComputeElGraph(elgraph,elgraphindex);
		int nel = elgraphindex.NElements()-1;
#ifdef USING_BOOST
		TPZBoostGraph boost(nel,nindep);
		boost.setGType(TPZBoostGraph::KMC);
		boost.SetElementGraph(elgraph, elgraphindex);
		boost.CompressedResequence(perm, iperm);
#else
		TPZSloan sloan(nel,nindep);
		sloan.SetElementGraph(elgraph, elgraphindex);
		sloan.Resequence(perm, iperm);
#endif
		fMesh->Permute(perm);
	}
	int nsub = NSubMesh(fCompMesh);
	int isub;
	
	for(isub=0; isub<nsub; isub++)
	{
		TPZSubCompMesh *submesh = SubMesh(fCompMesh, isub);
#ifdef PZDEBUG
		std::cout << '.'; std::cout.flush();
#endif
		if(!submesh) 
		{
			continue;
		}
		TPZVec<int64_t> perm,iperm;
		TPZStack<int64_t> elgraph,elgraphindex;
		int64_t nindep = submesh->NIndependentConnects();
		submesh->ComputeElGraph(elgraph,elgraphindex);
		int64_t nel = elgraphindex.NElements()-1;
#ifdef USING_BOOST
		TPZBoostGraph boost(nel,nindep);
		boost.setGType(TPZBoostGraph::KMC);
		boost.SetElementGraph(elgraph, elgraphindex);
		boost.CompressedResequence(perm, iperm);
#else
		TPZSloan sloan(nel,nindep);
		sloan.SetElementGraph(elgraph, elgraphindex);
		sloan.Resequence(perm, iperm);
#endif
		
		submesh->Permute(perm);
#ifdef PZDEBUG 
		std::stringstream filename;
		filename << "SubMatrix" << submesh->Index() << ".vtk";
		TPZFMatrix<REAL> fillin(50,50);
		submesh->ComputeFillIn(50, fillin);
		VisualMatrix(fillin, filename.str().c_str());
#endif
	}		
	
	tempo.ft1comput = timeforcopute.ReturnTimeDouble(); //end of time for compute
#ifdef PZDEBUG
	std::cout << tempo.ft1comput << std::endl;
	std::cout << "Identifying corner nodes\n";
	TPZfTime timefornodes; // init of timer
#endif
	
	IdentifyCornerNodes();
    
#ifdef PZDEBUG
	tempo.ft4identcorner = timefornodes.ReturnTimeDouble();
	std::cout << "Total for Identifying Corner Nodes: " << tempo.ft4identcorner << std::endl; // end of timer
#endif
    
	TPZDohrMatrix<STATE,TPZDohrSubstructCondense<STATE> > *dohr = new TPZDohrMatrix<STATE,TPZDohrSubstructCondense<STATE> >(assembly);

	int64_t neq = fMesh->NEquations();
	dohr->Resize(neq,neq);
	// fCornerEqs was initialized during the mesh generation process
	dohr->SetNumCornerEqs(this->fCornerEqs.size());
	
	assembly->fFineEqs.Resize(nsub);
	assembly->fCoarseEqs.Resize(nsub);
	for(isub=0; isub<nsub; isub++)
	{
		TPZSubCompMesh *submesh = SubMesh(fCompMesh, isub);
		if(!submesh) 
		{
			continue;
		}
		TPZAutoPointer<TPZDohrSubstructCondense<STATE> > substruct = new TPZDohrSubstructCondense<STATE>();
		submesh->ComputeNodElCon();
		int64_t neq = ((TPZCompMesh *)submesh)->NEquations();
		//    int neq = substruct->fStiffness->Rows();
		
		substruct->fNEquations = neq;
		
		
		// identify the equation numbers of the submesh
		std::map<int,int> globinv,globaleqs;
		// initialize the globaleqs data structure
		// the global eqs are ordered in the sequence the connects appear
		IdentifyEqNumbers(submesh, globaleqs ,globinv);
		int next = globaleqs.size();
		substruct->fNumExternalEquations = next;
        substruct->fNumInternalEquations = submesh->NumInternalEquations();
		assembly->fFineEqs[isub].Resize(next);
		std::map<int,int>::iterator it;
		int count = 0;
		for(it=globaleqs.begin(); it!=globaleqs.end(); it++)
		{
			assembly->fFineEqs[isub][count++] = it->second; 
		}
        
		
		// initialize the permutations from the mesh enumeration to the external enumeration
		typedef TPZDohrSubstructCondense<STATE>::ENumbering ENumbering;
		typedef std::pair<ENumbering,ENumbering> Numberingpair;
		ENumbering tsub,text,tint;
		tsub = TPZDohrSubstructCondense<STATE>::Submesh;
		text = TPZDohrSubstructCondense<STATE>::ExternalFirst;
		tint = TPZDohrSubstructCondense<STATE>::InternalFirst;
		
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
TPZMatrix<STATE> * TPZDohrStructMatrix::CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface,
                                                       unsigned numthreads_assemble, unsigned numthreads_decompose)
{
    TPZMatrix<STATE> *dohrgeneric = Create();
    Assemble(*dohrgeneric, rhs, guiInterface, numthreads_assemble, numthreads_decompose);
    return dohrgeneric;
}

template<class TVar>
class parallel_assemble_task_t
{
private:
    
    /** We divide the assemble procedure into N work items, which will
     be executed by one or several tasks. The TBB parallel_for
     construct automatically divide the work items in subsets and
     "creates" tasks to execute the work in each subset. Each task
     invokes the operator(blocked_range subset), which will be
     responsible for executing the work items in the subset. */
    template<class TTVar>
    struct work_item_t
    {
        work_item_t (unsigned submesh_idx, const TPZAutoPointer<TPZDohrSubstructCondense<TTVar> >& substruct) :
        fSubMeshIndex(submesh_idx), fSubstruct(substruct) {}
        
        unsigned fSubMeshIndex;
        TPZAutoPointer<TPZDohrSubstructCondense<TTVar> > fSubstruct;
    };
    
    /** Array of work items. */
    std::vector<work_item_t<TVar> > work_items;
    // TODO: Try the cache_aligned_allocator for improved performance.
    //std::vector<work_item_t<TVar>,cache_aligned_allocator<work_item_t<TVar> > > work_items;
    
    /* Pointers to shared data structures. */
    TPZAutoPointer<TPZDohrAssembly<TVar> > fAssembly;
    TPZAutoPointer<TPZCompMesh> fMesh;
    
public:
    
    parallel_assemble_task_t(TPZAutoPointer<TPZDohrAssembly<TVar> > assembly,
                             TPZAutoPointer<TPZCompMesh> mesh) :
    fAssembly(assembly), fMesh(mesh) {}
    
    /** Add a new work item to be list. */
    void push_work_item(unsigned submesh_idx, const TPZAutoPointer<TPZDohrSubstructCondense<TVar> >& substruct)
    {
        work_items.push_back(work_item_t<TVar>(submesh_idx, substruct));
    }
    
    /** Execute work items serially. */
    void run_serial()
    {
        typename std::vector<work_item_t<TVar> >::iterator it = work_items.begin();
        typename std::vector<work_item_t<TVar> >::iterator end = work_items.end();
        
        for (;it != end; it++)
        {
            work_item_t<TVar>& wi = *it;
            TPZSubCompMesh* submesh = SubMesh(fMesh, wi.fSubMeshIndex);
            ::AssembleMatrices(submesh, wi.fSubstruct, fAssembly,NULL);
            ::DecomposeBig(wi.fSubstruct, -2 /* Do not realloc */);
            ::DecomposeInternal(wi.fSubstruct, -2 /* Do not realloc */);
        }
    }
    
#ifdef USING_TBB
    /** Computing operator for the parallel for. */
    void operator()(const blocked_range<size_t>& range) const
    {
        for(size_t i=range.begin(); i!=range.end(); ++i ) {
            const work_item_t<TVar>& wi = work_items[i];
            TPZSubCompMesh* submesh = SubMesh(fMesh, wi.fSubMeshIndex);
            ::AssembleMatrices(submesh, wi.fSubstruct, fAssembly,NULL);
            ::DecomposeBig(wi.fSubstruct,-2 /* Do not realloc */);
            ::DecomposeInternal(wi.fSubstruct,-2 /* Do not realloc */);
        }
    }
    
    /** Execute work items in parallel. */
    void run_parallel_for()
    {
        /* TBB Parallel for. It will split the range into N sub-ranges and
         invoke the operator() for each sub-range.*/
        parallel_for(blocked_range<size_t>(0,work_items.size(), 1 /*IdealGrainSize*/), *this);
    }
#endif
    
};

void TPZDohrStructMatrix::AssembleTBB(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs,
                                      TPZAutoPointer<TPZGuiInterface> guiInterface)
{
    TPZMatrix<STATE> *dohrgeneric = &mat;
    TPZDohrMatrix<STATE,TPZDohrSubstructCondense<STATE> > *dohr =
    dynamic_cast<TPZDohrMatrix<STATE,TPZDohrSubstructCondense<STATE> > *> (dohrgeneric);
    
    const std::list<TPZAutoPointer<TPZDohrSubstructCondense<STATE> > > &sublist = dohr->SubStructures();
    unsigned isub;
    unsigned nsub = NSubMesh(fMesh);
    std::list<TPZAutoPointer<TPZDohrSubstructCondense<STATE> > >::const_iterator it = sublist.begin();
    parallel_assemble_task_t<STATE> parallel_tasks(fDohrAssembly, fMesh);
    
    /* Initialize work items. */
    std::cout << "Assembling " << nsub << " submeshes" << std::endl;
    for (isub=0; isub<nsub ; isub++) {
        TPZSubCompMesh *submesh = SubMesh(fMesh, isub);
        if(!submesh) continue;
        parallel_tasks.push_work_item(isub, *it);
        it++;
    }
    
    /* Run assemble and decompose. */
#ifdef USING_TBB
    parallel_tasks.run_parallel_for();
#else
    parallel_tasks.run_serial();
#endif
    
    /* Post processing. */
    for (isub=0, it=sublist.begin(); it != sublist.end(); it++, isub++) {
        TPZFMatrix<STATE> rhsloc((*it)->fNumExternalEquations,1,0.);
        (*it)->ContributeRhs(rhsloc);
        fDohrAssembly->Assemble(isub,rhsloc,rhs);
    }
    
    dohr->Initialize();
    
    TPZDohrPrecond<STATE,TPZDohrSubstructCondense<STATE> > *precond = new TPZDohrPrecond<STATE,TPZDohrSubstructCondense<STATE> > (*dohr,fDohrAssembly);
    
    precond->Initialize();
    
    fDohrPrecond = precond;
    
    return; // dohrgeneric;
}

template<class T>
struct ThreadDohrmanAssemblyList_ThreadArgs_t
{
    ThreadDohrmanAssemblyList_ThreadArgs_t() : thread_idx(-1), list(NULL) {}
    
    /* Thread index. */
    unsigned thread_idx;
    /* Thread descriptor. */
    pthread_t pthread;
    /* List of items to be assembled. */
    ThreadDohrmanAssemblyList<T>* list;
};

/**
 * @brief Assemble the global system of equations into the matrix which has already been created
 */

    

/* Run statistics. */
/** Jorge comments this code because is missing a file ARGLIB.CPP. 
 This file is in PERFUTIL directory and must to be added to solve linking problems.
*/
RunStatsTable dohr_ass   ("-tpz_dohr_ass", "Raw data table statistics for TPZDohrStructMatrix::Assemble assemble (first)");
RunStatsTable dohr_dec   ("-tpz_dohr_dec", "Raw data table statistics for TPZDohrStructMatrix::Assemble decompose (second)");


void TPZDohrStructMatrix::Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs,
                                   TPZAutoPointer<TPZGuiInterface> guiInterface,
                                   unsigned numthreads_assemble, unsigned numthreads_decompose)
{
#ifdef PERF_ANALYSIS
    ClockTimer timer;
    TimingAnalysis ta;
#endif
    
    TPZMatrix<STATE> *dohrgeneric = &mat;
    TPZDohrMatrix<STATE,TPZDohrSubstructCondense<STATE> > *dohr = dynamic_cast<TPZDohrMatrix<STATE,TPZDohrSubstructCondense<STATE> > *> (dohrgeneric);
    const std::list<TPZAutoPointer<TPZDohrSubstructCondense<STATE> > > &sublist = dohr->SubStructures();
    
    int nsub = NSubMesh(fCompMesh); // mod fMesh
    std::list<TPZAutoPointer<TPZDohrSubstructCondense<STATE> > >::const_iterator it = sublist.begin();
    
    /* Create a list of items to assemble. */
    ThreadDohrmanAssemblyList<STATE> worklist;
    for (int isub=0; isub<nsub ; isub++) {
        TPZSubCompMesh *submesh = SubMesh(fCompMesh, isub); // mod fMesh
        if(!submesh) continue;
        ThreadDohrmanAssembly<STATE> *work = new ThreadDohrmanAssembly<STATE>(fCompMesh,isub,*it,fDohrAssembly); //mod fMesh
        worklist.Append(work);
        it++;
    }
    
    if(guiInterface){
        if(guiInterface->AmIKilled()){
            return ;//0;
        }
    }
    
    // First pass : assembling the matrices
    ThreadDohrmanAssemblyList<STATE> worklistAssemble(worklist);
    std::list<TPZAutoPointer<ThreadDohrmanAssembly<STATE> > >::iterator itwork =
    worklistAssemble.fList.begin();
    while (itwork != worklistAssemble.fList.end()) {
        (*itwork)->fTask = ThreadDohrmanAssembly<STATE>::EComputeMatrix;
        itwork++;
    }
    
    
#ifdef USING_PAPI
    float rtime, ptime, mflops;
    int64_t flpops;
    PAPI_flops ( &rtime, &ptime, &flpops, &mflops );
#endif
    
    dohr_ass.start();
    if (numthreads_assemble == 0) {
        /* Put the main thread to work on all items. */
        ThreadDohrmanAssemblyList_ThreadArgs_t<STATE> targ;
        targ.thread_idx=0;
        targ.list = &worklistAssemble;
        ThreadDohrmanAssemblyList<STATE>::ThreadWork(&targ);
    }
    else {
        /* Threads arguments. */
        std::vector<ThreadDohrmanAssemblyList_ThreadArgs_t<STATE> > args(numthreads_assemble);
        
        /* Assemble multi-threaded */
        for(unsigned itr=0; itr<numthreads_assemble; itr++)
        {
            ThreadDohrmanAssemblyList_ThreadArgs_t<STATE>* targ = &(args[itr]);
            targ->thread_idx=itr;
            targ->list = &worklistAssemble;
            PZ_PTHREAD_CREATE(&targ->pthread, NULL,
                              ThreadDohrmanAssemblyList<STATE>::ThreadWork,
                              targ, __FUNCTION__);
        }
        /* Sync. */
        for(unsigned itr=0; itr<numthreads_assemble; itr++)
        {
            PZ_PTHREAD_JOIN(args[itr].pthread, NULL, __FUNCTION__);
        }
    }
    dohr_ass.stop();
    
#ifdef USING_PAPI
    float ltime;
    PAPI_flops ( &ltime, &ptime, &flpops, &mflops );
    
    printf("Assemble Time: %.2f \t", ltime-rtime);
    printf("Assemble Stiffness : %.2f seconds\n", stiff_sum);
    
#endif
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        int isub = 0;
        for (it=sublist.begin(); it!=sublist.end(); it++) {
            std::stringstream sout;
            sout << "Substructure number " << isub <<std::endl;
            isub++;
           // TPZDohrSubstructCondense<STATE> *ptr = (*it).operator->();
            (*it)->fMatRed->Print("Matred",sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
    }
#endif
    
    // Second  pass : decomposing
    // Philippe: this may be easier to adapt the code for NUMA.
    // Edson: TODO: Measure time again.
    ThreadDohrmanAssemblyList<STATE> worklistDecompose;
    itwork = worklist.fList.begin();
    while (itwork != worklist.fList.end()) {
        TPZAutoPointer<ThreadDohrmanAssembly<STATE> > pt1 = new ThreadDohrmanAssembly<STATE> (*itwork);
        pt1->fTask = ThreadDohrmanAssembly<STATE>::EDecomposeBig;
        worklistDecompose.Append(pt1);
        TPZAutoPointer<ThreadDohrmanAssembly<STATE> > pt2 = new ThreadDohrmanAssembly<STATE>(*itwork);
        pt2->fTask = ThreadDohrmanAssembly<STATE>::EDecomposeInternal;
        worklistDecompose.Append(pt2);
        itwork++;
    }
    
    dohr_dec.start();
    if (numthreads_decompose == 0) {
        /* Compute it sequentialy */
        ThreadDohrmanAssemblyList_ThreadArgs_t<STATE> targ;
        targ.thread_idx = 0;
        targ.list = &worklistDecompose;
        ThreadDohrmanAssemblyList<STATE>::ThreadWork(&targ);
    }
    else {
        /* Threads arguments. */
        std::vector<ThreadDohrmanAssemblyList_ThreadArgs_t<STATE> >
        args(numthreads_decompose);
        for(unsigned itr=0; itr<numthreads_decompose; itr++)
        {
            ThreadDohrmanAssemblyList_ThreadArgs_t<STATE>& targ = args[itr];
            targ.thread_idx=itr;
            targ.list = &worklistDecompose;
            PZ_PTHREAD_CREATE(&targ.pthread, NULL,
                              ThreadDohrmanAssemblyList<STATE>::ThreadWork,
                              &targ, __FUNCTION__);
        }
        for(unsigned itr=0; itr<numthreads_decompose; itr++)
        {
            PZ_PTHREAD_JOIN(args[itr].pthread, NULL, __FUNCTION__);
        }
    }
    dohr_dec.stop();
    
    // Post processing (TODO: check whethe it is time consuming
    int isub;
    for (it=sublist.begin(), isub=0; it != sublist.end(); it++,isub++) {
        TPZFMatrix<STATE> rhsloc((*it)->fNumExternalEquations,1,0.);
        (*it)->ContributeRhs(rhsloc);
        fDohrAssembly->Assemble(isub,rhsloc,rhs);
    }
    
    dohr->SetNumThreads(this->fNumThreads);

    dohr->Initialize();
    TPZDohrPrecond<STATE,TPZDohrSubstructCondense<STATE> > *precond = new TPZDohrPrecond<STATE,TPZDohrSubstructCondense<STATE> > (*dohr,fDohrAssembly);
    precond->Initialize();
    fDohrPrecond = precond;
    
    return; // dohrgeneric;
}

/**
 * @brief Assemble the global right hand side
 */
void TPZDohrStructMatrix::Assemble(TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface)
{
    int rows = fMesh->NEquations();
    rhs.Redim(rows,1);
    TPZDohrPrecond<STATE,TPZDohrSubstructCondense<STATE> > *precond = dynamic_cast<TPZDohrPrecond<STATE,TPZDohrSubstructCondense<STATE> > *>(fDohrPrecond.operator->());
    const std::list<TPZAutoPointer<TPZDohrSubstructCondense<STATE> > > &sublist = precond->Global();
    
    int nsub = NSubMesh(fMesh);
    std::list<TPZAutoPointer<TPZDohrSubstructCondense<STATE> > >::const_iterator it = sublist.begin();
    
    
    int isub;
    for (isub=0; isub<nsub ; isub++) {
        TPZSubCompMesh *submesh = SubMesh(fCompMesh, isub);
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
        TPZFMatrix<STATE> rhsloc((*it)->fNumExternalEquations,1,0.);
        (*it)->ContributeRhs(rhsloc);
        fDohrAssembly->Assemble(isub,rhsloc,rhs);
    }
    
    
}


// identify cornernodes
void TPZDohrStructMatrix::IdentifyCornerNodes()
{
    fCornerEqs.clear();
    TPZStack<int64_t> elementgraph,elementgraphindex;
    TPZStack<int64_t> expelementgraph,expelementgraphindex;
    std::set<int> subelindexes;
    int nelem = fMesh->NElements();
    int iel;
    for (iel=0; iel<nelem ; iel++) {
        TPZSubCompMesh *sub = dynamic_cast<TPZSubCompMesh *> (fMesh->ElementVec()[iel]);
        if (sub) {
            subelindexes.insert(iel);
        }
    }
    // Determine the eligible connect sequence numbers
    std::set<int64_t> cornerconnind;
	std::set<int> cornerconnseq;
    fMesh->BuildCornerConnectList(cornerconnind);
    std::set<int64_t>::iterator it;
    for (it=cornerconnind.begin(); it!=cornerconnind.end(); it++) {
        TPZConnect &c = fMesh->ConnectVec()[*it];
        int seqnum = c.SequenceNumber();
        cornerconnseq.insert(seqnum);
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
#ifdef LOG4CXX
    {
        std::stringstream sout;
        renum.Print(expelementgraph, expelementgraphindex,"Expanded graph",sout);
		if (logger->isDebugEnabled())
		{
			LOGPZ_DEBUG(logger, sout.str())
		}
    }
#endif
    //	renum.SetElementGraph(elementgraph, elementgraphindex);
    std::set<int> othercornereqs;
    renum.CornerEqs(3,nelprev,cornerconnseq,othercornereqs);
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
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
#ifdef PZDEBUG
    std::set<int> cornerseqnums;
#endif
    int nnodes = fMesh->Block().NBlocks();
    int in;
    for (in=0; in<nnodes; in++) {
        if (othercornereqs.find(in) != othercornereqs.end()) {
#ifdef PZDEBUG
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
#ifdef PZDEBUG
    std::cout << "Number cornereqs " << fCornerEqs.size() << std::endl;

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
                int fatherindex = submesh->NodeIndex(connectindex,fMesh);
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
    TPZManVector<int64_t> nodeindices(1,0);
    int ngeo = geonodeindices.NElements();
    int igeo;
    for (igeo=0; igeo<ngeo; igeo++) {
        nodeindices[0] = geonodeindices[igeo];
        int64_t index;
        pointgmesh->CreateGeoElement(EPoint,nodeindices,1,index);
    }
    pointgmesh->BuildConnectivity();
    std::ofstream arquivo("PointMesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(pointgmesh.operator->(),arquivo,true);
#endif
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
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
    int64_t ncon = sub->ConnectVec().NElements();
    // ncon is the number of connects of the subcompmesh
    TPZCompMesh *super = fMesh;
    int64_t ic;
#ifdef LOG4CXX_STOP
    std::stringstream sout;
    sout << "total submesh connects/glob/loc ";
#endif
    for(ic=0; ic<ncon; ic++)
    {
        int64_t glob = sub->NodeIndex(ic,super);
        // continue is the connect is internal
        if(glob == -1) continue;
        int64_t locseq = sub->ConnectVec()[ic].SequenceNumber();
        int64_t globseq = super->ConnectVec()[glob].SequenceNumber();
        int64_t locpos = sub->Block().Position(locseq);
        int64_t globpos = super->Block().Position(globseq);
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
	if (logger->isDebugEnabled())
	{
		LOGPZ_DEBUG(logger, sout.str())
	}
#endif
}

// return the number of submeshes
int64_t NSubMesh(TPZAutoPointer<TPZCompMesh> compmesh)
{
    int64_t nel = compmesh->NElements();
    TPZCompEl *cel;
    int64_t iel, count = 0;
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
    int64_t nel = compmesh->NElements();
    TPZCompEl *cel;
    int64_t iel, count = 0;
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
    TPZVec<int64_t> scatterpermuteblock;
    sub->ComputePermutationInternalFirst(scatterpermuteblock);
    TPZBlock<STATE> destblock = sub->Block();
    TPZBlock<STATE> &origblock = sub->Block();
    int64_t nblocks = origblock.NBlocks();
    if(scatterpermuteblock.NElements() != origblock.NBlocks())
    {
        std::cout << __PRETTY_FUNCTION__ << " something seriously wrong!!!\n";
    }
    int64_t ib;
    for(ib=0; ib<nblocks; ib++)
    {
        destblock.Set(scatterpermuteblock[ib],origblock.Size(ib));
    }
    destblock.Resequence();
    
    int64_t neq = ((TPZCompMesh *)sub)->NEquations();
    scatterpermute.Resize(neq);
    gatherpermute.Resize(neq);
    scatterpermute.Fill(-1);
    gatherpermute.Fill(-1);
    int64_t ncon = sub->ConnectVec().NElements();
#ifdef LOG4CXX_STOP
    std::stringstream sout;
    sout << "internal submesh connects/glob/loc ";
#endif
    int64_t ic;
    for(ic=0; ic<ncon; ic++)
    {
        // skip dependent connects
        TPZConnect &con = sub->ConnectVec()[ic];
        if(con.HasDependency() || con.IsCondensed() ) continue;
        int64_t locseq = sub->ConnectVec()[ic].SequenceNumber();
        // skip unused connects
        if(locseq < 0) continue;
        int destseq = scatterpermuteblock[locseq];
        int64_t locpos = origblock.Position(locseq);
        int64_t destpos = destblock.Position(destseq);
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
    int64_t ieq;
    for(ieq = 0; ieq < neq; ieq++)
    {
        gatherpermute[scatterpermute[ieq]] = ieq;
    }
#ifdef LOG4CXX_STOP
	if (logger->isDebugEnabled())
	{
		LOGPZ_DEBUG(logger, sout.str())
	}
#endif
    
}

// Identify the corner equations associated with a substructure
void TPZDohrStructMatrix::IdentifySubCornerEqs(std::map<int,int> &globaltolocal, TPZVec<int> &cornereqs,
                                               TPZVec<int> &coarseindex)
{
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
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
    int64_t count = 0;
    int64_t localcount = 0;
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
    
    int64_t nel = fMesh->NElements();
    int meshdim = fMesh->Dimension();
    int64_t nnodes = fMesh->NIndependentConnects();
    
    TPZMetis metis(nel,nnodes);
    TPZStack<int64_t> elgraph,elgraphindex;
    fMesh->ComputeElGraph(elgraph,elgraphindex);
    metis.SetElementGraph(elgraph, elgraphindex);
    TPZManVector<int> domain_index(nel,-1);
    metis.Subdivide(nsub, domain_index);
    CorrectNeighbourDomainIndex(fMesh, domain_index);
#ifdef PZDEBUG
    {
        TPZGeoMesh *gmesh = fMesh->Reference();
        int64_t nelgeo = gmesh->NElements();
        TPZVec<int> domaincolor(nelgeo,-999);
        int64_t cel;
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
    CorrectNeighbourDomainIndex(fMesh, domain_index);
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "Geometric mesh and domain indices\n";
        fMesh->Reference()->Print(sout);
        sout << "Domain indices : \n";
        int64_t nel = fMesh->NElements();
        for (int64_t el=0; el<nel; el++) {
            sout << "el " << el << " domain " << domain_index[el] << std::endl;
        }
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    
#ifdef PZDEBUG
    {
        TPZGeoMesh *gmesh = fMesh->Reference();
        int64_t nelgeo = gmesh->NElements();
        TPZVec<int> domaincolor(nelgeo,-999);
        int64_t cel;
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
        int64_t index;
#ifdef PZDEBUG
        std::cout << '^'; std::cout.flush();
#endif
        submeshes[isub] = new TPZSubCompMesh(*fMesh,index);
        if (index < domain_index.NElements()) {
            domain_index[index] = -1;
        }
    }
    int64_t iel;
    for (iel=0; iel<nel; iel++) {
        int domindex = domain_index[iel];
        if (domindex >= 0) {
            TPZCompEl *cel = fMesh->ElementVec()[iel];
            if (!cel) {
                continue;
            }
            submeshes[domindex]->TransferElement(fMesh,iel);
        }
    }
    for (isub = 0; isub<nsub; isub++) {
        int64_t nel = submeshes[isub]->NElements();
        if (nel == 0) {
            delete submeshes[isub];
            submeshes[isub] = 0;
        }
    }
    fMesh->ComputeNodElCon();
    for (isub=0; isub<nsub; isub++) {
        if (submeshes[isub])
        {
            submeshes[isub]->MakeAllInternal();
            submeshes[isub]->PermuteExternalConnects();
#ifdef PZDEBUG
            std::cout << '*'; std::cout.flush();
#endif
        }
    }
    
    fMesh->ComputeNodElCon();
    fMesh->CleanUpUnconnectedNodes();
}

// This is a lengthy process which should run on the remote processor assembling all
void AssembleMatrices(TPZSubCompMesh *submesh, TPZAutoPointer<TPZDohrSubstructCondense<STATE> > substruct, TPZAutoPointer<TPZDohrAssembly<STATE> > dohrassembly,
                      pthread_mutex_t* TestThread)
{
    //	static std::set<int> subindexes;
    //	int index = submesh->Index();
    //	if (subindexes.find(index) != subindexes.end()) {
    //		DebugStop();
    //	}
    //	subindexes.insert(index);
    
    
    {
        typedef TPZDohrSubstructCondense<STATE>::ENumbering ENumbering;
        typedef std::pair<ENumbering,ENumbering> pairnumbering;
        pairnumbering fromsub(TPZDohrSubstructCondense<STATE>::Submesh,TPZDohrSubstructCondense<STATE>::InternalFirst);
        TPZVec<int> &permutescatter = substruct->fPermutationsScatter[fromsub];
        // create a skyline matrix based on the current numbering of the mesh
        // put the stiffness matrix in a TPZMatRed object to facilitate the computation of phi and zi
        TPZSkylineStructMatrix skylstr(submesh);
        skylstr.EquationFilter().Reset();
        
        
        TPZAutoPointer<TPZMatrix<STATE> > Stiffness = skylstr.Create();
        
        
        TPZMatRed<STATE, TPZFMatrix<STATE> > *matredbig = new TPZMatRed<STATE,TPZFMatrix<STATE> >(Stiffness->Rows()+substruct->fCoarseNodes.NElements(),Stiffness->Rows());
        
        
        matredbig->SetK00(Stiffness);
        substruct->fMatRedComplete = matredbig;
        
        
        
        TPZVec<int64_t> permuteconnectscatter;
        
        substruct->fNumInternalEquations = submesh->NumInternalEquations();
        
#ifdef LOG4CXX
        if (logger->isDebugEnabled())
        {
            std::stringstream sout;
            sout << "SubMesh Index = " << submesh->Index() << " Before permutation sequence numbers ";
            int64_t i;
            int64_t ncon = submesh->ConnectVec().NElements();
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
        if (logger->isDebugEnabled())
        {
            std::stringstream sout;
            sout << "SubMesh Index = " << submesh->Index() << " After permutation sequence numbers ";
            int64_t i;
            int64_t ncon = submesh->ConnectVec().NElements();
            for (i=0; i<ncon; i++) {
                sout << i << '|' << submesh->ConnectVec()[i].SequenceNumber() << " ";
            }
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif
#ifdef LOG4CXX
        if (logger->isDebugEnabled())
        {
            std::stringstream sout;
            sout << "SubMesh Index = " << submesh->Index() << "\nComputed scatter vector ";
            sout << permuteconnectscatter;
            sout << "\nStored scatter vector " << permutescatter;
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif
        
        
        // create a "substructure matrix" based on the submesh using a skyline matrix structure as the internal matrix
        TPZMatRedStructMatrix<TPZSkylineStructMatrix,TPZVerySparseMatrix<STATE> > redstruct(submesh);
        TPZMatRed<STATE, TPZVerySparseMatrix<STATE> > *matredptr = dynamic_cast<TPZMatRed<STATE, TPZVerySparseMatrix<STATE> > *>(redstruct.Create());
        //TPZAutoPointer<TPZMatRed<TPZVerySparseMatrix> > matred = matredptr;
        
        // create a structural matrix which will assemble both stiffnesses simultaneously
        // permutescatter will reorder the equations to internal first
        TPZPairStructMatrix pairstructmatrix(submesh,permutescatter);
        
        // reorder the sequence numbering of the connects to reflect the original ordering
        TPZVec<int64_t> invpermuteconnectscatter(permuteconnectscatter.NElements());
        int64_t iel;
        for (iel=0; iel < permuteconnectscatter.NElements(); iel++) {
            invpermuteconnectscatter[permuteconnectscatter[iel]] = iel;
        }
        TPZAutoPointer<TPZMatrix<STATE> > InternalStiffness = matredptr->K00();
        
#ifdef PZDEBUG
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
#ifdef USING_PAPI
        float rtime, ptime, mflops, ltime;
        int64_t flpops;
        
        PAPI_flops ( &rtime, &ptime, &flpops, &mflops );
#endif
        pairstructmatrix.Assemble(Stiffness.operator->(), matredptr, substruct->fLocalLoad);
#ifdef USING_PAPI
        PAPI_flops ( &ltime, &ptime, &flpops, &mflops );
        //printf("Stiff: %.2f \t", ltime-rtime);
        
        stiff_sum += ltime-rtime;
#endif
        // fLocalLoad is in the original ordering of the submesh
        matredbig->SimetrizeMatRed();
        matredptr->SimetrizeMatRed();
        
        substruct->fWeights.Resize(Stiffness->Rows());
        int64_t i;
        for(i=0; i<substruct->fWeights.NElements(); i++)
        {
            substruct->fWeights[i] = Stiffness->GetVal(i,i);
        }
        // Desingularize the matrix without affecting the solution
        int64_t ncoarse = substruct->fCoarseNodes.NElements(), ic;
        int64_t neq = Stiffness->Rows();
        for(ic=0; ic<ncoarse; ic++)
        {
            int coarse = substruct->fCoarseNodes[ic];
            Stiffness->operator()(coarse,coarse) += 10.;
            //Philippe 7/6/2012
            //matredbig->operator()(coarse,coarse) += 10.;
            matredbig->operator()(neq+ic,coarse) = 1.;
            matredbig->operator()(coarse,neq+ic) = 1.;
        }
        //substruct->fStiffness = Stiffness;
        TPZStepSolver<STATE> *InvertedStiffness = new TPZStepSolver<STATE>(Stiffness);
        InvertedStiffness->SetMatrix(Stiffness);
        
        //EBORIN: Uncomment the following line to replace Cholesky by LDLt decomposition
        //#ifdef USE_LDLT_DECOMPOSITION
        
#ifdef USE_LDLT_DECOMPOSITION
        InvertedStiffness->SetDirect(ELDLt);
#else
        InvertedStiffness->SetDirect(ECholesky);
#endif
        matredbig->SetSolver(InvertedStiffness);
        
        
        TPZStepSolver<STATE> *InvertedInternalStiffness = new TPZStepSolver<STATE>(InternalStiffness);
        InvertedInternalStiffness->SetMatrix(InternalStiffness);
#ifdef DUMP_LDLT_MATRICES
        InvertedInternalStiffness->SetDirect(ELDLt);
#else
        InvertedInternalStiffness->SetDirect(ECholesky);
#endif
        matredptr->SetSolver(InvertedInternalStiffness);
        matredptr->SetReduced();
        TPZMatRed<STATE,TPZFMatrix<STATE> > *matredfull = new TPZMatRed<STATE,TPZFMatrix<STATE> >(*matredptr);
        substruct->fMatRed = matredfull;
        
        
    }
}

#ifdef DUMP_LDLT_MATRICES

#include "pzbfilestream.h"
pthread_mutex_t dump_matrix_mutex = PTHREAD_MUTEX_INITIALIZER;
unsigned matrix_unique_id = 0;

void dump_matrix(TPZAutoPointer<TPZMatrix<STATE> > Stiffness)
{
    PZ_PTHREAD_MUTEX_LOCK(&dump_matrix_mutex, "dump_matrix");
    std::cout << "Dump stiffness matrix at DecomposeBig..." << std::endl;
    std::stringstream fname;
    fname << "matrix_" << matrix_unique_id++ << ".bin";
    TPZBFileStream fs;
    fs.OpenWrite(fname.str());
    Stiffness->Write(fs, 0);
    std::cout << "Dump stiffness matrix at DecomposeBig... [Done]" << std::endl;
    PZ_PTHREAD_MUTEX_UNLOCK(&dump_matrix_mutex, "dump_matrix");
}

#endif

//EBORIN: consumes tasks from the ThreadDohrmanAssemblyList list. The tasks
//        are ThreadDohrmanAssembly::AssembleMatrices

#ifdef USING_LIBNUMA
#include<numa.h>
class NUMA_mng_t {
    
public:
    NUMA_mng_t() {
        max_node_id = numa_max_node();
        next_rr_node = 0;
    }
    /** Return the number of nodes on the system.
     *  Nodes are identified from 0 to num_nodes-1. */
    unsigned get_num_nodes() {return (max_node_id+1);}
    /** Return the next node on a round-robin fashion. */
    unsigned get_rr_node_id() {return (next_rr_node++);}
    
private:
    
    unsigned max_node_id;
    /** Next round-robin node. */
    unsigned next_rr_node;
};

NUMA_mng_t NUMA;
#endif

clarg::argBool naa("-naDALora", "NUMA aware Dohrman Assembly List thread work objects re-allocation.", false);
clarg::argInt  naat("-naDALorat", "NUMA aware Dohrman Assembly List thread work objects re-allocation threshold.", 0);

#ifdef USING_LIBNUMA
clarg::argBool nats("-naDALtws", "NUMA aware (node round-robin) Dohrman Assembly List thread work scheduling.", false);
#endif

void DecomposeBig(TPZAutoPointer<TPZDohrSubstructCondense<STATE> > substruct, int numa_node)
{
    TPZAutoPointer<TPZMatRed<STATE,TPZFMatrix<STATE> > > matredbig = substruct->fMatRedComplete;
    TPZAutoPointer<TPZMatrix<STATE> > Stiffness = matredbig->K00();

    if (Stiffness->MemoryFootprint() > naat.get_value()) {
      Stiffness.ReallocForNuma(numa_node);
    }
    
#ifdef USE_LDLT_DECOMPOSITION
    Stiffness->Decompose_LDLt();
#else
    Stiffness->Decompose_Cholesky();
#endif
    
    substruct->Initialize();
}

void DecomposeInternal(TPZAutoPointer<TPZDohrSubstructCondense<STATE> > substruct, int numa_node)
{
    TPZAutoPointer<TPZMatRed<STATE,TPZFMatrix<STATE> > > matred = substruct->fMatRed;
    TPZAutoPointer<TPZMatrix<STATE> > InternalStiffness = matred->K00();
    
    if (InternalStiffness->MemoryFootprint() > naat.get_value()) {
      InternalStiffness.ReallocForNuma(numa_node);
    }
    
#ifdef USE_LDLT_DECOMPOSITION
    InternalStiffness->Decompose_LDLt();
#else
    InternalStiffness->Decompose_Cholesky();
#endif
}

//EComputeMatrix, EDecomposeInternal, EDecomposeBig
template<class TVar>
void ThreadDohrmanAssembly<TVar>::AssembleMatrices(pthread_mutex_t &threadtest, int numa_node)
{
    ThreadDohrmanAssembly *threadData = this;
    TPZSubCompMesh *submesh = SubMesh(threadData->fMesh,threadData->fSubMeshIndex);
    switch (fTask) {
        case EComputeMatrix:
            ::AssembleMatrices(submesh,threadData->fSubstruct,threadData->fAssembly,&threadtest);
            break;
        case EDecomposeInternal:
            DecomposeInternal(threadData->fSubstruct, numa_node);
            break;
        case EDecomposeBig:
            DecomposeBig(threadData->fSubstruct, numa_node);
            break;
        default:
            DebugStop();
            break;
    }
#ifdef LOG4CXX
    if (fTask == EComputeMatrix)
        if (logger->isDebugEnabled())
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

template<class TVar>
ThreadDohrmanAssemblyList<TVar>::ThreadDohrmanAssemblyList()
{
    PZ_PTHREAD_MUTEX_INIT(&fAccessElement,NULL,"ThreadDohrmanAssemblyList::ThreadDohrmanAssemblyList()");
    PZ_PTHREAD_MUTEX_INIT(&fTestThreads,NULL,"ThreadDohrmanAssemblyList::ThreadDohrmanAssemblyList()");
}

template<class TVar>
ThreadDohrmanAssemblyList<TVar>::ThreadDohrmanAssemblyList(ThreadDohrmanAssemblyList<TVar> &cpy) : fList(cpy.fList)
{
    PZ_PTHREAD_MUTEX_INIT(&fAccessElement,NULL,"ThreadDohrmanAssemblyList::ThreadDohrmanAssemblyList()");
    PZ_PTHREAD_MUTEX_INIT(&fTestThreads,NULL,"ThreadDohrmanAssemblyList::ThreadDohrmanAssemblyList()");
}

template<class TVar>
ThreadDohrmanAssemblyList<TVar>::~ThreadDohrmanAssemblyList()
{
	PZ_PTHREAD_MUTEX_DESTROY(&fAccessElement,"ThreadDohrmanAssemblyList::~ThreadDohrmanAssemblyList()");
	PZ_PTHREAD_MUTEX_DESTROY(&fTestThreads,"ThreadDohrmanAssemblyList::~ThreadDohrmanAssemblyList()");
}

template<class TVar>
void ThreadDohrmanAssemblyList<TVar>::Append(TPZAutoPointer<ThreadDohrmanAssembly<TVar> > object)
{
    PZ_PTHREAD_MUTEX_LOCK(&fAccessElement, "ThreadDohrmanAssemblyList::Append()");
    fList.push_back(object);
    PZ_PTHREAD_MUTEX_UNLOCK(&fAccessElement, "ThreadDohrmanAssemblyList::Append()");
}

template<class TVar>
TPZAutoPointer<ThreadDohrmanAssembly<TVar> > ThreadDohrmanAssemblyList<TVar>::NextObject()
{
    TPZAutoPointer<ThreadDohrmanAssembly<TVar> > result;
    PZ_PTHREAD_MUTEX_LOCK(&fAccessElement, "ThreadDohrmanAssemblyList::NextObject()");
    if (fList.begin() != fList.end()) {
        result = *fList.begin();
        fList.pop_front();
    }
    PZ_PTHREAD_MUTEX_UNLOCK(&fAccessElement, "ThreadDohrmanAssemblyList::NextObject()");
    return result;
}

template<class TVar>
void *ThreadDohrmanAssemblyList<TVar>::ThreadWork(void *voidptr)
{
    ThreadDohrmanAssemblyList_ThreadArgs_t<STATE>* args =
    (ThreadDohrmanAssemblyList_ThreadArgs_t<STATE>*) (voidptr);
    
    /* bind thread and newlly allocated memory to node if -naDALtws is set. */
    int node_id = -2 /* Do not realloc */;
    
#ifdef USING_LIBNUMA
    if (nats.was_set()) {
        struct bitmask* nodemask = numa_allocate_nodemask();
        numa_bitmask_clearall(nodemask);
        numa_bitmask_setbit(nodemask,args->thread_idx%NUMA.get_num_nodes());
        numa_bind(nodemask);
        numa_free_nodemask(nodemask);
    }
    if (naa.was_set()) {
        node_id = args->thread_idx%NUMA.get_num_nodes();
    }
#else
    if (naa.was_set()) {
        node_id = -1; /* Realloc */
    }
#endif
    
    TPZAutoPointer<ThreadDohrmanAssembly<TVar> > runner = args->list->NextObject();
    
    while (runner) {
        runner->AssembleMatrices(args->list->fTestThreads,node_id);
        runner = args->list->NextObject();
    }
    
    return 0;
}

// Identify the external connects
void TPZDohrStructMatrix::IdentifyExternalConnectIndexes()
{
    // for each computational element
    std::set<int64_t> connectindexes;
    int64_t iel;
    int64_t nel = fMesh->NElements();
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
            int64_t ncomp = compneigh.NElements();
            int64_t ic;
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
                    int64_t ind = intel->SideConnectIndex(isc,is);
                    connectindexes.insert(ind);
                }
            }
        }
    }
    std::set<int64_t>::iterator it;
    fExternalConnectIndexes.Resize(connectindexes.size());
    int64_t i = 0;
    for (it=connectindexes.begin(); it != connectindexes.end(); it++,i++) {
        fExternalConnectIndexes[i] = *it;
    }
}

// Verifies if the subdomains are connected by sides of connectdimension and separate them if not
// returns the new number of subdomains
int TPZDohrStructMatrix::SeparateUnconnected(TPZVec<int> &domain_index, int nsub, int connectdimension)
{
    std::map<int,int> domain_index_count;
    int64_t iel;
    int64_t nel = fMesh->NElements();
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
                int64_t nneigh = elsidevec.NElements();
                int64_t neigh;
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
            domain_check.erase(mydomainindex);
            domain_check.insert(nsub);
            for (it=gelcluster.begin(); it!=gelcluster.end(); it++) {
                domain_index[(*it)->Index()]=nsub;
            }
            nsub++;
        }
    }
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "Number of elements per domain ";
        std::map<int,int>::iterator it;
        int64_t count = 0;
        for (it=domain_index_count.begin(); it != domain_index_count.end(); it++) {
            if (! (count++ %40)) {
                sout << std::endl;
            }
            sout << it->first << " " << it->second << " " << std::endl;
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
    int64_t nel = fMesh->NElements();
    int64_t mincount = nel/nsub/20;
    // contains for each subdomain the set of neighbouring domains
    TPZVec<std::set<int> > domain_neighbours(nsub);
    // contains for each domain the number of cells within that domain
    std::map<int,int> domain_index_count;
    int64_t iel;
    for (iel=0; iel<nel; iel++) {
        TPZCompEl *cel = fMesh->ElementVec()[iel];
        if (!cel) {
            continue;
        }
        int mydomainindex = domain_index[cel->Index()];
        //        if (mydomainindex == 0) {
        //            std::stringstream sout;
        //            cel->Print(sout);
        //            TPZGeoEl *gel = cel->Reference();
        //            if (gel) {
        //                gel->Print(sout);
        //            }
        //            LOGPZ_DEBUG(logger, sout.str())
        //        }
        domain_index_count[mydomainindex]++;
        TPZGeoEl *gel = cel->Reference();
        if (!gel) {
            continue;
        }
        int nsides = gel->NSides();
        int geldim = gel->Dimension();
        int is;
        for (is=0; is<nsides; is++) {
            int sidedim = gel->SideDimension(is);
            if (sidedim != connectdimension && geldim>=connectdimension) {
                continue;
            }
            if (geldim < connectdimension && is != nsides-1)
            {
                continue;
            }
            TPZGeoElSide gelside(gel,is);
            TPZStack<TPZCompElSide> elsidevec;
            gelside.ConnectedCompElementList(elsidevec, 0, 0);
            int64_t nneigh = elsidevec.NElements();
            int64_t neigh;
            int64_t nneighvalid = 0;
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
                // include the boundary as a ficticious neighbour index
                domain_neighbours[mydomainindex].insert(-1);
            }
        }
    }
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        for (int64_t i=0; i<domain_neighbours.size(); i++) {
            std::set<int>::const_iterator it;
            sout << "Domain index " << i << " neighbours ";
            for (it=domain_neighbours[i].begin(); it != domain_neighbours[i].end(); it++) {
                sout << *it << " ";
            }
            sout << std::endl;
        }
        std::map<int,int>::const_iterator it = domain_index_count.begin();
        while (it != domain_index_count.end()) {
            sout << "Domain index " << it->first << " number of elements " << it->second << std::endl;
            it++;
        }
        LOGPZ_DEBUG(logger, sout.str())
        
    }
#endif
    // compute a destination domain index for each domain (used for clustering domains)
    int isub;
    TPZManVector<int> domain_dest(nsub,-1);
    int64_t count = 0;
    for (isub=0; isub < nsub; isub++)
    {
        // if the subdomain is neighbour to only one subdomain
        // this means that the subdomain is isolated (only boundaries as neighbours) (not treated)
        // or that the domain is embedded in another domain
        if (domain_neighbours[isub].size() == 1 )
        {
            // merge both subdomains
            int target = *(domain_neighbours[isub].begin());
            // target == -1 is not treated here
            if (target == -1) {
                continue;
            }
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
            // the domain has very little elements
            // the domain has at least two neighbouring domains (may include the ficticious -1 domain)
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
    if (logger->isDebugEnabled())
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
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "Number of elements per domain ";
        std::map<int,int>::iterator it;
        int64_t count = 0;
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

void TPZDohrStructMatrix::Write( TPZStream &str, int withclassid ) const
{
    TPZPersistenceManager::WritePointer(fMesh, &str);
    int hasdohrassembly = 0;
    if (fDohrAssembly) {
        hasdohrassembly = 1;
    }
    str.Write(&hasdohrassembly);
    if (hasdohrassembly) {
        TPZPersistenceManager::WritePointer(fDohrAssembly.operator ->(), &str);
    }
    str.Write( fExternalConnectIndexes);
    str.Write(fCornerEqs);
}

void TPZDohrStructMatrix::Read(TPZStream &str, void *context )
{
    SetMesh(TPZAutoPointerDynamicCast<TPZCompMesh>(TPZPersistenceManager::GetAutoPointer(&str)));
    int hasdohrassembly;
    str.Read(&hasdohrassembly);
    if (hasdohrassembly) {
        fDohrAssembly = TPZAutoPointerDynamicCast<TPZDohrAssembly<STATE>>(TPZPersistenceManager::GetAutoPointer(&str));
    }
    str.Read( fExternalConnectIndexes);
    str.Read( fCornerEqs);
}

/** @brief Set the domain index of the lower dimension elements equal to the domain index of their neighbour */
void TPZDohrStructMatrix::CorrectNeighbourDomainIndex(TPZCompMesh *cmesh, TPZVec<int> &domainindex)
{
    int64_t nel = cmesh->NElements();
    TPZAdmChunkVector<TPZCompEl *> &elvec = cmesh->ElementVec();
    bool changed = true;
    while(changed)
    {
        changed = false;
        for (int64_t el=0; el<nel; el++) {
            TPZCompEl *cel = elvec[el];
            if (! cel) {
                continue;
            }
            TPZGeoEl *gel = cel->Reference();
            if (!gel) {
                continue;
            }
            int nsides = gel->NSides();
            TPZGeoElSide neighbour = gel->Neighbour(nsides-1);
            if (neighbour.Element() != gel) {
                TPZCompEl *neighcel = neighbour.Element()->Reference();
                if (! neighcel) {
                    continue;
                }
                int64_t neighindex = neighcel->Index();
                if (domainindex[el] != domainindex[neighindex]) {
                    domainindex[el] = domainindex[neighindex];
                    changed = true;
                }
            }
        }
    }
}

