/**
 * @file
 * @brief Contains the implementation of the TPZPairStructMatrix methods. 
 */

#include "tpzpairstructmatrix.h"
#include "TPZSSpStructMatrix.h"
#include "TPZTimer.h"
#include "TPZElementMatrixT.h"
#include "pzcompel.h"
#include "pzsubcmesh.h"
#include "TPZLinearAnalysis.h"
#include "TPZMaterial.h"

#include <thread>

using namespace std;

#include "pzlog.h"


#ifdef USING_TBB
#include <tbb/tbb.h>
using namespace tbb;
#endif

#ifdef PZ_LOG
static TPZLogger logger("pz.strmatrix.tpzpairstructmatrix");
static TPZLogger loggerel("pz.strmatrix.element");
#endif

int TPZPairStructMatrix::gNumThreads = 0;

TPZPairStructMatrix::TPZPairStructMatrix(TPZCompMesh *mesh, TPZVec<int> &permutescatter) : fStrMatrix( new TPZSSpStructMatrix<STATE>(mesh))
{
	fPermuteScatter = permutescatter;
#ifdef PERF_DEBUG
	std::cout << "fNumThreads e "<< fNumThreads << " gNumThreads esta setando para " << gNumThreads << endl; 
#endif
    fStrMatrix->SetNumThreads(gNumThreads);

}

#ifdef USING_TBB1

struct PipeItem_t
{
  TPZCompEl* el;
  TPZElementMatrix ek;
  TPZElementMatrix ef;
};

/** @brief First stage: select these TPZCompEl elements to be processed. */
class StageOne_t : public tbb::filter
{
  unsigned current_iel;
  unsigned nelem;
  TPZAdmChunkVector<TPZCompEl *> &elementvec;
  TPZCompMesh& mesh;

  std::vector<PipeItem_t> items;
  unsigned current_item;

  std::set<int>& fMaterialIds;

  bool ShouldCompute(int matid)
  {
    return fMaterialIds.size()==0 || fMaterialIds.find(matid) != fMaterialIds.end();
  }

public:

  StageOne_t(unsigned tokens, TPZCompMesh& _mesh, std::set<int>& materialIds) : 
    tbb::filter(/*is serial*/ true), current_iel(0), 
    nelem(_mesh.NElements()), elementvec(_mesh.ElementVec()), mesh(_mesh),
    items(tokens), current_item(0), fMaterialIds(materialIds)
  {}
  
  unsigned n_tokens() {return items.size();}

  void* operator()(void*)
  {
    TPZCompEl* el = NULL;
    bool found = false;
    
    while(current_iel < nelem) 
    {
      el = elementvec[current_iel++];

      if (!el) continue;

      if(fMaterialIds.size() == 0) {
        found = true;
        break;
      }

      //cout << "DEBUG: fMaterialIds.size() != 0 ????" << endl;

      TPZMaterial * mat = el->Material();
      TPZSubCompMesh *submesh = dynamic_cast<TPZSubCompMesh *> (el);
      if(!mat)
	  {
        if(!submesh) {
            continue;
        }
        else if(submesh->NeedsComputing(fMaterialIds) == false) {
          continue;
		}
      }
      else 
	  {
        int matid = mat->Id();
        if(this->ShouldCompute(matid) == false) continue;
      }
      found = true;
      break;
    }

    if (found) 
    {
      /* Ok, we found an item to process! */
      PipeItem_t& item = items[current_item];
      current_item = (current_item+1) % items.size();
      item.el = el;
      //cout << "Stage 1: Producing item: " << (current_iel-1) << endl;
      return &item;
    }
    else 
    {
      //cout << "Stage 1: Done" << endl;
      return NULL; // It's over!
    }
  }
};

/** @brief Second stage: compute the ek and ef matrices. */
class StageTwo_t: public tbb::filter 
{
  TPZCompMesh *mesh;
  int mineq;
  int maxeq;

public:
  StageTwo_t(TPZCompMesh* _mesh, int _mineq, int _maxeq) : 
    tbb::filter(/*is serial*/ false), mesh(_mesh), mineq(_mineq), maxeq(_maxeq)
  {}
  
  void* operator()(void* item_p)
  {
    PipeItem_t& item = *static_cast<PipeItem_t*>(item_p);

    item.ek.Reset(mesh, TPZElementMatrix::EK);
    item.ef.Reset(mesh, TPZElementMatrix::EF);

    item.el->CalcStiff(item.ek,item.ef);

    //cout << "Stage Two: processing item: " << item_p << endl;

    if (item.el->HasDependency())
    {
      // the element has dependent nodes
      item.ek.ApplyConstraints();
      item.ef.ApplyConstraints();
    }
      
    item.ek.ComputeDestinationIndices();

    if(mineq != -1 && maxeq != -1)
      TPZStructMatrix::FilterEquations(item.ek.fSourceIndex, item.ek.fDestinationIndex, 
                                       mineq, maxeq);

    return item_p;
  }
};

/** @brief Third stage: compute the first matrix and the rhs. */
template<class TVar>
class StageThree_t: public tbb::filter 
{
  TPZMatrix<TVar>& fGlobMatrix1;
  TPZFMatrix<TVar>& fGlobRhs;

public:
  StageThree_t(TPZMatrix<TVar>& _fGlobMatrix1, TPZFMatrix<TVar>& _fGlobRhs) : 
    tbb::filter(/*is serial*/ true), fGlobMatrix1(_fGlobMatrix1), fGlobRhs(_fGlobRhs)
  {}
  
  void* operator()(void* item_p)
  {
    PipeItem_t& item = *static_cast<PipeItem_t*>(item_p);
    TPZElementMatrix& ek = item.ek;
    TPZElementMatrix& ef = item.ef;

    // Assemble the matrix
    if(!ek.HasDependency())
    {
      fGlobMatrix1.AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
      fGlobRhs.AddFel(ef.fMat,ek.fSourceIndex,ek.fDestinationIndex);				
    }
    else
    {
      fGlobMatrix1.AddKel(ek.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
      fGlobRhs.AddFel(ef.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);				
    }

    return item_p;
  }
};

/** @brief Fourth stage: compute the second matrix. */
template<class TVar>
class StageFour_t: public tbb::filter 
{
  TPZMatrix<TVar>& fGlobMatrix2;
  const TPZVec<int>& fPermuteScatter;

  void PermuteScatter(TPZVec<int64_t> &index)
  {
	int nel = index.NElements();
	int iel;
	for(iel = 0; iel<nel; iel++)
		index[iel] = fPermuteScatter[index[iel]];
  }

public:
  StageFour_t(TPZMatrix<TVar>& _fGlobMatrix2, const TPZVec<int>& _fPermuteScatter) : 
    tbb::filter(/*is serial*/ true), fGlobMatrix2(_fGlobMatrix2), 
    fPermuteScatter(_fPermuteScatter)
  {}
  
  void* operator()(void* item_p)
  {
    PipeItem_t& item = *static_cast<PipeItem_t*>(item_p);
    TPZElementMatrix& ek = item.ek;

    PermuteScatter(ek.fDestinationIndex);
	
    // Assemble the matrix
    if(!ek.HasDependency())
      fGlobMatrix2.AddKel(ek.fMat, ek.fSourceIndex, ek.fDestinationIndex);
    else
      fGlobMatrix2.AddKel(ek.fConstrMat, ek.fSourceIndex, ek.fDestinationIndex);

    return NULL;
  }
};

#endif // USING TBB

/**
 * @brief Iterações do laço podem ser executadas em paralelo. \n
 * - 1: Utilizar um parallel_for: (simples, mas não explora paralelismo entre threadassembly 1 e 2. \n
 * - 2: Utilizar árvore de tarefas: Podemos explorar o paralelismo entre threadassembly 1 e 2, mas \n
 * O código pode ficar grande.
 */
template<class TVar>
class parallel_assemble_task_t
{
private:

  void PermuteScatter(TPZVec<int> &index)
  {
	int nel = index.NElements();
	int iel;
	for(iel = 0; iel<nel; iel++)
    {
      index[iel] = fPermuteScatter[index[iel]];
    }
  }

  /** Array of work items. */
  std::vector<TPZCompEl*> work_items;
  // TODO: Try the cache_aligned_allocator for improved performance.
  //std::vector<work_item_t<TVar>,cache_aligned_allocator<work_item_t<TVar> > > work_items;

  /* Pointers to shared data structures. */

  //TPZAutoPointer<TPZDohrAssembly<TVar> > fAssembly;
  //TPZAutoPointer<TPZCompMesh> fMesh;

public:

  /** Add a new work item to be list. */
  void push_work_item(TPZCompEl *el)
  {
    work_items.push_back(el);
  }

#ifdef USING_TBB1

protected:

  typedef spin_mutex MatrixMutex_t;
  MatrixMutex_t Matrix1Mutex;
  MatrixMutex_t Matrix2Mutex;

public:

  /** Execute work items in parallel. */
  void run_parallel_for()
  {
    /* TBB Parallel for. It will split the range into N sub-ranges and
       invoke the operator() for each sub-range.*/
    parallel_for(blocked_range<size_t>(0,work_items.size(), 1 /*IdealGrainSize*/), *this); 
  }

  /** Computing operator for the parallel for. */
  void operator()(const blocked_range<size_t>& range)
  { 
    MatrixMutex_t::scoped_lock lock;
    TPZElementMatrix ek(&mesh, TPZElementMatrix::EK);
    TPZElementMatrix ef(&mesh, TPZElementMatrix::EF);

    for(size_t i=range.begin(); i!=range.end(); ++i ) {

      TPZCompEl* el = work_items[i];

      el->CalcStiff(ek,ef);

      bool has_dependency = el->HasDependency();
      if(has_dependency) 
      {
        // the element has dependent nodes
        ek.ApplyConstraints();
        ef.ApplyConstraints();
      }
      
      ek.ComputeDestinationIndices();

      if(mineq != -1 && maxeq != -1)
        TPZStructMatrix::FilterEquations(ek.fSourceIndex,ek.fDestinationIndex,mineq,maxeq);
      
      // ThreadAssembly 1 -- Lock 1
      lock.acquire(Matrix1Mutex);
      if(!has_dependency) 
      {
        first->AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
        rhs.AddFel(ef.fMat,ek.fSourceIndex,ek.fDestinationIndex);
      } 
      else 
      {
        first->AddKel(ek.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
        rhs.AddFel(ef.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
      }
      lock.release();
      
      // ThreadAssembly 2 -- Lock 2
      // FIXME. Precisa de um lock
      lock.acquire(Matrix2Mutex);
      PermuteScatter(ek.fDestinationIndex);
      second->AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);      
      lock.release();
    }
  } 

#endif // USING TBB

  /** Execute work items serially. */
  void run_serial()
  {
    TPZElementMatrixT<STATE> ek(fStrMatrix->Mesh(), TPZElementMatrix::EK);
    TPZElementMatrixT<STATE> ef(fStrMatrix->Mesh(), TPZElementMatrix::EF);

    std::vector<TPZCompEl*>::iterator it = work_items.begin();
    for(; it != work_items.end(); it++)
    {   
      TPZCompEl* el = *it;

      el->CalcStiff(ek,ef);

      bool has_dependency = el->HasDependency();
      if(has_dependency) 
      {
        // the element has dependent nodes
        ek.ApplyConstraints();
        ef.ApplyConstraints();
      }
      
      ek.ComputeDestinationIndices();

    fStrMatrix->FilterEquations(ek.fSourceIndex,ek.fDestinationIndex);
      
      // ThreadAssembly 1
      if(!has_dependency) 
      {
        first->AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
        rhs.AddFel(ef.fMat,ek.fSourceIndex,ek.fDestinationIndex);
      } 
      else 
      {
        first->AddKel(ek.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
        rhs.AddFel(ef.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
      }
      
      // ThreadAssembly 2
      PermuteScatter(ek.fDestinationIndex);
      second->AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);      
    }
  }

  TPZMatrix<STATE>*  first;
  TPZMatrix<STATE>*  second;
  TPZFMatrix<STATE>& rhs;
  TPZVec<int>&       fPermuteScatter;
    TPZStructMatrix *fStrMatrix;

  parallel_assemble_task_t(TPZStructMatrix *strmat, int _mineq, int _maxeq, TPZMatrix<STATE> *_first, 
                           TPZMatrix<STATE> *_second, TPZFMatrix<STATE> &_rhs, TPZVec<int>& permuteScatter) :
    first(_first), second(_second), 
    rhs(_rhs), fPermuteScatter(permuteScatter), fStrMatrix(strmat)
  {
      fStrMatrix->SetEquationRange(_mineq, _maxeq);
  }

};

#include"arglib.h"
clarg::argInt tbb_pair_pipe_tokens("-tbb_pair_ntokens", "# tokens during assemble TPZPairStructMatrix (level 2) using TBB", 128);

void TPZPairStructMatrix::TBBAssemble(TPZMatrix<STATE> *first, 
                                      TPZMatrix<STATE> *second, TPZFMatrix<STATE> &rhs)
{
#ifndef USING_TBB1

  cerr << "TPZPairStructMatrix::TBBAssemble() invoked, but TBB "
       << "was not linked. Executing SerialAssemble!" << endl;
  return SerialAssemble(first, second, rhs);

#else // if USING_TBB
  int iel;
  TPZCompMesh &mesh = *fMesh;
  int nelem = mesh.NElements();
  TPZAdmChunkVector<TPZCompEl *> &elementvec = mesh.ElementVec();

  // Create the pipeline 
  tbb::pipeline pipeline; 
  
  // Create file-reading writing stage and add it to the pipeline 
  StageOne_t          filter1(tbb_pair_pipe_tokens.get_value() /* Number of tokens on the fly */, mesh, fMaterialIds);
  StageTwo_t          filter2(fMesh, mineq, maxeq);
  StageThree_t<STATE> filter3(*first, rhs);
  StageFour_t<STATE>  filter4(*second, fPermuteScatter);

  pipeline.add_filter(filter1);
  pipeline.add_filter(filter2);
  pipeline.add_filter(filter3);
  pipeline.add_filter(filter4);
  
  // Run the pipeline 
  pipeline.run(filter1.n_tokens() /* Max tokens on the fly: TODO - Tune this parameter */ ); 
  
  pipeline.clear(); 
#endif // USING_TBB
}

void TPZPairStructMatrix::SerialAssemble(TPZMatrix<STATE> *first, TPZMatrix<STATE> *second, TPZFMatrix<STATE> &rhs)
{
	int iel;
	TPZCompMesh &mesh = *fStrMatrix->Mesh();
	int nelem = mesh.NElements();
	TPZElementMatrixT<STATE> ek(&mesh, TPZElementMatrix::EK),ef(&mesh, TPZElementMatrix::EF);
    
		TPZTimer calcstiff("Computing the stiffness matrices");
	TPZTimer assemble("Assembling the stiffness matrices");
	TPZAdmChunkVector<TPZCompEl *> &elementvec = mesh.ElementVec();
	
	for(iel=0; iel < nelem; iel++) {
		TPZCompEl *el = elementvec[iel];
		if(!el) continue;
		calcstiff.start();
		
		el->CalcStiff(ek,ef);
		
		calcstiff.stop();
		assemble.start();
		
		if(!el->HasDependency()) {
			ek.ComputeDestinationIndices();
            fStrMatrix->FilterEquations(ek.fSourceIndex,ek.fDestinationIndex);
			
			first->AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
			rhs.AddFel(ef.fMat,ek.fSourceIndex,ek.fDestinationIndex);
			PermuteScatter(ek.fDestinationIndex);
			second->AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
#ifdef PZ_LOG
			if(loggerel.isDebugEnabled())
			{
				std::stringstream sout;
				ek.fMat.Print("Element stiffness matrix",sout);
				ef.fMat.Print("Element right hand side", sout);
				LOGPZ_DEBUG(loggerel,sout.str())
			}
#endif
		} else {
			// the element has dependent nodes
			ek.ApplyConstraints();
			ef.ApplyConstraints();
			ek.ComputeDestinationIndices();
            //FIXME: (Edson) - O operador é  && ou ||
            fStrMatrix->FilterEquations(ek.fSourceIndex,ek.fDestinationIndex);
			first->AddKel(ek.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
			rhs.AddFel(ef.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
			PermuteScatter(ek.fDestinationIndex);
            //FIXME: (Edson) - Não deveria utilizar ek.fConstrMat?
			second->AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
		}
		
		assemble.stop();
	}//fim for iel
#ifdef PZ_LOG
	{
		std::stringstream sout;
		sout << "Number of equations " << first->Rows() << std::endl;
		sout << calcstiff.processName() << " " << calcstiff << std::endl;
		sout << assemble.processName() << " " << assemble;
		/*    stiffness.Print("Matriz de Rigidez: ",sout);
		 rhs.Print("Vetor de Carga: ",sout);*/
		if (logger.isDebugEnabled())
		{
			LOGPZ_DEBUG(logger, sout.str().c_str());
		}
	}
#endif
	
}

void TPZPairStructMatrix::PermuteScatter(TPZVec<int64_t> &index)
{
	int nel = index.NElements();
	int iel;
	for(iel = 0; iel<nel; iel++)
	{
		index[iel] = ((int64_t)fPermuteScatter[index[iel]]);
	}
}
void TPZPairStructMatrix::PermuteScatter(TPZVec<int> &index)
{
	int nel = index.NElements();
	int iel;
	for(iel = 0; iel<nel; iel++)
	{
		index[iel] = fPermuteScatter[index[iel]];
	}
}

void TPZPairStructMatrix::Assemble(TPZMatrix<STATE> *first, TPZMatrix<STATE> *second, TPZFMatrix<STATE> &rhs)
{
#ifdef PERF_DEBUG
  std::cout << "TPZPairStructMatrix::Assemble = Assembling the system of equations with " << fNumThreads
	    << " threads (TPZPairStructMatrix::gNumThreads = " << TPZPairStructMatrix::gNumThreads  << ")\n";
#endif

    /// "Fixme!!!"   "Fixme!!!"   "Fixme!!!"   "Fixme!!!"
  /// @note Find a better way to select among TBB, std::thread or serial execution! in this TPZPairStructMatrix::Assemble

    int numthreads = fStrMatrix->GetNumThreads();
  if (numthreads < 0)
    TBBAssemble(first, second, rhs);
  else if(numthreads > 0)
    MultiThread_Assemble(first, second, rhs);
  else 
    SerialAssemble(first,second,rhs);
}

void TPZPairStructMatrix::MultiThread_Assemble(TPZMatrix<STATE> *first, TPZMatrix<STATE> *second, TPZFMatrix<STATE> &rhs)
{
    ThreadData threaddata(fStrMatrix.operator->(),*first,*second,rhs);
	threaddata.fPermuteScatter = fPermuteScatter;
	const int numthreads = fStrMatrix->GetNumThreads();
	std::vector<thread> allthreads(numthreads+1);
	int itr;
	for(itr=0; itr<numthreads; itr++)
	{
        allthreads[itr] = thread(ThreadData::ThreadWork, &threaddata);
	}

  // assemble the first matrix
  allthreads[itr] = thread(ThreadData::ThreadAssembly1, &threaddata);

	// assemble the second matrix
	ThreadData::ThreadAssembly2(&threaddata);
	
	for(itr=0; itr<numthreads+1; itr++)
	{
        allthreads[itr].join();
	}
}

TPZPairStructMatrix::ThreadData::ThreadData(TPZStructMatrix *strmat, TPZMatrix<STATE> &mat1, TPZMatrix<STATE> &mat2, 
											TPZFMatrix<STATE> &rhs)
: fStrMatrix(strmat), 
fGlobMatrix1(&mat1), fGlobMatrix2(&mat2), fGlobRhs(&rhs),fNextElement(0),fAccessElement()
{	
}

TPZPairStructMatrix::ThreadData::~ThreadData()
{
}

void *TPZPairStructMatrix::ThreadData::ThreadWork(void *datavoid)
{
	ThreadData *data = (ThreadData *) datavoid;
	// compute the next element (this method is threadsafe)
	int iel = data->NextElement();
	TPZCompMesh *cmesh = data->fStrMatrix->Mesh();
	int nel = cmesh->NElements();
	while(iel < nel)
	{
		
      TPZAutoPointer<TPZElementMatrixT<STATE>> ek =
        new TPZElementMatrixT<STATE>(cmesh,TPZElementMatrix::EK);
      TPZAutoPointer<TPZElementMatrixT<STATE>> ef =
        new TPZElementMatrixT<STATE>(cmesh,TPZElementMatrix::EF);
		
		TPZCompEl *el = cmesh->ElementVec()[iel];
		el->CalcStiff(ek,ef);
		
		
		if(!el->HasDependency()) {
			ek->ComputeDestinationIndices();
			
            data->fStrMatrix->FilterEquations(ek->fSourceIndex,ek->fDestinationIndex);
#ifdef PZ_LOG
			if(loggerel.isDebugEnabled())
			{
				std::stringstream sout;
				ek->fMat.Print("Element stiffness matrix",sout);
				ef->fMat.Print("Element right hand side", sout);
				LOGPZ_DEBUG(loggerel,sout.str())
			}
#endif
		} else {
			// the element has dependent nodes
			ek->ApplyConstraints();
			ef->ApplyConstraints();
			ek->ComputeDestinationIndices();
            data->fStrMatrix->FilterEquations(ek->fSourceIndex,ek->fDestinationIndex);

		}
		
		
		// put the elementmatrices on the stack to be assembled (threadsafe)
		data->ComputedElementMatrix(iel,ek,ef);
		// compute the next element (this method is threadsafe)
		iel = data->NextElement();
	}
	data->fAssembly1.Post();
	data->fAssembly2.Post();
	return 0;
}

void TPZPairStructMatrix::ThreadData::PermuteScatter(TPZVec<int> &index)
{
	int nel = index.NElements();
	int iel;
	for(iel = 0; iel<nel; iel++)
	{
		index[iel] = fPermuteScatter[index[iel]];
	}
}
void TPZPairStructMatrix::ThreadData::PermuteScatter(TPZVec<int64_t> &index)
{
	int nel = index.NElements();
	int iel;
	for(iel = 0; iel<nel; iel++)
	{
		index[iel] = ((int64_t)fPermuteScatter[index[iel]]);
	}
}


// The function which will compute the assembly
void *TPZPairStructMatrix::ThreadData::ThreadAssembly1(void *threaddata)
{
	ThreadData *data = (ThreadData *) threaddata;
	TPZCompMesh *cmesh = data->fStrMatrix->Mesh();
	int nel = cmesh->NElements();
    unique_lock<std::mutex> lock(data->fAccessElement);
	int nextel = data->fNextElement;
	int numprocessed = data->fProcessed1.size();
	while(nextel < nel || numprocessed)
	{
		std::set<int>::iterator itprocess;
		bool keeplooking = false;
		if(data->fSubmitted1.size() && data->fProcessed1.size())
		{
			auto itavail = data->fSubmitted1.begin();
			itprocess = data->fProcessed1.begin();
			if(itavail->first == *itprocess)
			{
				// make sure we come back to look for one more element
				keeplooking = true;
				// Get a hold of the data
				data->fProcessed1.erase(itprocess);
				TPZAutoPointer<TPZElementMatrixT<STATE>> ek = itavail->second.first;
				TPZAutoPointer<TPZElementMatrixT<STATE>> ef = itavail->second.second;
				data->fSubmitted1.erase(itavail);
#ifdef PZ_LOG
				int iel = *itprocess;
				std::stringstream sout;
				sout << "Assembling element " << iel;
				if (logger.isDebugEnabled())
				{
					LOGPZ_DEBUG(logger, sout.str())
				}
#endif
				// Release the mutex
                lock.unlock();
				// Assemble the matrix
				if(!ek->HasDependency())
				{
					data->fGlobMatrix1->AddKel(ek->fMat,ek->fSourceIndex,ek->fDestinationIndex);
					data->fGlobRhs->AddFel(ef->fMat,ek->fSourceIndex,ek->fDestinationIndex);				
				}
				else
				{
					data->fGlobMatrix1->AddKel(ek->fConstrMat,ek->fSourceIndex,ek->fDestinationIndex);
					data->fGlobRhs->AddFel(ef->fConstrMat,ek->fSourceIndex,ek->fDestinationIndex);				
				}
				// acquire the mutex
                lock.lock();
			}
		}
		if(!keeplooking)
		{
            lock.unlock();
#ifdef PZ_LOG
            if(logger.isDebugEnabled()) LOGPZ_DEBUG(logger, "Going to sleep within assembly")
#endif
			// wait for a signal
			data->fAssembly1.Wait();
#ifdef PZ_LOG
            if(logger.isDebugEnabled()) LOGPZ_DEBUG(logger, "Waking up for assembly")
#endif
            lock.lock();
		}
		nextel = data->fNextElement;
		numprocessed = data->fProcessed1.size();
		
	}
	return 0;
}		

// The function which will compute the assembly
void *TPZPairStructMatrix::ThreadData::ThreadAssembly2(void *threaddata)
{
	ThreadData *data = (ThreadData *) threaddata;
	TPZCompMesh *cmesh = data->fStrMatrix->Mesh();
	int nel = cmesh->NElements();
    unique_lock<std::mutex> lock(data->fAccessElement);
	int nextel = data->fNextElement;
	int numprocessed = data->fProcessed2.size();
	while(nextel < nel || numprocessed)
	{
		std::set<int>::iterator itprocess;
		bool keeplooking = false;
		if(data->fSubmitted2.size() && data->fProcessed2.size())
		{
			auto itavail = data->fSubmitted2.begin();
			itprocess = data->fProcessed2.begin();
			if(itavail->first == *itprocess)
			{
				// make sure we come back to look for one more element
				keeplooking = true;
				// Get a hold of the data
				data->fProcessed2.erase(itprocess);
				TPZAutoPointer<TPZElementMatrixT<STATE>> ek = itavail->second;
				data->fSubmitted2.erase(itavail);
#ifdef PZ_LOG
				int iel = *itprocess;
				std::stringstream sout;
				sout << "Assembling element " << iel;
				if (logger.isDebugEnabled())
				{
					LOGPZ_DEBUG(logger, sout.str())
				}
#endif
				// Release the mutex
                lock.unlock();
				TPZManVector<int64_t,300> destindex(ek->fDestinationIndex);
				data->PermuteScatter(destindex);
				
				// Assemble the matrix
				if(!ek->HasDependency())
				{
					data->fGlobMatrix2->AddKel(ek->fMat,ek->fSourceIndex,destindex);
				}
				else
				{
					data->fGlobMatrix2->AddKel(ek->fConstrMat,ek->fSourceIndex,destindex);
				}
				// acquire the mutex
                lock.lock();
			}
		}
		if(!keeplooking)
		  {
              lock.unlock();
#ifdef PZ_LOG
              if(logger.isDebugEnabled()) LOGPZ_DEBUG(logger, "Going to sleep within assembly")
#endif
              // wait for a signal
              data->fAssembly2.Wait();
#ifdef PZ_LOG
              if(logger.isDebugEnabled()) LOGPZ_DEBUG(logger, "Waking up for assembly")
#endif
              lock.lock();
		}
		nextel = data->fNextElement;
		numprocessed = data->fProcessed2.size();
		
	}
	return 0;
}		

int TPZPairStructMatrix::ThreadData::NextElement()
{
    std::unique_lock<std::mutex> lock(fAccessElement);
	int iel;
	int nextel = fNextElement;
	TPZCompMesh *cmesh = fStrMatrix->Mesh();
	TPZAdmChunkVector<TPZCompEl *> &elementvec = cmesh->ElementVec();
	int nel = elementvec.NElements();
	for(iel=fNextElement; iel < nel; iel++)
	{
		TPZCompEl *el = elementvec[iel];
		if(!el) continue;
		TPZMaterial * mat = el->Material();
		TPZSubCompMesh *submesh = dynamic_cast<TPZSubCompMesh *> (el);
		if(!mat)
		{
			if(!submesh)
			{
				continue;
			}
			else if(submesh->NeedsComputing(fStrMatrix->MaterialIds()) == false) continue;
		}
		else 
		{
			int matid = mat->Id();
			if(this->ShouldCompute(matid) == false) continue;
		}
		break;
	}
	fNextElement = iel+1;
	nextel = iel;
	if(iel<nel) 
	{
		fProcessed1.insert(iel);
		fProcessed2.insert(iel);
	}
#ifdef PZ_LOG
	{
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << " returning " << nextel << " fNextElement " << fNextElement;
		if (logger.isDebugEnabled())
		{
			LOGPZ_DEBUG(logger, sout.str())
		}
	}
#endif
	return nextel;
}

// put the computed element matrices in the map
void TPZPairStructMatrix::ThreadData::ComputedElementMatrix(int iel, TPZAutoPointer<TPZElementMatrixT<STATE>> &ek, TPZAutoPointer<TPZElementMatrixT<STATE>> &ef)
{
    unique_lock<std::mutex> lock(fAccessElement);
	std::pair< TPZAutoPointer<TPZElementMatrixT<STATE>>, TPZAutoPointer<TPZElementMatrixT<STATE>> > el(ek,ef);
	fSubmitted1[iel] = el;
	fSubmitted2[iel] = ek;
	fAssembly1.Post();
	fAssembly2.Post();
}

// Set the set of material ids which will be considered when assembling the system
void TPZPairStructMatrix::SetMaterialIds(const std::set<int> &materialids)
{
	fStrMatrix->SetMaterialIds(materialids);
#ifdef PZ_LOG
	{
		std::set<int>::const_iterator it;
		std::stringstream sout;
		sout << "setting input material ids ";
		for(it=materialids.begin(); it!= materialids.end(); it++)
		{
			sout << *it << " ";
		}
		if (logger.isDebugEnabled())
		{
			LOGPZ_DEBUG(logger, sout.str())
		}
	}
#endif
}
