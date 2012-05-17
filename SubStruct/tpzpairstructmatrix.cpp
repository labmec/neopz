/**
 * @file
 * @brief Contains the implementation of the TPZPairStructMatrix methods. 
 */

#include "tpzpairstructmatrix.h"
#include "TPZTimer.h"
#include "pzelmat.h"
#include "pzcompel.h"
#include "pzstrmatrix.h"
#include "pzsubcmesh.h"
#include "pzanalysis.h"
#include "pzmaterial.h"

using namespace std;

#include "pzlog.h"

#include "pz_pthread.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.strmatrix.tpzpairstructmatrix"));
static LoggerPtr loggerel(Logger::getLogger("pz.strmatrix.element"));
#endif

int TPZPairStructMatrix::gNumThreads = 0;

TPZPairStructMatrix::TPZPairStructMatrix(TPZCompMesh *mesh, TPZVec<int> &permutescatter)
{
	fMesh = mesh;
	fPermuteScatter = permutescatter;
#ifdef PERF_DEBUG
	std::cout << "fNumThreads e "<< fNumThreads << " gNumThreads esta setando para " << gNumThreads << endl; 
#endif
	fNumThreads = gNumThreads;
}

void TPZPairStructMatrix::SerialAssemble(int mineq, int maxeq, TPZMatrix<STATE> *first, TPZMatrix<STATE> *second, TPZFMatrix<STATE> &rhs)
{
	int iel;
	TPZCompMesh &mesh = *fMesh;
	int nelem = mesh.NElements();
	TPZElementMatrix ek(&mesh, TPZElementMatrix::EK),ef(&mesh, TPZElementMatrix::EF);
    
	
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
			if(mineq != -1 && maxeq != -1)
			{
				TPZStructMatrix::FilterEquations(ek.fSourceIndex,ek.fDestinationIndex,mineq,maxeq);
			}
			first->AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
			rhs.AddFel(ef.fMat,ek.fSourceIndex,ek.fDestinationIndex);
			PermuteScatter(ek.fDestinationIndex);
			second->AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
#ifdef LOG4CXX
			if(loggerel->isDebugEnabled())
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
			if(mineq != -1 & maxeq != -1)
			{
				TPZStructMatrix::FilterEquations(ek.fSourceIndex,ek.fDestinationIndex,mineq,maxeq);
			}
			first->AddKel(ek.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
			rhs.AddFel(ef.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
			PermuteScatter(ek.fDestinationIndex);
			second->AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
		}
		
		assemble.stop();
	}//fim for iel
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Number of equations " << first->Rows() << std::endl;
		sout << calcstiff.processName() << " " << calcstiff << std::endl;
		sout << assemble.processName() << " " << assemble;
		/*    stiffness.Print("Matriz de Rigidez: ",sout);
		 rhs.Print("Vetor de Carga: ",sout);*/
		LOGPZ_DEBUG(logger,sout.str().c_str());
	}
#endif
	
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

void TPZPairStructMatrix::Assemble(int mineq, int maxeq, TPZMatrix<STATE> *first, TPZMatrix<STATE> *second, TPZFMatrix<STATE> &rhs)
{
#ifdef PERF_DEBUG
  std::cout << "TPZPairStructMatrix::Assemble = Assembling the system of equations with " << fNumThreads
	    << " threads (TPZPairStructMatrix::gNumThreads = " << TPZPairStructMatrix::gNumThreads  << ")\n";
#endif
	if(fNumThreads <= 0)
	{
		SerialAssemble(mineq,maxeq,first,second,rhs);
		return;
	}
	else 
	{
		MultiThread_Assemble(mineq, maxeq, first, second, rhs);
	}
}

void TPZPairStructMatrix::MultiThread_Assemble(int mineq, int maxeq, TPZMatrix<STATE> *first, TPZMatrix<STATE> *second, TPZFMatrix<STATE> &rhs)
{
	ThreadData threaddata(*fMesh,*first,*second,rhs,mineq,maxeq);
	threaddata.fPermuteScatter = fPermuteScatter;
	const int numthreads = this->fNumThreads;
	std::vector<pthread_t> allthreads(numthreads+1);
	int itr;
	for(itr=0; itr<numthreads; itr++)
	{
	  PZ_PTHREAD_CREATE(&allthreads[itr], NULL,
			    ThreadData::ThreadWork, &threaddata, __FUNCTION__);
	}
	
	// assemble the first matrix
	PZ_PTHREAD_CREATE(&allthreads[itr], NULL,
			  ThreadData::ThreadAssembly1, &threaddata, __FUNCTION__);
	
	// assemble the second matrix
	ThreadData::ThreadAssembly2(&threaddata);
	
	for(itr=0; itr<numthreads+1; itr++)
	{
	  PZ_PTHREAD_JOIN(allthreads[itr], NULL, __FUNCTION__);
	}
}

TPZPairStructMatrix::ThreadData::ThreadData(TPZCompMesh &mesh, TPZMatrix<STATE> &mat1, TPZMatrix<STATE> &mat2, 
											TPZFMatrix<STATE> &rhs, int mineq, int maxeq)
: fMesh(&mesh), 
fGlobMatrix1(&mat1), fGlobMatrix2(&mat2), fGlobRhs(&rhs),
fMinEq(mineq), fMaxEq(maxeq),fNextElement(0)
{	
  PZ_PTHREAD_MUTEX_INIT(&fAccessElement,NULL,"TPZPairStructMatrix::ThreadData::ThreadData()");
}

TPZPairStructMatrix::ThreadData::~ThreadData()
{
  PZ_PTHREAD_MUTEX_DESTROY(&fAccessElement,"TPZPairStructMatrix::ThreadData::~ThreadData()");
}

void *TPZPairStructMatrix::ThreadData::ThreadWork(void *datavoid)
{
	ThreadData *data = (ThreadData *) datavoid;
	// compute the next element (this method is threadsafe)
	int iel = data->NextElement();
	TPZCompMesh *cmesh = data->fMesh;
	int nel = cmesh->NElements();
	while(iel < nel)
	{
		
		TPZAutoPointer<TPZElementMatrix> ek = new TPZElementMatrix(cmesh,TPZElementMatrix::EK);
		TPZAutoPointer<TPZElementMatrix> ef = new TPZElementMatrix(cmesh,TPZElementMatrix::EF);
		
		TPZCompEl *el = cmesh->ElementVec()[iel];
		TPZElementMatrix *ekp = ek.operator->();
		TPZElementMatrix *efp = ef.operator->();
		TPZElementMatrix &ekr = *ekp;
		TPZElementMatrix &efr = *efp;
		el->CalcStiff(ekr,efr);
		
		
		if(!el->HasDependency()) {
			ek->ComputeDestinationIndices();
			
			if(data->fMinEq != -1 || data->fMaxEq != -1)
			{
				TPZStructMatrix::FilterEquations(ek->fSourceIndex,ek->fDestinationIndex,data->fMinEq,data->fMaxEq);
			}
#ifdef LOG4CXX
			if(loggerel->isDebugEnabled())
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
			if(data->fMinEq != -1 || data->fMaxEq != -1)
			{
				TPZStructMatrix::FilterEquations(ek->fSourceIndex,ek->fDestinationIndex,data->fMinEq,data->fMaxEq);
			}
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


// The function which will compute the assembly
void *TPZPairStructMatrix::ThreadData::ThreadAssembly1(void *threaddata)
{
	ThreadData *data = (ThreadData *) threaddata;
	TPZCompMesh *cmesh = data->fMesh;
	int nel = cmesh->NElements();
	PZ_PTHREAD_MUTEX_LOCK(&(data->fAccessElement),"TPZPairStructMatrix::ThreadData::ThreadAssembly1()");
	int nextel = data->fNextElement;
	int numprocessed = data->fProcessed1.size();
	while(nextel < nel || numprocessed)
	{
		std::map<int, std::pair< TPZAutoPointer<TPZElementMatrix>, TPZAutoPointer<TPZElementMatrix> > >::iterator itavail;
		std::set<int>::iterator itprocess;
		bool keeplooking = false;
		if(data->fSubmitted1.size() && data->fProcessed1.size())
		{
			itavail = data->fSubmitted1.begin();
			itprocess = data->fProcessed1.begin();
			if(itavail->first == *itprocess)
			{
				// make sure we come back to look for one more element
				keeplooking = true;
				// Get a hold of the data
				int iel = *itprocess;
				data->fProcessed1.erase(itprocess);
				TPZAutoPointer<TPZElementMatrix> ek = itavail->second.first;
				TPZAutoPointer<TPZElementMatrix> ef = itavail->second.second;
				data->fSubmitted1.erase(itavail);
#ifdef LOG4CXX
				std::stringstream sout;
				sout << "Assembling element " << iel;
				LOGPZ_DEBUG(logger,sout.str())
#endif
				// Release the mutex
				PZ_PTHREAD_MUTEX_UNLOCK(&data->fAccessElement,"TPZPairStructMatrix::ThreadData::ThreadAssembly1()");
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
				PZ_PTHREAD_MUTEX_LOCK(&data->fAccessElement,"TPZPairStructMatrix::ThreadData::ThreadAssembly1()");
			}
		}
		if(!keeplooking)
		{
		        PZ_PTHREAD_MUTEX_UNLOCK(&data->fAccessElement,"TPZPairStructMatrix::ThreadData::ThreadAssembly1()");
			LOGPZ_DEBUG(logger,"Going to sleep within assembly")
			// wait for a signal
			data->fAssembly1.Wait();
			LOGPZ_DEBUG(logger,"Waking up for assembly")
		        PZ_PTHREAD_MUTEX_LOCK(&data->fAccessElement,"TPZPairStructMatrix::ThreadData::ThreadAssembly1()");
		}
		nextel = data->fNextElement;
		numprocessed = data->fProcessed1.size();
		
	}
	PZ_PTHREAD_MUTEX_UNLOCK(&data->fAccessElement,"TPZPairStructMatrix::ThreadData::ThreadAssembly1()");
	return 0;	
}		

// The function which will compute the assembly
void *TPZPairStructMatrix::ThreadData::ThreadAssembly2(void *threaddata)
{
	ThreadData *data = (ThreadData *) threaddata;
	TPZCompMesh *cmesh = data->fMesh;
	int nel = cmesh->NElements();
	PZ_PTHREAD_MUTEX_LOCK(&(data->fAccessElement),"TPZPairStructMatrix::ThreadData::ThreadAssembly2()");
	int nextel = data->fNextElement;
	int numprocessed = data->fProcessed2.size();
	while(nextel < nel || numprocessed)
	{
		std::map<int, TPZAutoPointer<TPZElementMatrix> >::iterator itavail;
		std::set<int>::iterator itprocess;
		bool keeplooking = false;
		if(data->fSubmitted2.size() && data->fProcessed2.size())
		{
			itavail = data->fSubmitted2.begin();
			itprocess = data->fProcessed2.begin();
			if(itavail->first == *itprocess)
			{
				// make sure we come back to look for one more element
				keeplooking = true;
				// Get a hold of the data
				int iel = *itprocess;
				data->fProcessed2.erase(itprocess);
				TPZAutoPointer<TPZElementMatrix> ek = itavail->second;
				data->fSubmitted2.erase(itavail);
#ifdef LOG4CXX
				std::stringstream sout;
				sout << "Assembling element " << iel;
				LOGPZ_DEBUG(logger,sout.str())
#endif
				// Release the mutex
				PZ_PTHREAD_MUTEX_UNLOCK(&(data->fAccessElement),"TPZPairStructMatrix::ThreadData::ThreadAssembly2()");
				TPZManVector<int,300> destindex(ek->fDestinationIndex);
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
				PZ_PTHREAD_MUTEX_LOCK(&(data->fAccessElement),"TPZPairStructMatrix::ThreadData::ThreadAssembly2()");
			}
		}
		if(!keeplooking)
		  {
		        PZ_PTHREAD_MUTEX_UNLOCK(&(data->fAccessElement),"TPZPairStructMatrix::ThreadData::ThreadAssembly2()");
			LOGPZ_DEBUG(logger,"Going to sleep within assembly")
			// wait for a signal
			data->fAssembly2.Wait();
			LOGPZ_DEBUG(logger,"Waking up for assembly")
			PZ_PTHREAD_MUTEX_LOCK(&data->fAccessElement,"TPZPairStructMatrix::ThreadData::ThreadAssembly2()");
		}
		nextel = data->fNextElement;
		numprocessed = data->fProcessed2.size();
		
	}
	PZ_PTHREAD_MUTEX_UNLOCK(&data->fAccessElement,"TPZPairStructMatrix::ThreadData::ThreadAssembly2()");
	return 0;	
}		

int TPZPairStructMatrix::ThreadData::NextElement()
{
        PZ_PTHREAD_MUTEX_LOCK(&fAccessElement,"TPZPairStructMatrix::ThreadData::NextElement()");
	int iel;
	int nextel = fNextElement;
	TPZCompMesh *cmesh = fMesh;
	TPZAdmChunkVector<TPZCompEl *> &elementvec = cmesh->ElementVec();
	int nel = elementvec.NElements();
	for(iel=fNextElement; iel < nel; iel++)
	{
		TPZCompEl *el = elementvec[iel];
		if(!el) continue;
		if(fMaterialIds.size() == 0) break;
		TPZMaterial * mat = el->Material();
		TPZSubCompMesh *submesh = dynamic_cast<TPZSubCompMesh *> (el);
		if(!mat)
		{
			if(!submesh)
			{
				continue;
			}
			else if(submesh->NeedsComputing(fMaterialIds) == false) continue;
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
        PZ_PTHREAD_MUTEX_UNLOCK(&fAccessElement,"TPZPairStructMatrix::ThreadData::NextElement()");
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << " returning " << nextel << " fNextElement " << fNextElement;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	return nextel;
}

// put the computed element matrices in the map
void TPZPairStructMatrix::ThreadData::ComputedElementMatrix(int iel, TPZAutoPointer<TPZElementMatrix> &ek, TPZAutoPointer<TPZElementMatrix> &ef)
{
  //FIXME: Edson, este metodo precisa de lock (exclusao mutua)?
        PZ_PTHREAD_MUTEX_LOCK(&fAccessElement,"TPZPairStructMatrix::ThreadData::ComputedElementMatrix()");
	std::pair< TPZAutoPointer<TPZElementMatrix>, TPZAutoPointer<TPZElementMatrix> > el(ek,ef);
	fSubmitted1[iel] = el;
	fSubmitted2[iel] = ek;
	fAssembly1.Post();
	fAssembly2.Post();
        PZ_PTHREAD_MUTEX_UNLOCK(&fAccessElement,"TPZPairStructMatrix::ThreadData::ComputedElementMatrix()");	
}

// Set the set of material ids which will be considered when assembling the system
void TPZPairStructMatrix::SetMaterialIds(const std::set<int> &materialids)
{
	fMaterialIds = materialids;
#ifdef LOG4CXX
	{
		std::set<int>::const_iterator it;
		std::stringstream sout;
		sout << "setting input material ids ";
		for(it=materialids.begin(); it!= materialids.end(); it++)
		{
			sout << *it << " ";
		}
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	if(!fMesh)
	{
		LOGPZ_WARN(logger,"SetMaterialIds called without mesh")
		return;
	}
	int iel;
	TPZAdmChunkVector<TPZCompEl*> &elvec = fMesh->ElementVec();
	int nel = elvec.NElements();
	for(iel=0; iel<nel; iel++)
	{
		TPZCompEl *cel = elvec[iel];
		if(!cel) continue;
		TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *> (cel);
		if(!subcmesh) continue;
		TPZAutoPointer<TPZAnalysis> anal = subcmesh->Analysis();
		if(!anal)
		{
			LOGPZ_ERROR(logger,"SetMaterialIds called for substructure without analysis object")
			DebugStop();
			//subcmesh->SetAnalysis();
			//anal=subcmesh->GetAnalysis();
		}
		TPZAutoPointer<TPZStructMatrix> str = anal->StructMatrix();
		if(!str)
		{
			LOGPZ_WARN(logger,"SetMaterialIds called for substructure without structural matrix")
			continue;
		}
		str->SetMaterialIds(materialids);
	}
}
