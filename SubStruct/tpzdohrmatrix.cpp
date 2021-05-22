/**
 * @file
 * @brief Contains the implementation of the TPZDohrMatrix methods. 
 */

#include "tpzdohrmatrix.h"

#include "tpzdohrassembly.h"
#include "pzlog.h"

#include "TPZSimpleTimer.h"
#include "TPZTimeTemp.h"

#include <thread>

#include "tpzparallelenviroment.h"

#ifdef PZ_LOG
static TPZLogger logger("substruct.dohrsubstruct");
#endif


template<class TVar, class TSubStruct>
int64_t TPZDohrMatrix<TVar,TSubStruct>::Size() const
{
  PZError<<__PRETTY_FUNCTION__;
  PZError<<"ERROR:Should not be called\n.Aborting...\n";
  DebugStop();
	return -1;
}
template<class TVar, class TSubStruct>
TVar* &TPZDohrMatrix<TVar,TSubStruct>::Elem()
{
  PZError<<__PRETTY_FUNCTION__;
  PZError<<"ERROR:Should not be called\n.Aborting...\n";
  DebugStop();
	static TVar* t{nullptr};
  return t;
}
template<class TVar, class TSubStruct>
const TVar* TPZDohrMatrix<TVar,TSubStruct>::Elem() const
{
  PZError<<__PRETTY_FUNCTION__;
  PZError<<"ERROR:Should not be called\n.Aborting...\n";
  DebugStop();
	return nullptr;
}

template<class TVar, class TSubStruct>
TPZDohrMatrix<TVar,TSubStruct>::TPZDohrMatrix(TPZAutoPointer<TPZDohrAssembly<TVar> > assembly)
: TPZMatrix<TVar>(), fNumThreads(0), fAssembly(assembly)
{
}

template<class TVar, class TSubStruct>
TPZDohrMatrix<TVar,TSubStruct>::~TPZDohrMatrix()
{
}

/** Threading Building Blocks */

#ifdef USING_TBB
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/partitioner.h>
using namespace tbb;
#endif

#include "arglib.h"

template<class TVar, class TSubStruct>
class ParallelAssembleTaskMatrix
{
private:
    /** @brief Array of work items */
    std::vector<TPZDohrThreadMultData<TSubStruct> > mWorkItems;
    
    /** @brief The vector with which we will multiply */
	const TPZFMatrix<TVar> *fInput;
	/** @brief Scalar multiplication factor */
	TVar fAlpha;
	/** @brief The data structure which defines the assemble destinations */
	TPZAutoPointer<TPZDohrAssembly<TVar> > fAssembly;
	/** @brief The list of data objects which need to treated by the threads */
	std::list<TPZDohrThreadMultData<TSubStruct> > fWork;
	/** @brief The local contribution to the v2 vector */
	TPZAutoPointer<TPZDohrAssembleList<TVar> > fAssemblyStructure;
	
	
    
public:
    
    /** @brief Constructor */
    ParallelAssembleTaskMatrix(const TPZFMatrix<TVar> &x, TVar alpha, TPZAutoPointer<TPZDohrAssembly<TVar> > assembly, TPZAutoPointer<TPZDohrAssembleList<TVar> > &assemblestruct) : fInput(&x), fAlpha(alpha), fAssembly(assembly), fAssemblyStructure(assemblestruct) {};
    
    /** @brief Add a new work item to the array */
    void addWorkItem(TPZDohrThreadMultData<TSubStruct> data) {
        mWorkItems.push_back(data);
    }
    
#ifdef USING_TBB
    /** @brief Computing operator for the parallel for. */
    void operator()(const blocked_range<size_t>& range) const
    {
        
        for(size_t i=range.begin(); i!=range.end(); ++i )
        {
            TPZDohrThreadMultData<TSubStruct> runner = mWorkItems[i];
            TPZFMatrix<TVar> xlocal;
            fAssembly->Extract(runner.fisub,*(fInput),xlocal);
            TPZAutoPointer<TPZDohrAssembleItem<TVar> > assembleItem = new TPZDohrAssembleItem<TVar>(runner.fisub,xlocal.Rows(),xlocal.Cols());
            runner.fSub->ContributeKULocal(fAlpha,xlocal,assembleItem->fAssembleData);
            fAssemblyStructure->AddItem(assembleItem);
        }
    }
    
    /** Execute work items in parallel. */
    void run_parallel_for(affinity_partitioner &ap)
    {
        /* TBB Parallel for. It will split the range
         * into N sub-ranges and
         * invoke the operator() for each sub-range.
         */
        parallel_for(blocked_range<size_t>(0, mWorkItems.size()), *this, ap);
    }
#endif
    
    
}; /* ParallelAssembleTask */

template<class TVar, class TSubStruct>
void TPZDohrMatrix<TVar,TSubStruct>::MultAddTBB(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
                                                const TVar alpha,const TVar beta,const int opt) const
{

#ifdef USING_TBB
    
	if ((!opt && this->Cols() != x.Rows()) || this->Rows() != x.Rows())
		this->Error( "Operator* <matrixs with incompatible dimensions>" );
	if(x.Cols() != y.Cols() || x.Cols() != z.Cols() || x.Rows() != y.Rows() || x.Rows() != z.Rows()) {
		this->Error ("TPZFMatrix::MultiplyAdd incompatible dimensions\n");
	}
	this->PrepareZ(y,z,beta,opt);
	
    
    unsigned int nglob = fGlobal.size();
    TPZAutoPointer<TPZDohrAssembleList<TVar> > assemblelist = new TPZDohrAssembleList<TVar>(nglob,z,this->fAssembly);
    
    ParallelAssembleTaskMatrix<TVar,TSubStruct> multwork(x,alpha,fAssembly,assemblelist);
    typename std::list<TPZAutoPointer<TSubStruct> >::const_iterator iter;
    int isub=0;
    for (iter=fGlobal.begin(); iter!=fGlobal.end(); iter++,isub++) {
        TPZDohrThreadMultData<TSubStruct> data(isub,*iter);
        
        multwork.addWorkItem(data);
    }

    multwork.run_parallel_for(pzenviroment.fSubstructurePartitioner);

    std::thread t(TPZDohrAssembleList<TVar>::Assemble, assemblelist.operator->());
    t.join();
#endif

}


template<class TVar, class TSubStruct>
void TPZDohrMatrix<TVar,TSubStruct>::MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
											 const TVar alpha,const TVar beta,const int opt) const
{
    
#ifdef USING_TBB 
        MultAddTBB(x, y, z, alpha, beta, opt);
        return;
#endif
        
	TPZSimpleTimer mult;
	if ((!opt && this->Cols() != x.Rows()) || this->Rows() != x.Rows())
		this->Error( "Operator* <matrixs with incompatible dimensions>" );
	if(x.Cols() != y.Cols() || x.Cols() != z.Cols() || x.Rows() != y.Rows() || x.Rows() != z.Rows()) {
		this->Error ("TPZFMatrix::MultiplyAdd incompatible dimensions\n");
	}
	this->PrepareZ(y,z,beta,opt);
	
	typename SubsList::const_iterator iter;
	int isub = 0;
	if (fNumThreads == 0) {
		for (iter=fGlobal.begin();iter!=fGlobal.end();iter++,isub++) {
            if(0)
            {
                TPZPersistenceManager::OpenWrite("dohr.txt");
                TPZPersistenceManager::WriteToFile(fAssembly.operator ->());
                TPZPersistenceManager::WriteToFile(&x);
                TPZAutoPointer<TSubStruct> point(*iter);
                TPZPersistenceManager::WriteToFile(point.operator ->());
                TPZPersistenceManager::CloseWrite();
                
            }
			TPZFMatrix<TVar> xlocal,zlocal;
			fAssembly->Extract(isub,x,xlocal);
			zlocal.Redim(xlocal.Rows(),xlocal.Cols());
			(*iter)->ContributeKULocal(alpha,xlocal,zlocal);
			fAssembly->Assemble(isub,zlocal,z);
			//         z.Print("Resultado intermediario");
		}		
	}
	else {
        unsigned int nglob = fGlobal.size();
		TPZAutoPointer<TPZDohrAssembleList<TVar> > assemblelist = new TPZDohrAssembleList<TVar>(nglob,z,this->fAssembly);
		
		TPZDohrThreadMultList<TVar,TSubStruct> multwork(x,alpha,fAssembly,assemblelist);
		typename std::list<TPZAutoPointer<TSubStruct> >::const_iterator iter;
		int isub=0;
		for (iter=fGlobal.begin(); iter!=fGlobal.end(); iter++,isub++) {
			TPZDohrThreadMultData<TSubStruct> data(isub,*iter);
            
            multwork.AddItem(data);
		}
		std::vector<std::thread> listThreads(fNumThreads);
        int i;
        for (i = 0; i < fNumThreads; i++) {
              listThreads[i] = std::thread(TPZDohrThreadMultList<TVar,TSubStruct>::ThreadWork, &multwork);
        }
        std::thread assembleThread(TPZDohrAssembleList<TVar>::Assemble, assemblelist.operator->());
        assembleThread.join();
        for (i = 0; i < fNumThreads; i++) {
          listThreads[i].join();
        }
	}
	tempo.fMultiply.Push(mult.ReturnTimeDouble());
}

template<class TVar, class TSubStruct>
void TPZDohrMatrix<TVar,TSubStruct>::Initialize() 
{
	std::cout << "Number of substructures " << fGlobal.size() << std::endl;
	tempo.fNumSub = fGlobal.size();																// alimenta timeTemp com o numero de substruturas
	TPZFMatrix<TVar> diag(this->Rows(),1,0.);
	typename SubsList::iterator iter;
	int isub = 0;
	for (iter=fGlobal.begin();iter!=fGlobal.end();iter++,isub++) {
        //Basic initialization for each substructure (compute the matrices)
        //(*iter)->Initialize();
		TPZFMatrix<TVar> diaglocal;
        (*iter)->ContributeDiagonalLocal(diaglocal);
#ifdef PZ_LOG
        if(logger.isDebugEnabled())
        {
            LOGPZ_DEBUG(logger,"Before assemble diagonal")
        }
#endif
		this->fAssembly->Assemble(isub,diaglocal,diag);
#ifdef PZ_LOG
        if(logger.isDebugEnabled())
		{
			std::stringstream sout;
			sout << "Substructure " << isub << " ";
			diag.Print("Global Diagonal matrix",sout);
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
        std::cout << '*';
        std::cout.flush();
	}
#ifdef PZ_LOG
    if(logger.isDebugEnabled())
	{
		std::stringstream sout;
		diag.Print("Global Diagonal matrix",sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	std::cout << std::endl;
	for (iter=fGlobal.begin(),isub=0;iter!=fGlobal.end();iter++,isub++) {
        //Computes the Weights for each substructure
		TPZFMatrix<TVar> diaglocal;
		this->fAssembly->Extract(isub,diag,diaglocal);
        (*iter)->ComputeWeightsLocal(diaglocal);
		
	}
}

/**
 * Adjust the residual to zero the residual of the internal connects
 */
template<class TVar, class TSubStruct>
void TPZDohrMatrix<TVar,TSubStruct>::AdjustResidual(TPZFMatrix<TVar> &res)
{
	typename SubsList::iterator iter;
	for (iter=fGlobal.begin();iter!=fGlobal.end();iter++) {
		(*iter)->AdjustResidual(res);
	}
}

/**
 * Add the solution corresponding to the internal residual
 */
template<class TVar, class TSubStruct>
void TPZDohrMatrix<TVar,TSubStruct>::AddInternalSolution(TPZFMatrix<TVar> &solution)
{
	typename SubsList::iterator iter;
	for (iter=fGlobal.begin();iter!=fGlobal.end();iter++) {
		(*iter)->AddInternalSolution(solution);
	}
}

template<class TVar, class TSubStruct>
void *TPZDohrThreadMultList<TVar,TSubStruct>::ThreadWork(void *ptr)
{
	TPZDohrThreadMultList<TVar,TSubStruct> *myptr = (TPZDohrThreadMultList<TVar,TSubStruct> *) ptr;
	TPZDohrThreadMultData<TSubStruct> runner = myptr->PopItem();
	while (runner.IsValid()) {
		TPZFMatrix<TVar> xlocal;
		myptr->fAssembly->Extract(runner.fisub,*(myptr->fInput),xlocal);
		TPZAutoPointer<TPZDohrAssembleItem<TVar> > assembleItem = new TPZDohrAssembleItem<TVar>(runner.fisub,xlocal.Rows(),xlocal.Cols());
		runner.fSub->ContributeKULocal(myptr->fAlpha,xlocal,assembleItem->fAssembleData);
		myptr->fAssemblyStructure->AddItem(assembleItem);
		runner = myptr->PopItem();
	}
	return ptr;
}

/**
 * @brief Unpacks the object structure from a stream of bytes
 * @param buf The buffer containing the object in a packed form
 * @param context 
 */
template <class TVar, class TSubStruct>
void TPZDohrMatrix<TVar, TSubStruct >::Read(TPZStream &buf, void *context )
{
    SAVEABLE_SKIP_NOTE(buf);
    TPZMatrix<TVar>::Read(buf, context);
    SAVEABLE_SKIP_NOTE(buf);
    fAssembly = TPZAutoPointerDynamicCast<TPZDohrAssembly<TVar>>(TPZPersistenceManager::GetAutoPointer(&buf));
    SAVEABLE_SKIP_NOTE(buf);
    buf.Read(&fNumCoarse);
    SAVEABLE_SKIP_NOTE(buf);
    buf.Read(&fNumThreads);
    int sz;
    SAVEABLE_SKIP_NOTE(buf);
    buf.Read(&sz);
    for (int i=0; i<sz; i++) {
        TPZAutoPointer<TSubStruct > sub = new TSubStruct;
        SAVEABLE_SKIP_NOTE(buf);
        sub->Read(buf,0);
        fGlobal.push_back(sub);
    }
    int classid;
    SAVEABLE_SKIP_NOTE(buf);
    buf.Read(&classid );
    if (classid != ClassId()) {
        DebugStop();
    }
}
/**
 * @brief Packs the object structure in a stream of bytes
 * @param buf Buffer which will receive the bytes
 * @param withclassid
 */
template <class TVar, class TSubStruct>
void TPZDohrMatrix<TVar,TSubStruct >::Write( TPZStream &buf, int withclassid ) const
{
    SAVEABLE_STR_NOTE(buf,"TPZMatrix<TVar>::Write ()");
    TPZMatrix<TVar>::Write(buf, withclassid);
    SAVEABLE_STR_NOTE(buf,"fAssembly->Write");
    TPZPersistenceManager::WritePointer(fAssembly.operator ->(), &buf);
    SAVEABLE_STR_NOTE(buf,"fNumCoarse");
    buf.Write(&fNumCoarse);
    SAVEABLE_STR_NOTE(buf,"fNumThreads");
    buf.Write(&fNumThreads);
    int size = fGlobal.size();
    SAVEABLE_STR_NOTE(buf,"fGlobal.size()");
    buf.Write(&size);
    for (auto it=fGlobal.begin(); it != fGlobal.end(); it++) {
        SAVEABLE_STR_NOTE(buf,"fGlobal[...]");
        (*it)->Write(buf,0);
    }
    size = 0;
    int classid = ClassId();
    SAVEABLE_STR_NOTE(buf,"ClassId");
    buf.Write(&classid );
}

template <>
void TPZDohrMatrix<long double, TPZDohrSubstructCondense<long double> >::Read(TPZStream &buf, void *context )
{
	DebugStop();
}
template <>
void TPZDohrMatrix<float, TPZDohrSubstruct<float> >::Read(TPZStream &buf, void *context )
{
    DebugStop();
}
template <>
void TPZDohrMatrix<double, TPZDohrSubstruct<double> >::Read(TPZStream &buf, void *context )
{
    DebugStop();
}
template <>
void TPZDohrMatrix<long double, TPZDohrSubstruct<long double> >::Read(TPZStream &buf, void *context )
{
    DebugStop();
}

template <>
void TPZDohrMatrix<long double,TPZDohrSubstructCondense<long double> >::Write( TPZStream &buf, int withclassid ) const
{
    DebugStop();
}
template <>
void TPZDohrMatrix<float, TPZDohrSubstruct<float> >::Write( TPZStream &buf, int withclassid ) const
{
    DebugStop();
}
template <>
void TPZDohrMatrix<double, TPZDohrSubstruct<double> >::Write( TPZStream &buf, int withclassid ) const
{
    DebugStop();
}
template <>
void TPZDohrMatrix<long double, TPZDohrSubstruct<long double> >::Write( TPZStream &buf, int withclassid ) const
{
    DebugStop();
}


template class TPZDohrMatrix<float, TPZDohrSubstruct<float> >;
template class TPZDohrMatrix<double, TPZDohrSubstruct<double> >;
template class TPZDohrMatrix<long double, TPZDohrSubstruct<long double> >;

template class TPZDohrMatrix<float, TPZDohrSubstructCondense<float> >;
template class TPZDohrMatrix<double, TPZDohrSubstructCondense<double> >;
template class TPZDohrMatrix<long double, TPZDohrSubstructCondense<long double> >;

//template class TPZDohrMatrix<std::complex<float>, TPZDohrSubstruct<std::complex<float> > >;
template class TPZDohrMatrix<std::complex<double >, TPZDohrSubstruct<std::complex<double> > >;
//template class TPZDohrMatrix<std::complex<long double>, TPZDohrSubstruct<std::complex<long double> > >;

//template class TPZDohrMatrix<std::complex<float>, TPZDohrSubstructCondense<std::complex<float> > >;
template class TPZDohrMatrix<std::complex<double>, TPZDohrSubstructCondense<std::complex<double> > >;
//template class TPZDohrMatrix<std::complex<long double>, TPZDohrSubstructCondense<std::complex<long double> > >;

#ifndef BORLAND
template class TPZRestoreClass< TPZDohrMatrix<double, TPZDohrSubstructCondense<double> > >;
template class TPZRestoreClass< TPZDohrMatrix<double, TPZDohrSubstruct<double> > >;
template class TPZRestoreClass< TPZDohrMatrix<float, TPZDohrSubstructCondense<float> > >;
template class TPZRestoreClass< TPZDohrMatrix<float, TPZDohrSubstruct<float> > >;
#endif
