/**
 * @file
 * @brief Contains the implementation of the TPZDohrPrecond methods.
 * @author Philippe Devloo
 * @since 2006
 */

#include "tpzdohrprecond.h"
#include "tpzdohrsubstructCondense.h"
#include "pzskylmat.h"

#include "pzvisualmatrix.h"
#include "tpzdohrassemblelist.h"

#include <sstream>
#include "pzlog.h"

#include "TPZSimpleTimer.h"
#include "TPZTimeTemp.h"


#include "arglib.h"

#include "tpzparallelenviroment.h"
#include "TPZPersistenceManager.h"

#include <thread>

#ifdef PZ_LOG
static TPZLogger logger("substruct.dohrprecond");
static TPZLogger loggerv1v2("substruct.v1v2");
#endif

using namespace std;

template<class TVar, class TSubStruct>
int64_t TPZDohrPrecond<TVar,TSubStruct>::Size() const
{
  PZError<<__PRETTY_FUNCTION__;
  PZError<<"ERROR:Should not be called\n.Aborting...\n";
  DebugStop();
	return -1;
}
template<class TVar, class TSubStruct>
TVar* &TPZDohrPrecond<TVar,TSubStruct>::Elem()
{
  PZError<<__PRETTY_FUNCTION__;
  PZError<<"ERROR:Should not be called\n.Aborting...\n";
  DebugStop();
	static TVar* t{nullptr};
  return t;
}
template<class TVar, class TSubStruct>
const TVar* TPZDohrPrecond<TVar,TSubStruct>::Elem() const
{
  PZError<<__PRETTY_FUNCTION__;
  PZError<<"ERROR:Should not be called\n.Aborting...\n";
  DebugStop();
	return nullptr;
}

template<class TVar, class TSubStruct>
TPZDohrPrecond<TVar, TSubStruct>::TPZDohrPrecond(TPZDohrMatrix<TVar, TSubStruct> &origin, TPZAutoPointer<TPZDohrAssembly<TVar> > assemble)
: TPZMatrix<TVar>(origin), fGlobal(origin.SubStructures()), fCoarse(0), fNumCoarse(origin.NumCoarse()), fNumThreads(0), fAssemble(assemble)
{
	fNumThreads = origin.NumThreads();
}

template<class TVar, class TSubStruct>
TPZDohrPrecond<TVar, TSubStruct>::TPZDohrPrecond(const TPZDohrPrecond<TVar, TSubStruct> &cp) : TPZMatrix<TVar>(cp), fGlobal(cp.fGlobal), fCoarse(0),
fNumCoarse(cp.fNumCoarse), fNumThreads(cp.fNumThreads), fAssemble(cp.fAssemble)
{
	if (cp.fCoarse) {
		fCoarse = (TPZStepSolver<TVar> *) cp.fCoarse->Clone();
	}
}

/** @brief Empty constructor for restoring */
template<class TVar, class TSubStruct>
TPZDohrPrecond<TVar, TSubStruct>::TPZDohrPrecond() : fCoarse(0), fNumCoarse(-1), fNumThreads(-1)
{
}

template<class TVar, class TSubStruct>
TPZDohrPrecond<TVar, TSubStruct>::~TPZDohrPrecond()
{
	if (fCoarse)
	{
		delete fCoarse;
		fCoarse = 0;
	}
}

/** Threading Building Blocks */

#ifdef USING_TBB
#include <tbb/blocked_range.h>
#include <tbb/partitioner.h>
#include <tbb/parallel_for.h>
#endif

template<class TVar, class TSubStruct>
class ParallelAssembleTask
{
private:
    /** @brief Array of work items */
    std::vector<TPZDohrPrecondV2SubData<TVar,TSubStruct> > mWorkItems;
    
    /** @brief The local contribution to the v2 vector */
	TPZAutoPointer<TPZDohrAssembleList<TVar> > fAssemblyStructure;
    
public:
    
    /** @brief Constructor */
    ParallelAssembleTask(TPZAutoPointer<TPZDohrAssembleList<TVar> > assemblyStruct) : fAssemblyStructure(assemblyStruct) {}
    
    /** @brief Add a new work item to the array */
    void addWorkItem(TPZDohrPrecondV2SubData<TVar, TSubStruct> data) {
        mWorkItems.push_back(data);
    }
    
#ifdef USING_TBB
    /** @brief Computing operator for the parallel for. */
    void operator()(const tbb::blocked_range<size_t>& range) const
    {
        
        for(size_t i=range.begin(); i!=range.end(); ++i )
        {
            
            TPZDohrPrecondV2SubData<TVar,TSubStruct> data  = mWorkItems[i];
            data.fInput_local.ReallocForNuma(0);
			data.fSubStructure->ReallocMatRed();
            data.fSubStructure->Contribute_v2_local(data.fInput_local, data.fv2_local->fAssembleData);
            fAssemblyStructure->AddItem(data.fv2_local);
        }
    }
    
    /** Execute work items in parallel. */
    void run_parallel_for(tbb::affinity_partitioner &ap)
    {
        /* TBB Parallel for. It will split the range
         * into N sub-ranges and
         * invoke the operator() for each sub-range.
         */
        parallel_for(tbb::blocked_range<size_t>(0, mWorkItems.size()), *this, ap);
    }
    
#endif
    
    
}; /* ParallelAssembleTask */


template<class TVar, class TSubStruct>
void TPZDohrPrecond<TVar, TSubStruct>::MultAddTBB(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z, const TVar alpha,const TVar beta,const int opt) const {

#ifdef USING_TBB
    if ((!opt && this->Cols() != x.Rows()) || this->Rows() != x.Rows())
		this->Error( "Operator* <matrices with incompatible dimensions>" );
	if(x.Cols() != y.Cols() || x.Cols() != z.Cols() || x.Rows() != y.Rows() || x.Rows() != z.Rows()) {
		this->Error ("TPZFMatrix::MultiplyAdd incompatible dimensions\n");
	}

	int64_t rows = this->Rows();
	int64_t cols = this->Cols();
	this->PrepareZ(y,z,beta,opt);
    
    TPZFMatrix<TVar> v1(rows,x.Cols(),0.);
	TPZFMatrix<TVar> v2(cols,x.Cols(),0.);

    std::vector<thread> AllThreads(2);
    TPZDohrPrecondThreadV1Data<TVar,TSubStruct> v1threaddata(this,x,v1);
    AllThreads[0] = thread((TPZDohrPrecondThreadV1Data<TVar,TSubStruct>::ComputeV1),
                           &v1threaddata);

    TPZAutoPointer<TPZDohrAssembleList<TVar> > assemblelist = new TPZDohrAssembleList<TVar>(fGlobal.size(),v2,this->fAssemble);
    
    ParallelAssembleTask<TVar, TSubStruct> tbb_work(assemblelist);
    
    typename std::list<TPZAutoPointer<TSubStruct> >::const_iterator it;
    
    int isub=0;
    /** Cria tarefa que execute a distribuicao de cada elemento do fGlobal */
    for(it= fGlobal.begin(); it != fGlobal.end(); it++,isub++)
    {
        TPZFMatrix<TVar> *Residual_local = new TPZFMatrix<TVar>;
        fAssemble->Extract(isub,x,*Residual_local);
        TPZDohrPrecondV2SubData<TVar,TSubStruct> data(isub,*it,Residual_local);
        tbb_work.addWorkItem(data);
    }
    

    tbb_work.run_parallel_for(pzenviroment.fSubstructurePartitioner);

    AllThreads[1] = thread(TPZDohrAssembleList<TVar>::Assemble,
                           assemblelist.operator->());

    for (int i=0; i<2; i++) {
        AllThreads[i].join();
    }
    
    v2 += v1;

	/** Soma v1+v2+v3 com z */
    int64_t xcols = x.Cols();
    for (int64_t ic=0; ic<xcols; ic++)
    {
        for (int64_t c=0; c<rows; c++) {
            z(c,ic) += v2(c,ic);
        }
    }

#endif
}

template<class TVar, class TSubStruct>
void TPZDohrPrecond<TVar, TSubStruct>::MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z, const TVar alpha,const TVar beta,const int opt) const {
    
#ifdef USING_TBB
	MultAddTBB(x, y, z, alpha, beta, opt);
	return;
#endif
    
	if ((!opt && this->Cols() != x.Rows()) || this->Rows() != x.Rows())
		this->Error( "Operator* <matrices with incompatible dimensions>" );
	if(x.Cols() != y.Cols() || x.Cols() != z.Cols() || x.Rows() != y.Rows() || x.Rows() != z.Rows()) {
		this->Error ("TPZFMatrix::MultiplyAdd incompatible dimensions\n");
	}
	TPZSimpleTimer precondi; // init of timer
	int64_t rows = this->Rows();
	int64_t cols = this->Cols();
	int64_t c;
	this->PrepareZ(y,z,beta,opt);
#ifdef PZ_LOG
	{
		std::stringstream sout;
		x.Print("x entry vector",sout);
		if (logger.isDebugEnabled())
		{
			LOGPZ_DEBUG(logger, sout.str());
		}
		if (loggerv1v2.isDebugEnabled())
		{
			LOGPZ_DEBUG(loggerv1v2, sout.str())
		}
	}
#endif
	TPZFMatrix<TVar> v1(rows,x.Cols(),0.);
	TPZFMatrix<TVar> v2(cols,x.Cols(),0.);
	if(fNumThreads <= 0)
	{
		ComputeV1(x,v1);
		ComputeV2(x,v2);
	}
	else
	{
        std::vector<std::thread> AllThreads(fNumThreads+2);
		TPZDohrPrecondThreadV1Data<TVar,TSubStruct> v1threaddata(this,x,v1);
		
        AllThreads[0] = thread((TPZDohrPrecondThreadV1Data<TVar,TSubStruct>::ComputeV1), &v1threaddata);
		
		TPZAutoPointer<TPZDohrAssembleList<TVar> > assemblelist = new TPZDohrAssembleList<TVar>(fGlobal.size(),v2,this->fAssemble);
		
		TPZDohrPrecondV2SubDataList<TVar,TSubStruct> v2work(assemblelist);
		typename std::list<TPZAutoPointer<TSubStruct> >::const_iterator it;
		
		int isub=0;
		//Criar tarefa que execute a distribuicao de cada elemento do fGlobal
		for(it= fGlobal.begin(); it != fGlobal.end(); it++,isub++)
		{
			TPZFMatrix<TVar> *Residual_local = new TPZFMatrix<TVar>;
			fAssemble->Extract(isub,x,*Residual_local);
			TPZDohrPrecondV2SubData<TVar,TSubStruct> data(isub,*it,Residual_local);
			v2work.AddItem(data);
		}
		
		int i;
		for (i=0; i<fNumThreads; i++) {
            AllThreads[i+2] = std::thread(TPZDohrPrecondV2SubDataList<TVar,TSubStruct>::ThreadWork, &v2work);
		}
		//		v2work.ThreadWork(&v2work);
		AllThreads[1] = thread(TPZDohrAssembleList<TVar>::Assemble,
                               assemblelist.operator->());
		
		for (i=0; i<fNumThreads+2; i++) {
            AllThreads[i].join();
		}
		//		ComputeV2(x,v2);
	}
	v2 += v1;
	
#ifndef MAKEINTERNAL
	isub=0;
	//Criar tarefa que execute a distribuicao de cada elemento do fGlobal
	for(it= fGlobal.begin(); it != fGlobal.end(); it++,isub++)
	{
		TPZFNMatrix<100> v2Expand((*it)->fNEquations,1,0.), v3Expand((*it)->fNEquations,1,0.);
		int64_t neqs = (*it)->fGlobalEqs.NElements();
		TPZFMatrix<TVar> v3_local(neqs,1,0.), v2_local(neqs,1,0.);
		fAssemble->Extract(isub,v2,v2_local);
		int64_t i;
		for (i=0;i<neqs;i++)
		{
			std::pair<int,int> ind = (*it)->fGlobalEqs[i];
			v2Expand(ind.first,0) += v2_local(i,0);
		}
		
		(*it)->Contribute_v3_local(v2Expand,v3Expand);
		for (i=0;i<neqs;i++)
		{
			std::pair<int,int> ind = (*it)->fGlobalEqs[i];
			v3_local(i,0) += v3Expand(ind.first,0);
		}
#ifdef PZ_LOG
		{
			std::stringstream sout;
			v2Expand.Print("v1+v2 Expand",sout);
			v3Expand.Print("v3 Expand", sout);
			v2_local.Print("v1+v2 local",sout);
			v3_local.Print("v3 local",sout);
			if (logger.isDebugEnabled())
			{
				LOGPZ_DEBUG(logger, sout.str());
			}
		}
#endif
		fAssemble->Assemble(isub,v3_local,v2);
	}
#endif
	// wait task para finalizacao da chamada
	// esperar a versao correta do v1
	/* Sum v1+v2+v3 with z */
    int64_t xcols = x.Cols();
    for (int64_t ic=0; ic<xcols; ic++)
    {
        for (c=0; c<rows; c++) {
            z(c,ic) += v2(c,ic);
        }
    }
	tempo.fPreCond.Push(precondi.ReturnTimeDouble()); // end of timer
}


template<class TVar, class TSubStruct>
void TPZDohrPrecond<TVar, TSubStruct>::Initialize()
{
	//Compute the skyline of the coarse equations
	TPZManVector<int64_t> skyline(fNumCoarse);
	int64_t ic;
	for (ic=0; ic<fNumCoarse; ic++) {
		skyline[ic] = ic;
	}
	int64_t nsub = fAssemble->fCoarseEqs.NElements();
	int64_t isub;
	for (isub=0; isub<nsub; isub++) {
		int64_t nc = fAssemble->fCoarseEqs[isub].NElements();
		int64_t ic;
		int64_t mineq = 0;
        if(nc != 0) mineq = fAssemble->fCoarseEqs[isub][0];
		for (ic=0; ic<nc; ic++) {
			int64_t eq = fAssemble->fCoarseEqs[isub][ic];
			mineq = mineq > eq ? eq : mineq;
		}
		for (ic=0; ic<nc; ic++) {
			int64_t eq = fAssemble->fCoarseEqs[isub][ic];
			if(skyline[eq] > mineq) skyline[eq] = mineq;
		}
	}
	/* Computing K(c) */
	TPZMatrix<TVar> *coarse = new TPZSkylMatrix<TVar>(fNumCoarse,skyline);
#ifdef PZDEBUG
	{
		TPZFMatrix<TVar> coarse2(*coarse);
		for (isub=0; isub<nsub; isub++) {
			int64_t nc = fAssemble->fCoarseEqs[isub].NElements();
			int64_t ic;
			for (ic=0; ic<nc; ic++) {
				int64_t ieq = fAssemble->fCoarseEqs[isub][ic];
				int64_t jc;
				for (jc=0; jc<nc; jc++) {
					int64_t jeq = fAssemble->fCoarseEqs[isub][jc];
					coarse2(ieq,jeq) = 1.;
				}
			}
			
		}
		VisualMatrix(coarse2,"CoarseMatrix.vtk");
	}
#endif
	typename std::list<TPZAutoPointer<TSubStruct> >::iterator it;
	int count = 0;
	for(it= fGlobal.begin(); it != fGlobal.end(); it++,count++)
	{
		(*it)->Contribute_Kc(*coarse,fAssemble->fCoarseEqs[count]);
	}
	fCoarse = new TPZStepSolver<TVar>(coarse);
}

template<class TVar, class TSubStruct>
void TPZDohrPrecond<TVar, TSubStruct>::ComputeV1(const TPZFMatrix<TVar> &x, TPZFMatrix<TVar> &v1) const
{
	/* Computing r(c) */
	TPZFMatrix<TVar> CoarseResidual(fNumCoarse,x.Cols());
	CoarseResidual.Zero();
	typename std::list<TPZAutoPointer<TSubStruct> >::const_iterator it;
	
	int isub = 0;
	for(it= fGlobal.begin(); it != fGlobal.end(); it++, isub++) {
		TPZFMatrix<TVar> xloc, CoarseResidual_local;
		fAssemble->Extract(isub,x,xloc);
		//		(*it)->LoadWeightedResidual(xloc);
		(*it)->Contribute_rc_local(xloc,CoarseResidual_local);
		fAssemble->AssembleCoarse(isub,CoarseResidual_local,CoarseResidual);
	}
#ifdef PZ_LOG
	{
		std::stringstream sout;
		CoarseResidual.Print("Coarse Residual",sout);
		if (logger.isDebugEnabled())
		{
			LOGPZ_DEBUG(logger, sout.str());
		}
	}
#endif
	/* Computing K(c)_inverted*r(c) and stores it in "product" */
	fCoarse->SetDirect(ELDLt);
	//Dado global
	TPZFMatrix<TVar> CoarseSolution(fNumCoarse,x.Cols());
	fCoarse->Solve(CoarseResidual,CoarseSolution);
#ifdef PZ_LOG
	{
		std::stringstream sout;
		CoarseSolution.Print("CoarseSolution",sout);
		if (logger.isDebugEnabled())
		{
			LOGPZ_DEBUG(logger, sout.str());
		}
	}
#endif
	isub=0;
	//Criar tarefa que execute a distribuicao de cada elemento do fGlobal
	for(it= fGlobal.begin(); it != fGlobal.end(); it++,isub++)
	{
		// Gerenciamento Global->Local sobre o product
		//product é administrado pelo DM mas permanece no processador 0
		// tarefa separada, expansao da solucao coarse
		TPZFMatrix<TVar> v1_local,CoarseSolution_local;
		fAssemble->ExtractCoarse(isub,CoarseSolution,CoarseSolution_local);
		(*it)->Contribute_v1_local(v1_local,CoarseSolution_local);
		
		fAssemble->Assemble(isub,v1_local,v1);
	}
#ifdef PZ_LOG
	{
		std::stringstream sout;
		v1.Print("v1 vector",sout);
		if (loggerv1v2.isDebugEnabled())
		{
			LOGPZ_DEBUG(loggerv1v2, sout.str())
		}
	}
#endif
}

template<class TVar, class TSubStruct>
void TPZDohrPrecond<TVar, TSubStruct>::ComputeV2(const TPZFMatrix<TVar> &x, TPZFMatrix<TVar> &v2) const
{
	
	typename std::list<TPZAutoPointer<TSubStruct> >::const_iterator it;
	
	int isub=0;
	//Criar tarefa que execute a distribuicao de cada elemento do fGlobal
	for(it= fGlobal.begin(); it != fGlobal.end(); it++,isub++)
	{
		// contribute v2 deve ser uma tarefa inicializada mais cedo
		TPZFNMatrix<100,TVar> Residual_local,v2_local;
		fAssemble->Extract(isub,x,Residual_local);
		(*it)->Contribute_v2_local(Residual_local,v2_local);
#ifdef PZ_LOG
		{
			std::stringstream sout;
			sout << "Substructure " << isub << std::endl;
			Residual_local.Print("Residual local",sout);
            v2_local.Print("v2_local",sout);
			if (logger.isDebugEnabled())
			{
				LOGPZ_DEBUG(logger, sout.str());
			}
		}
#endif
		//		v2_local += v1_local;
		fAssemble->Assemble(isub,v2_local,v2);
	}
#ifdef PZ_LOG
	{
		std::stringstream sout;
		v2.Print("v2 vector",sout);
		if (loggerv1v2.isDebugEnabled())
		{
			LOGPZ_DEBUG(loggerv1v2, sout.str())
		}
	}
#endif
}

template<class TVar, class TSubStruct>
void *TPZDohrPrecondV2SubDataList<TVar,TSubStruct>::ThreadWork(void *voidptr)
{
	TPZDohrPrecondV2SubDataList<TVar,TSubStruct>	*myptr = (TPZDohrPrecondV2SubDataList<TVar,TSubStruct> *) voidptr;
	TPZDohrPrecondV2SubData<TVar,TSubStruct> data = myptr->PopItem();
	while (data.IsValid()) {
		data.fSubStructure->Contribute_v2_local(data.fInput_local,data.fv2_local->fAssembleData);
		myptr->fAssemblyStructure->AddItem(data.fv2_local);
		data = myptr->PopItem();
	}
	return voidptr;
}

/**
 * @brief Unpacks the object structure from a stream of bytes
 * @param buf The buffer containing the object in a packed form
 * @param context
 */
template<class TVar, class TSubStruct>
void TPZDohrPrecond<TVar, TSubStruct>::Read(TPZStream &buf, void *context )
{
    TPZMatrix<TVar>::Read(buf,context);
    TPZDohrMatrix<TVar,TSubStruct> *ptr = (TPZDohrMatrix<TVar,TSubStruct> *)(context);
    fAssemble = ptr->fAssembly;
    fGlobal = ptr->SubStructures();
    buf.Read(&fNumCoarse);
    buf.Read(&fNumThreads);
    fCoarse = dynamic_cast<TPZStepSolver<TVar> *>(TPZPersistenceManager::GetInstance(&buf));
}
/**
 * @brief Packs the object structure in a stream of bytes
 * @param buf Buffer which will receive the bytes
 * @param withclassid
 */
template<class TVar, class TSubStruct>
void TPZDohrPrecond<TVar, TSubStruct>::Write( TPZStream &buf, int withclassid ) const
{
    TPZMatrix<TVar>::Write(buf, withclassid);
    buf.Write(&fNumCoarse);
    buf.Write(&fNumThreads);
    TPZPersistenceManager::WritePointer(fCoarse, &buf);
}


template class TPZDohrPrecond<float,TPZDohrSubstruct<float> >;
template class TPZDohrPrecond<double,TPZDohrSubstruct<double> >;
template class TPZDohrPrecond<long double,TPZDohrSubstruct<long double> >;

template class TPZDohrPrecond<float, TPZDohrSubstructCondense<float> >;
template class TPZDohrPrecond<double, TPZDohrSubstructCondense<double> >;
template class TPZDohrPrecond<long double, TPZDohrSubstructCondense<long double> >;

//template class TPZDohrPrecond<std::complex<float>,TPZDohrSubstruct<std::complex<float> > >;
template class TPZDohrPrecond<std::complex<double>,TPZDohrSubstruct<std::complex<double> > >;
//template class TPZDohrPrecond<std::complex<long double>,TPZDohrSubstruct<std::complex<long double> > >;

//template class TPZDohrPrecond<std::complex<float>, TPZDohrSubstructCondense<std::complex<float> > >;
template class TPZDohrPrecond<std::complex<double>, TPZDohrSubstructCondense<std::complex<double> > >;
//template class TPZDohrPrecond<std::complex<long double>, TPZDohrSubstructCondense<std::complex<long double> > >;


#ifndef BORLAND

template class TPZRestoreClass<TPZDohrPrecond<float, TPZDohrSubstruct<float> >>;
template class TPZRestoreClass<TPZDohrPrecond<double, TPZDohrSubstruct<double> >>;
template class TPZRestoreClass<TPZDohrPrecond<long double, TPZDohrSubstruct<long double> >>;

template class TPZRestoreClass<TPZDohrPrecond<std::complex<double>, TPZDohrSubstruct<std::complex<double>> >>;

template class TPZRestoreClass<TPZDohrPrecond<float, TPZDohrSubstructCondense<float> >>;
template class TPZRestoreClass<TPZDohrPrecond<double, TPZDohrSubstructCondense<double> >>;
template class TPZRestoreClass<TPZDohrPrecond<long double, TPZDohrSubstructCondense<long double> >>;

template class TPZRestoreClass<TPZDohrPrecond<std::complex<double>, TPZDohrSubstructCondense<std::complex<double>> >>;

#endif
