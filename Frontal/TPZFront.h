/**
 * @file
 * @brief Contains the TPZFront class which implements decomposition process of the frontal matrix.
 */

#ifndef TPZFRONT_H
#define TPZFRONT_H
#include <mutex>
#include <thread>
#include <vector>

#include "pzmatrix.h"
#include "pzstack.h"
#include "pzvec.h"
#include "TPZSemaphore.h"
#include "Hash/TPZHash.h"
template<class TVar>
class TPZEqnArray;

/** 
 * The Front matrix itself. \n
 * It is controled by TPZFrontMatrix.
 */
/**
 * @ingroup frontal
 * @brief Abstract class implements storage and decomposition process of the frontal matrix. \ref frontal "Frontal"
 */
template<class TVar>
class TPZFront : public TPZSavable {
	
public:
	
    ///** @brief Static main used for testing */	
	//static void main();
	
	int64_t NElements();
    /** @brief Simple destructor */
    virtual ~TPZFront();
    /** @brief Simple constructor */
    TPZFront();
    /** @brief Constructor with a initial size parameter */
	TPZFront(
			 int64_t GlobalSize //! Initial size of the Frontal Matrix
			 );
	
	TPZFront(const TPZFront<TVar> &cp);
	
    /** 
	 * @brief Decompose these equations in a symbolic way and store freed indexes in fFree 
	 * @param mineq Initial equation
	 * @param maxeq Final equation
	 */
    void SymbolicDecomposeEquations(int64_t mineq, int64_t maxeq);
	
	/** 
	 * @brief Add a contribution of a stiffness matrix using the indexes to compute the frontwidth 
	 * @param destinationindex Destination index of each element added
	 */
	void SymbolicAddKel(TPZVec < int64_t > & destinationindex);

	int Work() {
		return fWork;
	}
	
	/** @brief Indicate the first equation dedicated to rigid body modes */
	void SetNumRigidBodyModes(int nrigid)
	{
		fNextRigidBodyMode = fLocal.NElements()-nrigid;
		
		std::cout << " fNextRigidBody Mode neste ponto " << fNextRigidBodyMode<<std::endl;
	}
	
        public:
int ClassId() const override;
        
protected:
	int fWork;
private:    
	
    /**
     * @brief Sets the global equation as freed, allowing the space used by this equation to be used
	 * @param global Equation number to be freed
	 */
	/** By future assembly processes. 
     */
    void FreeGlobal(int64_t global);
    
	/** 
	 * @brief return a local index corresponding to a global equation number 
	 * @param global Global equation index which has a local indexation
     */
    int Local(int64_t global);
	
public:
	/** @brief Extracts the so far condensed matrix */
	virtual	void ExtractFrontMatrix(TPZFMatrix<TVar> &front) {
		std::cout << "TPZFront ExtractFrontMatrix should never be called\n";
		DebugStop();
	}
	
    /** @brief Returns decomposition type. \n Default LU*/
    DecomposeType GetDecomposeType() const
    {
        return fDecomposeType;
    }
    
    /// Set the decomposition type
    virtual void SetDecomposeType(DecomposeType dectype) = 0;

	/** @brief Return the number of equations in the condensed front matrix
	 * It would be equal to FrontSize if the front is compressed.
	 */
	int NonNullFrontSize() const{
		int maxeq = fLocal.NElements();
		int mineq = 0;
		for(mineq=0; mineq<maxeq; mineq++) if(fLocal[mineq] != -1) break;
		int numeq = maxeq-mineq;
		return numeq;
	}	
	
	/** @brief Returns the number of free equations */
	virtual int64_t NFree();
    /** Resets data structure */
	void Reset(int64_t GlobalSize=0);
	
    /** @brief It prints TPZFront data */
	void Print(const char *name, std::ostream& out) const;
	void PrintGlobal(const char *name, std::ostream& out = std::cout);
	
	/** @brief returns the actual front size */
	int FrontSize()
	{
		return fFront;
	}

protected:
    
    /** @brief Maximum size of the front */
    int fMaxFront;

    /**
     * @brief Global equation associated to each front equation.
	 */
	/**
	 * If we need a position in globalmatrix of a equation "i" in the frontmatrix \n
	 * then we can use fGlobal[i]. If the global equation "i" is not used \f$ then fGlobal[i]==-1 \f$
     */
    TPZManVector <int64_t> fGlobal;
	
    /** @brief Front equation to each global equation */
    /**
	 * If we need a position in frontmatrix of a global equation "i" \n
	 * then we can use fLocal[i]. If the global equation is not represented in the front then \f$ fLocal[i]==-1 \f$.
     */
    TPZVec<int64_t> fLocal;
	
    /** @brief Actual front size */
    int64_t fFront;
	
	/** @brief Equation where rigid body modes can be stored */
	int64_t fNextRigidBodyMode;
	
    /** @brief Colection of already decomposed equations still on the front */
    TPZStack <int> fFree;
	
    /** @brief Frontal matrix data */
    TPZVec<TVar> fData;
    
    /** @brief Expansion Ratio of frontal matrix */
    int fExpandRatio;

protected:
    
    /** @brief Used Decomposition method */
    DecomposeType fDecomposeType;
    

public:
	
	///struct para paralelizar a decomposicao da matriz
	struct STensorProductMTData{
		
		typedef std::pair<int,STensorProductMTData * > STensorProductThreadData;
		
	public:
		
		///dados para sincronizar o thread principal
		TPZSemaphore fWorkDoneSem;
		int fWorkDoneCount;
		std::mutex fMutexWorkDoneCS;
		
		///semaforos para sincronizar os threads de calculo
		TPZVec<TPZSemaphore> fWorkSem;
		
		///array de threads
		std::vector<std::thread> fThreads;//for now we cannot use TPZVec
		
		///vetores de operacao
		TPZVec<TVar> * fAuxVecCol, * fAuxVecRow;
        
        ///valor da diagonal
        TVar fDiagonal;
		
		///num threads
		int NThreads(){ return fThreads.size(); };
		
		//vec to storage
		TPZVec<STensorProductThreadData*> fThreadData;  
		
		///matriz TPZFront
		TPZFront<TVar> * fMat;
		
		bool fRunning;
		
		///construtor padrao
		STensorProductMTData(int nthreads, TPZFront<TVar> * frontMat){
			if(!frontMat) DebugStop();
			this->fMat = frontMat;
			this->fRunning = true;
			
			fWorkSem.Resize(nthreads);
			
      // threads must be initialised already with thread work
			// fThreads.resize(nthreads);
			fThreadData.Resize(nthreads);
			for(int i = 0; i < nthreads; i++){
				STensorProductThreadData * threadData = new STensorProductThreadData;
				threadData->first = i;
				threadData->second = this;
				fThreadData[i] = threadData;
        fThreads.push_back(std::thread(Execute, threadData));
			}
		}///construtor
		
		///destrutor
		~STensorProductMTData(){
			
			this->fRunning = false;
			const int nthreads = fWorkSem.NElements();
			for(int i = 0; i < nthreads; i++){
        fWorkSem[i].Post();
			}
			for(int i = 0; i < nthreads; i++){
        fThreads[i].join();
				delete fThreadData[i];
				fThreadData[i] = NULL;
			}
			
		}///destrutor
		
		static void *Execute(void *data){
			STensorProductThreadData * mypair = static_cast< STensorProductThreadData * >(data);
			if(!mypair) DebugStop();
			TPZFront<TVar> * mat = mypair->second->fMat;
			if(!mat) DebugStop();
			mat->TensorProductIJ(mypair->first,mypair->second);
			return NULL;
		}
		
		void WorkDone(){
			std::scoped_lock<std::mutex> lck(fMutexWorkDoneCS);
      
			fWorkDoneCount++;
			if(fWorkDoneCount == NThreads()){
        fWorkDoneSem.Post();
			}
		}
		
		void Run(TPZVec<TVar> &AuxVecCol, TPZVec<TVar> &AuxVecRow){
			this->fAuxVecCol = &AuxVecCol;
			this->fAuxVecRow = &AuxVecRow;
			this->fWorkDoneCount = 0;
			for(int i = 0; i < NThreads(); i++){
				fWorkSem[i].Post();
			}
      fWorkDoneSem.Wait();
		}
		
	};
	
protected:
	
	STensorProductMTData *fProductMTData;
	
public:
	
	void ProductTensorMTInitData(int nthreads){
		fProductMTData = new STensorProductMTData(nthreads, this);
	}///void
	
	void ProductTensorMTFinish(){
		delete this->fProductMTData;
	}///void
	
	void ProductTensorMT(TPZVec<TVar> &AuxVecCol, TPZVec<TVar> &AuxVecRow){
		this->fProductMTData->Run(AuxVecCol,AuxVecRow);
	}///void
	
public:
	
	///Faz o tensor product de fato
	virtual void TensorProductIJ(int ithread,typename TPZFront<TVar>::STensorProductMTData *data);
};

template<class TVar>
int TPZFront<TVar>::ClassId() const {
    return Hash("TPZFront") ^ ClassIdOrHash<TVar>() << 1;
}
#endif //TPZFRONT_H
