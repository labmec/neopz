/**
 * @file
 * @brief Contains the TPZDohrMatrix class which implements a matrix divided into substructures. \n
 * Also contains the TPZDohrThreadMultData and TPZDohrThreadMultList structs.
 * @author Philippe Devloo
 * @since 2006
 */

#ifndef TPZDOHRMATRIX_H
#define TPZDOHRMATRIX_H

#include "pzmatrix.h"
#include <list>
#include <sstream>
#include "tpzautopointer.h"
#include "tpzdohrsubstruct.h"
#include "tpzdohrsubstructCondense.h"
#include "tpzdohrassembly.h"

#include "tpzdohrassemblelist.h"


/**
 * @brief Implements a matrix divided into substructures. \ref matrix "Matrix" \ref substructure "Sub structure"
 * @ingroup substructure matrix
 * @author Philippe Devloo
 */
template <class TVar, class TSubStruct>
class TPZDohrMatrix : public TPZMatrix<TVar>
{
public:
	/** @brief The matrix class is a placeholder for a list of substructures */
	typedef typename std::list<TPZAutoPointer<TSubStruct> > SubsList;
private:
	SubsList fGlobal;
	
	int fNumCoarse; //n(c)
	
	/** @brief Number of threads that will be used during the matrix vector multiplication */
	int fNumThreads;

	int64_t Size() const override;
  TVar* &Elem() override;
  const TVar* Elem() const override;
	
public:
	
	TPZAutoPointer<TPZDohrAssembly<TVar> > fAssembly;
	
	TPZDohrMatrix(TPZAutoPointer<TPZDohrAssembly<TVar> > dohrassembly);
	
	TPZDohrMatrix() : TPZMatrix<TVar>()
    {
        
    }
	
	TPZDohrMatrix(const TPZDohrMatrix &cp) : fGlobal(cp.fGlobal), fNumCoarse(cp.fNumCoarse), fNumThreads(cp.fNumThreads),
	fAssembly(cp.fAssembly)
	{
	}
	
	//	CLONEDEF(TPZDohrMatrix)
	inline TPZDohrMatrix*NewMatrix() const override {return new TPZDohrMatrix{};}
	virtual TPZMatrix<TVar>*Clone() const  override { return new TPZDohrMatrix(*this); }
	
	~TPZDohrMatrix();

	void CopyFrom(const TPZMatrix<TVar> *  mat) override        
  {                                                           
    auto *from = dynamic_cast<const TPZDohrMatrix<TVar,TSubStruct> *>(mat);                
    if (from) {                                               
      *this = *from;                                          
    }                                                         
    else                                                      
      {                                                       
        PZError<<__PRETTY_FUNCTION__;                         
        PZError<<"\nERROR: Called with incompatible type\n."; 
        PZError<<"Aborting...\n";                             
        DebugStop();                                          
      }                                                       
  }
	const SubsList &SubStructures() const
	{
		return fGlobal;
	}
	
	int NumCoarse() const
	{
		return fNumCoarse;
	}
	
	int NumThreads() const
	{
		return fNumThreads;
	}
	
	void SetNumThreads(int numthreads)
	{
		fNumThreads = numthreads;
	}
	
	/** @brief Just a method for tests */
	TPZAutoPointer<TSubStruct> GetFirstSub() {
		return (*fGlobal.begin());
	}
	/** @brief Just a method for tests */
	void Print(const char *name, std::ostream& out,const MatrixOutputFormat form = EFormatted) const override
	{
		out << __PRETTY_FUNCTION__ << std::endl;
		out << name << std::endl;
		out << "Number of coarse equations " << fNumCoarse << std::endl;
		typename SubsList::const_iterator iter;
		for (iter=fGlobal.begin();iter!=fGlobal.end();iter++) {
			(*iter)->Print(out);
		}
	}
	
	void SetNumCornerEqs(int nc)
	{
		fNumCoarse = nc;
	}
	/** @brief Initialize the necessary datastructures */
	void Initialize();
	/** @brief It adds a substruct */
	void AddSubstruct(TPZAutoPointer<TSubStruct> substruct)
	{
		fGlobal.push_back(substruct);
	}
	
	/**
	 * @brief It computes z = beta * y + alpha * opt(this)*x but z and x can not overlap in memory.
	 * @param x Is x on the above operation
	 * @param y Is y on the above operation
	 * @param z Is z on the above operation
	 * @param alpha Is alpha on the above operation
	 * @param beta Is beta on the above operation
	 * @param opt Indicates if is Transpose or not
	 */
    
    /** The only method any matrix class needs to implement with TBB */
	virtual void MultAddTBB(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,const TVar alpha,const TVar beta,const int opt) const;
    
	/** The only method any matrix class needs to implement */
	virtual void MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,const TVar alpha,const TVar beta,const int opt) const override;
	
	/** @brief Adjust the residual to zero the residual of the internal connects */
	void AdjustResidual(TPZFMatrix<TVar> &res);
	
	/** @brief Add the solution corresponding to the internal residual */
	void AddInternalSolution(TPZFMatrix<TVar> &solution);
    
	/**
	 * @name TPZSavable
	 * Methods which would make TPZMatrix<TVar>compliant with TPZSavable
	 */
	//@{
	/**
	 * @brief Unpacks the object structure from a stream of bytes
	 * @param buf The buffer containing the object in a packed form
	 * @param context
	 */
	void Read(TPZStream &buf, void *context) override;
	/**
	 * @brief Packs the object structure in a stream of bytes
	 * @param buf Buffer which will receive the bytes
	 * @param withclassid
	 */
	void Write(TPZStream &buf, int withclassid) const override;
    
    /** @brief Routines to send and receive messages */
	public:
int ClassId() const override;

};

template <class TVar, class TSubStruct>
int TPZDohrMatrix<TVar, TSubStruct>::ClassId() const{
    return Hash("TPZDohrMatrix") ^ TPZMatrix<TVar>::ClassId() << 1 ^ TSubStruct().ClassId() << 2;
}

/**
 * @ingroup substructure
 * @brief .. . \ref substructure "Sub structure"
 */
template <class TSubStruct>
struct TPZDohrThreadMultData
{
	/** @brief Default constructor */
	TPZDohrThreadMultData() : fisub(-1), fSub(0)
	{
	}
	TPZDohrThreadMultData(int isub, TPZAutoPointer<TSubStruct> submesh) : fisub(isub), fSub(submesh)
	{
	}
	/** @brief Copy constructor */
	TPZDohrThreadMultData(const TPZDohrThreadMultData<TSubStruct> &cp) : fisub(cp.fisub), fSub(cp.fSub)
	{
	}
	/** @brief Implement the attribution operator. */
	TPZDohrThreadMultData<TSubStruct> &operator=(const TPZDohrThreadMultData<TSubStruct> &cp)
	{
		fisub = cp.fisub;
		fSub = cp.fSub;
		return *this;
	}
	int fisub;
	TPZAutoPointer<TSubStruct> fSub;
	
	bool IsValid()
	{
		return (fisub >= 0);
	}
};

/**
 * @ingroup substructure
 * @brief .. . \ref substructure "Sub structure"
 */
template <class TVar, class TSubStruct>
struct TPZDohrThreadMultList
{
	/** @brief The vector with which we will multiply */
	const TPZFMatrix<TVar> *fInput;
	/** @brief Scalar multiplication factor */
	TVar fAlpha;
	/** @brief Mutex which will enable the access protection of the list */
	std::mutex fAccessLock;
	/** @brief The data structure which defines the assemble destinations */
	TPZAutoPointer<TPZDohrAssembly<TVar> > fAssembly;
	/** @brief The list of data objects which need to treated by the threads */
	std::list<TPZDohrThreadMultData<TSubStruct> > fWork;
	/** @brief The local contribution to the v2 vector */
	TPZAutoPointer<TPZDohrAssembleList<TVar> > fAssemblyStructure;
	
	TPZDohrThreadMultList(const TPZFMatrix<TVar> &x, TVar alpha, TPZAutoPointer<TPZDohrAssembly<TVar> > assembly, TPZAutoPointer<TPZDohrAssembleList<TVar> > &assemblestruct) : fInput(&x), fAlpha(alpha),
	fAssembly(assembly), fAssemblyStructure(assemblestruct), fAccessLock()
	{
	}
	~TPZDohrThreadMultList()
	{
	}
	
	/** @brief The procedure which executes the lengthy process */
	static void *ThreadWork(void *voidptr);
	/** @brief Interface to add items in a thread safe way */
	void AddItem(TPZDohrThreadMultData<TSubStruct> &data)
	{
    std::scoped_lock<std::mutex> lock(fAccessLock);
		fWork.push_back(data);
	}
	/** @brief Interface to pop an item in a thread safe way */
	TPZDohrThreadMultData<TSubStruct> PopItem()
	{
		TPZDohrThreadMultData<TSubStruct> result;
    std::scoped_lock<std::mutex> lock(fAccessLock);
		if (fWork.size()) {
			result = *fWork.begin();
			fWork.pop_front();
		}
		return result;
	}
};

#endif
