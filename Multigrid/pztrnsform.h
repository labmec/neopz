/**
 * @file
 * @brief Contains the TPZTransform<> class which implements an affine transformation between points in parameter space.
 */

#ifndef TTRANSFORMH
#define TTRANSFORMH

#include "pzfmatrix.h"

#ifdef _AUTODIFF
#include "fad.h"
#endif

template<class T>
class TPZVec;

/**
 * @ingroup topologyutils
 * @brief Implements an affine transformation between points in parameter space. \ref topologyutils "Topology Utility"
 */
template<class T=REAL>
class TPZTransform {
	/** @brief Number of rows of the matrix associated with the transformation */
	int fRow;
	/** @brief Number of columns of the matrix associated with the transformation */
	int fCol;
	/** @brief Matrix used by multiplications */
	TPZFNMatrix<3,T> fMult;
	/** @brief Matrix used by sums */
	TPZFNMatrix<9,T> fSum;
	/** @brief Storage the matrix objects use to avoid */
	// dynamic memory allocation
public:
	/** @brief Constructor of the transformation with square matrix */
	TPZTransform(int dim);
	/** @brief Default constructor */
	TPZTransform();
	/** @brief Constructor of the transformation with rectangular matrix */
	TPZTransform(int fRow,int fCol);
	/** @brief Copy constructor */
	TPZTransform(const TPZTransform<T> &tr);
	/** @brief Default destructor */
	~TPZTransform();
	
    friend class TPZTransform<REAL>;
    
#ifdef _AUTODIFF
    friend class TPZTransform<Fad<REAL> >;
#endif
    
	/** @brief Overloading equal operator for transformation */
	TPZTransform<T> &operator=(const TPZTransform<T> &t);
	
    /** @brief Create a copy form a real transformation */
    void CopyFrom(const TPZTransform<REAL> &cp);
    
	const TPZFMatrix<T>  & Mult() const {return fMult;}
	
	const TPZFMatrix<T>  & Sum() const {return fSum;}
	
	TPZFMatrix<T>  & Mult() {return fMult;}
	
	TPZFMatrix<T>  & Sum()  {return fSum;}
	
	/** @brief Sets the transformation matrices */
	void SetMatrix(TPZFMatrix<T> &mult,TPZFMatrix<T> &sum);
	
	/** @brief Multiply the transformation object (to the right) with right (Multiplying matrices) */
	TPZTransform<T> Multiply(TPZTransform<T> &right);
	
	/** @brief Transforms the vector */
	void Apply(TPZVec<T> &vectorin,TPZVec<T> &vectorout);
	
	void PrintInputForm(std::ostream &out);
	
	/** @brief Compare the current transformation with t transformation considering a given tolerance */
	int Compare(TPZTransform<T> &t,REAL tol = 1.e-6);
	
	void Read(TPZStream &buf);
    
	void Write(TPZStream &buf);
	
};

/** @brief Overload operator << to write transform data */
template<class T>
inline std::ostream &operator<<(std::ostream &out, const TPZTransform<T> &tr)
{
	out << "mult = " << tr.Mult() << " sum = " << tr.Sum();
	return out;
}

#endif
