/**
 * @file
 * @brief Contains the TPZTransform class which implements an affine transformation between points in parameter space.
 */

#ifndef TTRANSFORMH
#define TTRANSFORMH

#include "pzfmatrix.h"

template<class T>
class TPZVec;

/**
 * @ingroup topologyutils
 * @brief Implements an affine transformation between points in parameter space. \ref topologyutils "Topology Utility"
 */
class TPZTransform {
	/** @brief Number of rows of the matrix associated with the transformation */
	int fRow;
	/** @brief Number of columns of the matrix associated with the transformation */
	int fCol;
	/** @brief Matrix used by multiplications */
	TPZFMatrix<REAL> fMult;
	/** @brief Matrix used by sums */
	TPZFMatrix<REAL> fSum;
	/** @brief Storage the matrix objects use to avoid */
	REAL fStore[12];
	// dynamic memory allocation
public:
	/** @brief Constructor of the transformation with square matrix */
	TPZTransform(int dim);
	/** @brief Default constructor */
	TPZTransform();
	/** @brief Constructor of the transformation with rectangular matrix */
	TPZTransform(int fRow,int fCol);
	/** @brief Copy constructor */
	TPZTransform(const TPZTransform &tr);
	/** @brief Default destructor */
	~TPZTransform();
	
	/** @brief Overloading equal operator for transformation */
	TPZTransform &operator=(const TPZTransform &t);
	
	const TPZFMatrix<REAL>  & Mult() const {return fMult;}
	
	const TPZFMatrix<REAL>  & Sum() const {return fSum;}
	
	TPZFMatrix<REAL>  & Mult() {return fMult;}
	
	TPZFMatrix<REAL>  & Sum()  {return fSum;}
	
	/** @brief Sets the transformation matrices */
	void SetMatrix(TPZFMatrix<REAL> &mult,TPZFMatrix<REAL> &sum);
	
	/** @brief Multiply the transformation object (to the right) with right (Multiplying matrices) */
	TPZTransform Multiply(TPZTransform &right);
	
	/** @brief Transforms the vector */
	void Apply(TPZVec<REAL> &vectorin,TPZVec<REAL> &vectorout);
	
	void PrintInputForm(std::ostream &out);
	
	/** @brief Compare the current transformation with t transformation considering a given tolerance */
	int Compare(TPZTransform &t,REAL tol = 1.e-6);
	
	void Read(TPZStream &buf);
    
	void Write(TPZStream &buf);
	
};

/** @brief Overload operator << to write transform data */
inline std::ostream &operator<<(std::ostream &out, const TPZTransform &tr)
{
	out << "mult = " << tr.Mult() << " sum = " << tr.Sum();
	return out;
}

#endif
