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
	
	int fRow,fCol;     //Matrix dimensions
	TPZFMatrix<REAL> fMult;	// multiplication matrix
	TPZFMatrix<REAL> fSum;		// matrix used to sum
	REAL fStore[12];	// storage the matrix objects use to avoid
	// dynamic memory allocation
public:
	
	TPZTransform(int dim);//square matrix
	
	TPZTransform();
	
	TPZTransform(int fRow,int fCol);
	
	TPZTransform(const TPZTransform &tr);
	
	~TPZTransform();
	
	TPZTransform &operator=(const TPZTransform &t);
	
	const TPZFMatrix<REAL>  & Mult() const {return fMult;}
	
	const TPZFMatrix<REAL>  & Sum() const {return fSum;}
	
	TPZFMatrix<REAL>  & Mult() {return fMult;}
	
	TPZFMatrix<REAL>  & Sum()  {return fSum;}
	
	/**Sets the transformation matrices*/
	void SetMatrix(TPZFMatrix<REAL> &mult,TPZFMatrix<REAL> &sum);
	
	/**Multiply the transformation object to the right with right*/
	TPZTransform Multiply(TPZTransform &right);
	
	/**Transforms the vector*/
	void Apply(TPZVec<REAL> &vectorin,TPZVec<REAL> &vectorout);
	
	void PrintInputForm(std::ostream &out);
	
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
