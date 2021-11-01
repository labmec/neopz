/**
 * @file
 * @brief Contains the TPZTransform<> class which implements an affine transformation between points in parameter space.
 */

#ifndef TTRANSFORMH
#define TTRANSFORMH

#include "pzfmatrix.h"

#include "fad.h"

template<class T>
class TPZVec;

/**
 * @ingroup topologyutils
 * @brief Implements an affine transformation between points in parameter space. \ref topologyutils "Topology Utility"
 */
template<class T=REAL>
class TPZTransform : public TPZSavable {
	/** @brief Number of rows of the matrix associated with the transformation */
	int fRow;
	/** @brief Number of columns of the matrix associated with the transformation */
	int fCol;
	/** @brief Matrix used by multiplications */
	TPZFNMatrix<9,T> fMult;
	/** @brief Matrix used by sums */
	TPZFNMatrix<3,T> fSum;
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
    
//    /** @brief Rows of transformation matrix */
//    int Rows();
//    
//    /** @brief Cols of transformation matrix */
//    int Cols();
	
    friend class TPZTransform<REAL>;
    
    friend class TPZTransform<Fad<REAL> >;

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
	void Apply(const TPZVec<T> &vectorin,TPZVec<T> &vectorout);
	
	void PrintInputForm(std::ostream &out);
	
	/** @brief Compare the current transformation with t transformation considering a given tolerance */
	int CompareTransform(TPZTransform<T> &t,REAL tol = 1.e-6);
	
        int ClassId() const override{
            return Hash("TPZTransform");
        }
        
        void Read(TPZStream& buf, void* context) override { //ok
            buf.Read(&fRow);
            buf.Read(&fCol);
            fMult.Read(buf,context);
            fSum.Read(buf,context);
        }
        
        void Write(TPZStream& buf, int withclassid) const override{ //ok
            buf.Write(&fRow);
            buf.Write(&fCol);
            fMult.Write(buf,withclassid);
            fSum.Write(buf,withclassid);
        }
	
};

/** @brief Overload operator << to write transform data */
template<class T>
inline std::ostream &operator<<(std::ostream &out, const TPZTransform<T> &tr)
{
	out << "mult = " << tr.Mult() << "\tsum = " << tr.Sum();
	return out;
}

#endif
