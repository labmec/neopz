/**
 * \file
 * @brief DEPRECATED CLASS. Contains the base class to define Coordinate Systems.
 */
//HEADER FILE FOR CLASS TPZCosys

#ifndef COSYSTHPP
#define COSYSTHPP


#include <math.h>
#include "pzreal.h"
#include "pzfmatrix.h"

enum gCosysType {reference,cartesian,cylindric,esferic};

class TPZCartsys;
/**
 * @deprecated DEPRECATED root coordinate system CLASS.
 * @brief Root class for the coordinate system. \n
 * Root class for the coordinate system
 */
class TPZCosys {
	
public:
	
	/** @brief Constructor without arguments */
	TPZCosys();
	
	/** @brief Copy constructor */
	TPZCosys(const TPZCosys &cp);
	
	
	/**
	 * @brief Constructor - Create one object from other (reference coordinate system) 
	 * @param num index of the coordinate system
	 * @param ref reference coordinate system
	 */
	TPZCosys(int num,TPZCartsys* ref = NULL);
	
	/** @brief Destructor - must be implemented at each child class */
	virtual  ~TPZCosys();
	
	/**
	 * @brief Resets all coordinate system data \n
	 * Sets-> tr=I ; ref=NULL ; origin=0.
	 */
	void Reset();
	
	/**
	 * @brief Verify if the parameters of the object have consistant values and adjusts \n
	 * their values to the range of the coordinate system
	 * @param co has the coordinate of the point to verify
	 */
	virtual void VerifyRange(TPZFMatrix & co)=0;
	
	/**
	 * @brief Calculates the transformation gradient Gradx given Gradient of  the point X 
	 * @note looking at the code, it doesnt seem to be correct, \n
	 * please send a mail message to ecrylo@uol.com.br to correct it
	 */
	virtual void TransformGradient(TPZVec<REAL> &X, TPZFMatrix &GradX, TPZVec<REAL> &x, TPZFMatrix &Gradx, TPZCosys *dest=0)=0;
	
	/** @brief Returns type of coordinate system */
	virtual int Type() = 0;
	
	/** @brief Returns the coordinate system identifier number */
	int Id(){return fNumber;}
    
	/**
	 * @brief Prints system coordinate data to especified file   "out " 
	 * @param out device to print
	 */
	void Print(std::ostream& out = std::cout);
	
	/**
	 * @brief The org parameter is expressed in reference coordinate system
	 * @param co new reference coordinate system
	 */
	void SetReference(TPZCartsys *co) {fReference = co;}
	
	/**
	 * @brief Transforms point from this system to the reference system
	 * @param x point to transfer to reference coordinate system
	 */
	//    virtual void ToReference(REAL point[3]) = 0;
	virtual void ToReference(TPZVec<REAL> &x) = 0;
	
	/**
	 * @brief Transforms point from reference system to this system
	 * @param x point to transfer from reference coordinate system
	 */
	//  virtual void FromReference(REAL point[3]) = 0;
	virtual void FromReference(TPZVec<REAL> &x) = 0;
    
	/**
	 * @brief Transforms point from current system to the reference "ref " system
	 * @param x new coordinate system origin
	 * @param ref new reference coordinate system
	 */
    void ToSpecific(TPZVec<REAL> &x, TPZCartsys *ref);
	//void ToSpecific(REAL x[3], TPZCosys *ref);
	
	/**
	 * @brief Transforms point from this system to the root reference
	 * @param point point to transfer to the global coordinate system
	 */
    //virtual void ToGlobal(REAL point[3]);
	virtual void ToGlobal(TPZVec<REAL> & point);
	
	/**
	 * @brief Transforms point from global to this system
	 * @param point point to transfer from the global coordinate system 
	 */
	//virtual void FromGlobal(REAL point[3]);
	virtual void FromGlobal(TPZVec<REAL> & point);

	/**
	 * @brief Returns the pointer to reference coordinate system 
	 */
	TPZCartsys * Reference(){return fReference;}
	
protected:
	
	/**
	 * Calculates the normalised vector 
	 */
	//void Normalise(REAL point[3]);
	
	/**
	 * Calculates the determinant of the matrix t 
	 */
	//REAL Determinant(REAL t[3][3]);
	
	/**coordinate system number*/
	int		fNumber;
	
	/**origin of the coordinate system, expressed in
	 * cartesian coordinates system*/
	//REAL	fOrigin[3];
	
	/**
	 * @brief Pointer to reference coordinate system if exist 
	 */
	TPZCartsys	*fReference;
	
	/**
	 * Rotation tensor -
	 * Composed by reference system cosine vectors
	 *  in cartesian coordinate system 
	 */
	//REAL	fTr[3][3];
	
};

#endif
