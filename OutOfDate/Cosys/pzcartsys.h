/**
 * @file
 * @brief DEPRECATED CLASS. Contains the class defining the Cartesian Coordinate System
 */
//HEADER FILE FOR CLASS COSYS
class TPZCosys;

#ifndef CARTSYSTHPP
#define CARTSYSTHPP

#include "pzcosys.h"
#include "pzvec.h"

/**
 * @deprecated DEPRECATED cartesian coordinate system CLASS.
 * @brief Defines the cartesian coordinate system class
 */
class  TPZCartsys : public TPZCosys {
	
public:
	
	/** @brief Default empty constructor */
	TPZCartsys();
	
	/** @brief Copy constructor */
	TPZCartsys(const TPZCartsys & cp);
	
	/**
	 * @brief Create object from one reference object and new origin 
	 * @param num coordinate system identificator
	 * @param ref reference coordinate system
	 * @param org pointer to the coordinate system origin - given in local coordinate system terms
	 * @param angles - rotation angles given in terms of euller angles
	 */
	TPZCartsys(int num, TPZCartsys * ref, TPZVec<REAL> * org, TPZVec<REAL> *angles);
	
	/**
	 * @brief Create object from one reference object and new origin 
	 * @param num coordinate system identificator
	 * @param ref reference coordinate system
	 * @param org pointer to the coordinate system origin - given in local coordinate system terms
	 * @param angles - rotation angles given in terms of rotation matrix 
	 */
	TPZCartsys(int num, TPZCartsys * ref, TPZVec<REAL> * org, TPZFMatrix *angles=NULL);
	
	/** Destructor */
	~TPZCartsys() ;
	
	/**
	 * @brief Returns the current coordinate system axes 
	 * @param t reference that will receive the coordinate system axes
	 */
	void GetAxes(TPZFMatrix &t);
	
	/**
	 * @brief Changes the current coordinate system axes giving two axes (vectors) 
	 * @param x new x axe
	 * @param z new z axe
	 */
	//void SetAxes(REAL x[3], REAL z[3]);
	void SetAxes(TPZVec<REAL>&x,TPZVec<REAL>&z);
	
	/**
	 * @brief Changes the current coordinate system axes 
	 * giving the Rotation Matrix 
	 * @param RotMat rotation matrix to set new coordinate system axes
	 */
	void SetAxes(TPZFMatrix &RotMat);
	
	/**
	 * @brief Changes the current coordinate system origin 
	 * @param org pointer to the new origin point
	 */
	void SetOrigin(TPZVec<REAL> *org);
	
	/** @brief Gives the object coordinate system type */
	int Type() {return cartesian;}
	
	/**
	 * @brief Calculates the transformation gradient Gradx given Gradient of the point X 
	 * @param X point that want to transform gradient
	 * @param GradX gradient of point X
	 * @param x will receive point X in terms of the destintation coordinate system
	 * @param Gradx will receive the tranformed gradient
	 * @param dest coordinate system to transform the gradient
	 */
	void TransformGradient(TPZVec<REAL> &X, TPZFMatrix &GradX, TPZVec<REAL> &x, TPZFMatrix &Gradx, TPZCosys *dest=0); 
	
	/** @brief Verify Range */
	/** None parameters must be verified in cartesian coordinates system */
	void VerifyRange(TPZFMatrix &){;}
	
	/**
	 * @brief Returns de the value of one coordinate given in reference system in actual system
	 * @param point point to transfer from reference coordinate system
	 */
	void FromReference (TPZVec<REAL> &point);
	
	/**
	 * @brief Returns de the value of one coordinate in reference coordinate system
	 * @param point point to transfer to reference coordinate system   
	 */
	void ToReference (TPZVec<REAL> &point) ;
	
	/**
	 * @brief Returns de the value of one coordinate given in reference system in actual system
	 * @param point point to transfer from reference coordinate system
	 */
	void RotateFromReference (TPZVec<REAL> &point);
	
	/**
	 * @brief Returns de the value of one coordinate in reference coordinate system
	 * @param point point to transfer to reference coordinate system
	 */
	void RotateToReference (TPZVec<REAL> &point) ;
	
	
	/** @brief Returns the rotation matrix */
	TPZFMatrix &Rotation () {return fTr;}
	
	void Print(std::ostream &out);
    double RotationAngle();
    void SetGlobal();
	
protected:
	/**
	 * @brief Rotation tensor \n
	 * Composed by reference system cosine vectors
	 * in cartesian coordinate system 
	 */
	TPZFMatrix fTr;
	
	/** @brief Origin of the coordinate system, expressed in cartesian reference coordinate system */
	TPZVec<REAL> fOrigin;
	
private:
	
	/**
	 * @brief Calculates the rotation matrix giving the  euller angles
	 * @param angles euller angles
	 */
	void CalcRotMatrix(TPZVec<REAL> &angles);
	
	/**
	 * @brief Calculates the normalised vector 
	 * @param point - point that defines the vector to be normalised
	 */
	void Normalise(TPZVec<REAL> &point);
	
	/**
	 * @brief Returns in the parameter norm the normal vector to the vectors vec1 and vec2 
	 * @param vec1 first vector that defines the plane whose want the normal
	 * @param vec2 second vector that defines the plane whose want the normal
	 * @param norm will receive the calculated normal vector
	 */
	void GetNormal(TPZVec<REAL> &vec1, TPZVec<REAL> &vec2, TPZVec<REAL> &norm);
	
};

#endif
