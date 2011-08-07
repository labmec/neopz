/**
 * @file
 * @brief Contains declaration of TPZReferredCompEl class which generates computational elements.
 */
//$Id: pzreferredcompel.h,v 1.11 2011-05-11 02:50:03 phil Exp $

#ifndef PZSPECIAL
#define PZSPECIAL

class TPZElementMatrix;
class TPZCompEl;
class TPZInterpolatedElement;
class TPZCompElDisc;
class TPZGeoEl;
class TPZCompMesh;
class TPZFMatrix;
#include "pzvec.h"
#include "pzmanvector.h"
#include <map>
#include <set>

/**
 * @ingroup CompElement
 * @brief Template to generate computational elements. \ref CompElement "Computational Element"
 */
template<class TCOMPEL>
class TPZReferredCompEl : public TCOMPEL {
public:
	
	/** @brief Class constructor */
	TPZReferredCompEl(TPZCompMesh &mesh, TPZGeoEl *gel, int &index);
	
	TPZReferredCompEl();
	
	TPZReferredCompEl(TPZCompMesh &mesh, const TPZReferredCompEl<TCOMPEL> &copy);
	
	TPZReferredCompEl(TPZCompMesh &mesh,
                      const TPZReferredCompEl<TCOMPEL> &copy,
                      std::map<int,int> & gl2lcConMap,
                      std::map<int,int> & gl2lcElMap);
	
	/** @brief Class destructor */
	~TPZReferredCompEl();
	
	/** @brief Set create function in TPZCompMesh to create elements of this type */
	virtual void SetCreateFunctions();
	
	/** @brief Returns referred element of this */
	TPZCompEl * ReferredElement();
	
	/**
	 * @brief Computes solution and its derivatives in local coordinate qsi
	 * @param qsi master element coordinate
	 * @param phi matrix containing shape functions compute in qsi point
	 * @param dphix matrix containing the derivatives of shape functions in the direction of the axes
	 * @param axes direction of the derivatives
	 * @param sol finite element solution
	 * @param dsol solution derivatives
	 */
	virtual void ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix &phi, TPZFMatrix &dphix,
								 const TPZFMatrix &axes, TPZVec<REAL> &sol, TPZFMatrix &dsol);
	
	/**
	 * @brief Computes solution and its derivatives in the local coordinate qsi.
	 * @param qsi master element coordinate of the interface element
	 * @param sol finite element solution
	 * @param dsol solution derivatives
	 * @param axes axes associated with the derivative of the solution
	 * @param leftsol finite element solution
	 * @param dleftsol solution derivatives
	 * @param leftaxes axes associated with the left solution
	 * @param rightsol finite element solution
	 * @param drightsol solution derivatives
	 * @param rightaxes axes associated with the right solution
	 */
	/**
	 * This method will function for both volumetric and interface elements
	 */
	virtual void ComputeSolution(TPZVec<REAL> &qsi,
								 TPZVec<REAL> &normal,
								 TPZVec<REAL> &leftsol, TPZFMatrix &dleftsol,TPZFMatrix &leftaxes,
								 TPZVec<REAL> &rightsol, TPZFMatrix &drightsol,TPZFMatrix &rightaxes);
	
	/**
	 * @brief Prints element data
	 * @param out indicates the device where the data will be printed
	 */
	virtual void Print(std::ostream & out = std::cout) const;
	
protected:
	
	/** @brief Append solution of the referred element. */
	void AppendOtherSolution(TPZVec<REAL> &qsi, TPZVec<REAL> &sol,
							 TPZFMatrix &dsol,  TPZFMatrix &axes);
	
	/** @brief Append solution of the referred element. */
	void AppendOtherSolution(TPZVec<REAL> &qsi, TPZVec<REAL> &sol,
							 TPZFMatrix &dsol,  const TPZFMatrix &axes);
	
	/** @brief Append solution of the referred element. */
	void AppendOtherSolution(TPZVec<REAL> &qsi, TPZVec<REAL> &normal,
							 TPZVec<REAL> &leftsol, TPZFMatrix &dleftsol, TPZFMatrix &leftaxes,
							 TPZVec<REAL> &rightsol, TPZFMatrix &drightsol,TPZFMatrix &rightaxes);
};

/** Adjust the derivatives from one system of axes to the other */
void AdjustSolutionDerivatives(TPZFMatrix &dsolfrom, TPZFMatrix &axesfrom,
                               TPZFMatrix &dsolto, const TPZFMatrix &axesto);

TPZCompEl *CreateReferredDisc(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
TPZCompEl *CreateReferredPointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
TPZCompEl *CreateReferredLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
TPZCompEl *CreateReferredQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
TPZCompEl *CreateReferredTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
TPZCompEl *CreateReferredCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
TPZCompEl *CreateReferredPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
TPZCompEl *CreateReferredPyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
TPZCompEl *CreateReferredTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);

void Append(TPZVec<REAL> &u1, TPZVec<REAL> &u2, TPZVec<REAL> &u12);
void Append(TPZFMatrix &u1, TPZFMatrix &u2, TPZFMatrix &u12);
bool AreEqual(const TPZVec<REAL> &A, const TPZVec<REAL> &B, REAL tol = 1e-10);

#endif
