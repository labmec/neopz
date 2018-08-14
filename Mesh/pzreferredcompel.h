/**
 * @file
 * @brief Contains declaration of TPZReferredCompEl class which generates computational elements.
 */

#ifndef PZSPECIAL
#define PZSPECIAL

struct TPZElementMatrix;
class TPZCompEl;
#include "pzintel.h"
class TPZGeoEl;
class TPZCompMesh;
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzmaterialdata.h"
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
	TPZReferredCompEl(TPZCompMesh &mesh, TPZGeoEl *gel, int64_t &index);
	
	TPZReferredCompEl();
	
	TPZReferredCompEl(TPZCompMesh &mesh, const TPZReferredCompEl<TCOMPEL> &copy);
	
	TPZReferredCompEl(TPZCompMesh &mesh,
                      const TPZReferredCompEl<TCOMPEL> &copy,
                      std::map<int64_t,int64_t> & gl2lcConMap,
                      std::map<int64_t,int64_t> & gl2lcElMap);
	
	/** @brief Class destructor */
	~TPZReferredCompEl();
	
	/** @brief Set create function in TPZCompMesh to create elements of this type */
	virtual void SetCreateFunctions(TPZCompMesh *mesh);
	
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
//    virtual void ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphix,
//                                 const TPZFMatrix<REAL> &axes, TPZSolVec &sol, TPZGradSolVec &dsol);
	
    /**
     * @brief Computes solution and its derivatives in local coordinate qsi
     * @param qsi master element coordinate
     * @param phi matrix containing shape functions compute in qsi point
     * @param dphix matrix containing the derivatives of shape functions in the direction of the axes
     * @param axes direction of the derivatives
     * @param sol finite element solution
     * @param dsol solution derivatives
     */
    virtual void ComputeSolution(TPZVec<REAL> &qsi, TPZMaterialData &data);
    
	/**
	 * @brief Computes solution and its derivatives in the local coordinate qsi.
	 * @param qsi master element coordinate of the interface element
	 * @param normal unit normal vector
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
								 TPZSolVec &leftsol, TPZGradSolVec &dleftsol,TPZFMatrix<REAL> &leftaxes,
								 TPZSolVec &rightsol, TPZGradSolVec &drightsol,TPZFMatrix<REAL> &rightaxes);
	
	/**
	 * @brief Prints element data
	 * @param out indicates the device where the data will be printed
	 */
	virtual void Print(std::ostream & out = std::cout) const;
	
    public:
virtual int ClassId() const;
    public:
protected:
	
	/** @brief Append solution of the referred element. */
	void AppendOtherSolution(TPZVec<REAL> &qsi, TPZSolVec &sol,
							 TPZGradSolVec &dsol,  TPZFMatrix<REAL> &axes);
	
    /** @brief Append solution of the referred element. */
    void AppendOtherSolution(TPZVec<REAL> &qsi, TPZSolVec &sol);
    
	/** @brief Append solution of the referred element. */
	void AppendOtherSolution(TPZVec<REAL> &qsi, TPZSolVec &sol,
							 TPZGradSolVec &dsol,  const TPZFMatrix<REAL> &axes);
	
	/** @brief Append solution of the referred element. */
	void AppendOtherSolution(TPZVec<REAL> &qsi, TPZVec<REAL> &normal,
							 TPZSolVec &leftsol, TPZGradSolVec &dleftsol, TPZFMatrix<REAL> &leftaxes,
							 TPZSolVec &rightsol, TPZGradSolVec &drightsol,TPZFMatrix<REAL> &rightaxes);
};

/** @brief Adjust the derivatives from one system of axes to the other */
void AdjustSolutionDerivatives(TPZFMatrix<STATE> &dsolfrom, TPZFMatrix<REAL> &axesfrom,
                               TPZFMatrix<STATE> &dsolto, const TPZFMatrix<REAL> &axesto);

/** @brief Creates discontinuous referred computational element related with geometric element gel */
TPZCompEl *CreateReferredDisc(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates point referred computational element related with geometric element gel */
TPZCompEl *CreateReferredPointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates linear referred computational element related with geometric element gel */
TPZCompEl *CreateReferredLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates quadrilateral referred computational element related with geometric element gel */
TPZCompEl *CreateReferredQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates triangular referred computational element related with geometric element gel */
TPZCompEl *CreateReferredTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates cube referred computational element related with geometric element gel */
TPZCompEl *CreateReferredCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates prismal referred computational element related with geometric element gel */
TPZCompEl *CreateReferredPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates pyramidal referred computational element related with geometric element gel */
TPZCompEl *CreateReferredPyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates tetrahedral referred computational element related with geometric element gel */
TPZCompEl *CreateReferredTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);

/** @brief Append u2 vector after u1 vector in u12 vector */
template<class TVar>
void Append(TPZVec<TVar> &u1, TPZVec<TVar> &u2, TPZVec<TVar> &u12);
/** @brief Append u2 matrix following u1 matrix in u12 matrix */
/** Returns u12 = [u1][u2]. Then: \f$ u12.Rows = max(u1.Rows, u2.Rows) \f$ and \f$ u12.Cols = u1.Cols + u2.Cols \f$ */
template<class TVar>
void Append(TPZFMatrix<TVar> &u1, TPZFMatrix<TVar> &u2, TPZFMatrix<TVar> &u12);
/** @brief Returns true whether \f$ |Aij - Bij| < tol \f$ for all the entries of the matrices */
bool AreEqual(const TPZVec<REAL> &A, const TPZVec<REAL> &B, REAL tol = 1e-10);

#endif
