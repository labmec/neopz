/**
 * @file
 * @brief Contains declaration of TPZInterfaceElement class which computes the contribution over an interface between two discontinuous elements.
 */
//$Id: TPZInterfaceEl.h,v 1.60 2010-06-17 17:46:56 phil Exp $

#ifndef ELEMINTERFACEHH
#define ELEMINTERFACEHH

#include "pzcompel.h"
#include "pzinterpolationspace.h"
#include "pzgeoel.h"
#include "pzcmesh.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzreal.h"

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include "TPZCompElDisc.h"
#include "pzmaterialdata.h"


/**
 @brief Computes the contribution over an interface between two discontinuous elements. \ref CompElement "Computational Element"
 @ingroup CompElement
 */
class TPZInterfaceElement : public TPZCompEl {
	
	private :
	
	/** @brief Element the left of the normal a interface */
	TPZCompElSide fLeftElSide;
	
	/** @brief Element the right of the normal a interface */
	TPZCompElSide fRightElSide;
	
	/** @brief Normal to the face element */
	TPZManVector<REAL,3> fCenterNormal;
	
	/** @brief Informs the connect that this element is no longer connected to it. */
	void DecreaseElConnected();
	
	/** @brief Informs the connect that this element is connected to it. */
	void IncrementElConnected();
	
protected:
	
	/** @brief Initialize a material data and its attributes based on element dimension, number of state variables and material definitions
	 */
	void InitMaterialData(TPZMaterialData &data, TPZInterpolationSpace *left);
    
    /** @brief Initialize the material data with the geometric data of the interface element */
    void InitMaterialData(TPZMaterialData &data);
	
	/** @brief Compute and fill data with requested attributes for neighbouring element */
	void ComputeRequiredData(TPZMaterialData &data,
							 TPZInterpolationSpace *elem,
							 TPZVec<REAL> &IntPoint);
    
    /** @brief Compute the required geometric data for the interface element */
    void ComputeRequiredData(TPZMaterialData &data);
	
public:
	
	/**
	 * @brief Extract connects from element el.
	 */
	void GetConnects(TPZCompElSide &elside, TPZVec<TPZConnect*> &connects, TPZVec<int> &connectindex);
	
	/** @brief Compute solution at neighbour element in a given master coordinate qsi. It returns the axes
	 * at which respect derivatives are computed.
	 * @param [in] Neighbor
	 * @param [in] qsi
	 * @param [out] sol
	 * @param [out] dsol
	 * @param [out] NeighborAxes
	 */
	void NeighbourSolution(TPZCompElSide & Neighbor, TPZVec<REAL> & qsi, TPZSolVec &sol, TPZGradSolVec &dsol, TPZFMatrix<REAL> &NeighborAxes);
	
protected:
	
	/** @brief Check consistency of mapped qsi performed by method TPZInterfaceElement::MapQsi by
	 * comparing the X coordinate of qsi and the correspondent NeighIntPoint.
	 * 
	 * It return true if everything is ok or false otherwise.
	 */
public:
	bool CheckConsistencyOfMappedQsi(TPZCompElSide &Neighbor, TPZVec<REAL> &qsi, TPZVec<REAL>&NeighIntPoint);
	
	void ComputeSideTransform(TPZCompElSide &Neighbor, TPZTransform &transf);
protected:
	/** @brief Computes normal at qsi point */
	void ComputeNormal(TPZVec<REAL>&qsi, TPZVec<REAL> &normal);
	
	/** @brief Computes normal based on already computed axes matrix. */
	/** Axes has been computed in the desired qsi coordinate */
	void ComputeNormal(TPZFMatrix<REAL> &axes, TPZVec<REAL> &normal);
	
	/** @brief Computes normal for linear geometric elements. */
	/** For linear geometry the normal vector is constant. */
	void ComputeCenterNormal(TPZVec<REAL> &normal);
	
public:
	
	void InitializeElementMatrix(TPZElementMatrix &ef);
	void InitializeElementMatrix(TPZElementMatrix &ek, TPZElementMatrix &ef);
	
	/**
	 * @brief Maps qsi coordinate at this master element to qsi coordinate at neighbor master element.
	 * @param Neighbor [in] may be this->LeftElementSide() or this->RightElementSide()
	 * @param qsi [in] is the point at this element master
	 * @param NeighIntPoint [out] is the point at neighbor element master. X[qsi] is equal to X[NeighIntPoint]
	 */
	void MapQsi(TPZCompElSide &Neighbor, TPZVec<REAL> &qsi, TPZVec<REAL> &NeighIntPoint);
	
	enum CalcStiffOptions{ENone = -1, EStandard /*Deprecated*/ = 0, EPenalty, EContDisc,EReferred};
	
	/** @brief Constuctor to continuous and/or discontinuous neighbours. */
	TPZInterfaceElement(TPZCompMesh &mesh,TPZGeoEl *geo,int &index,TPZCompElSide & left, TPZCompElSide &right);
	
	/** @brief Simple copy constructor. */
	TPZInterfaceElement(TPZCompMesh &mesh, const TPZInterfaceElement &copy);

	/**
	 * @brief Clone constructor to a patch mesh
	 * @param mesh reference to a clone mesh
	 * @param copy element to be copied
	 * @param gl2lcConIdx map with connects
	 * @param gl2lcElIdx map with computational elements
	 */
	TPZInterfaceElement(TPZCompMesh &mesh,
						const TPZInterfaceElement &copy,
						std::map<int,int> &gl2lcConIdx,
						std::map<int,int> &gl2lcElIdx);

	/** @brief Copy constructor with specified index */
	TPZInterfaceElement(TPZCompMesh &mesh, const TPZInterfaceElement &copy, int &index);
	
	/** @brief Empty constructor. */
	TPZInterfaceElement();
	
	/** @brief Default TPZCompEl constructor. SetLeftRightElements must be called before any computation. */
	TPZInterfaceElement(TPZCompMesh &mesh,TPZGeoEl *geo,int &index);
	
	/** @brief Destructor */
	~TPZInterfaceElement();
	
	virtual int IsInterface() { return 1; }
	
	
	/** @brief Set neighbors. */
	void SetLeftRightElements(TPZCompElSide & left, TPZCompElSide & right);
	
	/** @brief Makes a clone of this */
	virtual TPZCompEl *Clone(TPZCompMesh &mesh) const {
		return new TPZInterfaceElement(mesh, *this);
	}
	
	/** @see class TPZCompEl */
	virtual TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,std::map<int,int> &gl2lcConMap, std::map<int,int> &gl2lcElMap) const
	{
		return new TPZInterfaceElement(mesh, *this, gl2lcConMap,gl2lcElMap);
	}
	
	/** @brief Method used in TPZAgglomerateElement::CreateAgglomerateMesh */
	TPZCompEl * CloneInterface(TPZCompMesh &aggmesh,int &index, /*TPZCompElDisc **/TPZCompElSide & left, /*TPZCompElDisc **/ TPZCompElSide &right) const;
	
	/** @brief Identifies the elements of left and right volume of the interface */
	void VolumeEls(TPZCompEl &thirdel);
	
	/** @brief Returns the right element from the element interface */
	TPZCompEl *RightElement() const {
		return fRightElSide.Element();
	}
	
	/** @brief Returns the left element from the element interface */
	TPZCompEl *LeftElement() const {
		return fLeftElSide.Element();
	}
	
	/** @brief Returns left neighbor */
	TPZCompElSide &LeftElementSide(){ return this->fLeftElSide; }
	
	/** @brief Returns right neighbor */
	TPZCompElSide &RightElementSide(){ return this->fRightElSide; }
	
	/** @brief Returns the normal of this interface which goes from left to right neighbors */
	void CenterNormal(TPZVec<REAL> &CenterNormal) const;
	
	/** @brief Returns normal based on already computed axes matrix. */
	/**
	 * Axes has been computed in the desired qsi coordinate
	 * If geometric element has LinearMapping the CenterNormal is returned
	 */
	void Normal(TPZFMatrix<REAL> &axes, TPZVec<REAL> &normal);
	
	/** @brief Returns normal at qsi point */
	/** If geometric element has LinearMapping the CenterNormal is returned */
	void Normal(TPZVec<REAL>&qsi, TPZVec<REAL> &normal);
	
	/** @brief Returns the number from connectivities of the element */
	virtual int NConnects() const;
	
	/** @brief Returns the number from connectivities of the element related to right neighbour */
	int NRightConnects() const;
	
	/** @brief Returns the number from connectivities of the element related to left neighbour */
	int NLeftConnects() const;
	
	/** @brief Its return the connects of the left and right element associates */
	int ConnectIndex(int i) const;
	
	/** @brief This function should not be called */
	void SetConnectIndex(int node, int index);
	
	/** @brief Returns the dimension from the element interface */
	int Dimension() const {
		return this->Reference()->Dimension();
	}
	
	/** @brief Type of the element */
	MElementType Type() { return EInterface; }
	
	/**
	 * @brief Loads the solution within the internal data structure of the element
	 */
	/**
	 * Is used to initialize the solution of connect objects with dependency
	 * Is also used to load the solution within SuperElements
	 */
	virtual void LoadSolution(){
		//NOTHING TO BE DONE HERE
	}
	
	/**
	 * @brief CalcStiff computes the element stiffness matrix and right hand side
	 * @param ek element matrix
	 * @param ef element right hand side
	 */
	virtual void CalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef);
	
	/**
	 * @brief CalcResidual only computes the element residual
	 * @param ef element residual
	 */
	virtual void CalcResidual(TPZElementMatrix &ef);
	
	/**
	 * @brief Computes solution and its derivatives in the local coordinate qsi.
	 * @param qsi [in] master element coordinate
	 * @param normal normal vector
	 * @param leftsol [out] left finite element solution
	 * @param rightsol [out] right finite element solution
	 * @param dleftsol [out] left solution derivatives
	 * @param drightsol [out] right solution derivatives
	 * @param leftaxes [out] axes associated with the derivative of the left element
	 * @param rightaxes [out] axes associated with the derivative of the right element
	 */
	virtual void ComputeSolution(TPZVec<REAL> &qsi,
								 TPZVec<REAL> &normal,
								 TPZSolVec &leftsol, TPZGradSolVec &dleftsol,TPZFMatrix<REAL> &leftaxes,
								 TPZSolVec &rightsol, TPZGradSolVec &drightsol,TPZFMatrix<REAL> &rightaxes);
	
	/**
	 * @brief Computes solution and its derivatives in local coordinate qsi
	 * @param qsi master element coordinate
	 * @param phi matrix containing shape functions compute in qsi point
	 * @param dphix matrix containing the derivatives of shape functions in the direction of the axes
	 * @param axes axes indicating the direction of the derivatives
	 * @param sol finite element solution
	 * @param dsol solution derivatives
	 */
	virtual void ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphix,
								 const TPZFMatrix<REAL> &axes, TPZSolVec &sol, TPZGradSolVec &dsol);
	
	/**
	 * @brief Computes solution and its derivatives in the local coordinate qsi.
	 * @param qsi master element coordinate
	 * @param sol finite element solution
	 * @param dsol solution derivatives
	 * @param axes axes associated with the derivative of the solution
	 */
	virtual void ComputeSolution(TPZVec<REAL> &qsi,
								 TPZSolVec &sol, TPZGradSolVec &dsol,TPZFMatrix<REAL> &axes);
	
	void VetorialProd(TPZVec<REAL> &ivet,TPZVec<REAL> &jvet,TPZVec<REAL> &kvet);
	
	/** @brief Prints attributes of the object */
	void Print(std::ostream &out = std::cout) const;
	
	/**
	 * @see Base class for comments
	 * @brief Interface elements does not have graphical representation
	 */
	virtual void CreateGraphicalElement(TPZGraphMesh & graphmesh, int dimension){
		//Nothing to be done here
	}
	
	/**
	 * @brief Verifies the existence of interfaces associates
	 * with the side of an element
	 */
	/**
	 * Case to interface should exist and exists only a returns 1, \n
	 * case to interface should not exist and does not exist returns 1, \n
	 * otherwise returns 0
	 */
	static int ExistInterfaces(TPZCompElSide &comp);
	
	static int FreeInterface(TPZCompMesh &cmesh);
	
	/** @brief Make a clone of the fine mesh into clustered mesh.*/
	void CloneInterface(TPZCompMesh *aggmesh);
	
	static int main(TPZCompMesh &cmesh);
	
	void EvaluateError(void (*fp)(TPZVec<REAL> &loc,TPZVec<REAL> &val,TPZFMatrix<REAL> &deriv),
					   TPZVec<REAL> &errors, TPZBlock<REAL> * /*flux */);
	
	/** @brief ComputeError computes the element error estimator */
	virtual void ComputeErrorFace(int errorid,
								  TPZVec<REAL> &errorL,
								  TPZVec<REAL> &errorR);
	
	/** @brief Integrate a variable over the element. */
	virtual void Integrate(int variable, TPZVec<REAL> & value);
	
	void IntegrateInterface(int variable, TPZVec<REAL> & value);
	
	/** 
	 * \f$ opt = 0 \f$ ->  Evaluates \f$ \sqrt{ \int { (leftsol - rightsol)^2 } } \f$ \n
	 * \f$ opt = 1 \f$ ->  Evaluates \f$ Max { | leftsol - rightsol | } \f$
	 */
	void EvaluateInterfaceJump(TPZSolVec &jump, int opt);
	
	/** @brief Returns the unique identifier for reading/writing objects to streams */
	virtual int ClassId() const;
	/** @brief Saves the element data to a stream */
	virtual void Write(TPZStream &buf, int withclassid);
	
	/** @brief Reads the element data from a stream */
	virtual void Read(TPZStream &buf, void *context);
	
};

#endif

