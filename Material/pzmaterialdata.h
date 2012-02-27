/**
 * \file
 * @brief Contains the TPZMaterialData class which implements an interface between TPZCompEl::CalcStiff and TPZMaterial::Contribute methods.
 */
//$Id: pzmaterialdata.h,v 1.15 2011-05-30 20:19:12 denise Exp $

#ifndef PZMATERIALDATA_H
#define PZMATERIALDATA_H

#include "pzmanvector.h"
#include "pzfmatrix.h"
#include "pzmaterialid.h"


/**
 * @ingroup material
 * @brief This class implements an interface between TPZCompEl::CalcStiff and TPZMaterial::Contribute methods. \n
 * It request to the material which attributes must be computed by the computational element and trigger their computation.\n
 * Attributes are solution and its derivatives, X coordinate, etc.
 * @since April 10, 2007
 */

/// Represent the state variables of a finite element approximation
typedef TPZManVector<REAL, 10> TPZFemSol;
/// Represents the gradient of a state variable of a finite element approximation
typedef TPZFNMatrix<30> TPZFemGradSol;
typedef TPZVec<TPZFemSol > TPZSolVec;
typedef TPZVec<TPZFemGradSol > TPZGradSolVec;


class TPZMaterialData : public TPZSaveable {
	
public:
	
	/** @name Flags indicating whether some attributes shall be computed or not */
	/** @{ */
	bool fNeedsSol, fNeedsNeighborSol, fNeedsHSize, fNeedsNeighborCenter, fNeedsNormal;
	/** @} */
	
	/** @name Attributes to be computed in CalcStiff */
	/** @{ */
	TPZFNMatrix<220> phi;//, phil, phir;
	TPZFNMatrix<660> dphix;//, dphixl, dphixr;
	TPZFNMatrix<9> axes;//, axesleft, axesright;
	TPZFNMatrix<9> jacobian;//, leftjac, rightjac;
	TPZFNMatrix<9> jacinv;//, leftjacinv, rightjacinv;
	TPZManVector<REAL,3> normal;
	TPZManVector<REAL,3> x;
	int p;//, leftp, rightp;
	TPZManVector<TPZFemSol, 10> sol;//, soll, solr;
	TPZManVector<TPZFemGradSol, 10> dsol;//, dsoll, dsolr;
	REAL HSize;
	REAL detjac;//, leftdetjac, rightdetjac;
    TPZManVector<REAL,3> XCenter;
	//TPZManVector<REAL,3> XLeftElCenter, XRightElCenter;
	
	int numberdualfunctions;
	TPZManVector<std::pair<int,int> > fVecShapeIndex;
	TPZFNMatrix<100> fNormalVec;
	/** @} */


	/** @brief Index of the current integration point being evaluated **/
	/** Needed for materials with memory **/
	int intPtIndex;
	
	/** @brief Default constructor */
	TPZMaterialData();
	
	/** @brief Copy constructor */
	TPZMaterialData( const TPZMaterialData &cp );
	
	/** @brief Default destructor */
	~TPZMaterialData();
	
	/** @brief Set all flags at once */
	void SetAllRequirements(bool set);
	
	//void InvertLeftRightData();
	
	TPZMaterialData &operator= (const TPZMaterialData &cp );
	
	/** @brief Prints the data */
	void Print(std::ostream &out) const;
	/** @brief Prints the data in a format suitable for Mathematica */
	void PrintMathematica(std::ostream &out) const;
	/** @brief Saves the element data to a stream */
	virtual void Write(TPZStream &buf, int withclassid);
	
	/** @brief Reads the element data from a stream */
	virtual void Read(TPZStream &buf, void *context);
	
	/** @brief Compares the object for identity with the object pointed to, eventually copy the object */
	/**
	 * Compares both objects bitwise for identity. Put an entry in the log file if different
	 * overwrite the calling object if the override flag is true
	 */
	virtual bool Compare(TPZSaveable *copy, bool override = false);
	
	/** @brief Compares the object for identity with the object pointed to, eventually copy the object */
	/**
	 * Compares both objects bitwise for identity. Put an entry in the log file if different
	 * overwrite the calling object if the override flag is true
	 */
	virtual bool Compare(TPZSaveable *copy, bool override = false) const;

	virtual int ClassId() const
	{
		return TPZMATERIALDATAID;
	}
};

#endif
