//$Id: pzmaterialdata.h,v 1.13 2010-08-12 13:50:04 phil Exp $

#ifndef PZMATERIALDATA_H
#define PZMATERIALDATA_H

#include "pzmanvector.h"
#include "pzfmatrix.h"
#include "pzmaterialid.h"


/**
This class implements an interface between TPZCompEl::CalcStiff and TPZMaterial::Contribute methods. It request to the material which attributes must be computed by the computational element and trigger their computation. Attributes are solution and its derivatives, X coordinate, etc.

@since April 10, 2007
*/


class TPZMaterialData : public TPZSaveable {

public:

/** Flags indicating whether some attributes shall be computed or not */
  bool fNeedsSol, fNeedsNeighborSol, fNeedsHSize, fNeedsNeighborCenter, fNeedsNormal;

/** Attributes to be computed in CalcStiff */
  TPZFNMatrix<220> phi, phil, phir;
  TPZFNMatrix<660> dphix, dphixl, dphixr;
  TPZFNMatrix<9> axes, axesleft, axesright;
  TPZFNMatrix<9> jacobian, leftjac, rightjac;
  TPZFNMatrix<9> jacinv, leftjacinv, rightjacinv;
  TPZManVector<REAL,3> normal;
  TPZManVector<REAL,3> x;
  int p, leftp, rightp;
  TPZManVector<REAL,10> sol, soll, solr;
  TPZFNMatrix<30> dsol, dsoll, dsolr;
  REAL HSize;
  REAL detjac, leftdetjac, rightdetjac;
  TPZManVector<REAL,3> XLeftElCenter, XRightElCenter;

/** Index of the current integration point being evaluated **/
/** Needed for materials with memory **/

  int intPtIndex;

/** Class constructor */
  TPZMaterialData();

/** Copy constructor */
  TPZMaterialData( const TPZMaterialData &cp );

/** Class destructor */
  ~TPZMaterialData();

/** Set all flags at once */
  void SetAllRequirements(bool set);

  void InvertLeftRightData();

  TPZMaterialData &operator= (const TPZMaterialData &cp );

	/**
	 * Print the data
	 */
	void Print(std::ostream &out) const;
	/**
	 * Print the data in a format suitable for Mathematica
	 */
	void PrintMathematica(std::ostream &out) const;
	/**
	 * Save the element data to a stream
	 */
	virtual void Write(TPZStream &buf, int withclassid);
	
	/**
	 * Read the element data from a stream
	 */
	virtual void Read(TPZStream &buf, void *context);
	
	/// Compare the object for identity with the object pointed to, eventually copy the object
	/**
	 * compare both objects bitwise for identity. Put an entry in the log file if different
	 * overwrite the calling object if the override flag is true
	 */
	virtual bool Compare(TPZSaveable *copy, bool override = false);
	
	/// Compare the object for identity with the object pointed to, eventually copy the object
	/**
	 * compare both objects bitwise for identity. Put an entry in the log file if different
	 * overwrite the calling object if the override flag is true
	 */
	virtual bool Compare(TPZSaveable *copy, bool override = false) const;
	

	virtual int ClassId() const
	{
		return TPZMATERIALDATAID;
	}
};


#endif
