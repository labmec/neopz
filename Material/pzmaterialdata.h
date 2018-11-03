/**
 * @file
 * @brief Contains the TPZMaterialData class which implements an interface between TPZCompEl::CalcStiff and TPZMaterial::Contribute methods.
 */

#ifndef PZMATERIALDATA_H
#define PZMATERIALDATA_H

#include "pzmanvector.h"
#include "pzfmatrix.h"


/**
 * @ingroup material
 * @brief This class implements an interface between TPZCompEl::CalcStiff and TPZMaterial::Contribute methods. \n
 * It request to the material which attributes must be computed by the computational element and trigger their computation.\n
 * Attributes are solution and its derivatives, X coordinate, etc.
 * @since April 10, 2007
 */

/// Represent the state variables of a finite element approximation
typedef TPZManVector<STATE, 10> TPZFemSol;
/// Represents the gradient of a state variable of a finite element approximation
typedef TPZFNMatrix<15, STATE> TPZFemGradSol;
typedef TPZManVector<TPZFemSol,20> TPZSolVec;
typedef TPZManVector<TPZFemGradSol,20> TPZGradSolVec;


class TPZMaterialData : public TPZSavable {
    
public:
    
    
    enum MShapeFunctionType {EEmpty, EScalarShape, EVecandShape, EVecShape};
    // EScalarShape : regular shape functions: one shape function used for all state variables (default)
    // EVecandShape : HDiv type shape function and a scalar function
    // EVecShape : a vector valued shape function
    
    MShapeFunctionType fShapeType;
    /** @name Flags indicating whether some attributes shall be computed or not */
    /** @{ */
    bool fNeedsSol, fNeedsNeighborSol, fNeedsHSize, fNeedsNeighborCenter, fNeedsNormal;
    /** @} */
    
    /** @name Attributes to be computed in CalcStiff */
    /** @{ */
    
    /// vector of shapefunctions (format is dependent on the value of shapetype)
    TPZFNMatrix<220, REAL> phi;
    /// values of the derivative of the shape functions over the master element
    TPZFNMatrix<660, REAL> dphi;
    /// values of the derivative of the shape functions
    TPZFNMatrix<660, REAL> dphix;
    /// axes indicating the directions of the derivatives of the shapefunctions
    TPZFNMatrix<9,REAL> axes;
    /// value of the jacobian at the integration point
    TPZFNMatrix<9,REAL> jacobian;
    /// value of the inverse of the jacobian at the integration point
    TPZFNMatrix<9,REAL> jacinv;
    /// normal to the element at the integration point
    TPZManVector<REAL,3> normal;
    /// value of the coordinate at the integration point
    TPZManVector<REAL,3> x;
    /// value of the coordinate at the integration point corresponding to the x-parametric coordinate (master element)
    TPZManVector<REAL,3> xParametric;
    /// maximum polinomial order of the shape functions
    int p;
    /// vector of the solutions at the integration point
    TPZSolVec sol;
    /// vector of the derivatives of the solution at the integration point
    TPZGradSolVec dsol;
    /// measure of the size of the element
    REAL HSize;
    /// determinant of the jacobian
    REAL detjac;
    /// value of the coordinate at the center of the element
    TPZManVector<REAL,3> XCenter;
    
    /// number of dual function (e.g. pressure in HDiv approximations)
    int numberdualfunctions;
    
    //Id of associated geo element
    int gelElId;
    
    /// correspondence between normal vector index and index of the shape functions
    TPZManVector<std::pair<int,int64_t> > fVecShapeIndex;
    /// list of normal vectors
    TPZFNMatrix<180> fNormalVec;
    /** @} */
    
    
    
    /** @brief Index of the current integration point being evaluated **/
    /** Needed for materials with memory **/
    int intLocPtIndex;
    
    /** @brief global point index */
    int intGlobPtIndex;
    
    /** @brief amount of points in the integrstion rule */
    int NintPts;
    
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
    virtual void Write(TPZStream &buf, int withclassid) const;
    
    /** @brief Reads the element data from a stream */
    virtual void Read(TPZStream &buf, void *context);
    
    /** @brief Compares the object for identity with the object pointed to, eventually copy the object */
    /**
     * Compares both objects bitwise for identity. Put an entry in the log file if different
     * overwrite the calling object if the override flag is true
     */
    virtual bool Compare(TPZSavable *copy, bool override = false);
    
    /** @brief Compares the object for identity with the object pointed to, eventually copy the object */
    /**
     * Compares both objects bitwise for identity. Put an entry in the log file if different
     * overwrite the calling object if the override flag is true
     */
    virtual bool Compare(TPZSavable *copy, bool override = false) const;
    
    /** @brief Computes the flux values based on a Material of Hdiv approx space */
    void ComputeFluxValues(TPZFMatrix<REAL> & fluxes);
    
    public:
virtual int ClassId() const;

};

#endif
