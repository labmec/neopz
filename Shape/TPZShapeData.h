/**
 * @file
 * @brief Contains the TPZShapeData class which stores the shape function information in the master element coordinate system.
 */

#ifndef PZSHAPEDATA_H
#define PZSHAPEDATA_H

#include "pzmanvector.h"
#include "pzfmatrix.h"
#include "pztrnsform.h"

/**
 * @ingroup shape
 * @brief This class implements a structure to store data regarding shape functions in the master element coordinate system.
 *
 */
class TPZShapeData : virtual public TPZSavable {

public:
    //! Default constructor.
    TPZShapeData();
    //! Copy constructor.
    TPZShapeData(const TPZShapeData&) = default;
    //! Move constructor.
    TPZShapeData(TPZShapeData&&) = default;
    //! Destructor.
    virtual ~TPZShapeData() = default;
    //! Copy assignment operator
    TPZShapeData& operator=(const TPZShapeData&) = default;
    //! Move assignment operator
    TPZShapeData& operator=(TPZShapeData&&) = default;

    //@{
    //! Read and Write methods
    int ClassId() const override;
    
    void Write(TPZStream &buf, int withclassid) const override;
    
    void Read(TPZStream &buf, void *context) override;
    //@}
    //! Shape function type as a string.
    std::string ShapeFunctionType() const;
    /** @brief Prints the data */
    void Print(std::ostream &out) const;
    /** @brief Prints the data in a format suitable for Mathematica */
    void PrintMathematica(std::ostream &out) const;    

    static constexpr int MatDataNumPhi{60};
    static constexpr int MatDataNumDir{81};
    static constexpr int MatDataDimSol{10};//TODO:Remove?
    static constexpr int MatDataNumSol{20};//TODO:Remove?
    /*! Type of the shape function associated with an element*/
    enum MShapeFunctionType {EEmpty,
        EScalarShape,///< Scalar shape functions (H1, L2)
        EVecandShape,///< Composite shape function (scalar function and vector field, HDiv and HCurl spaces)
        EVecShape///< Vector shape function (e.g. vector H1 space)
    };//TODO:Remove?
    //! Type of shape function
    MShapeFunctionType fShapeType{EEmpty};
    //! Corner node ids determine the parameter transformations to the sides
    TPZManVector<int64_t,8> fCornerNodeIds;
    //! Connect orders determine the number of shape functions
    TPZManVector<int,27> fH1ConnectOrders;
    //! Number of shape functions by connect
    TPZManVector<int,27> fH1NumConnectShape;
    //! Parametric transforms from the interior to the side
    TPZManVector<TPZTransform<REAL>, 20> fSideTransforms;
    //! Transformation ids for each side
    TPZManVector<int, 20> fSideTransformationId;
    //! Vector of shapefunctions (format is dependent on the value of shapetype) over the master element
    TPZFNMatrix<MatDataNumPhi, REAL> fPhi;
    //! Values of the derivative of the shape functions over the master element
    TPZFNMatrix<3*MatDataNumPhi, REAL> fDPhi;

    // NEEDED BY HDIV OR HCURL
    //! Connect orders determine the number of shape functions
    TPZManVector<int,27> fHDivConnectOrders;
    //! Number of shape functions by connect
    TPZManVector<int,27> fHDivNumConnectShape;
    //! The orientation of the sides (either -1 or +1)
    TPZManVector<int,20> fSideOrient;
    //! Directions on the master element
    TPZFNMatrix<3*MatDataNumDir> fMasterDirections;
    //! Correspondence between direction vector index and index of the shape functions. Used for H(div) and H(curl) approximation spaces.
    TPZManVector<std::pair<int,int64_t>, MatDataNumPhi > fSDVecShapeIndex;
        
};
#endif

