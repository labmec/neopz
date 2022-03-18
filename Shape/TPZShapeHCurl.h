#ifndef TPZSHAPEHCURL_H
#define TPZSHAPEHCURL_H

#include "pzreal.h"
#include "pzvec.h"
#include "pztrnsform.h"
template<class T>
class TPZFMatrix;

#include "TPZShapeData.h"

/// Generates HCurl spaces in the reference element using constant vector fields and H1 shape functions.
template <class TSHAPE>
struct TPZShapeHCurl
{
    
    TPZShapeHCurl() = default;
    //! Should be called once per element. Initializes the data structure
    static void Initialize(const TPZVec<int64_t> &ids,
                    TPZVec<int> &connectorders,                    
                    TPZShapeData &data);
    //! Computes the pair (vec, h1 shape) that are used for creating the HCurl shape functions
    static void ComputeVecandShape(TPZShapeData &data);
    /** 
        @brief  Computes the number of shape functions for a given connect
        @param [in] connect which connect
        @param [in] order connect ordder
    */
    static int ComputeNConnectShapeF(const int connect, const int order);
    //! Number of HCurl shape functions for a given (filled) TPZShapeData
    static int NHCurlShapeF(const TPZShapeData &data);
    /**
     * @brief Calculates phi and curl phi at a given integration point.
     * @param[in] pt integration point in the reference element
     * @param[in] shapedata 
     */
    static void Shape(TPZVec<REAL> &pt, TPZShapeData &data,
                      TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &curlphi);

    /**
    * @brief Returns a matrix index of the shape and vector  associate to element
    * @param[in] connectOrder Order of the HCurl connects
    * @param[in] sirstH1ShapeFunc Index of the first H1 shape function associated with a given side
    * @param[in] sidesH1Ord H1 order of each side
    * @param[in] nodeIds Ids of the nodes. Used for adjusting order of HCurl functions.
    * @param[out] shapeCountVec Stores the number of shape function for each connect
    * @param[out] indexVecShape Indicates the pair vector/H1 shape function that will 
    be used for constructing the HCurl functions
    */
    static void StaticIndexShapeToVec(const TPZVec<int>& connectOrder,
                                      const TPZVec<int64_t>& firstH1ShapeFunc,
                                      const TPZVec<int> &sidesH1Ord,
                                      const TPZVec<int64_t>& nodeIds,
                                      TPZVec<unsigned int>& shapeCountVec,
                                      TPZVec<std::pair<int,int64_t>> & indexVecShape);

    /**
       @brief Internal method for calculating the needed order of H1 functions
       for each side given the HCurl connect orders listed in `ordHCurl` param
       @param[in] ordHCurl The effective order of each HCurl connect associated with
       a the side of shape TSIDESHAPE
       @param[out] the order of the h1 connects needed for computing the hcurl functions
       
       @note since the H1 vertex functions are always needed and their order
       does not change, `ord` has size `NSides-NCornerNodes`.
     */
    static void CalcH1ShapeOrders(const TPZVec<int> &ordHCurl,
                                  TPZVec<int> &ord);

    //! Maximum (actual) polynomial order. Useful for choosing integration rule
    [[nodiscard]] static int MaxOrder(const int ordh1);
};

#endif
