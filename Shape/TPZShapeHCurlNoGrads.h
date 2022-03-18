#ifndef TPZSHAPEHCURLNOGRADS_H
#define TPZSHAPEHCURLNOGRADS_H

#include "pzreal.h"
#include "pzvec.h"
#include "pztrnsform.h"
template<class T>
class TPZFMatrix;

#include "TPZShapeData.h"

/// Generates HCurl spaces in the reference element using constant vector fields and H1 shape functions.
template <class TSHAPE>
struct TPZShapeHCurlNoGrads
{
    
    TPZShapeHCurlNoGrads() = default;
    //! Should be called once per element. Initializes the data structure
    static void Initialize(const TPZVec<int64_t> &ids,
                    TPZVec<int> &connectorders,                    
                    TPZShapeData &data);
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


    //! Maximum (actual) polynomial order. Useful for choosing integration rule
    [[nodiscard]] static int MaxOrder(const int ordh1);

  /**
     @brief Given the original indices of functions of a a HCurl element,
     calculates the subset corresponding to
     the filtered higher-order (face, interior) functions.
     @param[in] firstHCurlFunc index of first Hurl function corresponding to a given connect.
     @param[in] conOrders connect orders.
     @param[in] filteredFuncs indices of desired of Hcurl functions.
  */
  static void HighOrderFunctionsFilter(
    const TPZVec<int> &firstHCurlFunc,
    const TPZVec<int> &conOrders,
    TPZVec<int> &filteredHCurlFuncs);
};

#endif
