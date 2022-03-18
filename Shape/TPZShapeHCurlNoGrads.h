#ifndef TPZSHAPEHCURLNOGRADS_H
#define TPZSHAPEHCURLNOGRADS_H

#include "pzreal.h"
#include "pzvec.h"
#include "pztrnsform.h"
template<class T>
class TPZFMatrix;

#include "TPZShapeData.h"

/**
   @brief This class generates a set of Hcurl conforming shape functions
   such that the high order gradient fields are removed.
   In order to generate a truly grad-free Hcurl-conforming approximation
   space, the TPZHCurlEquationFilter class should be used.
*/
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
     @brief Filters the high-order gradient fields out of a HCurl element.
     This filter acts only on face and volumetric functions.
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
