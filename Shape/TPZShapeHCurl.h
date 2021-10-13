#ifndef TPZSHAPEHCURL_H
#define TPZSHAPEHCURL_H

#include "pzreal.h"
#include "pzvec.h"
#include "pztrnsform.h"
template<class T>
class TPZFMatrix;

#include "TPZShapeData.h"

/// Traditional HDiv spaces, data structures that do not depend on the geometric map
template <class TSHAPE>
struct TPZShapeHCurl
{
    
    TPZShapeHCurl();
    
    static void Initialize(TPZVec<int64_t> &ids,
                    TPZVec<int> &connectorders,                    
                    TPZShapeData &data);
    
    static void ComputeMasterDirections(TPZShapeData &data);
    
    static void ComputeVecandShape(TPZShapeData &data);
    
    static int NConnectShapeF(int connect, TPZShapeData &data);

    static int NHCurlShapeF(TPZShapeData &data);
    
    static int NH1ShapeF(TPZShapeData &data)
    {
        return data.fPhi.Rows();
    }
    
    static void Shape(TPZVec<REAL> &pt, TPZShapeData &data, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &divphi);

    /**
    * @brief Returns a matrix index of the shape and vector  associate to element
    * @param[out] IndexVecShape Indicates the pair vector/shape function that will construct the approximation space
    * @param[in] connectOrder Order of the connects
    */

    template<class TSIDESHAPE=TSHAPE>
    static void StaticIndexShapeToVec(                                      TPZVec<std::pair<int,int64_t>> & indexVecShape, const TPZVec<int>& connectOrder,
                                      const TPZVec<int64_t>& firstH1ShapeFunc, const TPZVec<int> &sidesH1Ord, TPZVec<unsigned int>& shapeCountVec,
                                      const TPZVec<int64_t>& nodeIds);

    /**
       @brief Internal method for calculating the needed h1 order for each side given
       the HCurl connect orders listed in `ordHCurl` param
       @param[in] ordHCurl The effective order of each HCurl connect associated with
       a the side of shape TSIDESHAPE
       @param[out] the order of the h1 connects needed for computing the hcurl functions
       
       @note since the h1 vertex functions are always needed, ord has size
       `NSides-NCornerNodes`.
     */
    template<class TSIDESHAPE=TSHAPE>
    static void StaticCalcH1ShapeOrders(const TPZVec<int> &ordHCurl,
                                             TPZVec<int> &ord);

};

#endif
