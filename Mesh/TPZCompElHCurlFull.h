/**
 * @file
 * @brief Contains declaration of TPZCompElHCurl class which implements a full-order computational HCurl-conforming element
 */

#ifndef TPZCOMPELHCURLFULL_H
#define TPZCOMPELHCURLFULL_H

#include <TPZCompElHCurl.h>

/**
 * @brief This class implements a  full-order computational HCurl-conforming element. \ref CompElement "Computational Element"
 * @addtogroup CompElement
 * @{
 */
/** 
 * By varying the classes passed as template arguments, the complete family of computational elements are implemented
 */
template<class TSHAPE>
class TPZCompElHCurlFull : public TPZCompElHCurl<TSHAPE> {
public:
    TPZCompElHCurlFull(TPZCompMesh &mesh, TPZGeoEl *gel, int64_t &index);
	
	TPZCompElHCurlFull(TPZCompMesh &mesh, const TPZCompElHCurlFull<TSHAPE> &copy);
	
	/**
	 * @brief Constructor used to generate patch mesh... generates a map of connect index from
	 * global mesh to clone mesh
	 */
	TPZCompElHCurlFull(TPZCompMesh &mesh,
				  const TPZCompElHCurlFull<TSHAPE> &copy,
				  std::map<int64_t,int64_t> & gl2lcConMap,
				  std::map<int64_t,int64_t> & gl2lcElMap);
	
	TPZCompElHCurlFull();
	
	virtual ~TPZCompElHCurlFull() = default;

    /** @brief Returns the unique identifier for reading/writing objects to streams */
    int ClassId() const override;

    /**
    * @brief Number of shapefunctions of the connect associated
    * @param connect connect number
    * @return number of shape functions
    */
    int NConnectShapeF(int connect, int order) const override;

    void SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) override;
protected:
    /**
    * @brief Returns a matrix index of the shape and vector  associate to element
    * @param[out] IndexVecShape Indicates the pair vector/shape function that will construct the approximation space
    * @param[in] connectOrder Order of the connects
    */
    void IndexShapeToVec(TPZVec<std::pair<int,int64_t> > & indexVecShape, const TPZVec<int>& connectOrder) const override;

    template<class TSIDESHAPE=TSHAPE>
    static void StaticIndexShapeToVec(TPZVec<std::pair<int,int64_t>> & indexVecShape, const TPZVec<int>& connectOrder,
                                      const TPZVec<int64_t>& firstH1ShapeFunc, const TPZVec<int> &sidesH1Ord, TPZVec<unsigned int>& shapeCountVec,
                                      const TPZVec<int>& transformationIds);

    /**
     * @brief This method calculates the appropriate side orders for the correct calculation of the SCALAR shape functions.
     * @param[out] ord Vector that will be filled with the corresponding order of each connect to compute the desired h1 functions
     */
    void CalculateSideShapeOrders(TPZVec<int> &ord) const override;
};


/** @} */

#endif
