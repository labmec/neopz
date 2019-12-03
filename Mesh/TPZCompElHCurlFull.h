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
};


/** @} */

#endif
