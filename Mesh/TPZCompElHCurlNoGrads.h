#ifndef _TPZCOMPELHCURLNOGRADS_H_
#define _TPZCOMPELHCURLNOGRADS_H_

#include <TPZCompElHCurl.h>


/**
   @brief TPZCompElHCurlNoGrads<TSHAPE> generates a subset of HCurl functions
   by removing functions needed to represent the high-order gradient fields.
   
   Upon treatment of the mesh (i.e., filtering some of the edges out of the system),
   these elements can be used for generating an approximation space such that
   dim(curl(phi)) = dim(phi).
*/
template<class TSHAPE>
class TPZCompElHCurlNoGrads  : public TPZCompElHCurl<TSHAPE> {
public:
  //!Default constructor.
  TPZCompElHCurlNoGrads();
  //! Ctor taking mesh, geoel and returning index
  TPZCompElHCurlNoGrads(TPZCompMesh &mesh, TPZGeoEl *gel);
  /**
   * @brief Number of shapefunctions of the connect associated
   * @param connect connect number
   * @return number of shape functions
   */
  int NConnectShapeF(int connect, int order) const override;
  //! Fills data.fMasterDirections, data.fVecShapeIndex and data.fShapeType.
  void InitMaterialData(TPZMaterialData &data) override;
  //! Computes data.phi and data.curlphi to be used in the integration points.
  void ComputeRequiredData(TPZMaterialDataT<STATE> &data, TPZVec<REAL> &qsi) override{
    ComputeRequiredDataT(data,qsi);
  }
  //! Computes data.phi and data.curlphi to be used in the integration points.
  void ComputeRequiredData(TPZMaterialDataT<CSTATE> &data, TPZVec<REAL> &qsi) override{
    ComputeRequiredDataT(data,qsi);
  }
  //! Fills data.sol and data.curlsol.
  void ReallyComputeSolution(TPZMaterialDataT<STATE>& data) override{
    ReallyComputeSolutionT(data);
  }
  //! Fills data.sol and data.curlsol.
  void ReallyComputeSolution(TPZMaterialDataT<CSTATE>& data) override{
    ReallyComputeSolutionT(data);
  }
protected:
  //! Adjusts the number of shape functions and block size of the connects.
  void AdjustConnects();
  //! Combines/filters the HCurl functions from TPZCompElHCurlFull<TSHAPE>.
  void ComputeShape(TPZMaterialData &data, TPZFMatrix<REAL> &phiHCurl);
  template<int DIM>
  /** 
   * @brief Combines/filters the HCurl functions from TPZCompElHCurlFull<TSHAPE>
   and compute their curl, filling data.curlphi.*/
  void ComputeCurl(TPZMaterialData &data);
  //! Computes data.phi and data.curlphi to be used in the integration points.
  template<class TVar>
  void ComputeRequiredDataT(TPZMaterialDataT<TVar> &data, TPZVec<REAL> &qsi);
  //! Fills data.sol and data.curlsol.
  template<class TVar>
  void ReallyComputeSolutionT(TPZMaterialDataT<TVar> &data);
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
#endif /* _TPZCOMPELHCURLNOGRADS_H_ */
