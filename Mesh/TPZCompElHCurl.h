/**
 * @file
 * @brief Contains declaration of TPZCompElHCurl class which implements a generic computational HCurl-conforming element
 */

#ifndef TPZCOMPELHCURL_H
#define TPZCOMPELHCURL_H

#include <pzelctemp.h>
#include "TPZMaterialDataT.h"

/**
 * @brief This class implements a "generic" computational HCurl-conforming element. \ref CompElement "Computational Element"
 * @addtogroup CompElement
 * @{
 */
/** 
 * By varying the classes passed as template arguments, the complete family of computational elements are implemented
 */
template<class TSHAPE>
class TPZCompElHCurl : public TPZIntelGen<TSHAPE> {
protected:
    ///! Indexes of the connects associated with the elements
    TPZManVector<int64_t,TSHAPE::NSides - TSHAPE::NCornerNodes> fConnectIndexes =
        TPZManVector<int64_t,TSHAPE::NSides - TSHAPE::NCornerNodes>(TSHAPE::NSides-TSHAPE::NCornerNodes,-1);
    
    HCurlFamily fhcurlfam = DefaultFamily::fHCurlDefaultValue;
public:
    TPZCompElHCurl(TPZCompMesh &mesh, TPZGeoEl *gel, const HCurlFamily hcurlfam = DefaultFamily::fHCurlDefaultValue);
	
    TPZCompElHCurl(TPZCompMesh &mesh, const TPZCompElHCurl<TSHAPE> &copy);
	
    /**
     * @brief Constructor used to generate patch mesh... generates a map of connect index from
     * global mesh to clone mesh
     */
    TPZCompElHCurl(TPZCompMesh &mesh,
                   const TPZCompElHCurl<TSHAPE> &copy,
                   std::map<int64_t,int64_t> & gl2lcConMap,
                   std::map<int64_t,int64_t> & gl2lcElMap);
	
    TPZCompElHCurl();

    TPZCompEl *Clone(TPZCompMesh &mesh) const  override {
        return new TPZCompElHCurl<TSHAPE> (mesh, *this);
    }

    /**
     * @brief Create a copy of the given element. The clone copy have the connect indexes
     * mapped to the local clone connects by the given map
     * @param mesh Patch clone mesh
     * @param gl2lcConMap map the connects indexes from global element (original) to the local copy.
     * @param gl2lcElMap map the indexes of the elements between the original element and the patch element
     */
    virtual TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,std::map<int64_t,int64_t> & gl2lcConMap,std::map<int64_t,int64_t>&gl2lcElMap) const override
    {
        return new TPZCompElHCurl<TSHAPE> (mesh, *this, gl2lcConMap, gl2lcElMap);
    }
    ~TPZCompElHCurl();

    /** @brief Returns the unique identifier for reading/writing objects to streams */
    int ClassId() const override;
    /** @brief Set create function in TPZCompMesh to create elements of this type */
    void SetCreateFunctions(TPZCompMesh *mesh) override;

    MElementType Type() override;

    int NConnects() const override;
	/** @brief return the local index for connect */
    int SideConnectLocId(int con, int side) const override;

    int NSideConnects(int side) const override;

    int NCornerConnects() const override {
        return 0;
    }

    void GetInterpolationOrder(TPZVec<int> &ord) override;
        
    int64_t ConnectIndex(int con) const override;

    void SetConnectIndex(int i, int64_t connectindex) override;

    //! Reference to the connect vector
    const TPZVec<int64_t> & ConnectVec() const override{
        return fConnectIndexes;
    }

    /**
    * @brief Number of shapefunctions of the connect associated
    * @param connect connect number
    * @return number of shape functions
    */
    int NConnectShapeF(int connect, int order) const override;

    /**
    * @brief return the interpolation order of the polynomial for connect
    **/
    int ConnectOrder(int connect) const override;

    /** @brief Returns the actual interpolation order of the polynomial along the side*/
    int EffectiveSideOrder(int side) const override;

    /** @brief Sets the interpolation order of side to order*/
    void SetSideOrder(int side, int order) override;

    void RestrainSide(int side, TPZInterpolatedElement *large, int neighbourside) override;
    /** @brief Initialize a material data and its attributes based on element dimension, number
    * of state variables and material definitions */
    void InitMaterialData(TPZMaterialData &data) override;
    /**
     * @brief Computes the shape functions at the point qsi corresponding to x.
     * The values correspond to the function in the deformed element.    *
     * @param qsi point in master element coordinates
     * @param data fields regarding the geometric transformation and data.phi and data.curlphi will be filled
     */
    void ComputeShape(TPZVec<REAL> &qsi, TPZMaterialData &data) override;
    
    /**
     * @brief Computes the shape functions at the point qsi corresponding to x.
     * @param qsi point in master element coordinates
     * @param phi vector of values of shapefunctions, dimension (numshape,1)
     * @param dphi matrix of derivatives of shapefunctions in master element coordinates, dimension (dim,numshape)
     */
    void Shape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi) override;

    /**
     * @brief Computes the trace of the shape functions associated with a given side.
     * @param side side in which the traces will be calculated
     * @param point point in side's parametric coordinates in which the shape functions will be evaluated
     * @param phi trace of the shape functions
     * @param dphi trace of the curl of the shape functions
     */
    void SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) override;
    //@}
protected:

    /**
       @brief Applies the covariant Piola transform for obtaining the basis functions in the deformed element.
       @param [in] phiref HCurl shape functions in the reference element
       @param [in] jacinv inverse of the jacobian of the element mapping
       @param [in] axes relates the directions \$f (\xi,\eta,\zeta) \rightarrow (x,y,z) \$f
       @param [out] phi HCurl basis functions in the deformed element
    */
    static void TransformShape(const TPZFMatrix<REAL> &phiref, const REAL detjac,
                               const TPZFMatrix<REAL> &jacinv,
                               const TPZFMatrix<REAL> &axes,
                               TPZFMatrix<REAL> &phi);

    /**
      @brief Applies the appropriate transform for obtaining the curl of the basis functions in the deformed element
      @tparam TDIM Dimension of the element curl will have dim=2*TDIM-3 for TDIM=2,3 (and 1 for TDIM=1).
      @param [in] curlphiref curl in the reference element
      @param [in] detJac  determinant of the jacobian 
      @param [in] jacobian jacobian of the element mapping
      @param [out] curlphi curl of the vector shape functions in the deformed element
     */
    template<int TDIM=TSHAPE::Dimension>
    static void TransformCurl(const TPZFMatrix<REAL> &curlphiref, const REAL detjac,
                              const TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &curlphi);
    

    /**
     * @brief Creates the connects for a generic HCurl-conforming element.
     * This method must be called in the constructor of every child class of TPZCompElHCurl.
     * @param mesh Computational mesh in which the connects will be inserted
     */
    void CreateHCurlConnects(TPZCompMesh &mesh);

    void ReallyComputeSolution(TPZMaterialDataT<STATE>& data) override{
        ComputeSolutionHCurlT(data.phi, data.curlphi, data.sol, data.curlsol);
    }
    void ReallyComputeSolution(TPZMaterialDataT<CSTATE>& data) override{
        ComputeSolutionHCurlT(data.phi, data.curlphi, data.sol, data.curlsol);
    }

    /**
     * This method computes the solution of an HCurl-conforming element and its derivatives in the local coordinate qsi
     * @param phi hcurl shape function
     * @param curlphi curl of the shape function
     * @param sol finite element solution
     * @param curlsol curl of the finite element solution
     */
    template<class TVar>
    void ComputeSolutionHCurlT(const TPZFMatrix<REAL> &phi,
                              const TPZFMatrix<REAL> &curlphi,
                              TPZSolVec<TVar> &sol,
                              TPZSolVec<TVar> &curlsol);

    int MaxOrder() override; 
};


/** @brief Creates computational linear element for HCurl-conforming approximation space */
TPZCompEl *CreateHCurlLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HCurlFamily hcurlfam);
/** @brief Creates computational quadrilateral element for HCurl-conforming approximation space */
TPZCompEl *CreateHCurlQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HCurlFamily hcurlfam);
/** @brief Creates computational triangular element for HCurl-conforming approximation space */
TPZCompEl *CreateHCurlTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HCurlFamily hcurlfam);
/** @brief Creates computational cube element for HCurl-conforming approximation space */
TPZCompEl *CreateHCurlCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HCurlFamily hcurlfam);
/** @brief Creates computational prismal element for HCurl-conforming approximation space */
TPZCompEl *CreateHCurlPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HCurlFamily hcurlfam);
/** @brief Creates computational pyramidal element for HCurl-conforming approximation space */
TPZCompEl *CreateHCurlPyramEl(TPZGeoEl *gel,TPZCompMesh &mesh);
/** @brief Creates computational tetrahedral element for HCurl-conforming approximation space */
TPZCompEl *CreateHCurlTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HCurlFamily hcurlfam);
/** @brief Creates computational point element for HCurl-conforming approximation space */
TPZCompEl *CreateHCurlBoundPointEl(TPZGeoEl *gel,TPZCompMesh &mesh);
/** @brief Creates computational linear element for HCurl-conforming approximation space */
TPZCompEl *CreateHCurlBoundLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HCurlFamily hcurlfam);
/** @brief Creates computational quadrilateral element for HCurl-conforming approximation space */
TPZCompEl *CreateHCurlBoundQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HCurlFamily hcurlfam);


/** @} */

#endif
