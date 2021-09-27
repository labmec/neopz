/**
 * @file
 * @brief Contains declaration of TPZCompElHCurl class which implements a generic computational HCurl-conforming element
 */

#ifndef TPZCOMPELHCURL_H
#define TPZCOMPELHCURL_H

#include <pzelctemp.h>
#include "TPZMaterialDataT.h"
/**
 * @brief This class contains static methods pertinent to the implementation of different families of
 * HCurl-conforming approximation spaces.
 */
class TPZHCurlAuxClass{
public:
    enum EHCurlFamily{EFullOrder = 1};
private:
    static EHCurlFamily hCurlFamily;
public:
    static void SetHCurlFamily(EHCurlFamily fam){
        hCurlFamily = fam;
    }

    static EHCurlFamily GetHCurlFamily(){
        return hCurlFamily;
    }

    /**
     * @brief This computes the curl of a vector shape function implemented as \$f \phi\mathbf{v} \$f
     * @tparam TDIM Dimension of the element curl will have dim=2*TDIM-3 for TDIM=2,3.
     * @param vecShapeIndex associated which scalar function multiplies which vector
     * @param dphi derivatives of the scalar functions on the reference element
     * @param masterDirections vectors that are combined with the scalar functions to create the vector shape functions
     * @param jacobian jacobian of the element mapping
     * @param detJac determinant of the jacobian
     * @param axes relates the directions \$f (\xi,\eta,\zeta) \rightarrow (x,y,z) \$f
     * @param curlPhi curl of the vector shape functions
     */
    template<int TDIM>
    static void ComputeCurl(const TPZVec<std::pair<int, int64_t>> &vecShapeIndex, const TPZFMatrix<REAL> &dphi,
                            const TPZFMatrix<REAL> &masterDirections, const TPZFMatrix<REAL> &jacobian,
                            REAL detJac, const TPZFMatrix<REAL> &axes, TPZFMatrix<REAL> &curlPhi);

    /**
     * @brief This method multiplies the scalar shape function with the directions in order to compute the Hcurl basis.
     * It should be called at the beginning of every TPZMaterial::Contribute, in order to obtain the functions at the
     * integration point.
     * @param vecShapeIndex associated which scalar function multiplies which vector
     * @param phi scalar shape functions on the reference element
     * @param deformedDirections directions on the deformed element obtained by covariant piola transform
     * @param phiHCurl hcurl shape functions
     */
    static void ComputeShape(const TPZVec<std::pair<int, int64_t>> &vecShapeIndex, const TPZFMatrix<REAL> &phi,
                            const TPZMatrix<REAL> &deformedDirections, TPZMatrix<REAL> &phiHCurl);
};

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
    
    /// vector describing the permutation associated with each side
    TPZManVector<int, TSHAPE::NSides - TSHAPE::NCornerNodes> fSidePermutation;

    //the hcurl vectors already taking into account the sides transformation ids, but on the master element
    TPZFNMatrix<TSHAPE::Dimension * TSHAPE::Dimension * TSHAPE::NSides,REAL> fMasterDirections;

public:
	TPZCompElHCurl(TPZCompMesh &mesh, TPZGeoEl *gel, int64_t &index);
	
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
	
	virtual ~TPZCompElHCurl();

    /** @brief Returns the unique identifier for reading/writing objects to streams */
    int ClassId() const override;
    //since TPZCompElHCurl is an abstract class
    static int StaticClassId();
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
    int NConnectShapeF(int connect, int order) const override = 0;

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

    /** @brief Compute and fill data with requested attributes */
	void ComputeRequiredData(TPZMaterialDataT<STATE> &data, TPZVec<REAL> &qsi) override{
        ComputeRequiredDataT(data,qsi);
    }
    void ComputeRequiredData(TPZMaterialDataT<CSTATE> &data, TPZVec<REAL> &qsi) override{
        ComputeRequiredDataT(data,qsi);
    }
    

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
    void SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) override = 0;
    //@}
protected:

    /**
     * @brief Creates the connects for a generic HCurl-conforming element.
     * This method must be called in the constructor of every child class of TPZCompElHCurl.
     * @param mesh Computational mesh in which the connects will be inserted
     */
    void CreateHCurlConnects(TPZCompMesh &mesh);

    void ReallyComputeSolution(TPZMaterialDataT<STATE>&) override;
    void ReallyComputeSolution(TPZMaterialDataT<CSTATE>&) override;
    /**
     * @brief This computes the curl of a vector shape function implemented as \$f \phi\mathbf{v} \$f
     * @tparam TDIM Dimension of the element curl will have dim=2*TDIM-3 for TDIM=2,3.
     * @param vecShapeIndex associated which scalar function multiplies which vector
     * @param dphi derivatives of the scalar functions on the reference element
     * @param masterDirections vectors that are combined with the scalar functions to create the vector shape functions
     * @param jacobian jacobian of the element mapping
     * @param detJac determinant of the jacobian
     * @param axes relates the directions \$f (\xi,\eta,\zeta) \rightarrow (x,y,z) \$f
     * @param curlPhi curl of the vector shape functions
     */
    template<int TDIM=TSHAPE::Dimension>
    void ComputeCurl(const TPZVec<std::pair<int, int64_t>> &vecShapeIndex, const TPZFMatrix<REAL> &dphi,
                     const TPZFMatrix<REAL> &masterDirections, const TPZFMatrix<REAL> &jacobian,
                     REAL detJac, const TPZFMatrix<REAL> &axes, TPZFMatrix<REAL> &curlPhi){
        TPZHCurlAuxClass::ComputeCurl<TDIM>(vecShapeIndex,dphi,masterDirections,jacobian,detJac,axes,curlPhi);
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
    template<class TVar>
    void ComputeRequiredDataT(TPZMaterialDataT<TVar> &data, TPZVec<REAL> &qsi);


    /**
    * @brief Returns a matrix index of the shape and vector  associate to element
    * The directions are calculated based on the LOCAL side ids (and SideNodeLocId), such as the H1 shape functions.
    * For instance, for the triangle, the vectors are:
    * vea30 vea31 vea41 vea42 vea52 vea51 vet3 vet4 vet5 vfe63 vfe64 vfe65 vft1 vft2
    * and the shapes will be organized as follows:
    * phi0 phi1 phi2  phi31 phi32 ... phi3i   phi41 phi42 ... phi4j   phi51 phi52 ... phi5k   phi61 phi62 ... phi6n
    *^^^^^^^^^^^^^^^  ^^^^^^^^^^^^^^^^^^^^^   ^^^^^^^^^^^^^^^^^^^^^   ^^^^^^^^^^^^^^^^^^^^^   ^^^^^^^^^^^^^^^^^^^^^
    *  corner funcs       edge 3 funcs            edge 4 funcs            edge 5 funcs            internal funcs
    *
    * In order to ensure that the functions will coincide in two neighbouring elements, they will be then sorted by
    * their sides' GLOBAL ids instead of their LOCAL ids
    * @param[out] IndexVecShape Indicates the pair vector/shape function that will construct the approximation space
    * @param[in] connectOrder Order of the connects
    */
    virtual void IndexShapeToVec(TPZVec<std::pair<int,int64_t> > & indexVecShape, const TPZVec<int>& connectOrder) const = 0;

    /**
       @brief Computes the needed order of the H1 functions for creating the HCurl functions.
       @param[out] ord Order of the side connects different from the corner connects to compute the desired h1 functions
       @note This method skips the h1 connects associated with vertices
     */
    virtual void CalcH1ShapeOrders(TPZVec<int> &ord) const = 0;


    void CalcShapeSideTraces(const int side, const TPZFMatrix<REAL> &phi,
//                             const TPZFMatrix<REAL> &curlPhi,
                             TPZFMatrix<REAL> &phiTrace
//                             ,TPZFMatrix<REAL> &curlTrace
                             ) const;

    //Computes the deformed directions for TPZMaterialData
    void ComputeDeformedDirections(TPZMaterialData &data);

};


/** @brief Creates computational linear element for HCurl-conforming approximation space */
TPZCompEl *CreateHCurlLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational quadrilateral element for HCurl-conforming approximation space */
TPZCompEl *CreateHCurlQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational triangular element for HCurl-conforming approximation space */
TPZCompEl *CreateHCurlTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational cube element for HCurl-conforming approximation space */
TPZCompEl *CreateHCurlCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational prismal element for HCurl-conforming approximation space */
TPZCompEl *CreateHCurlPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational pyramidal element for HCurl-conforming approximation space */
TPZCompEl *CreateHCurlPyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational tetrahedral element for HCurl-conforming approximation space */
TPZCompEl *CreateHCurlTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational point element for HCurl-conforming approximation space */
TPZCompEl *CreateHCurlBoundPointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational linear element for HCurl-conforming approximation space */
TPZCompEl *CreateHCurlBoundLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational quadrilateral element for HCurl-conforming approximation space */
TPZCompEl *CreateHCurlBoundQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational triangular element for HCurl-conforming approximation space */
TPZCompEl *CreateHCurlBoundTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);

TPZCompEl * CreateRefHCurlLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
TPZCompEl * CreateRefHCurlQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
TPZCompEl * CreateRefHCurlTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
TPZCompEl * CreateRefHCurlCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
TPZCompEl * CreateRefHCurlPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
TPZCompEl * CreateRefHCurlPyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
TPZCompEl * CreateRefHCurlTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);


/** @} */

#endif
