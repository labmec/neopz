//
//  TPZSBFemVolume.hpp
//  PZ
//
//  Created by Philippe Devloo on 4/4/16.
//
//

#ifndef TPZSBFemVolume_hpp
#define TPZSBFemVolume_hpp

#include <stdio.h>
#include "pzcompel.h"
#include "TPZElementMatrixT.h"
#include "pzinterpolationspace.h"
#include "TPZMaterialDataT.h"

class TPZSBFemVolume : public TPZInterpolationSpace
{
protected:
    /// index of element group
    int64_t fElementGroupIndex = -1;
    
    /// pointer to the element group computational element
    TPZCompEl *fElementGroup = 0;
    
    /// index of the skeleton element
    int64_t fSkeleton = -1;
    
    /// pointer to the integration rule
    TPZIntPoints *fIntRule = 0;
    
    /// Section of the phi vector associated with this volume element
    TPZFNMatrix<30,std::complex<double> > fPhi;
    
    /// Section of the phi vector associated with this volume element
    TPZFNMatrix<30,std::complex<REAL> > fPhiBubble;
    
    /// Eigenvlues associated with the internal shape functions
    TPZManVector<std::complex<double> > fEigenvalues;
    
    /// Eigenvlues associated with the internal shape functions
    TPZManVector<std::complex<REAL> > fEigenvaluesBubble;
    
    /// Inverse of fPhi and bubble coefficients
    TPZFNMatrix<100,std::complex<REAL> > fPhiInv;
    
    /// Inverse of fPhi and bubble coefficients
    TPZFNMatrix<100,std::complex<REAL> > fPhiInvBubbles;
    
    /// Multiplier coeficients associated with the solution
    TPZFNMatrix<30,std::complex<double> > fCoeficients;
    
    /// vector of local indices of multipliers in the group
    TPZManVector<int64_t> fLocalIndices;
    
    /// extend the border shape functions for SBFem computations
    void ExtendShapeFunctions(TPZMaterialDataT<STATE> &data1d, TPZMaterialDataT<STATE> &data2d);
    
    /// Density associated with the mass matrix
    REAL fDensity = 1.;

    /// adjust the axes and jacobian of the 3D element
    void AdjustAxes3D(const TPZFMatrix<REAL> &axes2D, TPZFMatrix<REAL> &axes3D, TPZFMatrix<REAL> &jac3D, TPZFMatrix<REAL> &jacinv3D, REAL detjac);
public:
    
    TPZSBFemVolume(TPZCompMesh &mesh, TPZGeoEl *gel);
    
    virtual ~TPZSBFemVolume()
    {
        Reference()->ResetReference();
    }
    
    /// Compute the E0, E1 and E2 matrices
    void ComputeKMatrices(TPZElementMatrixT<STATE> &E0, TPZElementMatrixT<STATE> &E1, TPZElementMatrixT<STATE> &E2, TPZElementMatrixT<STATE> &M0);
    
    /// Data structure initialization
    void SetSkeleton(int64_t skeleton);
    
    int64_t SkeletonIndex()
    {
        return fSkeleton;
    }
    
    /**
     * @brief Initialize a material data and its attributes based on element dimension, number
     * of state variables and material definitions
     */
    virtual void InitMaterialData(TPZMaterialData &data) override;

    /// Initialize the data structure indicating the group index
    void SetElementGroupIndex(int64_t index);
    
    /** @brief Method for creating a copy of the element */
    virtual TPZCompEl *Clone(TPZCompMesh &mesh) const override
    {
        // till I remember how this works
        DebugStop();
        return 0;
    }
    
    int64_t ElementGroupIndex() const
    {
        return fElementGroupIndex;
    }
    /**
     * @brief Method for creating a copy of the element in a patch mesh
     * @param mesh Patch clone mesh
     * @param gl2lcConMap map the connects indexes from global element (original) to the local copy.
     * @param gl2lcElMap map the computational elements
     */
    /**
     * Otherwise of the previous clone function, this method don't
     * copy entire mesh. Therefore it needs to map the connect index
     * from the both meshes - original and patch
     */
    virtual TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,
                                    std::map<int64_t,int64_t> & gl2lcConMap,
                                    std::map<int64_t,int64_t> & gl2lcElMap) const override
    {
        // till I remember how this works
        DebugStop();
        return 0;
    }

    /** @brief Returns the number of nodes of the element */
    virtual int NConnects() const override
    {
        if (fElementGroup == 0) {
            return 0;
        }
        return fElementGroup->NConnects();
    }
    
    /**
     * @brief Returns the index of the ith connectivity of the element
     * @param i connectivity index who want knows
     */
    virtual int64_t ConnectIndex(int i) const override
    {
        if (fElementGroup == 0) {
            DebugStop();
        }
        return fElementGroup->ConnectIndex(i);
    }
    /** @brief Dimension of the element */
    virtual int Dimension() const override
    {
        TPZGeoEl *reference = Reference();
        return reference->Dimension();
    }
    
    /** @brief return the density associated with the element */
    REAL Density()
    {
        return fDensity;
    }
    
    /** @brief assign a different density */
    void SetDensity(REAL density)
    {
#ifdef PZDEBUG
        if(density <= 0.) DebugStop();
#endif
        fDensity = density;
    }
    
    /** @brief adds the connect indexes associated with base shape functions to the set */
    virtual void BuildCornerConnectList(std::set<int64_t> &connectindexes) const override;
    
    /**
     * @brief Set the index i to node inode
     * @param inode node to set index
     * @param index index to be set
     */
    virtual void SetConnectIndex(int inode, int64_t index) override
    {
        DebugStop();
    }
    /** @brief Returns the number of dof nodes along side iside*/
    virtual int NSideConnects(int iside) const override
    {
        DebugStop();
		return 0;
    }
    
    /**
     * @brief Returns the local node number of icon along is
     * @param icon connect number along side is
     * @param is side which is being queried
     */
    virtual int SideConnectLocId(int icon,int is) const override
    {
        DebugStop();
		return 0;
    }

    /** @brief It returns the shapes number of the element */
    virtual int NShapeF() const override
    {
        int nc = NConnects();
        int nshape = 0;
        for (int ic=0; ic<nc; ic++) {
            TPZConnect &c = Connect(ic);
            nshape += c.NShape();
        }
        return nshape;
    }
    
    /** @brief Returns the number of shapefunctions associated with a connect*/
    virtual int NConnectShapeF(int icon, int order) const override
    {
        TPZConnect &c = Connect(icon);
        return c.NShape();
    }
    
    void InitializeIntegrationRule() override 
    {
        if (fIntRule) {
            DebugStop();
        }
        int nsides = Reference()->NSides();
        fIntRule = Reference()->CreateSideIntegrationRule(nsides-1, 1);
    }
    
    virtual void SetIntegrationRule(int order) override;
    
    /** @brief Returns a reference to an integration rule suitable for integrating the interior of the element */
    virtual const TPZIntPoints &GetIntegrationRule() const override
    {
        if(!fIntRule) DebugStop();

		return *fIntRule;
    }
    
    /** @brief Returns a reference to an integration rule suitable for integrating the interior of the element */
    virtual TPZIntPoints &GetIntegrationRule() override
    {
        if(!fIntRule) InitializeIntegrationRule();
        return *fIntRule;
    }
    
    /**
     * @brief Change the preferred order for the element and proceed the
     * adjust of the aproximation space \n taking in acount the type
     * of formulation and the neighbours of the element
     */
    virtual void PRefine ( int order ) override;

    /**  @brief Defines the desired order for entire element. */
    virtual void SetPreferredOrder ( int order ) override;
    


    /// initialize the data structures of the eigenvectors and eigenvalues associated with this volume element
    void SetPhiEigVal(TPZFMatrix<std::complex<double> > &phi, TPZManVector<std::complex<double> > &eigval);
    
    TPZFMatrix<std::complex<double> > Phi()
    {
        return fPhi;
    }
    
    TPZManVector<std::complex<double> > Eigenvalues()
    {
        return fEigenvalues;
    }
    
    TPZFMatrix<std::complex<double> > Coeficients()
    {
        return fCoeficients;
    }
    
    TPZFMatrix<double> PhiReal()
    {
        int64_t rows = fPhi.Rows(),cols = fPhi.Cols();
        TPZFMatrix<double> phireal(rows,cols);
        for(int64_t i=0; i<rows; i++)
        {
            for(int64_t j=0; j<cols; j++)
            {
                phireal(i,j) = fPhi(i,j).real();
            }
        }
        return phireal;
    }
    
    TPZManVector<double> EigenvaluesReal()
    {
        int64_t nel = fEigenvalues.NElements();
        TPZManVector<double> eig(nel);
        for(int64_t el=0; el<nel; el++)
        {
            eig[el] = fEigenvalues[el].real();
        }
        return eig;
    }
    
    TPZFMatrix<double> CoeficientsReal()
    {
        int64_t rows = fCoeficients.Rows(),cols = fCoeficients.Cols();
        TPZFMatrix<double> coefreal(rows,cols);
        for(int64_t i=0; i<rows; i++)
        {
            for(int64_t j=0; j<cols; j++)
            {
                coefreal(i,j) = fCoeficients(i,j).real();
            }
        }
        return coefreal;
    }
    
    /** @brief Loads the solution within the internal data structure of the element */
    /**
     * Is used to initialize the solution of connect objects with dependency. \n
     * Is also used to load the solution within SuperElements
     */
    virtual void LoadCoef(TPZFMatrix<std::complex<double> > &coef);
protected:
    //! Compute solution based on a filled TPZMaterialData
    void ReallyComputeSolution(TPZMaterialDataT<STATE>&) override;
public:
    
    void EvaluateError(TPZVec<REAL> &errors,bool store_error) override;
    
    /**
     * @brief Computes the shape function set at the point x.
     * @param qsi point in master element coordinates
     * @param phi vector of values of shapefunctions, dimension (numshape,1)
     * @param dphi matrix of derivatives of shapefunctions in master element coordinates, dimension (dim,numshape)
     */
    /**
     * This method uses the order of interpolation
     * of the element along the sides to compute the number of shapefunctions
     */
    virtual void Shape(TPZVec<REAL> &qsi,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphidxi) override;
    
    /** @brief Compute shape functions based on master element in the classical FEM manner. */
    virtual void ComputeShape(TPZVec<REAL> &intpoint, TPZVec<REAL> &X, TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes,
                              REAL &detjac, TPZFMatrix<REAL> &jacinv, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi, TPZFMatrix<REAL> &dphidx) override;
    

    /** @brief Compute the solution at the integration point and store in the data structure
     */

    /**
     * @brief Calculates the solution - sol - for the variable var
     * at point qsi, where qsi is expressed in terms of the
     * master element coordinates
     * @param qsi master element coordinate
     * @param var variable name
     * @param sol vetor for the solution
     */
    virtual void Solution(TPZVec<REAL> &qsi,int var,TPZVec<STATE> &sol) override;
    

    /**
     * @brief Prints element data
     * @param out Indicates the device where the data will be printed
     */
    virtual void Print(std::ostream &out = std::cout) const override
    {
        out << "Printing " << __PRETTY_FUNCTION__ << std::endl;
        TPZCompEl::Print(out);
        out << "Group Element Index " << fElementGroupIndex << std::endl;
        out << "Skeleton Element Index " << fSkeleton << std::endl;
        out << "Local Indices " << fLocalIndices << std::endl;
        fCoeficients.Print("Coef =",out,EMathematicaInput);
        fPhi.Print("Phi = ",out,EMathematicaInput);
        if (fCoeficients.Rows())
        {
            TPZManVector<std::complex<double>,5> prod(fPhi.Rows(),0.);
            for (int i=0; i<fPhi.Rows(); i++) {
                for (int j=0; j<fPhi.Cols(); j++) {
                    prod[i] += fPhi.GetVal(i,j)*fCoeficients.GetVal(j,0);
                }
            }
            out << "Values at border " << prod << std::endl;
        }
        for (int i=0; i<NConnects(); i++) {
            Connect(i).Print(*Mesh(),out);
        }
    }
    
    void CreateGraphicalElement(TPZGraphMesh &, int) override;

    void LocalBodyForces(TPZFNMatrix<100,std::complex<REAL>> &f, TPZFNMatrix<100,std::complex<REAL>> &fbubble, TPZManVector<std::complex<REAL>> &eigval, TPZManVector<std::complex<REAL>> &eigvalbubble, int icon);

    void ComputeSolutionWithBubbles(TPZVec<REAL> &qsi,
                                    TPZSolVec<STATE> &sol, TPZGradSolVec<STATE> &dsol, TPZFMatrix<REAL> &axes);
    
    void SetCoefNonHomogeneous(TPZFNMatrix<100,std::complex<double>> &phi, TPZManVector<std::complex<double> > &eigval, TPZFNMatrix<100,std::complex<double> > &phiinv, TPZFNMatrix<100,std::complex<double> > &rot);


};


TPZCompEl * CreateSBFemCompEl(TPZGeoEl *gel,TPZCompMesh &mesh);



#endif /* TPZSBFemVolume_hpp */
