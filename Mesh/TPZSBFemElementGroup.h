//
//  TPZSBFemElementGroup.hpp
//  PZ
//
//  Created by Philippe Devloo on 4/4/16.
//
//

#ifndef TPZSBFemElementGroup_hpp
#define TPZSBFemElementGroup_hpp

#include <stdio.h>

#include "pzelementgroup.h"
#include "TPZSBFemVolume.h"
#include "pzcmesh.h"
#include "pzcondensedcompel.h"


class TPZSBFemElementGroup : public TPZElementGroup
{
    
public:
    enum EComputationMode {EStiff, EOnlyMass, EMass, EStiffBubble};

private:
    /// Default polynomial order for internal bubble functions
    // if its value is zero, there are no internal functions
    static int gDefaultPolynomialOrder;
    
    static bool gPolynomialShapeFunctions;
    
private:
    
    /// Matrix of eigenvectors which compose the stiffness matrix
    TPZFNMatrix<100,std::complex<double> > fPhi;

    /// Matrix of eigenvectors which compose the stiffness matrix
    TPZFNMatrix<100,std::complex<double> > fPhiBubble;
    
    /// Inverse of the eigenvector matrix (transfers eigenvector coeficients to side shape coeficients)
    TPZFNMatrix<100,std::complex<double> > fPhiInverse;

    /// Matrix that composes the bubble functions
    TPZFNMatrix<100,std::complex<double> > fMatBubble;
    
    /// Vector of eigenvalues of the SBFem analyis
    TPZManVector<std::complex<double> > fEigenvalues;
    
    /// Vector of eigenvalues of the SBFem analyis
    TPZManVector<std::complex<double> > fEigenvaluesBubble;

    /// Stiffness matrix related to the bubble functions
    TPZFNMatrix<100, CSTATE > fStiffBubble;
    
    /// Right hand side of the bubble functions
    TPZFNMatrix<100, CSTATE > fRhsBubble;

    /// Matrix of eigenvectors which compose the stiffness matrix
    TPZFNMatrix<100,std::complex<double> > fQVectors;
    
    /// Multiplying coefficients of each eigenvector
    TPZFNMatrix<100,std::complex<double> > fCoef;
    
    TPZFMatrix<STATE> fMassMatrix;
    
    EComputationMode fComputationMode = EStiff;
    
    /// multiplier to multiply the mass matrix
    REAL fMassDensity = 1.;
    
    /// timestep coeficient
    REAL fDelt = 1.;

    int fInternalPolynomialOrder = 0;
        
    bool fPolynomialShapeFunctions = false;

    int64_t fInternalConnectIndex = 0;

    /// @brief Flag to indicate the eigenmodes have been computed
    bool fEigenComputed = false;
    
    /// @brief Compute the eigenmodes of the element group
    void ComputeEigenmodes();

    /// Compute the mass matrix based on the value of M0 and the eigenvectors
    void ComputeMassMatrix(TPZElementMatrixT<STATE> &M0);

    /// @brief Computes the right hand side vector
    /// @param ef element load vector(s)
    /// this method is private because it is called by CalcStiff and assumes ef has been initialized to zero
    void ComputeRhs(TPZElementMatrixT<STATE> &ef);

    /// @brief Compute the fRhs datastructure
    void ComputeRhsBubble();
public:
    
    /// constructor
    TPZSBFemElementGroup() : TPZElementGroup()
    {
        
    }

    virtual ~TPZSBFemElementGroup()
    {
        
    }
    
    /// constructor
    TPZSBFemElementGroup(TPZCompMesh &mesh);
    
    /** @brief add an element to the element group
     */
    virtual void AddElement(TPZCompEl *cel) override;
    

    /// Compute the SBFem matrices
    /// method to assemble E0, E1, E2
    /// if withBC is false, only SBFemVolume elements of mesh dimension will be assembled
    void ComputeMatrices(TPZElementMatrixT<STATE> &E0, TPZElementMatrixT<STATE> &E1, TPZElementMatrixT<STATE> &E2, TPZElementMatrixT<STATE> &M0);
    
    /**
     * @brief Computes the element stifness matrix and right hand side
     * @param ek element stiffness matrix
     * @param ef element load vector
     */
    virtual void CalcStiff(TPZElementMatrixT<STATE> &ek,TPZElementMatrixT<STATE> &ef) override;

    void CalcStiffBlaze(TPZElementMatrixT<STATE> &ek,TPZElementMatrixT<STATE> &ef);

    /// copy the values of the stiffness matrix to the first block of stiff
    void ContributeStiffness(TPZFMatrix<STATE> &stiff);
    
    /// compute the integral of the shape functions
    void ComputeShapeFunctionIntegral(TPZVec<STATE> &phiint);
    
    /// Get the polynomial order of the internal bubble functions
    int InternalPolynomialOrder() const
    {
        return fInternalPolynomialOrder;
    }

    static int GetDefaultPolynomialOrder()
    {
        return gDefaultPolynomialOrder;
    }
    static void SetDefaultPolynomialOrder(int order)
    {
        gDefaultPolynomialOrder = order;
    }
    void SetInternalPolynomialOrder(int order)
    {
        fInternalPolynomialOrder = order;
    }
    /// set the density or specific heat of the material
    void SetDensity(REAL density)
    {
        fMassDensity = density;
    }
    /// Set the element to compute the mass matrix
    void SetComputeOnlyMassMatrix()
    {
        fComputationMode = EOnlyMass;
        SetEigenComputed(false);
    }
    
    /// Set the element to compute stiffness plus mass
    void SetComputeTimeDependent(REAL delt)
    {
        fDelt = delt;
        fComputationMode = EMass;
        SetEigenComputed(false);
    }
    
    void SetComputeStiff()
    {
        fComputationMode = EStiff;
        SetEigenComputed(false);
    }

    void SetComputeFullBubbleStiff()
    {
        fComputationMode = EStiffBubble;
        SetEigenComputed(false);
    }

    void SetEigenComputed(bool computed)
    {
        fEigenComputed = computed;
    }
    
    /// set the right hand side for the bubble functions
    void SetBubbleRhs(TPZFMatrix<CSTATE> &rhsbubble) {
        if(fRhsBubble.Rows() != rhsbubble.Rows()) {
            DebugStop();
        }
        fRhsBubble = rhsbubble;
    }
    /**
     * @brief Prints element data
     * @param out Indicates the device where the data will be printed
     */
    virtual void Print(std::ostream &out = std::cout) const override;

    
    /**
     * @brief Computes the element right hand side
     * @param ef element load vector(s)
     */
    virtual void CalcResidual(TPZElementMatrixT<STATE> &ef) override
    {
        TPZElementMatrixT<STATE> ek(Mesh(),TPZElementMatrixT<STATE>::EK);
        CalcStiff(ek,ef);
    }
    

    /** @brief Loads the solution within the internal data structure of the element */
    /**
     * Is used to initialize the solution of connect objects with dependency. \n
     * Is also used to load the solution within SuperElements
     */
    virtual void LoadSolution() override;
    
    /// @brief Load the solution using the given solution vector
    void LoadSolution(TPZFMatrix<CSTATE> &sol, TPZFMatrix<CSTATE> &rhsbubble);

    /** @brief Loads the geometric element referece */
    virtual void LoadElementReference() override
    {
        for (int64_t i = 0; i < fElGroup.size(); i++) {
            fElGroup[i]->LoadElementReference();
        }
    }
    
    
    int64_t NumEigenValues()
    {
        return fEigenvalues.size();
    }
    
    int64_t NumEigenValuesBubble() {
        return fEigenvaluesBubble.size();
    }
    
    /// Load the coeficients such that we visualize an eigenvector
    void LoadEigenVector(int64_t eig);
    /// method to compute the stiffness
    /// method to compute the solution
    
    TPZManVector<std::complex<REAL> > &EigenValues()
    {
        return fEigenvalues;
    }
    
    TPZManVector<std::complex<REAL> > &EigenValuesBubble()
    {
        return fEigenvaluesBubble;
    }
    
    TPZFMatrix<std::complex<REAL> > &Phi()
    {
        return fPhi;
    }
    
    TPZFMatrix<std::complex<REAL> > &PhiBubble()
    {
        return fPhiBubble;
    }
    
    TPZFMatrix<std::complex<REAL> > &PhiInverse()
    {
        return fPhiInverse;
    }
    
    TPZFMatrix<CSTATE> &MatBubble() {
        return fMatBubble;
    }
    
    TPZFMatrix<STATE> &MassMatrix()
    {
        return fMassMatrix;
    }
    
    TPZFMatrix<std::complex<REAL> > Coeficients()
    {
        return fCoef;
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
    
    TPZManVector<REAL> EigenvaluesReal()
    {
        int64_t nel = fEigenvalues.NElements();
        TPZManVector<double> eig(nel);
        for(int64_t el=0; el<nel; el++)
        {
            eig[el] = fEigenvalues[el].real();
        }
        return eig;
    }
    
    TPZManVector<REAL> EigenvaluesBubbleReal()
    {
        int64_t nel = fEigenvaluesBubble.NElements();
        TPZManVector<double> eig(nel);
        for(int64_t el=0; el<nel; el++)
        {
            eig[el] = fEigenvaluesBubble[el].real();
        }
        return eig;
    }
    
    TPZFMatrix<REAL> CoeficientsReal()
    {
        int64_t rows = fCoef.Rows(),cols = fCoef.Cols();
        TPZFMatrix<double> coefreal(rows,cols);
        for(int64_t i=0; i<rows; i++)
        {
            for(int64_t j=0; j<cols; j++)
            {
                coefreal(i,j) = fCoef(i,j).real();
            }
        }
        return coefreal;
    }

    // Initialize the connects related to the bubble functions
    void InitializeInternalConnect();

    // Update the fEigenvalues, fPhi and fPhiInverse with the values needed to construct the bubble functions
    void ComputeBubbleParameters();

    // Overwrite eigenvalues and eigenvectors to use a polynomial approx (i.e, uses a collapsed FE approx instead of the SBFEM approximation)
    void OverwritePhis(TPZElementMatrixT<STATE> &E0, TPZElementMatrixT<STATE> &E1, TPZElementMatrixT<STATE> &E2, TPZElementMatrixT<STATE> &ek,TPZElementMatrixT<STATE> &ef);

    // 
    void SolveEigenProblemSBFEM(TPZFMatrix<STATE> &globmatkeep, TPZManVector<std::complex<double> > &eigenvalues, TPZFNMatrix<100,std::complex<double> > &eigenvectors);
    
    /// @brief Compute the stiffness of Laplace equation using numerical integration
    /// this method is to test whether the derivative of the shape functions is computed correctly
    void ComputeStiffnessLaplace(TPZFMatrix<STATE> &ek);
};

#endif /* TPZSBFemElementGroup_hpp */
