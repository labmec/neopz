//
//  TPZSBFemElementGroup.hpp
//  PZ
//
//  Created by Karolinne Coelho on 22/01/21.
//
//

#pragma once

#include <stdio.h>
#include "TPZSBFemElementGroup.h"
#include "pzcondensedcompel.h"

using namespace std;

class TPZSBFemMultiphysicsElGroup : public TPZSBFemElementGroup
{
    
private:
    
    TPZElementGroup * fCondensedEls;

    TPZCondensedCompEl * fCondEl;
    
    /// Matrix of eigenvectors which compose the stiffness matrix
    TPZFMatrix<std::complex<double>> fPhi;
    
    /// Matrix of eigenvectors which compose the stiffness matrix
    TPZFMatrix<std::complex<double>> fPhiQ;
    
    /// Matrix of eigenvectors which compose the stiffness matrix for the flux
    TPZFMatrix<std::complex<double>> fPhiFlux;
    
    /// Inverse of the eigenvector matrix (transfers eigenvector coeficients to side shape coeficients)
    TPZFMatrix<std::complex<double>> fPhiInverse;
    
    /// Vector of eigenvalues of the SBFem analyis
    TPZManVector<std::complex<double>> fEigenvalues;
    
    /// Multiplying coefficients of each eigenvector
    TPZFMatrix<std::complex<double>> fCoef;
    
    /// Multiplying coefficients of each eigenvector
    TPZFMatrix<std::complex<double>> fCoefD;
    
    TPZManVector<int64_t> fLocalindices;

public:
    
    TPZSBFemMultiphysicsElGroup() : TPZSBFemElementGroup()
    {
        
    }
    
    /// constructor
    TPZSBFemMultiphysicsElGroup(TPZCompMesh &mesh) : TPZSBFemElementGroup(mesh)
    { 
    }

    TPZCondensedCompEl * CondensedCompEl()
    {
        return fCondEl;
    }

    void AddElement(TPZCompEl *cel) override;

    virtual void Print(std::ostream &out) const override;

    void GroupandCondense(set<int> & matidscondensed);


    void ComputeMatrices(TPZElementMatrixT<STATE> &E0, TPZElementMatrixT<STATE> &E1, TPZElementMatrixT<STATE> &E2);
    
    /**
     * @brief Computes the element stifness matrix and right hand side
     * @param ek element stiffness matrix
     * @param ef element load vector
     */
    virtual void CalcStiff(TPZElementMatrixT<STATE> &ek,TPZElementMatrixT<STATE> &ef) override;


    void ComputeMatrices(TPZElementMatrix &E0, TPZElementMatrix &E1, TPZElementMatrix &E2);


    virtual void BuildCornerConnectList(std::set<int64_t> &connectindexes) const override
    {
        for (auto sbfemvol : fElGroup)
        {
            sbfemvol->BuildCornerConnectList(connectindexes);
        }
    }

    void InitializeElementMatrix(TPZElementMatrixT<STATE> &ek, TPZElementMatrixT<STATE> &ef) const;

    void LoadSolution() override;

    void LoadFlux();
    
    TPZFMatrix<std::complex<double>> &PhiInverse()
    {
        return fPhiInverse;
    }

    void AdjustConnectivities();

    void SetLocalIndices(int64_t index);
    
    void SetLocalIndicesFlux(int64_t index);

    TPZManVector<double> EigenvaluesReal(TPZManVector<complex<double> > & eigenvalues)
    {
        int64_t nel = eigenvalues.NElements();
        TPZManVector<double> eig(nel);
        for(int64_t el=0; el<nel; el++)
        {
            eig[el] = eigenvalues[el].real();
        }
        return eig;
    }

    TPZFMatrix<double> EigenVectorsReal(TPZFMatrix<std::complex<double>> &eigvec)
    {
        auto nrows = eigvec.Rows();
        auto ncols = eigvec.Cols();
        TPZFMatrix<double> eig(nrows,ncols);

        for (auto ir = 0; ir < nrows; ir++)
        {
            for (auto ic = 0; ic < ncols; ic++)
            {
                eig(ir,ic) = eigvec(ir,ic).real();
            }
        }
        
        return eig;
    }
};
