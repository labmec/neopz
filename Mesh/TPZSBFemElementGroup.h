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


class TPZSBFemElementGroup : public TPZElementGroup
{
    
    /// Matrix of eigenvectors which compose the stiffness matrix
    TPZFMatrix<std::complex<double> > fPhi;
    
    /// Inverse of the eigenvector matrix (transfers eigenvector coeficients to side shape coeficients)
    TPZFMatrix<std::complex<double> > fPhiInverse;
    
    /// Vector of eigenvalues of the SBFem analyis
    TPZManVector<std::complex<double> > fEigenvalues;
    
    /// Multiplying coefficients of each eigenvector
    TPZFMatrix<std::complex<double> > fCoef;
    
    TPZFMatrix<STATE> fMassMatrix;
    
    /// Compute the mass matrix based on the value of M0 and the eigenvectors
    void ComputeMassMatrix(TPZElementMatrix &M0);
    
public:
    
    TPZSBFemElementGroup() : TPZElementGroup()
    {
        
    }
    
    /// constructor
    TPZSBFemElementGroup(TPZCompMesh &mesh, long &index) : TPZElementGroup(mesh,index)
    {
        
    }
    
    /// Compute the SBFem matrices
    /// method to assemble E0, E1, E2
    void ComputeMatrices(TPZElementMatrix &E0, TPZElementMatrix &E1, TPZElementMatrix &E2, TPZElementMatrix &M0);
    
    /**
     * @brief Computes the element stifness matrix and right hand side
     * @param ek element stiffness matrix
     * @param ef element load vector
     */
    virtual void CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef);
    
    /**
     * @brief Prints element data
     * @param out Indicates the device where the data will be printed
     */
    virtual void Print(std::ostream &out = std::cout) const
    {
        out << __PRETTY_FUNCTION__ << std::endl;
        TPZElementGroup::Print(out);
        int nel = fElGroup.size();
        out << "Element indexes of the volume elements ";
        for (int el=0; el<nel; el++) {
            out << fElGroup[el]->Index() << " ";
        }
        out << std::endl;
        out << "Indices of the associated skeleton elements\n";
        for (int el=0; el<nel; el++) {
            TPZCompEl *cel = fElGroup[el];
            TPZSBFemVolume *vol = dynamic_cast<TPZSBFemVolume *>(cel);
            if(!vol) DebugStop();
            out << vol->SkeletonIndex() << " ";
        }
        out << std::endl;
        out << "Connect indexes of the contained elements\n";
        for (int el=0; el<nel; el++) {
            TPZCompEl *cel = fElGroup[el];
            int nc = cel->NConnects();
            for (int ic=0; ic<nc; ic++) {
                out << cel->ConnectIndex(ic) << " ";
            }
            out << std::endl;
        }

/*
        out << "EigenVectors for displacement\n";
        fPhi.Print("Phi = ",out,EMathematicaInput);
        out << "Inverse EigenVectors\n";
        fPhiInverse.Print("PhiInv = ",out,EMathematicaInput);
        out << "EigenValues " << fEigenvalues << std::endl;
        out << "Mass Matrix\n";
        fMassMatrix.Print("Mass = ",out);
        out << "Solution Coeficients\n";
        fCoef.Print("Coef ",out);
        for (int el=0; el<nel; el++) {
            fElGroup[el]->Print(out);
        }
 */
        out << "End of " << __PRETTY_FUNCTION__ << std::endl;
    }
    

    
    /**
     * @brief Computes the element right hand side
     * @param ef element load vector(s)
     */
    virtual void CalcResidual(TPZElementMatrix &ef)
    {
        TPZElementMatrix ek(Mesh(),TPZElementMatrix::EK);
        CalcStiff(ek,ef);
    }
    

    /** @brief Loads the solution within the internal data structure of the element */
    /**
     * Is used to initialize the solution of connect objects with dependency. \n
     * Is also used to load the solution within SuperElements
     */
    virtual void LoadSolution();

    /** @brief Loads the geometric element referece */
    virtual void LoadElementReference()
    {
        for (long i = 0; i < fElGroup.size(); i++) {
            fElGroup[i]->LoadElementReference();
        }
    }
    


    long NumEigenValues()
    {
        return fEigenvalues.size();
    }
    
    /// Load the coeficients such that we visualize an eigenvector
    void LoadEigenVector(long eig);
    /// method to compute the stiffness
    /// method to compute the solution
    
    TPZVec<STATE> MultiplyingCoeficients()
    {
        int nel = fCoef.Rows();
        TPZVec<STATE> result(nel,0.);
        for (int ir=0; ir<nel; ir++) {
            result[ir] = fCoef(ir,0).real();
        }
        return result;
    }
    
    TPZVec<std::complex<double> > &EigenValues()
    {
        return fEigenvalues;
    }
    
    TPZFMatrix<std::complex<double> > Phi()
    {
        return fPhi;
    }
    
    
    TPZFMatrix<std::complex<double> > Coeficients()
    {
        return fCoef;
    }
    
    TPZFMatrix<double> PhiReal()
    {
        long rows = fPhi.Rows(),cols = fPhi.Cols();
        TPZFMatrix<double> phireal(rows,cols);
        for(long i=0; i<rows; i++)
        {
            for(long j=0; j<cols; j++)
            {
                phireal(i,j) = fPhi(i,j).real();
            }
        }
        return phireal;
    }
    
    TPZManVector<double> EigenvaluesReal()
    {
        long nel = fEigenvalues.NElements();
        TPZManVector<double> eig(nel);
        for(long el=0; el<nel; el++)
        {
            eig[el] = fEigenvalues[el].real();
        }
        return eig;
    }
    
    TPZFMatrix<double> CoeficientsReal()
    {
        long rows = fCoef.Rows(),cols = fCoef.Cols();
        TPZFMatrix<double> coefreal(rows,cols);
        for(long i=0; i<rows; i++)
        {
            for(long j=0; j<cols; j++)
            {
                coefreal(i,j) = fCoef(i,j).real();
            }
        }
        return coefreal;
    }

};

#endif /* TPZSBFemElementGroup_hpp */
