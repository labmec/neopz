//
//  TPZProjectEllipse.h
//  PZ
//
//  Created by Philippe Devloo on 11/4/13.
//
//

#ifndef __PZ__TPZProjectEllipse__
#define __PZ__TPZProjectEllipse__

#include <iostream>

#include "pzvec.h"
#include "pzfmatrix.h"

class TPZProjectEllipse
{
    /// Points which determine the projection
    TPZFMatrix<REAL> fPoints;
    /// Projection type
    int fType;
    /// coefficients determined by the least squares problem
    TPZManVector<REAL,6> fcoefficients;
    
    /// Tolerance to project the coefficients to zero value
    REAL fTol;

public:
    
    enum EType {EVerySimple2D, ESimple2D, E3D, E2D};
    
    TPZProjectEllipse(EType project, TPZFMatrix<REAL> &Points)
    {
        fType = project;
        fPoints = Points;
        int ncoef = 0;
        switch (project) {
            case EVerySimple2D:
                ncoef = 2;
                fcoefficients.resize(ncoef);
                LeastSquaresToGetVerySimpleEllipse();
                break;
            case ESimple2D:
                ncoef = 4;
                fcoefficients.resize(4);
                LeastSquaresToGetSimpleEllipse();
                break;
            default:
                DebugStop();
                break;
        }
    }
    
    TPZProjectEllipse(std::multimap<REAL, REAL> &Points)
    {
        fType = EVerySimple2D;
        fcoefficients.resize(2);
        fPoints.Resize(Points.size(), 2);
        std::multimap<REAL, REAL>::iterator it;
        int i=0;
        for (it = Points.begin(); it != Points.end(); it++, i++) {
            fPoints(i,0) = it->first;
            fPoints(i,1) = it->second;
        }
        LeastSquaresToGetVerySimpleEllipse();
    }

    
    void Getcoefficients(TPZVec<REAL> &coef)
    {
        coef = fcoefficients;
    }
    
    void PrintingAsSimpleEquation(TPZVec<REAL> &Center,TPZVec<REAL> &Ratios,std::ostream &out);

    bool StandardFormatForSimpleEllipse(TPZVec<REAL> &Center,TPZVec<REAL> &Ratios);

private:
    
    bool LeastSquaresToGetVerySimpleEllipse();
    // Least Squares Method to compute a ellipse nearest for a points in vector
    // The obtained ellipse has the axes parallels to rectangular axes
    // Format (x-x0)^2/a^2 + (y-y0)^2/b^2 = 1.
    bool AdjustingWithSimpleEllipse();
    // Least Squares Method to compute a ellipse nearest for a points in vector
    // The obtained ellipse has the axes parallels to rectangular axes
    // Format (x)^2/a^2 + (y)^2/b^2 = 1.
    bool AdjustingWithVerySimpleEllipse();
    
    // Least Squares Method to compute a ellipse nearest for a points in vector
    // The ellipse is a conic with second order equation as
    // y^2 = A*x^2 + B*xy + C*x + D*y + E  -> 2D case
    // z^2 = A*x^2 + B*xy + C*y^2 + D*yz + E*xz + F*x + G*y + H*z + I  -> 3D case
    // Then their axes could be rotated and translated
    bool AdjustingWithEllipse();
    
    // Diagonalize a quadratic form making avoid the quadratic term xy
    bool DiagonalizingQuadraticForm(int dim,TPZFMatrix<REAL> &Coeffs,TPZFMatrix<REAL> &NewCoeffs);
    
    // To print as zero all the values almost zero
    void AlmostZeroToZero(TPZFMatrix<REAL> &mat);
    void AlmostZeroToZero(TPZVec<REAL> &mat);
    
    //REAL Tol = 1.e-4;
    
    
    bool LeastSquaresToGetSimpleEllipse();
    
    /** Find coefficients of the quadratic equation (Quadratic form) with best adjust to points given
     * y^2 = a x^2 + b x + c xy + d y + e
     * z^2 = a x^2 + b x + c y^2 + d y + e xy + f z + g yz + h xz + i
     * using the least square method.
     */
    bool LeastSquaresToGetEllipse();
    

    /** Print the unitary vectors for rotated axes to ellipse */
    void PrintAxes(TPZFMatrix<REAL> &P,std::ostream &out);
    
   // bool DiagonalizingQuadraticForm(TPZVec<REAL> &NewCoeffs,std::ostream &out);
    

};

#endif /* defined(__PZ__TPZProjectEllipse__) */
