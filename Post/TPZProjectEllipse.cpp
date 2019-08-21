//
//  TPZProjectEllipse.cpp
//  PZ
//
//  Created by Philippe Devloo on 11/4/13.
//
//

#include "TPZProjectEllipse.h"

/** Find coefficients of the quadratic equation with best adjust to points given
 * y^2 = a x^2 + b x + c y + d
 * z^2 = a x^2 + b x + c y^2 + d y + e z + f
 * using the least square method.
 */
bool TPZProjectEllipse::LeastSquaresToGetSimpleEllipse() {
	int64_t npoints = fPoints.Rows();
    int64_t dim = fPoints.Cols();
	int nincog, i;
	nincog = 2*dim;
    
	if(npoints<nincog) return false;
    
	// Dimensioning vector of coefficients
	fcoefficients.Resize(nincog);
	fcoefficients.Fill(0.);
    
	// Will be solved y^2 = p*x^2 + q*x + r*y + s
	// Constructing matrix H and Transpose of H to compute by least squares method
	TPZFMatrix<REAL> DeltaH;
	TPZFMatrix<REAL> DeltaHTranspose;
	TPZFMatrix<REAL> DifSol;
	TPZFMatrix<REAL> A;
    
	// Redimensioning
	A.Redim(nincog,nincog);
    DeltaH.Redim(npoints,nincog);
    DeltaHTranspose.Redim(nincog,npoints);
    DifSol.Redim(npoints,1);
    
	if(dim == 2) {
		// Filling y^2 into Coefficients
		for(i=0;i<npoints;i++)
			DifSol.PutVal(i,0,fPoints(i,1)*fPoints(i,1));  // fill y*y
        
		// Filling elements for H matrix
		for(int i=0;i<npoints;i++) {
			DeltaH.PutVal(i,0,fPoints(i,0)*fPoints(i,0));      // fill x*x
			DeltaH.PutVal(i,1,fPoints(i,0));                  // fill x
			DeltaH.PutVal(i,2,fPoints(i,0));                // fill y
			DeltaH.PutVal(i,3,1.);                           // fill 1.
		}
	}
	else if(dim == 3) {
		// Filling y^2 into Coefficients
		for(i=0;i<npoints;i++)
			DifSol.PutVal(i,0,fPoints(i,2)*fPoints(i,2));   // fill z*z
        
		// Filling elements for H matrix
		for(int i=0;i<npoints;i++) {
			DeltaH.PutVal(i,0,fPoints(i,0)*fPoints(i,0));       // fill x*x
			DeltaH.PutVal(i,1,fPoints(i,0));                   // fill x
			DeltaH.PutVal(i,2,fPoints(i,1)*fPoints(i,1));   // fill y*y
			DeltaH.PutVal(i,3,fPoints(i,1));                 // fill y
			DeltaH.PutVal(i,4,fPoints(i,2));                 // fill z
			DeltaH.PutVal(i,5,1.);                            // fill 1.
		}
	}
	else
		return false;
    
    // Solving by least squares using product of matrix: DeltaH_t * DifSol = DeltaH_t * DeltaH * Coeffs(u)
    A.Zero();
    DeltaH.Transpose(&DeltaHTranspose);
    A = DeltaHTranspose*DeltaH;
    TPZFMatrix<REAL> Coefficients;
	Coefficients = DeltaHTranspose*DifSol;
	A.SolveDirect(Coefficients,ELU);
    for (int i=0; i<Coefficients.Rows(); i++) {
        fcoefficients[i] = Coefficients(i,0);
    }
    return true;
}

/** Find coefficients of the quadratic equation (Quadratic form) with best adjust to points given
 * y^2 = a x^2 + b x + c xy + d y + e
 * z^2 = a x^2 + b x + c y^2 + d y + e xy + f z + g yz + h xz + i
 * using the least square method.
 */
bool TPZProjectEllipse::LeastSquaresToGetEllipse() {
    int dim = fPoints.Cols();
	int64_t npoints = fPoints.Rows();
	int nincog, i;
	nincog = 2*dim+1+2*(dim-2);
    
	if(npoints<nincog) return false;
    
	// Dimensioning vector of coefficients
	fcoefficients.Resize(nincog);
	fcoefficients.Fill(0.);
    
	// Will be solved y^2 = p*x^2 + q*x + r*y + s
	// Constructing matrix H and Transpose of H to compute by least squares method
	TPZFMatrix<REAL> DeltaH;
	TPZFMatrix<REAL> DeltaHTranspose;
	TPZFMatrix<REAL> DifSol;
	TPZFMatrix<REAL> A;
    
	// Redimensioning
	A.Redim(nincog,nincog);
    DeltaH.Redim(npoints,nincog);
    DeltaHTranspose.Redim(nincog,npoints);
    DifSol.Redim(npoints,1);
    
	if(dim == 2) {
		// Filling y^2 into Coefficients
		for(i=0;i<npoints;i++)
			DifSol.PutVal(i,0,fPoints(i,1)*fPoints(i,1));  // fill y*y
        
		// Filling elements for H matrix
		for(int i=0;i<npoints;i++) {
			DeltaH.PutVal(i,0,fPoints(i,0)*fPoints(i,0));      // fill x*x
			DeltaH.PutVal(i,1,fPoints(i,0));                  // fill x
			DeltaH.PutVal(i,2,fPoints(i,0)*fPoints(i,1));    // fill x*y
			DeltaH.PutVal(i,3,fPoints(i,1));                // fill y
			DeltaH.PutVal(i,4,1.);                           // fill 1.
		}
	}
	else if(dim == 3) {
		// Filling y^2 into Coefficients
		for(i=0;i<npoints;i++)
			DifSol.PutVal(i,0,fPoints(i,2)*fPoints(i,2));   // fill z*z
        
		// Filling elements for H matrix
		for(int i=0;i<npoints;i++) {
			DeltaH.PutVal(i,0,fPoints(i,0)*fPoints(i,0));       // fill x*x
			DeltaH.PutVal(i,1,fPoints(i,0));                   // fill x
			DeltaH.PutVal(i,2,fPoints(i,1)*fPoints(i,1));   // fill y*y
			DeltaH.PutVal(i,3,fPoints(i,1));                 // fill y
			DeltaH.PutVal(i,4,fPoints(i,0)*fPoints(i,1));     // fill x*y
			DeltaH.PutVal(i,5,fPoints(i,2));                 // fill z
			DeltaH.PutVal(i,6,fPoints(i,1)*fPoints(i,2));   // fill y*z
			DeltaH.PutVal(i,7,fPoints(i,0)*fPoints(i,2));     // fill x*z
			DeltaH.PutVal(i,8,1.);                            // fill 1.
		}
	}
	else
		return false;
    
    // Solving by least squares using product of matrix: DeltaH_t * DifSol = DeltaH_t * DeltaH * Coeffs(u)
    A.Zero();
    DeltaH.Transpose(&DeltaHTranspose);
    A = DeltaHTranspose*DeltaH;
    TPZFMatrix<REAL> Coefficients;
	Coefficients = DeltaHTranspose*DifSol;
	A.SolveDirect(Coefficients,ELU);
    
	// Verifying that the equation corresponding a ellipse
	REAL temp;
	if(dim==2) {
		temp = Coefficients(1,0)*Coefficients(1,0)/Coefficients(0,0);
		temp -= (Coefficients(2,0)*Coefficients(2,0));
		if(Coefficients(0,0) > 0 || (4*Coefficients(3,0)) < temp)
			return false;
	}
	else {   // It is not complete
		temp = Coefficients(1,0)*Coefficients(1,0)/Coefficients(0,0);
		if(Coefficients(0,0) > 0)
			return false;
	}
    for (int i=0; i<Coefficients.Rows(); i++) {
        fcoefficients[i] = Coefficients(i,0);
    }
    return true;
}


bool TPZProjectEllipse::LeastSquaresToGetVerySimpleEllipse() {
	int64_t npoints = fPoints.Rows();
	int nincog = 2, i;
    
	if(npoints<nincog) return false;
    
	// Dimensioning vector of coefficients
    TPZFMatrix<REAL> Coefficients;
	Coefficients.Redim(nincog,1);
	Coefficients.Zero();
    
	// Will be solved y^2 = p*x^2 + q*x + r*y + s
	// Constructing matrix H and Transpose of H to compute by least squares method
	TPZFMatrix<REAL> DeltaH;
	TPZFMatrix<REAL> DeltaHTranspose;
	TPZFMatrix<REAL> DifSol;
	TPZFMatrix<REAL> A;
    
	// Redimensioning
	A.Redim(nincog,nincog);
    DeltaH.Redim(npoints,nincog);
    DeltaHTranspose.Redim(nincog,npoints);
    DifSol.Redim(npoints,1);
    //	DifSol.Zero();
    
    // Filling y^2 into Coefficients
    for(i=0;i<npoints;i++)
        DifSol.PutVal(i,0,fPoints(i,1)*fPoints(i,1));  // fill y*y
    
    // Filling elements for H matrix
    for(int i=0;i<npoints;i++) {
        DeltaH.PutVal(i,0,fPoints(i,0)*fPoints(i,0));      // fill x*x
        DeltaH.PutVal(i,1,1.);                           // fill 1.
    }
    
//    DeltaH.Print(std::cout);
    
    // Solving by least squares using product of matrix: DeltaH_t * DifSol = DeltaH_t * DeltaH * Coeffs(u)
    A.Zero();
    DeltaH.Transpose(&DeltaHTranspose);
    A = DeltaHTranspose*DeltaH;
	Coefficients = DeltaHTranspose*DifSol;
	A.SolveDirect(Coefficients,ELU);
    for (int i=0; i<Coefficients.Rows(); i++) {
        fcoefficients[i] = Coefficients(i,0);
    }
    return true;
}




/*********** Utilities *************/

// To print as zero all the values almost zero
void TPZProjectEllipse::AlmostZeroToZero(TPZFMatrix<REAL> &mat) {
	int nr, nc;
	for(nr=0;nr<mat.Rows();nr++) {
		for(nc=0;nc<mat.Cols();nc++) {
			if(fabs(mat.GetVal(nr,nc)) < fTol)
				mat.PutVal(nr,nc,0.);
		}
	}
}
void TPZProjectEllipse::AlmostZeroToZero(TPZVec<REAL> &mat) {
	int nr;
	for(nr=0;nr<mat.NElements();nr++) {
		if(fabs(mat[nr]) < fTol)
			mat[nr] = 0.;
	}
}

/** Print the unitary vectors for rotated axes to ellipse */
void TPZProjectEllipse::PrintAxes(TPZFMatrix<REAL> &P,std::ostream &out) {
	out << "\nPrinting axes of the ellipse:\n";
	AlmostZeroToZero(P);
	for(int i=0;i<P.Cols();i++) {
		out << "Vector" << i+1 << " = ( ";
		for(int j=0;j<P.Rows();j++) {
			if(j) out << " , ";
			out << P.GetVal(j,i);
		}
		out << " )\n";
	}
	out << std::endl;
}

/****************** END Utilities ************************/


/** From quadratic equation (Quadratic form)
 * y^2 = a x^2 + b x + c xy + d y + e
 * z^2 = a x^2 + b x + c y^2 + d y + e xy + f z + g yz + h xz + i
 * compute the coordinates of the center and axes of the ellipse. Print in the standard equation.
 */

bool TPZProjectEllipse::StandardFormatForSimpleEllipse(TPZVec<REAL> &Center,TPZVec<REAL> &Ratios) {
    //	int dim = Center.NElements();
	int ncoeffs = fcoefficients.size();
	REAL temp, temp1;
	if(ncoeffs ==4) {
		Center[0] = -(fcoefficients[1]/(2.*fcoefficients[0]));
		Center[1] = 0.5*fcoefficients[2];
		// Computing Ratios[1] in temp
		temp = fcoefficients[1]*fcoefficients[1]-fcoefficients[0]*fcoefficients[2]*fcoefficients[2]-4*fcoefficients[0]*fcoefficients[3];
		// Computing Ratios
		temp1 = temp/(4*fcoefficients[0]*fcoefficients[0]);
		if(temp1 < 0. || IsZero(temp1))
			return false;
		Ratios[0] = sqrt(temp1);
		temp1 = (-temp)/(4*fcoefficients[0]);
		Ratios[1] = sqrt(temp1);
	}
	else if(ncoeffs == 6) {
		Center[0] = -(fcoefficients[1]/(2.*fcoefficients[0]));
		Center[1] = -(fcoefficients[3]/(2.*fcoefficients[2]));
		Center[2] = 0.5*fcoefficients[4];
		// Computing Ratios[2]
		temp = fcoefficients[5]+(Center[2]*Center[2])-(Center[0]*Center[0]*fcoefficients[0])-(Center[1]*Center[1]*fcoefficients[2]);
		if(temp < 0.) return false;
		Ratios[2] = sqrt(temp);
		// Computing Ratios[0] in Ratios[1]
		temp = -(Ratios[2]*Ratios[2]/fcoefficients[2]);
		if(temp < 0.) return false;
		Ratios[1] = sqrt(temp);
		temp = -(Ratios[2]*Ratios[2]/fcoefficients[0]);
		if(temp < 0.) return false;
		Ratios[0] = sqrt(temp);
	}
    else if(ncoeffs == 2)
    {
        Center[0] = 0.;
		Center[1] = 0.;
		// Computing Ratios[1] in temp
		temp = fcoefficients[1]-(Center[0]*Center[0]*fcoefficients[0])+(Center[1]*Center[1]);
		// Computing Ratios[0] in Ratios[1]
		Ratios[1] = -temp/fcoefficients[0];
		if(temp < 0. || Ratios[1] < 0.)
			return false;
		Ratios[0] = sqrt(Ratios[1]);
		Ratios[1] = sqrt(temp);
    }
	return true;
}

//void PrintingAsSimpleEquation(TPZVec<REAL> &Center,TPZVec<REAL> &Ratios,std::ostream &out);

void TPZProjectEllipse::PrintingAsSimpleEquation(TPZVec<REAL> &Center,TPZVec<REAL> &Ratios, std::ostream &out) {
    if(fcoefficients.size() == 2)
    {
		out << std::endl << "y*y = " << fcoefficients[0] << "x*x + " << fcoefficients[1] << "\n";
		out << "\nElipse: (x - " << Center[0] << ")^2/" << Ratios[0]*Ratios[0] << " + (y - " << Center[1] << ")^2/" << Ratios[1]*Ratios[1] << " = 1.\n" << std::endl;
        
    }
	else if(fcoefficients.size() == 4) {
		out << std::endl << "y*y = " << fcoefficients[0] << "x*x + " << fcoefficients[1] << "x + " << fcoefficients[2] << "y + " << fcoefficients[3] << "\n";
		out << "\nElipse: (x - " << Center[0] << ")^2/" << Ratios[0]*Ratios[0] << " + (y - " << Center[1] << ")^2/" << Ratios[1]*Ratios[1] << " = 1.\n" << std::endl;
	}
	else if(fcoefficients.size() == 6)
    {
		out << std::endl << "z*z = " << fcoefficients[0] << "x*x + " << fcoefficients[1] << "x + " << fcoefficients[2] << "y*y + " << fcoefficients[3] << "y +";
		out << fcoefficients[4] << "z + " << fcoefficients[5] << std::endl;
		out << "\nElipse: (x - " << Center[0] << ")^2/" << Ratios[0]*Ratios[0] << " + (y - " << Center[1] << ")^2/" << Ratios[1]*Ratios[1];
		out << " + (z - " << Center[2] << ")^2/" << Ratios[2]*Ratios[2] << " = 1.\n" << std::endl;
	}
}

/**
 * Process to get the ellipse by least square method from a set of points.
 * First: Find coefficients of the simple quadratic form from a set of points using least square method.
 * Second: From coefficients get the center and axes of the ellipse adjusted
 * Third: Print the result and standard form into saida.
 */
bool TPZProjectEllipse::AdjustingWithSimpleEllipse() {
    
    int dim = fPoints.Cols();
	TPZManVector<REAL> Center(dim,0.);
	TPZManVector<REAL> Ratios(dim,0.);
    
	// Applying least squares for these five points
	if(!LeastSquaresToGetSimpleEllipse())
		return false;
    std::cout << "\n\nCoefficients after least square for quadratic equation:";
    
	// Making zero depending on Tolerance
	float Tol;
	ZeroTolerance(Tol);
    
	// Getting the values of the center and ratios for the ellipse from Coeffs values
	if(!StandardFormatForSimpleEllipse(Center,Ratios))
		return false;
	PrintingAsSimpleEquation(Center,Ratios,std::cout);
	return true;
}


bool TPZProjectEllipse::AdjustingWithVerySimpleEllipse() {
    
    int dim = 2;
	TPZManVector<REAL> Center(dim,0.);
	TPZManVector<REAL> Ratios(dim,0.);
    
	// Applying least squares for these five points
	LeastSquaresToGetVerySimpleEllipse();
	std::cout << "\n\nSolution:";
    
	// Making zero depending on Tolerance
	float Tol;
	ZeroTolerance(Tol);
    
	// Getting the values of the center and ratios for the ellipse from Coeffs values
	if(!StandardFormatForSimpleEllipse(Center,Ratios))
		return false;
	PrintingAsSimpleEquation(Center,Ratios,std::cout);
	return true;
}

/**
 * Process to get the ellipse by least square method from a set of points.
 * First: Find coefficients of the quadratic form (with mixed terms) from a set of points using least square method.
 * Second: Apply diagonalization to transform quadratic form with mixed terms to quadratic form with only quadratic terms of the variables.
 * Third: From new coefficients get the center and axes of the ellipse adjusted
 * Fourth: Print the result and standard form into saida.
 */
bool TPZProjectEllipse::AdjustingWithEllipse() {
    
    int dim = fPoints.Cols();
	TPZManVector<REAL> Center(dim,0.);
	TPZManVector<REAL> Ratios(dim,0.);
	TPZManVector<REAL> VectorI(dim,0.);
	TPZManVector<REAL> VectorJ(dim,0.);
	TPZManVector<REAL> VectorK(dim,0.);
	// Applying least squares for these five points
	if(!LeastSquaresToGetEllipse())
		return false;
    
	// In case a12 = 0, then don't exist rotation then it is a ellipse with parallel axes to cartesian axes
	if(dim==2 && IsZero(fcoefficients[2])) {
		TPZManVector<REAL> NewCoeffs(4,0.);
        NewCoeffs = fcoefficients;
		// Making zero depending on Tolerance
		float Tol;
		ZeroTolerance(Tol);
        
		// Getting the values of the center and ratios for the ellipse from Coeffs values
		if(!StandardFormatForSimpleEllipse(Center,Ratios))
			return false;
		PrintingAsSimpleEquation(Center,Ratios,std::cout);
		return true;
	}
    
	// Diagonalization of the quadratic form when Coeffs(2,0) is not null
	TPZManVector<REAL> NewCoeffs;
    DebugStop();
	//DiagonalizingQuadraticForm(NewCoeffs,std::cout);
    
	// The coefficients are multiplied by (-1/lambda2)
	if(dim==2) {
		fcoefficients[0] = -(NewCoeffs[0]/NewCoeffs[2]);
		fcoefficients[1] = -(NewCoeffs[1]/NewCoeffs[2]);
		fcoefficients[2] = -(NewCoeffs[3]/NewCoeffs[2]);
		fcoefficients[3] = -(NewCoeffs[4]/NewCoeffs[2]);
		// Getting the values of the center and ratios for the ellipse from Coeffs values
		if(!StandardFormatForSimpleEllipse(Center,Ratios))
			return false;
		PrintingAsSimpleEquation(Center,Ratios,std::cout);
	}
    
	return true;
}

/** Apply diagonalization process to transform a general quadratic form(with mixed terms) to simple quadratic form (without mixed terms)
 * that is, apply rotation on the axes.
 */
/*
 bool TPZProjectEllipse::DiagonalizingQuadraticForm(TPZVec<REAL> &NewCoeffs,std::ostream &out) {
 int dim = fPoints.Cols();
 NewCoeffs.resize(2*dim+1);
 // Changing signal
 fcoefficients *= (-1.);
 // Constructing a symmetric matrix
 TPZFMatrix<REAL> A(dim,dim);
 if(dim == 2) {
 A(0,0) = Coeffs(0,0);
 A(0,1) = A(1,0) = 0.5*Coeffs(2,0);
 A(1,1) = 1.;
 }
 else {
 A(0,0) = Coeffs(0,0);
 A(0,1) = A(1,0) = 0.5*Coeffs(4,0);
 A(1,1) = Coeffs(2,0);
 A(2,2) = 1.;
 A(0,2) = A(2,0) = 0.5*Coeffs(7,0);
 A(1,2) = A(2,1) = 0.5*Coeffs(6,0);
 }
 
 // Store the coefficients of the variables in the homogenous equation
 TPZFMatrix<REAL> B(dim,1);
 int nr, nc;
 REAL F = Coeffs.GetVal(4+(dim-2)*4,0);
 for(nr=0;nr<dim;nr++)
 B.PutVal(nr,0,Coeffs.GetVal(2*nr+1,0));
 
 // Computing eigenvalues and eigenvectors
 TPZVec<REAL> Eigenvalues(dim);
 Coeffs.Redim(dim,dim);
 REAL Tol, temp, norm;
 ZeroTolerance(Tol);
 int64_t niter = 1000;
 
 if(!A.SolveEigensystemJacobi(niter,Tol,Eigenvalues,Coeffs))
 return false;                            // Could be some eigenvector a null vector
 
 // Verifying Eigenvalues must to be positives to be ellipse or ellipsoid
 for(nr=0;nr<dim;nr++) {
 if(Eigenvalues[nr] > 0.) continue;
 return false;  // If some eigenvalue is zero it is parabole, if some of these are negative hyperbol
 }
 // Normalizing eigenvectors in matrix
 for(nc=0;nc<dim;nc++) {
 temp = 0.;
 for(nr=0;nr<dim;nr++)
 temp += Coeffs(nr,nc)*Coeffs(nr,nc);
 norm = sqrt(temp);
 for(nr=0;nr<dim;nr++)
 Coeffs(nr,nc) *= (1./norm);
 }
 // The transpose of the ortogonal matrix is not necessary, it exists, is enough
 
 // Print the direction vectors for axes of the ellipse
 PrintAxes(Coeffs,out);
 
 // Coefficients of the squares of the variables
 for(nr=0;nr<dim;nr++)
 NewCoeffs.PutVal(2*nr,0,Eigenvalues[nr]);
 
 // Coefficients of the variables (linear terms)
 if(dim==2) {
 temp = B(0,0)*Coeffs(0,0) + B(1,0)*Coeffs(0,1);
 NewCoeffs.PutVal(1,0,temp);
 temp = B(0,0)*Coeffs(1,0) + B(1,0)*Coeffs(1,1);
 NewCoeffs.PutVal(3,0,temp);
 NewCoeffs.PutVal(4,0,F);
 }
 else if(dim==3) {
 temp = B(0,0)*Coeffs(0,0) + B(1,0)*Coeffs(0,1) + B(2,0)*Coeffs(0,2);
 NewCoeffs.PutVal(1,0,temp);
 //It is not complete
 }
 
 return true;
 }
 */
