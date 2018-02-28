/**
 * @file
 * @brief Tutorial program showing a method to find a ellipse nearest for a set of points using least square method
 */

#include <iostream>

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "tpzcompmeshreferred.h"
#include "pzreferredcompel.h"
#include "pzcompel.h"
#include "pzbndcond.h"
#include "TPZInterfaceEl.h"
#include "pzinterpolationspace.h"
#include "pzpoisson3d.h"
#include "pzmat1dlin.h"

#include "pzgengrid.h"
#include "pzgradient.h"
#include "pzl2projection.h"

#include "tpzgeoelrefpattern.h"
#include "TPZGeoLinear.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "pzgeopoint.h"

#include "pzgeoelbc.h"

#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "pzstrmatrix.h"
#include "pzbstrmatrix.h"
#include "pzstepsolver.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"

#include "TPZVTKGeoMesh.h"

#include "TPZProjectEllipse.h"


#include <iostream>
#include <math.h>
using namespace std;

// Functions to fill coordinates for several points to test the adjusting by ellipse
void FillingPoints2D(TPZManVector<REAL> &Points);
void FillingPoints2D(TPZFMatrix<REAL> &Points);
void FillingPoints3D(TPZManVector<REAL> &Points);

// Least Squares Method to compute a ellipse nearest for a points in vector
// The obtained ellipse has the axes parallels to rectangular axes
// Format (x-x0)^2/a^2 + (y-y0)^2/b^2 = 1.
bool AdjustingWithSimpleEllipse(int dim,TPZManVector<REAL> &points,std::ostream &out=std::cout);
bool AdjustingWithVerySimpleEllipse(TPZManVector<REAL> &points);

// Least Squares Method to compute a ellipse nearest for a points in vector
// The ellipse is a conic with second order equation as
// y^2 = A*x^2 + B*xy + C*x + D*y + E  -> 2D case
// z^2 = A*x^2 + B*xy + C*y^2 + D*yz + E*xz + F*x + G*y + H*z + I  -> 3D case
// Then their axes could be rotated and translated
bool AdjustingWithEllipse(int dim,TPZManVector<REAL> &points,std::ostream &out=std::cout);

// Diagonalize a quadratic form making avoid the quadratic term xy
bool DiagonalizingQuadraticForm(int dim,TPZFMatrix<REAL> &Coeffs,TPZFMatrix<REAL> &NewCoeffs,std::ostream &out);

// To print as zero all the values almost zero
void AlmostZeroToZero(TPZFMatrix<REAL> &mat);
void AlmostZeroToZero(TPZVec<REAL> &mat);

REAL Tol = 1.e-4;
std::ofstream out("EllipseInfo.txt");

void LeastSquaresToGetVerySimpleEllipse(TPZManVector<REAL> &points,TPZFMatrix<REAL> &Coefficients);

int main(int argc, char *argv[]) {

#ifdef LOG4CXX
    InitializePZLOG();
#endif

    TPZFMatrix<REAL> matpoints;
    FillingPoints2D(matpoints);
    TPZProjectEllipse ellips(TPZProjectEllipse::EVerySimple2D,matpoints);
    TPZManVector<REAL,10> coefs;
    ellips.Getcoefficients(coefs);
    std::cout << "Coefficients " << coefs << std::endl;
    TPZManVector<REAL,3> center(2),ratios(2);
    ellips.StandardFormatForSimpleEllipse(center, ratios);
    std::cout << "center " << center << " ratios " << ratios << std::endl;
    ellips.PrintingAsSimpleEquation(center, ratios, std::cout);
    
	int dim = 2;
	TPZManVector<REAL> Points;
	if(dim==2)
		FillingPoints2D(Points);
	else
		FillingPoints3D(Points);

    if (dim == 2) {
        AdjustingWithVerySimpleEllipse(Points);
    }
	

	if(dim == 2) out << "Adjusting with ELLIPSE (axes parallels with cartesian axes):\n\n";
	else if(dim == 3) out << "Adjusting with ELLIPSOID (axes parallels with cartesian axes):\n\n";
	// Finding a ellipse nearest for all points
	if(!AdjustingWithSimpleEllipse(dim,Points,out)) {
		AdjustingWithEllipse(dim,Points,out);  // Whether don't exist ellipse with axes parallel to cartesian axes, the program compute arbitrary ellipse
		return 1;    // Return 1 means the simple ellipse couldn't be finded.
	}
    
	if(dim == 2) out << "\n\nAdjusting with ELLIPSE (Could be had rotation):\n\n";
	else if(dim == 3) out << "\n\nAdjusting with ELLIPSOID (Could be had rotation):\n\n";
	// Finding a ellipse nearest for all points
	if(dim == 2 && !AdjustingWithEllipse(dim,Points,out))
		return 2;   // Return 2 means any ellipse couldn't be finded.

	out.close();
	return 0;   // Success
}

/** Find coefficients of the quadratic equation with best adjust to points given
 * y^2 = a x^2 + b x + c y + d
 * z^2 = a x^2 + b x + c y^2 + d y + e z + f
 * using the least square method.
 */
bool LeastSquaresToGetSimpleEllipse(int dim,TPZManVector<REAL> &points,TPZFMatrix<REAL> &Coefficients) {
	if(dim < 2 || dim > 3) return false;
	int64_t npoints = points.NElements()/dim;
	int nincog, i;
	nincog = 2*dim;

	if(npoints<nincog) return false;

	// Dimensioning vector of coefficients
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

	if(dim == 2) {
		// Filling y^2 into Coefficients
		for(i=0;i<npoints;i++)
			DifSol.PutVal(i,0,points[2*i+1]*points[2*i+1]);  // fill y*y

		// Filling elements for H matrix
		for(int i=0;i<npoints;i++) {
			DeltaH.PutVal(i,0,points[2*i]*points[2*i]);      // fill x*x
			DeltaH.PutVal(i,1,points[2*i]);                  // fill x
			DeltaH.PutVal(i,2,points[2*i+1]);                // fill y
			DeltaH.PutVal(i,3,1.);                           // fill 1.
		}
	}
	else if(dim == 3) {
		// Filling y^2 into Coefficients
		for(i=0;i<npoints;i++)
			DifSol.PutVal(i,0,points[3*i+2]*points[3*i+2]);   // fill z*z

		// Filling elements for H matrix
		for(int i=0;i<npoints;i++) {
			DeltaH.PutVal(i,0,points[3*i]*points[3*i]);       // fill x*x
			DeltaH.PutVal(i,1,points[3*i]);                   // fill x
			DeltaH.PutVal(i,2,points[3*i+1]*points[3*i+1]);   // fill y*y
			DeltaH.PutVal(i,3,points[3*i+1]);                 // fill y
			DeltaH.PutVal(i,4,points[3*i+2]);                 // fill z
			DeltaH.PutVal(i,5,1.);                            // fill 1.
		}
	}
	else 
		return false;

    // Solving by least squares using product of matrix: DeltaH_t * DifSol = DeltaH_t * DeltaH * Coeffs(u)
    A.Zero();
    DeltaH.Transpose(&DeltaHTranspose);
    A = DeltaHTranspose*DeltaH;
	Coefficients = DeltaHTranspose*DifSol;
	A.SolveDirect(Coefficients,ELU);
	return true;
}

/** Find coefficients of the quadratic equation (Quadratic form) with best adjust to points given
 * y^2 = a x^2 + b x + c xy + d y + e
 * z^2 = a x^2 + b x + c y^2 + d y + e xy + f z + g yz + h xz + i
 * using the least square method.
 */
bool LeastSquaresToGetEllipse(int dim,TPZManVector<REAL> &points,TPZFMatrix<REAL> &Coefficients) {
	if(dim < 2 || dim > 3) return false;
	int64_t npoints = points.NElements()/dim;
	int nincog, i;
	nincog = 2*dim+1+2*(dim-2);

	if(npoints<nincog) return false;

	// Dimensioning vector of coefficients
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

	if(dim == 2) {
		// Filling y^2 into Coefficients
		for(i=0;i<npoints;i++)
			DifSol.PutVal(i,0,points[2*i+1]*points[2*i+1]);  // fill y*y

		// Filling elements for H matrix
		for(int i=0;i<npoints;i++) {
			DeltaH.PutVal(i,0,points[2*i]*points[2*i]);      // fill x*x
			DeltaH.PutVal(i,1,points[2*i]);                  // fill x
			DeltaH.PutVal(i,2,points[2*i]*points[2*i+1]);    // fill x*y
			DeltaH.PutVal(i,3,points[2*i+1]);                // fill y
			DeltaH.PutVal(i,4,1.);                           // fill 1.
		}
	}
	else if(dim == 3) {
		// Filling y^2 into Coefficients
		for(i=0;i<npoints;i++)
			DifSol.PutVal(i,0,points[3*i+2]*points[3*i+2]);   // fill z*z

		// Filling elements for H matrix
		for(int i=0;i<npoints;i++) {
			DeltaH.PutVal(i,0,points[3*i]*points[3*i]);       // fill x*x
			DeltaH.PutVal(i,1,points[3*i]);                   // fill x
			DeltaH.PutVal(i,2,points[3*i+1]*points[3*i+1]);   // fill y*y
			DeltaH.PutVal(i,3,points[3*i+1]);                 // fill y
			DeltaH.PutVal(i,4,points[3*i]*points[3*i+1]);     // fill x*y
			DeltaH.PutVal(i,5,points[3*i+2]);                 // fill z
			DeltaH.PutVal(i,6,points[3*i+1]*points[3*i+2]);   // fill y*z
			DeltaH.PutVal(i,7,points[3*i]*points[3*i+2]);     // fill x*z
			DeltaH.PutVal(i,8,1.);                            // fill 1.
		}
	}
	else 
		return false;

    // Solving by least squares using product of matrix: DeltaH_t * DifSol = DeltaH_t * DeltaH * Coeffs(u)
    A.Zero();
    DeltaH.Transpose(&DeltaHTranspose);
    A = DeltaHTranspose*DeltaH;
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
	return true;
}


void LeastSquaresToGetVerySimpleEllipse(TPZManVector<REAL> &points,TPZFMatrix<REAL> &Coefficients) {
	int64_t npoints = points.NElements()/2;
	int nincog = 2, i;
    
	if(npoints<nincog) DebugStop(); //return false;
    
	// Dimensioning vector of coefficients
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
        DifSol.PutVal(i,0,points[2*i+1]*points[2*i+1]);  // fill y*y
    
    // Filling elements for H matrix
    for(int i=0;i<npoints;i++) {
        DeltaH.PutVal(i,0,points[2*i]*points[2*i]);      // fill x*x
        DeltaH.PutVal(i,1,1.);                           // fill 1.
    }
    
	DeltaH.Print(std::cout);
    
    // Solving by least squares using product of matrix: DeltaH_t * DifSol = DeltaH_t * DeltaH * Coeffs(u)
    A.Zero();
    DeltaH.Transpose(&DeltaHTranspose);
    A = DeltaHTranspose*DeltaH;
	Coefficients = DeltaHTranspose*DifSol;
	A.SolveDirect(Coefficients,ELU);
	Coefficients.Print(std::cout);
}




/*********** Utilities *************/

// To print as zero all the values almost zero
void AlmostZeroToZero(TPZFMatrix<REAL> &mat) {
	int nr, nc;
	for(nr=0;nr<mat.Rows();nr++) {
		for(nc=0;nc<mat.Cols();nc++) {
			if(fabs(mat.GetVal(nr,nc)) < Tol)
				mat.PutVal(nr,nc,0.);
		}
	}
}
void AlmostZeroToZero(TPZVec<REAL> &mat) {
	int nr;
	for(nr=0;nr<mat.NElements();nr++) {
		if(fabs(mat[nr]) < Tol)
			mat[nr] = 0.;
	}
}
void PrintingAsSimpleEquation(TPZFMatrix<REAL> &Coeffs,TPZManVector<REAL> &Center,TPZManVector<REAL> &Ratios,std::ostream &out) {
	int dim = Center.NElements();
	AlmostZeroToZero(Coeffs);
	AlmostZeroToZero(Center);
	AlmostZeroToZero(Ratios);
	if(dim == 2) {
		out << std::endl << "y*y = " << Coeffs(0,0) << "x*x + " << Coeffs(1,0) << "x + " << Coeffs(2,0) << "y + " << Coeffs(3,0) << "\n";
		out << "\nElipse: (x - " << Center[0] << ")^2/" << Ratios[0]*Ratios[0] << " + (y - " << Center[1] << ")^2/" << Ratios[1]*Ratios[1] << " = 1.\n" << std::endl;
		out << "\nCenter = ( " << Center[0] << " , " << Center[1] << " ).\n" << "Axes =>   a = " << Ratios[0] << "\nand   b = " << Ratios[1] << std::endl;
	}
	else {
		out << std::endl << "z*z = " << Coeffs(0,0) << "x*x + " << Coeffs(1,0) << "x + " << Coeffs(2,0) << "y*y + " << Coeffs(3,0) << "y +";
		out << Coeffs(4,0) << "z + " << Coeffs(5,0) << std::endl;
		out << "\nElipse: (x - " << Center[0] << ")^2/" << Ratios[0]*Ratios[0] << " + (y - " << Center[1] << ")^2/" << Ratios[1]*Ratios[1];
		out << " + (z - " << Center[2] << ")^2/" << Ratios[2]*Ratios[2] << " = 1.\n" << std::endl;
		out << "\nCenter = ( " << Center[0] << " , " << Center[1] << " , " << Center[2] << " ).\n" << "Axes => a = " << Ratios[0] << "\nand   b = " << Ratios[1] << "\nand   c = " << Ratios[2] << std::endl;
	}
}

/** Print the unitary vectors for rotated axes to ellipse */
void PrintAxes(TPZFMatrix<REAL> &P,std::ostream &out) {
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

bool StandardFormatForSimpleEllipse(TPZFMatrix<REAL> &Coeffs,TPZManVector<REAL> &Center,TPZManVector<REAL> &Ratios) {
	int ncoeffs = Coeffs.Rows();
	REAL temp, temp1;
	if( Coeffs.Cols()!=1)
		return false;
	if(ncoeffs ==4) {
		Center[0] = -(Coeffs(1,0)/(2.*Coeffs(0,0)));
		Center[1] = 0.5*Coeffs(2,0);
		// Computing Ratios[1] in temp
		temp = Coeffs(1,0)*Coeffs(1,0)-Coeffs(0,0)*Coeffs(2,0)*Coeffs(2,0)-4*Coeffs(0,0)*Coeffs(3,0);
		// Computing Ratios
		temp1 = temp/(4*Coeffs(0,0)*Coeffs(0,0));
		if(temp1 < 0. || IsZero(temp1))
			return false;
		Ratios[0] = sqrt(temp1);
		temp1 = (-temp)/(4*Coeffs(0,0));
		Ratios[1] = sqrt(temp1);
	}
	else if(ncoeffs == 6) {
		Center[0] = -(Coeffs(1,0)/(2.*Coeffs(0,0)));
		Center[1] = -(Coeffs(3,0)/(2.*Coeffs(2,0)));
		Center[2] = 0.5*Coeffs(4,0);
		// Computing Ratios[2]
		temp = Coeffs(5,0)+(Center[2]*Center[2])-(Center[0]*Center[0]*Coeffs(0,0))-(Center[1]*Center[1]*Coeffs(2,0));
		if(temp < 0.) return false;
		Ratios[2] = sqrt(temp);
		// Computing Ratios[0] in Ratios[1]
		temp = -(Ratios[2]*Ratios[2]/Coeffs(2,0));
		if(temp < 0.) return false;
		Ratios[1] = sqrt(temp);
		temp = -(Ratios[2]*Ratios[2]/Coeffs(0,0));
		if(temp < 0.) return false;
		Ratios[0] = sqrt(temp);
	}
    else if(ncoeffs == 2)
    {
        Center[0] = 0.;
		Center[1] = 0.;
		// Computing Ratios[1] in temp
		temp = Coeffs(1,0)-(Center[0]*Center[0]*Coeffs(0,0))+(Center[1]*Center[1]);
		// Computing Ratios[0] in Ratios[1]
		Ratios[1] = -temp/Coeffs(0,0);
		if(temp < 0. || Ratios[1] < 0.)
			return false;
		Ratios[0] = sqrt(Ratios[1]);
		Ratios[1] = sqrt(temp);
    }
	return true;
}


void PrintingAsSimpleEquation(TPZFMatrix<REAL> &Coeffs,TPZManVector<REAL> &Center,TPZManVector<REAL> &Ratios) {
    if(Coeffs.Rows() == 2)
    {
		std::cout << std::endl << "y*y = " << Coeffs(0,0) << "x*x + " << Coeffs(1,0) << "\n";
		std::cout << "\nElipse: (x - " << Center[0] << ")^2/" << Ratios[0]*Ratios[0] << " + (y - " << Center[1] << ")^2/" << Ratios[1]*Ratios[1] << " = 1.\n" << std::endl;
        
    }
	else if(Coeffs.Rows() == 4) {
		std::cout << std::endl << "y*y = " << Coeffs(0,0) << "x*x + " << Coeffs(1,0) << "x + " << Coeffs(2,0) << "y + " << Coeffs(3,0) << "\n";
		std::cout << "\nElipse: (x - " << Center[0] << ")^2/" << Ratios[0]*Ratios[0] << " + (y - " << Center[1] << ")^2/" << Ratios[1]*Ratios[1] << " = 1.\n" << std::endl;
	}
	else if(Coeffs.Rows() == 6)
    {
		std::cout << std::endl << "z*z = " << Coeffs(0,0) << "x*x + " << Coeffs(1,0) << "x + " << Coeffs(2,0) << "y*y + " << Coeffs(3,0) << "y +";
		std::cout << Coeffs(4,0) << "z + " << Coeffs(5,0) << std::endl;
		std::cout << "\nElipse: (x - " << Center[0] << ")^2/" << Ratios[0]*Ratios[0] << " + (y - " << Center[1] << ")^2/" << Ratios[1]*Ratios[1];
		std::cout << " + (z - " << Center[2] << ")^2/" << Ratios[2]*Ratios[2] << " = 1.\n" << std::endl;
	}
}

/**
 * Process to get the ellipse by least square method from a set of points.
 * First: Find coefficients of the simple quadratic form from a set of points using least square method.
 * Second: From coefficients get the center and axes of the ellipse adjusted
 * Third: Print the result and standard form into saida.
 */
bool AdjustingWithSimpleEllipse(int dim,TPZManVector<REAL> &Points,std::ostream &out) {

	TPZFMatrix<REAL> Coeffs;
	TPZManVector<REAL> Center(dim,0.);
	TPZManVector<REAL> Ratios(dim,0.);

	// Applying least squares for these five points
	if(!LeastSquaresToGetSimpleEllipse(dim,Points,Coeffs))
		return false;
	out << "\n\nCoefficients after least square for quadratic equation:";

	// Making zero depending on Tolerance
	float Tol;
	ZeroTolerance(Tol);
	for(int i=0;i<Coeffs.Rows();i++)
		if(fabs(Coeffs(i)) < 100*Tol)
			Coeffs.PutVal(i,0,0.);

	// Getting the values of the center and ratios for the ellipse from Coeffs values
	if(!StandardFormatForSimpleEllipse(Coeffs,Center,Ratios))
		return false;
	PrintingAsSimpleEquation(Coeffs,Center,Ratios,out);
	return true;
}


bool AdjustingWithVerySimpleEllipse(TPZManVector<REAL> &Points) {
    
	TPZFMatrix<REAL> Coeffs;
    int dim = 2;
	TPZManVector<REAL> Center(dim,0.);
	TPZManVector<REAL> Ratios(dim,0.);
    
	// Applying least squares for these five points
	LeastSquaresToGetVerySimpleEllipse(Points,Coeffs);
	std::cout << "\n\nSolution:";
    
	// Making zero depending on Tolerance
	float Tol;
	ZeroTolerance(Tol);
	for(int i=0;i<Coeffs.Rows();i++)
		if(fabs(Coeffs(i)) < 100*Tol)
			Coeffs.PutVal(i,0,0.);
    
	// Getting the values of the center and ratios for the ellipse from Coeffs values
	if(!StandardFormatForSimpleEllipse(Coeffs,Center,Ratios))
		return false;
	PrintingAsSimpleEquation(Coeffs,Center,Ratios);
	return true;
}

/**
 * Process to get the ellipse by least square method from a set of points.
 * First: Find coefficients of the quadratic form (with mixed terms) from a set of points using least square method.
 * Second: Apply diagonalization to transform quadratic form with mixed terms to quadratic form with only quadratic terms of the variables.
 * Third: From new coefficients get the center and axes of the ellipse adjusted
 * Fourth: Print the result and standard form into saida.
 */
bool AdjustingWithEllipse(int dim,TPZManVector<REAL> &Points,std::ostream &out) {


	TPZFMatrix<REAL> Coeffs;
	TPZManVector<REAL> Center(dim,0.);
	TPZManVector<REAL> Ratios(dim,0.);
	TPZManVector<REAL> VectorI(dim,0.);
	TPZManVector<REAL> VectorJ(dim,0.);
	TPZManVector<REAL> VectorK(dim,0.);
	// Applying least squares for these five points
	if(!LeastSquaresToGetEllipse(dim,Points,Coeffs))
		return false;

	// In case a12 = 0, then don't exist rotation then it is a ellipse with parallel axes to cartesian axes
	if(dim==2 && IsZero(Coeffs(2,0))) {
		TPZFMatrix<REAL> NewCoeffs(4,1,0.);
		NewCoeffs.PutVal(0,0,Coeffs(0,0));
		NewCoeffs.PutVal(1,0,Coeffs(1,0));
		NewCoeffs.PutVal(2,0,Coeffs(3,0));
		NewCoeffs.PutVal(3,0,Coeffs(4,0));
		// Making zero depending on Tolerance
		float Tol;
		ZeroTolerance(Tol);
		for(int i=0;i<Coeffs.Rows();i++)
			if(fabs(Coeffs(i)) < 100*Tol)
				Coeffs.PutVal(i,0,0.);

		// Getting the values of the center and ratios for the ellipse from Coeffs values
		if(!StandardFormatForSimpleEllipse(Coeffs,Center,Ratios))
			return false;
		PrintingAsSimpleEquation(Coeffs,Center,Ratios,out);
		return true;
	}

	// Diagonalization of the quadratic form when Coeffs(2,0) is not null
	TPZFMatrix<REAL> NewCoeffs;
	DiagonalizingQuadraticForm(dim,Coeffs,NewCoeffs,out);
	Coeffs.Redim(2*dim,1);

	// The coefficients are multiplied by (-1/lambda2)
	if(dim==2) {
		Coeffs.PutVal(0,0,-(NewCoeffs(0,0)/NewCoeffs(2,0)));
		Coeffs.PutVal(1,0,-(NewCoeffs(1,0)/NewCoeffs(2,0)));
		Coeffs.PutVal(2,0,-(NewCoeffs(3,0)/NewCoeffs(2,0)));
		Coeffs.PutVal(3,0,-(NewCoeffs(4,0)/NewCoeffs(2,0)));
		// Getting the values of the center and ratios for the ellipse from Coeffs values
		if(!StandardFormatForSimpleEllipse(Coeffs,Center,Ratios))
			return false;
		PrintingAsSimpleEquation(Coeffs,Center,Ratios,out);
	}

	return true;
}

/** Apply diagonalization process to transform a general quadratic form(with mixed terms) to simple quadratic form (without mixed terms)
 * that is, apply rotation on the axes.
 */
bool DiagonalizingQuadraticForm(int dim,TPZFMatrix<REAL> &Coeffs,TPZFMatrix<REAL> &NewCoeffs,std::ostream &out) {
	NewCoeffs.Redim(2*dim+1,1);
	// Changing signal 
	Coeffs *= (-1.);
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

/************ For tests *****************************/
void FillingPoints2D(TPZFMatrix<REAL> &Points) {
    TPZManVector<REAL> vecpoints;
    FillingPoints2D(vecpoints);
    int np = vecpoints.size()/2;
    Points.Redim(np, 2);
    for (int i=0; i<np; i++) {
        Points(i,0) = vecpoints[2*i];
        Points(i,1) = vecpoints[2*i+1];
    }
    
}

/* Filling points 2D case */
void FillingPoints2D(TPZManVector<REAL> &Points) {
	int dim = 2;
	int npoints = 8;
	if(0) {
		Points.Resize(dim*npoints);

		// Filling points coordinates
		Points[0] = 1.;
		Points[1] = 3;
		Points[2] = 1.;
		Points[3] = -3;
		Points[4] = -1.;
		Points[5] = 0.;
		Points[6] = 3.;
		Points[7] = 0.;
		Points[8] = 0.;
		Points[9] = (3./2.)*sqrt(3.);
		Points[10] = 0.;
		Points[11] = -(3./2.)*sqrt(3.);
		Points[12] = 2.97203;
		Points[13] = 0.5;
		Points[14] = 2.49071;
		Points[15] = -2.;
	}
	else if(0) {
		// Consideremos a elipse (1./9.)*(x-2)^2 + (1./18.)*(y-1)^2 = 1
		// Podemos incrementar qualquer numero de pontos nas coordenadas x pertencente [-1.,5.]
		int i;
		double x[] = {-0.5, 0., 0.5, 1., 2., 3., 4., 4.5};
		double temp;
		Points.Resize(2*npoints*dim);
		for(i=0;i<npoints;i++) {
			Points[dim*i] = x[i];
			temp = sqrt(18 - 2*(x[i] - 2.)*(x[i] - 2.));
			Points[dim*i+1] = 1.+temp;
			Points[2*npoints+dim*i] = x[i];
			Points[2*npoints+dim*i+1] = 1.-temp;
		}
	}
	else if(0) {
		// Consideremos a elipse (1./9.)*(x-2)^2 + (1./18.)*(y-1)^2 = 1
		// Podemos incrementar qualquer numero de pontos nas coordenadas x pertencente [-1.,5.]
		// Verificando com valores deslocados muito pouco
		int i;
		double x[] = {-0.88, -0.5, 0.5, 1., 2., 3., 4., 4.5};
		double temp;
		Points.Resize(2*(npoints+1)*dim);
		for(i=0;i<npoints;i++) {
			Points[dim*i] = x[i];
			temp = sqrt(18 - 2*(x[i] - 2.)*(x[i] - 2.));
			Points[dim*i+1] = 1.+temp;
			Points[2*(npoints+1)+dim*i] = x[i];
			Points[2*(npoints+1)+dim*i+1] = 1.-temp;
		}
		Points[dim*i] = -1.1;
		Points[dim*i+1] = 1.03;
		Points[2*(npoints+1)+dim*i] = 2.;
		Points[2*(npoints+1)+dim*i+1] = 1.+sqrt(18.);
	}
	else if(0) {
		// Para um elipse com forma quadr�tica 2 x^2 + 2 xy + 2 y^2 = 9
		// with axes v1=(1/sqrt(2),-1/sqrt(2)) and v2 = (1/sqrt(2),1/sqrt(2))
		// Center = (0,0) and a = sqrt(3)  and  b = 3.  Vide Kolman (alg linear) pag 480.
		Points.Resize(26);
		Points[0] = Points[2] = 0.;
		Points[1] = 3./sqrt(2);
		Points[3] = -Points[1];
		Points[4] = 1.; Points[5] = 0.5*(-1.-sqrt(15));
		Points[6] = 1.; Points[7] = 0.5*(-1.+sqrt(15));
		Points[8] = .5; Points[9] = -2.32666;
		Points[10] = .5; Points[11] = 1.82666;
		Points[12] = -.5; Points[13] = 2.32666;
		Points[14] = -.5; Points[15] = -1.82666;
		Points[16] = 2.; Points[17] = 0.5*(-2.-sqrt(6));
		Points[18] = 2.; Points[19] = 0.5*(-2.+sqrt(6));
		Points[20] = -2.4; Points[21] = 0.775736;
		Points[22] = -2.4; Points[23] = 1.62426;
		Points[24] = -0.3; Points[25] = 2.25535;
	}
	else if(0) {
		// Para um elipse com forma quadr�tica 2 x^2 + 2 xy + 2 y^2 = 9
		// with axes v1=(1/sqrt(2),-1/sqrt(2)) and v2 = (1/sqrt(2),1/sqrt(2))
		// Center = (0,0) and a = sqrt(3)  and  b = 3.  Vide Kolman (alg linear) pag 480.
		// Nesse caso estamos fornecendo pontos sempre a um lado da elipse, para ver se d� a mesma, 
		// pois nos comentarios do Philippe a fratura (borda) tinha s� uma parte da elipse
		Points.Resize(26);
		Points[0] = 0.; Points[1] = 3./sqrt(2);
		Points[2] = 0.1; Points[3] = 2.06955;
		Points[4] = 1.; Points[5] = 0.5*(-1.+sqrt(15));
		Points[6] = 1.1; Points[7] = 1.34539;
		Points[8] = .45; Points[9] = 1.86022;
		Points[10] = .5; Points[11] = 1.82666;
		Points[12] = -.5; Points[13] = 2.32666;
		Points[14] = -.45; Points[15] = 2.31022;
		Points[16] = 1.8; Points[17] = 0.538749;
		Points[18] = 2.; Points[19] = 0.5*(-2.+sqrt(6));
		Points[20] = 0.3; Points[21] = 1.95535;
		Points[22] = -2.4; Points[23] = 1.62426;
		Points[24] = -0.3; Points[25] = 2.25535;
	}
	else if(0) {
		// Para um elipse com forma quadr�tica 2 x^2 + 2 xy + 2 y^2 = 9
		// with axes v1=(1/sqrt(2),-1/sqrt(2)) and v2 = (1/sqrt(2),1/sqrt(2))
		// Center = (sqrt(2),sqrt(2)) and a = sqrt(3)  and  b = 3.  Vide Kolman (alg linear) pag 480.
		// Nesse caso estamos fornecendo os mesmos pontos com a abscisa deslacada duas unidades
		Points.Resize(26);
		Points[0] = 2.; Points[1] = 3./sqrt(2);
		Points[2] = 2.1; Points[3] = 2.06955;
		Points[4] = 3.; Points[5] = 0.5*(-1.+sqrt(15));
		Points[6] = 3.1; Points[7] = 1.34539;
		Points[8] = 2.45; Points[9] = 1.86022;
		Points[10] = 2.5; Points[11] = 1.82666;
		Points[12] = 1.5; Points[13] = 2.32666;
		Points[14] = 1.55; Points[15] = 2.31022;
		Points[16] = 3.8; Points[17] = 0.538749;
		Points[18] = 4.; Points[19] = 0.5*(-2.+sqrt(6));
		Points[20] = 2.3; Points[21] = 1.95535;
		Points[22] = -.4; Points[23] = 1.62426;
		Points[24] = 1.7; Points[25] = 2.25535;
	}
	else {
		// Para um elipse com forma quadr�tica 2 x^2 + 2 xy + 2 y^2 = 9
		// with axes v1=(1/sqrt(2),-1/sqrt(2)) and v2 = (1/sqrt(2),1/sqrt(2))
		// Center = (sqrt(2),sqrt(2)) and a = sqrt(3)  and  b = 3.  Vide Kolman (alg linear) pag 480.
		// Nesse caso estamos fornecendo os mesmos pontos com a abscisa deslacada duas unidades
		Points.Resize(16);
		Points[0] = 2.; Points[1] = 3./sqrt(2);
		Points[2] = 2.1; Points[3] = 2.06955;
		Points[4] = 3.; Points[5] = 0.5*(-1.+sqrt(15));
		Points[6] = 3.1; Points[7] = 1.34539;
		Points[8] = 2.45; Points[9] = 1.86022;
		Points[10] = 2.5; Points[11] = 1.82666;
		Points[12] = 1.55; Points[13] = 2.31022;
		Points[14] = 2.3; Points[15] = 1.95535;
	}
}
/* Filling points 3D case */
void FillingPoints3D(TPZManVector<REAL> &Points) {
	Points.Resize(18);

	// Filling points coordinates
	Points[0] = 1.;
	Points[1] = 0.;
	Points[2] = 1.;

	Points[3] = 1.;
	Points[4] = 0.;
	Points[5] = -3.;

	Points[6] = -1.;
	Points[7] = 0.;
	Points[8] = -1.;

	Points[9] = 3.;
	Points[10] = 0.;
	Points[11] = -1.;

	Points[12] = 1.;
	Points[13] = 3.;
	Points[14] = -1.;

	Points[15] = 1.;
	Points[16] = -3.;
	Points[17] = -1.;
}
