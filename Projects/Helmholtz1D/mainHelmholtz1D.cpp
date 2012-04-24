/*
 * @file mainHelmholtz1D.h
 * @since 4/4/12.
 */

#include "pzhelmholtz1D.h"
#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "PZ_Process.h"
#include "pzlog.h"

// Right handside term of our Linear PDE
void ForcingFunction(const TPZVec<REAL> &pto, TPZVec<REAL> &xfloat) {
	
	if(!xfloat.NElements() || !pto.NElements())
		DebugStop();
	
	xfloat[0] = pto[0];
}

// Exact Solution u(x)
void ExactSolution(TPZVec<REAL> &pto, TPZVec<REAL> &u_exact,TPZFMatrix<REAL> &du_exact) {

}

// Variable coefficient of first derivative
void AlphaFunction(const TPZVec<REAL> &pto, TPZVec<REAL> &alpha) {
#ifdef DEBUG
        if(alpha.NElements() != 2) {            
            DebugStop();
        }            
#endif
        std::complex<REAL> m(2,-0.1);
        //std::complex<REAL> m(1.,0.);
        std::complex<REAL> eps = 4. + m * (1. - pto[0] / L) * (1. - pto[0] / L);
        //std::complex<REAL> eps(1., 0.);
        std::complex<REAL> alph = 1. / eps;   
        
        alpha[0] = alph.real();
        alpha[1] = alph.imag();
}

// Variable coefficient of the u
void BetaFunction(const TPZVec<REAL> &pto, TPZVec<REAL> &beta) {
#ifdef DEBUG
        if(beta.NElements() != 2) {            
            DebugStop();
        }            
#endif
        REAL k0 = 2. * M_PI / lambda;
        std::complex<REAL> m(2,-0.1);
        //std::complex<REAL> m(1.,0.);
        std::complex<REAL> eps = 4. + m * (1. - pto[0] / L) * (1. - pto[0] / L);
        //std::complex<REAL> eps(1., 0.);
        std::complex<REAL> alph = 1. / eps;
        std::complex<REAL> bet = -(k0 * k0) * (m - alph * (std::sin(theta) * std::sin(theta)));
        
        beta[0] = bet.real();
        beta[1] = bet.imag();
}

void PhiFunction(const TPZVec<REAL> &pto, TPZVec<REAL> &phi){
#ifdef DEBUG
        if(phi.NElements() != 2) {
            
            DebugStop();
        }            
#endif   
        phi[0] = 0;
        phi[1] = 0;
}


int main() {
	
#ifdef LOG4CXX
	// Initializing generation log comments as log4cxx done
	std::string logs("log4cxx.doubleprojection1d");
	InitializePZLOG();
#endif
	// File for computing error
	std::ofstream FileError("erros.txt");
	// File for Mathematica output format
	std::ofstream outMath("Mathematica.txt");
	
	// p -> interpolation order
	int p = 1;
	// h -> level of uniform refinement of the initial mesh
	int h = 4;
	// Last point of the one-dimensional domain
	//double L = 1.;
	
	// Creating main extremes and material for current project
	TPZManVector<REAL> x0(3,0.), x1(3,0.);  // Corners of the mesh. Coordinates are zeros.
	x1[0] = L;
	// id of material and boundary conditions
	TPZVec<int> matId(1,1);
	TPZVec<int> bc(2,-1);
	bc[1] = -2;
	// type of boundary conditions
	TPZVec<int> bcType(2,0);   // 0 - dirichlet -> For x = 0 the boundary condition is Dirichlet
	bcType[1] = 2;             // Boundary condition for x = 1.
	
	// Material data
	TPZHelmholtz1D *material;
	material = new TPZHelmholtz1D(matId[0],1);
        
        material->SetAlphaFunction(new TPZDummyFunction(AlphaFunction));
	material->SetBetaFunction(new TPZDummyFunction(BetaFunction));
        material->SetPhiFunction(new TPZDummyFunction(PhiFunction));
        
        TPZFMatrix<REAL> xkin(2, 2, 0.), xcin(2, 2, 0.), xbin(2, 2, 0.), xfin(2, 1, 0.); 

        material->SetMaterial(xkin,xcin,xbin,xfin);
	
	// inserting function force
	material->SetForcingFunction(new TPZDummyFunction(ForcingFunction));
	
	// FEM PROCESS
	// Creating a geometric mesh and printing geometric information. Note: The coordinates of the nodes are 3D
	TPZGeoMesh *gmesh = GeomMesh1D(h,matId,bc,x0,x1);
	
	// Creating a computational mesh and printing computational information.
	TPZCompMesh * cmesh = CompMesh1D(gmesh,p,material,bc,bcType);
	
	// Assembling and Solving linear system
	TPZAnalysis an(cmesh);

	SolveSist(an,cmesh);
	// Solution output for Mathematica
	OutputMathematica(outMath,1,10, cmesh);
        OutputMathematica(outMath,2,10, cmesh);
	
	// Computing error
	an.SetExact(ExactSolution);
	TPZVec<REAL> posproc;
	an.PostProcess(posproc,FileError); // Compute the errors
	
	return 0;
}