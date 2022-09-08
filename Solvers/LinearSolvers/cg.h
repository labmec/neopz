/**
 * @file
 * @brief Contains the implementation of the CG function which solves the symmetric positive definite linear system
 * using the Conjugate Gradient method.
 */

/**
 * @ingroup solver
 * @brief CG solves the symmetric positive definite linear system \f$ A x = b \f$ using the Conjugate Gradient method.
 * @return The return value indicates convergence within max_iter (input) iterations (0), or no convergence within max_iter iterations (1). \n
 * Upon successful return, output arguments have the following values: 
 * @param A  -- matrix of the system
 * @param b  -- vector of the system
 * @param M  -- preconditioner matrix
 * @param x  --  approximate solution to Ax = b
 * @param max_iter  --  the number of iterations performed before the tolerance was reached
 * @param tol  --  the residual after the final iteration
 * @param residual  -- residual vector (return)
 * @param FromCurrent  -- for type of operation (MultAdd)
*/
/**
 * Iterative template routine -- CG \n
 * CG follows the algorithm described on p. 15 in the SIAM Templates book.
 */
#define TEST
#ifdef TEST
#include <list> 
#include <fstream>
#include "TPZLapackEigenSolver.h"
#include "pzreal.h"
#endif
#include <iostream>
#include "pzreal.h"
#include "pzextractval.h"

template < class Matrix, class Vector, class Preconditioner, class Real >
int
CG( Matrix &A, Vector &x, const Vector &b,
   Preconditioner &M, Vector *residual, int64_t &max_iter, Real &tol,const int FromCurrent)
{
	Real resid;
	Vector p, z, q;
	REAL alpha, beta, rho, rho_1 = 0;
	
    REAL normb = TPZExtractVal::val(Norm(b));
	Vector resbackup;
	Vector *res = residual;
	
#ifdef TEST
	std::list< Vector > plist,qlist;
    std::list< Real > betalist, reslist;
	typename std::list<Vector>::iterator jt;
	typename std::list<Vector>::iterator kt;
	Vector Au;
#endif
	
	if(!res) res = &resbackup;
	Vector &r = *res;
	//  Vector r = b - A*x;
	if(FromCurrent) 
    {
        A.MultAdd(x,b,r,-1.,1.);
    }
	else {
		x.Zero();
		r = b;
	}
	
	if (normb == 0.0)
		normb = 1.0;
#ifdef TEST
        reslist.push_back(1./(TPZExtractVal::val(Norm(r))));
#endif

    if ((resid = (TPZExtractVal::val( Norm(r) ) ) / normb) <= tol) {
		tol = resid;
		max_iter = 0;
		return 0;
	}
	int64_t i;
	for (i = 1; i <= max_iter; i++) {
        M.Solve(r,z);
        rho = TPZExtractVal::val(Dot(r, z));
		
		if (i == 1)
			p = z;
		else {
			beta = rho / rho_1;
			p.TimesBetaPlusZ(beta,z);
		}
#ifdef TEST
		plist.push_back(p);
        betalist.push_back(beta);
#endif	 
		
		A.Multiply(p,q);
		alpha = rho / (TPZExtractVal::val(Dot(p, q)));
		
#ifdef TEST
		qlist.push_back(q);
#endif
		
		x.ZAXPY(alpha,p);
		r.ZAXPY(-alpha,q);
		
#ifdef TEST
		A.Multiply(x,Au);
		REAL energy = Dot(x,Au)/2.-Dot(x,b);
#endif
		
		if ((resid = (TPZExtractVal::val(Norm(r))) / normb) <= tol) {
			tol = resid;
			max_iter = i;
#ifdef PZDEBUG
			std::cout << "cg iter = " << i <<  " res = " << resid << std::endl;
#endif

#ifdef TEST    
    reslist.push_back(1./(TPZExtractVal::val(Norm(r))));
    // Real CN=ConditionNumber();
    //Structure to compute the condition number via Lanczos connection. See https://doi.org/10.1371/journal.pone.0130920 for further details
    TPZFMatrix<Real> B(i,i,0.),L(i,i,0.),D(i,i,0.);
    TPZAutoPointer<TPZFMatrix<Real>> T = new TPZFMatrix<Real>;
    T->AutoFill(i,i,true);
    Vector aux(i,1,0.);
    B.Identity();
    D.PutVal(0,0,*std::next(reslist.begin(),0));
    //Assemble the diagonal entries of Lambda matrix
    A.MultAdd(*plist.begin(),*plist.begin(),aux,1.,0.);
    L.PutVal(0,0,Dot(*plist.begin(),aux));
    for (int k = 1; k < i; k++)
    {
        D.PutVal(k,k,*std::next(reslist.begin(),k));
        B.PutVal(k-1,k,*std::next(betalist.begin(),k));
        aux.Zero();
        A.MultAdd(*std::next(plist.begin(),k),*std::next(plist.begin(),k),aux,1.,0.);
        L.PutVal(k,k,Dot(*std::next(plist.begin(),k),aux));
    }

    //Compute the tridiagonal Matrix
    B.Multiply(D,T);
    L.Multiply(T,T);
    B.Transpose();
    B.Multiply(T,T);
    D.Multiply(T,T);

    TPZAutoPointer<TPZFMatrix<Real>> T1,T2,T3;
    T1=T;
    T2=T;
    T3=T;
    
    std::ofstream rprint0,rprint1,rprint2,rprint3,rprint4;
    rprint0.open("condition_number0.txt",std::ios_base::app);
    rprint1.open("condition_number1.txt",std::ios_base::app);
    rprint2.open("condition_number2.txt",std::ios_base::app);
    rprint3.open("All_eigenvalues.txt",std::ios_base::app);
    rprint4.open("MaxMin_eigenvalues.txt",std::ios_base::app);
    rprint0 << T1->ConditionNumber(0) << "\n";
    rprint1 << T2->ConditionNumber(1) << "\n";
    rprint2 << T3->ConditionNumber(2) << "\n";

    TPZLapackEigenSolver<Real> eigSolver;
    
    TPZVec<std::complex<double>> eigenvalues;
    eigSolver.SetMatrixA(T);
    auto a1 = eigSolver.SolveEigenProblem(eigenvalues);

    Real maxEig = 0.;
    Real minEig = 1e3;
    Real minAbs = 1e3;
    Real maxAbs = 0.;
    for (int i = 0; i < eigenvalues.size(); i++)
    {
        rprint3 << eigenvalues[i].real() << std::endl;
        if (eigenvalues[i].real() > maxEig) maxEig = eigenvalues[i].real();
        if (eigenvalues[i].real() < minEig) minEig = eigenvalues[i].real();
        if (fabs(eigenvalues[i].real()) < minAbs) minAbs = fabs(eigenvalues[i].real());
        if (fabs(eigenvalues[i].real()) > maxAbs) maxAbs = fabs(eigenvalues[i].real());
    }
    rprint3 << std::endl;

    rprint4 << maxEig << " " << minEig << " " << maxAbs << " " << minAbs << std::endl;

    
    // std::cout << "Eigenvalues = " << eigenvalues << std::endl;

#endif
			return 0;
		}
#ifdef PZDEBUG
		std::cout << "cg iter = " << i <<  " res = " << resid /*<< " energy " << energy */ << std::endl;
#endif
#ifdef TEST
        
        reslist.push_back(1./(TPZExtractVal::val(Norm(r))));
		std::cout << " energy " << energy << std::endl;
		TPZFMatrix<REAL> inner(plist.size(),plist.size(),0.);
		{
			int64_t j,k;
			for(j=0, jt = plist.begin(); jt != plist.end(); jt++,j++)
			{
				for(k=0, kt = qlist.begin(); kt != qlist.end(); kt++,k++)
				{
					inner(j,k) = Dot((*jt),(*kt));
				}
			}
		}
		// inner.Print("Inner product of search directions");
#endif
		rho_1 = rho;

	}//CG Iterative process
	
	tol = resid;
	std::cout << "cg iter = " << i <<  " res = " << resid << std::endl;
    //END OF CG

	return 1;
}
