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
#ifdef TEST
#include <list> 
#endif

template<class TVar>
int
CG( TPZMatrix<TVar> &A, TPZFMatrix<TVar> &x, const TPZFMatrix<TVar> &b,
		TPZMatrixSolver<TVar> &M, TPZFMatrix<TVar> *residual, int64_t &max_iter, RTVar &tol,const int FromCurrent)
{
	RTVar resid;
	TPZFMatrix<TVar> p, z, q;
	TVar alpha, beta, rho, rho_1 = 0;
	
	RTVar normb = Norm(b);
	TPZFMatrix<TVar> resbackup;
	TPZFMatrix<TVar> *res = residual;
	
#ifdef TEST
	std::list< TPZFMatrix<TVar> > plist,qlist;
	std::list< TPZFMatrix<TVar> >::iterator jt;
	std::list< TPZFMatrix<TVar> >::iterator kt;
	TPZFMatrix<TVar> Au;
#endif
	
	if(!res) res = &resbackup;
	TPZFMatrix<TVar> &r = *res;
	//  TPZFMatrix<TVar> r = b - A*x;
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
	
	if ((resid = Norm(r) / normb) <= tol) {
		tol = resid;
		max_iter = 0;
		return 0;
	}
	int64_t iter;
	for (iter = 1; iter <= max_iter; iter++) {
		M.Solve(r,z);
		rho = Dot(r, z);
		
		if (iter == 1)
			p = z;
		else {
			beta = rho / rho_1;
			p.TimesBetaPlusZ(beta,z);
		}
#ifdef TEST
		plist.push_back(p);
#endif	 
		
		A.Multiply(p,q);
		alpha = rho / (Dot(p, q));
		
#ifdef TEST
		qlist.push_back(q);
#endif
		
		x.ZAXPY(alpha,p);
		r.ZAXPY(-alpha,q);
		
#ifdef TEST
		A.Multiply(x,Au);
		TVar energy = Dot(x,Au)/2.-Dot(x,b);
#endif
		
		if ((resid = (Norm(r)) / normb) <= tol) {
			tol = resid;
			max_iter = iter;
			std::cout << iter << "\t" << resid << std::endl;
			return 0;
		}
		std::cout << iter << "\t" << resid << std::endl;
#ifdef TEST
		std::cout << " energy " << energy << std::endl;
		TPZFMatrix<TVar> inner(plist.size(),plist.size(),0.);
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
		inner.Print("Inner product of search directions");
#endif
		rho_1 = rho;
	}
	
	tol = resid;
	return 1;
}
