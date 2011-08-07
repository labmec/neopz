/**
 * @file
 * @brief Contains the implementation of the BiCG function which solves the unsymmetric linear system \n
 * using Preconditioned BiConjugate Gradient method.
 */
/******************************************************************
 * @ingroup solver
 * @brief BiCG solves the unsymmetric linear system \f$ Ax = b \f$
 * using the Preconditioned BiConjugate Gradient method
 * @return The return value indicates convergence within max_iter (input)
 * iterations (0), or no convergence within max_iter iterations (1). \n
 * Upon successful return, output arguments have the following values:
 * @param A  -- matrix of the system
 * @param b  -- vector of the system
 * @param M  -- preconditioner matrix
 * @param x  --  approximate solution to Ax = b
 * @param max_iter  --  the number of iterations performed before the tolerance was reached
 * @param tol  --  the residual after the final iteration
 *****************************************************************/
/**
 * Iterative template routine -- BiCG \n
 * BiCG follows the algorithm described on p. 22 of the SIAM Templates book.
 */
template < class Matrix, class Vector, class Preconditioner, class Real >
int 
BiCG( Matrix &A, Vector &x, const Vector &b,
     Preconditioner &M, int &max_iter, Real &tol)
{
	Real resid;
	Vector rho_1(1), rho_2(1), alpha(1), beta(1);
	Vector z, ztilde, p, ptilde, q, qtilde;
	
	Real normb = Norm(b);
	Vector r = b - A * x;
	Vector rtilde = r;
	
	if (normb == 0.0)
		normb = 1;
	
	if ((resid = Norm(r) / normb) <= tol) {
		tol = resid;
		max_iter = 0;
		return 0;
	}
	
	for (int i = 1; i <= max_iter; i++) {
		M.Solve(r,z);
		//    ztilde = M.trans_solve(rtilde);
		//    ztilde = M.Solve(rtilde);
		M.Solve(rtilde,ztilde);
		rho_1(0) = Dot(z, rtilde);
		if (rho_1(0) == 0) { 
			tol = Norm(r) / normb;
			max_iter = i;
			return 2;
		}
		if (i == 1) {
			p = z;
			ptilde = ztilde;
		} else {
			beta(0) = rho_1(0) / rho_2(0);
			//      p = z + beta(0) * p;
			p.TimesBetaPlusZ(beta(0),z);
			//      ptilde = ztilde + beta(0) * ptilde;
			ptilde.TimesBetaPlusZ(beta(0),ztilde);
		}
		q = A * p;
		A.Multiply(ptilde,qtilde,1);
		//    qtilde = A.trans_mult(ptilde);
		alpha(0) = rho_1(0) / Dot(ptilde, q);
		x += alpha(0) * p;
		r -= alpha(0) * q;
		rtilde -= alpha(0) * qtilde;
		
		rho_2(0) = rho_1(0);
		if ((resid = Norm(r) / normb) < tol) {
			tol = resid;
			max_iter = i;
			return 0;
		}
	}
	
	tol = resid;
	return 1;
}

