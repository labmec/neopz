/**
 * @file
 * @brief Contains the implementation of the CHEBY function which solves the symmetric positive definite linear system using the 
 * Preconditioned Chebyshec method. 
 */

/**
 * @ingroup solver
 * @brief CHEBY solves the symmetric positive definite linear system \f$ Ax = b \f$ using the Preconditioned Chebyshev Method
 * @return The return value indicates convergence within max_iter (input) iterations (0), or no convergence within max_iter iterations (1). \n
 * Upon successful return, output arguments have the following values:
 * @param A  -- matrix of the system
 * @param b  -- vector of the system
 * @param M  -- preconditioner matrix
 * @param x  --  approximate solution to Ax = b
 * @param max_iter  --  the number of iterations performed before the tolerance was reached
 * @param tol  --  the residual after the final iteration
 * @param eigmin  -- minimun eigenvalue type
 * @param eigmax  -- maximun eigenvalue type
*/
/**
 * Iterative template routine -- CHEBY \n
 * CHEBY follows the algorithm described on p. 30 of the SIAM Templates book.
*/

template < class Matrix, class Vector, class Preconditioner, class Real,
class Type >
int 
CHEBY(const Matrix &A, Vector &x, const Vector &b,
      const Preconditioner &M, int64_t &max_iter, Real &tol,
      Type eigmin, Type eigmax)
{
	Real resid;
	Type alpha, beta, c, d;
	Vector p, q, z;
	
	Real normb = norm(b);
	Vector r = b - A * x;
	
	if (normb == 0.0)
		normb = 1;
	
	if ((resid = norm(r) / normb) <= tol) {
		tol = resid;
		max_iter = 0;
		return 0;
	}
	
	c = (eigmax - eigmin) / 2.0;
	d = (eigmax + eigmin) / 2.0;
	
	for (int64_t i = 1; i <= max_iter; i++) {
		z = M.solve(r);                 // apply preconditioner
		
		if (i == 1) {
			p = z;
			alpha = 2.0 / d;
		} else {
			beta = c * alpha / 2.0;       // calculate new beta
			beta = beta * beta;
			alpha = 1.0 / (d - beta);     // calculate new alpha
			p = z + beta * p;             // update search direction
		}
		
		q = A * p;
		x += alpha * p;                 // update approximation vector
		r -= alpha * q;                 // compute residual
		
		if ((resid = norm(r) / normb) <= tol) {
			tol = resid;
			max_iter = i;
			return 0;                     // convergence
		}
	}
	
	tol = resid;
	return 1;                         // no convergence
}
