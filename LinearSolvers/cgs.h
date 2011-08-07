/**
 * @file
 * @brief Contains the implementation of the CGS function which solves the unsymmetric linear system using the Conjudate Gradient Squared method. 
 */
/**
 * @ingroup solver
 * @brief CGS solves the unsymmetric linear system \f$ Ax = b \f$ using the Conjugate Gradient Squared method
 * @return The return value indicates convergence within max_iter (input) iterations (0), or no convergence within max_iter iterations (1). \n
 * Upon successful return, output arguments have the following values:
 * @param A  -- matrix of the system
 * @param b  -- vector of the system
 * @param M  -- preconditioner matrix
 * @param x  --  approximate solution to Ax = b
 * @param max_iter  --  the number of iterations performed before the tolerance was reached
 * @param tol  --  the residual after the final iteration
 */
/**
 * Iterative template routine -- CGS \n
 * CGS follows the algorithm described on p. 26 of the SIAM Templates book.
 */
template < class Matrix, class Vector, class Preconditioner, class Real >
int 
CGS(const Matrix &A, Vector &x, const Vector &b,
    const Preconditioner &M, int &max_iter, Real &tol)
{
	Real resid;
	Vector rho_1(1), rho_2(1), alpha(1), beta(1);
	Vector p, phat, q, qhat, vhat, u, uhat;
	
	Real normb = norm(b);
	Vector r = b - A*x;
	Vector rtilde = r;
	
	if (normb == 0.0)
		normb = 1;
	
	if ((resid = norm(r) / normb) <= tol) {
		tol = resid;
		max_iter = 0;
		return 0;
	}
	
	for (int i = 1; i <= max_iter; i++) {
		rho_1(0) = dot(rtilde, r);
		if (rho_1(0) == 0) {
			tol = norm(r) / normb;
			return 2;
		}
		if (i == 1) {
			u = r;
			p = u;
		} else {
			beta(0) = rho_1(0) / rho_2(0);
			u = r + beta(0) * q;
			p = u + beta(0) * (q + beta(0) * p);
		}
		phat = M.solve(p);
		vhat = A*phat;
		alpha(0) = rho_1(0) / dot(rtilde, vhat);
		q = u - alpha(0) * vhat;
		uhat = M.solve(u + q);
		x += alpha(0) * uhat;
		qhat = A * uhat;
		r -= alpha(0) * qhat;
		rho_2(0) = rho_1(0);
		if ((resid = norm(r) / normb) < tol) {
			tol = resid;
			max_iter = i;
			return 0;
		}
	}
	
	tol = resid;
	return 1;
}

