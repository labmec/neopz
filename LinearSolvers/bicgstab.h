/**
 * @file
 * @brief Contains the implementation of the BiCGSTAB function which solves the unsymmetric linear system \n
 * using BiConjugate Gradient Stabilized method.
 */
/**
 * @ingroup solvers
 * @brief BiCGSTAB solves the unsymmetric linear system \f$ Ax = b \f$ using the Preconditioned BiConjugate Gradient Stabilized method
 * @return The return value indicates convergence within max_iter (input) iterations (0), or no convergence within max_iter iterations (1). \n
 * Upon successful return, output arguments have the following values:
 * @param A  -- matrix of the system
 * @param b  -- vector of the system
 * @param M  -- preconditioner matrix
 * @param x  -- approximate solution to \f$ Ax = b \f$
 * @param max_iter  --  the number of iterations performed before the tolerance was reached
 * @param tol  --  the residual after the final iteration  
 * @param residual  -- residual vector (return)
 * @param FromCurrent  -- for type of operation (MultAdd)
 */
/**
 * Iterative template routine -- BiCGSTAB
 * BiCGSTAB follows the algorithm described on p. 27 of the SIAM Templates book.
 */
template < class Matrix, class Vector, class Preconditioner, class Real >
int 
BiCGSTAB(const Matrix &A, Vector &x, const Vector &b,
         /*const */Preconditioner &M, int &max_iter, Real &tol, Vector *residual,const int FromCurrent)
{
	Real resid;
	Vector rho_1(1), rho_2(1), alpha(1), beta(1), omega(1);
	Vector p, phat, s, shat, t, v;
	
	Vector resbackup;
	Vector *res = residual;
	if(!res) res = &resbackup;
	Vector &r = *res;
	
	Real normb = Norm(b);
	if(FromCurrent) A.MultAdd(x,b,r,-1.,1.);
	else {
		x.Zero();
		r = b;
	}
	
	Vector rtilde = r;
	
	if (normb == 0.0)
		normb = 1;
	
	//  if ((resid = norm(r) / normb) <= tol) {
	if ((resid = Norm(r) / normb) <= tol) {
		tol = resid;
		max_iter = 0;
		return 0;
	}
	
	for (int i = 1; i <= max_iter; i++) {
		//    rho_1(0) = dot(rtilde, r);
		rho_1(0) = Dot(rtilde, r);
		if (rho_1(0) == 0) {
			//      tol = norm(r) / normb;
			tol = Norm(r) / normb;
			return 2;
		}
		if (i == 1)
			p = r;
		else {
			beta(0) = (rho_1(0)/rho_2(0)) * (alpha(0)/omega(0));
			p = r + beta(0) * (p - omega(0) * v);
		}
		
		//    phat = M.solve(p);
		M.Solve(p, phat);
		//    v = A * phat;
		A.Multiply(phat, v);
		//    alpha(0) = rho_1(0) / dot(rtilde, v);
		alpha(0) = rho_1(0) / Dot(rtilde, v);
		s = r - alpha(0) * v;
		//    if ((resid = norm(s)/normb) < tol) {
		if ((resid = Norm(s)/normb) < tol) {
			x += alpha(0) * phat;
			tol = resid;
			return 0;
		}
		//    shat = M.solve(s);
		M.Solve(s, shat);
		//    t = A * shat;
		A.Multiply(shat, t);
		//    omega = dot(t,s) / dot(t,t);
		omega = Dot(t,s) / Dot(t,t);
		x += alpha(0) * phat + omega(0) * shat;
		r = s - omega(0) * t;
		
		rho_2(0) = rho_1(0);
		//    if ((resid = norm(r) / normb) < tol) {
		if ((resid = Norm(r) / normb) < tol) {
			tol = resid;
			max_iter = i;
			return 0;
		}
		if (omega(0) == 0) {
			//      tol = norm(r) / normb;
			tol = Norm(r) / normb;
			return 3;
		}
	}
	
	tol = resid;
	return 1;
}
