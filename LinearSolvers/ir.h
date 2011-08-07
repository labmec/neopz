/**
 * @file
 * @brief Contains the implementation of the IR function which solves the unsymmetric linear system using the Iterative Refinement method. 
 */
/**
 * @ingroup solvers
 * @brief IR solves the unsymmetric linear system Ax = b using Iterative Refinement (preconditioned Richardson iteration).
 * @return The return value indicates convergence within max_iter (input) iterations (0), or no convergence within max_iter iterations (1).\n
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
 * Iterative template routine -- Preconditioned Richardson
 */
template < class Matrix, class Vector, class Preconditioner, class Real >
int 
IR( Matrix &A, Vector &x,const Vector &b,
   Preconditioner &M, Vector *residual, int &max_iter, Real &tol,const int FromCurrent)
{
	Real resid;
	Vector z;
	
	Real normb = Norm(b);
	Vector resbackup;
	Vector *res = residual;
	if(!res) res = &resbackup;
	Vector &r = *res;
	//  Vector r = b - A*x;
	if(FromCurrent) A.MultAdd(x,b,r,-1.,1.);
	else {
		x.Zero();
		r = b;
	}
	//  Vector r = b - A*x;
	
	if (normb == 0.0)
		normb = 1;
	
	if ((resid = Norm(r) / normb) <= tol) {
		tol = resid;
		max_iter = 0;
		return 0;
	}
	
	for (int i = 1; i <= max_iter; i++) {
		//	 z = M.solve(r);
		M.Solve(r,z);
		x += z;
		A.MultAdd(x,b,r,-1.,1.);
		//	 r = b - A * x;
		
		if ((resid = Norm(r) / normb) <= tol) {
			tol = resid;
			max_iter = i;
			return 0;
		}
	}
	
	tol = resid;
	return 1;
}



