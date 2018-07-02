/**
 * @file
 * @brief Contains the implementation of the CMRES function which solves the unsymmetric linear system using the Generalized Minimum Residual method. 
 */

/** @ingroup util */
/** @brief Compute the values cs and sn parameters to rotation */
template<class Real>
void GeneratePlaneRotation(Real &dx, Real &dy, Real &cs, Real &sn);

/** @ingroup util */
/** @brief Makes rotation of the plane based on the cs and sn parameters */
template<class Real>
void ApplyPlaneRotation(Real &dx, Real &dy, Real &cs, Real &sn);

/**
 * @ingroup solver
 * @brief Computes back solve. Updates on v vector.
 */
template < class Matrix, class Vector >
void 
Update(Vector &x, int64_t k, Matrix &h, Vector &s, Vector v[])
{
	Vector y(s);
	
	// Backsolve:  
	for (int64_t i = k; i >= 0; i--) {
		y(i) /= h(i,i);
		for (int64_t j = i - 1; j >= 0; j--)
			y(j) -= h(j,i) * y(i);
	}
	
	for (int64_t j = 0; j <= k; j++)
		x.ZAXPY(y(j),v[j]);
}

/**
 * @ingroup solver
 * @brief Returns absolute value to real.
 */
template < class Real >
Real
abs(Real x)
{
	return (x > ((Real)0.) ? x : -x);
}

/**
 * @ingroup solver
 * @brief GMRES solves the unsymmetric linear system Ax = b using the Generalized Minimum Residual method
 * @return The return value indicates convergence within max_iter (input) iterations (0), or no convergence within max_iter iterations (1).\n
 * Upon successful return, output arguments have the following values:
 * @param A Matrix of the system
 * @param b Vector of the system
 * @param M Preconditioner matrix
 * @param H Auxiliar matrix (?)
 * @param m Size (?)
 * @param x Approximate solution to \f$ Ax = b \f$
 * @param max_iter The number of iterations performed before the tolerance was reached
 * @param tol The residual after the final iteration
 * @param residual Residual vector (return)
 * @param FromCurrent For type of operation (MultAdd)
 */
/**
 * Iterative template routine -- GMRES \n
 * GMRES follows the algorithm described on p. 20 of the SIAM Templates book.
 */
template < class Operator, class Vector, class Preconditioner,
class Matrix, class Real >
int 
GMRES( Operator &A, Vector &x, const Vector &b,
	  Preconditioner &M, Matrix &H, int &m, int64_t &max_iter,
	  Real &tol, Vector *residual,const int FromCurrent)
{
	Real resid;
	int64_t j = 1, k;
	Vector s(m+1), cs(m+1), sn(m+1), w1,w;
	
	//  Real normb = norm(M.Solve(b));
	Vector resbackup;
	Vector *res = residual;
	if(!res) res = &resbackup;
	Vector &r = *res;
    M.Solve(b,r);
	Real normb = TPZExtractVal::val(Norm(r));	//  Vector r = b - A*x;
	if(FromCurrent) 
    {
        A.MultAdd(x,b,r,-1.,1.);
    } 
    else 
    {
		x.Zero();
		r = b;
	}
	M.Solve(r,w);
	r=w;
	Real beta = TPZExtractVal::val(Norm(r));
	
	if (normb == 0.0)
		normb = 1;
	
	if ((resid = ((Real)TPZExtractVal::val(Norm(r))) / normb) <= tol) {
		tol = resid;
		x+=r;
		max_iter = 0;
		return 0;
	}
	
	Vector *v = new Vector[m+1];
	
	while (j <= max_iter) {
        int64_t rows = r.Rows();
        v[0].Resize(rows,1);
        for (int64_t i=0; i<rows; i++) {
            v[0](i) = TPZExtractVal::val(r(i)) * ((Real)(1.0/beta));
        }
        int64_t srows = s.Rows();
        for (int64_t i=0; i<srows; i++) {
            s(i) = REAL(0.);
        }
//		v[0] = r * (REAL(1.0 / beta));    // ??? r / beta
//		s = REAL(0.0);
		s(0) = beta;
		
		for (int64_t i = 0; i < m && j <= max_iter; i++, j++) {
			A.Multiply(v[i],w1);
			M.Solve(w1,w);
			for (k = 0; k <= i; k++) {
				H(k, i) = Dot(w, v[k]);
				w.ZAXPY(-H(k, i), v[k]);
			}
			H(i+1, i) = Norm(w);
			v[i+1] = w;
			v[i+1] *= (1.0)/TPZExtractVal::val(H(i+1,i));
			
			for (k = 0; k < i; k++)
				ApplyPlaneRotation(H(k,i), H(k+1,i), cs(k), sn(k));
			
			GeneratePlaneRotation(H(i,i), H(i+1,i), cs(i), sn(i));
			ApplyPlaneRotation(H(i,i), H(i+1,i), cs(i), sn(i));
			ApplyPlaneRotation(s(i), s(i+1), cs(i), sn(i));
			resid = ((Real)fabs(s(i+1)))/normb;
			if (resid < tol) {
				Update(x, i, H, s, v);
				tol = resid;
				max_iter = j;
				delete [] v;
				return 0;
			}
		}
		Update(x, m - 1, H, s, v);
		A.MultAdd(x,b,r,-1.,1.);
		M.Solve(r,r);
		beta = TPZExtractVal::val(Norm(r));
        resid = beta/normb;
		if (resid < tol) {
            std::cout << "iter " << j << " - " << resid << std::endl;
			tol = resid;
			max_iter = j;
			delete [] v;
			return 0;
		}
	}
	
	tol = resid;
	delete [] v;
	return 1;
}


#include <math.h>

/** @ingroup util */
/** @brief Compute the values cs and sn parameters to rotation */
template<class Real>
void GeneratePlaneRotation(Real &dx, Real &dy, Real &cs, Real &sn)
{
	if (dy == ((Real)0.0)) {
		cs = 1.0;
		sn = 0.0;
	} else if (abs(dy) > abs(dx)) {
		Real temp = dx / dy;
		sn = ((Real)1.0) / sqrt( ((Real)1.0) + temp*temp );
		cs = temp * sn;
	} else {
		Real temp = dy / dx;
		cs = ((Real)1.0) / sqrt( ((Real)1.0) + temp*temp );
		sn = temp * cs;
	}
}

/** @ingroup util */
/** @brief Makes rotation of the plane based on the cs and sn parameters */
template<class Real>
void ApplyPlaneRotation(Real &dx, Real &dy, Real &cs, Real &sn)
{
	Real temp  =  cs * dx + sn * dy;
	dy = -sn * dx + cs * dy;
	dx = temp;
}

