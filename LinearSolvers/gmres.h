//*****************************************************************
// Iterative template routine -- GMRES
//
// GMRES solves the unsymmetric linear system Ax = b using the 
// Generalized Minimum Residual method
//
// GMRES follows the algorithm described on p. 20 of the 
// SIAM Templates book.
//
// The return value indicates convergence within max_iter (input)
// iterations (0), or no convergence within max_iter iterations (1).
//
// Upon successful return, output arguments have the following values:
//  
//        x  --  approximate solution to Ax = b
// max_iter  --  the number of iterations performed before the
//               tolerance was reached
//      tol  --  the residual after the final iteration
//  
//*****************************************************************


template < class Matrix, class Vector >
void 
Update(Vector &x, int k, Matrix &h, Vector &s, Vector v[])
{
  Vector y(s);

  // Backsolve:  
  for (int i = k; i >= 0; i--) {
    y(i) /= h(i,i);
    for (int j = i - 1; j >= 0; j--)
      y(j) -= h(j,i) * y(i);
  }

  for (int j = 0; j <= k; j++)
		x.ZAXPY(y(j),v[j]);
//    x += v[j] * y(j);
}


template < class Real >
Real
abs(Real x)
{
  return (x > 0 ? x : -x);
}


template < class Operator, class Vector, class Preconditioner,
           class Matrix, class Real >
int 
GMRES(const Operator &A, Vector &x, const Vector &b,
		Preconditioner &M, Matrix &H, int &m, int &max_iter,
		Real &tol, Vector *residual,const int FromCurrent)
{
	Real resid;
	int i, j = 1, k;
	Vector s(m+1), cs(m+1), sn(m+1), w1,w;

//  Real normb = norm(M.Solve(b));
	Vector resbackup;
	Vector *res = residual;
	if(!res) res = &resbackup;
	Vector &r = *res;
	M.Solve(b,r);
	Real normb = Norm(r);
//  Vector r = b - A*x;
	if(FromCurrent) A.MultAdd(x,b,r,-1.,1.);
	else {
		x.Zero();
		r = b;
	}
//	Vector r = M.solve(b - A * x);
	M.Solve(r,w);
	r=w;
	Real beta = Norm(r);

  if (normb == 0.0)
	 normb = 1;

  if ((resid = Norm(r) / normb) <= tol) {
	 tol = resid;
	 x+=r;
	 max_iter = 0;
	 return 0;
  }

  Vector *v = new Vector[m+1];

  while (j <= max_iter) {
	 v[0] = r * (REAL(1.0) / beta);    // ??? r / beta
	 s = REAL(0.0);
	 s(0) = beta;

	 for (i = 0; i < m && j <= max_iter; i++, j++) {
//		w = M.solve(A * v[i]);
		A.Multiply(v[i],w1);
		M.Solve(w1,w);
		for (k = 0; k <= i; k++) {
		  H(k, i) = Dot(w, v[k]);
//		  w -= H(k, i) * v[k];
		  w.ZAXPY(-H(k, i), v[k]);
		}
		H(i+1, i) = Norm(w);
		v[i+1] = w;
		v[i+1] *= REAL(1.0)/H(i+1,i);
//		v[i+1] = w * (1.0 / H(i+1, i)); // ??? w / H(i+1, i)

		for (k = 0; k < i; k++)
		  ApplyPlaneRotation(H(k,i), H(k+1,i), cs(k), sn(k));

		GeneratePlaneRotation(H(i,i), H(i+1,i), cs(i), sn(i));
		ApplyPlaneRotation(H(i,i), H(i+1,i), cs(i), sn(i));
		ApplyPlaneRotation(s(i), s(i+1), cs(i), sn(i));

		if ((resid = abs(s(i+1)) / normb) < tol) {
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
//	 r = M.solve(b - A * x);
	 beta = Norm(r);
	 if ((resid = beta / normb) < tol) {
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


template<class Real>
void GeneratePlaneRotation(Real &dx, Real &dy, Real &cs, Real &sn)
{
  if (dy == 0.0) {
    cs = 1.0;
    sn = 0.0;
  } else if (abs(dy) > abs(dx)) {
	 Real temp = dx / dy;
	 sn = REAL(1.0) / sqrt( REAL(1.0) + temp*temp );
	 cs = temp * sn;
  } else {
	 Real temp = dy / dx;
	 cs = REAL(1.0) / sqrt( REAL(1.0) + temp*temp );
	 sn = temp * cs;
  }
}


template<class Real>
void ApplyPlaneRotation(Real &dx, Real &dy, Real &cs, Real &sn)
{
  Real temp  =  cs * dx + sn * dy;
  dy = -sn * dx + cs * dy;
  dx = temp;
}

