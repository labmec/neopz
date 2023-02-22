#include "pzreal.h"
#include "TPZSimpleTimer.h"
#include "tpzautopointer.h"
#include "pzvec.h"


template<class TVar>
void GeneratePlaneRotation(const TVar dx, const TVar dy, TVar &cs, TVar &sn);

template<class TVar>
void ApplyPlaneRotation(TVar &dx, TVar &dy, const TVar cs, const TVar sn);

template<class TVar>
void Update(TPZFMatrix<TVar>&x, const int k, const TPZFMatrix<TVar>&H,
            const TPZVec<TVar> &s, const TPZVec<TPZFMatrix<TVar>>&v);

template<class TVar>
int GMRES(const TPZMatrix<TVar> &A, TPZFMatrix<TVar> &x, const TPZFMatrix<TVar> &b,
					TPZMatrixSolver<TVar> &M, TPZFMatrix<TVar> &H, const int krylovdim,
					int64_t &max_iter, RTVar &tol,
					TPZFMatrix<TVar> *residual, const int fromcurrent){  
  TPZSimpleTimer gmres("GMRES");
  //allocating structures
  TPZVec<TPZFMatrix<TVar>> v(krylovdim+1);
  TPZVec<TVar> cs(krylovdim+1), sn(krylovdim+1), s(krylovdim+1);

  //compute rhs norm
	TPZFMatrix<TVar> r;
	M.Solve(b,r);
	const RTVar normb = [&r]{
		RTVar norm_b = Norm(r);
		if(IsZero(norm_b)) return (RTVar)1;
		return norm_b;
	}();
	
  if(fromcurrent){//we actually need to compute the residual
		//r = b-A*x without dynamic allocation
		A.MultAdd(x,b,r,-1,1);
    M.Solve(r,r);
  }else{//residual is equal to rhs for initial sol == 0   
    x.Zero();
  }
	
  RTVar beta = Norm(r);
  RTVar resid = beta;
  if ((resid / normb) <= tol) {
    tol = resid;
    max_iter = 0;
    return 0;
  }

	//we avoid allocating w at every step
	TPZFMatrix<TVar> w(r.Rows(),1,0), w1;
  int iter = 1;
  while (iter <= max_iter) {
    v[0] = r * (1.0 / beta);
    s = 0.0;
    s[0] = beta;

    int i = 0;
    for (; i < krylovdim && iter <= max_iter; i++, iter++) {
			A.Multiply(v[i],w1);
			M.Solve(w1,w);
      for (int k = 0; k <= i; k++) {
        const auto dotwv = Dot(w,v[k]);
        H(k, i) = dotwv;
        w -= dotwv * v[k];
      }
      const auto normw = Norm(w);
      if(IsZero(normw)){
        DebugStop();
      }
      H(i+1, i) = normw;
      w *= ((TVar)1.0 / H(i+1, i));
      v[i+1] = w;

      for (int k = 0; k < i; k++){
        ApplyPlaneRotation(H(k,i), H(k+1,i), cs[k], sn[k]);
      }
      
      GeneratePlaneRotation(H(i,i), H(i+1,i), cs[i], sn[i]);
      ApplyPlaneRotation(H(i,i), H(i+1,i), cs[i], sn[i]);
      ApplyPlaneRotation(s[i], s[i+1], cs[i], sn[i]);
      
      if ((resid = std::abs(s[i+1]) / normb) < tol) {
				i++;
        break;
      }
      std::cout << iter << "\t" << std::scientific << resid << std::endl;
			// if(iter % 5 == 0){
			// }
    }
    Update(x, i - 1, H, s, v);
		//r = b-A*x without dynamic allocation
		A.MultAdd(x,b,r,-1,1);
    M.Solve(r,r);
    beta = Norm(r);
    std::cout<<"beta "<<beta<<" resid "<<resid<<" normb "<<normb<<std::endl;
    if ((resid = beta / normb) < tol) {
      tol = resid;
      max_iter = iter;
      return 0;
    }
  }
  
  tol = resid;
  return 1;
}

/** @ingroup util */
/** @brief Compute the values cs and sn parameters to rotation
    @note These rotations are based on lapack Xlartg routines. (X=s,d,c,z)
    They will result in a plane rotation with real cosine and complex sine*/
template<class TVar>
void GeneratePlaneRotation(TVar dx, TVar dy, TVar &cs, TVar &sn)
{
  
  if(IsZero(dy)){
    cs = 1.;
    sn = 0;
  }else if (IsZero(dx)){
    cs = 0;
    if constexpr(is_complex<TVar>::value){
      sn = std::conj(dy)/abs(dy);
    }else{
      sn = dy/std::abs(dy);
    }
  }else{
    //note: std::norm is actually the squared magnitude
    const TVar norm = (dx / abs(dx) ) * sqrt(std::norm(dx) + std::norm(dy));
    sn = dy/norm;
    cs = dx/norm;
  }
  
}

/** @ingroup util */
/** @brief Makes rotation of the plane based on the cs and sn parameters
    @note These rotations are based on lapack Xlartg routines. (X=s,d,c,z)
    They will result in a plane rotation with real cosine and complex sine*/
template<class TVar>
void ApplyPlaneRotation(TVar &dx, TVar &dy, const TVar cs, const TVar sn)
{
  const TVar temp = cs * dx + sn * dy;
  if constexpr (is_complex<TVar>::value){
    dy = -std::conj(sn) * dx + cs * dy;
  }else{
    dy = -sn * dx + cs * dy;
  }
  dx = temp;
}

template<class TVar>
void Update(TPZFMatrix<TVar>&x, const int k, const TPZFMatrix<TVar>&H,
            const TPZVec<TVar> &s, const TPZVec<TPZFMatrix<TVar>>&v)
{
  TPZVec<TVar> y(s);

  // Backsolve:  
  for (int i = k; i >= 0; i--) {
    y[i] /= H(i,i);
    for (int j = i - 1; j >= 0; j--)
      y[j] -= H(j,i) * y[i];
  }

  for (int j = 0; j <= k; j++)
    x += v[j] * y[j];
}