#include "pzfmatrix.h"
//#include ".h"
#include "pzsolve.h"
void Orthogonalization_CGS(TPZFMatrix<REAL> &a, TPZFMatrix<REAL> &q,TPZFMatrix<REAL> &r);
int main() {

	TPZFMatrix<REAL> a(3,3,0.);
   TPZFMatrix<REAL> r(3,3,0.);
   TPZFMatrix<REAL> q(3,3,0.);
   a(0,0) = 1.;
   a(0,2) = 3.;
   a(1,0) = 2.;
   a(1,1) = 1.;
   a(2,0) = 1.;
   a(2,2) = 2.;
   Orthogonalization_CGS(a,q,r);
   q.Print();
   REAL resultstore[8];
   TPZFMatrix<REAL> qt(3,3,0.),result(0,0,resultstore,8);
   q.Multiply(q,result,1);
   //q.Transpose(&qt);
   result.Print();
/*   TPZFMatria(0,0) = 1;x c(5,5,1.0);
   TPZFMatrix<REAL> d(5,5);
      for(int i=0; i<b.Rows();i++)
    for(int j=0; j<b.Cols(); j++)
       b(i,j)=(i+1)*(j+2);

   TPZMatrixSolver f(&b);
   f.SetDirect(2);
   f.Solve(e,e,&h);
   e.Print();
     b.Transpose(&d);

     c=b*a;
     c=c+d;
     c.Print();
  // a.Resize(6,6);
  //  a.Redim(3,3);
   a.Print();  */
   return 0;
}


/**Orthogonalization process of Gram-Schmidt*/
void Orthogonalization_CGS(TPZFMatrix<REAL> &a, TPZFMatrix<REAL> &q,TPZFMatrix<REAL> &r) {

      if(a.Rows()!=q.Rows() && a.Cols()!=q.Cols()) {
      	std::cout << "TPZMatrixSolver::Orthogonalization_CGS isn't compatible dimension\n";
      }

      int rows = a.Rows();
      TPZFMatrix<REAL> sk(rows,1,0.);
      for(int k=0;k<rows;k++) {
         REAL sum;

         for(int i=0;i<k;i++) {
         	sum = 0.0;
            for(int l=0;l<rows;l++) sum += q(l,i)*a(l,k);
         	sk(i,0) = sum;
         }
//         int cols = q.Cols();
         TPZFMatrix<REAL> qq(rows,1,0.),zk(rows,1,0.);
         for(int i=0;i<k;i++) {
         	for(int l=0;l<rows;l++) qq(l,0) = sk(i,0)*q(l,i);
            zk += qq;
         }

         if(k==0) {
         	for(int i=0;i<rows;i++) zk(i,0) = a(i,0);
         } else {
         	for(int i=0;i<rows;i++) zk(i,0) = a(i,k) - zk(i,0);
         }
         TPZFMatrix<REAL> zaux(1,rows);
         zk.Transpose(&zaux);
         TPZFMatrix<REAL> rkk(1,1);
         rkk = zaux*zk;
         REAL Rkk = sqrt(rkk(0,0));
         REAL div = 1./Rkk;
         for(int l=0;l<rows;l++) q(l,k) = zk(l,0)*div;
         for(int i=0;i<k;i++) r(i,0) = sk(i,0)*div;
      }
}


