#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzartdiff.h"
#include "pzdiffmatrix.h"
#include "pzdiffmatrix.h"
#include "pzeulerconslaw.h"
#include "stdlib.h"
#include "pztempmat.h"
#include "pzvec.h"
#include "fadType.h"

typedef REAL Number;

void error(char * teste)
{

}

int main()
{
 {// Rotation Tests
   TPZDiffMatrix<Number> Rot, RotT, Temp;
   TPZVec<Number> sol;
   Number us;

   cout << "Testing Rotation 2d\n";

   sol.Resize(4);
   sol[0] = 1.;
   sol[1] = 3.3;
   sol[2] = 2.3;
   sol[3] = 5;
   TPZEulerConsLaw2::uRes(sol, us);

   TPZArtDiff::RotMatrix(sol, us, Rot, RotT);

   cout << "Rot" << Rot << endl;
   cout << "RotT" << RotT << endl;
   Rot.Multiply(RotT, Temp);
   cout << "Rot.RotT" << Temp << endl;


   cout << "Testing Rotation 3d\n";

   sol.Resize(5);
   sol[0] = 1.;
   sol[1] = 3.3;
   sol[2] = 2.3;
   sol[3] = 4.5;
   sol[4] = 5;
   TPZEulerConsLaw2::uRes(sol, us);

   TPZArtDiff::RotMatrix(sol, us, Rot, RotT);

   cout << "Rot" << Rot << endl;
   cout << "RotT" << RotT << endl;
   Rot.Multiply(RotT, Temp);
   cout << "Rot.RotT" << Temp << endl;
 }

 {
   TPZDiffMatrix<Number> M, Mi, Temp;
   TPZVec<Number> sol;
   Number us;

   cout << "Testing Matrix M 2d\n";

   sol.Resize(4);
   sol[0] = 1.;
   sol[1] = 3.3;
   sol[2] = 2.3;
   sol[3] = 5;
   TPZEulerConsLaw2::uRes(sol,us);

   TPZArtDiff::MMatrix(sol, us, 1.4, M, Mi);

   cout << "M" << M << endl;
   cout << "Mi" << Mi << endl;
   M.Multiply(Mi, Temp);
   cout << "M.Mi" << Temp << endl;


   cout << "Testing Matrix M 3d\n";

   sol.Resize(5);
   sol[0] = 1.;
   sol[1] = 3.3;
   sol[2] = 2.3;
   sol[3] = 4.5;
   sol[4] = 5;
   TPZEulerConsLaw2::uRes(sol,us);

   TPZArtDiff::MMatrix(sol, us, 1.4, M, Mi);

   cout << "M" << M << endl;
   cout << "Mi" << Mi << endl;
   M.Multiply(Mi, Temp);
   cout << "M.Mi" << Temp << endl;

  }


 {
   TPZDiffMatrix<Number> X, Xi, Lambda, Temp;
   TPZVec<Number> sol;
   Number us, c;

   cout << "Testing Matrix X 2d\n";

   sol.Resize(4);
   sol[0] = 1.;
   sol[1] = 3.3;
   sol[2] = 2.3;
   sol[3] = 10;
   TPZEulerConsLaw2::uRes(sol,us);

   TPZEulerConsLaw2::cSpeed(sol, 1.4, c);

   TPZArtDiff::EigenSystemSUPG(sol, us, c, 1.4, X, Xi, Lambda);

   cout << "X" << X << endl;
   cout << "Xi" << Xi << endl;
   X.Multiply(Xi, Temp);
   cout << "X.Xi" << Temp << endl;


   cout << "Testing Matrix X 3d\n";

   sol.Resize(5);
   sol[0] = 1.;
   sol[1] = 3.3;
   sol[2] = 2.3;
   sol[3] = 4.5;
   sol[4] = 20;
   TPZEulerConsLaw2::uRes(sol,us);

   TPZEulerConsLaw2::cSpeed(sol, 1.4, c);

   TPZArtDiff::EigenSystemSUPG(sol, us, c, 1.4, X, Xi, Lambda);

   cout << "X" << X << endl;
   cout << "Xi" << Xi << endl;
   X.Multiply(Xi, Temp);
   cout << "X.Xi" << Temp << endl;

  }


 {
   TPZDiffMatrix<Number> Y, Yi, Lambda, Temp;
   TPZVec<Number> sol, aaS;
   Number us, c;

   cout << "Testing Matrix X 2d\n";

   sol.Resize(4);
   sol[0] = 1.;
   sol[1] = 3.3;
   sol[2] = 2.3;
   sol[3] = 10;

   aaS.Resize(2);
   aaS[0] = 2.1;
   aaS[1] = 2.3;

   TPZEulerConsLaw2::uRes(sol,us);

   TPZEulerConsLaw2::cSpeed(sol, 1.4, c);

   TPZArtDiff::EigenSystemBornhaus(sol, us, c, 1.4, aaS, Y, Yi, Lambda);

   cout << "Y" << Y << endl;
   cout << "Yi" << Yi << endl;
   Y.Multiply(Yi, Temp);
   cout << "Y.Yi" << Temp << endl;

   cout << "Testing Matrix X 3d\n";

   sol.Resize(5);
   sol[0] = 1.;
   sol[1] = 3.3;
   sol[2] = 2.3;
   sol[3] = 4.5;
   sol[4] = 20;

   aaS.Resize(3);
   aaS[0] = 2.1;
   aaS[1] = 2.3;
   aaS[2] = 1.7;

   TPZEulerConsLaw2::uRes(sol,us);

   TPZEulerConsLaw2::cSpeed(sol, 1.4, c);

   TPZArtDiff::EigenSystemSUPG(sol, us, c, 1.4, Y, Yi, Lambda);

   cout << "Y" << Y << endl;
   cout << "Yi" << Yi << endl;
   Y.Multiply(Yi, Temp);
   cout << "Y.Yi" << Temp << endl;

  }
}

