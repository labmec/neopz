#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzartdiff.h"
#include "pzdiffmatrix.h"
#include "pzdiffmatrix.h"
#include "pzeulerconslaw.h"
#include "stdlib.h"
//#include "pztempmat.h"
#include "pzvec.h"
#include "fadType.h"

#define Using_FAD

#ifdef Using_FAD
typedef FADREAL Number;
#else
typedef REAL Number;
#endif

void error(char * teste)
{

}

int main()
{

   TPZDiffMatrix<Number> Rot, RotT,
			 X, Xi,
			 Y, Yi,
			 M, Mi,
			 RM, RMi,
			 Temp, Temp2,
			 LambdaSUPG, LambdaB,
			 BornhausCompl;

 {// 2d Tests

   TPZVec<Number> sol;
   TPZVec<REAL> aaS;
   Number us, c;

   //creating data
   sol.Resize(4);
   sol[0] = 1.1;
   sol[1] = 3.4;
   sol[2] = 2.7;
   sol[3] = 10.;

#ifdef Using_FAD
      sol[0].diff(0,4);
      sol[1].diff(1,4);
      sol[2].diff(2,4);
      sol[3].diff(3,4);
#endif

   TPZEulerConsLaw2::uRes(sol, us);
   TPZEulerConsLaw2::cSpeed(sol, 1.4, c);

   cout << "Testing Rotation 2d\n";

   TPZArtDiff::RotMatrix(sol, us, Rot, RotT);

   cout << "Rot" << Rot << endl;
   cout << "RotT" << RotT << endl;
   Rot.Multiply(RotT, Temp);
   cout << "Rot.RotT" << Temp << endl;

   cout << "Testing MMatrix 2d\n";

   TPZArtDiff::MMatrix(sol, us, 1.4, M, Mi);

   cout << "M" << M << endl;
   cout << "Mi" << Mi << endl;
   M.Multiply(Mi, Temp);
   cout << "M.Mi" << Temp << endl;

   cout << "Testing RMMatrix 2d\n";

   TPZArtDiff::RMMatrix(sol, us, 1.4, RM, RMi);

   cout << "RM" << RM << endl;
   cout << "RMi" << RMi << endl;
   RM.Multiply(RMi, Temp);
   cout << "RM.RMi" << Temp << endl;

   RotT.Multiply(M,Temp);
   Temp.Add(RM, -1);
   cout << "RM-RotT.M" << Temp << endl;

   Mi.Multiply(Rot,Temp);
   Temp.Add(RMi, -1);
   cout << "RMi-Mi.Rot" << Temp << endl;

   cout << "Testing SUPG 2d\n";

   TPZArtDiff::EigenSystemSUPG(sol, us, c, 1.4, X, Xi, LambdaSUPG);

   cout << "X" << X << endl;
   cout << "Xi" << Xi << endl;
   X.Multiply(Xi, Temp);
   cout << "X.Xi" << Temp << endl;

   cout << "LambdaSUPG" << LambdaSUPG;

   cout << "Testing whole SUPG decomposition\n";
   RotT. Multiply(M, Temp);
   Temp. Multiply(X, Temp2);
   Temp2.Multiply(LambdaSUPG, Temp);
   Temp. Multiply(Xi, Temp2);
   Temp2.Multiply(Mi, Temp);
   Temp. Multiply(Rot, Temp2);

   cout << "Result of: RotT.M.X.LambdaSUPG.Xi.Mi.Rot (should be 1/Sqr(A2+B2))\n";
   cout << Temp2;

   cout << "Testing Bornhaus 2d\n";

   aaS.Resize(2);
   aaS[0] = 2.1;
   aaS[1] = 12.3;

   TPZArtDiff::EigenSystemBornhaus(sol, us, c, 1.4, aaS, Y, Yi, LambdaB);

   cout << "Y" << Y << endl;
   cout << "Yi" << Yi << endl;
   Y.Multiply(Yi, Temp);
   cout << "Y.Yi" << Temp << endl;

   BornhausCompl.Redim(sol.NElements(), sol.NElements());
   TPZArtDiff::ContributeBornhaus(sol, us, c, 1.4, aaS, BornhausCompl);
   Y. Multiply(LambdaB, Temp);
   Temp.Multiply(Yi, Temp2);
   Temp2.Add(BornhausCompl, -1.);
   cout << "BornhausCompl" << endl;
   cout << BornhausCompl << endl;
   cout << "Y.Lambda.Yi - BornhausCompl" << endl;
   cout << Temp2<< endl;
 }

 {// 3d Tests

   TPZVec<Number> sol;
   TPZVec<REAL> aaS;
   Number us, c;

   //creating data
   sol.Resize(5);
   sol[0] = 1.1;
   sol[1] = 5.4;
   sol[2] = 2.7;
   sol[3] = 4.9;
   sol[4] = 2800.;


#ifdef Using_FAD
      sol[0].diff(0,5);
      sol[1].diff(1,5);
      sol[2].diff(2,5);
      sol[3].diff(3,5);
      sol[4].diff(4,5);
#endif

   TPZEulerConsLaw2::uRes(sol, us);
   TPZEulerConsLaw2::cSpeed(sol, 1.4, c);

   cout << "Testing Rotation 3d\n";

   TPZArtDiff::RotMatrix(sol, us, Rot, RotT);

   cout << "Rot" << Rot << endl;
   cout << "RotT" << RotT << endl;
   Rot.Multiply(RotT, Temp);
   cout << "Rot.RotT" << Temp << endl;

   cout << "Testing MMatrix 3d\n";

   TPZArtDiff::MMatrix(sol, us, 1.4, M, Mi);

   cout << "M" << M << endl;
   cout << "Mi" << Mi << endl;
   M.Multiply(Mi, Temp);
   cout << "M.Mi" << Temp << endl;

   cout << "Testing RMMatrix 3d\n";

   TPZArtDiff::RMMatrix(sol, us, 1.4, RM, RMi);

   cout << "RM" << RM << endl;
   cout << "RMi" << RMi << endl;
   RM.Multiply(RMi, Temp);
   cout << "RM.RMi" << Temp << endl;

   RotT.Multiply(M,Temp);
   Temp.Add(RM, -1);
   cout << "RTM-RotT.M" << Temp << endl;

   Mi.Multiply(Rot,Temp);
   Temp.Add(RMi, -1);
   cout << "RMi-Mi.Rot" << Temp << endl;

   cout << "Testing SUPG 3d\n";

   TPZArtDiff::EigenSystemSUPG(sol, us, c, 1.4, X, Xi, LambdaSUPG);

   cout << "X" << X << endl;
   cout << "Xi" << Xi << endl;
   X.Multiply(Xi, Temp);
   cout << "X.Xi" << Temp << endl;
   TPZEulerConsLaw2::cSpeed(sol, 1.4, c);

   cout << "LambdaSUPG" << LambdaSUPG;

   cout << "Testing whole SUPG decomposition\n";
   RotT. Multiply(M, Temp);
   Temp. Multiply(X, Temp2);
   Temp2.Multiply(LambdaSUPG, Temp);
   Temp. Multiply(Xi, Temp2);
   Temp2.Multiply(Mi, Temp);
   Temp. Multiply(Rot, Temp2);

   cout << Temp2;

   cout << "Testing Bornhaus 3d\n";

   aaS.Resize(3);
   aaS[0] = 2.1;
   aaS[1] = 2.3;
   aaS[2] = 1.7;

   TPZArtDiff::EigenSystemBornhaus(sol, us, c, 1.4, aaS, Y, Yi, LambdaB);

   cout << "Y" << Y << endl;
   cout << "Yi" << Yi << endl;
   Y.Multiply(Yi, Temp);
   cout << "Y.Yi" << Temp << endl;


   BornhausCompl.Redim(sol.NElements(), sol.NElements());
   TPZArtDiff::ContributeBornhaus(sol, us, c, 1.4, aaS, BornhausCompl);
   Y. Multiply(LambdaB, Temp);
   Temp.Multiply(Yi, Temp2);
   Temp2.Add(BornhausCompl, -1.);
   cout << "BornhausCompl" << endl;
   cout << BornhausCompl << endl;
   cout << "Y.Lambda.Yi - BornhausCompl" << endl;
   cout << Temp2<< endl;

 }
}

