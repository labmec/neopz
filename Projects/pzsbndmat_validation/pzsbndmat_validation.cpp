/**
 * @file
 * @brief Trying to validate new symmetric band matrix - implemented in order to
 * make use of LAPACK package
 * @author Francisco Orlandini
 * @since 2016
 */


#include <iostream>
#include <fstream>
#include "TPZSBMatrixLapack.h"
#include "pzfmatrix.cpp"
#include "pzlog.h"
#include "TPZTimer.h"
#include "pzreal.h"


int main(int argc, char *argv[])
{
  TPZTimer timer;

  timer.start();
	/*float*/
	TPZSBMatrixLapack < float > aFloat (3 , 2);
	aFloat.PutVal(0,0,4);
	aFloat.PutVal(1,1,37);
	aFloat.PutVal(2,2,98);
	
	aFloat.PutVal(0,1,12);
	aFloat.PutVal(1,2,-43);
	
	aFloat.PutVal(0,2,-16);
	
	aFloat.Print("a" , std::cout );
	
	aFloat.Decompose_Cholesky();
	TPZFMatrix< float > bFloat( 3 , 1);
	bFloat(0 , 0) = 28;
	bFloat(1 , 0) = 86;
	bFloat(2 , 0) = -102;
	bFloat.Print("b" , std::cout);
	aFloat.Print("a_decomposed" , std::cout );
	aFloat.Subst_Forward(&bFloat);
	bFloat.Print("b_f" , std::cout);
	aFloat.Subst_Backward(&bFloat);
	bFloat.Print("x" , std::cout);
	aFloat.Resize(10);
	aFloat.Print("a_decomposed" , std::cout );
	aFloat.SetBand(4);
	aFloat.Print("a_decomposed" , std::cout );
	/*double*/
  TPZSBMatrixLapack < double > aDouble (3 , 2);
	aDouble.PutVal(0,0,4);
	aDouble.PutVal(1,1,37);
	aDouble.PutVal(2,2,98);
	
	aDouble.PutVal(0,1,12);
	aDouble.PutVal(1,2,-43);
	
	aDouble.PutVal(0,2,-16);
  
  aDouble.Print("a" , std::cout );
	
	aDouble.Decompose_Cholesky();
	
	aDouble.Print("a" , std::cout );
	
	/*complex<double>*/
	TPZSBMatrixLapack < std::complex<double> > aCDouble (3 , 2);
	aCDouble.PutVal(0,0,4);
	aCDouble.PutVal(1,1,37);
	aCDouble.PutVal(2,2,98);
	
	aCDouble.PutVal(0,1,12);
	aCDouble.PutVal(1,2,-43);
	
	aCDouble.PutVal(0,2,-16);
	
	aCDouble.Print("a" , std::cout );
	
	aCDouble.Decompose_Cholesky();
	
	aCDouble.Print("a" , std::cout );
  timer.stop();
  return 0;
}