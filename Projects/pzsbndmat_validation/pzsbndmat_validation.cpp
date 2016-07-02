/**
 * @file
 * @brief Trying to validate new symmetric band matrix - implemented in order to
 * make use of LAPACK package
 * @author Francisco Orlandini
 * @since 2016
 */


#include <iostream>
#include <fstream>
#include "pzfmatrix.cpp"
#include "pzlog.h"
#include "TPZTimer.h"
#include "pzreal.h"


int main(int argc, char *argv[])
{
  TPZTimer timer;

	
	/*float*/
	TPZSBMatrixLapack < float > aFloat (3 , 2);
	aFloat.PutVal(0,0,4);
	aFloat.PutVal(1,1,37);
	aFloat.PutVal(2,2,98);
	
	aFloat.PutVal(0,1,12);
	aFloat.PutVal(1,2,-43);
	
	aFloat.PutVal(0,2,-16);
	
	aFloat.Print("a" , std::cout );
	
	
	TPZFMatrix< float > bFloat( 3 , 1);
	bFloat(0 , 0) = 28;
	bFloat(1 , 0) = 86;
	bFloat(2 , 0) = -102;
	timer.start();
	aFloat.Decompose_Cholesky();
	aFloat.Subst_Forward(&bFloat);
	aFloat.Subst_Backward(&bFloat);
	timer.stop();
	std::cout<<"time: "<<timer<<std::endl;
	timer.reset();
	bFloat.Print("x" , std::cout);
	/*double*/
  TPZSBMatrixLapack < double > aDouble (3 , 2);
	aDouble.PutVal(0,0,4);
	aDouble.PutVal(1,1,37);
	aDouble.PutVal(2,2,98);
	
	aDouble.PutVal(0,1,12);
	aDouble.PutVal(1,2,-43);
	
	aDouble.PutVal(0,2,-16);
  
  aDouble.Print("a" , std::cout );
	
	
	TPZFMatrix< double > bDouble( 3 , 1);
	bDouble(0 , 0) = 28;
	bDouble(1 , 0) = 86;
	bDouble(2 , 0) = -102;
	bDouble.Print("b" , std::cout);
	timer.start();
	aDouble.Solve_LinSys(bDouble);
	timer.stop();
	std::cout<<"time: "<<timer<<std::endl;
	timer.reset();
	bDouble.Print("sol" , std::cout);
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