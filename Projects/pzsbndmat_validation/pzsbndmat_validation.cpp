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
#include "pzlog.h"
#include "TPZTimer.h"
#include "pzreal.h"


int main(int argc, char *argv[])
{
  TPZTimer timer;
#ifdef LOG4CXX
  InitializePZLOG();
#endif
  timer.start();
  TPZSBMatrixLapack < double > a (4 , 1);
  TPZSBMatrixLapack < double > b (4 , 2);
  TPZSBMatrixLapack < double > c (4 , 2);
  a.PutVal(0,0,1);
  a.PutVal(1,1,1);
  a.PutVal(2,2,1);
  a.PutVal(3,3,1);
  
  a.PutVal(0,1,-1);
  a.PutVal(1,2,-1);
  a.PutVal(2,3,-1);
  
  b.PutVal(0,0,1);
  b.PutVal(1,1,1);
  b.PutVal(2,2,1);
  b.PutVal(3,3,1);
  
  b.PutVal(0,1,-1);
  b.PutVal(1,2,-1);
  b.PutVal(2,3,-1);
  
  b.PutVal(0,2,4);
  b.PutVal(1,3,4);
  
  a.Print("a" , std::cout );
  b.Print("b" , std::cout );
  
  c = a+b;
  c.Print("c" , std::cout);
  c = a-b;
  c.Print("c" , std::cout);
  timer.stop();
  return 0;
}