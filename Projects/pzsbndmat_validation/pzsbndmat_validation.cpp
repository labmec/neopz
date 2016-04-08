/**
 * @file
 * @brief Trying to validate new symmetric band matrix - implemented in order to
 * make use of LAPACK package
 * @author Francisco Orlandini
 * @since 2016
 */


#include <iostream>
#include <fstream>
#include "pzsbndmat.h"
#include "pzlog.h"
#include "TPZTimer.h"


int main(int argc, char *argv[])
{
  TPZTimer timer;
#ifdef LOG4CXX
  InitializePZLOG();
#endif
  timer.start();
  timer.stop();
  return 0;
}