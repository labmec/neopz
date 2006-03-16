// -*- c++ -*-

//$Id: performance.h,v 1.1 2006-03-16 13:25:52 tiago Exp $

#include "TPZTimer.h"
#include "pzreal.h"
#include <iostream>

class TPZPerformance{
private:
#ifdef contar
  TPZCounter fStartFlops, fStopFlops;
#else
  TPZTimer fTimer;
#endif
public:
  TPZPerformance();
  ~TPZPerformance();    
  void Start();  
  void Stop();  
  void Print(std::ostream &out);
};
