// -*- c++ -*-

//$Id: performance.cpp,v 1.1 2006-03-16 13:25:52 tiago Exp $

#include "performance.h"
#include "pzreal.h"
#include "TPZTimer.h"

TPZPerformance::TPZPerformance(){
#ifdef contar
//   this->fStartFlops = 0;
//   this->fStopFlops = 0;
#else
  //NOTHING TO BE DONE
#endif 
}

TPZPerformance::~TPZPerformance(){ 
  //NOTHING TO BE DONE
}

void TPZPerformance::Start(){
#ifdef contar
  this->fStartFlops = std::TPZFlopCounter::gCount;
#else
  this->fTimer.start();
#endif
}

void TPZPerformance::Stop(){
#ifdef contar
  this->fStopFlops = std::TPZFlopCounter::gCount;
#else
  this->fTimer.stop();
#endif
}

void TPZPerformance::Print(std::ostream &out){
#ifdef contar
  const TPZCounter noper = this->fStopFlops - this->fStartFlops;
  out << "Number of operations: " << noper << std::endl;
#else
  out << this->fTimer << std::endl;
#endif
}