#include "pzlog.h"
#include "TPZMatfracinterface.h"
#include "pzbndcond.h"
#include "pzfmatrix.h"
#include "pzaxestools.h"

#include <iostream>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.multiphase"));
#endif

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.material.multiphase.data"));
#endif

TPZMatfracinterface::TPZMatfracinterface(): TPZDiscontinuousGalerkin()
{
  fData = NULL;
}

TPZMatfracinterface::TPZMatfracinterface(int matid): TPZDiscontinuousGalerkin(matid)
{
  fData = NULL;
}

TPZMatfracinterface::~TPZMatfracinterface()
{
  
}

int TPZMatfracinterface::Dimension() const {return 1;};

int TPZMatfracinterface::NStateVariables() {return 1;}

void TPZMatfracinterface::Print(std::ostream &out) {
  out << "name of material : " << Name() << "\n";
  out << "Coeficient which multiplies the gradient operator "<< "my var" << std::endl;
  out << "Base Class properties :";
  TPZMaterial::Print(out);
  out << "\n";
}


void TPZMatfracinterface::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, TPZVec<TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
  int globPtIndex = dataright[0].intGlobPtIndex;
  
  
  const REAL ql=1.0; // AQUINATHAN prencher com o certo depois
  
}