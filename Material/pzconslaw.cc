$Id: pzconslaw.cc,v 1.2 2003-10-20 11:59:22 erick Exp $

#include "pzconslaw.h"
#include "pzelmat.h"
#include "pzbndcond.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include <math.h>
#include <stdlib.h>


TPZConservationLaw2::TPZConservationLaw2(int nummat,REAL timeStep,int dim) :
                                                      TPZMaterial(nummat),
fDim(dim),
fTimeStep(0),
fCFL(0),
//fDelta(0),
fGamma(1.4)
{
   fTimeStep = timeStep;
   if(timeStep < 0 || timeStep > 1)
   {
       PZError << "TPZConservationLaw2::TPZConservationLaw2 time step parameter > 1 , default 1.0\n";
       fTimeStep = 1.0;
   }

   if(dim < 1 || dim > 3)
   {
      PZError << "TPZConservationLaw2::TPZConservationLaw2 (abort) error dimension = " << dim << endl;
      exit(-1);
   }
   fDim = dim;
}


TPZConservationLaw2::~TPZConservationLaw2()
{
}

void TPZConservationLaw2::Print(ostream &out)
{
   out << "name of material : " << Name() << "\n";
   out << "properties : \n";
   out << "\tdimension: " << fDim << endl;
   out << "\ttime step: " << fTimeStep << endl;
   out << "\tCFL: " << fCFL << endl;
//   out << "\tDelta (diffusive term): " << fDelta << endl;
   out << "\tGamma: " << fGamma << endl;
   TPZMaterial::Print(out);
}

// void TPZConservationLaw2::Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &/*axes*/,int var,TPZVec<REAL> &Solout){

//   cout << "TPZConservationLaw2::Solution nao deve ser chamada\n";

// }

void TPZConservationLaw2::Flux(TPZVec<REAL> &/*x*/, TPZVec<REAL> &/*Sol*/, TPZFMatrix &/*DSol*/, TPZFMatrix &/*axes*/, TPZVec<REAL> &/*flux*/) {
  //Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix &DSol, TPZFMatrix &axes, TPZVec<REAL> &flux)
}

void TPZConservationLaw2::Errors(TPZVec<REAL> &/*x*/,TPZVec<REAL> &u,
			       TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &/*flux*/,
			       TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values) {

}
