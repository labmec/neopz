#ifndef CONVTEST_H
#define CONVTEST_H


#include "pzvec.h"
#include "pzgeoel.h"
#include "pzgeoelside.h"

  /**
  / Class made by Paulo Cesar de Alvarenga Lucci (Caju)
  / LabMeC - FEC - UNICAMP
  / 2007
 */

class ConvTest
{

public:

    ConvTest();
   ~ConvTest();

    void JacobianConv(TPZGeoEl &Object, TPZVec< REAL > QsiEta);
    void JacobianConv(TPZGeoElSide &Object, TPZVec< REAL > QsiEta);

};

#endif
