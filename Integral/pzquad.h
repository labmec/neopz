//HEADER FILE FOR THE INTEGRATION VECTOR

// _*_ c++ _*_

#ifndef INTQUADHPP
#define INTQUADHPP

#include <iostream>
using namespace std;
#include "pzreal.h"
template<class T>
class TPZVec;
template<class T, int N>
class TPZManVector;

class TPZIntRuleT3D;
class TPZIntRuleP3D;
class TPZIntRule;
#include "tpzintrulet.h"
#include "tpzintpoints.h"
#include "tpzprinteg.h"

//*******************************************************************
// Base Class TInt1D - 	which handles the integration
//								for 1D problems
//*******************************************************************

class TPZInt1d : public TPZIntPoints{
  int fOrdKsi;
  TPZIntRule *fIntP;
 public:
   enum {Dim = 1};
  TPZInt1d(int OrdK = 0);
  virtual ~TPZInt1d()
  {
  }
  virtual int NPoints();
  virtual void Point(int ip, TPZVec<REAL> &pos, REAL &w);
  virtual void SetOrder(TPZVec<int> &ord);
  virtual void GetOrder(TPZVec<int> &ord);
  virtual int GetMaxOrder();
  virtual int Dimension() 
  {
    return Dim;
  }
  virtual TPZIntPoints *PrismExtend(int order)
  {
    return new TPZPrInteg<TPZInt1d>(order);
  }
};
//*******************************************************************
// Base Class TPZIntTriang - which handles the integration
//									for 2D problems, triangle
//									elements
//*******************************************************************

class TPZIntTriang : public TPZIntPoints{
  int fOrdKsi;
  TPZIntRuleT *fIntKsi;
 public:
   enum {Dim = 2};
  TPZIntTriang(	int OrdK = 2);
  virtual ~TPZIntTriang()
  {
  }
  virtual void SetOrder(TPZVec<int> &ord);
  virtual int  NPoints();
  virtual void Point(int ip, TPZVec<REAL> &pos, REAL &w);
  virtual void GetOrder(TPZVec<int> &ord);
  virtual int GetMaxOrder();
  virtual int Dimension()
  {
    return Dim;
  }
  virtual TPZIntPoints *PrismExtend(int order)
  {
    return new TPZPrInteg<TPZIntTriang>(order);
  }

};
//*******************************************************************
// Base Class TPZIntQuad - 	which handles the integration
//									for 2D problems, quadrilaterals
//									elements
//*******************************************************************

class TPZIntQuad : public TPZIntPoints{
  int fOrdKsi;
  int fOrdEta;
  TPZIntRule *fIntKsi;
  TPZIntRule *fIntEta;
 public:
   enum {Dim = 2};
  TPZIntQuad(int OrdK = 2, int OrdE = 2);
  virtual ~TPZIntQuad()
  {
  }
  virtual void SetOrder(TPZVec<int> &ord);
  virtual int NPoints();
  virtual void Point(int ip, TPZVec<REAL> &pos, REAL &w);
  virtual void GetOrder(TPZVec<int> &ord);
  virtual int GetMaxOrder();  
  virtual int Dimension()
  {
    return Dim;
  }
  virtual TPZIntPoints *PrismExtend(int order)
  {
    return new TPZPrInteg<TPZIntQuad>(order);
  }
};

//#########################################################################
//          Cedric 24/04/98
//*******************************************************************
// Clase Base TPZIntCube3D - 	Manipula regra de integracao para
//			        problemas tridimensionais
//   				elemento hexaedro
//*******************************************************************
// Cedric

class TPZIntCube3D : public TPZIntPoints{
  int fOrdKsi;
  int fOrdEta;
  int fOrdZeta;
  TPZIntRule *fIntKsi;
  TPZIntRule *fIntEta;
  TPZIntRule *fIntZeta;
 public:
   enum {Dim = 3};
  TPZIntCube3D(int OrdK = 2, int OrdE = 2, int OrdZ = 2);
  virtual ~TPZIntCube3D()
  {
  }
  virtual void SetOrder(TPZVec<int> &ord);
  virtual int NPoints();
  virtual void Point(int ip, TPZVec<REAL> &pos, REAL &w);
  virtual void GetOrder(TPZVec<int> &ord);
  virtual int GetMaxOrder();
  virtual int Dimension()
  {
    return Dim;
  }
  virtual TPZIntPoints *PrismExtend(int order)
  {
    return new TPZPrInteg<TPZIntCube3D>(order);
  }
};

//#########################################################################
//          Cedric 24/05/98
//*******************************************************************
// Clase Base TPZIntTetra3D - Manipula regra de integracao para
//			      problemas tridimensionais
//  			      elemento Tetraedro
//*******************************************************************
// Cedric
class TPZIntTetra3D : public TPZIntPoints {
  int fOrdKsi;
  TPZIntRuleT3D *fIntKsi;
 public:
   enum {Dim = 3};
  TPZIntTetra3D(int OrdK = 2);
  virtual void SetOrder(TPZVec<int> &ord);
  virtual int NPoints();
  virtual void Point(int ip, TPZVec<REAL> &pos, REAL &w);
  virtual void GetOrder(TPZVec<int> &ord);
  virtual int GetMaxOrder();  
  virtual int Dimension()
  {
    return Dim;
  }
  virtual TPZIntPoints *PrismExtend(int order)
  {
    return new TPZPrInteg<TPZIntTetra3D>(order);
  }
};
//#########################################################################
//#########################################################################
//          Cedric 06/09/98
//*******************************************************************
// Clase Base TPZIntPyram3D - Manipula regra de integracao para
//			      problemas tridimensionais
//  			      elemento Pir�mide
//*******************************************************************
// Cedric
class TPZIntPyram3D : public TPZIntPoints {
  int fOrdKsi;
  TPZIntRuleP3D *fIntKsi;
 public:
   enum {Dim =3};
  TPZIntPyram3D(int OrdK = 2);
  virtual void SetOrder(TPZVec<int> &ord);
  virtual int NPoints();
  virtual void Point(int ip, TPZVec<REAL> &pos, REAL &w);
  virtual void GetOrder(TPZVec<int> &ord);
  virtual int GetMaxOrder();  
  virtual int Dimension()
  {
    return Dim;
  }
  virtual TPZIntPoints *PrismExtend(int order)
  {
    return new TPZPrInteg<TPZIntPyram3D>(order);
  }
};
//#########################################################################
//#########################################################################
//          Cedric 06/09/98
//*******************************************************************
// Clase Base TPZIntPrism3D - Manipula regra de integracao para
//			      problemas tridimensionais
//  			      elemento Prisma
//*******************************************************************
// Cedric
class TPZIntPrism3D  : public TPZIntPoints {
  int fOrdKsi,fOrdKti;
  TPZInt1d fIntRule1D;
  TPZIntTriang fIntTriang;
 public:

   enum {Dim = 3};
  TPZIntPrism3D(int OrdK = 2,int OrdL = 2);
  virtual ~TPZIntPrism3D();
  void SetOrder(TPZVec<int> &ord) ;
  int NPoints();
  void Point(int ip, TPZVec<REAL> &pos, REAL &w);
  void GetOrder(TPZVec<int> &ord);
  virtual int GetMaxOrder();  
  virtual int Dimension()
  {
    return Dim;
  }
  virtual TPZIntPoints *PrismExtend(int order)
  {
    return new TPZPrInteg<TPZIntPrism3D>(order);
  }
};

#endif
