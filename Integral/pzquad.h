//HEADER FILE FOR THE INTEGRATION VECTOR

// _*_ c++ _*_

#ifndef INTQUADHPP
#define INTQUADHPP

#include <iostream>

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
	TPZInt1d(const TPZInt1d &copy ) : TPZIntPoints(copy), fOrdKsi(copy.fOrdKsi), fIntP(copy.fIntP)
	{
	}
  virtual ~TPZInt1d()
  {
  }
  virtual int NPoints() const;
  virtual void Point(int ip, TPZVec<REAL> &pos, REAL &w) const;
  virtual void SetOrder(TPZVec<int> &ord);
  virtual void GetOrder(TPZVec<int> &ord) const;
  virtual int GetMaxOrder() const;
  virtual int Dimension() const
  {
    return Dim;
  }
  virtual TPZIntPoints *PrismExtend(int order)
  {
    return new TPZPrInteg<TPZInt1d>(order);
  }
	virtual TPZIntPoints *Clone() const
	{
		return new TPZInt1d(*this);
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
	TPZIntTriang(const TPZIntTriang &copy) : TPZIntPoints(copy), fOrdKsi(copy.fOrdKsi), fIntKsi(copy.fIntKsi)
	{
	}
  virtual void SetOrder(TPZVec<int> &ord);
  virtual int  NPoints() const;
  virtual void Point(int ip, TPZVec<REAL> &pos, REAL &w) const;
  virtual void GetOrder(TPZVec<int> &ord) const;
  virtual int GetMaxOrder() const;
  virtual int Dimension() const
  {
    return Dim;
  }
  virtual TPZIntPoints *PrismExtend(int order)
  {
    return new TPZPrInteg<TPZIntTriang>(order);
  }
	virtual TPZIntPoints *Clone() const
	{
		return new TPZIntTriang(*this);
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
	TPZIntQuad(const TPZIntQuad &copy) : TPZIntPoints(copy), fOrdKsi(copy.fOrdKsi), fOrdEta(copy.fOrdEta), fIntKsi(copy.fIntKsi),fIntEta(copy.fIntEta)
	{
	}
  virtual void SetOrder(TPZVec<int> &ord);
  virtual int NPoints() const;
  virtual void Point(int ip, TPZVec<REAL> &pos, REAL &w) const;
  virtual void GetOrder(TPZVec<int> &ord) const;
  virtual int GetMaxOrder() const;  
  virtual int Dimension() const
  {
    return Dim;
  }
  virtual TPZIntPoints *PrismExtend(int order)
  {
    return new TPZPrInteg<TPZIntQuad>(order);
  }
	virtual TPZIntPoints* Clone() const
	{
		return new TPZIntQuad(*this);
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
	TPZIntCube3D(const TPZIntCube3D &copy) : TPZIntPoints(copy), fOrdKsi(copy.fOrdKsi), fOrdEta(copy.fOrdEta), fOrdZeta(copy.fOrdZeta),
		fIntKsi(copy.fIntKsi), fIntEta(copy.fIntEta), fIntZeta(copy.fIntZeta)
	{
	}
  virtual void SetOrder(TPZVec<int> &ord);
  virtual int NPoints() const;
  virtual void Point(int ip, TPZVec<REAL> &pos, REAL &w) const;
  virtual void GetOrder(TPZVec<int> &ord) const;
  virtual int GetMaxOrder() const;
  virtual int Dimension() const
  {
    return Dim;
  }
  virtual TPZIntPoints *PrismExtend(int order)
  {
    return new TPZPrInteg<TPZIntCube3D>(order);
  }
	virtual TPZIntPoints *Clone() const
	{
		return new TPZIntCube3D(*this);
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
	TPZIntTetra3D(const TPZIntTetra3D &copy) : TPZIntPoints(copy), fOrdKsi(copy.fOrdKsi), fIntKsi(copy.fIntKsi)
	{
	}
  virtual void SetOrder(TPZVec<int> &ord);
  virtual int NPoints() const;
  virtual void Point(int ip, TPZVec<REAL> &pos, REAL &w) const;
  virtual void GetOrder(TPZVec<int> &ord) const;
  virtual int GetMaxOrder() const;  
  virtual int Dimension() const
  {
    return Dim;
  }
  virtual TPZIntPoints *PrismExtend(int order)
  {
    return new TPZPrInteg<TPZIntTetra3D>(order);
  }
	virtual TPZIntPoints* Clone() const
	{
		return new TPZIntTetra3D(*this);
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
	TPZIntPyram3D(const TPZIntPyram3D &copy) : TPZIntPoints(copy), fOrdKsi(copy.fOrdKsi), fIntKsi(copy.fIntKsi)
	{
	}
  virtual void SetOrder(TPZVec<int> &ord);
  virtual int NPoints() const;
  virtual void Point(int ip, TPZVec<REAL> &pos, REAL &w) const;
  virtual void GetOrder(TPZVec<int> &ord) const;
  virtual int GetMaxOrder() const;  
  virtual int Dimension() const
  {
    return Dim;
  }
  virtual TPZIntPoints *PrismExtend(int order)
  {
    return new TPZPrInteg<TPZIntPyram3D>(order);
  }
	TPZIntPoints *Clone() const
	{
		return new TPZIntPyram3D(*this);
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
	TPZIntPrism3D(const TPZIntPrism3D &copy) : TPZIntPoints(copy), fOrdKsi(copy.fOrdKsi), fOrdKti(copy.fOrdKti), fIntRule1D(copy.fIntRule1D),
			fIntTriang(copy.fIntTriang)
	{
	}
  virtual ~TPZIntPrism3D();
  void SetOrder(TPZVec<int> &ord) ;
  int NPoints() const;
  void Point(int ip, TPZVec<REAL> &pos, REAL &w) const;
  void GetOrder(TPZVec<int> &ord) const;
  virtual int GetMaxOrder() const;  
  virtual int Dimension() const
  {
    return Dim;
  }
  virtual TPZIntPoints *PrismExtend(int order)
  {
    return new TPZPrInteg<TPZIntPrism3D>(order);
  }
	TPZIntPoints *Clone() const
	{
		return new TPZIntPrism3D(*this);
	}
};

#endif
