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

/* **************************************************************************

CLASS NAME :

integv

PURPOSE OF THE CLASS :

TPZIntRule contains a vector of integration points and associated weights.
its constructor and destructor are private. Therefore, pointers to objects
of TPZIntRule can only be obtained through TPZIntRuleList

DEPENDENCIES :

TPZIntRule.h :

pzcondef.h

TPZIntRule.c :

TPZIntRule.h
pzerror.h
	pzcondef.h


ERRORS :

no known errors

ANOMALIES :

When an integration point or weight is selected which is out of range
0.0 is returned and a warning message is generated

WARNING :

check the size of the TPZIntRule after requesting the pointer from TPZIntRuleList
to make sure that the requested integration rule is available

IMPROVEMENTS :

TPZIntRule should be programmed with double instead of double

TPZIntRule could initialize itself based on a file instead of "hardwired"

PROGRAMMERS :

Jose Sergio Rodrigues Alves Filho
Philippe Remy Bernard Devloo

************************************************************************** */

//THIS IMPLEMENT THE INTEGRATION RULE AS A CLASS

class  TPZIntRule  {

  friend class TPZIntRuleList;
  short	   fNumInt;		// number of integration points for this object
  REAL	*fLocation;	// location of the integration point
  REAL	*fWeight;	// weight of the integration point

  TPZIntRule(int i);
  ~TPZIntRule();

 public:

  short NInt(){ return fNumInt;}	//return number of integration points

  REAL Loc(int i);						//return location of the ith pot

  REAL W(int i);						//return weight for the ith point
};

//***********************************************************************
// New Class TPZIntRuleT (Integration Rule for Triangles)
//***********************************************************************

class  TPZIntRuleT  {

  friend class TPZIntRuleList;
  short	   fNumInt;		// number of integration points for this object
  REAL	*fLocationKsi;	// location of the integration point Ksi
  REAL	*fLocationEta;	// location of the integration point Eta
  REAL	*fWeight;		// weight of the integration point

  TPZIntRuleT(int i);
  ~TPZIntRuleT();

 public:

  short NInt(){ return fNumInt;}	//return number of integration points

  void Loc(int i, TPZVec<REAL> &pos);			   //return location of the ith pot

  REAL W(int i);						//return weight for the ith point
};

/* **************************************************************************

CLS NAME :

TPZIntRuleList

PURPOSE OF THE CLASS :

TPZIntRulest is a class from which a unique object is created called
"integ". integ is then used to return pointers to objects of class
TPZIntRule.

DEPENDENCIES :

TPZIntRule.h :

pzcondef.h

TPZIntRule.c :

TPZIntRule.h
pzerror.h
	pzcondef.h

ERRORS :

none

ANOMALIES :

if getrule is called for an integration rule which is not available, a
NULL pointer is returned.

WARNING :

integ should not be deleted

IMPROVEMENTS :

a pointer to an TPZIntRule object could be returned which is nearest to the
one requested

TPZIntRuleList could be extended to contain rules other than the one-d gaussian
tegration rule

PROGRAMMERS :

Jose Sergio Rodrigues Alves Filho
Phippe Remy Bernd Devloo

************************************************************************** */

class TPZIntRuleList {

  int		 	   intavail;	// number of integration rules available
  int        	intavailT;  // number of integration rules available for triangles
  int        	intavailT3D;
  int        	intavailP3D;
  TPZIntRule	**intlist;	 	// pointer to an array of integration rules
  TPZIntRuleT   **intlistT; 	// pointer to an array of integration rules
  TPZIntRuleT3D **intlistT3D;
  TPZIntRuleP3D **intlistP3D;

  public :

    TPZIntRuleList();	// method which initializes all integration rules
  // should be called only once!

  ~TPZIntRuleList();

  TPZIntRule *GetRule(int numint);	// returns a pointer to an integration
  // rule with numint integration points

  TPZIntRuleT *GetRuleT(int numint); // returns a pointer to an integration
  // rule for a triangle
  TPZIntRuleT3D *GetRuleT3D(int numint);
  TPZIntRuleP3D *GetRuleP3D(int numint);
};

extern  TPZIntRuleList  gIntRuleList;

//*******************************************************************
//		Classes TPZIntPoints, TPZInt1d, TPZIntQuad, TPZIntTriang
//*******************************************************************
//*******************************************************************
// Base Class TIntPoint - 	which handles the integration
//*******************************************************************

class TPZIntPoints{
 public:
  virtual int NPoints() = 0;
  virtual void Point(int ip, TPZVec<REAL> &pos, REAL &w) = 0;
  virtual void SetOrder(TPZVec<int> &ord) = 0;
  virtual void GetOrder(TPZVec<int> &ord) = 0;
};

//*******************************************************************
// Base Class TInt1D - 	which handles the integration
//								for 1D problems
//*******************************************************************

class TPZInt1d : public TPZIntPoints{
  int fOrdKsi;
  TPZIntRule *fIntP;
 public:
  TPZInt1d(int OrdK = 0);
  virtual int NPoints();
  virtual void Point(int ip, TPZVec<REAL> &pos, REAL &w);
  virtual void SetOrder(TPZVec<int> &ord);
  virtual void GetOrder(TPZVec<int> &ord);
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
  TPZIntTriang(	int OrdK = 2);
  virtual void SetOrder(TPZVec<int> &ord);
  virtual int  NPoints();
  virtual void Point(int ip, TPZVec<REAL> &pos, REAL &w);
  virtual void GetOrder(TPZVec<int> &ord);
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
  TPZIntQuad(int OrdK = 2, int OrdE = 2);
  virtual void SetOrder(TPZVec<int> &ord);
  virtual int NPoints();
  virtual void Point(int ip, TPZVec<REAL> &pos, REAL &w);
  virtual void GetOrder(TPZVec<int> &ord);
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
  TPZIntCube3D(int OrdK = 2, int OrdE = 2, int OrdZ = 2);
  virtual void SetOrder(TPZVec<int> &ord);
  virtual int NPoints();
  virtual void Point(int ip, TPZVec<REAL> &pos, REAL &w);
  virtual void GetOrder(TPZVec<int> &ord);
};
//#########################################################################
//          Cedric 24/05/98
//***********************************************************************
// New Class TPZIntRuleT3D (Integration Rule for Tetrahedro)
//***********************************************************************

class  TPZIntRuleT3D  {

  friend class TPZIntRuleList;
  short	 fNumInt;		// number of integration points for this object
  REAL	*fLocationKsi;	// location of the integration point Ksi
  REAL	*fLocationEta;	// location of the integration point Eta
  REAL	*fLocationZeta;	// location of the integration point Eta
  REAL	*fWeight;		// weight of the integration point

  TPZIntRuleT3D(int i = 2);
  ~TPZIntRuleT3D();

 public:

  short NInt(){ return fNumInt;}	//return number of integration points

  void Loc(int i, TPZVec<REAL> &pos);			   //return location of the ith pot

  REAL W(int i);						//return weight for the ith point
};
//#########################################################################
//          Cedric 24/05/98
//***********************************************************************
// New Class TPZIntRuleP3D (Rule for Pyramid)
//***********************************************************************

class  TPZIntRuleP3D  {

  friend class TPZIntRuleList;
  short	 fNumInt;		// number of integration points for this object
  REAL	*fLocationKsi;	// location of the integration point Ksi
  REAL	*fLocationEta;	// location of the integration point Eta
  REAL	*fLocationZeta;	// location of the integration point Eta
  REAL	*fWeight;		// weight of the integration point

  TPZIntRuleP3D(int i = 2);
  ~TPZIntRuleP3D();

 public:

  short NInt(){ return fNumInt;}	//return number of integration points

  void Loc(int i, TPZVec<REAL> &pos);			   //return location of the ith pot

  REAL W(int i);						//return weight for the ith point
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
  TPZIntTetra3D(int OrdK = 2);
  virtual void SetOrder(TPZVec<int> &ord);
  virtual int NPoints();
  virtual void Point(int ip, TPZVec<REAL> &pos, REAL &w);
  virtual void GetOrder(TPZVec<int> &ord);
};
//#########################################################################
//#########################################################################
//          Cedric 06/09/98
//*******************************************************************
// Clase Base TPZIntPyram3D - Manipula regra de integracao para
//			      problemas tridimensionais
//  			      elemento Pirâmide
//*******************************************************************
// Cedric
class TPZIntPyram3D : public TPZIntPoints {
  int fOrdKsi;
  TPZIntRuleP3D *fIntKsi;
 public:
  TPZIntPyram3D(int OrdK = 2);
  virtual void SetOrder(TPZVec<int> &ord);
  virtual int NPoints();
  virtual void Point(int ip, TPZVec<REAL> &pos, REAL &w);
  virtual void GetOrder(TPZVec<int> &ord);
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

  TPZIntPrism3D(int OrdK = 2,int OrdL = 2);
  virtual ~TPZIntPrism3D();
  void SetOrder(TPZVec<int> &ord) ;
  int NPoints();
  void Point(int ip, TPZVec<REAL> &pos, REAL &w);
  void GetOrder(TPZVec<int> &ord);
};
//#########################################################################
//#########################################################################
//          Cedric 25/04/99
//*******************************************************************
// Clase Base TPZInt1Point - Regra de integracao para
//                            1 ponto
//*******************************************************************
// Cedric
class TPZInt1Point  : public TPZIntPoints {
  int fOrdKsi;
  
 public:

  TPZInt1Point();
  virtual ~TPZInt1Point();
  void SetOrder(TPZVec<int> &ord) ;
  int NPoints();
  void Point(int ip, TPZVec<REAL> &pos, REAL &w);
  void GetOrder(TPZVec<int> &ord);
};

#endif
