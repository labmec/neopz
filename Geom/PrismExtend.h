#pragma once


#include "pzstack.h"
#include "pztrnsform.h"
#include "pzquad.h"
#include "tpzprinteg.h"
#include "tpzpoint.h"

namespace pztopology 
{

/// This class defines the Prismatic extension of a topology
template<class TFather>
class Pr :
	public TFather
{
public:

	enum {NCornerNodes = 2*TFather::NCornerNodes, NSides = 3*TFather::NSides, Dimension = TFather::Dimension+1};

    Pr();

    virtual ~Pr();

static void LowerDimensionSides(int side,TPZStack<int> &smallsides);
static void LowerDimensionSides(int side,TPZStack<int> &smallsides, int DimTarget);

 /**
  * returns all sides whose closure contains side
  * @param side smaller dimension side
  * @param high vector which will contain all sides whose closure contain sidefrom
  */
static void HigherDimensionSides(int side, TPZStack<int> &high);
 /**
  * return the number of nodes (not connectivities) associated with a side
  */
static int NSideNodes(int side);
 /**
  * returns the local node number of the node "node" along side "side"
  */
static int SideNodeLocId(int side, int node);
/**
  * returns the barycentric coordinates in the master element space of the original element
  */
static void CenterPoint(int side, TPZVec<REAL> &center);
/**volume of the master element*/
static REAL RefElVolume();
/**
  * returns the dimension of the side
  */
static int SideDimension(int side);

/**
  * returns the transformation which takes a point from the side sidefrom ot
  * the side sideto
  * @param sidefrom side where the point resides
  * @param sideto side whose closure contains sidefrom
  */
static TPZTransform SideToSideTransform(int sidefrom, int sideto);
/**
 * Returns the transformation which transform a point from the interior of the element to the side
 * @param side side to which the point will be tranformed (0<=side<=2)
 * @return TPZTransform object
 * @see the class TPZTransform
 */
static TPZTransform TransformElementToSide(int side);
/**
 * Returns the transformation which transform a point from the side to the interior of the element
 * @param side side from which the point will be tranformed (0<=side<=2)
 * @return TPZTransform object
 * @see the class TPZTransform
 */
static TPZTransform TransformSideToElement(int side);

static TPZIntPoints *CreateSideIntegrationRule(int side, int order);

typedef typename TFather::IntruleType fatintruletype;
typedef TPZPrInteg<fatintruletype> IntruleType;
// typedef TPZGraphEl1dd GraphElType;

 /**
   * return the type of the element as specified in file pzeltype.h
   */
static std::string StrType() ;//{ return EOned;}

/**
  * return the type of the element as specified in file pzeltype.h
  */
static std::string StrType(int side);

/**
 * computes the linear map from an internal point to the parameter space of the side
 * returns the jacobian of the transformation
 */
static bool MapToSide(int side, TPZVec<REAL> &InternalPar, TPZVec<REAL> &SidePar, TPZFMatrix &JacToSide);

/**
 * Number of connects of the element (3)
 * @return number of connects of the element
 */
static int NumSides();

 /**
  * return the number of nodes (not connectivities) associated with a side
  */
static int NContainedSides(int side);
 /**
  * returns the local connect number of the connect "c" along side "side"
  */
static int ContainedSideLocId(int side, int c);

/**
 * uses log4cxx to print the results of all methods
 */
static void Diagnostic();
};

}

