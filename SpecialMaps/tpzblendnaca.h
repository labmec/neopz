/**
 * @file
 * @brief Contains the TPZBlendNACA class. It is a special map.
 */

#ifndef TPZBLENDNACA_H
#define TPZBLENDNACA_H

#include <math.h>
#include "pzreal.h"
#include "pzvec.h"
#include "pzeltype.h"

class TPZGeoEl;
class TPZGeoMesh;

/**
 * @author caju2008 < caju\@skol >
 * @ingroup geometry
 * @brief Special map to NACA. \ref geometry "Geometry"
 */
class TPZBlendNACA
{
public:
	/** @brief Default constructor */
    TPZBlendNACA();
    /** @brief Constructor */
    TPZBlendNACA(REAL cord, int FourDigits, REAL angle, TPZVec<REAL> &x0);
	/** @brief Default destructor */
    ~TPZBlendNACA();
	
    /// with attack angle
    /// superior profile
    template <class Type>
    Type xua(Type x);
    template <class Type>
    Type yua(Type x);
	
    /// inferior profile
    template <class Type>
    Type xla(Type x);
    template <class Type>
    Type yla(Type x);
    
    template <class Type>
    void ProjectPoint(TPZVec<Type> &pt, int maxPt = 1000);
	
public:
    
    virtual void ParametricDomainNodeCoord(int node, TPZVec<REAL> &nodeCoord)
    {
        //Dont have
    }
    
	/** @brief Creates a geometric element according to the type of the father element */
	// static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
	// 								  TPZVec<int64_t>& nodeindexes,
	// 								  int matid,
	// 								  int64_t& index);
private:
	
    template<class Type>
    Type P();
    template<class Type>
    Type M();
    template<class Type>
    Type TT();
	
    /** @brief Mean line for the wing */
    template <class Type>
    Type yc(Type x);
    template <class Type>
    Type dyc(Type x);
	
    /** @brief Thickness */
    template <class Type>
    Type yt(Type x);
	
    /** @brief Superior profile */
    template <class Type>
    Type xu(Type x);
    template <class Type>
    Type yu(Type x);
	
    /** @brief Inferior profile */
    template <class Type>
    Type xl(Type x);
    template <class Type>
    Type yl(Type x);
    
	template <class Type>
    void NearestParameter(TPZVec<REAL> &pt, int &uplow, int maxPt,Type &);

public:
    
    REAL fCord;
    int  fFourDigits;
    REAL fAngle;
    REAL fX0[3];
    REAL fP;
    REAL fM;
    REAL fTT;
};
#endif
