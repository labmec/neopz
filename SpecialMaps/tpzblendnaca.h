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
#include "pznoderep.h"
#include "tpzline.h"

class TPZGeoEl;
class TPZGeoMesh;

/**
 * @author caju2008 < caju\@skol >
 * @ingroup geometry
 * @brief Special map to NACA. \ref geometry "Geometry"
 */

namespace pzgeom
{

class TPZBlendNACA  : public pzgeom::TPZNodeRep<2,pztopology::TPZLine>
{
public:
    typedef pztopology::TPZLine Top;
    /** @brief Number of nodes (connects) */
    enum {NNodes = 2};

	/** @brief Default constructor */
    TPZBlendNACA();
    /** @brief Constructor */
    TPZBlendNACA(REAL cord, int FourDigits, REAL angle, TPZVec<REAL> &x0);

    TPZBlendNACA(TPZVec<int64_t> &nodeindexes) : TPZRegisterClassId(&TPZBlendNACA::ClassId),
        pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(nodeindexes) {

    }

    TPZBlendNACA(const TPZBlendNACA &cp, TPZGeoMesh &) : TPZRegisterClassId(&TPZBlendNACA::ClassId),
        pzgeom::TPZNodeRep<NNodes, pztopology::TPZLine>(cp){
            DebugStop();
        }
        /** @brief Copy constructor with map of nodes */
    TPZBlendNACA(const TPZBlendNACA &cp,std::map<int64_t,int64_t> & gl2lcNdMap) : 
        TPZRegisterClassId(&TPZBlendNACA::ClassId),
        pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(cp,gl2lcNdMap){
            DebugStop();
        
        }
	/** @brief Default destructor */
    ~TPZBlendNACA();
	
    /// with attack angle
    /// superior profile
    template <class Type>
    Type xua(Type x) const;
    template <class Type>
    Type yua(Type x) const;
	
    template <class Type>
    Type dxua(Type x) const;
    template <class Type>
    Type dyua(Type x) const;
	
    /// inferior profile
    template <class Type>
    Type xla(Type x) const;
    template <class Type>
    Type yla(Type x) const;
    
    /// inferior profile
    template <class Type>
    Type dxla(Type x) const;
    template <class Type>
    Type dyla(Type x) const;
    
    template <class Type>
    void ProjectPoint(TPZVec<Type> &pt, int maxPt = 1000);
	
public:
    
    // declare the methods needed to implement the interface of TPZGeoElRefPattern<TPZBlendNACA>
    

    virtual void ParametricDomainNodeCoord(int node, TPZVec<REAL> &nodeCoord)
    {
        //Dont have
    }
    
    int ClassId() const override;

	/** @brief Creates a geometric element according to the type of the father element */
	// static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
	// 								  TPZVec<int64_t>& nodeindexes,
	// 								  int matid,
	// 								  int64_t& index);
protected:
	
    template<class Type>
    Type P();
    template<class Type>
    Type M();
    template<class Type>
    Type TT();
	
    /** @brief Mean line for the wing */
    template <class Type>
    Type yc(Type x) const;

    /** @brief Derivative of the mean line */
    template <class Type>
    Type dyc(Type x) const;
	
    /** @brief Derivative of the mean line */
    template <class Type>
    Type tgphi(Type x) const{
        return dyc(x);
    }
	
    /** @brief Second derivative of the mean line */
    template <class Type>
    Type dtgphi(Type x) const;
	
    /** @brief Thickness */
    template <class Type>
    Type yt(Type x) const;
	
    /** @brief Derivative of the thickness */
    template <class toto>
    toto dyt(toto x) const;

    /** @brief Superior profile */
    template <class Type>
    Type xu(Type x) const;
    template <class Type>
    // Derivative of the superior profile
    Type dxu(Type x) const;

    template <class Type>
    Type yu(Type x) const;
    template <class Type>
    Type dyu(Type x) const;
	
    /** @brief Inferior profile */
    template <class Type>
    Type xl(Type x) const;
    template <class Type>
    Type yl(Type x) const;
    /** @brief Inferior profile */
    template <class Type>
    Type dxl(Type x) const;
    template <class Type>
    Type dyl(Type x) const;
    
public:
	template <class Type>
    void NearestParameter(TPZVec<REAL> &pt, int &uplow, int maxPt,Type &);

public:

    void Read(TPZStream& buf, void* context) override {
        DebugStop();
    }
    void Write(TPZStream &buf, int withclassid) const override
    {
        DebugStop();
    }

    void Initialize(TPZGeoEl *refel)
    {
        DebugStop();
    }
    
    template<class T>
    void X(TPZFMatrix<REAL> &coord,TPZVec<T> &loc,TPZVec<T> &result) const
    {
        DebugStop();
    }   

    template<class T>
    void GradX(TPZFMatrix<REAL> &coord, TPZVec<T> &par, TPZFMatrix<T> &gradx) const
    {
        DebugStop();
    }
    static std::string TypeName() { return "NACA";}
    
    static bool IsLinearMapping(int side)
    {
        return false;
    }

    REAL fCord;
    int  fFourDigits;
    REAL fAngle;
    REAL fX0[3];
    REAL fP;
    REAL fM;
    REAL fTT;
};
};

#endif
