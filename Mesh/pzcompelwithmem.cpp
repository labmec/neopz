/**
 * @file
 * @brief Contains the declaration of the TPZCompElWithMem class, it is as TPZCompEl with enable material memory feature.
 */

#include "pzcompelwithmem.h"
#include "pzshapepoint.h"
#include "pzshapelinear.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"
#include "pzshapecube.h"
#include "pzshapetetra.h"
#include "pzshapeprism.h"
#include "pzshapepiram.h"

template<class TBASE>
TPZCompElWithMem<TBASE>::~TPZCompElWithMem() {
	SetFreeIntPtIndices();
}

template <class TBASE>
inline void TPZCompElWithMem<TBASE>::PrepareIntPtIndices() {
    
    // This code was copied from TPZInterpolationSpace::CalcStiff with minor changes
    // regarding integration point index evaluation.
    // Inclusions were commented properly.
    
    if(fIntPtIndices.NElements())
    {
#ifdef LOG4CXX
        {
            std::stringstream sout;
            sout << __PRETTY_FUNCTION__ << " Attempting to add memory indices to an already configured TPZCompElWithMem";
            LOGPZ_ERROR(CompElWMemlogger,sout.str().c_str());
        }
#endif
        return;
    }
    
    TPZMaterial * material = TBASE::Material();
    if(!material){
        PZError << "Error at " << __PRETTY_FUNCTION__ << " this->Material() == NULL\n";
        return;
    }
    
    if (this->NumberOfCompElementsInsideThisCompEl() == 0) {
        // This is suposed to happen if in the constructor of a multiphysics element. The CompEl vector is only initialized after the autobuild
        return;
    }
    
    const TPZIntPoints &intrule = TBASE::GetIntegrationRule();
    
    
    int intrulepoints = intrule.NPoints();
    
    fIntPtIndices.Resize(intrulepoints);
    
    for(int int_ind = 0; int_ind < intrulepoints; ++int_ind){
        fIntPtIndices[int_ind] = this->Material()->PushMemItem();
        // Pushing a new entry in the material memory
    } //Loop over integratin points generating a reference vector of memory
    //entries in the related pzmatwithmem for further use.
    
}

#include "TPZInterfaceEl.h"

template <>
inline void TPZCompElWithMem<TPZInterfaceElement>::PrepareIntPtIndices() {
    
    // This code was copied from TPZInterpolationSpace::CalcStiff with minor changes
    // regarding integration point index evaluation.
    // Inclusions were commented properly.
    
    if(fIntPtIndices.NElements())
    {
#ifdef LOG4CXX
        {
            std::stringstream sout;
            sout << __PRETTY_FUNCTION__ << " Attempting to add memory indices to an already configured TPZCompElWithMem";
            LOGPZ_ERROR(CompElWMemlogger,sout.str().c_str());
        }
#endif
        return;
    }
    
    TPZMaterial * material = TPZInterfaceElement::Material();
    if(!material){
        PZError << "Error at " << __PRETTY_FUNCTION__ << " this->Material() == NULL\n";
        return;
    }
    if (this->LeftElementSide().Element() == 0 || this->RightElementSide().Element() == 0) {
        return;
    }
    
    if (this->NumberOfCompElementsInsideThisCompEl() == 0) {
        // This is suposed to happen if in the constructor of a multiphysics element. The CompEl vector is only initialized after the autobuild
        return;
    }
    
    const TPZIntPoints &intrule = TPZInterfaceElement::GetIntegrationRule();
    
    
    int intrulepoints = intrule.NPoints();
    
    fIntPtIndices.Resize(intrulepoints);
    
    for(int int_ind = 0; int_ind < intrulepoints; ++int_ind){
        fIntPtIndices[int_ind] = this->Material()->PushMemItem();
        // Pushing a new entry in the material memory
    } //Loop over integratin points generating a reference vector of memory
    //entries in the related pzmatwithmem for further use.
    
}

template class TPZRestoreClass<TPZCompElWithMem<TPZIntelGen<pzshape::TPZShapePoint> >>;
template class TPZRestoreClass<TPZCompElWithMem<TPZIntelGen<pzshape::TPZShapeLinear> >>;
template class TPZRestoreClass<TPZCompElWithMem<TPZIntelGen<pzshape::TPZShapeTriang> >>;
template class TPZRestoreClass<TPZCompElWithMem<TPZIntelGen<pzshape::TPZShapeQuad> >>;
template class TPZRestoreClass<TPZCompElWithMem<TPZIntelGen<pzshape::TPZShapeCube> >>;
template class TPZRestoreClass<TPZCompElWithMem<TPZIntelGen<pzshape::TPZShapeTetra> >>;
template class TPZRestoreClass<TPZCompElWithMem<TPZIntelGen<pzshape::TPZShapePrism> >>;
template class TPZRestoreClass<TPZCompElWithMem<TPZIntelGen<pzshape::TPZShapePiram> >>;

#include "pzgeopoint.h"
#include "pzgeoquad.h"
#include "pzgeotriangle.h"
#include "pzgeotetrahedra.h"
#include "pzgeoprism.h"
#include "TPZGeoCube.h"
#include "TPZGeoLinear.h"
#include "pzgeopyramid.h"
#include "TPZMultiphysicsInterfaceEl.h"


template class TPZRestoreClass<TPZCompElWithMem<TPZMultiphysicsCompEl<pzgeom::TPZGeoPoint> >>;
template class TPZRestoreClass<TPZCompElWithMem<TPZMultiphysicsCompEl<pzgeom::TPZGeoLinear> >>;
template class TPZRestoreClass<TPZCompElWithMem<TPZMultiphysicsCompEl<pzgeom::TPZGeoTriangle> >>;
template class TPZRestoreClass<TPZCompElWithMem<TPZMultiphysicsCompEl<pzgeom::TPZGeoQuad> >>;
template class TPZRestoreClass<TPZCompElWithMem<TPZMultiphysicsCompEl<pzgeom::TPZGeoCube> >>;
template class TPZRestoreClass<TPZCompElWithMem<TPZMultiphysicsCompEl<pzgeom::TPZGeoTetrahedra> >>;
template class TPZRestoreClass<TPZCompElWithMem<TPZMultiphysicsCompEl<pzgeom::TPZGeoPrism> >>;
template class TPZRestoreClass<TPZCompElWithMem<TPZMultiphysicsCompEl<pzgeom::TPZGeoPyramid> >>;
template class TPZRestoreClass<TPZCompElWithMem<TPZMultiphysicsInterfaceElement>>;
template class TPZRestoreClass<TPZCompElWithMem<TPZInterfaceElement>>;


/*
template class TPZMultiphysicsCompEl<pzgeom::TPZGeoPoint>;
template class TPZMultiphysicsCompEl<pzgeom::TPZGeoLinear>;
template class TPZMultiphysicsCompEl<pzgeom::TPZGeoTriangle>;
template class TPZMultiphysicsCompEl<pzgeom::TPZGeoQuad>;
template class TPZMultiphysicsCompEl<pzgeom::TPZGeoCube>;
template class TPZMultiphysicsCompEl<pzgeom::TPZGeoPrism>;
template class TPZMultiphysicsCompEl<pzgeom::TPZGeoTetrahedra>;
template class TPZMultiphysicsCompEl<pzgeom::TPZGeoPyramid>;*/
