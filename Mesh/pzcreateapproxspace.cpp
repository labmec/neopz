/**
 * @file
 * @brief Contains the implementation of the functions to creates computational elements.
 */

#include "pzcreateapproxspace.h"
#include "pzgeoel.h"
#include "pzcmesh.h"
#include "pzcondensedcompel.h"
#include "pzinterpolationspace.h" 

#include "pzshapecube.h"
#include "pzshapelinear.h"
#include "pzshapepiram.h"
#include "pzshapepoint.h"
#include "pzshapeprism.h"
#include "pzshapequad.h"
#include "pzshapetetra.h"
#include "pzshapetriang.h"

#include "pzgeopoint.h"
#include "pzgeoprism.h"
#include "pzgeopyramid.h"
#include "pzgeoquad.h"
#include "pzgeotetrahedra.h"
#include "pzgeotriangle.h"
#include "TPZGeoCube.h"
#include "TPZGeoLinear.h"

#include "TPZInterfaceEl.h"

#include "pzcompelwithmem.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.tpzcreateapproximationspace"));
#endif

/** @brief Creates computational point element */
TPZCompEl *CreatePointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational linear element */
TPZCompEl *CreateLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational quadrilateral element */
TPZCompEl *CreateQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational triangular element */
TPZCompEl *CreateTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational cube element */
TPZCompEl *CreateCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational prismal element */
TPZCompEl *CreatePrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational pyramidal element */
TPZCompEl *CreatePyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational tetrahedral element */
TPZCompEl *CreateTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);


using namespace pzshape;

TPZCompEl *CreateNoElement(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
#ifdef LOG4CXX
    if (logger->isWarnEnabled()) {
        std::stringstream sout;
        sout << "Mesh dimension " << mesh.Dimension() << " gel dimension " << gel->Dimension() << " will not create a computational element\n";
        LOGPZ_WARN(logger, sout.str())
    }
#endif
    index = -1;
	return NULL;
}


TPZCompEl *CreatePointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
    {
		return new TPZIntelGen<TPZShapePoint>(mesh,gel,index);
    }
    index = -1;
	return NULL;
}
TPZCompEl *CreateLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
    {
		TPZCompEl *result = new TPZIntelGen<TPZShapeLinear>(mesh,gel,index);
        return result;//new TPZCondensedCompel(result);
    }
    index = -1;
	return NULL;
}
TPZCompEl *CreateQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
    {
		return new TPZIntelGen<TPZShapeQuad>(mesh,gel,index);
    }
    index = -1;
	return NULL;
}

TPZCompEl *CreateTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZIntelGen<TPZShapeTriang>(mesh,gel,index);
    index = -1;
	return NULL;
}
TPZCompEl *CreateCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZIntelGen<TPZShapeCube>(mesh,gel,index);
    index = -1;
	return NULL;
}
TPZCompEl *CreatePrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZIntelGen<TPZShapePrism>(mesh,gel,index);
    index = -1;
	return NULL;
}
TPZCompEl *CreatePyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZIntelGen<TPZShapePiram>(mesh,gel,index);
    index = -1;
	return NULL;
}
TPZCompEl *CreateTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZIntelGen<TPZShapeTetra>(mesh,gel,index);
    index = -1;
	return NULL;
}


// with mem
TPZCompEl *CreatePointElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElWithMem <TPZIntelGen<TPZShapePoint> >(mesh,gel,index) ;
    index = -1;
	return NULL;
}
TPZCompEl *CreateLinearElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElWithMem <TPZIntelGen<TPZShapeLinear> >(mesh,gel,index);
    index = -1;
	return NULL;
}
TPZCompEl *CreateQuadElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElWithMem <TPZIntelGen<TPZShapeQuad> >(mesh,gel,index);
    index = -1;
	return NULL;
}
TPZCompEl *CreateTriangleElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElWithMem < TPZIntelGen<TPZShapeTriang> >(mesh,gel,index);
    index = -1;
	return NULL;
}
TPZCompEl *CreateCubeElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElWithMem <TPZIntelGen<TPZShapeCube> >(mesh,gel,index);
    index = -1;
	return NULL;
}
TPZCompEl *CreatePrismElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElWithMem <TPZIntelGen<TPZShapePrism> >(mesh,gel,index);
    index = -1;
	return NULL;
}
TPZCompEl *CreatePyramElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElWithMem <TPZIntelGen<TPZShapePiram> >(mesh,gel,index);
    index = -1;
	return NULL;
}
TPZCompEl *CreateTetraElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElWithMem <TPZIntelGen<TPZShapeTetra> >(mesh,gel,index);
    index = -1;
	return NULL;
}

/**
 * @brief Creates the computational elements, and the degree of freedom nodes
 */ 
/** Only element of material id in the set<int> will be created */
void TPZCreateApproximationSpace::BuildMesh(TPZCompMesh &cmesh, const std::set<int> &MaterialIDs) const {
	TPZAdmChunkVector<TPZGeoEl *> &elvec = cmesh.Reference()->ElementVec();
	int64_t i, nelem = elvec.NElements();
	int64_t neltocreate = 0;
	int64_t index;
	for(i=0; i<nelem; i++) {
		TPZGeoEl *gel = elvec[i];
		if(!gel) continue;
        if (gel->Reference()) {
            continue;
        }
		if(!gel->HasSubElement()) {
			neltocreate++;
		}
	}
	std::set<int> matnotfound;
	int64_t nbl = cmesh.Block().NBlocks();
    if(neltocreate > nbl) 
    {
        cmesh.Block().SetNBlocks(neltocreate);
    }
    cmesh.Block().SetNBlocks(nbl);
	
	for(i=0; i<nelem; i++) {
		TPZGeoEl *gel = elvec[i];
		if(!gel || gel->Reference()) continue;
		if(!gel->HasSubElement()) {
			int matid = gel->MaterialId();
			TPZMaterial * mat = cmesh.FindMaterial(matid);
			if(!mat)
			{
				matnotfound.insert(matid);
				continue;
			}
			int printing = 0;
			if (printing) {
				gel->Print(std::cout);
			}
			
			//checking material in MaterialIDs
            std::set<int>::const_iterator found = MaterialIDs.find(matid);
            if (found == MaterialIDs.end())
            {
                continue;
            }
			
			if(!gel->Reference() && gel->NumInterfaces() == 0)
			{
				CreateCompEl(gel,cmesh,index);
                if (fCreateHybridMesh) {
                    cmesh.ElementVec()[index]->Reference()->ResetReference();
                }
#ifdef LOG4CXX
                if (logger->isDebugEnabled())
                {

                    std::stringstream sout;
                    if (index < 0) {
                        if(gel->Dimension() == 0){
                            sout << "Zero dimensional element, is your approximation space Hdiv?. " << std::endl;
                            gel->Print(sout);
                            sout << "No computational element was created. " << std::endl;
                        }
                    }else{
                        cmesh.ElementVec()[index]->Print(sout);
                    }
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
			}
		}
	}
	cmesh.InitializeBlock();
}

/** @brief Creates the computational elements, and the degree of freedom nodes */
void TPZCreateApproximationSpace::BuildMesh(TPZCompMesh &cmesh) const {
    std::set<int> materialids;
    std::map<int,TPZMaterial *>::iterator it;
    for (it = cmesh.MaterialVec().begin(); it != cmesh.MaterialVec().end(); it++) {
        materialids.insert(it->first);
    }
    BuildMesh(cmesh,materialids);
}

void TPZCreateApproximationSpace::CreateInterfaces(TPZCompMesh &cmesh, const std::set<int> &MaterialIDs){
    for(int64_t el = 0; el < cmesh.ElementVec().NElements(); el++)
	{
		TPZCompEl * compEl = cmesh.ElementVec()[el];
		if(!compEl) continue;
        int matid = compEl->Reference()->MaterialId();
        
        //checking if matid is in MaterialIDs
        std::set<int>::iterator it;
        it = MaterialIDs.find(matid);
        if(it!=MaterialIDs.end())
        {
            TPZInterpolationSpace * interpEl = dynamic_cast<TPZInterpolationSpace *>(compEl);
            if(!interpEl) continue;
            interpEl->CreateInterfaces();
        }
    }
    
}

void TPZCreateApproximationSpace::CreateInterfaces(TPZCompMesh &cmesh){
    for(int64_t el = 0; el < cmesh.ElementVec().NElements(); el++)
	{
		TPZCompEl * compEl = cmesh.ElementVec()[el];
		if(!compEl) continue;
        
        TPZInterpolationSpace * interpEl = dynamic_cast<TPZInterpolationSpace *>(compEl);
        if(!interpEl) continue;
        interpEl->CreateInterfaces();
    }
    
}

void TPZCreateApproximationSpace::SetAllCreateFunctions(TPZCompEl &cel, TPZCompMesh *mesh){
	cel.SetCreateFunctions(mesh);
}



void TPZCreateApproximationSpace::SetAllCreateFunctionsDiscontinuous(){
	
    fp[EPoint] = TPZCompElDisc::CreateDisc;
    fp[EOned] = TPZCompElDisc::CreateDisc;
    fp[ETriangle] = TPZCompElDisc::CreateDisc;
    fp[EQuadrilateral] = TPZCompElDisc::CreateDisc;
    fp[ETetraedro] = TPZCompElDisc::CreateDisc;
    fp[EPiramide] = TPZCompElDisc::CreateDisc;
    fp[EPrisma] = TPZCompElDisc::CreateDisc;
    fp[ECube] = TPZCompElDisc::CreateDisc;
    
    /*
     pzgeom::TPZGeoPoint::fp = TPZCompElDisc::CreateDisc;
     pzgeom::TPZGeoLinear::fp = TPZCompElDisc::CreateDisc;
     pzgeom::TPZGeoQuad::fp = TPZCompElDisc::CreateDisc;
     pzgeom::TPZGeoTriangle::fp = TPZCompElDisc::CreateDisc;
     pzgeom::TPZGeoPrism::fp = TPZCompElDisc::CreateDisc;
     pzgeom::TPZGeoTetrahedra::fp = TPZCompElDisc::CreateDisc;
     pzgeom::TPZGeoPyramid::fp = TPZCompElDisc::CreateDisc;
     pzgeom::TPZGeoCube::fp = TPZCompElDisc::CreateDisc;
     */
}


void TPZCreateApproximationSpace::SetAllCreateFunctionsContinuous(){
    fp[EPoint] = CreatePointEl;
    fp[EOned] = CreateLinearEl;
    fp[EQuadrilateral] = CreateQuadEl;
    fp[ETriangle] = CreateTriangleEl;
    fp[ETetraedro] = CreateTetraEl;
    fp[EPiramide] = CreatePyramEl;
    fp[EPrisma] = CreatePrismEl;
    fp[ECube] = CreateCubeEl;
    
    /*
     pzgeom::TPZGeoPoint::fp =  CreatePointEl;
     pzgeom::TPZGeoLinear::fp =  CreateLinearEl;
     pzgeom::TPZGeoQuad::fp = CreateQuadEl;
     pzgeom::TPZGeoTriangle::fp =  CreateTriangleEl;
     pzgeom::TPZGeoPrism::fp = CreatePrismEl;
     pzgeom::TPZGeoTetrahedra::fp = CreateTetraEl;
     pzgeom::TPZGeoPyramid::fp = CreatePyramEl;
     pzgeom::TPZGeoCube::fp = CreateCubeEl;
     */
	
}

void TPZCreateApproximationSpace::SetAllCreateFunctionsContinuousWithMem()
{
	// These ones will be always continuous for viscoelasticity
	fp[EPoint] = CreatePointElWithMem;
	fp[EOned] = CreateLinearElWithMem;
	fp[ETriangle] = CreateTriangleElWithMem;
	fp[EQuadrilateral] = CreateQuadElWithMem;	
    fp[ETetraedro] = CreateTetraElWithMem;
    fp[EPiramide] = CreatePyramElWithMem;
    fp[EPrisma] = CreatePrismElWithMem;
    fp[ECube] = CreateCubeElWithMem;
    
}


#include "pzelchdiv.h"
#include "pzelchdivbound2.h"

void TPZCreateApproximationSpace::SetAllCreateFunctionsHDiv(int dimension){
	
    switch (dimension) {
        case 1:
            fp[EPoint] = CreateHDivBoundPointEl;
            fp[EOned] = CreateHDivLinearEl;
            fp[ETriangle] = CreateHDivTriangleEl;
            fp[EQuadrilateral] = CreateNoElement;
            fp[ETetraedro] = CreateNoElement;
            fp[EPiramide] = CreateNoElement;
            fp[EPrisma] = CreateNoElement;
            fp[ECube] = CreateNoElement;
            break;
        case 2:
            fp[EPoint] = CreateNoElement;
            fp[EOned] = CreateHDivBoundLinearEl;
            fp[ETriangle] = CreateHDivTriangleEl;
            fp[EQuadrilateral] = CreateHDivQuadEl;
            fp[ETetraedro] = CreateNoElement;
            fp[EPiramide] = CreateNoElement;
            fp[EPrisma] = CreateNoElement;
            fp[ECube] = CreateNoElement;
            break;
        case 3:
            fp[EPoint] = CreateNoElement;
            fp[EOned] = CreateNoElement;
            fp[ETriangle] = CreateHDivBoundTriangleEl;
            fp[EQuadrilateral] = CreateHDivBoundQuadEl;
            fp[ETetraedro] = CreateHDivTetraEl;
            fp[EPiramide] = CreateHDivPyramEl;
            fp[EPrisma] = CreateHDivPrismEl;
            fp[ECube] = CreateHDivCubeEl;
            break;
        default:
            DebugStop();
            break;
    }
    
    /*
     pzgeom::TPZGeoPoint::fp =  CreateHDivPointEl;
     pzgeom::TPZGeoLinear::fp =  CreateHDivLinearEl;
     pzgeom::TPZGeoQuad::fp = CreateHDivQuadEl;
     pzgeom::TPZGeoTriangle::fp =  CreateHDivTriangleEl;
     pzgeom::TPZGeoPrism::fp = CreateHDivPrismEl;
     pzgeom::TPZGeoTetrahedra::fp = CreateHDivTetraEl;
     pzgeom::TPZGeoPyramid::fp = CreateHDivPyramEl;
     pzgeom::TPZGeoCube::fp = CreateHDivCubeEl;
     */
}

/** @brief Create an approximation space with HDiv elements */
void TPZCreateApproximationSpace::SetAllCreateFunctionsHDivReferred(int dimension)
{
    switch (dimension) {
        case 1:
            fp[EPoint] = CreateRefHDivBoundPointEl;
            fp[EOned] = CreateRefHDivLinearEl;
            fp[ETriangle] = CreateRefHDivTriangleEl;
            fp[EQuadrilateral] = CreateNoElement;
            fp[ETetraedro] = CreateNoElement;
            fp[EPiramide] = CreateNoElement;
            fp[EPrisma] = CreateNoElement;
            fp[ECube] = CreateNoElement;
            break;
        case 2:
            fp[EPoint] = CreateNoElement;
            fp[EOned] = CreateRefHDivBoundLinearEl;
            fp[ETriangle] = CreateRefHDivTriangleEl;
            fp[EQuadrilateral] = CreateRefHDivQuadEl;
            fp[ETetraedro] = CreateNoElement;
            fp[EPiramide] = CreateNoElement;
            fp[EPrisma] = CreateNoElement;
            fp[ECube] = CreateNoElement;
            break;
        case 3:
            fp[EPoint] = CreateNoElement;
            fp[EOned] = CreateNoElement;
            fp[ETriangle] = CreateRefHDivBoundTriangleEl;
            fp[EQuadrilateral] = CreateRefHDivBoundQuadEl;
            fp[ETetraedro] = CreateRefHDivTetraEl;
            fp[EPiramide] = CreateRefHDivPyramEl;
            fp[EPrisma] = CreateRefHDivPrismEl;
            fp[ECube] = CreateRefHDivCubeEl;
            break;
        default:
            DebugStop();
            break;
    }

}

#if defined(USING_MKL) && defined(USING_LAPACK) && !defined(STATE_COMPLEX)

#include "TPZSBFemVolume.h"

void TPZCreateApproximationSpace::SetAllCreateFunctionsSBFem(int dimension){
    
    switch (dimension) {
        case 1:
            DebugStop();
            break;
        case 2:
            fp[EPoint] = CreatePointEl;
            fp[EOned] = CreateSBFemCompEl;
            fp[ETriangle] = CreateNoElement;
            fp[EQuadrilateral] = CreateSBFemCompEl;
            fp[ETetraedro] = CreateNoElement;
            fp[EPiramide] = CreateNoElement;
            fp[EPrisma] = CreateNoElement;
            fp[ECube] = CreateNoElement;
            break;
        case 3:
            fp[EPoint] = CreatePointEl;
            fp[EOned] = CreateSBFemCompEl;
            fp[ETriangle] = CreateTriangleEl;
            fp[EQuadrilateral] = CreateQuadEl;
            fp[ETetraedro] = CreateNoElement;
            fp[EPiramide] = CreateNoElement;
            fp[EPrisma] = CreateSBFemCompEl;
            fp[ECube] = CreateSBFemCompEl;
            break;
        default:
            DebugStop();
            break;
    }
}

#endif

#include "pzhdivfull.h"

void TPZCreateApproximationSpace::SetAllCreateFunctionsHDivFull(int dimension){
	
    switch (dimension) {
        case 1:
            fp[EPoint] = CreateHDivBoundPointEl;
            fp[EOned] = CreateHDivFullLinearEl;
            fp[ETriangle] = CreateNoElement;
            fp[EQuadrilateral] = CreateNoElement;
            fp[ETetraedro] = CreateNoElement;
            fp[EPiramide] = CreateNoElement;
            fp[EPrisma] = CreateNoElement;
            fp[ECube] = CreateNoElement;
            break;
        case 2:
            fp[EPoint] = CreateNoElement;
            fp[EOned] = CreateHDivBoundLinearEl;
            fp[ETriangle] = CreateHDivFullTriangleEl;
            fp[EQuadrilateral] = CreateHDivFullQuadEl;
            fp[ETetraedro] = CreateNoElement;
            fp[EPiramide] = CreateNoElement;
            fp[EPrisma] = CreateNoElement;
            fp[ECube] = CreateNoElement;
            break;
        case 3:
            fp[EPoint] = CreateNoElement;
            fp[EOned] = CreateNoElement;
            fp[ETriangle] = CreateHDivBoundTriangleEl;
            fp[EQuadrilateral] = CreateHDivBoundQuadEl;
            fp[ETetraedro] = CreateHDivFullTetraEl;
            fp[EPiramide] = CreateHDivFullPyramEl;
            fp[EPrisma] = CreateHDivFullPrismEl;
            fp[ECube] = CreateHDivFullCubeEl;
            break;
        default:
            DebugStop();
            break;
    }
    
}



#ifndef STATE_COMPLEX
#include "pzhdivpressure.h"

void TPZCreateApproximationSpace::SetAllCreateFunctionsHDivPressure(int dimension){
    switch (dimension) {
        case 1:
            fp[EPoint] = CreateHDivBoundPointEl;
            fp[EOned] = CreateHDivPressureLinearEl;
            fp[ETriangle] = CreateNoElement;
            fp[EQuadrilateral] = CreateNoElement;
            fp[ETetraedro] = CreateNoElement;
            fp[EPiramide] = CreateNoElement;
            fp[EPrisma] = CreateNoElement;
            fp[ECube] = CreateNoElement;
            break;
        case 2:
            fp[EPoint] = CreateNoElement;
            fp[EOned] = CreateHDivBoundLinearEl;
            fp[ETriangle] = CreateHDivPressureTriangleEl;
            fp[EQuadrilateral] = CreateHDivPressureQuadEl;
            fp[ETetraedro] = CreateNoElement;
            fp[EPiramide] = CreateNoElement;
            fp[EPrisma] = CreateNoElement;
            fp[ECube] = CreateNoElement;
            break;
        case 3:
            fp[EPoint] = CreateNoElement;
            fp[EOned] = CreateNoElement;
            fp[ETriangle] = CreateHDivBoundTriangleEl;
            fp[EQuadrilateral] = CreateHDivBoundQuadEl;
            fp[ETetraedro] = CreateHDivPressureTetraEl;
            fp[EPiramide] = CreateHDivPressurePyramEl;
            fp[EPrisma] = CreateHDivPressurePrismEl;
            fp[ECube] = CreateHDivPressureCubeEl;
            break;
        default:
            DebugStop();
            break;
    }
    
}
#endif


#include "pzreferredcompel.h"
#include "pzelctemp.h"
void TPZCreateApproximationSpace::SetAllCreateFunctionsDiscontinuousReferred(){
	
    fp[EPoint] = CreateReferredDisc;
    fp[EOned] = CreateReferredDisc;
    fp[ETriangle] = CreateReferredDisc;
    fp[EQuadrilateral] = CreateReferredDisc;
    fp[ETetraedro] = CreateReferredDisc;
    fp[EPiramide] = CreateReferredDisc;
    fp[EPrisma] = CreateReferredDisc;
    fp[ECube] = CreateReferredDisc;
    
    /*
     pzgeom::TPZGeoPoint::fp =  CreateReferredDisc;
     pzgeom::TPZGeoLinear::fp =  CreateReferredDisc;
     pzgeom::TPZGeoQuad::fp = CreateReferredDisc;
     pzgeom::TPZGeoTriangle::fp =  CreateReferredDisc;
     pzgeom::TPZGeoPrism::fp = CreateReferredDisc;
     pzgeom::TPZGeoTetrahedra::fp = CreateReferredDisc;
     pzgeom::TPZGeoPyramid::fp = CreateReferredDisc;
     pzgeom::TPZGeoCube::fp = CreateReferredDisc;
     */
}

void TPZCreateApproximationSpace::SetAllCreateFunctionsContinuousReferred(){
	
    fp[EPoint] = CreateReferredPointEl;
    fp[EOned] = CreateReferredLinearEl;
    fp[ETriangle] = CreateReferredTriangleEl;
    fp[EQuadrilateral] = CreateReferredQuadEl;
    fp[ETetraedro] = CreateReferredTetraEl;
    fp[EPiramide] = CreateReferredPyramEl;
    fp[EPrisma] = CreateReferredPrismEl;
    fp[ECube] = CreateReferredCubeEl;
    
    /*
     pzgeom::TPZGeoPoint::fp =  CreateReferredPointEl;
     pzgeom::TPZGeoLinear::fp =  CreateReferredLinearEl;
     pzgeom::TPZGeoQuad::fp = CreateReferredQuadEl;
     pzgeom::TPZGeoTriangle::fp =  CreateReferredTriangleEl;
     pzgeom::TPZGeoPrism::fp = CreateReferredPrismEl;
     pzgeom::TPZGeoTetrahedra::fp = CreateReferredTetraEl;
     pzgeom::TPZGeoPyramid::fp = CreateReferredPyramEl;
     pzgeom::TPZGeoCube::fp = CreateReferredCubeEl;
     */
	
}

#include "pzmultiphysicscompel.h"
void TPZCreateApproximationSpace::SetAllCreateFunctionsMultiphysicElem(){
	
    fp[EPoint] = CreateMultiphysicsPointEl;
    fp[EOned] = CreateMultiphysicsLinearEl;
    fp[ETriangle] = CreateMultiphysicsTriangleEl;
    fp[EQuadrilateral] = CreateMultiphysicsQuadEl;
    fp[ETetraedro] = CreateMultiphysicsTetraEl;
    fp[EPiramide] = CreateMultiphysicsPyramEl;
    fp[EPrisma] = CreateMultiphysicsPrismEl;
    fp[ECube] = CreateMultiphysicsCubeEl;
    
    /*
     pzgeom::TPZGeoPoint::fp =  CreateMultiphysicsPointEl;
     pzgeom::TPZGeoLinear::fp =  CreateMultiphysicsLinearEl;
     pzgeom::TPZGeoTriangle::fp =  CreateMultiphysicsTriangleEl;
     pzgeom::TPZGeoQuad::fp = CreateMultiphysicsQuadEl;
     pzgeom::TPZGeoCube::fp = CreateMultiphysicsCubeEl;
     pzgeom::TPZGeoPrism::fp = CreateMultiphysicsPrismEl;
     pzgeom::TPZGeoTetrahedra::fp = CreateMultiphysicsTetraEl;
     pzgeom::TPZGeoPyramid::fp = CreateMultiphysicsPyramEl;
     */
}

void TPZCreateApproximationSpace::SetAllCreateFunctionsMultiphysicElemWithMem()
{
	
	fp[EPoint] = CreateMultiphysicsPointElWithMem;
	fp[EOned] = CreateMultiphysicsLinearElWithMem;
	fp[ETriangle] = CreateMultiphysicsTriangleElWithMem;
	fp[EQuadrilateral] = CreateMultiphysicsQuadElWithMem;
	fp[ETetraedro] = CreateMultiphysicsTetraElWithMem;
	fp[EPiramide] = CreateMultiphysicsPyramElWithMem;
	fp[EPrisma] = CreateMultiphysicsPrismElWithMem;
	fp[ECube] = CreateMultiphysicsCubeElWithMem;
	
}

/*
 * @brief Create a computational element using the function pointer for the topology
 */
TPZCompEl *TPZCreateApproximationSpace::CreateCompEl(TPZGeoEl *gel, TPZCompMesh &mesh, int64_t &index) const
{
    switch (gel->Type()) {
        case EPoint:
            return fp[EPoint](gel,mesh,index);
            break;
        case EOned:
            return fp[EOned](gel,mesh,index);
            break;
        case EQuadrilateral:
            return fp[EQuadrilateral](gel,mesh,index);
            break;
        case ETriangle:
            return fp[ETriangle](gel,mesh,index);
            break;
        case EPiramide:
            return fp[EPiramide](gel,mesh,index);
            break;
        case EPrisma:
            return fp[EPrisma](gel,mesh,index);
            break;
        case ETetraedro:
            return fp[ETetraedro](gel,mesh,index);
            break;
        case ECube:
            return fp[ECube](gel,mesh,index);
            break;
        default:
            DebugStop();
            break;
    }
    return 0;
}

/**
 * @brief Set custom function pointers
 */

void TPZCreateApproximationSpace::SetCreateFunctions(TPZVec<TCreateFunction> &createfuncs)
{
    fp[EPoint] = createfuncs[EPoint];
    fp[EOned] = createfuncs[EOned];
    fp[ETriangle] = createfuncs[ETriangle];
    fp[EQuadrilateral] = createfuncs[EQuadrilateral];
    fp[ETetraedro] = createfuncs[ETetraedro];
    fp[EPiramide] = createfuncs[EPiramide];
    fp[EPrisma] = createfuncs[EPrisma];
    fp[ECube] = createfuncs[ECube];
    
}

#include "pzcondensedcompel.h"

/**
 * @brief Encapsulate the elements in condensed computational elements
 */
void TPZCreateApproximationSpace::CondenseLocalEquations(TPZCompMesh &cmesh)
{
    int64_t nel = cmesh.NElements();
    int64_t iel;
    for (iel=0; iel<nel; iel++) {
        TPZCompEl *cel = cmesh.ElementVec()[iel];
        if(!cel) {
            continue;
        }
        new TPZCondensedCompEl(cel);
    }
    
}

/**
 * @brief Undo the encapsulate elements
 */
void TPZCreateApproximationSpace::UndoCondenseLocalEquations(TPZCompMesh &cmesh)
{
    int64_t nel = cmesh.NElements();
    int64_t iel;
    for (iel=0; iel<nel; iel++) {
        TPZCompEl *cel = cmesh.ElementVec()[iel];
        TPZCondensedCompEl *condel = dynamic_cast<TPZCondensedCompEl *>(cel);
        if(!condel) {
            continue;
        }
        condel->Unwrap();
    }    
}

/**
 * @brief transform in low order Raviar Tomas
 */
void TPZCreateApproximationSpace::MakeRaviartThomas(TPZCompMesh &cmesh)
{
    int64_t numcell = cmesh.NElements();
    int64_t el;
    for (el = 0; el<numcell ; el++) {
        TPZCompEl *cel = cmesh.ElementVec()[el];
        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
        if (!intel) {
            continue;
        }
        intel->SetPreferredOrder(1);
    }
    cmesh.ExpandSolution();
    for (el = 0; el<numcell ; el++) {
        TPZCompEl *cel = cmesh.ElementVec()[el];
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) {
            continue;
        }
        TPZGeoEl *gel = intel->Reference();
        int geldim = gel->Dimension();
        int is;
        int nsides = gel->NSides();
        for (is=0; is<nsides; is++) {
            if (gel->SideDimension(is) != geldim-1) {
                continue;
            }
            int nsconnects = intel->NSideConnects(is);
            // only interested in HDiv elements
            if (nsconnects != 1) {
                continue;
            }
            int64_t cindex = intel->SideConnectIndex(0, is);
            TPZConnect &c = intel->Connect(intel->SideConnectLocId(0,is));
            if (c.HasDependency()) {
                continue;
            }
            int nshape = 1;
            int nstate = 1;
            int order = 0;
            int64_t cindex2 = cmesh.AllocateNewConnect(nshape, nstate, order);
            //            TPZConnect &c2 = cmesh.ConnectVec()[cindex];
            TPZFNMatrix<2,REAL> depmat(2,1,1.);
            c.AddDependency(cindex, cindex2, depmat, 0, 0, 2, 1);
        }
    }
    cmesh.ExpandSolution();
}

/**
 * @brief transform in low order Raviar Tomas
 */
void TPZCreateApproximationSpace::UndoMakeRaviartThomas(TPZCompMesh &cmesh)
{
    int64_t numcell = cmesh.NElements();
    int64_t el;
    for (el = 0; el<numcell ; el++) {
        TPZCompEl *cel = cmesh.ElementVec()[el];
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) {
            continue;
        }
        TPZGeoEl *gel = intel->Reference();
        int geldim = gel->Dimension();
        int is;
        int nsides = gel->NSides();
        for (is=0; is<nsides; is++) {
            if (gel->SideDimension(is) != geldim-1) {
                continue;
            }
            int nsconnects = intel->NSideConnects(is);
            // only interested in HDiv elements
            if (nsconnects != 1) {
                continue;
            }
            //            int cindex = intel->SideConnectIndex(0, is);
            TPZConnect &c = intel->Connect(intel->SideConnectLocId(0,is));
            if (c.HasDependency()) {
                c.RemoveDepend();
            }
        }
    }
    cmesh.ExpandSolution();
    cmesh.CleanUpUnconnectedNodes();
}

/** @brief Create interface elements between the computational elements */
void TPZCreateApproximationSpace::CreateInterfaceElements(TPZCompMesh *mesh, bool betweencontinuous, bool multiphysics)
{
    TPZChunkVector<TPZCompEl *> compelvec = mesh->ElementVec();
    int64_t nel = compelvec.NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = compelvec[el];
        if (!cel) {
            continue;
        }
        if(!multiphysics)
        {
            TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
            if (intel) {
                intel->CreateInterfaces(betweencontinuous);
            }
        }
        else {
            TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
            if (!mfcel) {
                continue;
            }
            mfcel->CreateInterfaces();
        }
        
    }
}

/// this method will substitute all interface elements with materialid within the set by three elements : one H1 element and two interface elements
void TPZCreateApproximationSpace::Hybridize(TPZCompMesh &cmesh,const std::set<int> &matids, bool isconnectedElem)
{
    cmesh.ApproxSpace().SetAllCreateFunctionsContinuous();
    cmesh.Reference()->ResetReference();
    int nel = cmesh.NElements();
    for (int el=0; el<nel; el++) {
        TPZCompEl *cel = cmesh.ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZInterfaceElement *face = dynamic_cast<TPZInterfaceElement *>(cel);
        if (!face) {
            continue;
        }
        int matid = face->Material()->Id();
        if (matids.find(matid) == matids.end()) {
            continue;
        }
        TPZCompElSide left, right;
        left = face->LeftElementSide();
        right = face->RightElementSide();
        int leftmatid = left.Element()->Material()->Id();
        int rightmatid = right.Element()->Material()->Id();
        TPZGeoEl *gelface = face->Reference();
        gelface = gelface->CreateBCGeoEl(gelface->NSides()-1, matid);
        delete face;
        int64_t index;
        
        //hp mesh
        TPZInterpolationSpace *leftint = dynamic_cast<TPZInterpolationSpace *>(left.Element());
        TPZInterpolationSpace *rightint = dynamic_cast<TPZInterpolationSpace *>(right.Element());
        int leftorder = leftint->MaxOrder();
        int rightorder = rightint->MaxOrder();
        int neworder =  MIN(leftorder,rightorder);
        if(neworder<0) neworder = MAX(leftorder,rightorder);
        cmesh.SetDefaultOrder(neworder);
//        int neworder = MAX(leftorder,rightorder);
//        cmesh.SetDefaultOrder(neworder);
       
        
        cmesh.ApproxSpace().CreateCompEl(gelface, cmesh, index);
        TPZCompEl *newcel = cmesh.ElementVec()[index];
        if(!isconnectedElem){
            gelface->ResetReference();
        }
        TPZCompElSide center(newcel,gelface->NSides()-1);
        //        
        TPZGeoEl *leftgelface = gelface->CreateBCGeoEl(gelface->NSides()-1, leftmatid);
        TPZGeoEl *rightgelface = gelface->CreateBCGeoEl(gelface->NSides()-1, rightmatid);

        new TPZInterfaceElement(cmesh,leftgelface,index,left,center);
        new TPZInterfaceElement(cmesh,rightgelface,index,right,center);
        
//        TPZInterfaceElement *faceleft = new TPZInterfaceElement(cmesh,leftgelface,index,left,center);
//        TPZInterfaceElement *faceright = new TPZInterfaceElement(cmesh,rightgelface,index,right,center);
        
    }
    cmesh.CleanUpUnconnectedNodes();
    cmesh.ExpandSolution();
}

int TPZCreateApproximationSpace::ClassId() const {
    return Hash("TPZCreateApproximationSpace");
}

void TPZCreateApproximationSpace::Read(TPZStream& buf, void* context) { //ok
    buf.Read(fCreateHybridMesh);
    buf.Read(fCreateLagrangeMultiplier);
    buf.Read(fCreateWithMemory);
}

void TPZCreateApproximationSpace::Write(TPZStream& buf, int withclassid) const { //ok
    buf.Write(fCreateHybridMesh);
    buf.Write(fCreateLagrangeMultiplier);
    buf.Write(fCreateWithMemory);
}

template class TPZRestoreClass<TPZCreateApproximationSpace>;
