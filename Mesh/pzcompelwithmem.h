/**
 * @file
 * @brief Contains the declaration of the TPZCompElWithMem class, it is as TPZCompEl with enable material memory feature.
 */

#ifndef PZCOMPELWITHMEM_H
#define PZCOMPELWITHMEM_H

class TPZMaterialData;

#include "pzinterpolationspace.h"
#include "pzstack.h"
#include "pzcmesh.h"
#include "pzquad.h"
#include "pzmaterial.h"
#include "pzelctemp.h"
#include "pzmultiphysicscompel.h"

//#include "tpzpoint.h"

#include "pzlog.h"



#ifdef LOG4CXX
static LoggerPtr CompElWMemlogger(Logger::getLogger("mesh.TPZCompElWithMem"));
#endif

/**
 * @brief This class implements the TPZCompEl structure to enable material memory feature. \n
 * It should be instantiated using one of TPZCompEl bottom classes as parent in the template parameter.
 * @ingroup CompElement
 * @since July, 23 2008
 */

template <class TBASE>
class TPZCompElWithMem : public TBASE
{
public:
	
	TPZCompElWithMem();
	
	virtual ~TPZCompElWithMem();
	
	TPZCompElWithMem(TPZCompMesh &mesh, TPZGeoEl *gel, long &index);
    
    TPZCompElWithMem(TPZCompMesh &mesh, TPZGeoEl *ref, long &index, TPZCompElSide left, TPZCompElSide right);
	
	TPZCompElWithMem(TPZCompMesh &mesh, const TPZCompElWithMem<TBASE> &copy);
	
	/** @brief used to generate patch mesh... generates a map of connect index from global mesh to clone mesh */
	TPZCompElWithMem(TPZCompMesh &mesh,
					 const TPZCompElWithMem<TBASE> &copy,
					 std::map<long,long> & gl2lcConMap,
					 std::map<long,long> & gl2lcElMap);
	
	virtual TPZCompEl *Clone(TPZCompMesh &mesh) const {
		return new TPZCompElWithMem<TBASE> (mesh, *this);
	}
	
    /** @brief Set create function in TPZCompMesh to create elements of this type */
	virtual void SetCreateFunctions(TPZCompMesh *mesh){
		mesh->SetAllCreateFunctionsContinuousWithMem();
	}

	/**
	 * @brief Create a copy of the given element. The clone copy have the connect indexes mapped to the local clone connects by the given map
	 * @param mesh Patch clone mesh
	 * @param gl2lcConMap map the connects indexes from global element (original) to the local copy.
	 * @param gl2lcElMap map the indexes of the elements between the original element and the patch element
	 */
	virtual TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,std::map<long,long> & gl2lcConMap,std::map<long,long>&gl2lcElMap) const
	{
		return new TPZCompElWithMem<TBASE> (mesh, *this, gl2lcConMap, gl2lcElMap);
	}
	
	void ComputeRequiredData(TPZMaterialData &data, TPZVec<REAL> &qsi);
	
	virtual void ComputeRequiredData(TPZVec<TPZMaterialData> &datavec, TPZVec<REAL> &intpointtemp, TPZManVector<TPZTransform> &trvec);

  long GetGlobalIntegrationPointIndex(TPZMaterialData &data);
	
protected:
	
	/** @brief PrepareIntPtIndices initializes the material damage varibles memory in the proper material class. */
	virtual void PrepareIntPtIndices();
  
	/** @brief PrepareIntPtIndices initializes the material damage varibles memory in the proper material class. */
	virtual void ForcePrepareIntPtIndices();
	
	/** @brief Frees the material damage varibles memory in the proper material class. */
	virtual void SetFreeIntPtIndices();
	
	void CopyIntPtIndicesFrom(const TPZCompElWithMem<TBASE> & copy);
	
public:
	
    /** @brief Get the indices of the vector of element memory associated with the integration points */
    /**
     * Will return an empty vector if no memory is associated with the integration point
     * Is implemented in TPZCompElWithMem
     */
  void GetMemoryIndices(TPZVec<long> &indices) const;

    /// Modify the maximum order an integration rule can integrate
	void SetIntegrationRule(int ord);
    
	/** @brief Saves the element data to a stream */
	virtual void Write(TPZStream &buf, int withclassid);
	
	/** @brief Reads the element data from a stream */
	virtual void Read(TPZStream &buf, void *context);
	
	/** @brief ClassId of the class. Is implemented for each type of compel in this .h */
	virtual int ClassId() const
	{
		return -1;
	}
	/**
	 * @name Print
	 * @brief Methods for print data structure
	 * @{
	 */
	
	/**
	 * @brief Prints element data
	 * @param out indicates the device where the data will be printed
	 */
	virtual void Print(std::ostream & out = std::cout) const
	{
		TBASE::Print(out);
		out << "Integration point indexes " << fIntPtIndices << std::endl;
	}
	
	/** @} */
	
private:
	
	TPZStack<long,128> fIntPtIndices;
	
};

template<class TBASE>
TPZCompElWithMem<TBASE>::TPZCompElWithMem() : TBASE() {
	//PrepareIntPtIndices();
}

template<class TBASE>
TPZCompElWithMem<TBASE>::TPZCompElWithMem(TPZCompMesh &mesh, TPZGeoEl *gel, long &index) :
TBASE(mesh, gel, index){
	PrepareIntPtIndices();
}

template<class TBASE>
TPZCompElWithMem<TBASE>::TPZCompElWithMem(TPZCompMesh &mesh, TPZGeoEl *ref, long &index, TPZCompElSide left, TPZCompElSide right) :
TBASE(mesh, ref, index, left, right){
    PrepareIntPtIndices();
}

template<class TBASE>
TPZCompElWithMem<TBASE>::TPZCompElWithMem(TPZCompMesh &mesh, const TPZCompElWithMem<TBASE> &copy) :
TBASE(mesh, copy) {
	CopyIntPtIndicesFrom(copy);
}

 
template<class TBASE>
TPZCompElWithMem<TBASE>::TPZCompElWithMem(TPZCompMesh &mesh,
										  const TPZCompElWithMem<TBASE> &copy,
										  std::map<long,long> & gl2lcConMap,
										  std::map<long,long> & gl2lcElMap) :
TBASE(mesh,copy,gl2lcConMap,gl2lcElMap)
{
	CopyIntPtIndicesFrom(copy);
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

template <class TBASE>
inline void TPZCompElWithMem<TBASE>::ForcePrepareIntPtIndices() {
	
//	if(fIntPtIndices.NElements())
//	{
//#ifdef LOG4CXX
//		{
//			std::stringstream sout;
//			sout << __PRETTY_FUNCTION__ << " Attempting to add memory indices to an already configured TPZCompElWithMem";
//			LOGPZ_ERROR(CompElWMemlogger,sout.str().c_str());
//		}
//#endif
//		return;
//	}
	
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



template <class TBASE>
inline void TPZCompElWithMem<TBASE>::SetFreeIntPtIndices() {
	
	TPZMaterial * material = TBASE::Material();
	if(!material){
		PZError << "Error at " << __PRETTY_FUNCTION__ << " this->Material() == NULL\n";
		return;
	}
	
	long n = fIntPtIndices.NElements();
	
	for(long i = 0; i < n; i++){
		this->Material()->FreeMemItem(fIntPtIndices[i]);
	}
	
	fIntPtIndices.Resize(0);
	
}

template <class TBASE>
inline void TPZCompElWithMem<TBASE>::SetIntegrationRule(int ord)
{
    TBASE::SetIntegrationRule(ord);
    // verify if the number of integration points changed
    const TPZIntPoints &intrule = TBASE::GetIntegrationRule();
	int intrulepoints = intrule.NPoints();
    if (intrulepoints != fIntPtIndices.size()) {
        SetFreeIntPtIndices();
        PrepareIntPtIndices();
    }

}


/** @brief Get the indices of the vector of element memory associated with the integration points */
/**
 * Will return an empty vector if no memory is associated with the integration point
 * Is implemented in TPZCompElWithMem
 */
template<class TBASE>
inline void TPZCompElWithMem<TBASE>::GetMemoryIndices(TPZVec<long> &indices) const
{
    indices = fIntPtIndices;
}


template <class TBASE>
void TPZCompElWithMem<TBASE>::CopyIntPtIndicesFrom(const TPZCompElWithMem<TBASE> & copy)
{
	
	TPZMaterial * material = TBASE::Material();
  	if(!material){
    	PZError << "Error at " << __PRETTY_FUNCTION__ << " this->Material() == NULL\n";
    	return;
  	}
	
	long i, n = copy.fIntPtIndices.NElements();
	fIntPtIndices.Resize(n);
	
	for(i = 0; i < n; i++)
    {
        fIntPtIndices[i] = copy.fIntPtIndices[i];
    }
}

/** Save the element data to a stream */
template <class TBASE>
inline void TPZCompElWithMem<TBASE>::ComputeRequiredData(TPZMaterialData &data,
														 TPZVec<REAL> &qsi){
	TBASE::ComputeRequiredData(data, qsi);
    data.intGlobPtIndex = GetGlobalIntegrationPointIndex(data);
	//material index for the n-th CompEl integration point
}

template <class TBASE>
inline void TPZCompElWithMem<TBASE>::ComputeRequiredData(TPZVec<TPZMaterialData> &datavec, TPZVec<REAL> &intpointtemp, TPZManVector<TPZTransform> &trvec)
{
	TBASE::ComputeRequiredData(datavec,intpointtemp,trvec);
	
	int nelofthismphysics = this->NumberOfCompElementsInsideThisCompEl();
	for (int icel = 0; icel < nelofthismphysics; icel++) {
		datavec[icel].intGlobPtIndex = GetGlobalIntegrationPointIndex(datavec[icel]);
	}
}

template <class TBASE>
inline long TPZCompElWithMem<TBASE>::GetGlobalIntegrationPointIndex(TPZMaterialData &data)
{
    long glIntegralPt = -1;
    if (data.intLocPtIndex >= 0) {
        glIntegralPt = fIntPtIndices[ data.intLocPtIndex ]; // returning the
    }
    return glIntegralPt;
}

template <class TBASE>
inline void TPZCompElWithMem<TBASE>::Write(TPZStream &buf, int withclassid)
{
	TBASE::Write(buf,withclassid);
	TPZSaveable::WriteObjects(buf, fIntPtIndices);
    int classid = ClassId();
    buf.Write(&classid);
}

/** Read the element data from a stream */
template <class TBASE>
inline void TPZCompElWithMem<TBASE>::Read(TPZStream &buf, void *context)
{
	TBASE::Read(buf,context);
	TPZSaveable::ReadObjects(buf, fIntPtIndices);
    int classid;
    buf.Read(&classid);
    if (classid != ClassId()) {
        DebugStop();
    }
}

#include "pzshapepoint.h"

/** @brief Routines to send and receive messages */
template<>
inline int TPZCompElWithMem<TPZIntelGen<pzshape::TPZShapePoint> >::ClassId() const
{
	return TPZCOMPELWITHMEMPOINTID;
}

#include "pzshapelinear.h"

template<>
inline int TPZCompElWithMem<TPZIntelGen<pzshape::TPZShapeLinear> >::ClassId() const
{
	return TPZCOMPELWITHMEMLINEARID;
}

#include "pzshapetriang.h"

template<>
inline int TPZCompElWithMem<TPZIntelGen<pzshape::TPZShapeTriang> >::ClassId() const
{
	return TPZCOMPELWITHMEMTRIANGID;
}

#include "pzshapequad.h"

template<>
inline int TPZCompElWithMem<TPZIntelGen<pzshape::TPZShapeQuad> >::ClassId() const
{
	return TPZCOMPELWITHMEMQUADID;
}

#include "pzshapecube.h"

template<>
inline int TPZCompElWithMem<TPZIntelGen<pzshape::TPZShapeCube> >::ClassId() const
{
	return TPZCOMPELWITHMEMCUBEID;
}

#include "pzshapetetra.h"

template<>
inline int TPZCompElWithMem<TPZIntelGen<pzshape::TPZShapeTetra> >::ClassId() const
{
	return TPZCOMPELWITHMEMTETRAID;
}

#include "pzshapeprism.h"

template<>
inline int TPZCompElWithMem<TPZIntelGen<pzshape::TPZShapePrism> >::ClassId() const
{
	return TPZCOMPELWITHMEMPRISMID;
}

#include "pzshapepiram.h"

template<>
inline int TPZCompElWithMem<TPZIntelGen<pzshape::TPZShapePiram> >::ClassId() const
{
	return TPZCOMPELWITHMEMPIRAMID;
}

#endif
