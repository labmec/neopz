/**
 * @file
 * @brief Contains the implementation of the TPZCompElHDivCollapsed methods.
 */


#include "TPZCompElHDivCollapsed.h"
#include "pzgeoel.h"
#include "TPZMaterial.h"
#include "pzlog.h"
#include "TPZShapeDisc.h"
#include "TPZCompElDisc.h"
#include "pzmaterialdata.h"
#include "pzelchdiv.h"


#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.TPZCompElHDivCollapsed"));
#endif

template<class TSHAPE>
TPZCompElHDivCollapsed<TSHAPE>::TPZCompElHDivCollapsed(TPZCompMesh &mesh, TPZGeoEl *gel, int64_t &index) :
TPZRegisterClassId(&TPZCompElHDivCollapsed::ClassId),
TPZCompElHDiv<TSHAPE>(mesh,gel,index), fBottom(mesh,gel,index), fTop(mesh,gel,index)
{
    index = this->fIndex;
    mesh.ElementVec().SetFree(fTop.Index());
    mesh.ElementVec().SetFree(fBottom.Index());
    int64_t bottom_c_index = fBottom.ConnectIndex(0);
    int64_t top_c_index = fTop.ConnectIndex(0);
    fBottom.SetIndex(-1);
    fTop.SetIndex(-1);
    

#ifdef LOG4CXX
    if (logger->isDebugEnabled())
	 {
         std::stringstream sout;
         sout << "Finalizando criacao do elemento ";
         this->Print(sout);
         LOGPZ_DEBUG(logger,sout.str())
	 }
#endif
	 
}

template<class TSHAPE>
TPZCompElHDivCollapsed<TSHAPE>::TPZCompElHDivCollapsed(TPZCompMesh &mesh, const TPZCompElHDivCollapsed<TSHAPE> &copy) :
TPZRegisterClassId(&TPZCompElHDivCollapsed::ClassId),
TPZCompElHDiv<TSHAPE>(mesh,copy),fBottom(copy.fBottom), fTop(copy.fTop)
{
}

// NAO TESTADO
template<class TSHAPE>
TPZCompElHDivCollapsed<TSHAPE>::TPZCompElHDivCollapsed(TPZCompMesh &mesh,
												 const TPZCompElHDivCollapsed<TSHAPE> &copy,
												 std::map<int64_t,int64_t> & gl2lcConMap,
												 std::map<int64_t,int64_t> & gl2lcElMap) :
TPZRegisterClassId(&TPZCompElHDivCollapsed::ClassId),
TPZCompElHDiv<TSHAPE>(mesh,copy,gl2lcConMap,gl2lcElMap),fBottom(mesh,copy.fBottom,gl2lcConMap,gl2lcElMap),
fTop(mesh,copy.fBottom,gl2lcConMap,gl2lcElMap)
{
	
	this-> fPreferredOrder = copy.fPreferredOrder;
	int i;
	for(i=0;i<TSHAPE::NFacets;i++)
	{
		int64_t lcIdx = -1;
		int64_t glIdx = copy.fConnectIndexes[i];
		if(glIdx == -1)
		{
			// nothing to clone
			this->fConnectIndexes[i] = -1;
			continue;
		}
		if (gl2lcConMap.find(glIdx) != gl2lcConMap.end())
		{
			lcIdx = gl2lcConMap[glIdx];
		}
		else
		{
			std::stringstream sout;
			sout << "ERROR in : " << __PRETTY_FUNCTION__
			<< " trying to clone the connect index: " << glIdx
			<< " wich is not in mapped connect indexes!";
			LOGPZ_ERROR(logger, sout.str().c_str());
			this-> fConnectIndexes[i] = -1;
			return;
		}
		this-> fConnectIndexes[i] = lcIdx;
	}
    // write the code when this constructor is called
    DebugStop();
}

// TESTADO
template<class TSHAPE>
TPZCompElHDivCollapsed<TSHAPE>::TPZCompElHDivCollapsed() :
TPZRegisterClassId(&TPZCompElHDivCollapsed::ClassId),
TPZCompElHDiv<TSHAPE>(),fBottom(),fTop()
{
	this->fPreferredOrder = -1;
	int i;
	for(i=0;i<TSHAPE::NFacets;i++) {
		this-> fConnectIndexes[i] = -1;
	}
}

// TESTADO
template<class TSHAPE>
TPZCompElHDivCollapsed<TSHAPE>::~TPZCompElHDivCollapsed(){
    TPZGeoEl *gel = this->Reference();
    if(!gel) return;
    if (gel && gel->Reference() != this) {
        return;
    }
    int side = TSHAPE::NSides-1;
    TPZGeoElSide gelside(this->Reference(),side);
    TPZStack<TPZCompElSide> celstack;
    TPZCompElSide largecel = gelside.LowerLevelCompElementList2(0);
    if (largecel) {
        int cindex = SideConnectLocId(0, side);
        TPZConnect &c = this->Connect(cindex);
        c.RemoveDepend();
    }
    gelside.HigherLevelCompElementList3(celstack, 0, 1);
    int64_t ncel = celstack.size();
    for (int64_t el=0; el<ncel; el++) {
        TPZCompElSide celsidesmall = celstack[el];
        TPZGeoElSide gelsidesmall = celsidesmall.Reference();
        if (gelsidesmall.Dimension() != gel->Dimension()) {
            continue;
        }
        TPZCompEl *cel = celsidesmall.Element();
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) {
            DebugStop();
        }
        int cindex = intel->SideConnectLocId(0, celsidesmall.Side());
        TPZConnect &c = intel->Connect(cindex);
        c.RemoveDepend();
    }
    if (gel){
        gel->ResetReference();
    }

}

// NAO TESTADO
template<class TSHAPE>
MElementType TPZCompElHDivCollapsed<TSHAPE>::Type() {
    DebugStop();
	return TSHAPE::Type();
}

// NAO TESTADO
template<class TSHAPE>
void TPZCompElHDivCollapsed<TSHAPE>::SetSideOrient(int side, int sideorient)
{
    if (side < TSHAPE::NFacets || side >= TSHAPE::NSides) {
        DebugStop();
    }
    if(side < TSHAPE::NSides-1) TPZCompElHDiv<TSHAPE>::SetSideOrient(side, sideorient);
    else fBottom.SetSideOrient(side, sideorient);
}

// NAO TESTADO
template<class TSHAPE>
int TPZCompElHDivCollapsed<TSHAPE>::GetSideOrient(int side)
{
    if (side >= TSHAPE::NFacets && side < TSHAPE::NSides - 1) {
        return TPZCompElHDiv<TSHAPE>::GetSideOrient(side);
    }
    else if(side == TSHAPE::NSides-1) return fBottom.GetSideOrient(side);
    DebugStop();
}

template<class TSHAPE>
int TPZCompElHDivCollapsed<TSHAPE>::NConnects() const {
	
	return TPZCompElHDiv<TSHAPE>::NConnects()+2;
}

template<class TSHAPE>
void TPZCompElHDivCollapsed<TSHAPE>::SetConnectIndex(int i, int64_t connectindex)
{
	if(i < TSHAPE::NFacets)
	{
        SetConnectIndex(i, connectindex);
	}
    else if(i == TSHAPE::NFacets)
    {
        fBottom.SetConnectIndex(0, connectindex);
    }
    else if(i == TSHAPE::NFacets+1)
    {
        fTop.SetConnectIndex(0, connectindex);
    }
    else
    {
        DebugStop();
    }
	
}

template<class TSHAPE>
int TPZCompElHDivCollapsed<TSHAPE>::NConnectShapeF(int connect, int connectorder) const
{
    if(connect <= TSHAPE::NFacets)
    {
        return TPZCompElHDiv<TSHAPE>::NConnectShapeF(connect, connectorder);
    }
    else if(connect <= TSHAPE::NFacets+2)
    {
        return fBottom.NConnectShapeF(0, connectorder);
    }
    else DebugStop();
    return -1;
}

template<class TSHAPE>
int TPZCompElHDivCollapsed<TSHAPE>::NSideConnects(int side) const
{
	if(side == TSHAPE::NSides-1)
	{
		return 1;
	}
    else if(side < TSHAPE::NSides-1)
    {
        return TPZCompElHDiv<TSHAPE>::NSideConnects(side);
    }
    DebugStop();
	return 0;
}

template<class TSHAPE>
int TPZCompElHDivCollapsed<TSHAPE>::SideConnectLocId(int node, int side) const
{
	if(side == TSHAPE::NSides-1 && node < 2)
	{
		return TSHAPE::NFacets+1+node;
	}
	else{
        return TPZCompElHDiv<TSHAPE>::SideConnectLocId(node,side);
	}
	
}

//Identifies the interpolation order on the connects of the element
template<class TSHAPE>
void TPZCompElHDivCollapsed<TSHAPE>::GetInterpolationOrder(TPZVec<int> &ord) {
    int myorder = ConnectOrder(0);
    ord.Resize(TSHAPE::NSides-TSHAPE::NCornerNodes, 0);
	int i;
	for(i=0; i<ord.size(); i++) {
		ord[i] = myorder;
	}
}

//return the preferred order of the polynomial along connect connect
template<class TSHAPE>
int TPZCompElHDivCollapsed<TSHAPE>::PreferredSideOrder(int side) {
	if(side != TSHAPE::NSides-1)
	{
		PZError << __PRETTY_FUNCTION__ << " side " << side << std::endl;
	}
	int connect= 0;
	int order =this->fPreferredOrder;
	return this->AdjustPreferredSideOrder(connect,order);
}

//sets the interpolation order of side to order
template<class TSHAPE>
void TPZCompElHDivCollapsed<TSHAPE>::SetSideOrder(int side, int order) {
	int connectaux= SideConnectLocId(0,side);
	if(connectaux<0 || connectaux > this-> NConnects()) {
		PZError << "TPZCompElHDiv::SetSideOrder. Bad paramenter side " << side << " order " << order << std::endl;
#ifdef LOG4CXX
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << " Bad side or order " << side << " order " << order;
		LOGPZ_DEBUG(logger,sout.str())
#endif
		return;
	}
	TPZConnect &c = this->Connect(connectaux);
    c.SetOrder(order,this->fConnectIndexes[connectaux]);
    int64_t seqnum = c.SequenceNumber();
    int nvar = 1;
    TPZMaterial * mat =this-> Material();
    if(mat) nvar = mat->NStateVariables();
    int nshape = NConnectShapeF(connectaux,order);
    c.SetNShape(nshape);
    c.SetNState(nvar);
    this-> Mesh()->Block().Set(seqnum,nshape*nvar);
    if(connectaux == NConnects()-1)
    {
		this->SetIntegrationRule(2*order);
    }
}

/**
 return the interpolation orderof the polynomial for connect
 **/
template<class TSHAPE>
int TPZCompElHDivCollapsed<TSHAPE>::ConnectOrder(int connect) const
{
	
	if (connect < 0 || connect >= this->NConnects())
	{
#ifdef LOG4CXX
		{
			std::stringstream sout;
			sout << "Connect index out of range connect " << connect <<
			" nconnects " << NConnects();
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
		return -1;
	}
	
	if (this->fConnectIndexes[connect] == -1) {
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << " connect " << connect
		<< " is not initialized" << std::endl;
#ifdef LOG4CXX
		LOGPZ_ERROR(logger,sout.str());
		DebugStop();
#else
		std::cout << sout.str() << std::endl;
#endif
		return -1;
	}
	const TPZConnect &c = this-> Connect(connect);
	return c.Order();
}

/** Initialize a material data and its attributes based on element dimension, number
 * of state variables and material definitions
 */

template<class TSHAPE>
void TPZCompElHDivCollapsed<TSHAPE>::InitMaterialData(TPZMaterialData &data)
{
	TPZCompElHDiv<TSHAPE>::InitMaterialData(data);
    data.fShapeType = TPZMaterialData::EVecShape;
    // expand the shape vector and normal vector
	
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
	{
        std::stringstream sout;
        sout << "After InitMaterialData\n";
        data.Print(sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	
}

template<class TSHAPE>
void TPZCompElHDivCollapsed<TSHAPE>::ComputeShapeIndex(TPZVec<int> &sides, TPZVec<int64_t> &shapeindex) {
	
	TPZManVector<int64_t> firstshapeindex;    // Para o que?
	TPZCompElHDiv<TSHAPE>::FirstShapeIndex(firstshapeindex);      // se foram calculados os indices mas n�o utilizados?
	int nshape = TPZIntelGen<TSHAPE>::NShapeF();
	shapeindex.Resize(nshape);
	int64_t nsides = sides.NElements();
	int64_t is, count=0;
	for(is=0 ; is<nsides; is++)
	{
		int side = sides[is];
		int sideorder= this->EffectiveSideOrder(side);
		int NShapeFace = TSHAPE::NConnectShapeF(side,sideorder);
		int ishapeface;
		for(ishapeface=0; ishapeface<NShapeFace; ishapeface++)
		{
			shapeindex[count++] = is;
		}
	}
	shapeindex.Resize(count);
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "count = " << count << " nshape " << nshape;
		sout << std::endl<<"sides associated with the normals "<< sides <<
		"\nnormal associated with each shape function : shape function indexes " << shapeindex;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	
}

template<class TSHAPE>
void TPZCompElHDivCollapsed<TSHAPE>::SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
	
    if(TSHAPE::SideDimension(side)!= TSHAPE::Dimension || point.size() != TSHAPE::Dimension ){
		DebugStop() ;
	}
    fBottom.SideShapeFunction(side, point, phi, dphi);
    
    return;
}

/** Compute the shape function at the integration point */
template<class TSHAPE>
void TPZCompElHDivCollapsed<TSHAPE>::Shape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi)
{
    fBottom.Shape(pt, phi, dphi);
}

/** Read the element data from a stream */
template<class TSHAPE>
void TPZCompElHDivCollapsed<TSHAPE>::Read(TPZStream &buf, void *context)
{
	TPZIntelGen<TSHAPE>::Read(buf,context);
    fBottom.Read(buf,context);
    fTop.Read(buf, context);
}

/** Save the element data to a stream */
template<class TSHAPE>
void TPZCompElHDivCollapsed<TSHAPE>::Write(TPZStream &buf, int withclassid) const
{
	TPZIntelGen<TSHAPE>::Write(buf,withclassid);
    fBottom.Write(buf, false);
    fTop.Write(buf, false);
}

/** @brief Prints the relevant data of the element to the output stream */
template<class TSHAPE>
void TPZCompElHDivCollapsed<TSHAPE>::Print(std::ostream &out) const
{
    out << __PRETTY_FUNCTION__ << std::endl;
    TPZCompElHDiv<TSHAPE>::Print(out);
    fBottom.Print(out);
    fTop.Print(out);
    
}

/** Returns the actual interpolation order of the polynomial along the side */
template<class TSHAPE>
int TPZCompElHDivCollapsed<TSHAPE>::EffectiveSideOrder(int side) const
{
    if(side >= 0 && side < TSHAPE::NSides-1)
    {
        return TPZCompElHDiv<TSHAPE>::EffectiveSideOrder(side);
    }
	if(side == TSHAPE::NSides-1)
	{
		return fBottom.EffectiveSideOrder(side);
	}
	else {
		return -1;
	}
}

template<class TSHAPE>
void TPZCompElHDivCollapsed<TSHAPE>::ComputeRequiredData(TPZMaterialData &data,
                                                TPZVec<REAL> &qsi){
    
//    TPZManVector<int,TSHAPE::NSides*TSHAPE::Dimension> normalsidesDG(TSHAPE::Dimension*TSHAPE::NSides);

    bool needsol = data.fNeedsSol;
    data.fNeedsSol = false;
    TPZIntelGen<TSHAPE>::ComputeRequiredData(data,qsi);
    data.fNeedsSol = needsol;
    
    int restrainedface = this->RestrainedFace();
    // Acerta o vetor data.fDeformedDirections para considerar a direcao do campo. fSideOrient diz se a orientacao e de entrada
    // no elemento (-1) ou de saida (+1), dependedo se aquele lado eh vizinho pela direita (-1) ou pela esquerda(+1)
    int firstface = TSHAPE::NSides - TSHAPE::NFacets - 1;
    int lastface = TSHAPE::NSides - 1;
    int cont = 0;
   
    TPZIntelGen<TSHAPE>::Reference()->HDivDirectionsMaster(data.fMasterDirections);

    for(int side = firstface; side < lastface; side++)
    {
        int nvec = TSHAPE::NContainedSides(side);
        for (int ivet = 0; ivet<nvec; ivet++)
        {
            for (int il = 0; il<3; il++)
            {
                data.fMasterDirections(il,ivet+cont) *= SideOrient(side-firstface);
            }

        }
        cont += nvec;
    }
    
    if(data.fNeedsDeformedDirectionsFad){
    #ifdef _AUTODIFF
        TPZIntelGen<TSHAPE>::Reference()->HDivDirections(qsi,data.fDeformedDirectionsFad,restrainedface);
        cont = 0;
        
        for(int side = firstface; side < lastface; side++)
        {
            int nvec = TSHAPE::NContainedSides(side);
            for (int ivet = 0; ivet<nvec; ivet++)
            {
                for (int il = 0; il<3; il++)
                {
                    data.fDeformedDirectionsFad(il,ivet+cont) *= SideOrient(side-firstface);
                }
                
            }
            cont += nvec;
        }
    #else
        DebugStop();
    #endif
    }else{
        TPZIntelGen<TSHAPE>::Reference()->HDivDirections(qsi,data.fDeformedDirections,restrainedface);
        cont = 0;
    
        for(int side = firstface; side < lastface; side++)
        {
            int nvec = TSHAPE::NContainedSides(side);
            for (int ivet = 0; ivet<nvec; ivet++)
            {
                for (int il = 0; il<3; il++)
                {
                    data.fDeformedDirections(il,ivet+cont) *= SideOrient(side-firstface);
                }
                
            }
            cont += nvec;
        }
    }
    
    data.ComputeFunctionDivergence();
    if (data.fNeedsSol) {
        TPZCompElHDiv<TSHAPE>::ComputeSolution(qsi, data);
    }


#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        data.fDeformedDirections.Print("Normal Vectors " , sout,EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    

}//void

template<class TSHAPE>
int64_t TPZCompElHDivCollapsed<TSHAPE>::ConnectIndex(int con) const
{
    if(con <= TSHAPE::NFacets) return TPZCompElHDiv<TSHAPE>::ConnectIndex(con);
    if(con > TSHAPE::NFacets + 2) DebugStop();
    if(con == TSHAPE::NFacets+1) return fBottom.ConnectIndex(0);
    if(con == TSHAPE::NFacets+2) return fTop.ConnectIndex(0);
    DebugStop();
    return -1;
}



template<class TSHAPE>
void TPZCompElHDivCollapsed<TSHAPE>::SetCreateFunctions(TPZCompMesh* mesh) {
    mesh->SetAllCreateFunctionsHDiv();
}

#include "pzshapetriang.h"
#include "pzshapepoint.h"
#include "pzshapelinear.h"
#include "pzshapequad.h"

using namespace pzshape;

template class TPZRestoreClass< TPZCompElHDivCollapsed<TPZShapeLinear>>;
template class TPZRestoreClass< TPZCompElHDivCollapsed<TPZShapeTriang>>;
template class TPZRestoreClass< TPZCompElHDivCollapsed<TPZShapeQuad>>;


template class TPZCompElHDivCollapsed<TPZShapeTriang>;
template class TPZCompElHDivCollapsed<TPZShapeLinear>;
template class TPZCompElHDivCollapsed<TPZShapeQuad>;

