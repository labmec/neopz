/**
 * @file
 * @brief Contains the implementation of the TPZCompElHDivBound2 methods.
 */


#include "pzelchdivbound2.h"
#include "pzgeoel.h"
#include "TPZMaterial.h"
#include "pzlog.h"
#include "TPZShapeDisc.h"
#include "TPZCompElDisc.h"
#include "TPZMaterialDataT.h"
#include "pzelchdiv.h"
#include "TPZShapeHDivBound.h"
#include "TPZShapeHDivConstantBound.h"
#include "TPZShapeHCurlNoGrads.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.TPZCompElHDivBound2");
#endif

template<class TSHAPE>
TPZCompElHDivBound2<TSHAPE>::TPZCompElHDivBound2(TPZCompMesh &mesh, TPZGeoEl *gel, const HDivFamily hdivfam) :
TPZRegisterClassId(&TPZCompElHDivBound2::ClassId),
TPZIntelGen<TSHAPE>(mesh,gel,1), fSideOrient(1), fhdivfam(hdivfam){
		
	//int i;
	this->TPZInterpolationSpace::fPreferredOrder = mesh.GetDefaultOrder();
	//for(i=0; i<TSHAPE::NSides; i++) this->fConnectIndexes[i]=-1;
  this->fConnectIndexes[0]=-1;
  gel->SetReference(this);
		
  this->fConnectIndexes[0] = this->CreateMidSideConnect(TSHAPE::NSides-1);
#ifdef PZ_LOG
  if (logger.isDebugEnabled())
		{
      std::stringstream sout;
      sout << "After creating boundary flux connect " << this->fConnectIndexes[0] << std::endl;
      //	this->Print(sout);
      LOGPZ_DEBUG(logger,sout.str())
        }
#endif
    
  mesh.ConnectVec()[this->fConnectIndexes[0]].IncrementElConnected();
    
		
	
#ifdef PZ_LOG
  if (logger.isDebugEnabled())
    {
      std::stringstream sout;
		
      sout << std::endl<<" Criando Connects: "<< std::endl;
      for(int j=0; j< NConnects();j++)
        {
          sout<<" "<< this->fConnectIndexes[j];
			
        }
      LOGPZ_DEBUG(logger,sout.str())
        }
#endif
	int sideorder = EffectiveSideOrder(TSHAPE::NSides-1);
	sideorder = 2*sideorder;
	if (sideorder > this->fIntRule.GetMaxOrder()) sideorder = this->fIntRule.GetMaxOrder();
	//  TPZManVector<int,3> order(3,2*sideorder+2);
	TPZManVector<int,3> order(3,sideorder);
	//TPZManVector<int,3> order(3,20);
	this->fIntRule.SetOrder(order);


    if (fhdivfam == HDivFamily::EHDivConstant) {
        // For HDiv constant, polynomial order was compatibilized in connectorders, 
        // see TPZShapeHDivConstantBound<TSHAPE>::Initialize. So now we need to update
        // the number of shape functions and also the integration rule
        if (TSHAPE::Type() == ETriangle || TSHAPE::Type() == EOned){
            for (int icon = 0; icon < this->NConnects(); icon++)
            {
                TPZConnect &c = this->Connect(icon);
                int nShapeF = NConnectShapeF(icon,c.Order());
                if (c.NShape() != nShapeF){
                    DebugStop();
                }
            }
        }
    }

#ifdef PZ_LOG
  if (logger.isDebugEnabled())
    {
      std::stringstream sout;
      sout << "Finalizando criacao do elemento ";
      this->Print(sout);
      LOGPZ_DEBUG(logger,sout.str())
        }
#endif
	 
}

template<class TSHAPE>
TPZCompElHDivBound2<TSHAPE>::TPZCompElHDivBound2(TPZCompMesh &mesh, const TPZCompElHDivBound2<TSHAPE> &copy) :
TPZRegisterClassId(&TPZCompElHDivBound2::ClassId),
TPZIntelGen<TSHAPE>(mesh,copy), fSideOrient(copy.fSideOrient), fConnectIndexes(copy.fConnectIndexes), fhdivfam(copy.fhdivfam)
{
#ifdef PZDEBUG
    if(fConnectIndexes[0] != copy.fConnectIndexes[0]) DebugStop();
#endif
}

// NAO TESTADO
template<class TSHAPE>
TPZCompElHDivBound2<TSHAPE>::TPZCompElHDivBound2(TPZCompMesh &mesh,
												 const TPZCompElHDivBound2<TSHAPE> &copy,
												 std::map<int64_t,int64_t> & gl2lcConMap,
												 std::map<int64_t,int64_t> & gl2lcElMap) :
TPZRegisterClassId(&TPZCompElHDivBound2::ClassId),
TPZIntelGen<TSHAPE>(mesh,copy,gl2lcConMap,gl2lcElMap), fSideOrient(copy.fSideOrient), fhdivfam(copy.fhdivfam)
{
	
	this-> fPreferredOrder = copy.fPreferredOrder;
	int i;
	for(i=0;i<TSHAPE::NSides;i++)
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
    
}

// TESTADO
template<class TSHAPE>
TPZCompElHDivBound2<TSHAPE>::TPZCompElHDivBound2() : 
TPZRegisterClassId(&TPZCompElHDivBound2::ClassId),
TPZIntelGen<TSHAPE>()
{
	this->fPreferredOrder = -1;
	int i;
	for(i=0;i<TSHAPE::NSides;i++) {
		this-> fConnectIndexes[i] = -1;
	}
}

// TESTADO
template<class TSHAPE>
TPZCompElHDivBound2<TSHAPE>::~TPZCompElHDivBound2(){
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
MElementType TPZCompElHDivBound2<TSHAPE>::Type() {
	return TSHAPE::Type();
}

// NAO TESTADO
template<class TSHAPE>
void TPZCompElHDivBound2<TSHAPE>::SetSideOrient(int side, int sideorient)
{
    if (side != TSHAPE::NSides - 1) {
        DebugStop();
    }
    fSideOrient = sideorient;
}

// NAO TESTADO
template<class TSHAPE>
int TPZCompElHDivBound2<TSHAPE>::GetSideOrient(int side)
{
    if (side != TSHAPE::NSides - 1) {
        DebugStop();
    }
    return fSideOrient;
}

template<class TSHAPE>
int TPZCompElHDivBound2<TSHAPE>::NConnects() const {
	
	return 1;
}

template<class TSHAPE>
void TPZCompElHDivBound2<TSHAPE>::SetConnectIndex(int i, int64_t connectindex)
{
	if(i)
	{
		DebugStop();
	}
	this->fConnectIndexes[i] = connectindex;
	
}

template<class TSHAPE>
int TPZCompElHDivBound2<TSHAPE>::NConnectShapeF(int connect, int connectorder) const
{
#ifdef DEBUG
    if (connect < 0 || connect > TSHAPE::NFacets) {
        DebugStop();
    }
#endif

    if(connect == 0)
    {
        if(connectorder == 0) return 1;
        TPZManVector<int,22> order(TSHAPE::NSides-TSHAPE::NCornerNodes,connectorder);
        return TSHAPE::NShapeF(order);
    }    
    
    return -1;
}

template<class TSHAPE>
int TPZCompElHDivBound2<TSHAPE>::NSideConnects(int side) const
{
	if(side == TSHAPE::NSides-1)
	{
		return 1;
	}
	return 0;
}

template<class TSHAPE>
int TPZCompElHDivBound2<TSHAPE>::SideConnectLocId(int node, int side) const
{
	if(side == TSHAPE::NSides-1 && node == 0)
	{
		return 0;
	}
	else{
	return -1;
	}
	
}

//Identifies the interpolation order on the connects of the element
template<class TSHAPE>
void TPZCompElHDivBound2<TSHAPE>::GetInterpolationOrder(TPZVec<int> &ord) {
    int myorder = ConnectOrder(0);
    ord.Resize(TSHAPE::NSides-TSHAPE::NCornerNodes, 0);
	int i;
	for(i=0; i<ord.size(); i++) {
		ord[i] = myorder;
	}
}

//return the preferred order of the polynomial along connect connect
template<class TSHAPE>
int TPZCompElHDivBound2<TSHAPE>::PreferredSideOrder(int side) {
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
void TPZCompElHDivBound2<TSHAPE>::SetSideOrder(int side, int order) {
	int connectaux= SideConnectLocId(0,side);
	if(connectaux<0 || connectaux > this-> NConnects()) {
		PZError << "TPZCompElHDiv::SetSideOrder. Bad paramenter side " << side << " order " << order << std::endl;
#ifdef PZ_LOG
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
int TPZCompElHDivBound2<TSHAPE>::ConnectOrder(int connect) const
{
	
	if (connect < 0 || connect >= this->NConnects())
	{
#ifdef PZ_LOG
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
#ifdef PZ_LOG
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
void TPZCompElHDivBound2<TSHAPE>::InitMaterialData(TPZMaterialData &data)
{
	TPZIntelGen<TSHAPE>::InitMaterialData(data);
    data.fShapeType = TPZMaterialData::EScalarShape;
	
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
	{
		LOGPZ_DEBUG(logger,"Initializing normal vectors")
	}
#endif
    TPZGeoEl *gel = this->Reference();
    int nc = gel->NCornerNodes();
    TPZManVector<int64_t,8> id(nc);
    for (int ic=0; ic<nc; ic++) {
        id[ic] = gel->Node(ic).Id();
    }
    int connectorder = this->Connect(0).Order();
    int sideorient = fSideOrient;

    // fill in the datastructures of shapedata
    switch (fhdivfam)
    {
    case HDivFamily::EHDivStandard:
        TPZShapeHDivBound<TSHAPE>::Initialize(id, connectorder, sideorient, data);
        break;
    case HDivFamily::EHDivConstant:
        TPZShapeHDivConstantBound<TSHAPE>::Initialize(id, connectorder, sideorient, data);
        break;
    
    default:
        DebugStop();
        break;
    }
}

template<class TSHAPE>
void TPZCompElHDivBound2<TSHAPE>::ComputeShapeIndex(TPZVec<int> &sides, TPZVec<int64_t> &shapeindex) {
	
	TPZManVector<int64_t> firstshapeindex;    // Para o que?
	FirstShapeIndex(firstshapeindex);      // se foram calculados os indices mas n�o utilizados?
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
    #ifdef PZ_LOGTPZCompElHDivBound2
	{
		std::stringstream sout;
		sout << "count = " << count << " nshape " << nshape;
		sout << std::endl<<"sides associated with the normals "<< sides <<
		"\nnormal associated with each shape function : shape function indexes " << shapeindex;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	
}

/**return the first shape associate to each side*/
template<class TSHAPE>
void TPZCompElHDivBound2<TSHAPE>::FirstShapeIndex(TPZVec<int64_t> &Index){
	
	Index.Resize(TSHAPE::NSides+1);
	Index[0]=0;
    int order = ConnectOrder(0);
		
    for(int iside=0;iside<TSHAPE::NSides;iside++)
    {
        
        if(TSHAPE::Type()==EQuadrilateral){
            Index[iside+1] = Index[iside] + TSHAPE::NConnectShapeF(iside,order);
        }
        else{
            Index[iside+1] = Index[iside] + TSHAPE::NConnectShapeF(iside,order);
        }
    }
    
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        sout << " FirsShapeIndex result " << Index;
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
	
	
}

template<class TSHAPE>
void TPZCompElHDivBound2<TSHAPE>::SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
	
    if(TSHAPE::SideDimension(side)!= TSHAPE::Dimension || point.size() != TSHAPE::Dimension ){
		DebugStop() ;
	}
    TPZGeoEl *gel = this->Reference();
    REAL detjac;
    {
        int dim = gel->SideDimension(side);
        TPZFNMatrix<9,REAL> jac(dim,dim),jacinv(dim,dim),axes(dim,3);
        gel->Jacobian(point, jac, axes, detjac, jacinv);
    }
    
    TPZShapeHDivBound<TSHAPE> shapehdiv;
    TPZShapeData shapedata;
    int nc = gel->NCornerNodes();
    TPZManVector<int64_t,8> id(nc);
    for (int ic=0; ic<nc; ic++) {
        id[ic] = gel->Node(ic).Id();
    }
    int connectorder = this->Connect(0).Order();
    int sideorient = 1;
    // fill in the datastructures of shapedata
    shapehdiv.Initialize(id, connectorder, sideorient, shapedata);
    // compute the shape functions at the integration point
    shapehdiv.Shape(point, shapedata, phi);
    phi *= 1./detjac;
        
}

/** Compute the shape function at the integration point */
template<class TSHAPE>
void TPZCompElHDivBound2<TSHAPE>::Shape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi)
{
    DebugStop(); // Is this ever used (Nov 2021). If yes, just uncomment me
    TPZManVector<int,TSHAPE::NSides> ordl;
    this->GetInterpolationOrder(ordl);
    TPZConnect &c = this->Connect(0);
    int nshape = NConnectShapeF(0,c.Order());
    phi.Resize(nshape, 1);
    dphi.Resize(TSHAPE::Dimension, nshape);
    SideShapeFunction(TSHAPE::NSides-1, pt, phi, dphi);

    if (fSideOrient == -1) {
        phi(0,0) *= -1.;
        dphi(0,0) *= -1.;
    }
    
}

template<class TSHAPE>
void TPZCompElHDivBound2<TSHAPE>::ComputeShape(TPZVec<REAL> &intpoint, TPZMaterialData &data){
    
    TPZShapeData shapedata(data);
    
    switch (fhdivfam)
    {
    case HDivFamily::EHDivStandard:
        {
            auto nShape = TPZShapeHDivBound<TSHAPE>::NShape(shapedata);
            data.phi.Resize(nShape, 1);
            TPZShapeHDivBound<TSHAPE>::Shape(intpoint, shapedata, data.phi);
        }
        break;
    case HDivFamily::EHDivConstant:
        {
            data.phi.Resize(this->NShapeF(), 1);
            TPZShapeHDivConstantBound<TSHAPE>::Shape(intpoint, shapedata, data.phi);
        }
        break;
       
    default:
        DebugStop();//You should chose an HDiv family space
        break;
    }
    
    data.phi *= 1./data.detjac;

}

template<class TSHAPE>
void TPZCompElHDivBound2<TSHAPE>::ComputeShape(TPZVec<REAL> &intpoint, TPZVec<REAL> &X, TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes,
                                         REAL &detjac, TPZFMatrix<REAL> &jacinv, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi, TPZFMatrix<REAL> &dphidx){
    
    std::cout << "Method not implement call the architec." << std::endl;
    DebugStop();
    
}

/** Read the element data from a stream */
template<class TSHAPE>
void TPZCompElHDivBound2<TSHAPE>::Read(TPZStream &buf, void *context)
{
	TPZIntelGen<TSHAPE>::Read(buf,context);
  buf.Read(&fSideOrient);
  buf.Read(fConnectIndexes.begin(),TSHAPE::NSides);
}

/** Save the element data to a stream */
template<class TSHAPE>
void TPZCompElHDivBound2<TSHAPE>::Write(TPZStream &buf, int withclassid) const
{
	TPZIntelGen<TSHAPE>::Write(buf,withclassid);
  buf.Write(&fSideOrient);
  buf.Write(fConnectIndexes.begin(),TSHAPE::NSides);
}

/** @brief Prints the relevant data of the element to the output stream */
template<class TSHAPE>
void TPZCompElHDivBound2<TSHAPE>::Print(std::ostream &out) const
{
    out << __PRETTY_FUNCTION__ << std::endl;
    out << "Side orientation " << fSideOrient << std::endl;
    TPZIntelGen<TSHAPE>::Print(out);
    
    
}

/** Returns the actual interpolation order of the polynomial along the side */
template<class TSHAPE>
int TPZCompElHDivBound2<TSHAPE>::EffectiveSideOrder(int side) const
{
	if(side == TSHAPE::NSides-1)
	{
		return ConnectOrder(0);
	}
	else {
		return -1;
	}
}

/** Return a matrix with index shape and vector associate to element */
template<class TSHAPE>
void TPZCompElHDivBound2<TSHAPE>::IndexShapeToVec(TPZVec<int> &VectorSide,TPZVec<std::pair<int,int64_t> > & ShapeAndVec){
	
	// VectorSide indicates the side associated with each vector entry
	TPZVec<int64_t> FirstIndex;
	// the first index of the shape functions
	FirstShapeIndex(FirstIndex);
#ifdef PZ_LOG
		{
				std::stringstream sout;
				sout << "FirstIndex of shape functions " << FirstIndex;
				LOGPZ_DEBUG(logger,sout.str())
		}
#endif
	
	int tamanho= this->NShapeF();
	
	ShapeAndVec.Resize(tamanho);
	int64_t count=0;
		//for(int jvec=0;jvec< VectorSide.NElements();jvec++)
	for(int64_t jvec=0;jvec< VectorSide.NElements();jvec++)//coloca-se -1 caso queira reduzir o espaco de fluxo
	{
		int lside=VectorSide[jvec];
		int64_t fshape1= FirstIndex[lside];
		int64_t fshape2= FirstIndex[lside+1];
		for (int64_t ishape=fshape1; ishape<fshape2; ishape++)
		{
			
#ifdef PZ_LOG
			std::stringstream sout;
			sout << " <vec,shape> " << "< "<<jvec << " * "<<ishape << "> "<<std::endl;
			LOGPZ_DEBUG(logger,sout.str())
#endif
			
			
			ShapeAndVec[count++]=std::pair<int,int64_t>(jvec,ishape);
		}
		
	}
	
#ifdef PZ_LOG
    std::stringstream sout;
    sout << "VecShapeIndex " << ShapeAndVec;
    LOGPZ_DEBUG(logger,sout.str())
#endif
	
	
}

template<class TSHAPE>
void TPZCompElHDivBound2<TSHAPE>::SetCreateFunctions(TPZCompMesh* mesh) {
    mesh->ApproxSpace().SetHDivFamily(fhdivfam);
    mesh->ApproxSpace().SetAllCreateFunctionsHDiv(TSHAPE::Dimension);
}


#include "pzshapetriang.h"
#include "pzshapepoint.h"
#include "pzshapelinear.h"
#include "pzshapequad.h"

using namespace pzshape;

template class TPZRestoreClass< TPZCompElHDivBound2<TPZShapePoint>>;
template class TPZRestoreClass< TPZCompElHDivBound2<TPZShapeLinear>>;
template class TPZRestoreClass< TPZCompElHDivBound2<TPZShapeTriang>>;
template class TPZRestoreClass< TPZCompElHDivBound2<TPZShapeQuad>>;


template class TPZCompElHDivBound2<TPZShapeTriang>;
template class TPZCompElHDivBound2<TPZShapePoint>;
template class TPZCompElHDivBound2<TPZShapeLinear>;
template class TPZCompElHDivBound2<TPZShapeQuad>;
