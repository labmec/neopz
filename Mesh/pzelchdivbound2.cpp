/**
 * @file
 * @brief Contains the implementation of the TPZCompElHDivBound2 methods.
 */
/*
 * pzelchdivbound.cpp
 *
 *  Created on: Sep 29, 2009
 *      Author: phil
 */

#include "pzelchdivbound2.h"
#include "pzgeoel.h"
#include "pzmaterial.h"
#include "pzlog.h"
#include "TPZShapeDisc.h"
#include "TPZCompElDisc.h"
#include "pzmaterialdata.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.TPZCompElHDivBound2"));
#endif

template<class TSHAPE>
TPZCompElHDivBound2<TSHAPE>::TPZCompElHDivBound2(TPZCompMesh &mesh, TPZGeoEl *gel, int &index) :
TPZIntelGen<TSHAPE>(mesh,gel,index,1) {
	int i;
	this->TPZInterpolationSpace::fPreferredOrder = mesh.GetDefaultOrder();
	for(i=0; i<TSHAPE::NSides; i++) this->fConnectIndexes[i]=-1;
	gel->SetReference(this);
	TPZGeoElSide neigh(gel->Neighbour(gel->NSides()-1));
	TPZCompElSide compneigh(neigh.Reference());
	int sideoffset = neigh.Element()->NSides()-neigh.Side();
	int neighnconnects = compneigh.Element()->NConnects();
	int connectnumber = neighnconnects-sideoffset;
	TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (compneigh.Element());
	connectnumber = intel->SideConnectLocId(0, compneigh.Side());
	this->fConnectIndexes[0] = compneigh.Element()->ConnectIndex(connectnumber);
	mesh.ConnectVec()[this->fConnectIndexes[0]].IncrementElConnected();
	
#ifdef LOG4CXX
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
	int sideorder = SideOrder(TSHAPE::NSides-1);
	sideorder = 2*sideorder;
	if (sideorder > this->fIntRule.GetMaxOrder()) sideorder = this->fIntRule.GetMaxOrder();
	//  TPZManVector<int,3> order(3,2*sideorder+2);
	TPZManVector<int,3> order(3,sideorder);
	//TPZManVector<int,3> order(3,20);
	this->fIntRule.SetOrder(order);
	/*
	 #ifdef LOG4CXX
	 {
	 std::stringstream sout;
	 sout << "Finalizando criacao do elemento ";
	 this->Print(sout);
	 LOGPZ_DEBUG(logger,sout.str())
	 }
	 #endif
	 */
}

template<class TSHAPE>
TPZCompElHDivBound2<TSHAPE>::TPZCompElHDivBound2(TPZCompMesh &mesh, const TPZCompElHDivBound2<TSHAPE> &copy) :
TPZIntelGen<TSHAPE>(mesh,copy)
{
	int i;
	for(i=0;i<TSHAPE::NSides;i++)
	{
		this-> fConnectIndexes[i] = copy.fConnectIndexes[i];
	}
}

// NAO TESTADO
template<class TSHAPE>
TPZCompElHDivBound2<TSHAPE>::TPZCompElHDivBound2(TPZCompMesh &mesh,
												 const TPZCompElHDivBound2<TSHAPE> &copy,
												 std::map<int,int> & gl2lcConMap,
												 std::map<int,int> & gl2lcElMap) :
TPZIntelGen<TSHAPE>(mesh,copy,gl2lcConMap,gl2lcElMap)
{
	
	this-> fPreferredOrder = copy.fPreferredOrder;
	int i;
	for(i=0;i<TSHAPE::NSides;i++)
	{
		int lcIdx = -1;
		int glIdx = copy.fConnectIndexes[i];
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
	//   gl2lcElMap[copy.fIndex] = this->Index();
}

// TESTADO
template<class TSHAPE>
TPZCompElHDivBound2<TSHAPE>::TPZCompElHDivBound2() : TPZIntelGen<TSHAPE>()
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
	
}

// NAO TESTADO
template<class TSHAPE>
MElementType TPZCompElHDivBound2<TSHAPE>::Type() {
	return TSHAPE::Type();
}


template<class TSHAPE>
int TPZCompElHDivBound2<TSHAPE>::NConnects() const {
	
	return 1;//acrescentando um connect mais pra variavel dual antes era apenas NumSides(dimension) + 1
}

template<class TSHAPE>
void TPZCompElHDivBound2<TSHAPE>::SetConnectIndex(int i, int connectindex)
{
	if(i)
	{
		DebugStop();
	}
	this->fConnectIndexes[i] = connectindex;
	
}

template<class TSHAPE>
int TPZCompElHDivBound2<TSHAPE>::NConnectShapeF(int connect) const
{
	if(connect == 0)
	{
		TPZManVector<int,1> order(1,ConnectOrder(connect));
		//se a ordem eh maior depedenra do tipo de TSHAPE
		
		if (order[0]==1) {
			return TSHAPE::NShapeF(order);
		}
		else{
			TPZGeoElSide gelside(this->Reference(),TSHAPE::NSides-1);
			TPZGeoElSide neighbour = gelside.Neighbour();
			while(gelside != neighbour)
			{	switch (neighbour.Element()->Type()) {
				case EQuadrilateral:
					return TSHAPE::NShapeF(order)-1;
					break;
				case ETriangle:
					return TSHAPE::NShapeF(order);
					break;
				default : DebugStop();return -1;
			}
				neighbour = neighbour.Neighbour();
				
				
			}
			
			
		}
		
	}
    return -1;
	/*if(connect == 0)
	{
		
		TPZManVector<int,1> order(1,ConnectOrder(connect));
		//se a ordem eh maior depedenra do tipo de TSHAPE
		
		if (order[0]>1) {
			return TSHAPE::NShapeF(order)-1;
		}
		else{
			return TSHAPE::NShapeF(order);
		}
		
	}
	
	
	else
	{
		DebugStop();
		return -1;
	}
	 */
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
	return -1;
	
}

//Identifies the interpolation order on the connects of the element
template<class TSHAPE>
void TPZCompElHDivBound2<TSHAPE>::GetInterpolationOrder(TPZVec<int> &ord) {
	ord.Resize(NConnects());
	int i;
	for(i=0; i<NConnects(); i++) {
		ord[i] = ConnectOrder(i);
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
#ifdef LOG4CXX
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << " Bad side or order " << side << " order " << order;
		LOGPZ_DEBUG(logger,sout.str())
#endif
		return;
 	}
	TPZConnect &c = this->Connect(connectaux);
    c.SetOrder(order);
    int seqnum = c.SequenceNumber();
    int nvar = 1;
    TPZAutoPointer<TPZMaterial> mat =this-> Material();
    if(mat) nvar = mat->NStateVariables();
    int nshape = NConnectShapeF(connectaux);
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
void TPZCompElHDivBound2<TSHAPE>::InitMaterialData(TPZMaterialData &data)
{
	TPZIntelGen<TSHAPE>::InitMaterialData(data);
	
#ifdef LOG4CXX
	{
		LOGPZ_DEBUG(logger,"Initializing normal vectors")
	}
#endif
	
	//data.fVecShapeIndex=true;
	
	TPZGeoElSide gelside(this->Reference(),TSHAPE::NSides-1);
	TPZGeoElSide neighbour = gelside.Neighbour();
	while(gelside != neighbour && neighbour.Element()->Dimension() != TSHAPE::Dimension+1)
	{
		neighbour = neighbour.Neighbour();
	}
	if(neighbour.Element()->Dimension() != TSHAPE::Dimension+1)
	{
		return;
	}
	TPZGeoEl *neighel = neighbour.Element();
	TPZManVector<int> normalsides;
	neighel->ComputeNormals(neighbour.Side(),data.fNormalVec, normalsides);
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "normal side depois do ComputeNormals " << normalsides << std::endl;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
		
	// relate the sides indicated in vecindex to the sides of the current element
	int nvec = normalsides.NElements();
	int ivec;
	for(ivec=0; ivec<nvec; ivec++)
	{
		TPZGeoElSide neigh(neighel,normalsides[ivec]);
#ifdef LOG4CXX
		{
			std::stringstream sout;
			sout << "normal side depois do TPZGeoElSide " << normalsides << std::endl;
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
		while(neigh.Element() != this->Reference())
		{
		
			neigh = neigh.Neighbour();
		}
		
			normalsides[ivec]=neigh.Side();

		
std::cout<< "neigh aqui ------"<<neigh.Side()<<std::endl;
		
#ifdef LOG4CXX
		{
			std::stringstream sout;
			sout << "ivec " << ivec<< "normal side " << neigh.Side() << std::endl;
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
	}
	IndexShapeToVec(normalsides,data.fVecShapeIndex);
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout<< "normal sides depois de IndexShapeToVec "<<normalsides<<std::endl;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	//	TPZInterpolationSpace *cel = dynamic_cast<TPZInterpolationSpace*> (neighel->Reference());
	data.numberdualfunctions = 0;
#ifdef LOG4CXX
	{
		std::stringstream sout;
		data.fNormalVec.Print("Normal vectors", sout);
		sout << "Vector/Shape indexes " << data.fVecShapeIndex << std::endl;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	 
}


template<class TSHAPE>
void TPZCompElHDivBound2<TSHAPE>::ComputeShapeIndex(TPZVec<int> &sides, TPZVec<int> &shapeindex){
	
	TPZManVector<int> firstshapeindex;
	FirstShapeIndex(firstshapeindex);
	int nshape = TPZIntelGen<TSHAPE>::NShapeF();
	shapeindex.Resize(nshape);
	int nsides = sides.NElements();
	int is, count=0;
	for(is=0 ; is<nsides; is++)
	{
		int side = sides[is];
		int sideorder= this->SideOrder(side);
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

/**return the first shape associate to each side*/
template<class TSHAPE>
void TPZCompElHDivBound2<TSHAPE>::FirstShapeIndex(TPZVec<int> &Index){
	
	Index.Resize(TSHAPE::NSides+1);
	Index[0]=0;
	
	for(int iside=0;iside<TSHAPE::NSides;iside++)
	{
		if(TSHAPE::Type()==EQuadrilateral){
		int order= SideOrder(iside)-1;
			Index[iside+1] = Index[iside] + TSHAPE::NConnectShapeF(iside,order);
		}
		else{
		int order= SideOrder(iside);
			Index[iside+1] = Index[iside] + TSHAPE::NConnectShapeF(iside,order);
		}
		
		
		
	}
	/*
	 Index.Resize(TSHAPE::NSides+1);
	 Index[0]=0;
	 
	 for(int iside=0;iside<TSHAPE::NSides;iside++)
	 {
	 #warning Alteri AQ estou verificando se o lado e de dimensao 1 e tirando a ultima function
	 if (TSHAPE::SideDimension(iside)==TSHAPE::Dimension) {
	 int order= SideOrder(iside)-1;//estava -1 para quadrilatero caso contrario tira
	 Index[iside+1] = Index[iside] + TSHAPE::NConnectShapeF(iside,order);
	 }
	 
	 else{
	 int order= SideOrder(iside);
	 Index[iside+1] = Index[iside] + TSHAPE::NConnectShapeF(iside,order);
	 
	 }
	 
	 }
	 */
	
#ifdef LOG4CXX
    std::stringstream sout;
    sout << " FirsShapeIndex result " << Index;
    LOGPZ_DEBUG(logger,sout.str())
#endif
	
	
}

template<class TSHAPE>
void TPZCompElHDivBound2<TSHAPE>::SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix &phi,TPZFMatrix &dphi) {
	
	if(TSHAPE::SideDimension(side)!= TSHAPE::Dimension ){
		return ;
	}
    phi.Zero();
    dphi.Zero();
    return;
//	TPZIntelGen<TSHAPE>::SideShapeFunction(side,point,phi,dphi);
	
}

/**
 * Compute the shape function at the integration point
 */
template<class TSHAPE>
void TPZCompElHDivBound2<TSHAPE>::Shape(TPZVec<REAL> &pt, TPZFMatrix &phi, TPZFMatrix &dphi)
{
	//TPZFNMatrix<300> phi1;
	//TPZFNMatrix<600> dphi1;
	
	TPZCompElSide thisside(this,TSHAPE::NSides-1);
	TPZGeoElSide thisgeoside(thisside.Reference());
	TPZGeoElSide neighgeo(thisgeoside.Neighbour());
	TPZCompElSide neigh(neighgeo.Reference());
	TPZInterpolatedElement *neighel = dynamic_cast<TPZInterpolatedElement *> (neigh.Element());
	if(!neighel)
	{
		LOGPZ_ERROR(logger,"Inconsistent neighbour")
		return;
	}
	int nshapeneigh = neighel->NSideShapeF(neighgeo.Side())+1;
	//	int nconnectneigh = neighel->NConnects();
	//	int ndiscshape = neighel->NConnectShapeF(nconnectneigh-1);
	phi.Redim(nshapeneigh, 1);
	dphi.Redim(neighgeo.Element()->Dimension(), nshapeneigh);
	TPZTransform tr(thisgeoside.Dimension()),tr2; 
	thisgeoside.SideTransform3(neighgeo, tr);
	TPZManVector<REAL,3> pt2(neighgeo.Dimension()),pt3(neighel->Dimension());
	tr.Apply(pt, pt2);
	neighel->SideShapeFunction(neigh.Side(), pt2, phi, dphi);
	//	TPZGeoEl *neighgeoel = neighgeo.Element();
	//	tr2 = neighgeoel->SideToSideTransform(neighgeo.Side(), neighgeoel->NSides()-1);
	//	tr2.Apply(pt2, pt3);
	//	neighel->SideShapeFunction(neighgeoel->NSides(), pt3, phi, dphi);
}

/**
 Read the element data from a stream
 */
template<class TSHAPE>
void TPZCompElHDivBound2<TSHAPE>::Read(TPZStream &buf, void *context)
{
	TPZIntelGen<TSHAPE>::Read(buf,context);
}

/**
 Save the element data to a stream
 */
template<class TSHAPE>
void TPZCompElHDivBound2<TSHAPE>::Write(TPZStream &buf, int withclassid)
{
	TPZIntelGen<TSHAPE>::Write(buf,withclassid);
}

/**returns the actual interpolation order of the polynomial along the side*/
template<class TSHAPE>
int TPZCompElHDivBound2<TSHAPE>::SideOrder(int side) const
{
	if(side == TSHAPE::NSides-1)
	{
		return ConnectOrder(0);
	}
	else {
		return -1;
	}
}

//return a matrix with index shape and vector associate to element
template<class TSHAPE>
void TPZCompElHDivBound2<TSHAPE>::IndexShapeToVec(TPZVec<int> &VectorSide,TPZVec<std::pair<int,int> > & ShapeAndVec){
	
	// VectorSide indicates the side associated with each vector entry
	TPZVec<int> FirstIndex;
	// the first index of the shape functions
	FirstShapeIndex(FirstIndex);
	/*
	 #ifdef LOG4CXX
	 {
	 std::stringstream sout;
	 sout << "FirstIndex of shape functions " << FirstIndex;
	 LOGPZ_DEBUG(logger,sout.str())
	 }
	 #endif
	 */
	int tamanho= this->NShapeF();
	
	ShapeAndVec.Resize(tamanho);
	int count=0;
	for(int jvec=0;jvec< VectorSide.NElements();jvec++)//coloquei -1
	{
		int lside=VectorSide[jvec];
		int fshape1= FirstIndex[lside];
		int fshape2= FirstIndex[lside+1];
		for (int ishape=fshape1; ishape<fshape2; ishape++)
		{
			
#ifdef LOG4CXX
			std::stringstream sout;
			sout << " <vec,shape> " << "< "<<jvec << " * "<<ishape << "> "<<std::endl;
			LOGPZ_DEBUG(logger,sout.str())
#endif
			
			
			ShapeAndVec[count++]=std::pair<int,int>(jvec,ishape);
		}
		
	}
	
#ifdef LOG4CXX
    std::stringstream sout;
    sout << "VecShapeIndex " << ShapeAndVec;
    LOGPZ_DEBUG(logger,sout.str())
#endif
	
	
}


#include "pzshapetriang.h"
#include "pzshapepoint.h"
#include "pzshapelinear.h"
#include "pzshapequad.h"

using namespace pzshape;
/**
 * returns the unique identifier for reading/writing objects to streams
 */
template<>
int TPZCompElHDivBound2<TPZShapePoint>::ClassId() const
{
	return TPZHDIVBOUND2POINTID;
}

template class
TPZRestoreClass< TPZCompElHDivBound2<TPZShapePoint>, TPZHDIVBOUNDPOINTID>;

template<>
int TPZCompElHDivBound2<TPZShapeLinear>::ClassId() const
{
	return TPZHDIVBOUND2LINEARID;
}

template class
TPZRestoreClass< TPZCompElHDivBound2<TPZShapeLinear>, TPZHDIVBOUNDLINEARID>;

template<>
int TPZCompElHDivBound2<TPZShapeTriang>::ClassId() const
{
	return TPZHDIVBOUND2TRIANGLEID;
}

template class
TPZRestoreClass< TPZCompElHDivBound2<TPZShapeTriang>, TPZHDIVBOUNDTRIANGLEID>;

template<>
int TPZCompElHDivBound2<TPZShapeQuad>::ClassId() const
{
	return TPZHDIVBOUND2QUADID;
}

template class
TPZRestoreClass< TPZCompElHDivBound2<TPZShapeQuad>, TPZHDIVBOUNDQUADID>;

/*
 const int TPZHDIVBOUNDPOINTID = 252;
 const int TPZHDIVBOUNDLINEARID = 252;
 const int TPZHDIVBOUNDTRIANGLEID = 252;
 const int TPZHDIVBOUNDQUADID = 252;
 */

template class TPZCompElHDivBound2<TPZShapeTriang>;
template class TPZCompElHDivBound2<TPZShapePoint>;
template class TPZCompElHDivBound2<TPZShapeLinear>;
template class TPZCompElHDivBound2<TPZShapeQuad>;
//template class TPZCompElHDiv<TPZShapeTetra>;
//template class TPZCompElHDiv<TPZShapePrism>;
//template class TPZCompElHDiv<TPZShapePiram>;
//template class TPZCompElHDiv<TPZShapeCube>;
