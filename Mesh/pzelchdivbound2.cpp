/**
 * @file
 * @brief Contains the implementation of the TPZCompElHDivBound2 methods.
 */


#include "pzelchdivbound2.h"
#include "pzgeoel.h"
#include "pzmaterial.h"
#include "pzlog.h"
#include "TPZShapeDisc.h"
#include "TPZCompElDisc.h"
#include "pzmaterialdata.h"
#include "pzelchdiv.h"


#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.TPZCompElHDivBound2"));
#endif

template<class TSHAPE>
TPZCompElHDivBound2<TSHAPE>::TPZCompElHDivBound2(TPZCompMesh &mesh, TPZGeoEl *gel, long &index) :
TPZIntelGen<TSHAPE>(mesh,gel,index,1), fSideOrient(1){
		
	//int i;
	this->TPZInterpolationSpace::fPreferredOrder = mesh.GetDefaultOrder();
	//for(i=0; i<TSHAPE::NSides; i++) this->fConnectIndexes[i]=-1;
		this->fConnectIndexes[0]=-1;
	gel->SetReference(this);
    TPZIntelGen<TSHAPE>::fConnectIndexes.resize(1);
		
    this->fConnectIndexes[0] = this->CreateMidSideConnect(TSHAPE::NSides-1);
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
		{
				std::stringstream sout;
				sout << "After creating boundary flux connect " << this->fConnectIndexes[0] << std::endl;
				//	this->Print(sout);
				LOGPZ_DEBUG(logger,sout.str())
		}
#endif
    
    mesh.ConnectVec()[this->fConnectIndexes[0]].IncrementElConnected();
    
		
	//TPZGeoElSide myInnerSide(gel,gel->NSides()-1);
//	TPZGeoElSide neigh = myInnerSide.Neighbour();
//	while(!neigh.Reference())
//	{
//		neigh = neigh.Neighbour();
//	}
//	if(neigh == myInnerSide)
//	{
//		/**
//		 O codigo pressupoe que os elementos computacionais 2D sao criados antes dos 1D.
//		 Quando serao criados os elementos computacionais 1D, os respectivos vizinhos 2D sao encontrados.
//		 Situacoes assim ocorrem (neste algoritmo) quando eh realizado refinamento uniforme, pois os primeiros elementos sem descendentes sao os 2D (e depois os descendentes 1D de contorno)
//		 
//		 Ocorreu o problema quando tentou-se realizar o refinamento do quadrilatero em 02 triangulos, em que o quadrilatero apresenta descendentes e as arestas nao.
//		 Neste caso a criacao de elementos computacionais eh iniciada pelos 1D, fazendo com que nao encontrem vizinhos computacionais 2D.
//		 Com isso a variavel int connectIndex0 eh setada com -1, dando o BUG observado.
//		 */
//		std::cout << "Nao foi encontrado elemento 2D com elemento computacional inicializado!!!\n"; 
//		DebugStop();
//	}
//	TPZCompElSide compneigh(neigh.Reference());
//    fneighbour = compneigh;
//	int sideoffset = neigh.Element()->NSides()-neigh.Side();
//	int neighnconnects = compneigh.Element()->NConnects();
//	int connectnumber = neighnconnects-sideoffset;
//	TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (compneigh.Element());
//	connectnumber = intel->SideConnectLocId(0,compneigh.Side());
//	int connectIndex0 = compneigh.Element()->ConnectIndex(connectnumber);
//	
//	this->fConnectIndexes[0] = connectIndex0;
//	mesh.ConnectVec()[connectIndex0].IncrementElConnected();
	
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
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
TPZCompElHDivBound2<TSHAPE>::TPZCompElHDivBound2(TPZCompMesh &mesh, const TPZCompElHDivBound2<TSHAPE> &copy) :
TPZIntelGen<TSHAPE>(mesh,copy), fSideOrient(copy.fSideOrient)
{
//	for(int i=0;i<TSHAPE::NSides;i++)
//	{
//		this-> fConnectIndexes[i] = copy.fConnectIndexes[i];
//	}
    long index = copy.fneighbour.Element()->Index();
    TPZCompEl *cel = this->Mesh()->ElementVec()[index];
    if (!cel) {
        DebugStop();
    }
    fneighbour = TPZCompElSide(cel,copy.fneighbour.Side());
}

// NAO TESTADO
template<class TSHAPE>
TPZCompElHDivBound2<TSHAPE>::TPZCompElHDivBound2(TPZCompMesh &mesh,
												 const TPZCompElHDivBound2<TSHAPE> &copy,
												 std::map<long,long> & gl2lcConMap,
												 std::map<long,long> & gl2lcElMap) :
TPZIntelGen<TSHAPE>(mesh,copy,gl2lcConMap,gl2lcElMap), fSideOrient(copy.fSideOrient)
{
	
	this-> fPreferredOrder = copy.fPreferredOrder;
	int i;
	for(i=0;i<TSHAPE::NSides;i++)
	{
		long lcIdx = -1;
		long glIdx = copy.fConnectIndexes[i];
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
    
    long neiIdx = copy.fneighbour.Element()->Index();
    if(gl2lcElMap.find(neiIdx)==gl2lcElMap.end())
    {
        DebugStop();
    }
    TPZCompEl *cel = mesh.ElementVec()[gl2lcElMap[neiIdx]];
    if (!cel) {
        DebugStop();
    }
    fneighbour = TPZCompElSide(cel,copy.fneighbour.Side());
}

// TESTADO
template<class TSHAPE>
TPZCompElHDivBound2<TSHAPE>::TPZCompElHDivBound2() : TPZIntelGen<TSHAPE>(),fneighbour()
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

// NAO TESTADO
template<class TSHAPE>
void TPZCompElHDivBound2<TSHAPE>::SetSideOrient(int sideorient)
{
    //fSideOrient = this->Reference()->NormalOrientation(side);
    fSideOrient = sideorient;
}

// NAO TESTADO
template<class TSHAPE>
int TPZCompElHDivBound2<TSHAPE>::GetSideOrient( )
{
    return fSideOrient;
}

template<class TSHAPE>
int TPZCompElHDivBound2<TSHAPE>::NConnects() const {
	
	return 1;
}

template<class TSHAPE>
void TPZCompElHDivBound2<TSHAPE>::SetConnectIndex(int i, long connectindex)
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
		
		TPZManVector<int,22> order(TSHAPE::NSides-TSHAPE::NCornerNodes,ConnectOrder(connect));
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
#ifdef LOG4CXX
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << " Bad side or order " << side << " order " << order;
		LOGPZ_DEBUG(logger,sout.str())
#endif
		return;
	}
	TPZConnect &c = this->Connect(connectaux);
    c.SetOrder(order);
    long seqnum = c.SequenceNumber();
    int nvar = 1;
    TPZMaterial * mat =this-> Material();
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
    data.fShapeType = TPZMaterialData::EScalarShape;
	
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
	{
		LOGPZ_DEBUG(logger,"Initializing normal vectors")
	}
#endif
	
	//data.fVecShapeIndex=true;
	/*
	TPZGeoElSide gelside(this->Reference(),TSHAPE::NSides-1);
	TPZGeoElSide neighbour = gelside.Neighbour();
	while(gelside != neighbour && neighbour.Element()->Dimension() != TSHAPE::Dimension+1)
	{
		neighbour = neighbour.Neighbour();
	}
	if(neighbour.Element()->Dimension() != TSHAPE::Dimension+1)
	{
		DebugStop();
	}
	TPZGeoEl *neighel = neighbour.Element();
	TPZManVector<int,9> normalsides;
//    TPZFNMatrix<100,REAL> normalvec;
	neighel->ComputeNormals(neighbour.Side(),data.fNormalVec, normalsides);
//#ifdef LOG4CXX
//	{
//		std::stringstream sout;
//		sout << "normal side depois do ComputeNormals " << normalsides << std::endl;
//		LOGPZ_DEBUG(logger,sout.str())
//	}
//#endif
	
	// relate the sides indicated in vecindex to the sides of the current element
	long nvec = normalsides.NElements();
	long ivec;
	for(ivec=0; ivec<nvec; ivec++)
	{
		TPZGeoElSide neigh(neighel,normalsides[ivec]);
//#ifdef LOG4CXX
//		{
//			std::stringstream sout;
//			sout << "normal side depois do TPZGeoElSide " << normalsides << std::endl;
//			LOGPZ_DEBUG(logger,sout.str())
//		}
//#endif
		while(neigh.Element() != this->Reference())
		{
			
			neigh = neigh.Neighbour();
		}
		
		normalsides[ivec]=neigh.Side();
	}
	IndexShapeToVec(normalsides,data.fVecShapeIndex);
	data.numberdualfunctions = 0;
#ifdef LOG4CXX
	{
		std::stringstream sout;
		data.fNormalVec.Print("Normal vectors", sout);
		sout << "Vector/Shape indexes " << data.fVecShapeIndex << std::endl;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	*/
}

template<class TSHAPE>
void TPZCompElHDivBound2<TSHAPE>::ComputeShapeIndex(TPZVec<int> &sides, TPZVec<long> &shapeindex) {
	
	TPZManVector<long> firstshapeindex;    // Para o que?
	FirstShapeIndex(firstshapeindex);      // se foram calculados os indices mas n�o utilizados?
	int nshape = TPZIntelGen<TSHAPE>::NShapeF();
	shapeindex.Resize(nshape);
	long nsides = sides.NElements();
	long is, count=0;
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
void TPZCompElHDivBound2<TSHAPE>::FirstShapeIndex(TPZVec<long> &Index){
	
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
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << " FirsShapeIndex result " << Index;
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
	
	
}

template<class TSHAPE>
void TPZCompElHDivBound2<TSHAPE>::SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
	
	if(TSHAPE::SideDimension(side)!= TSHAPE::Dimension ){
		DebugStop() ;
	}
    TPZGeoEl *gel = this->Reference();
    int nc = gel->NCornerNodes();
    TPZManVector<long,8> id(nc);
    for (int ic=0; ic<nc; ic++) {
        id[ic] = gel->Node(ic).Id();
    }
    TPZManVector<int,TSHAPE::NSides> ord;
    this->GetInterpolationOrder(ord);

    TPZFNMatrix<50,REAL> philoc(phi.Rows(),phi.Cols()),dphiloc(dphi.Rows(),dphi.Cols());
    TSHAPE::Shape(point,id,ord,philoc,dphiloc);
    
    //int idsize = id.size();
    TPZManVector<int,9> permutegather(TSHAPE::NSides);
    int transformid = TSHAPE::GetTransformId(id);
    TSHAPE::GetSideHDivPermutation(transformid, permutegather);
    
    TPZManVector<long,27> FirstIndex(TSHAPE::NSides+1);
    FirstShapeIndex(FirstIndex);
   

    int signQn = fSideOrient;
    int order = this->Connect(0).Order();
    for (int side=0; side < TSHAPE::NSides; side++) {
        int ifirst = FirstIndex[side];
        int kfirst = FirstIndex[permutegather[side]];
        int nshape = TSHAPE::NConnectShapeF(side,order);
        for (int i=0; i<nshape; i++) {
            phi(ifirst+i,0) = signQn*philoc(kfirst+i,0);
            for (int d=0; d< TSHAPE::Dimension; d++) {
                dphi(d,ifirst+i) = signQn*dphiloc(d,kfirst+i);
            }
        }
    }
    
    return;
}

/** Compute the shape function at the integration point */
template<class TSHAPE>
void TPZCompElHDivBound2<TSHAPE>::Shape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi)
{
    TPZManVector<int,TSHAPE::NSides> ordl;
    this->GetInterpolationOrder(ordl);
    int nshape = NConnectShapeF(0);
    phi.Resize(nshape, 1);
    dphi.Resize(TSHAPE::Dimension, nshape);
    SideShapeFunction(TSHAPE::NSides-1, pt, phi, dphi);
    return;
	/*
  TPZCompElSide thisside(this,TSHAPE::NSides-1);
	TPZGeoElSide thisgeoside(thisside.Reference());
	TPZGeoElSide neighgeo(thisgeoside.Neighbour());
	TPZCompElSide neigh(fneighbour);
	TPZInterpolatedElement *neighel = dynamic_cast<TPZInterpolatedElement *> (neigh.Element());
	if(!neighel)
	{
        LOGPZ_ERROR(logger,"Inconsistent neighbour")
        DebugStop();
		return;
	}
	int nshapeneigh = neighel->NSideShapeF(neighgeo.Side())+1;
	phi.Redim(nshapeneigh, 1);
	dphi.Redim(neighgeo.Element()->Dimension(), nshapeneigh);
	TPZTransform tr(thisgeoside.Dimension()),tr2; 
	thisgeoside.SideTransform3(neighgeo, tr);
	TPZManVector<REAL,3> pt2(neighgeo.Dimension()),pt3(neighel->Dimension());
	tr.Apply(pt, pt2);
	neighel->SideShapeFunction(neigh.Side(), pt2, phi, dphi);
	 */
		//tentando reimplementar 
		TPZManVector<long,TSHAPE::NSides-1> id(TSHAPE::NSides-1,0);
		
		TPZGeoEl *ref = this->Reference();
		int nnodes= ref->NNodes();
		for(int i=0; i<nnodes; i++) {
				id[i] = ref->NodePtr(i)->Id();
		}
		
#ifdef LOG4CXX
		{
				std::stringstream sout;
				sout<< "---Id local---"<<id<<std::endl;
				LOGPZ_DEBUG(logger,sout.str())
		}
#endif
		
		//-----ordenando os id's
//		int i, j, min, x;
//		for (i = 0; i < nnodes; ++i) {
//				min = i;
//				for (j = i+1; j < nnodes; ++j){
//						if (id[j] < id[min])  min = j;
//				x = id[i]; 
//				id[i] = id[min]; 
//				id[min] = x;
//				}
//		}	
		
#ifdef LOG4CXX
		{
				std::stringstream sout;
				sout<< "---Id local Reordenado---"<<id<<std::endl;
				LOGPZ_DEBUG(logger,sout.str())
		}
#endif
		
		//-----
		
	
		TPZCompElSide thisside(this,TSHAPE::NSides-1);
		TPZGeoElSide thisgeoside(thisside.Reference());
//		TPZGeoElSide neighgeo(thisgeoside.Neighbour());
//		TPZCompElSide neigh(fneighbour);
		TPZInterpolatedElement *neighel = dynamic_cast<TPZInterpolatedElement *> (thisside.Element());
//		if(!neighel)
//		{
//        LOGPZ_ERROR(logger,"Inconsistent neighbour")
//        DebugStop();
//				return;
//		}
		int nshapeneigh = neighel->NSideShapeF(thisgeoside.Side());
		phi.Redim(nshapeneigh, 1);
		dphi.Redim(thisgeoside.Element()->Dimension(), nshapeneigh);
		TPZVec<int> ord;
		neighel->GetInterpolationOrder(ord);
		TSHAPE::Shape(pt,id,ord,phi,dphi);
//    if(id[0] > id[1])
//    {
//        REAL phival = phi(0,0);
//        phi(0,0) = phi(1,0);
//        phi(1,0) = phival;
//        REAL dphival = dphi(0,0);
//        dphi(0,0) = dphi(0,1);
//        dphi(0,1) = dphival;
//    }
		
		
//		TPZTransform tr(thisgeoside.Dimension()),tr2; 
//		thisgeoside.SideTransform3(thisgeoside, tr);
//		TPZManVector<REAL,3> pt2(thisgeoside.Dimension()),pt3(neighel->Dimension());
//		tr.Apply(pt, pt2);
//		neighel->SideShapeFunction(thisside.Side(), pt2, phi, dphi);
    
#ifdef LOG4CXX
		//if (logger->isDebugEnabled())
		{
				std::stringstream sout;
				sout.precision(20);
				sout<< "---Phi Novo---"<<phi<<std::endl;
				sout<< "---Dphi Novo---"<<dphi<<std::endl;
				LOGPZ_DEBUG(logger,sout.str())
		}
#endif
}

/** Read the element data from a stream */
template<class TSHAPE>
void TPZCompElHDivBound2<TSHAPE>::Read(TPZStream &buf, void *context)
{
	TPZIntelGen<TSHAPE>::Read(buf,context);
    buf.Read(&fSideOrient);
}

/** Save the element data to a stream */
template<class TSHAPE>
void TPZCompElHDivBound2<TSHAPE>::Write(TPZStream &buf, int withclassid)
{
	TPZIntelGen<TSHAPE>::Write(buf,withclassid);
    buf.Write(&fSideOrient);
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

/** Return a matrix with index shape and vector associate to element */
template<class TSHAPE>
void TPZCompElHDivBound2<TSHAPE>::IndexShapeToVec(TPZVec<int> &VectorSide,TPZVec<std::pair<int,long> > & ShapeAndVec){
	
	// VectorSide indicates the side associated with each vector entry
	TPZVec<long> FirstIndex;
	// the first index of the shape functions
	FirstShapeIndex(FirstIndex);
#ifdef LOG4CXX
		{
				std::stringstream sout;
				sout << "FirstIndex of shape functions " << FirstIndex;
				LOGPZ_DEBUG(logger,sout.str())
		}
#endif
	
	int tamanho= this->NShapeF();
	
	ShapeAndVec.Resize(tamanho);
	long count=0;
		//for(int jvec=0;jvec< VectorSide.NElements();jvec++)
	for(long jvec=0;jvec< VectorSide.NElements();jvec++)//coloca-se -1 caso queira reduzir o espaco de fluxo
	{
		int lside=VectorSide[jvec];
		long fshape1= FirstIndex[lside];
		long fshape2= FirstIndex[lside+1];
		for (long ishape=fshape1; ishape<fshape2; ishape++)
		{
			
#ifdef LOG4CXX
			std::stringstream sout;
			sout << " <vec,shape> " << "< "<<jvec << " * "<<ishape << "> "<<std::endl;
			LOGPZ_DEBUG(logger,sout.str())
#endif
			
			
			ShapeAndVec[count++]=std::pair<int,long>(jvec,ishape);
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
/** returns the unique identifier for reading/writing objects to streams */
template<>
int TPZCompElHDivBound2<TPZShapePoint>::ClassId() const
{
	return TPZHDIVBOUND2POINTID;
}
template<>
int TPZCompElHDivBound2<TPZShapeLinear>::ClassId() const
{
	return TPZHDIVBOUND2LINEARID;
}
template<>
int TPZCompElHDivBound2<TPZShapeTriang>::ClassId() const
{
	return TPZHDIVBOUND2TRIANGLEID;
}
template<>
int TPZCompElHDivBound2<TPZShapeQuad>::ClassId() const
{
	return TPZHDIVBOUND2QUADID;
}


#ifndef BORLAND
template class
TPZRestoreClass< TPZCompElHDivBound2<TPZShapePoint>, TPZHDIVBOUND2POINTID>;

template class
TPZRestoreClass< TPZCompElHDivBound2<TPZShapeLinear>, TPZHDIVBOUND2LINEARID>;

template class
TPZRestoreClass< TPZCompElHDivBound2<TPZShapeTriang>, TPZHDIVBOUND2TRIANGLEID>;

template class
TPZRestoreClass< TPZCompElHDivBound2<TPZShapeQuad>, TPZHDIVBOUND2QUADID>;
#endif


template class TPZCompElHDivBound2<TPZShapeTriang>;
template class TPZCompElHDivBound2<TPZShapePoint>;
template class TPZCompElHDivBound2<TPZShapeLinear>;
template class TPZCompElHDivBound2<TPZShapeQuad>;
