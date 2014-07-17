//
//  pzcompelhdivcurved.cpp
//  PZ
//
//  Created by Douglas Castro on 6/25/14.
//
//

#include "pzcompelhdivcurved.h"

/**
 * @file
 * @brief Contains the implementation of the TPZCompElHdivCurved methods.
 */

#include "pzcmesh.h"
#include "pzelchdiv.h"
#include "pzquad.h"
#include "pzgeoel.h"
#include "pzmaterial.h"
#include "pzlog.h"
#include "pzgeoquad.h"
#include "TPZShapeDisc.h"
#include "TPZCompElDisc.h"
#include "pzmaterialdata.h"
#include "pzhdivpressure.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.TPZCompElHdivCurved"));
#endif

using namespace std;


template<class TSHAPE>
TPZCompElHdivCurved<TSHAPE>::TPZCompElHdivCurved(TPZCompMesh &mesh, TPZGeoEl *gel, long &index) :
TPZCompElHDiv<TSHAPE>(mesh,gel,index) {
	int i;
	this->TPZInterpolationSpace::fPreferredOrder = mesh.GetDefaultOrder();
	int nconflux= TPZCompElHDiv<TSHAPE>::NConnects();
    this->fConnectIndexes.Resize(nconflux);
	gel->SetReference(this);
	
	for(i=0;i< nconflux;i++)
	{
        int sideaux= i + TSHAPE::NCornerNodes;
		this->fConnectIndexes[i] = this->CreateMidSideConnect(sideaux);
#ifdef LOG4CXX
        {
            std::stringstream sout;
            sout << "After creating last flux connect " << i << std::endl;
            //	this->Print(sout);
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif
        
		mesh.ConnectVec()[this->fConnectIndexes[i]].IncrementElConnected();
		this->IdentifySideOrder(sideaux);
    }
    
	
    int sideorder = this->SideOrder(TSHAPE::NSides-1);
    if(TSHAPE::Type()==EQuadrilateral)
    {
        sideorder++;
    }
	sideorder = 2*sideorder;
	if (sideorder > this->fIntRule.GetMaxOrder()) sideorder = this->fIntRule.GetMaxOrder();
	TPZManVector<int,3> order(3,sideorder);
	this->fIntRule.SetOrder(order);
}

template<class TSHAPE>
TPZCompElHdivCurved<TSHAPE>::TPZCompElHdivCurved() :
TPZCompElHDiv<TSHAPE>()
{

}

template<class TSHAPE>
TPZCompElHdivCurved<TSHAPE>::~TPZCompElHdivCurved(){
	
}

template<class TSHAPE>
MElementType TPZCompElHdivCurved<TSHAPE>::Type() {
	return TSHAPE::Type();
}

template<class TSHAPE>
int TPZCompElHdivCurved<TSHAPE>::NConnectShapeF(int connect)const
{
	
    if(connect >= this->NConnects())
    {
        PZError << "TPZCompElHdivCurved<TSHAPE>::NConnectShapeF: there is not this connect " <<  endl;
        return -1;
    }
    
    if (TSHAPE::Type()==EQuadrilateral)
    {
        int iside = connect+TSHAPE::NCornerNodes;
        if(TSHAPE::SideDimension(iside)< this->Dimension()-2)
        {
            PZError << "TPZCompElHdivCurved<TSHAPE>::NConnectShapeF: no shape associate " <<  endl;
            return -1;
            
        }
        int order = this->ConnectOrder(connect);
        if(order < 0) return 0;
        
        TPZStack<int> smallsides;
        TSHAPE::LowerDimensionSides(iside,smallsides);
        //i.e., trata-se do lado de mesma dimensao que o elemento
        if(TSHAPE::SideDimension(iside) == TSHAPE::Dimension){
            
            ///------
            if (order==1) {
                int NShapeFace = 0;
                // funcoes para os lados menor que o proprio elemento
                for(int nside=TSHAPE::NCornerNodes; nside<smallsides.NElements();nside++)
                {
                    NShapeFace += TSHAPE::NConnectShapeF(nside,order+1);
                }
                return(NShapeFace);
            }
            else{
                int nshape=(order-1)*(order-1)+(order)*(order) + 4*(order)-1;
                return(nshape);
            }
            
            
            
        }
        
        if(TSHAPE::SideDimension(iside) == TSHAPE::Dimension-1)
        {
            int NShapeF = 0;
            for(int j=0;j< smallsides.NElements();j++)
            {
                NShapeF += TSHAPE::NConnectShapeF(j,order);
            }
            
            int result=NShapeF + TSHAPE::NConnectShapeF(iside,order);
            return(result);
        }
        
    }
    
    if (TSHAPE::Type()==ETriangle){
        
        if(connect< this->NConnects()){
            int iside = connect+TSHAPE::NCornerNodes;
            if(TSHAPE::SideDimension(iside)< this->Dimension()-2)
            {
                PZError << "TPZCompElHdivCurved<TSHAPE>::NConnectShapeF: no shape associate " <<  endl;
                return -1;
                
            }
            int order = this->ConnectOrder(connect);
            
            if(order < 0) return 0;
            
            TPZStack<int> smallsides;
            TSHAPE::LowerDimensionSides(iside,smallsides);
            if(TSHAPE::SideDimension(iside) == TSHAPE::Dimension)//i.e., trata-se do lado de mesma dimensao que o elemento
            {
                
                int result= order*order-1;
                
                
                return (result)  ;
                
            }
            else if(TSHAPE::SideDimension(iside) == TSHAPE::Dimension-1)//i.e., trata-se do lado de 1 dimensao a menos que a dimensao do elemento
            {
                int NShapeF = 0;
                for(int j=0;j< smallsides.NElements();j++)
                {
                    NShapeF += TSHAPE::NConnectShapeF(j,order);
                }
                int result=NShapeF + TSHAPE::NConnectShapeF(iside,order);
                return (result);
            }
        }
	}
    
    
    else {
        std::cout<< "No implemented yet"<<std::endl;
        DebugStop();
    }
    
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout <<__PRETTY_FUNCTION__<< "unhandled case ";
        LOGPZ_ERROR(logger,sout.str())
    }
#endif
    return -1;
    
}


template<class TSHAPE>
void TPZCompElHdivCurved<TSHAPE>::SetSideOrder(int side, int order){
	int connectaux= TPZCompElHDiv<TSHAPE>::SideConnectLocId(0,side);
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
    c.SetNState(nvar);
    int nshape =this-> NConnectShapeF(connectaux);
    c.SetNShape(nshape);
	this-> Mesh()->Block().Set(seqnum,nshape*nvar);
}


/**return the first shape associate to each side*/
template<class TSHAPE>
void TPZCompElHdivCurved<TSHAPE>::FirstShapeIndex(TPZVec<long> &Index)
{
	
    Index.Resize(TSHAPE::NSides+1);
    Index[0]=0;
    int maxorder=0;
	
    int ncon=this->NConnects();
    for (int icon=0; icon< ncon; icon++) {
        
        maxorder=(this->ConnectOrder(icon) > maxorder) ? this->ConnectOrder(icon) : maxorder;
        
    }
    
    
	for(int iside=0;iside<TSHAPE::NSides;iside++)
	{
		//int order= SideOrder(iside);
		//Index[iside+1] = Index[iside] + TSHAPE::NConnectShapeF(iside,order);
        if (TSHAPE::Type()==EQuadrilateral) {
            Index[iside+1] = Index[iside] + TSHAPE::NConnectShapeF(iside,maxorder+1);
        }
        else {
            Index[iside+1] = Index[iside] + TSHAPE::NConnectShapeF(iside,maxorder);
        }
        
        
		
	}
	
#ifdef LOG4CXX
    std::stringstream sout;
    sout << "First  Index " << Index;
    LOGPZ_DEBUG(logger,sout.str())
#endif
}


template<class TSHAPE>
int TPZCompElHdivCurved<TSHAPE>::NFluxShapeF() const{
    int in,result=0;
    int nn=this->NConnects();
    for(in=0;in<nn;in++){
        result += this->NConnectShapeF(in);
    }
    

    return result;
    
    
}


template<class TSHAPE>
void TPZCompElHdivCurved<TSHAPE>::IndexShapeToVec(TPZVec<int> &VectorSide,TPZVec<std::pair<int,long> > & ShapeAndVec, int pressureorder){
    //		{
    //		#ifdef LOG4CXX
    //												std::stringstream sout;
    //												sout << "VectorSide "<<VectorSide << std::endl;
    //												LOGPZ_DEBUG(logger,sout.str())
    //		#endif
    //		}
    
    // VectorSide indicates the side associated with each vector entry
    TPZManVector<long,27> FirstIndex;
    // the first index of the shape functions
    FirstShapeIndex(FirstIndex);
    
    long count=0;
    int nshapeflux= this->NFluxShapeF();
    
    ShapeAndVec.Resize(nshapeflux);
    if (TSHAPE::Type()==EQuadrilateral) {
        
        TPZManVector<long,4> ids(4,0);
        TPZGeoEl *gel = this->Reference();
        int id;
        for (id=0; id<4; id++) {
            ids[id] = gel->NodePtr(id)->Id();
        }
        
        for(long jvec=0;jvec< VectorSide.NElements();jvec++)
        {
            if (jvec==2||jvec==5||jvec==8||jvec==11)
            {
                int lside=VectorSide[jvec];
                int nconside=this->SideConnectLocId(0,lside);
                int nshapecon=NConnectShapeF(nconside);
                if (nshapecon > 2)
                {
                    long fshape1= FirstIndex[lside];
                    long fshape2= fshape1+ nshapecon-2;
                    //#ifdef LOG4CXX
                    //                    std::stringstream sout;
                    //                    sout << " fshape1 " <<fshape1 << " fshape2 "<<fshape2 << std::endl;
                    //                    LOGPZ_DEBUG(logger,sout.str())
                    //#endif
                    for (long ishape=fshape1; ishape<fshape2; ishape++)
                    {
                        //#ifdef LOG4CXX
                        //                        std::stringstream sout;
                        //                        sout << " <vec,shape> " << "< "<<jvec << " * "<<ishape << "> "<<std::endl;
                        //                        LOGPZ_DEBUG(logger,sout.str())
                        //#endif
                        ShapeAndVec[count++]=std::pair<int,long>(jvec,ishape);
                    }
                }
            }
            else if(jvec==16 || jvec ==17)
            {
                int lside = VectorSide[jvec];
                
                int conectside=TPZCompElHDiv<TSHAPE>::SideConnectLocId(0,lside);
                int order=this->ConnectOrder(conectside);
                int nshape = TSHAPE::NConnectShapeF(lside,order+1);//order-1);
                TPZFNMatrix<25> sideorders(2,nshape);
                int ksi,eta;
                int loccount = 0;
                int transid = TSHAPE::GetTransformId(ids);
                switch (transid)
                {
                    case 0:
                    case 3:
                    case 4:
                    case 7:
                        for (ksi=0; ksi<order/*-1*/; ksi++) {
                            for (eta=0; eta<order/*-1*/; eta++) {
                                sideorders(0,loccount) = ksi+2;
                                sideorders(1,loccount) = eta+2;
                                loccount++;
                            }
                        }
                        break;
                    case 2:
                    case 6:
                    case 1:
                    case 5:
                        for (eta=0; eta<order/*-1*/; eta++) {
                            for (ksi=0; ksi<order/*-1*/; ksi++) {
                                sideorders(0,loccount) = ksi+2;
                                sideorders(1,loccount) = eta+2;
                                loccount++;
                            }
                        }
                        
                        break;
                    default:
                        DebugStop();
                        break;
                }
                long ish=0;
                long fshape1 = FirstIndex[lside];
                //#ifdef LOG4CXX
                //                {
                //                    std::stringstream sout;
                //                    sideorders.Print("SideOrders= ", sout ,EFormatted);
                //                    LOGPZ_DEBUG(logger,sout.str())
                //                }
                //#endif
                
                
                for (ish=0; ish<nshape; ish++)
                {
                    int orderksi = sideorders(0,ish);
                    int ordereta = sideorders(1,ish);
                    
                    if (jvec ==16)
                    {
                        bool etacheck = ordereta <= pressureorder;
                        if (etacheck)
                        {
                            if (!(ordereta == pressureorder+1 && orderksi == pressureorder))
                            {
                                ShapeAndVec[count++]=std::pair<int,long>(jvec,fshape1+ish);
                                
                                //#ifdef LOG4CXX
                                //                                std::stringstream sout;
                                //                                sout << " <vec,shape> " << "< "<<jvec << " * "<<fshape1+ish << "> "<<std::endl;
                                //                                sout << " side order ksi " << sideorders(0,ish) << " side order eta " << sideorders(1,ish);
                                //                                LOGPZ_DEBUG(logger,sout.str())
                                //#endif
                            }
                        }
                    }
                    if (jvec ==17) {
                        if (orderksi<=pressureorder) {
                            if (!(orderksi == pressureorder+1 && ordereta == pressureorder))
                            {
                                //#ifdef LOG4CXX
                                //                                std::stringstream sout;
                                //                                sout << " <vec,shape> " << "< "<<jvec << " * "<<fshape1+ish << "> "<<std::endl;
                                //                                sout << " side order ksi " << sideorders(0,ish) << " side order eta " << sideorders(1,ish);
                                //                                LOGPZ_DEBUG(logger,sout.str())
                                //#endif
                                ShapeAndVec[count++]=std::pair<int,long>(jvec,fshape1+ish);
                            }
                        }
                    }
                }
            }
            else {
                int lside = VectorSide[jvec];
                long fshape1 = FirstIndex[lside];
                long fshape2 = FirstIndex[lside+1];
                //#ifdef LOG4CXX
                //                std::stringstream sout;
                //                sout << " lside "<< lside << " fshape1 " <<fshape1 << " fshape2 "<<fshape2 << std::endl;
                //                LOGPZ_DEBUG(logger,sout.str())
                //#endif
                for (long ishape=fshape1; ishape<fshape2; ishape++) {
                    //#ifdef LOG4CXX
                    //                    std::stringstream sout;
                    //                    sout << " <vec,shape> " << "< "<<jvec << " * "<<ishape << "> "<<std::endl;
                    //                    LOGPZ_DEBUG(logger,sout.str())
                    //#endif
                    ShapeAndVec[count++]=std::pair<int,long>(jvec,ishape);
                }
            }
        }
    }//end to EQuadrilateral
    else {
        long count=0;
        for(long jvec=0;jvec< VectorSide.NElements();jvec++)
        {
            if (jvec==2||jvec==5||jvec==8)
            {
                int lside=VectorSide[jvec];
                long fshape1= FirstIndex[lside];
                int nconside= this->SideConnectLocId(0,lside);
                int nshapecon=NConnectShapeF(nconside);
                long fshape2= fshape1+nshapecon-2;
                for (long ishape=fshape1; ishape<fshape2; ishape++)
                {
                    //#ifdef LOG4CXX
                    //                    std::stringstream sout;
                    //                    sout << " <vec,shape> " << "< "<<jvec << " * "<<ishape << "> "<<std::endl;
                    //                    LOGPZ_DEBUG(logger,sout.str())
                    //#endif
                    ShapeAndVec[count++]=std::pair<int,long>(jvec,ishape);
                }
            }
            else
            {
                int lside=VectorSide[jvec];
                long fshape1= FirstIndex[lside];
                long fshape2= FirstIndex[lside+1];
                for (long ishape=fshape1; ishape<fshape2; ishape++)
                {
                    //#ifdef LOG4CXX
                    //                   std::stringstream sout;
                    //                   sout << " <vec,shape> " << "< "<<jvec << " * "<<ishape << "> "<<std::endl;
                    //                    LOGPZ_DEBUG(logger,sout.str())
                    //#endif
                    ShapeAndVec[count++]=std::pair<int,long>(jvec,ishape);
                }
            }
        }
    }
    
    
#ifdef LOG4CXX
    std::stringstream sout;
    sout << " ShapeAndVec " << ShapeAndVec;
    LOGPZ_DEBUG(logger,sout.str())
#endif
    
}


template<class TSHAPE>

void TPZCompElHdivCurved<TSHAPE>::Shape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi) {
	TPZManVector<long,TSHAPE::NCornerNodes> id(TSHAPE::NCornerNodes,0);
	TPZManVector<int, TSHAPE::NSides-TSHAPE::NCornerNodes+1> ord(TSHAPE::NSides-TSHAPE::NCornerNodes,0);
    int i;
    TPZGeoEl *ref = this->Reference();
    for(i=0; i<TSHAPE::NCornerNodes; i++) {
		id[i] = ref->NodePtr(i)->Id();
    }
    
    int nconflux=this->NConnects();
    for(i=0; i< nconflux; i++)
    {
        ord[i] = this->ConnectOrder(i);
        
    }
    int dimension= TSHAPE::Dimension;
    
    
    int nshape=0;
    this->NShapeContinuous(ord, nshape );
    
    phi.Resize(nshape, 1);
    dphi.Resize(dimension, nshape);
    TSHAPE::Shape(pt,id,ord,phi,dphi);
    
}

template<class TSHAPE>
void TPZCompElHdivCurved<TSHAPE>::NShapeContinuous(TPZVec<int> &order, int &nshape ){
    
    int ncon=this->NConnects();
    order.Resize(ncon);
    int maxorder=0;
    //int ordercon=0;
    for (int icon=0; icon< ncon; icon++) {
        
        maxorder=(this->ConnectOrder(icon) > maxorder) ? this->ConnectOrder(icon) : maxorder;
        
    }
    
    std::cout << "necessario mudar?" << std::endl;
    
// Parece que temos que mudar isso //
    
    
    for (int i=0; i< ncon; i++) {
		if (TSHAPE::Type()==EQuadrilateral) {
            order[i]=maxorder+1;
		}
		else {
            order[i]=maxorder;
		}
        
		
    }
    
    
    nshape=TSHAPE::NShapeF(order);
    //#ifdef LOG4CXX
    //		{
    //				std::stringstream sout;
    //				sout << "ordem max "<<maxorder<< " vec order " << order<<" num func cont "<< nshape<<std::endl;
    //				LOGPZ_DEBUG(logger,sout.str())
    //		}
    //#endif
    
    
}




// Tentar usar isso par obter os vetores, normais ou tangenes, depende do caso


//template<class TSHAPE>
//TPZTransform TPZCompElHdivCurved<TSHAPE>::TransformSideToElement(int side){
//	return TSHAPE::TransformSideToElement(side);
//}


// Save the element data to a stream
template<class TSHAPE>
void TPZCompElHdivCurved<TSHAPE>::Write(TPZStream &buf, int withclassid)
{

}


// Read the element data from a stream
template<class TSHAPE>
void TPZCompElHdivCurved<TSHAPE>::Read(TPZStream &buf, void *context)
{
    
}

/** Initialize a material data and its attributes based on element dimension, number
 * of state variables and material definitions
 */

// Esse eu devo implementar
template<class TSHAPE>
void TPZCompElHdivCurved<TSHAPE>::InitMaterialData(TPZMaterialData &data)
{
	TPZIntelGen<TSHAPE>::InitMaterialData(data);
    //	if (TSHAPE::Type()==EQuadrilateral) {
    //        int maxorder = this->MaxOrder();
    //        data.p = maxorder+1;
    //    }
#ifdef LOG4CXX
    {
        LOGPZ_DEBUG(logger,"Initializing MaterialData of TPZCompElHDiv")
    }
#endif
	TPZVec<int> normalsides;
	TPZIntelGen<TSHAPE>::Reference()->ComputeNormals(data.fNormalVec, normalsides);
    
	// vecindex : lado associado a cada normal
	// vecindex indica apenas o numero do lado associado a cada normal
	// agora temos que expandir para formar pares : vecIndex e shapeindex
	//ComputeShapeIndex(data.fVecIndex,data.fVecShapeIndex);
	//data.numberdualfunctions = NConnectShapeF(NConnects()-1);
    
    int pressureorder=0;
    if (TSHAPE::Type()==EQuadrilateral) {
        pressureorder=this->fPreferredOrder;//ver como melhorar..?
    }
    else {
        pressureorder=this->fPreferredOrder-1;
    }
    
    IndexShapeToVec(normalsides,data.fVecShapeIndex,pressureorder);
    
#ifdef LOG4CXX
	{
		std::stringstream sout;
		data.fNormalVec.Print("Normal vector ", sout);
		sout << "NormalVector/Shape indexes " << data.fVecShapeIndex << std::endl;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
}

#include "pzshapecube.h"
#include "TPZRefCube.h"
#include "pzshapelinear.h"
#include "TPZRefLinear.h"
#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "pzshapetriang.h"
#include "pzreftriangle.h"
#include "pzgeotriangle.h"
#include "pzshapeprism.h"
#include "pzrefprism.h"
#include "pzgeoprism.h"
#include "pzshapetetra.h"
#include "pzreftetrahedra.h"
#include "pzgeotetrahedra.h"
#include "pzshapepiram.h"
#include "pzrefpyram.h"
#include "pzgeopyramid.h"
#include "pzrefpoint.h"
#include "pzgeopoint.h"
#include "pzshapepoint.h"
#include "pzgraphelq2dd.h"
#include "tpzgraphelt3d.h"
#include "pzgraphel1dd.h"
#include "pztrigraphd.h"
#include "pzgraphelq3dd.h"
#include "tpzgraphelprismmapped.h"
#include "tpzgraphelpyramidmapped.h"
#include "tpzgraphelt2dmapped.h"

using namespace pztopology;

#include "tpzpoint.h"
#include "tpzline.h"
#include "tpzquadrilateral.h"
#include "tpztriangle.h"
#include "tpzcube.h"
#include "tpztetrahedron.h"
#include "tpzprism.h"
#include "tpzpyramid.h"

#include "pzmeshid.h"

#include "pzelchdivbound2.h"

using namespace pzgeom;
using namespace pzshape;

template<>
void TPZCompElHdivCurved<TPZShapePoint>::CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) {
	if(dimension == 0) std::cout << "A point element has no graphical representation\n";
}

template<class TSHAPE>
void TPZCompElHdivCurved<TSHAPE>::CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) {
	if(dimension == TSHAPE::Dimension && this->Material()->Id() > 0) {
		new typename TSHAPE::GraphElType(this,&grafgrid);
	}
}

// Provavelmente te que mudar esses IDs

template<>
int TPZCompElHdivCurved<TPZShapePoint>::ClassId() const
{
	return TPZHDIVPOINTID;
}

template<>
int TPZCompElHdivCurved<TPZShapeLinear>::ClassId() const
{
	return TPZHDIVLINEARID;
}
template<>
int TPZCompElHdivCurved<TPZShapeTriang>::ClassId() const
{
	return TPZHDIVTRIANGLEID;
}
template<>
int TPZCompElHdivCurved<TPZShapeQuad>::ClassId() const
{
	return TPZHDIVQUADID;
}
template<>
int TPZCompElHdivCurved<TPZShapeCube>::ClassId() const
{
	return TPZHDIVCUBEID;
}
template<>
int TPZCompElHdivCurved<TPZShapeTetra>::ClassId() const
{
	return TPZHDIVTETRAID;
}
template<>
int TPZCompElHdivCurved<TPZShapePrism>::ClassId() const
{
	return TPZHDIVPRISMID;
}
template<>
int TPZCompElHdivCurved<TPZShapePiram>::ClassId() const
{
	return TPZHDIVPYRAMID;
}
/*
template class
TPZRestoreClass< TPZCompElHdivCurved<TPZShapePoint>, TPZHDIVPOINT_CURVED_ID>;

template class
TPZRestoreClass< TPZCompElHdivCurved<TPZShapeLinear>, TPZHDIVLINEAR_CURVED_ID>;

template class
TPZRestoreClass< TPZCompElHdivCurved<TPZShapeTriang>, TPZHDIVTRIANGLE_CURVED_ID>;

template class
TPZRestoreClass< TPZCompElHdivCurved<TPZShapeQuad>, TPZHDIVQUAD_CURVED_ID>;

template class
TPZRestoreClass< TPZCompElHdivCurved<TPZShapeCube>, TPZHDIVCUBE_CURVED_ID>;

template class
TPZRestoreClass< TPZCompElHdivCurved<TPZShapeTetra>, TPZHDIVTETRA_CURVED_ID>;

template class
TPZRestoreClass< TPZCompElHdivCurved<TPZShapePrism>, TPZHDIVPRISM_CURVED_ID>;

template class
TPZRestoreClass< TPZCompElHdivCurved<TPZShapePiram>, TPZHDIVPYRAM_CURVED_ID>;


*/

template class TPZCompElHdivCurved<TPZShapeTriang>;
template class TPZCompElHdivCurved<TPZShapePoint>;
template class TPZCompElHdivCurved<TPZShapeLinear>;
template class TPZCompElHdivCurved<TPZShapeQuad>;
template class TPZCompElHdivCurved<TPZShapeTetra>;
template class TPZCompElHdivCurved<TPZShapePrism>;
template class TPZCompElHdivCurved<TPZShapePiram>;
template class TPZCompElHdivCurved<TPZShapeCube>;


TPZCompEl * CreateHDivCurvedBoundPointEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
	return new TPZCompElHDivBound2< pzshape::TPZShapePoint>(mesh,gel,index);
}

TPZCompEl * CreateHDivCurvedLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
	return new TPZCompElHdivCurved< pzshape::TPZShapeLinear>(mesh,gel,index);
}

TPZCompEl * CreateHDivCurvedBoundLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
	return new TPZCompElHDivBound2< pzshape::TPZShapeLinear>(mesh,gel,index);
}

TPZCompEl * CreateHDivCurvedQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
	return new TPZCompElHdivCurved< pzshape::TPZShapeQuad>(mesh,gel,index);
}

TPZCompEl * CreateHDivCurvedTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
	return new TPZCompElHdivCurved< pzshape::TPZShapeTriang >(mesh,gel,index);
}

TPZCompEl * CreateHDivCurvedBoundQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
	return new TPZCompElHDivBound2< pzshape::TPZShapeQuad>(mesh,gel,index);
}

TPZCompEl * CreateHDivCurvedBoundTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
	return new TPZCompElHDivBound2< pzshape::TPZShapeTriang >(mesh,gel,index);
}

TPZCompEl * CreateHDivCurvedCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
	return new TPZCompElHdivCurved< pzshape::TPZShapeCube >(mesh,gel,index);
}

TPZCompEl * CreateHDivCurvedPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
	return new TPZCompElHdivCurved< pzshape::TPZShapePrism>(mesh,gel,index);
}

TPZCompEl * CreateHDivCurvedPyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
	return new TPZCompElHdivCurved< pzshape::TPZShapePiram >(mesh,gel,index);
}

TPZCompEl * CreateHDivCurvedTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
	return new TPZCompElHdivCurved< pzshape::TPZShapeTetra >(mesh,gel,index);
}

