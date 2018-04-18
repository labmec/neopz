/*
 *  pzhdivfull.cpp
 *  PZ
 *
 *  Created by labmec on 10/23/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#include "pzhdivfull.h"
#include "pzquad.h"
#include "pzgeoel.h"
#include "TPZMaterial.h"
#include "pzlog.h"
#include "pzgeoquad.h"


#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.TPZCompElHDivFull"));
#endif
template<class TSHAPE>

TPZCompElHDivFull<TSHAPE>::TPZCompElHDivFull(TPZCompMesh &mesh, TPZGeoEl *gel, int64_t &index) :
TPZRegisterClassId(&TPZCompElHDivFull::ClassId),
TPZCompElHDiv<TSHAPE>(mesh,gel,index) {
	
	int i;
	int nconflux= TPZCompElHDiv<TSHAPE>::NConnects();
    this->fConnectIndexes.Resize(nconflux);
	gel->SetReference(this);
	
	for(i=0;i< nconflux;i++)
	{
        int sideaux= i + TSHAPE::NCornerNodes;
		this->fConnectIndexes[i] = this->CreateMidSideConnect(sideaux);
#ifdef LOG4CXX
        if (logger->isDebugEnabled())
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
    
	
	
	
}



template<class TSHAPE>
TPZCompElHDivFull<TSHAPE>::TPZCompElHDivFull() :
TPZRegisterClassId(&TPZCompElHDivFull::ClassId),
TPZCompElHDiv<TSHAPE>()
{
	
}

template<class TSHAPE>
TPZCompElHDivFull<TSHAPE>::~TPZCompElHDivFull(){
	
}

template<class TSHAPE>
MElementType TPZCompElHDivFull<TSHAPE>::Type() {
	return TSHAPE::Type();
}

template<class TSHAPE>
void TPZCompElHDivFull<TSHAPE>::Write(TPZStream &buf, int withclassid) const
{
	
}


//Read the element data from a stream

template<class TSHAPE>
void TPZCompElHDivFull<TSHAPE>::Read(TPZStream &buf, void *context)
{
	
}

template<class TSHAPE>
void TPZCompElHDivFull<TSHAPE>::SetSideOrder(int side, int order){
	int connectaux= TPZCompElHDiv<TSHAPE>:: SideConnectLocId(0,side);
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
    c.SetNState(nvar);
    int nshape =this-> NConnectShapeF(connectaux,order);
    c.SetNShape(nshape);
	this-> Mesh()->Block().Set(seqnum,nshape*nvar);
}


template<class TSHAPE>
int TPZCompElHDivFull<TSHAPE>::NConnectShapeF(int connect, int order)const
{
	//TPZCompElHDiv<TSHAPE>::NConnectShapeF(connect);
	
	if(connect >= this->NConnects())
	{
		PZError << "TPZCompElHDivFull<TSHAPE>::NConnectShapeF: there is not this connect " <<  std::endl;
		return -1;
	}
	
	if (TSHAPE::Type()==EQuadrilateral)
	{
		int iside = connect+TSHAPE::NCornerNodes;
		if(TSHAPE::SideDimension(iside)< this->Dimension()-2)
		{
			PZError << "TPZCompElHDiv<TSHAPE>::NConnectShapeF: no shape associate " <<  std::endl;
			return -1;
			
		}
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
					NShapeFace += TSHAPE::NConnectShapeF(nside,order);
				}
				return(NShapeFace);
			}
			else{
				int NShapeFace = 0;
				int Nshape=0;
				for(int nside=TSHAPE::NCornerNodes; nside<smallsides.NElements();nside++)
				{
					NShapeFace += TSHAPE::NConnectShapeF(nside,order);
				}
				
				Nshape= NShapeFace + 2*(order-1)*(order-1);
				return(Nshape);

				
				//int nshape=(order-1)*(order-1)+(order)*(order) + 4*(order)-1;
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
		
		if(connect< TPZCompElHDiv<TSHAPE>::NConnects()){
			int iside = connect+TSHAPE::NCornerNodes;
			if(TSHAPE::SideDimension(iside)< this->Dimension()-2)
			{
				PZError << "TPZCompElHDivFull<TSHAPE>::NConnectShapeF: no shape associate " <<  std::endl;
				return -1;
				
			}
			
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
void TPZCompElHDivFull<TSHAPE>::IndexShapeToVec(TPZVec<int> &VectorSide,TPZVec<std::pair<int,int64_t> > & ShapeAndVec, int pressureorder){
	//		{	
	//		#ifdef LOG4CXX
	//												std::stringstream sout;
	//												sout << "VectorSide "<<VectorSide << std::endl;
	//												LOGPZ_DEBUG(logger,sout.str())
	//		#endif
	//		}
	
#ifdef LOG4CXX
	std::stringstream sout;
	sout << " ShapeAndVecQuadIn " << ShapeAndVec;
	LOGPZ_DEBUG(logger,sout.str())
#endif
	
    // VectorSide indicates the side associated with each vector entry
    TPZManVector<int64_t,27> FirstIndex;
    // the first index of the shape functions
    FirstShapeIndex(FirstIndex);
	
    int64_t count=0;
    int nshapeflux= this-> NFluxShapeF();
	
    ShapeAndVec.Resize(nshapeflux);    
   if (TSHAPE::Type()==EQuadrilateral) {
        
        TPZManVector<int64_t,4> ids(4,0);
        TPZGeoEl *gel = this->Reference();
       
        for (int id=0; id<4; id++) {
            ids[id] = gel->NodePtr(id)->Id();
        }
        
        for(int64_t jvec=0;jvec< VectorSide.NElements();jvec++)
        {
            if (jvec==2||jvec==5||jvec==8||jvec==11)
			{
				int lside=VectorSide[jvec];
                int64_t fshape1= FirstIndex[lside];
                int nconside=this-> SideConnectLocId(0,lside);
                TPZConnect &c = this->Connect(nconside);
                int nshapecon= c.NShape();
                int64_t fshape2= fshape1+nshapecon-2;
                for (int64_t ishape=fshape1; ishape<fshape2; ishape++)
                {
#ifdef LOG4CXX
                    std::stringstream sout;
                    sout << " <vec,shape> " << "< "<<jvec << " * "<<ishape << "> "<<std::endl;
                    LOGPZ_DEBUG(logger,sout.str())
#endif
                    ShapeAndVec[count++]=std::pair<int,int64_t>(jvec,ishape);
                }
            }
			else
            {
                int lside=VectorSide[jvec];
                int64_t fshape1= FirstIndex[lside];
                int64_t fshape2= FirstIndex[lside+1];
                for (int64_t ishape=fshape1; ishape<fshape2; ishape++)
                {
#ifdef LOG4CXX
                   std::stringstream sout;
                   sout << " <vec,shape> " << "< "<<jvec << " * "<<ishape << "> "<<std::endl;
                    LOGPZ_DEBUG(logger,sout.str())
#endif
                    ShapeAndVec[count++]=std::pair<int,int64_t>(jvec,ishape);
                }
            }
			
		}
	   
#ifdef LOG4CXX
	   std::stringstream sout;
	   sout << " ShapeAndVecQuad " << ShapeAndVec;
	   LOGPZ_DEBUG(logger,sout.str())
#endif
				
   }//end to EQuadrilateral 
	

	
    else {
        int64_t count=0;
        for(int64_t jvec=0;jvec< VectorSide.NElements();jvec++)
        {
            if (jvec==2||jvec==5||jvec==8)
            {
                int lside=VectorSide[jvec];
                int64_t fshape1= FirstIndex[lside];
                int nconside=this-> SideConnectLocId(0,lside);
                TPZConnect &c = this->Connect(nconside);
                int nshapecon=c.NShape();
                int64_t fshape2= fshape1+nshapecon-2;
                for (int64_t ishape=fshape1; ishape<fshape2; ishape++)
                {
					//#ifdef LOG4CXX
					//                    std::stringstream sout;
					//                    sout << " <vec,shape> " << "< "<<jvec << " * "<<ishape << "> "<<std::endl;
					//                    LOGPZ_DEBUG(logger,sout.str())
					//#endif
                    ShapeAndVec[count++]=std::pair<int,int64_t>(jvec,ishape);
                }
            }
            else
            {
                int lside=VectorSide[jvec];
                int64_t fshape1= FirstIndex[lside];
                int64_t fshape2= FirstIndex[lside+1];
                for (int64_t ishape=fshape1; ishape<fshape2; ishape++)
                {
					//#ifdef LOG4CXX
					//                   std::stringstream sout;
					//                   sout << " <vec,shape> " << "< "<<jvec << " * "<<ishape << "> "<<std::endl;
					//                    LOGPZ_DEBUG(logger,sout.str())
					//#endif
                    ShapeAndVec[count++]=std::pair<int,int64_t>(jvec,ishape);
                }
            }
        }
    }
	
    
//#ifdef LOG4CXX
//    std::stringstream sout;
//    sout << " ShapeAndVec " << ShapeAndVec;
//    LOGPZ_DEBUG(logger,sout.str())
//#endif
    
}

template<class TSHAPE>
int TPZCompElHDivFull<TSHAPE>::NFluxShapeF() const{
    int in,result=0;
    int nn=this->NConnects();
    for(in=0;in<nn;in++){
        TPZConnect &c = this->Connect(in);
        result += this-> NConnectShapeF(in,c.Order());
    }
	
	
    return result;
	
	
}


template<class TSHAPE>
void TPZCompElHDivFull<TSHAPE>::FirstShapeIndex(TPZVec<int64_t> &Index){
	
	Index.Resize(TSHAPE::NSides+1);
	Index[0]=0;
	int maxorder=0;
	
	int ncon=this->NConnects();
	for (int icon=0; icon< ncon; icon++) {
		
		maxorder=(this->ConnectOrder(icon) > maxorder) ? this->ConnectOrder(icon) : maxorder;
		
	}
	
	
	for(int iside=0;iside<TSHAPE::NSides;iside++)
	{
		Index[iside+1] = Index[iside] + TSHAPE::NConnectShapeF(iside,maxorder);		
		
		
	}
	
#ifdef LOG4CXX
    std::stringstream sout;
    sout << "First  Index " << Index;
    LOGPZ_DEBUG(logger,sout.str())
#endif
}

template<class TSHAPE>
void TPZCompElHDivFull<TSHAPE>::Shape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi) {
	
	TPZManVector<int64_t,TSHAPE::NCornerNodes> id(TSHAPE::NCornerNodes,0);
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
void TPZCompElHDivFull<TSHAPE>::NShapeContinuous(TPZVec<int> &order, int &nshape ){
	
	int ncon=this->NConnects();
	order.Resize(ncon);
	int maxorder=0;
	//int ordercon=0;
	for (int icon=0; icon< ncon; icon++) {
		
		maxorder=(this-> ConnectOrder(icon) > maxorder) ? this->ConnectOrder(icon) : maxorder;
		
	}
	
	for (int i=0; i< ncon; i++) {
		
			order[i]=maxorder;
		
		
	}
	
	
	nshape=TSHAPE::NShapeF(order);
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
		{
				std::stringstream sout;
				sout << "ordem max "<<maxorder<< " vec order " << order<<" num func cont "<< nshape<<std::endl;
				LOGPZ_DEBUG(logger,sout.str())
		}
#endif
	
	
}



template<class TSHAPE>
void TPZCompElHDivFull<TSHAPE>::InitMaterialData(TPZMaterialData &data)
{
	#ifdef LOG4CXX
	{
		LOGPZ_DEBUG(logger,"Initializing MaterialData of TPZCompElHDivFull")
	}
#endif
	
	TPZIntelGen<TSHAPE>::InitMaterialData(data);

	
#ifdef LOG4CXX
	{
		LOGPZ_DEBUG(logger,"Initializing MaterialData of TPZCompElHDivFull")
	}
#endif
	TPZVec<int> normalsides;
	TPZIntelGen<TSHAPE>::Reference()->ComputeNormals(data.fNormalVec, normalsides);
	
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

#include "pzelchdivbound2.h"

using namespace pzgeom;
using namespace pzshape;

//template<>
//void TPZCompElHDivFull<TPZShapePoint>::CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) {
//	if(dimension == 0) std::cout << "A point element has no graphical representation\n";
//}

template<class TSHAPE>
void TPZCompElHDivFull<TSHAPE>::CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) {
	if(dimension == TSHAPE::Dimension && this->Material()->Id() > 0) {
		new typename TSHAPE::GraphElType(this,&grafgrid);
	}
}



//------

//template class
//TPZRestoreClass< TPZCompElHDivFull<pzshape::TPZShapePoint>>;

#ifndef BORLAND
template class TPZRestoreClass< TPZCompElHDivFull<pzshape::TPZShapeLinear>>;
template class TPZRestoreClass< TPZCompElHDivFull<pzshape::TPZShapeTriang>>;
template class TPZRestoreClass< TPZCompElHDivFull<pzshape::TPZShapeQuad>>;
template class TPZRestoreClass< TPZCompElHDivFull<pzshape::TPZShapeCube>>;
template class TPZRestoreClass< TPZCompElHDivFull<pzshape::TPZShapeTetra>>;
template class TPZRestoreClass< TPZCompElHDivFull<pzshape::TPZShapePrism>>;
template class TPZRestoreClass< TPZCompElHDivFull<pzshape::TPZShapePiram>>;
#endif

template class TPZCompElHDivFull<pzshape::TPZShapeTriang>;
//template class TPZCompElHDivFull<pzshape::TPZShapePoint>;
template class TPZCompElHDivFull<pzshape::TPZShapeLinear>;
template class TPZCompElHDivFull<pzshape::TPZShapeQuad>;
template class TPZCompElHDivFull<pzshape::TPZShapeTetra>;
template class TPZCompElHDivFull<pzshape::TPZShapePrism>;
template class TPZCompElHDivFull<pzshape::TPZShapePiram>;
template class TPZCompElHDivFull<pzshape::TPZShapeCube>;


//TPZCompEl * CreateHDivFullPointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
//    return new TPZCompElHDivBound2<pzshape::TPZShapePoint>(mesh,gel,index);
//}


TPZCompEl * CreateHDivFullLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
    return new TPZCompElHDivFull< pzshape::TPZShapeLinear>(mesh,gel,index);
}

TPZCompEl * CreateHDivFullQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
    return new TPZCompElHDivFull< pzshape::TPZShapeQuad>(mesh,gel,index);
}

TPZCompEl * CreateHDivFullBoundTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
    return new TPZCompElHDivBound2< pzshape::TPZShapeTriang >(mesh,gel,index);
}

TPZCompEl * CreateHDivFullBoundLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
    return new TPZCompElHDivBound2< pzshape::TPZShapeLinear>(mesh,gel,index);
}

TPZCompEl * CreateHDivFullBoundQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
    return new TPZCompElHDivBound2< pzshape::TPZShapeQuad>(mesh,gel,index);
}

TPZCompEl * CreateHDivFullTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
    return new TPZCompElHDiv< pzshape::TPZShapeTriang >(mesh,gel,index);
}

TPZCompEl * CreateHDivFullCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
    return new TPZCompElHDiv< pzshape::TPZShapeCube >(mesh,gel,index);
}

TPZCompEl * CreateHDivFullPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
    return new TPZCompElHDiv< pzshape::TPZShapePrism>(mesh,gel,index);
}

TPZCompEl * CreateHDivFullPyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
    return new TPZCompElHDiv< pzshape::TPZShapePiram >(mesh,gel,index);
}

TPZCompEl * CreateHDivFullTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
    return new TPZCompElHDiv< pzshape::TPZShapeTetra >(mesh,gel,index);
}

