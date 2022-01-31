/**
 * @file
 * @brief Contains the implementation of the TPZCompElHDivPressure methods.
 */

#include "pzcmesh.h"
#include "pzhdivpressure.h"
#include "pzquad.h"
#include "pzgeoel.h"
#include "TPZMaterial.h"
#include "TPZMatSingleSpace.h"
#include "pzlog.h"
#include "pzgeoquad.h"
#include "TPZShapeDisc.h"
#include "TPZCompElDisc.h"
#include "TPZMaterialDataT.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.TPZCompElHDivPressure");
#endif

using namespace std;


// TESTADO
template<class TSHAPE>
TPZCompElHDivPressure<TSHAPE>::TPZCompElHDivPressure(TPZCompMesh &mesh, TPZGeoEl *gel, const HDivFamily hdivfam) :
TPZRegisterClassId(&TPZCompElHDivPressure::ClassId),
TPZCompElHDiv<TSHAPE>(mesh,gel,hdivfam) {
		
		if (TSHAPE::Type()==EQuadrilateral) {
				fPressureOrder = mesh.GetDefaultOrder();
		}
		else {
				fPressureOrder = mesh.GetDefaultOrder()-1;
		}

		
		//fPressureOrder = mesh.GetDefaultOrder()-1;
		
		
		//criando o connect da variavel dual
		int nshape = 0;
		
		if (TSHAPE::Type()==EQuadrilateral) {
				nshape =  pzshape::TPZShapeDisc::NShapeF(this->fPressureOrder, this->Dimension(), pzshape::TPZShapeDisc::  ETensorial);
		}
		if (TSHAPE::Type()==ETriangle) {
				nshape =  pzshape::TPZShapeDisc::NShapeF(this->fPressureOrder, this->Dimension(), pzshape::TPZShapeDisc::  EOrdemTotal);
		}
        int nstate = 1;
		this->fConnectIndexes.Resize(NConnects());
		int64_t newnodeindex = mesh.AllocateNewConnect(nshape,nstate,this->fPressureOrder);
		TPZConnect &newnod = mesh.ConnectVec()[newnodeindex];
        newnod.SetLagrangeMultiplier(1);
		this->fConnectIndexes[this->NConnects()-1]=newnodeindex;
		int64_t seqnum = newnod.SequenceNumber();
        newnod.SetLagrangeMultiplier(1);
        mesh.Block().Set(seqnum,nshape);
        mesh.ConnectVec()[this->fConnectIndexes[this->NConnects()-1]].IncrementElConnected();
		
		//		for (int i=0; i<this->NConnects(); i++) {
		//#ifdef PZ_LOG
		//				{
		//						std::stringstream sout;
		//						sout << "verificando  fConnectIndexes " <<  std::endl;
		//						sout<< " i "<<i<<"fConnectIndexes[i] "<<this->fConnectIndexes[i]<<std::endl;
		//						LOGPZ_DEBUG(logger,sout.str())
		//				}
		//#endif
		//		}
		
#ifdef PZ_LOG
    if(logger.isDebugEnabled())
		{
				std::stringstream sout;
				sout << "After creating last pressure connect " << newnodeindex << std::endl;
				this->Print(sout);
				LOGPZ_DEBUG(logger,sout.str())
		}
#endif
		
		//		for (int i=0; i<NConnects(); i++) {
		//				std::cout<< "fConnectIndexes[i] "<<this->fConnectIndexes[i]<<std::endl;
		//		}
		
}


template<class TSHAPE>
TPZCompElHDivPressure<TSHAPE>::TPZCompElHDivPressure(TPZCompMesh &mesh, const TPZCompElHDivPressure<TSHAPE> &copy) :
TPZRegisterClassId(&TPZCompElHDivPressure::ClassId),
TPZCompElHDiv<TSHAPE>(mesh,copy)
{
		fPressureOrder = copy.fPressureOrder;
		
}


template<class TSHAPE>
TPZCompElHDivPressure<TSHAPE>::TPZCompElHDivPressure(TPZCompMesh &mesh, const TPZCompElHDivPressure<TSHAPE> &copy,
													 std::map<int64_t,int64_t> & gl2lcConMap, std::map<int64_t,int64_t> & gl2lcElMap) :
		TPZRegisterClassId(&TPZCompElHDivPressure::ClassId),
TPZCompElHDiv<TSHAPE>(mesh,copy,gl2lcConMap,gl2lcElMap)
{
		
		fPressureOrder = copy.fPressureOrder;
}

template<class TSHAPE>
TPZCompElHDivPressure<TSHAPE>::TPZCompElHDivPressure() :
TPZRegisterClassId(&TPZCompElHDivPressure::ClassId),
TPZCompElHDiv<TSHAPE>()
{
		fPressureOrder = -1;
}

template<class TSHAPE>
TPZCompElHDivPressure<TSHAPE>::~TPZCompElHDivPressure() {
    TPZGeoEl *gel = this->Reference();
    if (gel) {
        TPZCompEl *cel = gel->Reference();
        if (cel) {
            this->RemoveSideRestraintsII(TPZInterpolatedElement::EDelete);
        }
        this->Reference()->ResetReference();
    }

}

template<class TSHAPE>
MElementType TPZCompElHDivPressure<TSHAPE>::Type() {
		return TSHAPE::Type();
}

template<class TSHAPE>
void TPZCompElHDivPressure<TSHAPE>::SetPressureOrder(int order){
		fPressureOrder = order;
#ifdef PZ_LOG
		{
				std::stringstream sout;
				sout << endl<<"Ordem da Variavel dual: "<< fPressureOrder<<std::endl;
				LOGPZ_DEBUG(logger,sout.str())
		}
#endif
		
}

template<class TSHAPE>
int TPZCompElHDivPressure<TSHAPE>::DualOrder() {
		return fPressureOrder;
}

template<class TSHAPE>
int TPZCompElHDivPressure<TSHAPE>::NConnects() const {
		
    return TPZCompElHDiv<TSHAPE>::NConnects()+1;//acrescentando um connect mais pra variavel dual
}

template<class TSHAPE>
void TPZCompElHDivPressure<TSHAPE>::SetConnectIndex(int i, int64_t connectindex){
#ifndef PZNODEBUG
		if(i<0 || i>= this->NConnects()) {
				std::cout << " TPZCompElHDivPressure<TSHAPE>::SetConnectIndex index " << i <<
				" out of range\n";
				DebugStop();
				return;
		}
#endif
		this-> fConnectIndexes[i] = connectindex;
#ifdef PZ_LOG
    if(logger.isDebugEnabled())
		{
				std::stringstream sout;
				sout << endl<<"Setting Connect : " << i << " to connectindex " << connectindex<<std::endl;
				LOGPZ_DEBUG(logger,sout.str())
		}
#endif
}

template<class TSHAPE>
int TPZCompElHDivPressure<TSHAPE>::NConnectShapeF(int connect, int order)const
{
    int dualorder=this->fPressureOrder;
    if(connect< NConnects()-1)
    {//tirando o connect da pressao
        int numshape=TPZCompElHDiv<TSHAPE>::NConnectShapeF(connect,order);
        return numshape;   
    }
		
		
    else {
        if(TSHAPE::Type()==EQuadrilateral)
        {
						
            int numshape=pzshape::TPZShapeDisc::NShapeF(dualorder, this->Dimension(), pzshape::TPZShapeDisc:: ETensorial);
            return(numshape);
        }
        else// (TSHAPE::Type()==ETriangle)
        {
            return pzshape::TPZShapeDisc::NShapeF(dualorder, this->Dimension(), pzshape::TPZShapeDisc:: EOrdemTotal);
        }
    }
}


//Identifies the interpolation order on the connects of the element
template<class TSHAPE>
void TPZCompElHDivPressure<TSHAPE>::GetInterpolationOrder(TPZVec<int> &ord) {
		ord.Resize(NConnects());
		int i;
		for(i=0; i<NConnects(); i++) {
				ord[i] = ConnectOrder(i);
		}
}

template<class TSHAPE>
int64_t TPZCompElHDivPressure<TSHAPE>::ConnectIndex(int con) const{
//#ifndef PZNODEBUG
//		if(con<0 || con>= this->NConnects()) {
//				std::cout << "TPZCompElHDivPressure::ConnectIndex wrong parameter connect " << con <<
//				" NConnects " << this-> NConnects() << std::endl;
//				DebugStop();
//				return -1;
//		}
//		
//#endif
		
				if(con<0 ) {
						std::cout << "TPZCompElHDivPressure::ConnectIndex wrong parameter connect " << con <<
						" NConnects " << this-> NConnects() << std::endl;
						DebugStop();
						return -1;
				}
				else{
						if ( con>= this->NConnects()) {
								int con2= con-TSHAPE::NCornerNodes;
								return this->fConnectIndexes[con2];
				
				
						}
		
						else{
		
		
								return this->fConnectIndexes[con];
				
						}
				}
}



//Sets the preferred interpolation order along a side
template<class TSHAPE>
void TPZCompElHDivPressure<TSHAPE>::SetPreferredOrder(int order)
{
		this->fPreferredOrder = order;
}


/**
 return the interpolation orderof the polynomial for connect
 **/
template<class TSHAPE>
int TPZCompElHDivPressure<TSHAPE>::ConnectOrder(int connect) const
{
		if (connect < NConnects() - 1) {
        return TPZCompElHDiv<TSHAPE>::ConnectOrder(connect);
    }
    else {
				return (this->fPressureOrder);//definindo ordem de interpolacao para o connect da pressao
		}
}


//compute the values of the shape function of the side
template<class TSHAPE>
void TPZCompElHDivPressure<TSHAPE>::SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
		
		// this is an exception
		// when the sides parameter is "out of bounds", the method appends the dualshape functions
		// and assumes the point coordinates are already referring to the interior of the element
		// I hate exceptions...
		if( side == TSHAPE::NSides)
		{
				int nshapedual = NConnectShapeF(NConnects()-1,fPressureOrder);
				TPZFNMatrix<300> phi1(nshapedual,1),phi2(phi);
				TPZFNMatrix<900> dphi1(TSHAPE::Dimension,nshapedual),dphi2(dphi);
				ShapeDual(point,phi1,dphi1);
				Append(phi2, phi1, phi);
				Append(dphi2, dphi1, dphi);
				return;
		}
    else TPZCompElHDiv<TSHAPE>::SideShapeFunction(side,point,phi,dphi);
		
}

template<class TSHAPE>
template<class TVar>
void TPZCompElHDivPressure<TSHAPE>::SolutionT(TPZVec<REAL> &qsi,int var,TPZVec<TVar> &sol){
		
        TPZMaterialDataT<TVar> data;
		InitMaterialData(data);
        data.p=this->MaxOrder();
        
        this->ComputeShape(qsi,data);
        this->ComputeSolutionHDivPressureT(data);
        
        data.x.Resize(3,0.);
        this->Reference()->X(qsi,data.x);
        auto *mat =
            dynamic_cast<TPZMatSingleSpaceT<TVar>*>(this->Material());
		mat->Solution(data,var,sol);
}

template<class TSHAPE>
template<class TVar>
void TPZCompElHDivPressure<TSHAPE>::ComputeSolutionHDivPressureT(TPZMaterialDataT<TVar> &data){
		
    const int numdof = this->Material()->NStateVariables();
    const int ncon = this->NConnects();
    
    
    TPZBlock &block =this->Mesh()->Block();
    TPZFMatrix<TVar> &MeshSol = this->Mesh()->Solution();
    int64_t numbersol = MeshSol.Cols();
    
    int nsol= this->Dimension()+2;
    data.sol.Resize(numbersol);
    data.dsol.Resize(numbersol);
    for (int64_t is=0; is<numbersol; is++) {
        data.sol[is].Resize(nsol,1);//2 componente para fluxo+ 1 para pressao +1 para div
        data.sol[is].Fill(0);
				
    }
		//solucao associada a fluxo
    int iv = 0,ishape=0,ivec=0,cols, jv=0;
    for(int in=0; in<ncon-1 ; in++) {//estou tirando o connect da pressao
				TPZConnect *df = &this->Connect(in);
				int64_t dfseq = df->SequenceNumber();
				int dfvar = block.Size(dfseq);
				int64_t pos = block.Position(dfseq);
				
				for(int jn=0; jn<dfvar; jn++) {
						ivec=data.fVecShapeIndex[jv ].first;
						ishape=data.fVecShapeIndex[jv].second;
						
						TPZFNMatrix<3> ivecDiv(3,1);
						ivecDiv(0,0) = data.fDeformedDirections(0,ivec);
						ivecDiv(1,0) = data.fDeformedDirections(1,ivec);
						ivecDiv(2,0) = data.fDeformedDirections(2,ivec);
						TPZFNMatrix<3> axesvec(3,1);
						data.axes.Multiply(ivecDiv,axesvec);
						
						for (int ilinha=0; ilinha<this->Dimension(); ilinha++) {
								cols=iv%numdof;
								
								//	 #ifdef PZ_LOG
								//	 std::stringstream sout;
								//	 sout << " vetor  " << ivec << " shape  " << ishape<<" coef "<< MeshSol(pos+jn,0)<<endl;
								//	 LOGPZ_DEBUG(logger,sout.str())
								//	 #endif
								for (int64_t is=0; is<numbersol; is++) {
                    data.sol[is][ilinha] += data.fDeformedDirections(ilinha,ivec)* data.phi(ishape,0)*MeshSol(pos+jn,is);
                    data.sol[is][nsol-1] +=  axesvec(ilinha,0)*data.dphix(ilinha,ishape)*MeshSol(pos+jn,is);//divergente
                    
                }
								
						}
						
						jv++;
				}
				
				iv++;
		}
		
		
    //colocando a solucao para o connect interno usando shape descontinua
    
    TPZConnect *df2 = &this->Connect(ncon-1);
    int64_t dfseq2 = df2->SequenceNumber();
    int64_t pos2 = block.Position(dfseq2);
    
    for (int64_t idesc=0; idesc<data.numberdualfunctions; idesc++) {
				int iphi= data.phi.Rows()-data.numberdualfunctions +idesc;
        for (int64_t is=0; is<numbersol; is++) {
            data.sol[is][nsol-2]+= data.phi(iphi,0)*MeshSol(pos2+idesc,is);            
        }
				
    }
    
}


template<class TSHAPE>
void TPZCompElHDivPressure<TSHAPE>::Append(TPZFMatrix<REAL> &u1, TPZFMatrix<REAL> &u2, TPZFMatrix<REAL> &u12)
{
		bool Is_u1PHI = (u1.Cols() == 1) ? true : false;
		bool Is_u2PHI = (u2.Cols() == 1) ? true : false;
		
		if(Is_u1PHI && Is_u2PHI)
		{
				int64_t nu1 = u1.Rows(),nu2 = u2.Rows();
				u12.Redim(nu1+nu2,1);
				int64_t i;
				for(i=0; i<nu1; i++) u12(i,0) = u1(i,0);
				for(i=0; i<nu2; i++) u12(i+nu1,0) = u2(i,0);
				
				
		}
		else if(!Is_u1PHI || !Is_u2PHI) // Se u1 e u2 nÃ£o sÃ£o Phi's, implica em serem dPhi's
		{
				int64_t ru1 = u1.Rows(), cu1 = u1.Cols(), ru2 = u2.Rows(), cu2 = u2.Cols();
				int64_t ru12 = ru1 < ru2 ? ru2 : ru1;
				int64_t cu12 = cu1+cu2;
				u12.Redim(ru12,cu12);
				int64_t i,j;
				for(i=0; i<ru1; i++) for(j=0; j<cu1; j++) u12(i,j) = u1(i,j);
				for(i=0; i<ru2; i++) for(j=0; j<cu2; j++) u12(i,j+cu1) = u2(i,j);//---modifiquei--
		}
		else
		{
				PZError << "TPZCompElHDivPressure::Append. Bad input parameters " << std::endl;//Este metodo so serve para u1 E u2 do mesmo tipo
				
				
		}
		
}

/** compute the shape functions corresponding to the dual space
 *
 */
template<class TSHAPE>
void TPZCompElHDivPressure<TSHAPE>::ShapeDual(TPZVec<REAL> &qsi, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi)
{
		int dimension= TSHAPE::Dimension;
		REAL C=1;//fator de escala utilizado neste metodo
		TPZManVector<REAL,3> X0(3,0.);//centro do elemento
		
    int degree= this->fPressureOrder;
		//	const int nshapedisc = pzshape::TPZShapeDisc::NShapeF(degree, dimension, pzshape::TPZShapeDisc:: ETensorial);
		pzshape::TPZShapeDisc::Shape(dimension,C,X0,qsi,degree,phi,dphi, pzshape::TPZShapeDisc:: ETensorial);
		
}

template<class TSHAPE>
void TPZCompElHDivPressure<TSHAPE>::Shape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi) {
    
    TPZFMatrix<REAL> phiCont;
  	TPZFMatrix<REAL> dphiCont;
    TPZCompElHDiv<TSHAPE>::Shape(pt,phiCont,dphiCont);
    
    // acrescentar as funcoes de pressao (via descontinuo)
    REAL C=1;//fator de escala utilizado neste metodo
    TPZManVector<REAL,3> X0(3,0.);//centro do elemento
		
    int dimension= TSHAPE::Dimension;
    int degree= this->fPressureOrder;
    
    int nshapedisc;
    nshapedisc = pzshape::TPZShapeDisc::NShapeF(degree, dimension, pzshape::TPZShapeDisc:: ETensorial);
    TPZFNMatrix<660> phiDisc(nshapedisc,1);
    TPZFNMatrix<660> dphiDisc(dimension,nshapedisc);
    
    if (TSHAPE::Type()==EQuadrilateral){
        pzshape::TPZShapeDisc::Shape(dimension,C,X0,pt,degree,phiDisc,dphiDisc, pzshape::TPZShapeDisc::ETensorial);
        
//#ifdef PZ_LOG
//        std::stringstream sout;
//        sout<<"\n ponto de integracao " << pt <<endl;
//        sout<< "\n vetor phiDisc "<<phiDisc<<endl;
//        LOGPZ_DEBUG(logger,sout.str())
//#endif

        
    }
    else if (TSHAPE::Type()==ETriangle) {
        pzshape::TPZShapeDisc::Shape(dimension,C,X0,pt,degree,phiDisc,dphiDisc, pzshape::TPZShapeDisc::EOrdemTotal);
    }
    else {
        DebugStop();
    }

    
    this->Append(phiCont,phiDisc,phi);
    this->Append(dphiCont,dphiDisc,dphi);
//	{	
//    #ifdef PZ_LOG
//    std::stringstream sout;
//   // sout<< "vetor phiCont"<<phiCont<<endl;
//    sout<< "\n vetor phi "<<phi<<endl;
//    LOGPZ_DEBUG(logger,sout.str())
//    #endif
//    }
}

template<class TSHAPE>
TPZTransform<> TPZCompElHDivPressure<TSHAPE>::TransformSideToElement(int side){
		return TSHAPE::TransformSideToElement(side);
}


/** Initialize a material data and its attributes based on element dimension, number
 * of state variables and material definitions
 */

template<class TSHAPE>
void TPZCompElHDivPressure<TSHAPE>::InitMaterialData(TPZMaterialData &data)
{
		
		//#ifdef PZ_LOG
		//	{
		//		LOGPZ_DEBUG(logger,"Initializing normal vectors")
		//	}
		//#endif
    data.fShapeType = TPZMaterialData::EVecandShape;
    TPZCompElHDiv<TSHAPE>::InitMaterialData(data);
    data.numberdualfunctions = NConnectShapeF(NConnects()-1,fPressureOrder);
}

///  Save the element data to a stream
template<class TSHAPE>
void TPZCompElHDivPressure<TSHAPE>::Write(TPZStream &buf, int withclassid) const
{
		TPZInterpolatedElement::Write(buf,withclassid);
		TPZManVector<int,3> order(3,0);
		this-> fIntRule.GetOrder(order);
		buf.Write(order);
		buf.Write(this->fConnectIndexes.begin(),TSHAPE::NSides);
		buf.Write(&this->fPreferredOrder,1);
        buf.Write(&this->fPressureOrder);
		int classid = this->ClassId();
		buf.Write ( &classid, 1 );
}


//Read the element data from a stream

template<class TSHAPE>
void TPZCompElHDivPressure<TSHAPE>::Read(TPZStream &buf, void *context)
{
		TPZInterpolatedElement::Read(buf,context);
		TPZManVector<int,3> order;
		buf.Read(order);
		this-> fIntRule.SetOrder(order);
		buf.Read(this->fConnectIndexes.begin(),TSHAPE::NSides);
		buf.Read(&this->fPreferredOrder,1);
        buf.Read(&this->fPressureOrder);
		int classid = -1;
		buf.Read( &classid, 1 );
		if ( classid != this->ClassId() )
		{
				std::stringstream sout;
				sout << "ERROR - " << __PRETTY_FUNCTION__
        << " trying to restore an object id " << this->ClassId() << " for an package of id = " << classid;
				LOGPZ_ERROR ( logger, sout.str().c_str() );
		}
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


template<class TSHAPE>
void TPZCompElHDivPressure<TSHAPE>::CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) {
		if(dimension == TSHAPE::Dimension && this->Material()->Id() > 0) {
				new typename TSHAPE::GraphElType(this,&grafgrid);
		}
}

//template class
//TPZRestoreClass< TPZCompElHDivPressure<TPZShapePoint>>;

#ifndef BORLAND
template class TPZRestoreClass< TPZCompElHDivPressure<TPZShapeLinear>>;
template class TPZRestoreClass< TPZCompElHDivPressure<TPZShapeTriang>>;
template class TPZRestoreClass< TPZCompElHDivPressure<TPZShapeQuad>>;
template class TPZRestoreClass< TPZCompElHDivPressure<TPZShapeCube>>;
template class TPZRestoreClass< TPZCompElHDivPressure<TPZShapeTetra>>;
template class TPZRestoreClass< TPZCompElHDivPressure<TPZShapePrism>>;
template class TPZRestoreClass< TPZCompElHDivPressure<TPZShapePiram>>;
#endif


template class TPZCompElHDivPressure<TPZShapeTriang>;
//template class TPZCompElHDivPressure<TPZShapePoint>;
template class TPZCompElHDivPressure<TPZShapeLinear>;
template class TPZCompElHDivPressure<TPZShapeQuad>;
template class TPZCompElHDivPressure<TPZShapeTetra>;
template class TPZCompElHDivPressure<TPZShapePrism>;
template class TPZCompElHDivPressure<TPZShapePiram>;
template class TPZCompElHDivPressure<TPZShapeCube>;


//TPZCompEl * CreateHDivPressurePointEl(TPZGeoEl *gel,TPZCompMesh &mesh) {
//		return new TPZCompElHDivPressure<TPZShapePoint>(mesh,gel);
//}


TPZCompEl * CreateHDivPressureLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
		return new TPZCompElHDivBound2< TPZShapeLinear>(mesh,gel,hdivfam);
}

TPZCompEl * CreateHDivPressureQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
		return new TPZCompElHDivPressure< TPZShapeQuad>(mesh,gel,hdivfam);
}

TPZCompEl * CreateHDivPressureTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
		return new TPZCompElHDivPressure< TPZShapeTriang >(mesh,gel,hdivfam);
}

TPZCompEl * CreateHDivPressureCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
		return new TPZCompElHDivPressure< TPZShapeCube >(mesh,gel,hdivfam);
}

TPZCompEl * CreateHDivPressurePrismEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
		return new TPZCompElHDivPressure< TPZShapePrism>(mesh,gel,hdivfam);
}

TPZCompEl * CreateHDivPressurePyramEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
		return new TPZCompElHDivPressure< TPZShapePiram >(mesh,gel,hdivfam);
}

TPZCompEl * CreateHDivPressureTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
		return new TPZCompElHDivPressure< TPZShapeTetra >(mesh,gel,hdivfam);
}

