// -*- c++ -*-

// $Id: pzelctemp.cpp,v 1.42 2008-11-20 23:30:41 phil Exp $

#include "pzelchdiv.h"
#include "pzquad.h"
#include "pzgeoel.h"
#include "pzmaterial.h"
#include "pzlog.h"
#include "pzgeoquad.h"
#include "TPZShapeDisc.h"
#include "TPZCompElDisc.h"
#include "pzmaterialdata.h"
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.TPZCompElHDiv"));
#endif

using namespace std;
//template<class TSHAPE>
//class TPZCompElHDiv : public TPZInterpolatedElement {

// TESTADO
template<class TSHAPE>
TPZCompElHDiv<TSHAPE>::TPZCompElHDiv(TPZCompMesh &mesh, TPZGeoEl *gel, int &index) :
TPZIntelGen<TSHAPE>(mesh,gel,index,1) {
	int i;
	this->TPZInterpolationSpace::fPreferredOrder = mesh.GetDefaultOrder();
//	fPressureOrder = mesh.GetDefaultOrder();
	fPressureOrder = mesh.GetDefaultOrder()-1;
	for(i=0; i<TSHAPE::NSides; i++) this->fConnectIndexes[i]=-1;
	gel->SetReference(this);
	
	for(i=0;i< this->NConnects()-1;i++) {//tirei o connect da pressao
		int sideaux= i + TSHAPE::NSides - this->NConnects() + 1;//tirei connect da pressao
		
		this->fConnectIndexes[i] = this->CreateMidSideConnect(sideaux);
		/* LOGPZ_DEBUG(logger,"before printing")
		 #ifdef LOG4CXX
		 {
		 std::stringstream sout;
		 sout << "After creating connect " << i << std::endl;
		 this->Print(sout);
		 LOGPZ_DEBUG(logger,sout.str())
		 }
		 #endif
		 */
		mesh.ConnectVec()[this->fConnectIndexes[i]].IncrementElConnected();
		this->IdentifySideOrder(sideaux);
		
	}
	
	//criando o connect da variavel dual
	
	int newnodeindex = mesh.AllocateNewConnect();
	TPZConnect &newnod = mesh.ConnectVec()[newnodeindex];
	this->fConnectIndexes[i]=newnodeindex;
	int seqnum = newnod.SequenceNumber();
	newnod.SetOrder(fPressureOrder);
	int nshape =  pzshape::TPZShapeDisc::NShapeF(this->fPressureOrder, this->Dimension(), pzshape::TPZShapeDisc:: ETensorial);
	mesh.Block().Set(seqnum,nshape);
	mesh.ConnectVec()[this->fConnectIndexes[i]].IncrementElConnected();
	
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "After creating last connect " << i << std::endl;
		this->Print(sout);
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
	
}

// NAO TESTADO
template<class TSHAPE>
TPZCompElHDiv<TSHAPE>::TPZCompElHDiv(TPZCompMesh &mesh, const TPZCompElHDiv<TSHAPE> &copy) :
TPZIntelGen<TSHAPE>(mesh,copy)
{
	this-> fPreferredOrder = copy.fPreferredOrder;
	fPressureOrder = copy.fPressureOrder;
	int i;
	for(i=0;i<TSHAPE::NSides;i++)
	{
		this-> fConnectIndexes[i] = copy.fConnectIndexes[i];
	}
}

// NAO TESTADO
template<class TSHAPE>
TPZCompElHDiv<TSHAPE>::TPZCompElHDiv(TPZCompMesh &mesh,
                                      const TPZCompElHDiv<TSHAPE> &copy,
                                      std::map<int,int> & gl2lcConMap,
                                      std::map<int,int> & gl2lcElMap) :
	TPZIntelGen<TSHAPE>(mesh,copy,gl2lcConMap,gl2lcElMap)
{
//   if (gl2lcElMap.find(copy.fIndex)!= gl2lcElMap.end())
//   {
//     std::stringstream sout;
//     sout << "ERROR in : " << __PRETTY_FUNCTION__
//         << " trying to clone an already cloned element index: " << copy.fIndex;
//     LOGPZ_ERROR(logger, sout.str().c_str());
//     exit(-1);
//  }

 this-> fPreferredOrder = copy.fPreferredOrder;
  fPressureOrder = copy.fPressureOrder;
  int i;
  for(i=0;i<TSHAPE::NSides;i++)
  {
    int lcIdx = -1;
    int glIdx = copy.fConnectIndexes[i];
    if (gl2lcConMap.find(glIdx) != gl2lcConMap.end()) lcIdx = gl2lcConMap[glIdx];
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
TPZCompElHDiv<TSHAPE>::TPZCompElHDiv() :
  TPZIntelGen<TSHAPE>()
{
  this->fPreferredOrder = -1;
  fPressureOrder = -1;
  int i;
  for(i=0;i<TSHAPE::NSides;i++) {
   this-> fConnectIndexes[i] = -1;
  }
}

// TESTADO
template<class TSHAPE>
TPZCompElHDiv<TSHAPE>::~TPZCompElHDiv(){

}

// NAO TESTADO
template<class TSHAPE>
MElementType TPZCompElHDiv<TSHAPE>::Type() {
  return TSHAPE::Type();
}

template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::SetPressureOrder(int order){
	fPressureOrder = order;
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << endl<<"Ordem da Variavel dual: "<< fPressureOrder<<std::endl;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif

}

template<class TSHAPE>
int TPZCompElHDiv<TSHAPE>::DualOrder() {
	return fPressureOrder;
}

template<class TSHAPE>
int TPZCompElHDiv<TSHAPE>::NConnects() const {
	int dimension = Dimension()-1;

	return TSHAPE::NumSides(dimension) + 2;//acrescentando um connect mais pra variavel dual antes era apenas NumSides(dimension) + 1


}

template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::SetConnectIndex(int i, int connectindex){
#ifndef NODEBUG
  if(i<0 || i>= this->NConnects()) {
    std::cout << " TPZCompElHDiv<TSHAPE>::SetConnectIndex index " << i <<
      " out of range\n";
    DebugStop();
    return;
  }
#endif
	 this-> fConnectIndexes[i] = connectindex;
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << endl<<"Setting Connect : " << i << " to connectindex " << connectindex<<std::endl;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
}

template<class TSHAPE>
int TPZCompElHDiv<TSHAPE>::NConnectShapeF(int connect) const
{
	if(connect< NConnects()-1){//tirando o connect da pressao
		int iside = TSHAPE::NSides - NConnects() + connect+1;//é o lado do elemento q estou, acrescentei um connect a mais
		if(TSHAPE::SideDimension(iside)< this->Dimension()-2)
		{
			PZError << "TPZCompElHDiv<TSHAPE>::NConnectShapeF: no shape associate " <<  endl;
			return -1;

		}
		int order = ConnectOrder(connect);

		if(order < 0) return 0;

		TPZStack<int> smallsides;
		TSHAPE::LowerDimensionSides(iside,smallsides);
		if(TSHAPE::SideDimension(iside) == TSHAPE::Dimension)//i.e., trata-se do lado de mesma dimensao que o elemento
		{
			int NShapeFace = 0;
			for(int nside=TSHAPE::NCornerNodes; nside<smallsides.NElements();nside++){
				int sideorder= this->SideOrder(nside);
				NShapeFace += TSHAPE::NConnectShapeF(nside,sideorder);

			}
			return (NShapeFace + TSHAPE::Dimension*TSHAPE::NConnectShapeF(iside,order))  ;

		}
			else if(TSHAPE::SideDimension(iside) == TSHAPE::Dimension-1)//i.e., trata-se do lado de 1 dimensao a menos que a dimensao do elemento
			{
				int NShapeF = 0;
				for(int j=0;j< smallsides.NElements();j++)
				{
					NShapeF += TSHAPE::NConnectShapeF(j,order);
				}
					int result=NShapeF + TSHAPE::NConnectShapeF(iside,order);//TSHAPE::Dimension*TSHAPE::NConnectShapeF(iside,order);
					return (result);
			}
	}

	else
	{
		int dualorder=this->fPressureOrder;
		return pzshape::TPZShapeDisc::NShapeF(dualorder, this->Dimension(), pzshape::TPZShapeDisc::ETensorial);

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
void TPZCompElHDiv<TSHAPE>::SetIntegrationRule(int ord) {
  TPZManVector<int,3> order(TSHAPE::Dimension,ord);
	this->fIntRule.SetOrder(order);
}



template<class TSHAPE>
int TPZCompElHDiv<TSHAPE>::NSideConnects(int side) const{
	if(TSHAPE::SideDimension(side)<= Dimension()-2) return 0;
	if(TSHAPE::SideDimension(side)==Dimension()-1) return 1;
	if(TSHAPE::SideDimension(side)== Dimension()) return this->NConnects()-1;//tirando o connect da varivael dual
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << "Side: " << side <<"unhandled case ";
		LOGPZ_ERROR(logger,sout.str())
	}
#endif
	return -1;

}
/** return the local index for side*/
template<class TSHAPE>
int TPZCompElHDiv<TSHAPE>::SideConnectLocId(int node,int side) const {
 	if(TSHAPE::SideDimension(side)<= TSHAPE::Dimension - 2 || node >= NSideConnects(side)) {
		PZError << "TPZCompElHDiv<TSHAPE>::SideConnectLocId no connect associate " <<  endl;
		return -1;
	}

	return side - TSHAPE::NSides + NConnects() - 1;//tirei o connect da pressao
}


//Identifies the interpolation order on the connects of the element
template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::GetInterpolationOrder(TPZVec<int> &ord) {
  ord.Resize(NConnects());
  int i;
  for(i=0; i<NConnects(); i++) {
    ord[i] = ConnectOrder(i);
  }
}

//return the preferred order of the polynomial along connect connect
template<class TSHAPE>
int TPZCompElHDiv<TSHAPE>::PreferredSideOrder(int side) {
	if(TSHAPE::SideDimension(side) < Dimension()-1)
	{
		PZError << __PRETTY_FUNCTION__ << " side " << side << std::endl;
	}
	int connect= SideConnectLocId(0,side);
	if(connect<0 || connect > NConnects()) {
		PZError << "TPZCompElHDiv<TSHAPE>::PreferredSideOrder no polynomial associate " <<  endl;
		return -1;
	}
	if(connect<NConnects()) {
	  int order =this->fPreferredOrder;
	  return this->AdjustPreferredSideOrder(connect,order);
  }
  PZError << "TPZCompElHDiv<TSHAPE>::PreferredSideOrder called for connect = " << connect << "\n";
  return 0;

}

template<class TSHAPE>
int TPZCompElHDiv<TSHAPE>::ConnectIndex(int con) const{
#ifndef NODEBUG
	if(con<0 || con>= this->NConnects()) {
    std::cout << "TPZCompElHDiv::ConnectIndex wrong parameter connect " << con <<
      " NConnects " << this-> NConnects() << std::endl;
    DebugStop();
    return -1;
	}

#endif

	 return this->fConnectIndexes[con];
 }



//Sets the preferred interpolation order along a side
template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::SetPreferredOrder(int order)
{
  this->fPreferredOrder = order;
}

//sets the interpolation order of side to order
template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::SetSideOrder(int side, int order) {
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
   this-> Mesh()->Block().Set(seqnum,NConnectShapeF(connectaux)*nvar);
    if(connectaux == NConnects()-1) {
      SetIntegrationRule(2*order);

  }
}

/**
 return the interpolation orderof the polynomial for connect
 **/
template<class TSHAPE>
int TPZCompElHDiv<TSHAPE>::ConnectOrder(int connect) const
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
		//DebugStop();
#else
		std::cout << sout.str() << std::endl;
#endif
		return 0;
	}
	if (connect < NConnects() - 1) {
		TPZConnect &c = this-> Connect(connect);
		return c.Order();
	}

	else {
		return (this->fPressureOrder);//definindo ordem de interpolacao para o connect da pressao
	}
}

template<class TSHAPE>
int TPZCompElHDiv<TSHAPE>::SideOrder(int side) const
{
	//se side=connect devolve ConnectOrder
	//int node;
	if(!NSideConnects(side)) return -1;
	int corder =SideConnectLocId(0, side);
	int maxorder = 0;
	int conectaux;
	if(corder>=0 || corder <= NConnects()) return ConnectOrder(corder);
	//caso contrario devolve a maior ordem de interpolacao
	TPZStack< int > high;
	TSHAPE::HigherDimensionSides(side, high);
	int highside= high.NElements();


	for(int j=0;j<highside;j++)//percorro todos os lados de dimensao maior
	{
		conectaux =SideConnectLocId(0, high[j]);
		maxorder = (ConnectOrder(conectaux) > maxorder) ? ConnectOrder(conectaux) : maxorder;
	}

	return maxorder;

}
/**return the first shape associate to each side*/
template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::FirstShapeIndex(TPZVec<int> &Index){

	Index.Resize(TSHAPE::NSides+1);
	Index[0]=0;

	for(int iside=0;iside<TSHAPE::NSides;iside++)
	{
		int order= SideOrder(iside);
		Index[iside+1] = Index[iside] + TSHAPE::NConnectShapeF(iside,order);

	}

#ifdef LOG4CXX
    std::stringstream sout;
    sout << " Index" << Index;
    LOGPZ_DEBUG(logger,sout.str())
#endif


}
//return a matrix with index shape and vector associate to element
template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::IndexShapeToVec(TPZVec<int> &VectorSide,TPZVec<std::pair<int,int> > & ShapeAndVec){

	// VectorSide indicates the side associated with each vector entry
	TPZVec<int> FirstIndex;
	// the first index of the shape functions
	FirstShapeIndex(FirstIndex);
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "FirstIndex of shape functions " << FirstIndex;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
//	int dualorder=this->fPressureOrder;
	int tamanho= this->NShapeF() - NConnectShapeF(NConnects()-1);//estou tirando as funcoes da variavel dual

	ShapeAndVec.Resize(tamanho);
	int count=0;
	for(int jvec=0;jvec< VectorSide.NElements();jvec++)
	{
		int lside=VectorSide[jvec];
		int fshape1= FirstIndex[lside];
		int fshape2= FirstIndex[lside+1];
		for (int ishape=fshape1; ishape<fshape2; ishape++)
		{
			ShapeAndVec[count++]=std::pair<int,int>(jvec,ishape);
		}

	}

#ifdef LOG4CXX
    std::stringstream sout;
    sout << " VecShapeIndex " << ShapeAndVec;
    LOGPZ_DEBUG(logger,sout.str())
#endif


}
//compute the values of the shape function of the side

template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix &phi,TPZFMatrix &dphi) {

	// this is an exception
	// when the sides parameter is "out of bounds", the method appends the dualshape functions
	// and assumes the point coordinates are already referring to the interior of the element
	// I hate exceptions...
	if( side == TSHAPE::NSides)
	{
		int nshapedual = NConnectShapeF(NConnects()-1);
		TPZFNMatrix<300> phi1(nshapedual,1),phi2(phi);
		TPZFNMatrix<900> dphi1(TSHAPE::Dimension,nshapedual),dphi2(dphi);
		ShapeDual(point,phi1,dphi1);
		Append(phi2, phi1, phi);
		Append(dphi2, dphi1, dphi);
		return;
	}
	if(TSHAPE::SideDimension(side)!= TSHAPE::Dimension -1 ){
			return ;
	}
	TPZIntelGen<TSHAPE>::SideShapeFunction(side,point,phi,dphi);
	/*
	TPZFNMatrix<300> phi1,dphi1,phi2,dphi2;
	TPZTransform tr = TSHAPE::SideToSideTransform(side,TSHAPE::NSides-1);
	TPZManVector<REAL,5> interior;
	tr.Apply(point,interior);
	ShapeDual(interior,phi2,dphi2);
	Append(phi1, phi2, phi);
	Append(dphi1, dphi2, dphi);
	 */

}

template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>:: Solution(TPZVec<REAL> &qsi,int var,TPZVec<REAL> &sol){
	TPZMaterialData data;
	InitMaterialData(data);
	this->ComputeShape(qsi, data.x,data.jacobian,data.axes, data.detjac,data.jacinv,data.phi, data.dphix);
	
	
	this->ComputeSolution(data);
	this->Material()->Solution(data,var,sol);
    
}


template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix &phi, TPZFMatrix &dphix,
                                            const TPZFMatrix &axes, TPZVec<REAL> &sol, TPZFMatrix &dsol){
    TPZMaterialData data;
    InitMaterialData(data);
    this->ComputeSolution(data);
    
}


template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::ComputeSolution(TPZMaterialData &data){
	//const int dim = this->Reference()->Dimension();
    const int numdof = this->Material()->NStateVariables();
    const int ncon = this->NConnects();
	
	
    
	TPZBlock &block =this->Mesh()->Block();
    TPZFMatrix &MeshSol = this->Mesh()->Solution();
    
	
	data.sol.Resize(4,1);
    /*
     
     #ifdef LOG4CXX
     {
     
     std::stringstream sout;
     sout << std::endl << "num.de connects " << ncon << std::endl;
     sout << std::endl << "num.de func. dual " << data.numberdualfunctions << std::endl;
     sout << std::endl << "Vec index and Shapes " << data.fVecShapeIndex << std::endl;
     LOGPZ_DEBUG(logger,sout.str());
     }
     #endif
     
     */
    data.sol[0]=0.;
    data.sol[1]=0.;
    data.sol[2]=0.;
    data.sol[3]=0.;
	
	int iv = 0,ishape=0,ivec=0,cols, jv=0;
    for(int in=0; in<ncon-1 ; in++) {//estou tirando o connect da pressao
		TPZConnect *df = &this->Connect(in);
		int dfseq = df->SequenceNumber();
		int dfvar = block.Size(dfseq);
		int pos = block.Position(dfseq);
        /*	#ifdef LOG4CXX
         {
         
         std::stringstream sout;
         sout << std::endl << "Sequence Number do connect " << in << " = " << dfseq << " tamanho do bloco  " << dfseq<< " = "<< dfvar<< std::endl;
         sout << std::endl <<" Numero de func shape do conect  " << in<< " =  " << NConnectShapeF(in)<<std::endl;
         LOGPZ_DEBUG(logger,sout.str());
         }
         #endif
         */
		
		for(int jn=0; jn<dfvar; jn++) {
			ivec=data.fVecShapeIndex[jv ].first;
			ishape=data.fVecShapeIndex[jv].second;
            /*
             #ifdef LOG4CXX
             {
             
             std::stringstream sout;
             sout << " Indice do vec = " << ivec<<" Indice da shape = " << ishape <<  std::endl;
             LOGPZ_DEBUG(logger,sout.str());
             }
             #endif
             */
			
			for (int ilinha=0; ilinha<3; ilinha++) {
				cols=iv%numdof;
                /*		#ifdef LOG4CXX
                 {
                 
                 std::stringstream sout;
                 sout << "solucao=  " << pos+jn <<" valor da solucao =  "<< MeshSol(pos+jn,0)<< std::endl;
                 sout<< " valor da phi = " <<  data.phi(ishape,0)<<" vetor na posicao = " << data.fNormalVec(ilinha,ivec)<< endl;
                 LOGPZ_DEBUG(logger,sout.str());
                 }
                 #endif
                 */
				data.sol[ilinha] += data.fNormalVec(ilinha,ivec)* data.phi(ishape,0)*MeshSol(pos+jn,0);
                /*			#ifdef LOG4CXX
                 {
                 
                 std::stringstream sout;
                 sout <<" produto = " <<  data.fNormalVec(ilinha,ivec)* data.phi(ishape,0)*MeshSol(pos+jn,0)<<endl;
                 LOGPZ_DEBUG(logger,sout.str());
                 }
                 #endif
                 */			
				
			}
			
			jv++;
		}
        
#ifdef LOG4CXX
		{
			std::stringstream sout;
			sout << std::endl << "Checando VecSol para este connect " <<data.sol<< std::endl;
			LOGPZ_DEBUG(logger,sout.str());
		}
#endif
        
		
		iv++;
	}
    //colocando a solucao para o connect interno usando shape descontinua
    
    TPZConnect *df2 = &this->Connect(ncon-1);
    int dfseq2 = df2->SequenceNumber();
    int dfvar2 = block.Size(dfseq2);
    int pos2 = block.Position(dfseq2);
    
    for (int idesc=0; idesc<data.numberdualfunctions; idesc++) {
        //	int indice=MeshSol.Rows()-data.numberdualfunctions +idesc;
        int iphi= data.phi.Rows()-data.numberdualfunctions +idesc;
        /*	#ifdef LOG4CXX
         {
         std::stringstream sout;
         sout << std::endl << "indice func descontinua = "<<iphi<<std::endl;
         sout << std::endl << "Indice para solucao da parte descontinua " <<indice<< std::endl;
         sout << std::endl << "coeficiente da parte descontinua " <<MeshSol(indice,0)<< std::endl;
         sout << std::endl << " valor da funcao descontinua " <<data.phi(iphi,0)<< std::endl;
         
         LOGPZ_DEBUG(logger,sout.str());
         }
         #endif*/
        data.sol[3]+= data.phi(iphi,0)*MeshSol(pos2+idesc,0);
		/*		#ifdef LOG4CXX
         {
         
         std::stringstream sout;
         sout << std::endl << "produto da func descon  " <<data.phi(iphi,0)*MeshSol(indice,0)<< std::endl;
         LOGPZ_DEBUG(logger,sout.str());
         }
         #endif*/
    }
    
    
	/*	#ifdef LOG4CXX
     {
     
     std::stringstream sout;
     sout << std::endl << "VecSol apos todas as iteracoes " <<data.sol<< std::endl;
     LOGPZ_DEBUG(logger,sout.str());
     }
     #endif
     
     */	
    
	
	
}


template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::Append(TPZFMatrix &u1, TPZFMatrix &u2, TPZFMatrix &u12)
{
	bool Is_u1PHI = (u1.Cols() == 1) ? true : false;
	bool Is_u2PHI = (u2.Cols() == 1) ? true : false;

	if(Is_u1PHI && Is_u2PHI)
	{
		int nu1 = u1.Rows(),nu2 = u2.Rows();
		u12.Redim(nu1+nu2,1);
		int i;
		for(i=0; i<nu1; i++) u12(i,0) = u1(i,0);
		for(i=0; i<nu2; i++) u12(i+nu1,0) = u2(i,0);


	}
	else if(!Is_u1PHI || !Is_u2PHI) // Se u1 e u2 não são Phi's, implica em serem dPhi's
	{
		int ru1 = u1.Rows(), cu1 = u1.Cols(), ru2 = u2.Rows(), cu2 = u2.Cols();
		int ru12 = ru1 < ru2 ? ru2 : ru1;
		int cu12 = cu1+cu2;
		u12.Redim(ru12,cu12);
		int i,j;
		for(i=0; i<ru1; i++) for(j=0; j<cu1; j++) u12(i,j) = u1(i,j);
		for(i=0; i<ru2; i++) for(j=0; j<cu2; j++) u12(i,j+cu1) = u2(i,j);//---modifiquei--
	}
	else
	{
		PZError << "TPZCompElHDiv::Append. Bad input parameters " << std::endl;//Este metodo so serve para u1 E u2 do mesmo tipo


	}

}

/** compute the shape functions corresponding to the dual space
 *
 */
template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::ShapeDual(TPZVec<REAL> &qsi, TPZFMatrix &phi, TPZFMatrix &dphi)
{
	  int dimension= TSHAPE::Dimension;
	REAL C=1;//fator de escala utilizado neste metodo
	TPZManVector<REAL,3> X0(3,0.);//centro do elemento

    int degree= this->fPressureOrder;
//	const int nshapedisc = pzshape::TPZShapeDisc::NShapeF(degree, dimension, pzshape::TPZShapeDisc:: ETensorial);
	pzshape::TPZShapeDisc::Shape(dimension,C,X0,qsi,degree,phi,dphi, pzshape::TPZShapeDisc:: ETensorial);

}

template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::Shape(TPZVec<REAL> &pt, TPZFMatrix &phi, TPZFMatrix &dphi) {
  TPZManVector<int,TSHAPE::NCornerNodes> id(TSHAPE::NCornerNodes,0);
  TPZManVector<int, TSHAPE::NSides-TSHAPE::NCornerNodes+1> ord(TSHAPE::NSides-TSHAPE::NCornerNodes,0); //nao vou acrescentar aq a ordem da variavel dual
//	TPZVec<int> ord(NConnects()-1);
  int i;
  TPZGeoEl *ref = this->Reference();
  for(i=0; i<TSHAPE::NCornerNodes; i++) {
    id[i] = ref->NodePtr(i)->Id();
  }


	for(i=0; i< NConnects()-1; i++) {//estou tirando o connect da variavel dual
		ord[i] = ConnectOrder(i);
		#ifdef LOG4CXX
		std::stringstream sout;
		sout << " ordem do connect " << i << " " << ord[i]<<endl;
		LOGPZ_DEBUG(logger,sout.str())
		#endif
  }

  int dimension= TSHAPE::Dimension;


	const int nshapecont = TSHAPE::NShapeF(ord);
	TPZFNMatrix<660> phiCont(nshapecont,1);
	TPZFNMatrix<660> dphiCont(dimension,nshapecont);
	TSHAPE::Shape(pt,id,ord,phiCont,dphiCont);


	// acrescentar as funcoes de pressao (via descontinuo)
	REAL C=1;//fator de escala utilizado neste metodo
	TPZManVector<REAL,3> X0(3,0.);//centro do elemento

    int degree= this->fPressureOrder;
	const int nshapedisc = pzshape::TPZShapeDisc::NShapeF(degree, dimension, pzshape::TPZShapeDisc:: ETensorial);
	TPZFNMatrix<660> phiDisc(nshapedisc,1);

	TPZFNMatrix<660> dphiDisc(dimension,nshapedisc);
	pzshape::TPZShapeDisc::Shape(dimension,C,X0,pt,degree,phiDisc,dphiDisc, pzshape::TPZShapeDisc:: ETensorial);
	this->Append(phiCont,phiDisc,phi);
	this->Append(dphiCont,dphiDisc,dphi);
	#ifdef LOG4CXX
	std::stringstream sout;
	sout << " Tamanho do vetor de Shape continuas " << phiCont.Rows()<<" Tamanho do vetor de Shape descontinuas " << phiDisc.Rows()<<" Tamanho do vetor de Shape Total " << phi.Rows()<<endl;
    sout << " Tamanho do vetor de DPHI continuas (linhas)" << dphiCont.Rows()<<" Tamanho do vetor de Shape descontinuas " << dphiDisc.Rows()<<" Tamanho do vetor de DPHI Total " << dphi.Rows()<<endl;
	sout << " Tamanho do vetor de DPHI continuas (colunas)" << dphiCont.Cols()<<" Tamanho do vetor de DPHI descontinuas " << dphiDisc.Cols()<<" Tamanho do vetor de DPHI Total " << dphi.Cols()<<endl;
	LOGPZ_DEBUG(logger,sout.str())
	#endif

}


// template<class TSHAPE>
// void TPZCompElHDiv<TSHAPE>::Solution(TPZVec<REAL> &qsi,int var,TPZManVector<REAL> &sol) {

// }

  /// Jorge 09/06/2001
   // Returns the transformation which transform a point from the side to the interior of the element
   //
template<class TSHAPE>
TPZTransform TPZCompElHDiv<TSHAPE>::TransformSideToElement(int side){
  return TSHAPE::TransformSideToElement(side);
}

template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::ComputeShapeIndex(TPZVec<int> &sides, TPZVec<int> &shapeindex){

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
		sout << endl<<"sides associated with the normals "<< sides <<
		"\nnormal associated with each shape function : shape function indexes " << shapeindex;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	/*
	int v=0; //Contador dos vetores
	int ind=0; //Contador de index ou funções de forma em Hdiv
	data.fVecIndex.Resize(0);
	data.fShapeIndex.Resize(0);

#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Elemento = "<< this->Reference()-> Id()<< std::endl;
		sout << "Loop sobre os lados do elemento"<< std::endl;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	int nSidelEl = this->Reference()->NSides();
	for(int isideEl=0; isideEl< nSidelEl; isideEl++)// loop sobre os lados do elemento
	{
#ifdef LOG4CXX
		{
			std::stringstream sout;
			sout << "	Lado "<< isideEl<< std::endl;
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif

		if(this->Reference()->SideDimension(isideEl)== this->Dimension()-1) //lados 1D -> Nós e Normais
		{
#ifdef	LOG4CXX
			{
				std::stringstream sout;
				sout << "		dim "<< this->Reference()->SideDimension(isideEl) << std::endl;
				sout << "		Loop sobre os nos do lado " << endl;
				LOGPZ_DEBUG(logger,sout.str())
			}
#endif
			//loop sob os nós do lado
			for(int jnode=0;jnode<this->Reference()->NSideNodes(isideEl);jnode++)
			{
#ifdef LOG4CXX
				{
					std::stringstream sout;
					sout << "		no: "<< jnode<< endl;
					LOGPZ_DEBUG(logger,sout.str())
				}
#endif
				int sideNod  = this->Reference()->SideNodeLocIndex(isideEl,jnode);
#ifdef LOG4CXX
				{
					std::stringstream sout;
					sout << endl;
					sout<< "			lado="<<sideNod<< endl;
					//sout<< "			NConnects----="<< NConnects()<<endl;
					sout<< "			numFuncoes="<< NConnectShapeF(sideNod)<<endl;
					sout<< "			firsIndex="<<FirstShapeIndex(sideNod)<< endl<< endl<<endl;
					LOGPZ_DEBUG(logger,sout.str())
				}
#endif
				int tamanho1=data.VecIndex.NElements()+NConnectShapeF(sideNod);
				data.VecIndex.Resize(tamanho1);
				data.ShapeIndex.Resize(tamanho1);
				int nada = NConnectShapeF(sideNod);
				for(int w=0; w < NConnectShapeF(sideNod);w++)
				{
					data.VecIndex[ind] = v ;
					data.ShapeIndex[ind] = FirstShapeIndex(sideNod)+w;
					ind++;
				}
				v++;
			}

#ifdef LOG4CXX
			{
				//				std::stringstream sout;
				sout << endl;
				sout<< "			lado="<< isideEl<< endl;
				sout<< "			numFuncoes="<< NConnectShapeF(isideEl)<<endl;
				sout<< "			firsIndex="<<FirstShapeIndex(isideEl)<< endl<< endl<<endl;
				//			LOGPZ_DEBUG(logger,sout.str())
			}
#endif
			//normal ao lado
			data.VecIndex.Resize(data.VecIndex.NElements()+NConnectShapeF(isideEl));
			data.ShapeIndex.Resize(data.ShapeIndex.NElements()+NConnectShapeF(isideEl));
			for(int w=0; w < NConnectShapeF(isideEl);w++)
			{
				data.VecIndex[ind] = v ;
				data.ShapeIndex[ind] = FirstShapeIndex(isideEl)+ w ;
				ind++;
			}
			v++;
		}
		if(this->Reference()->SideDimension(isideEl)== this->Dimension()) // 2D -> internas tangentes e do meio
		{
#ifdef LOG4CXX
			{
				//				std::stringstream sout;
				sout << "		dim "<< this->Reference()->SideDimension(isideEl) << std::endl;

				//internas tangentes
				sout << "		internas tangentes"<< endl;
				//			LOGPZ_DEBUG(logger,sout.str())
			}
#endif
			TPZStack<int> smallsides;
			this->Reference()-> LowerDimensionSides(isideEl,smallsides);
			for(int m=0;m<smallsides.NElements();m++)
			{
				if(this->Reference()->SideDimension(smallsides[m])==1)
				{
#ifdef LOG4CXX
					{
						//						std::stringstream sout;
						sout<< "				lado(internoTangente)="<<smallsides[m]<< endl;
						sout<< "				numFuncoes="<< NConnectShapeF(smallsides[m])<<endl;
						sout<< "				firsIndex="<<FirstShapeIndex(smallsides[m])<< endl<< endl<<endl;
						//					LOGPZ_DEBUG(logger,sout.str())
					}
#endif

					data.VecIndex.Resize(data.VecIndex.NElements()+NConnectShapeF(smallsides[m]));
					data.ShapeIndex.Resize(data.ShapeIndex.NElements()+NConnectShapeF(smallsides[m]));
					for(int w=0; w < NConnectShapeF(smallsides[m]);w++)
					{
						data.VecIndex[ind] = v ;
						data.ShapeIndex[ind] = FirstShapeIndex(smallsides[m])+ w ;
						ind++;
					}
					v++;
				}
			}
			//internas do meio
#ifdef LOG4CXX
			{
				//				std::stringstream sout;
				sout << "		internas meio"<< endl;
				//			LOGPZ_DEBUG(logger,sout.str())
			}
#endif

#ifdef LOG4CXX
			{
				//				std::stringstream sout;
				sout<< "				lado(internoMeio)="<<isideEl<< endl;
				sout<< "				numFuncoes="<< NConnectShapeF(isideEl)<<endl;
				sout<< "				firsIndex="<<FirstShapeIndex(isideEl)<< endl<< endl<<endl;
				//			LOGPZ_DEBUG(logger,sout.str())
			}
#endif
			data.VecIndex.Resize(data.VecIndex.NElements()+NConnectShapeF(isideEl));
			data.ShapeIndex.Resize(data.ShapeIndex.NElements()+NConnectShapeF(isideEl));
			cout<<"NConnectShapeF = "<<NConnectShapeF(isideEl)<<endl;
			for(int w=0; w < NConnectShapeF(isideEl);w++)
			{
				data.VecIndex[ind] = v ;
				data.ShapeIndex[ind] = FirstShapeIndex(isideEl)+ w ;
				ind++;
			}
			v++;

#ifdef LOG4CXX
			{
				//				std::stringstream sout;
				sout<< "				lado(internoMeio)="<<isideEl<< endl;
				sout<< "				numFuncoes="<< NConnectShapeF(isideEl)<<endl;
				sout<< "				firsIndex="<<FirstShapeIndex(isideEl)<< endl<< endl<<endl;
				//				LOGPZ_DEBUG(logger,sout.str())
			}
#endif
			data.VecIndex.Resize(data.VecIndex.NElements()+NConnectShapeF(isideEl));
			cout << "vendo o novo tamanho do vec index"<< data.VecIndex.NElements()+NConnectShapeF(isideEl)<<endl;
			data.ShapeIndex.Resize(data.ShapeIndex.NElements()+NConnectShapeF(isideEl));
			for(int w=0; w < NConnectShapeF(isideEl);w++)
			{
				data.VecIndex[ind] = v ;
				data.ShapeIndex[ind] = FirstShapeIndex(isideEl)+ w ;
				ind++;
			}
			v++;
		}//if
	}//for

#ifdef LOG4CXX
	{
		//std::stringstream sout;
		sout << "data.VecIndex  "<< "  data.ShapeIndex"<<endl;
		sout << data.VecIndex.NElements() <<endl;
		for(int z=0;z<data.VecIndex.NElements();z++)
		{
			sout<< data.VecIndex[z]<<"		"<<data.ShapeIndex[z]<<endl;
		}
		sout <<endl;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif

*/

}

/** Initialize a material data and its attributes based on element dimension, number
 * of state variables and material definitions
 */

template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::InitMaterialData(TPZMaterialData &data)
{
	TPZIntelGen<TSHAPE>::InitMaterialData(data);
#ifdef LOG4CXX
	{
		LOGPZ_DEBUG(logger,"Initializing normal vectors")
	}
#endif
	TPZManVector<int> normalsides;
	TPZIntelGen<TSHAPE>::Reference()->ComputeNormals(data.fNormalVec, normalsides);
	// vecindex : lado associado a cada normal
	// vecindex indica apenas o numero do lado associado a cada normal
	// agora temos que expandir para formar pares : vecIndex e shapeindex
//	ComputeShapeIndex(data.fVecIndex,data.fShapeIndex);
	data.numberdualfunctions = NConnectShapeF(NConnects()-1);
	IndexShapeToVec(normalsides,data.fVecShapeIndex);
#ifdef LOG4CXX
	{
		std::stringstream sout;
		data.fNormalVec.Print("Normal vectors", sout);
		sout << "NormalVector/Shape indexes " << data.fVecShapeIndex << std::endl;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
}

///  Save the element data to a stream

template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::Write(TPZStream &buf, int withclassid)
{
  TPZInterpolatedElement::Write(buf,withclassid);
  TPZManVector<int,3> order(3,0);
 this-> fIntRule.GetOrder(order);
this->  WriteObjects(buf,order);
  buf.Write(this->fConnectIndexes,TSHAPE::NSides);
  buf.Write(&this->fPreferredOrder,1);
  int classid = this->ClassId();
  buf.Write ( &classid, 1 );
}


  //Read the element data from a stream

template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::Read(TPZStream &buf, void *context)
{
  TPZInterpolatedElement::Read(buf,context);
  TPZManVector<int,3> order;
 this-> ReadObjects(buf,order);
 this-> fIntRule.SetOrder(order);
  buf.Read(this->fConnectIndexes,TSHAPE::NSides);
  buf.Read(&this->fPreferredOrder,1);
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


#include "tpzint1point.h"
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
void TPZCompElHDiv<TPZShapePoint>::CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) {
  if(dimension == 0) std::cout << "A point element has no graphical representation\n";
}



template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) {
  if(dimension == TSHAPE::Dimension && this->Material()->Id() > 0) {
    new typename TSHAPE::GraphElType(this,&grafgrid);
  }
}

template<>
int TPZCompElHDiv<TPZShapePoint>::ClassId() const
{
  return TPZHDIVPOINTID;
}

template class
		TPZRestoreClass< TPZCompElHDiv<TPZShapePoint>, TPZHDIVPOINTID>;


template<>
int TPZCompElHDiv<TPZShapeLinear>::ClassId() const
{
  return TPZHDIVLINEARID;
}

template class
		TPZRestoreClass< TPZCompElHDiv<TPZShapeLinear>, TPZINTELLINEARID>;

template<>
int TPZCompElHDiv<TPZShapeTriang>::ClassId() const
{
	return TPZHDIVTRIANGLEID;
}

template class
		TPZRestoreClass< TPZCompElHDiv<TPZShapeTriang>, TPZHDIVTRIANGLEID>;

template<>
int TPZCompElHDiv<TPZShapeQuad>::ClassId() const
{
	return TPZHDIVQUADID;
}

template class
		TPZRestoreClass< TPZCompElHDiv<TPZShapeQuad>, TPZHDIVQUADID>;

template<>
int TPZCompElHDiv<TPZShapeCube>::ClassId() const
{
	return TPZHDIVCUBEID;
}

template class
		TPZRestoreClass< TPZCompElHDiv<TPZShapeCube>, TPZHDIVCUBEID>;

template<>
int TPZCompElHDiv<TPZShapeTetra>::ClassId() const
{
	return TPZHDIVTETRAID;
}

template class
		TPZRestoreClass< TPZCompElHDiv<TPZShapeTetra>, TPZHDIVTETRAID>;

template<>
int TPZCompElHDiv<TPZShapePrism>::ClassId() const
{
	return TPZHDIVPRISMID;
}

template class
		TPZRestoreClass< TPZCompElHDiv<TPZShapePrism>, TPZHDIVPRISMID>;

template<>
int TPZCompElHDiv<TPZShapePiram>::ClassId() const
{
	return TPZHDIVPYRAMID;
}

template class
		TPZRestoreClass< TPZCompElHDiv<TPZShapePiram>, TPZHDIVPYRAMID>;


template class TPZCompElHDiv<TPZShapeTriang>;
template class TPZCompElHDiv<TPZShapePoint>;
template class TPZCompElHDiv<TPZShapeLinear>;
template class TPZCompElHDiv<TPZShapeQuad>;
template class TPZCompElHDiv<TPZShapeTetra>;
template class TPZCompElHDiv<TPZShapePrism>;
template class TPZCompElHDiv<TPZShapePiram>;
template class TPZCompElHDiv<TPZShapeCube>;


TPZCompEl * CreateHDivPointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	return new TPZCompElHDiv<TPZShapePoint>(mesh,gel,index);
}


TPZCompEl * CreateHDivLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	return new TPZCompElHDivBound2< TPZShapeLinear>(mesh,gel,index);
}

TPZCompEl * CreateHDivQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	return new TPZCompElHDiv< TPZShapeQuad>(mesh,gel,index);
}

TPZCompEl * CreateHDivTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	return new TPZCompElHDiv< TPZShapeTriang >(mesh,gel,index);
}

TPZCompEl * CreateHDivCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	return new TPZCompElHDiv< TPZShapeCube >(mesh,gel,index);
}

TPZCompEl * CreateHDivPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	return new TPZCompElHDiv< TPZShapePrism>(mesh,gel,index);
}

TPZCompEl * CreateHDivPyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	return new TPZCompElHDiv< TPZShapePiram >(mesh,gel,index);
}

TPZCompEl * CreateHDivTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	return new TPZCompElHDiv< TPZShapeTetra >(mesh,gel,index);
}



