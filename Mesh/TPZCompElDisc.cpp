// _*_ c++ _*_ 
#include "pzelgc3d.h"
#include "pzelgt3d.h"
#include "pzelgpi3d.h"
#include "pzelgpr3d.h"
#include "pzelgt2d.h"
#include "pzelgq2d.h"
//#include "pzelcq2d.h"
#include "pzmatrix.h"
#include "pzelmat.h"
#include "pzquad.h"
#include "pzgeoel.h"
#include "pzcmesh.h"
#include "pzerror.h"
#include "pzconnect.h"
#include "pzshapequad.h"
#include "pzmaterial.h"
#include "pzbndcond.h"
#include "pzmat1dlin.h"
#include "pztempmat.h"
#include "pzmanvector.h"
#include "TPZShapeDisc.h"
#include "TPZCompElDisc.h"
#include "TPZInterfaceEl.h"
//#include "TPZRefPattern.h"
#include "time.h"
#include "pzgeoel.h"
#include "pzcompel.h"
#include <math.h>
#include <stdio.h>

static int fGridDimension = 2;


TPZCompElDisc::TPZCompElDisc(TPZCompMesh &mesh,TPZGeoEl *ref,int &index) :
		TPZCompEl(mesh,index), fCenterPoint(3) {//dois pontos por eixo
  int i;
  fDegree = gOrder;
  fReference = ref;
  ref->SetReference(this);
  fMesh = &mesh;
  int materialid = ref->MaterialId();
  fMaterial = mesh.FindMaterial(materialid);
  fConnectIndex  = CreateMidSideConnect();
  mesh.ConnectVec()[fConnectIndex].IncrementElConnected();
  ref->CenterPoint(ref->NSides()-1,fCenterPoint);
  TPZVec<REAL> csi(fCenterPoint);
  ref->X(csi,fCenterPoint);
  fConstC = NormalizeConst();
  //criando os elementos interface
  CreateInterface();//argumento: número de faces
}

/*
//HEXAEDRO
TPZCompElDisc::TPZCompElDisc(TPZCompMesh &mesh,TPZGeoElC3d *ref,int &index) :
		TPZCompEl(mesh,index), fCenterPoint(3) {//dois pontos por eixo
  int i;
  fDegree = gOrder;
  fReference = ref;
  ref->SetReference(this);
  fMesh = &mesh;
  int materialid = ref->MaterialId();
  fMaterial = mesh.FindMaterial(materialid);
  fConnectIndex  = CreateMidSideConnect();
  mesh.ConnectVec()[fConnectIndex].IncrementElConnected();
  ref->CenterPoint(fCenterPoint,ref->NSides()-1);
  TPZVec<REAL> csi(fCenterPoint);
  ref->X(csi,fCenterPoint);
  fConstC = NormalizeConst();
  //criando os elementos interface
  CreateInterface();//argumento: número de faces
}

//TETRAEDRO
TPZCompElDisc::TPZCompElDisc(TPZCompMesh &mesh,TPZGeoElT3d *ref,int &index) :
		TPZCompEl(mesh,index), fCenterPoint(3) {//dois pontos por eixo
  int i;
  fDegree = gOrder;
  fReference = ref;
  ref->SetReference(this);
  fMesh = &mesh;
  int materialid = ref->MaterialId();
  fMaterial = mesh.FindMaterial(materialid);
  fConnectIndex  = CreateMidSideConnect();
  mesh.ConnectVec()[fConnectIndex].IncrementElConnected();
  ref->CenterPoint(fCenterPoint,ref->NSides()-1);
  TPZVec<REAL> csi(fCenterPoint);
  ref->X(csi,fCenterPoint);
  fConstC = NormalizeConst();
  //criando os elementos interface
  CreateInterface();//argumento: número de faces
}

//PIRÂMIDE
TPZCompElDisc::TPZCompElDisc(TPZCompMesh &mesh,TPZGeoElPi3d *ref,int &index) :
		TPZCompEl(mesh,index), fCenterPoint(3) {//dois pontos por eixo
  int i;
  fDegree = gOrder;
  fReference = ref;
  ref->SetReference(this);
  fMesh = &mesh;
  int materialid = ref->MaterialId();
  fMaterial = mesh.FindMaterial(materialid);
  fConnectIndex  = CreateMidSideConnect();
  mesh.ConnectVec()[fConnectIndex].IncrementElConnected();
  ref->CenterPoint(fCenterPoint,ref->NSides()-1);
  TPZVec<REAL> csi(fCenterPoint);
  ref->X(csi,fCenterPoint);
  fConstC = NormalizeConst();
  //criando os elementos interface
  CreateInterface();//argumento: número de faces
}


//PRISMA
TPZCompElDisc::TPZCompElDisc(TPZCompMesh &mesh,TPZGeoElPr3d *ref,int &index) :
		TPZCompEl(mesh,index), fCenterPoint(3) {//dois pontos por eixo
  int i;
  fDegree = gOrder;
  fReference = ref;
  ref->SetReference(this);
  fMesh = &mesh;
  int materialid = ref->MaterialId();
  fMaterial = mesh.FindMaterial(materialid);
  fConnectIndex  = CreateMidSideConnect();
  mesh.ConnectVec()[fConnectIndex].IncrementElConnected();
  ref->CenterPoint(fCenterPoint,ref->NSides()-1);
  TPZVec<REAL> csi(fCenterPoint);
  ref->X(csi,fCenterPoint);
  fConstC = NormalizeConst();
  //criando os elementos interface
  CreateInterface();//argumento: número de faces
}

//ELEMENTO DISCONTINUO TRIÂNGULO
TPZCompElDisc::TPZCompElDisc(TPZCompMesh &mesh,TPZGeoElT2d *ref,int &index) :
		TPZCompEl(mesh,index), fCenterPoint(3) {//dois pontos por eixo
  int i;
  fReference = ref;
  fDegree = gOrder;
  ref->SetReference(this);
  fMesh = &mesh;
  int materialid = ref->MaterialId();
  fMaterial = mesh.FindMaterial(materialid);
  fConnectIndex  = CreateMidSideConnect();
  mesh.ConnectVec()[fConnectIndex].IncrementElConnected();
  ref->CenterPoint(fCenterPoint,ref->NSides()-1);
  TPZVec<REAL> csi(fCenterPoint);
  ref->X(csi,fCenterPoint);
  fConstC = NormalizeConst();
  //criando os elementos interface
  CreateInterface();//para a face do elemento 2D
}

//ELEMENTO DISCONTINUO QUADRILÁTERO
TPZCompElDisc::TPZCompElDisc(TPZCompMesh &mesh,TPZGeoElQ2d *ref,int &index) :
		TPZCompEl(mesh,index), fCenterPoint(3) {//dois pontos por eixo
  int i;
  fReference = ref;
  fDegree = gOrder;
  ref->SetReference(this);
  fMesh = &mesh;
  int materialid = ref->MaterialId();
  fMaterial = mesh.FindMaterial(materialid);
  fConnectIndex  = CreateMidSideConnect();
  mesh.ConnectVec()[fConnectIndex].IncrementElConnected();
  ref->CenterPoint(fCenterPoint,ref->NSides()-1);
  TPZVec<REAL> csi(fCenterPoint);
  ref->X(csi,fCenterPoint);
  fConstC = NormalizeConst();
  //criando os elementos interface
  CreateInterface();//para a face do elemento 2D
}
*/

//INTERFACE
void TPZCompElDisc::CreateInterface(){

  TPZGeoElQ2d::SetCreateFunction(TPZInterfaceElement::CreateInterfaceQEl);//construtor do el. interface
  TPZGeoElT2d::SetCreateFunction(TPZInterfaceElement::CreateInterfaceTEl);//construtor do el. interface
  TPZStack<TPZCompElSide> list;
  int nsides = fReference->NSides();
  int vol = 0;
  //int dimension = Dimension();
  int side;
  nsides--;//last face
  for(side=nsides;side>=0;side--){
    if(fReference->SideDimension(side) != fGridDimension) continue;
    TPZCompElSide thisside(this,side);
    if(ExistsInterface(thisside)) {
      cout << "TPZCompElDisc::CreateInterface inconsistent\n";
      continue;
    }
    list.Resize(0);
    thisside.EqualLevelElementList(list,0,1);//retorna distinto ao atual ou nulo
    int size = list.NElements();
    if(size)//espera-se ter os elementos computacionais esquerdo e direito já criados antes de 
      fReference->CreateBCCompEl(side,fReference->MaterialId(),*Mesh());//criar o elemento interface
  }
  TPZGeoElQ2d::SetCreateFunction(CreateQ2Disc);//const. usual descontínuo
  TPZGeoElT2d::SetCreateFunction(CreateT2Disc);//const. usual descontínuo
}

REAL TPZCompElDisc::NormalizeConst(){
  //maior distancia entre o ponto interior e os vértices do elemento
  int nnodes = fReference->NNodes(),i;
  REAL maxdist,dist;
  TPZGeoNode node = fReference->Mesh()->NodeVec()[0];
  maxdist = pow(node.Coord(0)-fCenterPoint[0],2.)+pow(node.Coord(1)-fCenterPoint[1],2.);
  maxdist += pow(node.Coord(2)-fCenterPoint[2],2.);
  maxdist = sqrt(maxdist);
  for(i=1;i<nnodes;i++){
    node = fReference->Mesh()->NodeVec()[i];
    dist = pow(node.Coord(0)-fCenterPoint[0],2.)+pow(node.Coord(1)-fCenterPoint[1],2.);
    dist += pow(node.Coord(2)-fCenterPoint[2],2.);
    dist = sqrt(maxdist);
    if(maxdist < dist) maxdist = dist;
  }
  return maxdist;  
}

void TPZCompElDisc::Shape(TPZVec<REAL> X, TPZFMatrix &phi, TPZFMatrix &dphi) {

  if(Dimension()==1)
    TPZShapeDisc::Shape1D(fConstC,fCenterPoint,X,fDegree,phi,dphi);
  if(Dimension()==2)
    TPZShapeDisc::Shape2D(fConstC,fCenterPoint,X,fDegree,phi,dphi);
  if(Dimension()==3)
    TPZShapeDisc::Shape3D(fConstC,fCenterPoint,X,fDegree,phi,dphi);  
}

void TPZCompElDisc::Print(ostream &out) {

  out << "\nDiscontinous element : \n";
  out << "\tGeometric reference id : " << fReference->Id() << endl
      << "\tMaterial id : " << fReference->MaterialId() << endl
      << "\tDegrau of interpolation : " <<  fDegree << endl
      << "\tConnect index : " << fConnectIndex << endl
      << "\tNormalizing constant : " << fConstC << endl
      << "\tCenter point of the element : ";
  int size = fCenterPoint.NElements(),i;
  for(i=0;i<size-1;i++) out << fCenterPoint[i] << " , ";
  out << fCenterPoint[i] << endl;
}

/*
void TPZCompElDisc::NormalVector(int side,TPZVec<REAL> &int_point,
				    TPZVec<REAL> &normal, TPZFMatrix &axes, TPZFMatrix &norm_r3){
	Reference()->NormalVector(side,int_point,normal,axes,norm_r3);
}
*/

int TPZCompElDisc::ConnectIndex(int side) {

   return fConnectIndex;
}


int TPZCompElDisc::CreateMidSideConnect(){

  if(!Material())
    PZError << "\nTPZCompElDisc::CreateMidSideConnect Material nulo\n";
  int newnodeindex = Mesh()->AllocateNewConnect();
  return newnodeindex;
}

int TPZCompElDisc::NShapeF(){

  if(Dimension()==0) return 0;
  if(Dimension()==1) return (fDegree+1);
  if(Dimension()==2) return (fDegree+1)*(fDegree+2)/2;
  if(Dimension()==2){
    int i,sum=0;
    for(i=0;i<fDegree;i++) sum += (i+1)*(i+2)/2;;
    return sum;
  }
}


void TPZCompElDisc::InternalPoint(TPZVec<REAL> &point){

  point.Resize(3,0.);
  point[0] = fCenterPoint[0];
  point[1] = fCenterPoint[1];
  point[2] = fCenterPoint[2];
}


void TPZCompElDisc::CalcStiffDisc(TPZFMatrix &ek, TPZFMatrix &ef){

  if(fMaterial == NULL){
    cout << "TPZCompElDisc::CalcStiff : no material for this element\n";
    return;
  }

  int dim = Dimension();
  int nstate = fMaterial->NStateVariables();
  int nshape = NShapeF();
  int numeq = nshape * nstate;
  ek.Redim(numeq,numeq);
  ef.Redim(numeq,1);

  TPZFMatrix phix(nshape,1),dphix(dim,nshape);
  TPZFMatrix axes(3,3,0.);
    TPZFMatrix jacobian(dim,dim);
  TPZFMatrix jacinv(dim,dim);
  TPZVec<REAL> x(3,0.);
  TPZVec<REAL> intpoint(dim,0.);
  REAL detjac,weight,C = ConstC();
  TPZIntPoints *intrule = Reference()->CreateSideIntegrationRule(Reference()->NSides()-1,Degree());
  int npoints = intrule->NPoints(),ip;
  TPZVec<REAL> sol(nstate,0.);
  TPZFMatrix dsol(dim,nstate,0.);

  for(ip=0;ip<npoints;ip++){
    intrule->Point(ip,intpoint,weight);
    fReference->Jacobian( intpoint, jacobian, axes, detjac , jacinv);
    fReference->X(intpoint, x);
    weight *= fabs(detjac);
    Shape(x,phix,dphix);
    fMaterial->Contribute(x,jacinv,sol,dsol,weight,axes,phix,dphix,ek,ef);
  }

}

void TPZCompElDisc::Divide(int index,TPZVec<int> &subindex,int degree){

  int nsidesm2 = fReference->NSides()-2,face;
  TPZStack<TPZCompElSide> vec;
  int facedim =  fReference->Dimension()-1;
  for(face=nsidesm2;face>=0;face--){
    if(fReference->SideDimension(face) != facedim) break;//para no checar arestas e cantos
    TPZCompElSide thisside(this,face);
    thisside.EqualLevelElementList(vec,0,0);
    thisside.HigherLevelElementList(vec,0,0);    
  }
  int size = vec.NElements();
  for(face=0;face<size;face++){
    TPZCompEl *cel = vec[face].Element();
    if(cel->Reference()->Dimension() == facedim) delete cel;
  }
  fReference->ResetReference();
  //divide o elemento geométrico
  int nsubs = fReference->NSubElements();
  subindex.Resize(nsubs);
  TPZManVector<TPZGeoEl *> geosubs(nsubs);
  fReference->Divide(geosubs);
  if(!geosubs.NElements()) {
    subindex.Resize(0);
    return;
  }
  TPZCompElDisc *discel;
  int i,deg;
  if(degree) deg = degree;
  else deg = Degree();
  for (i=0;i<nsubs;i++)    {
    TPZGeoElC3d *ref = (TPZGeoElC3d *) geosubs[i];
    TPZCompElDisc son(*Mesh(),ref,subindex[i]);
    discel = (TPZCompElDisc *) fMesh->ElementVec()[subindex[i]];
    discel->SetDegree(deg);
  }
  Mesh()->ExpandSolution();//??????????????????????????????
  for(i=0; i<nsubs; i++) {
    discel = (TPZCompElDisc *) fMesh->ElementVec()[subindex[i]];
    discel->InterpolateSolution(*this);
  }
  //criando os elementos interface
  TPZCompEl *(*fpq)(TPZGeoElQ2d *geoel,TPZCompMesh &mesh,int &index);
  TPZCompEl *(*fpt)(TPZGeoElT2d *geoel,TPZCompMesh &mesh,int &index);
  fpq = TPZGeoElQ2d::fp;
  fpt = TPZGeoElT2d::fp;
  TPZGeoElQ2d::SetCreateFunction(TPZInterfaceElement::CreateInterfaceQEl);
  TPZGeoElT2d::SetCreateFunction(TPZInterfaceElement::CreateInterfaceTEl);
  for(i=0;i<nsubs;i++){
    TPZGeoEl *ref = geosubs[i];
    int matid = ref->MaterialId();
    int nsides = ref->NSides();
    int nnodes = ref->NNodes(),side;
    for(side=nnodes;side<nsides;side++)
      ref->CreateBCCompEl(side,matid,*Mesh());//partir em geométrico e computacional
  }
  //volta ao construtor usual
/*   TPZGeoElQ2d::SetCreateFunction(TPZGeoElQ2d::CreateEl); */
/*   TPZGeoElT2d::SetCreateFunction(TPZGeoElT2d::CreateEl); */
  TPZGeoElQ2d::SetCreateFunction(fpq);
  TPZGeoElT2d::SetCreateFunction(fpt);
}

void TPZCompElDisc::InterpolateSolution(TPZCompElDisc &coarsel){
  // accumulates the transfer coefficients between the current element and the
  // coarse element into the transfer matrix, using the transformation t
  TPZTransform t(Dimension());
  Reference()->BuildTransform2(fReference->NSides()-1,coarsel.Reference(),t);

  int locmatsize = NShapeF();
  int cormatsize = coarsel.NShapeF();
  int nvar = fMaterial->NStateVariables();
  int dimension = Dimension();

  TPZFMatrix loclocmat(locmatsize,locmatsize,0.);
  TPZFMatrix projectmat(locmatsize,nvar,0.);

  TPZVec<int> prevorder(1),order(1);
  TPZIntPoints *intrule = Reference()->CreateSideIntegrationRule(Reference()->NSides()-1,Degree());
  intrule->GetOrder(prevorder);
  order[0] = Degree()*2;
  intrule->SetOrder(order);

  TPZFMatrix locphi(locmatsize,1);
  TPZFMatrix locdphi(dimension,locmatsize);	// derivative of the shape function
  // in the master domain

  TPZFMatrix corphi(cormatsize,1);
  TPZFMatrix cordphi(dimension,cormatsize);	// derivative of the shape function
  // in the master domain

  TPZVec<REAL> int_point(dimension),coarse_int_point(dimension);
  TPZFMatrix jacobian(dimension,dimension),jacinv(dimension,dimension);
  TPZFMatrix axes(3,3,0.);
  REAL zero = 0.;
  TPZVec<REAL> x(3,zero);
  TPZVec<REAL> u(nvar);

  int numintpoints = intrule->NPoints();
  REAL weight;
  int lin,ljn,cjn;
  TPZConnect *df;
  TPZBlock &block = Mesh()->Block();

  for(int int_ind = 0; int_ind < numintpoints; ++int_ind) {

    intrule->Point(int_ind,int_point,weight);
    REAL jac_det = 1.;
    Reference()->Jacobian( int_point, jacobian , axes,jac_det,jacinv);
    Reference()->X(int_point, x);
    Shape(int_point,locphi,locdphi);
    weight *= jac_det;
    t.Apply(int_point,coarse_int_point);
    coarsel.Shape(coarse_int_point,corphi,cordphi);
    u.Fill(0.);
    int iv = 0;
    
    df = &coarsel.Connect(0);
    int dfseq = df->SequenceNumber();
    
    int dfvar = block.Size(dfseq);
    for(ljn=0; ljn<dfvar; ljn++) {
      u[iv%nvar] += corphi(iv/nvar,0)*block(dfseq,0,ljn,0);
      iv++;
    }

    for(lin=0; lin<locmatsize; lin++) {
      for(ljn=0; ljn<locmatsize; ljn++) {
	loclocmat(lin,ljn) += weight*locphi(lin,0)*locphi(ljn,0);
      }
      for(cjn=0; cjn<nvar; cjn++) {
	projectmat(lin,cjn) += weight*locphi(lin,0)*u[cjn];
      }
    }
    jacobian.Zero();
  }
  loclocmat.SolveDirect(projectmat,ELU);
  // identify the non-zero blocks for each row
  int iv=0;

  df = &Connect(0);
  int dfseq = df->SequenceNumber();
  int dfvar = block.Size(dfseq);
  for(ljn=0; ljn<dfvar; ljn++) {
    block(dfseq,0,ljn,0) = projectmat(iv/nvar,iv%nvar);
    iv++;
  }

  intrule->SetOrder(prevorder);
  
}

int TPZCompElDisc::ExistsInterface(TPZCompElSide compsd){

  TPZGeoElSide geosd = compsd.Reference();
  TPZGeoElSide  neighside = geosd.Neighbour();
  while(neighside.Element() && neighside.Element() != geosd.Element()){
    TPZCompElSide neighcompside = neighside.Reference();
    neighside = neighside.Neighbour();
    if(!neighcompside.Element()) continue;
    if(neighcompside.Element()->Type() == 17) 
      return 1;
  }
  return 0;
}

static int NFaces[7] = {0,1,1,4,5,6,8};//ponto,triângulo,quadrilatero,

void TPZCompElDisc::RemoveInterfaces(){

  int nsides = Reference()->NSides();
  int nfaces = NFaces[Type()],face;
  for(face=0;face<nfaces;face++){
    
  }

}

/* void TPZGeoEl::RemoveConectivities(){ */

/*   if(Reference() && Reference()->Type() != 17){ */
/*     PZError << "TPZGeoEl::RemoveConectivities warning: this element isn't type TPZInterfaceElement\n"; */
/*   } */
/*   int side,nsides = NSides()-1; */
/*   if(Dimension() == 3) nsides--;//tirando a connectividade associadas ao interior */
/*   for(side=0;side<nsides;side++){ */
/*     TPZGeoElSide thisside(this,side); */
/*     thisside.RemoveConnectivity(); */
/*   } */
/* } */

