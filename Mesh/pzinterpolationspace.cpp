//$Id: pzinterpolationspace.cpp,v 1.13 2007-06-08 00:02:28 cesar Exp $

#include "pzinterpolationspace.h"
#include "pzmaterialdata.h"
#include "pzbndcond.h"
#include "pzelmat.h"
#include "pzquad.h"
#include "TPZCompElDisc.h"
#include "TPZInterfaceEl.h"
#include "pztransfer.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.TPZInterpolationSpace"));
#endif

TPZInterpolationSpace::TPZInterpolationSpace()
 : TPZCompEl()
{
  fPreferredOrder = -1;
}

TPZInterpolationSpace::TPZInterpolationSpace(TPZCompMesh &mesh, const TPZInterpolationSpace &copy)
 : TPZCompEl(mesh, copy)
{
  fPreferredOrder = copy.fPreferredOrder;
}

TPZInterpolationSpace::TPZInterpolationSpace(TPZCompMesh &mesh, const TPZInterpolationSpace &copy, std::map<int,int> &gl2lcElMap)
  : TPZCompEl(mesh, copy, gl2lcElMap)
{
  fPreferredOrder = copy.fPreferredOrder;
}

TPZInterpolationSpace::TPZInterpolationSpace(TPZCompMesh &mesh, const TPZInterpolationSpace &copy, int &index)
  : TPZCompEl(mesh, copy, index)
{
  fPreferredOrder = copy.fPreferredOrder;
}

TPZInterpolationSpace::TPZInterpolationSpace(TPZCompMesh &mesh, TPZGeoEl *gel, int &index)
 : TPZCompEl(mesh,gel,index)
{
  fPreferredOrder = mesh.GetDefaultOrder();
}

TPZInterpolationSpace::~TPZInterpolationSpace(){}

int TPZInterpolationSpace::MaxOrder(){
  const int n = this->NConnects();
  int result = 0;
  int side;
  for(int i = 0; i < n; i++){
    side = this->Connect(i).Order();
    if (side > result) result = side;
  }//i
  return result;
}

void TPZInterpolationSpace::ComputeShape(TPZVec<REAL> &intpoint, TPZVec<REAL> &X,
                                         TPZFMatrix &jacobian, TPZFMatrix &axes,
                                         REAL &detjac, TPZFMatrix &jacinv,
                                         TPZFMatrix &phi, TPZFMatrix &dphix){
  TPZGeoEl * ref = this->Reference();
  if (!ref){
    PZError << "\nERROR AT " << __PRETTY_FUNCTION__ << " - this->Reference() == NULL\n";
    return;
  }//if
  TPZFNMatrix<660> dphi(dphix.Rows(), dphix.Cols(), 0.);
  int dim = this->Dimension();
  int nshape = NShapeF();

  ref->Jacobian( intpoint, jacobian, axes, detjac , jacinv);
  this->Shape(intpoint,phi,dphi);
  int ieq;
  switch(dim){
  case 0:
    break;
  case 1:
    dphix = dphi;
    dphix *= (1./detjac);
    break;
  case 2:
    for(ieq = 0; ieq < nshape; ieq++) {
      dphix(0,ieq) = jacinv(0,0)*dphi(0,ieq) + jacinv(1,0)*dphi(1,ieq);
      dphix(1,ieq) = jacinv(0,1)*dphi(0,ieq) + jacinv(1,1)*dphi(1,ieq);
    }
    break;
  case 3:
    for(ieq = 0; ieq < nshape; ieq++) {
      dphix(0,ieq) = jacinv(0,0)*dphi(0,ieq) + jacinv(1,0)*dphi(1,ieq) + jacinv(2,0)*dphi(2,ieq);
      dphix(1,ieq) = jacinv(0,1)*dphi(0,ieq) + jacinv(1,1)*dphi(1,ieq) + jacinv(2,1)*dphi(2,ieq);
      dphix(2,ieq) = jacinv(0,2)*dphi(0,ieq) + jacinv(1,2)*dphi(1,ieq) + jacinv(2,2)*dphi(2,ieq);
    }
    break;
  default:
    PZError << "Error at " << __PRETTY_FUNCTION__ << " please implement the " << dim << "d Jacobian and inverse\n";
  } //switch
  ref->X(intpoint, X);
}

REAL TPZInterpolationSpace::InnerRadius(){
  if (!this->Reference()){
    PZError << "Error at " << __PRETTY_FUNCTION__ << " this->Reference() == NULL\n";
    return 0.;
  }
  return this->Reference()->ElementRadius();
}

void TPZInterpolationSpace::InitMaterialData(TPZMaterialData &data){
  this->Material()->FillDataRequirements(data);
  const int dim = this->Dimension();
  const int nshape = this->NShapeF();
  const int nstate = this->Material()->NStateVariables();
  data.phi.Redim(nshape,1);
  data.dphix.Redim(dim,nshape);
  data.axes.Redim(3,3);
  data.jacobian.Redim(dim,dim);
  data.jacinv.Redim(dim,dim);
  data.x.Resize(3);
  if (data.fNeedsSol){
    data.sol.Resize(nstate);
    data.dsol.Redim(dim,nstate);
  }
}//void

void TPZInterpolationSpace::ComputeRequiredData(TPZMaterialData &data,
                                                TPZVec<REAL> &qsi){
  if (data.fNeedsNeighborSol){
    this->ComputeSolution(qsi, data.normal, data.soll, data.dsoll, data.axesleft, data.solr, data.dsolr, data.axesright);
  }//fNeedsNeighborSol

  if (data.fNeedsSol){
    if (data.phi.Rows()){//if shape functions are available
      this->ComputeSolution(qsi, data.phi, data.dphix, data.axes, data.sol, data.dsol);
    }
    else{//if shape functions are not available
      this->ComputeSolution(qsi, data.sol, data.dsol, data.axes);
    }
  }//fNeedsSol

  if (data.fNeedsHSize){
    data.HSize = 2.*this->InnerRadius();
  }//fNeedHSize
}//void

void TPZInterpolationSpace::CalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef){

  TPZAutoPointer<TPZMaterial> material = Material();
  if(!material){
    PZError << "Error at " << __PRETTY_FUNCTION__ << " this->Material() == NULL\n";
    ek.Reset();
    ef.Reset();
    return;
  }

  this->InitializeElementMatrix(ek,ef);

  TPZMaterialData data;
  this->InitMaterialData(data);
  data.p = this->MaxOrder();

  int dim = Dimension();
  TPZManVector<REAL,3> intpoint(dim,0.);
  REAL weight = 0.;

  TPZIntPoints &intrule = GetIntegrationRule();
  if(material->HasForcingFunction()) {
    TPZManVector<int,3> order(dim,intrule.GetMaxOrder());
    intrule.SetOrder(order);
  }

  int intrulepoints = intrule.NPoints();
  for(int int_ind = 0; int_ind < intrulepoints; ++int_ind){
    intrule.Point(int_ind,intpoint,weight);
    this->ComputeShape(intpoint, data.x, data.jacobian, data.axes, data.detjac, data.jacinv, data.phi, data.dphix);
    weight *= fabs(data.detjac);
    this->ComputeRequiredData(data, intpoint);
    material->Contribute(data,weight,ek.fMat,ef.fMat);
  }//loop over integratin points

}//CalcStiff

void TPZInterpolationSpace::CalcResidual(TPZElementMatrix &ef){

  TPZAutoPointer<TPZMaterial> material = Material();
  if(!material){
    PZError << "Error at " << __PRETTY_FUNCTION__ << " this->Material() == NULL\n";
    ef.Reset();
    return;
  }

  this->InitializeElementMatrix(ef);

  TPZMaterialData data;
  this->InitMaterialData(data);
  data.p = this->MaxOrder();

  int dim = Dimension();
  TPZManVector<REAL,3> intpoint(dim,0.);
  REAL weight = 0.;

  TPZIntPoints &intrule = GetIntegrationRule();
  if(material->HasForcingFunction()) {
    TPZManVector<int,3> order(dim,intrule.GetMaxOrder());
    intrule.SetOrder(order);
  }

  int intrulepoints = intrule.NPoints();
  for(int int_ind = 0; int_ind < intrulepoints; ++int_ind){
    intrule.Point(int_ind,intpoint,weight);
    this->ComputeShape(intpoint, data.x, data.jacobian, data.axes, data.detjac, data.jacinv, data.phi, data.dphix);
    weight *= fabs(data.detjac);
    this->ComputeRequiredData(data, intpoint);
    material->Contribute(data,weight,ef.fMat);
  }//loop over integratin points

}//CalcResidual

void TPZInterpolationSpace::InitializeElementMatrix(TPZElementMatrix &ek, TPZElementMatrix &ef){
  const int numdof = this->Material()->NStateVariables();
  const int ncon = this->NConnects();
  const int nshape = this->NShapeF();
  const int numeq = nshape*numdof;
  ek.fMat.Redim(numeq,numeq);
  ef.fMat.Redim(numeq,1);
  ek.fBlock.SetNBlocks(ncon);
  ef.fBlock.SetNBlocks(ncon);
  ek.fNumStateVars = numdof;
  ef.fNumStateVars = numdof;
  TPZManVector<REAL> sol(numdof,0.);
  int i;
  for (i = 0; i < ncon ; i++){
    ek.fBlock.Set(i,NConnectShapeF(i)*numdof);
    ef.fBlock.Set(i,NConnectShapeF(i)*numdof);
  }
  ek.fConnect.Resize(ncon);
  ef.fConnect.Resize(ncon);
  for(i=0; i<ncon; i++){
    (ef.fConnect)[i] = ConnectIndex(i);
    (ek.fConnect)[i] = ConnectIndex(i);
  }
}//void

void TPZInterpolationSpace::InitializeElementMatrix(TPZElementMatrix &ef){
  const int numdof = this->Material()->NStateVariables();
  const int ncon = this->NConnects();
  const int nshape = this->NShapeF();
  const int numeq = nshape*numdof;
  ef.fMat.Redim(numeq,1);
  ef.fBlock.SetNBlocks(ncon);
  ef.fNumStateVars = numdof;
  TPZManVector<REAL> sol(numdof,0.);
  int i;
  for (i = 0; i < ncon ; i++){
    ef.fBlock.Set(i,NConnectShapeF(i)*numdof);
  }
  ef.fConnect.Resize(ncon);
  for(i=0; i<ncon; i++){
    (ef.fConnect)[i] = ConnectIndex(i);
  }
}//void

void TPZInterpolationSpace::Solution(TPZVec<REAL> &qsi,int var,TPZVec<REAL> &sol) {

  if(var >= 100) {
    TPZCompEl::Solution(qsi,var,sol);
    return;
  }
  int dim = this->Dimension();
  if(var == 99) {
    sol[0] = this->MaxOrder();
    return;
  }

  TPZAutoPointer<TPZMaterial> material = this->Material();
  if(!material){
    sol.Resize(0);
    return;
  }

  int numdof = material->NStateVariables();
  TPZManVector<REAL,10> u(numdof);
  TPZFNMatrix<30> du(dim,numdof,0.);
  TPZFNMatrix<9> axes(3,3,0.);
  this->ComputeSolution(qsi, u, du, axes);
  material->Solution(u,du,axes,var,sol);
}

void TPZInterpolationSpace::InterpolateSolution(TPZInterpolationSpace &coarsel){
  // accumulates the transfer coefficients between the current element and the
  // coarse element into the transfer matrix, using the transformation t
  TPZAutoPointer<TPZMaterial> material = Material();
  if(!material){
    PZError << __PRETTY_FUNCTION__ << " this->Material() == NULL " << std::endl;
    return;
  }

  TPZTransform t(Dimension());
  TPZGeoEl *ref = Reference();

  //Cedric 16/03/99
  //  Reference()->BuildTransform(NConnects(),coarsel.Reference(),t);
  t = Reference()->BuildTransform2(ref->NSides()-1,coarsel.Reference(),t);

  int locnod = NConnects();
//   int cornod = coarsel.NConnects();
  int locmatsize = NShapeF();
  int cormatsize = coarsel.NShapeF();
  int nvar = material->NStateVariables();
  int dimension = Dimension();
  if (!dimension) {
    PZError << "\nExiting " << __PRETTY_FUNCTION__ << " - trying to interpolate a node solution.\n";
    return ;
  }

  TPZFMatrix loclocmat(locmatsize,locmatsize,0.);
  TPZFMatrix projectmat(locmatsize,nvar,0.);

  TPZManVector<int,3> prevorder(dimension);
  TPZIntPoints &intrule = GetIntegrationRule();
  intrule.GetOrder(prevorder);

  int thismaxorder = this->MaxOrder();
  int coarsemaxorder = coarsel.MaxOrder();
  int maxorder = (thismaxorder > coarsemaxorder) ? thismaxorder : coarsemaxorder;
  // Cesar 2003-11-25 -->> To avoid integration warnings...
  maxorder = (2*maxorder > intrule.GetMaxOrder() ) ? intrule.GetMaxOrder() : 2*maxorder;
  TPZManVector<int,3> order(dimension,maxorder);
  for(int dim = 0; dim < dimension; dim++){
    order[dim] = maxorder*2;
  }
  intrule.SetOrder(order);

  TPZFNMatrix<220> locphi(locmatsize,1);
  TPZFNMatrix<660> locdphi(dimension,locmatsize);//derivative of the shape function in the master domain

  TPZFNMatrix<220> corphi(cormatsize,1);
  TPZFNMatrix<660> cordphi(dimension,cormatsize);//derivative of the shape function in the master domain

  TPZManVector<REAL,3> int_point(dimension),coarse_int_point(dimension);
  TPZFNMatrix<9> jacobian(dimension,dimension),jacinv(dimension,dimension);
  TPZFNMatrix<9> axes(3,3,0.), coarseaxes(3,3,0.);
  REAL zero = 0.;
  TPZManVector<REAL,3> x(3,zero);
  TPZManVector<REAL,10> u(nvar);
  TPZFNMatrix<30> du(dimension,nvar);

  int numintpoints = intrule.NPoints();
  REAL weight;
  int lin,ljn,cjn;
  TPZConnect *df;
//   TPZBlock &coarseblock = coarsel.Mesh()->Block();

  for(int int_ind = 0; int_ind < numintpoints; ++int_ind) {
    intrule.Point(int_ind,int_point,weight);
    REAL jac_det = 1.;
    this->ComputeShape(int_point, x, jacobian, axes, jac_det, jacinv, locphi, locdphi);
    weight *= jac_det;
    t.Apply(int_point,coarse_int_point);
    coarsel.ComputeSolution(coarse_int_point, u, du, coarseaxes);

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
  TPZBlock &fineblock = Mesh()->Block();
  int iv=0,in;
  for(in=0; in<locnod; in++) {
    df = &Connect(in);
    int dfseq = df->SequenceNumber();
    int dfvar = fineblock.Size(dfseq);
    for(ljn=0; ljn<dfvar; ljn++) {
      fineblock(dfseq,0,ljn,0) = projectmat(iv/nvar,iv%nvar);
      iv++;
    }
  }
  intrule.SetOrder(prevorder);
}//InterpolateSolution

void TPZInterpolationSpace::CreateInterfaces(bool BetweenContinuous){
  //nao verifica-se caso o elemento de contorno
  //eh maior em tamanho que o interface associado
  //caso AdjustBoundaryElement nao for aplicado
  //a malha eh criada consistentemente
  TPZGeoEl *ref = Reference();
  int nsides = ref->NSides();
  int InterfaceDimension = this->Material()->Dimension() - 1;
  int side;
  nsides--;//last face
  for(side=nsides;side>=0;side--){
    if(ref->SideDimension(side) != InterfaceDimension) continue;
    TPZCompElSide thisside(this,side);
    if(this->ExistsInterface(thisside.Reference())) {
//      cout << "TPZCompElDisc::CreateInterface inconsistent: interface already exists\n";
      continue;
    }
    TPZStack<TPZCompElSide> highlist;
    thisside.HigherLevelElementList(highlist,0,1);
    //a interface se cria uma vez so quando existem ambos
    //elementos esquerdo e direito (computacionais)
    if(!highlist.NElements()) {
      this->CreateInterface(side, BetweenContinuous);//s�tem iguais ou grande => pode criar a interface
    } else {
      int ns = highlist.NElements();
      int is;
      for(is=0; is<ns; is++) {//existem pequenos ligados ao lado atual
        const int higheldim = highlist[is].Reference().Dimension();
	if(higheldim != InterfaceDimension) continue;
// 	TPZCompElDisc *del = dynamic_cast<TPZCompElDisc *> (highlist[is].Element());
// 	if(!del) continue;

	TPZCompEl *del = highlist[is].Element();
	if(!del) continue;

        TPZCompElSide delside( del, highlist[is].Side() );
        TPZInterpolationSpace * delSp = dynamic_cast<TPZInterpolationSpace*>(del);
        if (!delSp){
          PZError << "\nERROR AT " << __PRETTY_FUNCTION__ <<  " - CASE NOT AVAILABLE\n";
          return;
        }
        if ( delSp->ExistsInterface(delside.Reference()) ) {
//          cout << "TPZCompElDisc::CreateInterface inconsistent: interface already exists\n";
        }
        else{
	  delSp->CreateInterface(highlist[is].Side(), BetweenContinuous);
        }
      }
    }
  }
}

TPZInterfaceElement * TPZInterpolationSpace::CreateInterface(int side, bool BetweenContinuous)
{
//  LOGPZ_INFO(logger, "Entering CreateInterface");
  TPZInterfaceElement * newcreatedinterface = NULL;

  TPZGeoEl *ref = Reference();
  if(!ref) {
    LOGPZ_WARN(logger, "Exiting CreateInterface Null reference reached - NULL interface returned");
    return newcreatedinterface;
  }

  TPZCompElSide thisside(this,side);
  TPZStack<TPZCompElSide> list;
  list.Resize(0);
  thisside.EqualLevelElementList(list,0,1);//retorna distinto ao atual ou nulo
  int size = list.NElements();
  //espera-se ter os elementos computacionais esquerdo e direito
  //ja criados antes de criar o elemento interface
  if(size){
    //Interface has the same material of the neighbour with lesser dimension.
    //It makes the interface have the same material of boundary conditions (TPZCompElDisc with interface dimension)
    int matid;
    int thisdim = this->Dimension();
    int neighbourdim = list[0].Element()->Dimension();
    if (thisdim == neighbourdim){
//      matid = this->Material()->Id();
        matid = this->Mesh()->Reference()->InterfaceMaterial(this->Material()->Id(), list[0].Element()->Material()->Id() );
    }
    else { //one element is a boundary condition
      if (thisdim < neighbourdim) matid = this->Material()->Id();
      else matid = list[0].Element()->Material()->Id();
    }


    int index;

    TPZCompEl *list0 = list[0].Element();
    int list0side = list[0].Side();
    TPZCompElDisc * thisdisc  = dynamic_cast<TPZCompElDisc*>(this);
    TPZCompElDisc * neighdisc = dynamic_cast<TPZCompElDisc*>(list0);
    int thisside = side;
    int neighside = list0side;

    if (BetweenContinuous == false){
      //It means at least one element must be discontinuous
      if (!thisdisc && !neighdisc){
        return NULL;
      }
    }

    TPZGeoEl *gel = ref->CreateBCGeoEl(side,matid); //isto acertou as vizinhanas da interface geometrica com o atual


    if(Dimension() > list0->Dimension()){
       //o de volume eh o direito caso um deles seja BC
       //a normal aponta para fora do contorno
       TPZCompElSide thiscompelside(this, thisside);
       TPZCompElSide neighcompelside(list0, neighside);
       newcreatedinterface = new TPZInterfaceElement(*fMesh,gel,index,thiscompelside,neighcompelside);
    } else {
       //caso contrario ou caso ambos sejam de volume
       TPZCompElSide thiscompelside(this, thisside);
       TPZCompElSide neighcompelside(list0, neighside);
       newcreatedinterface = new TPZInterfaceElement(*fMesh,gel,index,neighcompelside,thiscompelside);
    }
    return newcreatedinterface;
  }

//If there is no equal level element, we try the lower elements.
//Higher elements will not be considered by this method. In that case the interface must be created by the neighbour.
  TPZCompElSide lower = thisside.LowerLevelElementList(0);
  if(lower.Exists()){
    //Interface has the same material of the neighbour with lesser dimension.
    //It makes the interface has the same material of boundary conditions (TPZCompElDisc with interface dimension)
    int matid;
    int thisdim = this->Dimension();
    int neighbourdim = lower.Element()->Dimension();

    if (thisdim == neighbourdim){
//      matid = this->Material()->Id();
        matid = this->Mesh()->Reference()->InterfaceMaterial(this->Material()->Id(), lower.Element()->Material()->Id() );
    }
    else { //one element is a boundary condition
      if (thisdim < neighbourdim) matid = this->Material()->Id();
      else matid = lower.Element()->Material()->Id();
    }


    TPZCompEl *lowcel = lower.Element();
    int lowside = lower.Side();
    TPZCompElDisc * thisdisc  = dynamic_cast<TPZCompElDisc*>(this);
    TPZCompElDisc * neighdisc = dynamic_cast<TPZCompElDisc*>(lowcel);
    int thisside = side;
    int neighside = lowside;

    if (BetweenContinuous == false){
      //It means at least one element must be discontinuous
      if (!thisdisc && !neighdisc){
        return NULL;
      }
    }

    //existem esquerdo e direito: this e lower
    TPZGeoEl *gel = ref->CreateBCGeoEl(side,matid);
    int index;


    if(Dimension() > lowcel->Dimension()){
       //para que o elemento esquerdo seja de volume
       TPZCompElSide thiscompelside(this, thisside);
       TPZCompElSide lowcelcompelside(lowcel, neighside);
       newcreatedinterface = new TPZInterfaceElement(*fMesh,gel,index,thiscompelside,lowcelcompelside);
    } else {
       TPZCompElSide thiscompelside(this, thisside);
       TPZCompElSide lowcelcompelside(lowcel, neighside);
       newcreatedinterface = new TPZInterfaceElement(*fMesh,gel,index,lowcelcompelside,thiscompelside);
    }
    return newcreatedinterface;
  }
  return newcreatedinterface;
}

int TPZInterpolationSpace::ExistsInterface(TPZGeoElSide geosd){

  TPZGeoElSide  neighside = geosd.Neighbour();
  while(neighside.Element() && neighside.Element() != geosd.Element()){
    TPZCompElSide neighcompside = neighside.Reference();
    neighside = neighside.Neighbour();
    if(!neighcompside.Element()) continue;
    if(neighcompside.Element()->Type() == EInterface)
      return 1;
  }
  return 0;
}

void TPZInterpolationSpace::RemoveInterfaces(){

  int nsides = Reference()->NSides();
  if (!this->Material()){
    std::stringstream mess;
    mess << __PRETTY_FUNCTION__ << " - this->Material() == NULL, I can't RemoveInterfaces()";
    PZError << mess.str() << std::endl;
    LOGPZ_ERROR(logger, mess.str());
    return;
  }
  int InterfaceDimension = this->Material()->Dimension() - 1;
  int is;
  TPZStack<TPZCompElSide> list,equal;
  for(is=0;is<nsides;is++){
    TPZCompElSide thisside(this,is);
    if(thisside.Reference().Dimension() != InterfaceDimension) continue;
    // procurar na lista de elementos iguais
    list.Resize(0);// o lado atual �uma face
    //thisside.EqualLevelElementList(list,0,0);// monta a lista de elementos iguais
    RemoveInterface(is);// chame remove interface do elemento atual (para o side atual)
    thisside.HigherLevelElementList(list,0,0);// procurar na lista de elementos menores (todos)
    int size = list.NElements(),i;            // 'isto pode incluir elementos interfaces'
    //tirando os elementos de interface da lista
    for(i=0;i<size;i++){
      if(list[i].Element()->Type() == EInterface) {
        LOGPZ_DEBUG(logger, "Removing interface element from the list of higher level elements");
        //This need to be done because otherwise list could be invalidated when an interface is removed.
        list[i] = TPZCompElSide();//tirando interface
      }
    }
    for(i=0;i<size;i++){// percorre os elementos menores
      if(!list[i].Element()) continue;
      TPZGeoElSide geolist = list[i].Reference();//TESTE
      if(geolist.Dimension() != InterfaceDimension) continue;
      equal.Resize(0);// para cada elemento menor e' preciso verificar a dimensao,
      list[i].EqualLevelElementList(equal,0,0);//montar a lista de elementos iguais (todos)
      equal.Push(list[i]);//n� �incorporado no m�odo anterior
      int neq = equal.NElements(),k=-1;
      while(++k < neq) if(equal[k].Element()->Type() != EInterface) break;//procurando elemento descont�uo cujo
      if(!neq || k == neq){                               //lado faz parte da parti� do lado side do this
	      LOGPZ_FATAL(logger, " Inconsistency of data");
	      exit(-1);//elemento descont�uo n� achado: ERRO
      }// chame removeinterface do elemento menor

      TPZInterpolationSpace * equalkSp = dynamic_cast<TPZInterpolationSpace*>(equal[k].Element());
      if (!equalkSp){
        PZError << "\nERROR AT " << __PRETTY_FUNCTION__ <<  " - CASE NOT AVAILABLE\n";
        return;
      }
      equalkSp->RemoveInterface(equal[k].Side());
    }
  }

}

void TPZInterpolationSpace::RemoveInterface(int side) {

  TPZStack<TPZCompElSide> list;
  list.Resize(0);
  TPZCompElSide thisside(this,side);
  thisside.EqualLevelElementList(list,0,0);// monta a lista de elementos iguais
  int size = list.NElements(),i=-1;
  while(++i < size) if(list[i].Element()->Type() == EInterface) break;// procura aquele que e derivado de TPZInterfaceEl
  if(!size || i == size){
    return;// nada a ser feito
  }
  // aqui existe a interface
  TPZCompEl *cel = list[i].Element();
  TPZGeoEl *gel = cel->Reference();
  gel->RemoveConnectivities();// deleta o elemento das vizinhancas
  TPZGeoMesh *gmesh = Mesh()->Reference();
  int index = gmesh->ElementIndex(gel);// identifica o index do elemento
  gmesh->ElementVec()[index] = NULL;
  delete cel;
  delete gel;// deleta o elemento
  gmesh->ElementVec().SetFree(index);// Chame SetFree do vetor de elementos da malha geometrica para o index

}

void TPZInterpolationSpace::EvaluateError(  void (*fp)(TPZVec<REAL> &loc,TPZVec<REAL> &val,TPZFMatrix &deriv),
                                            TPZVec<REAL> &errors,TPZBlock * /*flux */){
  int NErrors = this->Material()->NEvalErrors();
  errors.Resize(NErrors);
  errors.Fill(0.);
  TPZAutoPointer<TPZMaterial> material = Material();
  if(!material){
    PZError << "TPZInterpolatedElement::EvaluateError : no material for this element\n";
    Print(PZError);
    return;
  }
  if(dynamic_cast<TPZBndCond *>(material.operator ->())) {
    LOGPZ_INFO(logger,"Exiting EvaluateError - null error - boundary condition material.");
    return;
  }
  int problemdimension = Mesh()->Dimension();
  if(Reference()->Dimension() < problemdimension) return;

  // Adjust the order of the integration rule
  int dim = Dimension();
  TPZIntPoints &intrule = this->GetIntegrationRule();

  TPZManVector<int,3> prevorder(dim), maxorder(dim, this->MaxOrder());
  intrule.GetOrder(prevorder);
  intrule.SetOrder(maxorder);

  int ndof = material->NStateVariables();
  int nflux = material->NFluxes();
  TPZManVector<REAL,10> u_exact(ndof);
  TPZFNMatrix<90> du_exact(dim,ndof);
  TPZManVector<REAL,3> intpoint(3), values(NErrors);
  values.Fill(0.0);
  REAL weight;
  TPZManVector<REAL,9> flux_el(nflux,0.);

  TPZMaterialData data;
  this->InitMaterialData(data);

  for(int nint = 0; nint < intrule.NPoints(); nint++) {

    intrule.Point(nint,intpoint,weight);
    this->ComputeShape(intpoint, data.x, data.jacobian, data.axes, data.detjac, data.jacinv, data.phi, data.dphix);
    weight *= fabs(data.detjac);
    this->ComputeSolution(intpoint, data.phi, data.dphix, data.axes, data.sol, data.dsol);

    //contribuicoes dos erros
    if(fp) {
      fp(data.x,u_exact,du_exact);
      material->Errors(data.x,data.sol,data.dsol,data.axes,flux_el,u_exact,du_exact,values);
      for(int ier = 0; ier < NErrors; ier++)
        errors[ier] += values[ier]*weight;
    }
  }//fim for : integration rule
   //Norma sobre o elemento
  for(int ier = 0; ier < NErrors; ier++){
    errors[ier] = sqrt(errors[ier]);
  }//for ier

  intrule.SetOrder(prevorder);

}//method

void TPZInterpolationSpace::ComputeError(int errorid,
                                         TPZVec<REAL> &error){

  TPZAutoPointer<TPZMaterial> material = Material();
  if(!material){
    cout << "TPZCompElDisc::ComputeError : no material for this element\n";
    return;
  }

  TPZMaterialData data;
  this->InitMaterialData(data);

//   data.fPrimalExactSol = fp;
//   data.fDualExactSol = fd;

  REAL weight;
  int dim = Dimension();
  TPZVec<REAL> intpoint(dim,0.);

  TPZIntPoints &intrule = this->GetIntegrationRule();

  TPZManVector<int,3> prevorder(dim), maxorder(dim, this->MaxOrder());
  intrule.GetOrder(prevorder);
  intrule.SetOrder(maxorder);

  data.p = this->MaxOrder();
  data.HSize = 2.*this->InnerRadius();
  error.Fill(0.);
  int npoints = intrule.NPoints(), ip;
  for(ip=0;ip<npoints;ip++){
    intrule.Point(ip,intpoint,weight);
    this->ComputeShape(intpoint, data.x, data.jacobian, data.axes, data.detjac, data.jacinv, data.phi, data.dphix);
    weight *= fabs(data.detjac);
    this->ComputeSolution(intpoint, data.phi, data.dphix, data.axes, data.sol, data.dsol);
    material->ContributeErrors(data,weight,error,errorid);
  }
  intrule.SetOrder(prevorder);
}

void TPZInterpolationSpace::Integrate(int variable, TPZVec<REAL> & value){
  TPZAutoPointer<TPZMaterial> material = Material();
  if(!material){
    PZError << "Error at " << __PRETTY_FUNCTION__ << " : no material for this element\n";
    return;
  }
  if (!this->Reference()){
    PZError << "Error at " << __PRETTY_FUNCTION__ << " : no reference element\n";
    return;
  }
  const int dim = this->Dimension();
  REAL weight;
  TPZMaterialData data;
  this->InitMaterialData(data);

  TPZManVector<REAL, 3> intpoint(dim,0.);
  const int varsize = material->NSolutionVariables(variable);
  value.Resize(varsize);
  value.Fill(0.);

  TPZIntPoints &intrule = this->GetIntegrationRule();
  int npoints = intrule.NPoints(), ip, iv;
  for(ip=0;ip<npoints;ip++){
    intrule.Point(ip,intpoint,weight);
    data.sol.Fill(0.);
    this->Solution(intpoint, variable, data.sol);
    //Tiago: Next call is performet only for computing detcaj. The previous method (Solution) has already computed jacobian.
    //       It means that the next call would not be necessary if I write the whole code here.
    this->Reference()->Jacobian(intpoint, data.jacobian, data.axes, data.detjac, data.jacinv);
    weight *= fabs(data.detjac);
    for(iv = 0; iv < varsize; iv++){
      value[iv] += data.sol[iv]*weight;
    }//for iv
  }//for ip
}//method

void TPZInterpolationSpace::ProjectFlux(TPZElementMatrix &ek, TPZElementMatrix &ef) {

  TPZAutoPointer<TPZMaterial> material = Material();
  if(!material){
    stringstream sout;
    sout << "Exiting ProjectFlux: no material for this element\n";
    Print(sout);
    LOGPZ_ERROR(logger,sout.str());
    ek.Reset();
    ef.Reset();
    return;
  }

  int num_flux = material->NFluxes();
  int dim = Dimension();
  int nshape = NShapeF();
  int ncon = NConnects();
  TPZIntPoints &intrule = GetIntegrationRule();

  int numeq = nshape;
  ek.fMat.Resize(numeq,numeq);
  ek.fBlock.SetNBlocks(ncon);
  ef.fMat.Resize(numeq,num_flux);
  ef.fBlock.SetNBlocks(ncon);

  for(int i=0; i<ncon; ++i){
    (ef.fConnect)[i] = ConnectIndex(i);
    (ek.fConnect)[i] = ConnectIndex(i);
  }

  TPZMaterialData data;
  this->InitMaterialData(data);

  TPZManVector<REAL> flux(num_flux,1);
  TPZManVector<REAL,3> intpoint(dim);
  REAL weight = 0.;
  for(int int_ind = 0; int_ind < intrule.NPoints(); ++int_ind){

    intrule.Point(int_ind,intpoint,weight);
    this->ComputeShape(intpoint, data.x, data.jacobian, data.axes, data.detjac, data.jacinv, data.phi, data.dphix);
    weight *= fabs(data.detjac);
    this->ComputeSolution(intpoint, data.phi, data.dphix, data.axes, data.sol, data.dsol);

    material->Flux(data.x,data.sol,data.dsol,data.axes,flux);
    for(int in=0; in<nshape; in++){
      for(int ifl=0; ifl<num_flux; ifl++){
        (ef.fMat)(in,ifl) += flux[ifl]*data.phi(in,0)*weight;
      }//for ifl
      for(int jn = 0; jn<nshape; jn++){
        (ek.fMat)(in,jn) += data.phi(in,0)*data.phi(jn,0)*weight;
      }//for jn
    }//for in
  }//for int_ind
}//method

/**
 * Save the element data to a stream
 */
void TPZInterpolationSpace::Write(TPZStream &buf, int withclassid)
{
  TPZCompEl::Write(buf,withclassid);
  buf.Write(&fPreferredOrder,1);
}
void TPZInterpolationSpace::MinMaxSolutionValues(TPZVec<REAL> &min, TPZVec<REAL> &max){

  const int dim = Dimension();
  TPZManVector<REAL,3> intpoint(dim,0.);

  TPZIntPoints &intrule = GetIntegrationRule();
  TPZManVector<int,3> prevorder(dim,0);
  intrule.GetOrder(prevorder);

  TPZManVector<int,3> maxorder(dim,intrule.GetMaxOrder());
  intrule.SetOrder(maxorder);

  TPZManVector<REAL,10> sol;
  TPZFNMatrix<30> dsol;
  TPZFNMatrix<9> axes(3,3,0.);
  REAL weight;

  int intrulepoints = intrule.NPoints();
  intrule.Point(0,intpoint,weight);
  this->ComputeSolution(intpoint, sol, dsol, axes);
  min = sol;
  max = sol;
  const int nvars = sol.NElements();
  for(int int_ind = 1; int_ind < intrulepoints; int_ind++){
    intrule.Point(int_ind,intpoint,weight);
    this->ComputeSolution(intpoint, sol, dsol, axes);
    for(int iv = 0; iv < nvars; iv++){
      if (sol[iv] < min[iv]) min[iv] = sol[iv];
      if (sol[iv] > max[iv]) max[iv] = sol[iv];
    }//iv
  }//loop over integratin points

  intrule.SetOrder(prevorder);

}//void
void TPZInterpolationSpace::BuildTransferMatrix(TPZInterpolationSpace &coarsel, TPZTransform &t, TPZTransfer &transfer){
  // accumulates the transfer coefficients between the current element and the
  // coarse element into the transfer matrix, using the transformation t
  TPZGeoEl *ref = Reference();
  int locnod = NConnects();
  int cornod = coarsel.NConnects();
  int locmatsize = NShapeF();
  int cormatsize = coarsel.NShapeF();

  // compare interpolation orders
  // the minimum interpolation order of this needs to be larger than the maximum interpolation order of coarse

  int mymaxorder = MaxOrder();
  int ic;
  int coarsemaxorder = this->MaxOrder();
  if(coarsemaxorder > mymaxorder) {
    stringstream sout;
    sout << "Exiting BuildTransferMatrix - compute the transfer matrix coarse "
        << coarsemaxorder << " me " << mymaxorder << endl;
    LOGPZ_ERROR(logger,sout.str());
    return;
  }
  TPZStack<int> connectlistcoarse,dependencyordercoarse, corblocksize;
  connectlistcoarse.Resize(0);
  dependencyordercoarse.Resize(0);
  corblocksize.Resize(0);
  for(ic=0; ic<coarsel.NConnects(); ic++) connectlistcoarse.Push(coarsel.ConnectIndex(ic));
  coarsel.BuildConnectList(connectlistcoarse);
  TPZConnect::BuildDependencyOrder(connectlistcoarse,dependencyordercoarse,*coarsel.Mesh());

  // cornod = number of connects associated with the coarse element
  cornod = connectlistcoarse.NElements();
  int nvar = coarsel.Material()->NStateVariables();

  // number of blocks is cornod
  TPZBlock corblock(0,cornod);
  int in;

  cormatsize = 0;
  int c;
  for(in=0;in<cornod; in++) {
    c = connectlistcoarse[in];
    int blsize = coarsel.Mesh()->ConnectVec()[c].NDof(*(coarsel.Mesh()))/nvar;
    corblock.Set(in,blsize);
    corblocksize.Push(blsize);
    cormatsize += blsize;
  }
  corblock.Resequence();

//  REAL loclocmatstore[500] = {0.};
  // loclocmat is the inner product of the shape functions of the local element
  // loccormat is the inner product of the shape functions with the shape functions
  //    of the coarse element, both dependent and independent
  TPZFNMatrix<500> loclocmat(locmatsize,locmatsize);
  TPZFNMatrix<500> loccormat(locmatsize,cormatsize);
  loclocmat.Zero();
  loccormat.Zero();

  TPZIntPoints &intrule = GetIntegrationRule();
  int dimension = Dimension();

  TPZManVector<int> prevorder(dimension),order(dimension);
  intrule.GetOrder(prevorder);


  // compute the interpolation order of the shapefunctions squared
  int dim;
  for(dim=0; dim<dimension; dim++) {
    order[dim] = mymaxorder*2;
  }
  intrule.SetOrder(order);

  TPZBlock locblock(0,locnod);

  for(in = 0; in < locnod; in++) {
    locblock.Set(in,NConnectShapeF(in));
  }
  locblock.Resequence();

  REAL locphistore[50]={0.},locdphistore[150]={0.};
  TPZFMatrix locphi(locmatsize,1,locphistore,50);
  TPZFMatrix locdphi(dimension,locmatsize,locdphistore,150);
  locphi.Zero();
  locdphi.Zero();
  // derivative of the shape function
  // in the master domain

  TPZFMatrix corphi(cormatsize,1);
  TPZFMatrix cordphi(dimension,cormatsize);
  // derivative of the shape function
  // in the master domain

  REAL jacobianstore[9],
  axesstore[9];
  TPZManVector<REAL> int_point(dimension),
  coarse_int_point(dimension);
  TPZFMatrix jacobian(dimension,dimension,jacobianstore,9),jacinv(dimension,dimension);
  TPZFMatrix axes(3,3,axesstore,9);
  TPZManVector<REAL> x(3);
  int_point.Fill(0.,0);
  REAL jac_det = 1.;
  ref->Jacobian( int_point, jacobian , axes, jac_det, jacinv);
  REAL multiplier = 1./jac_det;

  int numintpoints = intrule.NPoints();
  REAL weight;
  int lin,ljn,cjn;

  for(int int_ind = 0; int_ind < numintpoints; ++int_ind) {

    intrule.Point(int_ind,int_point,weight);
    ref->Jacobian( int_point, jacobian , axes, jac_det, jacinv);
    ref->X(int_point, x);
    Shape(int_point,locphi,locdphi);
    weight *= jac_det;
    t.Apply(int_point,coarse_int_point);
    corphi.Zero();
    cordphi.Zero();
    coarsel.Shape(coarse_int_point,corphi,cordphi);

    coarsel.ExpandShapeFunctions(connectlistcoarse,dependencyordercoarse,corblocksize,corphi,cordphi);

    for(lin=0; lin<locmatsize; lin++) {
      for(ljn=0; ljn<locmatsize; ljn++) {
        loclocmat(lin,ljn) += weight*locphi(lin,0)*locphi(ljn,0)*multiplier;
      }
      for(cjn=0; cjn<cormatsize; cjn++) {
        loccormat(lin,cjn) += weight*locphi(lin,0)*corphi(cjn,0)*multiplier;
      }
    }
    jacobian.Zero();
  }
  loclocmat.SolveDirect(loccormat,ELDLt);


  for(in=0; in<locnod; in++) {
    //    int cind = connectlistcoarse[in];
    if(Connect(in).HasDependency()) continue;
    int locblocknumber = Connect(in).SequenceNumber();
    int locblocksize = locblock.Size(in);
    int locblockpos = locblock.Position(in);
    TPZStack<int> locblockvec;
    TPZStack<int> globblockvec;
    int numnonzero = 0,jn;
    //      if(transfer.HasRowDefinition(locblocknumber)) continue;

    for(jn = 0; jn<cornod; jn++) {
      int corblocksize = corblock.Size(jn);
      int corblockpos = corblock.Position(jn);
      int cind = connectlistcoarse[jn];
      TPZConnect &con = coarsel.Mesh()->ConnectVec()[cind];
      if(con.HasDependency()) continue;
      int corblocknumber = con.SequenceNumber();
      if(locblocksize == 0 || corblocksize == 0) continue;
      TPZFMatrix small(locblocksize,corblocksize,0.);
      loccormat.GetSub(locblockpos,corblockpos,
                       locblocksize,corblocksize,small);
      REAL tol = Norm(small);
      if(tol >= 1.e-10) {
        locblockvec.Push(jn);
        globblockvec.Push(corblocknumber);
        numnonzero++;
      }
    }
    if(transfer.HasRowDefinition(locblocknumber)) continue;
    transfer.AddBlockNumbers(locblocknumber,globblockvec);
    int jnn;
    for(jnn = 0; jnn<numnonzero; jnn++) {
      jn = locblockvec[jnn];
      int corblocksize = corblock.Size(jn);
      int corblockpos = corblock.Position(jn);
      if(corblocksize == 0 || locblocksize == 0) continue;
      TPZFMatrix small(locblocksize,corblocksize,0.);
      loccormat.GetSub(locblockpos,corblockpos,locblocksize,corblocksize,small);
      transfer.SetBlockMatrix(locblocknumber,globblockvec[jnn],small);
    }
  }
  intrule.SetOrder(prevorder);
}


void TPZInterpolationSpace::ExpandShapeFunctions(TPZVec<int> &connectlist, TPZVec<int> &dependencyorder, TPZVec<int> &blocksizes, TPZFMatrix &phi, TPZFMatrix &dphix) {
  int numblocks =  connectlist.NElements();
  TPZCompMesh &mesh = *Mesh();
  int nhandled=0;
  int current_order = 0;
  int current_block =0;
  while(nhandled < numblocks) {
    if(dependencyorder[current_block] == current_order) {
      nhandled++;
      int cind = connectlist[current_block];
      TPZConnect &con = mesh.ConnectVec()[cind];
      con.ExpandShape(cind,connectlist,blocksizes,phi,dphix);
    }
    current_block++;
    if(current_block == numblocks) {
      current_block = 0;
      current_order++;
    }
  }
}
/**
 * Read the element data from a stream
 */
void TPZInterpolationSpace::Read(TPZStream &buf, void *context)
{
  TPZCompEl::Read(buf,context);
  buf.Read(&fPreferredOrder,1);
}
