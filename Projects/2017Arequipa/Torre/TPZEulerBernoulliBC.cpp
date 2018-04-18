#include "TPZEulerBernoulliBC.h"
#include "TPZMaterial.h"
#include "pzelmat.h"

TPZEulerBernoulliBC::TPZEulerBernoulliBC()
	:TPZCompEl(),fPropertyData(),
    fConnectIndex(-1), fBCVal(), fMasses(){
  //nothing here
}

TPZEulerBernoulliBC::~TPZEulerBernoulliBC(){
  //nothing here
}

TPZEulerBernoulliBC::TPZEulerBernoulliBC(TPZCompMesh &mesh, TPZGeoEl *gel, int64_t &index)
     : TPZCompEl(mesh,gel,index),fPropertyData(),
       fConnectIndex(-1), fBCVal(), fMasses(){
  if(gel->Dimension() != 0) DebugStop();
  TPZCompElSide neigh = this->FindNeighbourCompEl();
  if(neigh.Element()){
    this->fConnectIndex = neigh.Element()->ConnectIndex( neigh.Side() );
  }
  else{//criar connect
	  int order = 1;
		int nshape = 1;
		int nstate = 6;

      const int newnodeindex = mesh.AllocateNewConnect(nshape,nstate,order);
      TPZConnect &newnod = mesh.ConnectVec()[newnodeindex];
      if(newnod.HasDependency()) DebugStop();
      const int seqnum = newnod.SequenceNumber();
      newnod.SetOrder(1,newnodeindex); //ordem p = 1 ? nao, mas uma funcao por no por variavel (3 deslocamentos e 3 giros)
      mesh.Block().Set(seqnum,6);//6 equacoes por no (e connect)
	  this->fConnectIndex = newnodeindex;
  }
  gel->SetReference(this);
}

TPZEulerBernoulliBC::TPZEulerBernoulliBC(const TPZEulerBernoulliBC &cp)
    : TPZCompEl(cp){
  this->fPropertyData = cp.fPropertyData;
  this->fConnectIndex = cp.fConnectIndex;
  this->fBCVal = cp.fBCVal;
  this->fMasses = cp.fMasses;
}

TPZCompElSide TPZEulerBernoulliBC::FindNeighbourCompEl() const{
  TPZGeoElSide thisGeo(this->Reference(),0);
  TPZGeoElSide neigh = thisGeo.Neighbour();
  while( thisGeo != neigh ){
    TPZCompEl * cel = neigh.Element()->Reference();
    if(cel) {
      return TPZCompElSide(cel, neigh.Side());
    }
    else {
      neigh = neigh.Neighbour();
    }
  }
  return TPZCompElSide(NULL,-1);
}

void TPZEulerBernoulliBC::StiffnessMatrix(TPZFMatrix<STATE> &K, TPZFMatrix<STATE> &F) {
  if(fBCVal.NElements() == 0) return;//no forces or supports
  if(fBCVal.NElements() != 6) DebugStop();
  for(int ibc = 0; ibc < 6; ibc++){
    BCVal bc = fBCVal[ibc];
    if(bc.fType == ENone){
      continue; //nothing to do
    }
    else if(bc.fType == EForce){
      F.PutVal(ibc,0, bc.fVal);
    }//force
    else if(bc.fType == ESupport){
      K.PutVal(ibc,ibc, TPZMaterial::gBigNumber);
      F.PutVal(ibc,0, bc.fVal*TPZMaterial::gBigNumber);
    }
    else {
      DebugStop();
    }
  }//for i

  //nonlinear context:
  /// F = F - Ku
  TPZFNMatrix<6> u(6,1), aux(6,1);
  this->GetSolutionVector(u);
  K.Multiply(u,aux,0);//aux = K.u
  F -= aux;
}

void TPZEulerBernoulliBC::MassMatrix(TPZFMatrix<STATE> &M) const{
  if(this->fMasses.NElements() == 0) return;//no nodal mass
  if(this->fMasses.NElements() != 3) DebugStop();//nodal mass needs 3 components x,y,z
  for(int i = 0; i < 3; i++){
    const REAL val = this->fMasses[i];
    M.PutVal(i,0, val );
  }
}

void TPZEulerBernoulliBC::GetSolutionVector(TPZFMatrix<STATE> &u) {
  const int whichSol = 0;

  TPZBlock<STATE> &block = Mesh()->Block();
  TPZFMatrix<STATE> &MeshSol = Mesh()->Solution();
  u.Redim( this->NEquations(), 1 );
  int count = 0;
  for(int ic = 0; ic < this->NConnects(); ic++){
    TPZConnect *df = &Connect(ic);
    int dfseq = df->SequenceNumber();
    int dfvar = block.Size(dfseq);
    int pos = block.Position(dfseq);
    for(int jn = 0; jn < dfvar; jn++) {
      u(count,0) = MeshSol(pos+jn,whichSol);
      count++;
    }
  }

  if(count != this->NEquations()) {
    DebugStop();
  }
}

void TPZEulerBernoulliBC::Print(std::ostream & out) const{
  out << "TPZEulerBernoulliBC\n";
  TPZCompEl::Print(out);
  out << "fConnectIndex: " << fConnectIndex << "\n";

  if(fBCVal.NElements()){
    if(fBCVal.NElements() != 6) DebugStop();
    out << "fBCVal:\n";
    for(int i = 0; i < 6; i++){
      out << fBCVal[i].fType << "\t" << fBCVal[i].fVal << "\n";
    }
  }
  else{
    out << "fBCVal.NElements = 0 \n";
  }

  if(fMasses.NElements()){
    if(fMasses.NElements() != 3) DebugStop();
    out << "fMasses: ";
    for(int i = 0; i < 3; i++){
      out << fMasses[i] << "\t" ;
    }
    out << "\n";
  }
  else{
    out << "fMasses.NElements = 0 \n";
  }

  out << "fPropertyData = " << fPropertyData.operator->() << "\n";
}

static int tamarindo = 0;
void TPZEulerBernoulliBC::CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef){

  const int ncon = this->NConnects();//1
  const int numdof = 6;
  const int numeq = this->NEquations();//6

  ek.fMesh = Mesh();
  ek.fType = TPZElementMatrix::EK;
  ef.fMesh = Mesh();
  ef.fType = TPZElementMatrix::EF;

  ek.fMat.Redim(numeq,numeq);
  ef.fMat.Redim(numeq,1);
  ek.fBlock.SetNBlocks(ncon);
  ef.fBlock.SetNBlocks(ncon);
  ek.fNumStateVars = numdof;
  ef.fNumStateVars = numdof;
  for(int i = 0; i < ncon ; i++){
    ek.fBlock.Set(i,6);
    ef.fBlock.Set(i,6);
  }
  ek.fConnect.Resize(ncon);
  ef.fConnect.Resize(ncon);
  for(int i = 0; i < ncon; i++){
    ef.fConnect[i] = this->ConnectIndex(i);
    ek.fConnect[i] = this->ConnectIndex(i);
  }

  if(tamarindo){
    return;//tamarindo
  }

  if(this->fPropertyData->ComputeStiffness()){
    this->StiffnessMatrix(ek.fMat,ef.fMat);
  }
  else if(this->fPropertyData->ComputeMass()){
    this->MassMatrix(ek.fMat);
  }
  else DebugStop();

}

void TPZEulerBernoulliBC::Solution(TPZVec<REAL> & /*qsi*/,int /*var*/,TPZVec<REAL> &/*sol*/){
  DebugStop(); //implement me
}

void TPZEulerBernoulliBC::Write(TPZStream &/*buf*/, int /*withclassid*/){
  DebugStop(); //implement me
}

void TPZEulerBernoulliBC::Read(TPZStream &/*buf*/, void * /*context*/){
  DebugStop(); //implement me
}



void TPZEulerBernoulliBC::GetEquationIndices(TPZVec<int> &indices) const{
  TPZBlock<STATE> &block = Mesh()->Block();
  indices.Resize( 6 );
  int count = 0;

  TPZConnect *df = &Connect(0);
  int dfseq = df->SequenceNumber();
  int dfvar = block.Size(dfseq);
  int pos = block.Position(dfseq);
  for(int jn = 0; jn < dfvar; jn++) {
    indices[count] = pos+jn;
    count++;
  }
  if(count != 6) {
    DebugStop();
  }
}

