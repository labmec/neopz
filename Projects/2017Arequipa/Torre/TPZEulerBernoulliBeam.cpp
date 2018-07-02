#define M_PI       3.14159265358979323846


#include "TPZEulerBernoulliBeam.h"
#include "pzelmat.h"

bool TPZEulerBernoulliBeam::SimplifiedCorotational() const{
  return false;
}

TPZEulerBernoulliBeam::TPZEulerBernoulliBeam()
 : TPZCompEl(), fMaterialId(-1),
   fSectionId(-1),
   fAlfa(0.), fConnectIndexes(2,-1), fFabricationErrorStrain(0.) {
  //nothing here
}

TPZEulerBernoulliBeam::~TPZEulerBernoulliBeam(){
  for(int i = 0; i < NConnects(); i++){
    this->Mesh()->ConnectVec()[this->ConnectIndex(i)].DecrementElConnected();
  }
  if(this->Reference()->Reference() == this){
    this->Reference()->ResetReference();
  }
}

TPZEulerBernoulliBeam::TPZEulerBernoulliBeam(TPZCompMesh &mesh, TPZGeoEl *gel, int64_t &index)
      : TPZCompEl(mesh,gel,index),
        fMaterialId(-1), fSectionId(-1), fAlfa(0.), fFabricationErrorStrain(0.){
  if(gel->Dimension() != 1) DebugStop();
  this->fConnectIndexes.Resize(2);
  for(int iside = 0; iside < 2; iside++){
    TPZCompElSide neigh = this->FindNeighbourCompEl(iside);
    if(neigh.Element()){
      this->fConnectIndexes[iside] = neigh.Element()->ConnectIndex( neigh.Side() );
    }
    else{//criar connect
		int nshape = 1;
		int nstate = 6;
		int order = 1;

      const int newnodeindex = mesh.AllocateNewConnect(nshape,nstate,order);
      TPZConnect &newnod = mesh.ConnectVec()[newnodeindex];
      if(newnod.HasDependency()) DebugStop();
      const int seqnum = newnod.SequenceNumber();
      newnod.SetOrder(1,newnodeindex); //ordem p = 1 ? nao, mas uma funcao por no por variavel (3 deslocamentos e 3 giros)
      mesh.Block().Set(seqnum,6);//6 equacoes por no (e connect)
      this->fConnectIndexes[iside] = newnodeindex;
    }
  }//for iside
  gel->SetReference(this);
}

TPZEulerBernoulliBeam::TPZEulerBernoulliBeam(const TPZEulerBernoulliBeam &cp)
  : TPZCompEl( cp ){
  this->fMaterialId = cp.fMaterialId;
  this->fSectionId = cp.fSectionId;
  this->fConnectIndexes = cp.fConnectIndexes;
  this->fAlfa = cp.fAlfa;
  this->fFabricationErrorStrain = cp.fFabricationErrorStrain;
  this->fPropertyData = cp.fPropertyData;
}

TPZCompElSide TPZEulerBernoulliBeam::FindNeighbourCompEl(int myside) const{
  TPZGeoElSide thisGeo(this->Reference(),myside);
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

void TPZEulerBernoulliBeam::SetPropertyData( TPZAutoPointer< TPZEulerBernoulliBeamData > PropertyData,
                                             int MaterialId, int SectionId, REAL alfa,
                                             REAL fabErrorStrain ){
  this->fPropertyData = PropertyData;
  this->fMaterialId = MaterialId;
  this->fSectionId = SectionId;
  this->fAlfa = alfa;
  this->fFabricationErrorStrain = fabErrorStrain;
}

void TPZEulerBernoulliBeam::Print(std::ostream & out) const{
  out << "TPZEulerBernoulliBeam\n";
  TPZCompEl::Print(out);
  out << "fMaterialId: " << fMaterialId << "\n";
  out << "fSectionId: " << fSectionId << "\n";
  out << "fAlfa: " << fAlfa << "\n";
  out << "fFabricationErrorStrain: " << fFabricationErrorStrain << "\n";
  out << "fConnectIndexes: ";
  for(int i = 0; i < fConnectIndexes.NElements(); i++) out << fConnectIndexes[i] << "\t";
  out << "\n";
  out << "fPropertyData = " << fPropertyData.operator->() << "\n";
}

void TPZEulerBernoulliBeam::CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef){
  const int ncon = this->NConnects();//2
  const int numdof = 6;
  const int numeq = this->NEquations();//12

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

  if(this->fPropertyData->ComputeStiffness()){
    this->StiffnessMatrix(ek.fMat,ef.fMat);
  }
  else if(this->fPropertyData->ComputeMass()){
    this->MassMatrix(ek.fMat);
  }
  else DebugStop();

}

REAL TPZEulerBernoulliBeam::GetStaticNormalForce(){
  TPZFNMatrix<12> NodalF;
  this->GetStaticNodalForces(NodalF);
  const REAL Na = -NodalF(0,0);
  const REAL Nb = NodalF(6,0);
  return (Na + Nb)/2.;
}

void TPZEulerBernoulliBeam::Solution(TPZVec<REAL> &/*qsi*/,int /*var*/,TPZVec<REAL> &/*sol*/){
  DebugStop();
//please implement me
}

void TPZEulerBernoulliBeam::Write(TPZStream &/*buf*/, int /*withclassid*/){
  DebugStop();
  //please, implement me
}

void TPZEulerBernoulliBeam::Read(TPZStream &/*buf*/, void * /*context*/){
  DebugStop();
  //please, implement me
}

REAL TPZEulerBernoulliBeam::LOriginal() const{
  TPZGeoEl * gel = this->Reference();
  TPZManVector<REAL,3> xi(3), xf(3);
  gel->NodePtr(0)->GetCoordinates(xi);
  gel->NodePtr(1)->GetCoordinates(xf);
  REAL sum = 0.;
  for(int i = 0; i < 3; i++){
    const REAL d = xf[i]-xi[i];
    sum += d*d;
  }
  const REAL Lorig = sqrt(sum);
  return Lorig;
}

void TPZEulerBernoulliBeam::L(const TPZFMatrix<STATE> &u, REAL & Lorig, REAL & Ldef) const{
  TPZGeoEl * gel = this->Reference();
  TPZManVector<REAL,3> xi(3), xf(3);
  gel->NodePtr(0)->GetCoordinates(xi);
  gel->NodePtr(1)->GetCoordinates(xf);
  REAL sum = 0.;
  for(int i = 0; i < 3; i++){
    const REAL d = xf[i]-xi[i];
    sum += d*d;
  }
  Lorig = sqrt(sum);

  if(this->SimplifiedCorotational()){
    sum = 0.;
    for(int i = 0; i < 3; i++){
      const REAL d = ( xf[i] + u.GetVal(i+6,0) ) - ( xi[i] + u.GetVal(i,0) );
      sum += d*d;
    }
    Ldef = sqrt(sum);
  }
  else{
    Ldef = Lorig;
  }
}

void TPZEulerBernoulliBeam::LocalStiffnessMatrix(const TPZFMatrix<STATE> & u,
	TPZFMatrix<STATE> &K, TPZFMatrix<STATE> &F, TPZFMatrix<STATE> &Rot) {
  K.Redim(12,12);//Redim faz Zero();

  //properties
  TPZEulerBernoulliBeamData::MaterialProperties mat;
  TPZEulerBernoulliBeamData::SectionProperties sec;
  this->fPropertyData->GetData(fMaterialId, fSectionId, mat, sec);
  const REAL E = mat.fE;
  const REAL G = mat.fG;
  const REAL A = sec.fA;
  const REAL J = sec.fJt;
  const REAL Iy = sec.fIy;
  const REAL Iz = sec.fIz;
  //atualizar com a solucao Laxial = original Lmomento = Latualizado
  REAL Lorig, Ldef;
  this->L(u, Lorig, Ldef);

  //local matrix
  REAL value;
  value = E * A / Lorig;
  K.PutVal(0,0,value);
  K.PutVal(0,6,-1.0*value);
  K.PutVal(6,0,-1.0*value);
  K.PutVal(6,6,value);

  value = 12.0 * E * Iz / ( Ldef*Ldef*Ldef );
  K.PutVal(1,1,value);
  K.PutVal(1,7,-1.0*value);
  K.PutVal(7,1,-1.0*value);
  K.PutVal(7,7,value);

  value = 6.0 * E * Iz / ( Ldef*Ldef );
  K.PutVal(1,5,value);
  K.PutVal(1,11,value);
  K.PutVal(7,5,-1.0*value);
  K.PutVal(7,11,-1.0*value);

  K.PutVal(5,1,value);
  K.PutVal(11,1,value);
  K.PutVal(5,7,-1.0*value);
  K.PutVal(11,7,-1.0*value);

  value = 12.0 * E * Iy / (Ldef*Ldef*Ldef);
  K.PutVal(2,2,value);
  K.PutVal(2,8,-1.0*value);
  K.PutVal(8,2,-1.0*value);
  K.PutVal(8,8,value);

  value = -6.0 * E * Iy / (Ldef*Ldef);
  K.PutVal(2,4,value);
  K.PutVal(2,10,value);
  K.PutVal(8,4,-1.0*value);
  K.PutVal(8,10,-1.0*value);

  K.PutVal(4,2,value);
  K.PutVal(10,2,value);
  K.PutVal(4,8,-1.0*value);
  K.PutVal(10,8,-1.0*value);

  value = G * J / Lorig;
  K.PutVal(3,3,value);
  K.PutVal(3,9,-1.0*value);
  K.PutVal(9,3,-1.0*value);
  K.PutVal(9,9,value);

  value = 4.0 * E * Iy / Ldef;
  K.PutVal(4,4,value);
  K.PutVal(4,10,value/2.0);
  K.PutVal(10,4,value/2.0);
  K.PutVal(10,10,value);

  value = 4.0 * E * Iz / Ldef;
  K.PutVal(5,5,value);
  K.PutVal(5,11,value/2.0);
  K.PutVal(11,5,value/2.0);
  K.PutVal(11,11,value);

  //////// Rotation matrix /////////////
  this->RotateMatrix(Rot, u, this->SimplifiedCorotational());

  ////////////// Load vector - self-weight /////////////////////////////////////////

  F.Redim(12,1);//Redim faz Zero()

  TPZFNMatrix<3> loadGlobal, loadLocal;
  this->fPropertyData->g(loadGlobal);
  loadGlobal *= mat.fRho * sec.fA;
  TPZFNMatrix<9> Rshort(3,3);
  for(int i = 0; i < 3; i++) for(int j = 0; j < 3; j++) Rshort(i,j) = Rot(i,j);
  Rshort.Multiply(loadGlobal,loadLocal);

  //forces in lX direction
  REAL q;
  q = loadLocal(0.0);
  value = q * Lorig / 2.0;
  F.PutVal(0,0, value);
  F.PutVal(6,0, value);

  //forces in lY direction
  q = loadLocal(1,0);
  value = q * Lorig / 2.0;
  F.PutVal(1,0, value);
  F.PutVal(7,0, value);

  value = q * Lorig * Lorig / 12.0;
  F.PutVal(5,0, value);
  F.PutVal(11,0, -1.*value);

  //forces in lZ direction
  q = loadLocal(2,0);
  value = q * Lorig / 2.0;
  F.PutVal(2,0, value);
  F.PutVal(8,0, value);

  value = -1.0 * (q * Lorig * Lorig ) / 12.0;
  F.PutVal(4,0, value);
  F.PutVal(10,0, -1.*value);

  //fabrication error
  value = E * A * this->fFabricationErrorStrain;
  F(0,0) += -value;
  F(6,0) += value;

}//local stiffness

void TPZEulerBernoulliBeam::StiffnessMatrix(TPZFMatrix<STATE> &K, TPZFMatrix<STATE> &F) {

  TPZFNMatrix<12> uGlobal(12,1);
  this->GetSolutionVector(uGlobal);

  TPZFNMatrix<144> localK(12,12,0.), localF(12,1,0.), Rot(12,12,0.);
  this->LocalStiffnessMatrix(uGlobal, localK, localF, Rot);

  ////////////// Load vector - self-weight /////////////////////////////////////////

  //girando para sistema global
  Rot.Multiply(localF,F,1);   //F = Rt . localF

  /// F = F - Ku
  if(this->SimplifiedCorotational() == true){//SimplifiedCorotational

    TPZFNMatrix<12> uLocal(12,1), ElasticForces(12,1);
    TPZFNMatrix<144> twistRot(12,12);
    REAL twistAngle;
    this->TransformDisplFromGlobalToLocalRefSystem(uGlobal, Rot, uLocal, twistAngle, twistRot);

    TPZFNMatrix<144> aux(12,12), aux2(12,12);

    //computing elastic forces in local ref system
    localK.Multiply(uLocal,aux);//aux = localK.deltaU

    TPZFNMatrix<12> auxTwisted(12,1);
    twistRot.Multiply(aux,auxTwisted,1);//auxTwisted = Transpose[Rtwist].aux

    //from local do global ref system
    Rot.Multiply(auxTwisted,ElasticForces,1);  //ElasticF = Rt . localK.deltaU
    F -= ElasticForces;

    //applying Rtwist to jacobian matrix
    //K = Transp[R] . Transp[Rtwist] . localK . Rtwist . R
    twistRot.Multiply(Rot,aux);//aux = Rtwist . R
    localK.Multiply(aux,aux2);//aux2 = localK . Rtwist . R
    twistRot.Multiply(aux2,aux,1);//aux = Transp[Rtwist] . localK . Rtwist . R
    Rot.Multiply(aux,K,1); //K = Transp[R] . Transp[Rtwist] . localK . Rtwist . R

  }
  else{//linear formulation
    TPZFNMatrix<144> aux(12,12);

  //   B = this * X
  //   If opt = 1 then B = Transpose[this] * X
  //  void ConstMultiply(const TPZFMatrix & x,TPZFMatrix & B,const int opt = 0) const;
    localK.Multiply(Rot,aux,0); //aux = localK.R
    Rot.Multiply(aux,K,1);   //K = Rt . localK . R

    K.Multiply(uGlobal,aux,0);//aux = K.u (everything in global ref system)
    F -= aux;
  }

}//local stiffness

void TPZEulerBernoulliBeam::GetInitialPosVec(TPZFMatrix<STATE> &InitialPos){
  InitialPos.Redim(12,1); //Redim zeroes the matrix
  TPZGeoEl * gel = this->Reference();
  TPZManVector<REAL,3> coord(3);
  gel->NodePtr(0)->GetCoordinates(coord);
  for(int i = 0; i < 3; i++) InitialPos(i,0) = coord[i];
  gel->NodePtr(1)->GetCoordinates(coord);
  for(int i = 0; i < 3; i++) InitialPos(6+i,0) = coord[i];
}

void TPZEulerBernoulliBeam::RotateMatrix(TPZFMatrix<STATE> &Rotate, const TPZFMatrix<STATE> & u, bool uComposeRotation) {

  TPZGeoEl * gel = this->Reference();
  TPZManVector<REAL,3> xi(3), xf(3);
  gel->NodePtr(0)->GetCoordinates(xi);
  gel->NodePtr(1)->GetCoordinates(xf);

  REAL deltaX = xf[0] - xi[0];
  REAL deltaY = xf[1] - xi[1];
  REAL deltaZ = xf[2] - xi[2];
  REAL L = sqrt(deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ);

  //aqui fica faltando a parcela da derivada da energia de deformacao pela rotacao. Quando EI=0, essa parcela vale zero. Pra cabos,
  //portanto, eh uma boa aproximacao. Para barras de inercia grande, entretanto, a formulacao esta incompleta.
  if(uComposeRotation){
    deltaX += u.GetVal(6,0) - u.GetVal(0,0);
    deltaY += u.GetVal(7,0) - u.GetVal(1,0);
    deltaZ += u.GetVal(8,0) - u.GetVal(2,0);
    L = sqrt(deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ);
  }

//R = Ralfa*Rgamma*Rbeta

  //Compute Ralfa
  REAL cosaux = cos(this->fAlfa);
  REAL sinaux = sin(this->fAlfa);
  Rotate.Resize(12,12);
  Rotate.Zero();
//        R 0 0 0
//Ralfa = 0 R 0 0
//        0 0 R 0
//        0 0 0 R
  for(int i = 0;  i < 4; i++){
    int j = i * 3;
    Rotate.PutVal(j,j,1.0);
    Rotate.PutVal(j+1,j+1,cosaux);
    Rotate.PutVal(j+1,j+2,sinaux);
    Rotate.PutVal(j+2,j+1,-1.0*sinaux);
    Rotate.PutVal(j+2,j+2,cosaux);
  }
//now Rotate = Ralfa

//Compute Rgamma
  cosaux = sqrt(deltaX * deltaX + deltaY * deltaY) / L;
  sinaux = deltaZ / L;

  TPZFNMatrix<144> Raux(12,12,0.0);
//         R 0 0 0
//Rgamma = 0 R 0 0
//         0 0 R 0
//         0 0 0 R
  for(int i = 0; i < 4; i++){
    int j = i * 3;
    Raux.PutVal(j,j,cosaux);
    Raux.PutVal(j,j+2,sinaux);
    Raux.PutVal(j+2,j,-1.0*sinaux);
    Raux.PutVal(j+2,j+2,cosaux);
    Raux.PutVal(j+1,j+1,1.0);
  }

  TPZFNMatrix<144> Raux2 = Rotate;
  Raux2.Multiply(Raux,Rotate);//now Rotate = Ralfa * Rgamma

  //Compute Rbeta
  Raux.Zero();
  if ( fabs(deltaX) < 1e-10 && fabs(deltaY) < 1e-10  ){
    cosaux = 1.0;
    sinaux = 0.0;
  }
  else{
    const REAL lproj = sqrt(deltaX * deltaX + deltaY * deltaY);
    cosaux = deltaX / lproj;
    sinaux = deltaY / lproj;
  }

  //        R 0 0 0
  //Rbeta = 0 R 0 0
  //        0 0 R 0
  //        0 0 0 R
  for(int i = 0; i < 4; i++){
    int j = i * 3;
    Raux.PutVal(j,j,cosaux);
    Raux.PutVal(j,j+1,sinaux);
    Raux.PutVal(j+1,j,-1.0*sinaux);
    Raux.PutVal(j+1,j+1,cosaux);
    Raux.PutVal(j+2,j+2,1.0);
  }

  Raux2 = Rotate;
  Raux2.Multiply(Raux,Rotate);//now Rotate = Ralfa * Rgamma * Rbeta

}


void TPZEulerBernoulliBeam::LocalMassMatrix(TPZFMatrix<STATE> &M) const{
  M.Redim(12,12);

  //properties
  TPZEulerBernoulliBeamData::MaterialProperties mat;
  TPZEulerBernoulliBeamData::SectionProperties sec;
  this->fPropertyData->GetData(fMaterialId, fSectionId, mat, sec);
  const REAL L = this->LOriginal();
  const REAL ml = mat.fRho * sec.fA * L;
  const REAL Ip = sec.fIp;
  const REAL A = sec.fA;

  //filling values
  REAL value = ml / 3.0;

  M.PutVal(0,0,value);
  M.PutVal(6, 6, value);

  value = ml * 13.0 / 35.0;
  M.PutVal(1,1,value);
  M.PutVal(2,2,value);
  M.PutVal(7,7,value);
  M.PutVal(8,8,value);

  value = ml * Ip / ( 3.0 * A );
  M.PutVal(3,3,value);
  M.PutVal(9,9,value);

  value = ml * L * L / 105.0;
  M.PutVal(4,4,value);
  M.PutVal(5,5,value);
  M.PutVal(10,10,value);
  M.PutVal(11,11,value);
  //acabou a diagonal

  value = ml * 11.0 * L / 210.0;
  M.PutVal(4,2,-value);
  M.PutVal(5,1,value);
  M.PutVal(2,4,-value);
  M.PutVal(1,5,value);

  value = ml / 6.0;
  M.PutVal(6, 0, value);
  M.PutVal(0, 6, value);

  value = ml * 9.0 / 70.0;
  M.PutVal(7,1,value);
  M.PutVal(1,7,value);

  value = ml * 13.0 * L / 420.0;
  M.PutVal(7,5,value);
  M.PutVal(5,7,value);

  value = ml * 9.0 / 70.0;
  M.PutVal(8,2,value);
  M.PutVal(2,8,value);

  value = - ml * 13.0 * L / 420.0;
  M.PutVal(8,4,value);
  M.PutVal(4,8,value);

  value = ml * Ip / ( 6.0 * A );
  M.PutVal(9,3,value);
  M.PutVal(3,9,value);

  value = ml * 13.0 * L / 420.0;
  M.PutVal(10,2,value);
  M.PutVal(2,10,value);

  value = - ml * L * L / 140.0;
  M.PutVal(10, 4, value);
  M.PutVal(4, 10, value);

  value = ml * 11.0 * L / 210.0;
  M.PutVal(10, 8, value);
  M.PutVal(8, 10, value);

  value = - ml * 13.0 * L / 420.0;
  M.PutVal(11, 1, value);
  M.PutVal(1, 11, value);

  value = - ml * L * L / 140.0;
  M.PutVal(11, 5, value);
  M.PutVal(5, 11, value);

  value = - ml * 11.0 * L / 210.0;
  M.PutVal(11, 7, value);
  M.PutVal(7, 11, value);

}

void TPZEulerBernoulliBeam::MassMatrix(TPZFMatrix<STATE> &M){
  TPZFNMatrix<144> localM(12,12,0.), Rot(12,12,0.);
  this->LocalMassMatrix(localM);

  TPZFNMatrix<12> u(12,1);
  this->GetSolutionVector(u);
  this->RotateMatrix(Rot,u,this->SimplifiedCorotational());

  TPZFNMatrix<144> aux(12,12);

//   B = this * X
//   If opt = 1 then B = Transpose[this] * X
//  void ConstMultiply(const TPZFMatrix & x,TPZFMatrix & B,const int opt = 0) const;
  localM.Multiply(Rot,aux,0); //aux = localM.R
  Rot.Multiply(aux,M,1);   //M = Rt . localM . R

}

void TPZEulerBernoulliBeam::GetStaticNodalForces(TPZFMatrix<STATE> & NodalF) {

  TPZFNMatrix< 144 > LocalStiff(12,12,0.), localF(12,1,0.), NodalDispl(12,1,0.),
                     Rot(12,12,0.), uLocal(12,1,0.);
  this->GetSolutionVector(NodalDispl);
  this->LocalStiffnessMatrix(NodalDispl, LocalStiff, localF, Rot);
  TPZFNMatrix<144> twistRot(12,12);
  REAL twistAngle;
  this->TransformDisplFromGlobalToLocalRefSystem(NodalDispl, Rot, uLocal, twistAngle, twistRot);

 if(this->SimplifiedCorotational() == true){//SimplifiedCorotational

    TPZFNMatrix<12> aux;
    LocalStiff.Multiply(uLocal,aux);//ElasticForces = localK.deltaU
    TPZFNMatrix<12> auxTwisted;
    twistRot.Multiply(aux,auxTwisted,1);//auxTwisted = Transpose[Rtwist].aux
    NodalF = auxTwisted;
    NodalF -= localF;
  }
  else{//linear formulation

    ////Multiply: Stiffness * Displacements = NodalForces (everything in local ref sys)
    LocalStiff.Multiply(uLocal, NodalF); //finding AM = NodalForces

    ////The nodal forces are given by: AM = SM.DM + EP (notation of Gere & Weaver)
    ////In our notation NodalForces = LocalStiff * LocalDisplacements - LoadVector
    NodalF -= localF;
  }

}

void TPZEulerBernoulliBeam::GetSolutionVector(TPZFMatrix<STATE> &u) {
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

void TPZEulerBernoulliBeam::GetEquationIndices(TPZVec<int64_t> &indices) {

	TPZBlock<STATE> &block = Mesh()->Block();
  indices.Resize( this->NEquations() );
  int count = 0;
  for(int ic = 0; ic < this->NConnects(); ic++){
    TPZConnect *df = &Connect(ic);
    int dfseq = df->SequenceNumber();
    int dfvar = block.Size(dfseq);
    int pos = block.Position(dfseq);
    for(int jn = 0; jn < dfvar; jn++) {
      indices[count] = pos+jn;
      count++;
    }
  }

  if(count != this->NEquations()) {
    DebugStop();
  }
}

void TPZEulerBernoulliBeam::Divide(int64_t index, TPZVec<int64_t> &subindex, int interpolate){
  if(interpolate) DebugStop();//not implemented yet
  if (this->Mesh()->ElementVec()[index] != this) DebugStop();

  TPZGeoEl * gel = this->Reference();
  TPZManVector<TPZGeoEl*> children;
  gel->Divide( children );
  const int nsons = children.NElements();
  subindex.Resize(nsons);
  for(int i = 0; i < nsons; i++){
    TPZEulerBernoulliBeam * cel = new TPZEulerBernoulliBeam(*Mesh(), children[i], index);
    subindex[i] = index;
    cel->fMaterialId = this->fMaterialId;
    cel->fSectionId = this->fSectionId;
    cel->fAlfa = this->fAlfa;
    cel->fFabricationErrorStrain = this->fFabricationErrorStrain;
    cel->fPropertyData = this->fPropertyData;
  }

  delete this;

}

void TPZEulerBernoulliBeam::Angles(const TPZVec<REAL> &v, TPZVec<REAL> &angles) const{
  const REAL deltaX = v[0];
  const REAL deltaY = v[1];
  const REAL deltaZ = v[2];
  angles.Resize(3);
  REAL tol = 1e-8;
//rotation over x axis
  REAL theta = 0.;
  if(fabs(deltaZ) > tol || fabs(deltaY) > tol) theta = +1.*atan2(deltaZ,deltaY);
  angles[0] = theta;

  //roatation over y axis
  theta = 0.;
  if(fabs(deltaZ) > tol || fabs(deltaX) > tol) theta = -1.*atan2(deltaZ,deltaX);
  angles[1] = theta;

  //roatation over z axis
  theta = 0.;
  if(fabs(deltaY) > tol || fabs(deltaX) > tol) theta = +1.*atan2(deltaY,deltaX);
  angles[2] = theta;
}//void

void TPZEulerBernoulliBeam::ComputeRigidBodyMotionRotations(const TPZFMatrix<STATE> & uGlobal, TPZVec<REAL> & angles) const{

  TPZManVector<REAL,3> vInitial(3), vFinal(3), anglesInitial(3), anglesFinal(3);
  TPZGeoEl * gel = this->Reference();
  TPZManVector<REAL,3> coordA(3), coordB(3);
  gel->NodePtr(0)->GetCoordinates(coordA);
  gel->NodePtr(1)->GetCoordinates(coordB);
  for(int i = 0; i < 3; i++){
    vInitial[i] = coordB[i]-coordA[i];
    vFinal[i] = (coordB[i]+uGlobal.GetVal(6+i,0)) - (coordA[i]+uGlobal.GetVal(i,0));
  }

  angles.Resize(3);
  this->Angles(vInitial, anglesInitial);
  this->Angles(vFinal, anglesFinal);
  for(int i = 0; i < 3; i++) angles[i] = anglesFinal[i]-anglesInitial[i];

  //quando o produto interno das projecoes no plano de rotacao for nulo, o angulo eh zero
  TPZFNMatrix<3> vIniProj(3,1), vFinProj(3,1);
  for(int axis = 0; axis < 3; axis++){
    //projecao dos vetores
    for(int i = 0; i < 3; i++){
      vIniProj(i,0) = vInitial[i];
      vFinProj(i,0) = vFinal[i];
    }
    vIniProj(axis,0) = 0.;
    vFinProj(axis,0) = 0.;

    //normalizando vetores
    REAL tol = 1e-8;
    if(Norm(vIniProj) > tol){
      vIniProj *= 1./Norm(vIniProj);
    }
    if(Norm(vFinProj) > tol){
      vFinProj *= 1./Norm(vFinProj);
    }

    //dot
    REAL IniFinDot = Dot(vIniProj,vFinProj);
    if( fabs(IniFinDot) < tol ){
      angles[ axis ] = 0.;
    }
  }//axis

}///void

void TPZEulerBernoulliBeam::ComputeRotationBetweenTwoVectors(const TPZVec<REAL> &a, const TPZVec<REAL> &b,
	TPZFMatrix<STATE> & RotBetweenVectors) const{
#ifdef DEBUG
  {//checking vectors a and b are normalized
    REAL normA = 0.;
    for(int i = 0; i < 3; i++) normA += a[i]*a[i];
    if( fabs( sqrt(normA) - 1.) > 1e-10 ) DebugStop();
    REAL normB = 0.;
    for(int i = 0; i < 3; i++) normB += b[i]*b[i];
    if( fabs( sqrt(normB) - 1.) > 1e-10 ) DebugStop();
  }
#endif
  REAL v[3] = { -a[2]*b[1] + a[1]*b[2], a[2]*b[0] - a[0]*b[2], -a[1]*b[0] + a[0]*b[1] };
  REAL vMatTemp[3][3] = {{0, -v[2], v[1]}, {v[2], 0, -v[0]}, {-v[1], v[0], 0}};
  TPZFNMatrix<9> vMat(3,3);
  for(int i = 0; i < 3; i++) for(int j = 0; j < 3; j++) vMat(i,j) = vMatTemp[i][j];

  REAL c = 0.; //c = a.b;
  for(int i = 0; i < 3; i++) c += a[i]*b[i];
  if( fabs( c - (-1.) ) < 1e-6 ) DebugStop();//I cannot compute for 180º

  //R = IdentityMatrix[3] + vMat + vMat.vMat 1/(1 + c);
  RotBetweenVectors.Resize(3,3);
  RotBetweenVectors.Identity();
  RotBetweenVectors += vMat;
  TPZFNMatrix<9> vMat2(3,3);
  vMat.Multiply(vMat,vMat2);
  vMat2 *= 1./(1.+c);
  RotBetweenVectors += vMat2;
}

REAL TPZEulerBernoulliBeam::ComputeTwistRigidBodyAngle(const TPZFMatrix<STATE> & RotOrigBig,
	const TPZFMatrix<STATE> & RotDeformedBig) const{
  if(this->SimplifiedCorotational() == false){//linear formulation
    DebugStop();
  }

//  return 0*M_PI/4.;//tamarindo

  TPZFNMatrix<9> RotOrig(3,3), RotDeformed(3,3);
  RotOrigBig.GetSub(0,0,3,3,RotOrig);
  RotDeformedBig.GetSub(0,0,3,3,RotDeformed);


  TPZFNMatrix<9> axisCanonical(3,3), axisOrig(3,3), axisDeformed(3,3);
  axisCanonical.Identity();

  RotOrig.Multiply(axisCanonical,axisOrig,1);
  RotDeformed.Multiply(axisCanonical,axisDeformed,1);
#ifdef DEBUG
   std::ofstream myfile("c:\\Temp\\CEATI\\axis.txt");
   myfile.precision(10);
   axisOrig.Print("orig=",myfile,EMathematicaInput);
   axisDeformed.Print("deformed=",myfile,EMathematicaInput);
   myfile.flush();
#endif


  TPZFNMatrix<3> RotBetweenVectors(3,3);
  {//computing rotation matrix between the two vectors
    TPZManVector<REAL,3> a(3), b(3);
    for(int i = 0; i < 3; i++){
      a[i] = axisOrig(i,0);
      b[i] = axisDeformed(i,0);
    }
    this->ComputeRotationBetweenTwoVectors(a, b, RotBetweenVectors);
#ifdef DEBUG
    myfile << "a= {";
    for(int i = 0; i < 3; i++) myfile << a[i] << ", ";
    myfile << "}\n";
    myfile << "b= {";
    for(int i = 0; i < 3; i++) myfile << b[i] << ", ";
    myfile << "}\n";
#endif
  }

  //rotating original axis to deformed coordinate system
  TPZFNMatrix<9> axisOrigRotated(3,3);
  RotBetweenVectors.Multiply(axisOrig, axisOrigRotated);

#ifdef DEBUG
  RotBetweenVectors.Print("RotBetweenVectors=",myfile,EMathematicaInput);  myfile.flush();
  axisOrigRotated.Print("axisOrigRotated=",myfile,EMathematicaInput); myfile.flush();
#endif

  //taking to master element reference system only to simplify twist angle computation
  TPZFNMatrix<9> OrigFinal(3,3);
  RotDeformed.Multiply(axisOrigRotated, OrigFinal);

#ifdef DEBUG
{
  TPZFNMatrix<9> DeformedFinal(3,3);
  RotDeformed.ConstMultiply(axisDeformed, DeformedFinal);
  OrigFinal.Print("OrigFinal=",myfile,EMathematicaInput); myfile.flush();
  DeformedFinal.Print("DeformedFinal=",myfile,EMathematicaInput); myfile.flush();

  //DeformedFinal must be identity
  DeformedFinal -= axisCanonical;
  REAL n = Norm(DeformedFinal);
  if(n > 1e-6) DebugStop();
}
#endif

  //angulo de rotacao do segundo eixo local (alinhado com y no elemento mestre)
  REAL dx = OrigFinal(0,1);
  REAL dy = OrigFinal(1,1);
  REAL dz = OrigFinal(2,1);
#ifdef DEBUG
  if(fabs(dx) > 1e-6) DebugStop();
#endif
  REAL twistAngle = atan2( dz,dy );

  return -1.*twistAngle;
}

static std::ofstream myfileTwist("myfileTwist.txt");  //tamarindo
void TPZEulerBernoulliBeam::TransformDisplFromGlobalToLocalRefSystem(const TPZFMatrix<STATE> & uGlobal,
	const TPZFMatrix<STATE> & Rot,	TPZFMatrix<STATE> &uLocal, REAL &twistAngle, TPZFMatrix<STATE> &twistRot){

  if(this->SimplifiedCorotational() == true){//SimplifiedCorotational

#ifdef DEBUG
  TPZFNMatrix<12> LocalDisplLinear(12,1,0.);
  Rot.ConstMultiply(uGlobal,LocalDisplLinear,0);
#endif

    //tamarindo
      REAL mytheta = -atan( ( uGlobal.GetVal(1+6,0)-uGlobal.GetVal(1,0) ) /  LOriginal() );
      if(mytheta > 0.)
        std::cout << "pronto";
    //tamarindo

    TPZManVector<REAL,3> RigidBodyAngles(3);
    this->ComputeRigidBodyMotionRotations(uGlobal, RigidBodyAngles);
    TPZFNMatrix<12> uGlobalPlusRigidBody = uGlobal;
    for(int i = 0; i < 3; i++){
      uGlobalPlusRigidBody(3+i,0) -= RigidBodyAngles[i];
      uGlobalPlusRigidBody(6+3+i,0) -= RigidBodyAngles[i];
    }

    //Transforming the nodal displacements to local ref sys
    //uLocal = R(u) . ( (xb+ub) -(xa+ua) ) - Rorig . (xb-xa)
    TPZFNMatrix<12> InitialPos(12,1,0.), FinalPos(12,1,0.);
    uLocal.Redim(12,1);
    this->GetInitialPosVec(InitialPos);
    FinalPos = InitialPos;
    for(int i = 0; i < 12; i++){
      FinalPos(i,0) += uGlobalPlusRigidBody.GetVal(i,0);//displacement and nodal rotations
    }
    Rot.Multiply(FinalPos, uLocal);
    TPZFNMatrix<144> RotOrig, aux;
    this->RotateMatrix(RotOrig, uGlobal, false);
    RotOrig.Multiply(InitialPos, aux);
    uLocal -= aux;

    twistAngle = this->ComputeTwistRigidBodyAngle(RotOrig, Rot);
    myfileTwist << 180.*twistAngle/M_PI << "\n"; myfileTwist.flush();
    //rotating uLocal over first local axes
    //         R 0 0 0
    //Rtwist = 0 R 0 0
    //         0 0 R 0
    //         0 0 0 R
    twistRot.Redim(12,12);
    REAL cosine = cos(twistAngle);
    REAL sine = sin(twistAngle);
    for(int i = 0; i < 4; i++){
      int pos = 3*i;
      twistRot(pos+0,pos+0) = 1.;
      twistRot(pos+1,pos+1) = cosine;
      twistRot(pos+1,pos+2) = -sine;
      twistRot(pos+2,pos+1) = +sine;
      twistRot(pos+2,pos+2) = cosine;
    }
    aux = uLocal;
    twistRot.Multiply(aux,uLocal);

#ifdef DEBUG
  {
    std::ofstream myfile("c:\\Temp\\rotacao.nb");
    Rot.Print("Rot=",myfile,EMathematicaInput);
    RotOrig.Print("Rorig=",myfile,EMathematicaInput);
  }
#endif
  }
  else{//linear formulation
    uLocal.Redim(12,1);
    ////Multiply R.D = DM. The nodal displacements will be represented in local ref sys (LocalDispl)
    Rot.Multiply(uGlobal,uLocal,0);
    twistAngle = 0.;
    twistRot.Resize(12,12);
    twistRot.Identity();
  }
}///void



