
#include "TPZPlacaOrthotropic.h"
#include "TPZMulticamadaOrtho.h"
#include "pzmatorthotropic.h"
#include "pzgmesh.h"
#include "pzintel.h"
#include "pzgeoel.h"
#include "pzbctension.h"
#include "pztempmat.h"
#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "pzskylstrmatrix.h"

TPZMulticamadaOrthotropic::TPZMulticamadaOrthotropic(REAL z,REAL dx,REAL dy, int nelx, int nely) : fDirx(3,0.), fDiry(3,0.) {

  fGeoMesh = new TPZGeoMesh();
  fCompMesh = new TPZCompMesh(fGeoMesh);
  fZMin  =  z;
  fZMax = z;
  fNelx = nelx;
  if(! fNelx%2) fNelx++;//se for par fica impar
  fNely = nely;
  if(! fNely%2) fNely++;
  fDx = dx/fNelx;
  fDy = dy/fNely;
  fGeoMesh->NodeVec().Resize((fNelx+1)*(fNely+1));
  int ix, iy;
  TPZManVector<REAL,3> coord(3,fZMax);
  for(ix=0; ix<= nelx; ix++) {
    for(iy=0; iy<= nely; iy++) {
      coord[0] = ix*fDx-dx/2.;
      coord[1] = iy*fDy-dy/2.;
      fGeoMesh->NodeVec()[ix+iy*(nelx+1)].Initialize(coord,*fGeoMesh);
    }
  }
  fLinearX = 0;
  fLinearY = 0;

  fDirx[0] = 0.25;
  fDiry[0] = 0.25;

  int i;
  for(i=0; i<3; i++) {
    fMX[i] = 0.;
    fMY[i] = 0.;
    fMXY[i] = 0.;
    fNX[i] = 0.;
    fNY[i] = 0.;
    fNXY[i] = 0.;
    fQX[i] = 0.;
    fQY[i] = 0.;
    fdMXdX[i]=0;
    fdMXdY[i]=0;
    fdMYdY[i]=0;
    fdMYdX[i]=0;
    fdMXYdX[i] = 0.;
    fdMXYdY[i] = 0;
    fdQXdX[i] = 0;
    fdQXdY[i] = 0;
    fdQYdX[i] = 0.;
    fdQYdY[i] = 0.;
    fdNXdX[i] = 0.;
    fdNXdY[i] = 0.;
    fdNYdX[i] = 0.;
    fdNYdY[i] = 0.;
    fdNXYdX[i] = 0.;
    fdNXYdY[i] = 0.;
  }
}


void TPZMulticamadaOrthotropic::GenerateMesh(){

  fGeoMesh->BuildConnectivity();
  TPZCompEl::gOrder = 3;
  fCompMesh->AutoBuild();
  int nplaca = fPlacaOrth.NElements();
  int ip;
  for(ip = 0; ip<nplaca; ip++) {
    fPlacaOrth[ip].IdentifyCompEl();
  }
}


void TPZMulticamadaOrthotropic::AddPlacaOrtho(TPZMatOrthotropic *material, REAL height){

  fCompMesh->InsertMaterialObject(material);
  TPZFNMatrix<9> val1(3,3,0.),val2(3,1,0.);
  TPZBCTension *bc = new TPZBCTension(material,-material->Id()*4,4,val1,val2, 1., this,fPlacaOrth.NElements());
  fCompMesh->InsertMaterialObject(bc);
  bc = new TPZBCTension(material,-material->Id()*4-1,4,val1,val2, 1., this,fPlacaOrth.NElements());
  fCompMesh->InsertMaterialObject(bc);
  bc = new TPZBCTension(material,-material->Id()*4-2,4,val1,val2, -1., this,fPlacaOrth.NElements());
  fCompMesh->InsertMaterialObject(bc);
  bc = new TPZBCTension(material,-material->Id()*4-3,4,val1,val2, -1.,this,fPlacaOrth.NElements());
  fCompMesh->InsertMaterialObject(bc);

  //  fPlacaOrth.Push(placa);
  int nnodes = fGeoMesh->NodeVec().NElements();
  fGeoMesh->NodeVec().Resize(nnodes+(fNelx+1)*(fNely+1));
  int ix, iy;
  TPZManVector<REAL,3> coord(3,fZMax+height);
  fZMax += height;
  for(ix=0; ix<= fNelx; ix++) {
    for(iy=0; iy<= fNely; iy++) {
      coord[0] = ix*fDx-fDx*fNelx/2.;
      coord[1] = iy*fDy-fDy*fNely/2.;
      fGeoMesh->NodeVec()[nnodes+ix+iy*(fNelx+1)].Initialize(coord,*fGeoMesh);
    }
  }
  int nodebase1 = nnodes - (fNelx+1)*(fNely+1);
  int elx, ely;
  TPZManVector<int,8> nodeindexes(8,-1);
  for(elx=0; elx<fNelx; elx++) {
    for(ely=0; ely<fNely; ely++) {
      nodeindexes[0] = nodebase1+elx+ely*(fNelx+1);
      nodeindexes[1] = nodebase1+elx+1+ely*(fNelx+1);
      nodeindexes[2] = nodebase1+elx+1+(ely+1)*(fNelx+1);
      nodeindexes[3] = nodebase1+elx+(ely+1)*(fNelx+1);
      int i;
      for(i=0; i<4; i++) nodeindexes[i+4] = nodeindexes[i]+(fNelx+1)*(fNely+1);
      int index;
      TPZGeoEl *gel = fGeoMesh->CreateGeoElement (ECube, nodeindexes, material->Id(), index);
      if(ely == 0) TPZGeoElBC gbc1(gel,21,-material->Id()*4,*fGeoMesh);
      if(elx == fNelx-1) TPZGeoElBC gbc2(gel,22,-material->Id()*4-1,*fGeoMesh);
      if(ely == fNely-1) TPZGeoElBC gbc3(gel,23,-material->Id()*4-2,*fGeoMesh);
      if(elx == 0) TPZGeoElBC gbc4(gel,24,-material->Id()*4-3,*fGeoMesh);
      int nplaca = fPlacaOrth.NElements();
      if(nplaca == 0 && elx == 0 && ely == 0) {
	TPZBndCond *bc2 = new TPZBndCond(material,-100,0,val1,val2);
	fCompMesh->InsertMaterialObject(bc2);
	TPZGeoElBC gbc5(gel,0,-100,*fGeoMesh);
	val1(0,0) = 1.e12;
	val1(2,2) = 1.e12;
	TPZBndCond *bc3 = new TPZBndCond(material,-101,2,val1,val2);
	fCompMesh->InsertMaterialObject(bc3);
	TPZGeoElBC gbc6(gel,3,-101,*fGeoMesh);
	val1.Zero();
	val1(0,0) = 0.;
	val1(2,2) = 1.e12;
	TPZBndCond *bc4 = new TPZBndCond(material,-102,2,val1,val2);
	fCompMesh->InsertMaterialObject(bc4);
	TPZGeoElBC gbc7(gel,2,-102,*fGeoMesh);
      }
      if(elx == fNelx/2 && ely == fNely/2) {
	TPZPlacaOrthotropic placa(gel,fZMax-height,fZMax);
	fPlacaOrth.Push(placa);
      }
    }
  }
}

void  TPZMulticamadaOrthotropic::Print(ostream &out){

  int i, nplaca=fPlacaOrth.NElements();
  out << "TPZMulticamadaOrthotropic::Print\n";
  out << nplaca << endl;
  for (i=0; i<nplaca; i++){
    out << "placa : " << i << endl;
    fPlacaOrth[i].Print(); 
  }
}
  
REAL TPZMulticamadaOrthotropic::Height(){

  return (fZMax - fZMin);
}

int TPZMulticamadaOrthotropic::NPlacas(){
  return fPlacaOrth.NElements();
}

void TPZMulticamadaOrthotropic::AnalyticTensor(TPZVec<REAL> &co, TPZFMatrix &tensor) {

  REAL z = co[2];
  REAL x = co[0];
  REAL y = co[1];
  tensor.Redim(3,3);
  REAL height = Height();
  REAL zrel = z-(fZMax-fZMin)/2.;
  if(height <= 0.) return;
  REAL height3 = height*height*height;
  tensor(0,0) += (fNX[2]+fLinearX*fdNXdX[2]*x)/height;
  tensor(1,1) += (fNY[2]+fLinearY*fdNYdY[2]*y)/height;
  tensor(1,0) += (fNXY[2]+fLinearX*fdNXYdX[2]*x+fLinearY*fdNXYdY[2]*y)/height;
  tensor(0,1) = tensor(1,0);
  tensor(0,0) += 12.*(fMX[2]+x*fLinearX*fdMXdX[2])*zrel/(height3);
  tensor(1,1) += 12.*(fMY[2]+y*fLinearY*fdMYdY[2])*zrel/height3;
  tensor(0,1) += 12.*(fMXY[2]+fLinearX*x*fdMXYdX[2]+fLinearY*y*fdMXYdY[2])*zrel/height3;
  tensor(1,0) += tensor(1,0);
  tensor(0,2) += -6.*(fQX[2]+fLinearX*x*fdQXdX[2])*(zrel*zrel-height*height/4.)/height3;
  tensor(2,0) += tensor(0,2);
  tensor(1,2) += -6.*(fQY[2]+fLinearY*y*fdQYdY[2])*(zrel*zrel-height*height/4.)/height3;
  tensor(2,1) += tensor(1,2);

}

void TPZMulticamadaOrthotropic::ComputeCenterForces() {

  int nplaca = fPlacaOrth.NElements();
  int ip;
  fMX[1] = 0.;
  fMY[1] = 0.;
  fMXY[1] = 0.;
  fNX[1] = 0.;
  fNY[1] = 0.;
  fNXY[1] = 0.;
  fQX[1] = 0.;
  fQY[1] = 0.;
  REAL zref = (fZMax+fZMin)/2.;
  TPZManVector<REAL,3> normal(3,0.),direction(3,0.);
  for(ip=0; ip<nplaca; ip++) {
    normal.Fill(0.);
    direction.Fill(0.);
    normal[0] = 1.;
    direction[0] = 1.;
    fMX[1] += fPlacaOrth[ip].Moment(zref,normal,direction);
    fNX[1] += fPlacaOrth[ip].Force(normal,direction);
    direction.Fill(0.);
    direction[1] = 1.;
    fMXY[1] += fPlacaOrth[ip].Moment(zref,normal,direction);
    fNXY[1] += fPlacaOrth[ip].Force(normal,direction);
    direction.Fill(0.);
    direction[2] = 1.;
    fQX[1] += fPlacaOrth[ip].Force(normal,direction);
    normal.Fill(0.);
    normal[1] = 1.;
    direction.Fill(0.);
    direction[1] = 1.;
    fMY[1] += fPlacaOrth[ip].Moment(zref,normal,direction);
    fNY[1] += fPlacaOrth[ip].Force(normal,direction);
    direction.Fill(0.);
    direction[2] = 1.;
    fQY[1] += fPlacaOrth[ip].Force(normal,direction);
    if(fLinearX) {
      normal.Fill(0.);
      direction.Fill(0.);
      normal[0] = 1.;
      direction[0] = 1.;
      fdMXdX[1] += fPlacaOrth[ip].GradMoment(zref,fDirx,normal,direction);
      fdNXdX[1] += fPlacaOrth[ip].GradForce(fDirx,normal,direction);
      direction.Fill(0.);
      direction[1] = 1.;
      fdMXYdX[1] += fPlacaOrth[ip].GradMoment(zref,fDirx,normal,direction);
      fdNXYdX[1] += fPlacaOrth[ip].GradForce(fDirx,normal,direction);
      direction.Fill(0.);
      direction[2] = 1.;
      fdQXdX[1] += fPlacaOrth[ip].GradForce(fDirx,normal,direction);
      normal.Fill(0.);
      normal[1] = 1.;
      direction.Fill(0.);
      direction[1] = 1.;
      fdMYdX[1] += fPlacaOrth[ip].GradMoment(zref,fDirx,normal,direction);
      fdNYdX[1] += fPlacaOrth[ip].GradForce(fDirx,normal,direction);
      direction.Fill(0.);
      direction[2] = 1.;
      fdQYdX[1] += fPlacaOrth[ip].GradForce(fDirx,normal,direction);
    } else {
      fdMXdX[1] = 0.;
      fdMYdX[1] = 0.;
      fdMXYdX[1] = 0.;
      fdNXdX[1] = 0.;
      fdNYdX[1] = 0.;
      fdNXYdX[1] = 0.;
      fdQXdX[1] = 0.;
      fdQYdX[1] = 0.;
    }
    if(fLinearY) {
      normal.Fill(0.);
      direction.Fill(0.);
      normal[0] = 1.;
      direction[0] = 1.;
      fdMXdY[1] += fPlacaOrth[ip].GradMoment(zref,fDiry,normal,direction);
      fdNXdY[1] += fPlacaOrth[ip].GradForce(fDiry,normal,direction);
      direction.Fill(0.);
      direction[1] = 1.;
      fdMXYdY[1] += fPlacaOrth[ip].GradMoment(zref,fDiry,normal,direction);
      fdNXYdY[1] += fPlacaOrth[ip].GradForce(fDiry,normal,direction);
      direction.Fill(0.);
      direction[2] = 1.;
      fdQXdY[1] += fPlacaOrth[ip].GradForce(fDiry,normal,direction);
      normal.Fill(0.);
      normal[1] = 1.;
      direction.Fill(0.);
      direction[1] = 1.;
      fdMYdY[1] += fPlacaOrth[ip].GradMoment(zref,fDiry,normal,direction);
      fdNYdY[1] += fPlacaOrth[ip].GradForce(fDiry,normal,direction);
      direction.Fill(0.);
      direction[2] = 1.;
      fdQYdY[1] += fPlacaOrth[ip].GradForce(fDiry,normal,direction);
    } else {
      fdMXdY[1] = 0.;
      fdMYdY[1] = 0.;
      fdMXYdY[1] = 0.;
      fdNXdY[1] = 0.;
      fdNYdY[1] = 0.;
      fdNXYdY[1] = 0.;
      fdQXdY[1] = 0.;
      fdQYdY[1] = 0.;
    }
  }
  fMX[2] = fMX[1]-fMX[0];
  fMY[2] = fMY[1]-fMY[0];
  fMXY[2] = fMXY[1]-fMXY[0];
  fNX[2] = fNX[1]-fNX[0];
  fNY[2] = fNY[1]-fNY[0];
  fNXY[2] = fNXY[1]-fNXY[0];
  fQX[2] = fQX[1]-fQX[0];
  fQY[2] = fQY[1]-fQY[0];

  fdMXdX[2] = fdMXdX[1]-fdMXdX[0];
  //  fdMYdX[2] = fdMYdX[1]-fdMYdX[0];
  fdMYdX[2] = 0.;
  fdMXYdX[2] = fdMXYdX[1]-fdMXYdX[0];
  fdNXdX[2] = fdNXdX[1]-fdNXdX[0];
  //  fdNYdX[2] = fdNYdX[1]-fdNYdX[0];
  fdNYdX[2] = 0.;
  fdNXYdX[2] = fdNXYdX[1]-fdNXYdX[0];
  fdQXdX[2] = fdQXdX[1]-fdQXdX[0];
  //  fdQYdX[2] = fdQYdX[1]-fdQYdX[0];
  fdQYdX[2] = 0.;

  //  fdMXdY[2] = fdMXdY[1]-fdMXdY[0];
  fdMXdY[2] = 0.;
  fdMYdY[2] = fdMYdY[1]-fdMYdY[0];
  fdMXYdY[2] = fdMXYdY[1]-fdMXYdY[0];
  // fdNXdY[2] = fdNXdY[1]-fdNXdY[0];
  fdNXdY[2] = 0.;
  fdNYdY[2] = fdNYdY[1]-fdNYdY[0];
  fdNXYdY[2] = fdNXYdY[1]-fdNXYdY[0];
  //  fdQXdY[2] = fdQXdY[1]-fdQXdY[0];
  fdQXdY[2] = 0.;
  fdQYdY[2] = fdQYdY[1]-fdQYdY[0];

}

void TPZMulticamadaOrthotropic::Tensor(TPZVec<REAL> &x, int placa, TPZFMatrix &tensor) {

  REAL zmin = fPlacaOrth[placa].ZMin();
  REAL zmax = fPlacaOrth[placa].ZMax();
  TPZManVector<REAL,3> ksi(3,0.);
  ksi[2] = -1.+2.*(x[2]-zmin)/(zmax-zmin);
  TPZFNMatrix<9> tensorcomp(3,3),tensoranalytic(3,3),gradtensorx(3,3),gradtensory(3,3);
  fPlacaOrth[placa].Tensor(ksi,tensorcomp);
  AnalyticTensor(x,tensoranalytic);
  tensor = tensorcomp+tensoranalytic;
  if(fLinearX) {
    fPlacaOrth[placa].GradTensor(fDirx,ksi,gradtensorx);
    gradtensorx *= x[0];
    tensor += gradtensorx;
  }
  if(fLinearY) {
    fPlacaOrth[placa].GradTensor(fDiry,ksi,gradtensory);
    gradtensory *= x[1];
    tensor += gradtensory;
  }

}

void TPZMulticamadaOrthotropic::ComputeSolution(ostream &out,int print){

  TPZAnalysis an(fCompMesh);
  TPZSkylineStructMatrix skyl(fCompMesh);
  an.SetStructuralMatrix(skyl);
  TPZStepSolver solve;
  solve.SetDirect(ELDLt);
  an.SetSolver(solve);
  an.Run();
  if(print) an.Print("* PRINT ANALISYS *",out);  
}

void TPZMulticamadaOrthotropic::PrintTensors(ostream &out) {

  int nplaca = fPlacaOrth.NElements();
  int ip;
  for(ip=0; ip< nplaca; ip++) {
    out << "Tensor para placa " << ip << endl;
    fPlacaOrth[ip].PrintTensors(out);
  }
}

void TPZMulticamadaOrthotropic::PrintCenterForces(ostream &out) {

  int i;
  out << "fMX ";
  for(i=0; i<3; i++) {
    out << 
    fMX[i] << " ";
  }
  out << "\nfMY ";
  for(i=0; i<3; i++) {
    out << 
    fMY[i] << " ";
  }
  out << "\nfMXY ";
  for(i=0; i<3; i++) {
    out << 
    fMXY[i] << " ";
  }
  out << "\nfNX ";
  for(i=0; i<3; i++) {
    out << 
    fNX[i] << " ";
  }
  out << "\nfNY ";
  for(i=0; i<3; i++) {
    out << 
    fNY[i] << " ";
  }
  out << "\nfNXY ";
  for(i=0; i<3; i++) {
    out << 
    fNXY[i] << " ";
  }
  out << "\nfQX ";
  for(i=0; i<3; i++) {
    out << 
    fQX[i] << " ";
  }
  out << "\nfQY ";
  for(i=0; i<3; i++) {
    out << 
    fQY[i] << " ";
  }
  out << "\nfdMXdX ";
  for(i=0; i<3; i++) {
    out << 
    fdMXdX[i] << " ";
  }
  out << "\nfdMXdY ";
  for(i=0; i<3; i++) {
    out << 
    fdMXdY[i] << " ";
  }
  out << "\nfdMYdY ";
  for(i=0; i<3; i++) {
    out << 
    fdMYdY[i] << " ";
  }
  out << "\nfMYdX ";
  for(i=0; i<3; i++) {
    out << 
    fdMYdX[i] << " ";
  }
  out << "\nfdMXYdX ";
  for(i=0; i<3; i++) {
    out << 
    fdMXYdX[i] << " ";
  }
  out << "\nfdMXYdY ";
  for(i=0; i<3; i++) {
    out << 
    fdMXYdY[i] << " ";
  }
  out << "\nfdQXdX ";
  for(i=0; i<3; i++) {
    out << 
    fdQXdX[i] << " ";
  }
  out << "\nfdQXdY ";
  for(i=0; i<3; i++) {
    out << 
    fdQXdY[i] << " ";
  }
  out << "\nfdQYdX ";
  for(i=0; i<3; i++) {
    out << 
    fdQYdX[i] << " ";
  }
  out << "\nfdQYdY ";
  for(i=0; i<3; i++) {
    out << 
    fdQYdY[i] << " ";
  }
  out << "\nfdNXdX ";
  for(i=0; i<3; i++) {
    out << 
    fdNXdX[i] << " ";
  }
  out << "\nfdNXdY ";
  for(i=0; i<3; i++) {
    out << 
    fdNXdY[i] << " ";
  }
  out << "\nfdNYdX ";
  for(i=0; i<3; i++) {
    out << 
    fdNYdX[i] << " ";
  }
  out << "\nfdNYdY ";
  for(i=0; i<3; i++) {
    out << 
    fdNYdY[i] << " ";
  }
  out << "\nfdNXYdX ";
  for(i=0; i<3; i++) {
    out << 
    fdNXYdX[i] << " ";
   }
  out << "\nfdNXYdY ";
  for(i=0; i<3; i++) {
    out << 
   fdNXYdY[i] << " ";
  }
  out << endl;
}
