// -*- c++ -*-

#include <stdlib.h>
#include "TPZNLMultGridAnalysis.h"
#include "TPZCompElDisc.h"
#include "TPZAgglomerateEl.h"

#include "pzcmesh.h"
#include "pzintel.h"
#include "pzgeoel.h"
#include "pztransfer.h"
#include "pzadmchunk.h"
#include "pzbdstrmatrix.h"
#include "pzblockdiag.h"
#include "pzskylmat.h"
#include "pzskylstrmatrix.h"
#include "TPZFrontSym.h"
#include "TPZFrontNonSym.h"
#include "TPZFrontStructMatrix.h"
#include "pzmgsolver.h"
#include "pzseqsolver.h"
#include "pzstepsolver.h"
#include "pzquad.h"

#include "pzonedref.h"

//class TPZTransfer;


TPZNonLinMultGridAnalysis::TPZNonLinMultGridAnalysis(TPZCompMesh *cmesh) : TPZAnalysis(cmesh) {
  fMeshes.Push(cmesh);
}

TPZNonLinMultGridAnalysis::~TPZNonLinMultGridAnalysis() {
  while (fMeshes.NElements()) delete fMeshes.Pop();
  while(fSolutions.NElements()) delete fSolutions.Pop();
  while(fSolvers.NElements()) delete fSolvers.Pop();
  while(fPrecondition.NElements()) delete fPrecondition.Pop();
}


void TPZNonLinMultGridAnalysis::AppendMesh(TPZCompMesh * mesh){

  if(fMeshes.NElements() != fSolvers.NElements() || fMeshes.NElements() != fSolutions.NElements() ||
     fPrecondition.NElements() != fMeshes.NElements()) {
    cout << "TPZNonLinMultGridAnalysis::AppendMesh can only be called after solving the coarse mesh\n";
    return;
  }
  fMeshes.Push(mesh);
}

TPZCompMesh *TPZNonLinMultGridAnalysis::PopMesh() {

  if(fMeshes.NElements() == 1) {
    cout << "TPZNonLinMultGridAnalysis cannot delete the root mesh, sorry\n";
    return 0;
  }
  if(fSolutions.NElements() == fMeshes.NElements()) delete fSolutions.Pop();
  delete fSolvers.Pop();
  SetSolver(*fSolvers[fSolvers.NElements()-1]);
  fCompMesh = fMeshes[fMeshes.NElements()-2];
  fSolution = *fSolutions[fSolutions.NElements()-1];
  return fMeshes.Pop();
}

void TPZNonLinMultGridAnalysis::Solve() {
//   if(fMeshes.NElements() == 1) {
//     TPZAnalysis::Solve();
//     if(fSolvers.NElements() == 0) {
//       fSolvers.Push((TPZMatrixSolver *) fSolver->Clone());
//     }
//     if(fPrecondition.NElements() == 0) {
//       fPrecondition.Push(0);
//     }
//     if(fSolutions.NElements() == 0) {
//       fSolutions.Push(new TPZFMatrix(fSolution));
//     } else {
//       int nsol = fSolutions.NElements();
//       *(fSolutions[nsol-1]) = fSolution;
//     }
//     return;
//   }
//   int numeq = fCompMesh->NEquations();
//   if(fRhs.Rows() != numeq ) return;
//   int nsolvers = fSolvers.NElements();
  
//   TPZFMatrix residual(fRhs);
//   TPZFMatrix delu(numeq,1);
//   TPZMatrixSolver *solve = dynamic_cast<TPZMatrixSolver *> (fSolvers[nsolvers-1]);
//   if(fSolution.Rows() != numeq) {
//     fSolution.Redim(numeq,1);
//   } else {
//     solve->Matrix()->Residual(fSolution,fRhs,residual);
//   }
  
//   REAL normrhs = Norm(fRhs);
//   REAL normres  = Norm(residual);
//   if(normrhs*1.e-6 >= normres) {
//     cout << "TPZNonLinMultGridAnalysis::Solve no need for iterations normrhs = " << normrhs << " normres = " << normres << endl;
//     if(fSolutions.NElements() < fMeshes.NElements()) {
//       fSolutions.Push(new TPZFMatrix(fSolution));
//     } else {
//       int nsol = fSolutions.NElements();
//       *(fSolutions[nsol-1]) = fSolution;
//     }
//     return ;
//   }
// //   REAL tol = 1.e-6*normrhs/normres;
// //   if(numeq > 1500) {
// //     fIterative->SetCG(200,*fPrecond,tol,0);
// //   } else {
// //     fIterative->SetDirect(ELDLt);
// //   }
//   TPZStepSolver *stepsolve = dynamic_cast<TPZStepSolver *> (solve);
//   if(stepsolve) stepsolve->SetTolerance(1.e-6*normrhs/normres);
//   cout << "TPZNonLinMultGridAnalysis::Run res : " << Norm(residual) << " neq " << numeq << endl;
//   solve->Solve(residual, delu);
//   fSolution += delu;
	
//   fCompMesh->LoadSolution(fSolution);
//   if(fSolutions.NElements() < fMeshes.NElements()) {
//     fSolutions.Push(new TPZFMatrix(fSolution));
//   } else {
//     int nsol = fSolutions.NElements();
//     *(fSolutions[nsol-1]) = fSolution;
//  
}

TPZCompMesh *TPZNonLinMultGridAnalysis::AgglomerateMesh (TPZCompMesh *finemesh,int levelnumbertorefine,int setdegree){

  TPZVec<int> accumlist;
  int numaggl,dim;
  TPZAgglomerateElement::ListOfGroupings(finemesh,accumlist,levelnumbertorefine,numaggl,dim);
  TPZCompMesh *aggmesh = new TPZCompMesh(finemesh->Reference());
  TPZCompElDisc::CreateAgglomerateMesh(finemesh,*aggmesh,accumlist,numaggl);
  return aggmesh;
}

TPZCompMesh  *TPZNonLinMultGridAnalysis::UniformlyRefineMesh(TPZCompMesh *coarcmesh,int levelnumbertorefine,int setdegree) {

  TPZGeoMesh *gmesh = coarcmesh->Reference();
  if(!gmesh) {
    cout << "TPZMGAnalysis::UniformlyRefineMesh mesh with null reference, cancelled method\n";
    return 0;
  }
  cout << "\nTPZNonLinMultGridAnalysis::UniformlyRefineMesh uniforme division of coarcmesh,"
       << " levels to be fine = " << levelnumbertorefine << endl;

  gmesh->ResetReference();
  TPZCompMesh *finemesh = new TPZCompMesh(gmesh);
  int nmat = coarcmesh->MaterialVec().NElements();
  int m;
  for(m=0; m<nmat; m++) {
    TPZMaterial *mat = coarcmesh->MaterialVec()[m];
    if(!mat) continue;
    mat->Clone(finemesh->MaterialVec());
  }
  TPZAdmChunkVector<TPZCompEl *> &elementvec = coarcmesh->ElementVec();
  int el,nelem = elementvec.NElements();
  for(el=0; el<nelem; el++) {
    TPZCompEl *cel = elementvec[el];
    if(!cel) continue;
    if(cel->Type() == EAgglomerate){
      PZError << "TPZNonLinMultGridAnalysis::UniformlyRefineMesh mesh error, mesh not refined\n";
      return new TPZCompMesh(NULL);
    }
    if(cel->Type() != EDiscontinuous) continue;
    TPZCompElDisc *disc = dynamic_cast<TPZCompElDisc *>(cel);
    int degree = disc->Degree();

    TPZGeoEl *gel = disc->Reference();
    if(!gel) {
      cout << "TPZMGAnalysis::UniformlyRefineMesh encountered an element without geometric reference\n";
      continue;
    }
    TPZStack<TPZGeoEl *> sub0,sub1,sub;
    //GetRefinedGeoEls(geo,sub);
    int lev = 0,k,nsons,i;
    gel->Divide(sub0);
    while(lev <  levelnumbertorefine){
      int nsubs = sub0.NElements();
      TPZVec<TPZGeoEl *> copy(sub0);
      for(i=0;i<nsubs;i++){
	copy[i]->Divide(sub1);
	nsons = sub1.NElements();
	if(lev == levelnumbertorefine){
	  for(k=0;k<nsons;k++) sub.Push(sub1[k]);
	} else {
	  for(k=0;k<nsons;k++) sub0.Push(sub1[k]);
	}
      }
      lev++;
    }
    int nsub = sub.NElements(),isub,index;
    //o construtor adequado ja deveria ter sido definido
    for(isub=0; isub<nsub; isub++) {
      disc = dynamic_cast<TPZCompElDisc *>(sub[isub]->CreateCompEl(*finemesh,index));
      if(setdegree != degree) disc->SetDegree(degree);
    }
  }
  return finemesh;
}

void TPZNonLinMultGridAnalysis::GetRefinedGeoEls(TPZGeoEl *geo,TPZStack<TPZGeoEl *> sub,int &levelnumbertorefine){

//   TPZVec<TPZGeoEl *> sub0;
//   int numref = 0,i;
//   if(!levelnumbertorefine){
//     sub.Push(geo);
//     return;
//   }
//   geo->Divide(sub0);
//   levelnumbertorefine--;
}

//   lev = 0;
//   geo->Divide(sub0);
//   do {
//     nsubs = sub0.NElements();
//     TPZVec<TPZGeoEl *> copy(sub0);
//     for(i=0;i<nsubs;i++){
//       copy[i]->Divide(sub1);
//       nsons = sub1.NElements();
//       for(k=0;k<nsons;k++) sub0.Push(sub1[k]);
//       if(lev == levelnumbertorefine) for(k=0;k<nsons;k++) sub.Push(sub1[k]);
//     }
//     lev++;
//   } while(lev <  levelnumbertorefine);


int TPZNonLinMultGridAnalysis::main() {

//   TPZCompMesh *cmesh = TPZCompElQ2d::CreateMesh();
//   cmesh->CleanUpUnconnectedNodes();
//   TPZFMatrix sol(cmesh->NEquations(),1);
//   ofstream out("output.txt");
//   int row = sol.Rows();
//   int r;
//    for(r=0; r<row; r++) {
// //     //    sol(r,0) = rand()/(RAND_MAX+1.);
//      sol(r,0) = 1.;
    
//    }
//   TPZNonLinMultGridAnalysis mgan(cmesh);
//   TPZSkylineStructMatrix strskyl(cmesh);
//   mgan.SetStructuralMatrix(strskyl);
//   TPZStepSolver direct;
//   direct.SetDirect(ELDLt);
//   mgan.SetSolver(direct);
//   mgan.Run();
//   TPZGeoMesh *gmesh = cmesh->Reference();
//   int nel = gmesh->ElementVec().NElements();
//   int el;
//   TPZVec<TPZGeoEl *> sub;
//   for(el=0; el<nel; el++) {
//     TPZGeoEl *gel = gmesh->ElementVec()[el];
//     if(!gel) continue;
//     //     if(!gel->HasSubElement(0)) {
//        gel->Divide(sub);
//        //       break;
//        //     }
//   }
//   gmesh->ResetReference();
//   TPZCompMesh *cmesh2 = new TPZCompMesh(gmesh);
//   //  cmesh2.LoadReferences();
//   TPZAdmChunkVector<TPZMaterial *> &matvec2 = cmesh2->MaterialVec();
//   TPZAdmChunkVector<TPZMaterial *> &matvec = cmesh->MaterialVec();
//   int nmat = matvec.NElements();
//   TPZMaterial *mat;
//   int im;
//   for(im=0; im<nmat; im++) {
//     mat = matvec[im];
//     if(!mat) continue;
//     mat->Clone(matvec2);
//   }
// //   TPZFMatrix xk(1,1,1.),xc(1,1,1.),xf(1,1,1.);
// //    //   xk(0,1) = xk(1,0) = xc(0,1) = xc(1,0) = 0.;
// //   mat->SetMaterial(xk,xc,xf);
// //   TPZFMatrix val1(1,1,0.),val2(1,1,0.);
// //   TPZBndCond *bc = mat->CreateBC(-1,0,val1,val2);
// //   cmesh2->InsertMaterialObject(mat);
// //   cmesh2->InsertMaterialObject(bc);

//   //  gmesh->Print();
//   cmesh2->AutoBuild();
//   cmesh2->AdjustBoundaryElements();
//   cmesh2->CleanUpUnconnectedNodes();
//   mgan.AppendMesh(cmesh2);

//   mgan.Run();
//   // TPZNonLinMultGridAnalysis an(cmesh2);
//   //  TPZTransfer *trf = TPZNonLinMultGridAnalysis::BuildTransferMatrix(cmesh2,cmesh);
//   cmesh->LoadSolution(sol);
//   TPZTransfer trf;
//   cmesh2->BuildTransferMatrix(*cmesh,trf);
//   trf.Print("Transfer Matrix",out);
//   TPZFMatrix sol2(cmesh2->NEquations(),1,0.);
//   trf.TransferSolution(sol,sol2);
//   cmesh2->LoadSolution(sol2);
//   gmesh->Print(out);
//   cmesh->Print(out);
//   cmesh2->Print(out);
//   cmesh->Solution().Print("Coarse mesh solution",out);
//   cmesh2->Solution().Print("Fine mesh solution",out);

//   TPZVec<REAL> ervec,truervec;
//   TPZNonLinMultGridAnalysis::MeshError(cmesh2,cmesh,ervec,mgan.fExact,truervec);
//   int i;
//   cout << "TPZNonLinMultGridAnalysis the error between both meshes is equal to \n";
//   for(i=0; i<ervec.NElements(); i++) cout << ervec[i] << ' ';
//   cout << endl;
  return 1;
}
