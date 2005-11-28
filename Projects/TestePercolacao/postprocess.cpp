//$Id: postprocess.cpp,v 1.1 2005-11-28 13:49:47 tiago Exp $

#include "meshes.h"

#include "TPZGeoElement.h"

#include "TPZGeoCube.h"
#include "pzshapecube.h"
#include "TPZRefCube.h"
#include "pzshapelinear.h"
#include "TPZGeoLinear.h"
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

#include "pzcompel.h"
#include "TPZCompElDisc.h"
#include "pzintel.h"
#include "pzbndcond.h"
#include "pzpoisson3d.h"
#include "TPZShapeDisc.h"
#include "pzgmesh.h"
#include "pzshapelinear.h"
#include "pzstack.h"
#include <set>
#include "pztrnsform.h"
#include "pzquad.h"

#include "postprocess.h"

using namespace pzshape;

void ComputeNormalFlux(std::ofstream &file, TPZCompMesh &cmesh, const int matid/*, const int npoints*/){
  std::set< std::pair<REAL,REAL> > SolutionSet;
  const int nel = cmesh.NElements();
  for(int i = 0; i < nel; i++){
    TPZCompEl * cel = cmesh.ElementVec()[i];  
    if (!cel) continue;
    if (cel->Material()->Id() != matid) continue;
    TPZCompElSide celside(cel, cel->Reference()->NSides()-1);
    TPZStack<TPZCompElSide> stack;
    celside.EqualLevelElementList(stack, 1, 1);
    if (stack.NElements() != 1){
      std::cout << "ERRO: Nao sei nada - nel = " << stack.NElements() << std::endl;
//       exit(-1);
    }//if
//     else std::cout << "nel = " << stack.NElements() << std::endl;
    
    int n = stack.NElements();
    for(int j = 0; j < n; j++){
        TPZCompEl * neig = stack[j].Element();
        int side = stack[j].Side();
        TPZIntPoints * intrule = neig->Reference()->CreateSideIntegrationRule( side, 10);
        
        TPZManVector<REAL,3> normal(2);
        normal[0] = 0.;
        normal[1] = 1.;
        
        int npt = intrule->NPoints();
        for(int k = 0; k < npt; k++){
          REAL w;
          TPZManVector<REAL,3> point(1), elpoint(2), gradU(2), X(3);
          intrule->Point(k, point, w);
          TPZTransform trans = neig->Reference()->SideToSideTransform(side, neig->Reference()->NSides()-1);
          trans.Apply(point, elpoint);
          ComputeGradient(dynamic_cast<TPZInterpolatedElement*>(neig), elpoint, gradU);
          REAL flux = 0.;
          for(int id = 0; id < normal.NElements(); id++) flux += gradU[id] * normal[id];
          neig->Reference()->X(elpoint, X);
          std::pair<REAL, REAL> mypair(X[0], flux);        
          SolutionSet.insert(mypair);           
        }//k        
    }//j    
  }//i
  
  std::set< std::pair<REAL,REAL> >::iterator w, e;
  e = SolutionSet.end();
  file << "Imprimindo pontos e solucao" << std::endl;
  file << "{" << std::endl;
  w = SolutionSet.begin();
  if (w != e){//first data
    file << "{ " << w->first << " , " << w->second << " }" << std::endl;
    w++;
  }
  for( ; w != e; w++){
    REAL val1 = w->first, 
         val2 = w->second;
    char * number1 = new char[32];
    char * number2 = new char[32];
    sprintf(number1, "%16.16lf", val1);
    sprintf(number2, "%16.16lf", val2);
     file << ",{ " << number1 << " , " << number2 << " }" << std::endl;
//     file << ",{ " << w->first << " , " << w->second << " }" << std::endl;
  }//for w  
  file << "}" << std::endl;
}//method

void ComputeGradient(TPZInterpolatedElement * cel, TPZVec<REAL> &intpoint, TPZVec<REAL> &gradU){
  if (!cel){
    //ERROR
    return;
  }
  int nstate = cel->Material()->NStateVariables();
  int dim = cel->Dimension();
  int nshape = cel->NShapeF();
  
  TPZFNMatrix<9> axes(3,3,0.);
  TPZFNMatrix<9> jacobian(dim,dim);
  TPZFNMatrix<9> jacinv(dim,dim);
  REAL detjac;
  TPZFNMatrix<220> phi(nshape,1);
  TPZFNMatrix<660> dphi(dim,nshape),dphix(dim,nshape);  

  TPZManVector<REAL> sol(nstate,0.);
  TPZFNMatrix<660> dsol(dim,nstate,0.);
  TPZGeoEl * ref = cel->Reference();
  if (!ref){
    //ERROR
    return;
  }

  ref->Jacobian( intpoint, jacobian, axes, detjac , jacinv);
  cel->Shape(intpoint,phi,dphi);
    
  int ieq;
  switch(dim) {
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
  }
  
  int iv=0,d;
  sol.Fill(0.);
  dsol.Zero();
  int ncon = cel->NConnects();
  TPZBlock &block = cel->Mesh()->Block();
  TPZFMatrix &MeshSol = cel->Mesh()->Solution();
  int numdof = nstate;
  for(int in=0; in<ncon; in++) {
    TPZConnect *df = &cel->Connect(in);
    int dfseq = df->SequenceNumber();
    int dfvar = block.Size(dfseq);
    int pos = block.Position(dfseq);
    for(int jn=0; jn<dfvar; jn++) {
      sol[iv%numdof] += phi(iv/numdof,0) * MeshSol(pos+jn,0);
      for(d=0; d<dim; d++)
        dsol(d,iv%numdof) += dphix(d,iv/numdof) * MeshSol(pos+jn,0);
      iv++;
    }
  }
  
  TPZMatPoisson3d * poisson = dynamic_cast<TPZMatPoisson3d*>(cel->Material());
  if (!poisson){
    //ERROR
    return;
  }
  poisson->Solution(sol, dsol, axes, 7, gradU);  
  
}

void ComputeSolution(std::ofstream &file, TPZCompMesh &cmesh, const int matid/*, const int npoints*/){
  std::set< std::pair<REAL,REAL> > SolutionSet;
  const int nel = cmesh.NElements();
  for(int i = 0; i < nel; i++){
    TPZCompEl * cel = cmesh.ElementVec()[i];  
    if (!cel) continue;
    if (!cel->Reference()) continue;
    if (cel->Material()->Id() != matid) continue;
    TPZCompElSide celside(cel, cel->Reference()->NSides()-1);
    TPZStack<TPZCompElSide> stack;
    celside.EqualLevelElementList(stack, 1, 1);
    if (stack.NElements() != 1){
      std::cout << "ERRO: Nao sei nada - nel = " << stack.NElements() << std::endl;
//       exit(-1);
    }//if
//     else std::cout << "nel = " << stack.NElements() << std::endl;    
    int n = stack.NElements();
    for(int j = 0; j < n; j++){
        TPZCompEl * neig = stack[j].Element();
        int side = stack[j].Side();
        TPZIntPoints * intrule = neig->Reference()->CreateSideIntegrationRule( side, 10);
         
        int npt = intrule->NPoints();
        for(int k = 0; k < npt; k++){
          REAL w;
          TPZManVector<REAL,3> point(1), elpoint(2), X(3);
          REAL sol;
          intrule->Point(k, point, w);
          TPZTransform trans = neig->Reference()->SideToSideTransform(side, neig->Reference()->NSides()-1);
          trans.Apply(point, elpoint);
          ComputePointSolution(dynamic_cast<TPZInterpolatedElement*>(neig), elpoint, sol);  
          neig->Reference()->X(elpoint, X);
          std::pair<REAL, REAL> mypair(X[0], sol);        
          SolutionSet.insert(mypair);           
        }//k        
    }//j    
  }//i
  
  std::set< std::pair<REAL,REAL> >::iterator w, e;
  e = SolutionSet.end();
  file << "Imprimindo pontos e solucao" << std::endl;
  file << "{" << std::endl;
  w = SolutionSet.begin();
  if (w != e){//first data
    file << "{ " << w->first << " , " << w->second << " }" << std::endl;
    w++;
  }
  for( ; w != e; w++){
    REAL val1 = w->first,
    val2 = w->second;
    char * number1 = new char[32];
    char * number2 = new char[32];
    sprintf(number1, "%16.16lf", val1);
    sprintf(number2, "%16.16lf", val2);
    file << ",{ " << number1 << " , " << number2 << " }" << std::endl;    
//     file << ",{ " << w->first << " , " << w->second << " }" << std::endl;
  }//for w  
  file << "}" << std::endl;
}//method


void ComputePointSolution(TPZInterpolatedElement * cel, TPZVec<REAL> &intpoint, REAL &solution){
  if (!cel){
    //ERROR
    return;
  }
  int nstate = cel->Material()->NStateVariables();
  int dim = cel->Dimension();
  int nshape = cel->NShapeF();
  
  TPZFNMatrix<9> axes(3,3,0.);
  TPZFNMatrix<9> jacobian(dim,dim);
  TPZFNMatrix<9> jacinv(dim,dim);
  REAL detjac;
  TPZFNMatrix<220> phi(nshape,1);
  TPZFNMatrix<660> dphi(dim,nshape),dphix(dim,nshape);  

  TPZManVector<REAL> sol(nstate,0.);
  TPZFNMatrix<660> dsol(dim,nstate,0.);
  TPZGeoEl * ref = cel->Reference();
  if (!ref){
    //ERROR
    return;
  }

  ref->Jacobian( intpoint, jacobian, axes, detjac , jacinv);
  cel->Shape(intpoint,phi,dphi);
    
  int ieq;
  switch(dim) {
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
  }
  
  int iv=0,d;
  sol.Fill(0.);
  dsol.Zero();
  int ncon = cel->NConnects();
  TPZBlock &block = cel->Mesh()->Block();
  TPZFMatrix &MeshSol = cel->Mesh()->Solution();
  int numdof = nstate;
  for(int in=0; in<ncon; in++) {
    TPZConnect *df = &cel->Connect(in);
    int dfseq = df->SequenceNumber();
    int dfvar = block.Size(dfseq);
    int pos = block.Position(dfseq);
    for(int jn=0; jn<dfvar; jn++) {
      sol[iv%numdof] += phi(iv/numdof,0) * MeshSol(pos+jn,0);
      for(d=0; d<dim; d++)
        dsol(d,iv%numdof) += dphix(d,iv/numdof) * MeshSol(pos+jn,0);
      iv++;
    }
  }
  
  TPZMatPoisson3d * poisson = dynamic_cast<TPZMatPoisson3d*>(cel->Material());
  if (!poisson){
    //ERROR
    return;
  }
  
  TPZManVector<REAL,1> SolutionVec(1);
  poisson->Solution(sol, dsol, axes, 1, SolutionVec);
  solution = SolutionVec[0];
  
}











