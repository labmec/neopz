/***************************************************************************
                          tmbxubinshih.cpp  -  description
                             -------------------
    begin                : Thu Oct 23 2003
    copyright            : (C) 2003 by cesar
    email                : cesar@becks.fec.unicamp.br
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "tmbxubinshih.h"
#include <pzcompel.h>
#include <pzgeoel.h>
//REAL TMBXubinShih::gcos45 = 0.70710678118655;

TMBXubinShih::TMBXubinShih(TPZCompMesh *cmesh, int nstate, TPZVec<int> &dimstate,
                            TPZVec<int> &statetoanalyse,TPZVec<REAL> &sol, int source) :
                            TPZErrorIndicator (nstate,dimstate,statetoanalyse,sol,cmesh){


   cos45 = sqrt(2.) / 2.;
   fSource = source;
   //nothing to do here...
   
                              
}

TMBXubinShih::~TMBXubinShih(){

}

/** @see TPZErrorIndicator class documentation */
void TMBXubinShih::MarkedElListH(TPZVec<int> &elindex, TPZVec<int> &side,int sidedim, int sidestate){

  int e,s;
  int nstate = fState.NElements();
  TPZFMatrix results(fNElements,nstate,0.);
  int nelem = fMesh->ElementVec().NElements();
  for (e=0;e<nelem;e++){
    TPZCompEl *cel = fMesh->ElementVec()[e];
    if (!cel) continue;
    TPZVec<REAL> elresults(nstate,0.);
    if (fSource) ElSourceIndicator(cel,fState,elresults);
    else ElLocationIndicator(cel,fState,elresults);
    for (s=0;s<nstate;s++) {
      results(e,s) = elresults[s];
      if (!fErAnType[s]){
        if (results(e,s) > fMaxError[s]) elindex[e] = 1;
        if (results(e,s) < fMinError[s]) elindex[e] = -1;
      }
    }
  }
  
  bool cont = false;
  for (s=0;s<fState.NElements();s++){
    if (fErAnType[s]){
      cont = true;
      break;
    }
  }

  if(cont){
    TPZFMatrix perm (fNElements,fState.NElements(),0.);
    Sort(results,perm);
    //cout << "Error " << results << endl;
    //cout << "Perm  " << perm << endl;
    TPZVec <REAL> sum (fState.NElements(),0.);
    TPZVec <REAL> higher (fState.NElements(),0.);
    TPZVec <REAL> lower (fState.NElements(),0.);
    for (e=0;e<fNElements;e++){
      for (s=0;s<fState.NElements();s++){
        higher[s] += results(e,s);
      }
    }
    for (s=0;s<fState.NElements();s++){
      lower[s] = higher[s] * (1. - fMinError[s]);
      higher[s] *= fMaxError[s];
    }
    //cout << "Higher " << higher << endl;
    //cout << "Lower  " << lower << endl;
    for (e=0;e<fNElements;e++){
      for (s=0;s<fState.NElements();s++){
        if (sum[s] < higher[s]) elindex[(int)perm(e,s)] = 1;
        else if (sum[s] > lower[s] && elindex[e] != 1) elindex[e] = -1;
        sum [s] += results((int)perm(e,s),s);
      }
    }
  }
  if (sidedim > -1 && sidestate > -1 && sidestate < fNState){
    side.Resize(fNElements);
    side.Fill(-1);
    for (e=0;e<fNElements;e++){
      TPZCompEl *cel = fMesh->ElementVec()[e];
      if (!cel || elindex[e] < 1) continue;
      side[e] = GetRefSide(cel,sidedim,sidestate,&results);
    }
  }  
}

/** Evaluates the source error indicator */
void TMBXubinShih::SourceIndicator(TPZVec<int> &statetoanalyse, TPZVec<REAL> &elresults){
//Evaluates the maximum angle between the element and its
//neighbours whose angle formed by the element center point
//and the neighbour center point are, in modulus, greater
//than 45 degrees.

//get the element
  int i,j;
  int nstate = statetoanalyse.NElements();
  TPZFMatrix results(fNElements,nstate,0.);
  int nelem = fMesh->ElementVec().NElements();  
  for (i=0;i<nelem;i++){
    TPZCompEl *cel = fMesh->ElementVec()[i];
    if (!cel) continue;
    TPZVec<REAL> elresults(nstate,0.);
    ElSourceIndicator(cel,statetoanalyse,elresults);
    for (j=0;j<nstate;j++) results(i,j) = elresults[j];
  }

  //sort or no?...

}


void TMBXubinShih::ElSourceIndicator(TPZCompEl *cel, TPZVec<int> &statetoanalyse, TPZVec<REAL> &results){
  int i,j,d;
  //take the element center point
//  TPZVec<REAL> cp (3,0.);
//  TPZVec<REAL> nc (3,0.);
  TPZGeoEl *gel = cel->Reference();
//  gel->CenterPoint(gel->NSides()-1,cp);
/*
 REAL elAngle = 0.;
// We will take all elements with the same dimension of the analysed element
  //take the neighbours whose angle is greater than Pi/4
  TPZStack <TPZGeoEl *> neighbours;
  for (j=0;j<gel->NSides();j++){
    TPZGeoElSide neighside = gel->Neighbour(j);
    if (neighside.Element() == gel) continue;
    TPZGeoEl *neigh = neighside.Element();
    // take the neighbour center point
    neigh->CenterPoint(neigh->NSides()-1,nc);
    elAngle = EvaluateCosAngle(cp,nc);
    //verify the angle between the element and the neighbour
    if (elAngle > -cos45 && elAngle < cos45 && neigh->Reference()) neighbours.Push(neigh);
    while (neigh != gel){
      neighside = neigh->Neighbour(neighside.Side());
      if (neighside.Element() == gel) break;
      neigh = neighside.Element();
      neigh->CenterPoint(neigh->NSides()-1,nc);
      elAngle = EvaluateCosAngle(cp,nc);
      if (elAngle > -cos45 && elAngle < cos45 && neigh->Reference()) neighbours.Push(neigh);
    }
  }
*/
  TPZStack <TPZGeoEl *> neighbours;
  int elemdim = gel->Dimension();
  for (j=0;j<gel->NSides();j++){
    TPZGeoElSide neighside = gel->Neighbour(j);
    if (neighside.Element() == gel) continue;
    if (neighside.Element()->Dimension() != elemdim) continue;
    TPZGeoEl *neigh = neighside.Element();
    neighbours.Push(neigh);
    
    // take the neighbour center point
    //neigh->CenterPoint(neigh->NSides()-1,nc);
    //elAngle = EvaluateCosAngle(cp,nc);
    //verify the angle between the element and the neighbour
    //if (elAngle > -cos45 && elAngle < cos45 && neigh->Reference()) neighbours.Push(neigh);
    while (neigh != gel){
      neighside = neigh->Neighbour(neighside.Side());
      if (neighside.Element() == gel) break;
      neigh = neighside.Element();
      neighbours.Push(neigh);
      //neigh->CenterPoint(neigh->NSides()-1,nc);
      //elAngle = EvaluateCosAngle(cp,nc);
      //if (elAngle > -cos45 && elAngle < cos45 && neigh->Reference()) neighbours.Push(neigh);
    }
  }
  
  int nstate = statetoanalyse.NElements();
  int elindex = cel->Index();
  for (i=0;i<neighbours.NElements();i++){
    TPZGeoEl *neigh = neighbours[i];
    TPZCompEl *cneigh = neigh->Reference();
    if(!cneigh) continue;
    int neighindex = neigh->Reference()->Index();  
    for (j=0;j<nstate;j++){
      int state = statetoanalyse[j];
      int dim = fDim[state];
      TPZVec<REAL> elsol(dim,0.);
      TPZVec<REAL> neigsol(dim,0.);
      for (d=0;d<dim;d++){
        int solelindex = Index(elindex,state,d);
        int neighsolindex = Index(neighindex,state,d);
        //take the element result for the selected state variables
        elsol[d] = fSolution[solelindex];
        //take the neighbour result for the selected state variables
        neigsol[d] = fSolution[neighsolindex];
      }
      //take the angle between the element and neighbour results
      REAL angle = EvaluateCosAngle(elsol,neigsol);
      angle = acos(angle);
      //verifies if the angle is greater than the previous result
      if (angle > results[j]) results[j] = angle;
    }
  }
}

/** Evaluates the Xubin-Shih location error indicator */
void TMBXubinShih::ElLocationIndicator(TPZCompEl *cel, TPZVec<int> &statetoanalyse, TPZVec<REAL> &results){
  int i,j,d;
  // Neighbours norm and delta evaluation
  int elindex = cel->Index();
  TPZGeoEl *gel = cel->Reference();
  TPZStack <TPZGeoEl *> neighbours;
  int elemdim = gel->Dimension();
  for (j=0;j<gel->NSides();j++){
    TPZGeoElSide neighside = gel->Neighbour(j);
    if (neighside.Element() == gel) continue;
    if (neighside.Element()->Dimension() != elemdim) continue;
    TPZGeoEl *neigh = neighside.Element();
    neighbours.Push(neigh);
    while (neigh != gel){
      neighside = neigh->Neighbour(neighside.Side());
      if (neighside.Element() == gel) break;
      neigh = neighside.Element();
      neighbours.Push(neigh);
    }
  }
  
  //This should be changed to neighbours by faces (3D) or ribs (2D)...
/*REAL elAngle = 0.;
  TPZVec<REAL> nc (3,0.);
  TPZVec<REAL> cp (3,0.);

  for (j=0;j<gel->NSides();j++){
    TPZGeoElSide neighside = gel->Neighbour(j);
    if (neighside.Element() == gel) continue;
    TPZGeoEl *neigh = neighside.Element();
    // take the neighbour center point
    neigh->CenterPoint(neigh->NSides()-1,nc);
    elAngle = EvaluateCosAngle(cp,nc);
    //verify the angle between the element and the neighbour
    if (elAngle > -cos45 && elAngle < cos45 && neigh->Reference()) neighbours.Push(neigh);
    while (neigh != gel){
      neighside = neigh->Neighbour(neighside.Side());
      if (neighside.Element() == gel) break;
      neigh = neighside.Element();
      neigh->CenterPoint(neigh->NSides()-1,nc);
      elAngle = EvaluateCosAngle(cp,nc);
      if (elAngle > -cos45 && elAngle < cos45 && neigh->Reference()) neighbours.Push(neigh);
    }
  }*/


  for (i=0;i<neighbours.NElements();i++){
    TPZGeoEl *neigh = neighbours[i];
    TPZCompEl *cneigh = neigh->Reference();
    if (!cneigh) continue;
    int neighindex = cneigh->Index();
    int nstate = statetoanalyse.NElements();
    for (j=0;j<nstate;j++){
      int state = statetoanalyse[j];
      int dim = fDim[state];
      REAL esolnorm = 0.;
      REAL nsolnorm = 0.;
      for (d=0;d<dim;d++){
        int solelindex = Index(elindex,state,d);
        int neighsolindex = Index(neighindex,state,d);
        //take the neighbour result for the selected state variables
        esolnorm += fSolution[solelindex]*fSolution[solelindex];
        nsolnorm += fSolution[neighsolindex]*fSolution[neighsolindex];
      }
      esolnorm = sqrt(esolnorm);
      nsolnorm = sqrt(nsolnorm);
      //take the difference between the element and neighbour results
      REAL delta = fabs(esolnorm - nsolnorm);
      //verifies if the angle is greater than the previous result
      if (delta > results[j]) results[j] = delta;
    }
  }
}

REAL TMBXubinShih::EvaluateCosAngle(TPZVec<REAL> u, TPZVec<REAL> v){
#ifndef NDEBUG
  if (u.NElements() != v.NElements()){
    cout << "TMBXubinShih::EvalutateCosAngle Error vectors have not the same dimension\n";
    return 0.;
  }
#endif
  int i;
  REAL cosin = 0.;
  int dim = u.NElements();
  for (i=0;i<dim;i++){
    cosin += u[i] * v[i];
  }
  cosin = sqrt(cosin);
  cosin /= L2Norm(u);
  cosin /= L2Norm(v);
  return cosin;
}

REAL TMBXubinShih::L2Norm(TPZVec<REAL> &vec){
  int i,size = vec.NElements();
  REAL norm = 0.;
  for (i=0;i<size;i++){
    norm += vec[i]*vec[i];
  }
  norm = sqrt(norm);
  return norm;
}
