#include "tmbheidi.h"

TMBHeidi::TMBHeidi(int nstate, TPZVec<int> &dimstate, TPZVec<int> &statetoanalyse, TPZVec<REAL> &sol) : TPZErrorIndicator(nstate,dimstate,statetoanalyse,sol,0){
}

TMBHeidi::~TMBHeidi(){
}

void TMBHeidi::MarkedElListH(TPZVec<int> &elindex, TPZVec<int> &side,int sidedim, int sidestate){
  int e,s,d;
  int state,index;
  bool cont = false;
  TPZVec <REAL> maxvec (fState.NElements(),0.);
  TPZVec <REAL> minvec (fState.NElements(),0.);
  TPZVec <REAL> delta_sol (fState.NElements(),0.);
  TPZVec <REAL> sum (fState.NElements(),0.);
  TPZVec <REAL> higher (fState.NElements(),0.);
  TPZVec <REAL> lower (fState.NElements(),0.);
  //matrix of error -->> rows = elements ; cols = state variables
  TPZFMatrix error (fNElements,fState.NElements(),0.);
  TPZFMatrix perm (fNElements,fState.NElements(),0.);

  FindMaxMin(maxvec,minvec,delta_sol);
  //cout << "fSolution " << fSolution << endl;
  //cout << "Maxvec " << maxvec << endl;
  //cout << "Minvec " << minvec << endl;
  //cout << "Delta " << delta_sol << endl;
  for (e=0;e<fNElements;e++){
    for (s=0;s<fState.NElements();s++){
      state = fState[s];
      for (d=0;d<fDim[state];d++){
        index = Index(e,state,d);
        REAL sol = fSolution[index];
        error(e,s) +=  sol * sol;
      }
      //cout << "( e , s ) = ( " << e << " , " << s << " )" << endl;
      //cout << "Error squared " << error (e,s) << endl;
      error(e,s) = sqrt (error(e,s));
      //cout << "Root  " << error(e,s) << endl;;
      error(e,s) = error(e,s) / delta_sol[s];
      //cout << "By delta " << error(e,s) << endl;
      //cout << error << endl;
      if (!fErAnType[s]){
        if (error(e,s) > fMaxError[s]) elindex[e] = 1;
        if (error(e,s) < fMinError[s]) elindex[e] = -1;
      }
    }
  }
  //cout << "Error " << error << endl;
  //cout << "elindex - 0 " << elindex << endl;
  for (s=0;s<fState.NElements();s++){
    if (fErAnType[s]){
      cont = true;
      break;
    }
  }
  if(cont){
    Sort(error,perm);
    //cout << "Error " << error << endl;
    //cout << "Perm  " << perm << endl;
    for (e=0;e<fNElements;e++){
      for (s=0;s<fState.NElements();s++){
        higher[s] += error(e,s);
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
        sum [s] += error((int)perm(e,s),s);
      }
    }
  }
  if (sidedim > -1 && sidestate > -1 && sidestate < fNState){
    side.Resize(fNElements);
    side.Fill(-1);
    for (e=0;e<fNElements;e++){
      TPZCompEl *cel = fMesh->ElementVec()[e];
      if (!cel || elindex[e] < 1) continue;
      side[e] = GetRefSide(cel,sidedim,sidestate,&error);
    }
  }
}

void TMBHeidi::FindMaxMin(TPZVec<REAL> &max, TPZVec<REAL> &min, TPZVec<REAL> &delta){
  int e,s,d,index,state ;
  REAL sol = 0.;
  for (s=0;s<fState.NElements();s++){
    max[s] = -1e300;
    min[s] = +1e300;
  }
  for (e=0;e<fNElements;e++){
    for(s=0;s<fState.NElements();s++){
      state = fState[s];
      sol = 0.;
      for(d=0;d<fDim[state];d++){
				index = Index(e,state,d);
				REAL aux = fSolution[index] ;
				sol += aux * aux;
      }
      sol = sqrt (sol);
      if (sol > max[s]) max[s] = sol;
      if (sol < min[s]) min[s] = sol;
    }
  }
  delta.Fill(0.);
  for (s=0;s<fState.NElements();s++){
    if (fabs(max[s]-min[s]) < 1e-12){
      cout << "TMBHeidi::FindMaxMin : ERROR : constant solution detected. Aborting !"<< endl;
      cout.flush();
      exit(-1);
    }
    delta[s] = max[s] - min[s];
  }
}
