#include "pzerror_ind.h"
#include "pzmaterial.h"
#include "pzcompel.h"
#include "pzgeoel.h"
#include "pzmatrix.h"

TPZErrorIndicator::TPZErrorIndicator(int nstate, TPZVec<int>& dimstate, TPZVec<int>& statetoanalyse,  TPZVec<REAL> &sol, TPZCompMesh *cmesh ){
  if (cmesh) SetMesh(cmesh);
  else fMesh = 0;
  SetSolution (sol,nstate,dimstate,statetoanalyse);
}

TPZErrorIndicator::~TPZErrorIndicator (){
  
}

void TPZErrorIndicator::SetMesh(TPZCompMesh *cmesh){
  if (!cmesh){
    cout << "TPZErrorIndicator::SetMesh - Error : Trying to set a null mesh in error analysis\n";
    cout << "If the selected indicator don't need mesh don't call SetMesh function\n";
    cout << "The adaptive modulus is going down...\n";
    exit (-1);
  }
  fMesh = cmesh;  
}

void TPZErrorIndicator::SetSolution(TPZVec<REAL> &sol, int nstate, TPZVec<int> &dimstate, TPZVec<int> &statetoanalyse){
  if (nstate <= 0){
    cout << "TPZErrorIndicator::SetSolution - Error : The number of state variables must be greater than zero\n";
    cout << "The adaptive modulus is going down...\n";
    exit (-1);
  }
  
  fNState = nstate;
  fState = statetoanalyse;
  fDim = dimstate;
  
  int i,dims  = 0;
  for (i=0;i<fNState;i++) dims += fDim[i];
  
  fNDataEl = dims;

  
  if (sol.NElements() % dims != 0){
    cout << "TPZErrorIndicator::SetSolution - Warning : The number of elements in solution vector is not multiple of state variables seted.\n";
  }
  
  fSolution = sol;
  fNElements = fSolution.NElements() / dims;
  if (fMesh && fMesh->ElementVec().NElements() != fNElements){
    cout << "TPZErrorIndicator::SetSolution - Warning : The number of elements in solution vector is not equal to the number of elements in the mesh.\n";
  }
  
}

void TPZErrorIndicator::SetError(TPZVec<REAL> &maxerror, TPZVec<REAL> &minerror, TPZVec<int> &erantype){
  if  (maxerror.NElements() != fState.NElements() || minerror.NElements() != fState.NElements() || erantype.NElements() != fState.NElements()){
    cout << "TPZErrorIndicator::SetError - Error : The error vector must have the dimension equal to the number of state variables\n";
    cout << "Try to set the number of state variables (by SetSolution) first...\n";
    exit (-1);
  }
  fMaxError = maxerror;
  fMinError = minerror;
  fErAnType = erantype;
}

int TPZErrorIndicator::Index(int elem, int state, int dim){
  int index = fNDataEl * elem;
  int i;
  for (i=0;i<state;i++) index+= fDim[i];
  index += dim;
  return index;
}

void TPZErrorIndicator::Sort(TPZFMatrix &error, TPZFMatrix &perm) {
  int i,j,k;
  int imin = 0;
  int imax = error.Rows();
  perm.Resize(imax,error.Cols());
  for (i=0;i<imax;i++)
    for (j=0;j<error.Cols();j++) perm(i,j) = i;
  	for(i=imin; i<imax; i++) {
    	for(j=i+1; j<imax; j++) {
      	for (k=0;k<error.Cols();k++){
					if(error((int)perm(i,k)) < error((int)perm(j,k))) {
	  				int kp = (int) perm(i,k);
	  				perm(i,k) = perm(j,k);
	  				perm(j,k) = kp;
					}
      	}
    	}
  }
}


/** Returns the side to refine.
@param cel element to analyse
@param sidedim dimension of the sides which will be analysed */
int TPZErrorIndicator::GetRefSide(TPZCompEl *cel, int sidedim, int sidestate, TPZMatrix *errormat){
  int e,s,nsides = cel->Reference()->NSides();
  REAL sum = 0.;
  int side = -1;
  for (s=0;s<nsides;s++){
    TPZCompElSide celside (cel,s);
    if (celside.Reference().Dimension() != sidedim) continue;
    TPZStack<TPZCompElSide> elsidevec;
    celside.EqualLevelElementList(elsidevec,0,0);
    REAL auxsum = 0;
    int neigbyside = elsidevec.NElements();
    for (e=0;e<neigbyside;e++){
      int index = elsidevec[e].Element()->Index();
      double val = 0.;
      val = errormat->Get(e,sidestate);
      auxsum += val;
    }
    auxsum /= ((REAL) neigbyside);
    if(auxsum > sum){
      sum = auxsum;
      side = s;
    }
  }
  return side;
}

