
#include "pzbndcond.h"
#include "pzadmchunk.h"

void TPZBndCond::Clone(TPZAdmChunkVector<TPZMaterial *> &matvec) {
  int matid = Id();
  int nmat = matvec.NElements();
  int m;

  TPZMaterial *refmaterial = Material();
  TPZMaterial *newrefmaterial = 0;
  int refmatid = 0;
  if(refmaterial) {
    refmaterial->Clone(matvec);
    nmat = matvec.NElements();
    refmatid = refmaterial->Id();
    for(m=0; m<nmat; m++) {
      TPZMaterial *mat = matvec[m];
      if(!mat) continue;
      if(mat->Id() == refmatid) {
	newrefmaterial = mat;
	break;
      }
    }
  }
  for(m=0; m<nmat; m++) {
    TPZMaterial *mat = matvec[m];
    if(!mat) continue;
    if(mat->Id() == matid) return;
  }
  int vecpos = matvec.AllocateNewElement();
  TPZMaterial *newmat = new TPZBndCond(*this, newrefmaterial);
  matvec[vecpos] = newmat;
}

#ifdef _AUTODIFF

void TPZBndCond::ContributeEnergy(TPZVec<REAL> &x,
	TPZVec<FADFADREAL> &sol, TPZVec<FADFADREAL> &dsol,
	FADFADREAL &U, REAL weight)
{

    int typetmp = fType;
    if (fType == 50){
	    int i;
	    for (i=0;i<sol.NElements();i++){
		fBCVal2(i,0) = gBigNumber*sol[i].val().val();
		fBCVal1(i,i) = gBigNumber;
	    }
	    fType = 2;
    }
    fMaterial->ContributeBCEnergy(x,sol,U,weight,*this);
    fType = typetmp;
  

}

#endif

void TPZBndCond::ContributeInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL,TPZVec<REAL> &solR,TPZFMatrix &dsolL,
				   TPZFMatrix &dsolR,REAL weight,TPZVec<REAL> &normal,TPZFMatrix &phiL,
				   TPZFMatrix &phiR,TPZFMatrix &dphiL,TPZFMatrix &dphiR,
				     TPZFMatrix &ek,TPZFMatrix &ef) {
  TPZDiscontinuousGalerkin *mat = dynamic_cast<TPZDiscontinuousGalerkin *>(fMaterial);
  if(!mat) return;
  if(phiL.Rows() == 0) {
    TPZManVector<REAL,3> nor(normal.NElements(),0.);
    for(int i=0; i<nor.NElements(); i++) nor[i] = -normal[i];
    mat->ContributeBCInterface(x,solR,dsolR,weight,nor,phiR,dphiR,ek,ef,*this);
    return;
  }
  if(phiR.Rows() == 0) {
    mat->ContributeBCInterface(x,solL,dsolL,weight,normal,phiL,dphiL,ek,ef,*this);
    return;
  }
}
