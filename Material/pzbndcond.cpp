
#include "pzbndcond.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"

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

//#ifdef _AUTODIFF

/*
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
*/
//#endif

void TPZBndCond::ContributeInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL,TPZVec<REAL> &solR,TPZFMatrix &dsolL,
				   TPZFMatrix &dsolR,REAL weight,TPZVec<REAL> &normal,TPZFMatrix &phiL,
				   TPZFMatrix &phiR,TPZFMatrix &dphiL,TPZFMatrix &dphiR,
				     TPZFMatrix &ek,TPZFMatrix &ef) {
  TPZDiscontinuousGalerkin *mat = dynamic_cast<TPZDiscontinuousGalerkin *>(fMaterial);
  if(!mat) return;
  if(fForcingFunction) {
      TPZManVector<REAL> result(fBCVal2.Rows());
      fForcingFunction(x,result);
      int i;
      for(i=0; i<fBCVal2.Rows(); i++) {
	fBCVal2(i,0) = result[i];
      }
  }
  
  if( this->fValFunction ) {
    TPZManVector<REAL> result(this->fBCVal2.Rows());
    this->fValFunction( x, this->fBCVal1, result, this->fType );
    int i;
    for(i = 0; i < this->fBCVal2.Rows(); i++) {
      this->fBCVal2(i,0) = result[i];
    }
  }//if     
  
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

void TPZBndCond::ContributeInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL,TPZVec<REAL> &solR,TPZFMatrix &dsolL,
				   TPZFMatrix &dsolR,REAL weight,TPZVec<REAL> &normal,TPZFMatrix &phiL,
				   TPZFMatrix &phiR,TPZFMatrix &dphiL,TPZFMatrix &dphiR,
				   TPZFMatrix &ef) {
  TPZDiscontinuousGalerkin *mat = dynamic_cast<TPZDiscontinuousGalerkin *>(fMaterial);
  if(!mat) return;
  if(fForcingFunction) {
      TPZManVector<REAL> result(fBCVal2.Rows());
      fForcingFunction(x,result);
      int i;
      for(i=0; i<fBCVal2.Rows(); i++) {
	fBCVal2(i,0) = result[i];
      }
  }
  
  if( this->fValFunction ) {
    TPZManVector<REAL> result(this->fBCVal2.Rows());
    this->fValFunction( x, this->fBCVal1, result, this->fType );
    int i;
    for(i = 0; i < this->fBCVal2.Rows(); i++) {
      this->fBCVal2(i,0) = result[i];
    }
  }//if
  
  if(phiL.Rows() == 0) {
    TPZManVector<REAL,3> nor(normal.NElements(),0.);
    for(int i=0; i<nor.NElements(); i++) nor[i] = -normal[i];
    mat->ContributeBCInterface(x,solR,dsolR,weight,nor,phiR,dphiR,ef,*this);
    return;
  }
  if(phiR.Rows() == 0) {
    mat->ContributeBCInterface(x,solL,dsolL,weight,normal,phiL,dphiL,ef,*this);
    return;
  }
}

void TPZBndCond::ContributeInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL,TPZVec<REAL> &solR,TPZFMatrix &dsolL,
			 TPZFMatrix &dsolR,REAL weight,TPZVec<REAL> &normal,TPZFMatrix &phiL,
			 TPZFMatrix &phiR,TPZFMatrix &dphiL,TPZFMatrix &dphiR,
			 TPZFMatrix &ek,TPZFMatrix &ef, int LeftPOrder, int RightPOrder, REAL faceSize){
   TPZDiscontinuousGalerkin *mat = dynamic_cast<TPZDiscontinuousGalerkin *>(this->fMaterial);

   if(!mat) return;
   if(fForcingFunction) {
      TPZManVector<REAL> result(fBCVal2.Rows());
      fForcingFunction(x,result);
      int i;
      for(i=0; i<fBCVal2.Rows(); i++) {
	fBCVal2(i,0) = result[i];
      }
   }
   
  if( this->fValFunction ) {
    TPZManVector<REAL> result(this->fBCVal2.Rows());
    this->fValFunction( x, this->fBCVal1, result, this->fType );
    int i;
    for(i = 0; i < this->fBCVal2.Rows(); i++) {
      this->fBCVal2(i,0) = result[i];
    }
  }//if

   if(phiL.Rows() == 0) {
      TPZManVector<REAL,3> nor(normal.NElements(),0.);
      for(int i=0; i<nor.NElements(); i++) nor[i] = -normal[i];
      mat->ContributeBCInterface(x,solR,dsolR,weight,nor,phiR,dphiR,ek,ef,*this, RightPOrder, faceSize );
      return;
   }

   if(phiR.Rows() == 0) {
      mat->ContributeBCInterface(x,solL,dsolL,weight,normal,phiL,dphiL,ek,ef,*this, LeftPOrder, faceSize);
      return;
   }

}

void TPZBndCond::InterfaceJumps(TPZVec<REAL> &x, TPZVec<REAL> &leftu, TPZVec<REAL> &leftNormalDeriv,
                                TPZVec<REAL> &rightu, TPZVec<REAL> &rightNormalDeriv,
                                TPZVec<REAL> &values){
   TPZDiscontinuousGalerkin *mat = dynamic_cast<TPZDiscontinuousGalerkin *>(this->fMaterial);

   if(!mat) return;
   if(fForcingFunction) {
      TPZManVector<REAL> result(fBCVal2.Rows());
      fForcingFunction(x,result);
      int i;
      for(i=0; i<fBCVal2.Rows(); i++) {
        fBCVal2(i,0) = result[i];
      }
   }

   if(leftu.NElements() == 0) {
      mat->BCInterfaceJumps(rightu, *this, values);
      return;
   }
   
   if(rightu.NElements() == 0) {
      mat->BCInterfaceJumps(leftu, *this, values);
      return;
   }   
   
   std::cout << __PRETTY_FUNCTION__ << " - Huge problem. Both leftu and rightu contain elements. Wich one is the actual element neighbour to the Dirichlet boundary ?" << std::endl;

}//InterfaceJumps

int TPZBndCond::ClassId() const
{
   return TPZBNDCONDID;
}
#ifndef BORLAND
template class TPZRestoreClass<TPZBndCond,TPZBNDCONDID>;
#endif

void TPZBndCond::Write(TPZStream &buf, int withclassid)
{
   TPZMaterial::Write(buf, withclassid);
   buf.Write(&fType, 1);
   fBCVal1.Write(buf, 0);
   fBCVal2.Write(buf, 0);
   int MatId =fMaterial->Id();
   buf.Write(&MatId, 1);
}

void TPZBndCond::Read(TPZStream &buf, void *context)
{
   TPZMaterial::Read(buf, context);
   buf.Read(&fType, 1);
   fBCVal1.Read(buf, 0);
   fBCVal2.Read(buf, 0);
   int MatId;
   buf.Read(&MatId,1);
   TPZCompMesh * pCM = (TPZCompMesh * )/*dynamic_cast<TPZCompMesh *>*/(context);
   fMaterial = pCM->FindMaterial(MatId);
}

void TPZBndCond::ContributeInterfaceErrors(TPZVec<REAL> &x,
                                       TPZVec<REAL> &solL,
                                       TPZVec<REAL> &solR,
                                       TPZFMatrix &dsolL, 
                                       TPZFMatrix &dsolR,
                                       REAL weight,
                                       TPZVec<REAL> &normal,
                                       TPZVec<REAL> &nkL, 
                                       TPZVec<REAL> &nkR,
                                       int LeftPOrder, 
                                       int RightPOrder, 
                                       REAL faceSize,
                                       int &errorid){

   TPZDiscontinuousGalerkin *mat = dynamic_cast<TPZDiscontinuousGalerkin *>(this->fMaterial);
   if(!mat) return;

   if(fForcingFunction) {
      TPZManVector<REAL> result(fBCVal2.Rows());
      fForcingFunction(x,result);
      int i;
      for(i=0; i<fBCVal2.Rows(); i++) {
        fBCVal2(i,0) = result[i];
      }
   }

   int POrder;
   TPZVec<REAL> sol;  
   TPZFMatrix dsol;
   TPZVec<REAL> nk;
   if(solL.NElements() < solR.NElements()){
      POrder=  RightPOrder;
      sol = solR;  
      dsol = dsolR;
      nk= nkR;
      for(int i=0; i<normal.NElements(); i++) normal[i] = -normal[i];
   }
   else {
     POrder=  LeftPOrder;
     sol = solL;
     dsol = dsolL;
     nk= nkL;
   }
 
   mat->ContributeInterfaceBCErrors(x,sol,dsol,weight,normal,nk,*this , POrder, faceSize, errorid);
}
