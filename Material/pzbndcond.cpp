
#include "pzbndcond.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"

void TPZBndCond::Clone(std::map<int, TPZAutoPointer<TPZMaterial> > &matvec) {
  int matid = Id();

  TPZAutoPointer<TPZMaterial> refmaterial = Material();
  TPZAutoPointer<TPZMaterial> newrefmaterial;
  int refmatid = 0;
  if(refmaterial) {
    refmaterial->Clone(matvec);
    refmatid = refmaterial->Id();
    newrefmaterial = matvec[refmatid];
  }
  std::map<int, TPZAutoPointer<TPZMaterial> >::iterator matit;
  matit = matvec.find(matid);
  if(matit != matvec.end())
  {
    TPZAutoPointer<TPZMaterial> newmat = TPZAutoPointer<TPZMaterial>(new TPZBndCond(*this, newrefmaterial));
    matvec[matid] = newmat;
  }
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
                                   TPZFMatrix &axesleft, TPZFMatrix &axesright,
				     TPZFMatrix &ek,TPZFMatrix &ef) {
//   TPZDiscontinuousGalerkin *mat = dynamic_cast<TPZDiscontinuousGalerkin *>(fMaterial.operator ->());
//   if(!mat) return;
//   if(fForcingFunction) {
//       TPZManVector<REAL> result(fBCVal2.Rows());
//       fForcingFunction(x,result);
//       int i;
//       for(i=0; i<fBCVal2.Rows(); i++) {
// 	fBCVal2(i,0) = result[i];
//       }
//   }
//   
//   if( this->fValFunction ) {
//     TPZManVector<REAL> result(this->fBCVal2.Rows());
//     this->fValFunction( x, this->fBCVal1, result, this->fType );
//     int i;
//     for(i = 0; i < this->fBCVal2.Rows(); i++) {
//       this->fBCVal2(i,0) = result[i];
//     }
//   }//if     
//   
//   if(phiL.Rows() == 0) {
//     TPZManVector<REAL,3> nor(normal.NElements(),0.);
//     for(int i=0; i<nor.NElements(); i++) nor[i] = -normal[i];
//     mat->ContributeBCInterface(x,solR,dsolR,weight,nor,phiR,dphiR,axesright, ek,ef,*this);
//     return;
//   }
//   if(phiR.Rows() == 0) {
//     mat->ContributeBCInterface(x,solL,dsolL,weight,normal,phiL,dphiL,axesleft,ek,ef,*this);
//     return;
//   }
}

void TPZBndCond::ContributeInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL,TPZVec<REAL> &solR,TPZFMatrix &dsolL,
				   TPZFMatrix &dsolR,REAL weight,TPZVec<REAL> &normal,TPZFMatrix &phiL,
				   TPZFMatrix &phiR,TPZFMatrix &dphiL,TPZFMatrix &dphiR,
                                   TPZFMatrix &axesleft, TPZFMatrix &axesright,
				   TPZFMatrix &ef) {
  TPZDiscontinuousGalerkin *mat = dynamic_cast<TPZDiscontinuousGalerkin *>(fMaterial.operator ->());
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
    mat->ContributeBCInterface(x,solR,dsolR,weight,nor,phiR,dphiR,axesright,ef,*this);
    return;
  }
  if(phiR.Rows() == 0) {
    mat->ContributeBCInterface(x,solL,dsolL,weight,normal,phiL,dphiL,axesleft,ef,*this);
    return;
  }
}

void TPZBndCond::ContributeInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL,TPZVec<REAL> &solR,TPZFMatrix &dsolL,
			 TPZFMatrix &dsolR,REAL weight,TPZVec<REAL> &normal,TPZFMatrix &phiL,
			 TPZFMatrix &phiR,TPZFMatrix &dphiL,TPZFMatrix &dphiR,
                         TPZFMatrix &axesleft, TPZFMatrix &axesright,
			 TPZFMatrix &ek,TPZFMatrix &ef, int LeftPOrder, int RightPOrder, REAL faceSize){
   TPZDiscontinuousGalerkin *mat = dynamic_cast<TPZDiscontinuousGalerkin *>(this->fMaterial.operator ->());

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
      mat->ContributeBCInterface(x,solR,dsolR,weight,nor,phiR,dphiR,axesright,ek,ef,*this, RightPOrder, faceSize );
      return;
   }

   if(phiR.Rows() == 0) {
      mat->ContributeBCInterface(x,solL,dsolL,weight,normal,phiL,dphiL,axesleft,ek,ef,*this, LeftPOrder, faceSize);
      return;
   }

}

void TPZBndCond::InterfaceJumps(TPZVec<REAL> &x, TPZVec<REAL> &leftu, TPZVec<REAL> &leftNormalDeriv,
                                TPZVec<REAL> &rightu, TPZVec<REAL> &rightNormalDeriv,
                                TPZVec<REAL> &values){
   TPZDiscontinuousGalerkin *mat = dynamic_cast<TPZDiscontinuousGalerkin *>(this->fMaterial.operator ->());

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

   TPZDiscontinuousGalerkin *mat = dynamic_cast<TPZDiscontinuousGalerkin *>(this->fMaterial.operator ->());
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

void TPZBndCond::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef){
  if( this->HasForcingFunction() ){
    const int n = this->fBCVal2.Rows();
    TPZManVector<REAL> result(n);
    this->fForcingFunction(x,result);
    for(int i = 0; i < n; i++){
      this->fBCVal2(i,0) = result[i];
    }//for i
  }//if (hasforcintfunction)

  if( this->fValFunction ){
    const int n = this->fBCVal2.Rows();
    TPZManVector<REAL> result(n);
    this->fValFunction( x, this->fBCVal1, result, this->fType );
    for(int i = 0; i < n; i++){
      this->fBCVal2(i,0) = result[i];
    }//for i
  }//if

  //clone meshes required analysis
  int typetmp = fType;
  if (fType == 50){
    int i;
    for (i=0;i<sol.NElements();i++){
        fBCVal2(i,0) = gBigNumber*sol[i];
        fBCVal1(i,i) = gBigNumber;
    }
    fType = 2;
  }

  this->fMaterial->ContributeBC(data,weight,ek,ef,*this);
  fType = typetmp;
}

void TPZBndCond::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix &ef){

}

void TPZBndCond::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef, TPZBndCond &bc){
///nothing to be done here
}

void TPZBndCond::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix &ef, TPZBndCond &bc){
///nothing to be done here
}

void TPZBndCond::ContributeInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef){
  TPZDiscontinuousGalerkin *mat = dynamic_cast<TPZDiscontinuousGalerkin *>(fMaterial.operator ->());
  if(!mat) return;
  if(fForcingFunction){
    TPZManVector<REAL> result(fBCVal2.Rows());
    fForcingFunction(x,result);
    int i;
    for(i=0; i<fBCVal2.Rows(); i++) {
      fBCVal2(i,0) = result[i];
    }
  }

  if( this->fValFunction ){
    TPZManVector<REAL> result(this->fBCVal2.Rows());
    this->fValFunction( x, this->fBCVal1, result, this->fType );
    int i;
    for(i = 0; i < this->fBCVal2.Rows(); i++) {
      this->fBCVal2(i,0) = result[i];
    }
  }//if

  if(phiL.Rows() == 0) {//it meanst right data has been filled
    //left data should be filled instead of right data
    data.InvertLeftRightData();
  
    TPZManVector<REAL,3> nor(normal.NElements(),0.);
    for(int i=0; i<nor.NElements(); i++) nor[i] = -normal[i];
    mat->ContributeBCInterface(x,solR,dsolR,weight,nor,phiR,dphiR,axesright, ek,ef,*this);
    return;
  }
  if(phiR.Rows() == 0) {
    mat->ContributeBCInterface(x,solL,dsolL,weight,normal,phiL,dphiL,axesleft,ek,ef,*this);
    return;
  }
}

void TPZBndCond::ContributeInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ef){

}

void TPZBndCond::ContributeBCInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc){
///nothing to be done here
}

void TPZBndCond::ContributeBCInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ef,TPZBndCond &bc){
///nothing to be done here
}
