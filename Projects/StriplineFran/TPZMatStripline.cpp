#include "TPZMatStripline.h"
#include "pzbndcond.h"

TPZMatStripline::TPZMatStripline(int id, REAL freq, STATE (& ur)( TPZVec<REAL>),STATE (& er)( TPZVec<REAL>)) : TPZVecL2(id), fUr(ur), fEr(er), fFreq(freq)
{
  
}

TPZMatStripline::TPZMatStripline(int id) : TPZVecL2(id), fUr(urDefault),
fEr(erDefault), fFreq(1.0e9)
{
  
}

/** @brief Default constructor */
TPZMatStripline::TPZMatStripline() : TPZVecL2(), fUr(urDefault),
fEr(erDefault), fFreq(1.0e9)
{
  
}


TPZMatStripline::TPZMatStripline(const TPZMatStripline &mat) : TPZVecL2(mat), fUr(mat.fUr),
fEr(mat.fEr), fFreq(mat.fFreq)
{
  
}

TPZMatStripline::~TPZMatStripline()
{
    
}

void TPZMatStripline::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
  TPZManVector<REAL,3> &x = data.x;
  TPZManVector<REAL,3> &normal = data.normal;
  const STATE mu = fUr(x);
  const STATE epsilon = fEr(x);
  
  TPZManVector<STATE,3> force(3);
  if(fForcingFunction) {
    fForcingFunction->Execute(x,force);
  }
  
  // Setting the phis
  TPZFMatrix<REAL> &phiQ = data.phi;
  //TPZFMatrix<REAL> &dphiQ = data.dphix; //AQUIFRAN
  TPZFMatrix<REAL> &dphiQdaxes = data.dphix;
  TPZFNMatrix<3,REAL> dphiQ;
  TPZAxesTools<REAL>::Axes2XYZ(dphiQdaxes, dphiQ, data.axes);
  
  TPZFNMatrix<9,REAL> auxMatrix(3,3,0.);
  auxMatrix(0,1) = -1.;
  auxMatrix(1,0) = 1.;
  auxMatrix(2,2) = 1.;
  const TPZFNMatrix<9,REAL> rotationMatrix = auxMatrix;
  int phrq;
  phrq = data.fVecShapeIndex.NElements();
  for(int iq=0; iq<phrq; iq++)
  {
    //ef(iq, 0) += 0.;
    int ivecind = data.fVecShapeIndex[iq].first;
    int ishapeind = data.fVecShapeIndex[iq].second;
    
    TPZFNMatrix<3,REAL> ivec(3,1,0.);
    for(int id=0; id<3; id++){
      ivec(id,0) = data.fNormalVec(id,ivecind);
    }
    //ROTATE FOR HCURL
    //ivec.Print("before rotation");
    ivec = rotationMatrix * ivec;
    //ivec.Print("after rotation");
    
    STATE ff = 0.;
    for (int i=0; i<3; i++) {
      ff += ivec(i,0)*force[i];
    }
    
    ef(iq,0) += weight*ff*phiQ(ishapeind,0);
    
    for (int jq=0; jq<phrq; jq++)
    {
      TPZFNMatrix<3,REAL> jvec(3,1,0.);
      int jvecind = data.fVecShapeIndex[jq].first;
      int jshapeind = data.fVecShapeIndex[jq].second;
      
      for(int id=0; id<3; id++){
        jvec(id,0) = data.fNormalVec(id,jvecind);
      }
      
      //ROTATE FOR HCURL
      jvec = rotationMatrix * jvec;
      STATE curlJz, curlIz, phiIdotphiJ;
      
      curlIz = dphiQ(0,ishapeind) * ivec(1,0) - dphiQ(1,ishapeind) * ivec(0,0);
      curlJz = dphiQ(0,jshapeind) * jvec(1,0) - dphiQ(1,jshapeind) * jvec(0,0);
      phiIdotphiJ = ( phiQ(ishapeind,0) * ivec(0,0) * phiQ(jshapeind,0) * jvec(0,0));
      phiIdotphiJ += ( phiQ(ishapeind,0) * ivec(1,0) * phiQ(jshapeind,0) * jvec(1,0));
      STATE stiff = (1./mu) * ( curlJz * curlIz );
      stiff += -1. * fW * fW * epsilon * phiIdotphiJ;
      STATE bTerm = -1. * (1./mu) * (( -1. * curlJz * normal[1] ) * phiQ(ishapeind,0) * ivec (0,0) + ( 1. * curlJz * normal[0] ) * phiQ(ishapeind,0) * ivec (1,0));
      stiff+=bTerm;
      ek(iq,jq) += weight * stiff;
    }
  }
}


void TPZMatStripline::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
  DebugStop();
}

void TPZMatStripline::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef)
{
  DebugStop();
}

void TPZMatStripline::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
  TPZManVector<REAL,3> &x = data.x;
  //TPZManVector<REAL,3> &normal = data.normal;
  
  // Setting the phis
  TPZFMatrix<REAL> &phiQ = data.phi;
  //  std::cout<<"bc0: "<<data.XCenter[0]<<" "<<data.XCenter[1]<<std::endl;
  //  TPZFMatrix<REAL> &dphiQ = data.dphix; //AQUIFRAN
  TPZFMatrix<REAL> &dphiQdaxes = data.dphi;
  TPZFNMatrix<3,REAL> dphiQ;
  TPZAxesTools<REAL>::Axes2XYZ(dphiQdaxes, dphiQ, data.axes);
  int nshape=phiQ.Rows();
  REAL BIG = TPZMaterial::gBigNumber;
  STATE v1 = bc.Val1()(0,0);//sera posto na matriz K no caso de condicao mista
  STATE v2 = bc.Val2()(0,0);//sera posto no vetor F
  switch ( bc.Type() )
  {
    case 0:
      for(int i = 0 ; i<nshape ; i++)
      {
        const STATE rhs = phiQ(i,0) * BIG * v2;
        ef(i,0) += rhs*weight;
        for(int j=0;j<nshape;j++)
        {
          const STATE stiff = phiQ(i,0) * phiQ(j,0) * BIG;
          ek(i,j) += stiff*weight;
        }
      }
      break;
    case 1:
			for(int i = 0 ; i<nshape ; i++)
			{
				const STATE rhs = phiQ(i,0) * v2;
				ef(i,0) += rhs*weight;
			}
      break;
    case 2:
			for(int i = 0 ; i<nshape ; i++)
			{
				const STATE rhs = phiQ(i,0) * v2;
				ef(i,0) += rhs*weight;
				for(int j=0;j<nshape;j++)
				{
					const STATE stiff = phiQ(i,0) * phiQ(j,0) * v1;
					ek(i,j) += stiff*weight;
				}
			}
      break;
  }

}

void TPZMatStripline::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
  DebugStop();
}

void TPZMatStripline::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef)
{
  DebugStop();
}

void TPZMatStripline::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
  DebugStop();
}

int TPZMatStripline::VariableIndex(const std::string &name)
{
  if(name == "E") {
    return 2;
  }
  else if ( name == "curlE")
  {
    return 3;
  }
  return TPZMaterial::VariableIndex(name);
}

/**
 * @brief Returns the number of variables associated with the variable indexed by var.
 * @param var Index variable into the solution, is obtained by calling VariableIndex
 */
int TPZMatStripline::NSolutionVariables(int var)
{
  int nVar = 0;
  switch (var) {
    case 2:
      nVar = 3;
      break;
    case 3:
      nVar = 1;
      break;
    default:
      nVar = TPZMaterial::NSolutionVariables(var);
      break;
  }
  return nVar;
}

/** @brief Returns the solution associated with the var index based on the finite element approximation */
void TPZMatStripline::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout)
{
  
  //  TPZFNMatrix<9,STATE> auxMatrix(3,3,0.);
  //  auxMatrix(0,1) = -1.;
  //  auxMatrix(1,0) = 1.;
  //  auxMatrix(2,2) = 1.;
  //
  //  const TPZFNMatrix<9,STATE> rotationMatrix = auxMatrix;//AQUIFRAN
  //  TPZFNMatrix<3,STATE> originalSol(3,1,0.);
  //  TPZManVector<STATE,3> auxVec = data.sol[0];
  //  originalSol(0,0) = auxVec[0];
  //  originalSol(1,0) = auxVec[1];
  //  originalSol(2,0) = auxVec[2];
  //  originalSol = rotationMatrix * originalSol;
  //  //originalSol.Print("sol");
  //  Solout[0] = originalSol(0,0);
  //  Solout[1] = originalSol(1,0);
  //  Solout[2] = originalSol(2,0);
  TPZFNMatrix<9,STATE> auxMatrix(3,3,0.);
  auxMatrix(0,1) = -1.;
  auxMatrix(1,0) = 1.;
  auxMatrix(2,2) = 1.;
  
  const TPZFNMatrix<9,STATE> rotationMatrix = auxMatrix;//AQUIFRAN
  TPZFNMatrix<3,STATE> originalSol(3,1,0.);
  TPZManVector<STATE,3> auxVec = data.sol[0];
  originalSol(0,0) = auxVec[0];
  originalSol(1,0) = auxVec[1];
  originalSol(2,0) = auxVec[2];
  originalSol = rotationMatrix * originalSol;
  //originalSol.Print("sol");
  Solout[0] = originalSol(0,0);
  Solout[1] = originalSol(1,0);
  Solout[2] = originalSol(2,0);
  switch (var) {
    case 2://E
#ifdef STATE_COMPLEX
      //      Solout[0] = std::norm(Solout[0]);
      //      Solout[1] = std::norm(Solout[1]);
      //      Solout[2] = std::norm(Solout[2]);
#endif
      break;
      
    default:
      break;
  }
}

STATE urDefault( TPZVec<REAL> x )
{
  return 1.0;
}

STATE erDefault( TPZVec<REAL> x )
{
  return 1.0;
}









