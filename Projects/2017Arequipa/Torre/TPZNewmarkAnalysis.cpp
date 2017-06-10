#include "TPZNewmarkAnalysis.h"
#include "pzsolve.h"
#include "pzstepsolver.h"
#include "pzskylmat.h"
#include "TPZEulerBernoulliBC.h"
#include "TPZEulerBernoulliBeam.h"
#include "TSWXExportFrameMesh.h"

TPZNewmarkAnalysis::TPZNewmarkAnalysis() : TPZAnalysis(),
  fPropertyData(), fLastSol(), fNextSol(),
  fMassMatrix(), fViscousMatrix(), fStiffnessMatrix(), fDeltaT(1.)
{

  this->SetViscousMatrixCoefficients();
  this->SetNewmarkParameters();
}

///////////////////////////////////////////////////////////////////////////////

TPZNewmarkAnalysis::TPZNewmarkAnalysis(TPZCompMesh *mesh,
                                       TPZAutoPointer< TPZEulerBernoulliBeamData > PropertyData)
                      : TPZAnalysis(mesh),
  fPropertyData(PropertyData), fLastSol(), fNextSol(),
  fMassMatrix(), fViscousMatrix(), fStiffnessMatrix(), fDeltaT(1.){

  this->SetViscousMatrixCoefficients();
  this->SetNewmarkParameters();
}

///////////////////////////////////////////////////////////////////////////////

TPZNewmarkAnalysis::~TPZNewmarkAnalysis(){

}

///////////////////////////////////////////////////////////////////////////////

void TPZNewmarkAnalysis::SetInitialSolution( const TPZFMatrix<STATE> & Displacement, const TPZFMatrix<STATE> & Velocity, const TPZFMatrix<STATE> & Acceleration ){
  this->fLastSol.fDisplacement = Displacement;
  this->fLastSol.fVelocity = Velocity;
  this->fLastSol.fAcceleration = Acceleration;
}

///////////////////////////////////////////////////////////////////////////////

void TPZNewmarkAnalysis::SetViscousMatrixCoefficients(REAL ViscousMassCoeff, REAL ViscousStiffCoef){
  this->fViscousMassCoeff = ViscousMassCoeff;
  this->fViscousStiffCoef = ViscousStiffCoef;
}

///////////////////////////////////////////////////////////////////////////////

void TPZNewmarkAnalysis::SetNewmarkParameters(REAL NewmarkBeta, REAL NewmarkGamma){
  this->fNewmarkBeta = NewmarkBeta;
  this->fNewmarkGamma = NewmarkGamma;
}

///////////////////////////////////////////////////////////////////////////////

void TPZNewmarkAnalysis::TimeSteps(std::ostream &out,
                 REAL DeltaT,
                 REAL SteadyStateTol, int MaxTimeSteps,
                 REAL NewtonTol, int NewtonMaxIter,
                 const std::string &filename){
  //time step
  this->fDeltaT = DeltaT;

  //error file
  std::stringstream errorFileName;
  errorFileName << filename << ".error";
  std::ofstream errorFile(errorFileName.str().c_str());

  ///initial state
  if(this->fLastSol.HasData() == false) DebugStop();
  this->fNextSol = this->fLastSol;
//    TPZManVector<std::string> scalarnames(0), vecnames(1);
 //   scalarnames[0] = "DisplacementX";
//    scalarnames[1] = "SigmaY";
//    scalarnames[2] = "Pressure";
//    vecnames[0] = "Displacement";
//    vecnames[1] = "NormalForce";
 //   vecnames[1] = "";
//    DefineGraphMesh(3,scalarnames,vecnames,"torre.vtk");
    
  this->PostProcess(0, 0., filename);
//    TPZAnalysis::PostProcess(1,3);

  ///loop over time steps
  for(int istep = 0; istep < MaxTimeSteps; istep++){
    const REAL currTime = istep*fDeltaT;
    out << "\nSolving istep / time : " << istep << " / " << currTime << "\n";
    bool linesearch = false;
    //std::pair< int, REAL>(iter,error)
    std::pair< int, REAL> NewtonError = this->SolveOneStep(istep, out, NewtonTol, NewtonMaxIter, linesearch);
    if(NewtonError.second > NewtonTol){
      errorFile << "Newton did not converge! istep = " << istep << " , Error = " << NewtonError << " , tolerance = " << NewtonTol << "\n";
      errorFile.flush();
    }

    //post processing solution
    this->PostProcess(istep+1, currTime, filename);
//      TPZAnalysis::PostProcess(1,3);

    //step in time
    this->fLastSol = this->fNextSol;
  }//for istep


}//void

///////////////////////////////////////////////////////////////////////////////

void TPZNewmarkAnalysis::AssembleMassMatrix(const TPZFMatrix<STATE> & u){
  this->fPropertyData->SetComputeMass();
  this->Solution() = u;
  this->LoadSolution();
  this->Assemble();
  this->fMassMatrix = this->Solver().Matrix()->Clone();
}

///////////////////////////////////////////////////////////////////////////////

void TPZNewmarkAnalysis::AssembleStiffnessMatrix(const TPZFMatrix<STATE> & u){
  this->fPropertyData->SetComputeStiffness();
  this->Solution() = u;
  this->LoadSolution();
  this->Assemble();
  this->fStiffnessMatrix = this->Solver().Matrix()->Clone();
}

///////////////////////////////////////////////////////////////////////////////

void TPZNewmarkAnalysis::AssembleElasticForces(const TPZFMatrix<STATE> & u){
  this->fPropertyData->SetComputeStiffness();
  this->Solution() = u;
  this->LoadSolution();
  this->AssembleResidual();
}

///////////////////////////////////////////////////////////////////////////////

void TPZNewmarkAnalysis::BuildViscousMatrix(){
  TPZSkylMatrix<STATE> * M = dynamic_cast< TPZSkylMatrix<STATE>* >(this->fMassMatrix.operator->());
  TPZSkylMatrix<STATE> * K = dynamic_cast< TPZSkylMatrix<STATE>* >(this->fStiffnessMatrix.operator->());
  if( !M || !K ) DebugStop();

  this->fViscousMatrix = M->Clone();//C = M
  TPZSkylMatrix<STATE> * C = dynamic_cast< TPZSkylMatrix<STATE>* >(this->fViscousMatrix.operator->());
  if(!C) DebugStop();

  C->operator*=( this->fViscousMassCoeff );//C = fViscousMassCoeff M
  C->AddSameStruct(*K, this->fViscousStiffCoef); //C = fViscousMassCoeff M + fViscousStiffCoef K
}

///////////////////////////////////////////////////////////////////////////////

void TPZNewmarkAnalysis::UpdateNextSol(const TPZFMatrix<STATE> &nextAcceleration){
  TPZFMatrix<STATE> aux;

  //acceleration
  this->fNextSol.fAcceleration = nextAcceleration;

  //velocity
  this->fNextSol.fVelocity = fLastSol.fVelocity;//v(n+1) = v(n)
  aux = this->fLastSol.fAcceleration;
  aux *= (1.-fNewmarkGamma)*fDeltaT;
  this->fNextSol.fVelocity += aux;//v(n+1) = v(n) + (1-gamma)deltaT a(n)

  aux = this->fNextSol.fAcceleration;
  aux *= fNewmarkGamma*fDeltaT;
  this->fNextSol.fVelocity += aux;//v(n+1) = v(n) + (1-gamma)dT a(n) + gamma dT a(n+1)

  //displacement
  this->fNextSol.fDisplacement = this->fLastSol.fDisplacement;//u(n+1) = u(n)

  aux = this->fLastSol.fVelocity;
  aux *= fDeltaT;
  this->fNextSol.fDisplacement += aux;//u(n+1) = u(n) + dT v(n)

  aux = this->fLastSol.fAcceleration;
  aux *= 1./2.*fDeltaT*fDeltaT*(1.-2.*fNewmarkBeta);
  this->fNextSol.fDisplacement += aux;//u(n+1) = u(n) + dT v(n) + 1/2 dT^2 (1-2 beta) a(n)

  aux = this->fNextSol.fAcceleration;
  aux *= 1./2.*fDeltaT*fDeltaT*(+2.*fNewmarkBeta);
  this->fNextSol.fDisplacement += aux;//u(n+1) = u(n) + dT v(n) + 1/2 dT^2 (1-2 beta) a(n) + 1/2 dT^2 2 beta a(n+1)

}

///////////////////////////////////////////////////////////////////////////////

void TPZNewmarkAnalysis::ResidualAndJacobian(bool updateJacobian){

  TPZFMatrix<STATE> aux;
  if(updateJacobian){
    this->AssembleStiffnessMatrix(this->fNextSol.fDisplacement);
  }
  else{
    this->AssembleElasticForces(this->fNextSol.fDisplacement);
  }

  this->fRhs *= -1.;// R = - (F-Ku) = Ku-F = K u(n+1) - F(n+1)

   fMassMatrix->Multiply(this->fNextSol.fAcceleration, aux);
  this->fRhs += aux; // R = M a(n+1) + K u(n+1) - F(n+1)
#ifdef apagar
  this->fViscousMatrix->Multiply(this->fNextSol.fAcceleration, aux);
  aux *= fNewmarkGamma * fDeltaT;
  this->fRhs += aux;// R = (M + C gamma dT) a(n+1) + K u(n+1) - F(n+1)

  this->fViscousMatrix->Multiply(this->fLastSol.fVelocity, aux);
  this->fRhs += aux;// R = (M + C gamma dT) a(n+1) + K u(n+1) - F(n+1) + C v(n)

  this->fViscousMatrix->Multiply(this->fLastSol.fAcceleration, aux);
  aux *= (1.-fNewmarkGamma)*fDeltaT;
  this->fRhs += aux;// R = (M + C gamma dT) a(n+1) + K u(n+1) - F(n+1) + C v(n) + C(1-gamma)dT a(n)
#endif
  if(updateJacobian){
    TPZSkylMatrix<STATE> * J = dynamic_cast< TPZSkylMatrix<STATE> * >(this->Solver().Matrix().operator->());
    TPZSkylMatrix<STATE> * M = dynamic_cast< TPZSkylMatrix<STATE>* >(this->fMassMatrix.operator->());
    TPZSkylMatrix<STATE> * K = dynamic_cast< TPZSkylMatrix<STATE>* >(this->fStiffnessMatrix.operator->());
    TPZSkylMatrix<STATE> * C = dynamic_cast< TPZSkylMatrix<STATE>* >(this->fViscousMatrix.operator->());
    if( !J || !M || !K || !C ) DebugStop();

    //J = K
    J->operator*=( fNewmarkBeta*fDeltaT*fDeltaT );//J = beta dt^2 K
#ifdef apagar
    J->AddSameStruct(*C, fNewmarkGamma*fDeltaT);//J = C gamma dT + beta dt^2 K
#endif
    J->AddSameStruct(*M,1.);//J = M + C gamma dT + beta dt^2 K
  }

}

///////////////////////////////////////////////////////////////////////////////

std::pair<int,REAL> TPZNewmarkAnalysis::SolveOneStep(int TimeStep, std::ostream &out,
                                                     REAL tol,int numiter, bool linesearch){
  const int numeq = fCompMesh->NEquations();

  //Creating mass and viscous matrices. Stiffness matrix will be updated in Newton's loop
  this->AssembleMassMatrix( this->fLastSol.fDisplacement );
  this->AssembleStiffnessMatrix( this->fLastSol.fDisplacement );
  this->BuildViscousMatrix();//requires Mass and Stiffness computed

  TPZFMatrix<STATE> NextAcceleration( this->fLastSol.fAcceleration.Rows(), 1, 0. );
  this->UpdateNextSol(NextAcceleration);
  int iter = 0;
  REAL error = 2.*tol + 1.e10;
  while(error > tol && iter < numiter) {

/////////////tamarindo
  if(0){
    this->UpdateNextSol(NextAcceleration);
    this->ResidualAndJacobian(true);
    TPZSkylMatrix<STATE> * J = dynamic_cast< TPZSkylMatrix<STATE> * >(this->Solver().Matrix().operator->());
    std::ofstream myfile("c:\\Temp\\CEATI\\jac.nb");
    J->Print("J=",myfile,EMathematicaInput);

    TPZFMatrix<STATE> deltaAcc(12,1,0.), NewRhs;
    for(int i = 0; i < 12; i++) deltaAcc(i,0) = 0.1 * 0.01*i;
    J->Multiply(deltaAcc,NewRhs);
    NewRhs += fRhs;

    NextAcceleration += deltaAcc;
    this->UpdateNextSol(NextAcceleration);
    this->ResidualAndJacobian(false);
    fRhs.Print("fRhs=",myfile,EMathematicaInput);
    NewRhs.Print("NewRhs=",myfile,EMathematicaInput);
    myfile << "Chop[Table[ 100 (fRhs[[i, 1]] - NewRhs[[i, 1]])/fRhs[[i, 1]], {i, 1, 12}]]" << std::endl;
    fRhs -= NewRhs;
   }
//////////////

      this->UpdateNextSol(NextAcceleration);
      if(iter % 10 == 0) this->ResidualAndJacobian(true);
      else this->ResidualAndJacobian(false);
      Solve();
      fSolution *= -1.;
      const REAL normDeltaSol = Norm(fSolution);

      REAL alphaLineSearch = 1.;
      if (linesearch){
        TPZFMatrix<STATE> nextSol;
        REAL LineSearchTol = 1e-3 * Norm(fSolution);
        const int niter = 10;
        alphaLineSearch = this->LineSearch(NextAcceleration, this->fSolution, LineSearchTol, niter, nextSol);
        fSolution = nextSol;
      }
      else{
        NextAcceleration += fSolution;
      }

      //Computing rhs
      this->UpdateNextSol(NextAcceleration);
      this->ResidualAndJacobian(false);
      std::pair<REAL,REAL> norm = this->RhsNorm(fRhs);
      const REAL NormRhs = sqrt(norm.first*norm.first + norm.second*norm.second);
      out << "Newton iteration : " << iter << ", alphaLS = " << alphaLineSearch
          << " : |Delta(acceleration)| , |Delta(rhs)| : " << normDeltaSol << " / " << NormRhs << std::endl;

      const REAL normAdopted = NormRhs; //normDeltaSol
      if(normAdopted < tol) {
         out << "Convergence reached\n";
      }
      else if( (normAdopted - error) > 1.e-3 ) {
         out << "Divergent Method\n";
      }
      error = normAdopted;
	    iter++;
	    out.flush();
   }

   return std::pair< int, REAL>(iter,error);
}

///////////////////////////////////////////////////////////////////////////////

void TPZNewmarkAnalysis::PostProcess(int istep, REAL time, const std::string &filename){
  TSWXExportFrameMesh graph;
  this->Solution() = this->fNextSol.fDisplacement;
  this->LoadSolution();
  graph.GenerateGraphMesh(*this->Mesh(),time);
  std::stringstream name;
  name << filename << istep << ".vtk";
  std::ofstream paraviewfile(name.str().c_str());
  graph.Mesh().ToParaview(paraviewfile);
}

///////////////////////////////////////////////////////////////////////////////

void TPZNewmarkAnalysis::SupportReaction(int istep, REAL time, int supportId, std::ostream &myfile){
//implement me
}

///////////////////////////////////////////////////////////////////////////////

void TPZNewmarkAnalysis::NodalDisplacement(int istep, REAL time, int nodeId, std::ostream &myfile){
//implement me
}

///////////////////////////////////////////////////////////////////////////////

std::pair<REAL,REAL> TPZNewmarkAnalysis::RhsNorm(const TPZFMatrix<STATE> &origrhs) const{
  //removing bignumbers
  TPZFMatrix<STATE> cpRhs = origrhs;
  TPZManVector<int> indices(6);
  int nel = this->Mesh()->NElements();
  for(int iel = 0; iel < nel; iel++){
    TPZEulerBernoulliBC * bc = dynamic_cast<TPZEulerBernoulliBC*>(this->Mesh()->ElementVec()[iel]);
    if(!bc) continue;
    bc->GetEquationIndices(indices);
    for(int i = 0; i < indices.NElements(); i++){
      if(bc->GetBCVal(i).fType == TPZEulerBernoulliBC::ESupport){
        cpRhs( indices[i], 0 ) = 0.;
      }
    }
  }
  return std::pair<REAL,REAL>(Norm(cpRhs),0.);//tamarindo melhor aqui
}

///////////////////////////////////////////////////////////////////////////////

//#define DEBUGLINESEARCH
#ifdef DEBUGLINESEARCH
static ofstream alphafile("c:\\Temp\\tmp\\alpha.txt");
#endif
REAL TPZNewmarkAnalysis::LineSearch(const TPZFMatrix<STATE> &Acceleration,
                                    const TPZFMatrix<STATE> DeltaAccel,
                                    REAL tol,
                                    int niter,
                                    TPZFMatrix<STATE> &NextAcc){
  REAL error = 2.*tol+1.;
  REAL A, B, L, M;
  TPZFMatrix<STATE> ak, bk, lambdak, muk, Interval;
  REAL NormResLambda, NormResMu;
  ///ak = Wn + 0.0001 * DeltaW
  ak = DeltaAccel;
  A = 0.0001; // preciso disso nos cabos
  ak *= A;
  ak += Acceleration;
  ///bk = Wn + DeltaW
  bk = DeltaAccel;
  B = 1.6181229773462784*(1.-0.382*A);//B is set such that M = 1
  bk *= B;
  bk += Acceleration;

  ///calculando residuo em A e B
  REAL NormResA, NormResB;
  {
    ///computing residual
    this->UpdateNextSol(ak);
    this->ResidualAndJacobian(false);
    NormResA = Norm(fRhs);
  }
  {
    ///computing residual
    this->UpdateNextSol(bk);
    this->ResidualAndJacobian(false);
    NormResB = Norm(fRhs);
  }

  ///Interval = (bk-ak)
  Interval = bk; Interval -= ak;
  int iter = 0;
  int KeptVal = -1;///0 means I have residual(labmda); 1 means I have residual(mu); -1 means I have nothing
  while(error > tol && iter < niter){
    iter++;

    if (KeptVal != 0){
      L = 0.382*(B-A)+A;
      ///lambdak = ak + 0.382*(bk-ak)
      lambdak = Interval; lambdak *= 0.382; lambdak += ak;
      this->UpdateNextSol(lambdak);
      this->ResidualAndJacobian(false);
      NormResLambda = Norm(fRhs);
    }

    if (KeptVal != 1){
      ///muk = ak + 0.618*(bk-ak)
      M = 0.618*(B-A)+A;
      muk = Interval; muk *= 0.618; muk += ak;
      this->UpdateNextSol(muk);
      this->ResidualAndJacobian(false);
      NormResMu = Norm(fRhs);
    }

    if (NormResLambda > NormResMu){
      A = L;
      NormResA = NormResLambda;
      L = M;
      ak = lambdak;
      lambdak = muk;
      NormResLambda = NormResMu;
      KeptVal = 0;
    }
    else{
      B = M;
      NormResB = NormResMu;
      M = L;
      bk = muk;
      muk = lambdak;
      NormResMu = NormResLambda;
      KeptVal = 1;
    }
    ///error = Norm(bk-ak)
    Interval = bk; Interval -= ak; error = Norm(Interval);

    ///alpha shall be alpha <= 1
    if(A > 1. && B > 1.) break;

  }///while

  ///procurando o menor valor de residuo entre A, B, Mu e lambda
  std::map<REAL,std::pair<REAL,TPZFMatrix<STATE>*> > mymap;
  mymap[NormResA] = std::make_pair(A,&ak);
  mymap[NormResB] = std::make_pair(B,&bk);
  mymap[NormResLambda] = std::make_pair(L,&lambdak);
  mymap[NormResMu] = std::make_pair(M,&muk);
  REAL ALPHA = mymap.begin()->second.first;
  NextAcc = *( mymap.begin()->second.second );

#ifdef DEBUGLINESEARCH
  ///debug: valor do alpha
  TPZFMatrix alpha;
  alpha = NextW;
  alpha -= Wn;
  REAL sum = 0.;
  int ncontrib = 0;
  for(int i = 0; i < alpha.Rows(); i++){
    if (DeltaW(i,0)){
      alpha(i,0) = alpha(i,0)/DeltaW(i,0);
      sum += alpha(i,0);
      ncontrib++;
    }
  }
  //REAL MeanAlpha = sum/ncontrib;
  alphafile << /*MeanAlpha << "\t" <<*/ "ALPHA = " << ALPHA << "\n";
  alphafile.flush();
#endif

  if(ALPHA > 1.){ ///alpha shall be alpha <= 1
    NextAcc = Acceleration;
    NextAcc += DeltaAccel;
#ifdef DEBUGLINESEARCH
 alphafile << "ALPHA LIMIT APPLIED. Alpha = 1.\n";
#endif
    return 1.;
  }

  return ALPHA;

}///void


