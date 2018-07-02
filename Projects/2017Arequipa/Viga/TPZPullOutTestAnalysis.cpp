//---------------------------------------------------------------------------

//#define _PRINT_DAMAGE_ //tamarindo

#include "TPZPullOutTestAnalysis.h"
//#include "TSWXConvertString.h"
#include "TSWXGraphMesh.h"
#include "TSWXGraphElement.h"
#include "pzelast3dGD.h"
//#include "System.hpp"
//#include "Dialogs.hpp"
#include "pzbndcond.h"
#include "pzintel.h"
#include "TPZInterfaceEl.h"
TPZPullOutTestAnalysis::TPZPullOutTestAnalysis():TPZNonLinearAnalysis(),fFullBCData(){

}

TPZPullOutTestAnalysis::TPZPullOutTestAnalysis(TPZCompMesh *mesh)
	: TPZNonLinearAnalysis(mesh, std::cout), fFullBCData() {

}

TPZPullOutTestAnalysis::~TPZPullOutTestAnalysis(){

}

void TPZPullOutTestAnalysis::TimeSteps(std::ostream &out, REAL tol,int numiter, bool linesearch, int nsubsteps,
                                       const std::string &filename){
#define FIBTEST
#ifdef PULLOUT
  std::ofstream reactionFile1("c:\\Temp\\Release\\reactionFile1.txt");
  std::ofstream reactionFile2("c:\\Temp\\Release\\reactionFile2.txt");
#endif
#ifdef FIBTEST
  std::ofstream reactionFile3("c:\\Temp\\FIBTest\\cargaFileFIBTest.txt");
  std::ofstream reactionFile4("c:\\Temp\\FIBTest\\displFileFIBTest.txt");
  std::ofstream reactionFile5("c:\\Temp\\FIBTest\\reactionFileFIBTestApoioMovel.txt");
  std::ofstream reactionFile7("c:\\Temp\\FIBTest\\reactionFileFIBTestApoioFixoVert.txt");
  std::ofstream reactionFile8("c:\\Temp\\FIBTest\\reactionFileFIBTestApoioFixoHoriz.txt");
  std::ofstream reactionFile6("c:\\Temp\\FIBTest\\displOpostoFileFIBTest.txt");
#endif

  {
    std::stringstream f;
    f << filename << "Start";
    this->PostProcess(0,0.,f.str().c_str());
  }

  TPZCompMesh * cmesh = this->Mesh();
  out << "NEquations = " << cmesh->NEquations() << "\n\n";
  out.flush();
  if(!cmesh) DebugStop();
  this->FillBCData();
  const REAL subStepSize = 1./(REAL)(nsubsteps);
  REAL lastStep = 0.;
  REAL subStep = 0.;
  int iStep = 0;
  while(lastStep < 1.){//enquanto for menor que todo o carregamento imposto (1.0)

    const TPZFMatrix<STATE> lastSol = this->Solution();
    this->LoadSolution();
    std::pair<int,REAL> convVals;

    //convergindo passo de carga
    subStep = lastStep + subStepSize;
    if(subStep > 1.){
      subStep = 1.;
    }

    out << "Starting istep = " << iStep << ",\t" << subStep << ",\t" << subStepSize << "\n";
      std::cout << "Starting istep = " << iStep << ",\t" << subStep << ",\t" << subStepSize << "\n";
    this->SetBCValue( subStep );
    int nMaxIter = -1;
    for(int iSubStep = 0; iSubStep < 30; iSubStep++){
      out << "\n\tStarting iSubStep = " << iSubStep << "\n";
      this->UpdateDamage(-1);
      convVals = this->SolveOneStep(iStep, out,tol,numiter,linesearch);
      nMaxIter = (nMaxIter > convVals.first) ? nMaxIter : convVals.first;//nMaxIter = max(nMaxIter,convVals.first)
      if(convVals.first == 1 && iSubStep >= 1){
        break;
      }
//#define LOGSUBSTEP
#ifdef LOGSUBSTEP
      {
        std::stringstream myfilenameSub;
        myfilenameSub << filename << "_debug_" << iStep << "_" << iSubStep;
        this->PostProcess(iStep,subStep,myfilenameSub.str());
      }
#endif
    }

    this->PostProcess(iStep,subStep,filename);
    const int Xdir = 0;
    const int Ydir = 1;
    const int Zdir = 2;
#ifdef PULLOUT
    this->SupportReaction(iStep,-1,Zdir,reactionFile1);
    this->SupportReaction(iStep,-2,Zdir,reactionFile2);
#endif
#ifdef FIBTEST
    bool simetriaFib = true;
    this->SupportReaction(iStep,+7,Ydir,reactionFile3);
    this->SupportReaction(iStep,-2,Ydir,reactionFile5);
    if(!simetriaFib){
      this->SupportReaction(iStep,-22,Ydir,reactionFile7);
      this->SupportReaction(iStep,-22,Xdir,reactionFile8);
      this->NodalDisplacement(iStep,+88,reactionFile6);
    }
    this->NodalDisplacement(iStep,+8,reactionFile4);

#endif
    this->UpdatePlasticDeformation();
    this->UpdateDamage(iStep);
    iStep++;
    lastStep = subStep;
    out << "done!\n\n";
    out.flush();

  }//tempo total de simulacao

}//void


void TPZPullOutTestAnalysis::NodalDisplacement(int istep, int nodeMatId, std::ostream &myfile){
  TPZGeoMesh * gmesh = this->Mesh()->Reference();
  for(int iel = 0; iel < gmesh->NElements(); iel++){
    TPZGeoEl * gelnode = gmesh->ElementVec()[iel];
    if(!gelnode) continue;
    if(gelnode->MaterialId() != nodeMatId) continue;
    TPZGeoElSide neigh = gelnode->Neighbour(0);
    if(neigh.Element()->Dimension() != 3) DebugStop();
    int side = neigh.Side();
    TPZGeoEl * gel = neigh.Element();
	TPZManVector<REAL> qsi(3);
	TPZSolVec sol(3);
	TPZGradSolVec dsol(9);
    TPZFNMatrix<9>  axes(3,3);
    gel->CenterPoint(side,qsi);
#ifdef DEBUG
{
    TPZManVector<REAL> x(3);
    gel->X(qsi,x);
    TPZManVector<REAL> nodeCoord(3);
    REAL diff = 0.;
    gelnode->NodePtr(0)->GetCoordinates(nodeCoord);
    for(int ii = 0; ii < 3; ii++) diff += pow(x[ii]-nodeCoord[ii],2);
    if(sqrt(diff) > 1e-6) DebugStop();
}
#endif
    TPZInterpolationSpace * sp = dynamic_cast< TPZInterpolationSpace* >(gel->Reference());
    if(!sp) DebugStop();
    sp->ComputeSolution(qsi, sol, dsol, axes);
    for(int iv = 0; iv < sol.NElements(); iv++){
      myfile << sol[iv] << "\t";
    }
    myfile << "\n";
    myfile.flush();
  }
}

void TPZPullOutTestAnalysis::SupportReaction(int istep, int supportMatId, int direction, std::ostream &myfile){
  TPZCompMesh * cmesh = this->Mesh();
  std::set<int> matids;
  std::map<int , TPZMaterial* >::const_iterator w;
  for(w = cmesh->MaterialVec().begin(); w != cmesh->MaterialVec().end(); w++){
    int id = w->first;
    if(id != supportMatId){
      matids.insert(id);
    }
  }
  this->StructMatrix()->SetMaterialIds(matids);
  this->AssembleResidual();

  //collecting connects
  std::set<int> eqset;
  for(int iel = 0; iel < cmesh->NElements(); iel++){
    TPZCompEl * cel = cmesh->ElementVec()[iel];
    if(!cel) continue;
    if(cel->Reference()->MaterialId() != supportMatId) continue;
    TPZInterpolatedElement * intel = dynamic_cast<TPZInterpolatedElement *>(cel);
    if(intel){
      for(int ic = 0; ic < intel->NCornerConnects(); ic++){
        int connectindex = intel->ConnectIndex(ic);
        TPZConnect &np = cmesh->ConnectVec()[connectindex];
        int blocknumber = np.SequenceNumber();
        int firsteq = cmesh->Block().Position(blocknumber);
        int ndf = cmesh->Block().Size(blocknumber);
        if(ndf != 3) DebugStop();
        eqset.insert(firsteq+direction);
      }
    }//if intel
    TPZInterfaceElement * face = dynamic_cast<TPZInterfaceElement *>(cel);
    if(face){
    //  TPZBndCond * bc = dynamic_cast< TPZBndCond * >(face->Material().operator->());
		TPZBndCond * bc = dynamic_cast< TPZBndCond * >(face->Material());
      if(!bc) DebugStop();//se eh support, tem que ser bndcond
      TPZCompElDisc * disc = dynamic_cast<TPZCompElDisc*>(face->LeftElement());
      if(!disc) DebugStop();//abaixo, tem que ser disc, nao pode ser intel
      int connectindex = disc->ConnectIndex(0);
      TPZConnect &np = cmesh->ConnectVec()[connectindex];
      int blocknumber = np.SequenceNumber();
      int firsteq = cmesh->Block().Position(blocknumber);
      int ndf = cmesh->Block().Size(blocknumber);
      int nshape = disc->NShapeF();
      //assumindo que a ultima shape eh a constante = 1, o rhs dela eh o valor do fluxo
      int equationindex = firsteq + 3 * (nshape-1) + direction;
      eqset.insert(equationindex);
    }//if face
  }

  //integrating
  std::pair<int,REAL> integral(0,0.0);
  for(std::set<int>::iterator w = eqset.begin(); w != eqset.end(); w++){
    const int index = *w;
    integral.second += fRhs(index,0);
    integral.first++;
  }

  if(integral.first == 0){
    DebugStop();
  }
  myfile << istep << "\t"
         << integral.first << "\t" << integral.second << "\n";
  myfile.flush();

  std::set<int> fakeEmpty;
  this->StructMatrix()->SetMaterialIds(fakeEmpty);

}

void TPZPullOutTestAnalysis::PostProcess(int istep, REAL subStep, const std::string &filename){
  this->PostProcessAverageValues();
  TSWXGraphMesh graphMesh;
  int resolution = 0;
  TSWXGraphElement graphEl(resolution);
  const int nvars = 23;
  TPZVec<std::string> nodalVarIndex(nvars), cellVarIndex(nvars);
  std::string varnames[nvars] = {"Displacement","DisplacementX","DisplacementY","DisplacementZ",
                                 "StressX","StressY","StressZ","PrincipalStress",
                                 "MatId","yielding","Damage","PrincipalStrain",
                                 "PrincipalPlasticStrain","StrainZ","PlasticStrainZ",
                                 "POrder","PlasticFunctionVal","ElementSize",
                                 "GeoElIndex",
                                 "AvDamage", "AvStressX", "AvStressY", "AvStressZ"};
  for(int i = 0; i < nvars; i++){
    nodalVarIndex[i] = varnames[i];
  }
  cellVarIndex = nodalVarIndex;
  graphEl.GenerateVTKData(this->Mesh(),this->Mesh()->Dimension(),subStep, nodalVarIndex, cellVarIndex, graphMesh);
  std::stringstream name;
  name << filename << istep << ".vtk";
  std::ofstream paraviewfile(name.str().c_str());
  paraviewfile.precision(8);
  graphMesh.ToParaview(paraviewfile);
  paraviewfile.close();
}

void TPZPullOutTestAnalysis::PostProcessAverageValues(){
  const int nel = this->Mesh()->NElements();
  const int nvars = 4;
  const std::string varnames[nvars] = { "Damage", "StressX", "StressY", "StressZ" };//que se tornarao "AvDamage", "AvStressX", "AvStressY", "AvStressZ"};
  TPZVec< TPZVec< REAL > > averageSols(nvars);
  for(int i = 0; i < nvars; i++){
    averageSols[i].Resize(nel);
    averageSols[i].Fill(-8511965.);
  }
  TPZManVector<REAL,3> local(1);
  for(int iel = 0; iel < nel; iel++){
    TPZInterpolationSpace * sp = dynamic_cast< TPZInterpolationSpace* >(this->Mesh()->ElementVec()[iel]);
    if(!sp) continue;
    if(sp->Dimension() != this->Mesh()->Dimension()){
      continue;
    }
    for(int iv = 0; iv < nvars; iv++){
      const int varindex = sp->Material()->VariableIndex( varnames[iv] );
//      REAL domain = 0.;
      sp->Integrate(varindex, local);
      if(local.NElements() != 1) DebugStop();
	  REAL value = local[0];   // / domain;
      averageSols[ iv ][ iel ] = value;
    }
  }//for iel

  for(int iv = 0; iv < nvars; iv++){
    this->Mesh()->SetElementSolution(iv,averageSols[iv]);
  }

}///void

#include "pzsolve.h"
#include "pzstepsolver.h"

REAL f(REAL val){
  if(fabs(val)<1e-10) return 0.;
  else return val;
}

std::pair<int,REAL> TPZPullOutTestAnalysis::SolveOneStep(int TimeStep, std::ostream &out,REAL tol,int numiter, bool linesearch) {

#ifdef DEBUG_LOG
  std::ofstream solFile("c:\\Temp\\solFile.txt");
  std::ofstream dsolFile("c:\\Temp\\dsolFile.txt");
  std::ofstream rhsFile("c:\\Temp\\rhsFile.txt");
#endif

   int iter = 0;
   int numeq = fCompMesh->NEquations();
   TPZFMatrix<STATE> prevsol(fSolution);
   if(prevsol.Rows() != numeq) prevsol.Redim(numeq,1);
   REAL normDeltaSol = 2*tol+1, NormResLambda = 2*tol+1;
   REAL NormResLambdaInfinity = 2*tol+1;
   while( (NormResLambdaInfinity > tol /*normDeltaSol > tol*/ /*|| NormResLambda > 1.*/ /*10.*/)
       && (iter < numiter)) {

      fSolution.Redim(0,0);
//      if(this->Solver().Matrix() == NULL || (iter >= 1)  // coesiva apenas, as vezes
      if(this->Solver().Matrix() == NULL || (iter == 1) ){ //quando usar a tangente elastica, assim basta
        out << "Assembling Jacobian matrix..."; out.flush();
        Assemble();
        out << "done!\n"; out.flush();
      }
      else{
        AssembleResidual();
      }

      Solve();

#ifdef DEBUG_LOG
{
      dsolFile << "\niter " << iter << std::endl;
      for(int inode = 0; inode < fSolution.Rows()/3; inode++){
        for(int ivar = 0; ivar < 3; ivar++) dsolFile << f(fSolution( 3*inode+ivar, 0 )) << "\t\t";
        dsolFile << "\n";
      }
      dsolFile.flush();
}
{
      rhsFile << "\niter " << iter << std::endl;
      for(int inode = 0; inode < fRhs.Rows()/3; inode++){
        for(int ivar = 0; ivar < 3; ivar++) rhsFile << f(fRhs( 3*inode+ivar, 0 )) << "\t\t";
        rhsFile << "\n";
      }
      rhsFile.flush();
}
#endif


      bool perturbacao = (iter % 100 == 0) ? 1 : 0;
      REAL alphaLineSearch = 1.;
      if (perturbacao == false && linesearch && iter > 2){
		  TPZFMatrix<STATE> nextSol;
        REAL LineSearchTol = 1e-3 * Norm(fSolution);
        int niter = 2; //aqui, nao adianta serem muitos passos. Torna a execucao lenta e muda pouco na convergencia
        alphaLineSearch = this->LineSearch(prevsol, fSolution, LineSearchTol, niter, nextSol);
        fSolution = nextSol;
      }
      else{
        fSolution += prevsol;
      }

      prevsol -= fSolution;
      normDeltaSol = Norm(prevsol);
      prevsol = fSolution;
      this->LoadSolution(fSolution);
      this->AssembleResidual();
      NormResLambda = this->RhsNorm(fRhs,true,false);
      NormResLambdaInfinity = this->RhsNorm(fRhs,false,true); //tamarindo this->RhsNorm(fRhs,true,true);

#ifdef DEBUG_LOG
{
      solFile << "\niter " << iter << std::endl;
      for(int inode = 0; inode < fSolution.Rows()/3; inode++){
        for(int ivar = 0; ivar < 3; ivar++) solFile << f(fSolution( 3*inode+ivar, 0 )) << "\t\t";
        solFile << "\n";
      }
      solFile.flush();
}
#endif

      if(linesearch){
        out << "Iteracao n : " << (iter+1) << ", alphaLS = " << alphaLineSearch
            << " : normas |Delta(Un)| e |Delta(rhs)| : " << normDeltaSol << " / "
            << NormResLambda << " / " << NormResLambdaInfinity << std::endl;
      }
      else{
        out << "Iteracao n : " << (iter+1) << " : normas |Delta(Un)| e |Delta(rhs)| : "
            << normDeltaSol << " / " << NormResLambda << " / " << NormResLambdaInfinity << std::endl;
      }

	    iter++;

	    out.flush();
   }

   if(iter >= numiter){
     out << "!!!!!! Não convergi passo " << TimeStep << " !!!!!!!!!\n";
     out.flush();
//    UnicodeString mess = "Não consegui convergir passo ";
//    mess += TimeStep;
//    ShowMessage(mess); tamarindo
   }

   return std::pair<int,REAL>(iter,NormResLambdaInfinity);
}


void TPZPullOutTestAnalysis::FillBCData(){
  TPZCompMesh * cmesh = this->Mesh();
  if(!cmesh) DebugStop();
  fFullBCData.clear();
  std::map<int , TPZMaterial* >::const_iterator w = cmesh->MaterialVec().begin();
  std::map<int , TPZMaterial* >::const_iterator e = cmesh->MaterialVec().end();
  for(; w != e; w++){
    TPZMaterial * mat = w->second;
    TPZBndCond *bc = dynamic_cast<TPZBndCond*>(mat);
    if(!bc) continue;
    fFullBCData[ bc ] = bc->Val2();
  }//for
}//void

void TPZPullOutTestAnalysis::SetBCValue(REAL subStep){
	std::map<TPZBndCond*, TPZFMatrix<STATE> >::iterator w;
  TPZBndCond *bc = NULL;
  for(w = fFullBCData.begin(); w != fFullBCData.end(); w++){
    bc = w->first;
    bc->Val2() = w->second;
    bc->Val2() *= subStep;
  }
}

//#define DEBUGLINESEARCH
#ifdef DEBUGLINESEARCH
static ofstream alphafile("c:\\Temp\\tmp\\alpha.txt");
#endif
REAL TPZPullOutTestAnalysis::LineSearch(const TPZFMatrix<STATE> &Wn,
	const TPZFMatrix<STATE> DeltaW, REAL tol, int niter, TPZFMatrix<STATE> &NextW) {

  const bool removeBigNumbers = false;

  REAL error = 2.*tol+1.;
  REAL A, B, L, M;
  TPZFMatrix<STATE> ak, bk, lambdak, muk, Interval;
  REAL NormResLambda, NormResMu;
  ///ak = Wn + 0.1 * DeltaW
  ak = DeltaW;
  A = 0.1;
  ak *= A;
  ak += Wn;
  ///bk = Wn + DeltaW
  bk = DeltaW;
//  B = 1.6181229773462784*(1.-0.382*A); //este valor eh exatamente pra que M seja igual a 1
  B = 2.6178010471204187*(1.-0.618*A); //este valor eh exatamente pra que L seja igual a 1
  bk *= B;
  bk += Wn;

  ///calculando residuo em A e B
  REAL NormResA, NormResB;
  {
    ///computing residual
    this->LoadSolution(ak);
    this->AssembleResidual();
    std::cout.flush();
    NormResA = RhsNorm(fRhs,removeBigNumbers,false);
  }
  {
    ///computing residual
    this->LoadSolution(bk);
    this->AssembleResidual();
    std::cout.flush();
    NormResB = RhsNorm(fRhs,removeBigNumbers,false);
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
      this->LoadSolution(lambdak);
      this->AssembleResidual();
      NormResLambda = RhsNorm(fRhs,removeBigNumbers,false);
    }

    if (KeptVal != 1){
      ///muk = ak + 0.618*(bk-ak)
      M = 0.618*(B-A)+A;
      muk = Interval; muk *= 0.618; muk += ak;
      this->LoadSolution(muk);
      this->AssembleResidual();
      NormResMu = RhsNorm(fRhs,removeBigNumbers,false);
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

    ///alpha shall be alpha <= 1 in base class. Not here.
  //  if(A > 1. && B > 1.) break;

  }///while

  ///procurando o menor valor de residuo entre A, B, Mu e lambda
  std::map<REAL, std::pair<REAL, TPZFMatrix<STATE>* > > mymap;
  mymap[NormResA] = std::make_pair(A,&ak);
  mymap[NormResB] = std::make_pair(B,&bk);
  mymap[NormResLambda] = std::make_pair(L,&lambdak);
  mymap[NormResMu] = std::make_pair(M,&muk);
  REAL ALPHA = mymap.begin()->second.first;
  NextW = *( mymap.begin()->second.second );

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

//In base class, but not here
/*  if(ALPHA > 1.){ ///alpha shall be alpha <= 1
    NextW = Wn;
    NextW += DeltaW;
#ifdef DEBUGLINESEARCH
 alphafile << "ALPHA LIMIT APPLIED. Alpha = 1.\n";
#endif
    return 1.;
  }         */

  return ALPHA;

}///void

REAL TPZPullOutTestAnalysis::LineSearch(const TPZFMatrix<STATE> &Wn,
	const TPZFMatrix<STATE> DeltaW,	TPZFMatrix<STATE> &NextW){

  const int ntrials = 10;
  const REAL table[ntrials] = {0.1, 0.25, 0.5, 0.75, 1., 1.5, 2., 5., 10.,20.};

  std::map< REAL, int > mymap;
  TPZFMatrix<STATE> ak;
  for(int i = 0; i < ntrials; i++){

    ///ak = Wn + A * DeltaW
    ak = DeltaW;
    const REAL A = table[i];
    ak *= A;
    ak += Wn;

    ///computing residual
    this->LoadSolution(ak);
    this->AssembleResidual();
    std::cout.flush();
    const REAL NormResA = RhsNorm(fRhs,false,false);

    mymap[ NormResA ] = i;

  }//for i

  const int bestI = mymap.begin()->second;
  const REAL ALPHA = table[ bestI ];
  NextW = DeltaW;
  NextW *= ALPHA;
  NextW += Wn;

  return ALPHA;

}///void

void TPZPullOutTestAnalysis::UpdatePlasticDeformation(){
  TPZCompMesh *cmesh = this->Mesh();
  if(!cmesh) DebugStop();
  std::map<int ,TPZMaterial* >::const_iterator w = cmesh->MaterialVec().begin();
  std::map<int ,TPZMaterial* >::const_iterator e = cmesh->MaterialVec().end();
  for(; w != e; w++){
    TPZMaterial * mat = w->second;
    TPZElasticity3DGD *elast = dynamic_cast<TPZElasticity3DGD*>(mat);
    if(!elast) continue;
//    if(elast->PlasticModel()) elast->PlasticModel()->UpdatePlasticDeformation();
  }//for
}//void

void TPZPullOutTestAnalysis::UpdateDamage(int istep){
  TPZCompMesh *cmesh = this->Mesh();
  if(!cmesh) DebugStop();
  std::map<int ,TPZMaterial* >::const_iterator w = cmesh->MaterialVec().begin();
  std::map<int ,TPZMaterial* >::const_iterator e = cmesh->MaterialVec().end();

#ifdef _PRINT_DAMAGE_
  std::ofstream myfile;
  if(istep >= 0){
    UnicodeString filename;
    filename = "C:\\Temp\\fHistoricalDamage";
    filename += istep;
    filename += ".txt";
    myfile.open(swx::Unicode2string(filename).c_str());
  }
#endif
  for(; w != e; w++){
    TPZMaterial * mat = w->second;
    TPZElasticity3DGD *elast = dynamic_cast<TPZElasticity3DGD*>(mat);
    if(!elast) continue;
//    elast->UpdateDamage();
#ifdef _PRINT_DAMAGE_
    if(istep >= 0){
      elast->DamageModel()->fHistoricalDamage.Print(myfile);
    }
#endif
  }//for

#ifdef _PRINT_DAMAGE_
  if(istep >= 0){
    myfile.close();
  }
#endif

}//void

void TPZPullOutTestAnalysis::ForgetDamageStep(){
  TPZCompMesh *cmesh = this->Mesh();
  if(!cmesh) DebugStop();
  std::map<int ,TPZMaterial* >::const_iterator w = cmesh->MaterialVec().begin();
  std::map<int ,TPZMaterial* >::const_iterator e = cmesh->MaterialVec().end();
  for(; w != e; w++){
    TPZMaterial * mat = w->second;
    TPZElasticity3DGD *elast = dynamic_cast<TPZElasticity3DGD*>(mat);
    if(!elast) continue;
  }//for
}//void

void GetEquationIndices(TPZCompEl* cel,TPZStack<int64_t> &eqs) {
	eqs.Resize(0);
	const int64_t ncon = cel->NConnects();
	TPZBlock<STATE> &block = cel->Mesh()->Block();
	for (int64_t in = 0; in<ncon; in++) {
		TPZConnect *df = &cel->Connect(in);
		int64_t dfseq = df->SequenceNumber();
		int dfvar = block.Size(dfseq);
		int pos = block.Position(dfseq);
		for (int jn = 0; jn<dfvar; jn++){
			const int64_t eqindex = pos + jn;
			eqs.Push(eqindex);
		}
	}
#ifdef DEBUG
	if (eqs.NElements() != cel->NEquations()){
		DebugStop();
	}
#endif
}

REAL TPZPullOutTestAnalysis::RhsNorm(const TPZFMatrix<STATE> &origrhs,
                                     bool removeBigNumbers,
                                     bool infinity) const{

  if(origrhs.Cols() != 1) DebugStop();//assumi que era uma unica coluna

//  removeBigNumbers = false;//tamarindo

  TPZFMatrix<STATE> rhs = origrhs;

  if(removeBigNumbers){
    TPZCompMesh *cmesh = this->Mesh();
    const int nel = cmesh->NElements();
    TPZStack<int64_t> eqs;
    for(int iel = 0; iel < nel; iel++){
      TPZCompEl * cel = cmesh->ElementVec()[iel];
      if(!cel) continue;
      TPZInterfaceElement * face = dynamic_cast<TPZInterfaceElement*>(cel);
      if(face) continue;//o saco eh que a interface retorna todas as equacoes do elemento descontinuo e nao apenas as da face restrita
      TPZBndCond * bnd = dynamic_cast< TPZBndCond * >(cel->Material());
      if(!bnd) continue;
      const int bcType = bnd->Type();
	  TPZElasticity3DGD* mat = (TPZElasticity3DGD *)(bnd->Material());
      if(!mat) DebugStop();
      bool usesbig = mat->UsesBigNumberForBC(bcType);
      if(usesbig == false) continue;
      eqs.Resize(0);
      GetEquationIndices(cel,eqs);
      for(int ieq = 0; ieq < eqs.NElements(); ieq++){
        const int index = eqs[ieq];
        rhs(index,0) = 0.;
      }
    }//iel
  }

  if(infinity){
    REAL max = 0.;
    int n = rhs.Rows();
    for(int i = 0; i < n; i++){
      const REAL val = fabs(rhs(i,0));
      max = max > val ? max : val;
    }
    return max;
  }
  else return Norm(rhs);
}

