#include "pzlog.h"
#include "TPZFracAnalysis.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "TPZVTKGeoMesh.h"
#include "TPZMatfrac1dhdiv.h"
#include "pzbndcond.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzanalysis.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstepsolver.h"
#include "pzmatred.h"

#ifdef PZ_LOG
static PZLogger logger("pz.frac");
#endif

TPZFracAnalysis::TPZFracAnalysis(TPZAutoPointer<TPZFracData> Data)
{
  fData = Data;
  fmeshvec.Resize(2);
  fgmesh = NULL;
  fcmeshMixed = NULL;
  for (int i = 0; i < 2; i++) {
    fmeshvec[i] = NULL;
  }
  fLastStepRhs.Redim(0, 0);
  fMatFrac = NULL;
  fmustStop = false;
}

TPZFracAnalysis::~TPZFracAnalysis()
{
  for (int i = 0; i < 2; i++) {
    delete fmeshvec[i];
  }
  delete fcmeshMixed;
  delete fgmesh;
  /*
   if (fMatFrac) {
   delete fMatFrac;
   }
   */
}

void TPZFracAnalysis::Run()
{
  
  // Creating empty gmesh
  fgmesh = new TPZGeoMesh;
  
  REAL vl = this->RunUntilOpen();
  
  const REAL lFrac = fData->ElSize();
  fData->SetLfrac(lFrac);
  
  // Criando primeira malha geometrica com um elemento
  this->CreateFirstGeoElWithBC();
  
  // Malhas computacionais - FluxoH1, PressaoL2, Multifisica para sistema misto quente
  fmeshvec[0] = CreateCMeshFluxH1();
  fmeshvec[1] = CreateCMeshPressureL2();
  TPZFMatrix<REAL> vlMatrix(1,1,vl);
  fcmeshMixed = CreateCMeshMixed(vlMatrix);
  
  // Analysis
  bool mustOptimizeBandwidth = false;
  TPZAnalysis *an = new TPZAnalysis(fcmeshMixed,mustOptimizeBandwidth);
  TPZSkylineNSymStructMatrix skyl(fcmeshMixed);
  TPZStepSolver<STATE> step;
  step.SetDirect(ELU);
  an->SetSolver(step);
  an->SetStructuralMatrix(skyl);
  
  // Plot of first solution
  this->PostProcessVTK(an);
  
  // Solving transiente system
  while (fmustStop == false) {
    bool propagate = SolveSistTransient(an);
    
    
    
    
    
    if (propagate)
    {
      
      // Novo comprimento de fratura
      REAL newLfrac = fData->Lfrac() + fData->ElSize();
      fData->SetLfrac(newLfrac);
      std::cout << "Lfrac = " << newLfrac << std::endl;
      
      // Apagando elemento de contorno que impoe pressao
      TPZGeoEl *gel = this->FindPressureBCElement();
      fgmesh->ResetReference();
      fmeshvec[0]->LoadReferences();
      delete gel->Reference();
      fgmesh->ResetReference();
      fmeshvec[1]->LoadReferences();
      delete gel->Reference();
      fgmesh->ResetReference();
      fcmeshMixed->LoadReferences();
      gel->Reference()->SetFreeIntPtIndices();
      delete gel->Reference();
      gel->RemoveConnectivities();
      delete gel;
      
      fmeshvec[0]->CleanUpUnconnectedNodes();
      fmeshvec[1]->CleanUpUnconnectedNodes();
      fcmeshMixed->CleanUpUnconnectedNodes();
      
      // Criando novo noh
      TPZAdmChunkVector<TPZGeoNode> &nodevec = fgmesh->NodeVec();
      const int nnodes = nodevec.NElements() + 1;
      nodevec.Resize(nnodes);
      const int lastnodeid = nnodes - 1;
      nodevec[lastnodeid].SetNodeId(lastnodeid);
      TPZVec<REAL> x(3,0.);
      x[0] = newLfrac;
      nodevec[lastnodeid].SetCoord(x);
      
      // Criando novo elemento geometrico
      TPZVec<int64_t> TopolLine(2);
      TopolLine[0] = nnodes - 2;
      TopolLine[1] = nnodes - 1;
      int64_t index;
      const int matid = 1;
      TPZGeoEl *newGel = fgmesh->CreateGeoElement(EOned, TopolLine, matid, index);
      fgmesh->BuildConnectivity();
      
      // Colocano a cond de contorno de pressao na malha geometrica
      const int bcpressureid = -2;
      TPZGeoEl *gelBCPress = newGel->CreateBCGeoEl(1, bcpressureid);
      fgmesh->BuildConnectivity();
      
      // Criando novo elemento computacional de fratura
      fgmesh->ResetReference();
      fmeshvec[0]->LoadReferences();
      fmeshvec[0]->CreateCompEl(newGel, index);
      fgmesh->ResetReference();
      fmeshvec[1]->LoadReferences();
      TPZCompEl *celPress = fmeshvec[1]->CreateCompEl(newGel, index);
      fgmesh->ResetReference();
      fcmeshMixed->LoadReferences();
      TPZCompEl *celMixed = fcmeshMixed->CreateCompEl(newGel, index);
      
      // Criando novo elemento computacional de BC de pressao
      fgmesh->ResetReference();
      fmeshvec[0]->LoadReferences();
      fmeshvec[0]->CreateCompEl(gelBCPress, index);
      fgmesh->ResetReference();
      fmeshvec[1]->LoadReferences();
      fmeshvec[1]->CreateCompEl(gelBCPress, index);
      fgmesh->ResetReference();
      fcmeshMixed->LoadReferences();
      TPZCompEl *celBCMixed = fcmeshMixed->CreateCompEl(gelBCPress, index);
      
      // Ajustando a estrutura de dados
      fmeshvec[0]->ComputeNodElCon();
      fmeshvec[1]->ComputeNodElCon();
      fcmeshMixed->ComputeNodElCon();
      fmeshvec[0]->CleanUpUnconnectedNodes();
      fmeshvec[1]->CleanUpUnconnectedNodes();
      fcmeshMixed->CleanUpUnconnectedNodes();
      fmeshvec[0]->ExpandSolution();
      fmeshvec[1]->ExpandSolution();
      fcmeshMixed->ExpandSolution();
      
      // Setando a pressao inicial no elemento novo
      TPZConnect &c = celPress->Connect(0);
      TPZBlock<STATE> &Block = fmeshvec[1]->Block();
      int sq = c.SequenceNumber();
      int sz = Block.Size(sq);
      if (sz > 1) {
        DebugStop(); // Only works for p = 0 in pressure space
      }
      int pos = Block.Position(sq);
      fmeshvec[1]->Solution()(pos,0) = fData->SigmaConf();
      
      // Transferindo para a multifisica
      TPZBuildMultiphysicsMesh::AddElements(fmeshvec, fcmeshMixed);
      TPZBuildMultiphysicsMesh::AddConnects(fmeshvec, fcmeshMixed);
      TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, fcmeshMixed);
      
#ifdef PZDEBUG
      std::ofstream out2("meshes2.txt");
      fgmesh->Print(out2);
      fmeshvec[0]->Print(out2);
      fmeshvec[1]->Print(out2);
      fcmeshMixed->Print(out2);
#endif
      
      // Preparando os index dos pontos de integracao.
      TPZFMatrix<> Vl(1,1,fData->AccumVl());
      fMatFrac->SetDefaultMem(Vl);
      celMixed->PrepareIntPtIndices();
      celBCMixed->PrepareIntPtIndices();
      // Vl eh resetado depois de inicializar o chute inicial de newton
      
#ifdef PZDEBUG
      std::ofstream outvtk("newfrac.vtk");
      TPZVTKGeoMesh::PrintGMeshVTK(fgmesh, outvtk,true);
#endif
    }
    
    TPZEquationFilter newEquationFilter(fcmeshMixed->NEquations());
    an->StructMatrix()->EquationFilter() = newEquationFilter; //AQUINATHAN
    
    //std::cout << "\nSolucao apos a propagacao:" << std::endl;
    //fcmeshMixed->Solution().Print("sol");
    
    an->LoadSolution(fcmeshMixed->Solution());
    this->PostProcessVTK(an);
  }
  
  fData->PrintDebugMapForMathematica("LfracO.nb");
  
  delete an;
}

void TPZFracAnalysis::RunTest()
{
  
  DebugStop();
  
  // Malha geometrica
  fgmesh = CreateGMesh();
  
  // Malhas computacionais - FluxoH1, PressaoL2, Multifisica para sistema misto quente
  fmeshvec[0] = CreateCMeshFluxH1();
  fmeshvec[1] = CreateCMeshPressureL2();
  TPZFMatrix<REAL> vlZero(1,1,0.);
  fcmeshMixed = CreateCMeshMixed(vlZero);
  
  // Analysis
  bool mustOptimizeBandwidth = false;
  TPZAnalysis *an = new TPZAnalysis(fcmeshMixed,mustOptimizeBandwidth);
  TPZSkylineNSymStructMatrix skyl(fcmeshMixed);
  TPZStepSolver<STATE> step;
  step.SetDirect(ELU);
  an->SetSolver(step);
  an->SetStructuralMatrix(skyl);
  
  // Plot of first solution
  this->PostProcessVTK(an);
  
  // Solving transiente system
  SolveSistTransient(an);
  
  //fData->PrintDebugMapForMathematica("DebugMapQl.nb");
  
  delete an;
}

TPZGeoMesh * TPZFracAnalysis::CreateGMesh()
{
  const int nel = fData->NelFrac();
  const REAL lfrac = fData->Lfrac();
  const int matIdFrac = 1, bcLeftId = -1, bcRightId = -2;
  TPZGeoMesh *gmesh = new TPZGeoMesh;
  
  // Nos
  const int nnodes = nel+1;
  const REAL elsize = lfrac / nel;
  gmesh->NodeVec().Resize(nnodes);
  TPZVec<REAL> coord(3,0.);
  for (int in = 0; in < nnodes; in++) {
    coord[0] = in * elsize;
    gmesh->NodeVec()[in].SetNodeId(in);
    gmesh->NodeVec()[in].SetCoord(coord);
  }
  
  // Elementos
  TPZVec<int64_t> TopolLine(2,0);
  int64_t index = 0;
  for (int iel = 0; iel < nel; iel++) {
    TopolLine[0] = iel;
    TopolLine[1] = iel+1;
    gmesh->CreateGeoElement(EOned, TopolLine, matIdFrac, index);
  }
  
  gmesh->BuildConnectivity();
  
  // Left
  TPZVec<int64_t> TopolPoint(1,0);
  gmesh->CreateGeoElement(EPoint, TopolPoint, bcLeftId, index);
  
  // Right
  TopolPoint[0] = nnodes-1;
  gmesh->CreateGeoElement(EPoint,TopolPoint,bcRightId,index);
  
  gmesh->BuildConnectivity();
  
#ifdef PZDEBUG
  std::ofstream out("GeoMesh.vtk");
  TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
#endif
  
  return gmesh;
}

TPZCompMesh * TPZFracAnalysis::CreateCMeshFluxH1()
{
  // Definicao de ids e tipos
  const int matIdFrac = 1, bcLeftId = -1, bcRightId = -2;
  const int typeFlux = 0, typePressure = 1;
  const int fluxorder = fData->PorderFlow();
  TPZFMatrix<STATE> val1(1,1,0.), val2(1,1,0.);
  
  // Malha computacional
  TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
  
  // Material da fratura
  TPZMatfrac1dhdiv *mat = new TPZMatfrac1dhdiv(matIdFrac);
  cmesh->InsertMaterialObject(mat);
  
  // Condicao de contorno na esquerda
  TPZBndCond * bcLeft = mat->CreateBC(mat, bcLeftId, typeFlux, val1, val2);
 	cmesh->InsertMaterialObject(bcLeft);
  
  // Condicao de contorno na direita
  TPZBndCond * bcRight = mat->CreateBC(mat, bcRightId, typePressure, val1, val2);
 	cmesh->InsertMaterialObject(bcRight);
  
  // Setando H1
  cmesh->SetDimModel(1);
  cmesh->SetDefaultOrder(fluxorder);
  cmesh->SetAllCreateFunctionsContinuous();
  cmesh->AutoBuild();
  
#ifdef PZDEBUG
  std::ofstream out("CMeshFluxH1.txt");
  cmesh->Print(out);
#endif
  
  return cmesh;
}

TPZCompMesh * TPZFracAnalysis::CreateCMeshPressureL2()
{
  // Definicao de ids e tipos
  const int matIdFrac = 1, bcLeftId = -1, bcRightId = -2;
  const int typeFlux = 0, typePressure = 1;
  const int pressureorder = fData->PorderPressure();
  TPZFMatrix<STATE> val1(1,1,0.), val2(1,1,0.);
  
  // Malha computacional
  TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
  
  // Material da fratura
  TPZMatfrac1dhdiv *mat = new TPZMatfrac1dhdiv(matIdFrac);
  cmesh->InsertMaterialObject(mat);
  
  // Condicao de contorno na esquerda
  TPZBndCond * bcLeft = mat->CreateBC(mat, bcLeftId, typeFlux, val1, val2);
 	cmesh->InsertMaterialObject(bcLeft);
  
  
  // Condicao de contorno na direita
  TPZBndCond * bcRight = mat->CreateBC(mat, bcRightId, typePressure, val1, val2);
 	cmesh->InsertMaterialObject(bcRight);
  
  // Setando L2
  cmesh->SetDimModel(1);
  cmesh->SetDefaultOrder(pressureorder);
  cmesh->SetAllCreateFunctionsDiscontinuous();
  
  cmesh->AutoBuild();
  
  //  int ncon = cmesh->NConnects();
  //  for(int i=0; i<ncon; i++)
  //  {
  //    TPZConnect &newnod = cmesh->ConnectVec()[i];
  //    newnod.SetLagrangeMultiplier(1);
  //  }
  
  for (int i = 0; i < cmesh->Solution().Rows(); i++) {
    cmesh->Solution()(i,0) = fData->SigmaConf();
  }
  
#ifdef PZDEBUG
  std::ofstream out("CMeshPressureL2.txt");
  cmesh->Print(out);
#endif
  
  return cmesh;
}

TPZCompMesh * TPZFracAnalysis::CreateCMeshMixed(TPZFMatrix<REAL> vlMatrix)
{
  // Definicao de ids e tipos
  const int matIdFrac = 1, bcLeftId = -1, bcRightId = -2;
  const int typeFlux = 0, typePressure = 1;
  TPZFMatrix<STATE> val1(1,1,0.), val2(1,1,0.);
  
  // Malha computacional
  TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
  
  // Material da fratura
  fMatFrac = new TPZMatfrac1dhdiv(matIdFrac);
  fMatFrac->SetSimulationData(fData);
  cmesh->InsertMaterialObject(fMatFrac);
  fMatFrac->SetDefaultMem(vlMatrix);
  
  // Condicao de contorno na esquerda
  val2(0,0) = fData->Q();
  TPZBndCond * bcLeft = fMatFrac->CreateBC(fMatFrac, bcLeftId, typeFlux, val1, val2);
  cmesh->InsertMaterialObject(bcLeft);
  
  // Condicao de contorno na direita
  val2(0,0) = fData->SigmaConf();
  TPZBndCond * bcRight = fMatFrac->CreateBC(fMatFrac, bcRightId, typePressure, val1, val2);
  cmesh->InsertMaterialObject(bcRight);
  
  // Setando Multifisico
  cmesh->SetDimModel(1);
  cmesh->SetAllCreateFunctionsMultiphysicElemWithMem();
  cmesh->AutoBuild();
  
  // Transferindo para a multifisica
  TPZBuildMultiphysicsMesh::AddElements(fmeshvec, cmesh);
  TPZBuildMultiphysicsMesh::AddConnects(fmeshvec, cmesh);
  TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, cmesh);
  
  // Preparando os index dos pontos de integracao.
  int64_t nel = cmesh->NElements();
  for (int64_t iel = 0; iel < nel; iel++) {
    TPZCompEl *cel = cmesh->ElementVec()[iel];
    cel->PrepareIntPtIndices();
  }
  
#ifdef PZDEBUG
  std::ofstream out("CMeshMultiPhysics.txt");
  cmesh->Print(out);
#endif
  
  return cmesh;
}

void TPZFracAnalysis::IterativeProcess(TPZAnalysis *an, std::ostream &out, int numiter)
{
  int iter = 0;
  REAL error = 1.e10, NormResLambdaLast = 1.e10;
  const REAL tol = 1.e-8;
  
  fData->SetCurrentState();
  int numeq = an->Mesh()->NEquations();
  
  TPZFMatrix<STATE> prevsol(an->Solution());
  TPZFMatrix<STATE> SoliterK(prevsol);
  if(prevsol.Rows() != numeq) prevsol.Redim(numeq,1);
  
  an->Assemble();
  an->Rhs() += fLastStepRhs;
  TPZAutoPointer< TPZMatrix<STATE> > matK;
  
  while(error > tol && iter < numiter) {
    
#ifdef PZ_LOG
    if(logger.isDebugEnabled())
    {
      std::stringstream sout;
      matK=an->Solver().Matrix();
      matK->Print("matK = ", sout,EMathematicaInput);
      an->Rhs().Print("Rhs = ",sout,EMathematicaInput);
      LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    
    an->Solve(); // o an->Solution() eh o deltaU aqui
    TPZFMatrix<STATE> DeltaU = an->Solution();
    SoliterK = prevsol - DeltaU;
    REAL normDeltaSol = Norm(DeltaU);
    
#ifdef PZ_LOG
    if(logger.isDebugEnabled())
    {
      std::stringstream sout;
      an->Solution().Print("DeltaX = ", sout,EMathematicaInput);
      SoliterK.Print("Xk = ", sout,EMathematicaInput);
      LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    
    an->LoadSolution(SoliterK); // Aqui o an->Solution() eh o U
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fcmeshMixed);
    an->Assemble();
    an->Rhs() += fLastStepRhs;
    
#ifdef PZ_LOG
    if(logger.isDebugEnabled())
    {
      std::stringstream sout;
      an->Rhs().Print("Res = ", sout,EMathematicaInput);
      LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    
    double NormResLambda = Norm(an->Rhs());
    double norm = normDeltaSol;
    out << "Iteracao n : " << (iter+1) << " : normas |Delta(Un)| e |Delta(rhs)| : " << normDeltaSol << " / " << NormResLambda << std::endl;
    
    //SoliterK.Print("Sol");
    //an->Rhs().Print("Rhs:");
    
    if(norm < tol) {
      out << "\nTolerancia do DELTAU atingida na iteracao : " << (iter+1) << std::endl;
      out << "\n\nNorma da solucao |Delta(Un)|  : " << norm << std::endl << std::endl;
      
    }
    else if( (NormResLambda - NormResLambdaLast) > 1.e-4 ) {
      out << "\nDivergent Method\n" << "Applying Line Search" << std::endl;
      STATE scale = 0.5;
      int itls = 0;
      while (NormResLambda > NormResLambdaLast) {
        SoliterK = prevsol - scale * DeltaU;
        //SoliterK.Print("Sol");
        an->LoadSolution(SoliterK);
        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fcmeshMixed);
        an->Assemble();
        an->Rhs() += fLastStepRhs;
        NormResLambda = Norm(an->Rhs());
        scale *= scale;
        out << " Line Search iter: " << itls << " : normas |Delta(Un)| e |Delta(rhs)| : " << Norm(scale*DeltaU) << " / " << NormResLambda << std::endl;
        itls++;
      }
      
    }
    
    NormResLambdaLast = NormResLambda;
    error = norm;
    iter++;
    prevsol = SoliterK;
    out.flush();
  }
  
  if (error > tol && numiter != 1) {
    DebugStop(); // Metodo nao convergiu!!
  }
  
}

void TPZFracAnalysis::AssembleLastStep(TPZAnalysis *an)
{
  fData->SetLastState();
  an->Assemble();
  fLastStepRhs = an->Rhs();
  
}

bool TPZFracAnalysis::SolveSistTransient(TPZAnalysis *an)
{
  
  bool propagate = false;
  int nfracel = this->HowManyFracElement();
  int it = 0;
  while (fmustStop == false && propagate == false) {
    
    AssembleLastStep(an);
    TPZFMatrix<STATE> lastSol = an->Solution();
    //fLastStepRhs.Print("Rhsatn: ");
    //lastSol.Print("SOlLast: ");
    //    TPZEquationFilter eqF(fcmeshMixed->NEquations());
    if (it == 0) { // aqui inicializo chutes iniciais para newton depois da propagacao
      
      if(nfracel == 1){ // esse caso eh para o primeiro elemento da simulacao
        //an->Solution().Print("anBefore");
        this->ComputeFirstSolForOneELement(an);
        //an->Solution().Print("anAfter");
        
        
      }
      else{ // aqui eh quando ha mais de 1 elemento
        //an->Solution().Print("anBefore");
        this->SetPressureOnLastElement(an);
        //an->Solution().Print("anAfter");
        
        //        // Fixo a pressao no ultimo elemento, e o resultado eh o chute inicial para o newton
        //        TPZCompEl *cel = fcmeshMixed->Element(fcmeshMixed->ElementVec().NElements()-2);
        //        TPZConnect &c = cel->Connect(3);
        //        const int sq = c.SequenceNumber();
        //        const int pos = fcmeshMixed->Block().Position(sq);
        //        if (fcmeshMixed->Block().Size(sq) != 1) {
        //          DebugStop();
        //        }
        //        std::set<int64_t> eqOut;
        //
        //        eqOut.insert(pos);
        //
        //        const int neq = fcmeshMixed->NEquations();
        //        TPZVec<int64_t> actEquations(neq-eqOut.size());
        //        int p = 0;
        //        for (int64_t eq = 0; eq < neq; eq++) {
        //          if (eqOut.find(eq) == eqOut.end()) {
        //            actEquations[p] = eq;
        //            p++;
        //          }
        //        }
        //
        //        eqF.SetActiveEquations(actEquations);
        //
        
      }
      
      fData->SetAccumVl(0.); // Zerando accumvl para os proximos
      
    }
    
//    const REAL tol = 1.e-8;
    //    do{
    //
    //      an->StructMatrix()->EquationFilter() = eqF;
    //      an->Solver().SetMatrix(NULL);
    //      IterativeProcess(an, std::cout);
    //      an->StructMatrix()->EquationFilter().Reset();
    //
    //      IterativeProcess(an, std::cout, 1);
    //      an->Rhs()(0,0) = 0.;
    //
    //    }while (Norm(an->Rhs()) > tol);

      
    IterativeProcess(an, std::cout);
    
    
    const REAL qtip = this->Qtip();
      
    if (qtip < 0.) {
      DebugStop();
      propagate = false;
    }
    else{
      propagate = this->VerifyIfPropagate(qtip);
    }
      
    if (propagate) {
      fData->DebugMap().insert(std::pair<REAL, REAL> (fData->Time(),fData->Lfrac()));
    }
    
    
    if (propagate) {
      fData->SetLastQtip(qtip);
      an->Solution() = lastSol;
      an->LoadSolution();
      TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fcmeshMixed);
      std::cout << "\n******************************** FRACTURE PROPAGATED ********************************" << std::endl;
    }
    else{
      fData->SetNextTime();
      if (qtip > 0.) {
        AcceptSolution(an); // updates leak off
      }
      this->PostProcessVTK(an);
    }
    
    REAL peteleco = 1.e-8;
    if( fData->Time() > (fData->TotalTime() - peteleco) )
    {
      fmustStop = true;
    }
    it++;
  }
  return propagate;
}

void TPZFracAnalysis::AcceptSolution(TPZAnalysis *an)
{
  fMatFrac->SetUpdateMem();
  an->AssembleResidual();
  fMatFrac->SetUpdateMem(false);
}

void TPZFracAnalysis::PostProcessVTK(TPZAnalysis *an)
{
  const int dim = 1;
  TPZStack<std::string> scalnames, vecnames;
  scalnames.Push("Pressure");
  scalnames.Push("Flow");
  scalnames.Push("Opening");
  an->DefineGraphMesh(dim, scalnames, vecnames, fData->PostProcessFileName());
  an->PostProcess(0);
}

REAL TPZFracAnalysis::Qtip()
{
  fgmesh->ResetReference();
  fcmeshMixed->LoadReferences();
  const int bcRightId = -2;
  const int64_t nel = fcmeshMixed->NElements();
  TPZCompEl *cel = NULL;
  for (int64_t iel = 0; iel < nel; iel++) {
    cel = fcmeshMixed->Element(iel);
    if (!cel) continue;
    TPZGeoEl *gel = cel->Reference();
    if (gel->Dimension() != 1) {
      continue;
    }
    TPZGeoElSide gelside(gel,1);
    TPZGeoElSide neigh = gelside.Neighbour();
    while (neigh != gelside) {
      if (neigh.Element()->MaterialId() == bcRightId) {
        break;
      }
      neigh = neigh.Neighbour();
    }
    if (neigh == gelside) {
      continue;
    }
    break;
  }
  
  TPZVec<REAL> qsi(3,1.);
  TPZVec<STATE> sol(1,0.);
  const int varQ = fMatFrac->VariableIndex("Flow");
  cel->Solution(qsi, varQ, sol);
  const REAL qTip = sol[0];
  std::cout << "\nqtip = " << qTip << std::endl;
  
  return qTip;
}



TPZGeoEl* TPZFracAnalysis::CreateFirstGeoElWithBC()
{
  const int nnodes = 2;
  const int nel = nnodes - 1;
  const int matIdFrac = 1, bcLeftId = -1, bcRightId = -2;
  
  // Nos
  const REAL elsize = fData->ElSize();
  fgmesh->NodeVec().Resize(nnodes);
  TPZVec<REAL> coord(3,0.);
  for (int in = 0; in < nnodes; in++) {
    coord[0] = in * elsize;
    fgmesh->NodeVec()[in].SetNodeId(in);
    fgmesh->NodeVec()[in].SetCoord(coord);
  }
  
  // Elementos
  TPZVec<int64_t> TopolLine(2,0);
  int64_t index = 0;
  TPZGeoEl *gel = NULL;
  for (int iel = 0; iel < nel; iel++) {
    TopolLine[0] = iel;
    TopolLine[1] = iel+1;
    gel = fgmesh->CreateGeoElement(EOned, TopolLine, matIdFrac, index);
  }
  
  fgmesh->BuildConnectivity();
  
  // Left
  TPZVec<int64_t> TopolPoint(1,0);
  fgmesh->CreateGeoElement(EPoint, TopolPoint, bcLeftId, index);
  
  // Right
  TopolPoint[0] = nnodes-1;
  fgmesh->CreateGeoElement(EPoint,TopolPoint,bcRightId,index);
  
  fgmesh->BuildConnectivity();
  
#ifdef PZDEBUG
  std::ofstream out("GeoMesh.vtk");
  TPZVTKGeoMesh::PrintGMeshVTK(fgmesh, out, true);
#endif
  
  return gel;
}

bool TPZFracAnalysis::VerifyIfPropagate(REAL qtip)
{
  const REAL dt = fData->TimeStep();
//  const REAL AccumVolThroughTip = fData->AccumVl() * fData->ElSize() * 2.;
//  const REAL volThroughTip = qtip * dt + AccumVolThroughTip;
  const REAL pfrac = fData->SigmaConf();
  const REAL tstar = fData->FictitiousTime(fData->AccumVl(), pfrac); // VlForNextPropag is vl from last propag here
  REAL vl = fData->VlFtau(pfrac, tstar+dt);
  const REAL totalLeakOff = 2. * fData->ElSize() * vl;
  const REAL totalLeakOffprev = 2. * fData->ElSize() * fData->AccumVl();
  const REAL ql = (totalLeakOff - totalLeakOffprev)/dt;
  
  const REAL qFreshNewEl = this->QOfAFreshNewElement();
  const REAL crit = this->PropagationFlowCriteria(qFreshNewEl,ql);
  if (qtip > crit) { // AQUINATHAN
    return true;
  }
  else{
    vl = fData->AccumVl() + qtip*dt/fData->ElSize()/2.;
    fData->SetAccumVl(vl);
    return false;
  }
}

REAL TPZFracAnalysis::RunUntilOpen()
{
  const int maxinitialit = 10000;
  const REAL qtip = fData->Q();
  int it = 0;
  for (it = 0; it < maxinitialit; it++) {
    const REAL dt = fData->TimeStep();
//    const REAL AccumVolThroughTip = fData->AccumVl() * fData->ElSize() * 2.;
//    const REAL volThroughTip = qtip * dt + AccumVolThroughTip;
    const REAL pfrac = fData->SigmaConf();
    const REAL tstar = fData->FictitiousTime(fData->AccumVl(), pfrac); // VlForNextPropag is vl from last propag here
    REAL vlnext = fData->VlFtau(pfrac, tstar+dt);
    const REAL totalLeakOff = 2. * fData->ElSize() * vlnext;
    const REAL totalLeakOffPrev = 2. * fData->ElSize() * fData->AccumVl();
    const REAL ql = (totalLeakOff - totalLeakOffPrev)/dt;
    
    const REAL qFreshNewEl = this->QOfAFreshNewElement();
    const REAL crit = this->PropagationFlowCriteria(qFreshNewEl,ql);
    if (qtip > crit) { //AQUINATHAN
      break;
    }
    
    vlnext = fData->AccumVl() + qtip*dt/fData->ElSize()/2.;
    fData->SetAccumVl(vlnext);
    fData->SetNextTime();
  }
  if (it == maxinitialit) {
    DebugStop();
  }
  
  std::cout << "#################### Opening of the fracture occured at time t = " << fData->Time() << " s ####################" << std::endl;
  std::cout << "Total vol injected = " << qtip *fData->Time() << std::endl;
  std::cout << "\nStarting First Simulation" << std::endl;
  
  return fData->AccumVl();
}

REAL TPZFracAnalysis::PropagationFlowCriteria(REAL qFreshNewEl, REAL ql){
  return MAX(qFreshNewEl/5.0, 3.*ql);
}

TPZGeoEl * TPZFracAnalysis::FindPressureBCElement()
{
  TPZGeoEl *gel = NULL;
  const int bcpressureid = -2;
  for (int64_t iel = 0; iel < fgmesh->NElements(); iel++) {
    gel = fgmesh->ElementVec()[iel];
    if (gel->MaterialId() == bcpressureid) {
      break;
    }
  }
#ifdef PZDEBUG
  if (gel == NULL) {
    DebugStop();
  }
#endif
  return gel;
}

REAL TPZFracAnalysis::QOfAFreshNewElement()
{
  const REAL dt = fData->TimeStep();
  const REAL pfrac = fData->SigmaConf();
  REAL vlnext = fData->VlFtau(pfrac, dt);
  const REAL totalLeakOff = 2. * fData->ElSize() * vlnext;
  const REAL qFresh = (totalLeakOff)/dt;
  
  return qFresh;
}

void TPZFracAnalysis::SetPressureOnLastElement(TPZAnalysis *an)
{
  fgmesh->ResetReference();
  fmeshvec[1]->LoadReferences();
  
  int64_t nfracel = this->HowManyFracElement();
  int nel = fmeshvec[1]->NElements();
  for (int iel = 0; iel < nel; iel++) {
    TPZCompEl * cel = fmeshvec[1]->Element(iel);
    if (!cel) continue;
    TPZGeoEl * gel = cel->Reference();
    
    if (gel->Dimension() != 1) {
      continue;
    }
    
    int side = 1;
    TPZGeoElSide gelside(gel,side);
    TPZGeoElSide neigh = gelside.Neighbour();
    
    // Seeking for condition 1
    while (gelside!= neigh) {
      if (neigh.Element()->Dimension()==0 && neigh.Element()->MaterialId() == -2) {
        break;
      }
      neigh = neigh.Neighbour();
    }
    if (gelside == neigh) {
      continue;
    }
    
    TPZBlock<STATE> & block = fmeshvec[1]->Block();
    TPZGeoElSide gelsideleft(gel,0);
    TPZGeoElSide neighTip = gelsideleft.Neighbour();
    
    // Seeking for condition 2
    while (gelsideleft != neighTip) {
      if (neighTip.Element()->Dimension() == 1 && neighTip.Element()->MaterialId() == 1 ) {
        break;
      }
      neighTip = neighTip.Neighbour();
    }
    if (gelsideleft == neighTip) {
      DebugStop();
    }
    TPZGeoEl * gelleft = neighTip.Element();
    TPZCompEl * celleft = gelleft->Reference();
    
    
    // Chanching value
#ifdef PZDEBUG
    if (celleft->NConnects() != 1 && cel->NConnects() != 1) {
      DebugStop();
    }
#endif
    TPZConnect &connectleft =  celleft->Connect(0);
    TPZConnect &connect =  cel->Connect(0);
    
    int seqleft = connectleft.SequenceNumber();
    int seq = connect.SequenceNumber();
#ifdef PZDEBUG
    if (block.Size(seqleft) != 1 && block.Size(seq) != 1) {
      DebugStop();
    }
#endif
    int posleft = block.Position(seqleft);
    int pos = block.Position(seq);
    fmeshvec[1]->Solution()(pos,0) = fmeshvec[1]->Solution()(posleft,0);
    
    if (nfracel < 5) {
      TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, fcmeshMixed);
      an->LoadSolution(fcmeshMixed->Solution());
      return;
    }
    
    // Seeking for condition 3
    TPZGeoElSide gelsidesecondleft(gelleft,0);
    TPZGeoElSide neighsecondeleft = gelsidesecondleft.Neighbour();
    while (gelsidesecondleft != neighsecondeleft) {
      if (neighsecondeleft.Element()->Dimension() == 1 && neighsecondeleft.Element()->MaterialId() == 1 ) {
        break;
      }
      neighsecondeleft = neighsecondeleft.Neighbour();
    }
    if (gelsidesecondleft == neighsecondeleft) {
      DebugStop();
    }
    TPZGeoEl * gelsecondleft = neighsecondeleft.Element();
    TPZCompEl * celsecondleft = gelsecondleft->Reference();
#ifdef PZDEBUG
    if (celsecondleft->NConnects() != 1) {
      DebugStop();
    }
#endif
    TPZConnect &connectsecondleft = celsecondleft->Connect(0);
    int seqsecondleft = connectsecondleft.SequenceNumber();
    
#ifdef PZDEBUG
    if (block.Size(seqsecondleft) != 1) {
      DebugStop();
    }
#endif
    int possecondleft = block.Position(seqsecondleft);
    
    fmeshvec[1]->Solution()(posleft,0) = fmeshvec[1]->Solution()(possecondleft,0);
    
    TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, fcmeshMixed);
    an->LoadSolution(fcmeshMixed->Solution());
    return;
  }
}

int TPZFracAnalysis::HowManyFracElement()
{
  int64_t nel = fgmesh->NElements();
  int64_t nfracel = 0;
  for (int64_t iel = 0; iel < nel; iel++) {
    TPZGeoEl *gel = fgmesh->Element(iel);
    if (gel->MaterialId() == 1) {
      nfracel++;
    }
  }
  return nfracel;
}

void TPZFracAnalysis::ComputeFirstSolForOneELement(TPZAnalysis * an)
{
  int64_t nfracel = this->HowManyFracElement();
  
  if (nfracel != 1) {
    PZError << "This method sould only be called when the mesh has a single frac element " << std::endl;
    DebugStop();
  }
  
  int nel = fgmesh->NElements();
  TPZGeoEl *gel = NULL;
  for (int iel = 0; iel < nel; iel++) {
    gel = fgmesh->Element(iel);
    if (gel->MaterialId() == 1) {
      break;
    }
  }
  
  TPZBlock<STATE> &blockQ = fmeshvec[0]->Block(), &blockP = fmeshvec[1]->Block();
  
  
  // Setando os valores dos fluxos na fratura
  fgmesh->ResetReference();
  fmeshvec[0]->LoadReferences();
  
  const REAL pfrac = fData->SigmaConf();
  const REAL tstar = fData->FictitiousTime(fData->AccumVl(), pfrac); // VlForNextPropag is vl from last propag here
  const REAL vlnext = fData->VlFtau(pfrac, tstar+fData->TimeStep());
  const REAL totalLeakOff = 2. * fData->ElSize() * vlnext;
  const REAL totalLeakOffPrev = 2. * fData->ElSize() * fData->AccumVl();
  const REAL ql = (totalLeakOff - totalLeakOffPrev)/fData->TimeStep();
  REAL qout = fData->Q() - ql;
  
  TPZCompEl *celQ = gel->Reference();
  if (celQ->NConnects() != 3) {
    DebugStop(); // Mesh H1 1D p = 1
  }
  TPZConnect &c1Q = celQ->Connect(0), &c2Q = celQ->Connect(1);
  int seq1Q = c1Q.SequenceNumber(), seq2Q = c2Q.SequenceNumber();
  int pos1Q = blockQ.Position(seq1Q), pos2Q = blockQ.Position(seq2Q);
  fmeshvec[0]->Solution()(pos1Q,0) = fData->Q();
  fmeshvec[0]->Solution()(pos2Q,0) = qout;
  
  
  // Setando a pressao
  fgmesh->ResetReference();
  fmeshvec[1]->LoadReferences();
  const REAL dwdp = fData->GetDwDp();
  const REAL pini = fData->SigmaConf() + pow(12. * fData->Viscosity() * qout * fData->ElSize() / (dwdp*dwdp*dwdp),1./4.);
  
  TPZCompEl *celP = gel->Reference();
  if (celP->NConnects() != 1) {
    DebugStop(); // Mesh L2 1D p = 0
  }
  TPZConnect &cP = celP->Connect(0);
  int seqP = cP.SequenceNumber();
  int posP = blockP.Position(seqP);
  
  fmeshvec[1]->Solution()(posP,0) = pini;
  
  TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, fcmeshMixed);
  an->LoadSolution(fcmeshMixed->Solution());
  
}

