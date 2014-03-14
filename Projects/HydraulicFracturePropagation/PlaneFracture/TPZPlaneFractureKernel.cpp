//
//  TPZPlaneFractureKernel.cpp
//  PZ
//
//  Created by Cesar Lucci on 18/11/13.
//
//

#include "TPZPlaneFractureKernel.h"


#include "pzbndcond.h"
#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzstepsolver.h"
#include "pzreducedspace.h"
#include "pzskylstrmatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzbuildmultiphysicsmesh.h"

#include "pzintel.h"
#include "tpzintpoints.h"
#include "TPZVTKGeoMesh.h"

#include "TSWXGraphMesh.h"
#include "TSWXGraphElement.h"

//Utilize 1. para output (Mathematica) em metros e 3.280829131 para output (Mathematica) em foot
const REAL feet = 1.;//3.280829131;

//Inicializando vetor de cores
const std::string TPZPlaneFractureKernel::color[12] = {"Red","Green","Blue","Black","Gray","Cyan","Magenta","Yellow","Brown","Orange","Pink","Purple"};


TPZPlaneFractureKernel::TPZPlaneFractureKernel() : actColor(0)
{
    DebugStop();//Use constructor below;
}
//------------------------------------------------------------------------------------------------------------

TPZPlaneFractureKernel::TPZPlaneFractureKernel(TPZVec<LayerProperties> & layerVec, REAL bulletTVDIni, REAL bulletTVDFin,
                                               REAL xLength, REAL yLength, REAL Lmax, bool just1Stripe, REAL Qinj_well, REAL visc,
                                               REAL Jradius,
                                               int pOrder,
                                               REAL MaxDispl,
                                               bool pressureIndependent) : actColor(0)
{
    this->fpoligonalChain.Resize(0);
    this->fstep = 0;
    
    this->fLmax = Lmax;
    
    globLayerStruct.SetLayerVec(layerVec);
    this->fPlaneFractureMesh = new TPZPlaneFractureMesh(bulletTVDIni, bulletTVDFin, xLength, yLength, Lmax, just1Stripe);
    
    this->fmeshVec.Resize(2);
    this->fmeshVec[0] = NULL;
    this->fmeshVec[1] = NULL;
    
    this->fmphysics = NULL;
    
    this->fHbullet = (bulletTVDFin - bulletTVDIni);
    
    REAL Qinj1wing = Qinj_well/2.;
    
    {
        globFractOutput3DData.SetQinj1wing(-Qinj1wing);
    }
    
    REAL Qinj1wing_Hbullet = Qinj1wing/this->fHbullet;
    
    this->fQinj1wing_Hbullet = Qinj1wing_Hbullet;
    
    this->fCenterTVD = (bulletTVDIni + bulletTVDFin)/2.;
    
    this->fvisc = visc;
    
    this->fJIntegralRadius = Jradius;
    
    this->fpOrder = pOrder;
    
    this->fMaxDispl = MaxDispl;
    
    this->fPath3D.Reset();
    
    if(pressureIndependent)
    {
        globLeakoffStorage.SetPressureIndependent();
    }
    else
    {
        globLeakoffStorage.SetPressureDependent();
    }
    
    this->fResTop = -layerVec[0].fTVDini;
    this->fResBottom = -layerVec[layerVec.NElements()-1].fTVDfin;
    this->fResRight = xLength;
}
//------------------------------------------------------------------------------------------------------------

TPZPlaneFractureKernel::~TPZPlaneFractureKernel()
{
    delete this->fPlaneFractureMesh;
    
    for(int m = 0; m < this->fmeshVec.NElements(); m++)
    {
        delete this->fmeshVec[m];
    }
    this->fmeshVec.Resize(0);
    
    delete this->fmphysics;
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::RunUncoupled()
{
    if(this->fPlaneFractureMesh->Just1Stripe() == false)
    {
        std::cout << "\n\n\nMétodo desacoplado só funciona para apenas 1 faixa (just1Stripe = true)\n\n\n";
        DebugStop();
    }
    
    this->InitializePoligonalChain();
    TPZCompMesh * lastPressureCMesh = NULL;
    
    REAL fractVolum = 0.;
    while(globTimeControl.ReachEndOftime() == false)
    {
        std::cout << "\n\n=============================================================================\n";
        std::cout << "STEP " << this->fstep << "\n";
        this->InitializeMeshes();
        
        if(this->fstep > 0)
        {
            this->TransferLastLeakoff(lastPressureCMesh);
        }
        
        REAL maxKI_KIc = 0.;
        std::set<int> whoPropagate;
        
        this->ComputeStressAppliedThatpropagatesWhithKI1(maxKI_KIc, whoPropagate);
        
        bool thereIsNegW = false;
        REAL negVol = 0.;
        fractVolum = this->IntegrateW(thereIsNegW, negVol);
        this->PredictActDeltaT(fractVolum);
        
        {
            this->CloseActualTimeStepUncoupled();
            this->DefinePropagatedPoligonalChain(maxKI_KIc, whoPropagate);
            this->fstep++;
        }
        
        lastPressureCMesh = this->fmeshVec[1];
    }//end of while(reachEndOfTime == false)
    
    std::ofstream outMath("GlobalProstProcess.txt");
    globFractOutput3DData.PrintMathematica(outMath);
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::RunCoupled()
{
    this->InitializePoligonalChain();
    
    REAL fractVolum = 0.;
    TPZCompMesh * lastPressureCMesh = NULL;
    
    while(globTimeControl.ReachEndOftime() == false)
    {
        std::cout << "\n\n=============================================================================\n";
        std::cout << "STEP " << this->fstep << "\n";
        this->InitializeMeshes();
        
        if(this->fstep > 0)
        {
            this->TransferElasticSolution(fractVolum);
            this->TransferLastLeakoff(lastPressureCMesh);
        }
        
        globTimeControl.RestartBissection();
        this->PredictActDeltaT(fractVolum);

        //Resolvendo o problema acoplado da nova geometria
        this->RunThisFractureGeometry(fractVolum, false);
        
        lastPressureCMesh = this->fmeshVec[1];
    }//end of while(reachEndOfTime == false)
    
    std::ofstream outMath("GlobalProstProcess.txt");
    globFractOutput3DData.PrintMathematica(outMath);
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::InitializePoligonalChain()
{
    /** Elipse */
    int nptsUp = 50;
    fpoligonalChain.Resize(0);
    
    REAL yc = -fCenterTVD;
    REAL sAx = this->fLmax/2.;
    REAL sAy = this->fHbullet/2. + this->fLmax/4.;
    
    REAL shiftX = this->fLmax;
    for(int p = 1; p <= nptsUp; p++)
    {
        REAL x  = MIN(sAx - 0.001,p * sAx/nptsUp);
        REAL fx = yc + (sAy*sqrt(sAx*sAx - x*x))/sAx;
        x += shiftX;
        if(x > this->fLmax)
        {
            int oldSize = fpoligonalChain.NElements();
            this->fpoligonalChain.Resize(oldSize+1);
            this->fpoligonalChain[oldSize] = std::make_pair(x,fx);
        }
    }
    for(int p = nptsUp-1; p > 0; p--)
    {
        REAL x  = MIN(sAx - 0.001,p * sAx/nptsUp);
        REAL fx = yc - (sAy*sqrt(sAx*sAx - x*x))/sAx;
        x += shiftX;
        if(x > this->fLmax)
        {
            int oldSize = fpoligonalChain.NElements();
            this->fpoligonalChain.Resize(oldSize+1);
            this->fpoligonalChain[oldSize] = std::make_pair(x,fx);
        }
    }
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::InitializeMeshes()
{
    globLayerStruct.ResetData();
    
    //GeoMesh
    this->fPlaneFractureMesh->InitializeFractureGeoMesh(fpoligonalChain);
    
    std::cout << "\n************** GERANDO MALHAS COMPUTACIONAIS\n";
    
    //Malha computacional elastica processada para newman 0 (serah referencia da this->fmeshVec[0])
    TPZCompMesh * fractureCMeshRef = this->fPlaneFractureMesh->GetFractureCompMesh(this->fpOrder);
    this->ProcessLinearElasticCMesh(fractureCMeshRef);
    //
    //Malha computacional do tipo CMeshReferred
    this->fmeshVec[0] = this->fPlaneFractureMesh->GetFractureCompMeshReferred(fractureCMeshRef, this->fpOrder);
    
    //Malha computacional de pressao
    this->fmeshVec[1] = this->fPlaneFractureMesh->GetPressureCompMesh(this->fQinj1wing_Hbullet, this->fpOrder);
    
    //Malha computacional de acoplamento (multifisica)
    this->fmphysics = this->fPlaneFractureMesh->GetMultiPhysicsCompMesh(this->fmeshVec,
                                                                        this->fQinj1wing_Hbullet,
                                                                        this->fvisc,
                                                                        this->fpOrder);
    this->ApplyInitialCondition();
    
    this->InitializePath3DVector();
    
    {//just4Debug
//        REAL sigYY0 = -(globLayerStruct.GetLayer(0).fSigYY);
//        REAL sigYY1 = -(globLayerStruct.GetLayer(1).fSigYY);
//        REAL sigYY2 = -(globLayerStruct.GetLayer(2).fSigYY);
//        
//        REAL stressApplied = sigYY1;
//        
//        TPZVec< std::pair<REAL,REAL> > poligonalChainIni = this->fpoligonalChain;
//        while (stressApplied < 3.75) {
//            this->fmeshVec[0]->Solution()(1,0) = stressApplied;
//            
//            ///TODOS OS CONTATOS
//            this->fmeshVec[0]->Solution()(2,0) = MAX(0.,sigYY0 - stressApplied);
//            this->fmeshVec[0]->Solution()(3,0) = MAX(0.,sigYY1 - stressApplied);
//            this->fmeshVec[0]->Solution()(4,0) = MAX(0.,sigYY2 - stressApplied);
//            
//            std::cout << "\n\nSTRESS APPLIED = " << stressApplied << "=======================================================\n";
//            REAL maxKI_KIc = 0.;
//            std::set<int> who;
//            int maxPos;
//            this->CheckPropagationCriteria(maxKI_KIc, who, maxPos, this->fmeshVec[0]->Solution());
//            
//            PostProcessElasticity();
//            PostProcessSolutions();
//            DefinePropagatedPoligonalChain(maxKI_KIc, who);
//            PostProcessFractGeometry();
//            this->fpoligonalChain = poligonalChainIni;
//            this->fstep++;
//            stressApplied += 0.1;
//        }
//        
//        std::cout << "\n\n\nCHEGUEI\n\n\n";
//        DebugStop();
    }
}
//------------------------------------------------------------------------------------------------------------

#include "TPZTimer.h"
void TPZPlaneFractureKernel::ProcessLinearElasticCMesh(TPZCompMesh * cmesh)
{
    std::cout << "\n************** CALCULANDO SOLUCOES ELASTICAS DE REFERENCIA\n";
    
    {
        int neq = cmesh->NEquations();
        std::cout << "\nNequacoes elastica 3D = " << neq << "\n";
    }
    TPZAnalysis * an = new TPZAnalysis(cmesh);

    TPZSkylineStructMatrix skyl(cmesh); //caso simetrico
    an->SetStructuralMatrix(skyl);
    TPZStepSolver<REAL> stepS;
    stepS.SetDirect(ECholesky);
    an->SetSolver(stepS);

    TPZFMatrix<STATE> solution0(cmesh->Solution().Rows(), 1);
    TPZFMatrix<STATE> solution1(cmesh->Solution().Rows(), 1);
    TPZFMatrix<STATE> solutionLayStripe(cmesh->Solution().Rows(), 1);
    
    TPZFMatrix<STATE> solutions(cmesh->Solution().Rows(), 1);
    
    ////////////////// Newman=0 em toda a fratura
    TPZTimer ts;
    an->Assemble();
    ts.start();
    an->Solve();
    ts.stop();
    std::cout << "Solve = " << ts.seconds() << " s\n\n";
    
    solution0 = cmesh->Solution();
    for(int r = 0; r < solution0.Rows(); r++)
    {
        solutions(r,0) = solution0(r,0);
    }
    
    
    int NStripes = this->fPlaneFractureMesh->NStripes();
    
    ////////////////// Pressao de fluido
    ////////////////// Newman=1.E7 * globStressScale em cada faixa da fratura
    for(int stripe = 0; stripe < NStripes; stripe++)
    {
        this->fPlaneFractureMesh->SetNewmanOnThisStripe(cmesh,stripe);
        an->Rhs().Zero();
        an->AssembleResidual();
        an->Solve();
        solution1 = cmesh->Solution();
        
        int oldSize = solutions.Cols();
        solutions.Resize(solutions.Rows(), oldSize+1);
        for(int r = 0; r < solution0.Rows(); r++)
        {
            solutions(r,oldSize) = (solution1(r,0) - solution0(r,0));
        }
    }
    
    ///////////////// Contato
    ///////////////// Newman=1.E7 * globStressScale em partes (cada subdominio formado por layer e stripe)
    bool newmanApplied = false;
    for(int lay = 0; lay < globLayerStruct.NLayers(); lay++)
    {
        for(int stripe = 0; stripe < NStripes; stripe++)
        {
            newmanApplied = this->fPlaneFractureMesh->SetNewmanOnThisLayerAndStripe(cmesh, lay, stripe);
            if(newmanApplied)
            {
                an->Rhs().Zero();
                an->AssembleResidual();
                an->Solve();
                
                solutionLayStripe = cmesh->Solution();
                
                int oldSize = solutions.Cols();
                solutions.Resize(solutions.Rows(), oldSize+1);
                for(int r = 0; r < solution0.Rows(); r++)
                {
                    solutions(r,oldSize) = (solutionLayStripe(r,0) - solution0(r,0));
                }
                
                globLayerStruct.SetLayerStripe_ContactSolutionRow(lay,stripe,oldSize);
            }
        }
    }
    
    cmesh->LoadSolution(solutions);
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::ComputeStressAppliedThatpropagatesWhithKI1(REAL & maxKI_KIc, std::set<int> & whoPropagate)
{
    ////////////////Stress Applied solution row (LowerPreStress layer share solution row with others layers that have lower prestress than stress applied)
    
    REAL minSigYY = globLayerStruct.GetLowerPreStress();
    REAL maxSigYY = globLayerStruct.GetHigherPreStress();
    
    REAL stressLeft = minSigYY;
    REAL stressRight = maxSigYY;
    
    REAL stressStep = (stressRight - stressLeft)/15.;
    
    stressLeft += stressStep;
    stressRight = -1.;//DO NOT DELETE!!! Will be used to know if stressApplied reach KI >= KIc once.
    
    bool couldStop = false;
    int maxKI_KIcPosition = -1;
    
    bool lastWasTurbo = false;
    REAL stressApplied = stressLeft;
    while(couldStop == false)
    {
        if(stressRight < 0.)
        {//Ainda nao atingiu KI >= KIc
            stressApplied = stressLeft;
        }
        else
        {//Jah atingiu KI >= KIc
            stressApplied = (stressLeft + stressRight) / 2.;
        }
        
        this->fmeshVec[0]->Solution()(0,0) = 1.;
        this->fmeshVec[0]->Solution()(1,0) = stressApplied;//Metodo desacoplado soh funciona para nstripes=1, portanto estah hardcode.
        this->ComputeContactStress();
        
        if( ( (maxKI_KIc > 0.85) && (maxKI_KIc < 1.5) ) == false )
        {
            std::cout << "Complete lap.................................";
            std::cout.flush();
            
            maxKI_KIc = 0.;
            whoPropagate.clear();
            this->CheckPropagationCriteria(maxKI_KIc, whoPropagate, maxKI_KIcPosition, this->fmeshVec[0]->Solution());
            lastWasTurbo = false;
        }
        else
        {
            std::cout << "Turbo.................................\n";
            std::cout.flush();
            globLayerStruct.SetElastSolutionMatrix(this->fmeshVec[0]->Solution());
            this->fPath3D.IntegratePath3D(maxKI_KIcPosition);
            
            maxKI_KIc = this->fPath3D.Path(maxKI_KIcPosition)->KI()/this->fPath3D.Path(maxKI_KIcPosition)->KIc();
            
            lastWasTurbo = true;
        }
        std::cout << "maxKI/respectiveKIc = " << maxKI_KIc << "\n\n";
        std::cout.flush();
        
        if(maxKI_KIc >= 1.)
        {
            if(stressRight < 0.)
            {//Chegou aqui pela primeira vez
                stressLeft -= stressStep;
            }
            stressRight = stressApplied;
        }
        else
        {
            if(stressRight <= 0.)
            {//Ainda nao atingiu KI >= KIc
                stressLeft += stressStep;
            }
            else
            {//Jah atingiu KI >= KIc
                stressLeft = stressApplied;
            }
        }
        if(maxKI_KIc >= 1. && maxKI_KIc <= 1.1)
        {
            couldStop = true;
        }
    }
    
    if(lastWasTurbo)
    {
        //Soh para calculo completo.
        std::cout << "Final lap.................................";
        std::cout.flush();
        this->CheckPropagationCriteria(maxKI_KIc, whoPropagate, maxKI_KIcPosition, this->fmeshVec[0]->Solution());
        std::cout << "maxKI/respectiveKIc = " << maxKI_KIc << "\n\n";
        std::cout.flush();
    }
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::ComputeContactStress()
{
    REAL stressApplied = this->fmeshVec[0]->Solution()(1,0);
    
    for(int lay = 0; lay < globLayerStruct.NLayers(); lay++)
    {
        for(int stripe = 0; stripe < this->fPlaneFractureMesh->NStripes(); stripe++)
        {
            int solutionRow = globLayerStruct.GetContactSolutionRow(lay, stripe);
            if(solutionRow >= 2)
            {
                REAL preStress = -(globLayerStruct.GetLayer(lay).fSigYY);
                this->fmeshVec[0]->Solution()(solutionRow,0) = MAX(0.,preStress - stressApplied);
            }
        }
    }
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::RunThisFractureGeometry(REAL & volAcum, bool justTransferingElasticSolution)
{
    TPZAnalysis * an = new TPZAnalysis(this->fmphysics);
    
    /** Convergence test */
    //CheckConv();
    /**********************/
    
    int nEq = this->fmphysics->NEquations();
    
    TPZFMatrix<REAL> matRes_total(nEq,1,0.);
    
    TPZFMatrix<REAL> Sol_0 = this->fmphysics->Solution();
    
    TPZAutoPointer< TPZMatrix<REAL> > matK;
    TPZFMatrix<REAL> matRes_partial(nEq,1,0.);
    TPZFMatrix<REAL> matMass(nEq,1,0.);
    
    /// Backup
    TPZFMatrix<REAL> backupSol_0 = Sol_0;
    
    std::cout << "\n\n\n************** CALCULANDO SOLUCAO ACOPLADA (KERNEL)\n\n";
    
    REAL maxKI_KIc = 0;
    std::set<int> whoPropagate;
    
    bool timeLimitsIsCloseEnough = false;
    bool propagate = false;
    
    while (timeLimitsIsCloseEnough == false)
    {
        maxKI_KIc = 0.;
        whoPropagate.clear();
        propagate = false;
        
        if(justTransferingElasticSolution == false)
        {
            globTimeControl.ComputeActDeltaT();
            
            //Calculo da matriz de massa para o deltaT atual
            if(this->fstep > 0)
            {
                this->MassMatrix(matMass);
            }
        }
        
        std::cout << "\n\ndtLeft = " << globTimeControl.LeftDeltaT() << "s " <<
        " : actDeltaT = " << globTimeControl.actDeltaT() << "s " <<
        " : dtRight = " << globTimeControl.RightDeltaT() << "s\n";
        
        matRes_partial.Zero();
        this->AssembleStiffMatrixLoadVec(an, matK, matRes_partial);
        
        matRes_total = matRes_partial + matMass;
        
        REAL normRes = Norm(matRes_total);
        REAL tol = 1.e-3;
        int maxit = 10;
        int nit = 0;
        
        // Metodo de Newton
        while(normRes > tol && nit < maxit)
        {
            an->Rhs() = matRes_total;
            an->Solve();
            
            TPZFMatrix<REAL> Sol_1_minus_Sol_0 = an->Solution();
            TPZFMatrix<REAL> Sol_1 = Sol_1_minus_Sol_0 + Sol_0;
            an->LoadSolution(Sol_1);
            
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(this->fmeshVec, this->fmphysics);
            
            Sol_0 = Sol_1;
            
            matRes_partial.Zero();
            this->AssembleStiffMatrixLoadVec(an, matK, matRes_partial);
            matRes_total = matRes_partial + matMass;
            
            normRes = Norm(matRes_total);
            std::cout << "||res|| = " << normRes << std::endl;
            nit++;
        }
        
        if(normRes > tol)
        {
            if(justTransferingElasticSolution)
            {
                std::cout << "\nNão convergiu na transferencia de solucao elastica...\n";
                PostProcessFractGeometry();
                DebugStop();
            }
            
            //Quando o actDeltaT leva a um instante em que Vleakoff = Vinj, nao converge, necessitando
            //trazer o limite esquerdo do deltaT da bisseccao para a direita
            globTimeControl.TimeisOnLeft();
            
            std::cout << "\nNao convergiu (provavelmente w muito pequeno resultando em grande gradiente de pressão)...\n";
            std::cout << "************* dt está a Esquerda de KI=KIc *************\n";
        }
        else
        {
            if(justTransferingElasticSolution)
            {
                return;
            }
            
            bool thereIsNegW;
            REAL negVol = 0.;
            volAcum = this->IntegrateW(thereIsNegW, negVol);
            
            if(thereIsNegW)
            {
                globTimeControl.TimeisOnLeft();
                std::cout << "\nNegativeW\n";
                std::cout << "************* dt está a Esquerda de KI=KIc *************\n";
            }
            else
            {
                int maxKI_KIcPosition = -1;
                propagate = this->CheckPropagationCriteria(maxKI_KIc, whoPropagate, maxKI_KIcPosition, this->fmeshVec[0]->Solution());
                
                if(propagate == false)
                {
                    globTimeControl.TimeisOnLeft();
                    std::cout << "KI < KIc\n";
                    std::cout << "************* dt está a Esquerda de KI=KIc *************\n";
                }
                else
                {
                    globTimeControl.TimeisOnRight();
                    std::cout << "\nPropagate: maxKI/respectiveKIc = " << maxKI_KIc << "\n";
                    std::cout << "************* t está a Direita de KI=KIc *************\n";
                }
            }
        }
        
        timeLimitsIsCloseEnough = globTimeControl.TimeLimitsIsCloseEnough();
        if(propagate && maxKI_KIc < 1.1)
        {
            timeLimitsIsCloseEnough = true;
        }
        if(timeLimitsIsCloseEnough == false)
        {
            ////Restoring backups for deltaT next try
            Sol_0 = backupSol_0;
            this->fmphysics->LoadSolution(backupSol_0);
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(this->fmeshVec, this->fmphysics);
        }
    }
    
    if(propagate == false)
    {
        std::cout << "\n\n\nConvergiu deltaT sem propagar???\n\n\n";
        DebugStop();
    }
    
    this->CloseActualTimeStepCoupled();
    
    this->DefinePropagatedPoligonalChain(maxKI_KIc, whoPropagate);
    
    this->fstep++;
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::PredictActDeltaT(REAL fractVolum)
{
    std::map<int,REAL> leakoffmap = globLeakoffStorage.GetLeakoffMap();
    
    REAL actTime = globTimeControl.actTime();
    REAL Qinj1wing = -(this->fQinj1wing_Hbullet * fHbullet);
    
    /////Metodo das secantes
    REAL dt0 = MAX(fractVolum/Qinj1wing , 0.7*globTimeControl.Ttot());
    REAL dt1 = globTimeControl.Ttot();
    
    REAL volInj0 = Qinj1wing * (actTime + dt0);
    REAL volInj1 = Qinj1wing * (actTime + dt1);
    
    globTimeControl.SetDeltaT(dt0);
    globLeakoffStorage.UpdateLeakoff(this->fmphysics,dt0);
    REAL volLeak0 = this->ComputeVlAcumLeakoff(this->fmeshVec[1]);
    globLeakoffStorage.SetLeakoffMap(leakoffmap);
    
    globTimeControl.SetDeltaT(dt1);
    globLeakoffStorage.UpdateLeakoff(this->fmphysics,dt1);
    REAL volLeak1 = this->ComputeVlAcumLeakoff(this->fmeshVec[1]);
    globLeakoffStorage.SetLeakoffMap(leakoffmap);
    
    REAL Res0 = (volInj0 - volLeak0) - fractVolum;
    REAL Res1 = (volInj1 - volLeak1) - fractVolum;
    
#ifdef DEBUG
    if(fabs(Res0 - Res1) < 1.E-10)
    {
        std::cout << "\n\nVai dar divisão por zero!\n";
        std::cout << "Ver " << __PRETTY_FUNCTION__ << "\n\n";
        DebugStop();
    }
#endif
    REAL dtNext = (dt1*Res0 - dt0*Res1)/(Res0 - Res1);
    
    while(fabs(dt1 - dtNext) > 1.E-2)
    {
        dt0 = dt1;
        dt1 = dtNext;
        
        volInj0 = Qinj1wing * (actTime + dt0);
        volInj1 = Qinj1wing * (actTime + dt1);
        
        globTimeControl.SetDeltaT(dt0);
        globLeakoffStorage.UpdateLeakoff(this->fmphysics,dt0);
        volLeak0 = this->ComputeVlAcumLeakoff(this->fmeshVec[1]);
        globLeakoffStorage.SetLeakoffMap(leakoffmap);
        
        globTimeControl.SetDeltaT(dt1);
        globLeakoffStorage.UpdateLeakoff(this->fmphysics,dt1);
        volLeak1 = this->ComputeVlAcumLeakoff(this->fmeshVec[1]);
        globLeakoffStorage.SetLeakoffMap(leakoffmap);
        
        Res0 = (volInj0 - volLeak0) - fractVolum;
        Res1 = (volInj1 - volLeak1) - fractVolum;
        
#ifdef DEBUG
        if(fabs(Res0 - Res1) < 1.E-10)
        {
            std::cout << "\n\nVai dar divisão por zero!\n";
            std::cout << "Ver " << __PRETTY_FUNCTION__ << "\n\n";
            DebugStop();
        }
#endif
        dtNext = (dt1*Res0 - dt0*Res1)/(Res0 - Res1);
    }
    
    if(globTimeControl.actDeltaT() < 1.)
    {
        globTimeControl.SetDeltaT(1.);
    }
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::CloseActualTimeStepUncoupled()
{
    globTimeControl.UpdateActTime();//atualizando esta rodada concluida para o tempo atual
    this->UpdateLeakoff();
    this->PostProcessAll_Uncoupled();
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::CloseActualTimeStepCoupled()
{
    globTimeControl.UpdateActTime();//atualizando esta rodada concluida para o tempo atual
    this->UpdateLeakoff();
    this->PostProcessAll_Coupled();
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::InitializePath3DVector()
{
    this->fPath3D.Reset();
    
    int nCrackTipElems = this->fPlaneFractureMesh->NCrackTipElements();
    
    for(int pos = 0; pos < nCrackTipElems; pos++)
    {
        TPZGeoEl * gel1D = this->fPlaneFractureMesh->GetCrackTipGeoElement(pos);
        
#ifdef DEBUG
        if(!gel1D || gel1D->Dimension() != 1 || gel1D->MaterialId() != globMaterialIdGen.CrackTipMatId())
        {
            DebugStop();
        }
#endif
        
        TPZVec<REAL> qsi(1,0.), center(3,0.), normal(3,0.);
        gel1D->X(qsi,center);
        
        REAL n0x = gel1D->NodePtr(0)->Coord(0);
        REAL n0z = gel1D->NodePtr(0)->Coord(2);
        REAL n1x = gel1D->NodePtr(1)->Coord(0);
        REAL n1z = gel1D->NodePtr(1)->Coord(2);
        
        /**
         * Da forma como o codigo foi construido, os elementos 1D do crackTip seguem em sentido anti-horario em relacao ao eixo Y.
         * O plano na integral-J eh definido por uma normal, a qual servirah como eixo da regra da mao direita.
         * Desta forma, a normal da integral-J corresponde aa direcao oposta do crackTip, para que o arco saia de y=0 (fora da fratura),
         * passe por y>0 (no interior do meio elastico 3D) e acabe em y=0 (dentro da fratura).
         */
        normal[0] = (n0x - n1x);//x
        normal[2] = (n0z - n1z);//z
        
        REAL zCoord = center[2];
        
        REAL KIc = 0.;
        this->fPlaneFractureMesh->GetKIc_FromLayerOfThisZcoord(zCoord, KIc);
        
        int myLayer = -1;
        int myStripe = -1;
        
        TPZGeoElSide gel1Dside(gel1D,2);
        TPZGeoElSide neighSide(gel1Dside.Neighbour());

        bool found = false;
        while(neighSide != gel1Dside && found == false)
        {
            int matId = neighSide.Element()->MaterialId();
            if(globMaterialIdGen.IsInsideFractMat(matId))
            {
                myLayer = globMaterialIdGen.WhatLayerFromInsideFracture(matId);
                myStripe = globMaterialIdGen.WhatStripe(matId);
                found = true;
                break;
            }
            neighSide = neighSide.Neighbour();
        }
        if(myLayer < 0 || myStripe < 0)
        {
            std::cout << "\n\n\nPath3D layer and/or stripe not found!!!\n";
            std::cout << "See " << __PRETTY_FUNCTION__ << ".\n\n\n";
            DebugStop();
        }
        this->fPath3D.PushBackPath3D( new Path3D(this->fmeshVec[0], center, KIc, myLayer, myStripe, normal, this->fJIntegralRadius) );
    }
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::AssembleStiffMatrixLoadVec(TPZAnalysis * an,
                                                        TPZAutoPointer< TPZMatrix<REAL> > & matK,
                                                        TPZFMatrix<REAL> & matRes)
{
	this->fPlaneFractureMesh->SetActualState();
    
    TPZFStructMatrix structMat(this->fmphysics);
	an->SetStructuralMatrix(structMat);
    
    this->ApplyEquationFilter(an);
    
	TPZStepSolver<REAL> stepS;
	stepS.SetDirect(ELU);
	an->SetSolver(stepS);
    
    an->Assemble();
    
    matK = an->Solver().Matrix();
	matRes = an->Rhs();
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::ApplyInitialCondition()
{
    this->fmeshVec[0]->Solution().Zero();
    this->fmeshVec[0]->Solution()(0,0) = 1.;
    for(int stripe = 0; stripe < this->fPlaneFractureMesh->NStripes(); stripe++)
    {
        this->fmeshVec[0]->Solution()(globLayerStruct.GetStressAppliedSolutionRow(stripe),0) = globLayerStruct.GetHigherPreStress();
    }

    this->PutConstantPressureOnFluidSolution();
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::PutConstantPressureOnFluidSolution()
{
    for(int r = 0; r < this->fmeshVec[1]->Solution().Rows(); r++)
    {
        this->fmeshVec[1]->Solution()(r,0) = globLayerStruct.GetHigherPreStress();
    }
    TPZBuildMultiphysicsMesh::TransferFromMeshes(this->fmeshVec, this->fmphysics);
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::ApplyEquationFilter(TPZAnalysis * an)
{
    std::set<long> eqOut;
    
    int neq = this->fmphysics->NEquations();
    long blockAlphaEslast = this->fmphysics->ConnectVec()[0].SequenceNumber();
    long posBlock = this->fmphysics->Block().Position(blockAlphaEslast);
    eqOut.insert(posBlock);

    int firstStripeRow = 1;
    int lastStripeRow = this->fPlaneFractureMesh->NStripes();
    for(int r = 0; r < this->fmeshVec[0]->Solution().Rows(); r++)
    {
        bool isPressAppliedRow = (r >= firstStripeRow && r <= lastStripeRow);
        if(isPressAppliedRow == false)
        {
            //Desligando os graus de liberdade que correspondem ao (newman zero) e (contatos).
            eqOut.insert(posBlock+r);
        }
    }
    
    TPZVec<long> actEquations(neq-eqOut.size());
    int p = 0;
    for(long eq = 0; eq < neq; eq++)
    {
        if(eqOut.find(eq) == eqOut.end())
        {
            actEquations[p] = eq;
            p++;
        }
    }
    TPZEquationFilter eqF(neq);
    eqF.SetActiveEquations(actEquations);
    
    an->StructMatrix()->EquationFilter() = eqF;
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::MassMatrix(TPZFMatrix<REAL> & massMat)
{
    massMat.Zero();
    
    this->fPlaneFractureMesh->SetPastState();
    
	TPZSpStructMatrix structMat(this->fmphysics);
	TPZAutoPointer<TPZGuiInterface> guiInterface;
    structMat.CreateAssemble(massMat,guiInterface);
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::CheckConv()
{
    std::cout << "\n\n\nInside CheckConv\n\n\n";
    
    long neq = this->fmphysics->NEquations();
    int nsteps = 10;
    
    TPZFMatrix<REAL> xIni(neq,1);
    for(long i = 0; i < xIni.Rows(); i++)
    {
        REAL val = (double)(rand())*(1.e-10);
        xIni(i,0) = val;
    }
    TPZAnalysis *an = new TPZAnalysis(this->fmphysics);
    an->LoadSolution(xIni);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(this->fmeshVec, this->fmphysics);
    
    TPZFMatrix<REAL> actX = xIni;
    
    TPZAutoPointer< TPZMatrix<REAL> > fL_xIni;
    TPZFMatrix<REAL> f_xIni(neq,1);
    
    AssembleStiffMatrixLoadVec(an, fL_xIni, f_xIni);
    
    TPZFMatrix<REAL> fAprox_x(neq,1);
    TPZFMatrix<REAL> fExato_x(neq,1);
    
    TPZFMatrix<REAL> errorVec(neq,1,0.);
	TPZFMatrix<REAL> errorNorm(nsteps,1,0.);
    
    
    TPZAutoPointer< TPZMatrix<REAL> > fLtemp;
    TPZFMatrix<REAL> dFx(neq,1);
    
    TPZVec<REAL> deltaX(neq,0.001), alphas(nsteps);
    double alpha;
    
    std::stringstream exatoSS, aproxSS;
    exatoSS << "exato={";
    aproxSS << "aprox={";
    for(int i = 0; i < nsteps; i++)
    {
        alpha = (i+1)/10.;
        alphas[i] = alpha;
        
        ///Fx aproximado
        dFx.Zero();
        for(long r = 0; r < neq; r++)
        {
            for(long c = 0; c < neq; c++)
            {
                dFx(r,0) +=  (-1.) * fL_xIni->GetVal(r,c) * (alpha * deltaX[c]); // (-1) porque fLini = -D[res,sol]
            }
        }
        fAprox_x = f_xIni + dFx;
        
//        int wantToSeeRow = 0;
        
//        {
//            REAL aproxSol = fAprox_x(wantToSeeRow,0);
//            aproxSS << aproxSol;
//            if(i < nsteps-1)
//            {
//                aproxSS << ",";
//            }
//        }
        
        ///Fx exato
        for(long r = 0; r < neq; r++)
        {
            actX(r,0) = xIni(r,0) + (alpha * deltaX[r]);
        }
        an->LoadSolution(actX);
        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(this->fmeshVec, this->fmphysics);
        
        fExato_x.Zero();
        if(fLtemp) fLtemp->Zero();
        this->AssembleStiffMatrixLoadVec(an, fLtemp, fExato_x);
        
//        {
//            REAL exatoSol = fExato_x(wantToSeeRow,0);
//            exatoSS << exatoSol;
//            if(i < nsteps-1)
//            {
//                exatoSS << ",";
//            }
//        }
        
        ///Erro
        errorVec.Zero();
        for(long r = 0; r < neq; r++)
        {
            errorVec(r,0) = fExato_x(r,0) - fAprox_x(r,0);
        }
        
        ///Norma do erro
        double XDiffNorm = Norm(errorVec);
        errorNorm(i,0) = XDiffNorm;
    }
    aproxSS << "};";
    exatoSS << "};";
    std::cout << aproxSS.str() << std::endl;
    std::cout << exatoSS.str() << std::endl;
    std::cout << "Show[ListPlot[aprox, Joined -> True, PlotStyle -> Red],ListPlot[exato, Joined -> True]]\n";
    
    std::cout << "Convergence Order:\n";
    for(int j = 1; j < nsteps; j++)
    {
        std::cout << ( log(errorNorm(j,0)) - log(errorNorm(j-1,0)) )/( log(alphas[j]) - log(alphas[j-1]) ) << "\n";
    }
    std::cout << "\n\n\nEnding CheckConv\n\n\n";
    DebugStop();
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::UpdateLeakoff()
{
    globLeakoffStorage.UpdateLeakoff(this->fmphysics,globTimeControl.actDeltaT());

    
    //O UpdateLeakoff atualiza o leakoff com base na pressao no centro do elemento (=volLeakoffComputed).
    //Jah o kernel considera o leakoff com as pressoes calculadas em cada ponto da regra de integracao (=volLeakoffExpected).
    //Por esta razao, alguma diferenca pode ocorrer, necessitando portanto corrigir o mapa de leakoff antes de prosseguir.
    {
        REAL volInjected = this->ComputeVolInjected();
        
        bool thereIsNegW;
        REAL negVol = 0.;
        REAL volW = this->IntegrateW(thereIsNegW, negVol);
        
        REAL volLeakoffExpected = volInjected - volW;
        REAL volLeakoffComputed = ComputeVlAcumLeakoff(this->fmeshVec[1]);
        
        if(fabs(volLeakoffComputed) > 1.E-10)
        {
            REAL alpha = volLeakoffExpected/volLeakoffComputed;
            std::map<int,REAL>::iterator itLk;
            for(itLk = globLeakoffStorage.GetLeakoffMap().begin(); itLk != globLeakoffStorage.GetLeakoffMap().end(); itLk++)
            {
                itLk->second *= alpha;
            }
        }
        
        //Verificacao
        if(fabs(volLeakoffComputed - volLeakoffExpected) > 0.1)
        {
            std::cout << "\n\n\nLeakoff não está sendo calculado corretamente!\n";
            std::cout << "Computed = " << volLeakoffComputed << "\n";
            std::cout << "volInjected = " << volInjected << "\n";
            std::cout << "volW = " << volW << "\n";
            std::cout << "Expected = volInjected - volW = " << volLeakoffExpected << "\n\n\n";
            DebugStop();
        }
    }
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::PostProcessAll_Uncoupled()
{
    std::cout << "\n\n*********** POS-PROCESSAMENTO ***********\n";
    
    this->PostProcessSolutions();
    this->PostProcessAcumVolW();
    this->PostProcessVolLeakoff();
    this->PostProcessElasticity();
    this->PostProcessFractGeometry();
    std::cout << "Tempo atual desta geometria = " << globTimeControl.actTime()/60. << " minuto(s)\n\n";
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::PostProcessAll_Coupled()
{
    std::cout << "\n\n*********** POS-PROCESSAMENTO ***********\n";
    
    this->PostProcessSolutions();
    this->PostProcessAcumVolW();
    this->PostProcessVolLeakoff();
    this->PostProcessElasticity();
    this->PostProcessPressure();
    this->PostProcessFractGeometry();
    std::cout << "Tempo atual desta geometria = " << globTimeControl.actTime()/60. << " minuto(s)\n\n";
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::PostProcessSolutions()
{
    std::stringstream nm;
    nm << "Solutions_Step" << this->fstep << ".txt";
    std::ofstream solutFile(nm.str().c_str());
    
    solutFile.precision(8);
    solutFile << "\nMeshvec[0]:\n";
    for(int r = 0; r < this->fmeshVec[0]->Solution().Rows(); r++)
    {
        solutFile << this->fmeshVec[0]->Solution()(r,0) << "\n";
    }
    solutFile << "\nMeshvec[1]:\n";
    for(int r = 0; r < this->fmeshVec[1]->Solution().Rows(); r++)
    {
        solutFile << this->fmeshVec[1]->Solution()(r,0) << "\n";
    }
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::PostProcessAcumVolW()
{
    bool thereIsNegW;
    REAL negVol = 0.;
    REAL wAcum = this->IntegrateW(thereIsNegW, negVol);
    globFractOutput3DData.InsertTAcumVolW(globTimeControl.actTime(), wAcum);
    
    std::cout.precision(10);
    std::cout << "wAcum = " << wAcum << ";\n";
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::PostProcessVolLeakoff()
{
    REAL leakAcum = this->ComputeVlAcumLeakoff(this->fmeshVec[1]);
    
    globFractOutput3DData.InsertTAcumVolLeakoff(globTimeControl.actTime(), leakAcum);
    
    {
        std::map<int,REAL>::iterator it;
        TPZVec<REAL> penetrationValues(this->fmeshVec[1]->NElements(), 0.);
        for(int el = 0; el < this->fmeshVec[1]->NElements(); el++)
        {
            TPZCompEl * cel = this->fmeshVec[1]->ElementVec()[el];
            int matId = cel->Reference()->MaterialId();
            if(globMaterialIdGen.IsInsideFractMat(matId))
            {
                REAL penetration = 0.;
                it = globLeakoffStorage.GetLeakoffMap().find(cel->Reference()->Id());
                if(it != globLeakoffStorage.GetLeakoffMap().end())
                {
                    penetration = it->second;
                }
                penetrationValues[el] = penetration;
            }
        }
        
        std::stringstream nm;
        nm << "LeakoffPenetration_Step" << this->fstep << ".vtk";
        std::ofstream file(nm.str().c_str());
        TPZVTKGeoMesh::PrintCMeshVTK(this->fmeshVec[1], file, penetrationValues);
    }
    
    std::cout.precision(10);
    std::cout << "leakAcum = " << leakAcum << ";\n";
    std::cout << "VlInjected = " << ComputeVolInjected() << ";\n";
    std::cout << "wAcum+leakAcum\n";
    std::cout << "Eficiencia = wAcum/VlInjected\n";
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::PostProcessElasticity()
{
    TSWXGraphMesh grMesh;
    TSWXGraphElement grEl(0);
    TPZVec<std::string> nodalSol(1), cellSol(0);
    nodalSol[0] = "Displacement";
    //nodalSol[1] = "StressY";
    
    grEl.GenerateVTKData(this->fmeshVec[0], 3, 0., nodalSol, cellSol, grMesh);
    
    std::stringstream nm;
    nm << "Elasticity_Step" << this->fstep << ".vtk";
    std::ofstream file(nm.str().c_str());
    grMesh.ToParaview(file);
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::PostProcessPressure()
{
    TSWXGraphMesh grMesh;
    TSWXGraphElement grEl(0);
    TPZVec<std::string> nodalSol(1), cellSol(0);
    nodalSol[0] = "Pressure";
    grEl.GenerateVTKData(this->fmeshVec[1], 2, 0., nodalSol, cellSol, grMesh);
    
    std::stringstream nm;
    nm << "Pressure_Step" << this->fstep << ".vtk";
    std::ofstream file(nm.str().c_str());
    grMesh.ToParaview(file);
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::PostProcessFractGeometry()
{
    REAL Lfrac = 0;
    REAL Hsup = 0.;
    REAL Hinf = 0.;
    
    for(int p = 0; p < fpoligonalChain.NElements(); p++)
    {
        REAL x = fpoligonalChain[p].first;
        REAL z = fpoligonalChain[p].second + fCenterTVD;
        
        Lfrac = MAX(Lfrac,x);
        Hsup = MAX(Hsup,z);
        Hinf = MIN(Hinf,z);
    }
    
    globFractOutput3DData.InsertTL(globTimeControl.actTime(), Lfrac);
    globFractOutput3DData.InsertTHsup(globTimeControl.actTime(), Hsup);
    globFractOutput3DData.InsertTHinf(globTimeControl.actTime(), Hinf);
    
    {   //Output
        std::ofstream outF("000FractContours.txt");
        {   //Preamble for PostProcessFractGeometry method
            outF << "(* colors = {Red,Green,Blue,Black,Gray,Cyan,Magenta,Yellow,Brown,Orange,Pink,Purple} *)\n";
            outF << "AllPolChains={};\n";
            outF << "Lgr={};\n";
            outF << "Hsupgr={};\n";
            outF << "Hinfgr={};\n\n";
        }
        
        std::stringstream nmMath, nmAux;
        nmAux << "pcm={";
        
        for(int p = 0; p < fpoligonalChain.NElements(); p++)
        {
            globFractOutput3DData.fFractContour << "fractureDots" << p << " = {" << feet * fpoligonalChain[p].first << ","
                                                          << feet * fpoligonalChain[p].second + fCenterTVD << "};\n";
            nmAux << "fractureDots" << p;
            if(p < fpoligonalChain.NElements()-1)
            {
                nmAux << ",";
            }
        }
        nmAux << "};\n";
        nmAux << "gr" << this->fstep << "=ListPlot[pcm,Joined->True,AxesOrigin->{0,0},AspectRatio->1,PlotStyle->" << color[actColor%12] << ",AxesLabel->{\"L (m)\", \"H (m)\"}];\n";
        nmAux << "L" << this->fstep << "=Max[Transpose[pcm][[1]]];\n";
        nmAux << "Hsup" << this->fstep << "=Max[Transpose[pcm][[2]]];\n";
        nmAux << "Hinf" << this->fstep << "=-Min[Transpose[pcm][[2]]];\n";
        nmAux << "AppendTo[AllPolChains,gr" << this->fstep << "];\n";
        nmAux << "AppendTo[Lgr,{" << globTimeControl.actTime()/60. << ",L" << this->fstep << "}];\n";
        nmAux << "AppendTo[Hsupgr,{" << globTimeControl.actTime()/60. << ",Hsup" << this->fstep << "}];\n";
        nmAux << "AppendTo[Hinfgr,{" << globTimeControl.actTime()/60. << ",Hinf" << this->fstep << "}];\n";
        nmAux << "Print[\"L" << this->fstep << " = \", L" << this->fstep << "]\n";
        nmAux << "Print[\"Hsup" << this->fstep << " = \", Hsup" << this->fstep << "]\n";
        nmAux << "Print[\"Hinf" << this->fstep << " = \", Hinf" << this->fstep << "]\n";
        nmAux << "Print[\"\"]\n\n";
        globFractOutput3DData.fFractContour << nmAux.str();
        actColor++;
        
        outF << globFractOutput3DData.fFractContour.str();
        
        outF << "\nShow[AllPolChains]\n";
        outF << "l={ListPlot[Lgr,AxesOrigin->{0,0},Filling->Axis,AxesLabel->{\"Time (min)\", \"L (green), Hsup (blue), Hinf (red) (m)\"}],";
        outF << "ListPlot[Lgr,Joined->True,AxesOrigin->{0,0},PlotStyle->Green]};\n";
        outF << "hs={ListPlot[Hsupgr,Filling->Axis,AxesOrigin->{0,0}],ListPlot[Hsupgr,Joined->True]};\n";
        outF << "hi={ListPlot[Hinfgr,PlotStyle->Red,Filling->Axis,AxesOrigin->{0,0}],ListPlot[Hinfgr,Joined->True,PlotStyle->Red]};\n";
        outF << "Show[l,hs,hi,PlotRange->All]\n";
        outF.close();
    }
}
//------------------------------------------------------------------------------------------------------------

REAL TPZPlaneFractureKernel::IntegrateW(bool & thereIsNegW, REAL & negVol)
{
    thereIsNegW = false;
    negVol = 0.;
    
    REAL integralW = 0.;
    for(int c = 0; c < fmeshVec[0]->NElements(); c++)
    {
        TPZCompEl * cel = fmeshVec[0]->ElementVec()[c];
        if(!cel || globMaterialIdGen.IsInsideFractMat(cel->Reference()->MaterialId()) == false)
        {
            continue;
        }
        TPZInterpolationSpace * intel = dynamic_cast<TPZInterpolationSpace*>(cel);
        if(!intel)
        {
            DebugStop();
        }
        TPZVec<REAL> value;
        intel->Integrate(0, value);
        
#ifdef DEBUG
        if(value.NElements() != 3)
        {
            DebugStop();
        }
#endif
        if(value[1] < 0.)
        {
            thereIsNegW = true;
            negVol += fabs(value[1]);
        }
        integralW += 2.*value[1];//Integrando w = (2*uy)
    }
    
    return integralW;
}
//------------------------------------------------------------------------------------------------------------

REAL TPZPlaneFractureKernel::Fracture1wing_Area()
{
    REAL Area = 0.;
    int nelemCompMesh = fmeshVec[1]->NElements();
    for(int i = 0;  i < nelemCompMesh; i++)
    {
        TPZCompEl * cel = fmeshVec[1]->ElementVec()[i];
        if(globMaterialIdGen.IsInsideFractMat(cel->Reference()->MaterialId()) == false)
        {
            continue;
        }
        
#ifdef DEBUG
        if(!cel || cel->Reference()->Dimension() != 2)
        {
            DebugStop();
        }
#endif
        
        TPZGeoEl * gel = cel->Reference();
        
#ifdef DEBUG
        if(!gel)
        {
            DebugStop();
        }
#endif
        
        Area += gel->SideArea(gel->NSides()-1);
    }
    
    return Area;
}
//------------------------------------------------------------------------------------------------------------

REAL TPZPlaneFractureKernel::ComputeVlAcumLeakoff(TPZCompMesh * fluidCMesh)
{
    REAL leakAcum = 0.;
    
    int nelemCompMesh = fluidCMesh->NElements();
    
    TPZVec<REAL> penetrationValues(nelemCompMesh, 0.);
    
    for(int i = 0;  i < nelemCompMesh; i++)
    {
        TPZCompEl * cel = fluidCMesh->ElementVec()[i];
        if(globMaterialIdGen.IsInsideFractMat(cel->Reference()->MaterialId()) == false)
        {
            continue;
        }
        
#ifdef DEBUG
        if(!cel || cel->Reference()->Dimension() != 2)
        {
            DebugStop();
        }
#endif
        
        TPZGeoEl * gel = cel->Reference();
        
#ifdef DEBUG
        if(!gel)
        {
            DebugStop();
        }
#endif
        
        REAL Area = gel->SideArea(gel->NSides()-1);
        
        //comprimento de penetracao do fluido
        REAL penetration = 0.;
        std::map<int,REAL>::iterator it = globLeakoffStorage.GetLeakoffMap().find(gel->Id());
        if(it != globLeakoffStorage.GetLeakoffMap().end())
        {
            penetration = it->second;
        }
        penetrationValues[i] = penetration;
        
        //volume acumulado
        leakAcum += 2. * (Area * penetration);//Aqui jah eh considerada a simetria
    }
    
    return leakAcum;
}
//------------------------------------------------------------------------------------------------------------

REAL TPZPlaneFractureKernel::ComputeVolInjected()
{
    REAL actTime = globTimeControl.actTime();
    REAL Qinj1wing = -(this->fQinj1wing_Hbullet * fHbullet);
    
    REAL volInjected = actTime * Qinj1wing;
    
    return volInjected;
}
//------------------------------------------------------------------------------------------------------------

bool TPZPlaneFractureKernel::CheckPropagationCriteria(REAL & maxKI_KIc, std::set<int> & whoPropagate,
                                                      int & maxKI_KIcPosition, TPZFMatrix<REAL> & ElastReducedSolution)
{
    globLayerStruct.SetElastSolutionMatrix(ElastReducedSolution);
    
    bool propagate = false;
    
    maxKI_KIc = 0.;
    whoPropagate.clear();
    
    //Calculo das integrais-J
    this->fPath3D.IntegratePath3D();
    
    //Calculo de KI
    for(int p = 0; p < fPath3D.NPaths(); p++)
    {
        REAL cracktipKI = fPath3D.Path(p)->KI();
        REAL cracktipKIc = fPath3D.Path(p)->KIc();
        
        TPZVec<REAL> originVec = this->fPath3D.Path(p)->Origin();
        TPZVec<REAL> Jdirection = this->fPath3D.Path(p)->JDirection();
        TPZVec<REAL> ptTemp(3,0.), qsiTemp(2,0.);
        ptTemp[0] = originVec[0] + 0.2 * Jdirection[0];
        ptTemp[2] = originVec[2] + 0.2 * Jdirection[2];
        
        TPZGeoEl * gel = fPlaneFractureMesh->Find2DElementNearCrackTip(p, ptTemp);
        if( !gel || globMaterialIdGen.IsInsideFractMat(gel->MaterialId()) )
        {
            //Nao serah propagado para pontos fora do dominio
            //ou no interior da fratura (i.e., fechamento).
            fPath3D.Path(p)->SetKI(0.);
            cracktipKI = 0.;
        }
        
        if(cracktipKI >= cracktipKIc)
        {
            propagate = true;
            whoPropagate.insert(p);
        }
        if(cracktipKI/cracktipKIc > maxKI_KIc)
        {
            maxKI_KIcPosition = p;
            maxKI_KIc = cracktipKI/cracktipKIc;
        }
    }
    
    return propagate;
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::DefinePropagatedPoligonalChain(REAL & maxKI_KIc, std::set<int> & whoPropagate)
{
    TPZVec< std::pair<REAL,REAL> > newPoligonalChain(0);
    
    std::set<int>::iterator it;
    
    REAL newX = 0.;
    REAL newZ = 0.;
    for(int p = 0; p < this->fPath3D.NPaths(); p++)
    {
        /** Todo mundo propagando */
//        it = whoPropagate.find(p);
//        if(it == whoPropagate.end())
//        {//Nao propagou, portanto serah mantido o ponto medio
//            TPZVec<REAL> originVec = this->fPath3D.Path(p)->Origin();
//            
//            newX = originVec[0];
//            newZ = originVec[2];
//        }
//        else
        {//Propagou, portanto serah definido novo ponto a partir de seu centro
            REAL KI = this->fPath3D.Path(p)->KI();
            REAL KIc = this->fPath3D.Path(p)->KIc();
            
            TPZVec<REAL> originVec = this->fPath3D.Path(p)->Origin();
            TPZVec<REAL> Jdirection = this->fPath3D.Path(p)->JDirection();
            
            //Variacao linear no decorrer do tempo de fMaxDispl para fMinDispl
            REAL alpha = 1.;//alphaMin = [1.0~2.0]
            
            REAL displacement = this->fMaxDispl * pow((KI/KIc)/(maxKI_KIc),alpha);
            
            newX = originVec[0] + displacement * Jdirection[0];
            newZ = originVec[2] + displacement * Jdirection[2];
        }
        if(newX >= fLmax/2.)
        {
            int oldSize = newPoligonalChain.NElements();
            newPoligonalChain.Resize(oldSize+1);
            newPoligonalChain[oldSize] = std::make_pair(newX,newZ);
        }
    }
    
    bool thereIsZigZag = true;
    while(thereIsZigZag)
    {
        thereIsZigZag = this->RemoveZigZag(newPoligonalChain);
    }
    RemoveLayerInvasion(newPoligonalChain);
    
    bool applyBezier = true;
    if(applyBezier)
    {
        BezierCurve bz(newPoligonalChain);
        
        int npts = 50;
        this->fpoligonalChain.Resize(npts);
        for(int p = 0; p < npts; p++)
        {
            std::pair<REAL,REAL> pt;
            REAL t = p/(double(npts-1));
            bz.F(t, pt);
            this->fpoligonalChain[p] = pt;
        }
    }
    else
    {
        this->fpoligonalChain = newPoligonalChain;
    }
    
    for(int p = 0; p < this->fpoligonalChain.NElements(); p++)
    {
        REAL x = this->fpoligonalChain[p].first;
        REAL z = this->fpoligonalChain[p].second;
        
        REAL tol = 1.E-3;
        if(z > this->fResTop - tol)
        {
            std::cout << "\n\nPoligonal chain reach the top limit of available domain!\n";
            std::cout << "Simulation should stop!\n\n";
            
            std::ofstream outMath("GlobalProstProcess.txt");
            globFractOutput3DData.PrintMathematica(outMath);
            DebugStop();
        }
        if(z < this->fResBottom + tol)
        {
            std::cout << "\n\nPoligonal chain reach the bottom limit of available domain!\n";
            std::cout << "Simulation should stop!\n\n";
            
            std::ofstream outMath("GlobalProstProcess.txt");
            globFractOutput3DData.PrintMathematica(outMath);
            DebugStop();
        }
        if(x > this->fResRight - tol)
        {
            std::cout << "\n\nPoligonal chain reach the right limit of available domain!\n";
            std::cout << "Simulation should stop!\n\n";
            
            std::ofstream outMath("GlobalProstProcess.txt");
            globFractOutput3DData.PrintMathematica(outMath);
            DebugStop();
        }
    }
    
    //    {//C++ tool
    //        std::stringstream nm;
    //        nm << "PoligonalChain_Cpp_Step" << this->fstep  << ".txt";
    //        std::ofstream outPolig(nm.str().c_str());
    //
    //        outPolig << "int npts = " << this->fpoligonalChain.NElements() << ";\n";
    //        outPolig << "fpoligonalChain.Resize(npts);\n";
    //        for(int p = 0; p < fpoligonalChain.NElements(); p++)
    //        {
    //            outPolig << "fpoligonalChain[" << p << "] = std::make_pair(" << this->fpoligonalChain[p].first << "," << this->fpoligonalChain[p].second << ");\n";
    //        }
    //        outPolig << "return;\n";
    //        outPolig.close();
    //    }
}
//------------------------------------------------------------------------------------------------------------

bool TPZPlaneFractureKernel::RemoveZigZag(TPZVec< std::pair<REAL,REAL> > &newPoligonalChain)
{
    bool thereWasZigZag = false;
    
    TPZVec< std::pair< REAL,REAL > > NOzigzagPoligonalChain(0);
    
    REAL v0x = 1.;
    REAL v0z = 0.;
    for(int p = 0; p < newPoligonalChain.NElements()-1; p++)
    {
        std::pair<REAL,REAL> p0 = newPoligonalChain[p];
        std::pair<REAL,REAL> p1 = newPoligonalChain[p+1];
        
        REAL v1x = p1.first - p0.first;
        REAL v1z = p1.second - p0.second;
        
        if(fabs(v1x) < 1.E-3 && fabs(v1z) < 1.E-3)
        {
            continue;
        }
        
        REAL innerProd = v0x*v1x + v0z*v1z;
        
        if(innerProd > 1.E-5)
        {
            int oldSize = NOzigzagPoligonalChain.NElements();
            NOzigzagPoligonalChain.Resize(oldSize+1);
            NOzigzagPoligonalChain[oldSize] = p0;
        }
        else
        {
            thereWasZigZag = true;
        }
        
        v0x = v1x;
        v0z = v1z;
    }
    int oldSize = NOzigzagPoligonalChain.NElements();
    NOzigzagPoligonalChain.Resize(oldSize+1);
    NOzigzagPoligonalChain[oldSize] = newPoligonalChain[newPoligonalChain.NElements()-1];
    
    newPoligonalChain = NOzigzagPoligonalChain;
    
    return thereWasZigZag;
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::RemoveLayerInvasion(TPZVec< std::pair<REAL,REAL> > &newPoligonalChain)
{
    int sz = newPoligonalChain.NElements();
    
    REAL maxZ = newPoligonalChain[0].second;
    REAL minZ = newPoligonalChain[sz-1].second;
    for(int p = 1; p < sz-1; p++)
    {
        if(newPoligonalChain[p].second > maxZ)
        {
            newPoligonalChain[p].second = maxZ;
        }
        if(newPoligonalChain[p].second < minZ)
        {
            newPoligonalChain[p].second = minZ;
        }
    }
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::TransferElasticSolution(REAL volAcum)
{
    std::cout << "\n\n\n************** TRANSFERINDO VOLUME PARA NOVA MALHA ELASTICA (PROPAGADA)\n";
    
    bool thereIsNegW;
    
    REAL newVol = 0.;
    REAL dtOrig = globTimeControl.actDeltaT();
    
    REAL Qequiv = -volAcum/60./fHbullet;
    std::map<int,TPZMaterial*>::iterator it;
    for(it = this->fmphysics->MaterialVec().begin(); it != this->fmphysics->MaterialVec().end(); it++)
    {
        if(globMaterialIdGen.IsBulletMaterial(it->first))
        {
            TPZBndCond * bnd = dynamic_cast<TPZBndCond*>(it->second);
            if(bnd)
            {
                bnd->Val2()(0,0) = Qequiv;
            }
        }
    }
    globTimeControl.SetDeltaT(60.);
    
    globLeakoffStorage.DisableLeakoff();
    this->RunThisFractureGeometry(volAcum, true);
    globLeakoffStorage.RestoreDefaultLeakoff();
    {
        std::map<int,TPZMaterial*>::iterator it;
        for(it = this->fmphysics->MaterialVec().begin(); it != this->fmphysics->MaterialVec().end(); it++)
        {
            if(globMaterialIdGen.IsBulletMaterial(it->first))
            {
                TPZBndCond * bnd = dynamic_cast<TPZBndCond*>(it->second);
                if(bnd)
                {
                    bnd->Val2()(0,0) = this->fQinj1wing_Hbullet;
                }
            }
        }
    }
    
    globTimeControl.SetDeltaT(dtOrig);
    
    REAL negVol = 0.;
    newVol = this->IntegrateW(thereIsNegW, negVol);
    
    if(fabs(volAcum - newVol) > 1.E-5)
    {
        std::cout << "\n\n\nW nao manteve volume na transferencia de solucao elastica!!!\n";
        std::cout << "volAntes = " << volAcum << std::endl;
        std::cout << "volDepois = " << newVol << std::endl;
        std::cout << "(newVol - volAcum) = " << (newVol - volAcum) << "\n\n\n";
        DebugStop();
    }
    else
    {
        std::cout << "\n\n\nTransferência OK!!!\n";
        PutConstantPressureOnFluidSolution();
    }
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::TransferLastLeakoff(TPZCompMesh * cmeshFrom)
{
#ifdef DEBUG
    REAL leakAcumBefore = this->ComputeVlAcumLeakoff(cmeshFrom);
#endif
    
    //Mapa que guardarah os volumes acumulados dos elementos de fratura da malha antiga que nao tem pais
    std::map<int,REAL> fatherGeoIndex_VolAcum, newLeakoffMap;

    std::map<int,REAL>::iterator it;
    
    //Guardando os volumes acumulados dos elementos de fratura da malha antiga que nao tem pais
    int nElem = cmeshFrom->NElements();
    for(int el = 0; el < nElem; el++)
    {
        TPZCompEl * cel = cmeshFrom->ElementVec()[el];
        TPZGeoEl * gel = cel->Reference();
        
#ifdef DEBUG
        if(!gel)
        {
            DebugStop();
        }
#endif
        
        if(globMaterialIdGen.IsInsideFractMat(gel->MaterialId()) == false)
        {
            continue;
        }
        
#ifdef DEBUG
        if(gel->Dimension() != 2)
        {
            DebugStop();
        }
#endif
        
        REAL Area = gel->SideArea(gel->NSides()-1);
        it = globLeakoffStorage.GetLeakoffMap().find(gel->Id());
        
#ifdef DEBUG
        if(it == globLeakoffStorage.GetLeakoffMap().end())
        {
            DebugStop();
        }
#endif
        
        REAL penetration = it->second;
        
        REAL volAcum = Area * penetration;
        
        //Procurando pelo elemento ancestral que pertence aa malha preservada.
        TPZGeoEl * preservedMeshGel = gel;
        while(preservedMeshGel->Father())
        {
            preservedMeshGel = preservedMeshGel->Father();
            if(this->fPlaneFractureMesh->GeoElementIsOnPreservedMesh(preservedMeshGel))
            {
                break;
            }
        }
        if(this->fPlaneFractureMesh->GeoElementIsOnPreservedMesh(preservedMeshGel) == false)
        {
            DebugStop();
        }
        
        it = fatherGeoIndex_VolAcum.find(preservedMeshGel->Index());
        if(it == fatherGeoIndex_VolAcum.end())
        {
            fatherGeoIndex_VolAcum[preservedMeshGel->Index()] = volAcum;
        }
        else
        {
            REAL volAcumOld = it->second;
            volAcum += volAcumOld;
            fatherGeoIndex_VolAcum[preservedMeshGel->Index()] = volAcum;
        }
    }
    
    //Distribuindo os volumes acumulados para os elementos de fratura da malha nova que nao tem filhos
    nElem = this->fmeshVec[1]->NElements();
    for(int el = 0; el < nElem; el++)
    {
        TPZCompEl * cel = this->fmeshVec[1]->ElementVec()[el];
        TPZGeoEl * gel = cel->Reference();
        
#ifdef DEBUG
        if(!gel)
        {
            DebugStop();
        }
#endif
        
        if(globMaterialIdGen.IsInsideFractMat(gel->MaterialId()) == false)
        {
            continue;
        }
        
#ifdef DEBUG
        if(gel->Dimension() != 2)
        {
            DebugStop();
        }
#endif
        
        //Procurando pelo elemento ancestral que pertence aa malha preservada.
        TPZGeoEl * preservedMeshGel = gel;
        while(preservedMeshGel->Father())
        {
            preservedMeshGel = preservedMeshGel->Father();
            if(this->fPlaneFractureMesh->GeoElementIsOnPreservedMesh(preservedMeshGel))
            {
                break;
            }
        }
        if(this->fPlaneFractureMesh->GeoElementIsOnPreservedMesh(preservedMeshGel) == false)
        {
            DebugStop();
        }
        
        it = fatherGeoIndex_VolAcum.find(preservedMeshGel->Index());
        if(it == fatherGeoIndex_VolAcum.end())
        {//Primeira vez que elemento pai vira interior de fratura
            continue;
        }
        
        REAL volAcum = it->second;
        
        if(preservedMeshGel == gel)
        {//Pelo fato do pai nao estar refinado, todo o volume acumulado vai para ele mesmo.
            REAL Area = gel->SideArea(gel->NSides()-1);
            REAL penetration = volAcum/Area;

            newLeakoffMap[gel->Id()] = penetration;
        }
        else
        {//Tem filhos, portanto precisa filtrar quais deles
         //eh interior de fratura e ratear o volume acumulado entre eles.
            std::vector<int> fractSons(0);
            REAL fractSonsArea = 0.;
            
            TPZVec<TPZGeoEl*> allSons(0);
            preservedMeshGel->GetHigherSubElements(allSons);
            for(int s = 0; s < allSons.NElements(); s++)
            {
                TPZGeoEl * son = allSons[s];
                if(globMaterialIdGen.IsInsideFractMat(son->MaterialId()))
                {
                    fractSons.push_back(s);
                    fractSonsArea += son->SideArea(son->NSides()-1);
                }
            }
            REAL penetration = volAcum/fractSonsArea;//penetracao igual a todos os filhos de dentro da fratura
            for(int s = 0; s < fractSons.size(); s++)
            {
                TPZGeoEl * fractSon = allSons[fractSons[s]];

                if(newLeakoffMap.find(fractSon->Id()) != newLeakoffMap.end())
                {
                    std::cout << "\nSobrescrevendo mapa newLeakoffMap!\n";
                    DebugStop();
                }

                newLeakoffMap[fractSon->Id()] = penetration;
            }
        }
        fatherGeoIndex_VolAcum.erase(it);
    }
    
    globLeakoffStorage.SetLeakoffMap(newLeakoffMap);
    
#ifdef DEBUG
    REAL leakAcumAfter = ComputeVlAcumLeakoff(this->fmeshVec[1]);
    if(fabs(leakAcumBefore - leakAcumAfter) > 1.E-5)
    {
        std::cout << "\n\nA transferencia de leakoff nao manteve volume!!!\n\n";
        std::cout << "leakAcumBefore = " << leakAcumBefore << "\n";
        std::cout << "leakAcumAfter = " << leakAcumAfter << "\n";
        DebugStop();
    }
#endif
}
//------------------------------------------------------------------------------------------------------------

//---------------------------------------------------------


BezierCurve::BezierCurve()
{
    DebugStop();//Utilize o outro construtor
}

BezierCurve::BezierCurve(TPZVec< std::pair< REAL,REAL > > &poligonalChain)
{
    this->falphaL = 1.;
    
    this->forder = poligonalChain.NElements()-1;
    this->fPoligonalChain = poligonalChain;
    REAL maxL = 0.;
    
    for(int i = 0; i < poligonalChain.NElements(); i++)
    {
        maxL = MAX(maxL,poligonalChain[i].first);
    }
    
    REAL funcMaxL = 0.;
    int npts = 100;
    std::pair< REAL,REAL > pt;
    for(int i = 0; i < npts; i++)
    {
        REAL t = i/(double(npts-1));
        this->F(t,pt);
        funcMaxL = MAX(funcMaxL,pt.first);
    }
    
    this->falphaL = maxL / funcMaxL;
}

BezierCurve::~BezierCurve()
{
    this->forder = 0;
    this->fPoligonalChain.Resize(0);
}

void BezierCurve::F(REAL t, std::pair< REAL,REAL > & ft)
{
    if(t < 1.E-5)
    {
        ft = this->fPoligonalChain[0];
        return;
    }
    else if(t > 0.99999)
    {
        ft = this->fPoligonalChain[forder];
        return;
    }
    else
    {
        REAL x = 0., z = 0.;
        for(int i = 0; i <= forder; i++)
        {
            REAL b = Bernstein(t,i);
            
            x += b * this->fPoligonalChain[i].first;
            z += b * this->fPoligonalChain[i].second;
        }
        ft = std::make_pair(this->falphaL * x,z);
    }
}

REAL BezierCurve::Bernstein(REAL t, REAL i)
{
    REAL b = Coef(forder,i,forder-i) * pow(1.-t,forder-i) * pow(t,i);
    return b;
}

REAL BezierCurve::Coef(int num1, int num2, int num3)
{
    REAL coefVal = 1.;
    while(num1 > 1 || num2 > 1 || num3 > 1)
    {
        num1 = MAX(num1,1);
        num2 = MAX(num2,1);
        num3 = MAX(num3,1);
        
        coefVal *= (double(num1)/(double(num2)*double(num3)));
        
        num1--;
        num2--;
        num3--;
    }
    
    return coefVal;
}