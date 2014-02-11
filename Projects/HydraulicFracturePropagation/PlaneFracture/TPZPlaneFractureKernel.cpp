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

#include "TPZVTKGeoMesh.h"


#define usingSWXGraphs

#ifdef usingSWXGraphs
#include "TSWXGraphMesh.h"
#include "TSWXGraphElement.h"
#endif

//Utilize 1. para output (Mathematica) em metros e 3.280829131 para output (Mathematica) em foot
const REAL feet = 1.;//3.280829131;

//Inicializando vetor de cores
const std::string TPZPlaneFractureKernel::color[12] = {"Red","Green","Blue","Black","Gray","Cyan","Magenta","Yellow","Brown","Orange","Pink","Purple"};


TPZPlaneFractureKernel::TPZPlaneFractureKernel() : actColor(0)
{
    DebugStop();//Use constructor below;
}
//------------------------------------------------------------------------------------------------------------

TPZPlaneFractureKernel::TPZPlaneFractureKernel(TPZVec<TPZLayerProperties> & layerVec, REAL bulletTVDIni, REAL bulletTVDFin,
                                               REAL xLength, REAL yLength, REAL Lmax, int nstripes, REAL Qinj_well, REAL visc,
                                               REAL Jradius,
                                               int pOrder,
                                               REAL MaxDispl_ini,
                                               REAL MaxDispl_fin,
                                               bool pressureIndependent,
                                               bool uncoupled) : actColor(0)
{
    if(nstripes < 1)
    {
        std::cout << "\nnstripes > 0 is needed\n";
        std::cout << "See " << __PRETTY_FUNCTION__ << ".\n";
        DebugStop();
    }
    
    this->fpoligonalChain.Resize(0);
    this->fstep = 0;
    
    this->fLmax = Lmax;
    
    this->fPlaneFractureMesh = new TPZPlaneFractureMesh(layerVec, bulletTVDIni, bulletTVDFin, xLength, yLength, Lmax, nstripes);
    
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
    this->fPoligonalChainInitialHeigh = 1.1 * (bulletTVDFin - fCenterTVD);
    
    this->fvisc = visc;
    
    this->fJIntegralRadius = Jradius;
    
    this->fpOrder = pOrder;
    
    this->fMaxDisplIni = MaxDispl_ini;
    this->fMaxDisplFin = MaxDispl_fin;
    
    this->fPath3D.Reset();
    
    if(pressureIndependent)
    {
        globLeakoffStorage.SetPressureIndependent();
    }
    else
    {
        globLeakoffStorage.SetPressureDependent();
    }
    this->fUncoupled = uncoupled;
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

void TPZPlaneFractureKernel::Run()
{
    {//Mathematica output (preamble for PostProcessFractGeometry method)
        std::ofstream outPoligMath("000MathPreamble.txt");
        outPoligMath << "(* colors = {Red,Green,Blue,Black,Gray,Cyan,Magenta,Yellow,Brown,Orange,Pink,Purple} *)\n";
        outPoligMath << "AllPolChains={};\n";
        outPoligMath << "Lgr={};\n";
        outPoligMath << "Hsupgr={};\n";
        outPoligMath << "Hinfgr={};";
    }
    
    this->InitializePoligonalChain();
    
    REAL fractVolum = 0.;
    TPZCompMesh * lastPressureCMesh = NULL;
    
    while(globTimeControl.ReachEndOftime() == false)
    {
        std::cout << "\n\n=============================================================================\n";
        std::cout << "STEP " << this->fstep << "\n";
        this->InitializeMeshes();

        if(this->fUncoupled)
        {
            REAL maxKI = 0.;
            REAL respectiveKIc = 0.;
            std::map< int, std::pair<REAL,REAL> > whoPropagate_KI;

            bool propagate = false;
            fractVolum = this->PredictFractVolume_WithNonNegativeUy();
            
            while(propagate == false)
            {
                maxKI = 0.;
                respectiveKIc = 0.;
                whoPropagate_KI.clear();

                propagate = this->CheckPropagationCriteria(maxKI, respectiveKIc, whoPropagate_KI);
                std::cout << "maxKI/respectiveKIc = " << maxKI/respectiveKIc << "\n\n";
                
                if(propagate == false)
                {
                    fractVolum *= 1.1;
                    this->TransferElasticSolution(fractVolum);
                }
            }

            fractVolum = IntegrateW();
            
            if(this->fstep > 0)
            {
                this->TransferLastLeakoff(lastPressureCMesh);
            }
            this->PutConstantPressureOnFluidSolution();
            this->PredictActDeltaT(fractVolum);
            
            {
                this->CloseActualTimeStepUncoupled();
                this->DefinePropagatedPoligonalChain(maxKI, respectiveKIc, whoPropagate_KI);
                this->fstep++;
            }
        }
        else
        {
            if(this->fstep > 0)
            {
                this->TransferElasticSolution(fractVolum);
                this->TransferLastLeakoff(lastPressureCMesh);
            }
            
            //Resolvendo o problema acoplado da nova geometria
            this->RunThisFractureGeometry(fractVolum);
        }
        
        lastPressureCMesh = this->fmeshVec[1];
    }//end of while(reachEndOfTime == false)
    
    std::ofstream outMath("ProstProcess.txt");
    globFractOutput3DData.PrintMathematica(outMath);
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::InitializePoligonalChain()
{
    //*********** TESTE ENGLAND-GREEN AQUICAJU
    /*
    int nptsUp = 2;
    fpoligonalChain.Resize(nptsUp);
    fpoligonalChain[0] = std::make_pair(0.5,-2110.);
    fpoligonalChain[1] = std::make_pair(69.5,-2110.);
     */
    
    /** Elipse */
    int nptsUp = 50;
    fpoligonalChain.Resize(0);
    
    REAL yc = -fCenterTVD;
    REAL sAx = 11.5;
    REAL sAy = sAx;
    
    for(int p = 1; p <= nptsUp; p++)
    {
        REAL x  = MIN(sAx - 0.001,p * sAx/nptsUp);
        REAL fx = yc + (sAy*sqrt(sAx*sAx - x*x))/sAx;
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
        if(x > this->fLmax)
        {
            int oldSize = fpoligonalChain.NElements();
            this->fpoligonalChain.Resize(oldSize+1);
            this->fpoligonalChain[oldSize] = std::make_pair(x,fx);
        }
    }
    
    
    /** CIRCULO */
    /*
    REAL Lmax = 0.5;
    int nsteps = M_PI * fPoligonalChainInitialHeigh / Lmax;
    if(nsteps < 10) nsteps = 10;
    REAL ang = M_PI / nsteps;
    
    int nnodes = nsteps + 1;
    fpoligonalChain.Resize(nnodes-2);
    
    for(int node = 1; node < nnodes-1; node++)
    {
        REAL vx = 0.1 + fPoligonalChainInitialHeigh*sin(node*ang);
        REAL vz = fPoligonalChainInitialHeigh*cos(node*ang) - fCenterTVD;
        fpoligonalChain[node-1] = std::make_pair(vx,vz);
    }
     */
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::InitializeMeshes()
{
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
    
    //Chute inicial
    this->fmeshVec[0]->Solution()(0,0) = 1.;
    for(int r = 1; r < this->fmeshVec[0]->Solution().Rows(); r++)
    {
        this->fmeshVec[0]->Solution()(r,0) = -(this->fPlaneFractureMesh->Max_MinCompressiveStress());
    }
    PutConstantPressureOnFluidSolution();
    
    this->InitializePath3DVector();
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

    {
        TPZSkylineStructMatrix skyl(cmesh); //caso simetrico
        an->SetStructuralMatrix(skyl);
        TPZStepSolver<REAL> stepS;
        stepS.SetDirect(ECholesky);
        an->SetSolver(stepS);
    }
    
//    {
//        TPZParFrontStructMatrix<TPZFrontSym<STATE> > skyl(cmesh);
//        skyl.SetQuiet(true);
//        skyl.SetNumberOfThreads(4);
//        an->SetStructuralMatrix(skyl);
//        TPZStepSolver<REAL> stepS;
//        stepS.SetDirect(ECholesky);
//        an->SetSolver(stepS);
//    }
    
    int NStripes = this->fPlaneFractureMesh->NStripes();

    TPZFMatrix<STATE> solution0(cmesh->Solution().Rows(), 1);
    TPZFMatrix<STATE> solutionStripe(cmesh->Solution().Rows(), 1);
    TPZFMatrix<STATE> solutions(cmesh->Solution().Rows(), 1+NStripes);
    
    TPZTimer ta, ts;
    ta.start();
    
    an->Assemble();
    
    ta.stop();
    ts.start();
    
    an->Solve();
    
    ts.stop();
    
    std::cout << "\nAssemble = " << ta.seconds() << " s\n";
    std::cout << "Solve = " << ts.seconds() << " s\n\n";
    
    solution0 = cmesh->Solution();
    for(int r = 0; r < solution0.Rows(); r++)
    {
        solutions(r,0) = solution0(r,0);
    }
    
    for(int stripe = 0; stripe < NStripes; stripe++)
    {
        this->fPlaneFractureMesh->SetSigmaNStripeNum(cmesh, stripe);
        
        an->Rhs().Zero();
        an->AssembleResidual();
        an->Solve();
        
        solutionStripe = cmesh->Solution();
        
        for(int r = 0; r < solution0.Rows(); r++)
        {
            solutions(r,stripe+1) = (solutionStripe(r,0) - solution0(r,0));
        }
    }
    cmesh->LoadSolution(solutions);
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
    
    REAL maxKI = 0.;
    REAL respectiveKIc = 0.;
    std::map< int, std::pair<REAL,REAL> > whoPropagate_KI;
    
    bool propagate = false;
    bool maxKIacceptable = false;
    
    long posBlock = -1;
    
    while (globTimeControl.TimeLimitsIsCloseEnough() == false && maxKIacceptable == false)
    {
        maxKI = 0.;
        respectiveKIc = 0.;
        whoPropagate_KI.clear();
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
        this->AssembleStiffMatrixLoadVec(an, matK, matRes_partial, posBlock);

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
            this->AssembleStiffMatrixLoadVec(an, matK, matRes_partial, posBlock);
            matRes_total = matRes_partial + matMass;

            normRes = Norm(matRes_total);
            std::cout << "||res|| = " << normRes << std::endl;
            nit++;
        }

        if(normRes > tol)
        {
            if(justTransferingElasticSolution)
            {
                std::cout << "\nNao convergiu na transferencia de solucao elastica...\n";
                DebugStop();
            }
            
            //Quando o actDeltaT leva a um instante em que Vleakoff = Vinj, nao converge, necessitando
            //trazer o limite esquerdo do deltaT da bisseccao para a direita
            globTimeControl.TimeisOnLeft();
            
            std::cout << "\nNao convergiu...\n";
            std::cout << "************* t está a Esquerda de KI=KIc *************\n";
        }
        else
        {
            if(justTransferingElasticSolution)
            {
                return;
            }
            
            bool thereWasNegativeW = false;
            REAL negVol = 0.;
            volAcum = this->IntegrateW(thereWasNegativeW, negVol);
            
            if(thereWasNegativeW)
            {
                globTimeControl.TimeisOnLeft();
                
                std::cout << "\nNegative W\n";
                std::cout << "************* t está a Esquerda de KI=KIc *************\n";
            }
            else
            {
                propagate = this->CheckPropagationCriteria(maxKI, respectiveKIc, whoPropagate_KI);
                
                if(propagate)
                {
                    globTimeControl.TimeisOnRight();
                    std::cout << "\nPropagate: maxKI/respectiveKIc = " << maxKI/respectiveKIc << "\n";
                    std::cout << "************* t está a Direita de KI=KIc *************\n";
                    REAL maxAlpha = 10.;//<<<<<<<<<<<<<<<<<<<<
                    maxKIacceptable = (maxKI/respectiveKIc < maxAlpha);
                }
                else
                {
                    globTimeControl.TimeisOnLeft();
                    std::cout << "KI < KIc\n";
                    std::cout << "************* t está a Esquerda de KI=KIc *************\n";
                }
            }
        }
        
        if(globTimeControl.TimeLimitsIsCloseEnough() == false && maxKIacceptable == false)
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
    
    this->DefinePropagatedPoligonalChain(maxKI, respectiveKIc, whoPropagate_KI);
    
    this->fstep++;
}
//------------------------------------------------------------------------------------------------------------

REAL TPZPlaneFractureKernel::PredictFractVolume_WithNonNegativeUy()
{
    bool thereIsNeg = true;
    
    this->TransferElasticSolution(0.);
    REAL negVol0 = 0.;
    REAL volFrac0 = this->IntegrateW(thereIsNeg, negVol0);
    
    this->TransferElasticSolution(0.01);
    REAL negVol1 = 0.;
    REAL volFrac1 = this->IntegrateW(thereIsNeg, negVol1);
    
    /////Metodo das secantes
    while(thereIsNeg)
    {
        REAL newVolFract = ((negVol0*volFrac1) - (negVol1*volFrac0))/(negVol0 - negVol1);
        
        volFrac0 = volFrac1;
        volFrac1 = newVolFract;
        
        negVol0 = negVol1;
        
        if(fabs(negVol0 - negVol1) < 1.E-5)
        {
            newVolFract *= 1.01;
            volFrac1 = newVolFract;
        }
        
        this->TransferElasticSolution(newVolFract);
        this->IntegrateW(thereIsNeg, negVol1);
    }
    
    return this->IntegrateW();
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
    
    while(fabs(dt1 - dtNext) > 1.E-1)
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
    
    REAL volumeBalance = MAX(fractVolum , volInj1 - volLeak1);
    this->TransferElasticSolution(volumeBalance);
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::CloseActualTimeStepUncoupled()
{
    globTimeControl.UpdateActTime();//atualizando esta rodada concluida para o tempo atual
    this->UpdateLeakoff();
    this->PostProcessAllPressureIndependent();
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::CloseActualTimeStepCoupled()
{
    globTimeControl.UpdateActTime();//atualizando esta rodada concluida para o tempo atual
    this->UpdateLeakoff();
    this->PostProcessAllPressureDependent();
    
    globTimeControl.RestartBissection();
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
        
        REAL Young = 0.;
        REAL KIc = 0.;
        this->fPlaneFractureMesh->GetYoung_and_KIc_FromLayerOfThisZcoord(zCoord, Young, KIc);
        
        this->fPath3D.PushBackPath3D( new Path3D(this->fmeshVec[0], this->fmeshVec[1],
                                                 center, Young, KIc, normal, this->fJIntegralRadius) );
    }
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::AssembleStiffMatrixLoadVec(TPZAnalysis *an,
                                                        TPZAutoPointer< TPZMatrix<REAL> > & matK, TPZFMatrix<REAL> & matRes,
                                                        long &posBlock)
{
	this->fPlaneFractureMesh->SetActualState();
    
    TPZFStructMatrix structMat(this->fmphysics);
	an->SetStructuralMatrix(structMat);
    
    //Filtrando equacao de alpha0
    bool mustFilter = true;
    if(mustFilter)
    {
        int neq = this->fmphysics->NEquations();
        long blockAlphaEslast = this->fmphysics->ConnectVec()[0].SequenceNumber();
        posBlock = this->fmphysics->Block().Position(blockAlphaEslast);
        
        TPZVec<long> actEquations(neq-1);
        int p = 0;
        for(long eq = 0; eq < neq; eq++)
        {
            if(eq != posBlock)
            {
                actEquations[p] = eq;
                p++;
            }
        }
        TPZEquationFilter eqF(neq);
        eqF.SetActiveEquations(actEquations);
        
        an->StructMatrix()->EquationFilter() = eqF;
        this->fmphysics->Solution()(posBlock,0) = 1.;//<<< Nao Mexer!!! Eh o alpha0
    }
    
	TPZStepSolver<REAL> stepS;
	stepS.SetDirect(ELU);
	an->SetSolver(stepS);
    
    an->Assemble();
    
    matK = an->Solver().Matrix();
	matRes = an->Rhs();
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
    
    globTimeControl.ComputeActDeltaT();
    
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
    
    long posBlock = -1;
    AssembleStiffMatrixLoadVec(an, fL_xIni, f_xIni, posBlock);
    
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
        this->AssembleStiffMatrixLoadVec(an, fLtemp, fExato_x, posBlock);
        
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
    REAL volInjected = this->ComputeVolInjected();
    
    REAL volW = this->IntegrateW();
    
    REAL volLeakoffExpected = volInjected - volW;
    REAL volLeakoffComputed = ComputeVlAcumLeakoff(this->fmeshVec[1]);
    
    if(fabs(volLeakoffComputed - volLeakoffExpected) > 0.1)
    {
        std::cout << "\n\n\nLeakoff não está sendo calculado corretamente!\n";
        std::cout << "Computed = " << volLeakoffComputed << "\n";
        std::cout << "volInjected = " << volInjected << "\n";
        std::cout << "volW = " << volW << "\n";
        std::cout << "Expected = volInjected - volW = " << volLeakoffExpected << "\n\n\n";
        DebugStop();
    }
    if(fabs(volLeakoffComputed) > 1.E-10)
    {
        REAL alpha = volLeakoffExpected/volLeakoffComputed;
        std::map<int,REAL>::iterator itLk;
        for(itLk = globLeakoffStorage.GetLeakoffMap().begin(); itLk != globLeakoffStorage.GetLeakoffMap().end(); itLk++)
        {
            itLk->second *= alpha;
        }
    }
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::PostProcessAllPressureIndependent()
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

void TPZPlaneFractureKernel::PostProcessAllPressureDependent()
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
    REAL wAcum = this->IntegrateW();
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
#ifdef usingSWXGraphs
    TSWXGraphMesh grMesh;
    TSWXGraphElement grEl(0);
    TPZVec<std::string> nodalSol(2), cellSol(0);
    nodalSol[0] = "Displacement";
    nodalSol[1] = "StressY";

    grEl.GenerateVTKData(this->fmeshVec[0], 3, 0., nodalSol, cellSol, grMesh);
    
    std::stringstream nm;
    nm << "Elasticity_Step" << this->fstep << ".vtk";
    std::ofstream file(nm.str().c_str());
    grMesh.ToParaview(file);
#endif
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::PostProcessPressure()
{
#ifdef usingSWXGraphs
    TSWXGraphMesh grMesh;
    TSWXGraphElement grEl(0);
    TPZVec<std::string> nodalSol(1), cellSol(0);
    nodalSol[0] = "Pressure";
    grEl.GenerateVTKData(this->fmeshVec[1], 2, 0., nodalSol, cellSol, grMesh);
    
    std::stringstream nm;
    nm << "Pressure_Step" << this->fstep << ".vtk";
    std::ofstream file(nm.str().c_str());
    grMesh.ToParaview(file);
#endif
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
    
    {//Mathematica output
        std::stringstream nmMath, nmAux;
        nmAux << "pcm={";
        nmMath << "PoligonalChainMath_Step" << this->fstep << ".txt";
        std::ofstream outPoligMath(nmMath.str().c_str());
        
        for(int p = 0; p < fpoligonalChain.NElements(); p++)
        {
            outPoligMath << "fractureDots" << p << " = {" << feet * fpoligonalChain[p].first << ","
                                                          << feet * fpoligonalChain[p].second + fCenterTVD << /*","
                                                          << globTimeControl.actTime()/60. << */"};\n";
            nmAux << "fractureDots" << p;
            if(p < fpoligonalChain.NElements()-1)
            {
                nmAux << ",";
            }
        }
        nmAux << "};\n";
        nmAux << "gr" << this->fstep << " = ListPlot[pcm, Joined -> True, AxesOrigin -> {0,0}, AspectRatio -> 2,PlotStyle->" << color[actColor%12] << "];\n";
        nmAux << "L" << this->fstep << " = Max[Transpose[pcm][[1]]];\n";
        nmAux << "Hsup" << this->fstep << " = Max[Transpose[pcm][[2]]];\n";
        nmAux << "Hinf" << this->fstep << " = -Min[Transpose[pcm][[2]]];\n";
        nmAux << "AppendTo[AllPolChains,gr" << this->fstep << "];\n";
        nmAux << "AppendTo[Lgr,{" << globTimeControl.actTime()/60. << ",L" << this->fstep << "}];\n";
        nmAux << "AppendTo[Hsupgr,{" << globTimeControl.actTime()/60. << ",Hsup" << this->fstep << "}];\n";
        nmAux << "AppendTo[Hinfgr,{" << globTimeControl.actTime()/60. << ",Hinf" << this->fstep << "}];\n";
        nmAux << "Print[\"L" << this->fstep << " = \", L" << this->fstep << "]\n";
        nmAux << "Print[\"Hsup" << this->fstep << " = \", Hsup" << this->fstep << "]\n";
        nmAux << "Print[\"Hinf" << this->fstep << " = \", Hinf" << this->fstep << "]\n";
        outPoligMath << nmAux.str();
        actColor++;
    }
}
//------------------------------------------------------------------------------------------------------------

REAL TPZPlaneFractureKernel::IntegrateW()
{
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
        
        integralW += 2.*value[1];//Integrando w = (2*uy)
    }
    
    return integralW;
}
//------------------------------------------------------------------------------------------------------------

REAL TPZPlaneFractureKernel::IntegrateW(bool & thereWasNegativeW, REAL &negVol)
{
    thereWasNegativeW = false;
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
            negVol += 2.*value[1];
            thereWasNegativeW = true;
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

bool TPZPlaneFractureKernel::CheckPropagationCriteria(REAL &maxKI, REAL &respectiveKIc,
                                                      std::map< int, std::pair<REAL,REAL> > &whoPropagate_KI)
{
    bool propagate = false;
    
    maxKI = 0.;
    respectiveKIc = 0.;
    whoPropagate_KI.clear();
    
    //Calculo das integrais-J
    this->fPath3D.IntegratePath3D();
    
    //Calculo de KI
    for(int p = 0; p < fPath3D.NPaths(); p++)
    {
        REAL cracktipKI = fPath3D.Path(p)->KI();
        REAL cracktipKIc = fPath3D.Path(p)->KIc();
        
        {
            TPZVec<REAL> originVec = this->fPath3D.Path(p)->Origin();
            TPZVec<REAL> Jdirection = this->fPath3D.Path(p)->JDirection();
            TPZVec<REAL> ptTemp(3,0.), qsiTemp(2,0.);
            ptTemp[0] = originVec[0] + 0.1 * Jdirection[0];
            ptTemp[2] = originVec[2] + 0.1 * Jdirection[2];

            TPZGeoEl * gel = fPlaneFractureMesh->Find2DElementNearCrackTip(p, ptTemp);
            if( !gel || globMaterialIdGen.IsInsideFractMat(gel->MaterialId()) )
            {
                //Nao serah propagado para pontos fora do dominio
                //ou no interior da fratura (i.e., fechamento).
                cracktipKI = 0.;
            }
        }
        if(cracktipKI >= cracktipKIc)
        {
            propagate = true;
            whoPropagate_KI[p] = std::make_pair(cracktipKI,cracktipKIc);
            
            if(cracktipKI > maxKI)
            {
                maxKI = cracktipKI;
                respectiveKIc = cracktipKIc;
            }
        }
    }
    
    return propagate;
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::DefinePropagatedPoligonalChain(REAL maxKI, REAL respectiveKIc,
                                                            std::map< int, std::pair<REAL,REAL> > &whoPropagate_KI)
{
    TPZVec< std::pair<REAL,REAL> > newPoligonalChain(0);
    
    std::map< int, std::pair<REAL,REAL> >::iterator it;
    
    REAL newX = 0.;
    REAL newZ = 0.;
    for(int p = 0; p < this->fPath3D.NPaths(); p++)
    {
        it = whoPropagate_KI.find(p);
        if(it == whoPropagate_KI.end())
        {//Nao propagou, portanto serah mantido o ponto medio
            TPZVec<REAL> originVec = this->fPath3D.Path(p)->Origin();
            
            newX = originVec[0];
            newZ = originVec[2];
        }
        else
        {//Propagou, portanto serah definido novo ponto a partir de seu centro
            REAL KI = it->second.first;
            REAL KIc = it->second.second;
            
            TPZVec<REAL> originVec = this->fPath3D.Path(p)->Origin();
            TPZVec<REAL> Jdirection = this->fPath3D.Path(p)->JDirection();
            
            //Variacao linear no decorrer do tempo de fMaxDispl para fMinDispl
            REAL dLmax = this->fMaxDisplIni + (globTimeControl.actTime()*(this->fMaxDisplFin - this->fMaxDisplIni))/globTimeControl.Ttot();
            REAL alpha = 1.;//alphaMin = [1.0~2.0]
            
            REAL displacement = dLmax * pow((KI - KIc)/(maxKI - KIc),alpha);
            
            newX = originVec[0] + displacement * Jdirection[0];
            newZ = originVec[2] + displacement * Jdirection[2];
        }
        if(newX > fLmax/2.)
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
    
//    {//C++ tool : AQUICAJU
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
        
        if(innerProd > -1.E-5)
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

void TPZPlaneFractureKernel::TransferElasticSolution(REAL volAcum)
{
    if(this->fUncoupled == false)
    {
        std::cout << "\n\n\n************** TRANSFERINDO VOLUME PARA NOVA MALHA ELASTICA (PROPAGADA)\n";
    }
    
    if(this->fPlaneFractureMesh->NStripes() == 1)
    {
        int rows = this->fmeshVec[0]->Solution().Rows();
        int cols = this->fmeshVec[0]->Solution().Cols();
        if(cols != 1)
        {
            std::cout << "\n\n\nSoh era pra ter 1 coluna\n";
            std::cout << "Ver " << __PRETTY_FUNCTION__ << ".\n";
            DebugStop();
        }
        
        //Volume para alpha0 = 1 e alphas1 = 0
        TPZFMatrix<REAL> wSol(rows, cols);
        wSol.Zero();
        wSol(0,0) = 1.;
        this->fmeshVec[0]->LoadSolution(wSol);
        REAL volAlpha1_ini = this->IntegrateW();
        
        //Volume para alpha0 = 1 e alphas1 != 0
        for(int r = 1; r < rows; r++)
        {
            wSol(r,0) = 1.;
        }
        this->fmeshVec[0]->LoadSolution(wSol);
        REAL volAlpha1_fin = this->IntegrateW();
        
        //Calculando alpha para as stripes
        REAL alpha = (volAcum - volAlpha1_ini)/(volAlpha1_fin - volAlpha1_ini);
        for(int r = 1; r < rows; r++)
        {
            wSol(r,0) *= alpha;
        }
        this->fmeshVec[0]->LoadSolution(wSol);
    }
    else
    {
        REAL dtOrig = globTimeControl.actDeltaT();
     
        //DeltaT que, sem leakoff, deixa a fratura atual com o mesmo volume fornecido como parametro
        REAL dtEquiv = -volAcum/(this->fQinj1wing_Hbullet*this->fHbullet);
        globTimeControl.SetDeltaT(dtEquiv);
        
        globLeakoffStorage.DisableLeakoff();
        this->RunThisFractureGeometry(volAcum, true);
        globLeakoffStorage.RestoreDefaultLeakoff();
        
        globTimeControl.SetDeltaT(dtOrig);
        
//        REAL newVolAcum = this->IntegrateW();
//        REAL alpha = volAcum/newVolAcum;
//        
//        //Sintonia fina (corrigindo eveltuais discrepancias entre volume antes e depois)
//        if(fabs(volAcum - newVolAcum) > 1.E-10)
//        {
//            TPZFMatrix<REAL> newSolution = fmeshVec[0]->Solution();
//            for(int r = 0; r < newSolution.Rows(); r++)
//            {
//                for(int c = 0; c < newSolution.Cols(); c++)
//                {
//                    newSolution(r,c) *= alpha;
//                }
//            }
//            fmeshVec[0]->LoadSolution(newSolution);
//        }
    }
    
    //Verificando se a correcao deu certo
    REAL newVolAcum = this->IntegrateW();
    if(fabs(volAcum - newVolAcum) > 1.E-5)
    {
        std::cout << "\n\n\nW nao manteve volume na transferencia de solucao elastica!!!\n";
        std::cout << "volAntes = " << volAcum << std::endl;
        std::cout << "volDepois = " << newVolAcum << std::endl;
        std::cout << "(newVolAcum - volAcum) = " << (newVolAcum - volAcum) << "\n\n\n";
        DebugStop();
    }
    else
    {
        if(this->fUncoupled == false)
        {
            std::cout << "\n\n\nTransferência OK!!!\n";
        }
        TPZBuildMultiphysicsMesh::TransferFromMeshes(this->fmeshVec, this->fmphysics);
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

void TPZPlaneFractureKernel::PutConstantPressureOnFluidSolution()
{
    for(int r = 0; r < this->fmeshVec[1]->Solution().Rows(); r++)
    {
        this->fmeshVec[1]->Solution()(r,0) = this->fmeshVec[0]->Solution()(1,0);
    }
    TPZBuildMultiphysicsMesh::TransferFromMeshes(this->fmeshVec, this->fmphysics);
}


//---------------------------------------------------------


BezierCurve::BezierCurve()
{
    DebugStop();//Utilize o outro construtor
}

BezierCurve::BezierCurve(TPZVec< std::pair< REAL,REAL > > &poligonalChain)
{
    this->forder = poligonalChain.NElements()-1;
    this->fPoligonalChain = poligonalChain;
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
        REAL alpha = 1.;
        ft = std::make_pair(alpha * x,z);
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