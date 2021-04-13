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
#include "TPZGeoLinear.h"
#include "TPZVTKGeoMesh.h"

#include "TSWXGraphMesh.h"
#include "TSWXGraphElement.h"
#include "pzvtkmesh.h"

#include <algorithm>



TPZPlaneFractureKernel::TPZPlaneFractureKernel()
{
    DebugStop();//Use constructor below;
}
//------------------------------------------------------------------------------------------------------------

TPZPlaneFractureKernel::TPZPlaneFractureKernel(TPZVec<LayerProperties> & layerVec, REAL bulletTVDIni, REAL bulletTVDFin,
                                               REAL xLength, REAL yLength, REAL Lmax, REAL Qinj_well, REAL visc,
                                               REAL Jradius,
                                               REAL MaxDispl,
                                               bool pressureIndependent,
                                               bool just1stripe,
                                               bool layersStripesToo)
{
    if(Jradius > 0.7 * Lmax)
    {
        std::cout << "\n\n\nJradius > 0.7 * Lmax. Its value was changed from " << Jradius << " to " << 0.7 * Lmax << "\n\n\n";
        Jradius = 0.7 * Lmax;
    }
    
    this->fpoligonalChain.Resize(0);
    this->fstep = 0;
    
    this->fLmax = Lmax;
    this->fHbullet = (bulletTVDFin - bulletTVDIni);
    
    globLayerStruct.SetData(layerVec, bulletTVDFin, this->fHbullet);
    
    this->fPlaneFractureMesh = new TPZPlaneFractureMesh(bulletTVDIni, bulletTVDFin, xLength, yLength, Lmax, just1stripe, layersStripesToo);
    
    this->fmeshVec.Resize(2);
    this->fmeshVec[0] = NULL;
    this->fmeshVec[1] = NULL;
    
    this->fmphysics = NULL;
    
    REAL Qinj1wing = Qinj_well/2.;
    
    {
        globFractOutput3DData.SetQinj1wing(-Qinj1wing);
    }
    
    this->fQinj1wing_Hbullet = Qinj1wing/this->fHbullet;
    
    this->fCenterTVD = (bulletTVDIni + bulletTVDFin)/2.;
    
    this->fvisc = visc;
    
    this->fJIntegralRadius = Jradius;
    
    this->fMaxDispl = MaxDispl;
    
    this->fPath3D.Reset();
    
    if(pressureIndependent)
    {
        globLeakoffStorage.SetPressureIndependent();
        globLayerStruct.SetStaticPressureAsPrestress();
    }
    else
    {
        globLeakoffStorage.SetPressureDependent();
    }
    
    this->fResTop = -layerVec[0].fTVDini;
    this->fResBottom = -layerVec[layerVec.NElements()-1].fTVDfin;
    this->fResRight = xLength;
    
    this->fjust1stripe = just1stripe;
    this->flayersStripesToo = layersStripesToo;
    
    this->fZbulletTop = -bulletTVDIni;
    this->fZbulletBottom = -bulletTVDFin;
}
//------------------------------------------------------------------------------------------------------------

TPZPlaneFractureKernel::~TPZPlaneFractureKernel()
{
    delete this->fPlaneFractureMesh;
    this->fPlaneFractureMesh = NULL;
    
    for(int m = 0; m < this->fmeshVec.NElements(); m++)
    {
        delete this->fmeshVec[m];
        this->fmeshVec[m] = NULL;
    }
    this->fmeshVec.Resize(0);
    
    delete this->fmphysics;
    this->fmphysics = NULL;
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::Run()
{
    this->InitializePoligonalChain();
    
    TPZCompMesh * lastElasticCMesh = NULL;
    TPZCompMesh * lastPressureCMesh = NULL;
    
    while(globTimeControl.ReachEndOftime() == false)
    {
        std::cout << "\n\n=============================================================================\n";
        std::cout << "STEP " << this->fstep << "\n";
        this->InitializeMeshes(lastElasticCMesh);
        
        if(this->fstep > 0)
        {
            this->TransferLastLeakoff(lastPressureCMesh);
        }
        
        //Resolvendo o problema acoplado da nova geometria
        this->RunThisFractureGeometry();
        
        if(lastElasticCMesh)
        {
            TPZCompMeshReferred * cmeshRef = dynamic_cast<TPZCompMeshReferred*>(lastElasticCMesh);
            if(cmeshRef)
            {
                delete cmeshRef->ReferredMesh();
            }
            else
            {
                DebugStop();
            }
            delete lastElasticCMesh;
        }
        lastElasticCMesh = this->fmeshVec[0];
        
        if(lastPressureCMesh)
        {
            delete lastPressureCMesh;
        }
        lastPressureCMesh = this->fmeshVec[1];
        
        this->fPlaneFractureMesh->DeleteLastRefinedMesh();
        
        globTimeControl.RestartBissection();
    }
    std::cout << "\n\n\n\n\n********** THE END **********\n\n\n";
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::InitializePoligonalChain()
{
    ///Elipse
    int nptsUp = 50;
    fpoligonalChain.Resize(0);
    
    REAL yc = -fCenterTVD;
    
    REAL shiftX = 2.;
    REAL sAx = 4.5 - shiftX;
    REAL sAy = this->fHbullet/2. + 0.2;//Nao posso ter altura de fratura menor que hbullet, senao a BC de Qinj nao contribuira totalmente.
    
    for(int p = 1; p <= nptsUp; p++)
    {
        REAL x  = MIN(sAx - 0.001,p * sAx/nptsUp);
        REAL fx = yc + (sAy*sqrt(sAx*sAx - x*x))/sAx;
        x += shiftX;
        
        int oldSize = fpoligonalChain.NElements();
        this->fpoligonalChain.Resize(oldSize+1);
        this->fpoligonalChain[oldSize] = std::make_pair(x,fx);
    }
    for(int p = nptsUp-1; p > 0; p--)
    {
        REAL x  = MIN(sAx - 0.001,p * sAx/nptsUp);
        REAL fx = yc - (sAy*sqrt(sAx*sAx - x*x))/sAx;
        x += shiftX;
        
        int oldSize = fpoligonalChain.NElements();
        this->fpoligonalChain.Resize(oldSize+1);
        this->fpoligonalChain[oldSize] = std::make_pair(x,fx);
    }
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::InitializeMeshes(TPZCompMesh * lastElastCMesh)
{
    //GeoMesh
    this->fPlaneFractureMesh->InitializeFractureGeoMesh(fpoligonalChain);
    
    std::cout << "\n************** GERANDO MALHAS COMPUTACIONAIS\n";
    
    //Malha computacional elastica, processada com newman aplicado em toda a fratura e por faixas
    TPZCompMesh * fractureCMeshRef = this->fPlaneFractureMesh->GetFractureCompMesh();
    this->ProcessLinearElasticCMesh(fractureCMeshRef);
    
    //Malha computacional do tipo CMeshReferred
    this->fmeshVec[0] = this->fPlaneFractureMesh->GetFractureCompMeshReferred(fractureCMeshRef);
    
    //Malha computacional de pressao
    this->fmeshVec[1] = this->fPlaneFractureMesh->GetPressureCompMesh();
    
    //Malha computacional de acoplamento (multifisica)
    if(this->fmphysics)
    {
        delete this->fmphysics;
    }
    this->fmphysics = this->fPlaneFractureMesh->GetMultiPhysicsCompMesh(this->fmeshVec, lastElastCMesh,
                                                                        this->fQinj1wing_Hbullet,
                                                                        this->fvisc);
    
    this->ApplyInitialCondition(1.*globLayerStruct.GetHigherPreStress());
    
    this->InitializePath3DVector();
    
    //    {
    //        this->fmeshVec[0]->Solution().Zero();
    //        this->fmeshVec[0]->Solution()(0,0) = 1.;
    //
    //        this->PostProcessElasticity(0);
    //
    //        for(int s = 0; s < this->fPlaneFractureMesh->NStripes(); s++)
    //        {
    //            this->fmeshVec[0]->Solution().Zero();
    //            this->fmeshVec[0]->Solution()(0,0) = 1.;
    //            this->fmeshVec[0]->Solution()(s+1,0) = 1.;
    //
    //            this->PostProcessElasticity(s);
    //        }
    //        DebugStop();
    //    }
}
//------------------------------------------------------------------------------------------------------------

#include "TPZTimer.h"
#define Skyline
void TPZPlaneFractureKernel::ProcessLinearElasticCMesh(TPZCompMesh * cmesh)
{
    std::cout << "\n************** CALCULANDO SOLUCOES ELASTICAS DE REFERENCIA\n";
    
    {
        int neq = cmesh->NEquations();
        std::cout << "\nNequacoes elastica 3D = " << neq << "\n";
    }
    
    TPZAnalysis * an = new TPZAnalysis(cmesh, true);
    
#ifdef Skyline
    TPZSkylineStructMatrix skyl(cmesh); //caso simetrico
#else
    TPZParFrontStructMatrix< TPZFrontSym<STATE> > skyl(cmesh);
    skyl.SetNumThreads(4);//<<<<< NAO ESTA FUNCIONANDO!!!
#endif
    
    an->SetStructuralMatrix(skyl);
    TPZStepSolver<STATE> stepS;
    stepS.SetDirect(ECholesky);
    an->SetSolver(stepS);
    
    TPZFMatrix<STATE> solution0(cmesh->Solution().Rows(), 1);
    TPZFMatrix<STATE> solution1(cmesh->Solution().Rows(), 1);
    //
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
    
    ////////////////// Newman=1. em toda a fratura
    this->fPlaneFractureMesh->SetNeumannOnEntireFracture(cmesh);
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
    
    if(this->fjust1stripe == false && this->flayersStripesToo == false)
    {
        ////////////////// Newman=1. em cada faixa da fratura quando nao temos layersStripesToo
        int NStripes = this->fPlaneFractureMesh->NStripes();
        for(int stripe = 0; stripe < NStripes; stripe++)
        {
            this->fPlaneFractureMesh->SetNeumannOnThisStripe(cmesh,stripe);
            
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
    }
    if(this->flayersStripesToo == true)
    {
        ////////////////// Newman=1. em cada layer da fratura, podendo tb ser separado por faixas
        int NStripes = this->fPlaneFractureMesh->NStripes();
        if(this->fjust1stripe == true && NStripes != 1)
        {
            DebugStop();
        }
        for(int lay = 0; lay < globLayerStruct.NLayers(); lay++)
        {
            for(int stripe = 0; stripe < NStripes; stripe++)
            {
                bool neumannApplied = this->fPlaneFractureMesh->SetNeumannOnThisLayerAndStripe(cmesh,lay,stripe);
                if(neumannApplied)
                {
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
            }
        }
    }
    
    cmesh->LoadSolution(solutions);
    
    delete an;
}
//------------------------------------------------------------------------------------------------------------


void TPZPlaneFractureKernel::RunThisFractureGeometry()
{
    globFractOutput3DData.fKI_KI_history << "STEP " << this->fstep << ":\n";
    globFractOutput3DData.fKI_KI_history << "actTime = " << globTimeControl.actTime() << "s\n\n";
    std::map<REAL,std::string> deltaT_historyMap;
    
    TPZAnalysis * an = new TPZAnalysis(this->fmphysics);
    
    /** Convergence test */
//    CheckConv();
    /**********************/
    
    int nEq = this->fmphysics->NEquations();
    
    TPZFMatrix<STATE> matRes_total(nEq,1,0.);
    
    TPZFMatrix<STATE> Sol_0 = this->fmphysics->Solution();
    
    TPZAutoPointer< TPZMatrix<STATE> > matK;
    TPZFMatrix<STATE> matRes_partial(nEq,1,0.);
    TPZFMatrix<STATE> matMass(nEq,1,0.);
    
    /// Backup
    TPZFMatrix<STATE> backupSol_0 = Sol_0;
    
    std::cout << "\n\n\n************** CALCULANDO SOLUCAO ACOPLADA (KERNEL)\n\n";
    
    REAL maxKI_KIc = 0.;
    std::set<int> whoPropagate;
    
    bool timeLimitsIsCloseEnough = false;
    bool propagate = false;
    
    while(timeLimitsIsCloseEnough == false)
    {
        maxKI_KIc = 0.;
        whoPropagate.clear();
        propagate = false;
        
        globTimeControl.ComputeActDeltaT();
        
        std::cout << "\n\ndtLeft = " << globTimeControl.LeftDeltaT() << "s " <<
        " : actDeltaT = " << globTimeControl.actDeltaT() << "s " <<
        " : dtRight = " << globTimeControl.RightDeltaT() << "s\n";
        
        this->MassMatrix(matMass);
        
        EWhoBlock whoBlock = EBlockStripes;//set what degree of freedom is locked on equation filter
        
        this->AssembleStiffMatrixLoadVec(an, matK, matRes_partial, whoBlock);
        
        matRes_total = matRes_partial + matMass;
        
        ///Metodo de Newton
        std::cout << "\n-> Método de Newton\n";
        
        REAL normRes = 1.;
        REAL stripeTol = 1.E-1;
        REAL tolRes = 1.E-1;
        int maxit = 15;
        int nit = 0;
        while(normRes > tolRes && nit < maxit)
        {
            if(nit == 0 && whoBlock == EBlockStripes)
            {
                std::cout << "---> Entire fracture:\n";
            }
            else if(nit == 0 && whoBlock == EBlockEntireFracure)
            {
                if(this->fjust1stripe == true && this->flayersStripesToo == true)
                {
                    std::cout << "\n---> Layers:\n";
                }
                else if(this->fjust1stripe == false && this->flayersStripesToo == true)
                {
                    std::cout << "\n---> Layers and Stripes:\n";
                }
                else if(this->fjust1stripe == false && this->flayersStripesToo == false)
                {
                    std::cout << "\n---> Stripes:\n";
                }
                else
                {
                    std::cout << "EBlockEntireFracure but (this->fjust1stripe == true && this->flayersStripesToo == false)!!!\n\n";
                    DebugStop();
                }
            }
            
            an->Rhs() = matRes_total;
            an->Solve();
            
            TPZFMatrix<STATE> Sol_1_minus_Sol_0 = an->Solution();
            TPZFMatrix<STATE> Sol_1 = Sol_1_minus_Sol_0 + Sol_0;
            
            an->LoadSolution(Sol_1);
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(this->fmeshVec, this->fmphysics);
            
            Sol_0 = an->Solution();
            
            this->AssembleStiffMatrixLoadVec(an, matK, matRes_partial, whoBlock);
            
            matRes_total = matRes_partial + matMass;
            
            normRes = Norm(matRes_total);
            
            if((this->fjust1stripe == false || this->flayersStripesToo == true) && (normRes <= stripeTol && whoBlock == EBlockStripes))
            {
                maxit = 10;
                nit = 0;
                whoBlock = EBlockEntireFracure;
                this->AssembleStiffMatrixLoadVec(an, matK, matRes_partial, whoBlock);
                matRes_total = matRes_partial + matMass;
                
                std::cout << "||res|| = " << normRes << std::endl;
                normRes = 1.;
            }
            else
            {
                std::cout << "||res|| = " << normRes << std::endl;
                nit++;
            }
        }///end of Newton
        if(normRes > tolRes)
        {
            globTimeControl.TimeisOnLeft();
            
            std::cout << "\nNao convergiu (provavelmente w muito pequeno resultando em grande gradiente de pressão)...\n";
            std::cout << "************* dt está a Esquerda de KI=KIc *************\n";
            
            {
                std::ostringstream strVal;
                strVal << globTimeControl.actDeltaT();
                std::string tempstr = "deltaT ";
                tempstr += strVal.str();
                tempstr += " s : Não convergiu\n";
                deltaT_historyMap[globTimeControl.actDeltaT()] = tempstr;
            }
        }
        else
        {
            bool thereIsnegNetPress = this->ThereIsNegativeNetPressure();
            
            bool thereIsNegW;
            REAL negVol = 0.;
            this->IntegrateW(thereIsNegW, negVol);
            
            if(thereIsnegNetPress)
            {
                globTimeControl.TimeisOnLeft();
                std::cout << "\nNegative Net Pressure\n";
                std::cout << "************* dt está a Esquerda de KI=KIc *************\n";
                {
                    std::ostringstream strVal, strVol;
                    strVal << globTimeControl.actDeltaT();
                    strVol << negVol;
                    std::string tempstr = "deltaT ";
                    tempstr += strVal.str();
                    tempstr += " s : Net pressure negativa\n";
                    deltaT_historyMap[globTimeControl.actDeltaT()] = tempstr;
                }
            }
            else if(thereIsNegW)
            {
                globTimeControl.TimeisOnLeft();
                std::cout << "\nNegativeW = " << negVol << "\n";
                std::cout << "************* dt está a Esquerda de KI=KIc *************\n";
                {
                    std::ostringstream strVal, strVol;
                    strVal << globTimeControl.actDeltaT();
                    strVol << negVol;
                    std::string tempstr = "deltaT ";
                    tempstr += strVal.str();
                    tempstr += " s : W negativo = ";
                    tempstr += strVol.str();
                    tempstr += " m3\n";
                    deltaT_historyMap[globTimeControl.actDeltaT()] = tempstr;
                }
            }
            else
            {
                int maxKI_KIcPosition = -1;
                propagate = this->CheckPropagationCriteria(maxKI_KIc, whoPropagate, maxKI_KIcPosition);
                if(propagate == false)
                {
                    globTimeControl.TimeisOnLeft();
                    std::cout << "KI < KIc\n";
                    std::cout << "************* dt está a Esquerda de KI=KIc *************\n";
                    
                    {
                        std::ostringstream strVal, strKI;
                        strVal << globTimeControl.actDeltaT();
                        strKI << maxKI_KIc;
                        std::string tempstr = "deltaT ";
                        tempstr += strVal.str();
                        tempstr += " s : KI < KIc (relativo = ";
                        tempstr += strKI.str();
                        tempstr += ")\n";
                        deltaT_historyMap[globTimeControl.actDeltaT()] = tempstr;
                    }
                }
                else
                {
                    globTimeControl.TimeisOnRight();
                    std::cout << "\nPropagate: maxKI/respectiveKIc = " << maxKI_KIc << "\n";
                    std::cout << "************* t está a Direita de KI=KIc *************\n";
                    
                    {
                        std::ostringstream strVal, strKI;
                        strVal << globTimeControl.actDeltaT();
                        strKI << maxKI_KIc;
                        std::string tempstr = "deltaT ";
                        tempstr += strVal.str();
                        tempstr += " s : KI >= KIc (relativo = ";
                        tempstr += strKI.str();
                        tempstr += ")\n";
                        deltaT_historyMap[globTimeControl.actDeltaT()] = tempstr;
                    }
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
    
    this->CloseActualTimeStep();
    
    this->DefinePropagatedPoligonalChain(maxKI_KIc, whoPropagate);
    
    std::map<REAL,std::string>::iterator itdeltaTmap;
    for(itdeltaTmap = deltaT_historyMap.begin(); itdeltaTmap != deltaT_historyMap.end(); itdeltaTmap++)
    {
        globFractOutput3DData.fKI_KI_history << (itdeltaTmap->second);
    }
    globFractOutput3DData.fKI_KI_history << "\n----------------------------------------------------------\n";
    std::ofstream outKIhistory("000stepLOG.txt");
    outKIhistory << globFractOutput3DData.fKI_KI_history.str();
    outKIhistory.close();
    
    this->fstep++;
    
    delete an;
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::CloseActualTimeStep()
{
    globTimeControl.UpdateActTime();//atualizando esta rodada concluida para o tempo atual
    this->UpdateLeakoffKernel();
    this->PostProcessAll();
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::InitializePath3DVector()
{
    this->fPath3D.Reset();
    
    int nCrackTipElems = this->fPlaneFractureMesh->NCrackTipElements();
    
    for(int pos = 0; pos < nCrackTipElems; pos++)
    {
        TPZGeoEl * gel1D = this->fPlaneFractureMesh->GetCrackTipGeoElement(pos);
        
#ifdef PZDEBUG
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
        
        this->fPath3D.PushBackPath3D( new Path3D(center, KIc, normal, this->fJIntegralRadius, this->fmeshVec[0], this->fmeshVec[1]) );
    }
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::AssembleStiffMatrixLoadVec(TPZAnalysis * an,
                                                        TPZAutoPointer< TPZMatrix<STATE> > & matK,
                                                        TPZFMatrix<STATE> & matRes,
                                                        EWhoBlock whoBlock)
{
    this->fPlaneFractureMesh->SetActualState();
    
    TPZFStructMatrix structMat(this->fmphysics);
    an->SetStructuralMatrix(structMat);
    
    this->ApplyEquationFilter(an,whoBlock);
    
    TPZStepSolver<STATE> stepS;
    stepS.SetDirect(ELU);
    an->SetSolver(stepS);
    
    an->Assemble();
    
    matK = an->MatrixSolver<STATE>().Matrix();
    matRes = an->Rhs();
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::ApplyInitialCondition(STATE val)
{
    this->fmeshVec[0]->Solution().Zero();
    TPZFMatrix<STATE> &fmesh0sol = fmeshVec[0]->Solution();
    ///Newman=0
    fmesh0sol(0,0) = 1.;
    
    ///Newman=1 on entire fracture
    fmesh0sol(1,0) = val;
    
    ///Newman=0 on stripes
    for(int r = 2; r < fmesh0sol.Rows(); r++)
    {
        fmesh0sol(r,0) = 0.;
    }
    
    this->PutConstantPressureOnFluidSolution(val);
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::PutConstantPressureOnFluidSolution(STATE val)
{
    TPZFMatrix<STATE> &fmesh1sol = fmeshVec[1]->Solution();
    for(int r = 0; r < fmesh1sol.Rows(); r++)
    {
        fmesh1sol(r,0) = val;
    }
    TPZBuildMultiphysicsMesh::TransferFromMeshes(this->fmeshVec, this->fmphysics);
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::ApplyEquationFilter(TPZAnalysis * an, EWhoBlock whoBlock)
{
    std::set<int64_t> eqOut;
    
    int neq = this->fmphysics->NEquations();
    int64_t blockAlphaEslast = this->fmphysics->ConnectVec()[0].SequenceNumber();
    int64_t posBlock = this->fmphysics->Block().Position(blockAlphaEslast);
    eqOut.insert(posBlock);
    
    if(whoBlock == EAllBlock || whoBlock == EBlockEntireFracure)
    {
        eqOut.insert(posBlock+1);
    }
    if(whoBlock == EAllBlock || whoBlock == EBlockStripes)
    {
        for(int r = 2; r < this->fmeshVec[0]->Solution().Rows(); r++)
        {
            eqOut.insert(posBlock+r);
        }
    }
    
    TPZVec<int64_t> actEquations(neq-eqOut.size());
    int p = 0;
    for(int64_t eq = 0; eq < neq; eq++)
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

void TPZPlaneFractureKernel::MassMatrix(TPZFMatrix<STATE> & massMat)
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
    
    int64_t blockAlphaEslast = this->fmphysics->ConnectVec()[0].SequenceNumber();
    int64_t posBlock = this->fmphysics->Block().Position(blockAlphaEslast);
    
    globTimeControl.SetDeltaT(120.);
    int64_t neq = this->fmphysics->NEquations();
    int nsteps = 10;
    
    this->ApplyInitialCondition(globLayerStruct.GetHigherPreStress());

    TPZFMatrix<STATE> xIni = this->fmphysics->Solution();
    for(int64_t i = 0; i < xIni.Rows(); i++)
    {
        STATE val = (STATE)(rand())*(1.e-8);
        xIni(i,0) = val;
    }
    xIni(posBlock) = 1.;
    
    TPZAnalysis *an = new TPZAnalysis(this->fmphysics);
    an->LoadSolution(xIni);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(this->fmeshVec, this->fmphysics);
    
    TPZFMatrix<STATE> actX = xIni;
    
    TPZAutoPointer< TPZMatrix<STATE> > fL_xIni;
    TPZFMatrix<STATE> f_xIni(neq,1);
    
    AssembleStiffMatrixLoadVec(an, fL_xIni, f_xIni, ENoBlock);
    
    TPZFMatrix<STATE> fAprox_x(neq,1);
    TPZFMatrix<STATE> fExato_x(neq,1);
    
    TPZFMatrix<REAL> errorVec(neq,1,0.);
    TPZFMatrix<REAL> errorNorm(nsteps,1,0.);
    
    TPZAutoPointer< TPZMatrix<STATE> > fLtemp;
    TPZFMatrix<STATE> dFx(neq,1);
    
    TPZVec<REAL> deltaX(neq,0.001), alphas(nsteps);
    deltaX[posBlock] = 0.;
    REAL alpha;
    
    std::stringstream exatoSS, aproxSS;
    exatoSS << "exato={";
    aproxSS << "aprox={";
    for(int i = 0; i < nsteps; i++)
    {
        alpha = (i+1)/10.;
        alphas[i] = alpha;
        
        ///Fx aproximado
        dFx.Zero();
        for(int64_t r = 0; r < neq; r++)
        {
            int activeRow = r;
            if(activeRow == posBlock)
            {
                continue;
            }
            else if(activeRow > posBlock)
            {
                activeRow--;
            }
            for(int64_t c = 0; c < neq; c++)
            {
                int activeCol = c;
                if(activeCol == posBlock)
                {
                    continue;
                }
                if(activeCol >= posBlock)
                {
                    activeCol--;
                }
                dFx(r,0) +=  (-1.) * fL_xIni->GetVal(activeRow,activeCol) * (alpha * deltaX[c]); // (-1) porque fLini = -D[res,sol]
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
        for(int64_t r = 0; r < neq; r++)
        {
            if(r != posBlock)
            {
                actX(r,0) = xIni(r,0) + (alpha * deltaX[r]);
            }
        }
        an->LoadSolution(actX);
        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(this->fmeshVec, this->fmphysics);
        
        fExato_x.Zero();
        if(fLtemp) fLtemp->Zero();
        this->AssembleStiffMatrixLoadVec(an, fLtemp, fExato_x, ENoBlock);
        
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
        for(int64_t r = 0; r < neq; r++)
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
    
    delete an;
    
    DebugStop();
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::UpdateLeakoffKernel()
{
    globLeakoffStorage.UpdateLeakoff(this->fmphysics, globTimeControl.actDeltaT());
    
    
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
        
        REAL alpha = 0.;
        if(volLeakoffComputed > 0.)
        {
            alpha = volLeakoffExpected/volLeakoffComputed;
        }
        std::map<int,REAL>::iterator itLk;
        for(itLk = globLeakoffStorage.GetLeakoffMap().begin(); itLk != globLeakoffStorage.GetLeakoffMap().end(); itLk++)
        {
            itLk->second *= alpha;
        }
        std::cout << "\n-----------------------\nalpha leak off = " << alpha << "\n\n";
        
        //Verificacao
        volLeakoffComputed = ComputeVlAcumLeakoff(this->fmeshVec[1]);
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

void TPZPlaneFractureKernel::PostProcessAll()
{
    std::cout << "\n\n*********** POS-PROCESSAMENTO ***********\n";
    
    this->PostProcessSolutions();
    this->PostProcessAcumVolW();
    
    this->PostProcessLeakoff();
    this->PostProcessElasticity();
    //this->PostProcessInvariants();
    this->PostProcessPressure();
    
    this->PostProcessOutputData();
    
    std::cout << "Tempo atual desta geometria = " << globTimeControl.actTime()/60. << " minuto(s)\n\n";
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::PostProcessSolutions(int num)
{
    if(num < 0)
    {
        num = this->fstep;
    }
    std::stringstream nm;
    nm << "Solutions_Step" << num << ".txt";
    std::ofstream solutFile(nm.str().c_str());
    TPZFMatrix<STATE> &fmesh0sol = fmeshVec[0]->Solution();
    TPZFMatrix<STATE> &fmesh1sol = fmeshVec[1]->Solution();
    solutFile.precision(10);
    solutFile << "\nMeshvec[0]:\n";
    solutFile << "(Entire fracture : " << fmesh0sol(1,0) << ")\n";
    for(int r = 2; r < fmesh0sol.Rows(); r++)
    {
        solutFile << "Stripe" << (r-2) << " : " << fmesh0sol(1,0) + fmesh0sol(r,0) << ")\n";
    }
    solutFile << "\nMeshvec[1]:\n";
    for(int r = 0; r < fmesh1sol.Rows(); r++)
    {
        solutFile << fmesh1sol(r,0) << "\n";
    }
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::PostProcessAcumVolW()
{
    bool thereIsNegW;
    REAL negVol = 0.;
    REAL wAcum = this->IntegrateW(thereIsNegW, negVol);
    globFractOutput3DData.InsertTAcumVolW(globTimeControl.actTime()/60., wAcum);
    
    REAL FractArea = this->Fracture1wing_Area();
    REAL meanWm = wAcum/FractArea;//metros
    REAL meanWmm = meanWm*1000.;//milimetros
    globFractOutput3DData.InsertTmeanW(globTimeControl.actTime()/60., meanWmm);
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::PostProcessLeakoff(int num)
{
    if(num < 0)
    {
        num = this->fstep;
    }
    REAL leakAcum = this->ComputeVlAcumLeakoff(this->fmeshVec[1]);
    
    globFractOutput3DData.InsertTAcumVolLeakoff(globTimeControl.actTime()/60., leakAcum);
    
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
        nm << "LeakoffPenetration_Step" << num << ".vtk";
        std::ofstream file(nm.str().c_str());
        TPZVTKGeoMesh::PrintCMeshVTK(this->fmeshVec[1], file, penetrationValues, "PenetrationLength");
    }
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::PostProcessElasticity(int num)
{
    if(num < 0)
    {
        num = this->fstep;
    }
    TSWXGraphMesh grMesh;
    TSWXGraphElement grEl(0);
    TPZVec<std::string> nodalSol(1), cellSol(0);
    nodalSol[0] = "Displacement";
    //nodalSol[1] = "StressY";
    
    grEl.GenerateVTKData(this->fmeshVec[0], 3, 0., nodalSol, cellSol, grMesh);
    
    std::stringstream nm;
    nm << "Elasticity_Step" << num << ".vtk";
    std::ofstream file(nm.str().c_str());
    grMesh.ToParaview(file);
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::PostProcessInvariants(int num)
{
    if(num < 0)
    {
        num = this->fstep;
    }
    
    TPZMaterial * mat = this->fmeshVec[0]->FindMaterial(globMaterialIdGen.RockMatId(1));
    
    TPZManVector<std::string> scalnames(3), vecnames(0), tennames(0);
    scalnames[0] = "I1";
    scalnames[1] = "I2";
    scalnames[2] = "I3";
    
    std::stringstream nm;
    nm << "Invariants_Step" << num << ".vtk";
    
    int dim = 3;
    if(!mat){
        DebugStop();
    }
    std::set<int> matids;
    int matid = mat->Id();
    matids.insert(matid);
    TPZVTKGraphMesh vtkmesh(this->fmeshVec[0],dim,matids,scalnames,vecnames,tennames);
    vtkmesh.SetFileName(nm.str());
    vtkmesh.SetResolution(0);
    int numcases = 1;
    
    // Iteracoes de tempo
    int istep = 0;
    vtkmesh.DrawMesh(numcases);
    vtkmesh.DrawSolution(istep, 1.);
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::PostProcessPressure(int num)
{
    if(num < 0)
    {
        num = this->fstep;
    }
    TSWXGraphMesh grMesh;
    TSWXGraphElement grEl(0);
    TPZVec<std::string> nodalSol(1), cellSol(0);
    nodalSol[0] = "Pressure";
    grEl.GenerateVTKData(this->fmeshVec[1], 2, 0., nodalSol, cellSol, grMesh);
    
    std::stringstream nm;
    nm << "Pressure_Step" << num << ".vtk";
    std::ofstream file(nm.str().c_str());
    grMesh.ToParaview(file);
    
    {
        REAL integralP = 0.;
        for(int c = 0; c < this->fmeshVec[1]->NElements(); c++)
        {
            TPZCompEl * cel = this->fmeshVec[1]->ElementVec()[c];
            if(!cel || globMaterialIdGen.IsInsideFractMat(cel->Reference()->MaterialId()) == false)
            {
                continue;
            }
            TPZInterpolationSpace * intel = dynamic_cast<TPZInterpolationSpace*>(cel);
            if(!intel)
            {
                DebugStop();
            }
            TPZVec<STATE> value;
            intel->Integrate(0, value);
            
            integralP += Max(STATE(0),value[0]);
        }
        
        REAL fractA = this->Fracture1wing_Area();
        REAL pmedTot = integralP/fractA;
        REAL netpressure = pmedTot - globLayerStruct.GetLowerPreStress();
        
        globFractOutput3DData.InsertTNetPressure(globTimeControl.actTime()/60., netpressure);
    }
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::PostProcessOutputData(int num)
{
    if(num < 0)
    {
        num = this->fstep;
    }
    globFractOutput3DData.PrintFractureGeometry(num, this->fpoligonalChain, this->fCenterTVD);
    globFractOutput3DData.PrintConservationMass();
}
//------------------------------------------------------------------------------------------------------------

bool TPZPlaneFractureKernel::ThereIsNegativeNetPressure()
{
    return this->fPath3D.ThereIsNegativeNetPressure();
}
//------------------------------------------------------------------------------------------------------------

REAL TPZPlaneFractureKernel::IntegrateW(bool & thereIsNegW, REAL & negVol)
{
    thereIsNegW = false;
    negVol = 0.;
    
    this->fmeshVec[0]->LoadReferences();
    
    REAL integralW = 0.;
    for(int c = 0; c < this->fmeshVec[0]->NElements(); c++)
    {
        TPZCompEl * cel = this->fmeshVec[0]->ElementVec()[c];
        if(!cel || globMaterialIdGen.IsInsideFractMat(cel->Reference()->MaterialId()) == false)
        {
            continue;
        }
        TPZInterpolationSpace * intel = dynamic_cast<TPZInterpolationSpace*>(cel);
        if(!intel)
        {
            DebugStop();
        }
        TPZManVector<STATE> value(3,0.);
        intel->Integrate(0,value);
        
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
        
#ifdef PZDEBUG
        if(!cel || cel->Reference()->Dimension() != 2)
        {
            DebugStop();
        }
#endif
        
        TPZGeoEl * gel = cel->Reference();
        
#ifdef PZDEBUG
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
    if(globLeakoffStorage.DefaultLeakoffEnabled() == false)
    {
        return 0.;
    }
    
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
        
#ifdef PZDEBUG
        if(!cel || cel->Reference()->Dimension() != 2)
        {
            DebugStop();
        }
#endif
        
        TPZGeoEl * gel = cel->Reference();
        
#ifdef PZDEBUG
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
    REAL Qinj1wing_Hbullet = -(this->fQinj1wing_Hbullet);
    REAL Qinj1wing = Qinj1wing_Hbullet * this->fHbullet;
    
    REAL volInjected = actTime * Qinj1wing;
    
    return volInjected;
}
//------------------------------------------------------------------------------------------------------------

bool TPZPlaneFractureKernel::CheckPropagationCriteria(REAL & maxKI_KIc, std::set<int> & whoPropagate,
                                                      int & maxKI_KIcPosition)
{
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
        
        if(cracktipKI >= 0.9 * cracktipKIc)
        {
            if(cracktipKI >= 1.0 * cracktipKIc)
            {
                propagate = true;
            }
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
        it = whoPropagate.find(p);
        if(it == whoPropagate.end())
        {//Nao propagou, portanto serah mantido o ponto medio
            TPZVec<REAL> originVec = this->fPath3D.Path(p)->Origin();
            
            newX = originVec[0];
            newZ = originVec[2];
        }
        else
        {//Propagou, portanto serah definido novo ponto a partir de seu centro
            REAL KI = this->fPath3D.Path(p)->KI();
            REAL KIc = this->fPath3D.Path(p)->KIc();
            
            TPZVec<REAL> originVec = this->fPath3D.Path(p)->Origin();
            TPZVec<REAL> Jdirection = this->fPath3D.Path(p)->JDirection();
            
            //Variacao linear no decorrer do tempo de fMaxDispl para fMinDispl
            REAL alpha = 2.;
            REAL displacement = this->fMaxDispl * pow((KI/KIc)/(maxKI_KIc),alpha);
            
            newX = originVec[0] + displacement * Jdirection[0];
            newZ = originVec[2] + displacement * Jdirection[2];
            
            //Verificacoes se a fratura terh altura menor do que o canhoneado
            if(p == 0 && newZ < this->fZbulletTop)
            {
                newZ = this->fZbulletTop + 0.2;
            }
            if(p == this->fPath3D.NPaths()-1 && newZ > this->fZbulletBottom)
            {
                newZ = this->fZbulletBottom - 0.2;
            }
            
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
    this->RemoveLayerInvasion(newPoligonalChain);
    
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
            DebugStop();
        }
        if(z < this->fResBottom + tol)
        {
            std::cout << "\n\nPoligonal chain reach the bottom limit of available domain!\n";
            std::cout << "Simulation should stop!\n\n";
            DebugStop();
        }
        if(x > this->fResRight - tol)
        {
            std::cout << "\n\nPoligonal chain reach the right limit of available domain!\n";
            std::cout << "Simulation should stop!\n\n";
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
        else if(newPoligonalChain[p].second < minZ)
        {
            newPoligonalChain[p].second = minZ;
        }
    }
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::TransferLastLeakoff(TPZCompMesh * cmeshFrom)
{
    REAL leakAcumBefore = this->ComputeVlAcumLeakoff(cmeshFrom);
    
    //Mapa que guardarah os volumes acumulados dos elementos de fratura da malha antiga que nao tem pais
    std::map<int,REAL> fatherGeoIndex_VolAcum, newLeakoffMap;
    
    std::map<int,REAL>::iterator it;
    
    //Guardando os volumes acumulados dos elementos de fratura da malha antiga que nao tem pais
    int nElem = cmeshFrom->NElements();
    for(int el = 0; el < nElem; el++)
    {
        TPZCompEl * cel = cmeshFrom->ElementVec()[el];
        TPZGeoEl * gel = cel->Reference();
        
#ifdef PZDEBUG
        if(!gel)
        {
            DebugStop();
        }
#endif
        
        if(globMaterialIdGen.IsInsideFractMat(gel->MaterialId()) == false)
        {
            continue;
        }
        
#ifdef PZDEBUG
        if(gel->Dimension() != 2)
        {
            DebugStop();
        }
#endif
        
        REAL Area = gel->SideArea(gel->NSides()-1);
        it = globLeakoffStorage.GetLeakoffMap().find(gel->Id());
        
#ifdef PZDEBUG
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
        
#ifdef PZDEBUG
        if(!gel)
        {
            DebugStop();
        }
#endif
        
        if(globMaterialIdGen.IsInsideFractMat(gel->MaterialId()) == false)
        {
            continue;
        }
        
#ifdef PZDEBUG
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
    
    REAL leakAcumAfter = ComputeVlAcumLeakoff(this->fmeshVec[1]);
    if(fabs(leakAcumBefore - leakAcumAfter) > 1.E-2)
    {
        std::cout << "\n\nA transferencia de leakoff nao manteve volume!!!\n\n";
        std::cout << "leakAcumBefore = " << leakAcumBefore << "\n";
        std::cout << "leakAcumAfter = " << leakAcumAfter << "\n";
    }
}
//------------------------------------------------------------------------------------------------------------


//------------------------------------------------------------------------------------------------------------
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
