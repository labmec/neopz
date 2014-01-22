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
#include "TPZSpStructMatrix.h"
#include "pzbuildmultiphysicsmesh.h"

#include "pzintel.h"
#include "tpzintpoints.h"
#include "TPZVTKGeoMesh.h"

#include "TPZVTKGeoMesh.h"


//#define usingSWXGraphs

#ifdef usingSWXGraphs
#include "TSWXGraphMesh.h"
#include "TSWXGraphElement.h"
#endif

//Utilize 1. para output (Mathematica) em metros e 3.280829131 para output (Mathematica) em foot
const REAL feet = 1.;//3.280829131;


TPZPlaneFractureKernel::TPZPlaneFractureKernel()
{
    DebugStop();//Use constructor below;
}
//------------------------------------------------------------------------------------------------------------

TPZPlaneFractureKernel::TPZPlaneFractureKernel(TPZVec<TPZLayerProperties> & layerVec, REAL bulletTVDIni, REAL bulletTVDFin,
                                               REAL xLength, REAL yLength, REAL Lmax, int nstripes, REAL Qinj_well, REAL visc,
                                               REAL Jradius,
                                               int pOrder,
                                               REAL MaxDispl,
                                               REAL MinDispl)
{
    this->fLmax = Lmax;
    
    this->fPlaneFractureMesh = new TPZPlaneFractureMesh(layerVec, bulletTVDIni, bulletTVDFin, xLength, yLength, Lmax, nstripes);
    this->fmeshVec.Resize(2);
    this->fmphysics = NULL;
    
    this->fHbullet = (bulletTVDFin - bulletTVDIni);
    
    REAL Qinj1wing = Qinj_well/2.;
    
    {
        globFractOutput3DData.SetQinj1wing(-Qinj1wing);
    }
    
    REAL Qinj1wing_Hbullet = Qinj1wing/this->fHbullet;
    this->fQinj1wing_Hbullet = Qinj1wing_Hbullet;
    
    this->fCenterTVD = (bulletTVDIni + bulletTVDFin)/2.;
    this->fPoligonalChainInitialRadius = 1.1 * (bulletTVDFin - fCenterTVD);
    
    this->fvisc = visc;
    
    this->fJIntegralRadius = Jradius;
    
    this->fpOrder = pOrder;
    
    this->fMaxDispl = MaxDispl;
    this->fMinDispl = MinDispl;
    
    this->fPath3D.Reset();
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
        std::ofstream outPoligMath("PoligonalChainMath0.txt");
        outPoligMath << "AllPolChains={};\n";
        outPoligMath << "Lgr={};\n";
        outPoligMath << "Hsupgr={};\n";
        outPoligMath << "Hinfgr={};";
    }
    
    TPZVec< std::pair<REAL,REAL> > poligonalChain(0);
    FillInitialFractureDots(poligonalChain);
    PostProcessFractGeometry(0, poligonalChain);
    
    REAL volAcum = 0.;
    TPZCompMesh * lastPressureCMesh = NULL;
    
    bool firstGeometry = true;
    bool reachEndOfTime = false;
    
    int step = 1;
    while(reachEndOfTime == false)
    {
        //GeoMesh
        this->fPlaneFractureMesh->InitializeFractureGeoMesh(poligonalChain);
        
        std::cout << "\n\n************** GERANDO MALHAS COMPUTACIONAIS\n";
        
        //Malha computacional elastica processada por faixas (serah referencia da fmeshVec[0])
        TPZCompMesh * fractureCMeshRef = this->fPlaneFractureMesh->GetFractureCompMesh(this->fpOrder);
        this->ProcessElasticCMeshByStripes(fractureCMeshRef);
        
        //Malha computacional do tipo CMeshReferred
        this->fmeshVec[0] = this->fPlaneFractureMesh->GetFractureCompMeshReferred(fractureCMeshRef, this->fpOrder);
        
        //Malha computacional de pressao
        this->fmeshVec[1] = this->fPlaneFractureMesh->GetPressureCompMesh(this->fQinj1wing_Hbullet, this->fpOrder);
        
        //Malha computacional de acoplamento (multifisica)
        this->fmphysics = this->fPlaneFractureMesh->GetMultiPhysicsCompMesh(this->fmeshVec, this->fQinj1wing_Hbullet, this->fvisc, this->fpOrder);
        
        this->InitializePath3DVector();

        if(firstGeometry == false)
        {
            this->TransferElasticSolution(volAcum);
            this->TransferLeakoff(lastPressureCMesh);
        }
        
        //Resolvendo o problema acoplado marchando no tempo
        this->RunThisFractureGeometry(poligonalChain, firstGeometry, step);
        
        volAcum = this->IntegrateW(fmeshVec[0]);
        lastPressureCMesh = fmeshVec[1];
        
        firstGeometry = false;
        
        REAL timeTol = 1.E-1;//Just avoiding float comparisson
        reachEndOfTime = globTimeControl.actTime() >= globTimeControl.Ttot() - timeTol;
    }//end of while(reachEndOfTime == false)
    
    std::ofstream outMath("ProstProcess.txt");
    globFractOutput3DData.PrintMathematica(outMath);
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::FillInitialFractureDots(TPZVec< std::pair<REAL,REAL> > &fractureDots)
{
    //*********** TESTE ENGLAND-GREEN AQUICAJU
//    int nptsUp = 2;
//    fractureDots.Resize(nptsUp);
//    fractureDots[0] = std::make_pair(0.5,-2110.);
//    fractureDots[1] = std::make_pair(69.5,-2110.);

    /** Elipse */
    int nptsUp = 50;
    fractureDots.Resize(0);
    
    REAL yc = -fCenterTVD;
    REAL sAx = 11.;
    REAL sAy = 10.;
    
    for(int p = 1; p <= nptsUp; p++)
    {
        REAL x  = MIN(sAx - 0.001,p * sAx/nptsUp);
        REAL fx = yc + (sAy*sqrt(sAx*sAx - x*x))/sAx;
        if(x > this->fLmax)
        {
            int oldSize = fractureDots.NElements();
            fractureDots.Resize(oldSize+1);
            fractureDots[oldSize] = std::make_pair(x,fx);
        }
    }
    for(int p = nptsUp-1; p > 0; p--)
    {
        REAL x  = MIN(sAx - 0.001,p * sAx/nptsUp);
        REAL fx = yc - (sAy*sqrt(sAx*sAx - x*x))/sAx;
        if(x > this->fLmax)
        {
            int oldSize = fractureDots.NElements();
            fractureDots.Resize(oldSize+1);
            fractureDots[oldSize] = std::make_pair(x,fx);
        }
    }
    
    
    /** CIRCULO */
//    REAL Lmax = 0.5;
//    int nsteps = M_PI * fPoligonalChainInitialRadius / Lmax;
//    if(nsteps < 10) nsteps = 10;
//    REAL ang = M_PI / nsteps;
//    
//    int nnodes = nsteps + 1;
//    fractureDots.Resize(nnodes-2);
//    
//    for(int node = 1; node < nnodes-1; node++)
//    {
//        REAL vx = 0.1 + fPoligonalChainInitialRadius*sin(node*ang);
//        REAL vz = fPoligonalChainInitialRadius*cos(node*ang) - fCenterTVD;
//        fractureDots[node-1] = std::make_pair(vx,vz);
//    }
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::RunThisFractureGeometry(TPZVec<std::pair<REAL,REAL> > &poligonalChain,
                                                     bool initialElasticKickIsNeeded,
                                                     int &step)
{
    TPZAnalysis * an = new TPZAnalysis(this->fmphysics);
    
    /** Convergence test */
    //CheckConv();
    /**********************/
    
    int nrows = an->Solution().Rows();
    TPZFMatrix<REAL> res_total(nrows,1,0.);

    TPZFMatrix<REAL> SolIterK = this->fmphysics->Solution();
    TPZAutoPointer< TPZMatrix<REAL> > matK;
    TPZFMatrix<REAL> fres(this->fmphysics->NEquations(),1);
    TPZFMatrix<REAL> fmat(this->fmphysics->NEquations(),1);
    fres.Zero();
    fmat.Zero();

    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(this->fmeshVec, this->fmphysics);

    MassMatrix(fmat);

    if(initialElasticKickIsNeeded)
    {
        TPZFMatrix<REAL> chutenewton(this->fmeshVec[0]->Solution().Rows(), this->fmeshVec[0]->Solution().Cols(), 1.);
        this->fmeshVec[0]->LoadSolution(chutenewton);
        TPZBuildMultiphysicsMesh::TransferFromMeshes(this->fmeshVec, this->fmphysics);
    }

    bool reachEndOfTime = false;
    bool propagate = false;
    
    while(reachEndOfTime == false && propagate == false)
    {
        std::cout << "\n\n\n************** CALCULANDO SOLUCAO ACOPLADA (KERNEL)\n";
        
        fres.Zero();
        StiffMatrixLoadVec(an, matK, fres);

        res_total = fres + fmat;

        REAL res = Norm(res_total);
        REAL tol = 1.e-4;
        int maxit = 15;
        int nit = 0;

        // Metodo de Newton
        while(res > tol && nit < maxit)
        {
            an->Rhs() = res_total;
            an->Solve();
            an->LoadSolution(SolIterK + an->Solution());

            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(this->fmeshVec, this->fmphysics);

            SolIterK = an->Solution();

            fres.Zero();
            StiffMatrixLoadVec(an, matK, fres);
            res_total = fres + fmat;

            res = Norm(res_total);
            std::cout << "||res|| = " << res << std::endl;
            nit++;
        }

        if(res > 10.*tol)
        {
            std::cout << "\n\nAtingido o numero maximo de iteracoes, nao convergindo portanto!!!\n";
            std::cout << "||Res|| = " << res << std::endl;
            DebugStop();
        }

        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(this->fmeshVec, this->fmphysics);
        
        globTimeControl.UpdateActTime();//atualizando esta rodada concluida para o tempo atual
        REAL timeTol = 1.E-1;//Just avoiding float comparisson
        reachEndOfTime = globTimeControl.actTime() >= globTimeControl.Ttot() - timeTol;
        
        globLeakoffStorage.UpdateLeakoff(fmphysics,globTimeControl.actDeltaT());
        {
            //O UpdateLeakoff atualiza o leakoff com base na pressao no centro do elemento (=volLeakoffComputed).
            //Jah o kernel considera o leakoff com as pressoes calculadas em cada ponto da regra de integracao (=volLeakoffExpected).
            //Por esta razao, alguma diferenca pode ocorrer, necessitando portanto corrigir o mapa de leakoff antes de prosseguir.
            REAL volInjected = ComputeVolInjected();
            REAL volW = IntegrateW(fmeshVec[0]);
            
            REAL volLeakoffExpected = volInjected - volW;
            REAL volLeakoffComputed = ComputeVlAcumLeakoff(fmeshVec[1]);
            
            if(volW < 0.)
            {
                std::cout << "\n\n\nO leakoff deve ter subtraido mais volume do que entrou\nocasionando a inversao da abertura da fratura!\n";
                std::cout << "volInjected = " << volInjected << "\n";
                std::cout << "volLeakoffComputed = " << volLeakoffComputed << "\n\n\n";
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
        
        propagate = CheckPropagationCriteria(poligonalChain);
        if(propagate)
        {//reducao do deltaT da proxima rodada para o valor minimo.
            globTimeControl.SetMinDeltaT();
            std::cout << "\n\n************* PROPAGOU! *************\n\n\n";
        }
        else
        {//avancar no valor de deltaT para a proxima rodada.
            globTimeControl.SetNextDeltaT();
            fmat.Zero();
            MassMatrix(fmat);
            std::cout << "\n\n************* NAO PROPAGOU! *************\n\n";
        }
        
        PostProcessAll(step, poligonalChain);
    }
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
        
        this->fPath3D.PushBackPath3D( new Path3D(this->fmeshVec[0], this->fmeshVec[1], center, normal, this->fJIntegralRadius) );
    }
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::ProcessElasticCMeshByStripes(TPZCompMesh * cmesh)
{
    std::cout << "\n\n************** CALCULANDO SOLUCAO ELASTICA DE REFERENCIA POR FAIXAS\n";
    
    TPZAnalysis * an = new TPZAnalysis;
    
    int NStripes = this->fPlaneFractureMesh->NStripes();
    TPZFMatrix<STATE> solutions(cmesh->Solution().Rows(), NStripes);
    
    for(int stripe = 0; stripe < NStripes; stripe++)
    {
        std::cout << "\nFaixa " << stripe+1 << " de " << NStripes << " *********************************\n";
        
        // Resolvendo um problema modelo de elastica linear para utilizar a
        // solucao como espaco de aproximacao do problema nao linear acoplado
        this->fPlaneFractureMesh->SetSigmaNStripeNum(cmesh, stripe); 
        
        bool mustOptimizeBandwidth = (stripe == 0);
        an->SetCompMesh(cmesh, mustOptimizeBandwidth);
        this->SolveElasticity(*an, cmesh);
        
        for(int r = 0; r < cmesh->Solution().Rows(); r++)
        {
            solutions(r,stripe) = cmesh->Solution()(r,0);
        }
        
//        {
//            std::stringstream nm;
//            nm << "PreElastic" << stripe << ".vtk";
//
//            TPZManVector<std::string,10> scalnames(0), vecnames(1);
//            vecnames[0] = "Displacement";
//            
//            const int dim = 3;
//            int div = 0;
//            an->SetStep(stripe);
//            an->DefineGraphMesh(dim,scalnames,vecnames,nm.str());
//            an->PostProcess(div);
//        }
    }

    cmesh->LoadSolution(solutions);
    
    //*********** TESTE ENGLAND-GREEN AQUICAJU
//    {
//        this->fPath3D.Reset();
//        
//        int nCrackTipElems = this->fPlaneFractureMesh->NCrackTipElements();
//        
//        for(int pos = 0; pos < nCrackTipElems; pos++)
//        {
//            TPZGeoEl * gel1D = this->fPlaneFractureMesh->GetCrackTipGeoElement(pos);
//            
//#ifdef DEBUG
//            if(!gel1D || gel1D->Dimension() != 1 || gel1D->MaterialId() != globMaterialIdGen.CrackTipMatId())
//            {
//                DebugStop();
//            }
//#endif
//            
//            TPZVec<REAL> qsi(1,0.), center(3,0.), normal(3,0.);
//            gel1D->X(qsi,center);
//            
//            REAL n0x = gel1D->NodePtr(0)->Coord(0);
//            REAL n0z = gel1D->NodePtr(0)->Coord(2);
//            REAL n1x = gel1D->NodePtr(1)->Coord(0);
//            REAL n1z = gel1D->NodePtr(1)->Coord(2);
//            
//            /**
//             * Da forma como o codigo foi construido, os elementos 1D do crackTip seguem em sentido anti-horario em relacao ao eixo Y.
//             * O plano na integral-J eh definido por uma normal, a qual servirah como eixo da regra da mao direita.
//             * Desta forma, a normal da integral-J corresponde aa direcao oposta do crackTip, para que o arco saia de y=0 (fora da fratura),
//             * passe por y>0 (no interior do meio elastico 3D) e acabe em y=0 (dentro da fratura).
//             */
//            normal[0] = (n0x - n1x);//x
//            normal[2] = (n0z - n1z);//z
//
//            this->fPath3D.PushBackPath3D( new Path3D(cmesh, this->fmeshVec[1], center, normal, this->fJIntegralRadius) );
//        }
//    }
    //*******************************
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::SolveElasticity(TPZAnalysis &an, TPZCompMesh * cmesh)
{
	TPZSkylineStructMatrix full(cmesh); //caso simetrico
	an.SetStructuralMatrix(full);
	TPZStepSolver<REAL> step;
	step.SetDirect(ELDLt); //caso simetrico
	an.SetSolver(step);
	an.Run();
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::StiffMatrixLoadVec(TPZAnalysis *an, TPZAutoPointer< TPZMatrix<REAL> > & matK1, TPZFMatrix<REAL> &fvec)
{
	this->fPlaneFractureMesh->SetActualState();
    
    TPZFStructMatrix matsk(this->fmphysics);
    
	an->SetStructuralMatrix(matsk);
	TPZStepSolver<REAL> step;
    
	step.SetDirect(ELU);
	an->SetSolver(step);
    
    an->Assemble();
	
    matK1 = an->Solver().Matrix();
    
	fvec = an->Rhs();
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::MassMatrix(TPZFMatrix<REAL> & Un)
{
    this->fPlaneFractureMesh->SetPastState();
    
	TPZSpStructMatrix matsp(this->fmphysics);
	TPZAutoPointer<TPZGuiInterface> guiInterface;
    matsp.CreateAssemble(Un,guiInterface);
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::CheckConv()
{
    long neq = fmphysics->NEquations();
    int nsteps = 10;
    
    TPZFMatrix<REAL> xIni(neq,1);
    for(long i = 0; i < xIni.Rows(); i++)
    {
        REAL val = (double)(rand())*(1.e-10);
        xIni(i,0) = val;
    }
    TPZAnalysis *an = new TPZAnalysis(fmphysics);
    an->LoadSolution(xIni);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshVec, fmphysics);
    
    TPZFMatrix<REAL> actX = xIni;
    
    TPZAutoPointer< TPZMatrix<REAL> > fL_xIni;
    TPZFMatrix<REAL> f_xIni(neq,1);
    
    StiffMatrixLoadVec(an, fL_xIni, f_xIni);
    if(fL_xIni->Rows() != neq || fL_xIni->Cols() != neq || fL_xIni->IsDecomposed())
    {
        DebugStop();
    }
    
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
        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshVec, fmphysics);
        
        fExato_x.Zero();
        if(fLtemp) fLtemp->Zero();
        StiffMatrixLoadVec(an, fLtemp, fExato_x);
        
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
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::PostProcessAll(int & step, TPZVec<std::pair<REAL,REAL> > & PoligonalChain)
{
    std::cout << "\n\n\n*********** POS-PROCESSAMENTO ***********\n";
    
    PostProcessAcumVolW();
    PostProcessVolLeakoff(step);
    PostProcessElasticity(step);
    PostProcessPressure(step);
    PostProcessFractGeometry(step, PoligonalChain);
    std::cout << "Final do processamento " << globTimeControl.actTime()/60. << " minuto(s)\n";
    step++;
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::PostProcessAcumVolW()
{
    REAL wAcum = IntegrateW(fmeshVec[0]);
    globFractOutput3DData.InsertTAcumVolW(globTimeControl.actTime(), wAcum);
    
    std::cout.precision(10);
    std::cout << "wAcum = " << wAcum << ";\n";
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::PostProcessVolLeakoff(int step)
{
    REAL vlAcum = ComputeVlAcumLeakoff(fmeshVec[1]);
    
    globFractOutput3DData.InsertTAcumVolLeakoff(globTimeControl.actTime(), vlAcum);
    
    {
        std::map<int,REAL>::iterator it;
        TPZVec<REAL> penetrationValues(fmeshVec[1]->NElements(), 0.);
        for(int el = 0; el < fmeshVec[1]->NElements(); el++)
        {
            TPZCompEl * cel = fmeshVec[1]->ElementVec()[el];
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
        nm << "LeakoffPenetration" << step << ".vtk";
        std::ofstream file(nm.str().c_str());
        TPZVTKGeoMesh::PrintCMeshVTK(fmeshVec[1], file, penetrationValues);
    }
    
    std::cout.precision(10);
    std::cout << "vlAcum = " << vlAcum << ";\n";
    std::cout << "VlInjected = " << ComputeVolInjected() << ";\n";
    std::cout << "wAcum+vlAcum\n";
    std::cout << "Eficiencia=wAcum/VlInjected\n";
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::PostProcessElasticity(int step)
{
#ifdef usingSWXGraphs
    TSWXGraphMesh grMesh;
    TSWXGraphElement grEl(0);
    TPZVec<std::string> nodalSol(1), cellSol(0);
    nodalSol[0] = "Displacement";
    //nodalSol[1] = "StressY";

    grEl.GenerateVTKData(fmeshVec[0], 3, 0., nodalSol, cellSol, grMesh);
    
    std::stringstream nm;
    nm << "Elasticity" << step << ".vtk";
    std::ofstream file(nm.str().c_str());
    grMesh.ToParaview(file);
#endif
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::PostProcessPressure(int step)
{
#ifdef usingSWXGraphs
    TSWXGraphMesh grMesh;
    TSWXGraphElement grEl(0);
    TPZVec<std::string> nodalSol(1), cellSol(0);
    nodalSol[0] = "Pressure";
    grEl.GenerateVTKData(fmeshVec[1], 2, 0., nodalSol, cellSol, grMesh);
    
    std::stringstream nm;
    nm << "Pressure" << step << ".vtk";
    std::ofstream file(nm.str().c_str());
    grMesh.ToParaview(file);
#endif
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::PostProcessFractGeometry(int step, TPZVec<std::pair<REAL,REAL> > & PoligonalChain)
{
    REAL Lfrac = 0;
    REAL Hsup = 0.;
    REAL Hinf = 0.;
    
    for(int p = 0; p < PoligonalChain.NElements(); p++)
    {
        REAL x = PoligonalChain[p].first;
        REAL z = PoligonalChain[p].second + fCenterTVD;
        
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
        nmMath << "PoligonalChainMath" << step+1 << ".txt";
        std::ofstream outPoligMath(nmMath.str().c_str());
        
        for(int p = 0; p < PoligonalChain.NElements(); p++)
        {
            outPoligMath << "fractureDots" << p << " = {" << feet * PoligonalChain[p].first << "," << feet * PoligonalChain[p].second + fCenterTVD << "};\n";
            nmAux << "fractureDots" << p;
            if(p < PoligonalChain.NElements()-1)
            {
                nmAux << ",";
            }
        }
        nmAux << "};\n";
        nmAux << "gr" << step+1 << " = ListPlot[pcm, Joined -> True, AxesOrigin -> {0,0}, AspectRatio -> 2];\n";
        nmAux << "L" << step+1 << " = Max[Transpose[pcm][[1]]];\n";
        nmAux << "Hsup" << step+1 << " = Max[Transpose[pcm][[2]]];\n";
        nmAux << "Hinf" << step+1 << " = -Min[Transpose[pcm][[2]]];\n";
        nmAux << "AppendTo[AllPolChains,gr" << step+1 << "];\n";
        nmAux << "AppendTo[Lgr,{" << globTimeControl.actTime()/60. << ",L" << step+1 << "}];\n";
        nmAux << "AppendTo[Hsupgr,{" << globTimeControl.actTime()/60. << ",Hsup" << step+1 << "}];\n";
        nmAux << "AppendTo[Hinfgr,{" << globTimeControl.actTime()/60. << ",Hinf" << step+1 << "}];\n";
        nmAux << "Print[\"L" << step+1 << " = \", L" << step+1 << "]\n";
        nmAux << "Print[\"Hsup" << step+1 << " = \", Hsup" << step+1 << "]\n";
        nmAux << "Print[\"Hinf" << step+1 << " = \", Hinf" << step+1 << "]\n";
        outPoligMath << nmAux.str();
    }
}
//------------------------------------------------------------------------------------------------------------

REAL TPZPlaneFractureKernel::IntegrateW(TPZCompMesh * elasticCMesh)
{
    REAL integral = 0.;
    for(int c = 0; c < elasticCMesh->NElements(); c++)
    {
        TPZCompEl * cel = elasticCMesh->ElementVec()[c];
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
        
        integral += value[1];//Integrando w (uy)
    }
    
    return 2.*integral;//Aqui jah eh considerada a simetria
}
//------------------------------------------------------------------------------------------------------------

REAL TPZPlaneFractureKernel::ComputeVlAcumLeakoff(TPZCompMesh * fluidCMesh)
{
    REAL vlAcum = 0.;
    
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
        vlAcum += 2. * (Area * penetration);//Aqui jah eh considerada a simetria
    }
    
    return vlAcum;
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

bool TPZPlaneFractureKernel::CheckPropagationCriteria(TPZVec< std::pair<REAL,REAL> > &newPoligonalChain)
{
    bool propagate = false;
    
    std::map< int, std::pair<REAL,REAL> > whoPropagate_KI;
    
    //Calculo das integrais-J
    fPath3D.IntegratePath3D();
    
    REAL maxKI = 0., respectiveKIc = 0.;
    //Calculo de KI
    for(int p = 1; p < fPath3D.NPaths() - 1; p++)
    {
        REAL Jintegral = fPath3D.Path(p)->Jintegral();
        /**
         * Obs.:    Precisa verificar a validade da situacao em que, quando o vetor direcao da integral-J estiver para
         *          dentro da fratura, a nova frente da fratura terah retrocedido, indicando um estado compressivo,
         *          ocasionando o fechamento da fratura neste ponto.
         */
    
        REAL originZcoord = fPath3D.Path(p)->OriginZcoord();
        REAL layerKIc = fPlaneFractureMesh->GetKIcFromLayerOfThisZcoord(originZcoord);

        REAL young = 0., poisson = 0.;
        fPlaneFractureMesh->GetYoung_and_PoissonFromLayerOfThisZcoord(originZcoord, young, poisson);
        REAL cracktipKI = sqrt( Jintegral*young );
        
        {
            TPZVec<REAL> originVec = fPath3D.Path(p)->Origin();
            TPZVec<REAL> Jdirection = fPath3D.Path(p)->JDirection();
            TPZVec<REAL> ptTemp(3,0.), qsiTemp(2,0.);
            ptTemp[0] = originVec[0] + 0.1 * Jdirection[0];
            ptTemp[2] = originVec[2] + 0.1 * Jdirection[2];

            TPZGeoEl * gel = fPlaneFractureMesh->Find2DElementNearCrackTip(p, ptTemp);
            if( !gel || globMaterialIdGen.IsInsideFractMat(gel->MaterialId()) )
            {
                //Nao sera propagado para pontos fora do dominio
                //ou no interior da fratura (i.e., fechamento).
                cracktipKI = 0.;
            }
        }
        if(cracktipKI >= layerKIc)
        {
            propagate = true;
            whoPropagate_KI[p] = std::make_pair(cracktipKI,layerKIc);
            
            if(cracktipKI > maxKI)
            {
                maxKI = cracktipKI;
                respectiveKIc = layerKIc;
            }
        }
    }
    
    //Definicao da newPoligonalChain (propagada)
    if(propagate)
    {
        DefinePropagatedPoligonalChain(maxKI, respectiveKIc, whoPropagate_KI, newPoligonalChain);
    }
    
    return propagate;
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::DefinePropagatedPoligonalChain(REAL maxKI, REAL respectiveKIc,
                                                            std::map< int, std::pair<REAL,REAL> > &whoPropagate_KI,
                                                            TPZVec< std::pair<REAL,REAL> > &poligonalChain)
{
    TPZVec< std::pair<REAL,REAL> > newPoligonalChain(0);
    
    std::map< int, std::pair<REAL,REAL> >::iterator it;
    
    REAL newX = 0.;
    REAL newZ = 0.;
    for(int p = 1; p < poligonalChain.NElements() - 2; p++)
    {
        it = whoPropagate_KI.find(p);
        if(it == whoPropagate_KI.end())
        {//Nao propagou, portanto serah mantido o ponto medio
            TPZVec<REAL> originVec = fPath3D.Path(p)->Origin();
            
            newX = originVec[0];
            newZ = originVec[2];
        }
        else
        {//Propagou, portanto serah definido novo ponto a partir de seu centro
            REAL KI = it->second.first;
            REAL KIc = it->second.second;
            
            TPZVec<REAL> originVec = fPath3D.Path(p)->Origin();
            TPZVec<REAL> Jdirection = fPath3D.Path(p)->JDirection();
            
            //Variacao linear no decorrer do tempo de fMaxDispl para fMinDispl
            REAL dLmax = this->fMaxDispl + (globTimeControl.actTime()*(this->fMinDispl - this->fMaxDispl))/globTimeControl.Ttot();
            REAL alpha = 1.;//alphaMin = [1.0~2.0]
            
            REAL displacement = dLmax * pow((KI - KIc)/(maxKI - KIc),alpha);
            
            newX = originVec[0] + displacement * Jdirection[0];
            newZ = originVec[2] + displacement * Jdirection[2];
        }
        if(newX > fLmax)
        {
            int oldSize = newPoligonalChain.NElements();
            newPoligonalChain.Resize(oldSize+1);
            newPoligonalChain[oldSize] = std::make_pair(newX,newZ);
        }
    }
    
    RemoveZigZag(newPoligonalChain);
    RemoveZigZag(newPoligonalChain);
    
    bool applyBezier = true;
    if(applyBezier)
    {
        BezierCurve bz(newPoligonalChain);
        
        int npts = 100;
        poligonalChain.Resize(npts);
        for(int p = 0; p < npts; p++)
        {
            std::pair<REAL,REAL> pt;
            REAL t = p/(double(npts-1));
            bz.F(t, pt);
            poligonalChain[p] = pt;
        }
    }
    else
    {
        poligonalChain = newPoligonalChain;
    }
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::RemoveZigZag(TPZVec< std::pair< REAL,REAL > > &poligonalChain)
{
    TPZVec< std::pair< REAL,REAL > > NOzigzagPoligonalChain(0);
    
    REAL v0x = 1.;
    REAL v0z = 0.;
    for(int p = 0; p < poligonalChain.NElements()-1; p++)
    {
        std::pair<REAL,REAL> p0 = poligonalChain[p];
        std::pair<REAL,REAL> p1 = poligonalChain[p+1];
        
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
        
        v0x = v1x;
        v0z = v1z;
    }
    int oldSize = NOzigzagPoligonalChain.NElements();
    NOzigzagPoligonalChain.Resize(oldSize+1);
    NOzigzagPoligonalChain[oldSize] = poligonalChain[poligonalChain.NElements()-1];
    
    poligonalChain = NOzigzagPoligonalChain;
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::TransferElasticSolution(REAL volAcum)
{
    std::cout << "\n\n\n************** TRANSFERINDO VOLUME PARA NOVA MALHA ELASTICA (PROPAGADA)\n";
    
    TPZFMatrix<REAL> wSol(this->fmeshVec[0]->Solution().Rows(), this->fmeshVec[0]->Solution().Cols(), 1.);
    this->fmeshVec[0]->LoadSolution(wSol);
    
    REAL newVolAcum = this->IntegrateW(fmeshVec[0]);
    REAL alpha = volAcum/newVolAcum;
    
    //Sintonia fina (corrigindo os volumes de antes e depois)
    if(fabs(volAcum - newVolAcum) > 1.E-10)
    {
        TPZFMatrix<REAL> newSolution = fmeshVec[0]->Solution();
        for(int r = 0; r < newSolution.Rows(); r++)
        {
            for(int c = 0; c < newSolution.Cols(); c++)
            {
                newSolution(r,c) *= alpha;
            }
        }
        fmeshVec[0]->LoadSolution(newSolution);
    }
    
    //Verificando se a correcao deu certo
    newVolAcum = this->IntegrateW(fmeshVec[0]);
    if(fabs(volAcum - newVolAcum) > 1.E-10)
    {
        std::cout << "\n\n\nW nao manteve volume na transferencia de solucao elastica!!!\n\n\n";
        std::cout << "volAntes = " << volAcum << std::endl;
        std::cout << "volDepois = " << newVolAcum << std::endl;
        std::cout << "Dif = " << (newVolAcum - volAcum) << std::endl;
        DebugStop();
    }
    else
    {
        TPZBuildMultiphysicsMesh::TransferFromMeshes(this->fmeshVec, this->fmphysics);
    }
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::TransferLeakoff(TPZCompMesh * cmeshFrom)
{
#ifdef DEBUG
    REAL vlAcumBefore = ComputeVlAcumLeakoff(cmeshFrom);
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
    nElem = fmeshVec[1]->NElements();
    for(int el = 0; el < nElem; el++)
    {
        TPZCompEl * cel = fmeshVec[1]->ElementVec()[el];
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
    REAL vlAcumAfter = ComputeVlAcumLeakoff(fmeshVec[1]);
    if(fabs(vlAcumBefore - vlAcumAfter) > 1.E-5)
    {
        std::cout << "\n\nA transferencia de leakoff nao manteve volume!!!\n\n";
        std::cout << "vlAcumBefore = " << vlAcumBefore << "\n";
        std::cout << "vlAcumAfter = " << vlAcumAfter << "\n";
        DebugStop();
    }
#endif
}


//---------------------------------------------------------


BezierCurve::BezierCurve()
{
    DebugStop();//Utilize o outro construtor
}

BezierCurve::BezierCurve(TPZVec< std::pair< REAL,REAL > > &poligonalChain)
{
    forder = poligonalChain.NElements()-1;
    fPoligonalChain = poligonalChain;
}

BezierCurve::~BezierCurve()
{
    forder = 0;
    fPoligonalChain.Resize(0);
}

void BezierCurve::F(REAL t, std::pair< REAL,REAL > & ft)
{
    if(t < 1.E-5)
    {
        ft = fPoligonalChain[0];
        return;
    }
    else if(t > 0.99999)
    {
        ft = fPoligonalChain[forder];
        return;
    }
    else
    {
        REAL x = 0., z = 0.;
        for(int i = 0; i <= forder; i++)
        {
            REAL b = Bernstein(t,i);
            
            x += b * fPoligonalChain[i].first;
            z += b * fPoligonalChain[i].second;
        }
        ft = std::make_pair(x,z);
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