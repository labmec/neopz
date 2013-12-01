//
//  TPZPlaneFractureKernel.cpp
//  PZ
//
//  Created by Cesar Lucci on 18/11/13.
//
//

#include "TPZPlaneFractureKernel.h"

#include "pzanalysis.h"
#include "pzbndcond.h"
#include "pzfstrmatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "pzbuildmultiphysicsmesh.h"

#include "TPZVTKGeoMesh.h"

TPZPlaneFractureKernel::TPZPlaneFractureKernel()
{
    this->fPlaneFractureMesh = NULL;
    this->fmeshVec.Resize(2);
    this->fmphysics = NULL;
    fpOrder = 2;
}
//------------------------------------------------------------------------------------------------------------

TPZPlaneFractureKernel::TPZPlaneFractureKernel(TPZVec<TPZLayerProperties> & layerVec, REAL bulletTVDIni, REAL bulletTVDFin,
                                               REAL xLength, REAL yLength, REAL Lmax, int nstripes, int pOrder)
{
    this->fPlaneFractureMesh = new TPZPlaneFractureMesh(layerVec, bulletTVDIni, bulletTVDFin, xLength, yLength, Lmax, nstripes);
    this->fmeshVec.Resize(2);
    this->fmphysics = NULL;
    
    fpOrder = pOrder;
}
//------------------------------------------------------------------------------------------------------------

TPZPlaneFractureKernel::~TPZPlaneFractureKernel()
{
//    delete this->fPlaneFractureMesh;
//    
//    for(int m = 0; m < this->fmeshVec.NElements(); m++)
//    {
//        delete this->fmeshVec[m];
//    }
//    this->fmeshVec.Resize(0);
//    
//    delete this->fmphysics;
}
//------------------------------------------------------------------------------------------------------------

bool TPZPlaneFractureKernel::RunThisFractureGeometry(const TPZVec<std::pair<REAL,REAL> > &poligonalChain,
                                                     std::string vtkFile,
                                                     bool printVTKfile)
{
    REAL Qinj = -1.;//AQUICAJU
    REAL visc = 0.001E-6;//AQUICAJU

    //----------------------------------------------------------------------------------------
    this->fPlaneFractureMesh->InitializeRefinedMesh(poligonalChain);
    
    //Malha computacional elastica processada por faixas (serah referencia da fmeshVec[0])
    TPZCompMesh * fractureCMeshRef = this->fPlaneFractureMesh->GetFractureCompMesh(fpOrder);
    ProcessElasticCMeshByStripes(fractureCMeshRef);
    
    //Malha computacional do tipo CMeshReferred
    this->fmeshVec[0] = this->fPlaneFractureMesh->GetFractureCompMeshReferred(fractureCMeshRef, fpOrder);
    
    //Malha computacional de pressao
    this->fmeshVec[1] = this->fPlaneFractureMesh->GetPressureCompMesh(Qinj, fpOrder);
    
    //Malha computacional de acoplamento (multifisica)
    this->fmphysics = this->fPlaneFractureMesh->GetMultiPhysicsCompMesh(this->fmeshVec, Qinj, visc, fpOrder);
    
    //----------------------------------------------------------------------------------------
    TPZAnalysis * an = new TPZAnalysis(this->fmphysics);
    
    {
        CheckConv();
    }
    
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

    bool initialElasticKickIsNeeded = true;//AQUICAJU
    if(initialElasticKickIsNeeded)
    {
        TPZFMatrix<REAL> chutenewton(this->fmeshVec[0]->Solution().Rows(), this->fmeshVec[0]->Solution().Cols(), 1.);
        this->fmeshVec[0]->LoadSolution(chutenewton);
        TPZBuildMultiphysicsMesh::TransferFromMeshes(this->fmeshVec, this->fmphysics);
    }

    bool propagate = false;
    bool fMustStop = false;//AQUICAJU
    while( fMustStop == false && propagate == false )
    {
        fres.Zero();
        StiffMatrixLoadVec(an, matK, fres);

        res_total = fres + fmat;
        REAL res = Norm(res_total);
        REAL tol = 1.e-8;
        int maxit = 15;
        int nit = 0;

        while(res > tol && nit < maxit) //itercao de Newton
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

        if(res >= tol)
        {
            std::cout << "\nAtingido o numero maximo de iteracoes, nao convergindo portanto!!!\n";
            std::cout << "||Res|| = " << res << std::endl;
            DebugStop();
        }

        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(this->fmeshVec, this->fmphysics);
//        globFractInputData.UpdateLeakoff(fmeshVec[1]);
//        globFractInputData.UpdateActTime();
//
//        PostprocessPressure();
//
//        REAL KI = ComputeKIPlaneStrain();
//        if(KI > globFractInputData.KIc())
//        {//propagou!!!
//            globFractInputData.SetMinDeltaT();
//            propagate = true;
//        }
//        else
//        {//nao propagou!!!
//            globFractInputData.SetNextDeltaT();
//            fmat.Zero();
//            MassMatrix(fmat);
//
//            globFractOutputData.PlotElasticVTK(an);
//            PostProcessAcumVolW();
//            PostProcessVolLeakoff();
//        }
//        globFractOutputData.InsertTKI(globFractInputData.actTime(), KI);//its for output data to txt (Mathematica format)
//        REAL peteleco = 1.E-8;
//        if( globFractInputData.actTime() > (globFractInputData.Ttot() - peteleco) )
//        {
//            fMustStop = true;
//        }
        
        fMustStop = true;//AQUICAJU
    }

    return propagate;
    
//    {//Solve caso fosse nao linear!!!
//        REAL chute = 10.;
//        for(int r = 0; r < anRef.Solution().Rows(); r++)
//            for(int c = 0; c < anRef.Solution().Cols(); c++)
//                anRef.Solution()(r,c) = chute;
//        anRef.Run();
//        for(int r = 0; r < anRef.Solution().Rows(); r++)
//            for(int c = 0; c < anRef.Solution().Cols(); c++)
//                anRef.Solution()(r,c) += chute;
//    }
    
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::ProcessElasticCMeshByStripes(TPZCompMesh * cmesh)
{
    TPZAnalysis * an = new TPZAnalysis;
    
    int NStripes = this->fPlaneFractureMesh->NStripes();
    TPZFMatrix<STATE> solutions(cmesh->Solution().Rows(), NStripes);
    
    for(int stripe = 0; stripe < NStripes; stripe++)
    {
        // Resolvendo um problema modelo de elastica linear para utilizar a
        // solucao como espaco de aproximacao do problema nao linear acoplado
        this->fPlaneFractureMesh->SetSigmaNStripeNum(cmesh, stripe);
        
        bool mustOptimizeBandwidth = (stripe == 0);
        an->SetCompMesh(cmesh, mustOptimizeBandwidth);
        this->SolveInitialElasticity(*an, cmesh);
        
        for(int r = 0; r < cmesh->Solution().Rows(); r++)
        {
            solutions(r,stripe) = cmesh->Solution()(r,0);
        }
        
        {
            std::stringstream nm;
            nm << "SolStripeLinear.vtk";

            TPZManVector<std::string,10> scalnames(0), vecnames(1);
            
            vecnames[0] = "Displacement";
            
            const int dim = 3;
            int div = 0;
            an->SetStep(stripe);
            an->DefineGraphMesh(dim,scalnames,vecnames,nm.str());
            an->PostProcess(div);
        }
    }
    cmesh->LoadSolution(solutions);
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFractureKernel::SolveInitialElasticity(TPZAnalysis &an, TPZCompMesh * cmesh)
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
        
        int wantToSeeRow = 0;
        {
            REAL aproxSol = fAprox_x(wantToSeeRow,0);
            aproxSS << aproxSol;
            if(i < nsteps-1)
            {
                aproxSS << ",";
            }
        }
        
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
        
        {
            REAL exatoSol = fExato_x(wantToSeeRow,0);
            exatoSS << exatoSol;
            if(i < nsteps-1)
            {
                exatoSS << ",";
            }
        }
        
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

