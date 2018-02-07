#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "Common.h"

#include "TPZSBFemElementGroup.h"
#include "pzbndcond.h"

#include "pzskylstrmatrix.h"
#include "TPZSSpStructMatrix.h"
#include "pzstepsolver.h"


#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
#endif

struct locconfig
{
    int porder;
    int refskeleton;
    int nelxcount;
    int nelx;
    int neq;
    REAL delt;
    int postprocfreq;
    int64_t nsteps;
};

locconfig LocalConfig;

void SubstituteMaterialObjects(TPZCompMesh *cmesh);

void InitializeSolution(TPZCompMesh *cmesh);

//    Compute a number of timesteps in parabolic analysis
void SolveParabolicProblem(TPZAnalysis *an, REAL delt, int nsteps, int numthreads);



int main(int argc, char *argv[])
{
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    bool scalarproblem = true;

    int maxnelxcount = 5;
    int maxrefskeleton = 3;
    int maxporder = 4;
    int counter = 1;
    
#ifdef _AUTODIFF
    TimeLaplaceExact.fProblemType = TLaplaceExampleTimeDependent::ECos;
#endif
    for ( int POrder = 2; POrder < maxporder; POrder += 1)
    {
        for (int irefskeleton = 0; irefskeleton < maxrefskeleton; irefskeleton++)
        {
            if (POrder == 3 && !scalarproblem) {
                maxnelxcount = 3;
            }
            for(int nelxcount = 1; nelxcount < maxnelxcount; nelxcount += 1)
            {
                int nelx = 2 << (nelxcount-1);
                bool useexact = true;
                

                TPZCompMesh *SBFem = SetupSquareMesh(nelx,irefskeleton,POrder, scalarproblem,useexact);
                {
                    TPZMaterial *mat = SBFem->FindMaterial(Ebc1);
                    TPZBndCond *bnd = dynamic_cast<TPZBndCond *>(mat);
                    bnd->SetType(0);
#ifdef _AUTODIFF
                    bnd->TPZMaterial::SetForcingFunction(TimeLaplaceExact.TensorFunction());
#endif
                }
                {
                    TPZMaterial *mat = SBFem->FindMaterial(Ebc2);
                    TPZBndCond *bnd = dynamic_cast<TPZBndCond *>(mat);
                    bnd->SetType(0);
#ifdef _AUTODIFF
                    bnd->TPZMaterial::SetForcingFunction(TimeLaplaceExact.TensorFunction());
#endif
                }
                {
                    TPZMaterial *mat = SBFem->FindMaterial(Ebc3);
                    TPZBndCond *bnd = dynamic_cast<TPZBndCond *>(mat);
                    bnd->SetType(0);
#ifdef _AUTODIFF
                    bnd->TPZMaterial::SetForcingFunction(TimeLaplaceExact.TensorFunction());
#endif
                }
                {
                    TPZMaterial *mat = SBFem->FindMaterial(Ebc4);
                    TPZBndCond *bnd = dynamic_cast<TPZBndCond *>(mat);
                    bnd->SetType(0);
#ifdef _AUTODIFF
                    bnd->TPZMaterial::SetForcingFunction(TimeLaplaceExact.TensorFunction());
#endif
                }
                {
                    TPZMaterial *mat = SBFem->FindMaterial(ESkeleton);
                    TPZBndCond *bnd = dynamic_cast<TPZBndCond *>(mat);
                    bnd->SetType(0);
#ifdef _AUTODIFF
                    bnd->TPZMaterial::SetForcingFunction(TimeLaplaceExact.TensorFunction());
#endif
                }

                
#ifdef LOG4CXX
                if(logger->isDebugEnabled())
                {
                    std::stringstream sout;
                    SBFem->Print(sout);
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
                
                std::cout << "nelx = " << nelx << std::endl;
                std::cout << "irefskeleton = " << irefskeleton << std::endl;
                std::cout << "POrder = " << POrder << std::endl;
                LocalConfig.porder = POrder;
                LocalConfig.refskeleton = irefskeleton;
                LocalConfig.nelxcount = nelxcount;
                LocalConfig.nelx = nelx;
                LocalConfig.neq = SBFem->NEquations();
                LocalConfig.delt = 1./20000.;
                LocalConfig.postprocfreq = 2000;
                LocalConfig.nsteps = 20001;
#ifdef _AUTODIFF
                TimeLaplaceExact.fTime = 0.;
                TimeLaplaceExact.fDelt = LocalConfig.delt;
#endif

                // Visualization of computational meshes
                bool mustOptimizeBandwidth = true;
                TPZAnalysis * Analysis = new TPZAnalysis(SBFem,mustOptimizeBandwidth);
//                Analysis->SetStep(counter++);
                std::cout << "neq = " << LocalConfig.neq << std::endl;
                int numthreads = 0;
#ifdef _AUTODIFF
                Analysis->SetExact(TimeLaplace_exact);
#endif
                SolveParabolicProblem(Analysis, LocalConfig.delt, LocalConfig.nsteps, numthreads);
                
                
                delete Analysis;
                delete SBFem;
                //                exit(-1);
            }
            //            exit(-1);
        }
    }
    std::cout << "Check:: Calculation finished successfully" << std::endl;
    return EXIT_SUCCESS;
}






void UniformRefinement(TPZGeoMesh *gMesh, int nh)
{
    for ( int ref = 0; ref < nh; ref++ ){
        TPZVec<TPZGeoEl *> filhos;
        int64_t n = gMesh->NElements();
        for ( int64_t i = 0; i < n; i++ ){
            TPZGeoEl * gel = gMesh->ElementVec() [i];
            if (gel->Dimension() == 2 || gel->Dimension() == 1) gel->Divide (filhos);
        }//for i
    }//ref
}

void InitializeSolution(TPZCompMesh *cmesh)
{
    // create a H1 projection material object
    // set the exact solution
    // project the solution
    
}

void SwitchComputationMode(TPZCompMesh *cmesh, TPZSBFemElementGroup::EComputationMode mode, REAL delt)
{
    int64_t nel = cmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZSBFemElementGroup *elgr = dynamic_cast<TPZSBFemElementGroup *>(cel);
        if(!elgr) continue;
        switch (mode) {
            case TPZSBFemElementGroup::EStiff:
                elgr->SetComputeStiff();
                break;
            case TPZSBFemElementGroup::EMass:
                elgr->SetComputeTimeDependent(delt);
                break;
            case TPZSBFemElementGroup::EOnlyMass:
                elgr->SetComputeOnlyMassMatrix();
                break;
            default:
                DebugStop();
                break;
        }
    }
}

/// set the timestep of all SBFem Element groups
void SetSBFemTimestep(TPZCompMesh *CMesh, REAL delt)
{
    int64_t nel = CMesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = CMesh->Element(el);
        TPZSBFemElementGroup *elgr = dynamic_cast<TPZSBFemElementGroup *>(cel);
        if (!elgr) {
            continue;
        }
        if (delt > 0.) {
            elgr->SetComputeTimeDependent(delt);
        } else
        {
            elgr->SetComputeOnlyMassMatrix();
        }
    }
}


void PostProcess(TPZAnalysis *Analysis, int step)
{
    TPZManVector<REAL,10> errors;
    std::cout << "Compute errors\n";

    Analysis->PostProcessError(errors);



    std::stringstream sout;
    sout << "../ParabolicSolutionErrors.txt";

    std::ofstream results(sout.str(),std::ios::app);
    results.precision(15);
    results << "(* nx " << LocalConfig.nelx << " numrefskel " << LocalConfig.refskeleton << " " << " POrder " << LocalConfig.porder << " neq " << LocalConfig.neq << "*)" << std::endl;
    TPZFMatrix<double> errmat(1,6);
    for(int i=0;i<3;i++) errmat(0,i) = errors[i]*1.e6;
    errmat(0,3) = LocalConfig.neq;
    errmat(0,4) = LocalConfig.nelx;
#ifdef _AUTODIFF
    errmat(0,5) = TimeLaplaceExact.fTime;
#endif
    std::stringstream varname;
    varname << "Errmat[[" << step << "]][[" << LocalConfig.nelxcount << "]][[" << LocalConfig.refskeleton+1 << "]][[" << LocalConfig.porder << "]] = (1/1000000)*";
    errmat.Print(varname.str().c_str(),results,EMathematicaInput);

}

//    Compute a number of timesteps in parabolic analysis
void SolveParabolicProblem(TPZAnalysis *an, REAL delt, int nsteps, int numthreads)
{
    TPZCompMesh *Cmesh = an->Mesh();
    
#ifdef _AUTODIFF
    TimeLaplaceExact.fDelt = delt;
#endif
    
    SetSBFemTimestep(Cmesh, delt);
    //    TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(Cmesh);
#ifndef USING_MKL
    TPZSkylineStructMatrix strmat(Cmesh);
#else
    TPZSymetricSpStructMatrix strmat(Cmesh);
#endif
    
    int64_t neq = Cmesh->NEquations();
    
    if(neq > 20000)
    {
        std::cout << "Entering Assemble Equations\n";
        std::cout.flush();
    }
#ifdef USING_BOOST
    boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
#endif
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);
    an->SetSolver(step);
    
    std::set<int> matids;
    matids.insert(Emat1);
    matids.insert(Ebc1);
    matids.insert(Ebc2);
    matids.insert(Ebc3);
    matids.insert(Ebc4);
    strmat.SetMaterialIds(matids);
    strmat.SetNumThreads(numthreads);
    an->SetStructuralMatrix(strmat);
    an->Assemble();
    if(0)
    {
        std::ofstream andrade("KM.nb");
        andrade.precision(16);
        an->Solver().Matrix()->Print("KM = ",andrade,EMathematicaInput);
        an->Rhs().Print("Rhs = ",andrade,EMathematicaInput);
        std::cout << "KM printed\n";
    }
    TPZAutoPointer<TPZMatrix<STATE> > stiff = an->Solver().Matrix();
    TPZFMatrix<STATE> rhs = an->Rhs();
    
    matids.clear();
    matids.insert(Emat1);
    strmat.SetMaterialIds(matids);
    an->SetStructuralMatrix(strmat);
    
    SetSBFemTimestep(Cmesh, 0.);
    
    an->Solver().ResetMatrix();
    an->Assemble();
    TPZAutoPointer<TPZMatrix<STATE> > mass =  an->Solver().Matrix();
    
    if(0)
    {
        std::ofstream andrade("Mass.nb");
        andrade.precision(16);
        an->Solver().Matrix()->Print("Mass =",andrade,EMathematicaInput);
        std::cout << "Mass printed\n";
    }

    an->Solver().ResetMatrix();
    
    
    an->Solver().SetMatrix(stiff);
    // project the initial solution
    if(an->GetStep() == 0)
    {
        TPZAnalysis an2(an->Mesh(),false);
#ifndef USING_MKL
        TPZSkylineStructMatrix strmat(Cmesh);
#else
        TPZSymetricSpStructMatrix strmat(Cmesh);
#endif
        strmat.SetNumThreads(numthreads);
        std::set<int> matids;
        matids.insert(ESkeleton);
        matids.insert(Ebc1);
        matids.insert(Ebc2);
        matids.insert(Ebc3);
        matids.insert(Ebc4);
        strmat.SetMaterialIds(matids);
        an2.SetStructuralMatrix(strmat);
        TPZStepSolver<STATE> step;
        step.SetDirect(ECholesky);
        an2.SetSolver(step);
        an2.Run();
        if(1)
        {
            TPZStack<std::string> vecnames,scalnames;
            // scalar
            scalnames.Push("State");
            an2.DefineGraphMesh(2, scalnames, vecnames, "InitialSolution.vtk");
            an2.PostProcess(2);
        }
        an->LoadSolution(an2.Solution());
        std::cout << "compmesh solution norm " << Norm(Cmesh->Solution()) << std::endl;
    }
    
    
#ifdef USING_BOOST
    boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
#endif
    
    
#ifdef USING_BOOST
    std::cout << "Time for assembly " << t2-t1 << std::endl;
#endif
    
    for (int istep = 0; istep < nsteps; istep++)
    {
        
        if(istep%LocalConfig.postprocfreq == 0)
        {
            std::cout << "\n";
            int postprocindex = istep/LocalConfig.postprocfreq + 1;
            PostProcess(an, postprocindex);
        }
        /**
         * @brief It computes z = beta * y + alpha * opt(this)*x but z and x can not overlap in memory.
         * @param x Is x on the above operation
         * @param y Is y on the above operation
         * @param z Is z on the above operation
         * @param alpha Is alpha on the above operation
         * @param beta Is beta on the above operation
         * @param opt Indicates if is Transpose or not
         */
        //        virtual void MultAdd(const TPZFMatrix<TVar> & x,const TPZFMatrix<TVar>& y, TPZFMatrix<TVar>& z,
        //                             const TVar alpha=1., const TVar beta = 0., const int opt = 0) const;
        TPZFMatrix<STATE> rhstimestep;
        mass->MultAdd(an->Solution(), rhs, rhstimestep, 1./delt, 1.);
        an->Rhs() = rhstimestep;
        an->Solve();
        if(istep%LocalConfig.postprocfreq == 0)
        {
            std::stringstream sout;
            sout << "../Parabolic_nelx" << LocalConfig.nelx << "_p" << LocalConfig.porder << "_refsk" <<
                LocalConfig.refskeleton << ".vtk";
            TPZStack<std::string> vecnames,scalnames;
            // scalar
            scalnames.Push("State");
            an->DefineGraphMesh(2, scalnames, vecnames, sout.str());
            an->PostProcess(2);
        }
#ifdef _AUTODIFF
        TimeLaplaceExact.fTime += delt;
#endif
//        std::cout << "*";
    }
}

