#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "Common.h"

#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzgeoelside.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"

#include "pzgengrid.h"
#include "TPZMatElasticity2D.h"
#include "TPZMatLaplacian.h"
#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "TPZLagrangeMultiplier.h"
#include "pzmat1dlin.h"

#include "pzanalysis.h"

#include "pzelasmat.h"
#include "pzmultiphysicselement.h"
#include "pzmultiphysicscompel.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzlog.h"
#include <iostream>
#include <string>
#include "TPZVTKGeoMesh.h"
#include "pzfunction.h"
#include "pzmultiphysicselement.h"
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"
#include "TPZInterfaceEl.h"
#include "TPZSBFemVolume.h"
#include "TPZSBFemElementGroup.h"

#include "TPZSSpStructMatrix.h"

#include "TPZBuildSBFem.h"

#include "tpzarc3d.h"
#include "tpzgeoblend.h"

#include "tpzgeoelrefpattern.h"

#include "JSON.hpp"
void rect_mesh();

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

#include <cmath>
#include <set>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
#endif


REAL mult[] = {10./45.,9./45.,8./45.,7./45.,6./45.,5./45.};

void SingularNeumann(const TPZVec<REAL> &x, TPZVec<STATE> &val)
{
    REAL Lambda0 = 2./3.;
    REAL r = sqrt(x[0]*x[0]+x[1]*x[1]);
    REAL theta = atan2(x[1],x[0]);
    if(theta < 0.) theta += 2.*M_PI;
    val[0] = 0;
    for (int i=0; i<6; i++) {
        REAL Lambda = Lambda0*(i+1);
        val[0] += mult[i]*Lambda*pow(r,Lambda-1.)*cos(Lambda*theta);
    }
    std::cout << " x " << x << " theta " << theta << " val " << val[0] << std::endl;
}

void Singular_exact(const TPZVec<REAL> &x, TPZVec<STATE> &val, TPZFMatrix<STATE> &deriv)
{
    REAL Lambda0 = 2./3.;
    REAL r = sqrt(x[0]*x[0]+x[1]*x[1]);
    REAL theta = atan2(x[1],x[0]);
    if (theta<0.) {
        theta += 2.*M_PI;
    }
    
    val[0] = 0;
    deriv.Zero();
    for (int i=0; i<6; i++) {
        REAL Lambda = Lambda0*(i+1);
        val[0] += mult[i]*pow(r,Lambda)*cos(Lambda*theta);
        deriv(0,0) += mult[i]*(Lambda*pow(r,Lambda-2.)*x[0]*cos(Lambda*theta)+pow(r,Lambda-2)*(Lambda)*sin(Lambda*theta)*(x[1]));
        deriv(1,0) += mult[i]*(Lambda*pow(r,Lambda-2.)*x[1]*cos(Lambda*theta)-pow(r,Lambda-2)*(Lambda)*sin(Lambda*theta)*(x[0]));
    }
    
    
}




int main(int argc, char *argv[])
{
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    int minrefskeleton = 2;
    int maxrefskeleton = 3;
    int minporder = 2;
    int maxporder = 9;
    int counter = 1;
    for (int irefskeleton = minrefskeleton; irefskeleton < maxrefskeleton; irefskeleton++)
    {
        for ( int POrder = minporder; POrder < maxporder; POrder += 1)
        {

            REAL angle = M_PI*6./4.;
            TPZCompMesh *SBFem = SetupOneArc(irefskeleton,POrder,angle);
            {
                TPZMaterial *BCond2 = SBFem->FindMaterial(Ebc2);
                TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(SingularNeumann);
                TPZAutoPointer<TPZFunction<STATE> > autodummy = dummy;
                BCond2->SetForcingFunction(autodummy);
                TPZBndCond *BC1 = dynamic_cast<TPZBndCond *>(SBFem->FindMaterial(Ebc1));
                BC1->Val2()(0,0) = 1;
            }
            TPZSBFemElementGroup *celgrp = 0;
            long nel = SBFem->NElements();
            for (long el=0; el<nel; el++) {
                TPZSBFemElementGroup *cel = dynamic_cast<TPZSBFemElementGroup *>(SBFem->Element(el));
                if(cel)
                {
                    celgrp = cel;
                    break;
                }
            }
            
            
            std::cout << "irefskeleton = " << irefskeleton << std::endl;
            std::cout << "POrder = " << POrder << std::endl;
            
            // Visualization of computational meshes
            bool mustOptimizeBandwidth = true;
            TPZAnalysis * Analysis = new TPZAnalysis(SBFem,mustOptimizeBandwidth);
            Analysis->SetStep(counter++);
            std::cout << "neq = " << SBFem->NEquations() << std::endl;
            SolveSist(Analysis, SBFem);
            
            
            
            
            std::cout << "Post processing\n";
            Analysis->SetExact(Singular_exact);
            
            
            long neq = SBFem->Solution().Rows();
            
            if(1)
            {
                TPZStack<std::string> vecnames,scalnames;
                // scalar
                scalnames.Push("State");
                Analysis->DefineGraphMesh(2, scalnames, vecnames, "../SingularSolution.vtk");
                int res = POrder+1;
                if (res >5) {
                    res = 5;
                }
                Analysis->PostProcess(res);
            }
            
            if(0)
            {
                std::ofstream out("../CompMeshWithSol.txt");
                SBFem->Print(out);
            }
            
            TPZManVector<REAL> errors(3,0.);
            Analysis->PostProcessError(errors);
            
            
            
            std::stringstream sout;
            sout << "../SingularSolution.txt";
            
            std::ofstream results(sout.str(),std::ios::app);
            results.precision(15);
            results << "(* numrefskel " << irefskeleton << " " << " POrder " << POrder << " neq " << neq << "*)" << std::endl;
            TPZFMatrix<double> errmat(1,3);
            for(int i=0;i<3;i++) errmat(0,i) = errors[i]*1.e6;
            std::stringstream varname;
            varname << "Errmat[[" << POrder << "][" << irefskeleton+1 << "]] = (1/1000000)*";
            errmat.Print(varname.str().c_str(),results,EMathematicaInput);
            
            if(0)
            {
                std::multimap<REAL,REAL> eigmap;
                TPZManVector<double> eigval = celgrp->EigenvaluesReal();
                TPZFMatrix<double> coef = celgrp->CoeficientsReal();
                for (int i=0; i<eigval.size(); i++) {
                    eigmap.insert(std::pair<double,double>(eigval[i],coef(i,0)));
                }
                for (std::multimap<double, double>::reverse_iterator it = eigmap.rbegin(); it!=eigmap.rend(); it++) {
                    results << it->first << "|" << it->second << " ";
                }
            }
            //results << std::endl;
            //results << celgrp->EigenValues() << std::endl;
            
            std::cout << "Plotting shape functions\n";
            if(0 && irefskeleton == 0)
            {
                int numshape = 25;
                if (numshape > SBFem->NEquations()) {
                    numshape = SBFem->NEquations();
                }
                TPZVec<long> eqindex(numshape);
                for (int i=0; i<numshape; i++) {
                    eqindex[i] = i;
                }
                Analysis->ShowShape("Singular.vtk", eqindex);
            }
            
            delete Analysis;
            delete SBFem;
            //                exit(-1);
        }
        //            exit(-1);
    }
    std::cout << "Check:: Calculation finished successfully" << std::endl;
    return EXIT_SUCCESS;
}





void UniformRefinement(TPZGeoMesh *gMesh, int nh)
{
    for ( int ref = 0; ref < nh; ref++ ){
        TPZVec<TPZGeoEl *> filhos;
        long n = gMesh->NElements();
        for ( long i = 0; i < n; i++ ){
            TPZGeoEl * gel = gMesh->ElementVec() [i];
            if (gel->Dimension() == 2 || gel->Dimension() == 1) gel->Divide (filhos);
        }//for i
    }//ref
}


