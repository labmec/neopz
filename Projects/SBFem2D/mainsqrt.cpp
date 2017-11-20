#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "Common.h"

#include "TPZGeoLinear.h"
#include "pzgeopoint.h"
#include "tpzarc3d.h"
#include "tpzgeoblend.h"

#include "TPZVTKGeoMesh.h"

#include "pzbndcond.h"
#include "TPZMatLaplacian.h"

#include "pzfunction.h"
#include "TPZSBFemElementGroup.h"
#include "TPZBuildSBFem.h"

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif


#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
#endif


REAL mult[] = {1.,10./45.,9./45.,8./45.,7./45.,6./45.,5./45.};

void SingularNeumann(const TPZVec<REAL> &x, TPZVec<STATE> &val)
{
    REAL Lambda0 = 1./2.;
    REAL r = sqrt(x[0]*x[0]+x[1]*x[1]);
    REAL theta = atan2(x[1],x[0]);
    if(theta < 0.) theta += 2.*M_PI;
    val[0] = 0;
    for (int i=0; i<1; i++) {
        REAL Lambda = Lambda0+i;
        val[0] += mult[i]*Lambda*pow(r,Lambda-1.)*sin(Lambda*theta);
    }
    std::cout << " x " << x << " theta " << theta << " val " << val[0] << std::endl;
}

void Singular_exact(const TPZVec<REAL> &x, TPZVec<STATE> &val, TPZFMatrix<STATE> &deriv)
{
    REAL Lambda0 = 1./2;
    REAL r = sqrt(x[0]*x[0]+x[1]*x[1]);
    REAL theta = atan2(x[1],x[0]);
    if (theta<0.) {
        theta += 2.*M_PI;
    }
    
    val[0] = 0;
    deriv.Zero();
    for (int i=0; i<1; i++) {
        REAL Lambda = Lambda0+i;
        val[0] += mult[i]*pow(r,Lambda)*sin(Lambda*theta);
        deriv(0,0) += mult[i]*(Lambda*pow(r,Lambda-2.)*x[0]*sin(Lambda*theta)-pow(r,Lambda-2)*(Lambda)*cos(Lambda*theta)*(x[1]));
        deriv(1,0) += mult[i]*(Lambda*pow(r,Lambda-2.)*x[1]*sin(Lambda*theta)+pow(r,Lambda-2)*(Lambda)*cos(Lambda*theta)*(x[0]));
    }
    
    
}

TPZCompMesh *SetupOneArcWithRestraint(int numrefskeleton, int porder, REAL angle);



int main(int argc, char *argv[])
{
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    int minrefskeleton = 0;
    int maxrefskeleton = 5;
    int minporder = 1;
    int maxporder = 9;
    int counter = 1;
    for (int irefskeleton = minrefskeleton; irefskeleton < maxrefskeleton; irefskeleton++)
    {
        for ( int POrder = minporder; POrder < maxporder; POrder += 1)
        {

            REAL angle = M_PI;
            TPZCompMesh *SBFem = SetupOneArcWithRestraint(irefskeleton,POrder, angle);
            {
                TPZMaterial *BCond2 = SBFem->FindMaterial(Ebc2);
                TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(SingularNeumann);
                TPZAutoPointer<TPZFunction<STATE> > autodummy = dummy;
                BCond2->SetForcingFunction(autodummy);
                TPZBndCond *BC1 = dynamic_cast<TPZBndCond *>(SBFem->FindMaterial(Ebc1));
                BC1->Val2()(0,0) = 0;
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
            
            TPZManVector<STATE> errors(3,0.);
            
            long neq = SBFem->Solution().Rows();
            
            if(1)
            {
                TPZStack<std::string> vecnames,scalnames;
                // scalar
                scalnames.Push("State");
                Analysis->DefineGraphMesh(2, scalnames, vecnames, "../SqrtSolution.vtk");
                Analysis->PostProcess(3);
            }
            
            if(0)
            {
                std::ofstream out("../CompMeshWithSol.txt");
                SBFem->Print(out);
            }
            
            Analysis->PostProcessError(errors);
            
            
            
            std::stringstream sout;
            sout << "../SqrtSolution.txt";
            
            std::ofstream results(sout.str(),std::ios::app);
            results.precision(15);
            results << "(* numrefskel " << irefskeleton << " " << " POrder " << POrder << " neq " << neq << "*)" << std::endl;
            TPZFMatrix<double> errmat(1,3);
            for(int i=0;i<3;i++) errmat(0,i) = errors[i]*1.e6;
            std::stringstream varname;
            varname << "Errmat" << POrder <<  irefskeleton << " = (1/1000000)*";
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
            if(POrder == maxporder-1 && irefskeleton == 0)
            {
                int numshape = 25;
                if (numshape > SBFem->NEquations()) {
                    numshape = SBFem->NEquations();
                }
                TPZVec<long> eqindex(numshape);
                for (int i=0; i<numshape; i++) {
                    eqindex[i] = i;
                }
                Analysis->SetStep(0);
                Analysis->ShowShape("Sqrt.vtk", eqindex);
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


TPZCompMesh *SetupOneArcWithRestraint(int numrefskeleton, int porder, REAL angle)
{
    TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh;
    gmesh->SetDimension(2);
    gmesh->NodeVec().Resize(4);
    TPZManVector<REAL,3> co(3,0.);
    gmesh->NodeVec()[0].Initialize(co, gmesh);
    co[0] = 1.;
    gmesh->NodeVec()[1].Initialize(co, gmesh);
    co.Fill(0.);
    co[0] = cos(angle);
    co[1] = sin(angle);
    gmesh->NodeVec()[2].Initialize(co, gmesh);
    co.Fill(0.);
    co[0] = cos(angle/2.);
    co[1] = sin(angle/2.);
    gmesh->NodeVec()[3].Initialize(co, gmesh);
    co.Fill(0.);
    TPZManVector<long,4> nodeindex(1,0);
    
    nodeindex[0] = 1;
    long elementid = 1;
    gmesh->CreateGeoElement(EPoint, nodeindex, Ebc1, elementid);
    
    nodeindex.Resize(3);
    // Definition of Arc coordenates
    // Create Geometrical Arc #1
    nodeindex[0] = 1;
    nodeindex[1] = 2;
    nodeindex[2] = 3;
    elementid = 1;
    TPZGeoEl *arc = new TPZGeoElRefPattern < pzgeom::TPZArc3D > (nodeindex, Ebc2, gmesh,elementid);
    
    nodeindex.Resize(4);
    nodeindex[0] = 1;
    nodeindex[1] = 2;
    nodeindex[2] = 0;
    nodeindex[3] = 0;
    elementid = 2;
    TPZGeoEl *gblend = new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> > (nodeindex, EGroup, gmesh,elementid);
    
    nodeindex.Resize(2);
    nodeindex[0] = 1;
    nodeindex[1] = 0;
    TPZGeoEl *oned = new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> (nodeindex, EGroup+1, gmesh, elementid);
    
    gmesh->BuildConnectivity();
    
    //    gmesh->Print(std::cout);
    TPZManVector<REAL,3> xi(1),x(3);
    for (REAL s=-1.; s<=1.; s+= 1./10.) {
        xi[0] = s;
        arc->X(xi, x);
        std::cout << "xi " << xi << " x " << x << std::endl;
    }
    std::map<int,int> matmap;
    matmap[EGroup] = Emat1;
    matmap[EGroup+1] = Emat2;
    TPZBuildSBFem build(gmesh,ESkeleton,matmap);
    
    TPZStack<long> elids;
    elids.Push(gblend->Index());
    elids.Push(oned->Index());
    build.AddPartition(elids, 0);
    
    //        AddSkeletonElements(gmesh);
    /// generate the SBFem elementgroups
    
    /// put sbfem pyramids into the element groups
    TPZCompMesh *SBFem = new TPZCompMesh(gmesh);
    SBFem->SetDefaultOrder(porder);
    
    // problemtype - 1 laplace equation
    int problemtype  = 1;
    InsertMaterialObjects(SBFem,problemtype);
    
    
    TPZMaterial *mat1 = SBFem->FindMaterial(Emat1);
    TPZMaterial *mat2 = mat1->NewMaterial();
    TPZMatLaplacian *mat2lapl = dynamic_cast<TPZMatLaplacian *>(mat2);
    mat2->SetId(Emat2);
    mat2lapl->SetDimension(1);
    mat2lapl->SetParameters(1.e9, 0);
    SBFem->InsertMaterialObject(mat2);
    
    std::set<int> volmatids,boundmatids;
    volmatids.insert(Emat1);
    volmatids.insert(Emat2);
    boundmatids.insert(Ebc1);
    boundmatids.insert(Ebc2);
    boundmatids.insert(ESkeleton);
    build.DivideSkeleton(numrefskeleton, volmatids);

    build.BuildComputationMesh(*SBFem,volmatids,boundmatids);
    
    {
        std::ofstream outg("GMesh.txt");
        gmesh->Print(outg);
        std::ofstream outc("CMesh.txt");
        SBFem->Print(outc);
        std::ofstream out("Geometry.vtk");
        TPZVTKGeoMesh vtk;
        vtk.PrintGMeshVTK(gmesh, out,true);
    }
    return SBFem;
    
}
