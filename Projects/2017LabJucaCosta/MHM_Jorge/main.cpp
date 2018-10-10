#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "TPZInterfaceEl.h"
#include "pzgeoelside.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "tpzcompmeshreferred.h"
#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "pzanalysis.h"

#include "TPZParSkylineStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "pzstepsolver.h"
#include "pzstrmatrix.h"
#include "pzfstrmatrix.h"
#include "TPZFrontNonSym.h"
#include "TPZFrontSym.h"
#include "TPBSpStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzbstrmatrix.h"

#include "pzpoisson3d.h"
#include "pzpoisson3dreferred.h"
#include "mixedpoisson.h"
#include "pzelasmat.h"
#include "pzelasthybrid.h"
#include "pzmat1dlin.h"
#include "TPZVecL2.h"
#include "TPZMatLaplacianHybrid.h"
#include "TPZLagrangeMultiplier.h"
#include "pzmatmixedpoisson3d.h"

#include "pzbuildmultiphysicsmesh.h"
#include "pzelementgroup.h"
#include "TPZCompMeshTools.h"
#include "pzcondensedcompel.h"
#include "pzfunction.h"
#include "pzgraphmesh.h"
#include "pzfmatrix.h"

#include "pzlog.h"

#include "TPZVTKGeoMesh.h"
#include "pzvisualmatrix.h"
#include "pzgengrid.h"
#include "TPZExtendGridDimension.h"
#include "pzcheckgeom.h"

#include "TPZMHMeshControl.h"
#include "TPZMHMixedMeshControl.h"
#include "TPZMHMixedHybridMeshControl.h"

#include <iostream>
#include <string>

#include <math.h>
#include <set>

using namespace std;

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mainskeleton"));
#endif


/** PHILIPPE IMPLEMENTATIONS */

/// Insert material objects for the MHM-H(div) solution
void InsertMaterialObjects(TPZMHMixedMeshControl &control);
/// Insert material objects for the MHM-H(div) solution
//void InsertMaterialObjects(TPZMHMixedHybridMeshControl &control);

/// unwrap de TPZCondensedCompel and TPZElementGroup elements
void UnwrapMesh(TPZCompMesh *cmesh);
/// function that returns the permeability for a given coordinate
void Permeability(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE> &diff);

struct TRunConfig
{
    int nelxcoarse = -1;
    int nelycoarse = -1;
    int numHDivisions = 0;
    int pOrderInternal = 1;
    int numDivSkeleton = 0;
    int pOrderSkeleton = 1;
    int Hybridize = 0;
    int Condensed = 1;
    int LagrangeMult = 0;
    int newline = 0;
    int n_threads = 0;
    
    /// number of equations when not condensing anything
    int64_t fGlobalSystemSize = -1;
    /// number of equations considering local condensation
    int64_t fGlobalSystemWithLocalCondensationSize = -1;
    /// number of equations of the global system
    int64_t fNumeq = -1;
    
    REAL fDeltaT = 1.;
    
    /// number of timesteps
    int64_t nTimeSteps = 10;
    
    std::ostream &ConfigPrint(std::ostream &out)
    {
        out << nelxcoarse << "x" << nelycoarse << "_HSkel" << numDivSkeleton << "_pSkel" << pOrderSkeleton << "_HDiv" << numHDivisions << "_pInt" << pOrderInternal;
        return out;
    }
};

void ErrorHDiv2(TPZCompMesh *hdivmesh, std::ostream &out, TPZVec<STATE> &errorHDiv);
void ErrorH1(TPZCompMesh *l2mesh, std::ostream &out, STATE &errorL2, STATE &errordu);

void SolveProblem(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<TPZAutoPointer<TPZCompMesh> > compmeshes, std::string prefix, TRunConfig config);

/**
 JORGE IMPLEMENTATIONS
 */

bool UseTriangle = true;

//To create quadratic geometric mesh in quadrilaterals or triangles from x0 until x1
TPZGeoMesh *MalhaGeom(int nelx, int nely,TPZVec<REAL> &x0,TPZVec<REAL> &x1,bool ftriang, bool zigzag, int ndiv,TPZVec<int64_t> &coarseindices);
void ReadPorous(TPZFMatrix<REAL> &porous);

// Forcing function and exact solution of the problem
void SolProblema_ParedesThesis(const TPZVec<REAL> &pt, TPZVec<STATE> &u, TPZFMatrix<STATE> &flux);
void ForcingF_ParedesThesis(const TPZVec<REAL> &pt, TPZVec<STATE> &res);

//Compute pressure jump on all faces from right pressure and left pressure on center point of the face
bool ComputePressureJumpOnFaces(TPZCompMesh *cmesh,int matid,STATE &Error,STATE &ErrorNi);

// Constants to estimation error "a posteriori" (Tese Paredes - LNCC
REAL Cl = 3.;   //Depend only on l, (espace Lambda_l). The value Cl = 3 was taked from Paredes's thesis to l = 0. (pag 82)
REAL Cmin = 1.;   //Cmin and Cmax. From eliptic problem. The value 1 is for K = I (identity matrix)
REAL Cmax = 1.;

REAL alpha = 1.;   // Factor to K tensor

const int matId = 1;

/*
 END JORGE IMPLEMENTATIONS
 **/

const int neumann = 1;
int const bc1=-1;

bool permeabilityisfuncion = false;

const int PorousSize = 200;
TPZFMatrix<REAL> gPorous(PorousSize,PorousSize,0.);


// PROBLEM TO MIXED POISSON FORMULATION
int main(int argc, char *argv[]) {

#ifdef LOG4CXX
    InitializePZLOG();
#endif
    // to uniform refinement
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);

    // CONFIGURATION OF THE MHM METHOD
    TRunConfig Config;
    /// PolynomialOrder - p-order
    Config.pOrderInternal = 1;
    Config.Condensed = 1;
    Config.n_threads = 2;
    Config.pOrderSkeleton = 1;
    Config.numDivSkeleton = 0;
    // to avoid singular internal matrices
    if (Config.numHDivisions == 0 && Config.pOrderInternal < Config.pOrderSkeleton) {
        Config.pOrderInternal = Config.pOrderSkeleton+1;
    }
    Config.Hybridize = 0;
    
    REAL hsize = 1.;   /// FALTA CALCULAR O H PARA CADA MALHA APOS REFINAMENTO H
    TPZManVector<REAL,3> x0(3,0.),x1(3,1.);

    // To print errors calculated
    ofstream saidaerros;
    if(!Config.Hybridize)
        saidaerros.open("../Erro-Misto.txt");
    else
        saidaerros.open("../Erro-Hybridized.txt");

    ReadPorous(gPorous);
    
    // P maximun and h-refinement maximum
    int maxorder = 5;
    int maxhref = 7;
    
    // To store errors computed
    TPZFMatrix<STATE> L2ErrorPrimal(maxhref,maxorder-1,0.);
    TPZFMatrix<STATE> L2ErrorDual(maxhref,maxorder-1,0.);
    TPZFMatrix<STATE> L2ErrorDiv(maxhref,maxorder-1,0.);
    TPZFMatrix<STATE> HDivErrorDual(maxhref,maxorder-1,0.);
    TPZFMatrix<STATE> JumpPressure(maxhref,maxorder-1,0.);
    TPZFMatrix<STATE> JumpPressureErrorNi(maxhref,maxorder-1,0.);
    TPZFMatrix<STATE> EfectivityIndex(maxhref,maxorder-1,0.);
    // To store the convergence rate
    TPZFMatrix<STATE> L2ConvergPrimal(maxhref-1,maxorder-1,0.);
    TPZFMatrix<STATE> L2ConvergDual(maxhref-1,maxorder-1,0.);
    TPZFMatrix<STATE> L2ConvergDiv(maxhref-1,maxorder-1,0.);
    TPZFMatrix<STATE> HDivConverg(maxhref-1,maxorder-1,0.);
    TPZFMatrix<STATE> JumpPressureConverg(maxhref-1,maxorder-1,0.);
    TPZFMatrix<STATE> JumpPressureNiConverg(maxhref-1,maxorder-1,0.);
    // To store auxiliar data for iterations
    TPZFMatrix<int> porders(maxhref,maxorder-1,0.);
    TPZFMatrix<int> numhref(maxhref,maxorder-1,0.);
    TPZFMatrix<int> DofTotal(maxhref,maxorder-1,0.);
    TPZFMatrix<int> DofCond(maxhref,maxorder-1,0.);
	Config.nelxcoarse = 1 << 2;
	Config.nelycoarse = 1 << 2;

    for ( int order=1;order<maxorder;order++) {

		/// numhdiv - number of h-refinements
		Config.numHDivisions = 0;
		Config.pOrderInternal = order;
        Config.pOrderSkeleton = ((order/2 > 0)? (order/2): 1);
        hsize = 0.5;

        for(int href=2;href<maxhref;href++) {
            
			/// numhdiv - number of h-refinements
			Config.numHDivisions++;
			hsize = hsize/2.;
            Config.numDivSkeleton = 0;

            TPZGeoMesh *gmesh = 0;
            // To store indices of the subdomains (initial coarse elements)
            TPZManVector<int64_t> coarseindices;
            //creating geometric mesh with coarse elements, the geometric elements created will be a subdomains
            gmesh = MalhaGeom(Config.nelxcoarse, Config.nelycoarse,x0,x1,UseTriangle,false,Config.numHDivisions,coarseindices);

            TPZAutoPointer<TPZGeoMesh> gmeshpointer(gmesh);

            TPZMHMeshControl *control;
            if(!Config.Hybridize)
                control = new TPZMHMixedMeshControl(gmeshpointer);
            else
                control = new TPZMHMixedHybridMeshControl(gmeshpointer);

            
            // Fill fGeoToMHMDomains into control and creating the skeleton elements from subdomains
            control->DefinePartitionbyCoarseIndices(coarseindices);
            std::set<int> matids;
            matids.insert(matId);
            control->fMaterialIds = matids;
            matids.clear();
            matids.insert(bc1);
            control->fMaterialBCIds = matids;
        
            InsertMaterialObjects(*((TPZMHMixedMeshControl*)control));
        
            control->SetInternalPOrder(Config.pOrderInternal);
            control->SetSkeletonPOrder(Config.pOrderSkeleton);
        
            control->DivideSkeletonElements(Config.numDivSkeleton);

            control->SetHybridize(Config.Hybridize);
        
            bool substructure = (bool) Config.Condensed;
            control->BuildComputationalMesh(substructure);
            // Making material as needed by the problem
            TPZMatMixedPoisson3D *novomat = new TPZMatMixedPoisson3D(matId,2);
            if(!Config.Hybridize)
                ((TPZMHMixedMeshControl*)control)->FluxMesh()->MaterialVec()[matId] = novomat;
            else
                ((TPZMHMixedHybridMeshControl *)control)->FluxMesh()->MaterialVec()[matId] = novomat;

            std::cout << "MHM Hdiv Computational meshes created\n";
            std::cout << "Number of equations MHMixed " << control->CMesh()->NEquations() << std::endl;
            std::string configuration;
    
            int nDofTotal;
            if(!Config.Hybridize)
                nDofTotal = ((TPZMHMixedMeshControl*)control)->FluxMesh()->NEquations() + control->PressureMesh()->NEquations();
            else
                nDofTotal = ((TPZMHMixedHybridMeshControl*)control)->FluxMesh()->NEquations() + control->PressureMesh()->NEquations();

            {
                std::stringstream sout;
                sout << "H" << Config.numHDivisions << "-P" << Config.pOrderInternal;
                configuration = sout.str();
            }

            std::stringstream MHMMixedPref;
            MHMMixedPref << "MHMixed";
            if(Config.LagrangeMult) {
                MHMMixedPref << "_Lagr";
            }
            if (Config.Hybridize) {
                MHMMixedPref << "_Hybr";
            }
            
            //Variables para calculo de erros
            STATE errorPrimalL2;
            STATE errorDuL2;
            TPZVec<STATE> errorsHDiv;
            STATE JumpAsError = 0.;
            STATE ErrorNi = 0.;
        
            // compute the MHM H(div) solution
            Config.fGlobalSystemWithLocalCondensationSize = control->fGlobalSystemWithLocalCondensationSize;
            Config.fGlobalSystemSize = control->fGlobalSystemSize;
            Config.fNumeq = control->fNumeq;
            SolveProblem(control->CMesh(), control->GetMeshes(), MHMMixedPref.str(), Config);

            saidaerros << "H Size = " << hsize << std::endl;
            if(!Config.Hybridize)
                ErrorHDiv2(((TPZMHMixedMeshControl*)control)->FluxMesh().operator->(), saidaerros,errorsHDiv);
            else
                ErrorHDiv2(((TPZMHMixedHybridMeshControl*)control)->FluxMesh().operator->(), saidaerros,errorsHDiv);

            ErrorH1(control->PressureMesh().operator->(),saidaerros,errorPrimalL2,errorDuL2);

            if(ComputePressureJumpOnFaces(control->PressureMesh().operator->(), matId, JumpAsError, ErrorNi)) {
                saidaerros << "Jump of pressure = "    << JumpAsError << "\nError Ni = " << ErrorNi << std::endl;
                saidaerros << std::endl;
            }
            else
                saidaerros << "Jump of pressure couldn't to be computed."    << std::endl << std::endl;

            L2ErrorPrimal(href,order-1) = errorPrimalL2;
            L2ErrorDual(href,order-1) = errorsHDiv[0];
            L2ErrorDiv(href,order-1) = errorsHDiv[1];
            HDivErrorDual(href,order-1) = errorsHDiv[2];
                
            JumpPressure(href,order-1) = JumpAsError;   /// Jorge
            JumpPressureErrorNi(href,order-1) = ErrorNi;
            porders(href,order-1) = order;
            numhref(href,order-1) = href;
            DofTotal(href,order-1) = nDofTotal;
                
            // Computing efectivity index
            EfectivityIndex(href,order-1) = ErrorNi/(Cmin*((1/hsize)*errorPrimalL2+errorDuL2));   // FALTA Norm(KV(p-pa)div
                
            if(control->CMesh().operator->()) delete control->CMesh().operator->();
            if(!Config.Hybridize) {
                if(((TPZMHMixedMeshControl*)control)->FluxMesh().operator->())
                    delete ((TPZMHMixedMeshControl*)control)->FluxMesh().operator->();
            } else {
                if(((TPZMHMixedHybridMeshControl*)control)->FluxMesh().operator->())
                    delete ((TPZMHMixedHybridMeshControl*)control)->FluxMesh().operator->();
            }
            if(control->PressureMesh().operator->()) delete control->PressureMesh().operator->();
            if(control) delete gmesh;
        }
    }
    saidaerros.close();

    return 0;
}

void InsertMaterialObjects(TPZMHMixedMeshControl &control)
{
    TPZCompMesh &cmesh = control.CMesh();
    TPZGeoMesh &gmesh = control.GMesh();
    TPZFMatrix<STATE> val1(2,2,0.), val2Flux(2,1,0.);

    int dim = gmesh.Dimension();
    cmesh.SetDimModel(dim);
    
    TPZCompMesh *MixedFluxPressureCmesh = &cmesh;
    
    // Material medio poroso
    TPZMatMixedPoisson3D * mat = new TPZMatMixedPoisson3D(matId,dim);
    mat->SetPermeability(1.);
    if(permeabilityisfuncion) {
        TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(Permeability, 5);
        dummy->SetPolynomialOrder(0);
        TPZAutoPointer<TPZFunction<STATE> > func(dummy);
        mat->SetPermeabilityFunction(func);
    }
    TPZFMatrix<REAL> Ktensor(3,3,0.);
    TPZFMatrix<REAL> InvK(3,3,0.);
    for(int i=0; i<3; i++) {
        Ktensor(i,i) = alpha;
        InvK(i,i) = 1./alpha;
    }
    
    mat->SetPermeabilityTensor(Ktensor,InvK);

    // exact solution of the problem
    TPZAutoPointer<TPZFunction<STATE> > solexata;
    solexata = new TPZDummyFunction<STATE>(SolProblema_ParedesThesis, 5);
    mat->SetForcingFunctionExact(solexata);
    
    int int_order = 10;
    // setting forcing function of the equation
    TPZAutoPointer<TPZFunction<STATE> > fforce;
    TPZDummyFunction<STATE> *dum;
    dum = new TPZDummyFunction<STATE>(ForcingF_ParedesThesis, 5);
    dum->SetPolynomialOrder(int_order);
    fforce = dum;
    mat->SetForcingFunction(fforce);

    MixedFluxPressureCmesh->InsertMaterialObject(mat);
    
    // one boundary condition
    TPZBndCond * bc = mat->CreateBC(mat, bc1, neumann, val1, val2Flux);
    MixedFluxPressureCmesh->InsertMaterialObject(bc);
}

void Permeability(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE> &diff) {
    int factor = 300;
    int64_t ix = x[0]*factor;
    int64_t iy = x[1]*factor;
    static int count = 0;
    if((fabs(ix-x[0]*factor) < 1.e-6 || fabs(ix-x[1]*factor) < 1.e-6) && count < 100)
    {
        count++;
        std::cout << "probing for a permeability at the interface of two regions\n";
        std::cout << "x = " << x << std::endl;
    }
    if (IsZero(x[1]-1.0)) {
        iy = factor-1;
    }
    if (IsZero(x[0]-1.)) {
        ix = factor-1;
    }

    REAL valporous = gPorous(ix,iy);
    for (int i=0; i<2; i++) {
        diff(i,i) = valporous;
        diff(i,1-i)=0.;
        diff(2+i,i) = 1./valporous;
        diff(2+i,1-i) = 0.;
    }
}

void UnwrapMesh(TPZCompMesh *cmesh) {
    int64_t nel = cmesh->NElements();
    bool change = true;
    while(change) {
        change = false;
        for (int64_t el=0; el<nel; el++) {
            
            TPZCompEl *cel = cmesh->Element(el);
            TPZCondensedCompEl *condense = dynamic_cast<TPZCondensedCompEl *>(cel);
            if (condense) {
                condense->Unwrap();
                change = true;
            }
            cel = cmesh->Element(el);
            TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *>(cel);
            if (elgr) {
                elgr->Unwrap();
                change = true;
            }
        }
    }
}

/** PHILIPPE IMPLEMENTATIONS */

void ReadPorous(TPZFMatrix<REAL> &porous) {
#ifdef MACOSX
    std::ifstream pores("../porous_scaled.txt");
#else
    std::ifstream pores("porous_scaled.txt");
#endif
    for (int j=0; j<PorousSize; j++) {
        for (int i=0; i<PorousSize; i++) {
            pores >> gPorous(i,j);
            if (!pores) {
                DebugStop();
            }
        }
    }
}
void SolProblema_ParedesThesis(const TPZVec<REAL> &pt, TPZVec<STATE> &u, TPZFMatrix<STATE> &flux) {
    
    REAL x = pt[0];
    REAL y = pt[1];
    
    u.Resize(1, 0.);
    flux.Resize(3, 1);
    flux(0,0)=0., flux(1,0)=0., flux(2,0)=0.;
    
    //solucao u
    u[0] = cos(2.*M_PI*x)*cos(2.*M_PI*y);
    
    //fluxo em x (-K(p_x,p_y))                                    /// JORGE ??
    flux(0,0) = 2.*M_PI*sin(2.*M_PI*x)*cos(2.*M_PI*y);
    
    //fluxo em y
    flux(1,0) = 2.*M_PI*cos(2.*M_PI*x)*sin(2.*M_PI*y);
    
    //Solucao do divergente de u
    flux(2,0) = 8.*M_PI*M_PI*u[0]; //valor do divergente
    
    //------- Solucao da pressão: p(x,y) = cos(2*pi*x)*cos(2*pi*y) -----------
    //        u[0] = cos(2*pi*x)*cos(2*pi*y);
    //        flux(0,0) = -K*dp/dx = -2.*pi*sin(2*pi*x)*cos(2*pi*y);
    //        flux(1,0) = -K*dp/dy = -2.*pi*cos(2*pi*x)*sin(2*pi*y);
    //        flux(2,0) = V.(-K*V) = 8.*pi*pi*cos(2*pi*x)*cos(2*pi*y); //valor do divergente -Div ?
}

void ForcingF_ParedesThesis(const TPZVec<REAL> &pt, TPZVec<STATE> &res) {
    
    double x = pt[0];
    double y = pt[1];
    res[0] = 0.;
    
    REAL solp =  cos(2*M_PI*x)*cos(2*M_PI*y);
    
    res[0] = 8.*M_PI*M_PI*solp;
}

void ErrorHDiv2(TPZCompMesh *hdivmesh, std::ostream &out, TPZVec<STATE> &errorHDiv)
{
    hdivmesh->Reference()->ResetReference();
    hdivmesh->LoadReferences();
    
    int64_t nel = hdivmesh->NElements();
    int dim = hdivmesh->Dimension();
    TPZManVector<STATE,10> globerrors(10,0.);
    TPZStack<REAL> vech;
    
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = hdivmesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        
        TPZGeoEl *gel = cel->Reference();
        if (!gel || gel->Dimension() != dim) {
            continue;
        }
        
        TPZManVector<REAL,10> elerror(10,0.);
        cel->EvaluateError(SolProblema_ParedesThesis, elerror, 0);
        int nerr = elerror.size();
        for (int i=0; i<nerr; i++) {
            globerrors[i] += elerror[i]*elerror[i];
        }
    }
    
    out << "L2 Error Norm for flux = "    << sqrt(globerrors[1]) << std::endl;
    errorHDiv.Resize(3,0.);
    errorHDiv[0] = sqrt(globerrors[1]);
    errorHDiv[1] = sqrt(globerrors[2]);
    errorHDiv[2] = sqrt(globerrors[3]);
    out << "L2 Norm for divergence = "    << sqrt(globerrors[2])  <<std::endl;
    out << "Hdiv Norm for flux = "    << sqrt(globerrors[3])  <<std::endl;
    
}
void ErrorH1(TPZCompMesh *l2mesh, std::ostream &out, STATE &errorL2, STATE &errordu)
{
    l2mesh->Reference()->ResetReference();
    l2mesh->LoadReferences();
    
    int64_t nel = l2mesh->NElements();
    TPZManVector<STATE,10> globerrors(10,0.);
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = l2mesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZManVector<REAL,10> elerror(10,0.);
        elerror.Fill(0.);
        cel->EvaluateError(SolProblema_ParedesThesis, elerror, 0);
        
        int nerr = elerror.size();
        globerrors.resize(nerr);
        for (int i=0; i<nerr; i++) {
            globerrors[i] += elerror[i]*elerror[i];
        }
    }
    //out << "Errors associated with L2 or H1 space\n";
    out << "H1 Error Norm = "    << sqrt(globerrors[0]) << std::endl;
    out << "L2 Error Norm = "    << sqrt(globerrors[1]) << std::endl;
    out << "Semi H1 Norm = "    << sqrt(globerrors[2]) << std::endl;
    out << "=============================\n"<<std::endl;
    errorL2 = sqrt(globerrors[1]);
    errordu = sqrt(globerrors[2]);
}
void SolveProblem(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<TPZAutoPointer<TPZCompMesh> > compmeshes, std::string prefix, TRunConfig config)
{
    //calculo solution
    bool shouldrenumber = true;
    TPZAnalysis an(cmesh,shouldrenumber);
#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(cmesh.operator->());
    strmat.SetNumThreads(config.n_threads);
    
#else
    TPZSkylineStructMatrix strmat(cmesh.operator->());
    strmat.SetNumThreads(config.n_threads);
#endif
    
    
    an.SetStructuralMatrix(strmat);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    std::cout << "Assembling\n";
    an.Assemble();
    if(0)
    {
        std::string filename = prefix;
        filename += "_Global.nb";
        std::ofstream global(filename.c_str());
        TPZAutoPointer<TPZStructMatrix> strmat = an.StructMatrix();
        an.Solver().Matrix()->Print("Kg = ",global,EMathematicaInput);
        an.Rhs().Print("Fg = ",global,EMathematicaInput);
    }
    std::cout << "Solving\n";
    an.Solve();
    std::cout << "Finished\n";
    an.LoadSolution(); // compute internal dofs
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(compmeshes, cmesh);
    
#ifdef PZDEBUG
    {
        std::ofstream out(prefix+"_MeshWithSol.txt");
        cmesh->Print(out);
    }
#endif
    
    std::string plotfile;
    std::stringstream sout_geo;
    std::stringstream sout;
    {
        sout << prefix << "Approx_";
        config.ConfigPrint(sout);
        plotfile = sout.str() + "_dim2.vtk";
    }
    {
        sout_geo << prefix << "Geo_";
        config.ConfigPrint(sout_geo) << "_dim2.vtk";
    }
    
    std::ofstream plotfile3(sout_geo.str());
    TPZVTKGeoMesh::PrintGMeshVTK(cmesh.operator->()->Reference(), plotfile3, true);
    
    std::cout << "plotfiles " << " " << plotfile.c_str() << std::endl;
    TPZStack<std::string> scalnames,vecnames;
    TPZMaterial *mat = cmesh->FindMaterial(1);
    if (!mat) {
        DebugStop();
    }
    if (mat->NStateVariables() == 2)
    {
        scalnames.Push("SigmaX");
        scalnames.Push("SigmaY");
        scalnames.Push("TauXY");
        vecnames.Push("Displacement");
    }
    else if(mat->NStateVariables() == 1)
    {
        scalnames.Push("Pressure");
 //       scalnames.Push("Permeability");
        vecnames.Push("Flux");
 //       vecnames.Push("Derivative");
    }
    an.DefineGraphMesh(cmesh->Dimension(), scalnames, vecnames, plotfile);
    int resolution = 0;
    an.PostProcess(resolution,cmesh->Dimension());
}


/** JORGE IMPLEMENTATIONS */

// Auxiliar function to know the maxime index of the comp elements of a mesh
long MaxCompElementsIndex(TPZCompMesh *cmesh) {
    long i, nel = cmesh->NElements();
    long maxelindex = 0;
    for(i=0L;i<nel;i++) {
        long index = cmesh->ElementVec()[i]->Index();
        if(index>maxelindex)
            maxelindex = index;
    }
    return maxelindex;
}

TPZGeoMesh *MalhaGeom(int nelx, int nely, TPZVec<REAL> &x0,TPZVec<REAL> &x1,bool ftriang, bool zigzag, int ndiv,TPZVec<int64_t> &coarseindices){
    
    int dimension = 2;
    if(x0.NElements()<dimension || x1.NElements()<dimension)
        return NULL;
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    
    TPZManVector<REAL,3> xcorner(3,0.);
    TPZManVector<int,3> nx(2);
    nx[0] = nelx;
    nx[1] = nely;
    TPZGenGrid gengrid(nx,x0,x1);
    gengrid.SetRefpatternElements(true);
    
    if (ftriang && zigzag) {
        std::cout << "Zigzag meshes cannot be created with triangular meshes\n";
        DebugStop();
    }
    
    if (ftriang) {
        gengrid.SetElementType(ETriangle);
    }
    if (zigzag) {
        gengrid.SetDistortion(0.25);
    }
    //    gengrid.SetDistortion(0.75);
    gengrid.Read(gmesh,matId);
    xcorner[0] = x1[1];
    xcorner[1] = x0[0];
    gengrid.SetBC(gmesh, x0, xcorner, bc1);
    gengrid.SetBC(gmesh, xcorner, x1, bc1);
    xcorner[0] = x0[0];
    xcorner[1] = x1[1];
    gengrid.SetBC(gmesh, x1, xcorner, bc1);
    gengrid.SetBC(gmesh, xcorner, x0, bc1);
    
    int64_t nel = gmesh->NElements();
    
    coarseindices.resize(nel);
    int64_t elcount = 0;
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (gel->HasSubElement() ||  gel->Dimension() != dimension) {
            continue;
        }
        coarseindices[elcount] = el;
        elcount++;
    }
    coarseindices.resize(elcount);
    
    TPZCheckGeom geom(gmesh);
    geom.UniformRefine(ndiv);

    return gmesh;
}

bool IdentifyingFaces(TPZCompMesh *cmesh,TPZStack<TPZCompElSide> &Faces, TPZStack<TPZCompElSide> &AnotherSideFaces) {
    
    if(!cmesh) return false;
    cmesh->Reference()->ResetReference();
    cmesh->LoadReferences();
    
    // To know the dimension of the computational elements over we search
    int ModelDimension = cmesh->Dimension();
    int MaxIndex = MaxCompElementsIndex(cmesh);
    
    TPZAdmChunkVector<TPZCompEl *> elvec = cmesh->ElementVec();
    long i, nel = elvec.NElements();
    // To check if the Face was identified already
    TPZFMatrix<int> FoundedFaces(MaxIndex+1,27,0);
    
    /** Finding faces at mesh: Face is boundary whether it has no neighboard, Face is inner if it has (only one) neighboard. */
    for (i = 0L; i<nel; i++) {
        TPZCompEl *el = (TPZCompEl *)elvec[i];
        // elements with model dimension
        if (!el || el->Dimension() != ModelDimension) continue;
        
        int j, nsides = el->Reference()->NSides();
        for(j=0;j<nsides;j++) {
            TPZCompElSide celside(el,j);
            TPZStack<TPZCompElSide> neigh;
            TPZCompElSide el_neigh;
            // Only over the sides of the codimension 1, if it is not founded
            if(celside.Reference().Dimension() != ModelDimension-1 || FoundedFaces(el->Index(),j))
                continue;
            celside.EqualLevelElementList(neigh, 1, 1);
            if(!neigh.NElements()) {
                el_neigh = celside.LowerLevelElementList(1);
                // existing or not el_neigh we can to register a face
                Faces.push_back(celside);
                AnotherSideFaces.push_back(el_neigh);
                FoundedFaces(el->Index(),j) = 1;
                if(el_neigh.Element()) {
                    FoundedFaces(el_neigh.Element()->Index(),el_neigh.Side()) = 1;
                }
            }
            else if(neigh.NElements() == 1) {
                Faces.push_back(celside);
                AnotherSideFaces.push_back(neigh[0]);
                FoundedFaces(el->Index(),j) = 1;
                FoundedFaces(neigh[0].Element()->Index(),neigh[0].Side()) = 1;
            }
            
        }
    }
    return true;
}

/** To compute Cmin and Cmax of the elliptic equation based on tensor K */
bool ComputeCMinAndCMaxFromTensorK(void (*fp)(TPZVec<REAL> &loc,TPZFMatrix<REAL> &K),TPZCompMesh* cmesh,REAL &Cmin, REAL &Cmax) {
    TPZVec<REAL> pt(3,0.);
    REAL norm, prod;
    int dim = cmesh->Dimension();
    long nel = cmesh->NElements();
    TPZFMatrix<REAL> TensorK(dim,dim,0.);
    TPZFMatrix<REAL> Qsi(dim,1,0.), Psi(dim,1,0.);
    Cmin = Cmax = 1.;
    int i, j;
    
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = cmesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZMaterial * material = cel->Material();
        
        if (!material) continue;
        int dimcel = cel->Dimension();
        TPZAutoPointer<TPZIntPoints> intrule = ((TPZInterpolationSpace*)cel)->GetIntegrationRule().Clone();
        
        TPZManVector<REAL,3> intpoint(dimcel);
        TPZManVector<REAL,3> pt(3,0.);
        REAL weight;
        
        int nintpoints = intrule->NPoints();
        
        for(int nint = 0; nint < nintpoints; nint++) {
            norm = 0.; Qsi.Zero(); Psi.Zero(); prod = 0.;
            intrule->Point(nint,intpoint,weight);
            
            TPZGeoEl * ref = cel->Reference();
            if(!ref) continue;
            ref->X(intpoint, pt);
            
            fp(pt,TensorK);
            //            TensorKFunction(pt,TensorK);
            for(i=0; i<dim; i++) {
                norm += pt[i]*pt[i];
                for(j=0;j<dim;j++)
                    Qsi(i,0) += TensorK(i,j)*pt[j];
                prod += pt[i]*Qsi(i,0);
            }
            if(IsZero(norm)) continue;
            prod /= norm;
            
            Cmin = (Cmin < prod) ? Cmin : prod;
            Cmax = (Cmax < prod) ? prod : Cmax;
            //            if(fp) {
            // if exist function to calculate Tensor K
        }
        
    }
    return true;
    for(long i=0;i < cmesh->Reference()->NodeVec().NElements(); i++) {
        cmesh->Reference()->NodeVec()[i].GetCoordinates(pt);
        norm = 0.;
        for(j=0;j<dim;j++) {
            norm += pt[j]*pt[j];
        }
        if(dim==1) prod = pt[0]*pt[0]*TensorK(0,0);
        else if(dim==2) prod = pt[0]*pt[0]*TensorK(0,0)+pt[0]*pt[1]*(TensorK(0,1)+TensorK(1,0))+pt[1]*pt[1]*TensorK(1,1);
        else return false;
        if(!IsZero(norm)) {
            prod /= norm;
            Cmin = (Cmin < prod) ? Cmin : prod;
            Cmax = (Cmax < prod) ? prod : Cmax;
        }
    }
    return true;
}

void TensorKFunction(TPZVec<REAL> &x,TPZFMatrix<REAL> &K) {
    REAL alpha = 1.;
    int p = K.Rows();
    K.Zero();
    for(int i=0;i<p;i++)
        K(i,i) = alpha;
}
void Dirichlet(const TPZVec<REAL> &loc, TPZVec<STATE> &result)
{
    TPZFMatrix<STATE> du(3,1);
    SolProblema_ParedesThesis(loc,result,du);
}

/**
 * Consider faces, side with codimension 1 on the boundary of the domain and no boundary side with codimension 1 to higher level elements.
 */
bool ComputePressureJumpOnFaces(TPZCompMesh *cmesh,int matid,STATE &Error,STATE &ErrorNi) {
    
    if(!cmesh) return false;
    int ModelDimension = cmesh->Dimension();
    
    TPZStack<TPZCompElSide> Faces;
    TPZStack<TPZCompElSide> AnotherSideFaces;
    IdentifyingFaces(cmesh,Faces,AnotherSideFaces);
    
    //    TPZAdmChunkVector<TPZCompEl *> elvec = cmesh->ElementVec();
    long i, nfaces = Faces.NElements();   //elvec.NElements();
    
    // Identifying material and variable as pressure
    TPZMaterial *mat = cmesh->FindMaterial(matid);
    int varpress = mat->VariableIndex("Pressure");
    if(varpress < 0) return false;
    
    int dimvar = mat->NSolutionVariables(varpress);
    if(dimvar < 1) return false;
    TPZVec<STATE> sol(dimvar,0.);
    TPZVec<STATE> solneigh(dimvar,0.);
    
    // Initializing
    REAL volEl = 0.;
    Error = 0.;
    ErrorNi = 0.;
    STATE val;
    
    // Determining max index of the comp element in mesh
    long maxelindex = MaxCompElementsIndex(cmesh);
    TPZVec<REAL> PressureJump(maxelindex+1,0.);
    
    ComputeCMinAndCMaxFromTensorK(&TensorKFunction,cmesh,Cmin,Cmax);
    
    /** Computing error for all elements with same dimension of the model */
    for (i = 0L; i<nfaces; i++) {
        TPZCompElSide celside = Faces.Pop();
        TPZCompElSide neighcelside = AnotherSideFaces.Pop();
        
        // element with higher level
        TPZCompEl *el = (TPZCompEl *)celside.Element();
        // elements with model dimension
        if (!el || el->Dimension() != ModelDimension) continue;
        TPZGeoEl *gel = el->Reference();
        if(!gel) DebugStop();
        //       int j = celside.Side();
        
        // Working on the face
        //Computing pressure on center of the side
        TPZManVector<REAL,3> pt(3,0.);
        TPZManVector<REAL,3> pt_el(3,0.);
        volEl = celside.Reference().Area();
        // pt - point on face with codimension 1
        celside.Reference().CenterPoint(pt);
        TPZTransform<> tr;
        TPZGeoElSide geosideh(gel,gel->NSides()-1);
        tr = celside.Reference().SideToSideTransform(geosideh);
        // pt_el - point on faces with dimension ModelDimension
        tr.Apply(pt, pt_el);
        // Solution over computational side element+
        el->Solution(pt_el,varpress,sol);
        STATE Err;
        
        // If faces is boundary we need identify if it is dirichlet or Neumann
        if(!neighcelside.Element()) {
            int matidbc = gel->MaterialId();
            if(matidbc > 0) continue;
            TPZMaterial *mat = cmesh->FindMaterial(matid);
            TPZBndCond *bc = dynamic_cast<TPZBndCond *>(mat);
            if(!bc->Type()) {
                TPZVec<STATE> pressbc;
                if(!bc->HasBCForcingFunction())
                    pressbc[0] = bc->Val1()(0,0);
                else
                    Dirichlet(pt,pressbc);
                // Computing (pressure - g).n
                Err = fabs(sol[0]-pressbc[0]);                  // Don't need fabs but only to match formule
                PressureJump[el->Index()] += Cl*Cl*Cmin*(Err*Err);
                Err = volEl*Err;
                Error += Err*Err;
            }
            continue;
        }
        
        //if side has codimension 1 and it is inner then compute the jump pressure
        TPZManVector<REAL,3> pt_n(3,0.);
        TPZManVector<REAL,3> pt_el_n(3,0.);
        
        // working on faces from neighboard element with commom face
        // TPZCompElSide celside_n = AnotherSideFaces.Pop();
        TPZGeoElSide gelside_n = neighcelside.Reference();
        gelside_n.CenterPoint(pt_n);
        TPZGeoElSide gelsideh_n(gelside_n.Element(),gelside_n.Element()->NSides()-1);
        tr = gelside_n.SideToSideTransform(gelsideh_n);
        tr.Apply(pt_n,pt_el_n);
        neighcelside.Element()->Solution(pt_el_n,varpress,solneigh);
        Err = (-0.5)*(sol[0] - solneigh[0]);
        val = Cl*Cl*Cmin*(1./volEl)*(Err*Err);
        PressureJump[el->Index()] += val;
        PressureJump[neighcelside.Element()->Index()] += val;
        Err = volEl*Err;
        Error += Err*Err;
    }
    Error = sqrt(Error);
    
    // Computing ErrorNi
    for(i=0;i<maxelindex+1;i++) {
        val = PressureJump[i];
        ErrorNi += val*val;
    }
    ErrorNi = sqrt(ErrorNi);
    
    return true;
}

