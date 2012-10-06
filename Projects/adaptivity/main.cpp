#include "pzgclonemesh.h"
#include "pzcclonemesh.h"
#include "pzonedref.h"

#include "pzvec.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzvec_extras.h"
#include "pzdebug.h"
#include "pzcheckgeom.h"

#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgeoelside.h"
#include "pzgeoelbc.h"

#include "pzintel.h"
#include "pzcompel.h"

#include "pzmatrix.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "pzadaptmesh.h"

#include "pzmaterial.h"
#include "pzbndcond.h"
#include "pzelasmat.h"
#include "pzplaca.h"
#include "pzmat2dlin.h"
#include "pzmathyperelastic.h"
#include "pzmattest3d.h"
#include "pzmatplaca2.h"

#include "pzgengrid.h"
#include "TPZGenSpecialGrid.h"

#include "TPZRefPatternDataBase.h"

#include "pzfunction.h"
#include "TPZVTKGeoMesh.h"
#include "pzvtkmesh.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.adaptivity"));
static LoggerPtr loggerconv(Logger::getLogger("pz.adaptivity.conv"));
static LoggerPtr loggerpoint(Logger::getLogger("pz.adaptivity.points"));
#endif

#include <time.h>
#include <stdio.h>
#include <fstream>

/** VARIABLES */
/** Printing level */
int gPrintLevel = 0;
/** angle ??? */
static REAL angle = 0.2;
/** Whether the user must to be input data through keyboard */
bool user = true;

/** FUNCTION DEFINITIONS */

/** Forcing function: disp = Forcing1(x)  */
void Forcing1(TPZVec<REAL> &x, TPZVec<REAL> &disp);

/** Functions to create particular mesh */
/** Allows user to choose an option for creating a mesh */
TPZCompMesh *ReadCase(int &nref, int &dim, int &opt,bool usr);
/** Creating selected mesh */
TPZCompMesh *CreateMesh();
TPZCompMesh *CreateSillyMesh();
TPZCompMesh *CreateTriangularMesh();
TPZCompMesh *Create3DMesh();
TPZCompMesh *CreateSimple3DMesh();
TPZCompMesh *Create3DTetraMesh();
TPZCompMesh *Create3DPrismMesh();
TPZCompMesh *CreatePlanMesh();
TPZCompMesh *CreateTestMesh();
TPZCompMesh *CreatePyramTetraMesh();
TPZCompMesh *CreateAleatorioMesh();
TPZCompMesh *Create3DDiscMesh();
TPZCompMesh *Create3DExpMesh();

/** Read mesh from file and create a computational mesh */
TPZCompMesh *ReadKumar(char *filename);
/** Identify maxime level of refinement for all geometric elements referenced by computational elements */
int MaxLevel(TPZCompMesh *mesh);


/** VARIABLES ??? */
static int nstate = 1;

static std::ofstream MALHAG("malhageometrica");
static int mygorder = 1;
int gDebug = 1;

/** Rotation data */
bool rotating = false;
TPZFNMatrix<16,REAL> Rot(4,4,0.),RotInv(4,4,0.);
REAL alfa = M_PI/6.;
REAL transx = 3.;
REAL transy = 0.;

/** FUNCTIONS TO TRANSFORM (ROTATION) MESH */
void InitializeRotation(REAL alfa,REAL transx,REAL transy,TPZMatrix<REAL> &rot,TPZMatrix<REAL> &rotinv);
void TransformMesh(TPZGeoMesh *gmesh,TPZMatrix<REAL> &rot);
void TransformX(TPZVec<REAL> &x,TPZMatrix<REAL> &rot);
void TransformInvX(TPZVec<REAL> &x,TPZMatrix<REAL> &rotinv);

/** Functions to apply boundary conditions */
void NeumannExp(const TPZVec<REAL> &x, TPZVec<REAL> &force);
void Neumann2(const TPZVec<REAL> &x, TPZVec<REAL> &force);
void Neumann3(const TPZVec<REAL> &x, TPZVec<REAL> &force);
void Neumann4(const TPZVec<REAL> &x, TPZVec<REAL> &force);
void Neumann5(const TPZVec<REAL> &x, TPZVec<REAL> &force);

/** Functions with the values of the solutions */
void Exact(const TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<REAL> &dsol);
void ExactSimple3D(const TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<REAL> &dsol);
void Exact3D(const TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<REAL> &dsol);
void Exact3DExp(const TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<REAL> &dsol);
void CompareNeighbours(TPZGeoMesh *mesh);
void BCSolution(const TPZVec<REAL> &x,TPZVec<REAL> &result);
void Solution(const TPZVec<REAL> &x,TPZVec<REAL> &result,TPZFMatrix<REAL> &deriv);
void LoadSolution(TPZFMatrix<REAL> &axes,const TPZVec<REAL> &X,TPZFMatrix<REAL> &u,TPZFMatrix<REAL> &du);

/** MAIN FUNCTION */
int main() {

#ifdef LOG4CXX
    InitializePZLOG();
#endif
	
	// Initializing a ref patterns
	//gRefDBase.InitializeAllUniformRefPatterns();
	//gRefDBase.InitializeRefPatterns();
	gRefDBase.InitializeUniformRefPattern(EOned);
	gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
	gRefDBase.InitializeUniformRefPattern(ETriangle);

	// To visualization of the geometric mesh
	std::ofstream fgeom("GMesh.vtk");
	std::ofstream fgeomfromcomp("GMeshFromComp.vtk");

	// Initializing variables
    int nref = 0;
    int dim = 0;
    int opt = 0;

	// output files
    std::ofstream convergence("conv3d.txt");
    std::ofstream out("output.txt");

	/** Set polynomial order */
    TPZCompEl::SetgOrder(mygorder);
    gDebug = 0;
    
	/** Choosing mesh */
    TPZCompMesh *cmesh = ReadCase(nref,dim,opt,user);
	
	if(rotating) {
		InitializeRotation(alfa,transx,transy,Rot,RotInv);
		TransformMesh(cmesh->Reference(),Rot);
	}

    cmesh->Reference()->SetName("Malha Geometrica original");
    cmesh->SetName("Malha Computacional Original");
    
    cmesh->CleanUpUnconnectedNodes();
		
	// Printing geo mesh to check
	TPZVTKGeoMesh::PrintGMeshVTK(cmesh->Reference(),fgeom);
	TPZVTKGeoMesh::PrintCMeshVTK(cmesh->Reference(),fgeomfromcomp);

	/** Variable names for post processing */
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("POrder");
    scalnames.Push("Error");
    
    if(nstate == 1) {
        scalnames.Push("TrueError");
        scalnames.Push("EffectivityIndex");
        vecnames.Push("state");
    }else if(nstate == 2) {
        scalnames.Push("sig_x");
        scalnames.Push("sig_y");
        scalnames.Push("tau_xy");
        vecnames.Push("state");
    }
    if(dim < 3 && nstate == 2){
        vecnames.Push("displacement");
    }
    
    //Multigrid======================
    //   TPZMGAnalysis mgan (cmesh);
    //   mgan.SetStructuralMatrix(strskyl);
    //   TPZStepSolver *direct = new TPZStepSolver;
    //   direct->SetDirect(ELDLt);
    //   mgan.SetSolver(*direct);
    //   delete direct;
    //   direct = 0;
    //   mgan.Run();
    //   TPZCompMesh *finemesh = cmesh;
    // ===================================
    
	// Solving adaptive process
    {
        int r;
        for(r=0; r<nref; r++) {
            
			// Introduzing exact solution depending on the case
            TPZAnalysis an (cmesh);
            if (opt == 4 || opt == 7 || opt == 8){
                an.SetExact(ExactSimple3D);
            }
            else if(opt==1 || opt==2){
                an.SetExact(Exact);
            }
            else if(opt==12){
                an.SetExact(Exact3D);
            }
            {
                std::stringstream sout;
                int angle = (int) (alfa*180./M_PI + 0.5);
                sout << "hptestAngo" << angle << "." << r << ".vtk";
                an.DefineGraphMesh(dim,scalnames,vecnames,sout.str());
            }
            std::string MeshFileName;
            {
                std::stringstream sout;
                int angle = (int) (alfa*180./M_PI + 0.5);
                sout << "meshAngle" << angle << "." << r << ".vtk";
                MeshFileName = sout.str();
            }
            
            cmesh->SetName("Malha computacional adaptada");
            // Printing geometric and computational mesh
            if (gDebug == 1){
                cmesh->Reference()->Print(std::cout);
                cmesh->Print(std::cout);
            }
            
			// Solve using symmetric matrix then using Cholesky (direct method)
            TPZSkylineStructMatrix strskyl(cmesh);
            an.SetStructuralMatrix(strskyl);
            
            TPZStepSolver<REAL> *direct = new TPZStepSolver<REAL>;
            direct->SetDirect(ECholesky);
            an.SetSolver(*direct);
            delete direct;
            direct = 0;
            
            an.Run();
			
			// Post processing
            an.PostProcess(4,dim);
            {
                std::ofstream out(MeshFileName.c_str());
                cmesh->LoadReferences();
                TPZVTKGeoMesh::PrintCMeshVTK(cmesh->Reference(), out, false);
                
            }
            
            REAL valerror =0.;
            REAL valtruerror=0.;
            TPZVec<REAL> ervec,truervec,effect;
            
            //Multigrid==========================================
            //       finemesh = mgan.UniformlyRefineMesh(finemesh);
            //       mgan.AppendMesh(finemesh);
            //       mgan.Run();
            //       TPZCompMesh *adaptive = mgan.RefinementPattern(finemesh,cmesh,error,truerror,effect);
            //===================================================
            
            TPZAdaptMesh adapt;
            adapt.SetCompMesh (cmesh);
            
            std::cout << "\n\n\n\nEntering Auto Adaptive Methods... step " << r << "\n\n\n\n";
            
            //if(r==4) gPrintLevel = 1;
            time_t sttime;
            time (& sttime);
            TPZCompMesh *adptmesh;
            
            switch (opt){
                case (1) :{
                    adptmesh = adapt.GetAdaptedMesh(valerror,valtruerror,ervec,Exact,truervec,effect,0);
                    break;
                }
                case (2) :{
                    adptmesh = adapt.GetAdaptedMesh(valerror,valtruerror, ervec,Exact,truervec,effect,0);
                    break;
                }
                case (4) :{
                    adptmesh = adapt.GetAdaptedMesh(valerror,valtruerror, ervec,ExactSimple3D,truervec,effect,1);
                    break;
                }
                case (12) : {
                    adptmesh = adapt.GetAdaptedMesh(valerror,valtruerror, ervec,Exact3D,truervec,effect,0);
                    break;
                }
                case (13):{
                    adptmesh = adapt.GetAdaptedMesh(valerror,valtruerror, ervec,Exact3DExp,truervec,effect,1);
                    break;
                }
                default:
                    adptmesh = adapt.GetAdaptedMesh(valerror,valtruerror, ervec,0,truervec,effect,0);
            }
            
            time_t endtime;
            time (& endtime);
            
            int time_elapsed = endtime - sttime;
            std::cout << "\n\n\n\nExiting Auto Adaptive Methods....step " << r
            << "time elapsed " << time_elapsed << "\n\n\n\n";
            
            int prt;
            std::cout << "neq = " << cmesh->NEquations() << " error estimate = " << valerror
            << " true error " << valtruerror <<  " effect " << valerror/valtruerror << std::endl;
            
#ifdef LOG4CXX
            if (loggerconv->isDebugEnabled())
            {
                std::stringstream sout;
                sout << "neq = " << cmesh->NEquations() << " error estimate = " << valerror
                << " true error " << valtruerror <<  " effect " << valerror/valtruerror << std::endl;
                LOGPZ_DEBUG(loggerconv, sout.str())
            }
#endif
            
            convergence  << cmesh->NEquations() << "\t"
            << valerror << "\t"
            << valtruerror << "\t"
            << ( valtruerror / valerror ) <<  "\t"
            << sttime <<std::endl;
            for (prt=0;prt<ervec.NElements();prt++){
                std::cout <<"error " << ervec[prt] << "  truerror = " << truervec[prt] << "  Effect " << effect[prt] << std::endl;
                // convergence << '\t' << ervec[prt] << '\t' << truervec[prt] << "  Effect " << effect[prt] <<  std::endl;
                //  adptmesh->Print(cout);
            }
            
            std::cout.flush();
            cmesh->Reference()->ResetReference();
            cmesh->LoadReferences();
            adapt.DeleteElements(cmesh);
            delete cmesh;
            cmesh = adptmesh;
            
            cmesh->CleanUpUnconnectedNodes();
        }
    }
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        cmesh->ExpandSolution();
        cmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    CompareNeighbours(cmesh->Reference());
    delete cmesh;
    return 0;
}

/** Let be the user choose the type to created mesh, returning this option (opt), dimension of the problem and number of required refinements */
TPZCompMesh *ReadCase(int &nref, int &dim, int &opt,bool user){
    
    std::cout << "**************************************" << std::endl;
    std::cout << "******Auto adaptive test code********" << std::endl;
    std::cout << "**************************************" << std::endl;
    
    std::cout << "Select the analysis type: \n0 - Simple quadrilateral 2D \n1 - L Shape Quadrilateral\n"
    << "2 - Triangular Simples \n3 - Plane mesh (quadrilateral and triangular elements)"
    << "\n4 - 3D Simples \n5 - 3D Canto\n" <<"6 - Tetraedro\n7 - Prisma\n8 - All elements\n9 - All topologies\n10 Aleatorio\n"
    << "11 Pyramid and Tetrahedre\n12Exact 3d Poisson\n"
    << "13 Cube Exp\n";
    opt = 0;
	if(user)
		std::cin >> opt;
	else
		std::cout << "Option " << opt << std::endl;

    
    TPZCompMesh *cmesh;
    
    switch (opt){
        case (0) :{
            cmesh = CreateSillyMesh();
            break;
        }
        case (1) :{
            cmesh = CreateMesh();
            break;
        }
        case (2) :{
            cmesh = CreateTriangularMesh();
            break;
        }
        case (3):{
            cmesh = CreatePlanMesh();
            break;
        }
        case (4) :{
            cmesh = CreateSimple3DMesh();
            break;
        }
        case (5) :{
            cmesh = Create3DMesh();
            break;
        }
        case (6) :{
            cmesh = Create3DTetraMesh();
            break;
        }
        case (7) :{
            cmesh = Create3DPrismMesh();
            break;
        }
        case (8) :{
            cmesh = CreateTestMesh();
            break;
        }
        case (10) :{
            cmesh = CreateAleatorioMesh();
            break;
        }
        case (11) :{
            cmesh = CreatePyramTetraMesh();
            break;
        }
        case (12) :{
            cmesh = Create3DDiscMesh();
            break;
        }
        case (13) :{
            cmesh = Create3DExpMesh();
            break;
        }
            
        default:
            cmesh = CreateMesh();
    }
    
    dim = 2;
    opt > 3 ? dim=3 : dim = 2;
    
    std::cout << "number of refinement steps : ";
    nref = 0;
	if(user)
		std::cin >> nref;
	else
		std::cout << "N Ref " << nref << std::endl;
    
    std::cout << "Maximum p order:    ";
    int p = 1;
	if(user)
		std::cin >> p;
	else
		std::cout << "Order " << p << std::endl;
    std::cout << std::endl;
    
    TPZOneDRef::gMaxP = p;
    
    return cmesh;
}


//*************************************
//************Option 0*****************
//*******L Shape Quadrilateral*********
//*************************************
TPZCompMesh *CreateSillyMesh(){
	    
    //malha quadrada de numrel x numcel
    const	int numrel = 2;
    const	int numcel = 2;

//    int numel = numrel*numcel;
    
    // criar um objeto tipo malha geometrica
    TPZGeoMesh *geomesh = new TPZGeoMesh();
    
    // criar nos
    int i,j;
    for(i=0; i<(numrel+1); i++) {
        for (j=0; j<(numcel+1); j++) {
            int nodind = geomesh->NodeVec().AllocateNewElement();
            TPZVec<REAL> coord(2);
            coord[0] = j;
            coord[1] = i;
            geomesh->NodeVec()[nodind] = TPZGeoNode(i*(numcel+1)+j,coord,*geomesh);
        }
    }
    
    // criação dos elementos
    int elc, elr;
    TPZGeoEl *gel[numrel*numcel];
    TPZVec<int> indices(4);
    for(elr=0; elr<numrel; elr++) {
        for(elc=0; elc<numcel; elc++) {
            indices[0] = (numcel+1)*elr+elc;
            indices[1] = indices[0]+1;
            indices[3] = indices[0]+numcel+1;
            indices[2] = indices[1]+numcel+1;
            // O proprio construtor vai inserir o elemento na malha
            //      gel[elr*numcel+elc] = new TPZGeoElQ2d(elr*numcel+elc,indices,1,*geomesh);
            int index;
            gel[elr*numcel+elc] = geomesh->CreateGeoElement(EQuadrilateral,indices,1,index);
        }
    }

    // Descomentar o trecho abaixo para habilitar a
    // divisão dos elementos geométricos criado
	//Divisão dos elementos
/*	TPZVec<TPZGeoEl *> sub;
	TPZGeoEl *gel = geomesh->ElementVec()[0];
	gel->Divide(sub);
	//     for(int i=0;i<(sub.NElements()-1);i++){
	//		 TPZVec<TPZGeoEl *> subsub;
	//		 sub[i]->Divide(subsub);
	//   } */

    geomesh->BuildConnectivity();
    TPZGeoElBC t3(gel[0],4,-1);
    TPZGeoElBC t4(gel[0],6,-2);
    geomesh->Print(std::cout);
	    
    // Criação da malha computacional
    TPZCompMesh *comp = new TPZCompMesh(geomesh);
    
    // Criar e inserir os materiais na malha
    TPZMaterial * mat = new TPZElasticityMaterial(1,1.e5,0.2,0,0);
    comp->InsertMaterialObject(mat);
    
    TPZMaterial * meumat = mat;
    
    // Condições de contorno
    // Dirichlet
    TPZFMatrix<REAL> val1(3,3,0.),val2(3,1,0.);
    TPZMaterial * bnd = meumat->CreateBC (meumat,-1,0,val1,val2);
    comp->InsertMaterialObject(bnd);
    bnd = meumat->CreateBC (meumat,-2,0,val1,val2);
    
    // bnd->SetForcingFunction(Forcing1);
    // comp->InsertMaterialObject(bnd);
    
    // Neumann
    TPZFMatrix<REAL> val3(3,3,1);
    val2(0,0)=1.;
    bnd = meumat->CreateBC (meumat,-2,1,val1,val2);
    comp->InsertMaterialObject(bnd);
    
    // comp->Print(cout);
    
    // Ajuste da estrutura de dados computacional
    comp->AutoBuild();
    //  comp->Print(cout);
    comp->AdjustBoundaryElements();
    //  comp->Print(cout);
    comp->CleanUpUnconnectedNodes();
    //  comp->Print(cout);
//	geomesh->Print(std::cout);

    /*  //	comp->Print(output); */
    /*  TPZInterpolatedElement *intel = dynamic_cast <TPZInterpolatedElement *> (comp->ElementVec()[1]); */
    /*  TPZVec<int> subelindex; */
    /*  intel->Divide(1,subelindex,1); */
    /*  int isub; */
    /*  int nsides = intel->NConnects(); */
    /*  int porder = intel->PreferredSideOrder(nsides-1); */
    /*   for (isub=0; isub<subelindex.NElements();isub++){ */
    /*     TPZInterpolatedElement *cintel = dynamic_cast<TPZInterpolatedElement *> (comp->ElementVec()[subelindex[isub]]); */
    /*     std::cintel->PRefine(porder+1); */
    /*   } */
    /*comp->ExpandSolution(); */
    
    comp->SetName("Malha Computacional Original");
    //   comp->Print(cout);
    //    std::cout << std::endl << "Number of equations: " << comp->NEquations() << std::endl;
    // std::cout.flush();
    return comp;
}

//*************************************
//************Option 1*****************
//*******L Shape Quadrilateral*********
//*************************************
TPZCompMesh *CreateMesh() {
    REAL co[8][2] = {{0.,0.},{0.,-1.},{1.,-1.},{1.,0.},{1.,1.},{0.,1.},{-1.,1.},{-1.,0.}};
    int indices[3][4] = {{0,1,2,3},{0,3,4,5},{0,5,6,7}};
    TPZGeoEl *elvec[3];
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    int nnode = 8;
    int nod;
    for(nod=0; nod<nnode; nod++) {
        int nodind = gmesh->NodeVec().AllocateNewElement();
        TPZVec<REAL> coord(2);
        coord[0] = co[nod][0];
        coord[1] = co[nod][1];
        gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
    }
    
    int el;
    int nelem = 3;
    for(el=0; el<nelem; el++) {
        TPZVec<int> nodind(4);
        for(nod=0; nod<4; nod++) nodind[nod]=indices[el][nod];
        //    elvec[el] = new TPZGeoElQ2d(el,nodind,1);
        int index;
        elvec[el] = gmesh->CreateGeoElement(EQuadrilateral,nodind,1,index);
    }
    
    gmesh->BuildConnectivity();
    
    TPZVec<TPZGeoEl *> sub;

    // bc -1 -> Dirichlet
    TPZGeoElBC gbc1(elvec[0],4,-1);
    // bc -2 -> Neumann at the bottom y==-1
    TPZGeoElBC gbc2(elvec[0],5,-2);
    // bc -3 -> Neumann at the right x==1
    TPZGeoElBC gbc3(elvec[0],6,-3);
    
    // bc -3 -> Neumann at the right x==1
    TPZGeoElBC gbc4(elvec[1],5,-3);
    
    // bc -4 -> Neumann at the top y==1
    TPZGeoElBC gbc5(elvec[1],6,-4);
    
    // bc -4 -> Neumann at the top y==1
    TPZGeoElBC gbc6(elvec[2],5,-4);
    
    // bc -5 -> Neumann at the left x==-1
    TPZGeoElBC gbc7(elvec[2],6,-5);
    
    // bc -6 -> Homogeneous Neumann
    TPZGeoElBC gbc8(elvec[2],7,-6);
    
    int nel = gmesh->NElements();
    for (int iel=0; iel<nel; iel++) {
        TPZManVector<TPZGeoEl *> subels;
        TPZGeoEl *gel = gmesh->ElementVec()[iel];
        gel->Divide(subels);
    }
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(2);
    
    TPZMaterial * mat;
    if(nstate == 2) {
        mat = new TPZElasticityMaterial(1,2.,0.3,1.,1.);
    } else {
        TPZMat2dLin *mat2d = new TPZMat2dLin(1);
        int ist,jst;
        TPZFMatrix<REAL> xk(nstate,nstate,1.),xc(nstate,nstate,0.),xf(nstate,1,0.);
        for(ist=0; ist<nstate; ist++) {
            if(nstate != 1) xf(ist,0) = 1.;
            for(jst=0; jst<nstate; jst++) {
                if(ist != jst) xk(ist,jst) = 0.;
            }
        }
        mat2d->SetMaterial(xk,xc,xf);
        mat = mat2d;
    }
    TPZFMatrix<REAL> val1(nstate,nstate,0.),val2(nstate,1,0.);
    TPZManVector<TPZMaterial *,6> bc(6);
    bc[0] = mat->CreateBC(mat,-1,0,val1,val2);
    int i;
    if(nstate == 1) {
        for(i=1; i<6; i++) {
            bc[i] = mat->CreateBC(mat,-i-1,1,val1,val2);
        }
        bc[1]->SetForcingFunction(new TPZDummyFunction<STATE>(Neumann2));
        bc[2]->SetForcingFunction(new TPZDummyFunction<STATE>(Neumann3));
        bc[3]->SetForcingFunction(new TPZDummyFunction<STATE>(Neumann4));
        bc[4]->SetForcingFunction(new TPZDummyFunction<STATE>(Neumann5));
    } else {
        for(i=1; i<6; i++) {
            bc[i] = mat->CreateBC(mat,-i-1,0,val1,val2);
        }
    }
    
    cmesh->InsertMaterialObject(mat);
    for(i=0; i<6; i++) cmesh->InsertMaterialObject(bc[i]);
    
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    //  cmesh->ExpandSolution();
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        cmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    return cmesh;
}


//*************************************
//************Option 2*****************
//******Simple Triangular Mesh*********
//*************************************
TPZCompMesh *CreateTriangularMesh(){
    
    const int nelem = 1;
    const int nnode = 3;
    
    REAL co[nnode][3] = {{0.,0.,0.},{1.,0.,0.},{0.,1.,0.}};
    int indices[nelem][nnode] = {{0,1,2}};
    
    TPZGeoEl *elvec[nelem];
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    int nod;
    for(nod=0; nod<nnode; nod++) {
        int nodind = gmesh->NodeVec().AllocateNewElement();
        TPZVec<REAL> coord(3);
        coord[0] = co[nod][0];
        coord[1] = co[nod][1];
        coord[2] = co[nod][2];
        gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
    }
    
    int el;
    for(el=0; el<nelem; el++) {
        TPZVec<int> nodind(3);
        for(nod=0; nod<3; nod++) nodind[nod]=indices[el][nod];
        //    elvec[el] = new TPZGeoElT2d(el,nodind,1);
        int index;
        elvec[el] = gmesh->CreateGeoElement(ETriangle,nodind,1,index);
    }
    
    gmesh->BuildConnectivity();
      TPZStack<TPZGeoEl*> subel;
      elvec[0]->Divide(subel);
    
    
    // bc -1 -> Dirichlet
    TPZGeoElBC gbc1(elvec[0],3,-1);
    
    // bc -2 -> Neumann at the right x==1
    TPZGeoElBC gbc2(elvec[0],5,-2);
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    
    // inserir os materiais
    TPZMaterial * mat = new TPZElasticityMaterial(1,1.e5,0.2,0,0);
    cmesh->InsertMaterialObject(mat);
    
    TPZMaterial * meumat = mat;
    
    // inserir a condicao de contorno
    TPZFMatrix<REAL> val1(3,3,0.),val2(3,1,0.);
    
    TPZMaterial * bnd = meumat->CreateBC (meumat,-1,0,val1,val2);
    cmesh->InsertMaterialObject(bnd);
    
    val2(0,0)=1.;
    bnd = meumat->CreateBC (meumat,-2,1,val1,val2);
    cmesh->InsertMaterialObject(bnd);
    
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
}


//*************************************
//************Option 3*****************
//*****Plane Quad & Triang Mesh********
//*************************************
TPZCompMesh *CreatePlanMesh() {
    
    REAL co[5][2] = {
        {0.,0.},
        {1.,0.},
        {2.,0.},
        {0.,1.},
        {1.,1.}
    };
    
    int indices[2][4] = {
        {0,1,4,3},
        {1,2,4,-1}
    };
    
    TPZGeoEl *elvec[2];
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    int nnode = 5;
    int nelem = 2;
    
    int nod;
    for(nod=0; nod<nnode; nod++) {
        int nodind = gmesh->NodeVec().AllocateNewElement();
        TPZVec<REAL> coord(2);
        coord[0] = co[nod][0];
        coord[1] = co[nod][1];
        gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
    }
    
    int el;
    
    for(el=0; el<nelem; el++) {
        int ncnodes = el > 0 ? 3 : 4;
        TPZVec<int> nodind(ncnodes);
        for(nod=0; nod<ncnodes; nod++) nodind[nod]=indices[el][nod];
        int index;
        if (el == 0){
            //      elvec[el] = new TPZGeoElQ2d(el,nodind,1);
            elvec[el] = gmesh->CreateGeoElement(EQuadrilateral,nodind,1,index);
        }else{
            //      elvec[el] = new TPZGeoElT2d(el,nodind,1);
            elvec[el] = gmesh->CreateGeoElement(ETriangle,nodind,1,index);
        }
    }
    
    gmesh->BuildConnectivity();
    
    // bc -1 -> Dirichlet
    TPZGeoElBC gbc1(elvec[0],4,-1);
    TPZGeoElBC gbc2(elvec[1],3,-2);
    
    // bc -4 -> Neumann at the top y==1
    TPZGeoElBC gbc3(elvec[0],6,-3);
    
    gmesh->SetName ("Original Geometric Mesh");
    gmesh->Print(std::cout);
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetName ("Original Computational Mesh");
    
    // Criar e inserir os materiais na malha
    TPZMaterial * mat = new TPZElasticityMaterial(1,1.e5,0.2,0,0);
    cmesh->InsertMaterialObject(mat);
    
    TPZMaterial * meumat = mat;
    
    // Condições de contorno
    // Dirichlet
    TPZFMatrix<REAL> val1(3,3,0.),val2(3,1,0.);
    TPZMaterial * bnd = meumat->CreateBC (meumat,-1,0,val1,val2);
    cmesh->InsertMaterialObject(bnd);
    
    bnd = meumat->CreateBC (meumat,-2,0,val1,val2);
    cmesh->InsertMaterialObject(bnd);
    
    // Neumann
    val2(0,0) = 1.;
    bnd = meumat->CreateBC (meumat,-3,1,val1,val2);
    cmesh->InsertMaterialObject(bnd);
    
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    //cmesh->ExpandSolution();
    
    //cmesh->Print(cout);
    return cmesh;
}

//*************************************
//************Option 4*****************
//********Simple 3D Cube Mesh**********
//*************************************
TPZCompMesh *CreateSimple3DMesh() {
    
    REAL co[8][3] = {
        {0.,0.,0.},
        {1.,0.,0.},
        {1.,1.,0.},
        {0.,1.,0.},
        {0.,0.,1.},
        {1.,0.,1.},
        {1.,1.,1.},
        {0.,1.,1.}
    };
    int indices[1][8] = {{0,1,2,3,4,5,6,7}};
    
    const int nelem = 1;
    int nnode = 8;
    
    TPZGeoEl *elvec[nelem];
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    int nod;
    for(nod=0; nod<nnode; nod++) {
        int nodind = gmesh->NodeVec().AllocateNewElement();
        TPZVec<REAL> coord(3);
        coord[0] = co[nod][0];
        coord[1] = co[nod][1];
        coord[2] = co[nod][2];
        gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
    }
    
    int el;
    for(el=0; el<nelem; el++) {
        TPZManVector<int> nodind(8);
        for(nod=0; nod<8; nod++) nodind[nod]=indices[el][nod];
        int index;
        elvec[el] = gmesh->CreateGeoElement(ECube,nodind,1,index);
        //    elvec[el] = new TPZGeoElC3d(el,nodind,1);
    }
    
    gmesh->BuildConnectivity();
    
    //TPZVec<TPZGeoEl *> sub;
    //elvec[0]->Divide(sub);
    //   	elvec[1]->Divide(sub);
    //   	elvec[2]->Divide(sub);
    
    
    // bc -1 -> Dirichlet at the bottom face of the cube
    TPZGeoElBC gbc1(elvec[0],20,-1);
    
    // bc -2 -> Neumann at the top face of the cube
    TPZGeoElBC gbc2(elvec[0],25,-2);
    
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    
    TPZMaterial * mat = new TPZMaterialTest3D(1);
    if(nstate == 3) {
        //		mat = new TPZMatHyperElastic(1,2.,400);
        // mat = new TPZMaterialTest3D(1);
        TPZFMatrix<REAL> mp (3,1,0.);
        TPZMaterialTest3D * mataux = dynamic_cast<TPZMaterialTest3D *> (mat);
        TPZMaterialTest3D::geq3=1;
        mataux->SetMaterial(mp);
    } else {
        TPZMat2dLin *mat2d = new TPZMat2dLin(1);
        int ist,jst;
        TPZFMatrix<REAL> xk(nstate,nstate,1.),xc(nstate,nstate,0.),xf(nstate,1,0.);
        for(ist=0; ist<nstate; ist++) {
            if(nstate != 1) xf(ist,0) = 1.;
            for(jst=0; jst<nstate; jst++) {
                if(ist != jst) xk(ist,jst) = 0.;
            }
        }
        mat2d->SetMaterial(xk,xc,xf);
        // mat = mat2d;
    }
    TPZFMatrix<REAL> val1(1,1,0.),val2(1,1,0.);
    TPZMaterial * bc[2];
    bc[0] = mat->CreateBC(mat,-1,0,val1,val2);
    int i;
    
    val2(0,0)=1.;
    bc[1] = mat->CreateBC(mat,-2,1,val1,val2);
    
    cmesh->InsertMaterialObject(mat);
    for(i=0; i<2; i++) cmesh->InsertMaterialObject(bc[i]);
    
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    //cmesh->ExpandSolution();
    
    return cmesh;
}

//*************************************
//************Option 5*****************
//********Corner 3D Cube Mesh**********
//*************************************
TPZCompMesh *Create3DMesh() {
    
    REAL co[26][3] = {
        {0.,0.,0.},{1.,0.,0.},{2.,0.,0.},
        {0.,1.,0.},{1.,1.,0.},{2.,1.,0.},
        {0.,2.,0.},{1.,2.,0.},{2.,2.,0.},
        {0.,0.,1.},{1.,0.,1.},{2.,0.,1.},
        {0.,1.,1.},{1.,1.,1.},{2.,1.,1.},
        {0.,2.,1.},{1.,2.,1.},{2.,2.,1.},
        {0.,0.,2.},{1.,0.,2.},{2.,0.,2.},
        {0.,1.,2.},{1.,1.,2.},{2.,1.,2.},
        {0.,2.,2.},{1.,2.,2.}
    };
    
    int indices[7][8] = {
        {0,1,4,3,9,10,13,12},
        {1,2,5,4,10,11,14,13},
        {3,4,7,6,12,13,16,15},
        {4,5,8,7,13,14,17,16},
        {9,10,13,12,18,19,22,21},
        {10,11,14,13,19,20,23,22},
        {12,13,16,15,21,22,25,24}
    };
    
    const int nelem = 7;
    int nnode = 26;
    
    TPZGeoEl *elvec[nelem];
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    int nod;
    for(nod=0; nod<nnode; nod++) {
        int nodind = gmesh->NodeVec().AllocateNewElement();
        TPZVec<REAL> coord(3);
        coord[0] = co[nod][0];
        coord[1] = co[nod][1];
        coord[2] = co[nod][2];
        gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
    }
    
    int el;
    for(el=0; el<nelem; el++) {
        TPZVec<int> nodind(8);
        for(nod=0; nod<8; nod++) nodind[nod]=indices[el][nod];
        int index;
        elvec[el] = gmesh->CreateGeoElement(ECube,nodind,1,index);
        //    elvec[el] = new TPZGeoElC3d(el,nodind,1);
    }
    
    gmesh->BuildConnectivity();
    
    //  	TPZVec<TPZGeoEl *> sub;
    //   	elvec[0]->Divide(sub);
    //   	elvec[1]->Divide(sub);
    //   	elvec[2]->Divide(sub);
    
    
    // bc -1 -> Dirichlet
    TPZGeoElBC gbc1(elvec[0],21,-1);
    
    TPZGeoElBC gbc2(elvec[0],22,-1);
    
    TPZGeoElBC gbc3(elvec[0],25,-1);
    
    // bc -3 -> Neumann at the right x==1
    TPZGeoElBC gbc4(elvec[3],25,-2);
    
    // bc -4 -> Neumann at the top y==1
    TPZGeoElBC gbc5(elvec[5],23,-3);
    
    // bc -4 -> Neumann at the top y==1
    TPZGeoElBC gbc6(elvec[6],22,-4);
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    
    TPZMaterial * mat;
    if(nstate == 2) {
        mat = new TPZMatHyperElastic(1,2.,400);
    } else {
        TPZMat2dLin *mat2d = new TPZMat2dLin(1);
        int ist,jst;
        TPZFMatrix<REAL> xk(nstate,nstate,1.),xc(nstate,nstate,0.),xf(nstate,1,0.);
        for(ist=0; ist<nstate; ist++) {
            if(nstate != 1) xf(ist,0) = 1.;
            for(jst=0; jst<nstate; jst++) {
                if(ist != jst) xk(ist,jst) = 0.;
            }
        }
        mat2d->SetMaterial(xk,xc,xf);
        mat = mat2d;
    }
    TPZFMatrix<REAL> val1(3,3,0.),val2(3,1,0.);
    TPZMaterial * bc[4];
    bc[0] = mat->CreateBC(mat,-1,0,val1,val2);
    int i;
    val2(2,0)=-1.;
    bc[1] = mat->CreateBC(mat,-2,1,val1,val2);
    val2(2,0)=0.;
    val2(1,0)=-1.;
    bc[2] = mat->CreateBC(mat,-3,1,val1,val2);
    val2(1,0)=0.;
    val2(0,0)=-1.;
    bc[3] = mat->CreateBC(mat,-4,1,val1,val2);
    
    cmesh->InsertMaterialObject(mat);
    for(i=0; i<4; i++) cmesh->InsertMaterialObject(bc[i]);
    
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    //cmesh->ExpandSolution();
    
    return cmesh;
}

//*************************************
//************Option 6*****************
//*********3D Tetrahedra Mesh**********
//*************************************
TPZCompMesh *Create3DTetraMesh() {
    
    REAL co[4][3] = {
        {0.,0.,0.},
        {1.,0.,0.},
        {0.,1.,0.},
        {0.,0.,1.}
    };
    
    int indices[1][4] = {{0,1,2,3}};
    const int nelem = 1;
    int nnode = 4;
    
    TPZGeoEl *elvec[nelem];
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    int nod;
    for(nod=0; nod<nnode; nod++) {
        int nodind = gmesh->NodeVec().AllocateNewElement();
        TPZVec<REAL> coord(3);
        coord[0] = co[nod][0];
        coord[1] = co[nod][1];
        coord[2] = co[nod][2];
        gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
    }
    
    int el;
    for(el=0; el<nelem; el++) {
        TPZVec<int> nodind(4);
        for(nod=0; nod<4; nod++) nodind[nod]=indices[el][nod];
        int index;
        elvec[el] = gmesh->CreateGeoElement(ETetraedro,nodind,1,index);
        
        //    elvec[el] = new TPZGeoElT3d(el,nodind,1);
    }
    
    
    //TPZStack<TPZGeoEl*> subel;
    //elvec[0]->Divide(subel);
    
    
    // bc -1 -> Dirichlet
    TPZGeoElBC gbc1(elvec[0],10,-1);
    
    // bc -2 -> Neumann at the right x==1
    TPZGeoElBC gbc2(elvec[0],13,-2);
    
    gmesh->BuildConnectivity();
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    
    TPZMaterial * mat;
    if(nstate == 3) {
        //		mat = new TPZMatHyperElastic(1,2.,400);
        mat = new TPZMaterialTest3D(1);
        TPZFMatrix<REAL> mp (3,1,1.);
        
        TPZMaterialTest3D * mataux = dynamic_cast<TPZMaterialTest3D *> (mat);
        TPZMaterialTest3D::geq3=1;
        mataux->SetMaterial(mp);
    } else {
        TPZMat2dLin *mat2d = new TPZMat2dLin(1);
        int ist,jst;
        TPZFMatrix<REAL> xk(nstate,nstate,1.),xc(nstate,nstate,0.),xf(nstate,1,0.);
        for(ist=0; ist<nstate; ist++) {
            if(nstate != 1) xf(ist,0) = 1.;
            for(jst=0; jst<nstate; jst++) {
                if(ist != jst) xk(ist,jst) = 0.;
            }
        }
        mat2d->SetMaterial(xk,xc,xf);
        mat = mat2d;
    }
    TPZFMatrix<REAL> val1(3,3,0.),val2(3,1,0.);
    TPZMaterial * bc[2];
    bc[0] = mat->CreateBC(mat,-1,0,val1,val2);
    int i;
    
    val2(2,0)=-1.;
    bc[1] = mat->CreateBC(mat,-2,1,val1,val2);
    
    cmesh->InsertMaterialObject(mat);
    for(i=0; i<2; i++) cmesh->InsertMaterialObject(bc[i]);
    
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    //cmesh->ExpandSolution();
    
    return cmesh;
}

//*************************************
//************Option 7*****************
//**********3D Prism Mesh**************
//*************************************
TPZCompMesh *Create3DPrismMesh() {
    
    REAL co[6][3] = {
        {0.,0.,0.},
        {1.,0.,0.},
        {0.,1.,0.},
        {0.,0.,1.},
        {1.,0.,1.},
        {0.,1.,1.}
    };
    
    int indices[1][6] = {{0,1,2,3,4,5}};
    
    const int nelem = 1;
    int nnode = 6;
    
    TPZGeoEl *elvec[nelem];
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    int nod;
    for(nod=0; nod<nnode; nod++) {
        int nodind = gmesh->NodeVec().AllocateNewElement();
        TPZVec<REAL> coord(3);
        coord[0] = co[nod][0];
        coord[1] = co[nod][1];
        coord[2] = co[nod][2];
        gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
    }
    
    int el;
    for(el=0; el<nelem; el++) {
        TPZVec<int> nodind(6);
        for(nod=0; nod<6; nod++) nodind[nod]=indices[el][nod];
        int index;
        elvec[el] = gmesh->CreateGeoElement(EPrisma,nodind,1,index);
        //    elvec[el] = new TPZGeoElPr3d(el,nodind,1);
    }
    
    //  TPZStack<TPZGeoEl*> subel;
    //  elvec[0]->Divide(subel);
    
    
    // bc -1 -> Neumann
    TPZGeoElBC gbc1(elvec[0],15,-1);
    
    // bc -2 -> Neumann at the right x==1
    TPZGeoElBC gbc2(elvec[0],19,-2);
    
    // bc -2 -> Dirichlet at point 0
    TPZGeoElBC gbc3(elvec[0],0,-3);
    
    
    gmesh->BuildConnectivity();
    //ofstream MALHAG("malhageometrica");
    gmesh->Print(MALHAG);
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    
    TPZMaterial * mat;
    if(nstate == 3) {
        //		mat = new TPZMatHyperElastic(1,2.,400);
        mat = new TPZMaterialTest3D(1);
        TPZFMatrix<REAL> mp (1,1,0.);
        
        TPZMaterialTest3D * mataux = dynamic_cast<TPZMaterialTest3D *> (mat);
        TPZMaterialTest3D::geq3=1;
        mataux->SetMaterial(mp);
    } else {
        TPZMat2dLin *mat2d = new TPZMat2dLin(1);
        int ist,jst;
        TPZFMatrix<REAL> xk(nstate,nstate,1.),xc(nstate,nstate,0.),xf(nstate,1,0.);
        for(ist=0; ist<nstate; ist++) {
            if(nstate != 1) xf(ist,0) = 1.;
            for(jst=0; jst<nstate; jst++) {
                if(ist != jst) xk(ist,jst) = 0.;
            }
        }
        mat2d->SetMaterial(xk,xc,xf);
        mat = mat2d;
    }
    TPZFMatrix<REAL> val1(1,1,0.),val2(1,1,0.);
    TPZMaterial * bc[3];
    bc[0] = mat->CreateBC(mat,-3,0,val1,val2);
    int i;
    val2(0,0)=-1.;
    bc[1] = mat->CreateBC(mat,-2,1,val1,val2);
    val2(0,0)=1.;
    bc[2] = mat->CreateBC(mat,-1,1,val1,val2);
    
    cmesh->InsertMaterialObject(mat);
    for(i=0; i<3; i++) cmesh->InsertMaterialObject(bc[i]);
    
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    cmesh->Print(std::cout);
    
    /*   int o; */
    /*   for (o=0;o<10;o++){ */
    /*     TPZIntPrism3D prismrule(2*o+2,2*o+2); */
    /*     TPZVec<int> ord (3,2*o+2); */
    /*     prismrule.SetOrder(ord); */
    /*     int np = prismrule.NPoints(); */
    /*     int p; */
    /*     std::cout << std::endl << std::endl <<"Ordem o = " << 2*o+2 << std::endl; */
    /*     for (p=0;p<np;p++){ */
    /*       TPZVec<REAL> loc(3,0.); */
    /*       REAL weight = -1.; */
    /*       prismrule.Point(p,loc,weight); */
    /*       std::cout << "Point " << p << "  (x,y,z) = " << loc[0] << " , "   */
    /* 	   << loc[1] << " , "  << loc[2] << "  weight = "  << weight << std::endl; */
    /*     } */
    /*   } */
    
    return cmesh;
}

//*************************************
//************Option 8*****************
//*****All element types Mesh**********
//*************************************
TPZCompMesh * CreateTestMesh() {
    
    REAL nodeco[12][3] = {
        {0.,0.,0.},
        {1.,0.,0.},
        {2.,0.,0.},
        {0.,1.,0.},
        {1.,1.,0.},
        {2.,1.,0.},
        {0.,0.,1.},
        {1.,0.,1.},
        {2.,0.,1.},
        {0.,1.,1.},
        {1.,1.,1.},
        {2.,1.,1.}
    };
    
    int nodind[7][8] = {
        {0,1,4,3,6,7,10,9},
        {2,4,10,8,5},
        {8,10,11,5},
        {2,4,1,8,10,7},
        {0,1},
        {0,1,7,6},
        {1,2,7}
    };
    
    int numnos[7] = {8,5,4,6,2,4,3};
    
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    int noind[12];
    int no;
    for(no=0; no<12; no++) {
        noind[no] = gmesh->NodeVec().AllocateNewElement();
        TPZVec<REAL> coord(3);
        coord[0] = nodeco[no][0];
        coord[1] = nodeco[no][1];
        coord[2] = nodeco[no][2];
        gmesh->NodeVec()[noind[no]].Initialize(coord,*gmesh);
    }
    int matid = 1;
    TPZVec<int> nodeindex;
    int nel;
    TPZVec<TPZGeoEl *> gelvec;
    gelvec.Resize(4);
    for(nel=0; nel<4; nel++) {
        int in;
        nodeindex.Resize(numnos[nel]);
        for(in=0; in<numnos[nel]; in++) {
            nodeindex[in] = nodind[nel][in];
        }
        int index;
        switch(nel) {
            case 0:
                //      elvec[el] = gmesh->CreateGeoElement(ECube,nodeindex,1,index);
                //      gelvec[nel]=new TPZGeoElC3d(nodeindex,matid);
                break;
            case 1:
                gelvec[nel] = gmesh->CreateGeoElement(EPiramide,nodeindex,matid,index);
                //       gelvec[nel]=new TPZGeoElPi3d(nodeindex,matid);
                break;
            case 2:
                gelvec[nel] = gmesh->CreateGeoElement(ETetraedro,nodeindex,matid,index);
                //       gelvec[nel]=new TPZGeoElT3d(nodeindex,matid);
                break;
            case 3:
                //       gelvec[nel]=new TPZGeoElPr3d(nodeindex,matid);
                //      gelvec[nel] = gmesh->CreateGeoElement(EPrisma,nodeindex,matid,index);
                break;
            case 4:
                //      gelvec[nel]=new TPZGeoEl1d(nodeindex,2);
                break;
            case 5:
                //      gelvec[nel]=new TPZGeoElQ2d(nodeindex,3);
                break;
            case 6:
                //      gelvec[nel]=new TPZGeoElT2d(nodeindex,3);
                break;
            default:
                break;
        }
    }
    gmesh->BuildConnectivity();
    
    //TPZVec<TPZGeoEl *> sub;
    //elvec[0]->Divide(sub);
    //   	elvec[1]->Divide(sub);
    //   	elvec[2]->Divide(sub);
    
    
    // bc -1 -> Dirichlet
    //  TPZGeoElBC gbc1(gelvec[0],20,-1);
    TPZGeoElBC gbc11(gelvec[1],14,-1);
    //  TPZGeoElBC gbc12(gelvec[3],15,-1);
    
    
    
    // bc -2 -> Neumann at the right x==1
    //  TPZGeoElBC gbc2(gelvec[0],25,-2);
    //  TPZGeoElBC gbc21(gelvec[3],19,-2);
    TPZGeoElBC gbc22(gelvec[2],10,-2);
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    
    TPZMaterial * mat;
    if(nstate == 3) {
        mat = new TPZMaterialTest3D(1);
        TPZFMatrix<REAL> mp (3,1,0.);
        TPZMaterialTest3D * mataux = dynamic_cast<TPZMaterialTest3D *> (mat);
        TPZMaterialTest3D::geq3=1;
        mataux->SetMaterial(mp);
    } else {
        TPZMat2dLin *mat2d = new TPZMat2dLin(1);
        int ist,jst;
        TPZFMatrix<REAL> xk(nstate,nstate,1.),xc(nstate,nstate,0.),xf(nstate,1,0.);
        for(ist=0; ist<nstate; ist++) {
            if(nstate != 1) xf(ist,0) = 1.;
            for(jst=0; jst<nstate; jst++) {
                if(ist != jst) xk(ist,jst) = 0.;
            }
        }
        mat2d->SetMaterial(xk,xc,xf);
        mat = mat2d;
    }
    
    TPZFMatrix<REAL> val1(3,3,0.),val2(3,1,0.);
    TPZMaterial * bc[2];
    
    bc[0] = mat->CreateBC(mat,-1,0,val1,val2);
    val2(0,0) = 1.;
    bc[1] = mat->CreateBC(mat,-2,1,val1,val2);
    cmesh->InsertMaterialObject(mat);
    
    int i;
    for(i=0; i<2; i++) cmesh->InsertMaterialObject(bc[i]);
    
    gmesh->Print(std::cout);
    
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    gmesh->Print(std::cout);
    return cmesh;
}

/** Forcing function: disp = Forcing1(x)  */
void Forcing1(TPZVec<REAL> &x, TPZVec<REAL> &disp){
    TPZManVector<REAL,3> x2(x);
    TransformInvX(x2,RotInv);
    TPZFNMatrix<3,REAL> grad(3,1,0.),grad2(3,1,0.);
	
    grad(0,0) = -(x2[1]-0.5)*sin(angle)+(x2[0]-0.5)*cos(angle)-(x2[0]-0.5);
    grad(1,0) = (x2[1]-0.5)*cos(angle)+(x2[0]-0.5)*sin(angle)-(x2[1]-0.5);
    grad(2,0) = 0.;
    Rot.Multiply(grad, grad2);
    disp[0] = grad2(0,0);
    disp[1] = grad2(1,0);
    disp[2] = grad2(2,0);
}

/** Exact solutions to calculate the rate of convergence */

static REAL onethird = 0.33333333333333333;
static REAL PI = 3.141592654;

void Exact(const TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<REAL> &dsol) {
    TPZManVector<REAL,3> x2(x);
    TransformInvX(x2,RotInv);
  	REAL r = sqrt(x2[0]*x2[0]+x2[1]*x2[1]);
  	REAL theta = atan2(x2[1],x2[0]);
#ifdef LOG4CXX
    if (loggerpoint->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "Point " << x2;
        LOGPZ_DEBUG(loggerpoint, sout.str())
    }
#endif
  	REAL rexp = pow(r,onethird);
  	sol[0] = rexp*sin(onethird*(theta+PI/2));
    TPZFNMatrix<3,REAL> grad(4,1,0.),grad2(4,1,0.);
  	grad(0,0) = onethird*sin(onethird*(PI/2.-2.*theta))/(rexp*rexp);
  	grad(1,0) = onethird*cos(onethird*(PI/2.-2.*theta))/(rexp*rexp);
    Rot.Multiply(grad, grad2);
    dsol(0,0) = grad2(0,0);
    dsol(1,0) = grad2(1,0);
}

void Exact3D(const TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<REAL> &dsol) {
    TPZManVector<REAL,3> x2(x);
    TransformInvX(x2,RotInv);
  	REAL r = sqrt(x2[0]*x2[0]+x2[1]*x2[1]+x2[2]*x2[2]);
  	REAL theta = atan2(x2[1],x2[0]);
  	REAL rexp = pow(r,onethird);
  	sol[0] = rexp*sin(onethird*(theta+PI/2));
    TPZFNMatrix<3,REAL> grad(4,1,0.),grad2(4,1,0.);
  	grad(0,0) = onethird*sin(onethird*(PI/2.-2.*theta))/(rexp*rexp);
  	grad(1,0) = onethird*cos(onethird*(PI/2.-2.*theta))/(rexp*rexp);
    Rot.Multiply(grad, grad2);
    dsol(0,0) = grad2(0,0);
    dsol(1,0) = grad2(1,0);
}

void ExactSimple3D(const TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<REAL> &dsol) {
    TPZManVector<REAL,3> x2(x);
    TransformInvX(x2,RotInv);
	sol[0] = x2[2];
    TPZFNMatrix<3,REAL> grad(4,1,0.),grad2(4,1,0.);
  	grad(0,0) = 0.;
  	grad(1,0) = 0.;
    grad(2,0) = 1.;
    Rot.Multiply(grad, grad2);
    dsol(0,0) = grad2(0,0);
    dsol(1,0) = grad2(1,0);
    dsol(2,0) = grad2(2,0);
}

void Exact3DExp(const TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<REAL> &dsol) {
    TPZManVector<REAL,3> x2(x);
    TransformInvX(x2,RotInv);
	REAL one = (REAL)1.;

    sol[0] = pow(exp(one),1/(0.1 + pow(-0.5 + x2[0],2) + pow(-0.5 + x2[1],2)))*(1 - x2[0])*x2[0]*(1 - x2[1])*x2[1]*x2[2];
    TPZFNMatrix<3,REAL> grad(4,1,0.),grad2(4,1,0.);
    grad(0,0) = pow(exp(one),1/(0.1 + pow(-0.5 + x2[0],2) + pow(-0.5 + x2[1],2)))*(1 - x2[0])*(1 - x2[1])*x2[1]*x2[2] -
    pow(exp(one),1/(0.1 + pow(-0.5 + x2[0],2) + pow(-0.5 + x2[1],2)))*x2[0]*(1 - x2[1])*x2[1]*x2[2] -
    (2*pow(exp(one),1/(0.1 + pow(-0.5 + x2[0],2) + pow(-0.5 + x2[1],2)))*(1 - x2[0])*(-0.5 + x2[0])*x2[0]*
     (1 - x2[1])*x2[1]*x2[2])/pow(0.1 + pow(-0.5 + x2[0],2) + pow(-0.5 + x2[1],2),2);
    grad(1,0) = pow(exp(one),1/(0.1 + pow(-0.5 + x2[0],2) + pow(-0.5 + x2[1],2)))*(1 - x2[0])*x2[0]*(1 - x2[1])*x2[2] -
    pow(exp(one),1/(0.1 + pow(-0.5 + x2[0],2) + pow(-0.5 + x2[1],2)))*(1 - x2[0])*x2[0]*x2[1]*x2[2] -
    (2*pow(exp(one),1/(0.1 + pow(-0.5 + x2[0],2) + pow(-0.5 + x2[1],2)))*(1 - x2[0])*x2[0]*(1 - x2[1])*
     (-0.5 + x2[1])*x2[1]*x2[2])/pow(0.1 + pow(-0.5 + x2[0],2) + pow(-0.5 + x2[1],2),2);
    grad(2,0) = pow(exp(one),1/(0.1 + pow(-0.5 + x2[0],2) + pow(-0.5 + x2[1],2)))*(1 - x2[0])*x2[0]*(1 - x2[1])*x2[1];
    Rot.Multiply(grad, grad2);
    dsol(0,0) = grad2(0,0);
    dsol(1,0) = grad2(1,0);
    dsol(2,0) = grad2(2,0);
}


/** Boundary conditions */
void NeumannExp(const TPZVec<REAL> &x, TPZVec<REAL> &f) {
    TPZManVector<REAL,3> x2(x);
    TransformInvX(x2,RotInv);
	REAL one = (REAL) 1.0;
    f[0] = pow(exp(one),1/(0.1 + pow(-0.5 + x2[0],2) + pow(-0.5 + x2[1],2)))*(1 - x2[0])*x2[0]*(1 - x2[1])*x2[1];
}

void Neumann2(const TPZVec<REAL> &x, TPZVec<REAL> &f) {
    TPZManVector<REAL,3> x2(x);
    TransformInvX(x2,RotInv);
  	REAL r = sqrt(x2[0]*x2[0]+x2[1]*x2[1]);
  	REAL theta = atan2(x2[1],x2[0]);
  	REAL rexp = pow(r,onethird);
  	f[0] = -onethird*cos(onethird*(PI/2.-2.*theta))/(rexp*rexp);
}

void Neumann3(const TPZVec<REAL> &x, TPZVec<REAL> &f) {
    TPZManVector<REAL,3> x2(x);
    TransformInvX(x2,RotInv);
  	REAL r = sqrt(x2[0]*x2[0]+x2[1]*x2[1]);
  	REAL theta = atan2(x2[1],x2[0]);
  	REAL rexp = pow(r,onethird);
  	f[0] = onethird*sin(onethird*(PI/2.-2.*theta))/(rexp*rexp);
}

void Neumann4(const TPZVec<REAL> &x, TPZVec<REAL> &f) {
    TPZManVector<REAL,3> x2(x);
    TransformInvX(x2,RotInv);
    /*
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "Point " << x2;
        LOGPZ_DEBUG(loggerpoint, sout.str())
    }
#endif
     */
  	REAL r = sqrt(x2[0]*x2[0]+x2[1]*x2[1]);
  	REAL theta = atan2(x2[1],x2[0]);
  	REAL rexp = pow(r,onethird);
  	f[0] = onethird*cos(onethird*(PI/2.-2.*theta))/(rexp*rexp);
}

void Neumann5(const TPZVec<REAL> &x, TPZVec<REAL> &f) {
    TPZManVector<REAL,3> x2(x);
    TransformInvX(x2,RotInv);
    /*
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "Point " << x2;
        LOGPZ_DEBUG(loggerpoint, sout.str())
    }
#endif
     */
  	REAL r = sqrt(x2[0]*x2[0]+x2[1]*x2[1]);
  	REAL theta = atan2(x2[1],x2[0]);
  	REAL rexp = pow(r,onethird);
  	f[0] = -onethird*sin(onethird*(PI/2.-2.*theta))/(rexp*rexp);
}

/** Identify maxime level of refinement for all geometric elements referenced by computational elements */
int MaxLevel(TPZCompMesh *mesh) {
  	int nel = mesh->NElements();
  	int el;
  	int level = 0;
  	for(el=0; el<nel; el++) {
        TPZCompEl *cel = mesh->ElementVec()[el];
        if(!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if(!gel) continue;
        int gellev = gel->Level();
        level = (level <gellev) ? gellev : level;
  	}
  	return level;
}

/** Read mesh from file and create a computational mesh */
TPZCompMesh *ReadKumar(char *filename) {
    
  	int nnodes,nelem,nmat,nbcd,nbc;
  	TPZGeoMesh *gmesh = new TPZGeoMesh();
    std::ifstream input(filename);
  	if(!input) {
        std::cout << "file could not be opened " << filename << std::endl;
        return 0;
  	}
  	char buf[256];
  	input.getline(buf,256);
  	input >> nnodes >> nelem >> nmat >> nbcd >> nbc;
  	gmesh->NodeVec().Resize(nnodes);
  	input.getline(buf,256);
  	while(buf[0] != '#') input.getline(buf,256);
  	int id;
  	TPZVec<REAL> coord(2);
  	int nod;
  	for(nod=0; nod< nnodes; nod++) {
        input >> id >> coord[0] >> coord[1];
        gmesh->NodeVec()[id] = TPZGeoNode(id,coord,*gmesh);
  	}
    
  	input.getline(buf,256);
  	while(buf[0] != '#') input.getline(buf,256);
  	char c = input.peek();
  	while(c == '#') {
        input.getline(buf,256);
        c = input.peek();
  	}
    
  	int nel;
  	TPZVec<int> elvertices(4);
  	int elnodes[9],matindex;
  	for(nel=0; nel<nelem; nel++) {
        for(nod=0; nod<9; nod++) input>>elnodes[nod];
        input >> matindex;
        for(nod=0; nod<4; nod++) elvertices[nod] = elnodes[nod];
		//    		new TPZGeoElQ2d(nel,elvertices,matindex);
		int index;
		gmesh->CreateGeoElement(EQuadrilateral,elvertices,matindex,index);
  	}
    
  	gmesh->BuildConnectivity();
	input.getline(buf,256);
	char *compare = strstr(buf,"# BC records");
    
  	while(!compare) {
        input.getline(buf,256);
        compare = strstr(buf,"# BC records");
  	}
    
  	int elnum, side, bcnum;
  	int bc;
  	for(bc=0; bc<nbc; bc++) {
        input >> elnum >> side >> bcnum;
        TPZGeoElBC(gmesh->ElementVec()[elnum-1],side+4-1,-bcnum-1);
  	}
    
 	while(!strstr(buf,"real value)")) input.getline(buf,256);
  	REAL e1111,e1122,e2222,e1212;
  	input >> buf >> e1111 >> buf >> e1122 >> buf >> e2222 >> buf >> e1212;
    
  	REAL E,nu;
  	nu = e1122/e1111;
  	E = e1212*(1+nu);
    
  	TPZMaterial * mat = new TPZElasticityMaterial(3,E,nu,0.,0.);
  	TPZFMatrix<REAL> val1(2,2,0.),val2(2,1,0.);
  	TPZMaterial * bc1 = mat->CreateBC(mat,-1,0,val1,val2);
  	val2(1,0) = -1.;
  	TPZMaterial * bc2 = mat->CreateBC(mat,-2,1,val1,val2);
  	TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  	cmesh->InsertMaterialObject(mat);
  	cmesh->InsertMaterialObject(bc1);
  	cmesh->InsertMaterialObject(bc2);
    
  	cmesh->AutoBuild();
  	cmesh->AdjustBoundaryElements();
  	cmesh->CleanUpUnconnectedNodes();
  	//  TPZCompMesh *fine = TPZMGAnalysis::UniformlyRefineMesh(cmesh);
  	//  delete cmesh;
  	return cmesh;
}

//*************************************
//************Option 10*****************
//**********Aleatorio Mesh**************
//*************************************
TPZCompMesh *CreateAleatorioMesh() {
    
    REAL co[6][3] = {
        {0.,0.,0.},
        {1.,0.,0.},
        {0.,1.,0.},
        {0.,0.,1.},
        {1.,0.,1.},
        {0.,1.,1.}
    };
    
    int indices[1][6] = {{0,1,2,3,4,5}};
    
    const int nelem = 1;
    int nnode = 6;
    
    TPZGeoEl *elvec[nelem];
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    int nod;
    for(nod=0; nod<nnode; nod++) {
        int nodind = gmesh->NodeVec().AllocateNewElement();
        TPZVec<REAL> coord(3);
        coord[0] = co[nod][0];
        coord[1] = co[nod][1];
        coord[2] = co[nod][2];
        gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
    }
    
    int el;
    for(el=0; el<nelem; el++) {
        TPZVec<int> nodind(6);
        for(nod=0; nod<6; nod++) nodind[nod]=indices[el][nod];
        int index;
        elvec[el] = gmesh->CreateGeoElement(EPrisma,nodind,1,index);
        //    elvec[el] = new TPZGeoElPr3d(el,nodind,1);
    }
    
    //  TPZStack<TPZGeoEl*> subel;
    //  elvec[0]->Divide(subel);
    
    gmesh->BuildConnectivity();
    
    // bc -1 -> Dirichlet
    TPZGeoElBC gbc1(elvec[0],15,-1);
    
    // bc -2 -> Neumann at the right x==1
    TPZGeoElBC gbc2(elvec[0],19,-2);
    
    gmesh->BuildConnectivity();
    //ofstream MALHAG("malhageometrica");
    gmesh->Print(MALHAG);
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    
    TPZMaterial * mat;
    if(nstate == 3) {
        //		mat = new TPZMatHyperElastic(1,2.,400);
        mat = new TPZMaterialTest3D(1);
        TPZFMatrix<REAL> mp (3,1,0.);
        
        TPZMaterialTest3D * mataux = dynamic_cast<TPZMaterialTest3D *> (mat);
        TPZMaterialTest3D::geq3=1;
        mataux->SetMaterial(mp);
    } else {
        TPZMat2dLin *mat2d = new TPZMat2dLin(1);
        int ist,jst;
        TPZFMatrix<REAL> xk(nstate,nstate,1.),xc(nstate,nstate,0.),xf(nstate,1,0.);
        for(ist=0; ist<nstate; ist++) {
            if(nstate != 1) xf(ist,0) = 1.;
            for(jst=0; jst<nstate; jst++) {
                if(ist != jst) xk(ist,jst) = 0.;
            }
        }
        mat2d->SetMaterial(xk,xc,xf);
        mat = mat2d;
    }
    TPZFMatrix<REAL> val1(3,3,0.),val2(3,1,0.);
    TPZMaterial * bc[2];
    bc[0] = mat->CreateBC(mat,-1,0,val1,val2);
    int i;
    val2(0,0)=-1.;
    bc[1] = mat->CreateBC(mat,-2,0,val1,val2);
    
    cmesh->InsertMaterialObject(mat);
    for(i=0; i<2; i++) cmesh->InsertMaterialObject(bc[i]);
    
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    
    TPZVec <int> subelvec;
    cmesh->ElementVec()[0]->Divide(0,subelvec,1);
    
    cmesh->ElementVec()[subelvec [7]]->Divide(subelvec [7],subelvec,1);
    TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cmesh->ElementVec()[subelvec [7]]);
    intel->PRefine(2);
    cmesh->ElementVec()[subelvec [7]]->Divide(subelvec [7],subelvec,1);
    TPZInterpolatedElement *intel1 = dynamic_cast<TPZInterpolatedElement *> (cmesh->ElementVec()[subelvec [7]]);
    intel1->PRefine(3);
    
    //  cmesh->ElementVec()[subelvec [7]]->PRefine(4);
    //  cmesh->ElementVec()[0]->Divide(,subelvec,1);
    
    //  cmesh->Print(cout);
    return cmesh;
}


//**************************************
//************Option 11*****************
//******Pyramid and Tetrahedre**********
//**************************************
TPZCompMesh *CreatePyramTetraMesh() {
    
    REAL co[8][3] = {
        {0.,0.,0.},
        {1.,0.,0.},
        {1.,1.,0.},
        {0.,1.,0.},
        {0.,0.,1.},
        {1.,0.,1.},
        {1.,1.,1.},
        {0.,1.,1.}
        
    };
    
    int noel [4] = {5,5,4,4};
    
    int indices[4][5] = {
        {0,1,2,3,4},
        {2,3,7,6,4},
        {1,2,6,4},
        {1,6,5,4},
    };
    
    const int nelem = 4;
    int nnode = 8;
    
    TPZGeoEl *elvec[nelem];
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    int nod;
    for(nod=0; nod<nnode; nod++) {
        int nodind = gmesh->NodeVec().AllocateNewElement();
        TPZVec<REAL> coord(3);
        coord[0] = co[nod][0];
        coord[1] = co[nod][1];
        coord[2] = co[nod][2];
        gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
    }
    
    int el;
    for(el=0; el<nelem; el++) {
        TPZVec<int> nodind(6);
        for(nod=0; nod<noel[el]; nod++) nodind[nod]=indices[el][nod];
        int index;
        //     if (noel[el] == 5) elvec[el] = new TPZGeoElPr3d(el,nodind,1);
        //     if (noel[el] == 4) elvec[el] = new TPZGeoElT3d(el,nodind,1);
        if (noel[el] == 5){
            nodind.Resize(5);
            elvec[el] = elvec[el] = gmesh->CreateGeoElement(EPrisma,nodind,1,index);
        }
        if (noel[el] == 4) {
            nodind.Resize(4);
            elvec[el] = elvec[el] = gmesh->CreateGeoElement(ETetraedro,nodind,1,index);
        }
    }
    
    //  TPZStack<TPZGeoEl*> subel;
    //  elvec[0]->Divide(subel);
    
    //  TPZGeoElBC gbc;
    
    // bc -1 -> Neumann
    TPZGeoElBC gbc1(elvec[0],13,-1);
    
    // bc -2 -> Neumann
    TPZGeoElBC gbc2(elvec[1],15,-2);
    
    // bc -3 -> Neumann
    TPZGeoElBC gbc3(elvec[3],12,-3);
    
    // bc -3 -> Dirichlet
    TPZGeoElBC gbc4(elvec[0],0,-4);
    
    
    //  gmesh->BuildConnectivity2();
    gmesh->BuildConnectivity();
    std::ofstream MALHAG("malhageometrica");
    gmesh->Print(MALHAG);
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    
    TPZMaterial * mat;
    //  if(nstate == 3) {
    //		mat = new TPZMatHyperElastic(1,2.,400);
    mat = new TPZMaterialTest3D(1);
    TPZFMatrix<REAL> mp (3,1,0.);
    
    TPZMaterialTest3D * mataux = dynamic_cast<TPZMaterialTest3D *> (mat);
    TPZMaterialTest3D::geq3=1;
    mataux->SetMaterial(mp);
    //   } else {
    //     TPZMat2dLin *mat2d = new TPZMat2dLin(1);
    //     int ist,jst;
    //     TPZFMatrix<REAL> xk(nstate,nstate,1.),xc(nstate,nstate,0.),xf(nstate,1,0.);
    //     for(ist=0; ist<nstate; ist++) {
    //       if(nstate != 1) xf(ist,0) = 1.;
    //       for(jst=0; jst<nstate; jst++) {
    // 	if(ist != jst) xk(ist,jst) = 0.;
    //       }
    //     }
    //     mat2d->SetMaterial(xk,xc,xf);
    //     mat = mat2d;
    //   }
    TPZFMatrix<REAL> val1(3,3,0.),val2(3,1,0.);
    TPZMaterial * bc[2];
    bc[0] = mat->CreateBC(mat,-1,0,val1,val2);
    int i;
    val2(0,0)=-1.;
    bc[1] = mat->CreateBC(mat,-2,0,val1,val2);
    
    cmesh->InsertMaterialObject(mat);
    for(i=0; i<2; i++) cmesh->InsertMaterialObject(bc[i]);
    
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    
    TPZVec <int> subelvec;
    cmesh->ElementVec()[0]->Divide(0,subelvec,1);
    
    //   int  pord = 3;
    //   cmesh->ElementVec()[subelvec [7]]->Divide(subelvec [7],subelvec,1);
    //   TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cmesh->ElementVec()[subelvec [7]]);
    //   intel->PRefine(2);
    //   cmesh->ElementVec()[subelvec [7]]->Divide(subelvec [7],subelvec,1);
    //   TPZInterpolatedElement *intel1 = dynamic_cast<TPZInterpolatedElement *> (cmesh->ElementVec()[subelvec [7]]);
    //   intel1->PRefine(3);
    
    //  cmesh->ElementVec()[subelvec [7]]->PRefine(4);
    //  cmesh->ElementVec()[0]->Divide(,subelvec,1);
    
    
    //  cmesh->Print(cout);
    return cmesh;
}


//*************************************
//************Option 12*****************
//*******3D Discontinuous Mesh**********
//*************************************
TPZCompMesh *Create3DDiscMesh() {
    
    REAL co[26][3] = {
        {0.,0.,0.},{0.5,0.,0.},{1,0.,0.},{0.,0.5,0.},{0.5,0.5,0.},{1,0.5,0.},{0.,1.,0.},{0.5,1.,0.},{1,1.,0.},
        {0.,0.,0.5},{0.5,0.,0.5},{1,0.,0.5},{0.,0.5,0.5},{0.5,0.5,0.5},{1,0.5,0.5},{0.,1.,0.5},{0.5,1.,0.5},{1,1.,0.5},
        {0.,0.,1.},{0.5,0.,1.},{1,0.,1.},{0.,0.5,1.},{0.5,0.5,1.},{1,0.5,1.},{0.,1.,1.},{0.5,1.,1.}
    };
    
    int indices[7][8] = {
        {0,1,4,3,9,10,13,12},
        {1,2,5,4,10,11,14,13},
        {3,4,7,6,12,13,16,15},
        {4,5,8,7,13,14,17,16},
        {9,10,13,12,18,18,22,21},
        {10,11,14,13,19,20,23,22},
        {12,13,16,15,21,22,25,24}
    };
    
    const int nelem = 7;
    int nnode = 26;
    const int nodeperel = 8;
    
    TPZGeoEl *elvec[nelem];
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    int nod;
    for(nod=0; nod<nnode; nod++) {
        int nodind = gmesh->NodeVec().AllocateNewElement();
        TPZVec<REAL> coord(3);
        coord[0] = co[nod][0];
        coord[1] = co[nod][1];
        coord[2] = co[nod][2];
        gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
    }
    
    int el;
    for(el=0; el<nelem; el++) {
        TPZVec<int> nodind(nodeperel);
        for(nod=0; nod<nodeperel; nod++) nodind[nod]=indices[el][nod];
        int index;
        elvec[el] = gmesh->CreateGeoElement(ECube,nodind,1,index);
        //    elvec[el] = new TPZGeoElPr3d(el,nodind,1);
    }
    
    
    
    //Condicoes de Neumann
    // bc -1 -> Face inferior
    TPZGeoElBC gbc1(elvec[0],20,-1);
    TPZGeoElBC gbc2(elvec[1],20,-1);
    TPZGeoElBC gbc3(elvec[2],20,-1);
    TPZGeoElBC gbc4(elvec[3],20,-1);
    
    // bc -2 -> Face lateral esquerda
    TPZGeoElBC gbc5(elvec[0],24,-2);
    TPZGeoElBC gbc6(elvec[2],24,-2);
    TPZGeoElBC gbc7(elvec[4],24,-2);
    TPZGeoElBC gbc8(elvec[6],24,-2);
    
    // bc -3 -> Face frontal
    TPZGeoElBC gbc9(elvec[0],21,-3);
    TPZGeoElBC gbc10(elvec[1],21,-3);
    TPZGeoElBC gbc11(elvec[4],21,-3);
    TPZGeoElBC gbc12(elvec[5],21,-3);
    
    // bc -3 -> Face lateral direita
    TPZGeoElBC gbc13(elvec[1],22,-4);
    TPZGeoElBC gbc14(elvec[3],22,-4);
    TPZGeoElBC gbc15(elvec[5],22,-4);
    
    // bc -3 -> Face posterior
    TPZGeoElBC gbc16(elvec[2],23,-5);
    TPZGeoElBC gbc17(elvec[3],23,-5);
    TPZGeoElBC gbc18(elvec[6],23,-5);
    
    //Condicoes Dirichlet
    TPZGeoElBC gbc19(elvec[3],25,-6);
    TPZGeoElBC gbc20(elvec[5],23,-6);
    TPZGeoElBC gbc21(elvec[6],22,-6);
    
    gmesh->BuildConnectivity();
    //ofstream MALHAG("malhageometrica");
    //gmesh->Print(MALHAG);
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    
    TPZMaterial *mat;
    if(nstate == 3) {
        mat = new TPZMaterialTest3D(1);
        TPZFMatrix<REAL> mp (3,1,0.);
        TPZMaterialTest3D * mataux = dynamic_cast<TPZMaterialTest3D *> (mat);
        TPZMaterialTest3D::geq3=1;
        mataux->SetMaterial(mp);
    } else {
        //   TPZPoison3D *mat3d = new TPZPoison3D();
        int ist,jst;
        TPZFMatrix<REAL> xk(nstate,nstate,1.),xc(nstate,nstate,0.),xf(nstate,1,0.);
        for(ist=0; ist<nstate; ist++) {
            if(nstate != 1) xf(ist,0) = 1.;
            for(jst=0; jst<nstate; jst++) {
                if(ist != jst) xk(ist,jst) = 0.;
            }
        }
        //    mat2d->SetMaterial(xk,xc,xf);
        //    mat = mat2d;
    }
    
    
    //CreateBC
    //....
    
    
    cmesh->InsertMaterialObject(mat);
    // for(int i=0; i<2; i++) cmesh->InsertMaterialObject(bc[i]);
    
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    
    TPZVec <int> subelvec;
    cmesh->ElementVec()[0]->Divide(0,subelvec,1);
    
    cmesh->ElementVec()[subelvec [7]]->Divide(subelvec [7],subelvec,1);
    TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cmesh->ElementVec()[subelvec [7]]);
    intel->PRefine(2);
    cmesh->ElementVec()[subelvec [7]]->Divide(subelvec [7],subelvec,1);
    TPZInterpolatedElement *intel1 = dynamic_cast<TPZInterpolatedElement *> (cmesh->ElementVec()[subelvec [7]]);
    intel1->PRefine(3);
    
    //  cmesh->ElementVec()[subelvec [7]]->PRefine(4);
    //  cmesh->ElementVec()[0]->Divide(,subelvec,1);
        
    //  cmesh->Print(cout);
    return cmesh;
}


void CompareNeighbours(TPZGeoMesh *mesh) {
    
    TPZAdmChunkVector<TPZGeoEl *> &geovec = mesh->ElementVec();
    int nel = geovec.NElements();
    int iel;
    for(iel=0; iel<nel; iel++) {
        TPZGeoEl *gel = geovec[iel];
        if(!gel) continue;
        int nsides = gel->NSides();
        int is;
        for(is=0; is<nsides; is++) {
            TPZStack<TPZGeoElSide> st1,st2;
            TPZGeoElSide gelside = TPZGeoElSide(gel,is);
            gelside.AllNeighbours(st1);
            gelside.ComputeNeighbours(st2);
            Sort<TPZGeoElSide>(st1);
            Sort<TPZGeoElSide>(st2);
            int nlist1 = st1.NElements();
            int nlist2 = st2.NElements();
            if(nlist1 != nlist2) {
                std::cout << "AllNeighbours is different form ComputeNeighbours\n";
                continue;
            }
            int il;
            for(il=0; il<nlist1; il++) {
                if(st1[il] != st2[il]) {
                    std::cout << "Different neighbours\n";
                }
            }
        }
    }
}

/** Placa do Cedric **/
TPZCompMesh *PlateMesh () {
	TPZCompMesh *compmesh = 0;
    /*
	 static double Teste2Nos[5][3]  = { {0.,0.,0.},{2.,0.,0.},{2.,2.,0.},{0.,2.,0.},{1.,1.,0.}};
	 static int    Teste2Elems[4][8] = {{0,1,4},{1,2,4},{2,3,4},{3,0,4}};

	 //inserindo a ordem de interpola¢ão dos elementos e do espa¢o
     int ord;
     std::cout << "Entre ordem 1,2,3,4,5 : -> 1";
     //cin >> ord;
     ord = 1; std::cout << std::endl;
     TPZCompEl::gOrder = ord;
     
     //malha geometrica
     TPZGeoMesh *geomesh = new TPZGeoMesh;
     
     geomesh->Reference();
     
     //cria nós geométricos
     geomesh->NodeVec().Resize(5);
     TPZVec<REAL> coord(3);
     int i;
     for(i=0;i<num;i++) {
     coord[0]= Teste2Nos[i][0];
     coord[1]= Teste2Nos[i][1];
     coord[2]= Teste2Nos[i][2];
     geomesh.NodeVec()[i].Initialize(coord,geomesh);
     }
     
     //elementos geométricos
     TPZVec<int> index(3);
     int i,j,indice;
     for(i=0;i<4;i++) {
     for(j=0;j<3;j++) {
     index[j] =  Teste2Elems[i][j];
     }
     geomesh.CreateGeoElement(ETriangle,index,1,indice);
     }
     
     //montagem de conectividades entre elementos
     geomesh->BuildConnectivity();
     
     //cria malha computacional
     compmesh = new TPZCompMesh(geomesh);
     
     //cria material do problema
     TPZMaterial *placa = LerMaterial("placa.in",*compmesh);
     
     //transferindo para o material uma carga conhecida
     //placa->SetForcingFunction(LoadSolution);
     (dynamic_cast<TPZPlaca *>(placa))->SetExactFunction(LoadSolution);
     
     //cria condi¢ões de contorno
     if(1) CriaCondContTeste4(*geomesh);
     
     //cria elementos computacionais
     compmesh->AutoBuild();
     cmesh->AdjustBoundaryElements();
     cmesh->CleanUpUnconnectedNodes();
     
     
     //calcula o erro energia da solu¢ão
     SolutionError(compmesh);
     
     //arquivo de saida de dados
     ofstream data("mesh.out");
     geomesh->Print(data);
     compmesh->Print(data);
     data.flush();
     data.close();
     erros.close();
     
     delete compmesh;
     delete geomesh;
	 */
     return compmesh;
}

TPZMaterial *LerMaterial(char *filename, TPZCompMesh &cmesh) {
    std::ifstream input(filename);
    TPZFMatrix<REAL> naxes(3,3);
    REAL ni1,ni2,h,E1,E2,G12,G13,G23,f;
    REAL n00,n01,n02,n10,n11,n12,n20,n21,n22;
    TPZVec<REAL> xf(6);
    int matindex;
    input >> matindex;
    input >> f   >>  h  >>
    E1  >> E2  >>
    G12 >> G13 >> G23 >>
    ni1 >> ni2;
    input >> n00 >> n01 >> n02;
    input >> n10 >> n11 >> n12;
    input >> n20 >> n21 >> n22;
    input >> xf[0] >> xf[1] >> xf[2] >> xf[3] >> xf[4] >> xf[5];
    naxes(0,0) =  n00;    naxes(0,1) =  n01;    naxes(0,2) =  n02;
    naxes(1,0) =  n10;    naxes(1,1) =  n11;    naxes(1,2) =  n12;
    naxes(2,0) =  n20;    naxes(2,1) =  n21;    naxes(2,2) =  n22;
    TPZMaterial *placa = new TPZPlaca(matindex,h,f,E1,E2,ni1,ni2,G12,G13,G23,naxes,xf);
    cmesh.InsertMaterialObject(placa);
    return placa;
}

void CriaCondContTeste4(TPZGeoMesh &gmesh){
    int indicematerial = 1;
    TPZMaterial * placa = gmesh.Reference()->FindMaterial(indicematerial);
    if(!placa){
        std::cout << "main::CriaCond material nao existe, CC nao criadas\n";
        std::cout << "\t\tindice material pedido : " << indicematerial << std::endl;
        return;
    }
    TPZGeoEl *elg0 = gmesh.FindElement(0);
    TPZGeoEl *elg1 = gmesh.FindElement(1);
    TPZGeoEl *elg2 = gmesh.FindElement(2);
    TPZGeoEl *elg3 = gmesh.FindElement(3);
    //malha computacional
    TPZCompMesh *cmesh = gmesh.Reference();
    //BIG number
    //REAL big = 1.e12;
    //valor das CC
    TPZFMatrix<REAL> val1(6,6,0.),val2(6,1,0.);
    TPZGeoElBC(elg0,3,-1);
    TPZGeoElBC(elg1,3,-1);
    TPZGeoElBC(elg2,3,-1);
    TPZGeoElBC(elg3,3,-1);
    TPZMaterial * bc = placa->CreateBC(placa,-1,0,val1,val2);
    bc->SetForcingFunction(new TPZDummyFunction<STATE>(BCSolution));
    cmesh->InsertMaterialObject(bc);
}

void BCSolution(const TPZVec<REAL> &x,TPZVec<REAL> &result){
    TPZFMatrix<REAL> deriv(2,6);
    Solution(x,result,deriv);
}

void Solution(const TPZVec<REAL> &x,TPZVec<REAL> &result,TPZFMatrix<REAL> &deriv){
    TPZFMatrix<REAL> eixos(3,3,0.);
    eixos(0,0) = 1.0;
    eixos(1,1) = 1.0;
    eixos(2,2) = 1.0;
    TPZFMatrix<REAL> u(6,1);
    LoadSolution(eixos,x,u,deriv);
    for(int i=0;i<6;i++) result[i] = u(i,0);
}

void LoadSolution(TPZFMatrix<REAL> &axes,const TPZVec<REAL> &X,TPZFMatrix<REAL> &u,TPZFMatrix<REAL> &du){
    
    REAL k,eps,div;
    static int key = 1;
    static REAL val;
    
    k = 0.0;//0.00001;
    eps = 0.4;//-1.0 é interessante
    div = 100.0;
    
    if(key){
        std::cout << "main::LoadSolution valor de div = ";
        std::cin >> val;
        std::cout << std::endl;
        key = 0;
    }
    div = val;
    
    REAL x = X[0],y=X[1];
    
    REAL x2 = (x-1.0)*(x-1.0);
    REAL y2 = (y-1.0)*(y-1.0);
    REAL eps2 = eps*eps;
    
    REAL expx = exp(-x2/eps2);
    REAL expy = exp(-y2/eps2);
    
    REAL exdivmk = expx/div+k;
    REAL eydivmk = expy/div+k;
    
    REAL exy = expx*expy;
    
    REAL exmeydivmk = 2.*expx*eydivmk / div / eps2;
    REAL eymexdivmk = 2.*expy*exdivmk / div / eps2;
    //u exact
    u(0,0) =  0.0;//u
    u(1,0) =  0.0;//v
    u(2,0) =  (expx/div+k)*(expy/div+k);//w
    u(3,0) = -( 2.*expy*(expx/div+k)*(y-1.) ) / div / eps2;//øx
    u(4,0) =  ( 2.*expx*(expy/div+k)*(x-1.) ) / div / eps2;//øy
    u(5,0) =  0.0;//øz
    
    //du exact
    TPZFMatrix<REAL> dur(2,6);
    dur(0,0)  =  0.0;//du/dx
    dur(1,0)  =  0.0;//du/dy
    
    dur(0,1)  =  0.0;//dv/dx
    dur(1,1)  =  0.0;//dv/dx
    
    dur(0,2)  =  -exmeydivmk*(x-1.);//dw/dx
    dur(1,2)  =  -eymexdivmk*(y-1.);//dw/dy
    
    dur(0,3)  =  4.*exy*(x-1.)*(y-1.) / (div*div) / (eps2*eps2);//døx/dx
    dur(1,3)  = -eymexdivmk + 2.*eymexdivmk*(y-1.)*(y-1.) / eps2;//døx/dy
    
    dur(0,4)  =  exmeydivmk - 2.*exmeydivmk*(x-1.)*(x-1.) / eps2;//døy/dx
    dur(1,4)  = -4.*exy*(x-1.)*(y-1.) / (div*div) / (eps2*eps2);//døy/dy
    
    dur(0,5)  = 0.0;//døz/dx
    dur(1,5)  = 0.0;//døz/dy
    
    du(0,0) = axes(0,0)*dur(0,0) + axes(1,0)*dur(1,0);
    du(1,0) = axes(0,1)*dur(0,0) + axes(1,1)*dur(1,0);
    
    du(0,1) = axes(0,0)*dur(0,1) + axes(1,0)*dur(1,1);
    du(1,1) = axes(0,1)*dur(0,1) + axes(1,1)*dur(1,1);
    
    du(0,2) = axes(0,0)*dur(0,2) + axes(1,0)*dur(1,2);
    du(1,2) = axes(0,1)*dur(0,2) + axes(1,1)*dur(1,2);
    
    du(0,3) = axes(0,0)*dur(0,3) + axes(1,0)*dur(1,3);
    du(1,3) = axes(0,1)*dur(0,3) + axes(1,1)*dur(1,3);
    
    du(0,4) = axes(0,0)*dur(0,4) + axes(1,0)*dur(1,4);
    du(1,4) = axes(0,1)*dur(0,4) + axes(1,1)*dur(1,4);
    
    du(0,5) = axes(0,0)*dur(0,5) + axes(1,0)*dur(1,5);
    du(1,5) = axes(0,1)*dur(0,5) + axes(1,1)*dur(1,5);
    
}


//*************************************
//************Option 13****************
//**********Exp 3D Cube Mesh***********
//*************************************
TPZCompMesh *Create3DExpMesh() {
    
    REAL co[8][3] = {
        {0.,0.,0.},
        {1.,0.,0.},
        {1.,1.,0.},
        {0.,1.,0.},
        {0.,0.,1.},
        {1.,0.,1.},
        {1.,1.,1.},
        {0.,1.,1.}
    };
    int indices[1][8] = {{0,1,2,3,4,5,6,7}};
    
    const int nelem = 1;
    int nnode = 8;
    
    TPZGeoEl *elvec[nelem];
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    int nod;
    for(nod=0; nod<nnode; nod++) {
        int nodind = gmesh->NodeVec().AllocateNewElement();
        TPZVec<REAL> coord(3);
        coord[0] = co[nod][0];
        coord[1] = co[nod][1];
        coord[2] = co[nod][2];
        gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
    }
    
    int el;
    for(el=0; el<nelem; el++) {
        TPZVec<int> nodind(8);
        for(nod=0; nod<8; nod++) nodind[nod]=indices[el][nod];
        int index;
        elvec[el] = gmesh->CreateGeoElement(ECube,nodind,1,index);
        //    elvec[el] = new TPZGeoElC3d(el,nodind,1);
    }
    
    gmesh->BuildConnectivity();
    
    //TPZVec<TPZGeoEl *> sub;
    //elvec[0]->Divide(sub);
    //   	elvec[1]->Divide(sub);
    //   	elvec[2]->Divide(sub);
    
    //  TPZGeoElBC gbc;
    
    // bc -1 -> Dirichlet at the bottom face of the cube
    TPZGeoElBC bc1(elvec[0],20,-1);
    
    // bc -2 -> Neumann at the top face of the cube
    TPZGeoElBC bc2(elvec[0],25,-2);
    
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    
    TPZMaterial * mat;
    //  if(nstate == 3) {
    //		mat = new TPZMatHyperElastic(1,2.,400);
    mat = new TPZMaterialTest3D(1);
    TPZFMatrix<REAL> mp (3,1,0.);
    TPZMaterialTest3D * mataux = dynamic_cast<TPZMaterialTest3D *> (mat);
    TPZMaterialTest3D::geq3=1;
    mataux->SetMaterial(mp);
    //    mat->SetForcingFunction(NeumannExp);
    //   } else {
    //     TPZMat2dLin *mat2d = new TPZMat2dLin(1);
    //     int ist,jst;
    //     TPZFMatrix<REAL> xk(nstate,nstate,1.),xc(nstate,nstate,0.),xf(nstate,1,0.);
    //     for(ist=0; ist<nstate; ist++) {
    //       if(nstate != 1) xf(ist,0) = 1.;
    //       for(jst=0; jst<nstate; jst++) {
    // 	if(ist != jst) xk(ist,jst) = 0.;
    //       }
    //     }
    //     mat2d->SetMaterial(xk,xc,xf);
    //     mat = mat2d;
    //   }
    TPZFMatrix<REAL> val1(1,1,0.),val2(1,1,0.);
    TPZMaterial * bc[2];
    bc[0] = mat->CreateBC(mat,-1,0,val1,val2);
    int i;
    
    val2(0,0)=1.;
    bc[1] = mat->CreateBC(mat,-2,0,val1,val2);
    bc[1]->SetForcingFunction(new TPZDummyFunction<STATE>(NeumannExp));
    
    cmesh->InsertMaterialObject(mat);
    for(i=0; i<2; i++) cmesh->InsertMaterialObject(bc[i]);
    
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    //cmesh->ExpandSolution();
    
    return cmesh;
}

void InitializeRotation(REAL alfa,REAL transx,REAL transy,TPZMatrix<REAL> &rot,TPZMatrix<REAL> &rotinv)
{
	if(rot.Rows()<3 || rot.Cols() != rot.Rows() || rot.Rows()!=rotinv.Rows()) {
		rot.Zero();
		rotinv.Zero();
		return;
	}
    for (int i=0; i<rot.Rows(); i++) {
        rot(i,i) = 1.;
        rotinv(i,i) = 1.;
    }

    REAL cosa = cos(alfa);
    REAL sina = sin(alfa);
    rot(0,0) = cosa;
    rot(1,0) = sina;
    rot(0,1) = -sina;
    rot(1,1) = cosa;
    rotinv(0,0) = cosa;
    rotinv(1,0) = -sina;
    rotinv(0,1) = sina;
    rotinv(1,1) = cosa;
	
    rot(0,3) = transx;
    rotinv(0,3) = -transx;
    rot(1,3) = transy;
    rotinv(1,3) = -transy;
}
void TransformMesh(TPZGeoMesh *gmesh,TPZMatrix<REAL> &rot)
{
    int nnod = gmesh->NNodes();
    for (int inod=0; inod<nnod; inod++) {
        TPZGeoNode &gnod = gmesh->NodeVec()[inod];
        TPZManVector<REAL,3> x(3,0.);
        gnod.GetCoordinates(x);
        TransformX(x,rot);
        gnod.SetCoord(x);
    }
}
void TransformX(TPZVec<REAL> &x,TPZMatrix<REAL> &rot)
{
	if(rot.Rows() != 4)
		return;
    TPZFNMatrix<4,REAL> xtr(4,1,1.), xtr2(4,1);
    for (int i=0; i<3; i++) {
        xtr(i,0) = x[i];
    }
    rot.Multiply(xtr, xtr2);
    for (int i=0; i<3; i++) {
        x[i] = xtr2(i,0);
    }
}
void TransformInvX(TPZVec<REAL> &x,TPZMatrix<REAL> &rotinv)
{
	if(rotinv.Rows()!=4)
		return;
    TPZFNMatrix<4,REAL> xtr(4,1,1.), xtr2(4,1);
    for (int i=0; i<3; i++) {
        xtr(i,0) = x[i];
    }
    rotinv.Multiply(xtr, xtr2);
    for (int i=0; i<3; i++) {
        x[i] = xtr2(i,0);
    }
}
