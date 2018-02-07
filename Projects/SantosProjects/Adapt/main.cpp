#include <iostream>
#include <fstream>
#include <string>
#include "pzgmesh.h"
#include "TPZVTKGeoMesh.h"
#include "pzanalysis.h"
#include "pzbndcond.h"
#include "TPZRefPatternTools.h"
#include "TPZRefPatternDataBase.h"
#include "TPZRefPattern.h"
#include "TISSMaterial.h"
#include "pzgraphmesh.h"

#include "tpzchangeel.h"
#include "TPZGeoElement.h"
#include "pzreftriangle.h"
#include "tpzgeoelrefpattern.h"
#include "TPZRefPattern.h"

using namespace pzgeom;

#include "/Applications/MATLAB_R2014a.app/extern/include/engine.h"

#include "/usr/local/include/libiomp/omp.h"


//arquivos que irão entrar no ISSM
#include "/Users/santos/Documents/issm/trunk-jpl/src/c/classes/AdaptiveMeshRefinement.h"

struct ElemNode{
    
    ElemNode(){
        fNodeId.clear();
    }
    
    std::vector<int64_t> fNodeId;
};

//----------------------------------------

struct Segment{
    
    Segment(){
        fId.clear();
    }
    
    //node i, node i+1, element j
    std::vector<int64_t> fId;
};

//----------------------------------------

bool ReadMesh(std::ifstream &MeshFile,
              std::vector<double> &x,
              std::vector<double> &y,
              std::vector<ElemNode> &ElemNodeVec,
              std::vector<Segment> &SegmentVec,
              std::vector<int64_t> &ElemIdVec);

bool ReadMesh(std::ifstream &MeshFile,int &nvertices,int &nelements,int &nsegments,double *x,double *y,double *z,int *elements,int *segments);

//----------------------------------------
TPZGeoMesh *CreateGMesh(std::vector<double> &x,
                        std::vector<double> &y,
                        std::vector<ElemNode> &ElemNodeVec,
                        std::vector<Segment> &SegmentVec);

TPZGeoMesh *CreateGMesh(std::string &MeshName);

//----------------------------------------
TPZCompMesh *CreateCompMesh(TPZGeoMesh *gmesh);

//----------------------------------------
TPZCompMesh *CopyCompMesh(TPZCompMesh *cmesh, TPZGeoMesh *gmesh);

//----------------------------------------
void Permute(TPZCompMesh * cmesh);

//----------------------------------------
void Permute2(TPZCompMesh * cmesh);

//----------------------------------------
void RefineMeshCopy(TPZGeoMesh *mesh,
                std::vector<std::pair<int,int> > &FlagVec,
                int meshID);

//----------------------------------------
TPZGeoEl * ChangeToGeoElRefPattern(TPZGeoMesh *Mesh, int64_t ElemIndex);

//----------------------------------------
bool SidesToRefine(TPZGeoEl *gel, TPZVec<int> &sidestorefine);

//----------------------------------------
void RefineClosureElements(TPZCompMesh *cmesh, std::map<int64_t,TPZVec<int> > geoelindex);

//----------------------------------------
void RefineMesh(TPZGeoMesh *gmesh,
                std::vector<int> &ElemVec,
                int nlevel);

//----------------------------------------
void RefineUniformGeoMesh(TPZGeoMesh *gmesh, std::vector<TPZVec<REAL> > &GLvec, const int &MaxLevel);

//----------------------------------------
void RefineGeoMesh(TPZGeoMesh *gmesh, std::vector<TPZVec<REAL> > &GLvec, const int &MaxLevel);

//----------------------------------------
void RefineMesh(TPZCompMesh *cmesh,
                /*std::vector<int64_t>*//*std::set<int64_t>*/std::map<int64_t, int>  &ElementIndex,
                std::vector<int64_t> &ClosureElementIndex);

//----------------------------------------
void RefineMesh2(TPZCompMesh *cmesh, std::vector<TPZVec<REAL> > &GLvec, const int &MaxLevel);

//----------------------------------------
void RefineMesh(TPZCompMesh *cmesh, std::vector<TPZVec<REAL> > &GLvec, const int &MaxLevel,REAL &alpha);


//----------------------------------------
void PrintNewMesh(std::string &MeshName,
                  TPZGeoMesh *mesh);

//----------------------------------------
void PrintNewMesh(std::ofstream &MeshFile, int &nvertices, int &nelements, int &nsegments, int &elementswidth, double *x, double *y, double *z, int * elements, int * segments);


//----------------------------------------
void PrintInitialData(std::string &DataName,
                      TPZGeoMesh *mesh);


//----------------------------------------
void FlagElements(TPZGeoMesh *gmesh, std::vector<std::pair<int,int> > &FlagVec, int step, int dstep);

//----------------------------------------
void SetElementsToRefine(TPZCompMesh *cmesh,               // mesh with solution
                         /*std::vector<int64_t>*//*std::set<int64_t>*/std::map<int64_t, int> &ElementIndex,
                         std::vector<TPZVec<REAL> > &GLvec,
                         REAL &alpha);   // elementos of cmesh to refine

//----------------------------------------
void SetRefPatterns();


//----------------------------------------
std::string toStr(int intVal);

//----------------------------------------
std::string toStr(const char *charVal);


//----------------------------------------
void LoadSolution(std::string &SolutionName, TPZCompMesh *cmesh);

//results / data for each mesh node
bool LoadSolution(std::ifstream &SolutionFile,
                  std::vector<double> &surface,
                  std::vector<double> &base,
                  std::vector<double> &bed,
                  std::vector<double> &pressure,
                  std::vector<double> &temperature,
                  std::vector<double> &vx,
                  std::vector<double> &vy,
                  std::vector<double> &masklevelset);

bool LoadSolution(std::ifstream &SolutionFile, double **vx, double **vy, double **masklevelset);

//----------------------------------------
void TransferSolution(TPZCompMesh *cmesh, TPZCompMesh *cmesh2);

//----------------------------------------
void TransferSolution(TPZFMatrix<STATE> &sol2, TPZFMatrix<STATE> &sol1,
                      int64_t &nodeIndex2, int64_t &nodeIndex1);

//----------------------------------------
void GetVariableNames(TPZStack<std::string> &scalnames, TPZStack<std::string> &vecnames);


//----------------------------------------
void SetFileNames(const int &h,
                  const int &step,
                  bool ChangeSolutionFile,
                  std::string &strSolutionFile,
                  std::string &strMeshFile,
                  std::string &strDataFile,
                  std::string &strSolutionMatlab,
                  std::string &strMeshMatlab,
                  std::string &strDataMatlab);

//----------------------------------------
bool XDif(TPZVec<REAL> &X1, TPZVec<REAL> &X2, REAL &Tol);

//----------------------------------------
bool FindNodeIndex( TPZVec<REAL> &NodeX, int64_t &NodeIndexOnSons, TPZVec<TPZGeoEl *> unrefinedSons );

//----------------------------------------
void PrintMesh(const int &step,
              TPZGeoMesh *gmesh);

//----------------------------------------
void BuildPermuteVector(TPZCompMesh *cmesh);

//----------------------------------------
void GeoToConnect(TPZCompMesh * cmesh, TPZVec<int64_t> &geotoconnect);

//----------------------------------------
void CheckSequenceNumber(TPZCompMesh * cmesh);

//----------------------------------------
void SetModelNum(const int &level, std::string &strModelNum);

//----------------------------------------
std::string GetModelNum(const int &level);

//----------------------------------------
std::string GetExperimentPrefix(const int &Steps);

//----------------------------------------
std::string GetRepository(const int &level);

//----------------------------------------
std::string GetClusterName();


/// path do projeto a ser executado
#define WORKPATH "/Users/santos/Documents/projects/Misomip/"

#define MY_REFPATTERNDIR "/Users/santos/Documents/neopz/Projects/SantosProjects/Adapt/rpt/"

#define BUFSIZE 1024

#define CALLMATLAB

#define IS_UNIFORM false

#define EXP_APLHA 1.0

#define RUN_3D 0

#define MAX_GL_DISTANCE 40000.0

//se the excel file (step)
#define MAX_TIME_STEP 180

//see the excel file (yr)
#define FINAL_TIME 5

#define InitLevel 5

#define LastLevel 5

// use 4, 8 or 12
#define VertLayer 12

//see the excel file
#define VTKPrintRatio 10 //100

//use true to run Exp1 and Exp2
//use false to run SteadyState
#define IsRestart true

// cluster 0 = local
// cluster 1 = capitao
// cluster 2 = labmec-mac-pro
// cluster 3 = whopperdooper

#define Cluster 1

// this is used to set the files name (ex. mesh(step).txt, solution(step).txt etc)
// use 0 to run since the beginning
// use the last TimeStep to run restart (each experiment for example)
#define InitStep 51

// set the variable steps
// 'Transient_Steadystate' = 3
// 'Transient_Steadystate_3D_HO' = 5
// 'Experiment0' = 6
// 'Experiment1' = 7
// 'Experiment1a' = 8
// 'Experiment1a_short' = 9 NÃO RODAR!
// 'Experiment1r' = 10
// 'Experiment2' = 11
// 'Experiment2a' = 12
// 'Experiment2a_short' = 13 NÃO RODAR!
// 'Experiment2r' = 14

#define STEPS_SIM 14

#define IsViscous true


int mainMISOMIP(int argc, char *argv[]);

int mainMISMIP(int argc, char *argv[]);

int mainAMR(int argc, char *argv[]);

//just to debug anything
int mainDebug(int argc, char *argv[]);

TPZGeoMesh * CreateNewGMesh(TPZGeoMesh *cp);

int main(int argc, char *argv[]){
    
    mainDebug(argc,argv);
    
    int value = 1;
    

  //  printf("Hello from thread %d, nthreads %d\n", omp_get_thread_num(), omp_get_num_threads());
    
    //para rodar o MISOMIP
    //int value = mainMISOMIP(argc, argv);
    
    
    //para testar o AMR
   // int value = mainAMR(argc, argv);
    
    //para debugar
    //int value = mainDebug(argc, argv);
    
    
 //   int value = omp_get_num_threads();

//#pragma omp parallel num_threads(4)
//#pragma omp for
 //   for(int n=0; n<10; ++n){
 //       double dvalue = n*n;
 //       std::cout << dvalue << std::endl;
  //  }
    
    return value;
}
//###################################################
int mainDebug(int argc, char *argv[]){

    TPZGeoMesh *gmesh   = new TPZGeoMesh();
    int nnodes          = 4;
    const int mat       = 1;
    const int reftype   = 1; //1 is refpatern
    const int hmax      = 2;
    int64_t index;
    TPZAutoPointer<TPZRefPattern> refp = NULL;
    TPZManVector<REAL,3> coord(3,0.);
    TPZManVector<int64_t,3> tria(3,0);
    TPZVec<TPZGeoEl *> sons;
    SetRefPatterns();
    
    gmesh->NodeVec().Resize( nnodes );
    coord[0] = 0.; coord[1] = 0.; //nó 0
    gmesh->NodeVec()[0].SetCoord(coord); gmesh->NodeVec()[0].SetNodeId(0);
    coord[0] = 10.; coord[1] = 0.; //nó 1
    gmesh->NodeVec()[1].SetCoord(coord); gmesh->NodeVec()[1].SetNodeId(1);
    coord[0] = 10.; coord[1] = 10.; //nó 2
    gmesh->NodeVec()[2].SetCoord(coord); gmesh->NodeVec()[2].SetNodeId(2);
    coord[0] = 0.; coord[1] = 10.; //nó 3
    gmesh->NodeVec()[3].SetCoord(coord); gmesh->NodeVec()[3].SetNodeId(3);
    
    tria[0] = 3; tria[1] = 1; tria[2] = 2; //element 0
    gmesh->CreateGeoElement(ETriangle,tria,mat,index,reftype);
    gmesh->ElementVec()[index]->SetId(index);
    tria[0] = 0; tria[1] = 1; tria[2] = 3; //element 1
    gmesh->CreateGeoElement(ETriangle,tria,mat,index,reftype);
    gmesh->ElementVec()[index]->SetId(index);

    gmesh->BuildConnectivity();

    std::ofstream file0("/Users/santos/Desktop/mesh0.txt");     gmesh->Print(file0);
    std::ofstream filevtk0("/Users/santos/Desktop/mesh0.vtk");  TPZVTKGeoMesh::PrintGMeshVTK(gmesh,filevtk0 );
    
    /* Refinement process */
    for(int l=0;l<hmax;l++){
        const int64_t nelem = gmesh->NElements();
        for(int64_t i=0;i<nelem;i++){
            if(gmesh->Element(i)->HasSubElement()) continue;
            gmesh->Element(i)->Divide(sons); sons.clear();
        }
    }
    if(1){
        refp = TPZRefPatternTools::PerfectMatchRefPattern(gmesh->Element(1));
        if(refp){
            gmesh->Element(1)->SetRefPattern(refp);
            gmesh->Element(1)->Divide(sons);sons.clear();
        }else{
            DebugStop();
        }
    }
    gmesh->BuildConnectivity();

    std::ofstream file1("/Users/santos/Desktop/mesh1.txt");     gmesh->Print(file1);
    std::ofstream filevtk1("/Users/santos/Desktop/mesh1.vtk");  TPZVTKGeoMesh::PrintGMeshVTK(gmesh,filevtk1 );

    /* Delete process */
    index=2;
    gmesh->Element(index)->GetHigherSubElements(sons);
    gmesh->Element(index)->ResetSubElements();
    for (int i=0;i<sons.size();i++){
        gmesh->DeleteElement(sons[i],sons[i]->Index());
    }
    sons.clear();
    gmesh->BuildConnectivity();

    std::ofstream file2("/Users/santos/Desktop/mesh2.txt");     gmesh->Print(file2);
    std::ofstream filevtk2("/Users/santos/Desktop/mesh2.vtk");  TPZVTKGeoMesh::PrintGMeshVTK(gmesh,filevtk2 );
    
    /* Refine again*/
    if(1)
    { //Divide again, uniform
        gmesh->Element(index)->Divide(sons); sons.clear();
    }
    else
    { //Divide again, not uniform
        if(gmesh->Element(index)->HasSubElement()) DebugStop();
        refp = TPZRefPatternTools::PerfectMatchRefPattern(gmesh->Element(index));
        if(refp){
            gmesh->Element(index)->SetRefPattern(refp);
            gmesh->Element(index)->Divide(sons); sons.clear();
        }else{
            DebugStop();
        }
    }
    gmesh->BuildConnectivity();

    std::ofstream file3("/Users/santos/Desktop/mesh3.txt");     gmesh->Print(file3);
    std::ofstream filevtk3("/Users/santos/Desktop/mesh3.vtk");  TPZVTKGeoMesh::PrintGMeshVTK(gmesh,filevtk3 );
    
    return 0;
}
TPZGeoMesh * CreateNewGMesh(TPZGeoMesh *cp)
{
    TPZGeoMesh *newgmesh = new TPZGeoMesh();
    newgmesh->CleanUp();
    
    int nnodes  = cp->NNodes();
    int nelem   = cp->NElements();
    int mat     = 1;
    int reftype = 1;
    int64_t index;
    
    //nodes
    newgmesh->NodeVec().Resize(nnodes);
    for(int i=0;i<nnodes;i++) newgmesh->NodeVec()[i] = cp->NodeVec()[i];
    
    //elements
    for(int i=0;i<nelem;i++){
        TPZGeoEl * geoel = cp->Element(i);
        TPZManVector<int64_t> elem(3,0);
        for(int j=0;j<3;j++) elem[j] = geoel->NodeIndex(j);
     
        newgmesh->CreateGeoElement(ETriangle,elem,mat,index,reftype);
        newgmesh->ElementVec()[index]->SetId(geoel->Id());
        
        TPZGeoElRefPattern<TPZGeoTriangle>* newgeoel = dynamic_cast<TPZGeoElRefPattern<TPZGeoTriangle>*>(newgmesh->ElementVec()[index]);
        
        //old neighbourhood
        const int nsides = TPZGeoTriangle::NSides;
        TPZVec< std::vector<TPZGeoElSide> > neighbourhood(nsides);
        TPZVec<int64_t> NodesSequence(0);
        for(int s = 0; s < nsides; s++)
        {
            neighbourhood[s].resize(0);
            TPZGeoElSide mySide(geoel,s);
            TPZGeoElSide neighS = mySide.Neighbour();
            if(mySide.Dimension() == 0)
            {
                int64_t oldSz = NodesSequence.NElements();
                NodesSequence.resize(oldSz+1);
                NodesSequence[oldSz] = geoel->NodeIndex(s);
            }
            //if(TPZChangeEl::CreateMiddleNodeAtEdge(Mesh, ElemIndex, s, midN))
            //{
            //    int64_t oldSz = NodesSequence.NElements();
            //    NodesSequence.resize(oldSz+1);
            //    NodesSequence[oldSz] = midN;
            //}
            while(mySide != neighS)
            {
                neighbourhood[s].push_back(neighS);
                neighS = neighS.Neighbour();
            }
        }
        
        //inserting in new element
        for(int s = 0; s < nsides; s++)
        {
            TPZGeoEl * tempEl = newgeoel;
            TPZGeoElSide tempSide(newgeoel,s);
            int byside = s;
            for(uint64_t n = 0; n < neighbourhood[s].size(); n++)
            {
                TPZGeoElSide neighS = neighbourhood[s][n];
                tempEl->SetNeighbour(byside, neighS);
                tempEl = neighS.Element();
                byside = neighS.Side();
            }
            tempEl->SetNeighbour(byside, tempSide);
        }
        
        int64_t fatherindex = geoel->FatherIndex();
        if(fatherindex>-1) newgeoel->SetFather(fatherindex);
        
        if(!geoel->HasSubElement()) continue;
        
        int nsons = geoel->NSubElements();

        TPZAutoPointer<TPZRefPattern> ref = gRefDBase.GetUniformRefPattern(ETriangle);
        newgeoel->SetRefPattern(ref);
        
        for(int j=0;j<nsons;j++){
            TPZGeoEl* son = geoel->SubElement(j);
            if(!son){
                DebugStop();
            }
            newgeoel->SetSubElement(j,son);
        }
        
    }

    newgmesh->BuildConnectivity();
    
    return newgmesh;
    
}
//###################################################
int mainAMR(int argc, char *argv[]){
  
//    int isrestart   = atoi(argv[1]);
//    std::ifstream SolutionFile(argv[2]);
//    std::ofstream NewMeshFile(argv[3]);
//    std::string AMRfile = argv[4];
//    
//    /* starting AMR */
//    int elementswidth = 3; //itapopo malha 2D tringular
//    AdaptiveMeshRefinement *AMR = new AdaptiveMeshRefinement();
//    
//    if(isrestart){
//        
//        /* type of process. 0: refine the same mesh; 1 refine the mesh 0 (unrefine process) */
//        int type_process = atoi(argv[5]); /* 0=refinement; 1=unrefinement*/
//        
//        /* read refpattern*/
//        SetRefPatterns();
//        
//        /* start from AdaptiveMeshRefinement file */
//        TPZFileStream fstr;
//        fstr.OpenRead(AMRfile.c_str());
//        TPZSavable *sv = TPZSavable::Restore(fstr,0);
//        AMR = dynamic_cast<AdaptiveMeshRefinement*>(sv);
//        
//        /* read solution of the mesh i */
//        double *vx;
//        double *vy;
//        double *masklevelset;
//        bool IsOk = LoadSolution(SolutionFile, &vx, &vy, &masklevelset);
//        if(!IsOk) DebugStop();
//        
//        /* Refine the mesh */
//        double *newx;
//        double *newy;
//        double *newz;
//        int *newelements;
//        int *newsegments;
//        int newnumberofvertices, newnumberofelements, newnumberofsegments;
//        AMR->ExecuteRefinement(type_process,vx,vy,masklevelset,newnumberofvertices,newnumberofelements,newnumberofsegments,&newx,&newy,&newz,&newelements,&newsegments);
//    
//        /* Printing new mesh */
//        //PrintNewMesh(NewMeshFile, newnumberofvertices, newnumberofelements, newnumberofsegments,elementswidth, newx, newy, newz, newelements,newsegments);
//        
//        /* delete data*/
//        if(newx)    delete newx;
//        if(newy)    delete newy;
//        if(newz)    delete newz;
//        if(vx)      delete vx;
//        if(vy)      delete vy;
//        if(masklevelset) delete masklevelset;
//       // for(int64_t i=0;i<newnumberofelements;i++) delete newelements[i];
//        if(newelements) delete newelements;
//        //for(int64_t i=0;i<newnumberofsegments;i++) delete newsegments[i];
//        if(newsegments) delete newsegments;
//        
//    } else {
//        /* read initial mesh (mesh 0) */
//        double *x;
//        double *y;
//        double *z;
//        int *elements;
//        int *segments;
//        int nvertices, nelements, nsegments;
//        int hmax        = atoi(argv[6]);
//        std::ifstream InitialMeshFile(argv[7]);
//        //bool IsOk = ReadMesh(InitialMeshFile,nvertices,nelements,nsegments,&x,&y,&z,&elements,&segments);
//        //if(!IsOk) DebugStop();
//    
//        /* create initial mesh (father mesh and previous mesh. Previous mesh is equal to father mesh in this initialization*/
//        AMR->SetHMax(hmax);
//        //AMR->CreateInitialMesh(nvertices, nelements, nsegments, elementswidth, x, y, z, elements, segments);
//        
//        /* delete data*/
//        if(x) delete x;
//        if(y) delete y;
//        if(z) delete z;
//       // for(int64_t i=0;i<nelements;i++) delete elements[i];
//        if(elements) delete elements;
//       // for(int64_t i=0;i<nsegments;i++) delete segments[i];
//        if(segments) delete segments;
//    }
//    
//    /* save AMR in hard disc */
//    TPZFileStream fstr;
//    fstr.OpenWrite(AMRfile.c_str());
//    AMR->Write(fstr,1);
//    if(AMR) delete AMR;
    
    return 0;
}

//###################################################
int mainMISOMIP(int argc, char *argv[]){
    
    /// set the Engine MatLab object
    Engine *ep;
    char buffer[BUFSIZE+1];
    buffer[BUFSIZE] = '\0';
    
    /// set the auxiliar strings
    int step;
    std::string strSolutionMatlab; //solution from ISSM
    std::string strMeshMatlab; // mesh from ISSM
    std::string strDataMatlab; //data from ISSM
    
    std::string strDataFile;
    std::string strSolutionFile;
    std::string strMeshFile;
    
    /*
     * Call engOpen with a NULL string. This starts a MATLAB process
     * on the current host using the command "matlab".
     * engOpen use csh; so, command "matlab" needs be into the csh PATH.
     * A simple way is: sudo ln -s /Applications/MATLAB_R2014a.app/bin/matlab /usr/bin/matlab
     */
//    if (!(ep = engOpen(NULL))) {
//        fprintf(stderr, "\nCan't start MATLAB engine\n");
//        return EXIT_FAILURE;
//    }
//    
    // Set the ISSM_DIR, system variable - ITAPOPO isso deveria ser lido pelo processo que dispara o matlab
//    engEvalString(ep, "setenv('ISSM_DIR', '/Users/santos/Documents/issm/trunk-jpl');");
    
    // cd to workpath and run the initial setup for matlab
 //   engEvalString(ep, "cd /Users/santos/Documents/projects/Misomip");
  //  engEvalString(ep, "run ../scripts/init.m");
    
    /**
     * Looping for each h level. hlevel is the level of refinement.
     * This method was used to run automticaly severel levels of refinement.
     */
    for(int hlevel=InitLevel; hlevel <= LastLevel; hlevel++){
        
        //clear the workspace
    //    engEvalString(ep, "clear");
        
        //set the buffer
    //    engOutputBuffer(ep, buffer, BUFSIZE);
        
        //set the cluster
        std::string strClusterName = GetClusterName();
    //    engEvalString(ep,strClusterName.c_str());
        
        //set the FinalTime
        std::string strFinalTime = "FinalTime=" + toStr(FINAL_TIME) + ";";
     //   engEvalString(ep, strFinalTime.c_str());
        
        //set the Model Number. See runme.m
        std::string strModelNum;
        SetModelNum(hlevel, strModelNum);
   //     engEvalString(ep, strModelNum.c_str());
        
        std::string strRun3D = "run3D=" + toStr(RUN_3D);
     //   engEvalString(ep, strRun3D.c_str());
        
        std::string strVertLayer = "VertLayer=" + toStr(VertLayer);
       // engEvalString(ep, strVertLayer.c_str());
        
        //set the parfile
        //engEvalString(ep, "parfile='./Exp_Par/Mismip.par';");
        
        //set the variable to run with adapted meshes
        //engEvalString(ep, "RefineMesh=1;");
        
        //execute the runme.m to generate the first mesh, mesh0
        step=InitStep;
        SetFileNames(hlevel, step, IsRestart,strSolutionFile, strMeshFile,strDataFile,
                     strSolutionMatlab,strMeshMatlab,strDataMatlab);

        // mesh file and solution file to use in C++
        std::cout << " ----- Executing runme -----" << std::endl;
        //engEvalString(ep, strDataMatlab.c_str());
        //engEvalString(ep, strMeshMatlab.c_str());
        //engEvalString(ep, strSolutionMatlab.c_str());

        if(IsRestart){
            //engEvalString(ep, "steps=[];");
            //engEvalString(ep, "runme2");
            //engEvalString(ep, "RestartSimulations");
        }else{
            //engEvalString(ep, "steps=[1, 2];");
            //engEvalString(ep, "runme2");
        }
        
        /// it is necessary to use the mesh0 to generate the FatherMesh
        if(IsRestart){
            strMeshFile = GetRepository(hlevel) + "/mesh/mesh" + toStr(0) + ".txt";
        }
        
        ///generating the father mesh
        TPZGeoMesh *FatherGMesh = CreateGMesh(strMeshFile);
        
        ///generating comp mesh  and loading results
        TPZGeoMesh *gmesh = new TPZGeoMesh(*FatherGMesh);
        TPZCompMesh *cmesh = CreateCompMesh(gmesh);
        LoadSolution(strSolutionFile, cmesh);
        
        ///setting refinement patterns
        SetRefPatterns();
        
        ///print solution to Paraview
        TPZStack<std::string> scalnames, vecnames;
        bool optimizeBandwidth = false; //impede a renumeracao das equacoes do problema
        const int postProcessResolution = 0;//define resolucao do pos processamento
        const int dim = 2;
        TPZAnalysis an(cmesh, optimizeBandwidth); //cria objeto de analise que gerenciaria a analise do problema
        GetVariableNames(scalnames, vecnames);
        std::string plotfile;
        std::string vtkfile, vtkfile_restart;
        vtkfile = toStr(WORKPATH) + GetRepository(hlevel) + "/vtk/" + GetExperimentPrefix(STEPS_SIM) + "_solution.vtk";
        vtkfile_restart = toStr(WORKPATH) + GetRepository(hlevel) + "/vtk/" + GetExperimentPrefix(STEPS_SIM) + "_solution_restart.vtk";
        
        if(IsRestart) {
            plotfile = vtkfile_restart;
        } else {
            plotfile = vtkfile;
        }
        
        an.DefineGraphMesh(dim, scalnames, vecnames, plotfile);//define malha grafica
        an.SetStep(step);
        an.PostProcess(postProcessResolution);//realiza pos processamento
        
        plotfile = vtkfile;
        
        // set elements to refine
        std::map<int64_t, int> ElementIndex;
        
        // refine mesh
        TPZGeoMesh *gmesh1 = new TPZGeoMesh(*FatherGMesh);
        TPZCompMesh *cmesh1 = CreateCompMesh(gmesh1);
        
        const int MaxLevel = hlevel; //máximo nível de refinamento
        std::vector<TPZVec<REAL> > GLvec;
        REAL alpha = 1.0;
        if(IsRestart) alpha = 2.0;
        SetElementsToRefine(cmesh, ElementIndex, GLvec, alpha);
        RefineMesh(cmesh1, GLvec, MaxLevel, alpha);
        
        #ifdef PZDEBUG
        plotfile = toStr(WORKPATH) + GetRepository(hlevel) + "/vtk/solution_restart_PosRefine.vtk";
        an.SetCompMesh(cmesh1, optimizeBandwidth);
        an.CloseGraphMesh();
        an.DefineGraphMesh(dim, scalnames, vecnames, plotfile);//define malha grafica
        an.SetStep(step);
        an.PostProcess(postProcessResolution);//realiza pos processamento
        plotfile = vtkfile;
        #endif
        
        #ifdef PZDEBUG
        //std::ofstream fileDebugMesh("/Users/santos/Documents/projects/Misomip/MeshDebug.txt");
        //cmesh1->Print(fileDebugMesh);
        #endif
        
        Permute2(cmesh1);
        
        TransferSolution(cmesh, cmesh1);
        
        /// print the new mesh and the initial data
        step++;
        SetFileNames(hlevel, step, false,strSolutionFile, strMeshFile,strDataFile,
                     strSolutionMatlab,strMeshMatlab,strDataMatlab);
        
        PrintNewMesh(strMeshFile, gmesh1);
        PrintInitialData(strDataFile, gmesh1);
        
        delete cmesh;
        delete gmesh;
        
        gmesh = gmesh1;
        cmesh = cmesh1;
        
        an.SetCompMesh(cmesh, optimizeBandwidth);
        an.SetStep(step);
        
        std::string strSetSteps = "steps=" + toStr(STEPS_SIM);
        //engEvalString(ep, strSetSteps.c_str());
    
        //std::string strRun3D = "run3D=" + toStr(RUN_3D);
        //engEvalString(ep, strRun3D.c_str());
        
        /// running time step. In each TimeStep, the mesh is refined
        const int MaxTimeStep = MAX_TIME_STEP;
        
        int RealTime;;
        
        for(int TimeStep = 1; TimeStep <= MaxTimeStep; TimeStep++ ){
            
            std::cout << std::endl;
            std::cout << "RUNNING TIME STEP = " << TimeStep << "/" << MaxTimeStep << ", HLevel = " << hlevel << std::endl;
            
            //setting the real time, used to save the .mat files
            RealTime = TimeStep*FINAL_TIME; //in yrs
            std::string strRealTime = "time="+toStr(RealTime);
            //engEvalString(ep, strRealTime.c_str());
            
            //execute the main to run
            std::cout << std::endl;
            std::cout << " ----- Executing runme -----" << std::endl;
            //engEvalString(ep, strDataMatlab.c_str());
            //engEvalString(ep, strMeshMatlab.c_str());
            //engEvalString(ep, strSolutionMatlab.c_str());
            //engEvalString(ep, "runme2");
            
            //load the solution
            LoadSolution(strSolutionFile, cmesh);
            
            //print solution to vtk
            if(TimeStep%VTKPrintRatio == 0){
                an.CloseGraphMesh();
                an.DefineGraphMesh(dim, scalnames, vecnames, plotfile);//define malha grafica
                an.PostProcess(postProcessResolution);
            }
            
            //if( TimeStep == MaxTimeStep ) break;
            
            // generating the new mesh
            TPZGeoMesh *gmesh2 = new TPZGeoMesh(*FatherGMesh);
            TPZCompMesh *cmesh2 = CreateCompMesh(gmesh2);
            
            //refine the mesh
            alpha =1.0;
            SetElementsToRefine(cmesh, ElementIndex, GLvec, alpha);
            RefineMesh(cmesh2, GLvec, MaxLevel, alpha);
            
#ifdef PZDEBUG
            //std::ofstream fileDebugMesh("/Users/santos/Documents/projects/Misomip/MeshDebug.txt");
            //cmesh2->Print(fileDebugMesh);
#endif
            
            Permute2(cmesh2);
            
            /// transfering the solution to the new mesh
            TransferSolution(cmesh, cmesh2);
            
            //print the new mesh
            step++;
#ifdef PZDEBUG
            PrintMesh(step, gmesh2);
#endif
            SetFileNames(hlevel, step, false,strSolutionFile, strMeshFile,strDataFile,
                         strSolutionMatlab,strMeshMatlab,strDataMatlab);
            PrintNewMesh(strMeshFile, gmesh2);
            PrintInitialData(strDataFile, gmesh2);
            
            ///setting the new mesh pointer
            delete cmesh;
            delete gmesh;
            
            gmesh = gmesh2;
            cmesh = cmesh2;
            
            an.SetCompMesh(cmesh, optimizeBandwidth);
            an.SetStep(step);
            
        }
        
        delete cmesh;
        delete gmesh;
        
    }//for hlevel
    
    //close the workspace and close the engine
    //engEvalString(ep, "close;");
    //engClose(ep);
    
    std::cout << "FINISHED!" << std::endl;
    
    return 0;
}
//###################################################
std::string GetExperimentPrefix(const int &Steps){
    
    
    std::string experiment;
    if( Steps==3 | Steps==5 ){
        
        experiment="SteadyState";

    } else if(Steps==7 | Steps==8 | Steps==10){
        
        experiment="Exp1";
        
    } else if(Steps==11){
        
        experiment = "Ice2r";
        
    } else if(Steps==12){
        
        experiment = "Ice2ra";
    
    } else if(Steps==14){
        
        experiment="Ice2rr";
        
    } else {
        DebugStop();
    }
    
    return experiment;
}

//###################################################
std::string GetModelNum(const int &level){
    
    int prefix = IsViscous ? 1 : 2;
    
    std::string ModelNum = toStr(prefix) + toStr(level);
    
    return ModelNum;
    
}

//###################################################
std::string GetClusterName(){
    
    std::string ClusterName;
    switch (Cluster) {
        case 0:
            ClusterName = "'recruta'";
            break;
        
        case 1:
            ClusterName = "'capitao'";
            break;
            
        case 2:
            ClusterName = "'macpro'";
            break;
            
        case 3:
            ClusterName = "'whopper'";
            break;
            
        default:
            break;
    }
    
    std::string TotalName;
    
    TotalName = "clustername="+ClusterName;
    
    return TotalName;

}

//###################################################

void SetModelNum(const int &level, std::string &strModelNum){
    
    //int prefix = IsViscous ? 1 : 2;
    
//    strModelNum = "modelnum=" + toStr(prefix) + toStr(level) + ";";

    strModelNum = "modelnum=" + GetModelNum(level) + ";";
}

//###################################################

void CheckSequenceNumber(TPZCompMesh * cmesh){
    
    TPZVec<int64_t> geotoconnect;
    GeoToConnect(cmesh, geotoconnect);
    
    int error = 0;
    {
        std::ofstream out("../geoseq.txt");
        
        int64_t ngeonodes = geotoconnect.size();
        
        for(int i = 0; i < ngeonodes;i++){
            int64_t connectIndex = geotoconnect[i];
            TPZConnect &c= cmesh->ConnectVec()[connectIndex];
            int64_t seqnum = c.SequenceNumber();
            out << "i " << i << " seqnum " << seqnum << std::endl;
            if(seqnum != i) {
                error = 1;
            }
        }
    }
    if (error) {
        DebugStop();
    }
}

////////////////
int mainMISMIP(int argc, char *argv[])
{
    
    DebugStop();
    
    /// set the Engine MatLab object
    Engine *ep;
    char buffer[BUFSIZE+1];
    buffer[BUFSIZE] = '\0';
    int step;
    std::string strSolutionMatlab;
    std::string strMeshMatlab;
    std::string strDataMatlab;
    
    std::string strDataFile;
    std::string strSolutionFile;
    std::string strMeshFile;
    
 
#ifdef CALLMATLAB
    /*
     * Call engOpen with a NULL string. This starts a MATLAB process
     * on the current host using the command "matlab".
     */
//    if (!(ep = engOpen(NULL))) {
//        fprintf(stderr, "\nCan't start MATLAB engine\n");
//        return EXIT_FAILURE;
//    }
//    
    /**
     * Looping for each h level. hlevel is the level of refinement.
     * This method was used to run automticaly severel levels of refinement.
     */
for(int hlevel=4; hlevel < 5; hlevel++){//ITAPOPO
    
    std::cout << "#####     RUNNING H LEVEL = " << hlevel << "     #####" << std::endl;
    
    /// generating the files in the WORKPATH folder
   /* 
    itapopo arrumar a geração dos arquivos
    char teste;
    const char nfile = hlevel;
    std::strcat(&teste,MKSOL);
    std::strcat(&teste,&nfile);
    
    system(&teste);
    
    system(MKDEBUG);
    system(MKSOL);
    system(MKDATA);
    system(MKMESH);
    system(MKVTK);*/
    
    //clear the workspace
//    engEvalString(ep, "clear");
    
    //set the buffer
//    engOutputBuffer(ep, buffer, BUFSIZE);
    
    // Set the ISSM_DIR, system variable - ITAPOPO isso deveria ser lido pelo processo que dispara o matlab
//    engEvalString(ep, "setenv('ISSM_DIR', '/Users/santos/Documents/issm/trunk');");
    
    // add the path of ISSM and MatLab codes
//    engEvalString(ep, "addpath /Users/santos/Documents/issm/trunk/bin /Users/santos/Documents/issm/trunk/lib /Users/santos/Documents/NeoPZ/neopz/Projects/SantosProjects/Adapt/matlab");
    
    //change to the WORKPATH
    std::string cdWORKPATH = "cd " + toStr(WORKPATH);
//    engEvalString(ep, cdWORKPATH.c_str());//"cd /Users/santos/Documents/_PROJETOS/Criosfera/Adapt2D/");

    //set the variable isRun: 0, print mesh 0; 1, run the ISSM with new mesh
 //   engEvalString(ep, "isRun=0;");
    
    //set the expfile
//    engEvalString(ep, "expfile='Domain.exp';");
    
    //set the parfile
//    engEvalString(ep, "parfile='BC.par';");
    
    //set the resolution to generate the mesh0
//    engEvalString(ep, "resolution=20000;");//5000
    
#endif
    //execute the main to generate the first mesh, mesh0
    step=0;
    SetFileNames(hlevel, step, false,strSolutionFile, strMeshFile,strDataFile,
                 strSolutionMatlab,strMeshMatlab,strDataMatlab);
#ifdef CALLMATLAB
   // engEvalString(ep, strMeshMatlab.c_str());
    //engEvalString(ep, strSolutionMatlab.c_str());
    //engEvalString(ep, "main");
    //print the buffer
    //printf("%s", buffer);
#endif
    
    ///generating the father mesh
    TPZGeoMesh *FatherGMesh = CreateGMesh(strMeshFile);
    
    ///generating comp mesh  and loading results
    TPZGeoMesh *gmesh = new TPZGeoMesh(*FatherGMesh);
    TPZCompMesh *cmesh = CreateCompMesh(gmesh);
    LoadSolution(strSolutionFile, cmesh);
    
#ifdef PZDEBUG
    PrintMesh(step, gmesh);
#endif
    
    ///setting refinement patterns
    SetRefPatterns();
    
    ///print solution to Paraview
    TPZStack<std::string> scalnames, vecnames;
    bool optimizeBandwidth = false; //impede a renumeracao das equacoes do problema
    const int postProcessResolution = 0;//define resolucao do pos processamento
    const int dim = 2;
    TPZAnalysis an(cmesh, optimizeBandwidth); //cria objeto de analise que gerenciaria a analise do problema
    GetVariableNames(scalnames, vecnames);
    std::string plotfile = toStr(WORKPATH) + "vtk" + toStr(hlevel) + "/solution.vtk";

    an.DefineGraphMesh(dim, scalnames, vecnames, plotfile);//define malha grafica
    an.PostProcess(postProcessResolution);//realiza pos processamento
   
    // set elements to refine
    /*std::vector<int64_t>*/ /*std::set<int64_t>*/std::map<int64_t, int> ElementIndex;
    std::vector<int64_t> ClosureElementIndex;//elements to close (avoid hanging nodes)
    
    // refine mesh
    /// itapopo
    TPZGeoMesh *gmeshTeste = new TPZGeoMesh(*FatherGMesh);
    TPZCompMesh *cmeshTeste = CreateCompMesh(gmeshTeste);
    
    ///itapopo
    const int MaxLevel = hlevel;//1;///máximo nível de refinamento
    std::vector<TPZVec<REAL> > GLvec;
    REAL alpha =1.0;
    SetElementsToRefine(cmesh, ElementIndex,GLvec, alpha);
    //RefineMesh(cmeshTeste, ElementIndex, ClosureElementIndex);
    RefineMesh(cmeshTeste, GLvec, MaxLevel, alpha);
    Permute(cmeshTeste);
    
    
    ///itapopo
    TransferSolution(cmesh, cmeshTeste);
    ///itapopo
    
    
    /// print the new mesh and the initial data
    step++;
#ifdef PZDEBUG
    PrintMesh(step, gmeshTeste);
#endif

    SetFileNames(hlevel, step, false,strSolutionFile, strMeshFile,strDataFile,
                 strSolutionMatlab,strMeshMatlab,strDataMatlab);
    
    PrintNewMesh(strMeshFile, gmeshTeste);
    PrintInitialData(strDataFile, gmeshTeste);
    
    /// Run with new mesh
    
    //set the variable isRun: 0, print mesh 0; 1, run the ISSM with new mesh
   // engEvalString(ep, "isRun=1;");
    const int runmax = 10;///5000 * 10 = 50000
    
    
    ///itapopo
    delete cmesh;
    delete gmesh;
    
    gmesh = gmeshTeste;
    cmesh = cmeshTeste;
    
    an.SetCompMesh(cmesh, optimizeBandwidth);
    an.SetStep(step);
    
    ///itapopo apagar
    ///testando o TransferSolution
    //an.CloseGraphMesh();
    //an.DefineGraphMesh(dim, scalnames, vecnames, plotfile);//define malha grafica
    //an.PostProcess(postProcessResolution);
    
    ///itapopo
    
    
    
    for(int irun = 0; irun < runmax; irun++ ){
        
        std::cout << "RUNNING STEP = " << step << std::endl;
        
        //execute the main to run
//        engEvalString(ep, strDataMatlab.c_str());
//        engEvalString(ep, strMeshMatlab.c_str());
//        engEvalString(ep, strSolutionMatlab.c_str());
//        engEvalString(ep, "main");
//    
        //print the buffer
        //printf("%s", buffer);
    
        //load the solution
        LoadSolution(strSolutionFile, cmesh);
        
        //print solution to vtk
        an.CloseGraphMesh();
        an.DefineGraphMesh(dim, scalnames, vecnames, plotfile);//define malha grafica
        an.PostProcess(postProcessResolution);
    
        // generating the new mesh
        TPZGeoMesh *gmesh2 = new TPZGeoMesh(*FatherGMesh);
        TPZCompMesh *cmesh2 = CreateCompMesh(gmesh2);
        
        //refine the mesh
        SetElementsToRefine(cmesh, ElementIndex, GLvec, alpha);
        //RefineMesh(cmesh2, ElementIndex, ClosureElementIndex);//cmesh
        RefineMesh(cmesh2, GLvec, MaxLevel,alpha);
        Permute(cmesh2);//cmesh
        
        /// transfering the solution to the new mesh
        TransferSolution(cmesh, cmesh2);
        
        //print the new mesh
        step++;
#ifdef PZDEBUG
        PrintMesh(step, gmesh2);
#endif
        SetFileNames(hlevel, step, false,strSolutionFile, strMeshFile,strDataFile,
                     strSolutionMatlab,strMeshMatlab,strDataMatlab);
        PrintNewMesh(strMeshFile, gmesh2);
        PrintInitialData(strDataFile, gmesh2);
        
        ///setting the new mesh pointer
        delete cmesh;
        delete gmesh;
        
        gmesh = gmesh2;
        cmesh = cmesh2;
        
        an.SetCompMesh(cmesh, optimizeBandwidth);
        an.SetStep(step);
        
        ///itapopo apagar
        ///testando o TransferSolution
        //an.CloseGraphMesh();
        //an.DefineGraphMesh(dim, scalnames, vecnames, plotfile);//define malha grafica
        //an.PostProcess(postProcessResolution);
        
    }
    
    //save the model
    std::string SaveModel = "save model" + toStr(hlevel) + " md;";
//    engEvalString( ep, SaveModel.c_str() );//"save model md;"); /// alterar de acordo com a simulação
    
    //itapopo
    delete cmesh;
    delete gmesh;
    //itapopo
    
}//for hlevel
    
    //close the workspace and close the engine
  //  engEvalString(ep, "close;");
    //engClose(ep);
   
    //delete gmesh;
   // delete cmesh;
    //delete gmesh;
    
	std::cout << "FINISHED!" << std::endl;
	
	return 0;
}

///##############################################################
TPZCompMesh *CopyCompMesh(TPZCompMesh *cmesh, TPZGeoMesh *gmesh){
    
    
    TPZCompMesh * copycmesh = new TPZCompMesh(*cmesh);
    
    int matID = 1;
    TISSMaterial *material = new TISSMaterial(matID);
    copycmesh->InsertMaterialObject(material);

    copycmesh->SetReference(gmesh);
    copycmesh->LoadReferences();
    
    return copycmesh;


}
///##############################################################

bool XDif(TPZVec<REAL> &X1, TPZVec<REAL> &X2, REAL &Tol){
    
#ifdef PZDEBUG
    if(X1.NElements() != X2.NElements()) DebugStop();
#endif
    
    REAL DelX = 0.;
    
    for(int n = 0; n < X1.NElements(); n++){
        DelX += (X1[n]-X2[n]) * (X1[n]-X2[n]);
    }
    
    REAL error = std::sqrt(DelX);

    if(error <= Tol)
    {
        return true;
    }
    
    return false;
}
///##############################################################
bool FindNodeIndex( TPZVec<REAL> &NodeX, int64_t &NodeIndexOnSons, TPZVec<TPZGeoEl *> unrefinedSons ){
    
    /// NodeX is the target coordinates (node coordinates)
    
    std::set<int64_t> nodes;
    
    REAL Tol;
    ZeroTolerance(Tol);
    Tol *= 1.e3;
    
    for(int64_t i = 0; i < unrefinedSons.NElements();i++){
    
        for(int n = 0; n < unrefinedSons[i]->NNodes(); n++){
            
            int64_t NodeIndex = unrefinedSons[i]->NodeIndex(n);
            
            std::set<int64_t>::iterator it;
            it = nodes.find(NodeIndex);
            
            if( it == nodes.end() ){
            
                /// NodeIndex is not in the nodes set
                nodes.insert(NodeIndex);
                TPZVec<REAL> X(3,0);
                unrefinedSons[i]->NodePtr(n)->GetCoordinates(X);
                
                bool IsThisNode = XDif(NodeX, X, Tol);

                if(IsThisNode) {
                    NodeIndexOnSons = NodeIndex;
                    return true;
                }//if
    
            }//if
            
        }//for
    
    }//for
    
    NodeIndexOnSons = -1;
    
    return false;
}

///##############################################################
void TransferSolution(TPZFMatrix<STATE> &sol2, TPZFMatrix<STATE> &sol1,
                      int64_t &nodeIndex2, int64_t &nodeIndex1){
    
    /// 8 state variables
    const int nstate = 8;
    
    for(int var = 1; var <= nstate; var++ ){
        
        sol2(var-1 + nodeIndex2*nstate, 0) = sol1(var-1 + nodeIndex1*nstate, 0);
        
    }///for var
    
}


///##############################################################

void TransferSolution(TPZCompMesh *cmesh1, TPZCompMesh *cmesh2){
    
    /// set to avoid repeat each node
    std::set<int64_t> nodeset;
    std::set<int64_t>::iterator it;
    
    /// 1 solution
    const int nsol = 1;
    
    /// 8 state variables
    const int nstate = 8;
    
    /// number of nodes of the mesh2
    if(!cmesh2->Reference()) DebugStop();
    if(!cmesh1->Reference()) DebugStop();
    
    const int64_t nnodes = cmesh2->Reference()->NNodes();
    
    /// solution of the mesh2
    TPZFMatrix<STATE> sol2(nnodes*nstate, nsol, 0);
    
    /// solution of the mesh1
    TPZFMatrix<STATE> sol1 = cmesh1->Solution();
    
    REAL Tol;
    ZeroTolerance(Tol);
    
    std::cout << " ----- Starting TransferSolution ----- " << std::endl;
    
    for(int64_t i = 0; i < cmesh2->NElements(); i++){
        
        TPZCompEl *compel2 = cmesh2->Element(i);
        
        if( !compel2 ) continue;
        
#ifdef PZDEBUG
        if( !compel2->Reference() ) DebugStop();
#endif
        
        TPZGeoEl *geoel2 = compel2->Reference();
        
#ifdef PZDEBUG
        if( geoel2->HasSubElement() ) DebugStop();
#endif
        
        /// this ElIndex is the same in the gmesh1
        int64_t ElIndex;
        
        /// geoel2 has compel attached
        if( geoel2->LowestFather() ){
            ElIndex = geoel2->LowestFather()->Index();
        } else {
            ElIndex = geoel2->Index();//element i is not refined
        }
        
        /// geoel1 can be refined or not
        TPZGeoEl *geoel1 = cmesh1->Reference()->Element(ElIndex);
        
        /// geoel1 is refined
        if( geoel1->HasSubElement() ){
            
            /// the higher subelements. They have compel attached
            TPZVec<TPZGeoEl *> unrefinedSons;
            geoel1->GetHigherSubElements(unrefinedSons);
            
            for(int n = 0; n < geoel2->NNodes();n++){
                
                TPZVec<REAL> X2(3,0);
                geoel2->NodePtr(n)->GetCoordinates(X2);
                
                int64_t nodeIndex2 = geoel2->NodeIndex(n);
                int64_t nodeIndex1;
                
                it = nodeset.find(nodeIndex2);
                if( it != nodeset.end() ) continue; /// nodeIndex2 was already used
                nodeset.insert(nodeIndex2);
                
                if( FindNodeIndex( X2, nodeIndex1, unrefinedSons ) ){
                    
                    TransferSolution(sol2, sol1, nodeIndex2, nodeIndex1);
                    
                } else {

                    for(int s = 0; s < unrefinedSons.size(); s++){
                    
                        TPZVec<REAL> qsi(2,0);
                        TPZVec<STATE> sol(1,0);
                    
                        if(!unrefinedSons[s]->ComputeXInverse(X2, qsi, Tol) ) continue;
                    
                        #ifdef PZDEBUG
                        if(!unrefinedSons[s]->Reference() ) DebugStop();
                        #endif
                    
                        for(int var = 1; var <= nstate; var++ ){
                     
                            unrefinedSons[s]->Reference()->Solution(qsi, var, sol);
                            sol2(var-1 + nodeIndex2*nstate, 0) = sol[0];
                     
                        }///for var
                        
                    }/// for s
                    
                } /// else
                
            }/// for n
            
        } else { ///geoel1 is not refined (level = 0)
            
            for(int n = 0; n < geoel2->NNodes();n++){
                
                TPZVec<REAL> X2(3,0);
                geoel2->NodePtr(n)->GetCoordinates(X2);
                
                int64_t nodeIndex2 = geoel2->NodeIndex(n);
                int64_t nodeIndex1;
                
                it = nodeset.find(nodeIndex2);
                if( it != nodeset.end() ) continue; /// nodeIndex2 was already used
                nodeset.insert(nodeIndex2);
                
                TPZVec<TPZGeoEl *> OneEl(1,NULL);
                OneEl[0] = geoel1;
                
                if(FindNodeIndex(X2, nodeIndex1, OneEl)){
                
                    TransferSolution(sol2, sol1, nodeIndex2, nodeIndex1);
                
                } else {
                    
                    TPZVec<REAL> qsi(2,0);
                    TPZVec<STATE> sol(1,0);
                
                    if(geoel1->ComputeXInverse(X2, qsi, Tol)){
                    
                        for(int var = 1; var <= nstate; var++ ){
                            
                            geoel1->Reference()->Solution(qsi, var, sol);
                            sol2(var-1 + nodeIndex2*nstate, 0) = sol[0];
                            
                        }///for var
                    
                    } else {
                        
                        DebugStop();
                        
                    } /// else
                    
                
                } /// else
        
                
            }/// for n
            
            
        }/// if else
        
    }/// for i
    
    #ifdef PZDEBUG
    if(nnodes != nodeset.size()) DebugStop();
    #endif
    
    cmesh2->LoadSolution(sol2);
    
}

///##############################################################

void TransferSolution_bkp(TPZCompMesh *cmesh1, TPZCompMesh *cmesh2){
    
    /// 1 solution
    const int nsol = 1;
    
    /// 8 state variables
    const int nstate = 8;
    
    /// number of nodes of the mesh2
    if(!cmesh2->Reference()) DebugStop();
    if(!cmesh1->Reference()) DebugStop();
    
    const int64_t nnodes = cmesh2->Reference()->NNodes();
    
    /// solution of the mesh2
    TPZFMatrix<STATE> sol2(nnodes*nstate, nsol, 0);
    
    /// solution of the mesh1
    TPZFMatrix<STATE> sol1 = cmesh1->Solution();
    
    REAL Tol;
    ZeroTolerance(Tol);
    
    for(int64_t i = 0; i < cmesh2->NElements(); i++){
        
        TPZCompEl *compel2 = cmesh2->Element(i);
        
        if( !compel2 ) continue;
        
        #ifdef PZDEBUG
        if( !compel2->Reference() ) DebugStop();
        #endif
        
        TPZGeoEl *geoel2 = compel2->Reference();
        
        #ifdef PZDEBUG
        if( geoel2->HasSubElement() ) DebugStop();
        #endif
        
        /// this ElIndex is the same in the gmesh1
        int64_t ElIndex;
        
        /// geoel2 has compel attached
        if( geoel2->LowestFather() ){
            ElIndex = geoel2->LowestFather()->Index();
        } else {
            ElIndex = geoel2->Index();//element i is not refined
        }
        
        /// geoel1 can be refined or not
        TPZGeoEl *geoel1 = cmesh1->Reference()->Element(ElIndex);
        
        /// geoel1 is refined
        if( geoel1->HasSubElement() ){
            
            /// the higher subelements. They have compel attached
            TPZVec<TPZGeoEl *> unrefinedSons;
            geoel1->GetHigherSubElements(unrefinedSons);
            
            for(int n = 0; n < geoel2->NNodes();n++){
                
                TPZVec<REAL> X2(3,0);
                geoel2->NodePtr(n)->GetCoordinates(X2);
                const int64_t nodeIndex2 = geoel2->NodeIndex(n);
                TPZVec<REAL> qsi(2,0);
                TPZVec<STATE> sol(1,0);
                
                for(int s = 0; s < unrefinedSons.size(); s++){
                
                    //if(!unrefinedSons[s]->ComputeXInverse(X2, qsi, Tol) ) continue;
                
                    #ifdef PZDEBUG
                    if(!unrefinedSons[s]->Reference() ) DebugStop();
                    #endif
                    
                    /// testando
                    int n1;
                    for(n1 = 0; n1 < unrefinedSons[s]->NNodes(); n1++){
                        
                        TPZVec<REAL> X1(3,0);
                        unrefinedSons[s]->NodePtr(n1)->GetCoordinates(X1);
                        
                        bool IsThisNode = XDif(X1, X2, Tol);
                        
                        if (IsThisNode) {
                            
                            int64_t nodeIndex1 = unrefinedSons[s]->NodeIndex(n1);
                            
                            for(int var = 1; var <= nstate; var++ ){
                            
                                sol2(var-1 + nodeIndex2*nstate, 0) = sol1(var-1 + nodeIndex1*nstate, 0);
                                
                            }///for var
                            
                            break;
                        } //if
                    } //for
                    /*
                    if( n1==unrefinedSons[s]->NNodes() ) {
                        
                        for(int var = 1; var <= nstate; var++ ){
                            
                            unrefinedSons[s]->Reference()->Solution(qsi, var, sol);
                            sol2(var-1 + nodeIndex2*nstate, 0) = sol[0];
                            
                        }///for var
                        
                    }*/
                    
                    /// testando
                    
                    /** era original
                    
                     for(int var = 1; var <= nstate; var++ ){
                    
                        unrefinedSons[s]->Reference()->Solution(qsi, var, sol);
                        sol2(var-1 + nodeIndex2*nstate, 0) = sol[0];
                        
                    }///for var
                     
                     era original
                     */
                
                }/// for s
                
            }/// for n
            
        } else { ///geoel1 is not refined (level = 0)
            
            for(int n = 0; n < geoel2->NNodes();n++){
                
                TPZVec<REAL> X2(3,0);
                geoel2->NodePtr(n)->GetCoordinates(X2);
                const int64_t nodeIndex2 = geoel2->NodeIndex(n);
                
                TPZVec<REAL> qsi(2,0);
                TPZVec<STATE> sol(1,0);
                
                //if(geoel1->ComputeXInverse(X2, qsi, Tol)){
                    
                /// testando
                int n1;
                for(n1 = 0; n1 < geoel1->NNodes(); n1++){
                        
                    TPZVec<REAL> X1(3,0);
                    geoel1->NodePtr(n1)->GetCoordinates(X1);
                    bool IsThisNode = XDif(X1, X2, Tol);
                        
                    if (IsThisNode) {
                            
                        int64_t NodeIndex1 = geoel1->NodeIndex(n1);
                            
                        for(int var = 1; var <= nstate; var++ ){
                                
                            sol2(var-1 + nodeIndex2*nstate, 0) = sol1(var-1 + NodeIndex1*nstate, 0);
                                
                        }///for var
                            
                        break;
                    } //if
                } //for
                    
                if( n1==geoel1->NNodes() ) {
                        
                    for(int var = 1; var <= nstate; var++ ){
                        
                        bool IsInverse = geoel1->ComputeXInverse(X2, qsi, Tol);
                        
                        geoel1->Reference()->Solution(qsi, var, sol);
                        sol2(var-1 + nodeIndex2*nstate, 0) = sol[0];
                    
                        #ifdef PZDEBUG
                        if(!IsInverse) {
                            DebugStop();
                        }
                        #endif
                        
                    }///for var
                        
                }
                    
                    /// testando
                    
                    
                    /** era original
                    for(int var = 1; var <= nstate; var++ ){
                        
                        geoel1->Reference()->Solution(qsi, var, sol);
                        sol2(var-1 + nodeIndex2*nstate, 0) = sol[0];
                
                    }///for var
                     era original
                     */
                     
               // } else {/// if
                 
                 //   DebugStop();
                
               // }///else
        
            }/// for n
            
            
        }/// if else
        
    }/// for i
    
    cmesh2->LoadSolution(sol2);
  
}

///##############################################################

void PrintMesh(const int &step,
                TPZGeoMesh *gmesh){
    
    std::string strGeoMesh = toStr(WORKPATH) + "debug/geoMesh" + toStr(step) + ".txt";
    std::ofstream fileGeoMesh(strGeoMesh.c_str());
    gmesh->Print(fileGeoMesh);
    
    std::string strCompMesh = toStr(WORKPATH) + "debug/compMesh" + toStr(step) + ".txt";
    std::ofstream fileCompMesh(strCompMesh.c_str());
    if(gmesh->Reference()){
        gmesh->Reference()->Print(fileCompMesh);
    }
    
    std::string strGeoMeshVTK = toStr(WORKPATH) + "debug/geoMesh" + toStr(step) + ".vtk";
    std::ofstream fileGeoMeshVTK(strGeoMeshVTK.c_str());
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, fileGeoMeshVTK,true);
    
}

///##############################################################

std::string GetRepository(const int &level){
    
    std::string ModelName;
   
    if(IsViscous){
        
        ModelName = "Level" + toStr(level) + "_viscous";
        
    } else {
        
        ModelName = "Level" + toStr(level) + "_coulomb";
        
    }
    
    std::string MainFolder;
    
    if(RUN_3D == 0){
        
        MainFolder = "2D_SSA";
        
    }else if(RUN_3D==1){
        
        MainFolder = "3D_HO";
        
    }else {
        
        DebugStop();
    }
    
    std::string Repository = MainFolder + "/Models_" + ModelName;
    
    return Repository;
    
}


///##############################################################
void SetFileNames(const int &h,
                  const int &step,
                  bool ChangeSolutionFile,
                  std::string &strSolutionFile,
                  std::string &strMeshFile,
                  std::string &strDataFile,
                  std::string &strSolutionMatlab,
                  std::string &strMeshMatlab,
                  std::string &strDataMatlab){
    
    if(!ChangeSolutionFile){
        // the file name follows the step number
        strSolutionFile = GetRepository(h) + "/sol/" + GetExperimentPrefix(STEPS_SIM) + "_solution" + toStr(step) + ".txt";
        
        strMeshFile = GetRepository(h) + "/mesh/" + GetExperimentPrefix(STEPS_SIM) + "_mesh" + toStr(step) + ".txt";
        
        strDataFile = GetRepository(h) + "/data/" + GetExperimentPrefix(STEPS_SIM) + "_data" + toStr(step) + ".txt";
        
    }else{
        // in Restart, it is necessary to print the restart file
        
        //THIS NEEDS TO BE THE FIRST TIME WHEN RESTARTING.
        
        const int sim_step = 11; //Ice2r itapopo
        
        strSolutionFile = GetRepository(h) + "/sol/" + GetExperimentPrefix(sim_step) + "_solution" + toStr(step) + "restart.txt";
        
        strMeshFile = GetRepository(h) + "/mesh/" + GetExperimentPrefix(sim_step) + "_mesh" + toStr(step) + ".txt";
        
        strDataFile = GetRepository(h) + "/data/" + GetExperimentPrefix(sim_step) + "_data" + toStr(step) + ".txt";
        
    }
    
    strSolutionMatlab = "solutionfile = '" + strSolutionFile + "';";
    
    //strMeshFile = GetRepository(h) + "/mesh/" + GetExperimentPrefix() + "_mesh" + toStr(step) + ".txt";
    strMeshMatlab = "meshfile = '" + strMeshFile + "';";
    
    //strDataFile = GetRepository(h) + "/data/" + GetExperimentPrefix() + "_data" + toStr(step) + ".txt";
    strDataMatlab = "datafile = '" + strDataFile + "';";

    /** bkp antes do respositótio
    strSolutionFile = "sol" + toStr(h) + "/solution" + toStr(step) + ".txt";
    strSolutionMatlab = "solutionfile = '" + strSolutionFile + "';";
    
    strMeshFile = "mesh" + toStr(h) + "/mesh" + toStr(step) + ".txt";
    strMeshMatlab = "meshfile = '" + strMeshFile + "';";
    
    strDataFile = "data" + toStr(h) + "/data" + toStr(step) + ".txt";
    strDataMatlab = "datafile = '" + strDataFile + "';";
    */
}


///##############################################################

void GetVariableNames(TPZStack<std::string> &scalnames, TPZStack<std::string> &vecnames){
    
    scalnames.clear();
    vecnames.clear();
    
    //setando os nomes das variáveis
    scalnames.Push("Surface");
    scalnames.Push("Base");
    scalnames.Push("Bed");
    scalnames.Push("Pressure");
    scalnames.Push("Temperature");
    scalnames.Push("Vx");
    scalnames.Push("Vy");
    scalnames.Push("MaskLevelSet");
    
}

///##############################################################
void SetElementsToRefine(TPZCompMesh *cmesh,
                         /*std::vector<int64_t>*//*std::set<int64_t>*/std::map<int64_t, int> &ElementIndex,
                         std::vector<TPZVec<REAL> > &GLvec,
                         REAL &alpha){
    
    ElementIndex.clear();///index of the element in the geometric mesh index
    
    const int64_t nelem = cmesh->NElements();
    const int var = 8;//masklevelset
    
    const REAL MaxLevelSet = 200.;
    
    const bool IsUniform = IS_UNIFORM;
    
    const bool IsUsingLevelSetValue = false;
    
    const bool IsRadiusValue = !IsUniform;
    
    const REAL MaxDistance = MAX_GL_DISTANCE*alpha;
    
    //std::vector<TPZVec<REAL> > GLvec;
    GLvec.clear();
    
    for(int64_t i = 0; i < nelem; i++){
        
        TPZCompEl *compEl = cmesh->Element(i);
        
        if(!compEl) continue;
        
        TPZGeoEl *geoel = compEl->Reference();
        
        if( geoel->HasSubElement() ) continue;
        
        int64_t FatherIndex;
        
        ///uniform refinement
        if(IsUniform){
            
            if( geoel->LowestFather() ) {
                FatherIndex = geoel->LowestFather()->Index();
            } else {
                FatherIndex = geoel->Index();
            }
            //ElementIndex.insert(std::make_pair<int64_t,int>(FatherIndex,1));//insert(FatherIndex);//push_back(i);
            continue;
        }
   
        TPZVec<REAL> qsi(3,0);
        TPZVec<STATE> sol1(1,0), sol2(1,0), sol3(1,0) ;
        
        /// qsi = (0, 0)
        compEl->Solution(qsi, var, sol1);
        
        /// qsi = (1, 0)
        qsi[0] = 1;
        compEl->Solution(qsi, var, sol2);
        
        /// qsi = (0, 1)
        qsi[0] = 0;
        qsi[1] = 1;
        compEl->Solution(qsi, var, sol3);
        
        ///refining using the level set value
        if(IsUsingLevelSetValue){
            
            if( sol1[0]*sol2[0] > 0. && sol1[0]*sol3[0] > 0. ){
            
                STATE mean = ( sol1[0] + sol2[0] + sol3[0] ) / 3.;
                if( std::fabs(mean) > MaxLevelSet ) continue; //it is not necessary to refine
            
            }
        
            /// it is necessary to refine
            if( geoel->LowestFather() ) {
                FatherIndex = geoel->LowestFather()->Index();
            } else {
                FatherIndex = geoel->Index();
            }
        
            //ElementIndex.insert(std::make_pair<int64_t,int>(FatherIndex,1));//insert(FatherIndex);//push_back(i);
            
        } /// if IsUsingLevelSetValue
        
        if(IsRadiusValue){
            
            if( sol1[0]*sol2[0] < 0. || sol1[0]*sol3[0] < 0. ){
            
                
                const int side = 6;
                TPZVec<REAL> qsi(2,0.);
                TPZVec<REAL> X(3,0.);
                geoel->CenterPoint(side, qsi);
                geoel->X(qsi, X);
                
                GLvec.push_back(X);
            
            }
            
            
        }/// if IsRadiusValue
        
        
    }///for i / nelem
    
    if(IsRadiusValue){
        
        for(int64_t i = 0; i < nelem; i++){
            
            TPZCompEl *compEl = cmesh->Element(i);
            
            if(!compEl) continue;
            
            TPZGeoEl *geoel = compEl->Reference();
            
            if( geoel->HasSubElement() ) continue;
            
            int64_t FatherIndex;
            
            TPZGeoEl *LowestFather;
            
            if( geoel->LowestFather() ) {
                LowestFather = geoel->LowestFather();
            } else {
                LowestFather = geoel;
            }
            
            FatherIndex = LowestFather->Index();
            
            const int side = 6;
            TPZVec<REAL> qsi(2,0.);
            TPZVec<REAL> centerPoint(3,0.);
            LowestFather->CenterPoint(side, qsi);
            LowestFather->X(qsi, centerPoint);
            
            REAL distance = MaxDistance;
            
            for (int64_t j = 0; j < GLvec.size(); j++) {
            
                REAL value = ( GLvec[j][0] - centerPoint[0] ) * ( GLvec[j][0] - centerPoint[0] ); // (x2-x1)^2
                value += ( GLvec[j][1] - centerPoint[1] ) * ( GLvec[j][1] - centerPoint[1] );// (y2-y1)^2
                value = std::sqrt(value); ///Radius
                
                if(value < distance){
                    distance = value; //finding the min distance to the Grounding line
                } ///if
                
            } ///for j / GLvec.size()
            
            //if(distance < MaxDistance/8){ /// it is necessary to refine
                
                /*if( geoel->LowestFather() ) {
                    FatherIndex = geoel->LowestFather()->Index();
                } else {
                    FatherIndex = geoel->Index();
                }*/
         /*
                ElementIndex.insert(std::make_pair<int64_t,int>(FatherIndex,2));//insert(FatherIndex);//push_back(i);
                continue;
          
            }/// if
            */
            if(distance < MaxDistance){ /// it is necessary to refine
                /*
                if( geoel->LowestFather() ) {
                    FatherIndex = geoel->LowestFather()->Index();
                } else {
                    FatherIndex = geoel->Index();
                }
                */
                //ElementIndex.insert(std::make_pair<int64_t,int>(FatherIndex,1));//insert(FatherIndex);//push_back(i);
                continue;
                
            }/// if
            
            
        }/// for i / nelem
        
    } /// if IsRadiusValue
    
    
}

///##############################################################

void LoadSolution(std::string &SolutionName, TPZCompMesh *cmesh){
    
    std::string FullName = toStr(WORKPATH) + SolutionName;
    
    std::ifstream SolutionFile(FullName.c_str());
    
    bool IsOk = false;
    
    std::vector<double> surface;
    std::vector<double> base;
    std::vector<double> bed;
    std::vector<double> pressure;
    std::vector<double> temperature;
    std::vector<double> vx;
    std::vector<double> vy;
    std::vector<double> masklevelset;
    
    // reading from file
    IsOk = LoadSolution(SolutionFile,
                 surface,
                 base,
                 bed,
                 pressure,
                 temperature,
                 vx,
                 vy,
                 masklevelset);
    
    if(!IsOk) DebugStop();
    
    /// 1 solution
    const int nsol = 1;
    
    /// 8 state variables
    const int nstate = 8;
    
    /// number of nodes
    const int nnodes = surface.size();
    
    TPZFMatrix<STATE> sol(nnodes*nstate, nsol, 0);
    
    for(int i = 0; i < nnodes; i ++){
        
        sol(0+i*nstate,0) = surface[i];
        sol(1+i*nstate,0) = base[i];
        sol(2+i*nstate,0) = bed[i];
        sol(3+i*nstate,0) = pressure[i];
        sol(4+i*nstate,0) = temperature[i];
        sol(5+i*nstate,0) = vx[i];
        sol(6+i*nstate,0) = vy[i];
        sol(7+i*nstate,0) = masklevelset[i];
        
    }
    
    cmesh->LoadSolution(sol);
    
}

///##############################################################
TPZCompMesh *CreateCompMesh(TPZGeoMesh *gmesh){
    
    /// domain dimention (2D)
    const int dim = 2;
    
    /// polynomial order
    const int pOrder = 1;
    
    /// material ID
    const int matID = 1;
    
    /// creating any material
    TISSMaterial *material = new TISSMaterial(matID);

    /// creating the computational mesh
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    cmesh->InsertMaterialObject(material);
    cmesh->SetAllCreateFunctionsContinuous();
    
    /// adjust the data structure
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    Permute2(cmesh);
    
    return cmesh;
    
    
}

///##############################################################

void GeoToConnect(TPZCompMesh * cmesh, TPZVec<int64_t> &geotoconnect){
    
    TPZGeoMesh *gmesh = cmesh->Reference();
    geotoconnect.resize(gmesh->NNodes());
    geotoconnect.Fill(-1);
    
    int64_t nel=cmesh->ElementVec().NElements();
    
    for (int64_t jel=0; jel<nel; jel++) {
        
        TPZCompEl *compEl = cmesh->ElementVec()[jel];
        if(!compEl) continue;
        
        TPZGeoEl *geoEl = compEl->Reference();
        if(!geoEl) DebugStop();
        
        const int64_t nnodes = geoEl->NCornerNodes();
        
        for(int64_t node = 0; node < nnodes; node++){
            
            int64_t nodeindex = geoEl->NodeIndex(node);
            int64_t connectindex = compEl->ConnectIndex(node);
            
            if (geotoconnect[nodeindex] != -1 && geotoconnect[nodeindex] != connectindex) {
                DebugStop();
            }
            
            geotoconnect[nodeindex] = connectindex;
            
        }///for node
        
    }/// for jel
    
}

///##############################################################


void Permute2(TPZCompMesh * cmesh){

    BuildPermuteVector(cmesh);
    
}

///##############################################################


void BuildPermuteVector(TPZCompMesh *cmesh)
{
#ifdef PZDEBUG
    TPZGeoMesh *gmesh = cmesh->Reference();
    if (gmesh->Reference() != cmesh) {
        DebugStop();
    }
#endif
    cmesh->LoadReferences();
    
    TPZVec<int64_t> geotoconnect;

    GeoToConnect(cmesh, geotoconnect);
    
    //must be the internal number of connects
    const int64_t nblocks = cmesh->NIndependentConnects();//cmesh->NConnects();
    
#ifdef PZDEBUG
    if (cmesh->Block().NBlocks() != nblocks) {
        DebugStop();
    }
#endif
    
    TPZVec<int64_t> permutegather(nblocks,-1);
    TPZVec<int64_t> used(nblocks,-1);
    
    // fill in initial part of permute
    int64_t ngeo = geotoconnect.size();
    for (int64_t inode=0; inode < ngeo; inode++) {
        int64_t connectindex = geotoconnect[inode];
        int64_t seqnum = cmesh->ConnectVec()[connectindex].SequenceNumber();
        // it verifies if seqnum is >= 0
#ifdef PZDEBUG
        if(seqnum < 0){
            DebugStop();
        }
#endif
        permutegather[seqnum] = inode;
        used[seqnum] = 1;
    }
    
    // fill the other part of permute
    int64_t usedcounter = 0;
    for (int64_t index = ngeo; index < nblocks; index++) {
        
        while (used[usedcounter] == 1)  {
            //permutegather[index] = usedcounter;//using a seqnum nerver used before
            //used[usedcounter] = 1;
            usedcounter++;
        }///while
    
        permutegather[usedcounter] = index;//using a seqnum nerver used before
        used[usedcounter] = 1;
        
        //if (usedcounter < nconnects-1) {
         //   usedcounter++;
        //};
    }
    
#ifdef PZDEBUG

    for(int64_t i = 0; i < nblocks;i++){
        if(used[i] == -1 || permutegather[i] == -1){
            DebugStop();
        }
    }
    
    // verifica se permutegahter eh permutacao
    std::set<int64_t> permset;
    for (int64_t index = 0; index<nblocks; index++) {
        permset.insert(permutegather[index]);
    }
    if (permset.size() != nblocks) {
        DebugStop();
    }
    {
        std::ofstream out("../permutedata.txt");
        for (int64_t ig=0; ig<ngeo; ig++) {
            int64_t connectindex = geotoconnect[ig];
            TPZConnect &c = cmesh->ConnectVec()[connectindex];
            out << "geonode " << ig << " cindex " << connectindex << " seqnum " << c.SequenceNumber() << " permute " << permutegather[ig] << std::endl;
        }
    }
#endif
    
    cmesh->Permute(permutegather);
    
#ifdef PZDEBUG

    {
        std::ofstream out("../permutedataAfter.txt");
        for (int64_t ig=0; ig<ngeo; ig++) {
            int64_t connectindex = geotoconnect[ig];
            TPZConnect &c = cmesh->ConnectVec()[connectindex];
            out << "geonode " << ig << " cindex " << connectindex << " seqnum " << c.SequenceNumber() << " permute " << permutegather[ig] << std::endl;
        }
    }
    CheckSequenceNumber(cmesh);
#endif
    
}

///##############################################################
void Permute(TPZCompMesh * cmesh){
    
    TPZVec<int64_t> permute;
    int64_t numinternalconnects = cmesh->NIndependentConnects();
    permute.Resize(numinternalconnects,0);
    
    int64_t nel=cmesh->ElementVec().NElements();
    
    for (int64_t jel=0; jel<nel; jel++) {
        
        TPZCompEl *compEl = cmesh->ElementVec()[jel];
        if(!compEl) continue;
        
        TPZGeoEl *geoEl = compEl->Reference();
        if(!geoEl) continue;
        
        const int64_t nnodes = geoEl->NNodes();
        
        for(int64_t node = 0; node < nnodes; node++){
            
            for (int64_t ip=0; ip<permute.NElements(); ip++) {
                permute[ip]=ip;
            }
            
            int64_t nodeindex = geoEl->NodeIndex(node);
            int64_t seqnum = compEl->Connect(node).SequenceNumber();
            
            int64_t v1 = permute[seqnum];
            permute[nodeindex] = v1;
            permute[seqnum] = nodeindex;
            
            cmesh->Permute(permute);
            
        }///for node
        
    }/// for jel
    
}

///##############################################################

std::string toStr(int intVal){
    
    std::stringstream ss;
    ss << intVal;
    std::string strVal = ss.str();
    
    return strVal;
}

///##############################################################

std::string toStr(const char *charVal){
    
    std::stringstream ss;
    ss << charVal;
    std::string strVal = ss.str();
    
    return strVal;
}

///##############################################################

bool LoadSolution(std::ifstream &SolutionFile, double **vx, double **vy, double **masklevelset){
    
    std::vector<double> surface;
    std::vector<double> base;
    std::vector<double> bed;
    std::vector<double> pressure;
    std::vector<double> temperature;
    std::vector<double> mvx;
    std::vector<double> mvy;
    std::vector<double> mmasklevelset;
    
    double *vxptr;
    double *vyptr;
    double *masklevelsetptr;
    
    // reading from file
    bool IsOk = LoadSolution(SolutionFile,surface,base,bed,pressure,temperature,mvx,mvy,mmasklevelset);
    if(!IsOk) DebugStop();
    
    int64_t nvertices = mmasklevelset.size();
    masklevelsetptr = new double[nvertices];
    for(int64_t i=0;i<nvertices;i++) masklevelsetptr[i]=mmasklevelset[i];
    
    vxptr = NULL;
    vyptr = NULL;
    
    *masklevelset = masklevelsetptr;
    *vx = vxptr;
    *vy = vyptr;
    
    return IsOk;
}

///##############################################################

bool LoadSolution(std::ifstream &SolutionFile,
                 std::vector<double> &surface,
                 std::vector<double> &base,
                 std::vector<double> &bed,
                 std::vector<double> &pressure,
                 std::vector<double> &temperature,
                 std::vector<double> &vx,
                 std::vector<double> &vy,
                 std::vector<double> &masklevelset){
    
    if( !SolutionFile.is_open() ) return false;
    
    /// going to the begining of the file
    SolutionFile.seekg(0);
    
    /// reading number of nodes
    int64_t nnodes;
    
    SolutionFile >> nnodes;
    
    surface.clear();
    base.clear();
    bed.clear();
    pressure.clear();
    temperature.clear();
    vx.clear();
    vy.clear();
    masklevelset.clear();
    
    /// reading surface
    for(int64_t i = 0; i < nnodes; i++){
        
        double value;
        
        SolutionFile >> value;
   
        surface.push_back(value);
        
    }
    
    /// reading base
    for(int64_t i = 0; i < nnodes; i++){
        
        double value;
        
        SolutionFile >> value;
        
        base.push_back(value);
        
    }
    
    /// reading bed
    for(int64_t i = 0; i < nnodes; i++){
        
        double value;
        
        SolutionFile >> value;
        
        bed.push_back(value);
        
    }
    
    /// reading pressure
    for(int64_t i = 0; i < nnodes; i++){
        
        double value;
        
        SolutionFile >> value;
        
        pressure.push_back(value);
        
    }
    
    /// reading temperature
    for(int64_t i = 0; i < nnodes; i++){
        
        double value;
        
        SolutionFile >> value;
        
        temperature.push_back(value);
        
    }
    
    /// reading vx
    for(int64_t i = 0; i < nnodes; i++){
        
        double value;
        
        SolutionFile >> value;
        
        vx.push_back(value);
        
    }
    
    /// reading vy
    for(int64_t i = 0; i < nnodes; i++){
        
        double value;
        
        SolutionFile >> value;
        
        vy.push_back(value);
        
    }
    
    /// reading masklevelset
    for(int64_t i = 0; i < nnodes; i++){
        
        double value;
        
        SolutionFile >> value;
        
        masklevelset.push_back(value);
        
    }
    
    return true;
    
}

///##############################################################

void SetRefPatterns(){
    
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    
    //gRefDBase.InitializeRefPatterns();
    std::string filepath = REFPATTERNDIR;
    std::string filename1 = filepath + "/2D_Triang_Rib_3.rpt";
    std::string filename2 = filepath + "/2D_Triang_Rib_4.rpt";
    std::string filename3 = filepath + "/2D_Triang_Rib_OnlyTriang_Side_3_4.rpt";
    std::string filename4 = filepath + "/2D_Triang_Rib_OnlyTriang_Side_3_4_permuted.rpt";
    std::string filename5 = filepath + "/2D_Triang_Rib_OnlyTriang_Side_3_5.rpt";
    std::string filename6 = filepath + "/2D_Triang_Rib_OnlyTriang_Side_3_5_permuted.rpt";
    std::string filename7 = filepath + "/2D_Triang_Rib_5.rpt";
    
    TPZAutoPointer<TPZRefPattern> refpat1 = new TPZRefPattern(filename1);
    TPZAutoPointer<TPZRefPattern> refpat2 = new TPZRefPattern(filename2);
    TPZAutoPointer<TPZRefPattern> refpat3 = new TPZRefPattern(filename3);
    TPZAutoPointer<TPZRefPattern> refpat4 = new TPZRefPattern(filename4);
    TPZAutoPointer<TPZRefPattern> refpat5 = new TPZRefPattern(filename5);
    TPZAutoPointer<TPZRefPattern> refpat6 = new TPZRefPattern(filename6);
    TPZAutoPointer<TPZRefPattern> refpat7 = new TPZRefPattern(filename7);
    
    if(!gRefDBase.FindRefPattern(refpat1))
    {
        gRefDBase.InsertRefPattern(refpat1);
    }
    //refpat1->InsertPermuted();
    
    if(!gRefDBase.FindRefPattern(refpat2))
    {
        gRefDBase.InsertRefPattern(refpat2);
    }
    //refpat2->InsertPermuted();
    
    if(!gRefDBase.FindRefPattern(refpat3))
    {
        gRefDBase.InsertRefPattern(refpat3);
    }
    //refpat3->InsertPermuted();
    
    if(!gRefDBase.FindRefPattern(refpat4))
    {
        gRefDBase.InsertRefPattern(refpat4);
    }
    //refpat4->InsertPermuted();
    
    if(!gRefDBase.FindRefPattern(refpat5))
    {
        gRefDBase.InsertRefPattern(refpat5);
    }
    //refpat5->InsertPermuted();
    
    if(!gRefDBase.FindRefPattern(refpat6))
    {
        gRefDBase.InsertRefPattern(refpat6);
    }
    //refpat6->InsertPermuted();
    
    if(!gRefDBase.FindRefPattern(refpat7))
    {
        gRefDBase.InsertRefPattern(refpat7);
    }
    //refpat7->InsertPermuted();
    
    
}

///##############################################################

TPZGeoMesh *CreateGMesh(std::string &MeshName){
    
    std::string FullName = toStr(WORKPATH) + MeshName;
    
    std::ifstream MeshFile(FullName.c_str());
    
    bool IsOk = false;
    
    std::vector<double> x, y;
    std::vector<ElemNode> ElemNodeVec;
    std::vector<Segment> SegmentVec;
    std::vector<int64_t> ElemIdVec;
    
    /// reading mesh
    IsOk = ReadMesh(MeshFile, x, y, ElemNodeVec, SegmentVec, ElemIdVec);
    if(!IsOk) DebugStop();
    
    /// creating gmesh
    TPZGeoMesh *gmesh = CreateGMesh(x, y, ElemNodeVec, SegmentVec);
    
    return gmesh;
    
}


///##############################################################
bool ReadMesh(std::ifstream &MeshFile,int64_t &nvertices,int64_t &nelements,int64_t &nsegments,double **x,double **y,double **z,int64_t ***elements,int64_t ***segments){
    
    std::vector<double> mx,my;
    std::vector<ElemNode> ElemNodeVec;
    std::vector<Segment> SegmentVec;
    std::vector<int64_t> ElemIdVec;
    
    bool IsOk = ReadMesh(MeshFile, mx, my, ElemNodeVec, SegmentVec, ElemIdVec);
    if(!IsOk) DebugStop();
    
    double * xptr;
    double * yptr;
    double * zptr;
    int64_t **elementsptr;
    int64_t **segmentsptr;
    
    nvertices = mx.size();
    nelements = ElemNodeVec.size();
    nsegments = SegmentVec.size();
    
    xptr = new double[nvertices];
    yptr = new double[nvertices];
    zptr = new double[nvertices];
    
    for(int64_t i = 0; i < nvertices; i++){
        xptr[i] = mx[i];
        yptr[i] = my[i];
        zptr[i] = 0;
    }
    
    elementsptr = new int64_t*[nelements];
    for(int64_t i = 0; i < nelements; i++){
        elementsptr[i] = new int64_t[3];
        for(int j = 0; j < 3; j++){
            elementsptr[i][j] = ElemNodeVec[i].fNodeId[j] -1; //0 -1 e para inicializar de 0
        }
    }
    
    segmentsptr = new int64_t*[nsegments];
    for(int64_t i = 0; i < nsegments; i++){
        segmentsptr[i] = new int64_t[2];
        for(int j = 0; j < 2; j++){
            segmentsptr[i][j] = SegmentVec[i].fId[j] -1; //0 -1 e para inicializar de 0
        }
    }
    
    *x          = xptr;
    *y          = yptr;
    *z          = zptr;
    *elements   = elementsptr;
    *segments   = segmentsptr;
    
    return IsOk;
    
}
///##############################################################
bool ReadMesh(std::ifstream &MeshFile,
              std::vector<double> &x,
              std::vector<double> &y,
              std::vector<ElemNode> &ElemNodeVec,
              std::vector<Segment> &SegmentVec,
              std::vector<int64_t> &ElemIdVec){
    
    if( !MeshFile.is_open() ) return false;
    
/// reading nodes
    int64_t nnodes;
    
    MeshFile >> nnodes;
    
    x.clear();
    y.clear();
    
    for(int64_t i = 0; i < nnodes; i++){
        
        double xvalue, yvalue;
        
        MeshFile >> xvalue;
        MeshFile >> yvalue;
        
        x.push_back(xvalue);
        y.push_back(yvalue);
        
    }

///reading elements
    int64_t nelements;
    
    MeshFile >> nelements;
    
    ElemNodeVec.clear();
    
    for(int64_t i = 0; i < nelements; i++){
        
        int64_t n1, n2, n3;
        
        MeshFile >> n1;
        MeshFile >> n2;
        MeshFile >> n3;
        
        ElemNode Nodes;
        
        Nodes.fNodeId.push_back(n1);
        Nodes.fNodeId.push_back(n2);
        Nodes.fNodeId.push_back(n3);
        
        ElemNodeVec.push_back(Nodes);
    }
    
///reading segments
    int64_t nsegments;
    
    MeshFile >> nsegments;
    
    SegmentVec.clear();
    
    for(int64_t i = 0; i < nsegments; i++){
        
        int64_t nId1, nId2, eId;
        
        MeshFile >> nId1;
        MeshFile >> nId2;
        MeshFile >> eId;
        
        Segment OneSegment;
        
        OneSegment.fId.push_back(nId1);
        OneSegment.fId.push_back(nId2);
        OneSegment.fId.push_back(eId);
        
        SegmentVec.push_back(OneSegment);
    }
    
///reading element's ID to refine
 /*   int nids;
    
    MeshFile >> nids;
    
    ElemIdVec.clear();
    
    for(int i = 0; i < nids; i++){
        
        int64_t ID;
        
        MeshFile >> ID;
        
        ElemIdVec.push_back(ID);
        
    }
    
    */
    return true;
    
}


///##############################################################
void FlagElements(TPZGeoMesh *gmesh, std::vector<std::pair<int,int> > &FlagVec, int step, int dstep){

    FlagVec.clear();
    
    //region
    double dS1 = 200000.; //200 km
    double dS2 = 50000.; //50 km
    
    double Si = 300000.; //300 km
    double Sf = 700000.; //700 km
    double vel = (Sf-Si)/dstep;
    
    double S1f = vel*step + Si;
    double S1i = S1f - dS1;
    
    double S2f = vel*step + Si - (dS1-dS2)*0.5;
    double S2i = S2f - dS2;
    
    for(int64_t i = 0; i < gmesh->NElements(); i++){
        //apenas elementos triangulares
        if(gmesh->Element(i)->MaterialId() != 1) continue;
        
        TPZManVector<REAL,3> qsi(gmesh->Element(i)->Dimension()), xCenter(3,0.);
        gmesh->Element(i)->CenterPoint(gmesh->Element(i)->NSides()-1, qsi);
        gmesh->Element(i)->X(qsi,xCenter);
        
        double Xtria = xCenter[0];
        
        
        if( Xtria > S1i && Xtria < S1f){
            
            if( Xtria > S2i && Xtria < S2f){
                
                FlagVec.push_back(std::make_pair<int,int>(i+1, 2));
            
            } else {
                
                FlagVec.push_back(std::make_pair<int,int>(i+1, 1));
            }
            
        }
        
    }
    
    
    
    
}

///##############################################################
TPZGeoMesh *CreateGMesh(std::vector<double> &x,
                        std::vector<double> &y,
                        std::vector<ElemNode> &ElemNodeVec,
                        std::vector<Segment> &SegmentVec)
{
	TPZGeoMesh * gmesh = new TPZGeoMesh;
	
	gmesh->NodeVec().Resize( x.size() );

    //definindo os nós da malha
	for(int64_t i = 0 ; i < x.size(); i++){
        
        TPZManVector<REAL,3> coord(3,0.);
        
        coord[0]= x[i];
        coord[1]= y[i];
		
        gmesh->NodeVec()[i].SetCoord(coord);
		gmesh->NodeVec()[i].SetNodeId(i);
	}
	
	//materials ID
    int64_t id;
    const int matTria = 1;
    const int matBoundary = 2;
    
    // Criando Elementos Triangulares
    TPZManVector<int64_t,3> tria(3,0.);
    
	for(int64_t iel = 0; iel < ElemNodeVec.size(); iel++){

        //tira-se "-1" dos IDs pois o MatLab começa a numeração dos nós em "1" e não em "0"
        tria[0] = ElemNodeVec[iel].fNodeId[0] - 1;
        tria[1] = ElemNodeVec[iel].fNodeId[1] - 1;
        tria[2] = ElemNodeVec[iel].fNodeId[2] - 1;

        //definindo reftyoe = 0, ou seja, será utilizado padrão uniforme por default (mais rápido)
        const int reftype = 1;
        gmesh->CreateGeoElement(ETriangle, tria, matTria, id,reftype);//cria elemento triangular
        
        gmesh->ElementVec()[id]->SetId(iel);
        
	}
    
    //Criando os elementos unidimensionais que compõem todo o contorno
    TPZManVector<int64_t,2> boundary(2,0.);
    
    for(int64_t iel = ElemNodeVec.size(); iel < ElemNodeVec.size() + SegmentVec.size(); iel++){
        
        //tira-se "-1" dos IDs pois o MatLab começa a numeração dos nós em "1" e não em "0"
        boundary[0] = SegmentVec[iel-ElemNodeVec.size()].fId[0] - 1;
        boundary[1] = SegmentVec[iel-ElemNodeVec.size()].fId[1] - 1;
        
        //definindo reftyoe = 0, ou seja, será utilizado padrão uniforme por default (mais rápido)
        const int reftype = 0;
        gmesh->CreateGeoElement(EOned, boundary, matBoundary, id, reftype);//cria elemento unidimensional
        
        gmesh->ElementVec()[id]->SetId(iel);
        
    }
    
    gmesh->BuildConnectivity();
    
	return gmesh;
}

///##############################################################
bool SidesToRefine(TPZGeoEl *gel, TPZVec<int> &sidestorefine)
{
    
#ifdef PZDEBUG
    if(!gel){
        DebugStop();
    }
#endif
    
    bool thereIsAnyNeighbourRefined = false;
    
    int ncorners = gel->NCornerNodes();
    int nsides = gel->NSides();
    
    sidestorefine.Resize(nsides,0);
        
    for(int s = ncorners; s < nsides; s++)
    {
        TPZGeoElSide gelside(gel, s);
        TPZGeoElSide neighside = gelside.Neighbour();
        if(!neighside.Exists())
        {
            break;
        }
        while(neighside != gelside)
        {
            if(neighside.Element()->HasSubElement() && neighside.Element()->NSideSubElements(neighside.Side()) > 1)
            {
                thereIsAnyNeighbourRefined = true;
                sidestorefine[s] = 1;
                
                //break;
            } 
            neighside = neighside.Neighbour();
        }
    }
    
    return thereIsAnyNeighbourRefined;
}


///##############################################################
#include "TPZGeoElement.h"
#include "TPZGeoCube.h"
#include "TPZGeoLinear.h"
#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "pzshapetriang.h"
#include "pzreftriangle.h"
#include "pzgeotriangle.h"
#include "pzshapeprism.h"
#include "pzrefprism.h"
#include "pzgeoprism.h"
#include "pzshapetetra.h"
#include "pzreftetrahedra.h"
#include "pzgeotetrahedra.h"
#include "pzshapepiram.h"
#include "pzrefpyram.h"
#include "pzgeopyramid.h"
#include "pzgeopoint.h"
#include "pzrefpoint.h"
#include "pzshapepoint.h"


TPZGeoEl * ChangeToGeoElRefPattern(TPZGeoMesh *Mesh, int64_t ElemIndex)
{
    TPZGeoEl * OldElem = Mesh->ElementVec()[ElemIndex];
    
    TPZCompEl * OldCompEl = OldElem->Reference();
    
    /////////////////////////
#ifdef verifyNeighbourhood
    std::ofstream before("before.txt");
    for(int s = 0; s < OldElem->NSides(); s++)
    {
        TPZGeoElSide oldSide(OldElem,s);
        TPZGeoElSide neighSide(oldSide.Neighbour());
        while(oldSide != neighSide)
        {
            before << s << "\t" << neighSide.Element()->Id() << "\t" << neighSide.Side() << "\n";
            neighSide = neighSide.Neighbour();
        }
    }
    before.close();
    TPZGeoEl * oldFather = OldElem->Father();
    int oldMePosition = -1;
    if(oldFather)
    {
        for(int s = 0; s < oldFather->NSubElements(); s++)
        {
            if(oldFather->SubElement(s) == OldElem)
            {
                oldMePosition = s;
                break;
            }
        }
    }
#endif
    /////////////////////////
    
#ifdef PZDEBUG
    if(!OldElem)
    {
        std::cout << "Null geoel on " << __PRETTY_FUNCTION__ << std::endl;
        DebugStop();
    }
#endif
    
    TPZGeoEl * father = OldElem->Father();
    
    int64_t midN;
    int nsides = OldElem->NSides();
    
    //backingup oldElem neighbourhood
    TPZVec<REAL> Coord(3);
    TPZVec< std::vector<TPZGeoElSide> > neighbourhood(nsides);
    TPZVec<int64_t> NodesSequence(0);
    for(int s = 0; s < nsides; s++)
    {
        neighbourhood[s].resize(0);
        TPZGeoElSide mySide(OldElem,s);
        TPZGeoElSide neighS = mySide.Neighbour();
        if(mySide.Dimension() == 0)
        {
            int64_t oldSz = NodesSequence.NElements();
            NodesSequence.resize(oldSz+1);
            NodesSequence[oldSz] = OldElem->NodeIndex(s);
        }
        if(TPZChangeEl::CreateMiddleNodeAtEdge(Mesh, ElemIndex, s, midN))
        {
            int64_t oldSz = NodesSequence.NElements();
            NodesSequence.resize(oldSz+1);
            NodesSequence[oldSz] = midN;
        }
        while(mySide != neighS)
        {
            neighbourhood[s].push_back(neighS);
            neighS = neighS.Neighbour();
        }
    }
    
    MElementType elType = OldElem->Type();
    int64_t oldId = OldElem->Id();
    int64_t oldMatId = OldElem->MaterialId();
    
    TPZGeoEl * NewElem = NULL;
    
    /** Deleting OldElem */
    Mesh->DeleteElement(OldElem);
    
    switch(elType) /** Inserting New Element in Mesh */
    {
        case(EOned) :
        {
            NewElem = new TPZGeoElRefPattern< TPZGeoPoint >(oldId,NodesSequence,oldMatId,*Mesh);
            break;
        }
        case(ETriangle) :
        {
            NewElem = new TPZGeoElRefPattern< TPZGeoTriangle >(oldId,NodesSequence,oldMatId,*Mesh);
            break;
        }
        case(EQuadrilateral) :
        {
            NewElem = new TPZGeoElRefPattern< TPZGeoQuad >(oldId,NodesSequence,oldMatId,*Mesh);
            break;
        }
        case(ETetraedro) :
        {
            NewElem = new TPZGeoElRefPattern< TPZGeoTetrahedra >(oldId,NodesSequence,oldMatId,*Mesh);
            break;
        }
        case(EPiramide) :
        {
            NewElem = new TPZGeoElRefPattern< TPZGeoPyramid >(oldId,NodesSequence,oldMatId,*Mesh);
            break;
        }
        case(EPrisma) :
        {
            NewElem = new TPZGeoElRefPattern< TPZGeoPrism >(oldId,NodesSequence,oldMatId,*Mesh);
            break;
        }
        case(ECube) :
        {
            NewElem = new TPZGeoElRefPattern< TPZGeoCube >(oldId,NodesSequence,oldMatId,*Mesh);
            break;
        }
        default :
        {
            DebugStop();
            break;
        }
    }
    
    if(father)
    {
        NewElem->SetFather(father);
    }
    
    // melhor utilizar neigh.SetConnectivity...
    for(int s = 0; s < nsides; s++)
    {
        TPZGeoEl * tempEl = NewElem;
        TPZGeoElSide tempSide(NewElem,s);
        int byside = s;
        for(uint64_t n = 0; n < neighbourhood[s].size(); n++)
        {
            TPZGeoElSide neighS = neighbourhood[s][n];
            tempEl->SetNeighbour(byside, neighS);
            tempEl = neighS.Element();
            byside = neighS.Side();
        }
        tempEl->SetNeighbour(byside, tempSide);
    }
    
    if(NewElem->HasSubElement())
    {
        DebugStop();
        //Mudar subelementos para TPZGeoElMapped
    }
    
    NewElem->SetReference(OldCompEl);
    
    /////////////////////////
#ifdef verifyNeighbourhood
    std::ofstream after("after.txt");
    for(int s = 0; s < NewElem->NSides(); s++)
    {
        TPZGeoElSide newSide(NewElem,s);
        TPZGeoElSide neighSide(newSide.Neighbour());
        while(newSide != neighSide)
        {
            after << s << "\t" << neighSide.Element()->Id() << "\t" << neighSide.Side() << "\n";
            neighSide = neighSide.Neighbour();
        }
    }
    after.close();
    TPZGeoEl * newFather = NewElem->Father();
    int newMePosition = -1;
    if(newFather)
    {
        for(int s = 0; s < newFather->NSubElements(); s++)
        {
            if(newFather->SubElement(s) == NewElem)
            {
                newMePosition = s;
                break;
            }
        }
    }
    if(oldFather != newFather || oldMePosition != newMePosition)
    {
        DebugStop();
    }
#endif
    /////////////////////////
    
    return NewElem;
}

///##############################################################

void RefineClosureElements(TPZCompMesh *cmesh, std::map<int64_t,TPZVec<int> > geoelindex){
    
    const int Interpolate = 1; //a solução será interpolada para os filhos
    
    std::map<int64_t,TPZVec<int> >::iterator it;

    TPZGeoMesh *gmesh = cmesh->Reference();
    
    for(it = geoelindex.begin(); it != geoelindex.end(); it++){
        
        int64_t elid = it->first;
        
        TPZVec<int> sides = it->second;
        
        TPZGeoEl * NewGeoEl = ChangeToGeoElRefPattern(gmesh, elid);

        TPZAutoPointer<TPZRefPattern> refp = TPZRefPatternTools::PerfectMatchRefPattern(NewGeoEl, sides);
        
        if(refp)
        {
            NewGeoEl->SetRefPattern(refp);
            TPZVec<int64_t> SonsIndex;
            int64_t CompElIndex = NewGeoEl->Reference()->Index();
            cmesh->Divide(CompElIndex, SonsIndex, Interpolate);
            /**for(int64_t j = 0; j < SonsIndex.NElements();j++){
             ClosureElementIndex.push_back(SonsIndex[j]);
             }*/
        }
        
        
    }
    
}

///##############################################################

// este método altera o elemento para TPZRefPattern

void RefineMesh2(TPZCompMesh *cmesh, std::vector<TPZVec<REAL> > &GLvec, const int &MaxLevel){
    
    DebugStop();
    
    const int Interpolate = 1; //a solução será interpolada para os filhos
    std::vector<int> sides(3);
    sides[0] = 3;
    sides[1] = 4;
    sides[2] = 5;
    
    int level = 0;
    
    ///distance around grounding line
    const REAL Region = MAX_GL_DISTANCE;
    
    const bool IsUniform = IS_UNIFORM;
    
    while(level < MaxLevel){
    
        level++;
        
        cmesh->Reference()->ResetReference();
        cmesh->LoadReferences();
    
        ///distance around grounding line
        const REAL MaxDistance = Region / std::exp(level-1);//(level*level); //itapopo teste
    
        ///refinando os elementos triangulares necessários
        const int64_t nelem = cmesh->NElements();

        for(int64_t i = 0; i < nelem; i++){
        
            TPZCompEl *compel = cmesh->Element(i);
            if(!compel) continue;
        
            TPZGeoEl *geoel = compel->Reference();
            if(!geoel) DebugStop();
        
            if(IsUniform){
                
                TPZVec<int64_t> SonsIndex;
                cmesh->Divide(i, SonsIndex, Interpolate);
            
            } else {
            
                const int side2D = 6;
                TPZVec<REAL> qsi(2,0.);
                TPZVec<REAL> centerPoint(3,0.);
                geoel->CenterPoint(side2D, qsi);
                geoel->X(qsi, centerPoint);
        
                REAL distance = MaxDistance;
        
                for (int64_t j = 0; j < GLvec.size(); j++) {
            
                    REAL value = ( GLvec[j][0] - centerPoint[0] ) * ( GLvec[j][0] - centerPoint[0] ); // (x2-x1)^2
                    value += ( GLvec[j][1] - centerPoint[1] ) * ( GLvec[j][1] - centerPoint[1] );// (y2-y1)^2
                    value = std::sqrt(value); ///Radius
            
                    if(value < distance){
                        distance = value; //finding the min distance to the Grounding line
                    } ///if
            
                } ///for j / GLvec.size()

                if(distance < MaxDistance){
                    TPZVec<int64_t> SonsIndex;
                    cmesh->Divide(i, SonsIndex, Interpolate);
                } else {
                    continue;
                }
                
            }// if / IsUniform
            
            //refinando os elementos unidimensionais vizinhos ao elemento refinado
            // esses elementos não tem malha computacional
            for(int j = 0; j < sides.size(); j++ ){
            
                TPZGeoElSide Neighbour = geoel->Neighbour(sides[j]);
            
                if( Neighbour.Element()->MaterialId() == 2 && !Neighbour.Element()->HasSubElement() ){
                    ///original
                    TPZVec<TPZGeoEl *> pv2;
                    Neighbour.Element()->Divide(pv2);
                    ///original
               
                }///if
            }/// for j
        
        }///for i
    
        ///refinando os elementos triangulares para tirar os hanging nodes.
        /// se for uniforme, não será preciso refinar
        if(!IsUniform){
            
            TPZGeoMesh *gmesh = cmesh->Reference();
            int64_t NElem = gmesh->NElements();
            std::map<int64_t, TPZVec<int> > geoelindex;
            
            for(int64_t i = 0; i < NElem; i++){
        
                TPZGeoEl * geoel = gmesh->Element(i);
        
                if(!geoel->HasSubElement()){ //não pode ter sido refinado antes
                    
                    TPZVec<int> sides;
                    if( SidesToRefine(geoel, sides) ){
                        
                        int64_t elindex = i;
                        geoelindex.insert( std::pair<int64_t,TPZVec<int> >(elindex, sides) );
                    
                    }/// if
                }/// if
            }/// for
            
            /// refine the elements
            RefineClosureElements(cmesh, geoelindex);

        }/// IsUniform
    
        cmesh->AdjustBoundaryElements();
        cmesh->ExpandSolution();
        cmesh->CleanUpUnconnectedNodes();
    
        cmesh->Reference()->BuildConnectivity();
        
    } ///while
}

///##############################################################
void RefineUniformGeoMesh(TPZGeoMesh *gmesh, std::vector<TPZVec<REAL> > &GLvec, const int &MaxLevel){

    std::vector<int> sides(3);
    sides[0] = 3;
    sides[1] = 4;
    sides[2] = 5;
    
    int level = 0;
    
    while(level < MaxLevel){
        
        level++;
        
        ///refinando os elementos triangulares necessários
        const int64_t nelem = gmesh->NElements();
        
        for(int64_t i = 0; i < nelem; i++){
            
            TPZGeoEl *geoel = gmesh->Element(i);
            
#ifdef PZDEBUG
            if(!geoel) DebugStop();
#endif
            
            if(geoel->MaterialId() != 1) continue;
            
            TPZVec<TPZGeoEl *> Sons;
            geoel->Divide(Sons);
            
            //refinando os elementos unidimensionais vizinhos ao elemento refinado
            // esses elementos não tem malha computacional
            for(int j = 0; j < sides.size(); j++ ){
                
                TPZGeoElSide Neighbour = geoel->Neighbour(sides[j]);
                
                if( Neighbour.Element()->MaterialId() == 2 && !Neighbour.Element()->HasSubElement() ){
                    
                    TPZVec<TPZGeoEl *> pv2;
                    Neighbour.Element()->Divide(pv2);
                    
                    
                }///if
            }/// for j
            
        }///for i
        
        gmesh->BuildConnectivity();
        
    } ///while
    
    
}

///##############################################################

void RefineGeoMesh(TPZGeoMesh *gmesh, std::vector<TPZVec<REAL> > &GLvec, const int &MaxLevel){
    
    const bool IsUniform = IS_UNIFORM;

    if(IsUniform){
        
        RefineUniformGeoMesh(gmesh, GLvec, MaxLevel);
        return;
        
    }
    
    std::vector<int> sides(3);
    sides[0] = 3;
    sides[1] = 4;
    sides[2] = 5;
    
    int level = 0;
    
    std::vector< std::vector<TPZGeoEl *> >GeoElToRefine;
    
    GeoElToRefine.resize(MaxLevel);
    
    ///distance around grounding line
    const REAL Region = MAX_GL_DISTANCE;
    
    while(level < MaxLevel){
        
        level++;
        
        ///distance around grounding line
        const REAL MaxDistance = Region / std::exp(level-1);//(level*level); //itapopo teste
        
        ///refinando os elementos triangulares necessários
        const int64_t nelem = (level == 1) ? gmesh->NElements() : GeoElToRefine[level-2].size();
        
        #ifdef PZDEBUG
        //if(nelem == 0) DebugStop();
        #endif
        
        for(int64_t i = 0; i < nelem; i++){
            
            TPZGeoEl *geoel = (level == 1) ? gmesh->Element(i) : GeoElToRefine[level-2][i];
            
            #ifdef PZDEBUG
            if(!geoel) DebugStop();
            if( geoel->HasSubElement() && geoel->MaterialId() == 1) DebugStop();
            #endif
            
            if(geoel->MaterialId() != 1) continue;
            
                
            const int side2D = 6;
            TPZVec<REAL> qsi(2,0.);
            TPZVec<REAL> centerPoint(3,0.);
            geoel->CenterPoint(side2D, qsi);
            geoel->X(qsi, centerPoint);
                
            REAL distance = MaxDistance;
                
            for (int64_t j = 0; j < GLvec.size(); j++) {
                    
                REAL value = ( GLvec[j][0] - centerPoint[0] ) * ( GLvec[j][0] - centerPoint[0] ); // (x2-x1)^2
                value += ( GLvec[j][1] - centerPoint[1] ) * ( GLvec[j][1] - centerPoint[1] );// (y2-y1)^2
                value = std::sqrt(value); ///Radius
                    
                if(value < distance){
                    distance = value; //finding the min distance to the Grounding line
                } ///if
                    
            } ///for j / GLvec.size()
                
            if(distance < MaxDistance){
                    
                TPZVec<TPZGeoEl *> Sons;
                geoel->Divide(Sons);
                for(int j = 0; j < Sons.NElements(); j++) GeoElToRefine[level-1].push_back(Sons[j]);
                
            } else {
                continue;
            }
            
            //refinando os elementos unidimensionais vizinhos ao elemento refinado
            // esses elementos não tem malha computacional
            for(int j = 0; j < sides.size(); j++ ){
                
                TPZGeoElSide Neighbour = geoel->Neighbour(sides[j]);
                
                if( Neighbour.Element()->MaterialId() == 2 && !Neighbour.Element()->HasSubElement() ){
                    
                    TPZVec<TPZGeoEl *> pv2;
                    Neighbour.Element()->Divide(pv2);
                    
                    
                }///if
            }/// for j
            
        }///for i
        
        ///refinando os elementos triangulares para tirar os hanging nodes.
        /// se for uniforme, não será preciso refinar
            
        int64_t NElem = gmesh->NElements();
            
        for(int64_t i = 0; i < NElem; i++){
                
            TPZGeoEl * geoel = gmesh->Element(i);
                
            if(!geoel->HasSubElement()){ //não pode ter sido refinado antes
                    
                TPZAutoPointer<TPZRefPattern> refp = TPZRefPatternTools::PerfectMatchRefPattern(geoel);
                if(refp)
                {
                    geoel->SetRefPattern(refp);
                    TPZVec<TPZGeoEl *> Sons;
                    geoel->Divide(Sons);
                        
                }
                    
            }
        }
        
        //cmesh->AdjustBoundaryElements();
        //cmesh->ExpandSolution();
        //cmesh->CleanUpUnconnectedNodes();
        
        gmesh->BuildConnectivity();
        
    } ///while
}

///##############################################################

void RefineMesh(TPZCompMesh *cmesh, std::vector<TPZVec<REAL> > &GLvec, const int &MaxLevel, REAL &alpha){
    
    const int Interpolate = 0; //a solução NÃO será interpolada para os filhos. Cmesh não tem solution ainda
    std::vector<int> sides(3);
    sides[0] = 3;
    sides[1] = 4;
    sides[2] = 5;
    
    int level = 0;
    
    ///distance around grounding line
    const REAL Region = MAX_GL_DISTANCE*alpha;
    
    const bool IsUniform = IS_UNIFORM;
    
    std::cout << " ----- Starting RefineMesh ----- " << std::endl;
    
    while(level < MaxLevel){
        
        level++;
        std::cout << "      Level = " << level << " " << std::endl;
        
        
        cmesh->Reference()->ResetReference();
        cmesh->LoadReferences();

        ///distance around grounding line
        const REAL MaxDistance = Region / std::exp(EXP_APLHA*(level-1));
        
        ///refinando os elementos triangulares necessários
        const int64_t nelem = cmesh->NElements();
        
        for(int64_t i = 0; i < nelem; i++){
            
            TPZCompEl *compel = cmesh->Element(i);
            if(!compel) continue;
            
            TPZGeoEl *geoel = compel->Reference();
            if(!geoel) DebugStop();
            
            if(IsUniform){
                
                TPZVec<int64_t> SonsIndex;
                cmesh->Divide(i, SonsIndex, Interpolate);
                
            } else {
                
                const int side2D = 6;
                TPZVec<REAL> qsi(2,0.);
                TPZVec<REAL> centerPoint(3,0.);
                geoel->CenterPoint(side2D, qsi);
                geoel->X(qsi, centerPoint);
                
                REAL distance = MaxDistance;
                
                for (int64_t j = 0; j < GLvec.size(); j++) {
                    
                    REAL value = ( GLvec[j][0] - centerPoint[0] ) * ( GLvec[j][0] - centerPoint[0] ); // (x2-x1)^2
                    value += ( GLvec[j][1] - centerPoint[1] ) * ( GLvec[j][1] - centerPoint[1] );// (y2-y1)^2
                    value = std::sqrt(value); ///Radius
                    
                    if(value < distance){
                        distance = value; //finding the min distance to the Grounding line
                    } ///if
                    
                } ///for j / GLvec.size()
                
                if(distance < MaxDistance){
                    TPZVec<int64_t> SonsIndex;
                    cmesh->Divide(i, SonsIndex, Interpolate);
                } else {
                    continue;
                }
                
            }// if / IsUniform
            
            //refinando os elementos unidimensionais vizinhos ao elemento refinado
            // esses elementos não tem malha computacional
            for(int j = 0; j < sides.size(); j++ ){
                
                TPZGeoElSide Neighbour = geoel->Neighbour(sides[j]);
                
                if( Neighbour.Element()->MaterialId() == 2 && !Neighbour.Element()->HasSubElement() ){
                    ///original
                    TPZVec<TPZGeoEl *> pv2;
                    Neighbour.Element()->Divide(pv2);
                    ///original
                    
                }///if
            }/// for j
            
        }///for i
        
        ///refinando os elementos triangulares para tirar os hanging nodes.
        /// se for uniforme, não será preciso refinar
        if(!IsUniform){
            
            TPZGeoMesh *gmesh = cmesh->Reference();
            int64_t NElem = gmesh->NElements();
            
            for(int64_t i = 0; i < NElem; i++){
                
                TPZGeoEl * geoel = gmesh->Element(i);
                
                if(!geoel->HasSubElement()){ //não pode ter sido refinado antes
                    
                    TPZAutoPointer<TPZRefPattern> refp = TPZRefPatternTools::PerfectMatchRefPattern(geoel);
                    if(refp)
                    {
                        geoel->SetRefPattern(refp);
                        TPZVec<int64_t> SonsIndex;
                        int64_t CompElIndex = geoel->Reference()->Index();
                        cmesh->Divide(CompElIndex, SonsIndex, Interpolate);
                        /**for(int64_t j = 0; j < SonsIndex.NElements();j++){
                         ClosureElementIndex.push_back(SonsIndex[j]);
                         }*/
                    }
                    
                }
            }
        }/// IsUniform
        
        cmesh->AdjustBoundaryElements();
        cmesh->ExpandSolution();
        cmesh->CleanUpUnconnectedNodes();
        
        cmesh->Reference()->BuildConnectivity();
        
    } ///while
}


///##############################################################
void RefineMesh(TPZCompMesh *cmesh,
                /*std::vector<int64_t>*//*std::set<int64_t>*/std::map<int64_t, int>  &ElementIndex,///element index in the geometric mesh and level of refinement
                std::vector<int64_t> &ClosureElementIndex){
    
    DebugStop();
    
    ClosureElementIndex.clear();
    
    const int Interpolate = 1; //a solução será interpolada para os filhos
    std::vector<int> sides(3);
    sides[0] = 3;
    sides[1] = 4;
    sides[2] = 5;
    
    cmesh->Reference()->ResetReference();
    cmesh->LoadReferences();
    
    
    ///refinando os elementos triangulares necessários
    //std::set<int64_t>::iterator it;
    std::map<int64_t,int>::iterator it;
    //for(int64_t i = 0; i < ElementIndex.size(); i++){
    for(it=ElementIndex.begin(); it != ElementIndex.end(); it++){
        
        const int64_t geoElIndex = it->first;//*it;//ElementIndex[i];
        const int maxLevel = it->second;
        
        TPZGeoEl *GeoEl = cmesh->Reference()->Element(geoElIndex);
        if(GeoEl->MaterialId() != 1) DebugStop();
        
        
        std::vector<int64_t> index;
        int level = 0;
        
        if(GeoEl->Reference()){
            ///original
            int64_t compElIndex = GeoEl->Reference()->Index();
            //cmesh->Divide(compElIndex, SonsIndex, Interpolate);
            ///original
            
            index.push_back(compElIndex);
            
            while(level < maxLevel){
                
                std::vector<int64_t> SonsIndex;
                SonsIndex.clear();
                
                for(int i = 0; i < index.size(); i++){
                    
                    TPZVec<int64_t> Sons;
                    Sons.clear();
                    
                    cmesh->Divide(index[i], Sons, Interpolate);
                    
                    for(int s=0; s < Sons.size();s++){
                        SonsIndex.push_back(Sons[s]);
                    }///for s
                    
                }///for i
                
                index.clear();
                index = SonsIndex;
                
                level++;
            
            }///while
            
        } else {
            DebugStop();
        }
        
        
        //refinando os elementos unidimensionais vizinhos ao elemento refinado
        // esses elementos não tem malha computacional
        for(int j = 0; j < sides.size(); j++ ){
            
            TPZGeoElSide Neighbour = GeoEl->Neighbour(sides[j]);
            
            if( Neighbour.Element()->MaterialId() == 2 && !Neighbour.Element()->HasSubElement() ){
                ///original
                TPZVec<TPZGeoEl *> pv2;
                //Neighbour.Element()->Divide(pv2);
                ///original
                
                std::vector<TPZGeoEl *> UniElemVec;
                UniElemVec.clear();
                
                UniElemVec.push_back(Neighbour.Element());
                
                int level = 0;
                while(level<maxLevel){
                    
                    std::vector<TPZGeoEl *> vecAux;
                    vecAux.clear();
                    
                    for(int i = 0; i < UniElemVec.size(); i++){
                        
                        UniElemVec[i]->Divide(pv2);
                        
                        for(int s=0; s < pv2.size(); s++){
                            vecAux.push_back(pv2[s]);
                        }
                    
                    }///for i
                    
                    UniElemVec.clear();
                    
                    UniElemVec = vecAux;
                    
                    level++;
                }///while
            }///if
        }/// for j
        
    }///for iterator
    
    ///refinando os elementos triangulares para tirar os hanging nodes.
    TPZGeoMesh *gmesh = cmesh->Reference();
    int64_t NElem = gmesh->NElements();
    
    for(int64_t i = 0; i < NElem; i++){
        
        TPZGeoEl * GeoEl = gmesh->Element(i);
        
        if(!GeoEl->HasSubElement()){ //não pode ter sido refinado antes
            
            TPZAutoPointer<TPZRefPattern> refp = TPZRefPatternTools::PerfectMatchRefPattern(GeoEl);
            if(refp)
            {
                GeoEl->SetRefPattern(refp);
                TPZVec<int64_t> SonsIndex;
                int64_t CompElIndex = GeoEl->Reference()->Index();
                cmesh->Divide(CompElIndex, SonsIndex, Interpolate);
                for(int64_t j = 0; j < SonsIndex.NElements();j++){
                    ClosureElementIndex.push_back(SonsIndex[j]);
                }
            }
            
        }
    }
   
    cmesh->AdjustBoundaryElements();
    cmesh->ExpandSolution();
    cmesh->CleanUpUnconnectedNodes();
   
    cmesh->Reference()->BuildConnectivity();
}


///##############################################################
void RefineMesh(TPZGeoMesh *gmesh,
                std::vector<int> &ElemVec,
                int nlevel){ //itapopo nlevel não é usado por enquanto, default = 1
    
    DebugStop();
    
    ///refinando os elementos triangulares necessários
    std::vector<int> sides(3);
    sides[0] = 3;
    sides[1] = 4;
    sides[2] = 5;
    
    for(int64_t i = 0; i < ElemVec.size(); i++){
        
        int64_t ID = ElemVec[i];
        
        TPZGeoEl *GeoEl = gmesh->FindElement(ID);
        
        if(GeoEl->MaterialId() != 1) DebugStop();
        
        TPZVec<TPZGeoEl *> sons;
        GeoEl->Divide(sons);
        
        //itapopo arrumar para permitir mais refinamentos
       /** if(nlevel == 2){
            TPZVec<TPZGeoEl *> pv;
            for(int s = 0; s < sons.size(); s++){
                sons[s]->Divide(pv);
                for(int t = 0; t < pv.size(); t++){
                    TPZVec<TPZGeoEl *> pvTeste;
                    pv[t]->Divide(pvTeste);
                    //pv[t]->Level()
                }
            }
        }*/
        
        //refinando os elementos unidimensionais vizinhos ao elemento refinado
        for(int j = 0; j < sides.size(); j++ ){
            
            TPZGeoElSide Neighbour = GeoEl->Neighbour(sides[j]);
            
            if( Neighbour.Element()->MaterialId() == 2 && !Neighbour.Element()->HasSubElement() ){
                TPZVec<TPZGeoEl *> pv2;
                Neighbour.Element()->Divide(pv2);
            }
        }
        
    }
    
    ///refinando os elementos triangulares para tirar os hanging nodes.
    
    int NElem = gmesh->NElements();
    
    for(int64_t i = 0; i < NElem; i++){
        
        TPZGeoEl * GeoEl = gmesh->Element(i);
        
        if(!GeoEl->HasSubElement()){ //não pode ter sido refinado antes
            
            TPZAutoPointer<TPZRefPattern> refp = TPZRefPatternTools::PerfectMatchRefPattern(GeoEl);
            if(refp)
            {
                GeoEl->SetRefPattern(refp);
                TPZVec<TPZGeoEl*> sons;
                GeoEl->Divide(sons);
            }
            
        }
    }
    
    for(int i = 0; i < gmesh->NNodes(); i++){
        
        gmesh->NodeVec()[i].SetNodeId(i);
        
    }
    
    gmesh->BuildConnectivity();
    
}

///##############################################################

void RefineMeshCopy(TPZGeoMesh *mesh,
                std::vector<std::pair<int,int> > &FlagVec,
                int meshID){
    DebugStop();
    
///refinando os elementos triangulares necessários (vindos do MatLab)
    std::vector<int> sides(3);
    sides[0] = 3;
    sides[1] = 4;
    sides[2] = 5;
    
    for(int i = 0; i < FlagVec.size(); i++){
        
        int64_t ID = FlagVec[i].first;
        
        TPZGeoEl *GeoEl = mesh->FindElement(ID-1);//o -1 eh devido a numeracao que vem do MatLab (comeca em 1)
        
        if(GeoEl->MaterialId() != 1) DebugStop();
        
        const int nlevel = FlagVec[i].second;
        
        TPZVec<TPZGeoEl *> sons;
        GeoEl->Divide(sons);
        
        //itapopo arrumar para permitir mais refinamentos
        if(nlevel == 2){
            TPZVec<TPZGeoEl *> pv;
            for(int s = 0; s < sons.size(); s++){
                sons[s]->Divide(pv);
                for(int t = 0; t < pv.size(); t++){
                    TPZVec<TPZGeoEl *> pvTeste;
                    pv[t]->Divide(pvTeste);
                    //pv[t]->Level()
                }
            }
        }
        
        for(int j = 0; j < sides.size(); j++ ){
            
            TPZGeoElSide Neighbour = GeoEl->Neighbour(sides[j]);
            
            if( Neighbour.Element()->MaterialId() == 2 && !Neighbour.Element()->HasSubElement() ){
                TPZVec<TPZGeoEl *> pv2;
                Neighbour.Element()->Divide(pv2);
            }
        }
        
    }
    
///refinando os elementos triangulares para tirar os hanging nodes.
    
    int NElem = mesh->NElements();
    
    for(int64_t i = 0; i < NElem; i++){
        
        TPZGeoEl * GeoEl = mesh->Element(i);
        
        if(!GeoEl->HasSubElement()){ //não pode ter sido refinado antes
        
            TPZAutoPointer<TPZRefPattern> refp = TPZRefPatternTools::PerfectMatchRefPattern(GeoEl);
            if(refp)
            {
                GeoEl->SetRefPattern(refp);
                TPZVec<TPZGeoEl*> sons;
                GeoEl->Divide(sons);
            }
        
        }
    }
    
///refinando os elementos unidimensionais de contorno
    ///ITAPOPO por enquanto todos elementos estão sendo refinados

 /*   NElem = mesh->NElements();
    
    for(int i = 0; i < NElem; i++){
        
        if( mesh->ElementVec()[i]->MaterialId() == 2 ){
            
            TPZVec<TPZGeoEl *> pv;
            
            mesh->ElementVec()[i]->Divide(pv);
        }
        
    }*/
    
    
    
    for(int i = 0; i < mesh->NNodes(); i++){
        
        mesh->NodeVec()[i].SetNodeId(i);
        
    }
    
    mesh->BuildConnectivity();
    

}

///##############################################################
void PrintNewMesh(std::ofstream &MeshFile, int64_t &nvertices, int64_t &nelements, int64_t &nsegments,int &elementswidth,double *x, double *y, double *z, int64_t ** elements, int64_t **segments){
    
    if(!MeshFile.is_open()) DebugStop();
    
    MeshFile.precision(12);
    MeshFile << std::scientific;
    
    MeshFile << nvertices << "\n";
    for(int64_t i = 0; i < nvertices; i++ ) MeshFile << x[i] << "\t" << y[i] << "\t" << "\n";
    
    MeshFile << nelements << "\n";
    for(int64_t i = 0; i < nelements; i++ ){
        for(int64_t j = 0; j < 3; j++){
            MeshFile << elements[i][j] +1 << "\t";//soma-se 1 quando escreve-se para o Matlab
        }
        MeshFile << "\n";
    }
    
    MeshFile << nsegments << "\n";
    for(int64_t i=0;i<nsegments;i++){
        for(int64_t j=0;j<3;j++){
            MeshFile << segments[i][j] +1 << "\t";//soma-se 1 quando escreve-se para o MatLab
        }
        MeshFile << "\n";
    }
    
    MeshFile.flush();
    MeshFile.close();

}

///##############################################################
void PrintNewMesh(std::string &MeshName,
                  TPZGeoMesh *mesh){
    
    std::string FullName = toStr(WORKPATH) + MeshName;
    
    std::ofstream MeshFile(FullName.c_str());
    
    if(!MeshFile.is_open()) DebugStop();
    
    MeshFile << std::setprecision(18);
    
///Printing nodes
    int64_t nnodes = mesh->NNodes();
    
    MeshFile << nnodes << "\n";
    
    for(int64_t i = 0; i < nnodes; i++ ){
        TPZVec<REAL> coords(3,0.);
        mesh->NodeVec()[i].GetCoordinates(coords);
        MeshFile << coords[0] << "\t" << coords[1] << "\t" << "\n";
        
    }
    
///Printing elements and Father Index
    int64_t nTotalElements = mesh->NElements();
    
    std::vector<TPZGeoEl *> TriaVec;
    std::vector<TPZGeoEl *> SegmentVec;
    TriaVec.clear();
    SegmentVec.clear();
    
    
    for(int64_t i = 0; i < nTotalElements; i++){
    
        if( mesh->ElementVec()[i]->HasSubElement() ) continue;
        
        if( mesh->ElementVec()[i]->MaterialId() == 1 ){
        //separando apenas os elementos triangulares que não tem filhos
            //if( mesh->ElementVec()[i]->HasSubElement() ) continue;
        
            TriaVec.push_back( mesh->ElementVec()[i]);
        }
        
        if( mesh->ElementVec()[i]->MaterialId() == 2 ){
        //separando apenas os elementos unidimensionais que não tem filhos
            //if( mesh->ElementVec()[i]->HasSubElement() ) continue;
            
            SegmentVec.push_back( mesh->ElementVec()[i]);
        }
    }
    
    MeshFile << TriaVec.size() << "\n";
    
    for(int64_t i = 0; i < TriaVec.size(); i++){
        
        //int RefID = TriaVec[i]->FatherIndex();
        //if(RefID < 0) RefID = TriaVec[i]->Index();
        
        MeshFile << TriaVec[i]->NodeIndex(0)+1 << "\t" << TriaVec[i]->NodeIndex(1)+1 << "\t" << TriaVec[i]->NodeIndex(2)+1 << /*"\t" << RefID+1 <<*/ "\n"; //o 1 eh devido a numercao do MatLab
    }
    
//Printing segments
    MeshFile << SegmentVec.size() << "\n";
    
    for(int64_t i = 0; i < SegmentVec.size(); i++){
        
        TPZGeoElSide Neighbour = SegmentVec[i]->Neighbour(2);
        int NeighbourID = -1;
        
        for(int j = 0; j < TriaVec.size(); j++){
            
            // este if deve encontrar o elemento triangulo refinado uniformemente
            if( TriaVec[j]->Index() == Neighbour.Element()->Index() ){
                NeighbourID = j;
                break;
            
            // este if deve verificar se o elemento unidimensional é vizinho de algum pai.
            // Isso ocorre quando elementos de fechamento que estão na borda são refinados de forma a fechar novamente.
            // Caso muito raro de acontecer = refinarmento 2 vezes nesses elementos
            } else if(TriaVec[j]->FatherIndex() == Neighbour.Element()->Index() ){
               
                NeighbourID = j;
                break;
                
               /* TPZGeoEl * Father = TriaVec[j]->Father();
                
                while (Father){
                    if(Father->Index() == Neighbour.Element()->Index() ){
                        NeighbourID = j;
                        break;
                    } else {
                        Father = Father->Father();
                    }
                }
                
                delete Father;*/
            
            }//else
            
        }
        
        if(NeighbourID == -1) DebugStop();
        
        MeshFile << SegmentVec[i]->NodeIndex(0)+1 << "\t" << SegmentVec[i]->NodeIndex(1)+1 << "\t" << NeighbourID+1 <<"\n"; //o 1 eh devido a numercao do MatLab
    }
    
    MeshFile.flush();
    MeshFile.close();
    
}

///##############################################################

void PrintInitialData(std::string &DataName,
                      TPZGeoMesh *gmesh){
    
    std::string FullName = toStr(WORKPATH) + DataName;
    
    std::ofstream DataFile(FullName.c_str());
    
    if(!DataFile.is_open()) DebugStop();
    
    TPZCompMesh *cmesh = gmesh->Reference();
    TPZFMatrix<STATE> sol = cmesh->Solution();
    
    std::map<int, TPZMaterial*>::iterator matit;
    int nstate = -1;
    
    for(matit = cmesh->MaterialVec().begin(); matit != cmesh->MaterialVec().end(); matit++){
        
        int matid = matit->first;
    
        if(matid == 1){
            nstate = matit->second->NStateVariables();
            break;
        }
    }

    if(nstate == -1) DebugStop();
    
    //std::scientific;
    //MeshFile.precision(18);
    
    ///Printing nnodes
    ///Importante: deve-se garantir de que o numero de nós da malha é igual ao número de connects com funções de forma associado.
    int64_t nnodes = gmesh->NNodes();
    
    DataFile << nnodes << "\n";
    
    DataFile << std::setprecision(18);
    
    ///Printing surface data
    for(int64_t i = 0; i < nnodes; i++ ){
        
        STATE value = sol(0+i*nstate,0);
        
        DataFile << value << "\n";
        
    }
    
    ///Printing base data
    for(int64_t i = 0; i < nnodes; i++ ){
        
        STATE value = sol(1+i*nstate,0);
        
        DataFile << value << "\n";
        
    }
    
    ///Printing bed data
    for(int64_t i = 0; i < nnodes; i++ ){
        
        STATE value = sol(2+i*nstate,0);
        
        DataFile << value << "\n";
        
    }
    
    ///Printing pressure data
    for(int64_t i = 0; i < nnodes; i++ ){
        
        STATE value = sol(3+i*nstate,0);
        
        DataFile << value << "\n";
        
    }
    
    ///Printing temperature data
    for(int64_t i = 0; i < nnodes; i++ ){
        
        STATE value = sol(4+i*nstate,0);
        
        DataFile << value << "\n";
        
    }
    
    ///Printing vx data
    for(int64_t i = 0; i < nnodes; i++ ){
        
        STATE value = sol(5+i*nstate,0);
        
        DataFile << value << "\n";
        
    }
    
    ///Printing vy data
    for(int64_t i = 0; i < nnodes; i++ ){
        
        STATE value = sol(6+i*nstate,0);
        
        DataFile << value << "\n";
        
    }
    
    ///Printing masklevelset data
    for(int64_t i = 0; i < nnodes; i++ ){
        
        STATE value = sol(7+i*nstate,0);
        
        DataFile << value << "\n";
        
    }
    
    DataFile.flush();
    DataFile.close();
    
}
