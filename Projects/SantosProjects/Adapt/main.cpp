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

#include "/Applications/MATLAB_R2014a.app/extern/include/engine.h"

struct ElemNode{
    
    ElemNode(){
        fNodeId.clear();
    }
    
    std::vector<long> fNodeId;
};

//----------------------------------------

struct Segment{
    
    Segment(){
        fId.clear();
    }
    
    //node i, node i+1, element j
    std::vector<long> fId;
};

//----------------------------------------

bool ReadMesh(std::ifstream &MeshFile,
              std::vector<double> &x,
              std::vector<double> &y,
              std::vector<ElemNode> &ElemNodeVec,
              std::vector<Segment> &SegmentVec,
              std::vector<long> &ElemIdVec);

//----------------------------------------
TPZGeoMesh *CreateGMesh(std::vector<double> &x,
                        std::vector<double> &y,
                        std::vector<ElemNode> &ElemNodeVec,
                        std::vector<Segment> &SegmentVec);

TPZGeoMesh *CreateGMesh(std::string &MeshName);

//----------------------------------------
TPZCompMesh *CreateCompMesh(TPZGeoMesh *gmesh);

//----------------------------------------
void Permute(TPZCompMesh * cmesh);

//----------------------------------------
void RefineMeshCopy(TPZGeoMesh *mesh,
                std::vector<std::pair<int,int> > &FlagVec,
                int meshID);

//----------------------------------------
void RefineMesh(TPZGeoMesh *gmesh,
                std::vector<int> &ElemVec,
                int nlevel);

//----------------------------------------
void RefineMesh(TPZCompMesh *cmesh,
                /*std::vector<long>*//*std::set<long>*/std::map<long, int>  &ElementIndex,
                std::vector<long> &ClosureElementIndex);

//----------------------------------------
void RefineMesh(TPZCompMesh *cmesh, std::vector<TPZVec<REAL> > &GLvec, const int &MaxLevel);


//----------------------------------------
void PrintNewMesh(std::string &MeshName,
                  TPZGeoMesh *mesh);

//----------------------------------------
void PrintInitialData(std::string &DataName,
                      TPZGeoMesh *mesh);


//----------------------------------------
void FlagElements(TPZGeoMesh *gmesh, std::vector<std::pair<int,int> > &FlagVec, int step, int dstep);

//----------------------------------------
void SetElementsToRefine(TPZCompMesh *cmesh,               // mesh with solution
                         /*std::vector<long>*//*std::set<long>*/std::map<long, int> &ElementIndex,
                         std::vector<TPZVec<REAL> > &GLvec);   // elementos of cmesh to refine

//----------------------------------------
void DeleteClosureElements(TPZCompMesh *cmesh,               // mesh with solution
                          std::vector<long> &ClosureElementIndex);   // elementos of cmesh to delete

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


//----------------------------------------
void TransferSolution(TPZCompMesh *cmesh, TPZCompMesh *cmesh2);

//----------------------------------------
void GetVariableNames(TPZStack<std::string> &scalnames, TPZStack<std::string> &vecnames);


//----------------------------------------
void SetFileNames(const int &run,
                  const int &step,
                  std::string &strSolutionFile,
                  std::string &strMeshFile,
                  std::string &strDataFile,
                  std::string &strSolutionMatlab,
                  std::string &strMeshMatlab,
                  std::string &strDataMatlab);


//----------------------------------------
void PrintMesh(const int &step,
              TPZGeoMesh *gmesh);

//----------------------------------------

/// path do projeto a ser executado
#define WORKPATH "/Users/santos/Documents/_PROJETOS/Criosfera/Mismip2D/"

#define MY_REFPATTERNDIR "/Users/santos/Documents/NeoPZ/neopz/Projects/SantosProjects/Adapt/rpt/"

//pastas a serem geradas no WORKPATH
#define MKDEBUG "mkdir /Users/santos/Documents/_PROJETOS/Criosfera/Mismip2D/debug"
#define MKSOL "mkdir /Users/santos/Documents/_PROJETOS/Criosfera/Mismip2D/sol"
#define MKDATA "mkdir /Users/santos/Documents/_PROJETOS/Criosfera/Mismip2D/data"
#define MKMESH "mkdir /Users/santos/Documents/_PROJETOS/Criosfera/Mismip2D/mesh"
#define MKVTK "mkdir /Users/santos/Documents/_PROJETOS/Criosfera/Mismip2D/vtk"


#define BUFSIZE 1024

#define CALLMATLAB

int main(int argc, char *argv[])
{
    
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
    if (!(ep = engOpen(NULL))) {
        fprintf(stderr, "\nCan't start MATLAB engine\n");
        return EXIT_FAILURE;
    }
    
for(int nrun=4; nrun < 5; nrun++){//ITAPOPO
    
    std::cout << "#####     RUNNING H = " << nrun << "     #####" << std::endl;
    
    /// generating the files in the WORKPATH folder
   /* 
    itapopo arrumar a geração dos arquivos
    char teste;
    const char nfile = nrun;
    std::strcat(&teste,MKSOL);
    std::strcat(&teste,&nfile);
    
    system(&teste);
    
    system(MKDEBUG);
    system(MKSOL);
    system(MKDATA);
    system(MKMESH);
    system(MKVTK);*/
    
    //clear the workspace
    engEvalString(ep, "clear");
    
    //set the buffer
    engOutputBuffer(ep, buffer, BUFSIZE);
    
    // Set the ISSM_DIR, system variable - ITAPOPO isso deveria ser lido pelo processo que dispara o matlab
    engEvalString(ep, "setenv('ISSM_DIR', '/Users/santos/Documents/issm/trunk');");
    
    // add the path of ISSM and MatLab codes
    engEvalString(ep, "addpath /Users/santos/Documents/issm/trunk/bin /Users/santos/Documents/issm/trunk/lib /Users/santos/Documents/NeoPZ/neopz/Projects/SantosProjects/Adapt/matlab");
    
    //change to the WORKPATH
    std::string cdWORKPATH = "cd " + toStr(WORKPATH);
    engEvalString(ep, cdWORKPATH.c_str());//"cd /Users/santos/Documents/_PROJETOS/Criosfera/Adapt2D/");

    //set the variable isRun: 0, print mesh 0; 1, run the ISSM with new mesh
    engEvalString(ep, "isRun=0;");
    
    //set the expfile
    engEvalString(ep, "expfile='Domain.exp';");
    
    //set the parfile
    engEvalString(ep, "parfile='BC.par';");
    
    //set the resolution to generate the mesh0
    engEvalString(ep, "resolution=20000;");//5000
    
#endif
    //execute the main to generate the first mesh, mesh0
    step=0;
    SetFileNames(nrun, step, strSolutionFile, strMeshFile,strDataFile,
                 strSolutionMatlab,strMeshMatlab,strDataMatlab);
#ifdef CALLMATLAB
    engEvalString(ep, strMeshMatlab.c_str());
    engEvalString(ep, strSolutionMatlab.c_str());
    engEvalString(ep, "main");
    //print the buffer
    //printf("%s", buffer);
#endif
    
    ///generating the father mesh
    TPZGeoMesh *FatherGMesh = CreateGMesh(strMeshFile);
    
    ///generating comp mesh  and loading results
    TPZGeoMesh *gmesh = new TPZGeoMesh(*FatherGMesh);
    TPZCompMesh *cmesh = CreateCompMesh(gmesh);
    LoadSolution(strSolutionFile, cmesh);
    
#ifdef DEBUG
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
    std::string plotfile = toStr(WORKPATH) + "vtk" + toStr(nrun) + "/solution.vtk";

    an.DefineGraphMesh(dim, scalnames, vecnames, plotfile);//define malha grafica
    an.PostProcess(postProcessResolution);//realiza pos processamento
   
    // set elements to refine
    /*std::vector<long>*/ /*std::set<long>*/std::map<long, int> ElementIndex;
    std::vector<long> ClosureElementIndex;//elements to close (avoid hanging nodes)
    
    // refine mesh
    /// itapopo
    TPZGeoMesh *gmeshTeste = new TPZGeoMesh(*FatherGMesh);
    TPZCompMesh *cmeshTeste = CreateCompMesh(gmeshTeste);
    
    ///itapopo
    const int MaxLevel = nrun;//1;///máximo nível de refinamento
    std::vector<TPZVec<REAL> > GLvec;
    SetElementsToRefine(cmesh, ElementIndex,GLvec);
    //RefineMesh(cmeshTeste, ElementIndex, ClosureElementIndex);
    RefineMesh(cmeshTeste, GLvec, MaxLevel);
    Permute(cmeshTeste);
    
    
    ///itapopo
    TransferSolution(cmesh, cmeshTeste);
    ///itapopo
    
    
    /// print the new mesh and the initial data
    step++;
#ifdef DEBUG
    PrintMesh(step, gmeshTeste);
#endif

    SetFileNames(nrun, step, strSolutionFile, strMeshFile,strDataFile,
                 strSolutionMatlab,strMeshMatlab,strDataMatlab);
    
    PrintNewMesh(strMeshFile, gmeshTeste);
    PrintInitialData(strDataFile, gmeshTeste);
    
    /// Run with new mesh
    
    //set the variable isRun: 0, print mesh 0; 1, run the ISSM with new mesh
    engEvalString(ep, "isRun=1;");
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
        engEvalString(ep, strDataMatlab.c_str());
        engEvalString(ep, strMeshMatlab.c_str());
        engEvalString(ep, strSolutionMatlab.c_str());
        engEvalString(ep, "main");
    
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
        SetElementsToRefine(cmesh, ElementIndex, GLvec);
        //RefineMesh(cmesh2, ElementIndex, ClosureElementIndex);//cmesh
        RefineMesh(cmesh2, GLvec, MaxLevel);
        Permute(cmesh2);//cmesh
        
        /// transfering the solution to the new mesh
        TransferSolution(cmesh, cmesh2);
        
        //print the new mesh
        step++;
#ifdef DEBUG
        PrintMesh(step, gmesh2);
#endif
        SetFileNames(nrun, step, strSolutionFile, strMeshFile,strDataFile,
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
    std::string SaveModel = "save model" + toStr(nrun) + " md;";
    engEvalString( ep, SaveModel.c_str() );//"save model md;"); /// alterar de acordo com a simulação
    
    //itapopo
    delete cmesh;
    delete gmesh;
    //itapopo
    
}//for nrun
    
    //close the workspace and close the engine
    engEvalString(ep, "close;");
    engClose(ep);
   
    //delete gmesh;
   // delete cmesh;
    //delete gmesh;
    
	std::cout << "FINISHED!" << std::endl;
	
	return 0;
}

///##############################################################
void TransferSolution(TPZCompMesh *cmesh1, TPZCompMesh *cmesh2){
    
    /// 1 solution
    const int nsol = 1;
    
    /// 8 state variables
    const int nstate = 8;
    
    /// number of nodes of the mesh2
    if(!cmesh2->Reference()) DebugStop();
    if(!cmesh1->Reference()) DebugStop();
    
    const long nnodes = cmesh2->Reference()->NNodes();
    
    /// solution of the mesh2
    TPZFMatrix<STATE> sol2(nnodes*nstate, nsol, 0);
    
    /// solution of the mesh1
    TPZFMatrix<STATE> sol1 = cmesh1->Solution();
    
    REAL Tol;
    ZeroTolerance(Tol);
    //const REAL Tol = 10E-6;
    
    for(long i = 0; i < cmesh2->NElements(); i++){
        
        TPZCompEl *compel2 = cmesh2->Element(i);
        if( !compel2 ) continue;
        
        if( !compel2->Reference() ) DebugStop();
        TPZGeoEl *geoel2 = compel2->Reference();
        if( geoel2->HasSubElement() ) DebugStop();
        
        /// this ElIndex is the same in the gmesh1
        long ElIndex;
        if( geoel2->LowestFather() ){
            ElIndex = geoel2->LowestFather()->Index();
        } else {
            ElIndex = geoel2->Index();//element i is not refined
        }
        
        TPZGeoEl *geoel1 = cmesh1->Reference()->Element(ElIndex);
        
        if( geoel1->HasSubElement() ){
            
            TPZVec<TPZGeoEl *> unrefinedSons;
            geoel1->GetHigherSubElements(unrefinedSons);
            
            for(int n = 0; n < geoel2->NNodes();n++){
                
                TPZVec<REAL> X(3,0);
                geoel2->NodePtr(n)->GetCoordinates(X);
                const long nodeIndex = geoel2->NodeIndex(n);
                TPZVec<REAL> qsi(2,0);
                TPZVec<STATE> sol(1,0);
                
                for(int s = 0; s < unrefinedSons.size(); s++){
                
                    if(!unrefinedSons[s]->ComputeXInverse(X, qsi, Tol) ) continue;
                
                    if(!unrefinedSons[s]->Reference() ) DebugStop();
                    
                    for(int var = 1; var <= nstate; var++ ){
                    
                        unrefinedSons[s]->Reference()->Solution(qsi, var, sol);
                        sol2(var-1 + nodeIndex*nstate, 0) = sol[0];
                        
                    }///for var
                
                }/// for s
           
            }/// for n
            
        } else { ///geoel1 is not refined (level = 0)
            
            for(int n = 0; n < geoel2->NNodes();n++){
                
                TPZVec<REAL> X(3,0);
                geoel2->NodePtr(n)->GetCoordinates(X);
                const long nodeIndex = geoel2->NodeIndex(n);
                
                TPZVec<REAL> qsi(2,0);
                TPZVec<STATE> sol(1,0);
                
                if(geoel1->ComputeXInverse(X, qsi, Tol)){
                
                    for(int var = 1; var <= nstate; var++ ){
                        
                        geoel1->Reference()->Solution(qsi, var, sol);
                        sol2(var-1 + nodeIndex*nstate, 0) = sol[0];
                
                    }///for var
                
                } else {/// if
                 
                    DebugStop();
                }///else
        
            }/// for n
            
            
        }/// if else
        
    }/// for i
    
    cmesh2->LoadSolution(sol2);
    
   /*
    
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
    
    */
    
    
    
    
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

void SetFileNames(const int &run,
                  const int &step,
                  std::string &strSolutionFile,
                  std::string &strMeshFile,
                  std::string &strDataFile,
                  std::string &strSolutionMatlab,
                  std::string &strMeshMatlab,
                  std::string &strDataMatlab){
    
    
    strSolutionFile = "sol" + toStr(run) + "/solution" + toStr(step) + ".txt";
    strSolutionMatlab = "solutionfile = '" + strSolutionFile + "';";
    
    strMeshFile = "mesh" + toStr(run) + "/mesh" + toStr(step) + ".txt";
    strMeshMatlab = "meshfile = '" + strMeshFile + "';";
    
    strDataFile = "data" + toStr(run) + "/data" + toStr(step) + ".txt";
    strDataMatlab = "datafile = '" + strDataFile + "';";
    
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
void DeleteClosureElements(TPZCompMesh *cmesh,
                           std::vector<long> &ClosureElementIndex){
    
    
    // Apagando elemento de contorno que impoe pressao
  /**  TPZGeoEl *gel = this->FindPressureBCElement();
    fgmesh->ResetReference();
    fmeshvec[0]->LoadReferences();
    delete gel->Reference();
    fgmesh->ResetReference();
    fmeshvec[1]->LoadReferences();
    delete gel->Reference();
    fgmesh->ResetReference();
    fcmeshMixed->LoadReferences();
    gel->Reference()->SetFreeIntPtIndices();
    delete gel->Reference();
    gel->RemoveConnectivities();
    delete gel;
    
    fmeshvec[0]->CleanUpUnconnectedNodes();
    fmeshvec[1]->CleanUpUnconnectedNodes();
    fcmeshMixed->CleanUpUnconnectedNodes();
    
   */
    for(long i = 0; i , ClosureElementIndex.size(); i++){
        
        TPZGeoEl *geoel = cmesh->Element(i)->Reference();
        //cmesh->Element(i)->
        //geoel->
        
        
    }

}
///##############################################################
void SetElementsToRefine(TPZCompMesh *cmesh,
                         /*std::vector<long>*//*std::set<long>*/std::map<long, int> &ElementIndex,
                         std::vector<TPZVec<REAL> > &GLvec){
    
    ElementIndex.clear();///index of the element in the geometric mesh index
    
    const long nelem = cmesh->NElements();
    const int var = 8;//masklevelset
    
    const REAL MaxLevelSet = 200.;
    
    const bool IsUniform = true;//itapopo colocar false
    
    const bool IsUsingLevelSetValue = false;
    
    const bool IsRadiusValue = !IsUniform;
    
    const REAL MaxDistance = 150000.; ///150 km
    
    //std::vector<TPZVec<REAL> > GLvec;
    GLvec.clear();
    
    for(long i = 0; i < nelem; i++){
        
        TPZCompEl *compEl = cmesh->Element(i);
        
        if(!compEl) continue;
        
        TPZGeoEl *geoel = compEl->Reference();
        
        if( geoel->HasSubElement() ) continue;
        
        long FatherIndex;
        
        ///uniform refinement
        if(IsUniform){
            
            if( geoel->LowestFather() ) {
                FatherIndex = geoel->LowestFather()->Index();
            } else {
                FatherIndex = geoel->Index();
            }
            ElementIndex.insert(std::make_pair<long,int>(FatherIndex,1));//insert(FatherIndex);//push_back(i);
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
        
            ElementIndex.insert(std::make_pair<long,int>(FatherIndex,1));//insert(FatherIndex);//push_back(i);
            
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
        
        for(long i = 0; i < nelem; i++){
            
            TPZCompEl *compEl = cmesh->Element(i);
            
            if(!compEl) continue;
            
            TPZGeoEl *geoel = compEl->Reference();
            
            if( geoel->HasSubElement() ) continue;
            
            long FatherIndex;
            
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
            
            for (long j = 0; j < GLvec.size(); j++) {
            
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
                ElementIndex.insert(std::make_pair<long,int>(FatherIndex,2));//insert(FatherIndex);//push_back(i);
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
                ElementIndex.insert(std::make_pair<long,int>(FatherIndex,1));//insert(FatherIndex);//push_back(i);
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
    
    Permute(cmesh);
    
    return cmesh;
    
    
}

///##############################################################

void Permute(TPZCompMesh * cmesh){
    

    TPZVec<long> permute;
    long numinternalconnects = cmesh->NIndependentConnects();
    permute.Resize(numinternalconnects,0);
    
    long nel=cmesh->ElementVec().NElements();
    for (long jel=0; jel<nel; jel++) {
        
        /**for (long ip=0; ip<permute.NElements(); ip++) {
            permute[ip]=ip;
        }*/
        
        TPZCompEl *compEl = cmesh->ElementVec()[jel];
        if(!compEl) continue;
        
        TPZGeoEl *geoEl = compEl->Reference();
        if(!geoEl) continue;
        
        const long nnodes = geoEl->NNodes();
        
        for(long node = 0; node < nnodes; node++){
            
            for (long ip=0; ip<permute.NElements(); ip++) {
                permute[ip]=ip;
            }
            
            long nodeindex = geoEl->NodeIndex(node);
            long seqnum = compEl->Connect(node).SequenceNumber();
            
            long v1 = permute[seqnum];
            permute[nodeindex] = v1;
            permute[seqnum] = nodeindex;
            
            /**for(int p=0; p < permute.NElements();p++){
                std::cout << permute[p] << std::endl;
            }*/
            
            cmesh->Permute(permute);
            
        }

        //cmesh->Permute(permute);
        
    }		
    
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
    long nnodes;
    
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
    for(long i = 0; i < nnodes; i++){
        
        double value;
        
        SolutionFile >> value;
   
        surface.push_back(value);
        
    }
    
    /// reading base
    for(long i = 0; i < nnodes; i++){
        
        double value;
        
        SolutionFile >> value;
        
        base.push_back(value);
        
    }
    
    /// reading bed
    for(long i = 0; i < nnodes; i++){
        
        double value;
        
        SolutionFile >> value;
        
        bed.push_back(value);
        
    }
    
    /// reading pressure
    for(long i = 0; i < nnodes; i++){
        
        double value;
        
        SolutionFile >> value;
        
        pressure.push_back(value);
        
    }
    
    /// reading temperature
    for(long i = 0; i < nnodes; i++){
        
        double value;
        
        SolutionFile >> value;
        
        temperature.push_back(value);
        
    }
    
    /// reading vx
    for(long i = 0; i < nnodes; i++){
        
        double value;
        
        SolutionFile >> value;
        
        vx.push_back(value);
        
    }
    
    /// reading vy
    for(long i = 0; i < nnodes; i++){
        
        double value;
        
        SolutionFile >> value;
        
        vy.push_back(value);
        
    }
    
    /// reading masklevelset
    for(long i = 0; i < nnodes; i++){
        
        double value;
        
        SolutionFile >> value;
        
        masklevelset.push_back(value);
        
    }
    
    return true;
    
}

///##############################################################

void SetRefPatterns(){
    
    //gRefDBase.InitializeUniformRefPattern(ETriangle);
    
    //gRefDBase.InitializeRefPatterns();
    std::string filepath = REFPATTERNDIR;
    std::string filename1 = filepath + "/2D_Triang_Rib_3.rpt";
    std::string filename2 = filepath + "/2D_Triang_Rib_4.rpt";
    std::string filename3 = toStr(MY_REFPATTERNDIR) + "2D_Triang_Rib2_Side_3_4.rpt";
    std::string filename4 = toStr(MY_REFPATTERNDIR)  + "2D_Triang_Rib2_Side_3_4permuted.rpt";
    std::string filename5 = toStr(MY_REFPATTERNDIR)  + "2D_Triang_Rib2_Side_3_5.rpt";
    std::string filename6 = toStr(MY_REFPATTERNDIR)  + "2D_Triang_Rib2_Side_3_5permuted.rpt";
    std::string filename7 = toStr(MY_REFPATTERNDIR)  + "2D_Triang_Rib_5.rpt";
    
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
    std::vector<long> ElemIdVec;
    
    /// reading mesh
    IsOk = ReadMesh(MeshFile, x, y, ElemNodeVec, SegmentVec, ElemIdVec);
    if(!IsOk) DebugStop();
    
    /// creating gmesh
    TPZGeoMesh *gmesh = CreateGMesh(x, y, ElemNodeVec, SegmentVec);
    
    return gmesh;
    
}


///##############################################################

bool ReadMesh(std::ifstream &MeshFile,
              std::vector<double> &x,
              std::vector<double> &y,
              std::vector<ElemNode> &ElemNodeVec,
              std::vector<Segment> &SegmentVec,
              std::vector<long> &ElemIdVec){
    
    if( !MeshFile.is_open() ) return false;
    
/// reading nodes
    long nnodes;
    
    MeshFile >> nnodes;
    
    x.clear();
    y.clear();
    
    for(long i = 0; i < nnodes; i++){
        
        double xvalue, yvalue;
        
        MeshFile >> xvalue;
        MeshFile >> yvalue;
        
        x.push_back(xvalue);
        y.push_back(yvalue);
        
    }

///reading elements
    long nelements;
    
    MeshFile >> nelements;
    
    ElemNodeVec.clear();
    
    for(long i = 0; i < nelements; i++){
        
        long n1, n2, n3;
        
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
    long nsegments;
    
    MeshFile >> nsegments;
    
    SegmentVec.clear();
    
    for(long i = 0; i < nsegments; i++){
        
        long nId1, nId2, eId;
        
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
        
        long ID;
        
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
    
    for(long i = 0; i < gmesh->NElements(); i++){
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
	for(long i = 0 ; i < x.size(); i++){
        
        TPZManVector<REAL,3> coord(3,0.);
        
        coord[0]= x[i];
        coord[1]= y[i];
		
        gmesh->NodeVec()[i].SetCoord(coord);
		gmesh->NodeVec()[i].SetNodeId(i);
	}
	
	//materials ID
    long id;
    const int matTria = 1;
    const int matBoundary = 2;
    
    // Criando Elementos Triangulares
    TPZManVector<long,3> tria(3,0.);
    
	for(long iel = 0; iel < ElemNodeVec.size(); iel++){

        //tira-se "-1" dos IDs pois o MatLab começa a numeração dos nós em "1" e não em "0"
        tria[0] = ElemNodeVec[iel].fNodeId[0] - 1;
        tria[1] = ElemNodeVec[iel].fNodeId[1] - 1;
        tria[2] = ElemNodeVec[iel].fNodeId[2] - 1;

        gmesh->CreateGeoElement(ETriangle, tria, matTria, id);//cria elemento triangular
        
        gmesh->ElementVec()[id]->SetId(iel);
 
	}
    
    //Criando os elementos unidimensionais que compõem todo o contorno
    TPZManVector<long,2> boundary(2,0.);
    
    for(long iel = ElemNodeVec.size(); iel < ElemNodeVec.size() + SegmentVec.size(); iel++){
        
        //tira-se "-1" dos IDs pois o MatLab começa a numeração dos nós em "1" e não em "0"
        boundary[0] = SegmentVec[iel-ElemNodeVec.size()].fId[0] - 1;
        boundary[1] = SegmentVec[iel-ElemNodeVec.size()].fId[1] - 1;
        
        gmesh->CreateGeoElement(EOned, boundary, matBoundary, id);//cria elemento triangular
        
        gmesh->ElementVec()[id]->SetId(iel);
        
    }
    
    gmesh->BuildConnectivity();
    
	return gmesh;
}

///##############################################################
void RefineMesh(TPZCompMesh *cmesh, std::vector<TPZVec<REAL> > &GLvec, const int &MaxLevel){
    
    const int Interpolate = 1; //a solução será interpolada para os filhos
    std::vector<int> sides(3);
    sides[0] = 3;
    sides[1] = 4;
    sides[2] = 5;
    
    int level = 0;
    
    ///distance around grounding line
    const REAL Region = 200000.; ///160 km
    
    const bool IsUniform = true; //itapopo deixar false
    
    while(level < MaxLevel){
    
        level++;
        
        cmesh->Reference()->ResetReference();
        cmesh->LoadReferences();
    
        ///distance around grounding line
        const REAL MaxDistance = Region / std::exp(level-1);//(level*level); //itapopo teste
    
        ///refinando os elementos triangulares necessários
        const long nelem = cmesh->NElements();
        for(long i = 0; i < nelem; i++){
        
            TPZCompEl *compel = cmesh->Element(i);
            if(!compel) continue;
        
            TPZGeoEl *geoel = compel->Reference();
            if(!geoel) DebugStop();
        
            if(IsUniform){
                
                TPZVec<long> SonsIndex;
                cmesh->Divide(i, SonsIndex, Interpolate);
            
            } else {
            
                const int side2D = 6;
                TPZVec<REAL> qsi(2,0.);
                TPZVec<REAL> centerPoint(3,0.);
                geoel->CenterPoint(side2D, qsi);
                geoel->X(qsi, centerPoint);
        
                REAL distance = MaxDistance;
        
                for (long j = 0; j < GLvec.size(); j++) {
            
                    REAL value = ( GLvec[j][0] - centerPoint[0] ) * ( GLvec[j][0] - centerPoint[0] ); // (x2-x1)^2
                    value += ( GLvec[j][1] - centerPoint[1] ) * ( GLvec[j][1] - centerPoint[1] );// (y2-y1)^2
                    value = std::sqrt(value); ///Radius
            
                    if(value < distance){
                        distance = value; //finding the min distance to the Grounding line
                    } ///if
            
                } ///for j / GLvec.size()

                if(distance < MaxDistance){
                    TPZVec<long> SonsIndex;
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
        TPZGeoMesh *gmesh = cmesh->Reference();
        long NElem = gmesh->NElements();
    
        for(long i = 0; i < NElem; i++){
        
            TPZGeoEl * geoel = gmesh->Element(i);
        
            if(!geoel->HasSubElement()){ //não pode ter sido refinado antes
            
                TPZAutoPointer<TPZRefPattern> refp = TPZRefPatternTools::PerfectMatchRefPattern(geoel);
                if(refp)
                {
                    geoel->SetRefPattern(refp);
                    TPZVec<long> SonsIndex;
                    long CompElIndex = geoel->Reference()->Index();
                    cmesh->Divide(CompElIndex, SonsIndex, Interpolate);
                    /**for(long j = 0; j < SonsIndex.NElements();j++){
                     ClosureElementIndex.push_back(SonsIndex[j]);
                     }*/
                }
            
            }
        }
    
        cmesh->AdjustBoundaryElements();
        cmesh->ExpandSolution();
        cmesh->CleanUpUnconnectedNodes();
    
        cmesh->Reference()->BuildConnectivity();
        
    } ///while
}

///##############################################################
void RefineMesh(TPZCompMesh *cmesh,
                /*std::vector<long>*//*std::set<long>*/std::map<long, int>  &ElementIndex,///element index in the geometric mesh and level of refinement
                std::vector<long> &ClosureElementIndex){
    
    ClosureElementIndex.clear();
    
    const int Interpolate = 1; //a solução será interpolada para os filhos
    std::vector<int> sides(3);
    sides[0] = 3;
    sides[1] = 4;
    sides[2] = 5;
    
    cmesh->Reference()->ResetReference();
    cmesh->LoadReferences();
    
    
    ///refinando os elementos triangulares necessários
    //std::set<long>::iterator it;
    std::map<long,int>::iterator it;
    //for(long i = 0; i < ElementIndex.size(); i++){
    for(it=ElementIndex.begin(); it != ElementIndex.end(); it++){
        
        const long geoElIndex = it->first;//*it;//ElementIndex[i];
        const int maxLevel = it->second;
        
        TPZGeoEl *GeoEl = cmesh->Reference()->Element(geoElIndex);
        if(GeoEl->MaterialId() != 1) DebugStop();
        
        
        std::vector<long> index;
        int level = 0;
        
        if(GeoEl->Reference()){
            ///original
            long compElIndex = GeoEl->Reference()->Index();
            //cmesh->Divide(compElIndex, SonsIndex, Interpolate);
            ///original
            
            index.push_back(compElIndex);
            
            while(level < maxLevel){
                
                std::vector<long> SonsIndex;
                SonsIndex.clear();
                
                for(int i = 0; i < index.size(); i++){
                    
                    TPZVec<long> Sons;
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
    long NElem = gmesh->NElements();
    
    for(long i = 0; i < NElem; i++){
        
        TPZGeoEl * GeoEl = gmesh->Element(i);
        
        if(!GeoEl->HasSubElement()){ //não pode ter sido refinado antes
            
            TPZAutoPointer<TPZRefPattern> refp = TPZRefPatternTools::PerfectMatchRefPattern(GeoEl);
            if(refp)
            {
                GeoEl->SetRefPattern(refp);
                TPZVec<long> SonsIndex;
                long CompElIndex = GeoEl->Reference()->Index();
                cmesh->Divide(CompElIndex, SonsIndex, Interpolate);
                for(long j = 0; j < SonsIndex.NElements();j++){
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
    
    ///refinando os elementos triangulares necessários
    std::vector<int> sides(3);
    sides[0] = 3;
    sides[1] = 4;
    sides[2] = 5;
    
    for(long i = 0; i < ElemVec.size(); i++){
        
        long ID = ElemVec[i];
        
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
    
    for(long i = 0; i < NElem; i++){
        
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
    
///refinando os elementos triangulares necessários (vindos do MatLab)
    std::vector<int> sides(3);
    sides[0] = 3;
    sides[1] = 4;
    sides[2] = 5;
    
    for(int i = 0; i < FlagVec.size(); i++){
        
        long ID = FlagVec[i].first;
        
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
    
    for(long i = 0; i < NElem; i++){
        
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

void PrintNewMesh(std::string &MeshName,
                  TPZGeoMesh *mesh){
    
    std::string FullName = toStr(WORKPATH) + MeshName;
    
    std::ofstream MeshFile(FullName.c_str());
    
    if(!MeshFile.is_open()) DebugStop();
    
    MeshFile << std::setprecision(18);
    
///Printing nodes
    long nnodes = mesh->NNodes();
    
    MeshFile << nnodes << "\n";
    
    for(long i = 0; i < nnodes; i++ ){
        TPZVec<REAL> coords(3,0.);
        mesh->NodeVec()[i].GetCoordinates(coords);
        MeshFile << coords[0] << "\t" << coords[1] << "\t" << "\n";
        
    }
    
///Printing elements and Father Index
    long nTotalElements = mesh->NElements();
    
    std::vector<TPZGeoEl *> TriaVec;
    std::vector<TPZGeoEl *> SegmentVec;
    TriaVec.clear();
    SegmentVec.clear();
    
    
    for(long i = 0; i < nTotalElements; i++){
    
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
    
    for(long i = 0; i < TriaVec.size(); i++){
        
        //int RefID = TriaVec[i]->FatherIndex();
        //if(RefID < 0) RefID = TriaVec[i]->Index();
        
        MeshFile << TriaVec[i]->NodeIndex(0)+1 << "\t" << TriaVec[i]->NodeIndex(1)+1 << "\t" << TriaVec[i]->NodeIndex(2)+1 << /*"\t" << RefID+1 <<*/ "\n"; //o 1 eh devido a numercao do MatLab
    }
    
//Printing segments
    MeshFile << SegmentVec.size() << "\n";
    
    for(long i = 0; i < SegmentVec.size(); i++){
        
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
    ///itapopo deve-se garantir de que o numero de nós da malha é igual ao número de connects com funções de forma associado.
    long nnodes = gmesh->NNodes();
    
    DataFile << nnodes << "\n";
    
    DataFile << std::setprecision(18);
    
    ///Printing surface data
    for(long i = 0; i < nnodes; i++ ){
        
        STATE value = sol(0+i*nstate,0);
        
        DataFile << value << "\n";
        
    }
    
    ///Printing base data
    for(long i = 0; i < nnodes; i++ ){
        
        STATE value = sol(1+i*nstate,0);
        
        DataFile << value << "\n";
        
    }
    
    ///Printing bed data
    for(long i = 0; i < nnodes; i++ ){
        
        STATE value = sol(2+i*nstate,0);
        
        DataFile << value << "\n";
        
    }
    
    ///Printing pressure data
    for(long i = 0; i < nnodes; i++ ){
        
        STATE value = sol(3+i*nstate,0);
        
        DataFile << value << "\n";
        
    }
    
    ///Printing temperature data
    for(long i = 0; i < nnodes; i++ ){
        
        STATE value = sol(4+i*nstate,0);
        
        DataFile << value << "\n";
        
    }
    
    ///Printing vx data
    for(long i = 0; i < nnodes; i++ ){
        
        STATE value = sol(5+i*nstate,0);
        
        DataFile << value << "\n";
        
    }
    
    ///Printing vy data
    for(long i = 0; i < nnodes; i++ ){
        
        STATE value = sol(6+i*nstate,0);
        
        DataFile << value << "\n";
        
    }
    
    ///Printing masklevelset data
    for(long i = 0; i < nnodes; i++ ){
        
        STATE value = sol(7+i*nstate,0);
        
        DataFile << value << "\n";
        
    }
    
    DataFile.flush();
    DataFile.close();
    
}
