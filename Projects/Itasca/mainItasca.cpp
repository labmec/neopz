
#include <stdio.h>
#include <pzreal.h>
#include <pzvec.h>
#include <fstream>
#include <list>
#include <map>
#include "TPZCopySolve.h"
#include "pzstack.h"
#include "tpzautopointer.h"
#include "time.h"
#include "TPZGuiInterface.h"

//#include "threadExecuteProdIndex.h"

#include "pzskylstrmatrix.h"

//#include "TSWXPropertyData.h"
//#include "TSWXGraphMesh.h"
//#include "TPBRCreateFractureDefinition.h"
//#include "TPBRFractureVerticalWell.h"
//#include "TPBRFractureHorizontalWell.h"

#include "TPZRefPatternDataBase.h"

//#include "TSWXIP3DGuiIfaceMsg.h"

//#include "TSocket.h"
//#include "TSocketGuiInterface.h"
#include "pzvisualmatrix.h"

//#include "DShow.h"
//#include "winuser.h"
#ifdef SWX_BUILDER_XE3
#endif

//#include "TSocket.h"
//#include "TSocketGuiInterface.h"

#include "TSWXIMPGRIDData.h"
#include "TSWXIMPGRIDParser.h"
#include "TSWXGridDataToPZ.h"
#include "TPZVTKGeoMesh.h"

#include "pzgeoelside.h"
#include "pzgeoelbc.h"
#include "pzmatrix.h"
#include "pzelast3d.h"
#include "pzbndcond.h"

#include "pzfunction.h"
//#include "TSWXPropertyData.h"

#include "TPZParFrontStructMatrix.h"
#include "TPZFrontStructMatrix.h"
#include "TPZFrontSym.h"
#include "pzstepsolver.h"
//#include "TSWXGraphMesh.h"
//#include "TSWXGraphElement.h"
#include "TPZSpStructMatrix.h"
//#include "TSWXGraphElement.h"

#include "SmallerGMeshGenerator.h"

#include "pzinterpolationspace.h"

void CheckHangNodes(TPZGeoMesh * gmesh,
                    std::multimap<TPZGeoElSide,TPZGeoElSide> & restrictedFaceNode_dependencyNodes,
                    std::multimap<TPZGeoElSide,TPZGeoElSide> & restrictedEdgeNode_dependencyNodes,
                    std::set<TPZGeoElSide> & boundaryFaces);

void CheckFace(TPZGeoElSide & gelFace,
               std::multimap<TPZGeoElSide,TPZGeoElSide> & restrictedFaceNode_dependencyNodes,
               std::multimap<TPZGeoElSide,TPZGeoElSide> & restrictedEdgeNode_dependencyNodes,
               std::set<TPZGeoElSide> & boundaryFaces,
               std::set<TPZGeoElSide> & hangFaces);

void EstablishingFaceNodeDependencies(TPZGeoElSide & cornerNode,
                                      TPZGeoElSide & centerNode, TPZGeoElSide & faceSide,
                                      std::multimap<TPZGeoElSide,TPZGeoElSide> & restrictedFaceNode_dependencyNodes,
                                      std::multimap<TPZGeoElSide,TPZGeoElSide> & restrictedEdgeNode_dependencyNodes);

void CheckEdges(std::multimap<TPZGeoElSide,TPZGeoElSide> & restrictedEdgeNode_dependencyNodes,
                std::set<TPZGeoElSide> & boundaryFaces);

void EstablishingEdgeNodeDependencies(TPZGeoElSide & faceSide,
                                      std::set<TPZGeoElSide> & edgeNodesSet,
                                      std::set<TPZGeoElSide> & cornerNodesSet,
                                      std::multimap<TPZGeoElSide,TPZGeoElSide> & restrictedEdgeNode_dependencyNodes);

void GenerateBCsElements(std::set<TPZGeoElSide> & elFacesNoNeigh);

void GenerateItascaInputFile(TPZCompMesh *cmesh, TPZVec<int> &PZGeoElIndex2ItascaID,
                             const std::string &filename,
                             const std::string &displfilename);
void GenerateItascaInputFile(TPZCompMesh *cmesh, TPZVec<int> &PZGeoElIndex2ItascaID,
                             const std::string &filename);

TPZCompMesh * GetCompMesh(TPZGeoMesh * gmesh,
                          std::multimap<TPZGeoElSide,TPZGeoElSide> & restrictedFaceNode_dependencyNodes,
                          std::multimap<TPZGeoElSide,TPZGeoElSide> & restrictedEdgeNode_dependencyNodes);

/// este metodo vai colocar uma camada de faces vizinhos do seed
void GenerateBoundaryFaces(TPZGeoMesh *gmesh, TPZGeoElSide seed, int matid);

/// cria uma malha superficial com matid
void BuildSurfaceMesh(TPZGeoMesh *gmesh, int matid);

/// find an element without a 2d neighbour
TPZGeoElSide FindSeed(TPZGeoMesh *gmesh);

class BodyForceFunction : public TPZFunction<REAL>
{
public:
    
    BodyForceFunction()
    {
        //Utilize o outro construtor
        DebugStop();
    }
    
    BodyForceFunction(REAL rho1, REAL rho2, REAL TVDrhoInterface,
                      REAL porosity, std::map<REAL,REAL> & waterTable_x_z) :
    fRhoAparente1(rho1), fRhoAparente2(rho2), fTVD_Rho_interface(TVDrhoInterface),
    fPorosity(porosity), fWaterTable_x_z(waterTable_x_z)
    {
    }
    
    virtual int NFunctions()
    {
        return 1;
    }
    
    virtual int PolynomialOrder()
    {
        return 1;
    }
    
    virtual void Execute(const TPZVec<REAL> &x, TPZVec<REAL> &f, TPZFMatrix<REAL> &df)
    {
        DebugStop();
    }
    
    
    virtual void Execute(const TPZVec<REAL> &x, TPZVec<REAL> &f)
    {
#ifdef DEBUG
        if(f.NElements() != 3)
        {
            DebugStop();
        }
#endif
        
        f[0] = 0.;
        f[1] = 0.;
        
        REAL gravity = 9.8061999; //m/(s2)
        bool rhoWater = 0.;//kg/(m3)
        
        REAL waterZ = 0.;
        std::map<REAL,REAL>::iterator it = fWaterTable_x_z.lower_bound(x[0] + 1.E-3);
        if(it != fWaterTable_x_z.end())
        {
            waterZ = it->second;
        }
        else
        {
            it--;
            waterZ = it->second;
        }
        
        bool saturated = (x[2] <= waterZ);
        if(saturated)
        {
            rhoWater = 1.;//ton/(m3)
        }
        
        REAL pointTVD = -x[2];
        if(pointTVD <= fTVD_Rho_interface)
        {//1st layer, i.e., rho1
            f[2] = - (fRhoAparente1 + fPorosity*rhoWater) * gravity;
        }
        else
        {//2nd layer, i.e., rho2
            f[2] = - (fRhoAparente2 + fPorosity*rhoWater) * gravity;
        }
        if(fabs(f[2]) < 1.)
        {
            DebugStop();
        }
    }
    
    
private:
    
    REAL fRhoAparente1;
    REAL fRhoAparente2;
    REAL fTVD_Rho_interface;
    REAL fPorosity;
    std::map<REAL,REAL> fWaterTable_x_z;
};

//////////////////METODOS AUXILIARES
double DistBetweenCoords(TPZManVector<REAL> & A, TPZManVector<REAL> & B);

bool GetFaceThatContainsThisNodes(TPZGeoElSide & node0, TPZGeoElSide & node2, TPZGeoElSide & faceSide);

void GetNormal(TPZGeoElSide & face, TPZManVector<REAL> & normal);

bool IsX_direction(TPZManVector<REAL> & normal);
bool IsY_direction(TPZManVector<REAL> & normal);
bool IsZ_direction(TPZManVector<REAL> & normal);

void VTKPointsOnHangNodes_and_Dependencies(TPZGeoMesh * gmesh,
                                           std::multimap<TPZGeoElSide,TPZGeoElSide> & restrictedFaceNode_dependencyNodes,
                                           std::multimap<TPZGeoElSide,TPZGeoElSide> & restrictedEdgeNode_dependencyNodes);

void PostProcess(TPZCompMesh * cmesh);

void GetDomainEnvelopment(TPZGeoMesh * gmesh,
                          double & xRange, double & yRange, double & zRange, TPZVec<REAL> & XYZanchor);

//#define smallerGMesh
int main(int argc, char* argv[])
{
    
    
    //  double vec1[]={160256410.2564103,
    //                 160256410.2564103,
    //                 -160256410.2564102,
    //                 -160256410.2564103};
    //
    //  double sum1 = 0., sum2 = 0.;
    //  for (int i = 0 ; i < 4 ; i++){
    //    sum1 += vec1[i];
    //    sum2 += vec1[3-i];
    //  }
    //  double dif1 = sum1 - sum2;
    //  std::cout.precision(30);
    //  std::cout << std::setprecision(30) << "dif = " << dif1 << std::endl;
    //  getchar();
    //  return 0;
    
    //  int neq = 2;
    //  TPZFYsmpMatrix spM(neq,neq);
    //  double * aa = new double[4];
    //  aa[0] = 1.;
    //  aa[1] = 2.;
    //  aa[2] = 3.;
    //  aa[3] = 4.;
    //
    //  int * ia = new int[3];
    //  ia[0] = 0;
    //  ia[1] = 2;
    //  ia[2] = 4;
    //
    //  int * ja = new int[4];
    //  ja[0] = 0;
    //  ja[1] = 1;
    //  ja[2] = 0;
    //  ja[3] = 1;
    //  spM.SetData(ia,ja,aa);
    //
    //  TPZFMatrix fM(2,1), multM(2,1,0.);
    //  fM.PutVal(0,0,1.);
    //  fM.PutVal(1,0,2.);
    ////  fM.PutVal(1,0,3.);
    ////  fM.PutVal(1,1,4.);
    //
    //  spM.Multiply(fM,multM);
    //
    //  std::ofstream oouu("000sdbvdv.txt");
    //  multM.Print("000chvsdhcv",oouu);
    //  return 0;
    
    
    std::ofstream logfile("../tempo.txt");
    time_t inicio = time(NULL);
    
    //  std::string filename = "hangnodes.dat";
    std::string filename = "../example1 BIG.dat";
    logfile << "1 " << time(NULL) << "\n";
    
#ifdef smallerGMesh
    TPZVec<REAL> XYZanchor(3,0.);
    double Lx = 40.;//52.;
    double Ly = 40.;//32.;
    double Lz = 40.;//20.;
    double Lcharact = 10.;
    TPZGeoMesh * gmesh = SmallerGMeshGenerator::GetGMesh(XYZanchor, Lx, Ly, Lz, Lcharact);
    TPZVec<int> PZGeoElIndex2ItascaID(gmesh->NElements(),-1);
    for(int i = 0; i < gmesh->NElements(); i++)
    {
        PZGeoElIndex2ItascaID[i] = i;
    }
#else
    TSWXGridData * mygriddata = swx::ReadGridDataFile(filename);
    logfile << "2 " << time(NULL) << "\n";
    
    //GeoMesh///////////////////////////////////////
    std::cout << "\n\nGerando malha geometrica...\n";
    std::map<int, std::string> matId2GroupName;
    TPZVec<int> PZGeoElIndex2ItascaID;
    TPZGeoMesh * gmesh = swx::ConvertGridDataToPZ(*mygriddata, matId2GroupName,PZGeoElIndex2ItascaID);
    logfile << "3 " << time(NULL) << "\n";
    std::map<int, std::string>::const_iterator it;
    for (it = matId2GroupName.begin(); it != matId2GroupName.end(); it++) {
        std::cout << "Matid " << it->first << " group name " << it->second.c_str() << std::endl;
    }
#endif
    
    int matid = matId2GroupName.rbegin()->first+10;
    BuildSurfaceMesh(gmesh, matid);
    /*
    int matseed = matId2GroupName.rbegin()->first+10;
    TPZGeoElSide seed = FindSeed(gmesh);
    while (seed) {
        GenerateBoundaryFaces(gmesh, seed, matseed);
        matseed++;
        seed = FindSeed(gmesh);
    }
    */
    
    if(1){
        std::cout << "Gerando VTK...\n";
        std::stringstream vtkfilename;
        vtkfilename << "malhaGeometrica.vtk";
        std::ofstream myfile(vtkfilename.str().c_str());
//        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, myfile);
        std::set<int> myMaterial;
        myMaterial.insert(matid);
        TPZVTKGeoMesh::PrintGMeshVTKmy_material(gmesh, myfile, myMaterial);
    }
    
    std::cout << "Procurando hang nodes...\n";
    std::multimap<TPZGeoElSide, TPZGeoElSide> restrictedFaceNode_dependencyNodes;
    std::multimap<TPZGeoElSide, TPZGeoElSide> restrictedEdgeNode_dependencyNodes;
    std::set<TPZGeoElSide> boundaryFaces;
    CheckHangNodes(gmesh,
                   restrictedFaceNode_dependencyNodes,
                   restrictedEdgeNode_dependencyNodes,
                   boundaryFaces);
    
    {
        VTKPointsOnHangNodes_and_Dependencies(gmesh,
                                              restrictedFaceNode_dependencyNodes,
                                              restrictedEdgeNode_dependencyNodes);
    }
    
    logfile << "4 " << time(NULL) << "\n";
    
    std::cout << "Gerando faces...\n";
    GenerateBCsElements(boundaryFaces);
    logfile << "5 " << time(NULL) << "\n";
    
    /*
     //CompMesh////////////////////////////////////////////////////////////////////
     std::cout << "\nGerando cmesh...\n";
     TPZCompMesh * cmesh = GetCompMesh(gmesh,
     restrictedFaceNode_dependencyNodes,
     restrictedEdgeNode_dependencyNodes);
     
     logfile << "6 " << time(NULL) << "\n";
     
     std::cout << "Rodando o problema... (be patient)\n";
     std::cout << "NEquations = " << cmesh->NEquations() << std::endl;
     TPZAnalysis an(cmesh);
     std::cout << "Criei analysis" << std::endl;
     
     const bool DIRECTSOLVER = false;
     const int nthreads = 32;
     if(DIRECTSOLVER){
     TPZParFrontStructMatrix<TPZFrontSym>  frontMatrix(cmesh);
     an.SetStructuralMatrix(frontMatrix);
     an.StructMatrix()->SetNumThreads(nthreads);
     TPZStepSolver stepS;
     stepS.SetDirect(ECholesky);
     an.SetSolver(stepS);
     }
     else{
     TPZSpStructMatrix StrMatrix(cmesh);
     //TPZSkylineStructMatrix StrMatrix(cmesh);
     an.SetStructuralMatrix(StrMatrix);
     an.StructMatrix()->SetNumThreads(nthreads);
     TPZStepSolver stepS;
     TPZStepSolver precond(NULL);
     //  TPZCopySolve precond(NULL);
     #ifdef smallerGMesh
     REAL tol = 1e-13;
     #else
     REAL tol = 1e-3;
     #endif
     
     stepS.SetCG(cmesh->NEquations(), precond, tol, 0);
     //   stepS.SetGMRES(cmesh->NEquations(), 1000, precond, tol, 0);
     an.SetSolver(stepS);
     std::cout << "------------ Criando a matriz esparsa ------------" << std::endl;
     an.Solver().SetMatrix( an.StructMatrix()->Create() );
     
     TPZStepSolver * cast = dynamic_cast<TPZStepSolver*>(&an.Solver());
     if(cast){
     TPZStepSolver *castPrecond = dynamic_cast<TPZStepSolver*>(cast->PreConditioner());
     if(castPrecond) castPrecond->ShareMatrix( an.Solver() );
     }
     }
     
     if(0){
     TPZFMatrix matrix;
     cmesh->ComputeFillIn(200,matrix);
     VisualMatrixVTK(matrix, "matrix.vtk");
     std::cout << "matrix visual criada\n";
     }
     std::cout << "------------ NEle " << cmesh->NElements();
     
     std::cout << "------------ Entrando no Assemble ------------" << std::endl;
     an.Assemble();
     std::cout << "------------ Saindo do Assemble ------------" << std::endl;
     
     if(0){
     std::cout << "Imprimindo a matriz" << std::endl;
     if (nthreads == 0){
     std::ofstream toto("MatSerialSkylMyOrder.nb");
     an.Solver().Matrix()->Print("serSkylMyOrder=",toto,EMathematicaInput);
     }
     else{
     std::ofstream toto("Mat10bThreadSparse.nb");
     an.Solver().Matrix()->Print("Mat10bThreadSparse=",toto,EMathematicaInput);
     }
     }
     
     
     std::cout << "------------ Entrando no Solve ------------" << std::endl;
     an.Solve();
     std::cout << "------------ Saindo do Solve ------------" << std::endl;
     logfile << "7 " << time(NULL) << "\n";
     
     time_t fim = time(NULL);
     std::cout << "Tempo antes do postproc = " << difftime(fim,inicio) << " s" << std::endl;
     
     bool tamarindo = false;//NAO GERA VTK
     if(tamarindo){
     std::cout << "Pos-processamentos... (be patient)\n";
     PostProcess(cmesh);
     }
     
     GenerateItascaInputFile(cmesh,PZGeoElIndex2ItascaID, "StressVector.txt");
     //  GenerateItascaInputFile(cmesh,PZGeoElIndex2ItascaID,"StressVector.txt","DisplacementVector.txt");
     
     delete cmesh;
     */
    delete gmesh;
    
    std::cout << "\n\nDONE!!!\n\n";
    
    time_t fim = time(NULL);
    double diferenca = difftime(fim,inicio);
    std::cout << "tempo : " << diferenca << " s" << "\n";
    std::cout << "inicio : " << inicio << "\n";
    std::cout << "final : " << fim << "\n";
    std::cout.flush();
    
    //  int stop;
    //  std::cin >> stop;
    return 0;
}

/*
 void GenerateItascaInputFile(TPZCompMesh *cmesh,
 TPZVec<int> &PZGeoElIndex2ItascaID,
 const std::string &filename,
 const std::string &displfilename){
 //stress
 std::ofstream myfile(filename.c_str());
 myfile << "def SimworxStresses\n";
 const int nel = cmesh->NElements();
 TPZManVector<REAL,3> qsi(3);
 TPZManVector<REAL,6> stressVec(6);
 std::ofstream myfileTracaoZ("tracaoZZ.txt");
 int count = 0;
 for(int iel = 0; iel < nel; iel++){
 TPZCompEl *cel = cmesh->ElementVec()[iel];
 if(!cel) continue;
 TPZGeoEl *gel = cel->Reference();
 if(gel->Dimension() != 3) continue;
 gel->CenterPoint(gel->NSides()-1,qsi);
 cel->Solution(qsi,TPZElasticity3D::EStressVector,stressVec);
 const int ID = PZGeoElIndex2ItascaID[ gel->Index() ];
 for(int is = 0; is < 6; is++){
 const int pos = is+1;
 myfile << "sig(" << ID << "," << pos << ")=" << stressVec[is] << "\n";
 }
 if(stressVec[3] > 0) myfileTracaoZ << ID << "\t" << stressVec[3] << "\n";
 count++;
 }//for iel
 myfile << "end\n\n";
 myfile << "@SimworxStresses\n";
 
 std::cout << "\nPZGeoElIndex2ItascaID.NElements() = " <<
 PZGeoElIndex2ItascaID.NElements() << "\n";
 std::cout <<  "GenerateItascaInputFile::count = " << count << "\n";
 std::cout.flush();
 
 //Displacement
 std::cout << "\nDisplacement file started!!!\n";
 
 std::ofstream displfile(displfilename);
 displfile << "def SimworxDisplacement\n";
 
 int nelements = cmesh->NElements();
 int InitialElIndex = 0;
 
 std::set<int> nodesProcessed;
 
 for(int el = 0; el < nelements; el++)
 {
 TPZCompEl * cel = cmesh->ElementVec()[el];
 TPZGeoEl * gel = cel->Reference();
 int nnodes = gel->NNodes();
 
 for(int n = 0; n < nnodes; n++)
 {
 TPZGeoNode * nd = gel->NodePtr(n);
 if( nodesProcessed.find(nd->Id()) != nodesProcessed.end())
 {
 continue;
 }
 else
 {
 nodesProcessed.insert(nd->Id());
 }
 
 TPZManVector<REAL,3> ndcoord(3);
 TPZManVector<REAL,3> ndqsi(3);
 nd->GetCoordinates(ndcoord);
 gel->ComputeXInverse(ndcoord,ndqsi);
 TPZManVector<REAL,3> displVec(3);
 cel->Solution(qsi,TPZElasticity3D::EDisplacement,displVec);
 
 displfile << "desloc(" << nd->Id() << ")=vector(" << displVec[0] << "," << displVec[1] << "," << displVec[2] << ")\n";
 }
 }
 displfile << "end\n\n";
 displfile << "@SimworxDisplacement\n";
 
 std::cout << "\nDisplacement file finished!!!\n";
 }//void
 void GenerateItascaInputFile(TPZCompMesh *cmesh,
 TPZVec<int> &PZGeoElIndex2ItascaID,
 const std::string &filename)
 {
 //stress
 std::ofstream myfile(filename);
 myfile << "def SimworxStresses\n";
 const int nel = cmesh->NElements();
 TPZManVector<REAL,3> qsi(3);
 TPZManVector<REAL,6> stressVec(6);
 std::ofstream myfileTracaoZ("tracaoZZ.txt");
 int count = 0;
 for(int iel = 0; iel < nel; iel++){
 TPZCompEl *cel = cmesh->ElementVec()[iel];
 if(!cel) continue;
 TPZGeoEl *gel = cel->Reference();
 if(gel->Dimension() != 3) continue;
 gel->CenterPoint(gel->NSides()-1,qsi);
 cel->Solution(qsi,TPZElasticity3D::EStressVector,stressVec);
 const int ID = PZGeoElIndex2ItascaID[ gel->Index() ];
 for(int is = 0; is < 6; is++){
 const int pos = is+1;
 myfile << "sig(" << ID << "," << pos << ")=" << stressVec[is] << "\n";
 }
 if(stressVec[3] > 0) myfileTracaoZ << ID << "\t" << stressVec[3] << "\n";
 count++;
 }//for iel
 myfile << "end\n\n";
 myfile << "@SimworxStresses\n";
 
 std::cout << "\nPZGeoElIndex2ItascaID.NElements() = " <<
 PZGeoElIndex2ItascaID.NElements() << "\n";
 std::cout <<  "GenerateItascaInputFile::count = " << count << "\n";
 std::cout.flush();
 
 }//void
 
 void PostProcess(TPZCompMesh * cmesh)
 {
 TPZVec<std::string> nodalVarIndex(4);
 nodalVarIndex[0] = "Displacement";
 nodalVarIndex[1] = "StressX";
 nodalVarIndex[2] = "StressY";
 nodalVarIndex[3] = "StressZ";
 
 TPZVec<std::string> cellVarIndex(0);
 
 TSWXGraphMesh graphMesh;
 TSWXGraphElement graphEl(0);
 std::cout << "######### Entering Generate VTK Data ############" << std::endl;
 graphEl.GenerateVTKData(cmesh,3,0.,nodalVarIndex,cellVarIndex,graphMesh);
 
 std::ofstream file("ITASCA.vtk");
 std::cout << "---------- Generating VTK file --------------" << std::endl;
 graphMesh.ToParaview(file);
 }
 */

void CheckHangNodes(TPZGeoMesh * gmesh,
                    std::multimap<TPZGeoElSide,TPZGeoElSide> & restrictedFaceNode_dependencyNodes,
                    std::multimap<TPZGeoElSide,TPZGeoElSide> & restrictedEdgeNode_dependencyNodes,
                    std::set<TPZGeoElSide> & boundaryFaces)
{
    std::set<TPZGeoElSide> hangFaces;
    
    restrictedFaceNode_dependencyNodes.clear();
    hangFaces.clear();
    
    // percorrendo os elementos 3D da malha
    int nelem = gmesh->NElements();
    for(int el = 0; el < nelem; el++)
    {
        TPZGeoEl * gel = gmesh->ElementVec()[el];
        if(!gel || gel->Dimension() != 3)
        {
            continue;
        }
        
        // percorrendo as faces do elemento
        for(int s = gel->NNodes(); s < gel->NSides() - 1; s++)
        {
            TPZGeoElSide gelFace(gel, s);
            if(gelFace.Dimension() != 2)
            {
                continue;
            }
            if(gelFace.Neighbour() != gelFace)
            {
                // estamos interessados somente em faces sem vizinho,
                // pois seria Contorno ou HangFace.
                continue;
            }
            CheckFace(gelFace,
                      restrictedFaceNode_dependencyNodes,
                      restrictedEdgeNode_dependencyNodes,
                      boundaryFaces,
                      hangFaces);
        }
    }
    CheckEdges(restrictedEdgeNode_dependencyNodes, boundaryFaces);
}

void CheckFace(TPZGeoElSide & gelFace,
               std::multimap<TPZGeoElSide,TPZGeoElSide> & restrictedFaceNode_dependencyNodes,
               std::multimap<TPZGeoElSide,TPZGeoElSide> & restrictedEdgeNode_dependencyNodes,
               std::set<TPZGeoElSide> & boundaryFaces,
               std::set<TPZGeoElSide> & hangFaces)
{
    if(hangFaces.find(gelFace) == hangFaces.end())
    {
        // se logo abaixo for detectado que eh hang face, serah tirada desta estrutura!
        boundaryFaces.insert(gelFace);
    }
    else
    {
        return;
    }
    
    // coordenadas do centro da face, em (x,y,z)
    TPZManVector<REAL> XfaceCenter(3,0.);
    gelFace.CenterX(XfaceCenter);
    
    // capturando os lados desta face.
    TPZStack<int> smallsides;
    gelFace.Element()->LowerDimensionSides(gelFace.Side(), smallsides);
    
    // percorrendo os nohs desta face
    for(int ssVecPos = 0; ssVecPos < smallsides.NElements(); ssVecPos++)
    {
        int smallS = smallsides[ssVecPos];
        TPZGeoElSide nodeSide(gelFace.Element(), smallS);
        if(nodeSide.Dimension() != 0)
        {
            continue;
        }
        
        // Percorrendo vizinhos pelos nohs da face
        TPZGeoElSide neighNodeSide(nodeSide.Neighbour());
        while (neighNodeSide != nodeSide)
        {
            // verificando se o vizinho tem
            // algum noh que coincide com XfaceCenter.
            // percorrendo os vizinhos pelos nos da face
            for(int nn = 0; nn < neighNodeSide.Element()->NNodes(); nn++)
            {
                TPZGeoElSide centerNeighNode(neighNodeSide.Element(), nn);
                
                TPZGeoNode * neighNodePtr = neighNodeSide.Element()->NodePtr(nn);
                TPZManVector<REAL> coords(3, 0.);
                neighNodePtr->GetCoordinates(coords);
                double difToCenter = DistBetweenCoords(XfaceCenter, coords);
                if(difToCenter < 1.e-3)
                {
                    // eh hangFace, portanto serah removida da lista de candidatos
                    {
                        std::set<TPZGeoElSide>::iterator myIt = boundaryFaces.find(gelFace);
                        if(myIt != boundaryFaces.end())
                        {
                            boundaryFaces.erase(myIt);
                        }
                        
                        hangFaces.insert(gelFace);
                        
                        TPZGeoElSide whichFace;
                        if(GetFaceThatContainsThisNodes(neighNodeSide, centerNeighNode, whichFace))
                        {
                            std::set<TPZGeoElSide>::iterator itSubFace = boundaryFaces.find(whichFace);
                            if(itSubFace != boundaryFaces.end())
                            {
                                boundaryFaces.erase(itSubFace);
                            }
                            
                            hangFaces.insert(whichFace);
                            EstablishingFaceNodeDependencies(neighNodeSide, centerNeighNode, whichFace,
                                                             restrictedFaceNode_dependencyNodes,
                                                             restrictedEdgeNode_dependencyNodes);
                        }
                        else
                        {
                            DebugStop();
                        }
                    }
                }
            }
            neighNodeSide = neighNodeSide.Neighbour();
        }
    }
}

void EstablishingFaceNodeDependencies(TPZGeoElSide & cornerNode,
                                      TPZGeoElSide & centerNode, TPZGeoElSide & faceSide,
                                      std::multimap<TPZGeoElSide,TPZGeoElSide> & restrictedFaceNode_dependencyNodes,
                                      std::multimap<TPZGeoElSide,TPZGeoElSide> & restrictedEdgeNode_dependencyNodes)
{
    TPZStack<int> smallsides;
    faceSide.Element()->LowerDimensionSides(faceSide.Side(), smallsides);
    
    int edgesFound = 0;
    
    for(int s = 0; s < smallsides.NElements(); s++)
    {
        TPZGeoElSide faceSmallSide(faceSide.Element(), smallsides[s]);
        if(faceSmallSide.Dimension() == 0)
        { // eh noh
            if(faceSmallSide != cornerNode && faceSmallSide != centerNode)
            { // o noh eh de aresta
                edgesFound++;
                restrictedFaceNode_dependencyNodes.insert(std::make_pair(centerNode,faceSmallSide));
                restrictedEdgeNode_dependencyNodes.insert(std::make_pair(faceSmallSide,cornerNode));
            }
        }
    }
    
    // #ifdef DEBUG
    if(edgesFound != 2)
    {
        DebugStop();
    }
    // #endif
}

void CheckEdges(std::multimap<TPZGeoElSide,TPZGeoElSide> & restrictedEdgeNode_dependencyNodes,
                std::set<TPZGeoElSide> & boundaryFaces)
{
    std::set<TPZGeoElSide> hangFaces;
    
    std::set<TPZGeoElSide>::iterator itfor;
    for(itfor = boundaryFaces.begin(); itfor != boundaryFaces.end(); itfor++)
    {
        TPZGeoElSide gelFace(*itfor);
        
        TPZStack<int> smallsides;
        gelFace.Element()->LowerDimensionSides(gelFace.Side(), smallsides);
        
        ///Preenchendo vetor de centros das arestas desta face
        std::vector<TPZManVector<REAL> > edgesCenters(0);
        for(int ssVecPos = 0; ssVecPos < smallsides.NElements(); ssVecPos++)
        {
            int smallS = smallsides[ssVecPos];
            TPZGeoElSide edgeSide(gelFace.Element(), smallS);
            if(edgeSide.Dimension() != 1)
            {
                continue;
            }
            // Coordenadas em X do centro da aresta
            TPZManVector<REAL> XedgeCenter(3, 0.);
            edgeSide.CenterX(XedgeCenter);
            edgesCenters.push_back(XedgeCenter);
        }
        
        ///Percorrendo vizinhos por aresta
        for(int ssVecPos = 0; ssVecPos < smallsides.NElements(); ssVecPos++)
        {
            int smallS = smallsides[ssVecPos];
            TPZGeoElSide edgeSide(gelFace.Element(), smallS);
            if(edgeSide.Dimension() != 1)
            {
                continue;
            }
            TPZGeoElSide neighEdgeSide(edgeSide.Neighbour());
            while (neighEdgeSide != edgeSide)
            {
                TPZGeoEl * neighGel = neighEdgeSide.Element();
                
                ///percorrendo as faces do vizinho, procurando por
                // aquelas que aparecem na lista de faces candidatas
                for(int ns = neighGel->NNodes(); ns < neighGel->NSides() - 1; ns++)
                {
                    TPZGeoElSide neighFaceSide(neighGel, ns);
                    if(neighFaceSide.Dimension() < 2)
                    {
                        continue;
                    }
                    else if(neighFaceSide.Dimension() > 2)
                    {
                        break;
                    }
                    
                    // Se esta face jah foi excluida, vamos para a proxima!
                    std::set<TPZGeoElSide>::iterator it = hangFaces.find(neighFaceSide);
                    if(it != hangFaces.end())
                    {
                        continue;
                    }
                    
                    it = boundaryFaces.find(neighFaceSide);
                    if(it != boundaryFaces.end())
                    { // esta face estah na lista de faces candidatas
                        TPZStack<int> neighsmallsides;
                        neighGel->LowerDimensionSides(ns, neighsmallsides);
                        
                        int nnodesOnEdgeCenter = 0;
                        std::set<TPZGeoElSide> cornerNodesSet, edgeNodesSet;
                        for(int p = 0; p < neighsmallsides.NElements(); p++)
                        {
                            ///percorrendo nohs desta face, verificando se algum
                            // tem coordenada em alguns dos edgesCenters
                            TPZGeoElSide neighNodeSide(neighGel, neighsmallsides[p]);
                            if(neighNodeSide.Dimension() != 0)
                            {
                                continue;
                            }
                            
                            TPZGeoNode * neighNodePtr = neighNodeSide.Element()->NodePtr(neighsmallsides[p]);
                            TPZManVector<REAL> coords(3, 0.);
                            neighNodePtr->GetCoordinates(coords);
                            
                            bool thisNodeIsOnNeighEdge = false;
                            for(int edgeC = 0; edgeC < edgesCenters.size(); edgeC++)
                            {
                                double difToCenter = DistBetweenCoords(edgesCenters[edgeC],coords);
                                if(difToCenter < 1.e-3)
                                {
                                    thisNodeIsOnNeighEdge = true;
                                    
                                    nnodesOnEdgeCenter++;
                                    
                                    if(nnodesOnEdgeCenter > 1)
                                    {
                                        // eh hangFace
                                        hangFaces.insert(gelFace);
                                        hangFaces.insert(neighFaceSide);
                                    }
                                    
                                    edgeNodesSet.insert(neighNodeSide);
                                }
                            }
                            if(thisNodeIsOnNeighEdge == false)
                            {
                                cornerNodesSet.insert(neighNodeSide);
                            }
                        }
                        if(edgeNodesSet.size() == 2 && cornerNodesSet.size() == 2)
                        {
                            EstablishingEdgeNodeDependencies(neighFaceSide,
                                                             edgeNodesSet, cornerNodesSet,
                                                             restrictedEdgeNode_dependencyNodes);
                        }
                    }
                }
                neighEdgeSide = neighEdgeSide.Neighbour();
            }
        }
    }
    
    ///excluindo as hangFaces
    for(itfor = hangFaces.begin(); itfor != hangFaces.end(); itfor++)
    {
        std::set<TPZGeoElSide>::iterator itToExclude = boundaryFaces.find(*itfor);
        if(itToExclude != boundaryFaces.end())
        {
            boundaryFaces.erase(itToExclude);
        }
    }
}

void EstablishingEdgeNodeDependencies(TPZGeoElSide & faceSide,
                                      std::set<TPZGeoElSide> & edgeNodesSet,
                                      std::set<TPZGeoElSide> & cornerNodesSet,
                                      std::multimap<TPZGeoElSide,TPZGeoElSide> & restrictedEdgeNode_dependencyNodes)
{
    // #ifdef DEBUG
    if(edgeNodesSet.size() != 2 || cornerNodesSet.size() != 2)
    {
        DebugStop();
    }
    // #endif
    
    int edgesPaired = 0;
    
    TPZStack<int> smallsides;
    faceSide.Element()->LowerDimensionSides(faceSide.Side(), smallsides);
    for(int s = 0; s < smallsides.NElements(); s++)
    {
        TPZGeoElSide faceSmallSide(faceSide.Element(), smallsides[s]);
        if(faceSmallSide.Dimension() == 1)
        { // eh aresta
            int node0LocIndex = faceSmallSide.GelLocIndex(0);
            TPZGeoElSide node0Side(faceSide.Element(), node0LocIndex);
            bool node0IsOn_edgeNodesSet = (edgeNodesSet.find(node0Side) != edgeNodesSet.end());
            bool node0IsOn_cornerNodesSet = (cornerNodesSet.find(node0Side) != cornerNodesSet.end());
            
            int node1LocIndex = faceSmallSide.GelLocIndex(1); TPZGeoElSide node1Side(faceSide.Element(), node1LocIndex);
            bool node1IsOn_edgeNodesSet = (edgeNodesSet.find(node1Side) != edgeNodesSet.end());
            bool node1IsOn_cornerNodesSet = (cornerNodesSet.find(node1Side) != cornerNodesSet.end());
            
            if(node0IsOn_edgeNodesSet && node1IsOn_cornerNodesSet)
            { // o noh eh de aresta
                edgesPaired++;
                restrictedEdgeNode_dependencyNodes.insert(std::make_pair(node0Side,node1Side));
            }
            else if(node0IsOn_cornerNodesSet && node1IsOn_edgeNodesSet)
            { // o noh eh de aresta
                edgesPaired++;
                restrictedEdgeNode_dependencyNodes.insert(std::make_pair(node1Side,node0Side));
            }
        }
    }
    
    // #ifdef DEBUG
    if(edgesPaired != 2)
    {
        DebugStop();
    }
    // endif
}

void GenerateBCsElements(std::set<TPZGeoElSide> & elFacesNoNeigh)
{
    int normalX_bcMatId = -10;
    int normalY_bcMatId = -20;
    int normalZ_bottom_bcMatId = -30;
    int top_bcMatId = -40;
    
    std::set<TPZGeoElSide>::iterator it;
    for(it = elFacesNoNeigh.begin(); it != elFacesNoNeigh.end(); it++)
    {
        TPZGeoElSide gelFace(*it);
        if(gelFace.Dimension() != 2)
        {
            continue;
        }
        
        TPZGeoElSide neighFace(gelFace.Neighbour());
        if(neighFace != gelFace)
        {
            // era pra ser face sem vizinho
            DebugStop();
        }
        else
        {
            int s = gelFace.Side();
            TPZGeoEl * gel = gelFace.Element();
            
            TPZManVector<REAL> normalVec(3, 0);
            GetNormal(gelFace, normalVec);
            
            if(IsX_direction(normalVec))
            {
                gel->CreateBCGeoEl(s, normalX_bcMatId);
            }
            else if(IsY_direction(normalVec))
            {
                gel->CreateBCGeoEl(s, normalY_bcMatId);
            }
            else if(IsZ_direction(normalVec))
            {
                // verificando se eh face da base ou topo
                TPZManVector<REAL> QSIcentrVolum(3, 0.), QSIcentrFace(3, 0.);
                gel->CenterPoint(gel->NSides() - 1, QSIcentrVolum);
                gel->CenterPoint(s, QSIcentrFace);
                
                TPZManVector<REAL> XcentrVolum(3, 0.), XcentrFace(3, 0.);
                gel->X(QSIcentrVolum, XcentrVolum);
                gel->X(QSIcentrFace, XcentrFace);
                
                if(XcentrVolum[2] > XcentrFace[2])
                { // eh face do fundo
                    gel->CreateBCGeoEl(s, normalZ_bottom_bcMatId);
                }
                else
                { // eh face do topo com normal Z
                    gel->CreateBCGeoEl(s, top_bcMatId);
                }
            }
            else
            { // eh face inclinada do topo
                gel->CreateBCGeoEl(s, top_bcMatId);
            }
        }
    }
}

TPZCompMesh * GetCompMesh(TPZGeoMesh * gmesh,
                          std::multimap<TPZGeoElSide,TPZGeoElSide> & restrictedFaceNode_dependencyNodes,
                          std::multimap<TPZGeoElSide,TPZGeoElSide> & restrictedEdgeNode_dependencyNodes)
{
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    
    cmesh->SetDimModel(3);
    cmesh->SetDefaultOrder(1);
    
    TPZFNMatrix<9> k(3,3,0.), f(3,1,0.);
    int dirich = 0;//fixed at the bottom of domain
    int dirichDir = 3;//roller condition at the side faces
    
    REAL rho1 = 1.6E3; //(kg)/(m3) (TVDini = 0 , TVDfin = 27)
    REAL rho2 = 2.0E3;  //(kg)/(m3) (TVD > 27)
    REAL TVD_RhoInterface = 27.;//m
    REAL porosity = 0.3;
    //
    TPZManVector<REAL> force(3,0.);
    std::map<REAL,REAL> waterTable_x_z;
    waterTable_x_z[0.0]   =  -1.8;
    waterTable_x_z[66.8]  =  -2.0;
    waterTable_x_z[80.0]  =  -8.6;
    waterTable_x_z[94.6]  =  -8.6;
    waterTable_x_z[110.4] =  -9.3;
    waterTable_x_z[125.0] =  -9.6;
    waterTable_x_z[135.7] = -14.0;
    waterTable_x_z[147.5] = -17.0;
    waterTable_x_z[150.5] = -17.0;
    waterTable_x_z[165.3] = -15.6;
    waterTable_x_z[178.7] = -13.6;
    waterTable_x_z[190.0] = -12.5;
    
    TPZAutoPointer<TPZFunction<REAL> > bodyForceFunc =
    new BodyForceFunction(rho1,rho2,TVD_RhoInterface,porosity,waterTable_x_z);
    
    REAL young = 1.E8;//Pa
    REAL poisson = 0.3;
    REAL prestressXX = 0.;// <<<<<<<<< Effect. horiz. stress/Effect. vert. stress = 0.8
    REAL prestressYY = 0.;// <<<<<<<<< Effect. horiz. stress/Effect. vert. stress = 0.8
    REAL prestressZZ = 0.;// <<<<<<<<< Effect. horiz. stress/Effect. vert. stress = 0.8
    
    TPZMaterial  *mat3Dlin = NULL;
    int Nelast3Dmat = 24;
    for(int m = 1; m <= Nelast3Dmat ; m++)
    {
        ///Elast3D Materials
        force[2] = -rho2 * 9.81;
        mat3Dlin = new TPZElasticity3D(m, young, poisson, force, prestressXX, prestressYY, prestressZZ);
        //    mat3Dlin->SetForcingFunction(bodyForceFunc);  //AQUICAJUX
        cmesh->InsertMaterialObject(mat3Dlin);
    }
    std::cout << "Criando materiais ---- " << std::endl;
    ///BCs
    k.Zero();
    f.Zero();
    f(0,0) = 1.;//blocked X displacement
    TPZBndCond * blockedXmat = new TPZBndCond(mat3Dlin,-10, dirichDir, k, f);
    cmesh->InsertMaterialObject(blockedXmat);
    
    k.Zero();
    f.Zero();
    f(1,0) = 1.;//blocked Y displacement
    TPZBndCond * blockedYmat = new TPZBndCond(mat3Dlin,-20, dirichDir, k, f);
    cmesh->InsertMaterialObject(blockedYmat);
    
    k.Zero();
    f.Zero();//null dirichlet at bottom face
    TPZBndCond * blockedZmat_bottom = new TPZBndCond(mat3Dlin,-30, dirich, k, f);
    cmesh->InsertMaterialObject(blockedZmat_bottom);
    
    cmesh->SetAllCreateFunctionsContinuous();
    
    const bool tamarindo = false;
    if(tamarindo){
        std::set<int> targetMatId;
        targetMatId.insert(13);
        cmesh->AutoBuild(targetMatId); ///tamarindo
    }
    else{
        std::cout << "###### Realizando AutoBuild ######" << std::endl;
        cmesh->AutoBuild();
    }
    cmesh->InitializeBlock();
    
    {// soh escopo
        cmesh->LoadReferences();
        
        std::multimap<TPZGeoElSide,TPZGeoElSide>::iterator it;
        
        ///Connects da aresta para os vertices
        for (it = restrictedEdgeNode_dependencyNodes.begin();
             it != restrictedEdgeNode_dependencyNodes.end();
             it++)
        {
            int connectfrom = (it->first).Element()->Reference()->ConnectIndex((it->first).Side());
            int connectto = (it->second).Element()->Reference()->ConnectIndex((it->second).Side());
            TPZConnect &from = cmesh->ConnectVec()[connectfrom];
            TPZConnect &to = cmesh->ConnectVec()[connectto];
            int seqnumfrom = from.SequenceNumber();
            int seqnumto = to.SequenceNumber();
            int blocksize = cmesh->Block().Size(seqnumfrom);
            TPZFNMatrix<200> depmat(blocksize, blocksize);
            depmat.Identity();
            depmat *= 0.5;
            from.AddDependency(connectfrom, connectto, depmat, 0, 0, blocksize, blocksize);
        }
        
        ///Connects da face para as arestas
        for (it = restrictedFaceNode_dependencyNodes.begin();
             it != restrictedFaceNode_dependencyNodes.end();
             it++)
        {
            int connectfrom = (it->first).Element()->Reference()->ConnectIndex((it->first).Side());
            int connectto = (it->second).Element()->Reference()->ConnectIndex((it->second).Side());
            TPZConnect &from = cmesh->ConnectVec()[connectfrom];
            TPZConnect &to = cmesh->ConnectVec()[connectto];
            int seqnumfrom = from.SequenceNumber();
            int seqnumto = to.SequenceNumber();
            int blocksize = cmesh->Block().Size(seqnumfrom);
            TPZFNMatrix<200> depmat(blocksize, blocksize);
            depmat.Identity();
            depmat *= 0.25;
            from.AddDependency(connectfrom, connectto, depmat, 0, 0, blocksize, blocksize);
        }
        
        
        
        cmesh->InitializeBlock();
    }
    
    return cmesh;
}




//////////////////METODOS AUXILIARES

double DistBetweenCoords(TPZManVector<REAL> & A, TPZManVector<REAL> & B)
{
    if(A.NElements() != B.NElements())
    {
        DebugStop();
    }
    
    double answ = 0.;
    for(int c = 0; c < A.NElements(); c++)
    {
        answ += (B[c] - A[c])*(B[c] - A[c]);
    }
    answ = sqrt(answ);
    
    return answ;
}

bool GetFaceThatContainsThisNodes(TPZGeoElSide & node0,
                                  TPZGeoElSide & node2,
                                  TPZGeoElSide & faceSide)
{
    TPZGeoEl * gel = node0.Element();
    if(gel->Dimension() != 3)
    {
        return false;
    }
    
#ifdef DEBUG
    if(gel != node2.Element())
    {
        DebugStop();
    }
#endif
    
    for(int s = gel->NNodes(); s < gel->NSides() - 1; s++)
    {
        TPZGeoElSide gelSide(gel, s);
        
        TPZStack<int> smallsides;
        gel->LowerDimensionSides(s, smallsides);
        
        bool foundNode0 = false;
        bool foundNode2 = false;
        
        for(int n = 0; n < smallsides.NElements(); n++)
        {
            TPZGeoElSide smallSide(gel, smallsides[n]);
            if(node0 == smallSide)
            {
                foundNode0 = true;
            }
            if(node2 == smallSide)
            {
                foundNode2 = true;
            }
            
            if(foundNode0 && foundNode2)
            {
                if(gelSide.Dimension() == 1)
                {
                    // Este metodo foi concebido para achar a UNICA face
                    // que contem os nohs fornecidos, que deveriam ser nohs diagonais.
                    // Aqui foi encontrado uma aresta que contem os 2 nohs, o que
                    // implica em 2 faces que contem os respectivos nohs, sendo portanto
                    // um uso indevido deste metodo!
                    DebugStop();
                }
                else if(gelSide.Dimension() == 2)
                {
                    faceSide = gelSide;
                    return true; // face encontrada
                }
            }
        }
    }
    
    return false; // face nao encontrada
}

void GetNormal(TPZGeoElSide & face, TPZManVector<REAL> & normal)
{
    if(face.Dimension() != 2)
    {
        DebugStop();
    }
    
    TPZManVector<REAL> centr(2);
    TPZFNMatrix<9> jac(2, 2), ax(2, 3), jacinv(2, 2);
    REAL det;
    face.Element()->CenterPoint(face.Side(), centr);
    face.Jacobian(centr, jac, ax, det, jacinv);
    
    normal.Resize(3);
    normal[0] = -(ax(0,2) * ax(1,1)) + (ax(0,1) * ax(1,2));
    normal[1] = +(ax(0,2) * ax(1,0)) - (ax(0,0) * ax(1,2));
    normal[2] = -(ax(0,1) * ax(1,0)) + (ax(0,0) * ax(1,1));
    
    REAL norm = 0.;
    for(int i = 0; i < 3; i++)
    {
        norm += normal[i] * normal[i];
    }
    norm = sqrt(norm);
    for(int i = 0; i < 3; i++)
    {
        normal[i] = normal[i] / norm;
    }
}

bool IsX_direction(TPZManVector<REAL> & normal)
{
    if(fabs(normal[0]) > 0.99)
    {
        return true;
    }
    return false;
}

bool IsY_direction(TPZManVector<REAL> & normal)
{
    if(fabs(normal[1]) > 0.99)
    {
        return true;
    }
    return false;
}

bool IsZ_direction(TPZManVector<REAL> & normal)
{
    if(fabs(normal[2]) > 0.99)
    {
        return true;
    }
    return false;
}



#include "pzgnode.h"
void VTKPointsOnHangNodes_and_Dependencies(TPZGeoMesh * gmesh,
                                           std::multimap<TPZGeoElSide,TPZGeoElSide> & restrictedFaceNode_dependencyNodes,
                                           std::multimap<TPZGeoElSide,TPZGeoElSide> & restrictedEdgeNode_dependencyNodes)
{
    std::cout << "\n\nVTKPointsOnHangNodes_and_Dependencies:\n\n";
    
    std::multimap<TPZGeoElSide,TPZGeoElSide>::iterator it;
    
    int faceNodeMatId   = -500;
    int edgeNodeMatId1  = -501;
    int edgeNodeMatId2  = -502;
    int cornerNodeMatId = -503;
    
    
    for(it  = restrictedFaceNode_dependencyNodes.begin();
        it != restrictedFaceNode_dependencyNodes.end();
        it++)
    {
        TPZGeoElSide faceNode = it->first;
        TPZGeoElSide edgeNode = it->second;
        
        faceNode.Element()->CreateBCGeoEl(faceNode.Side(), faceNodeMatId);
        edgeNode.Element()->CreateBCGeoEl(edgeNode.Side(), edgeNodeMatId1);
    }
    
    for(it  = restrictedEdgeNode_dependencyNodes.begin();
        it != restrictedEdgeNode_dependencyNodes.end();
        it++)
    {
        TPZGeoElSide edgeNode = it->first;
        TPZGeoElSide cornerNode = it->second;
        
        edgeNode.Element()->CreateBCGeoEl(edgeNode.Side(), edgeNodeMatId2);
        cornerNode.Element()->CreateBCGeoEl(cornerNode.Side(), cornerNodeMatId);
    }
    
    std::ofstream out("Dependencies.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
}


void GetDomainEnvelopment(TPZGeoMesh * gmesh,
                          double & xRange, double & yRange, double & zRange, TPZVec<REAL> & XYZanchor)
{
    double minX, minY, minZ;
    minX = minY = minZ = +1.E20;
    
    double maxX, maxY, maxZ;
    maxX = maxY = maxZ = -1.E20;
    
    int nnodes = gmesh->NNodes();
    for(int n = 0; n < nnodes; n++)
    {
        double nX = gmesh->NodeVec()[n].Coord(0);
        double nY = gmesh->NodeVec()[n].Coord(1);
        double nZ = gmesh->NodeVec()[n].Coord(2);
        
        minX = std::min(minX,nX);
        minY = std::min(minY,nY);
        minZ = std::min(minZ,nZ);
        
        maxX = std::max(maxX,nX);
        maxY = std::max(maxY,nY);
        maxZ = std::max(maxZ,nZ);
    }
    
    XYZanchor.Resize(3);
    XYZanchor[0] = minX;
    XYZanchor[1] = minY;
    XYZanchor[2] = minZ;
    
    xRange = maxX - minX;
    yRange = maxY - minY;
    zRange = maxZ - minZ;
}

/// este metodo vai colocar uma camada de faces vizinhos do seed
void GenerateBoundaryFaces(TPZGeoMesh *gmesh, TPZGeoElSide seed, int matid)
{
    if (seed.Dimension() != 2) {
        DebugStop();
    }
    long numelcreated = 0;
    std::set<TPZGeoElSide> workingset;
    workingset.insert(seed);
    while (workingset.size()) {
        std::set<TPZGeoElSide>::iterator it = workingset.begin();
        TPZGeoElSide active = *it;
        if (active.Neighbour() != active || active.Dimension() != 2) {
            DebugStop();
        }
        workingset.erase(it);
        TPZGeoElBC gelbc(active,matid);
        numelcreated++;
        TPZStack<int> lowsides;
        TPZGeoEl *gelactive = active.Element();
        active.Element()->LowerDimensionSides(active.Side(), lowsides);
        int nsides = lowsides.size();
        // nos vamos percorrerer os vizinhos ao longo dos lados de dimensao um para achar elementos/lados de contorno nao contabilizados
        for (int i=0; i<nsides; i++) {
            int side = lowsides[i];
            if (gelactive->SideDimension(side) != 1) {
                continue;
            }
            TPZGeoElSide hinge(gelactive, side);
            TPZGeoElSide neighbour = hinge;
            while (neighbour != hinge) {
                if (neighbour.Element()->MaterialId() == matid || neighbour.Element()->Dimension() == 2) {
                    break;
                }
                neighbour = neighbour.Neighbour();
            }
            if (neighbour != hinge) {
                continue;
            }
            do {
                TPZStack<TPZGeoElSide> highsides;
                TPZGeoEl *gelloc = neighbour.Element();
                int side = neighbour.Side();
                gelloc->AllHigherDimensionSides(side,2,highsides);
                int numhigh = highsides.size();
                for (int ih=0; ih<numhigh; ih++) {
                    if (highsides[ih].Dimension() != 2) {
                        DebugStop();
                    }
                    if (highsides[ih].Neighbour() == highsides[ih]) {
                        workingset.insert(highsides[ih]);
                        break;
                    }
                }
                neighbour = neighbour.Neighbour();
            } while (neighbour != hinge);
        }
    }
    std::cout << "Number of boundary elements created " << numelcreated << " matid " << matid << std::endl;
}

/// find an element without a 2d neighbour
TPZGeoElSide FindSeed(TPZGeoMesh *gmesh)
{
    long nel = gmesh->NElements();
    for (long el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (!gel || gel->Dimension() !=3) {
            continue;
        }
        int ns = gel->NSides();
        for (int is = 0; is<ns; is++) {
            if (gel->SideDimension(is) != 2) {
                continue;
            }
            TPZGeoElSide gelside(gel, is);
            if (gelside.Neighbour() == gelside ) {
                return gelside;
            }
        }
    }
    return TPZGeoElSide();
}

void BuildSurfaceMesh(TPZGeoMesh *gmesh, int matid)
{
    long nel = gmesh->NElements();
    for (long el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (!gel || gel->Dimension() !=3) {
            continue;
        }
        int ns = gel->NSides();
        for (int is = 0; is<ns; is++) {
            if (gel->SideDimension(is) != 2) {
                continue;
            }
            TPZGeoElSide gelside(gel, is);
            if (gelside.Neighbour() == gelside ) {
                TPZGeoElBC(gelside, matid);
            }
        }
    }
    TPZGeoElSide gelside;
    nel = gmesh->NElements();
    for (long el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(!gel) continue;
        int ns = gel->NSides();
        for(int is=0; is<ns; is++) gel->SetNeighbour(is, gelside);
    }
    TPZAdmChunkVector<TPZGeoEl *> gelvec = gmesh->ElementVec();
    nel = gmesh->NElements();
    for (long el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (gel && gel->Dimension() == 3) {
            gmesh->ElementVec()[el]=0;
        }
    }
    gmesh->BuildConnectivity();
    for (long el=0; el<nel; el++) {
        TPZGeoEl *gel = gelvec[el];
        if (gel && gel->Dimension() == 3) {
            gmesh->ElementVec()[el]=gel;
        }
    }
}