//$Id: main.cc,v 1.19 2009-06-24 20:14:56 phil Exp $
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgmesh.h"

#include "pzintel.h"
#include "pzcompel.h"

#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzstack.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"

#include <TPZRefPattern.h>

#include "pzmaterial.h"
#include "pzelasmat.h"
#include "pzlog.h"

#include <time.h>
#include <stdio.h>

#include <math.h>

#include <iostream>
#include <fstream>
#include <string.h>
#include <sstream>

#include <set>
#include <map>
//#include "c0-simplequad.cpp"
//#include "c0-simpleline.cpp"
#include "c0-simpletetra.cpp"

/*#include <dl_entities.h>
#include <dl_dxf.h>
#include <dl_attributes.h>
*/
#include <vector>

// #ifdef LOG4CXX
// #include <log4cxx/logger.h>
// #include <log4cxx/basicconfigurator.h>
// #include <log4cxx/propertyconfigurator.h>
// #endif

using namespace std ;

TPZGeoMesh *ReadGeoMesh(ifstream &arq);
void TriangleRefine1(TPZGeoMesh *gmesh);
void TriangleRefine2(TPZGeoMesh *gmesh);
void QuadRefine(TPZGeoMesh *gmesh);
TPZAutoPointer<TPZRefPattern> GetBestRefPattern(TPZVec<int> &sides, std::list<TPZAutoPointer<TPZRefPattern> > &patlist);

void RefineDirectional(TPZGeoEl *gel,std::set<int> &matids, int destmatid);


//#include "MeshReader.cpp"
void removeComents (std::istream &file, std::string &numberOf);
void ReadNodes (std::istream & file, int nnos, TPZGeoMesh & gMesh);
void ReadElements (std::istream & file, int nelem, TPZGeoMesh & gMesh);
void ReadMaterials (std::istream & file, int nmat, TPZCompMesh & cMesh);
void ReadBCs (std::istream & file, int nmat, TPZCompMesh & cMesh);
TPZCompMesh * ReadMesh(std::istream &file);
void WriteMesh(TPZGeoMesh *mesh,std::ofstream &arq, int materialindex);
void WriteElement (TPZGeoEl *el,int elindex, std::ofstream &arq,TPZVec<int> &elementtype);
void WriteElementMesh(TPZGeoMesh *mesh,std::ofstream &arq, int matindex,int eltype);

void QuadOneRibRefine(TPZGeoMesh *gmesh);
void QuadTwoAdjacentRibRefine(TPZGeoMesh *gmesh);

void InsertNewPattern(std::ifstream &arquivo, TPZGeoMesh &geomesh, std::ofstream &padroes);

void LoadRefPattern(std::string &path);

void  FilterBoundingBox(TPZGeoMesh *geomesh);


 struct coord{
   public:
     double X;
     double Y;
     double Z;
     coord():X(-9.),Y(-9.),Z(-9.){}
 };

 void ReadGidMesh (std::ifstream &fonte, std::ofstream &saida){

   int i;
   double x;
   int d;
   int k;
   for (i=0; i<55; i++){
     fonte >> d;
     saida << d-1 << "  " ;

     for (k=0; k<3; k++){
       fonte >> x ;
       saida << x << "  ";
     }

     saida << std::endl;
   }

   for (i=0; i<107; i++){

     fonte >> d;
     saida << d-1 << "  1  ";

     for (k=0; k<4; k++){
       fonte >> x;
       saida << x << "  ";
     }
     saida << std::endl;
   }
 }


void InsertNewPattern(std::ifstream &arquivo, TPZGeoMesh &geomesh, std::ofstream &padroes)
{
   TPZRefPattern *pattern = new TPZRefPattern (&geomesh, arquivo);
   pattern->InsertPermuted(/*geomesh*/);
   delete pattern;
   geomesh.RefPatternFile(padroes);
 }


void ReadF17(TPZGeoMesh *geomesh)
{
// // Beggining of the work with the f17.
// // Beggining of the work with the F17's mesh.
//
//   // First step: nodes reading and initialization.
  std::string path;
#ifdef HAVE_CONFIG_H
  path = PZSOURCEDIR;
  path += "/Projects/Heman/Files/Meshes/";
#else
  path = "";
#endif

  std::cout << "===============================================================\n"
            << "Reading F17 mesh\n";


  std::string nodFile = path;
  nodFile += "yf17.xyz";
  std::cout << "\t\tInput file for nodes = " << nodFile.c_str()
            << "\n\t\t\tprocessing nodes...\n";
  std::ifstream malha_nos (nodFile.c_str());

//  double  x,y,z ;

  TPZVec< REAL > nodes(3);
  int i;
  while(malha_nos) {
    malha_nos >> nodes[0];
    malha_nos >> nodes[1];
    malha_nos >> nodes[2];
    if(!malha_nos) continue;
    int nodind = geomesh->NodeVec().AllocateNewElement();

    geomesh->NodeVec()[nodind] = TPZGeoNode(i,nodes,*geomesh);
  }

  // Second step: Generation of the triangles.

  std::string triFile = path;
  triFile += "yf17.tri";
  std::cout << "\t\tInput file for faces = " << triFile.c_str()
            << "\n\t\t\tprocessing faces...\n";
  std::ifstream malha_triangulos (triFile.c_str());
  TPZVec<int> indices(3);
  TPZGeoEl *gel;

//  int k;
  while (malha_triangulos){
    malha_triangulos >> indices[0] ;
    malha_triangulos >> indices[1] ;
    malha_triangulos >> indices[2] ;
    int index;
    indices[0]--;
    indices[1]--;
    indices[2]--;
    if(!malha_triangulos) continue;
    gel = geomesh->CreateGeoElement(ETriangle,indices,-1,index,1);

  }

  // Third step: Generation of the thetra.
  std::string tetraFile = path;
  tetraFile += "yf17.tet";
  std::cout << "\t\tInput file for volumes = " << tetraFile.c_str()
            << "\n\t\t\tprocessing volumes...\n";
  std::ifstream malha_tetraedros (tetraFile.c_str());
  TPZVec<int> indices2(4);

//  int p;
  while (malha_tetraedros){
    malha_tetraedros >> indices2[0] ;
    malha_tetraedros >> indices2[1] ;
    malha_tetraedros >> indices2[2] ;
    malha_tetraedros >> indices2[3] ;
    indices2[0]--;
    indices2[1]--;
    indices2[2]--;
    indices2[3]--;
    if(!malha_tetraedros) continue;
    int index;
    gel = geomesh->CreateGeoElement(ETetraedro,indices2,1,index,1);
  }

  // Fourth step: Building the mesh.

  std::cout << "\n\tGenerating the connectivities and neighborhood information...\n";
  geomesh->BuildConnectivity();
  std::cout << "F17 mesh read\n"
            << "===============================================================\n";
//   // End of work with the F17's mesh.
//
//   // End of the work with f17.
}

void LoadRefPatternDataBase(TPZGeoMesh *geomesh)
{
  std::string path;
#ifdef HAVE_CONFIG_H
  path = PZSOURCEDIR;
#else
  path = "";
#endif

  std::string prefix = path;
  prefix += "/Refine/RefPatterns/";
  std::string allfiles = prefix + "Patternlist";

  std::string allpatterns = prefix + "ListedPatterns10.txt";
  std::cout << "Reading refinement patterns: " << allpatterns.c_str() <<  std::endl;
  std::ifstream fonte (allpatterns.c_str(),ios::in);
  if(!fonte.fail())
  {
    try {
      geomesh->PatternFileLoad(fonte);
    }
    catch (...)
    {
      std::cout << __PRETTY_FUNCTION__ << std::endl;
      exit (-1);
    }
  }
  else
  {
    std::ifstream filelist(allfiles.c_str(),ios::in);
    while(!filelist.fail())
    {
      std::string filename;
      getline(filelist,filename);
      std::string fullname = prefix+filename;
      if(!filelist.fail())
      {
        TPZRefPattern *refp = new TPZRefPattern(geomesh, fullname);
        refp->InsertPermuted(/**geomesh*/);
        delete refp;
      }
    }
    std::ofstream out(allpatterns.c_str());
    geomesh->RefPatternFile(out);

  }

  std::ofstream shortlist("shortlist.txt");
  geomesh->PatternSidesFile(shortlist);

  {
    std::ofstream arquivo ("NotListedPatterns4.txt");
  }
  {
    std::ofstream padroes ("ListedPatterns.txt");
  }
  {
    std::ofstream padroes_lados ("PatternsandSides3.txt");
  }
}


void InitializeLOG()
{
  std::string path;
  std::string configfile;
#ifdef HAVE_CONFIG_H
  path = PZSOURCEDIR;
#else
  path = "";
#endif
#ifdef LOG4CXX
  configfile = path;
  configfile += "/Util/log4cxx.cfg";
  log4cxx::PropertyConfigurator::configure(configfile);

    //  log4cxx::BasicConfigurator::configure();
  {
    log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("pz.mesh.tpzcompel"));
 //   logger->setLevel(log4cxx::Level::WARNING);
  }
  {
    log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("pz.mesh.tpzcompelside"));
 //   logger->setLevel(log4cxx::Level::WARNING);
  }
  {
    log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("pz.mesh.tpzinterpolatedelement"));
 //   logger->setLevel(log4cxx::Level::WARNING);
  }
  {
    log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("pz.mesh.tpzgeoelrefpattern"));
    logger->setAdditivity(false);
//    logger->setLevel(log4cxx::Level::DEBUG);
  }
 {
    log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("pz.mesh.refpattern"));
    logger->setAdditivity(false);
//    logger->setLevel(log4cxx::Level::DEBUG);
  }
#endif
}


TPZGeoMesh * choiceMesh(std::string &meshName)
{
  TPZGeoMesh *geomesh;
  int opt = -1;
  cout << "Choose the mesh to be processed:\n"
      << "\t1 - plane plate\n"
      << "\t2 - yf17\n";
  cin >> opt;

  switch (opt)
  {
    case (1) :
    {
      std::string file_path = PZSOURCEDIR;
      file_path += "/Projects/Heman/Files/Meshes/placa_plana.txt";
      std::ifstream file (file_path.c_str());
      TPZCompMesh *compmesh = ReadMesh (file);
      geomesh = compmesh->Reference();
      meshName = "placa_plana";
      break;
    }
    case (2) :
    {
      geomesh = new TPZGeoMesh;
      ReadF17(geomesh);
      meshName = "yf17";
      break;
    }
    default:
      cout << "\nOption: " << opt << " is not a valid choice!\n"
          << "Exiting without process any mesh!\n";
      exit (-1);
      break;
  }
  return geomesh;
}

void LoadPatternDB(TPZGeoMesh *gMesh)
{
  cout << "\tinitialize uniform refinement patterns...\n";
  gMesh->InitializeRefPatterns();

  cout << "\tinitialize refinement patterns from file data base...\n";
  LoadRefPatternDataBase(gMesh);
/*  cout << "...0...\n";
  std::string file_path = PZSOURCEDIR;
  file_path += "/Projects/Heman/Files/RefPatterns/Hexa_Half.rpt";
  std::ifstream inparq (file_path.c_str());
  std::string outfile = PZSOURCEDIR;
  outfile += "/Projects/Heman/Files/RefPatterns/hexa_plus_tetra.rpt";
  std::ofstream outarq (outfile.c_str());

  cout << "...1...\n";
  InsertNewPattern(inparq, *gMesh, outarq);

  file_path = PZSOURCEDIR;
  file_path += "/Projects/Heman/Files/RefPatterns/Hexa_2AdjRibs.rpt";
  std::ifstream inparq2 (file_path.c_str());
  cout << "...2...\n";
  InsertNewPattern(inparq2, *gMesh, outarq);

  file_path = PZSOURCEDIR;
  file_path += "/Projects/Heman/Files/RefPatterns/Prism_2AdjRibs.rpt";
  std::ifstream inparq3 (file_path.c_str());
  cout << "...3...\n";
  InsertNewPattern(inparq3, *gMesh, outarq);*/
}

int main ()
{
  cout << "Initilizing log system...\n";
  InitializePZLOG();

  std::string meshname;
  TPZGeoMesh *geomesh = choiceMesh(meshname);

  cout << "Loading refinement patterns...\n";
  LoadPatternDB(geomesh);

  cout << "Filtering bounding box...\n";
  FilterBoundingBox(geomesh);

  ofstream dx_arq ("surface.dx");
  WriteMesh(geomesh,dx_arq,1);

  std::set<int> matids;
  matids.insert(-1);
  int count = 0;
  int destmatid = 2;

  int el, nref = -1;
  cout << "number of refinement steps:\n";
  cin >> nref;
  cout << endl;

  int i;
  for(i=0; i<nref; i++)
  {
    int nelements = geomesh->NElements();
    cout << "\tEntering refinement step: " << i << " number of elements = " << nelements << endl;
    cout << "\tRefining...\n";

    //map<set<int>, TPZRefPattern*> MyMap;
    TPZStack<int> refinesides;
    for (el=0; el<nelements; el++)
    {
      cout << ".";
      if (!(el+1)%40) cout << endl;
      else cout.flush();
      TPZGeoEl *elemento = geomesh->ElementVec()[el];
      RefineDirectional(elemento,matids,destmatid);
    }
    cout << "\nWriting dx files...\n";
    {
      std::stringstream nome(meshname);
      nome << count << "_boundary_a.dx";
      ofstream dx_arq_ref (nome.str().c_str());
      WriteMesh(geomesh,dx_arq_ref,destmatid);
    }
    {
      cout << "\tfor all boundaries\n";
      std::stringstream nome(meshname);
      nome << count++ << "_boundary_ref";
      std::string elMesh = nome.str().c_str();
      nome << ".dx";
      ofstream dx_arq_ref (nome.str().c_str());
      WriteMesh(geomesh,dx_arq_ref,destmatid+1);

      cout << "\tfor tetrahedres\n";
      std::string aux = elMesh;
      aux += "-tetra.dx";
      std::ofstream dx_tetra_arq (aux.c_str());
      WriteElementMesh(geomesh,dx_tetra_arq,destmatid+1,4);

      cout << "\tfor pyramides\n";
      aux = elMesh;
      aux += "-pyramid.dx";
      std::ofstream dx_pyra_arq (aux.c_str());
      WriteElementMesh(geomesh,dx_pyra_arq,destmatid+1,5);

      cout << "\tfor prims\n";
      aux = elMesh;
      aux += "-prism.dx";
      std::ofstream dx_prism_arq (aux.c_str());
      WriteElementMesh(geomesh,dx_prism_arq,destmatid+1,6);

      cout << "\tfor hexahedres\n";
      aux = elMesh;
      aux += "-hexa.dx";
      std::ofstream dx_hexa_arq (aux.c_str());
      WriteElementMesh(geomesh,dx_hexa_arq,destmatid+1,7);
    }
    destmatid += 2;
  }
  cout << "Finalizing...\n";
  destmatid--;
  //geomesh->Print(cout);
  //string ref = "refpattern.txt" ;
  /*TPZRefPattern *patt = */ //new TPZRefPattern(ref) ;
  //geomesh->Print();
  return 0 ;
}


void RefineDirectional(TPZGeoEl *gel,std::set<int> &matids,int destmatid)
{
  int matid = gel->MaterialId();
  if(matids.count(matid)) return;
  TPZManVector<int,27> sidestorefine(gel->NSides(),0);
  TPZManVector<int,27> cornerstorefine(gel->NSides(),0);
  // look for corners which are on the boundary
  int in;
  int numrefribs = 0;
  for(in=0; in<gel->NCornerNodes(); in++)
  {
    TPZGeoElSide gels(gel,in);
    TPZGeoElSide neigh(gels.Neighbour());
    while(gels != neigh)
    {
      if(matids.count(neigh.Element()->MaterialId()))
      {
        cornerstorefine[in] = 1;
        break;
      }
      neigh = neigh.Neighbour();
    }
  }
  // look for ribs which touch the boundary but which do no lay on the boundary
  int is;
  for(is=gel->NCornerNodes(); is<gel->NSides(); is++)
  {
    // we are only interested in ribs
    if(gel->SideDimension(is) != 1) continue;

    // the side is a candidate if it contains a corner which is neighbour of the boundary condition
    if(cornerstorefine[gel->SideNodeLocIndex(is,0)] || cornerstorefine[gel->SideNodeLocIndex(is,1)])
    {
      sidestorefine[is] = 1;
      numrefribs++;
      TPZGeoElSide gels(gel,is);
      TPZGeoElSide neigh(gels.Neighbour());
      while(neigh != gels)
      {
        // if this condition is true the rib lies on the boundary
        if(matids.count(neigh.Element()->MaterialId()))
        {
          sidestorefine[is] = 0;
          numrefribs--;
          break;
        }
        neigh = neigh.Neighbour();
      }
    }
  }
  if(!numrefribs)
  {
    return;
  }
//  TPZGeoMesh *gmesh = gel->Mesh();
  std::list<TPZAutoPointer<TPZRefPattern> > patlist;
  TPZRefPattern::GetCompatibleRefinementPatterns(gel, patlist);
  TPZAutoPointer<TPZRefPattern> patt = GetBestRefPattern(sidestorefine,patlist);
  static int count = 1;
  if(patt)
  {
    gel->SetMaterialId(destmatid+1);
    gel->SetRefPattern(patt);
    TPZManVector<TPZGeoEl *> subel;
    gel->Divide(subel);
    gel->SetMaterialId(destmatid);
    std::cout << "-";
  }
  else
  {
    std::ofstream arquivo ("NotListedPatterns4.txt",std::ios::app);
    if(count++ == 1) std::cout << "couldnt find a suitable refinement pattern\n";
    std::cout << "|";
    std::list<TPZAutoPointer<TPZRefPattern> >::iterator it;
    arquivo << "Compatible refinement patterns\n";
    for(it=patlist.begin(); it!=patlist.end(); it++)
    {
      (*it)->ShortPrint(arquivo); arquivo << (void*) (it->operator->()); arquivo << endl;
    }
    arquivo << std::endl;
    arquivo << "Element Type :" << gel->Type() << std::endl;
    arquivo << "Sides selected for refinement :" << std::endl;
    int i;
    for (i=0 ; i<gel->NSides() ; i++){
      if(cornerstorefine[i] == 1)
      {
        arquivo << " " << i << " ";
      }
      if (sidestorefine[i] == 1) {
        arquivo << " " << i << " " ;
      }
    }
    gel->Print(arquivo);
    int in;
    arquivo << std::endl;
    arquivo << "Neighbouring information \n";
    for(in=0; in<gel->NSides(); in++)
    {
      arquivo << "Side : " << in << " ";
      TPZGeoElSide gels(gel,in);
      arquivo << "Dim " << gels.Dimension() << " ";
      TPZGeoElSide neigh(gels.Neighbour());
      while(gels != neigh)
      {
        if(matids.count(neigh.Element()->MaterialId()))
        {
          arquivo << neigh.Element()->Id() << "-l-" << neigh.Side() << " ";
          if (neigh.Side() == 9 && gel->Type() == ETetraedro) {
            arquivo << "Teje pego meliante..." << std::endl;
            neigh.Element()->Print(arquivo);
          }
        }
        neigh = neigh.Neighbour();
      }
      arquivo << std::endl;
    }
    arquivo << std::endl;

    arquivo << "Element information : " << gel->Index() << std::endl;
    arquivo << "Vizinhos dos lados marcados para refinamento:" << std::endl;
    for (i=0 ; i<gel->NSides() ; i++){
      if(cornerstorefine[i] == 1 || sidestorefine[i] == 1)
      {
        TPZGeoElSide gelside (gel,i);
        TPZGeoElSide neigh = gelside.Neighbour();
        while (neigh != gelside)
        {
          arquivo << "*********** my side = " << i << " neighside " << neigh.Side() << std::endl;
          neigh.Element()->Print(arquivo);
          neigh = neigh.Neighbour();
        }
      }
    }
    arquivo << std::endl << std::endl << std::endl << std::endl;
    // Here we will provide the necessary information to develop a new ref. patt.
  }
  count++;
  if(!(count%20))
  {
    std::cout << count << std::endl;
  }

  return;
}

TPZAutoPointer<TPZRefPattern> GetBestRefPattern(TPZVec<int> &sides, std::list<TPZAutoPointer<TPZRefPattern> > &patlist)
{
  std::list<TPZAutoPointer<TPZRefPattern> >::iterator it;
  for(it = patlist.begin(); it != patlist.end(); it++)
  {
    if(! (*it)) continue;
    TPZGeoEl *gel = (*it)->Element(0);
    int ncorners = gel->NCornerNodes();
    int nsides = gel->NSides();
    int is;
    for(is = ncorners; is<nsides; is++)
    {
      if(sides[is] != (*it)->NSideNodes(is)) break;
    }
    if(is == nsides) return (*it);
  }
  return 0;
}

void FilterBoundingBox(TPZGeoMesh *geomesh)
{
  int no,nnode,idf;
  nnode = geomesh->NodeVec().NElements();
  TPZFNMatrix<6> xminmax(3,2,0.);
  for(no=0; no<nnode; no++)
  {
    for(idf=0; idf<3; idf++)
    {
      if(geomesh->NodeVec()[no].Coord(idf) < xminmax(idf,0)) xminmax(idf,0) = geomesh->NodeVec()[no].Coord(idf);
      if(geomesh->NodeVec()[no].Coord(idf) > xminmax(idf,1)) xminmax(idf,1) = geomesh->NodeVec()[no].Coord(idf);
    }
  }
  int el,nelem;
  nelem = geomesh->ElementVec().NElements();
  for(el=0; el<nelem; el++)
  {
      TPZGeoEl *gel = geomesh->ElementVec()[el];
      if(!geomesh->ElementVec()[el] || ! (geomesh->ElementVec()[el]->Type()== ETriangle) || gel->MaterialId() != -1) continue;
      TPZFNMatrix<9> axes(3,3),jac(2,2),jacinv(2,2),xminmaxloc(3,2,0.);
      for(no=0; no<gel->NNodes(); no++)
      {
        TPZGeoNode *gno = gel->NodePtr(no);
        for(idf=0; idf<3; idf++)
        {
          if(no == 0)
          {
            xminmaxloc(idf,0) = gno->Coord(idf);
            xminmaxloc(idf,1) = gno->Coord(idf);
          }
          if(gno->Coord(idf) < xminmaxloc(idf,0)) xminmaxloc(idf,0) = gno->Coord(idf);
          if(gno->Coord(idf) > xminmaxloc(idf,1)) xminmaxloc(idf,1) = gno->Coord(idf);
        }
      }
      xminmaxloc -= xminmax;
      REAL mindif = 1.;
      for(no=0; no<2; no++) for(idf=0; idf<3; idf++) mindif = (mindif < fabs(xminmaxloc(idf,no)))? mindif : fabs(xminmaxloc(idf,no));
      if(mindif < 0.001)
      {
        REAL detjac;
        TPZManVector<REAL,3> coor(2,0.3333);
        gel->Jacobian(coor,jac,axes,detjac,jacinv);
        if((fabs(axes(2,0))-1.) < 1.e-3 || (fabs(axes(2,1))-1.) < 1.e-3 || (fabs(axes(2,2))-1.) < 1.e-3)
        {
          gel->SetMaterialId(-4);
        }
        else
        {
          gel->SetMaterialId(-1);
        }
      }
  }
}
