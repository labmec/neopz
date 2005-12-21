#include "pzgeoel.h"
#include "pzgnode.h"

#include "pzintel.h"
#include "pzcompel.h"
#include "pzelcq2d.h"

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
#include "c0-simplequad.cpp"

/*#include <dl_entities.h>
#include <dl_dxf.h>
#include <dl_attributes.h>
*/
#include <vector>

#ifdef LOG4CXX
#include <log4cxx/logger.h>
#include <log4cxx/basicconfigurator.h>
#include <log4cxx/propertyconfigurator.h>
#endif

using namespace std ;

TPZGeoMesh *ReadGeoMesh(ifstream &arq);
void TriangleRefine1(TPZGeoMesh *gmesh);
void TriangleRefine2(TPZGeoMesh *gmesh);
void QuadRefine(TPZGeoMesh *gmesh);
TPZRefPattern *GetBestRefPattern(TPZVec<int> &sides, std::list<TPZRefPattern *> &patlist);

void RefineDirectional(TPZGeoEl *gel,std::set<int> &matids);


//#include "MeshReader.cpp"
void removeComents (std::istream &file, std::string &numberOf);
void ReadNodes (std::istream & file, int nnos, TPZGeoMesh & gMesh);
void ReadElements (std::istream & file, int nelem, TPZGeoMesh & gMesh);
void ReadMaterials (std::istream & file, int nmat, TPZCompMesh & cMesh);
void ReadBCs (std::istream & file, int nmat, TPZCompMesh & cMesh);
TPZCompMesh * ReadMesh(std::istream &file);
void WriteMesh(TPZGeoMesh *mesh,std::ofstream &arq);
void WriteElement (TPZGeoEl *el,int elindex, std::ofstream &arq,TPZVec<int> &elementtype);

void QuadOneRibRefine(TPZGeoMesh *gmesh);
void QuadTwoAdjacentRibRefine(TPZGeoMesh *gmesh);


void LoadRefPattern(std::string &path);

 struct coord{
   public:
     double X;
     double Y;
     double Z;
     coord():X(-9.),Y(-9.),Z(-9.){}
 };

int main ()
{
#ifdef LOG4CXX
    log4cxx::PropertyConfigurator::configure("/home/ic/luis/NeoPZ/Util/config.h");

    //  log4cxx::BasicConfigurator::configure();
  {
    log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("pz.mesh.tpzcompel"));
    logger->setLevel(log4cxx::Level::WARN);
  }
  {
    log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("pz.mesh.tpzcompelside"));
    logger->setLevel(log4cxx::Level::WARN);
  }
  {
    log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("pz.mesh.tpzinterpolatedelement"));
    logger->setLevel(log4cxx::Level::WARN);
  }
  {
    log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("pz.mesh.tpzgeoelrefpattern"));
    logger->setAdditivity(false);
    logger->setLevel(log4cxx::Level::DEBUG);
  }
 {
    log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("pz.mesh.refpattern"));
    logger->setAdditivity(false);
    logger->setLevel(log4cxx::Level::DEBUG);
  }
#endif

// Beggining of the work with the f17.
// Beggining of the work with the F17's mesh.

  // First step: nodes reading and initialization.
  TPZGeoMesh *geomesh = new TPZGeoMesh();
  std::ifstream malha_nos ("yf17.xyz");

  double  x,y,z ;
   
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

  std::ifstream malha_triangulos ("yf17.tri");
  TPZVec<int> indices(3);
  TPZGeoEl *gel;

  int k;
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

  std::ifstream malha_tetraedros ("yf17.tet");
  TPZVec<int> indices2(4);

  int p;
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

  geomesh->BuildConnectivity();

    
  // End of work with the F17's mesh.
  
  // End of the work with f17. 




  //ifstream arq ("Files/Meshes/cruz_hexaedros.txt");
  //ifstream arq ("Files/Meshes/cubo_direcional_3faces.txt");
  //TPZCompMesh *ccmesh = ReadMesh (arq);
  // ccmesh->Reference()->Print(cout);
  // ccmesh->Print(cout);

  //TPZCompMesh *cmesh = CreateSillyMesh();
  //TPZGeoMesh *Gmesh = ccmesh->Reference();

  ofstream dx_arq ("yf17.dx");
  WriteMesh(geomesh,dx_arq);

//   Gmesh->Print(cout);
//   TPZRefPattern * rpt = new TPZRefPattern ("Files/RefPatterns/Triang_Rib_Side_4_5.rpt");
//   rpt->Print();
//   rpt->InsertPermuted(*geomesh);
//   
//   Gmesh->Print(cout);
//   cout << endl;

//   const int nref = 23;
//   TPZRefPattern *refpat [nref];
// 
//   //Line - just uniform refpattern
//   std::cout << "\n\n\n\n\ninsert permuted : line " << std::endl;
//   refpat[0] = new TPZRefPattern ("Files/RefPatterns/Line_Unif.rpt");
// 
// 
//   //Triangle: uniform; one rib and two ribs refinement
//   std::cout << "\n\n\n\n\ninsert permuted : tri unif " << std::endl;
//   refpat[1] = new TPZRefPattern ("Files/RefPatterns/Triang_Unif.rpt");
//   std::cout << "\n\n\n\n\ninsert permuted : tri 1 " << std::endl;
//   refpat[2] = new TPZRefPattern ("Files/RefPatterns/Triang_Rib_Side_3.rpt");
//   std::cout << "\n\n\n\n\ninsert permuted : tri 2 " << std::endl;
//   refpat[3] = new TPZRefPattern ("Files/RefPatterns/Triang_Rib_Side_4_5.rpt");
// 
//   //Quadrilateral: uniform, one rib, two parallel ribs, two adjacent ribs and three ribs
//   std::cout << "\n\n\n\n\ninsert permuted : quad unif " << std::endl;
//   refpat[4] = new TPZRefPattern ("Files/RefPatterns/Quad_Unif.rpt");
//   std::cout << "\n\n\n\n\ninsert permuted : quad 1 " << std::endl;
//   refpat[5] = new TPZRefPattern ("Files/RefPatterns/Quad_Rib_Side_4.rpt");
//   std::cout << "\n\n\n\n\ninsert permuted : quad 2 par " << std::endl;
//   refpat[6] = new TPZRefPattern ("Files/RefPatterns/Quad_Rib_Side_4_6.rpt");
//   std::cout << "\n\n\n\n\ninsert permuted : quad 2 adj " << std::endl;
//   refpat[7] = new TPZRefPattern ("Files/RefPatterns/Quad_Rib_Side_4_5.rpt");
//   std::cout << "\n\n\n\n\ninsert permuted : quad 3 " << std::endl;
//   refpat[8] = new TPZRefPattern ("Files/RefPatterns/Quad_Rib_Side_4_5_6.rpt");
// 
//   //Hexaedros: Um, dois e três lados.
//   refpat[9] = new TPZRefPattern ("Files/RefPatterns/cubo_direcional_3faces_cruz.txt");
//   refpat[10] = new TPZRefPattern ("Files/RefPatterns/cubo_direcional_1.txt");
//   refpat[11] = new TPZRefPattern ("Files/RefPatterns/cubo_direcional_2faces.txt");
//   refpat[12] = new TPZRefPattern ("Files/RefPatterns/cubo_direcional_3faces.txt");
//   refpat[13] = new TPZRefPattern ("Files/RefPatterns/Piram_1lado_base.txt");
//   refpat[14] = new TPZRefPattern ("Files/RefPatterns/Piram_2lados_base.txt");
//   refpat[15] = new TPZRefPattern ("Files/RefPatterns/Piram_2lados_lateral.txt");
//   refpat[16] = new TPZRefPattern ("Files/RefPatterns/Piram_3lados_base.txt");
//   refpat[17] = new TPZRefPattern ("Files/RefPatterns/Piram_3lados_lateral.txt");
//   refpat[18] = new TPZRefPattern ("Files/RefPatterns/Tetra_1lado.txt");
//   refpat[19] = new TPZRefPattern ("Files/RefPatterns/Tetra_2lados.txt");
//   refpat[20] = new TPZRefPattern ("Files/RefPatterns/Tetra_2nodes_001.txt");
//   refpat[21] = new TPZRefPattern ("Files/RefPatterns/Tetra_2nodes_002.txt");
//   refpat[22] = new TPZRefPattern ("Files/RefPatterns/Tetra_2nodes_003.txt");
// 
// 
//   for (int r=0;r<nref;r++)
//   {
//     std::cout << "insert permuted : " << r << std::endl;
//     refpat[r]->InsertPermuted(*geomesh);
//   }
// 


    std::ifstream fonte ("ListedPatterns2.txt");
    geomesh->PatternFileLoad(fonte);

   
  //Gmesh->InsertRefPattern(refpat[10]);
  //Gmesh->InsertRefPattern(refpat[11]);
  //OBS: QUANDO INSIRO O PAT[10] PELO INSERTPERMUTED TENHO ERRO, MAS PELO INSERT
  // REFPATTERN FUNCIONA! 
 
  
  //TPZRefPattern *reftetra = new TPZRefPattern ("/home/pos/cesar/RefPattern/Hexa_Unif.rpt");
  //TPZRefPattern *reftetra = new TPZRefPattern ("/home/pos/cesar/RefPattern/Hexa_Unif.rpt");
  //TPZRefPattern *reftetra = new TPZRefPattern ("/home/ic/luis/Documents/ic/catalogacao/meuRefPattern/Hexaedro_5");
  //TPZRefPattern *reftetra = new TPZRefPattern ("/home/pos/cesar/RefPattern/Piram_Rib_Side_7.rpt");
  //TPZGeoMesh *Gmesh = reftetra->Mesh();
  
  //dxfdraw(Gmesh);  
  //dxfdrawsep(Gmesh);

  //Teste da fun�o de reconhecimento de lados
  //siderecog(Gmesh);
  //Fim do teste

  //TriangleRefine1(Gmesh);
  //TriangleRefine2(Gmesh);
  //QuadRefine(Gmesh);
  //QuadOneRibRefine(Gmesh);
  //QuadTwoAdjacentRibRefine(Gmesh);
  {
    std::ofstream arquivo ("NotListedPatterns.txt");
  }
  {
    std::ofstream padroes ("ListedPatterns.txt");
  }
  {
    std::ofstream padroes_lados ("PatternsandSides.txt");
  }

//   std::ofstream padroes ("ListedPatterns2.txt");
//   geomesh->RefPatternFile(padroes);

  //std::ofstream padroes_lados ("PatternsandSides.txt");
  //geomesh->PatternSidesFile(padroes_lados);

  
  std::set<int> matids;
  matids.insert(-1);
  int count = 0;
  for(i=0; i<10; i++)
  {
    int nelements = geomesh->NElements();
    int el;  
    map<set<int>, TPZRefPattern*> MyMap;
    TPZStack<int> refinesides;
    for (el=0; el<nelements; el++)
    {
      TPZGeoEl *elemento = geomesh->ElementVec()[el];
      RefineDirectional(elemento,matids);
    }
    std::stringstream nome;
    nome << count++ << "_cruz_ref.dx";
    ofstream dx_arq_ref (nome.str().c_str());
    WriteMesh(geomesh,dx_arq_ref);
  }
  ofstream dx_arq_ref ("cruz_ref.dx");
  WriteMesh(geomesh,dx_arq_ref);


  //geomesh->Print(cout);
  //string ref = "refpattern.txt" ;
  /*TPZRefPattern *patt = */ //new TPZRefPattern(ref) ;
  return 0 ;
}




void RefineDirectional(TPZGeoEl *gel,std::set<int> &matids)
{
  int matid = gel->MaterialId();
  if(matids.count(matid)) return;
  TPZVec<int> sidestorefine(gel->NSides(),0);
  TPZVec<int> cornerstorefine(gel->NSides(),0);
  // look for corners which are on the boundary
  int in;
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
      TPZGeoElSide gels(gel,is);
      TPZGeoElSide neigh(gels.Neighbour());
      while(neigh != gels)
      {
        // if this condition is true the rib lies on the boundary
        if(matids.count(neigh.Element()->MaterialId()))
        {
          sidestorefine[is] = 0;
          break;
        }
        neigh = neigh.Neighbour();
      }
    }    
  }
//  TPZGeoMesh *gmesh = gel->Mesh();
  std::list<TPZRefPattern *> patlist;
  TPZRefPattern::GetCompatibleRefinementPatterns(gel, patlist);
  TPZRefPattern *patt = GetBestRefPattern(sidestorefine,patlist);
  static int count = 1;
  if(patt) 
  {
    gel->SetRefPattern(patt);
    TPZManVector<TPZGeoEl *> subel;
    gel->Divide(subel);
    std::cout << "S";
  }
  else
  {
    std::ofstream arquivo ("NotListedPatterns.txt",std::ios::app);
    if(count++ == 1) std::cout << "couldnt find a suitable refinement pattern\n";
    std::cout << "*";
    if(!count%100) std::cout << std::endl;
    arquivo << std::endl;
    arquivo << "Element Type :" << gel->Type() << std::endl;
    arquivo << "Sides selected for refinement :" << std::endl;
    for (int i=0 ; i<gel->NSides() ; i++){
      if(cornerstorefine[i] == 1)
      {
        arquivo << " " << i << " ";
      }
      if (sidestorefine[i] == 1) {
        arquivo << " " << i << " " ;
      }      
    }
    int in;
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
        }
        neigh = neigh.Neighbour();
      }
      arquivo << std::endl;
    }
    arquivo << std::endl;
    // Here we will provide the necessary information to develop a new ref. patt.
  }
  return;
}

TPZRefPattern *GetBestRefPattern(TPZVec<int> &sides, std::list<TPZRefPattern *> &patlist)
{
  std::list<TPZRefPattern *>::iterator it;
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

