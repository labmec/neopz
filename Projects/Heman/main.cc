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


  ifstream arq ("/Files/Meshes/cruz_hexaedros.txt");
  TPZCompMesh *ccmesh = ReadMesh (arq);
  // ccmesh->Reference()->Print(cout);
  // ccmesh->Print(cout);

  //TPZCompMesh *cmesh = CreateSillyMesh();
  TPZGeoMesh *Gmesh = ccmesh->Reference();

  ofstream dx_arq ("cruz.dx");
  WriteMesh(Gmesh,dx_arq);

  TPZRefPattern * rpt = new TPZRefPattern ("Files/RefPatterns/Triang_Rib_Side_4_5.rpt");
   rpt->Print();
   rpt->InsertPermuted(*Gmesh);
  
  //Gmesh->Print(cout);
  //cout << endl;

  const int nref = 11;
  TPZRefPattern *refpat [nref];

  //Line - just uniform refpattern
  std::cout << "\n\n\n\n\ninsert permuted : line " << std::endl;
  refpat[0] = new TPZRefPattern ("Files/RefPatterns/Line_Unif.rpt");


  //Triangle: uniform; one rib and two ribs refinement
  std::cout << "\n\n\n\n\ninsert permuted : tri unif " << std::endl;
  refpat[1] = new TPZRefPattern ("Files/RefPatterns/Triang_Unif.rpt");
  std::cout << "\n\n\n\n\ninsert permuted : tri 1 " << std::endl;
  refpat[2] = new TPZRefPattern ("Files/RefPatterns/Triang_Rib_Side_3.rpt");
  std::cout << "\n\n\n\n\ninsert permuted : tri 2 " << std::endl;
  refpat[3] = new TPZRefPattern ("Files/RefPatterns/Triang_Rib_Side_4_5.rpt");

  //Quadrilateral: uniform, one rib, two parallel ribs, two adjacent ribs and three ribs
  std::cout << "\n\n\n\n\ninsert permuted : quad unif " << std::endl;
  refpat[4] = new TPZRefPattern ("Files/RefPatterns/Quad_Unif.rpt");
  std::cout << "\n\n\n\n\ninsert permuted : quad 1 " << std::endl;
  refpat[5] = new TPZRefPattern ("Files/RefPatterns/Quad_Rib_Side_4.rpt");
  std::cout << "\n\n\n\n\ninsert permuted : quad 2 par " << std::endl;
  refpat[6] = new TPZRefPattern ("Files/RefPatterns/Quad_Rib_Side_4_6.rpt");
  std::cout << "\n\n\n\n\ninsert permuted : quad 2 adj " << std::endl;
  refpat[7] = new TPZRefPattern ("Files/RefPatterns/Quad_Rib_Side_4_5.rpt");
  std::cout << "\n\n\n\n\ninsert permuted : quad 3 " << std::endl;
  refpat[8] = new TPZRefPattern ("Files/RefPatterns/Quad_Rib_Side_4_5_6.rpt");
  
  //Hexaedros: Um, dois e três lados.
  refpat[9] = new TPZRefPattern ("Files/RefPatterns/cubo_direcional_1.txt");
  refpat[10] = new TPZRefPattern ("Files/RefPatterns/cubo_direcional_2faces.txt");
  //refpat[11] = new TPZRefPattern ("Files/RefPatterns/cubo_direcional_3faces.txt");
  
  for (int r=0;r<nref;r++)
  {
    std::cout << "insert permuted : " << r << std::endl;
    refpat[r]->InsertPermuted(*Gmesh);
  }
  
  //Gmesh->InsertRefPattern(refpat[10]);
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

  std::set<int> matids;
  matids.insert(-1);
  int i,count = 0;
  for(i=0; i<10; i++)
  {
    int nelements = Gmesh->NElements();
    int el;  
    map<set<int>, TPZRefPattern*> MyMap;
    TPZStack<int> refinesides;
    for (el=0; el<nelements; el++)
    {
      TPZGeoEl *elemento = Gmesh->ElementVec()[el];
      RefineDirectional(elemento,matids);
    }
    std::stringstream nome;
    nome << count++ << "_cruz_ref.dx";
    ofstream dx_arq_ref (nome.str().c_str());
    WriteMesh(Gmesh,dx_arq_ref);
  }
  ofstream dx_arq_ref ("cruz_ref.dx");
  WriteMesh(Gmesh,dx_arq_ref);
  //string ref = "refpattern.txt" ;
  /*TPZRefPattern *patt = */ //new TPZRefPattern(ref) ;
  return 0 ;
}




void RefineDirectional(TPZGeoEl *gel,std::set<int> &matids)
{
  int matid = gel->MaterialId();
  if(matids.count(matid)) return;
  TPZVec<int> sidestorefine(gel->NSides(),0);
  TPZVec<int> cornerstorefine(gel->NCornerNodes(),0);
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
  if(patt) 
  {
    gel->SetRefPattern(patt);
    TPZManVector<TPZGeoEl *> subel;
    gel->Divide(subel);
  }
  else
  {
    std::cout << "couldnt find a suitable refinement pattern\n";
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

