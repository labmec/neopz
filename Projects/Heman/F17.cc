//$Id: main.cc,v 1.21 2010-07-22 17:43:43 caju Exp $
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "pzbfilestream.h"
#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgmesh.h"
#include "pzbiharmonic.h"
#include "pzcmesh.h"
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
#include "pzgeopyramid.h"
#include "TPZGeoLinear.h"

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
#include "TPZRefPatternDataBase.h"
#include "TPZRefPatternTools.h"
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

void PrintGMeshVTKcopy(TPZGeoMesh * gmesh, std::ofstream &file);
void PrintGMeshVTKneighbourcopy(TPZGeoMesh * gmesh, std::ofstream &file, int neighMaterial);
void PrintGMeshVTKmaterialcopy(TPZGeoMesh * gmesh, std::ofstream &file, std::set<int> Material);

void RefineDirectional(TPZGeoEl *gel,std::set<int> &matids, int destmatid);

TPZCompMesh *CreateCompMesh ( TPZGeoMesh &gmesh, int porder );


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

void InsertNewPattern(std::ifstream &arquivo, TPZGeoMesh &geomesh, std::ofstream &padroes)
{
   TPZRefPattern *pattern = new TPZRefPattern(arquivo);
   pattern->InsertPermuted();
	
   delete pattern;
	
	//AQUI
   gRefDBase.WriteRefPatternDBase(padroes);
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
  path = "/Users/Cesar/Documents/Projects/NeoPZ/Projects/Heman/Files/Meshes/";
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
      //AQUI
      gRefDBase.ReadRefPatternDBase(fonte);
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
        TPZRefPattern *refp = new TPZRefPattern(fullname);
        refp->InsertPermuted();
        delete refp;
      }
    }
    std::ofstream out(allpatterns.c_str());
	  
	//AQUI
    gRefDBase.WriteRefPatternDBase(out);
  }

		//AQUI
	/*
  std::ofstream shortlist("shortlist.txt");
  geomesh->PatternSidesFile(shortlist);
*/
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
	int opt = 2;//-1;
  cout << "Choose the mesh to be processed:\n"
      << "\t1 - plane plate\n"
      << "\t2 - yf17\n";
  //cin >> opt;

  switch (opt)
  {
    case (1) :
    {
			//AQUI
		std::string file_path;// = PZSOURCEDIR;
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

int main/*F17*/()
{
  cout << "Initilizing log system...\n";
  InitializePZLOG("log4cxx.cfg");
	
  gRefDBase.InitializeRefPatterns();
	
  std::string meshname;
  TPZGeoMesh *geomesh = choiceMesh(meshname);

  cout << "Filtering bounding box...\n";
  FilterBoundingBox(geomesh);

  ofstream dx_arq ("surface.dx");
  WriteMesh(geomesh,dx_arq,1);

  std::set<int> matids;
  matids.insert(-1);
  int count = 0;
  int destmatid = 2;

  int el, nref = 3;//-1;
  cout << "number of refinement steps:\n";
  //cin >> nref;
  cout << endl;

	ofstream outt("F17casca0ref_FINAL.vtk");
	std::set<int> mat1;
	mat1.insert(-1);
	
	TPZRefPatternTools::PrintGMeshVTKneighbour_material(geomesh, outt, -1);	
	
  int i;
  for(i=0; i<nref; i++)
  {
    int nelements = geomesh->NElements();
    cout << "\tEntering refinement step: " << i << " number of elements = " << nelements << endl;
    cout << "\tRefining...\n";

    for (el=0; el<nelements; el++)
    {
      TPZGeoEl *elemento = geomesh->ElementVec()[el];
	  TPZRefPatternTools::RefineDirectional(elemento, matids, destmatid);
    }
	destmatid += 2;
	  
	  std::stringstream sout;
	  int ii = i + 1;
	  sout << "F17casca" << ii << "ref_FINAL.vtk";
	  
	  ofstream out(sout.str().c_str());
	  std::set<int> mat2;
	  mat2.insert(-1);
	  
	  TPZRefPatternTools::PrintGMeshVTKneighbour_material(geomesh, out, -1);	  
  }
  cout << "Finalizing...\n";
  destmatid--;
	
	
  return 0 ;
}

void RefineDirectional(TPZGeoEl *gel, std::set<int> &matids, int destmatid)
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
	
	std::list<TPZAutoPointer<TPZRefPattern> > compatible;
		//TPZAutoPointer<TPZRefPattern> patt = TPZRefPatternTools::PerfectMatchRefPattern(gel);
	TPZRefPatternTools::GetCompatibleRefPatterns(gel,patlist);
	TPZAutoPointer<TPZRefPattern> patt = GetBestRefPattern(sidestorefine, patlist);
	
  static int count = 1;
  if(patt)
  {
    gel->SetMaterialId(destmatid+1);
    gel->SetRefPattern(patt);
    TPZManVector<TPZGeoEl *> subel;
    gel->Divide(subel);
    gel->SetMaterialId(destmatid);
  }
  else
  {
	cout << "Nao encontrei o padrao!!!" << std::endl;
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
	int num4 = 0;
	int num1 = 0;
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
		  
		  TPZManVector<REAL,3> vectorialProd(3,0.);
		  vectorialProd[0] = -axes(0,2)*axes(1,1) + axes(0,1)*axes(1,2);
		  vectorialProd[1] = axes(0,2)*axes(1,0) - axes(0,0)*axes(1,2);
		  vectorialProd[2] = -axes(0,1)*axes(1,0) + axes(0,0)*axes(1,1);
        if( fabs(fabs(vectorialProd[0]) -1.) < 1.e-3 || fabs (fabs(vectorialProd[1])-1.) < 1.e-3 || fabs (fabs(vectorialProd[2])-1.) < 1.e-3 )
        {
          gel->SetMaterialId(-4);
			num4++;
        }
        else
        {
          gel->SetMaterialId(-1);
			num1++;
        }
      }
  }
	std::cout << "Number of elements with -4 condition " << num4 << " with -1 condition " << num1 << std::endl;
}

TPZCompMesh *CreateCompMesh( TPZGeoMesh &gmesh, int porder )
{
	TPZCompEl::SetgOrder ( porder );
	TPZCompMesh *result = new TPZCompMesh ( &gmesh );
	result->SetDimModel ( 2 );
	result->SetAllCreateFunctionsDiscontinuousReferred();
		//result->SetAllCreateFunctionsContinuousReferred();
	TPZBiharmonic *material ;
	material = new TPZBiharmonic( 1,0. );
	TPZAutoPointer<TPZMaterial> material9 = new TPZBiharmonic( 9,0. );
	TPZAutoPointer<TPZMaterial> mat ( material );
	
	result->InsertMaterialObject ( mat );
	result->InsertMaterialObject(material9);
	
	TPZFMatrix val1 ( 2,1,0. ), val2 ( 2,1,0. ); // (2,1,0.) pq cada cond. tem Dirichlet e Neumann
	TPZAutoPointer<TPZMaterial> bnd1 = material->CreateBC ( mat,-1,0, val1, val2 );
	TPZAutoPointer<TPZMaterial> bnd2 = material->CreateBC ( mat,-2,0, val1, val2 );
	TPZAutoPointer<TPZMaterial> bnd3 = material->CreateBC ( mat,-3,0, val1, val2 );
	TPZAutoPointer<TPZMaterial> bnd4 = material->CreateBC ( mat,-4,0, val1, val2 );
	
	result->InsertMaterialObject ( bnd1 );
	result->InsertMaterialObject ( bnd2 );
	result->InsertMaterialObject ( bnd3 );
	result->InsertMaterialObject ( bnd4 );
	
	result->AutoBuild();
	ofstream arc ( "cmesh.txt" );
	result->Print ( arc );
	return result;
}

void PrintGMeshVTKcopy(TPZGeoMesh * gmesh, std::ofstream &file)
{
	file.clear();
	int nelements = gmesh->NElements();
	
	stringstream node, connectivity, type;
	
	//Header
	file << "# vtk DataFile Version 3.0" << std::endl;
	file << "TPZGeoMesh VTK Visualization" << std::endl;
	file << "ASCII" << std::endl << std::endl;
	
	file << "DATASET UNSTRUCTURED_GRID" << std::endl;
	file << "POINTS ";
	
	int actualNode = -1, size = 0, nVALIDelements = 0;
	
	for(int el = 0; el < nelements; el++)
	{				
		if(gmesh->ElementVec()[el]->Type() == EOned)//Exclude Lines and Arc3D
		{
			continue;
		}
		if(gmesh->ElementVec()[el]->HasSubElement())
		{
			continue;
		}
		
		int elNnodes = gmesh->ElementVec()[el]->NNodes();
		size += (1+elNnodes);
		connectivity << elNnodes;
		
		for(int t = 0; t < elNnodes; t++)
		{
			for(int c = 0; c < 3; c++)
			{
				double coord = gmesh->NodeVec()[gmesh->ElementVec()[el]->NodeIndex(t)].Coord(c);
				node << coord << " ";
			}			
			node << std::endl;
			
			actualNode++;
			connectivity << " " << actualNode;
		}
		connectivity << std::endl;
		
		int elType;
		switch (gmesh->ElementVec()[el]->Type())
		{
			case (ETriangle):
			{
				elType = 5;
				break;				
			}
			case (EQuadrilateral):
			{
				elType = 9;
				break;				
			}
			case (ETetraedro):
			{
				elType = 10;
				break;				
			}
			case (EPiramide):
			{
				elType = 14;
				break;				
			}
			case (EPrisma):
			{
				elType = 13;
				break;				
			}
			case (ECube):
			{
				elType = 12;
				break;				
			}
			default:
			{
				elType = -1;//ElementType NOT Found!!!
				DebugStop();
				break;	
			}
		}
		type << elType << std::endl;
		nVALIDelements++;
	}
	node << std::endl;
	actualNode++;
	file << actualNode << " float" << std::endl << node.str();
	
	file << "CELLS " << nVALIDelements << " ";
	
	file << size << std::endl;
	file << connectivity.str() << std::endl;
	
	file << "CELL_TYPES " << nVALIDelements << std::endl;
	file << type.str();
}

void PrintGMeshVTKneighbourcopy(TPZGeoMesh * gmesh, std::ofstream &file, int neighMaterial)
{
	file.clear();
	int nelements = gmesh->NElements();
	
	stringstream node, connectivity, type;
	
	//Header
	file << "# vtk DataFile Version 3.0" << std::endl;
	file << "TPZGeoMesh VTK Visualization" << std::endl;
	file << "ASCII" << std::endl << std::endl;
	
	file << "DATASET UNSTRUCTURED_GRID" << std::endl;
	file << "POINTS ";
	
	int actualNode = -1, size = 0, nVALIDelements = 0;
	
	for(int el = 0; el < nelements; el++)
	{				
		if(gmesh->ElementVec()[el]->Type() == EOned)//Exclude Lines and Arc3D
		{
			continue;
		}
		if(gmesh->ElementVec()[el]->HasSubElement())
		{
			continue;
		}
		
		if(neighMaterial != -1E-15)
		{
			bool matFound = false;
			for(int s = 0; s < gmesh->ElementVec()[el]->NSides(); s++)
			{
				TPZGeoElSide thisSide(gmesh->ElementVec()[el], s);
				TPZGeoElSide neighSide = thisSide.Neighbour();
				
				while(thisSide != neighSide)
				{
					if(neighSide.Element()->MaterialId() == neighMaterial)
					{
						matFound = true;
						break;
					}
					neighSide = neighSide.Neighbour();
				}
				if(matFound)
				{
					break;
				}
			}
			if(!matFound)
			{
				continue;
			}
		}
		
		int elNnodes = gmesh->ElementVec()[el]->NNodes();
		size += (1+elNnodes);
		connectivity << elNnodes;
		
		for(int t = 0; t < elNnodes; t++)
		{
			for(int c = 0; c < 3; c++)
			{
				double coord = gmesh->NodeVec()[gmesh->ElementVec()[el]->NodeIndex(t)].Coord(c);
				node << coord << " ";
			}			
			node << std::endl;
			
			actualNode++;
			connectivity << " " << actualNode;
		}
		connectivity << std::endl;
		
		int elType;
		switch (gmesh->ElementVec()[el]->Type())
		{
			case (ETriangle):
			{
				elType = 5;
				break;				
			}
			case (EQuadrilateral):
			{
				elType = 9;
				break;				
			}
			case (ETetraedro):
			{
				elType = 10;
				break;				
			}
			case (EPiramide):
			{
				elType = 14;
				break;				
			}
			case (EPrisma):
			{
				elType = 13;
				break;				
			}
			case (ECube):
			{
				elType = 12;
				break;				
			}
			default:
			{
				elType = -1;//ElementType NOT Found!!!
				DebugStop();
				break;	
			}
		}
		type << elType << std::endl;
		nVALIDelements++;
	}
	node << std::endl;
	actualNode++;
	file << actualNode << " float" << std::endl << node.str();
	
	file << "CELLS " << nVALIDelements << " ";
	
	file << size << std::endl;
	file << connectivity.str() << std::endl;
	
	file << "CELL_TYPES " << nVALIDelements << std::endl;
	file << type.str();
}

void PrintGMeshVTKmaterialcopy(TPZGeoMesh * gmesh, std::ofstream &file, std::set<int> Material)
{
	file.clear();
	int nelements = gmesh->NElements();
	
	stringstream node, connectivity, type;
	
	//Header
	file << "# vtk DataFile Version 3.0" << std::endl;
	file << "TPZGeoMesh VTK Visualization" << std::endl;
	file << "ASCII" << std::endl << std::endl;
	
	file << "DATASET UNSTRUCTURED_GRID" << std::endl;
	file << "POINTS ";
	
	int actualNode = -1, size = 0, nVALIDelements = 0;
	
	for(int el = 0; el < nelements; el++)
	{				
		if(gmesh->ElementVec()[el]->Type() == EOned)//Exclude Lines and Arc3D
		{
			continue;
		}
		int mat = gmesh->ElementVec()[el]->MaterialId();
		bool found = !(Material.find(mat) == Material.end() );
		if(gmesh->ElementVec()[el]->HasSubElement() || !found)
		{
			continue;
		}
		
		int elNnodes = gmesh->ElementVec()[el]->NNodes();
		size += (1+elNnodes);
		connectivity << elNnodes;
		
		for(int t = 0; t < elNnodes; t++)
		{
			for(int c = 0; c < 3; c++)
			{
				double coord = gmesh->NodeVec()[gmesh->ElementVec()[el]->NodeIndex(t)].Coord(c);
				node << coord << " ";
			}			
			node << std::endl;
			
			actualNode++;
			connectivity << " " << actualNode;
		}
		connectivity << std::endl;
		
		int elType;
		switch (gmesh->ElementVec()[el]->Type())
		{
			case (ETriangle):
			{
				elType = 5;
				break;				
			}
			case (EQuadrilateral):
			{
				elType = 9;
				break;				
			}
			case (ETetraedro):
			{
				elType = 10;
				break;				
			}
			case (EPiramide):
			{
				elType = 14;
				break;				
			}
			case (EPrisma):
			{
				elType = 13;
				break;				
			}
			case (ECube):
			{
				elType = 12;
				break;				
			}
			default:
			{
				elType = -1;//ElementType NOT Found!!!
				DebugStop();
				break;	
			}
		}
		type << elType << std::endl;
		nVALIDelements++;
	}
	node << std::endl;
	actualNode++;
	file << actualNode << " float" << std::endl << node.str();
	
	file << "CELLS " << nVALIDelements << " ";
	
	file << size << std::endl;
	file << connectivity.str() << std::endl;
	
	file << "CELL_TYPES " << nVALIDelements << std::endl;
	file << type.str();
}
