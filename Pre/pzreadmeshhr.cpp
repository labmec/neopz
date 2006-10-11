/*****************************************************************************
 * O conteúdo desse arquivo é de propriedade do LabMeC-DES-FEC-UNICAMP e do
 * CENPES-Petrobras. 
 * O uso de qualquer parte ou do todo está condicionado à expressa autorização
 * dos proprietários.
 *****************************************************************************/
#include "pzreadmeshhr.h"

#include <pzgeoel.h>
#include <pzgmesh.h>
#include <pzcmesh.h>
#include <pzelasmat.h>
#include <pzmat2dlin.h>
#include <pzbndcond.h>
#include <pzcompel.h>
#include <pzlog.h>

#include <iostream>
#include <sstream>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.pre.tpzmeshreaderhr"));
#endif


TPZReadMeshHR::TPZReadMeshHR(const char* infInputFile): TPZReadMesh(infInputFile)
{}


TPZReadMeshHR::~TPZReadMeshHR()
{}


TPZCompMesh* TPZReadMeshHR::ReadMesh()
{
  TPZGeoMesh *gmesh = new TPZGeoMesh;
  std::string numberOf;
  removeComents (numberOf);
  int nnos = atoi (numberOf.c_str());
  ReadNodes (nnos, *gmesh);
  removeComents (numberOf);
  int nelem = atoi (numberOf.c_str());
  ReadElements(nelem, *gmesh);
  gmesh->BuildConnectivity();
  TPZCompMesh *cmesh = new TPZCompMesh (gmesh);
  removeComents (numberOf);
  int nmat = atoi (numberOf.c_str());
  ReadMaterials(nmat,*cmesh);
  removeComents (numberOf);
  int nbcs = atoi (numberOf.c_str());
  ReadBCs(nbcs,*cmesh);
  cmesh->AutoBuild();
  return cmesh;
}

void TPZReadMeshHR::removeComents (std::string &NumberOf)
{
  if (!fInputFile) return;
  while (fInputFile)
  {
    char buf[512];
    fInputFile.getline(buf,512);
    std::string aux (buf);
    int pos = aux.find (":");
    if ((pos > -1 && pos < (int)aux.size()) || !aux.size()) continue; //encontrou um comentário na linha
    NumberOf = aux;
    break;
  }
}


void TPZReadMeshHR::ReadNodes (int NNos, TPZGeoMesh & GMesh)
{
  int i,c;
  int id;
  TPZVec <REAL> coord(3,0.);
  for (i=0;i<NNos;i++)
  {
    fInputFile >> id;
    for (c=0;c<3;c++) fInputFile >> coord[c];
    int index = GMesh.NodeVec().AllocateNewElement();
    GMesh.NodeVec()[index] = TPZGeoNode(id,coord,GMesh);
  }
}


void TPZReadMeshHR::ReadElements (int NElem, TPZGeoMesh & GMesh)
{
  int i,c;
  int id,type,matid;
  for (i=0;i<NElem;i++)
  {
    fInputFile >> id >> type >> matid;
    int nnos = -1;
    MElementType etype;
    switch (type)
    {
      case (0):
        nnos = 1;
        etype = EPoint;
        break;
      case (1):
        nnos = 2;
        etype = EOned;
        break;
      case (2):
        nnos = 3;
        etype = ETriangle;
        break;
      case (3):
        nnos = 4;
        etype = EQuadrilateral;
        break;
      case (4):
        nnos = 4;
        etype = ETetraedro;
        break;
      case (5):
        nnos = 5;
        etype = EPiramide;
        break;
      case (6):
        nnos = 6;
        etype = EPrisma;
        break;
      case (7):
        nnos = 8;
        etype = ECube;
        break;
      default:
        std::stringstream sout;
#ifndef WINDOWS
        sout << __PRETTY_FUNCTION__;
#endif
        sout << "Não sei que elemento " << type << " é esse indicado para o elemento " << id;
#ifdef LOG4CXX
        LOGPZ_WARN (logger,sout.str().c_str());
#else
        std::cout << sout.str().c_str() << std::endl;
#endif
        continue;
    }
    TPZVec<int> corneridx(nnos,-1);
    int id, idx;
    for (c=0;c<nnos;c++)
    {
      fInputFile >> id;
      idx = GetNodeIndex(&GMesh,id);
      corneridx [c] = idx;
    }
    GMesh.CreateGeoElement(etype,corneridx,matid,id,0);
  }
}


void TPZReadMeshHR::ReadMaterials (int NMat, TPZCompMesh & CMesh)
{
  int i;
  int id;
  double e, nu, px, py;
  for (i=0;i<NMat;i++)
  {
    fInputFile >> id >> e >> nu >> px >> py;
    TPZMaterial *mat = new TPZElasticityMaterial(id,e,nu,px,py,0);
    CMesh.InsertMaterialObject(mat);
  }
}


void TPZReadMeshHR::ReadBCs (int NMat, TPZCompMesh & CMesh)
{
  int i;
  int id, type;
  for (i=0;i<NMat;i++)
  {
    fInputFile >> id >> type;
    TPZMaterial *mat = CMesh.MaterialVec()[0];//cMesh.FindMaterial(gel->MaterialId());
    if (!mat)
    {
      std::stringstream sout;
#ifndef WINDOWS
      sout << __PRETTY_FUNCTION__;
#endif
      sout << "\tNão encontrei material na malha!";
#ifdef LOG4CXX
      LOGPZ_ERROR (logger, sout.str().c_str());
#else
      cout << sout.str().c_str() << endl;
#endif
      continue;
    }
    TPZFMatrix val1(3,3,0.),val2(3,1,0.);
    fInputFile >> val1 (0,0) >> val1(0,1) >> val1(0,2)
        >> val1 (1,0) >> val1(1,1) >> val1(1,2)
        >> val1 (2,0) >> val1(2,1) >> val1(2,2);
    fInputFile >> val2(0,0) >> val2(1,0) >> val2(2,0);
    TPZMaterial *bnd;
    bnd = mat->CreateBC (id,type,val1,val2);
/*    switch (type)
    {
      case (0) :
      { // Dirichlet
        bnd = mat->CreateBC (id,0,val1,val2);
        break;
      }
      case (1) :
      { // Neumann
        bnd = mat->CreateBC (id,1,val1,val2);
        break;
      }
      case (2) :
      {// Mista
        bnd = mat->CreateBC (id,2,val1,val2);
        break;
      }
      default:
      {
        std::stringstream sout;
#ifndef WINDOWS
        sout << __PRETTY_FUNCTION__;
#endif
        sout << "\tBC tipo " << type << " não identificada!";
#ifdef LOG4CXX
        LOGPZ_WARN (logger, sout.str().c_str());
#else
        cout << sout.str().c_str() << endl;
#endif
        continue;
      }
    }*/
    CMesh.InsertMaterialObject(bnd);
  }

  //Materiais da Fratura
  TPZFMatrix val1(3,3,0.),val2(3,1,0.);
  //val1(1,1) = 10000000000000.0;
  //val2(1,0) = 1.;
  TPZBndCond *bndFrac = new TPZBndCond (CMesh.MaterialVec()[0],-100,1,val1,val2);
  CMesh.InsertMaterialObject(bndFrac);
//    val2(1,0) = -1.;
  bndFrac = new TPZBndCond (CMesh.MaterialVec()[0],-101,1,val1,val2);
  CMesh.InsertMaterialObject(bndFrac);
}


int TPZReadMeshHR::GetNodeIndex(TPZGeoMesh *GMesh,int Id)
{
  TPZAdmChunkVector<TPZGeoNode> nodeVec = GMesh->NodeVec();
  int vid,index,size = nodeVec.NElements();
  if (Id < size) index = Id;
  else index = size - 1;
  vid = nodeVec[index].Id();
  if (vid == Id) return index;

  int i;
  for (i=index-1;i>-1;i--)
  {
    vid = nodeVec[i].Id();
    if (vid == Id) return i;
  }
  std::stringstream sout;
#ifndef WINDOWS
        sout << __PRETTY_FUNCTION__;
#endif
        sout << " Nó " << Id << " não encontrado!";
#ifdef LOG4CXX
        LOGPZ_WARN (logger, sout.str().c_str());
#else
        cout << sout.str().c_str() << endl;
#endif

  return -1;
}
