//
// C++ Implementation: %{MODULE}
//
// Description:
//
//
// Author: %{AUTHOR} <%{EMAIL}>, (C) %{YEAR}
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include <fstream>
#include <iostream>
#include <string>

#include <pzgeoel.h>
#include <pzgmesh.h>
#include <pzcmesh.h>
#include <pzelasmat.h>
#include <pzmat2dlin.h>
#include <pzbndcond.h>

void removeComents (std::istream &file, std::string &numberOf)
{
  if (!file) return;
  while (file)
  {
    char buf[512];
    file.getline(buf,512);
    std::string aux (buf);
    int pos = aux.find (":");
    if ((pos > -1 && pos < (int)aux.size()) || !aux.size()) continue; //encontrou um comentário na linha
    numberOf = aux;
    break;
  }
}

void ReadNodes (std::istream & file, int nnos, TPZGeoMesh & gMesh)
{
  int i,c;
  int id;
  TPZVec <REAL> coord(3,0.);
  for (i=0;i<nnos;i++)
  {
    file >> id;
    for (c=0;c<3;c++) file >> coord[c];
    int index = gMesh.NodeVec().AllocateNewElement();
    gMesh.NodeVec()[index] = TPZGeoNode(id,coord,gMesh);
  }
}

void ReadElements (std::istream & file, int nelem, TPZGeoMesh & gMesh)
{
  int i,c;
  int id,type,matid;
  for (i=0;i<nelem;i++)
  {
    file >> id >> type >> matid;
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
        std::cout << "Não sei que elemento " << type << " é esse indicado para o elemento " << id << std::endl;
        return;
    }
    TPZVec<int> corneridx(nnos,-1);
    for (c=0;c<nnos;c++) file >> corneridx [c];
    gMesh.CreateGeoElement(etype,corneridx,matid,id,1);
  }
}

void ReadMaterials (std::istream & file, int nmat, TPZCompMesh & cMesh)
{
  int i;
  int id;
  double e, nu, px, py;
  for (i=0;i<nmat;i++)
  {
    file >> id >> e >> nu >> px >> py;
    TPZAutoPointer<TPZMaterial> mat = new TPZElasticityMaterial(id,e,nu,px,py);
    cMesh.InsertMaterialObject(mat);
  }
}

void ReadBCs (std::istream & file, int nmat, TPZCompMesh & cMesh)
{
  int i;
  int id, type, elid, side;
  for (i=0;i<nmat;i++)
  {
    file >> id >> type >> elid >> side;
    TPZGeoEl *gel = cMesh.Reference()->FindElement(elid);
    if (!gel)
    {
      std::cout << "Não encontrei o elemento cujo id é " << elid << std::endl;
      continue;
    }
    TPZGeoElBC heman_1(gel,side,type,*(cMesh.Reference()));
    TPZAutoPointer<TPZMaterial> mat = cMesh.FindMaterial(gel->MaterialId());
    if (!mat)
    {
      std::cout << "Não encontrei o material cujo id é " << gel->MaterialId() << std::endl;
      continue;
    }
    if (type == -1)
    { // Dirichlet
      TPZFMatrix val1(3,3,0.),val2(3,1,0.);
      TPZAutoPointer<TPZMaterial> bnd = mat->CreateBC (mat,type,0,val1,val2);
      cMesh.InsertMaterialObject(bnd);
    }
    else
    { // Neumann
      TPZFMatrix val1(3,3,0.),val2(3,1,1.);
      TPZAutoPointer<TPZMaterial> bnd = mat->CreateBC (mat,type,1,val1,val2);
      cMesh.InsertMaterialObject(bnd);
    }
  }
}

TPZCompMesh * ReadMesh(std::istream &file)
{
  TPZGeoMesh *gmesh = new TPZGeoMesh;
  
  std::string numberOf;
  removeComents (file,numberOf);
  int nnos = atoi (numberOf.c_str());
  ReadNodes (file,nnos, *gmesh);
  removeComents (file,numberOf);
  int nelem = atoi (numberOf.c_str());
  ReadElements(file,nelem, *gmesh);
  
  gmesh->BuildConnectivity();
  TPZCompMesh *cmesh = new TPZCompMesh (gmesh);

  removeComents (file,numberOf);
  int nmat = atoi (numberOf.c_str());
  ReadMaterials(file,nmat,*cmesh);
  removeComents (file,numberOf);
  int nbcs = atoi (numberOf.c_str());
  ReadBCs(file,nbcs,*cmesh);

  cmesh->AutoBuild();
  return cmesh;
}


//===================================================================
//==================Prints the Mesh in DX format=====================
//===================================================================
void WriteElement (TPZGeoEl *el,int elindex, std::ofstream &arq,TPZVec<int> &elementtype){
  int ncon = el->NNodes();
  elementtype[elindex] = ncon;
  switch (ncon)  {
  case (2) : {
    //rib
    
    int ni = el->NodeIndex(0);
    int nf = el->NodeIndex(1);
    arq << ni << "\t" << nf << "\t"
        << ni << "\t" << nf << "\t"
        << ni << "\t" << nf << "\t"
        << ni << "\t" << nf << std::endl;
    break;
  }
  case (3) : {
    //triangle
    int n0 = el->NodeIndex(0);
    int n1 = el->NodeIndex(1);
    int n2 = el->NodeIndex(2);
    arq << n0 << "\t" << n1 << "\t"
        << n2 << "\t" << n2 << "\t"
        << n0 << "\t" << n1 << "\t"
        << n2 << "\t" << n2 << std::endl;
    break;
  }
  case (4) : {
    if (el->Dimension() == 2){
      //quad
      int n0 = el->NodeIndex(0);
      int n1 = el->NodeIndex(1);
      int n2 = el->NodeIndex(3);
      int n3 = el->NodeIndex(2);
      arq << n0 << "\t" << n1 << "\t"
          << n2 << "\t" << n3 << "\t"
          << n0 << "\t" << n1 << "\t"
          << n2 << "\t" << n3 << std::endl;
    }else{
      //tetrahedre
      int n0 = el->NodeIndex(0);
      int n1 = el->NodeIndex(1);
      int n2 = el->NodeIndex(2);
      int n3 = el->NodeIndex(3);
      arq << n0 << "\t" << n1 << "\t"
          << n2 << "\t" << n2 << "\t"
          << n3 << "\t" << n3 << "\t"
          << n3 << "\t" << n3 << std::endl;
    }
    break;
  }
  case (5) : {
    //pyramid
    int n0 = el->NodeIndex(0);
    int n1 = el->NodeIndex(1);
    int n2 = el->NodeIndex(3);
    int n3 = el->NodeIndex(2);
    int n4 = el->NodeIndex(4);
    arq << n0 << "\t" << n1 << "\t"
        << n2 << "\t" << n3 << "\t"
        << n4 << "\t" << n4 << "\t"
        << n4 << "\t" << n4 << std::endl;
    break;
  }
  case (6) : {
    //pyramid
    int n0 = el->NodeIndex(0);
    int n1 = el->NodeIndex(1);
    int n2 = el->NodeIndex(2);
    int n3 = el->NodeIndex(3);
    int n4 = el->NodeIndex(4);
    int n5 = el->NodeIndex(5);
    arq << n0 << "\t" << n1 << "\t"
        << n2 << "\t" << n2 << "\t"
        << n3 << "\t" << n4 << "\t"
        << n5 << "\t" << n5 << std::endl;
    break;
  }
  case (8) : {
    int n0 = el->NodeIndex(0);
    int n1 = el->NodeIndex(1);
    int n2 = el->NodeIndex(3);
    int n3 = el->NodeIndex(2);
    int n4 = el->NodeIndex(4);
    int n5 = el->NodeIndex(5);
    int n6 = el->NodeIndex(7);
    int n7 = el->NodeIndex(6);
    arq << n0 << "\t" << n1 << "\t"
        << n2 << "\t" << n3 << "\t"
        << n4 << "\t" << n5 << "\t"
        << n6 << "\t" << n7 << std::endl;
    break;
  }
  default:
    std::cout << "Erro..." << std::endl;
  }
  return;
}


void WriteMesh(TPZGeoMesh *mesh,std::ofstream &arq, int matindex){
  arq << "object 1 class array type float rank 1 shape 3 items ";
  arq << mesh->NodeVec().NElements() << " data follows" << std::endl;
  int i;
  //Print Nodes
  for (i=0;i<mesh->NodeVec().NElements(); i++){
    TPZGeoNode *node = &mesh->NodeVec()[i];
    arq /*<< node->Id() << "\t"*/
      << node->Coord(0) << "\t"
      << node->Coord(1) << "\t"
      << node->Coord(2) << std::endl;
  }
  int numelements = 0;
  for (i=0;i<mesh->ElementVec().NElements();i++){
    TPZGeoEl *el = mesh->ElementVec()[i];
    if ( !el || el->Dimension() != 3) continue;
    if(el->MaterialId() == matindex) numelements++;
  }
  arq << "object 2 class array type integer rank 1 shape 8 items ";
  arq << numelements << " data follows" << std::endl;
  TPZVec<int> elementtype(mesh->ElementVec().NElements(),0);
  for (i=0;i<mesh->ElementVec().NElements();i++){
    TPZGeoEl *el = mesh->ElementVec()[i];
    if ( !el || el->MaterialId() != matindex || el->Dimension() != 3) continue;
    WriteElement (el,i,arq,elementtype);
  }
  arq << "attribute \"element type\" string \"cubes\"" << std::endl
    << "attribute \"ref\" string \"positions\"" << std::endl;
  arq << "object 3 class array type integer rank 0 items ";
  arq << numelements << " data follows" << std::endl;
  for (i=0;i<mesh->ElementVec().NElements();i++){
    TPZGeoEl *el = mesh->ElementVec()[i];
    if ( !el || el->MaterialId() != matindex || el->Dimension()!= 3) continue;
    arq << elementtype[i] << std::endl;
  }
  arq << "attribute \"dep\" string \"connections\"" << std::endl;
  arq << "object 4 class field" << std::endl
    << "component \"positions\" value 1" << std::endl
    << "component \"connections\" value 2" << std::endl
    << "component \"data\" value 3" << std::endl;
}

void WriteElementMesh(TPZGeoMesh *mesh,std::ofstream &arq, int matindex,int eltype){
  arq << "object 1 class array type float rank 1 shape 3 items ";
  arq << mesh->NodeVec().NElements() << " data follows" << std::endl;
  int i;
  //Print Nodes
  for (i=0;i<mesh->NodeVec().NElements(); i++){
    TPZGeoNode *node = &mesh->NodeVec()[i];
    arq /*<< node->Id() << "\t"*/
      << node->Coord(0) << "\t"
      << node->Coord(1) << "\t"
      << node->Coord(2) << std::endl;
  }
  int numelements = 0;
  for (i=0;i<mesh->ElementVec().NElements();i++){
    TPZGeoEl *el = mesh->ElementVec()[i];
    if ( !el || el->Type() != eltype) continue;
    if(el->MaterialId() == matindex) numelements++;
  }
  arq << "object 2 class array type integer rank 1 shape 8 items ";
  arq << numelements << " data follows" << std::endl;
  TPZVec<int> elementtype(mesh->ElementVec().NElements(),0);
  for (i=0;i<mesh->ElementVec().NElements();i++){
    TPZGeoEl *el = mesh->ElementVec()[i];
    if ( !el || el->MaterialId() != matindex || el->Type() != eltype || el->HasSubElement()) continue;
    WriteElement (el,i,arq,elementtype);
  }
  arq << "attribute \"element type\" string \"cubes\"" << std::endl
    << "attribute \"ref\" string \"positions\"" << std::endl;
  arq << "object 3 class array type integer rank 0 items ";
  arq << numelements << " data follows" << std::endl;
  for (i=0;i<mesh->ElementVec().NElements();i++){
    TPZGeoEl *el = mesh->ElementVec()[i];
    if ( !el || el->MaterialId() != matindex || el->Type() != eltype  || el->HasSubElement()) continue;
    arq << elementtype[i] << std::endl;
  }
  arq << "attribute \"dep\" string \"connections\"" << std::endl;
  arq << "object 4 class field" << std::endl
    << "component \"positions\" value 1" << std::endl
    << "component \"connections\" value 2" << std::endl
    << "component \"data\" value 3" << std::endl;
}
