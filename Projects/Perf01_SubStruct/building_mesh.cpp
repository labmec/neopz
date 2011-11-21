#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "building_mesh.h"
#include "pzelast3d.h"
#include "pzbndcond.h"
#include <string>

using namespace std;

void InsertElasticity(TPZAutoPointer<TPZCompMesh> mesh)
{
  mesh->SetDimModel(3);
  int nummat = 1;
  REAL E = 1.e6;
  REAL poisson = 0.3;
  TPZManVector<REAL> force(3,0.);
  force[1] = 20.;
  TPZElasticity3D *elast = new TPZElasticity3D(nummat,E,poisson,force);
  TPZAutoPointer<TPZMaterial> elastauto(elast);
  TPZFMatrix val1(3,3,0.),val2(3,1,0.);
  TPZBndCond *bc = elast->CreateBC(elastauto, -1, 0, val1, val2);
  TPZAutoPointer<TPZMaterial> bcauto(bc);
  mesh->InsertMaterialObject(elastauto);
  mesh->InsertMaterialObject(bcauto);
}

TPZGeoMesh* BuildBuildingMesh(const char* input_filename)
{
  //int nBCs = 1;
  int numnodes=-1;
  int numelements=-1;
	
  {
    bool countnodes = false;
    bool countelements = false;
		
    ifstream read (input_filename);		
    while(read) {
      char buf[1024];
      read.getline(buf, 1024);
      std::string str(buf);

      if(str == "Coordinates") 
	countnodes = true;
      if(str == "end coordinates") 
	countnodes = false;
      if(countnodes) 
	numnodes++;
      if(str == "Elements") 
	countelements = true;
      if(str == "end elements") 
	countelements = false;
      if(countelements) 
	numelements++;
    }
  }
	
  TPZGeoMesh * gMesh = new TPZGeoMesh;
  gMesh -> NodeVec().Resize(numnodes);
	
  TPZVec <int> TopolTetra(4);
	
  const int Qnodes = numnodes;
  TPZVec <TPZGeoNode> Node(Qnodes);
	
  //setting nodes coords
  int nodeId = 0, elementId = 0, matElId = 1;
  double nodecoordX , nodecoordY , nodecoordZ ;
  char buf[1024];
  ifstream read;
  
  read.open(input_filename);
  read.getline(buf, 1024);
  read.getline(buf, 1024);
  std::string str(buf);
  for(int in=0; in<numnodes; in++) { 
    read >> nodeId;
    read >> nodecoordX;
    read >> nodecoordY;
    read >> nodecoordZ;
    Node[nodeId-1].SetNodeId(nodeId);
    Node[nodeId-1].SetCoord(0,nodecoordX);
    Node[nodeId-1].SetCoord(1,nodecoordY);
    Node[nodeId-1].SetCoord(2,nodecoordZ);
    gMesh->NodeVec()[nodeId-1] = Node[nodeId-1];
  }
  read.close();
  

  read.open(input_filename);
		
  int m = numnodes+5;
  for(int l=0; l<m; l++) {
    read.getline(buf, 1024);
  }
		
  int el;
  int matBCid = -1;
  //std::set<int> ncoordz; //jeitoCaju
  for(el=0; el<numelements; el++) {
    read >> elementId;
    read >> TopolTetra[0]; //node 1
    read >> TopolTetra[1]; //node 2
    read >> TopolTetra[2]; //node 3
    read >> TopolTetra[3]; //node 4
    
    // O GID comeca com 1 na contagem dos nodes, e nao zero como no PZ, assim o node 1 na verdade é o node 0
    TopolTetra[0]--;
    TopolTetra[1]--;
    TopolTetra[2]--;
    TopolTetra[3]--;
    
    int index;
    TPZGeoEl * tetra = gMesh->CreateGeoElement(ETetraedro, TopolTetra, matElId, index);
    
    // Colocando as condicoes de contorno
    TPZVec <TPZGeoNode> Nodefinder(4);
    TPZVec <REAL> nodecoord(3);
    TPZVec<int> ncoordzVec(0); int sizeOfVec = 0;
    //ncoordz.clear(); //jeitoCaju
    for (int i = 0; i < 4; i++) {
      Nodefinder[i] = gMesh->NodeVec()[TopolTetra[i]];
      Nodefinder[i].GetCoordinates(nodecoord);
      if (nodecoord[2] == 0.) {
	//ncoordz.insert(TopolTetra[i]); //jeitoCaju
	sizeOfVec++;
	ncoordzVec.Resize(sizeOfVec);
	ncoordzVec[sizeOfVec-1] = TopolTetra[i];
      }
    }
    //if(ncoordz.size() == 3) //jeitoCaju
    if(ncoordzVec.NElements() == 3) {
      int lado = tetra->WhichSide(ncoordzVec);
      TPZGeoElSide tetraSide(tetra, lado);
      TPZGeoElBC(tetraSide,matBCid);		
    }
  }
		
  gMesh->BuildConnectivity();
  	
  return gMesh;
}


