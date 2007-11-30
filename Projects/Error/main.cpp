#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzgeoel.h"
#include "pzcompel.h"
#include "pzanalysis.h"
#include "pzvec.h"
#include "pzerror.h"
#include "tmbadapinterface.h"
#include <pzmat1dlin.h>
#include <pzmattest3d.h>
#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"

//Shape Classes
#include <pzshapepoint.h>
#include <pzshapelinear.h>
#include <pzshapetriang.h>
#include <pzshapequad.h>
#include <pzshapetetra.h> 
#include <pzshapepiram.h>
#include <pzshapeprism.h>
#include <pzshapecube.h>

//Geometry Classes
#include <pzgeopoint.h>
#include <TPZGeoLinear.h>
#include <pzgeotriangle.h>
#include <pzgeoquad.h>
#include <pzgeotetrahedra.h>
#include <pzgeopyramid.h>
#include <pzgeoprism.h>
#include <TPZGeoCube.h> 

#include <time.h>
#include <stdio.h>
#include <fstream>
using namespace std;

int gDebug = 2;

void ReadControl (ifstream &arq, TPZVec<int> &meshsize, int &nstate,
                    TPZVec<int> &dimstate, TPZVec<int> &usestate,
                    TPZVec<int> &erind, TPZVec<int> &erantype,
                    TPZVec<REAL> &maxerror, TPZVec<REAL> &minerr,int &nsteps);

TPZCompMesh *ReadMesh(ifstream &arq,
                      TPZVec<int> &meshside);

void ReadSolution(ifstream &arq,
                  TPZVec<REAL> &sol,
                  TPZCompMesh *cmesh,
                  int &nstate,
                  TPZVec<int> &statedim);

void WriteCompMesh(TPZCompMesh *cmesh,ofstream &arq);

void EvaluateSolution(TPZVec<REAL> &pt, TPZVec<REAL> &sol);

TPZGeoEl *GeoElement(TPZGeoMesh *gmesh, int ncon,
                      TPZVec<int> &nodind, int material,
                      int index);

TPZCompMesh *ReadElementsMesh();
TPZCompMesh *TetraMesh();

void WriteMesh(TPZGeoMesh *mesh,ofstream &arq);

void WriteElement (TPZGeoEl *el,ofstream &arq,TPZVec<int> &elementtype);

TPZGeoEl *GeoElementRefPattern(TPZGeoMesh *gmesh,
                               int ncon,
                               TPZVec<int> &nodind,
                               int material,
                               int index,
                               TPZVec<TPZRefPattern *> refpattern);


int main(){
#define HUGE_DEBUG
  int i;
  TPZCompMesh *cmesh;
  TPZVec <REAL> solution;
  int nstate, nsteps;
  TPZVec<int> meshsize(3,0);
  TPZVec<int> dimstate(3,0);
  TPZVec<int> usestate(3,0);
  TPZVec<int> erind(3,0);
  TPZVec<int> erantype(3,0);
  TPZVec <REAL> maxerror;
  TPZVec <REAL> minerror;
  ifstream arq_control ("/compile/cesar/NeoPZ/Projects/Error/control.txt");
  ReadControl(arq_control,meshsize,nstate,dimstate,usestate,erind,erantype,maxerror,minerror,nsteps);
  //  cout << "Reading Mesh\n";
  ifstream arq_mesh ("/compile/cesar/NeoPZ/Projects/Error/Mesh.data");
  //cmesh = ReadMesh(arq_mesh,meshsize);
  //cmesh = ReadElementsMesh();
  cmesh = TetraMesh();
  cout << "Reading Solution\n";
  ifstream arq_solution ("/compile/cesar/NeoPZ/Projects/Error/Solution.data");
  TPZStack<char *> scalnames,vecnames;
  scalnames.Push("POrder");
  //cmesh->Print(cout);
  //cmesh->Reference()->Print(cout);
  for (i=0;i<nsteps;i++){
    cout << "\n\nEntering step...: " << i << endl;
    //Just for visualization purposes...
    TPZAnalysis an (cmesh);
    char buf [256];
    sprintf(buf,"htest%d.dx",i);
    //an.DefineGraphMesh(3,scalnames,vecnames,buf);
    //an.PostProcess(0,3);
    ofstream arq(buf);
    WriteMesh(cmesh->Reference(),arq);

    ReadSolution(arq_solution,solution,cmesh,nstate,dimstate);
    //  cout << solution << endl;
    TPZCompMesh *adaptmesh;
    TMBAdaptInterface adapt(cmesh,nstate,dimstate,usestate,solution);
    adapt.SetMaxMinError(maxerror,minerror);
    adaptmesh = adapt.GetAdaptedMesh(erind,erantype,true,1,1,0);
    //cmesh->Reference()->ResetReference();
    //   cmesh->LoadReferences();
    //delete cmesh;
    cmesh = adaptmesh;
    //cmesh->Reference()->Print(cout);
  }
  //Just for visualization purposes...
  TPZAnalysis an (cmesh);
  char buf [256];
  sprintf(buf,"htest%d.dx",i);
  ofstream arq (buf);
  WriteMesh(cmesh->Reference(),arq);
//  an.DefineGraphMesh(3,scalnames,vecnames,buf);
//  an.PostProcess(0,3);

//  WriteCompMesh(cmesh,cout);
  delete cmesh;
  cout << "End..." << endl;
  return 0;
}

void ReadControl (ifstream &arq, TPZVec<int> &meshsize, int &nstate,
                  TPZVec<int> &dimstate, TPZVec<int> &usestate,
                  TPZVec<int> &erind, TPZVec<int> &erantype,
                  TPZVec<REAL> &maxerror, TPZVec<REAL> &minerror,int &nsteps){
  cout << "================================================\n";
  cout << "CFD Group - Embraer / Fapesp\nAdaptive M�dulus\n";
  cout << "================================================\n\n";
  cout << endl << endl;
  cout << "Reading control data\n";
  arq >> meshsize[0] >> meshsize[1] >> meshsize[2];
  arq >> nstate;
  dimstate.Resize(nstate);
  int i,use;
  for (i=0;i<nstate;i++)	arq >> dimstate[i];
  arq >> use;
  usestate.Resize(use);
  erind.Resize(use);
  erantype.Resize(use);
  maxerror.Resize(use);
  minerror.Resize(use);
  
  for (i=0;i<use;i++)	arq >> usestate[i];
  for (i=0;i<use;i++)	arq >> erind[i];
  for (i=0;i<use;i++) arq >> erantype[i];
  for (i=0;i<use;i++) arq >> maxerror[i];
  for (i=0;i<use;i++) arq >> minerror[i];
  arq >> nsteps;

  //Print the control data
  cout << "Data Control Information: \n";
  cout << "Mesh Size: - x = " <<  meshsize[0]
        << "  y = "	<< meshsize[1]
        << "  z = " << meshsize[2] << endl;
  cout << "Number of state variables in solution vector = " << nstate << endl;
  cout << "Dimension of each state variable = " << dimstate << endl;
  cout << "State variable to use in analysis = " << usestate << endl;
  cout << "Error indicator for each state variable = " << erind << endl;
  cout << "Error analysis type for each state variable  = " << erantype << endl;
}

TPZCompMesh *ReadMesh(ifstream &arq, TPZVec<int> &meshsize){
  int nx = meshsize[0];
  int ny = meshsize[1];
  int nz = meshsize[2];
  int i,j,k;
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  TPZGeoEl * elvec[(const int)((nx-1)*(ny-1)*(nz-1))];

  TPZVec <REAL> coord (3,0.);
  TPZVec <int> connect(8,0);
  REAL lx = 1.;
  REAL ly = 1.;
  REAL lz = 1.;
  int id, index;
  //Nodes initialization
  for(i = 0; i < nx; i++){
    for(j = 0; j < ny; j++){
      for(k = 0; k < nz; k++){
        id = (i)*nz*ny + (j)*nz + k;
        coord[0] = (i)*lx/(nx - 1);
        coord[1] = (j)*ly/(ny - 1);
        coord[2] = (k)*lz/(nz - 1);
        //cout << coord << endl;
        index = gmesh->NodeVec().AllocateNewElement();
        gmesh->NodeVec()[index] = TPZGeoNode(id,coord,*gmesh);
      }
    }
  }

  //Element connectivities
  TPZRefPattern *unifcube = new TPZRefPattern ("/home/pos/cesar/RefPattern/Hexa_Unif.rpt");
  for(i = 0; i < (nx - 1); i++){
    for(j = 0; j < (ny - 1); j++){
      for(k = 0; k < (nz - 1); k++){
        index = (i)*(nz - 1)*(ny - 1) + (j)*(nz - 1) + k;
        connect[0] = (i)*nz*ny + (j)*nz + k;
        connect[1] = connect[0]+(ny)*(nz);
        connect[2] = connect[1]+(nz);
        connect[3] = connect[0]+(nz);
        connect[4] = connect[0] + 1;
        connect[5] = connect[1] + 1;
        connect[6] = connect[2] + 1;
        connect[7] = connect[3] + 1;
        //cout << connect << endl;
//        elvec[index] = gmesh->CreateGeoElement(ECube,connect,1,id);
        TPZGeoElRefPattern <TPZShapeCube,TPZGeoCube> *gel =
             new TPZGeoElRefPattern <TPZShapeCube,TPZGeoCube> (index,connect,1,*gmesh,unifcube);
        elvec[index] = gel;
      }
    }
  }
  //Generate neighborhod information
  gmesh->BuildConnectivity();
  gmesh->Print(cout);
  //Create computational mesh

  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  TPZMaterial *mat;
  mat = new TPZMaterialTest3D (1);
  cmesh->InsertMaterialObject(mat);

  cmesh->AutoBuild();
  return cmesh;
}

void ReadSolution(ifstream &arq, TPZVec<REAL> &sol, TPZCompMesh *cmesh, int &nstate, TPZVec<int> &dimstate){

  int i,j,totaldim=0;
  for(i=0;i<nstate;i++) totaldim += dimstate[i];
  TPZVec<REAL> pt(3,0.);
  TPZVec<REAL> coord(3,0.);
  TPZVec<REAL> auxsol(totaldim,0.);

  int iter = 0;
  int nel = cmesh->NElements();
  int solsize = totaldim * nel;
  sol.Resize(solsize);
  sol.Fill(0.);

  for(i=0; i<nel; i++){
    TPZCompEl *el = cmesh->ElementVec()[i];
    if (!el) continue;
    el->Reference()->CenterPoint(el->Reference()->NSides()-1,pt);
    el->Reference()->X(pt,coord);
    EvaluateSolution(coord,auxsol);
    for (j=0;j<totaldim;j++){
      sol[iter] = auxsol[j];
      iter++;
    }
  }
}

void EvaluateSolution(TPZVec<REAL> &pt, TPZVec<REAL> &sol){
  REAL x = pt[0];
  REAL y = pt[1];
  REAL z = pt[2];

  REAL x2 = x * x;
  REAL y2 = y * y;
  REAL z2 = z * z;
  REAL delta = 0.1;
  REAL aux = delta + x2 + y2 + z2;

  //Velocity;
  sol[0] = sqrt((4.*x2)/(aux*aux*aux*aux) + (4.*y2)/(aux*aux*aux*aux) + (4.*z2)/(aux*aux*aux*aux));
  //Vx
  sol[1] = (-2.*x)/(aux*aux);
  //Vy
  sol[2] = (-2.*y)/(aux*aux);
  //Vz
  sol[3] = (-2.*z)/(aux*aux);
  //dV/dx
  sol[4] = ((-32.* x2 *x) / (aux*aux*aux*aux*aux) - (32.*x*y2)/(aux*aux*aux*aux*aux) -
     (32.*x*z2)/(aux*aux*aux*aux*aux) + (8.*x)/(aux*aux*aux*aux))/
   (2.*sqrt((4.*x2)/(aux*aux*aux*aux) + (4.*y2)/(aux*aux*aux*aux) + (4.*z2)/(aux*aux*aux*aux)));
  //dV/dy
  sol[5] = ((-32.*x2*y)/(aux*aux*aux*aux*aux) - (32.* y2 * y)/(aux*aux*aux*aux*aux) -
     (32.*y*z2)/(aux*aux*aux*aux*aux) + (8.*y)/(aux*aux*aux*aux))/
   (2.*sqrt((4.*x2)/(aux*aux*aux*aux) + (4.*y2)/(aux*aux*aux*aux) + (4.*z2)/(aux*aux*aux*aux)));
  //dV/dz
  sol[6] = ((-32.*x2*z)/(aux*aux*aux*aux*aux) - (32.*y2*z)/(aux*aux*aux*aux*aux) -
     (32.*z2*z)/(aux*aux*aux*aux*aux) + (8.*z)/(aux*aux*aux*aux))/
   (2.*sqrt((4.*x2)/(aux*aux*aux*aux) + (4.*y2)/(aux*aux*aux*aux) + (4.*z2)/(aux*aux*aux*aux)));
}

void WriteCompMesh(TPZCompMesh *cmesh,ostream &out){
  cmesh->Print(cout);
  return;
}


TPZCompMesh *ReadElementsMesh(){

  REAL Coord [18][3] = {
    {0.,0.,0.},{0.5,0.,0.},{1.,0.,0.},
    {0.,1.,0.},{0.5,1.,0.},{1.,1.,0.},

    {0.,0.,0.5},{0.5,0.,0.5},{1.,0.,0.5},
    {0.,1.,0.5},{0.5,1.,0.5},{1.,1.,0.5},

    {0.,0.,1.},{0.5,0.,1.},{1.,0.,1.},
    {0.,1.,1.},{0.5,1.,1.},{1.,1.,1.}
  };

  int NodesPerEl [11] = { 8,
                          5,5,5,
                          6,6,
                          4,4,4,4,4};

  int Connects [11][8] = {
    {0,1,4,3,6,7,10,9},

    {6,7,10,9,16,-1,-1,-1},
    {6,7,13,12,16,-1,-1,-1},
    {6,9,15,12,16,-1,-1,-1},

    {1,2,5,7,8,11,-1,-1},
    {1,5,4,7,11,10,-1,-1},

    {7,8,11,14,-1,-1,-1,-1},
    {7,11,10,16,-1,-1,-1,-1},
    {16,17,14,11,-1,-1,-1,-1},
    {16,14,13,7,-1,-1,-1,-1},
    {7,11,16,14,-1,-1,-1,-1}
  };

  int i,j;
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  TPZGeoEl * elvec[11];

  TPZVec <REAL> coord (3,0.);

  int index;

  //Nodes initialization
  for(i = 0; i < 18; i++){
    for(j=0;j<3;j++){
      coord[j] = Coord[i][j];
    }
    index = gmesh->NodeVec().AllocateNewElement();
    gmesh->NodeVec()[index] = TPZGeoNode(i,coord,*gmesh);
  }

  TPZVec<TPZRefPattern *> refinement_Patterns;
  refinement_Patterns.Resize(8);
  refinement_Patterns[0] = 0;
  refinement_Patterns[1] = new TPZRefPattern("/home/pos/cesar/RefPattern/Line_Unif_Side_2.rpt");
  refinement_Patterns[2] = new TPZRefPattern("/home/pos/cesar/RefPattern/Triang_Unif.rpt");
  refinement_Patterns[3] = new TPZRefPattern("/home/pos/cesar/RefPattern/Quad_Unif.rpt");
  refinement_Patterns[4] = new TPZRefPattern("/home/pos/cesar/RefPattern/Tetra_Unif.rpt");
  refinement_Patterns[5] = new TPZRefPattern("/home/pos/cesar/RefPattern/Piram_Unif.rpt");
  refinement_Patterns[6] = new TPZRefPattern("/home/pos/cesar/RefPattern/Prism_Unif.rpt");
  refinement_Patterns[7] = new TPZRefPattern("/home/pos/cesar/RefPattern/Hexa_Unif.rpt");
                             
  
  for (i=0;i<11;i++){
    int ncon = NodesPerEl[i];
    TPZVec <int> connect(ncon,0);
    for(j=0; j<ncon;j++){
      connect[j] = Connects[i][j];
    }
    //elvec[i] = GeoElement(gmesh,ncon,connect,1,i);
    if (ncon ==4 ) ncon = 7;
    elvec[i] = GeoElementRefPattern(gmesh,ncon,connect,1,i,refinement_Patterns);
  }
  //Generate neighborhod information
  gmesh->Print(cout);
  gmesh->BuildConnectivity();
  gmesh->Print(cout);
  //Create computational mesh
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  TPZMaterial *mat;
  mat = new TPZMaterialTest3D (1);
  cmesh->InsertMaterialObject(mat);

  cmesh->AutoBuild();
  return cmesh;
}

TPZCompMesh *TetraMesh(){
  REAL Coord [8][3] = {
    {0.,0.,0.},{1.,0.,0.},{1.,1.,0.},{0.,1.,0.},
    {0.,0.,1.},{1.,0.,1.},{1.,1.,1.},{0.,1.,1.}
  };

  int Connects [5][4] = {
    {0,1,3,4},
    {1,2,3,6},
    {5,6,4,1},
    {7,6,4,3},
    {1,3,4,6}
  };

  int i,j;
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  TPZGeoEl * elvec[5];
  TPZVec <REAL> coord (3,0.);
  int index;
  //Nodes initialization
  for(i = 0; i < 8; i++){
    for(j=0;j<3;j++){
      coord[j] = Coord[i][j];
    }
    index = gmesh->NodeVec().AllocateNewElement();
    gmesh->NodeVec()[index] = TPZGeoNode(i,coord,*gmesh);
  }

  TPZVec<TPZRefPattern *> refinement_Patterns(6,0);
//  refinement_Patterns.Resize(6);
  refinement_Patterns[0] = new TPZRefPattern("/home/pos/cesar/RefPattern/Tetra_Rib_Side_4.rpt");
  refinement_Patterns[1] = new TPZRefPattern("/home/pos/cesar/RefPattern/Tetra_Rib_Side_5.rpt");
  refinement_Patterns[2] = new TPZRefPattern("/home/pos/cesar/RefPattern/Tetra_Rib_Side_6.rpt");
  refinement_Patterns[3] = new TPZRefPattern("/home/pos/cesar/RefPattern/Tetra_Rib_Side_7.rpt");
  refinement_Patterns[4] = new TPZRefPattern("/home/pos/cesar/RefPattern/Tetra_Rib_Side_8.rpt");
  refinement_Patterns[5] = new TPZRefPattern("/home/pos/cesar/RefPattern/Tetra_Rib_Side_9.rpt");

  for (i=0;i<6;i++) gmesh->InsertRefPattern(refinement_Patterns[i]);

  for (i=0;i<5;i++){
    int ncon = 4;
    TPZVec <int> connect(ncon,0);
    for(j=0; j<ncon;j++){
      connect[j] = Connects[i][j];
    }
    elvec[i] = GeoElementRefPattern(gmesh,7,connect,1,i,refinement_Patterns);
  }
  //Generate neighborhod information
//  gmesh->Print(cout);
  gmesh->BuildConnectivity();
//  gmesh->Print(cout);
  //Create computational mesh
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  TPZMaterial *mat;
  mat = new TPZMaterialTest3D (1);
  cmesh->InsertMaterialObject(mat);

  cmesh->AutoBuild();
  return cmesh;
}

TPZGeoEl *GeoElement(TPZGeoMesh *gmesh, int ncon, TPZVec<int> &nodind, int material, int index){
  switch (ncon)  {
  case (2) : {
    return gmesh->CreateGeoElement(EOned,nodind,material,index);
    break;
  }
  case (3) : {
    return gmesh->CreateGeoElement(ETriangle,nodind,material,index);
    break;
  }
  case (4) : {
    return gmesh->CreateGeoElement(ETetraedro,nodind,material,index);
    break;
  }
  case (5) : {
    return gmesh->CreateGeoElement(EPiramide,nodind,material,index);
    break;
  }
  case (6) : {
    return gmesh->CreateGeoElement(EPrisma,nodind,material,index);
    break;
  }
  case (8) : {
    return gmesh->CreateGeoElement(ECube,nodind,material,index);
    break;
  }
  default:
    cout << "Erro..." << endl;
  }
  return 0;
}

TPZGeoEl *GeoElementRefPattern(TPZGeoMesh *gmesh,
                               int ncon,
                               TPZVec<int> &nodind,
                               int material,
                               int index,
                               TPZVec<TPZRefPattern *> refpattern){
  switch (ncon)  {
  case (1) : {
    TPZGeoElRefPattern <TPZShapePoint,TPZGeoPoint> *gel =
              new TPZGeoElRefPattern <TPZShapePoint,TPZGeoPoint> (index,nodind,1,*gmesh,refpattern[0]);
    return gel;
    break;
  }
  case (2) : {
    TPZGeoElRefPattern <TPZShapeLinear,TPZGeoLinear> *gel =
              new TPZGeoElRefPattern <TPZShapeLinear,TPZGeoLinear> (index,nodind,1,*gmesh,refpattern[1]);
    return gel;
    break;
  }
  case (3) : {
    TPZGeoElRefPattern <TPZShapeTriang,TPZGeoTriangle> *gel =
              new TPZGeoElRefPattern <TPZShapeTriang,TPZGeoTriangle> (index,nodind,1,*gmesh,refpattern[2]);
    return gel;
    break;
  }
  case (4) : {
    TPZGeoElRefPattern <TPZShapeQuad,TPZGeoQuad> *gel =
              new TPZGeoElRefPattern <TPZShapeQuad,TPZGeoQuad> (index,nodind,1,*gmesh,refpattern[3]);
    return gel;

   // return gmesh->CreateGeoElement(ETetraedro,nodind,material,index);
    break;
  }
  case (7) : { ///usado para o tetrahedro...
    TPZGeoElRefPattern <TPZShapeTetra,TPZGeoTetrahedra> *gel =
              new TPZGeoElRefPattern <TPZShapeTetra,TPZGeoTetrahedra> (index,nodind,1,*gmesh,refpattern[4]);
    return gel;
    break;
  }
  case (5) : {
    TPZGeoElRefPattern <TPZShapePiram,TPZGeoPyramid> *gel =
              new TPZGeoElRefPattern <TPZShapePiram,TPZGeoPyramid> (index,nodind,1,*gmesh,refpattern[5]);
    return gel;
    break;
  }
  case (6) : {
    TPZGeoElRefPattern <TPZShapePrism,TPZGeoPrism> *gel =
              new TPZGeoElRefPattern <TPZShapePrism,TPZGeoPrism> (index,nodind,1,*gmesh,refpattern[6]);
    return gel;
    break;
  }
  case (8) : {
    TPZGeoElRefPattern <TPZShapeCube,TPZGeoCube> *gel =
              new TPZGeoElRefPattern <TPZShapeCube,TPZGeoCube> (index,nodind,1,*gmesh,refpattern[7]);
    return gel;
    break;
  }
  default:
    cout << "Erro..." << endl;
  }
  return 0;
}


void WriteMesh(TPZGeoMesh *mesh,ofstream &arq){

  //mesh->Print(cout);

  arq << "object 1 class array type float rank 1 shape 3 items ";
  arq << mesh->NodeVec().NElements() << " data follows" << endl;
  int i;
  //Print Nodes
  for (i=0;i<mesh->NodeVec().NElements(); i++){
    TPZGeoNode *node = &mesh->NodeVec()[i];
    arq /*<< node->Id() << "\t"*/
      << node->Coord(0) << "\t"
      << node->Coord(1) << "\t"
      << node->Coord(2) << endl;
  }
  arq << "object 2 class array type integer rank 1 shape 8 items ";
  arq << mesh->ElementVec().NElements() << " data follows" << endl;

  TPZVec<int> elementtype(mesh->ElementVec().NElements(),0);
  for (i=0;i<mesh->ElementVec().NElements();i++){
    TPZGeoEl *el = mesh->ElementVec()[i];
    WriteElement (el,arq,elementtype);
  }

  arq << "attribute \"element type\" string \"cubes\"" << endl
    << "attribute \"ref\" string \"positions\"" << endl;

  arq << "object 3 class array type integer rank 0 items ";
  arq << mesh->ElementVec().NElements() << " data follows" << endl;
  for (i=0;i<mesh->ElementVec().NElements();i++){
    arq << elementtype[i] << endl;
  }

  arq << "attribute \"dep\" string \"connections\"" << endl;

  arq << "object 4 class field" << endl
    << "component \"positions\" value 1" << endl
    << "component \"connections\" value 2" << endl
    << "component \"data\" value 3" << endl;
}


void WriteElement (TPZGeoEl *el,ofstream &arq,TPZVec<int> &elementtype){
  int ncon = el->NNodes();
  elementtype[el->Id()] = ncon;
  switch (ncon)  {
  case (2) : {
    //rib
    int ni = el->NodePtr(0)->Id();
    int nf = el->NodePtr(1)->Id();
    arq << ni << "\t" << nf << "\t"
        << ni << "\t" << nf << "\t"
        << ni << "\t" << nf << "\t"
        << ni << "\t" << nf << endl;
    break;
  }
  case (3) : {
    //triangle
    int n0 = el->NodePtr(0)->Id();
    int n1 = el->NodePtr(1)->Id();
    int n2 = el->NodePtr(2)->Id();
    arq << n0 << "\t" << n1 << "\t"
        << n2 << "\t" << n2 << "\t"
        << n0 << "\t" << n1 << "\t"
        << n2 << "\t" << n2 << endl;
    break;
  }
  case (4) : {
    if (el->Dimension() == 2){
      //quad
      int n0 = el->NodePtr(0)->Id();
      int n1 = el->NodePtr(1)->Id();
      int n2 = el->NodePtr(3)->Id();
      int n3 = el->NodePtr(2)->Id();
      arq << n0 << "\t" << n1 << "\t"
          << n2 << "\t" << n3 << "\t"
          << n0 << "\t" << n1 << "\t"
          << n2 << "\t" << n3 << endl;
    }else{
      //tetrahedre
      int n0 = el->NodePtr(0)->Id();
      int n1 = el->NodePtr(1)->Id();
      int n2 = el->NodePtr(2)->Id();
      int n3 = el->NodePtr(3)->Id();
      arq << n0 << "\t" << n1 << "\t"
          << n2 << "\t" << n2 << "\t"
          << n3 << "\t" << n3 << "\t"
          << n3 << "\t" << n3 << endl;
    }
    break;
  }
  case (5) : {
    //pyramid
    int n0 = el->NodePtr(0)->Id();
    int n1 = el->NodePtr(1)->Id();
    int n2 = el->NodePtr(3)->Id();
    int n3 = el->NodePtr(2)->Id();
    int n4 = el->NodePtr(4)->Id();
    arq << n0 << "\t" << n1 << "\t"
        << n2 << "\t" << n3 << "\t"
        << n4 << "\t" << n4 << "\t"
        << n4 << "\t" << n4 << endl;
    break;
  }
  case (6) : {
    //pyramid
    int n0 = el->NodePtr(0)->Id();
    int n1 = el->NodePtr(1)->Id();
    int n2 = el->NodePtr(2)->Id();
    int n3 = el->NodePtr(3)->Id();
    int n4 = el->NodePtr(4)->Id();
    int n5 = el->NodePtr(5)->Id();
    arq << n0 << "\t" << n1 << "\t"
        << n2 << "\t" << n2 << "\t"
        << n3 << "\t" << n4 << "\t"
        << n5 << "\t" << n5 << endl;
    break;
  }
  case (8) : {
    int n0 = el->NodePtr(0)->Id();
    int n1 = el->NodePtr(1)->Id();
    int n2 = el->NodePtr(3)->Id();
    int n3 = el->NodePtr(2)->Id();
    int n4 = el->NodePtr(4)->Id();
    int n5 = el->NodePtr(5)->Id();
    int n6 = el->NodePtr(7)->Id();
    int n7 = el->NodePtr(6)->Id();
    arq << n0 << "\t" << n1 << "\t"
        << n2 << "\t" << n3 << "\t"
        << n4 << "\t" << n5 << "\t"
        << n6 << "\t" << n7 << endl;
    break;
  }
  default:
    cout << "Erro..." << endl;
  }
  return;
}
