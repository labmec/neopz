
/**
 * Gera arquivo de entrada para curso de Programacao de Alto Desempenho
 * 21/09/2004
 */

#include "pzvec.h"

#include "pzcmesh.h"

#include "pzdebug.h"
#include "pzcheckgeom.h"
//#include "pzerror.h"

#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgeoelside.h"

#include "pzintel.h"
#include "pzcompel.h"
#include "TPZCompElDisc.h"

#include "pzpoisson3d.h"

#include "pzfilebuffer.h"
#include "pzbfilestream.h"
#include "pzstring.h"
#include "pzelmat.h"

#include "TPZSpStructMatrix.h"
#include "pzysmp.h"
#include "pzstepsolver.h"

#include "pzanalysis.h"
#include "pzmaterial.h"
#include "pzbndcond.h"

#include <vector>
#include <map>

#define CONVERGENCE_OUTPUT

using namespace std;

//#define Comentarios
//#define BinFile
int gDivide[6] = {0,0,0,0,0,0};
//int gDivide[6] = {1,1,1,1,1,1};
TPZCompMesh * BuildMesh();
void Forcing1(TPZVec<REAL> &x, TPZVec<REAL> &disp);
void WriteMatrix();
TPZFYsmpMatrix * ReadMatrix (istream &arq, TPZVec<REAL> &B);
void PZProceed();


int main(){
  //PZProceed();
  ifstream arq ("teste.mat");
  TPZVec<REAL> vec;

  TPZFYsmpMatrix *matrix =  ReadMatrix(arq,vec);


  return 0;
}

void PZProceed(){
  TPZCompMesh *cMesh = BuildMesh();
  TPZAnalysis an (cMesh);
  TPZSpStructMatrix spStrMatrix (cMesh);
  an.SetStructuralMatrix(spStrMatrix);

  int solId = -1;
  cout << "Which solver?:\n\t( 0 ) - Jacobi\n\t( 1 ) - GausseSeidel\n\t ( 2 ) - SOR\n";
  cout << "\tConjugate Gradient (3)\n";
  cin >> solId;

  double tol = 1e-10;
  int maxIters = 100;
  double overRelax = 1.79;

  int neq = cMesh->NEquations();
  TPZFMatrix  cheia(neq,neq, 0.);
  TPZStepSolver step(&cheia);
  TPZStepSolver precond(step);

  TPZStepSolver iterSolver;
  switch (solId){
    case  ( 0 ) : {
      iterSolver.SetJacobi(maxIters,tol,0);
      break;
    }
    case  ( 1 ) : {
      iterSolver.SetSOR(maxIters,overRelax,tol,0);
      break;
    }
    case  ( 2 ) : {
      iterSolver.SetSOR(maxIters,1.,tol,0);
    }
    case  ( 3 ) : {
      iterSolver.SetCG(maxIters,precond,tol,0);
    }
    default :
      cout << "Error... undefined solver\n";
      exit (-1);
  }

  an.SetSolver(iterSolver);
  an.Run();

  //an.Rhs().Print();
  cout << Norm(an.Rhs());


}

TPZCompMesh * BuildMesh(){

  REAL co[9][2] = {{0.,0.},{0.,-1.},{1.,-1.},{1.,0.},{1.,1.},{0.,1.},{-1.,1.},{-1.,0.},{-1,-1}};
  int indices[4][4] = {{0,1,2,3},{0,3,4,5},{0,5,6,7},{0,7,8,1}};
  TPZGeoEl *elvec[4];
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  int nnode = 9;
  int nelem = 4;
  int nod;
  for(nod=0; nod<nnode; nod++) {
    int nodind = gmesh->NodeVec().AllocateNewElement();
    TPZManVector<REAL,2> coord(2);
    coord[0] = co[nod][0];
    coord[1] = co[nod][1];
    gmesh->NodeVec()[nodind].Initialize(nod,coord,*gmesh);
  }

  int el;
  for(el=0; el<nelem; el++) {
    TPZManVector<int,4> nodind(4);
    for(nod=0; nod<4; nod++) nodind[nod]=indices[el][nod];
    int index;
    elvec[el] = gmesh->CreateGeoElement(EQuadrilateral,nodind,1,index);
  }

  gmesh->BuildConnectivity();


  cout << "Iniciando divisoes de elementos" << endl;
  for(int i = 0; i < nelem; i++){
    TPZVec<TPZGeoEl *> children, netos, bisnetos, tata1, tata2, tata3;
    //    cout << "\ngDivide[0] = \n" << gDivide[0];
    if(gDivide[0] == 1) {
      elvec[i]->Divide(children);
      //      cout <<  "\n Primeira divisao \n" ;
      if (gDivide[1] == 1) {
        for(int j = 0; j < children.NElements(); j++) {
          children[j]->Divide(netos);
	  //          cout <<  "\n Segunda divisao \n" ;
          if(gDivide[2] == 1) {
            for(int k = 0; k < netos.NElements(); k++) {
              netos[k]->Divide(bisnetos);
	      //              cout <<  "\n Terceira divisao \n" ;
              if(gDivide[3] == 1) {
                for(int k = 0; k < bisnetos.NElements(); k++) {
                  bisnetos[k]->Divide(tata1);
		  //                  cout <<  "\n Quarta divisao \n" ;
                  if(gDivide[4] == 1) {
                    for(int k = 0; k < tata1.NElements(); k++) {
                      tata1[k]->Divide(tata2);
		      //                      cout <<  "\n Quinta divisao \n" ;
                      if(gDivide[5] == 1) {
                        for(int k = 0; k < tata2.NElements(); k++) tata2[k]->Divide(tata3);
			//                        cout <<  "\n Sexta divisao \n" ;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  cout << "Elementos divididos" << endl;



  // bc -1 -> Dirichlet homogeneo
  TPZGeoElBC gbc1(elvec[0],5,-2,*gmesh); // bottom
  TPZGeoElBC gbc2(elvec[0],6,-2,*gmesh); // right
  TPZGeoElBC gbc3(elvec[1],5,-1,*gmesh); // right
  TPZGeoElBC gbc4(elvec[1],6,-1,*gmesh); // top
  TPZGeoElBC gbc5(elvec[2],5,-2,*gmesh); // top
  TPZGeoElBC gbc6(elvec[2],6,-2,*gmesh); // left
  TPZGeoElBC gbc7(elvec[3],5,-1,*gmesh); // left
  TPZGeoElBC gbc8(elvec[3],6,-1,*gmesh); // bottom

  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(2);

  TPZMatPoisson3d *mat;
  mat = new TPZMatPoisson3d(1, 2);
  mat->SetForcingFunction(Forcing1);

  TPZManVector<REAL,2> convdir(2,0.);
  convdir[0] = 1.;
  convdir[1] = 1.;
  mat->SetParameters(1., 0., convdir);
  int nstate = 1;
  TPZFMatrix val1(nstate,nstate,0.),val2(nstate,1,0.);
  TPZBndCond *bc[2];

  val2(0,0) = 0.618; //Dirichlet nao homogeneo
  //val2.Zero(); //Dirichlet homogeneo
  bc[0] = mat->CreateBC(-1,0,val1,val2);
  //  bc[0]->SetForcingFunction(Dirichlet1);

  //val2.Zero(); //Dirichlet homogeneo
  val2(0,0) = 3.14159;//Dirichlet nao homogeneo
  bc[1] = mat->CreateBC(-2,0,val1,val2);
  //  bc[1]->SetForcingFunction(Dirichlet2);

  cmesh->InsertMaterialObject(mat);
  int i;
  for(i=0; i<2; i++) cmesh->InsertMaterialObject(bc[i]);


//   TPZGeoElement<TPZShapeCube,TPZGeoCube,TPZRefCube>::SetCreateFunction(TPZCompElDisc::CreateDisc);
//   TPZGeoElement<TPZShapeLinear,TPZGeoLinear,TPZRefLinear>::SetCreateFunction(TPZCompElDisc::CreateDisc);
//   TPZGeoElement<TPZShapeQuad,TPZGeoQuad,TPZRefQuad>::SetCreateFunction(TPZCompElDisc::CreateDisc);
//   TPZGeoElement<TPZShapeTriang,TPZGeoTriangle,TPZRefTriangle>::SetCreateFunction(TPZCompElDisc::CreateDisc);
//   TPZGeoElement<TPZShapePrism,TPZGeoPrism,TPZRefPrism>::SetCreateFunction(TPZCompElDisc::CreateDisc);
//   TPZGeoElement<TPZShapeTetra,TPZGeoTetrahedra,TPZRefTetrahedra>::SetCreateFunction(TPZCompElDisc::CreateDisc);
//   TPZGeoElement<TPZShapePiram,TPZGeoPyramid,TPZRefPyramid>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  //template class TPZGeoElement<TPZShapePoint,TPZGeoPoint,TPZRefPoint>;

  cmesh->AutoBuild();
  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();

  return cmesh;
}

void Forcing1(TPZVec<REAL> &x, TPZVec<REAL> &disp) {
  disp[0] = - exp(0.75 * (x[0] + x[1])) * (8. * (1. - x[1] * x[1]) + 12. * x[0] * (1. - x[1] * x[1]) -4.5 * (1. - x[0] * x[0])* (1. - x[1] * x[1] )+
					  8. * (1. - x[0] * x[0]) + 12. * x[1] * (1. - x[0] * x[0]));
}


void WriteMatrix(){
  int aux;
  cout << "Binary Files: yes = 1 ; false = 0" << endl;
  cin >> aux ;

  TPZCompEl::SetgOrder(1);

#ifdef BinFile
  TPZBFileStream Out;
#endif
#ifndef BinFile
  TPZFileStream Out;
#endif

  Out.OpenWrite("/compile/cesar/NeoPZ/Projects/AltoDes/CursoAltoDes.in");

  cout << "Gerando a malha" << endl;
  TPZCompMesh * cmesh = BuildMesh();
  cout << "Malha pronta" << endl;

  int neq = cmesh->NEquations();

#ifdef Comentarios
  TPZString messg;
  messg = "Arquivo de entrada do curso de alto desempenho";
  Out.WriteSameLine(&messg[0], messg.Length() );

  messg = "Numero de equacoes:";
  Out.WriteSameLine(&messg[0], messg.Length() );
#endif

  Out.Write(&neq, 1);

  int nel = 0;
  for(int i = 0; i < cmesh->NElements(); i++){
    TPZCompEl* cel = cmesh->ElementVec()[i];
    if (!cel) continue;
    nel++;
  }

#ifdef Comentarios
  messg = "Numero de elementos:";
  Out.WriteSameLine(&messg[0], messg.Length() );
#endif
  Out.Write(&nel, 1);

#ifdef Comentarios
  messg = "Eh simetrica:";
  Out.WriteSameLine(&messg[0], messg.Length() );
#endif

//acho que eh simetrico.
#warning Verificar se eh simetrico
  int issimetric = 1;
  Out.Write(&issimetric, 1);

  TPZElementMatrix ek, ef;

//nos dos elementos
//incidencia
  cout << "Escrevendo incidencia" << endl;

#ifdef Comentarios
  messg = "nnos";
  Out.WriteSameLine(&messg[0], messg.Length() );

  messg = "no_0, no_1 ... no_n";
  Out.WriteSameLine(&messg[0], messg.Length() );
#endif

  for(int i = 0; i < cmesh->ElementVec().NElements(); i++){
    TPZCompEl * cel = cmesh->ElementVec()[i];
    if (!cel) continue;
    cel->CalcStiff (ek, ef);
    //numero de nos
    int cel_neq = ek.fMat.Rows();
    Out.Write( &cel_neq, 1 );
    int nconnects = cel->NConnects();
    int destindex = 0;
    TPZManVector<int> destinationindex(cel_neq);
    for(int in = 0; in < nconnects; in++){
      int connect_index = ek.ConnectIndex(in);
      TPZConnect &Connect = cmesh->ConnectVec()[connect_index];
      int blocknumber = Connect.SequenceNumber();
      int firsteq = cmesh->Block().Position(blocknumber);
      int ndof = cmesh->Block().Size(blocknumber);
      for(int idof = 0; idof < ndof; idof++)
        destinationindex[destindex++] = firsteq+idof;
    }
      //nos:
      Out.Write(&destinationindex[0], cel_neq);
      //cout << i << endl;
  }

  // matriz de rigidez e vetor de carga
  cout << "Escrevendo matrizes de rigidez e vetores de carga" << endl;
  for(int i = 0; i < cmesh->ElementVec().NElements(); i++){
    TPZCompEl * cel = cmesh->ElementVec()[i];
    //só vou pegar elementos bi-dimensionais. As condicoes de contorno eu deixo pra la,
    // porque a matriz de rigidez dos elementos de contorno é diferente da matriz padrao.
    if (!cel) continue;
    cel->CalcStiff (ek, ef);

#ifdef Comentarios
    messg = "Matriz de rigidez:";
    Out.WriteSameLine(&messg[0], messg.Length() );
#endif
    int cel_neq = ek.fMat.Rows();

#ifdef Comentarios
    messg = "Matriz - elementos em coluna:";
    Out.WriteSameLine(&messg[0], messg.Length() );
#endif

    Out.Write( &( ek.fMat(0,0) ), cel_neq * cel_neq );

#ifdef Comentarios
    messg = "Vetor de carga";
    Out.WriteSameLine(&messg[0], messg.Length() );
#endif

    Out.Write( &( ef.fMat(0,0) ), cel_neq );

    cout << i << endl;
  }
}



TPZFYsmpMatrix * ReadMatrix (istream &arq, TPZVec<REAL> &B){

  //Number of equations
  int size;
  arq >> size;

  //multiset structure...
  //              row               col     val
  //std::multimap < int , std::pair < int , double > > values;
  std::vector< std::map <int,double> > sparsemat (size);
  B.Resize(size);
  //std::vector<double> B(size);

  //Number of elements
  int nel;
  arq >> nel;
  std::vector< std::vector<int> > nodes (nel);

  //Is symetric?
  int issym;
  arq >> issym;

  //Read the elements connectivities
  register int e,n,auxn;
  int nnodes;
  for (e=0;e<nel;e++){
    arq >> nnodes;
    nodes[e].resize(nnodes);
    for (n=0;n<nnodes;n++){
      arq >> (nodes[e])[n];
    }
    for (n=0;n<nnodes;n++){
      //each node represents one line...
      std::map < int,double > &linemap = sparsemat[ (nodes[e])[n] ];
      for (auxn=0;auxn<nnodes;auxn++){
        linemap [ (nodes[e])[auxn] ] = 0.;
      }
    }
  }//the vector of maps is created and all values is zero!

  double value;
  for (e=0;e<nel;e++){
    nnodes = nodes[e].size();
    for (n=0;n<nnodes;n++){
      std::map < int,double > &linemap = sparsemat[ (nodes[e])[n] ];
      for (auxn=0;auxn<nnodes;auxn++){
        arq >> value;
        linemap [ (nodes[e])[auxn] ] += value;
      }
    }
    for (n=0;n<nnodes;n++){
      arq >> value;
      B[(nodes[e])[n]] += value;
    }
  }

  //View the sparse matrix
  int i,j;
  std::map<int,double>::iterator it;
  for (i=0;i<sparsemat.size();i++){
    std::map<int,double> &linemap = sparsemat[i];
    for (it=linemap.begin();it!=linemap.end();it++){
      std::cout << " \t\t( " << i << " , " << (*it).first << " )  = " << (*it).second;
    }
    std::cout << " \t\t f = " << B[i] << std::endl;
  }



  int count_row = 0;
  int nz = 0;

  int m = size;
  n = size;
  std::vector<int> *rowptr = new std::vector<int>(size+1);
  std::vector<int> *colptr = new std::vector<int>(nz);
  std::vector<double> *val = new std::vector<double>(nz);

  int ncols =  0;
  int countcols = 0;
  int line;
  int counter = 0;

  for (i=0;i<size;i++){
    (*rowptr)[i] = countcols;
    std::map<int,double> &linemap = sparsemat[i];
    ncols = linemap.size();
    countcols += ncols;
    it = linemap.begin();

    for (j=0;j<ncols;j++,it++){
      (*colptr)[counter] = (*it).first;
      (*val)[counter] = (*it).second;
/*      std::cout << "col " << col << " value " << v << "colptr " << (*colptr)[counter] << std::endl;*/
      counter++;
    }
  }
  (*rowptr)[size] = nz;

//   double val[] = {1.0, 4.0, 5.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0};
//   int colind[] = {0,   0,   1,   1,   2,   2,   3,    0,    4   };
//   int rowptr[] = {0,   1,        3,        5,         7,           9};

/*  TNT::Sparse_Matrix_CompRow<double> *A =
  new TNT::Sparse_Matrix_CompRow<double>(m,n,nz,&(*val->begin()), &(*rowptr->begin()),&(*colptr->begin()));*/

  TPZFYsmpMatrix *A = new TPZFYsmpMatrix (size,size);
  A->SetData(&(*rowptr->begin()),&(*colptr->begin()),&(*val->begin()));
  return A;
}
/*
void GenerateVector(int size,double valmin, double valmax,TNT::Array1D<double> &vec)
{
  int i;
  //vec.resize(size);
  for (i=0;i<size;i++)
  {
    vec[i] = GenerateValue(valmin,valmax);
  }
}*/
