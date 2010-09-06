#include <sstream>

#include "pzvec.h"
#include "pzcmesh.h"
#include "pzdebug.h"
#include "pzcheckgeom.h"

#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgeoelside.h"

#include "pzgeoel.h"
#include "pzcompel.h"
#include "TPZCompElDisc.h"
#include "pzmatrix.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZParFrontMatrix.h"
#include "TPZFrontNonSym.h"
#include "pzbdstrmatrix.h"
#include "pzblockdiag.h"
#include "TPZSpStructMatrix.h"
#include "TPZCopySolve.h"
#include "TPZStackEqnStorage.h"

#include "pzbstrmatrix.h"
#include "pzstepsolver.h"

#include "pzbndcond.h"
#include "pzpoisson3d.h"

#include "pzvisualmatrix.h"

#include "TPZGeoElement.h"
#include "pzgeoel.h"

#include "TPZGeoCube.h"
#include "pzshapecube.h"
#include "TPZRefCube.h"
#include "pzshapelinear.h"
#include "TPZGeoLinear.h"
#include "TPZRefLinear.h"
#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "pzshapetriang.h"
#include "pzreftriangle.h"
#include "pzgeotriangle.h"
#include "pzshapeprism.h"
#include "pzrefprism.h"
#include "pzgeoprism.h"
#include "pzshapetetra.h"
#include "pzreftetrahedra.h"
#include "pzgeotetrahedra.h"
#include "pzshapepiram.h"
#include "pzrefpyram.h"
#include "pzgeopyramid.h"
#include "pzrefpoint.h"
#include "pzgeopoint.h"
#include "pzshapepoint.h"
#include "pzskylstrmatrix.h"

#include <time.h>
#include <stdio.h>
#include "pzl2projection.h"
#include "tpzgeoelmapped.h"

/*******************************************TPZGeoMesh
 *   Made by Caju                          *
 *   LabMeC - 2008                         *
 *******************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include <cstdlib>

#include "tpzmathtools.h"
#include "pzgeotriangle.h"
#include "tpzcurvedtriangle.h"
#include "tpzquadratictrig.h"
//#include "tpzquadratictetra.h"
#include "tpzquadraticquad.h"
#include "tpzarc3d.h"
#include "tpzgeoelrefpattern.h"
#include "tpzgeoelrefpattern.h.h"
#include "pzgeoelside.h"
#include "tpzgeoblend.h"
#include "tpzchangeel.h"
#include "pzgeoprism.h"
#include "pzgeopyramid.h"
#include "TPZGeoCube.h"
#include "pzcompel.h"
//#include "tpzellipse.h"
#include "tpzblendnaca.h"
#include "pzelasAXImat.h"
#include "pzmaterialdata.h"
#include "pzmganalysis.h"
#include "TPZNLMultGridAnalysis.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("main"));
#endif

#include <sstream>
using namespace std;
using namespace pzgeom;
using namespace pzshape;
using namespace pzrefine;
struct MatIds;

void AuxEspiral(TPZVec <REAL> &vecIn, TPZVec <REAL> &vecOut, TPZVec <REAL> &C, double R);
TPZGeoMesh * CxEspiral2D(double Bb, double Hr, double Bt, double Hl, double Cx, double Cy,
                         double R, double h, double b, double Dx, double Dy, double e1,
                         double e2, double e3, double X1, double X2, double X3, double X4,
                         double X5, double h1, double h2, double Py, double es, MatIds &mat);
TPZCompMesh * SquareMesh();
TPZCompMesh * SpiralMesh(MatIds &MatId);
void ApplyDisplacements(TPZCompMesh & cmesh);

void WriteMesh(TPZGeoMesh *mesh,std::ofstream &arq);
void WriteElement(TPZGeoEl *el,int elindex, std::ofstream &arq,TPZVec<int> &elementtype);

void IntegralForces(TPZCompMesh * cmesh, int MatId, TPZVec<REAL> &forces);

REAL CalcArea(TPZCompMesh *mesh);

void MathOutUP(TPZBlendNACA * naca);
void MathOutDW(TPZBlendNACA * naca);


//   int ElMat = 1;
//   int PreDistribMat = 1;
//   int PreDistribAnelMat = 1;

//   int ElMat = 1;
//   int PreDistribMat = 2;
//   int PreDistribAnelMat = 3;

//   int Nonemat = -8;
//   int Dotmat = -9;
//   int PNmat = -8;
//   int PD2mat = PNmat;
//   int TTmat = Nonemat;
//   int PDmat = Nonemat;
//   int FFmat = Nonemat;
//   int FEmat = Nonemat;
//   int PGmat = Nonemat;
//   int Basemat = Nonemat;
//   int RRBc = Nonemat;

//   int TTmat = -2;
//   int PDmat = -3;
//   int FFmat = -4;
//   int FEmat = -5;
//   int PGmat = -6;
//   int Basemat = -1;
//   int RRBc = -7;
//   int Dotmat = -9;

struct MatIds
 {
   enum {ElMat =1, PreDistribMat = 2, PreDistribAnelMat = 3, Basemat = -1, TTmat = -2, PDmat = -3, FFmat = -4, FEmat = -5, PGmat = -6, RRBc = -7, PNmat = -8, Dotmat = -9, OutBound = -10, PD2mat = -11, Nonemat = -100};
  std::map<int, int> fMat;
  REAL f_rho;
  MatIds()
  {
    SetFullLoad();
  }
  void SetFullLoad()
  {
    fMat[ElMat] = ElMat;
    fMat[PreDistribMat] = PreDistribMat;
    fMat[PreDistribAnelMat] = PreDistribAnelMat;
    fMat[TTmat] = TTmat;
    fMat[PDmat] = PDmat;
    fMat[PD2mat] = PD2mat;
    fMat[FFmat] = FFmat;
    fMat[FEmat] = FEmat;
    fMat[PGmat] = PGmat;
    fMat[Basemat] = Basemat;
    fMat[RRBc] = RRBc;
    fMat[PNmat] = PNmat;
    fMat[OutBound] = Nonemat;
    fMat[Dotmat] = Dotmat;
    f_rho = -0.025;
  }
  void SetInternalPressure()
  {
    fMat[ElMat] = ElMat;
    fMat[PreDistribMat] = PreDistribMat;
    fMat[PreDistribAnelMat] = PreDistribAnelMat;
    fMat[TTmat] = Nonemat;
    fMat[PDmat] = Nonemat;
    fMat[PD2mat] = PNmat;
    fMat[FFmat] = Nonemat;
    fMat[FEmat] = Nonemat;
    fMat[PGmat] = Nonemat;
    fMat[Basemat] = Nonemat;
    fMat[RRBc] = Nonemat;
    fMat[PNmat] = PNmat;
    fMat[OutBound] = Nonemat;
    fMat[Dotmat] = Dotmat;
    f_rho = 0.;
  }
  void SetHydrostatic()
  {
    fMat[ElMat] = ElMat;
    fMat[PreDistribMat] = ElMat;
    fMat[PreDistribAnelMat] = ElMat;
    fMat[TTmat] = PNmat;
    fMat[PDmat] = Nonemat;
    fMat[PD2mat] = PNmat;
    fMat[FFmat] = PNmat;
    fMat[FEmat] = PNmat;
    fMat[PGmat] = PNmat;
    fMat[Basemat] = PNmat;
    fMat[RRBc] = PNmat;
    fMat[PNmat] = PNmat;
    fMat[OutBound] = PNmat;
    fMat[Dotmat] = Dotmat;
    f_rho = 0.;
  }
 };
/**
 *
 * @param argc
 * @param argv[]
 * @return
 */
int main(int argc, char *argv[])
{
#ifdef LOG4CXX
  InitializePZLOG("log4cxx.cfg");
#endif

//     TPZCompMesh * cmesh = SquareMesh();
    MatIds start;
//    start.SetInternalPressure();
    TPZCompMesh * cmesh = SpiralMesh(start);

//     TPZCompMesh *refmesh = TPZNonLinMultGridAnalysis::UniformlyRefineMesh (cmesh,1,4);
//     delete cmesh;
//     cmesh = refmesh;
//     cmesh->LoadReferences();

//////////////////////////////////////////////////////////////////////
#ifdef LOG4CXX
{//Verificando se foi esquecido algum elemento (linear) de contorno
    std::stringstream sout;
    sout << endl;

    TPZGeoMesh * gmesh = cmesh->Reference();
    int nelements = gmesh->NElements();

    for(int el = 0; el < nelements; el++)
    {
          TPZGeoEl * El = gmesh->ElementVec()[el];

          for(int side = 0; side < El->NSides(); side++)
          {
              TPZGeoElSide SideEl(El,side);
              if(SideEl.Dimension() == 1)
              {
                  TPZStack<TPZGeoElSide> allneigh;
                  SideEl.AllNeighbours(allneigh);

                  if(allneigh.NElements() != 1)
                  {
                          sout << endl << "Elemento: " << el << " dim " << El->Dimension() << " Lado: " << side << " -> NOT Ok!" <<
                              " numero de elementos vizinhos " << allneigh.NElements() <<
                              " sidenodes " << SideEl.SideNodeIndex(0) << " " << SideEl.SideNodeIndex(1) <<   endl;
                  }
              }
          }
    }
    LOGPZ_DEBUG(logger,sout.str());
}
#endif

	TPZCompMesh *fine = TPZMGAnalysis::UniformlyRefineMesh (cmesh,false);
	delete cmesh;
	cmesh = fine;
	cmesh->LoadReferences();
//////////////////////////////////////////////////////////////////////

    ///Inserindo os deslocamentos "na unha!"
//     ApplyDisplacements(*cmesh);

#ifdef LOG4CXX
{
  std::stringstream sout;
  cmesh->Reference()->Print(sout);
  LOGPZ_DEBUG(logger,sout.str());
}


    {
        std::stringstream sout;
        cmesh->Print(sout);
        LOGPZ_DEBUG(logger,sout.str());
    }
#endif



    TPZAnalysis an(cmesh);
    TPZSkylineStructMatrix full(cmesh);
    an.SetStructuralMatrix(full);
    TPZStepSolver step;
    step.SetDirect(ECholesky);
    an.SetSolver(step);
    an.Run();

    TPZVec<std::string> scalnames(5), vecnames(3);
    vecnames[0] = "Eigenvector1";
    vecnames[1] = "Eigenvector2";
    vecnames[2] = "Eigenvector3";
    scalnames[0] = "Sigmarr";
    scalnames[1] = "Sigmazz";
    scalnames[2] = "Sigmatt";
    scalnames[3] = "Taurz";
    scalnames[4] = "MohrCoulomb";
    std::string plotfile("saida.dx");
    an.DefineGraphMesh(2,scalnames,vecnames,plotfile);
    an.PostProcess(3);

    std::ofstream out("malha.txt");
    an.Print("nothing",out);

#ifdef LOG4CXX
{
  std::stringstream sout;
  TPZVec<REAL> forces;
  IntegralForces(cmesh, MatIds::Basemat, forces);
  sout << "\n===========================================================\n" << forces << " CalcArea " << CalcArea(cmesh) <<  "\n===========================================================\n";
  LOGPZ_DEBUG(logger,sout.str());
}
#endif

///End Malha Axi-simetrica

    return EXIT_SUCCESS;
}

void AuxEspiral(TPZVec <REAL> &vecIn, TPZVec <REAL> &vecOut, TPZVec <REAL> &C, double R)
{
    double norm = sqrt(vecIn[0]*vecIn[0] + vecIn[1]*vecIn[1] + vecIn[2]*vecIn[2]);
    for(int i = 0; i < 3; i++)
    {
        vecOut[i] = vecIn[i]/norm*R + C[i];
    }
}

TPZGeoMesh * CxEspiral2D(double Bb, double Hr, double Bt, double Hl, double Cx, double Cy,
                         double R, double h, double b, double Dx, double Dy, double e1,
                         double e2, double e3, double X1, double X2, double X3, double X4,
                         double X5, double h1, double h2, double Py, double es, MatIds &mat)
{
    TPZGeoMesh * Mesh = new TPZGeoMesh;
    //Mesh->InitializeRefPatterns();
    int Qnodes = 62;
    TPZVec < TPZVec <REAL> > NodeCoord(Qnodes);
    for(int i = 0; i < Qnodes; i++) NodeCoord[i].Resize(3);

    NodeCoord[0][0] = 0.;           NodeCoord[1][0] = Cx;       NodeCoord[2][0] = Bb;   NodeCoord[3][0] = Bb + 2.*Dx/3.;
    NodeCoord[0][1] = 0.;           NodeCoord[1][1] = 0.;       NodeCoord[2][1] = 0.;   NodeCoord[3][1] = 2.*Dy/3.;
    NodeCoord[0][2] = 0.;           NodeCoord[1][2] = 0.;       NodeCoord[2][2] = 0.;   NodeCoord[3][2] = 0.;

    NodeCoord[4][0] = Bb+11.*Dx/12.;  NodeCoord[5][0] = Bb + Dx;  NodeCoord[6][0] = Bb + Dx - e1;
    NodeCoord[4][1] = 11.*Dy/12.;     NodeCoord[5][1] = Dy;       NodeCoord[6][1] = Dy;
    NodeCoord[4][2] = 0.;           NodeCoord[5][2] = 0.;       NodeCoord[6][2] = 0.;

    NodeCoord[7][0] = Bt + e2;      NodeCoord[8][0] = Bt;       NodeCoord[9][0] = Bt;               NodeCoord[10][0] = Bt;
    NodeCoord[7][1] = Py;           NodeCoord[8][1] = Py;       NodeCoord[9][1] = Cy - h/2.;        NodeCoord[10][1] = Cy + h/2.;
    NodeCoord[7][2] = 0.;           NodeCoord[8][2] = 0.;       NodeCoord[9][2] = 0.;               NodeCoord[10][2] = 0.;

    NodeCoord[11][0] = Bt;          NodeCoord[12][0] = Bt;      NodeCoord[13][0] = X1+X2+X3+X4+X5;  NodeCoord[14][0] = X1 + X2 + X3 + X4;
    NodeCoord[11][1] = Hl;          NodeCoord[12][1] = Hr;      NodeCoord[13][1] = Hr;              NodeCoord[14][1] = Hl + h2;
    NodeCoord[11][2] = 0.;          NodeCoord[12][2] = 0.;      NodeCoord[13][2] = 0.;              NodeCoord[14][2] = 0.;

    NodeCoord[15][0] = X1+X2+X3;    NodeCoord[16][0] = X1+X2+X3;NodeCoord[17][0] = X1+X2;           NodeCoord[18][0] = X1 + X2;
    NodeCoord[15][1] = Hl + h2;     NodeCoord[16][1] = Hl;      NodeCoord[17][1] = Hl;              NodeCoord[18][1] = Hl + h1;
    NodeCoord[15][2] = 0.;          NodeCoord[16][2] = 0.;      NodeCoord[17][2] = 0.;              NodeCoord[18][2] = 0.;

    NodeCoord[19][0] = X1;          NodeCoord[20][0] = X1;      NodeCoord[21][0] = 0.;              NodeCoord[22][0] = 0.;
    NodeCoord[19][1] = Hl + h1;     NodeCoord[20][1] = Hl;      NodeCoord[21][1] = Hl;              NodeCoord[22][1] = Hl/2.;
    NodeCoord[19][2] = 0.;          NodeCoord[20][2] = 0.;      NodeCoord[21][2] = 0.;              NodeCoord[22][2] = 0.;

    NodeCoord[23][0] = Bt + e3;     NodeCoord[24][0] = Bt+e3-b; NodeCoord[25][0] = Bt + e3 - b;     NodeCoord[26][0] = Bt + e3;
    NodeCoord[23][1] = Cy - h/2.;   NodeCoord[24][1] = Cy-h/2.; NodeCoord[25][1] = Cy + h/2.;       NodeCoord[26][1] = Cy + h/2.;
    NodeCoord[23][2] = 0.;          NodeCoord[24][2] = 0.;      NodeCoord[25][2] = 0.;              NodeCoord[26][2] = 0.;

    NodeCoord[27][0] = X1+X2+X3+X4; NodeCoord[28][0] = X1+X2+X3+X4+X5;
    NodeCoord[27][1] = Hl;          NodeCoord[28][1] = Hl;
    NodeCoord[27][2] = 0.;          NodeCoord[28][2] = 0.;

    TPZVec <REAL> C(3,0.); C[0] = Cx; C[1] = Cy;
    TPZVec <REAL> In(3,0.); TPZVec <REAL> Out(3,0.);

    NodeCoord[29][0] = (2.*Cx + sqrt(-h*h + 4.*R*R))/2.;
    NodeCoord[29][1] = Cy + h/2.;
    NodeCoord[29][2] = 0.;

    for(int i = 0; i < 3; i++) In[i] = (NodeCoord[10][i] + NodeCoord[11][i])/2. - C[i];
    AuxEspiral(In, Out, C, R);
    NodeCoord[30][0] = Out[0];
    NodeCoord[30][1] = Out[1];
    NodeCoord[30][2] = Out[2];

    for(int i = 0; i < 3; i++) In[i] = NodeCoord[11][i] - C[i];
    AuxEspiral(In, Out, C, R);
    NodeCoord[31][0] = Out[0];
    NodeCoord[31][1] = Out[1];
    NodeCoord[31][2] = Out[2];

    for(int i = 0; i < 3; i++) In[i] = (NodeCoord[11][i] + NodeCoord[28][i])/2. - C[i];
    AuxEspiral(In, Out, C, R);
    NodeCoord[32][0] = Out[0];
    NodeCoord[32][1] = Out[1];
    NodeCoord[32][2] = Out[2];

    for(int i = 0; i < 3; i++) In[i] = NodeCoord[28][i] - C[i];
    AuxEspiral(In, Out, C, R);
    NodeCoord[33][0] = Out[0];
    NodeCoord[33][1] = Out[1];
    NodeCoord[33][2] = Out[2];

    for(int i = 0; i < 3; i++) In[i] = (NodeCoord[27][i] + NodeCoord[28][i])/2. - C[i];
    AuxEspiral(In, Out, C, R);
    NodeCoord[34][0] = Out[0];
    NodeCoord[34][1] = Out[1];
    NodeCoord[34][2] = Out[2];

    for(int i = 0; i < 3; i++) In[i] = NodeCoord[27][i] - C[i];
    AuxEspiral(In, Out, C, R);
    NodeCoord[35][0] = Out[0];
    NodeCoord[35][1] = Out[1];
    NodeCoord[35][2] = Out[2];

    for(int i = 0; i < 3; i++) In[i] = (NodeCoord[27][i] + NodeCoord[16][i])/2. - C[i];
    AuxEspiral(In, Out, C, R);
    NodeCoord[36][0] = Out[0];
    NodeCoord[36][1] = Out[1];
    NodeCoord[36][2] = Out[2];

    for(int i = 0; i < 3; i++) In[i] = NodeCoord[16][i] - C[i];
    AuxEspiral(In, Out, C, R);
    NodeCoord[37][0] = Out[0];
    NodeCoord[37][1] = Out[1];
    NodeCoord[37][2] = Out[2];

    for(int i = 0; i < 3; i++) In[i] = (NodeCoord[17][i] + NodeCoord[16][i])/2. - C[i];
    AuxEspiral(In, Out, C, R);
    NodeCoord[38][0] = Out[0];
    NodeCoord[38][1] = Out[1];
    NodeCoord[38][2] = Out[2];

    for(int i = 0; i < 3; i++) In[i] = NodeCoord[17][i] - C[i];
    AuxEspiral(In, Out, C, R);
    NodeCoord[39][0] = Out[0];
    NodeCoord[39][1] = Out[1];
    NodeCoord[39][2] = Out[2];

    for(int i = 0; i < 3; i++) In[i] = (NodeCoord[17][i] + NodeCoord[20][i])/2. - C[i];
    AuxEspiral(In, Out, C, R);
    NodeCoord[40][0] = Out[0];
    NodeCoord[40][1] = Out[1];
    NodeCoord[40][2] = Out[2];

    for(int i = 0; i < 3; i++) In[i] = NodeCoord[20][i] - C[i];
    AuxEspiral(In, Out, C, R);
    NodeCoord[41][0] = Out[0];
    NodeCoord[41][1] = Out[1];
    NodeCoord[41][2] = Out[2];

    for(int i = 0; i < 3; i++) In[i] = (NodeCoord[20][i] + NodeCoord[21][i])/2. - C[i];
    AuxEspiral(In, Out, C, R);
    NodeCoord[42][0] = Out[0];
    NodeCoord[42][1] = Out[1];
    NodeCoord[42][2] = Out[2];

    for(int i = 0; i < 3; i++) In[i] = NodeCoord[21][i] - C[i];
    AuxEspiral(In, Out, C, R);
    NodeCoord[43][0] = Out[0];
    NodeCoord[43][1] = Out[1];
    NodeCoord[43][2] = Out[2];

    for(int i = 0; i < 3; i++) In[i] = (NodeCoord[21][i] + NodeCoord[22][i])/2. - C[i];
    AuxEspiral(In, Out, C, R);
    NodeCoord[44][0] = Out[0];
    NodeCoord[44][1] = Out[1];
    NodeCoord[44][2] = Out[2];

    for(int i = 0; i < 3; i++) In[i] = NodeCoord[22][i] - C[i];
    AuxEspiral(In, Out, C, R);
    NodeCoord[45][0] = Out[0];
    NodeCoord[45][1] = Out[1];
    NodeCoord[45][2] = Out[2];

    for(int i = 0; i < 3; i++) In[i] = (NodeCoord[22][i] + NodeCoord[0][i])/2. - C[i];
    AuxEspiral(In, Out, C, R);
    NodeCoord[46][0] = Out[0];
    NodeCoord[46][1] = Out[1];
    NodeCoord[46][2] = Out[2];

    for(int i = 0; i < 3; i++) In[i] = NodeCoord[0][i] - C[i];
    AuxEspiral(In, Out, C, R);
    NodeCoord[47][0] = Out[0];
    NodeCoord[47][1] = Out[1];
    NodeCoord[47][2] = Out[2];

    for(int i = 0; i < 3; i++) In[i] = (NodeCoord[0][i] + NodeCoord[1][i])/2. - C[i];
    AuxEspiral(In, Out, C, R);
    NodeCoord[48][0] = Out[0];
    NodeCoord[48][1] = Out[1];
    NodeCoord[48][2] = Out[2];

    for(int i = 0; i < 3; i++) In[i] = NodeCoord[1][i] - C[i];
    AuxEspiral(In, Out, C, R);
    NodeCoord[49][0] = Out[0];
    NodeCoord[49][1] = Out[1];
    NodeCoord[49][2] = Out[2];

    for(int i = 0; i < 3; i++) In[i] = (NodeCoord[1][i] + NodeCoord[2][i])/2. - C[i];
    AuxEspiral(In, Out, C, R);
    NodeCoord[50][0] = Out[0];
    NodeCoord[50][1] = Out[1];
    NodeCoord[50][2] = Out[2];

    for(int i = 0; i < 3; i++) In[i] = NodeCoord[2][i] - C[i];
    AuxEspiral(In, Out, C, R);
    NodeCoord[51][0] = Out[0];
    NodeCoord[51][1] = Out[1];
    NodeCoord[51][2] = Out[2];

    for(int i = 0; i < 3; i++) In[i] = NodeCoord[8][i] - C[i];
    AuxEspiral(In, Out, C, R);
    NodeCoord[52][0] = Out[0];
    NodeCoord[52][1] = Out[1];
    NodeCoord[52][2] = Out[2];

    NodeCoord[53][0] = (2.*Cx + sqrt(-h*h + 4.*R*R))/2.;    NodeCoord[54][0] = Bt + e3 - b;
    NodeCoord[53][1] = Cy - h/2.;                           NodeCoord[54][1] = Cy - h/2. + es;
    NodeCoord[53][2] = 0.;                                  NodeCoord[54][2] = 0.;

    NodeCoord[55][0] = (2.*Cx + sqrt(-h*h + 4.*R*R))/2.;    NodeCoord[56][0] = Bt;
    NodeCoord[55][1] = Cy - h/2. + es;                      NodeCoord[56][1] = Cy - h/2. + es;
    NodeCoord[55][2] = 0.;                                  NodeCoord[56][2] = 0.;

    NodeCoord[57][0] = Bt + e3;                             NodeCoord[58][0] = Bt + e3 - b;
    NodeCoord[57][1] = Cy - h/2. + es;                      NodeCoord[58][1] = Cy + h/2. - es;
    NodeCoord[57][2] = 0.;                                  NodeCoord[58][2] = 0.;

    NodeCoord[59][0] = (2.*Cx + sqrt(-h*h + 4.*R*R))/2.;    NodeCoord[60][0] = Bt;
    NodeCoord[59][1] = Cy + h/2. - es;                      NodeCoord[60][1] = Cy + h/2. - es;
    NodeCoord[59][2] = 0.;                                  NodeCoord[60][2] = 0.;

    NodeCoord[61][0] = Bt + e3;                             NodeCoord[58][0] = Bt + e3 - b;
    NodeCoord[61][1] = Cy + h/2. - es;                      NodeCoord[58][1] = Cy + h/2. - es;
    NodeCoord[61][2] = 0.;
    //Ufa!!!!!!!!!!!!!!

    Mesh->NodeVec().Resize(Qnodes);
    TPZVec <TPZGeoNode> Node(Qnodes);
    for(int n = 0; n < Qnodes; n++)
    {
        Node[n].SetNodeId(n);
        Node[n].SetCoord(&NodeCoord[n][0]);
        Mesh->NodeVec()[n] = Node[n];
    }

    TPZVec <int> Topol(4);
    /// *********** CHOOSE MATERIAL FOR EACH ELEMENT ***********

    int id = 0;
    ///Quadrilaterals
    {
        //El0
        Topol[0] = 22; Topol[1] = 0; Topol[2] = 47; Topol[3] = 45;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,mat.fMat[MatIds::ElMat],*Mesh);
        id++;

        //El1
        Topol[0] = 0; Topol[1] = 1; Topol[2] = 49; Topol[3] = 47;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,mat.fMat[MatIds::ElMat],*Mesh);
        id++;

        //El2
        Topol[0] = 1; Topol[1] = 2; Topol[2] = 51; Topol[3] = 49;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,mat.fMat[MatIds::ElMat],*Mesh);
        id++;

        //El3
        Topol[0] = 2; Topol[1] = 3; Topol[2] = 8; Topol[3] = 51;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,mat.fMat[MatIds::ElMat],*Mesh);
        id++;

        //El4
        Topol[0] = 3; Topol[1] = 4; Topol[2] = 7; Topol[3] = 8;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,mat.fMat[MatIds::ElMat],*Mesh);
        id++;

        //El5
        Topol[0] = 4; Topol[1] = 5; Topol[2] = 6; Topol[3] = 7;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,mat.fMat[MatIds::ElMat],*Mesh);
        id++;

        //El6
        Topol[0] = 8; Topol[1] = 9; Topol[2] = 53; Topol[3] = 51;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,mat.fMat[MatIds::ElMat],*Mesh);
        id++;

        //El7
        Topol[0] = 24; Topol[1] = 53; Topol[2] = 55; Topol[3] = 54;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,mat.fMat[MatIds::PreDistribAnelMat],*Mesh);
        id++;

        //El8
        Topol[0] = 53; Topol[1] = 9; Topol[2] = 56; Topol[3] = 55;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,mat.fMat[MatIds::PreDistribAnelMat],*Mesh);
        id++;

        //El9
        Topol[0] = 9; Topol[1] = 23; Topol[2] = 57; Topol[3] = 56;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,mat.fMat[MatIds::PreDistribAnelMat],*Mesh);
        id++;

        //El10
        Topol[0] = 54; Topol[1] = 55; Topol[2] = 59; Topol[3] = 58;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,mat.fMat[MatIds::PreDistribAnelMat],*Mesh);
        id++;

        //El11
        Topol[0] = 55; Topol[1] = 56; Topol[2] = 60; Topol[3] = 59;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,mat.fMat[MatIds::PreDistribMat],*Mesh);
        id++;

        //El12
        Topol[0] = 56; Topol[1] = 57; Topol[2] = 61; Topol[3] = 60;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,mat.fMat[MatIds::PreDistribMat],*Mesh);
        id++;

        //El13
        Topol[0] = 58; Topol[1] = 59; Topol[2] = 29; Topol[3] = 25;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,mat.fMat[MatIds::PreDistribAnelMat],*Mesh);
        id++;

        //El14
        Topol[0] = 59; Topol[1] = 60; Topol[2] = 10; Topol[3] = 29;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,mat.fMat[MatIds::PreDistribAnelMat],*Mesh);
        id++;

        //El15
        Topol[0] = 60; Topol[1] = 61; Topol[2] = 26; Topol[3] = 10;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,mat.fMat[MatIds::PreDistribAnelMat],*Mesh);
        id++;

        //El16
        Topol[0] = 10; Topol[1] = 11; Topol[2] = 31; Topol[3] = 29;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,mat.fMat[MatIds::ElMat],*Mesh);
        id++;

        //El17
        Topol[0] = 11; Topol[1] = 28; Topol[2] = 33; Topol[3] = 31;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,mat.fMat[MatIds::ElMat],*Mesh);
        id++;

        //El18
        Topol[0] = 11; Topol[1] = 12; Topol[2] = 13; Topol[3] = 28;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,mat.fMat[MatIds::ElMat],*Mesh);
        id++;

        //El19
        Topol[0] = 28; Topol[1] = 27; Topol[2] = 35; Topol[3] = 33;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,mat.fMat[MatIds::ElMat],*Mesh);
        id++;

        //El20
        Topol[0] = 13; Topol[1] = 14; Topol[2] = 27; Topol[3] = 28;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,mat.fMat[MatIds::ElMat],*Mesh);
        id++;

        //El21
        Topol[0] = 27; Topol[1] = 16; Topol[2] = 37; Topol[3] = 35;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,mat.fMat[MatIds::ElMat],*Mesh);
        id++;

        //El22
        Topol[0] = 14; Topol[1] = 15; Topol[2] = 16; Topol[3] = 27;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,mat.fMat[MatIds::ElMat],*Mesh);
        id++;

        //El23
        Topol[0] = 16; Topol[1] = 17; Topol[2] = 39; Topol[3] = 37;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,mat.fMat[MatIds::ElMat],*Mesh);
        id++;

        //El24
        Topol[0] = 17; Topol[1] = 20; Topol[2] = 41; Topol[3] = 39;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,mat.fMat[MatIds::ElMat],*Mesh);
        id++;

        //El25
        Topol[0] = 17; Topol[1] = 18; Topol[2] = 19; Topol[3] = 20;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,mat.fMat[MatIds::ElMat],*Mesh);
        id++;

        //El26
        Topol[0] = 20; Topol[1] = 21; Topol[2] = 43; Topol[3] = 41;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,mat.fMat[MatIds::ElMat],*Mesh);
        id++;

        //El27
        Topol[0] = 21; Topol[1] = 22; Topol[2] = 45; Topol[3] = 43;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,mat.fMat[MatIds::ElMat],*Mesh);
        id++;
    }

#ifndef Arc3D
    Topol.Resize(3);

    Topol[0] = 31;
    Topol[1] = 29;
    Topol[2] = 30;
    new TPZGeoElRefPattern<TPZArc3D> (id,Topol,mat.fMat[MatIds::PNmat],*Mesh);
    id++;

    Topol[0] = 33;
    Topol[1] = 31;
    Topol[2] = 32;
    new TPZGeoElRefPattern<TPZArc3D> (id,Topol,mat.fMat[MatIds::PNmat],*Mesh);
    id++;

    Topol[0] = 35;
    Topol[1] = 33;
    Topol[2] = 34;
    new TPZGeoElRefPattern<TPZArc3D> (id,Topol,mat.fMat[MatIds::PNmat],*Mesh);
    id++;

    Topol[0] = 37;
    Topol[1] = 35;
    Topol[2] = 36;
    new TPZGeoElRefPattern<TPZArc3D> (id,Topol,mat.fMat[MatIds::PNmat],*Mesh);
    id++;

    Topol[0] = 39;
    Topol[1] = 37;
    Topol[2] = 38;
    new TPZGeoElRefPattern<TPZArc3D> (id,Topol,mat.fMat[MatIds::PNmat],*Mesh);
    id++;

    Topol[0] = 41;
    Topol[1] = 39;
    Topol[2] = 40;
    new TPZGeoElRefPattern<TPZArc3D> (id,Topol,mat.fMat[MatIds::PNmat],*Mesh);
    id++;

    Topol[0] = 43;
    Topol[1] = 41;
    Topol[2] = 42;
    new TPZGeoElRefPattern<TPZArc3D> (id,Topol,mat.fMat[MatIds::PNmat],*Mesh);
    id++;

    Topol[0] = 45;
    Topol[1] = 43;
    Topol[2] = 44;
    new TPZGeoElRefPattern<TPZArc3D> (id,Topol,mat.fMat[MatIds::PNmat],*Mesh);
    id++;

    Topol[0] = 47;
    Topol[1] = 45;
    Topol[2] = 46;
    new TPZGeoElRefPattern<TPZArc3D> (id,Topol,mat.fMat[MatIds::PNmat],*Mesh);
    id++;

    Topol[0] = 49;
    Topol[1] = 47;
    Topol[2] = 48;
    new TPZGeoElRefPattern<TPZArc3D> (id,Topol,mat.fMat[MatIds::PNmat],*Mesh);
    id++;

    Topol[0] = 51;
    Topol[1] = 49;
    Topol[2] = 50;
    new TPZGeoElRefPattern<TPZArc3D> (id,Topol,mat.fMat[MatIds::PNmat],*Mesh);
    id++;

    Topol[0] = 53;
    Topol[1] = 51;
    Topol[2] = 52;
    new TPZGeoElRefPattern<TPZArc3D> (id,Topol,mat.fMat[MatIds::PNmat],*Mesh);
    id++;

#else

    Topol.Resize(2);

    Topol[0] = 31;
    Topol[1] = 29;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::PNmat],*Mesh);
    id++;

    Topol[0] = 33;
    Topol[1] = 31;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::PNmat],*Mesh);
    id++;

    Topol[0] = 35;
    Topol[1] = 33;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::PNmat],*Mesh);
    id++;

    Topol[0] = 37;
    Topol[1] = 35;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::PNmat],*Mesh);
    id++;

    Topol[0] = 39;
    Topol[1] = 37;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::PNmat],*Mesh);
    id++;

    Topol[0] = 41;
    Topol[1] = 39;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::PNmat],*Mesh);
    id++;

    Topol[0] = 43;
    Topol[1] = 41;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::PNmat],*Mesh);
    id++;

    Topol[0] = 45;
    Topol[1] = 43;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::PNmat],*Mesh);
    id++;

    Topol[0] = 47;
    Topol[1] = 45;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::PNmat],*Mesh);
    id++;

    Topol[0] = 49;
    Topol[1] = 47;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::PNmat],*Mesh);
    id++;

    Topol[0] = 51;
    Topol[1] = 49;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::PNmat],*Mesh);
    id++;

    Topol[0] = 53;
    Topol[1] = 51;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::PNmat],*Mesh);
    id++;

#endif

    Topol.Resize(2);

    Topol[0] = 0;
    Topol[1] = 1;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::Basemat],*Mesh);
    id++;

    Topol[0] = 1;
    Topol[1] = 2;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::Basemat],*Mesh);
    id++;

    Topol[0] = 2;
    Topol[1] = 3;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::OutBound],*Mesh);
    id++;

    Topol[0] = 3;
    Topol[1] = 4;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::OutBound],*Mesh);
    id++;

    Topol[0] = 4;
    Topol[1] = 5;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::OutBound],*Mesh);
    id++;

    Topol[0] = 5;
    Topol[1] = 6;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::TTmat],*Mesh);
    id++;

    Topol[0] = 6;
    Topol[1] = 7;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::OutBound],*Mesh);
    id++;

    Topol[0] = 7;
    Topol[1] = 8;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::OutBound],*Mesh);
    id++;

    Topol[0] = 8;
    Topol[1] = 9;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::OutBound],*Mesh);
    id++;

    Topol[0] = 9;
    Topol[1] = 23;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::OutBound],*Mesh);
    id++;

    Topol[0] = 23;
    Topol[1] = 57;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::OutBound],*Mesh);
    id++;

    Topol[0] = 57;
    Topol[1] = 61;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::OutBound],*Mesh);
    id++;

    Topol[0] = 61;
    Topol[1] = 26;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::OutBound],*Mesh);
    id++;

    Topol[0] = 26;
    Topol[1] = 10;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::OutBound],*Mesh);
    id++;

    Topol[0] = 10;
    Topol[1] = 11;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::OutBound],*Mesh);
    id++;

    Topol[0] = 11;
    Topol[1] = 12;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::OutBound],*Mesh);
    id++;

    Topol[0] = 12;
    Topol[1] = 13;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::OutBound],*Mesh);
    id++;

    Topol[0] = 13;
    Topol[1] = 14;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::OutBound],*Mesh);
    id++;

    Topol[0] = 14;
    Topol[1] = 15;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::FFmat],*Mesh);
    id++;

    Topol[0] = 15;
    Topol[1] = 16;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::OutBound],*Mesh);
    id++;

    Topol[0] = 16;
    Topol[1] = 17;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::FEmat],*Mesh);
    id++;

    Topol[0] = 17;
    Topol[1] = 18;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::OutBound],*Mesh);
    id++;

    Topol[0] = 18;
    Topol[1] = 19;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::PGmat],*Mesh);
    id++;

    Topol[0] = 19;
    Topol[1] = 20;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::OutBound],*Mesh);
    id++;

    Topol[0] = 20;
    Topol[1] = 21;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::OutBound],*Mesh);
    id++;

    Topol[0] = 21;
    Topol[1] = 22;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::RRBc],*Mesh);
    id++;

    Topol[0] = 22;
    Topol[1] = 0;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::RRBc],*Mesh);
    id++;

    Topol[0] = 29;
    Topol[1] = 25;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::PD2mat],*Mesh);
    id++;

    Topol[0] = 25;
    Topol[1] = 58;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::PD2mat],*Mesh);
    id++;

    Topol[0] = 58;
    Topol[1] = 54;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::PD2mat],*Mesh);
    id++;

    Topol[0] = 54;
    Topol[1] = 24;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::PD2mat],*Mesh);
    id++;

    Topol[0] = 24;
    Topol[1] = 53;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::PD2mat],*Mesh);
    id++;

     Topol[0] = 53;
     Topol[1] = 9;
     new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,mat.fMat[MatIds::PDmat],*Mesh);
     id++;

    Topol.Resize(1);
    Topol[0] = 0;
    new TPZGeoElRefPattern<TPZGeoPoint> (id,Topol,mat.fMat[MatIds::Dotmat],*Mesh);

    Mesh->BuildConnectivity();


#ifdef LOG4CXX
    {
        for(int el = 0; el < Mesh->NElements(); el++)
        {
            if(!Mesh->ElementVec()[el] || Mesh->ElementVec()[el]->Dimension() != 2) continue;
            for(int byside = 0; byside < Mesh->ElementVec()[el]->NSides(); byside++)
            {
                TPZGeoElSide NeighSide(Mesh->ElementVec()[el],byside);
                if(NeighSide.Neighbour().Exists() && NeighSide.Element()->Dimension() == 1)
                {
                    //============================= Setting Vector 1
                    TPZVec<REAL> CP_el(3,0.), XCP_el(3,0.);
                    Mesh->ElementVec()[el]->CenterPoint(byside, CP_el);
                    Mesh->ElementVec()[el]->X(CP_el,XCP_el);

                    TPZVec<REAL> CP_neigh(3,0.), XCP_neigh(3,0.);
                    NeighSide.Element()->CenterPoint(byside, CP_neigh);
                    Mesh->ElementVec()[el]->X(CP_neigh,XCP_neigh);

                    TPZVec<REAL> v1(2,0.);
                    for(int i = 0; i < 2; i++)
                    {
                        v1[i] = XCP_neigh[i] - XCP_el[i];
                    }
                    //============================= Setting Vector 2
                    TPZVec<REAL> v2(2,0.);
                    int nodeIni = NeighSide.Element()->NodeIndex(0);
                    int nodeFin = NeighSide.Element()->NodeIndex(1);
                    for(int j = 0; j < 2; j++)
                    {
                        v2[j] = Mesh->NodeVec()[nodeFin].Coord(j) - Mesh->NodeVec()[nodeIni].Coord(j);
                    }
                    //============================= Computing Vectorial Product (v1 x v2), where the interest target is just Z component
                    //============================= because GeoMesh is already in xy plane!
                    REAL compZ = v1[0]*v2[1] - v1[1]*v2[0];
                    if(compZ <= 1.E-5)
                    {
                      std::stringstream sout;
                      sout << endl << "Elemento(s) de contorno(s) 1D declarado(s) com sentido(s) errado(s)!" << endl;
                      sout << endl << "IDelement = " << NeighSide.Element()->Id() << endl;
                      sout << "Os BCs externos devem ser anti-horarios e os BCs internos horarios!" << endl;
                      sout << "Verify CxEspiral2D() Method..." << endl;
                      LOGPZ_DEBUG(logger,sout.str());
                      exit(-1);
                    }
                }
            }
        }
    }
#endif

    return Mesh;
}

TPZCompMesh * SquareMesh()
{
  TPZGeoMesh * gmesh = new TPZGeoMesh;
  //gmesh->InitializeRefPatterns();

  int Qnodes = 4;
  TPZVec < TPZVec <REAL> > NodeCoord(Qnodes);
  for(int i = 0; i < Qnodes; i++) NodeCoord[i].Resize(3,0.);
  NodeCoord[0][0] = 1.0; NodeCoord[0][1] = 0.5;
  NodeCoord[1][0] = 1.5; NodeCoord[1][1] = 0.5;
  NodeCoord[2][0] = 1.5; NodeCoord[2][1] = 4.5;
  NodeCoord[3][0] = 1.0; NodeCoord[3][1] = 4.5; // *** o eixo de revolucao serah o eixo Y !!!

  gmesh->NodeVec().Resize(Qnodes);
  TPZVec <TPZGeoNode> Node(Qnodes);
  for(int n = 0; n < Qnodes; n++)
  {
    Node[n].SetNodeId(n);
    Node[n].SetCoord(&NodeCoord[n][0]);
    gmesh->NodeVec()[n] = Node[n];
  }

    //TOPOLOGIES
  TPZVec <int> QTopol(4), DownTopol(2), ExtTopol(2), UpTopol(2), IntTopol(2);
  for(int t = 0; t < 4; t++) QTopol[t] = t;
  DownTopol[0] = 0; DownTopol[1] = 1;
  ExtTopol[0] = 1;  ExtTopol[1] = 2;
  UpTopol[0] = 2;   UpTopol[1] = 3;
  IntTopol[0] = 3;  IntTopol[1] = 0;

    //MATERIALS Id's
  int QuadMat =  1;
  int DownMat = -1;
  int ExtMat  = -2;
  int UpMat   = -3;
  int IntMat  = -4;
  int DotMat = -5;

    //Quadrilateral
  int id = 0;
  TPZGeoEl * QuadEl = new TPZGeoElRefPattern<TPZGeoQuad> (id,QTopol,QuadMat,*gmesh);
  id++;

    //BC Bottom Side
  TPZGeoElRefLess<TPZGeoBlend<pzgeom::TPZGeoLinear> > * BCdownEl = new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoLinear> > (id,DownTopol,DownMat,*gmesh);
  id++;

    //BC External Side
  TPZGeoElRefLess<TPZGeoBlend<pzgeom::TPZGeoLinear> > * BCextEl = new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoLinear> > (id,ExtTopol,ExtMat,*gmesh);
  id++;

    //BC Upper Side
  TPZGeoElRefLess<TPZGeoBlend<pzgeom::TPZGeoLinear> > * BCupEl = new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoLinear> > (id,UpTopol,UpMat,*gmesh);
  id++;

    //BC Internal Side
  TPZGeoElRefLess<TPZGeoBlend<pzgeom::TPZGeoLinear> > * BCintEl = new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoLinear> > (id,IntTopol,IntMat,*gmesh);
  id++;

  TPZGeoElBC(QuadEl,0,DotMat,*gmesh);
  gmesh->BuildConnectivity();
//**************************

///Materials
  REAL fx = 0., fy = 0.;
  TPZAutoPointer<TPZMaterial> mat = new TPZElasticityAxiMaterial(QuadMat, 2500., 0.0, fx, fy);

  TPZManVector<REAL> Orig(3);  Orig[0] = 2.5; Orig[1] = 0.; Orig[2] = 0.;
  TPZManVector<REAL> AxisZ(3); AxisZ[0] = 0.; AxisZ[1] = 1.; AxisZ[2] = 0.;
  TPZManVector<REAL> AxisR(3); AxisR[0] = 1.; AxisR[1] = 0.; AxisR[2] = 0.;

  (dynamic_cast<TPZElasticityAxiMaterial*>(mat.operator->()))->SetOrigin(Orig, AxisZ, AxisR);

  TPZFMatrix Bval1(2,2,0.), Bval2(2,1,0.);
  Bval2(0,0) = 1.;
  TPZAutoPointer<TPZMaterial> bcmatB = mat->CreateBC(mat, DownMat, 3, Bval1, Bval2);

  TPZFMatrix Uval1(2,2,0.), Uval2(2,1,0.);
  Uval2(0,0) = 1.;
  TPZAutoPointer<TPZMaterial> bcmatE = mat->CreateBC(mat, ExtMat, 3, Uval1, Uval2);

  TPZFMatrix Eval1(2,2,0.), Eval2(2,1,0.);
  Eval2(0,0) = 1.;
  TPZAutoPointer<TPZMaterial> bcmatU = mat->CreateBC(mat, UpMat, 3, Eval1, Eval2);

  TPZFMatrix Ival1(2,2,0.), Ival2(2,1,0.);
  Ival2(0,0) = 1.;
  TPZAutoPointer<TPZMaterial> bcmatI = mat->CreateBC(mat, IntMat, 3, Ival1, Ival2);

  TPZFMatrix Dval1(2,2,0.), Dval2(2,1,0.);
  Dval2(1,1) = 1.e-10;
  TPZAutoPointer<TPZMaterial> Dotmat = mat->CreateBC(mat, DotMat, 2, Dval1, Dval2);
//**************************

///Computacional Mesh
  const int dim = 2;
  TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
  cmesh->InsertMaterialObject(mat);
  cmesh->InsertMaterialObject(bcmatB);
  cmesh->InsertMaterialObject(bcmatE);
  cmesh->InsertMaterialObject(bcmatU);
  cmesh->InsertMaterialObject(bcmatI);
//**************************

  TPZCompEl::SetgOrder(2);
  cmesh->AutoBuild();

  return cmesh;
}

TPZCompMesh * SpiralMesh(MatIds &MatId)
{
/*  int ElMat = 1;
  int PreDistribMat = 2;
  int PreDistribAnelMat = 3;

  int Nonemat = -100;
  int Dotmat = -9;
  int PNmat = -8;
  int TTmat = -2;
  int PDmat = -3;
  int PD2mat = -8;
  int FFmat = -4;
  int FEmat = -5;
  int PGmat = -6;
  int Basemat = -1;
  int RRBc = -7;*/


  ///Geometria Caixa Espiral 2D

  double Bb = 8.6; double Hr = 10.8; double Bt = 8.4; double Hl = 10.15; double Cx = 4.9; double Cy = 5.9; double R = 3.4;
  double h = 2.1; double b = 0.55; double Dx = 0.65; double Dy = 3.5; double e1 = 0.25; double e2 = 0.45; double e3 = 0.1;
  double X1 = 3.8; double X2 = 1.; double X3 = 1.9; double X4 = 0.6; double X5 = 0.3; double h1 = 0.45; double h2 = 0.45; double Py = 3.3; double es = 0.1;
  TPZGeoMesh * gmesh = CxEspiral2D(Bb, Hr, Bt, Hl, Cx, Cy, R, h, b, Dx, Dy, e1, e2, e3, X1, X2, X3, X4, X5, h1, h2, Py, es, MatId);
///Materials
  REAL fx = 0., fy = MatId.f_rho; //fy = -0.025;

  double E = 18.5*1000.; // Young modulus = 18.5 GPa = 18.5*10^3 MPa
  double nu = 0.2; // poisson coefficient


  TPZAutoPointer<TPZMaterial> mat = new TPZElasticityAxiMaterial(MatIds::ElMat, E, nu, fx, fy);
  TPZElasticityAxiMaterial *aximat = dynamic_cast<TPZElasticityAxiMaterial*>(mat.operator->());
  TPZManVector<REAL> Orig(3);  Orig[0] = 15.;  Orig[1] = 0.;  Orig[2] = 0.;
  TPZManVector<REAL> AxisZ(3); AxisZ[0] = 0.; AxisZ[1] = 1.; AxisZ[2] = 0.;
  TPZManVector<REAL> AxisR(3); AxisR[0] = -1.; AxisR[1] = 0.; AxisR[2] = 0.;
  aximat->SetOrigin(Orig, AxisZ, AxisR);
  REAL EPre = 0.116 * 210.e3;
  REAL nuPre = 0.2;
  TPZAutoPointer<TPZMaterial> mat2 = new TPZElasticityAxiMaterial(MatIds::PreDistribMat, EPre, nuPre, fx, fy);
  TPZElasticityAxiMaterial *aximat2 = dynamic_cast<TPZElasticityAxiMaterial*>(mat2.operator->());
  aximat2->SetOrigin(Orig, AxisZ, AxisR);
  REAL EPreAnel = 210.e3;
  REAL nuPreAnel = 0.2;
  TPZAutoPointer<TPZMaterial> mat3 = new TPZElasticityAxiMaterial(MatIds::PreDistribAnelMat, EPreAnel, nuPreAnel, fx, fy);
  TPZElasticityAxiMaterial *aximat3 = dynamic_cast<TPZElasticityAxiMaterial*>(mat3.operator->());
  aximat3->SetOrigin(Orig, AxisZ, AxisR);

  //MohrCoulomb data
  double phi = M_PI/6.;
  double fc = 15./1.4;//15./1.4; //fck = 15 MPa, fc = fck/1.4
  double c = fc*(1. - sin(phi))/(2.*cos(phi));
  aximat->SetMohrCoulomb(c,phi);
  aximat2->SetMohrCoulomb(0.,0.);
  aximat3->SetMohrCoulomb(0.,0.);

  TPZFMatrix Base1(2,2,0.), Base2(2,1,0.);
  Base1(1,1) = 1.E9;
  TPZAutoPointer<TPZMaterial> bcmatBase = mat->CreateBC(mat, MatIds::Basemat, 2, Base1, Base2);

  TPZFMatrix PGForce1(2,2,0.), PGForce2(2,1,0.);
  PGForce2(1,0) = -700.e-3;
  TPZAutoPointer<TPZMaterial> bcmatPGForce = mat->CreateBC(mat, MatIds::PGmat, 1, PGForce1, PGForce2);

  TPZFMatrix FEForce1(2,2,0.), FEForce2(2,1,0.);
  FEForce2(1,0) = -53.41e-3;
  TPZAutoPointer<TPZMaterial> bcmatFEForce = mat->CreateBC(mat, MatIds::FEmat, 1, FEForce1, FEForce2);

  TPZFMatrix FFForce1(2,2,0.), FFForce2(2,1,0.);
  FFForce2(1,0) = -254.94e-3;
  TPZAutoPointer<TPZMaterial> bcmatFFForce = mat->CreateBC(mat, MatIds::FFmat, 1, FFForce1, FFForce2);

  TPZFMatrix PDForce1(2,2,0.), PDForce2(2,1,0.);
  PDForce2(1,0) = -2759.68e-3;
  TPZAutoPointer<TPZMaterial> bcmatPDForce = mat->CreateBC(mat, MatIds::PDmat, 1, PDForce1, PDForce2);

  TPZFMatrix TTForce1(2,2,0.), TTForce2(2,1,0.);
  TTForce2(1,0) = -61.04e-3;
  TPZAutoPointer<TPZMaterial> bcmatTTForce = mat->CreateBC(mat, MatIds::TTmat, 1, TTForce1, TTForce2);

  TPZFMatrix RRForce1(2,2,0.), RRForce2(2,1,0.);
  RRForce1(0,0) = 1.e9;
  TPZAutoPointer<TPZMaterial> bcmatRR = mat->CreateBC(mat, MatIds::RRBc, 1, RRForce1, RRForce2);

  TPZFMatrix PNForce1(2,2,0.), PNForce2(2,1,0.);
  PNForce2(0,0) = 0.0;//Sem Pressao Interna
//   PNForce2(0,0) = 0.286;//Sem Golpe com Alivio
//   PNForce2(0,0) = 0.434;//Com Golpe com Alivio
//   PNForce2(0,0) = 0.565;//Sem Golpe sem Alivio
//   PNForce2(0,0) = 0.750;//Com Golpe sem Alivio

  TPZAutoPointer<TPZMaterial> bcmatPNForce = mat->CreateBC(mat, MatIds::PNmat, 3, PNForce1, PNForce2);

  TPZFMatrix Dot1(2,2,0.), Dot2(2,1,0.);
  Dot1(1,1) = 0.01;
  TPZAutoPointer<TPZMaterial> bcmatDot = mat->CreateBC(mat, MatIds::Dotmat, 2, Dot1, Dot2);
//**************************

///Computational Mesh
  const int dim = 2;
  TPZCompEl::SetgOrder(6);
  TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
  cmesh->InsertMaterialObject(mat);
  cmesh->InsertMaterialObject(mat2);
  cmesh->InsertMaterialObject(mat3);
  cmesh->InsertMaterialObject(bcmatBase);
  cmesh->InsertMaterialObject(bcmatPGForce);
  cmesh->InsertMaterialObject(bcmatPNForce);
  cmesh->InsertMaterialObject(bcmatFEForce);
  cmesh->InsertMaterialObject(bcmatFFForce);
  cmesh->InsertMaterialObject(bcmatPDForce);
  cmesh->InsertMaterialObject(bcmatTTForce);
  cmesh->InsertMaterialObject(bcmatRR);
  cmesh->InsertMaterialObject(bcmatDot);
//**************************

  cmesh->AutoBuild();

  return cmesh;
}

void ApplyDisplacements(TPZCompMesh & cmesh)
{
  TPZElasticityAxiMaterial * AxiMat = (dynamic_cast<TPZElasticityAxiMaterial*>(cmesh.MaterialVec()[1].operator->()));
  TPZManVector<REAL> AxisR = AxiMat->GetAxisR();
  TPZManVector<REAL> AxisZ = AxiMat->GetAxisZ();
  TPZManVector<REAL> Orig = AxiMat->GetOrigin();

  double alpha = 0.1;

  for(int compel = 0; compel < cmesh.NElements(); compel++)
  {
    for(int Conect = 0; Conect < cmesh.ElementVec()[compel]->NConnects(); Conect++)
    {
      TPZCompEl * cel = cmesh.ElementVec()[compel];
      TPZGeoEl * gel = cel->Reference();

      if(gel->SideDimension(Conect) == 0)
      {
        int blockIndex = cel->Connect(Conect).SequenceNumber();
        int nsides = gel->NSides();
        TPZTransform Tr = gel->SideToSideTransform(Conect,nsides-1);
        TPZVec<REAL> From(0), To(gel->Dimension()), xco(3);
        Tr.Apply(From,To);
        gel->X(To,xco);

        double R = (xco[0] - Orig[0])*AxisR[0] + (xco[1] - Orig[1])*AxisR[1] + (xco[2] - Orig[2])*AxisR[2];
        double Z = (xco[0] - Orig[0])*AxisZ[0] + (xco[1] - Orig[1])*AxisZ[1] + (xco[2] - Orig[2])*AxisZ[2];

        cmesh.Block()(blockIndex,0,0,0) = alpha*R;
        cmesh.Block()(blockIndex,0,1,0) = alpha*Z;
      }
    }
  }
}

void WriteMesh(TPZGeoMesh *mesh,std::ofstream &arq)
{
  arq << "object 1 class array type float rank 1 shape 3 items ";
  arq << mesh->NodeVec().NElements() << " data follows" << std::endl;
  int i;

    //Print Nodes
  for (i=0;i<mesh->NodeVec().NElements(); i++)
  {
    TPZGeoNode *node = &mesh->NodeVec()[i];
    arq << node->Coord(0) << "\t" << node->Coord(1) << "\t" << node->Coord(2) << std::endl;
  }

  int numelements = 0;
  for (i=0;i<mesh->ElementVec().NElements();i++)
  {
    TPZGeoEl *el = mesh->ElementVec()[i];
    if ( !el /*|| el->Dimension() != 2*/) continue;
    numelements++;
  }
  arq << "object 2 class array type integer rank 1 shape 8 items ";
  arq << numelements << " data follows" << std::endl;

  TPZVec<int> elementtype(mesh->ElementVec().NElements(),0);
  for (i=0;i<mesh->ElementVec().NElements();i++)
  {
    TPZGeoEl *el = mesh->ElementVec()[i];
    if ( !el /*|| el->Dimension() != 2*/) continue;
    WriteElement (el,i,arq,elementtype);
  }
  arq << "attribute \"element type\" string \"cubes\"" << std::endl
      << "attribute \"ref\" string \"positions\"" << std::endl;
  arq << "object 3 class array type integer rank 0 items ";
  arq << numelements << " data follows" << std::endl;

  for (i=0;i<mesh->ElementVec().NElements();i++)
  {
    TPZGeoEl *el = mesh->ElementVec()[i];
    if ( !el /*|| el->Dimension()!= 2*/) continue;
    arq << elementtype[i] << std::endl;
  }
  arq << "attribute \"dep\" string \"connections\"" << std::endl;
  arq << "object 4 class field" << std::endl
      << "component \"positions\" value 1" << std::endl
      << "component \"connections\" value 2" << std::endl
      << "component \"data\" value 3" << std::endl;
}

void WriteElement (TPZGeoEl *el,int elindex, std::ofstream &arq,TPZVec<int> &elementtype)
{
  int ncon = el->NNodes();
  elementtype[elindex] = ncon;
  switch (ncon)
  {
    case (2):
    {
            //rib
      int ni = el->NodeIndex(0);
      int nf = el->NodeIndex(1);
      arq << ni << "\t" << nf << "\t" << ni << "\t" << nf << "\t" << ni << "\t" << nf << "\t" << ni << "\t" << nf << std::endl;
      break;
    }

    case (3):
    {
            //triangle
      int n0 = el->NodeIndex(0);
      int n1 = el->NodeIndex(1);
      int n2 = el->NodeIndex(2);
      arq << n0 << "\t" << n1 << "\t" << n2 << "\t" << n2 << "\t" << n0 << "\t" << n1 << "\t" << n2 << "\t" << n2 << std::endl;
      break;
    }

    case (4):
    {
      if (el->Dimension() == 2)
      {
                //quad
        int n0 = el->NodeIndex(0);
        int n1 = el->NodeIndex(1);
        int n2 = el->NodeIndex(3);
        int n3 = el->NodeIndex(2);
        arq << n0 << "\t" << n1 << "\t" << n2 << "\t" << n3 << "\t" << n0 << "\t" << n1 << "\t" << n2 << "\t" << n3 << std::endl;
      }
      else
      {
                //tetrahedre
        int n0 = el->NodeIndex(0);
        int n1 = el->NodeIndex(1);
        int n2 = el->NodeIndex(2);
        int n3 = el->NodeIndex(3);
        arq << n0 << "\t" << n1 << "\t" << n2 << "\t" << n2 << "\t" << n3 << "\t" << n3 << "\t" << n3 << "\t" << n3 << std::endl;
      }
      break;
    }

    case (5):
    {
            //pyramid
      int n0 = el->NodeIndex(0);
      int n1 = el->NodeIndex(1);
      int n2 = el->NodeIndex(3);
      int n3 = el->NodeIndex(2);
      int n4 = el->NodeIndex(4);
      arq << n0 << "\t" << n1 << "\t" << n2 << "\t" << n3 << "\t" << n4 << "\t" << n4 << "\t" << n4 << "\t" << n4 << std::endl;
      break;
    }

    case (6):
    {
            //pyramid
      int n0 = el->NodeIndex(0);
      int n1 = el->NodeIndex(1);
      int n2 = el->NodeIndex(2);
      int n3 = el->NodeIndex(3);
      int n4 = el->NodeIndex(4);
      int n5 = el->NodeIndex(5);
      arq << n0 << "\t" << n1 << "\t" << n2 << "\t" << n2 << "\t" << n3 << "\t" << n4 << "\t" << n5 << "\t" << n5 << std::endl;
      break;
    }

    case (8):
    {
      int n0 = el->NodeIndex(0);
      int n1 = el->NodeIndex(1);
      int n2 = el->NodeIndex(3);
      int n3 = el->NodeIndex(2);
      int n4 = el->NodeIndex(4);
      int n5 = el->NodeIndex(5);
      int n6 = el->NodeIndex(7);
      int n7 = el->NodeIndex(6);
      arq << n0 << "\t" << n1 << "\t" << n2 << "\t" << n3 << "\t" << n4 << "\t" << n5 << "\t" << n6 << "\t" << n7 << std::endl;
      break;
    }
    default:
      std::cout << "Erro..." << std::endl;
  }
  return;
}

void IntegralForces(TPZCompMesh * cmesh, int MatId, TPZVec<REAL> &forces)
{
  forces.Resize(2);
  forces[0] = 0.;
  forces[1] = 0.;
  TPZElasticityAxiMaterial * Mat = dynamic_cast<TPZElasticityAxiMaterial*>(cmesh->MaterialVec()[1].operator->());
  for(int el = 0; el < cmesh->ElementVec().NElements(); el++)
  {
    TPZCompEl * cel = cmesh->ElementVec()[el];
    if(cel->Material()->Id() == MatId)
    {
      TPZGeoElSide geoside(cel->Reference(),2);
      TPZGeoElSide neigh = geoside.Neighbour();
      TPZGeoElSide internalSide = TPZGeoElSide(neigh.Element(),neigh.Element()->NSides() - 1);

      TPZTransform transf1(1,1), transf2(1,2), transf(1,2);
      transf1 = geoside.NeighbourSideTransform(neigh);
      transf2 = neigh.SideToSideTransform(internalSide);
      transf = transf2.Multiply(transf1);

      TPZIntPoints * rule = cel->Reference()->CreateSideIntegrationRule(geoside.Side(),15);

      TPZVec<REAL> location(1);
      REAL w;
      int npts = rule->NPoints();
      for(int pt = 0; pt < npts; pt++)
      {
        rule->Point(pt, location,w);
        TPZFMatrix jac(1,1), jacinv(1,1), axes(1,3);
        REAL detJac;
        geoside.Jacobian(location,jac, axes, detJac, jacinv);

        TPZVec<REAL> locout(2);
        transf.Apply(location,locout);

        TPZVec<REAL> x(3);
        geoside.X(location,x);
        REAL R = Mat->ComputeR(x);
        R = fabs(R);

        TPZVec<REAL> sol;
        neigh.Element()->Reference()->Solution(locout,4,sol);
        forces[0] += w * (2.*M_PI*R) * detJac * sol[0];
        neigh.Element()->Reference()->Solution(locout,0,sol);
        forces[1] += w * (2.*M_PI*R) * detJac * sol[1];
      }

      delete rule;
    }
  }
  int bye=0;
}

REAL CalcArea(TPZCompMesh *mesh)
{
  REAL area = 0.;
  TPZElasticityAxiMaterial * Mat = dynamic_cast<TPZElasticityAxiMaterial*>(mesh->MaterialVec()[MatIds::ElMat].operator->());
  for(int el = 0; el < mesh->ElementVec().NElements(); el++)
  {
    TPZCompEl * cel = mesh->ElementVec()[el];
    if(cel->Material()->Id() == MatIds::ElMat)
    {

      TPZGeoEl *gel = cel->Reference();
      TPZIntPoints * rule = cel->Reference()->CreateSideIntegrationRule(gel->NSides()-1,15);

      TPZVec<REAL> location(2);
      REAL w;
      int npts = rule->NPoints();
      for(int pt = 0; pt < npts; pt++)
      {
        rule->Point(pt, location,w);
        TPZFMatrix jac(2,2), jacinv(2,2), axes(2,3);
        REAL detJac;
        cel->Reference()->Jacobian(location,jac, axes, detJac, jacinv);
        TPZManVector<REAL> x(3,0.);
        cel->Reference()->X(location,x);
        REAL R = fabs(Mat->ComputeR(x));

        area += w*detJac*2.*M_PI*R;
      }

      delete rule;
    }
  }
  return area;

}


void MathOutUP(TPZBlendNACA * naca)
{
  cout << endl << endl << "NacaPtosUp = {";
  for(int pos = int(naca->fX0[0]); pos <= int(naca->fX0[0] + naca->fCord); pos++)
  {
    cout << "{" << naca->xua(pos) << "," << naca->yua(pos) << "}";
    if(pos != int(naca->fX0[0] + naca->fCord)) cout << ",";
  }
  cout << "};" << endl;
  cout << "A=ListPlot[NacaPtosUp, Joined -> True,AspectRatio -> 0.2,PlotRange->{{";
  cout << naca->fX0[0] << "," << naca->fX0[0] + int(naca->fCord) << "},{-5,5}}];" << endl << endl;
}

void MathOutDW(TPZBlendNACA * naca)
{
  cout << "NacaPtosDw = {";
  for(int pos = int(naca->fX0[0]); pos <= int(naca->fX0[0] + naca->fCord); pos++)
  {
    cout << "{" << naca->xla(pos) << "," << naca->yla(pos) << "}";
    if(pos != int(naca->fX0[0] + naca->fCord)) cout << ",";
  }
  cout << "};" << endl;
  cout << "B=ListPlot[NacaPtosDw, Joined -> True,AspectRatio -> 0.2];" << endl;
}




