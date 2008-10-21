#include <sstream>

#include "pzvec.h"
#include "pzcmesh.h"
#include "pzdebug.h"
#include "pzcheckgeom.h"

#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgeoelside.h"

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
#include "tpzquadratictetra.h"
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
#include <pzcompel.h>
#include "tpzellipse.h"
#include "tpzblendnaca.h"
#include "pzelasAXImat.h"
#include "pzmaterialdata.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("main"));
#endif

#include <sstream>
using namespace std;
using namespace pzgeom;
using namespace pzshape;
using namespace pzrefine;

void WriteElement(TPZGeoEl *el,int elindex, std::ofstream &arq,TPZVec<int> &elementtype);
void WriteMesh(TPZGeoMesh *mesh,std::ofstream &arq);

TPZGeoMesh * HorizWellGenE(double Wr, double Hr, double Lr, double Lw, double Dw, int part, double excentr);

TPZGeoMesh * DefaultMesh(double R);
TPZGeoMesh * CurvedTriangleMesh(double R);
TPZGeoMesh * BlendTriangleMesh(double R);
TPZGeoMesh * CxEspiral2D(double Bb, double Hr, double Bt, double Hl, double Cx, double Cy,
                         double R, double h, double b, double Dx, double Dy, double e1,
                         double e2, double e3, double X1, double X2, double X3, double X4,
                         double X5, double h1, double h2, double Py);
void AuxEspiral(TPZVec <REAL> &vecIn, TPZVec <REAL> &vecOut, TPZVec <REAL> &C, double R);

/** by Caju 2007 */
/** Input:  gMesh.NodeVec(), IdIni, IdFin, Tol */
/** Output: NodesHunted */
/** This method computes the nodes that belongs to the IdIni->IdFin line */
void NodesHunter(TPZGeoMesh &gMesh, vector<int>& NodesHunted, int IdIni, int IdFin, double Tol = 1.E-3);

void IntegralForces(TPZCompMesh * cmesh, int MatId, TPZVec<REAL> &forces);

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


void ApplyDisplacements(TPZCompMesh & cmesh);

TPZCompMesh * SquareMesh();

TPZCompMesh * SpiralMesh();

/**
 * 
 * @param argc 
 * @param argv[] 
 * @return 
 */
int main(int argc, char *argv[])
{
#ifdef LOG4CXX
  InitializePZLOG();
#endif

//     TPZCompMesh * cmesh = SquareMesh();
    TPZCompMesh * cmesh = SpiralMesh();

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
                          sout << endl << "Elemento: " << el << " Lado: " << side << " -> NOT Ok!" << endl;
                  }
              }
          }
    }
    LOGPZ_DEBUG(logger,sout.str());
}
#endif


//////////////////////////////////////////////////////////////////////

    ///Inserindo os deslocamentos "na unha!"
//     ApplyDisplacements(*cmesh);

// #ifdef LOG4CXX
// {
//   std::stringstream sout;
//   cmesh->Reference()->Print(sout);
//   LOGPZ_DEBUG(logger,sout.str());
// }
// 
//     {
//         std::stringstream sout;
//         cmesh->Print(sout);
//         LOGPZ_DEBUG(logger,sout.str());
//     }
// #endif

    TPZAnalysis an(cmesh);
    TPZSkylineStructMatrix full(cmesh);
    an.SetStructuralMatrix(full);
    TPZStepSolver step;
    step.SetDirect(ECholesky);
    an.SetSolver(step);
    an.Run();

    TPZVec<std::string> scalnames(4), vecnames(3);
    vecnames[0] = "Eigenvector1";
    vecnames[1] = "Eigenvector2";
    vecnames[2] = "Eigenvector3";
    scalnames[0] = "Sigmarr";
    scalnames[1] = "Sigmazz";
    scalnames[2] = "Sigmatt";
    scalnames[3] = "Taurz";
    std::string plotfile("saida.dx");
    an.DefineGraphMesh(2,scalnames,vecnames,plotfile);
    an.PostProcess(2);

    std::ofstream out("malha.txt");
    an.Print("nothing",out);

#ifdef LOG4CXX
{
  std::stringstream sout;
  TPZVec<REAL> forces;
  IntegralForces(cmesh, -8, forces);
  sout << "\n===========================================================\n" << forces << "\n===========================================================\n";
  LOGPZ_DEBUG(logger,sout.str());
}
#endif

///End Malha Axi-simetrica

    return EXIT_SUCCESS;
}


void ApplyDisplacements(TPZCompMesh & cmesh)
{
      TPZElasticityAxiMaterial * AxiMat = (dynamic_cast<TPZElasticityAxiMaterial*>(cmesh.MaterialVec()[1].operator->()));
      std::vector<REAL> AxisR = AxiMat->GetAxisR();
      std::vector<REAL> AxisZ = AxiMat->GetAxisZ();
      std::vector<REAL> Orig = AxiMat->GetOrigin();

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










TPZGeoMesh * DefaultMesh(double R)
{
    TPZGeoMesh * Mesh = new TPZGeoMesh;

    int Qnodes = 15;
    TPZVec < TPZVec <REAL> > NodeCoord(Qnodes);
    for(int i = 0; i < Qnodes; i++) NodeCoord[i].Resize(3);
    TPZFMatrix Rot(3,3,0.);
    TPZFMatrix IniVec1(3,1), IniVec2(3,1), IniVec3(3,1);
    TPZFMatrix FinalVec(3,1);

    double r00, r01, r10, r11, th;
//     double Pi = atan(1.)*4.;

    IniVec1(0,0) = 1.*R; IniVec1(1,0) = 0.; IniVec1(2,0) = 0.;
    IniVec2(0,0) = 2.*R; IniVec2(1,0) = 0.; IniVec2(2,0) = 0.;
    IniVec3(0,0) = 3.*R; IniVec3(1,0) = 0.; IniVec3(2,0) = 0.;

    for(int i = 0; i < 5; i++)
    {
        th = i*Pi/8.;
        r00 = cos(th); r01 = -sin(th);
        r10 = sin(th); r11 = cos(th);

        Rot(0,0) = r00; Rot(0,1) = r01;
        Rot(1,0) = r10; Rot(1,1) = r11; Rot(2,2) = 1.;

        Rot.Multiply(IniVec1,FinalVec);
        NodeCoord[3*i][0] = FinalVec(0,0);
        NodeCoord[3*i][1] = FinalVec(1,0);
        NodeCoord[3*i][2] = FinalVec(2,0);

        Rot.Multiply(IniVec2,FinalVec);
        NodeCoord[3*i+1][0] = FinalVec(0,0);
        NodeCoord[3*i+1][1] = FinalVec(1,0);
        NodeCoord[3*i+1][2] = FinalVec(2,0);

        Rot.Multiply(IniVec3,FinalVec);
        NodeCoord[3*i+2][0] = FinalVec(0,0);
        NodeCoord[3*i+2][1] = FinalVec(1,0);
        NodeCoord[3*i+2][2] = FinalVec(2,0);
    }

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
    int QuadMat = 1;
    int arcMat = -1;
    int id = 0;

    ///Quadrilaterals
    TPZVec< TPZGeoElRefLess<TPZGeoBlend<pzgeom::TPZGeoQuad> > *> QuadEl(4);
    for(int layer1 = 0; layer1 < 2; layer1++)
    {
        for(int layer2 = 0; layer2 < 2; layer2++)
        {
            Topol[0] =  6*layer1 + layer2; Topol[1] = Topol[0] + 1; Topol[2] = Topol[1] + 6; Topol[3] = Topol[0] + 6;
            QuadEl[2*layer1 + layer2] = new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,QuadMat,*Mesh);
            id++;
        }
    }

    ///Arcs
    Topol.Resize(3);
    TPZVec <TPZGeoEl*> ArcEl(6);
    for(int layer1 = 0; layer1 < 2; layer1++)
    {
        for(int layer2 = 0; layer2 < 3; layer2++)
        {
            Topol[0] = 6*layer1 + layer2;
            Topol[1] = Topol[0] + 6;
            Topol[2] = Topol[0] + 3;
            ArcEl[3*layer1 + layer2] = new TPZGeoElRefPattern<TPZArc3D> (id,Topol,arcMat,*Mesh);
            id++;
        }
    }

    Mesh->BuildBlendConnectivity();

    return Mesh;
}

TPZGeoMesh * CurvedTriangleMesh(double R)
{
    TPZGeoMesh * Mesh = new TPZGeoMesh;

    int Qnodes = 8;
    int Allnodes = 15; 
    TPZVec < TPZVec <REAL> > NodeCoord(Allnodes);
    for(int i = 0; i < Allnodes; i++) NodeCoord[i].Resize(3);
    TPZFMatrix Rot(3,3,0.);
    TPZFMatrix IniVec1(3,1), IniVec2(3,1), IniVec3(3,1);
    TPZFMatrix FinalVec(3,1);

    double r00, r01, r10, r11, th;
//     double Pi = atan(1.)*4.;

    IniVec1(0,0) = 1.*R; IniVec1(1,0) = 0.; IniVec1(2,0) = 0.;
    IniVec2(0,0) = 2.*R; IniVec2(1,0) = 0.; IniVec2(2,0) = 0.;
    IniVec3(0,0) = 3.*R; IniVec3(1,0) = 0.; IniVec3(2,0) = 0.;

    for(int i = 0; i < 5; i++)
    {
        th = i*Pi/8.;
        r00 = cos(th); r01 = -sin(th);
        r10 = sin(th); r11 = cos(th);

        Rot(0,0) = r00; Rot(0,1) = r01;
        Rot(1,0) = r10; Rot(1,1) = r11; Rot(2,2) = 1.;

        Rot.Multiply(IniVec1,FinalVec);
        NodeCoord[3*i][0] = FinalVec(0,0);
        NodeCoord[3*i][1] = FinalVec(1,0);
        NodeCoord[3*i][2] = FinalVec(2,0);

        Rot.Multiply(IniVec2,FinalVec);
        NodeCoord[3*i+1][0] = FinalVec(0,0);
        NodeCoord[3*i+1][1] = FinalVec(1,0);
        NodeCoord[3*i+1][2] = FinalVec(2,0);

        Rot.Multiply(IniVec3,FinalVec);
        NodeCoord[3*i+2][0] = FinalVec(0,0);
        NodeCoord[3*i+2][1] = FinalVec(1,0);
        NodeCoord[3*i+2][2] = FinalVec(2,0);
    }

    Mesh->NodeVec().Resize(Qnodes);
    TPZVec <TPZGeoNode> Node(Qnodes);

    Node[0].SetNodeId(0);
    Node[0].SetCoord(&NodeCoord[0][0]);
    Mesh->NodeVec()[0] = Node[0];

    Node[1].SetNodeId(1);
    Node[1].SetCoord(&NodeCoord[2][0]);
    Mesh->NodeVec()[1] = Node[1];

    Node[2].SetNodeId(2);
    Node[2].SetCoord(&NodeCoord[5][0]);
    Mesh->NodeVec()[2] = Node[2];

    NodeCoord[7][0] = R;
    NodeCoord[7][1] = R;
    NodeCoord[7][2] = 0.;
    Node[3].SetNodeId(3);
    Node[3].SetCoord(&NodeCoord[7][0]);
    Mesh->NodeVec()[3] = Node[3];

    Node[4].SetNodeId(4);
    Node[4].SetCoord(&NodeCoord[8][0]);
    Mesh->NodeVec()[4] = Node[4];

    Node[5].SetNodeId(5);
    Node[5].SetCoord(&NodeCoord[11][0]);
    Mesh->NodeVec()[5] = Node[5];

    Node[6].SetNodeId(6);
    Node[6].SetCoord(&NodeCoord[12][0]);
    Mesh->NodeVec()[6] = Node[6];

    Node[7].SetNodeId(7);
    Node[7].SetCoord(&NodeCoord[14][0]);
    Mesh->NodeVec()[7] = Node[7];


    TPZVec <int> Topol(3);

    /// *********** CHOOSE MATERIAL FOR EACH ELEMENT ***********
    int TrigMat = 1;
    int CTrigMat = 2;
    int arcMat = -1;
    int id = 0;

    ///Triangles
    TPZVec< TPZGeoElRefLess<pzgeom::TPZGeoTriangle > *> TrigEl(2);
    {
        Topol[0] = 0; Topol[1] = 1; Topol[2] = 3;
        TrigEl[0] = new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle > (id,Topol,TrigMat,*Mesh);
        id++;

        Topol[0] = 7; Topol[1] = 6; Topol[2] = 3;
        TrigEl[1] = new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle > (id,Topol,TrigMat,*Mesh);
        id++;
    }

    TPZVec< TPZGeoElRefLess<TPZGeoBlend<pzgeom::TPZGeoTriangle> > *> BTrigEl(2);
    {
        Topol[0] = 1; Topol[1] = 4; Topol[2] = 3;
        BTrigEl[0] = new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoTriangle> > (id,Topol,TrigMat,*Mesh);
        id++;

        Topol[0] = 4; Topol[1] = 7; Topol[2] = 3;
        BTrigEl[1] = new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoTriangle> > (id,Topol,TrigMat,*Mesh);
        id++;
    }

    TPZVec< TPZGeoElRefLess<TPZCurvedTriangle> *> CTrigEl(1);
    {
        Topol[0] = 3; Topol[1] = 6; Topol[2] = 0;
        CTrigEl[0] = new TPZGeoElRefPattern<TPZCurvedTriangle > (id,Topol,CTrigMat,*Mesh);
        id++;
    }

    ///Arcs
    Topol.Resize(3);
    TPZVec <TPZGeoEl*> ArcEl(2);

    Topol[0] = 1;
    Topol[1] = 4;
    Topol[2] = 2;
    ArcEl[0] = new TPZGeoElRefPattern<TPZArc3D> (id,Topol,arcMat,*Mesh);
    id++;

    Topol[0] = 4;
    Topol[1] = 7;
    Topol[2] = 5;
    ArcEl[1] = new TPZGeoElRefPattern<TPZArc3D> (id,Topol,arcMat,*Mesh);

    Mesh->BuildBlendConnectivity();

    return Mesh;
}

TPZGeoMesh * BlendTriangleMesh(double R)
{
    TPZGeoMesh * Mesh = new TPZGeoMesh;

    int Qnodes = 9;
    int Allnodes = 15; 
    TPZVec < TPZVec <REAL> > NodeCoord(Allnodes);
    for(int i = 0; i < Allnodes; i++) NodeCoord[i].Resize(3);
    TPZFMatrix Rot(3,3,0.);
    TPZFMatrix IniVec1(3,1), IniVec2(3,1), IniVec3(3,1);
    TPZFMatrix FinalVec(3,1);

    double r00, r01, r10, r11, th;
//     double Pi = atan(1.)*4.;

    IniVec1(0,0) = 1.*R; IniVec1(1,0) = 0.; IniVec1(2,0) = 0.;
    IniVec2(0,0) = 2.*R; IniVec2(1,0) = 0.; IniVec2(2,0) = 0.;
    IniVec3(0,0) = 3.*R; IniVec3(1,0) = 0.; IniVec3(2,0) = 0.;

    for(int i = 0; i < 5; i++)
    {
        th = i*Pi/8.;
        r00 = cos(th); r01 = -sin(th);
        r10 = sin(th); r11 = cos(th);

        Rot(0,0) = r00; Rot(0,1) = r01;
        Rot(1,0) = r10; Rot(1,1) = r11; Rot(2,2) = 1.;

        Rot.Multiply(IniVec1,FinalVec);
        NodeCoord[3*i][0] = FinalVec(0,0);
        NodeCoord[3*i][1] = FinalVec(1,0);
        NodeCoord[3*i][2] = FinalVec(2,0);

        Rot.Multiply(IniVec2,FinalVec);
        NodeCoord[3*i+1][0] = FinalVec(0,0);
        NodeCoord[3*i+1][1] = FinalVec(1,0);
        NodeCoord[3*i+1][2] = FinalVec(2,0);

        Rot.Multiply(IniVec3,FinalVec);
        NodeCoord[3*i+2][0] = FinalVec(0,0);
        NodeCoord[3*i+2][1] = FinalVec(1,0);
        NodeCoord[3*i+2][2] = FinalVec(2,0);
    }

    Mesh->NodeVec().Resize(Qnodes);
    TPZVec <TPZGeoNode> Node(Qnodes);

    Node[0].SetNodeId(0);
    Node[0].SetCoord(&NodeCoord[0][0]);
    Mesh->NodeVec()[0] = Node[0];

    Node[1].SetNodeId(1);
    Node[1].SetCoord(&NodeCoord[2][0]);
    Mesh->NodeVec()[1] = Node[1];

    Node[2].SetNodeId(2);
    Node[2].SetCoord(&NodeCoord[5][0]);
    Mesh->NodeVec()[2] = Node[2];

    NodeCoord[7][0] = R;
    NodeCoord[7][1] = R;
    NodeCoord[7][2] = 0.;
    Node[3].SetNodeId(3);
    Node[3].SetCoord(&NodeCoord[7][0]);
    Mesh->NodeVec()[3] = Node[3];

    Node[4].SetNodeId(4);
    Node[4].SetCoord(&NodeCoord[8][0]);
    Mesh->NodeVec()[4] = Node[4];

    Node[5].SetNodeId(5);
    Node[5].SetCoord(&NodeCoord[11][0]);
    Mesh->NodeVec()[5] = Node[5];

    Node[6].SetNodeId(6);
    Node[6].SetCoord(&NodeCoord[12][0]);
    Mesh->NodeVec()[6] = Node[6];

    Node[7].SetNodeId(7);
    Node[7].SetCoord(&NodeCoord[14][0]);
    Mesh->NodeVec()[7] = Node[7];

    Node[8].SetNodeId(8);
    Node[8].SetCoord(&NodeCoord[6][0]);
    Mesh->NodeVec()[8] = Node[8];


    TPZVec <int> Topol(3);

    /// *********** CHOOSE MATERIAL FOR EACH ELEMENT ***********
    int TrigMat = 1;
    int arcMat = -1;
    int id = 0;

    ///Triangles
    TPZVec< TPZGeoElRefLess<pzgeom::TPZGeoTriangle > *> TrigEl(2);
    {
        Topol[0] = 0; Topol[1] = 1; Topol[2] = 3;
        TrigEl[0] = new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle > (id,Topol,TrigMat,*Mesh);
        id++;

        Topol[0] = 7; Topol[1] = 6; Topol[2] = 3;
        TrigEl[1] = new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle > (id,Topol,TrigMat,*Mesh);
        id++;
    }

    TPZVec< TPZGeoElRefLess<TPZGeoBlend<pzgeom::TPZGeoTriangle> > *> BTrigEl(3);
    {
        Topol[0] = 1; Topol[1] = 4; Topol[2] = 3;
        BTrigEl[0] = new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoTriangle> > (id,Topol,TrigMat,*Mesh);
        id++;

        Topol[0] = 4; Topol[1] = 7; Topol[2] = 3;
        BTrigEl[1] = new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoTriangle> > (id,Topol,TrigMat,*Mesh);
        id++;

        Topol[0] = 0; Topol[1] = 3; Topol[2] = 6;
        BTrigEl[2] = new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoTriangle> > (id,Topol,TrigMat,*Mesh);
        id++;
    }

    ///Arcs
    Topol.Resize(3);
    TPZVec <TPZGeoEl*> ArcEl(3);

    Topol[0] = 1;
    Topol[1] = 4;
    Topol[2] = 2;
    ArcEl[0] = new TPZGeoElRefPattern<TPZArc3D> (id,Topol,arcMat,*Mesh);
    id++;

    Topol[0] = 4;
    Topol[1] = 7;
    Topol[2] = 5;
    ArcEl[1] = new TPZGeoElRefPattern<TPZArc3D> (id,Topol,arcMat,*Mesh);

    Topol[0] = 0;
    Topol[1] = 6;
    Topol[2] = 8;
    ArcEl[2] = new TPZGeoElRefPattern<TPZArc3D> (id,Topol,arcMat,*Mesh);

    Mesh->BuildBlendConnectivity();

    return Mesh;
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
                         double X5, double h1, double h2, double Py)
{
    TPZGeoMesh * Mesh = new TPZGeoMesh;
    int Qnodes = 58;
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
    NodeCoord[53][1] = Cy - h/2.;                           NodeCoord[54][1] = Cy;
    NodeCoord[53][2] = 0.;                                  NodeCoord[54][2] = 0.;

    NodeCoord[55][0] = (2.*Cx + sqrt(-h*h + 4.*R*R))/2.;    NodeCoord[56][0] = Bt;
    NodeCoord[55][1] = Cy;                                  NodeCoord[56][1] = Cy;
    NodeCoord[55][2] = 0.;                                  NodeCoord[56][2] = 0.;

    NodeCoord[57][0] = Bt + e3;
    NodeCoord[57][1] = Cy;
    NodeCoord[57][2] = 0.;
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
    int ElMat = 1;
    int PreDistribMat = 1;

    int Nonemat = -1;
    int PNmat = -1;
    int TTmat = -1;
    int PDmat = -1;
    int FFmat = -1;
    int FEmat = -1;
    int PGmat = -7;
    int Basemat = -8;
    int Dotmat = -1;

    int id = 0;
    ///Quadrilaterals
    {
        //El0
        Topol[0] = 22; Topol[1] = 0; Topol[2] = 47; Topol[3] = 45;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,ElMat,*Mesh);
        id++;

        //El1
        Topol[0] = 0; Topol[1] = 1; Topol[2] = 49; Topol[3] = 47;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,ElMat,*Mesh);
        id++;

        //El2
        Topol[0] = 1; Topol[1] = 2; Topol[2] = 51; Topol[3] = 49;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,ElMat,*Mesh);
        id++;

        //El3
        Topol[0] = 2; Topol[1] = 3; Topol[2] = 8; Topol[3] = 51;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,ElMat,*Mesh);
        id++;

        //El4
        Topol[0] = 3; Topol[1] = 4; Topol[2] = 7; Topol[3] = 8;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,ElMat,*Mesh);
        id++;

        //El5
        Topol[0] = 4; Topol[1] = 5; Topol[2] = 6; Topol[3] = 7;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,ElMat,*Mesh);
        id++;

        //El6
        Topol[0] = 8; Topol[1] = 9; Topol[2] = 53; Topol[3] = 51;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,ElMat,*Mesh);
        id++;

        //El7
        Topol[0] = 24; Topol[1] = 53; Topol[2] = 55; Topol[3] = 54;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,PreDistribMat,*Mesh);
        id++;

        //El8
        Topol[0] = 53; Topol[1] = 9; Topol[2] = 56; Topol[3] = 55;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,PreDistribMat,*Mesh);
        id++;

        //El9
        Topol[0] = 9; Topol[1] = 23; Topol[2] = 57; Topol[3] = 56;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,PreDistribMat,*Mesh);
        id++;

        //El10
        Topol[0] = 54; Topol[1] = 55; Topol[2] = 29; Topol[3] = 25;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,PreDistribMat,*Mesh);
        id++;

        //El11
        Topol[0] = 55; Topol[1] = 56; Topol[2] = 10; Topol[3] = 29;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,PreDistribMat,*Mesh);
        id++;

        //El12
        Topol[0] = 56; Topol[1] = 57; Topol[2] = 26; Topol[3] = 10;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,PreDistribMat,*Mesh);
        id++;

        //El13
        Topol[0] = 10; Topol[1] = 11; Topol[2] = 31; Topol[3] = 29;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,ElMat,*Mesh);
        id++;

        //El14
        Topol[0] = 11; Topol[1] = 28; Topol[2] = 33; Topol[3] = 31;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,ElMat,*Mesh);
        id++;

        //El15
        Topol[0] = 11; Topol[1] = 12; Topol[2] = 13; Topol[3] = 28;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,ElMat,*Mesh);
        id++;

        //El16
        Topol[0] = 28; Topol[1] = 27; Topol[2] = 35; Topol[3] = 33;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,ElMat,*Mesh);
        id++;

        //El17
        Topol[0] = 13; Topol[1] = 14; Topol[2] = 27; Topol[3] = 28;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,ElMat,*Mesh);
        id++;

        //El18
        Topol[0] = 27; Topol[1] = 16; Topol[2] = 37; Topol[3] = 35;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,ElMat,*Mesh);
        id++;

        //El19
        Topol[0] = 14; Topol[1] = 15; Topol[2] = 16; Topol[3] = 27;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,ElMat,*Mesh);
        id++;

        //El20
        Topol[0] = 16; Topol[1] = 17; Topol[2] = 39; Topol[3] = 37;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,ElMat,*Mesh);
        id++;

        //El21
        Topol[0] = 17; Topol[1] = 20; Topol[2] = 41; Topol[3] = 39;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,ElMat,*Mesh);
        id++;

        //El22
        Topol[0] = 17; Topol[1] = 18; Topol[2] = 19; Topol[3] = 20;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,ElMat,*Mesh);
        id++;

        //El23
        Topol[0] = 20; Topol[1] = 21; Topol[2] = 43; Topol[3] = 41;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,ElMat,*Mesh);
        id++;

        //El24
        Topol[0] = 21; Topol[1] = 22; Topol[2] = 45; Topol[3] = 43;
        new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,ElMat,*Mesh);
        id++;
    }

#ifndef Arc3D
    Topol.Resize(3);

    Topol[0] = 31;
    Topol[1] = 29;
    Topol[2] = 30;
    new TPZGeoElRefPattern<TPZArc3D> (id,Topol,PNmat,*Mesh);
    id++;

    Topol[0] = 33;
    Topol[1] = 31;
    Topol[2] = 32;
    new TPZGeoElRefPattern<TPZArc3D> (id,Topol,PNmat,*Mesh);
    id++;

    Topol[0] = 35;
    Topol[1] = 33;
    Topol[2] = 34;
    new TPZGeoElRefPattern<TPZArc3D> (id,Topol,PNmat,*Mesh);
    id++;

    Topol[0] = 37;
    Topol[1] = 35;
    Topol[2] = 36;
    new TPZGeoElRefPattern<TPZArc3D> (id,Topol,PNmat,*Mesh);
    id++;

    Topol[0] = 39;
    Topol[1] = 37;
    Topol[2] = 38;
    new TPZGeoElRefPattern<TPZArc3D> (id,Topol,PNmat,*Mesh);
    id++;

    Topol[0] = 41;
    Topol[1] = 39;
    Topol[2] = 40;
    new TPZGeoElRefPattern<TPZArc3D> (id,Topol,PNmat,*Mesh);
    id++;

    Topol[0] = 43;
    Topol[1] = 41;
    Topol[2] = 42;
    new TPZGeoElRefPattern<TPZArc3D> (id,Topol,PNmat,*Mesh);
    id++;

    Topol[0] = 45;
    Topol[1] = 43;
    Topol[2] = 44;
    new TPZGeoElRefPattern<TPZArc3D> (id,Topol,PNmat,*Mesh);
    id++;

    Topol[0] = 47;
    Topol[1] = 45;
    Topol[2] = 46;
    new TPZGeoElRefPattern<TPZArc3D> (id,Topol,PNmat,*Mesh);
    id++;

    Topol[0] = 49;
    Topol[1] = 47;
    Topol[2] = 48;
    new TPZGeoElRefPattern<TPZArc3D> (id,Topol,PNmat,*Mesh);
    id++;

    Topol[0] = 51;
    Topol[1] = 49;
    Topol[2] = 50;
    new TPZGeoElRefPattern<TPZArc3D> (id,Topol,PNmat,*Mesh);
    id++;

    Topol[0] = 53;
    Topol[1] = 51;
    Topol[2] = 52;
    new TPZGeoElRefPattern<TPZArc3D> (id,Topol,PNmat,*Mesh);
    id++;

#else

    Topol.Resize(2);

    Topol[0] = 31;
    Topol[1] = 29;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,PNmat,*Mesh);
    id++;

    Topol[0] = 33;
    Topol[1] = 31;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,PNmat,*Mesh);
    id++;

    Topol[0] = 35;
    Topol[1] = 33;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,PNmat,*Mesh);
    id++;

    Topol[0] = 37;
    Topol[1] = 35;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,PNmat,*Mesh);
    id++;

    Topol[0] = 39;
    Topol[1] = 37;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,PNmat,*Mesh);
    id++;

    Topol[0] = 41;
    Topol[1] = 39;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,PNmat,*Mesh);
    id++;

    Topol[0] = 43;
    Topol[1] = 41;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,PNmat,*Mesh);
    id++;

    Topol[0] = 45;
    Topol[1] = 43;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,PNmat,*Mesh);
    id++;

    Topol[0] = 47;
    Topol[1] = 45;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,PNmat,*Mesh);
    id++;

    Topol[0] = 49;
    Topol[1] = 47;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,PNmat,*Mesh);
    id++;

    Topol[0] = 51;
    Topol[1] = 49;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,PNmat,*Mesh);
    id++;

    Topol[0] = 53;
    Topol[1] = 51;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,PNmat,*Mesh);
    id++;

#endif

    Topol.Resize(2);

    Topol[0] = 0;
    Topol[1] = 1;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,Basemat,*Mesh);
    id++;

    Topol[0] = 1;
    Topol[1] = 2;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,Basemat,*Mesh);
    id++;

    Topol[0] = 2;
    Topol[1] = 3;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,Nonemat,*Mesh);
    id++;

    Topol[0] = 3;
    Topol[1] = 4;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,Nonemat,*Mesh);
    id++;

    Topol[0] = 4;
    Topol[1] = 5;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,Nonemat,*Mesh);
    id++;

    Topol[0] = 5;
    Topol[1] = 6;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,TTmat,*Mesh);
    id++;

    Topol[0] = 6;
    Topol[1] = 7;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,Nonemat,*Mesh);
    id++;

    Topol[0] = 7;
    Topol[1] = 8;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,Nonemat,*Mesh);
    id++;

    Topol[0] = 8;
    Topol[1] = 9;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,Nonemat,*Mesh);
    id++;

    Topol[0] = 9;
    Topol[1] = 23;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,Nonemat,*Mesh);
    id++;

    Topol[0] = 23;
    Topol[1] = 57;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,Nonemat,*Mesh);
    id++;

    Topol[0] = 57;
    Topol[1] = 26;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,Nonemat,*Mesh);
    id++;

    Topol[0] = 26;
    Topol[1] = 10;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,Nonemat,*Mesh);
    id++;

    Topol[0] = 10;
    Topol[1] = 11;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,Nonemat,*Mesh);
    id++;

    Topol[0] = 11;
    Topol[1] = 12;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,Nonemat,*Mesh);
    id++;

    Topol[0] = 12;
    Topol[1] = 13;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,Nonemat,*Mesh);
    id++;

    Topol[0] = 13;
    Topol[1] = 14;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,Nonemat,*Mesh);
    id++;

    Topol[0] = 14;
    Topol[1] = 15;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,FFmat,*Mesh);
    id++;

    Topol[0] = 15;
    Topol[1] = 16;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,Nonemat,*Mesh);
    id++;

    Topol[0] = 16;
    Topol[1] = 17;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,FEmat,*Mesh);
    id++;

    Topol[0] = 17;
    Topol[1] = 18;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,Nonemat,*Mesh);
    id++;

    Topol[0] = 18;
    Topol[1] = 19;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,PGmat,*Mesh);
    id++;

    Topol[0] = 19;
    Topol[1] = 20;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,Nonemat,*Mesh);
    id++;

    Topol[0] = 20;
    Topol[1] = 21;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,Nonemat,*Mesh);
    id++;

    Topol[0] = 21;
    Topol[1] = 22;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,Nonemat,*Mesh);
    id++;

    Topol[0] = 22;
    Topol[1] = 0;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,Nonemat,*Mesh);
    id++;

    Topol[0] = 29;
    Topol[1] = 25;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,Nonemat,*Mesh);
    id++;

    Topol[0] = 25;
    Topol[1] = 54;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,Nonemat,*Mesh);
    id++;

    Topol[0] = 54;
    Topol[1] = 24;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,Nonemat,*Mesh);
    id++;

    Topol[0] = 24;
    Topol[1] = 53;
    new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,Nonemat,*Mesh);
    id++;

//     Topol[0] = 53;
//     Topol[1] = 9;
//     new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,PDmat,*Mesh);
//     id++;

    Topol.Resize(1);
    Topol[0] = 0;
    new TPZGeoElRefPattern<TPZGeoPoint> (id,Topol,Dotmat,*Mesh);

    Mesh->BuildBlendConnectivity();


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

#include "pzgeoelrefless.h.h"
///CreateGeoElement -> TPZArc3D
template< >
TPZGeoEl *TPZGeoElRefLess<TPZArc3D >::CreateGeoElement(MElementType type, TPZVec<int>& nodeindexes, int matid, int& index)
{
    TPZGeoMesh &mesh = *(this->Mesh());
    if(!&mesh) return 0;
    return CreateGeoElementMapped(mesh,type,nodeindexes,matid,index);
}

    #define TPZGEOELEMENTARC3DID 300
    template<>
    int TPZGeoElRefPattern<TPZArc3D>::ClassId() const {
        return TPZGEOELEMENTARC3DID;
    }
    template class 
    TPZRestoreClass< TPZGeoElRefPattern<TPZArc3D>, TPZGEOELEMENTARC3DID>;

    template<>
    TPZCompEl *(*TPZGeoElRefLess<TPZArc3D>::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreateLinearEl;

    template class TPZGeoElRefLess<TPZArc3D>;


///CreateGeoElement -> TPZEllipse
template< >
TPZGeoEl *TPZGeoElRefLess<TPZEllipse >::CreateGeoElement(MElementType type, TPZVec<int>& nodeindexes, int matid, int& index)
{
    TPZGeoMesh &mesh = *(this->Mesh());
    if(!&mesh) return 0;
    return CreateGeoElementMapped(mesh,type,nodeindexes,matid,index);
}

    #define TPZGEOELEMENTELLIPSEID 301
    template<>
    int TPZGeoElRefPattern<TPZEllipse>::ClassId() const
    {
        return TPZGEOELEMENTELLIPSEID;
    }
    template class 
    TPZRestoreClass< TPZGeoElRefPattern<TPZEllipse>, TPZGEOELEMENTELLIPSEID>;

    template<>
    TPZCompEl *(*TPZGeoElRefLess<TPZEllipse>::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreateLinearEl;

    template class TPZGeoElRefLess<TPZEllipse>;


///CreateGeoElement -> TPZCurvedTriangle
template< >
TPZGeoEl *TPZGeoElRefLess<TPZCurvedTriangle >::CreateGeoElement(MElementType type, TPZVec<int>& nodeindexes, int matid, int& index)
{
    TPZGeoMesh &mesh = *(this->Mesh());
    if(!&mesh) return 0;
    return CreateGeoElementMapped(mesh,type,nodeindexes,matid,index);
}

    #define TPZGEOELEMENTCURVEDTRIANGLEID 302
    template<>
    int TPZGeoElRefPattern<TPZCurvedTriangle>::ClassId() const {
        return TPZGEOELEMENTCURVEDTRIANGLEID;
    }
    template class 
    TPZRestoreClass< TPZGeoElRefPattern<TPZCurvedTriangle>, TPZGEOELEMENTCURVEDTRIANGLEID>;

    template<>
    TPZCompEl *(*TPZGeoElRefLess<TPZCurvedTriangle>::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreateTriangleEl;

    template class TPZGeoElRefLess<TPZCurvedTriangle>;



///CreateGeoElement -> TPZGeoBlend
#define IMPLEMENTBLEND(TGEO,CLASSID,CREATEFUNCTION) \
template< > \
TPZGeoEl *TPZGeoElRefLess<TPZGeoBlend<TGEO> >::CreateGeoElement(MElementType type, TPZVec<int>& nodeindexes, int matid, int& index) \
{ \
    TPZGeoMesh &mesh = *(this->Mesh()); \
    if(!&mesh) return 0; \
    return CreateGeoElementMapped(mesh,type,nodeindexes,matid,index); \
} \
\
    template<> \
    int TPZGeoElRefPattern<TPZGeoBlend<TGEO>  >::ClassId() const { \
        return CLASSID; \
    } \
    template class \
    TPZRestoreClass< TPZGeoElRefPattern<TPZGeoBlend<TGEO> >, CLASSID>; \
\
    template<> \
    TPZCompEl *(*TPZGeoElRefLess<TPZGeoBlend<TGEO> >::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CREATEFUNCTION; \
\
    template class TPZGeoElRefLess<TPZGeoBlend<TGEO> >;

#define TPZGEOBLENDPOINTID 303
#define TPZGEOBLENDLINEARID 304
#define TPZGEOBLENDQUADID 305
#define TPZGEOBLENDTRIANGLEID 306
#define TPZGEOBLENDCUBEID 307
#define TPZGEOBLENDPRISMID 308
#define TPZGEOBLENDPYRAMIDID 309
#define TPZGEOBLENDTETRAHEDRAID 310

IMPLEMENTBLEND(pzgeom::TPZGeoPoint,TPZGEOBLENDPOINTID,CreatePointEl)
IMPLEMENTBLEND(pzgeom::TPZGeoLinear,TPZGEOBLENDLINEARID,CreateLinearEl)
IMPLEMENTBLEND(pzgeom::TPZGeoQuad,TPZGEOBLENDQUADID,CreateQuadEl)
IMPLEMENTBLEND(pzgeom::TPZGeoTriangle,TPZGEOBLENDTRIANGLEID,CreateTriangleEl)
IMPLEMENTBLEND(pzgeom::TPZGeoCube,TPZGEOBLENDCUBEID,CreateCubeEl)
IMPLEMENTBLEND(pzgeom::TPZGeoPrism,TPZGEOBLENDPRISMID,CreatePrismEl)
IMPLEMENTBLEND(pzgeom::TPZGeoPyramid,TPZGEOBLENDPYRAMIDID,CreatePyramEl)
IMPLEMENTBLEND(pzgeom::TPZGeoTetrahedra,TPZGEOBLENDTETRAHEDRAID,CreateTetraEl)


///CreateGeoElement -> TPZQuadraticQuad
template< >
TPZGeoEl *TPZGeoElRefLess<TPZQuadraticQuad >::CreateGeoElement(MElementType type, TPZVec<int>& nodeindexes, int matid, int& index)
{
    TPZGeoMesh &mesh = *(this->Mesh());
    if(!&mesh) return 0;
    return CreateGeoElementMapped(mesh,type,nodeindexes,matid,index);
}

    #define TPZGEOELEMENTQUADRATICQUADID 311
    template<>
    int TPZGeoElRefPattern<TPZQuadraticQuad>::ClassId() const {
        return TPZGEOELEMENTQUADRATICQUADID;
    }
    template class 
    TPZRestoreClass< TPZGeoElRefPattern<TPZQuadraticQuad>, TPZGEOELEMENTQUADRATICQUADID>;

    template<>
    TPZCompEl *(*TPZGeoElRefLess<TPZQuadraticQuad>::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreateQuadEl;

    template class TPZGeoElRefLess<TPZQuadraticQuad>;


///CreateGeoElement -> TPZQuadraticTetra
template< >
TPZGeoEl *TPZGeoElRefLess<TPZQuadraticTetra >::CreateGeoElement(MElementType type, TPZVec<int>& nodeindexes, int matid, int& index)
{
    TPZGeoMesh &mesh = *(this->Mesh());
    if(!&mesh) return 0;
    return CreateGeoElementMapped(mesh,type,nodeindexes,matid,index);
}

    #define TPZGEOELEMENTQUADRATICTETRAID 312
    template<>
    int TPZGeoElRefPattern<TPZQuadraticTetra>::ClassId() const {
        return TPZGEOELEMENTQUADRATICTETRAID;
    }
    template class 
    TPZRestoreClass< TPZGeoElRefPattern<TPZQuadraticTetra>, TPZGEOELEMENTQUADRATICTETRAID>;

    template<>
    TPZCompEl *(*TPZGeoElRefLess<TPZQuadraticTetra>::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreateTetraEl;

    template class TPZGeoElRefLess<TPZQuadraticTetra>;


///CreateGeoElement -> TPZQuadraticTrig
template< >
TPZGeoEl *TPZGeoElRefLess<TPZQuadraticTrig >::CreateGeoElement(MElementType type, TPZVec<int>& nodeindexes, int matid, int& index)
{
    TPZGeoMesh &mesh = *(this->Mesh());
    if(!&mesh) return 0;
    return CreateGeoElementMapped(mesh,type,nodeindexes,matid,index);
}

    #define TPZGEOELEMENTQUADRATICTRIANGLEID 313
    template<>
    int TPZGeoElRefPattern<TPZQuadraticTrig>::ClassId() const {
        return TPZGEOELEMENTQUADRATICTRIANGLEID;
    }
    template class 
    TPZRestoreClass< TPZGeoElRefPattern<TPZQuadraticTrig>, TPZGEOELEMENTQUADRATICTRIANGLEID>;

    template<>
    TPZCompEl *(*TPZGeoElRefLess<TPZQuadraticTrig>::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreateTriangleEl;

    template class TPZGeoElRefLess<TPZQuadraticTrig>;


#include "pznoderep.h.h"
template class pzgeom::TPZNodeRep<2,TPZEllipse>;
template class pzgeom::TPZNodeRep<3,TPZArc3D>;
template class pzgeom::TPZNodeRep<3,TPZCurvedTriangle>;
template class pzgeom::TPZNodeRep<8,TPZGeoBlend<TPZGeoCube> >;
template class pzgeom::TPZNodeRep<6,TPZGeoBlend<TPZGeoPrism> >;
template class pzgeom::TPZNodeRep<8,TPZQuadraticQuad>;
template class pzgeom::TPZNodeRep<10,TPZQuadraticTetra>;
template class pzgeom::TPZNodeRep<6,TPZQuadraticTrig>;



#include <pzgeoel.h>


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

 /**
  * Method that generate an 3D Horizontal Well_Mesh - Elliptical Reservoir
  * Caju - 2007
  */
TPZGeoMesh * HorizWellGenE(double Wr, double Hr, double Lr, double Lw, double Dw, int part, double excentr)
{
    TPZGeoMesh * Mesh = new TPZGeoMesh;
    int id = 0;

    #ifndef Debug
    if(Dw > Wr || Dw > Hr)
    {
        cout << "Well Diameter (Dw) Bigger than Reservoir Dimensions!\n";
        exit(-1);
    }
    #endif

    double Yneg    = -Lr/2.;
    double Ya      = -Lw/2.;
    double Yb      =  Lw/2.;
    double Ypos    =  Lr/2.;
    double h       =  Lw/double(part);

    int Qnodes = 48 + (part-1)*16;
    TPZVec < TPZVec <REAL> > NodesCoords(Qnodes);
    for(int i = 0; i < Qnodes; i++) NodesCoords[i].Resize(3);

    for(int i = 0; i < part + 1; i++)
    {
        NodesCoords[0 + i*16][0] = -Wr/2.;
        NodesCoords[0 + i*16][1] =  Ya + i*h;
        NodesCoords[0 + i*16][2] = -Hr/2.;

        NodesCoords[1 + i*16][0] =  Wr/2.;
        NodesCoords[1 + i*16][1] =  Ya + i*h;
        NodesCoords[1 + i*16][2] = -Hr/2.;

        NodesCoords[2 + i*16][0] =  Wr/2.;
        NodesCoords[2 + i*16][1] =  Ya + i*h;
        NodesCoords[2 + i*16][2] =  Hr/2.;

        NodesCoords[3 + i*16][0] = -Wr/2.;
        NodesCoords[3 + i*16][1] =  Ya + i*h;
        NodesCoords[3 + i*16][2] =  Hr/2.;

        NodesCoords[4 + i*16][0] = -Hr/2.;
        NodesCoords[4 + i*16][1] =  Ya + i*h;
        NodesCoords[4 + i*16][2] = -Hr/2.;

        NodesCoords[5 + i*16][0] =  Hr/2.;
        NodesCoords[5 + i*16][1] =  Ya + i*h;
        NodesCoords[5 + i*16][2] = -Hr/2.;

        NodesCoords[6 + i*16][0] =  Hr/2.;
        NodesCoords[6 + i*16][1] =  Ya + i*h;
        NodesCoords[6 + i*16][2] =  Hr/2.;

        NodesCoords[7 + i*16][0] = -Hr/2.;
        NodesCoords[7 + i*16][1] =  Ya + i*h;
        NodesCoords[7 + i*16][2] =  Hr/2.;

        NodesCoords[8 + i*16][0] = -Dw*sqrt(2.)/4.;
        NodesCoords[8 + i*16][1] =  Ya + i*h;
        NodesCoords[8 + i*16][2] = -Dw*sqrt(2.)/4. + excentr;

        NodesCoords[9 + i*16][0] =  0.;
        NodesCoords[9 + i*16][1] =  Ya + i*h;
        NodesCoords[9 + i*16][2] = -Dw/2. + excentr;

        NodesCoords[10 + i*16][0] =  Dw*sqrt(2.)/4.;
        NodesCoords[10 + i*16][1] =  Ya + i*h;
        NodesCoords[10 + i*16][2] = -Dw*sqrt(2.)/4. + excentr;

        NodesCoords[11 + i*16][0] =  Dw/2.;
        NodesCoords[11 + i*16][1] =  Ya + i*h;
        NodesCoords[11 + i*16][2] =  0. + excentr;

        NodesCoords[12 + i*16][0] =  Dw*sqrt(2.)/4.;
        NodesCoords[12 + i*16][1] =  Ya + i*h;
        NodesCoords[12 + i*16][2] =  Dw*sqrt(2.)/4. + excentr;

        NodesCoords[13 + i*16][0] =  0.;
        NodesCoords[13 + i*16][1] =  Ya + i*h;
        NodesCoords[13 + i*16][2] =  Dw/2. + excentr;

        NodesCoords[14 + i*16][0] = -Dw*sqrt(2.)/4.;
        NodesCoords[14 + i*16][1] =  Ya + i*h;
        NodesCoords[14 + i*16][2] =  Dw*sqrt(2.)/4. + excentr;

        NodesCoords[15 + i*16][0] = -Dw/2.;
        NodesCoords[15 + i*16][1] =  Ya + i*h;
        NodesCoords[15 + i*16][2] =  0. + excentr;

        if(i == 0)
        {
            NodesCoords[0 + i*16][0] = -Wr/2.;
            NodesCoords[0 + i*16][1] =  Ya - Hr/2.;
            NodesCoords[0 + i*16][2] = -Hr/2.;

            NodesCoords[1 + i*16][0] =  Wr/2.;
            NodesCoords[1 + i*16][1] =  Ya - Hr/2.;
            NodesCoords[1 + i*16][2] = -Hr/2.;

            NodesCoords[2 + i*16][0] =  Wr/2.;
            NodesCoords[2 + i*16][1] =  Ya - Hr/2.;
            NodesCoords[2 + i*16][2] =  Hr/2.;

            NodesCoords[3 + i*16][0] = -Wr/2.;
            NodesCoords[3 + i*16][1] =  Ya - Hr/2.;
            NodesCoords[3 + i*16][2] =  Hr/2.;

            NodesCoords[4 + i*16][0] = -Hr/2.;
            NodesCoords[4 + i*16][1] =  Ya - Hr/2.;
            NodesCoords[4 + i*16][2] = -Hr/2.;

            NodesCoords[4 + i*16][0] = -Hr/2.;
            NodesCoords[4 + i*16][1] =  Ya - Hr/2.;
            NodesCoords[4 + i*16][2] = -Hr/2.;

            NodesCoords[5 + i*16][0] =  Hr/2.;
            NodesCoords[5 + i*16][1] =  Ya - Hr/2.;
            NodesCoords[5 + i*16][2] = -Hr/2.;

            NodesCoords[6 + i*16][0] =  Hr/2.;
            NodesCoords[6 + i*16][1] =  Ya - Hr/2.;
            NodesCoords[6 + i*16][2] =  Hr/2.;

            NodesCoords[7 + i*16][0] = -Hr/2.;
            NodesCoords[7 + i*16][1] =  Ya - Hr/2.;
            NodesCoords[7 + i*16][2] =  Hr/2.;
        }

        if(i == part)
        {
            NodesCoords[0 + i*16][0] = -Wr/2.;
            NodesCoords[0 + i*16][1] =  Yb + Hr/2.;
            NodesCoords[0 + i*16][2] = -Hr/2.;

            NodesCoords[1 + i*16][0] =  Wr/2.;
            NodesCoords[1 + i*16][1] =  Yb + Hr/2.;
            NodesCoords[1 + i*16][2] = -Hr/2.;

            NodesCoords[2 + i*16][0] =  Wr/2.;
            NodesCoords[2 + i*16][1] =  Yb + Hr/2.;
            NodesCoords[2 + i*16][2] =  Hr/2.;

            NodesCoords[3 + i*16][0] = -Wr/2.;
            NodesCoords[3 + i*16][1] =  Yb + Hr/2.;
            NodesCoords[3 + i*16][2] =  Hr/2.;

            NodesCoords[4 + i*16][0] = -Hr/2.;
            NodesCoords[4 + i*16][1] =  Yb + Hr/2.;
            NodesCoords[4 + i*16][2] = -Hr/2.;

            NodesCoords[5 + i*16][0] =  Hr/2.;
            NodesCoords[5 + i*16][1] =  Yb + Hr/2.;
            NodesCoords[5 + i*16][2] = -Hr/2.;

            NodesCoords[6 + i*16][0] =  Hr/2.;
            NodesCoords[6 + i*16][1] =  Yb + Hr/2.;
            NodesCoords[6 + i*16][2] =  Hr/2.;

            NodesCoords[7 + i*16][0] = -Hr/2.;
            NodesCoords[7 + i*16][1] =  Yb + Hr/2.;
            NodesCoords[7 + i*16][2] =  Hr/2.;
        }
    }
    NodesCoords[Qnodes-1][0] = -Wr/2.;
    NodesCoords[Qnodes-1][1] =  Ypos;
    NodesCoords[Qnodes-1][2] =  Hr/2.;

    NodesCoords[Qnodes-2][0] = -Hr/2.;
    NodesCoords[Qnodes-2][1] =  Ypos;
    NodesCoords[Qnodes-2][2] =  Hr/2.;

    NodesCoords[Qnodes-3][0] =  Hr/2.;
    NodesCoords[Qnodes-3][1] =  Ypos;
    NodesCoords[Qnodes-3][2] =  Hr/2.;

    NodesCoords[Qnodes-4][0] =  Wr/2.;
    NodesCoords[Qnodes-4][1] =  Ypos;
    NodesCoords[Qnodes-4][2] =  Hr/2.;

    NodesCoords[Qnodes-5][0] =  Wr/2.;
    NodesCoords[Qnodes-5][1] =  Ypos;
    NodesCoords[Qnodes-5][2] = -Hr/2.;

    NodesCoords[Qnodes-6][0] =  Hr/2.;
    NodesCoords[Qnodes-6][1] =  Ypos;
    NodesCoords[Qnodes-6][2] = -Hr/2.;

    NodesCoords[Qnodes-7][0] = -Hr/2.;
    NodesCoords[Qnodes-7][1] =  Ypos;
    NodesCoords[Qnodes-7][2] = -Hr/2.;

    NodesCoords[Qnodes-8][0] = -Wr/2.;
    NodesCoords[Qnodes-8][1] =  Ypos;
    NodesCoords[Qnodes-8][2] = -Hr/2.;

    NodesCoords[Qnodes-9][0] = -Wr/2.;
    NodesCoords[Qnodes-9][1] =  Yneg;
    NodesCoords[Qnodes-9][2] =  Hr/2.;

    NodesCoords[Qnodes-10][0] = -Hr/2.;
    NodesCoords[Qnodes-10][1] =  Yneg;
    NodesCoords[Qnodes-10][2] =  Hr/2.;

    NodesCoords[Qnodes-11][0] =  Hr/2.;
    NodesCoords[Qnodes-11][1] =  Yneg;
    NodesCoords[Qnodes-11][2] =  Hr/2.;

    NodesCoords[Qnodes-12][0] =  Wr/2.;
    NodesCoords[Qnodes-12][1] =  Yneg;
    NodesCoords[Qnodes-12][2] =  Hr/2.;

    NodesCoords[Qnodes-13][0] =  Wr/2.;
    NodesCoords[Qnodes-13][1] =  Yneg;
    NodesCoords[Qnodes-13][2] = -Hr/2.;

    NodesCoords[Qnodes-14][0] =  Hr/2.;
    NodesCoords[Qnodes-14][1] =  Yneg;
    NodesCoords[Qnodes-14][2] = -Hr/2.;

    NodesCoords[Qnodes-15][0] = -Hr/2.;
    NodesCoords[Qnodes-15][1] =  Yneg;
    NodesCoords[Qnodes-15][2] = -Hr/2.;

    NodesCoords[Qnodes-16][0] = -Wr/2.;
    NodesCoords[Qnodes-16][1] =  Yneg;
    NodesCoords[Qnodes-16][2] = -Hr/2.;

    Mesh->NodeVec().Resize(Qnodes);
    TPZVec <TPZGeoNode> Node(Qnodes);
    for(int n = 0; n < Qnodes; n++)
    {
        Node[n].SetNodeId(n);
        Node[n].SetCoord(&NodesCoords[n][0]);
        Mesh->NodeVec()[n] = Node[n];
    }

    TPZVec <int> Topol(8);


    /// *********** CHOOSE MATERIAL FOR EACH ELEMENT ***********
    int reservMat = 1;
    int wellMat = 2;
    int arcMat = 10;


    ///Horizontal Well
    TPZVec< TPZGeoElRefLess<TPZGeoBlend<pzgeom::TPZGeoCube> > *> WellEl(part);
    for(int w = 0; w < part; w++)
    {
        Topol[0] =  8 + w*16; Topol[1] = 10 + w*16; Topol[2] = 26 + w*16; Topol[3] = 24 + w*16;
        Topol[4] = 14 + w*16; Topol[5] = 12 + w*16; Topol[6] = 28 + w*16; Topol[7] = 30 + w*16;
        WellEl[w] = new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoCube> > (id,Topol,wellMat,*Mesh);
        id++;
    }

    ///Arcs in Well x_z Planes
    Topol.Resize(3);
    TPZVec <TPZGeoEl*> ArcEl( 4*(part + 1) );
    for(int interface = 0; interface < part + 1; interface++)
    {
        for(int edge = 0; edge < 4; edge++)
        {
            if(edge == 3)
            {
                Topol[0] = 14 + interface*16;
                Topol[1] = Topol[0] - 6;
                Topol[2] = Topol[0] + 1;
            }
            else
            {
                Topol[0] = 8 + 2*edge + interface*16;
                Topol[1] = Topol[0] + 2;
                Topol[2] = Topol[1] - 1;
            }
            ArcEl[4*interface + edge] = new TPZGeoElRefPattern<TPZArc3D> (id,Topol,arcMat,*Mesh);
            id++;
        }
    }

    ///Hexaedrons that Surround Horizontal Well
    Topol.Resize(8);
    TPZVec < TPZGeoElRefLess<TPZGeoBlend<pzgeom::TPZGeoCube> >* > HexaEl(4*part);
    for(int w = 0; w < part; w++)
    {
        for(int face = 0; face < 4; face++)
        {
            if(face == 3)
            {
                Topol[0] = 7 + w*16;   Topol[1] = Topol[0] - 3; Topol[2] = Topol[0] + 13;  Topol[3] = 23 + w*16;
                Topol[4] = 14 + w*16; Topol[5] = Topol[0] + 1; Topol[6] = Topol[0] + 17; Topol[7] = 30 + w*16;
            }
            else
            {
                Topol[0] = 4 + face + w*16;   Topol[1] =  5 + face + w*16;   Topol[2] = 21 + face + w*16;   Topol[3] = 20 + face + w*16;
                Topol[4] = 8 + 2*face + w*16; Topol[5] = 10 + 2*face + w*16; Topol[6] = 26 + 2*face + w*16; Topol[7] = 24 + 2*face + w*16;
            }
            HexaEl[4*w + face] = new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoCube> >  (id, Topol,reservMat,*Mesh);
            id++;
        }
    }

    ///Lateral Hexaedrons
    TPZVec < TPZGeoElRefLess<TPZGeoBlend<pzgeom::TPZGeoCube> >* > HexaElLateral(2*part);
    for(int w = 0; w < part; w++)
    {
        Topol[0] = 0 + 16*w; Topol[1] = 4 + 16*w; Topol[2] = 20 + 16*w; Topol[3] = 16 + 16*w;
        Topol[4] = 3 + 16*w; Topol[5] = 7 + 16*w; Topol[6] = 23 + 16*w; Topol[7] = 19 + 16*w;
        HexaElLateral[2*w + 0] = new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoCube> >  (id, Topol,reservMat,*Mesh);
        id++;

        Topol[0] = 5 + 16*w; Topol[1] = 1 + 16*w; Topol[2] = 17 + 16*w; Topol[3] = 21 + 16*w;
        Topol[4] = 6 + 16*w; Topol[5] = 2 + 16*w; Topol[6] = 18 + 16*w; Topol[7] = 22 + 16*w;
        HexaElLateral[2*w + 1] = new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoCube> >  (id, Topol,reservMat,*Mesh);
        id++;
    }

    ///Hexahedron WellTips
    TPZVec < TPZGeoElRefLess<TPZGeoBlend<pzgeom::TPZGeoCube> >* > HexaElTips(2);
    {
        Topol[0] = 4; Topol[1] = 5; Topol[2] = 10; Topol[3] = 8;
        Topol[4] = 7; Topol[5] = 6; Topol[6] = 12; Topol[7] = 14;
        HexaElTips[0] = new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoCube> >  (id, Topol,reservMat,*Mesh);
        id++;

        Topol[0] = Qnodes - 24; Topol[1] = Qnodes - 22; Topol[2] = Qnodes - 27; Topol[3] = Qnodes - 28;
        Topol[4] = Qnodes - 18; Topol[5] = Qnodes - 20; Topol[6] = Qnodes - 26; Topol[7] = Qnodes - 25;
        HexaElTips[1] = new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoCube> >  (id, Topol,reservMat,*Mesh);
        id++;
    }

    ///Hexahedron ReservoirTips
    TPZVec < TPZGeoElRefLess<TPZGeoBlend<pzgeom::TPZGeoCube> >* > HexaResTips(2);
    {
        Topol[0] = Qnodes-15; Topol[1] = Qnodes-14; Topol[2] = 5; Topol[3] = 4;
        Topol[4] = Qnodes-10; Topol[5] = Qnodes-11; Topol[6] = 6; Topol[7] = 7;
        HexaResTips[0] = new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoCube> >  (id, Topol,reservMat,*Mesh);
        id++;

        Topol[0] = Qnodes-28; Topol[1] = Qnodes-27; Topol[2] = Qnodes-6; Topol[3] = Qnodes-7;
        Topol[4] = Qnodes-25; Topol[5] = Qnodes-26; Topol[6] = Qnodes-3; Topol[7] = Qnodes-2;
        HexaResTips[1] = new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoCube> >  (id, Topol,reservMat,*Mesh);
        id++;
    }

    ///Prism ReservoirTips
    Topol.Resize(6);
    TPZVec < TPZGeoElRefLess<TPZGeoBlend<pzgeom::TPZGeoPrism> >* > PrismResTips(8);
    {
        Topol[0] = Qnodes-10; Topol[1] = Qnodes- 9; Topol[2] = 7;
        Topol[3] = Qnodes-15; Topol[4] = Qnodes-16; Topol[5] = 4;
        PrismResTips[0] = new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoPrism> >  (id, Topol,reservMat,*Mesh);
        id++;

        Topol[0] = 0; Topol[1] = Qnodes-16; Topol[2] = 4;
        Topol[3] = 3; Topol[4] = Qnodes- 9; Topol[5] = 7;
        PrismResTips[1] = new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoPrism> >  (id, Topol,reservMat,*Mesh);
        id++;

        Topol[0] = Qnodes-14; Topol[1] = Qnodes-13; Topol[2] = 5;
        Topol[3] = Qnodes-11; Topol[4] = Qnodes-12; Topol[5] = 6;
        PrismResTips[2] = new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoPrism> >  (id, Topol,reservMat,*Mesh);
        id++;

        Topol[0] = 2; Topol[1] = Qnodes-12; Topol[2] = 6;
        Topol[3] = 1; Topol[4] = Qnodes-13; Topol[5] = 5;
        PrismResTips[3] = new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoPrism> >  (id, Topol,reservMat,*Mesh);
        id++;

        Topol[0] = Qnodes-29; Topol[1] = Qnodes-1; Topol[2] = Qnodes-25;
        Topol[3] = Qnodes-32; Topol[4] = Qnodes-8; Topol[5] = Qnodes-28;
        PrismResTips[4] = new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoPrism> >  (id, Topol,reservMat,*Mesh);
        id++;

        Topol[0] = Qnodes-7; Topol[1] = Qnodes-8; Topol[2] = Qnodes-28;
        Topol[3] = Qnodes-2; Topol[4] = Qnodes-1; Topol[5] = Qnodes-25;
        PrismResTips[5] = new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoPrism> >  (id, Topol,reservMat,*Mesh);
        id++;

        Topol[0] = Qnodes-31; Topol[1] = Qnodes-5; Topol[2] = Qnodes-27;
        Topol[3] = Qnodes-30; Topol[4] = Qnodes-4; Topol[5] = Qnodes-26;
        PrismResTips[6] = new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoPrism> >  (id, Topol,reservMat,*Mesh);
        id++;

        Topol[0] = Qnodes-3; Topol[1] = Qnodes-4; Topol[2] = Qnodes-26;
        Topol[3] = Qnodes-6; Topol[4] = Qnodes-5; Topol[5] = Qnodes-27;
        PrismResTips[7] = new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoPrism> >  (id, Topol,reservMat,*Mesh);
        id++;
    }

    Topol.Resize(2);
    vector<int> node;
    NodesHunter(*Mesh,node,Qnodes-16,Qnodes-13);
    for(unsigned int n = 0; n < node.size() - 1; n++)
    {
        Topol[0] = node[n];
        Topol[1] = node[n+1];
        TPZGeoElRefPattern<TPZEllipse> *Geo = new TPZGeoElRefPattern<TPZEllipse> (id,Topol,arcMat,*Mesh);
        Geo->Geom().SetAxes(Wr/2.,Lr/2.);
        id++;
    }
    NodesHunter(*Mesh,node,Qnodes-13,Qnodes-5);
    for(int n = 0; n < node.size() - 1; n++)
    {
        Topol[0] = node[n];
        Topol[1] = node[n+1];
        TPZGeoElRefPattern<TPZEllipse> *Geo = new TPZGeoElRefPattern<TPZEllipse> (id,Topol,arcMat,*Mesh);
        Geo->Geom().SetAxes(Wr/2.,Lr/2.);
        id++;
    }
    NodesHunter(*Mesh,node,Qnodes-5,Qnodes-8);
    for(int n = 0; n < node.size() - 1; n++)
    {
        Topol[0] = node[n];
        Topol[1] = node[n+1];
        TPZGeoElRefPattern<TPZEllipse> *Geo = new TPZGeoElRefPattern<TPZEllipse> (id,Topol,arcMat,*Mesh);
        Geo->Geom().SetAxes(Wr/2.,Lr/2.);
        id++;
    }
    NodesHunter(*Mesh,node,Qnodes-8,Qnodes-16);
    for(int n = 0; n < node.size() - 1; n++)
    {
        Topol[0] = node[n];
        Topol[1] = node[n+1];
        TPZGeoElRefPattern<TPZEllipse> *Geo = new TPZGeoElRefPattern<TPZEllipse> (id,Topol,arcMat,*Mesh);
        Geo->Geom().SetAxes(Wr/2.,Lr/2.);
        id++;
    }
    NodesHunter(*Mesh,node,Qnodes-9,Qnodes-12);
    for(int n = 0; n < node.size() - 1; n++)
    {
        Topol[0] = node[n];
        Topol[1] = node[n+1];
        TPZGeoElRefPattern<TPZEllipse> *Geo = new TPZGeoElRefPattern<TPZEllipse> (id,Topol,arcMat,*Mesh);
        Geo->Geom().SetAxes(Wr/2.,Lr/2.);
        id++;
    }
    NodesHunter(*Mesh,node,Qnodes-12,Qnodes-4);
    for(int n = 0; n < node.size() - 1; n++)
    {
        Topol[0] = node[n];
        Topol[1] = node[n+1];
        TPZGeoElRefPattern<TPZEllipse> *Geo = new TPZGeoElRefPattern<TPZEllipse> (id,Topol,arcMat,*Mesh);
        Geo->Geom().SetAxes(Wr/2.,Lr/2.);
        id++;
    }
    NodesHunter(*Mesh,node,Qnodes-4,Qnodes-1);
    for(int n = 0; n < node.size() - 1; n++)
    {
        Topol[0] = node[n];
        Topol[1] = node[n+1];
        TPZGeoElRefPattern<TPZEllipse> *Geo = new TPZGeoElRefPattern<TPZEllipse> (id,Topol,arcMat,*Mesh);
        Geo->Geom().SetAxes(Wr/2.,Lr/2.);
        id++;
    }
    NodesHunter(*Mesh,node,Qnodes-1,Qnodes-9);
    for(int n = 0; n < node.size() - 1; n++)
    {
        Topol[0] = node[n];
        Topol[1] = node[n+1];
        TPZGeoElRefPattern<TPZEllipse> *Geo = new TPZGeoElRefPattern<TPZEllipse> (id,Topol,arcMat,*Mesh);
        Geo->Geom().SetAxes(Wr/2.,Lr/2.);
        id++;
    }

    Mesh->BuildBlendConnectivity();

    return Mesh;
}

void NodesHunter(TPZGeoMesh &gMesh, std::vector<int>& NodesHunted, int IdIni, int IdFin, double Tol)
{
    /** Although this method considers coordinates in R3, Coord-Z must be constant for all nodes */
    /** i.e., nodes coordinates belong to XY parallel plane */
    int dim = 3;

    int VecSize = gMesh.NodeVec().NElements();
    int posIni = -1;
    int posFin = -1;

    /** Changing NodesCoords from Cartesian Notation to Vectorial Notation */
    /** with respect to IniNode coordinate */
    TPZFMatrix VectorialNotation(VecSize,dim);
    for(int pos = 0; pos < VecSize; pos++)
    {
        if(gMesh.NodeVec()[pos].Id() == IdIni) posIni = pos;
        if(gMesh.NodeVec()[pos].Id() == IdFin) posFin = pos;

        for(int coord = 0; coord < dim; coord++)
        {
            double val = gMesh.NodeVec()[pos].Coord(coord) - gMesh.NodeVec()[IdIni].Coord(coord);
            VectorialNotation.PutVal(pos,coord,val);
        }
    }

    if(posIni == -1 || posFin == -1)
    {
        cout << "Initial Node index or Final Node index doesn't belong to the given TPZGeoNode vector!\n";
        cout << "See NodesHunter method!\n"; exit(-1);
    }

    /// Computing BasisChange Matrix
    /// Where NewBase X_axis is defined by IniNode->FinNode orientation
    /// and   NewBase Y_axis is perpendicular to X_axis in XY plane counter-clockwise
    TPZFMatrix IfromCntoBase(dim,dim,0.);
    IfromCntoBase(dim-1,dim-1) = 1.;
    double norm = 0.;
    //computing X_axis
    for(int coord = 0; coord < dim; coord++)
    {
        double val = VectorialNotation.GetVal(posFin,coord) - VectorialNotation.GetVal(posIni,coord);
        IfromCntoBase.PutVal(0,coord,val);
        norm += val*val;
    }
    //normalizing X_axis
    for(int coord = 0; coord < dim; coord++)
    {
        IfromCntoBase.PutVal(0,coord,IfromCntoBase.GetVal(0,coord)/sqrt(norm));
    }
    //computing Y_axis, i.e., if X_axis=(x,y) then Y_axis=(-y,x)
    IfromCntoBase.PutVal(1,0,-IfromCntoBase.GetVal(0,1));
    IfromCntoBase.PutVal(1,1, IfromCntoBase.GetVal(0,0));

    /** Changing VectorialNotation from Canonic Base to NewBase */
    for(int i = 0; i < VecSize; i++)
    {
            vector <double> temp; temp.resize(3);
            for(int j = 0; j < dim; j++)
            {
                double val = 0.;
                for(int k = 0; k < dim; k++)
                {
                    val += IfromCntoBase.GetVal(j,k)*VectorialNotation.GetVal(i,k);
                }
                temp[j] = val;
            }
            for(int p = 0; p < dim; p++) VectorialNotation.PutVal(i,p,temp[p]);
    }

    /** Hunting Nodes */
    std::map <double,int> mymap;
    for(int h = 0; h < VecSize; h++)
    {
        if( VectorialNotation.GetVal(h,0) >= VectorialNotation.GetVal(posIni,0) &&
            VectorialNotation.GetVal(h,0) <= VectorialNotation.GetVal(posFin,0) &&
            (fabs(VectorialNotation.GetVal(h,1)) + fabs(VectorialNotation.GetVal(h,2))) < Tol )
        {
            std::pair< double , int> Item(VectorialNotation.GetVal(h,0), gMesh.NodeVec()[h].Id());
            mymap.insert(Item);
        }
    }
    NodesHunted.resize(mymap.size());
    std::map<double, int>::iterator it = mymap.begin();
    int i = 0;
    for(it = mymap.begin(); it!= mymap.end(); it++, i++)
    {
        NodesHunted[i] = it->second;
    }
}

void IntegralForces(TPZCompMesh * cmesh, int MatId, TPZVec<REAL> &forces)
{
  forces.Resize(1);
  forces[0] = 0.;
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

        TPZIntPoints * rule = cel->Reference()->CreateSideIntegrationRule(geoside.Side(),5);

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
        }

        delete rule;
    }
  }
  int bye=0;
}

TPZCompMesh * SquareMesh()
{
    TPZGeoMesh * gmesh = new TPZGeoMesh;

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

    vector<REAL> Orig(3);  Orig[0] = 2.5; Orig[1] = 0.; Orig[2] = 0.;
    vector<REAL> AxisZ(3); AxisZ[0] = 0.; AxisZ[1] = 1.; AxisZ[2] = 0.;
    vector<REAL> AxisR(3); AxisR[0] = 1.; AxisR[1] = 0.; AxisR[2] = 0.;

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

TPZCompMesh * SpiralMesh()
{
  ///Geometria Caixa Espiral 2D

  double Bb = 8.6; double Hr = 10.8; double Bt = 8.4; double Hl = 10.15; double Cx = 4.9; double Cy = 5.9; double R = 3.4;
  double h = 2.1; double b = 0.55; double Dx = 0.65; double Dy = 3.5; double e1 = 0.25; double e2 = 0.45; double e3 = 0.1;
  double X1 = 3.8; double X2 = 1.; double X3 = 1.9; double X4 = 0.6; double X5 = 0.3; double h1 = 0.45; double h2 = 0.45; double Py = 3.3;
  TPZGeoMesh * gmesh = CxEspiral2D(Bb, Hr, Bt, Hl, Cx, Cy, R, h, b, Dx, Dy, e1, e2, e3, X1, X2, X3, X4, X5, h1, h2, Py);

///Materials
  REAL fx = 0., fy = 0.;//-0.025;

  double E = 18.5*1000.; // Young modulus = 18.5 GPa = 18.5*10^3 MPa
  double nu = 0.2; // poisson coefficient

  TPZAutoPointer<TPZMaterial> mat = new TPZElasticityAxiMaterial(1, E, nu, fx, fy);
  vector<REAL> Orig(3);  Orig[0] = 15.;  Orig[1] = 0.;  Orig[2] = 0.;
  vector<REAL> AxisZ(3); AxisZ[0] = 0.; AxisZ[1] = 1.; AxisZ[2] = 0.;
  vector<REAL> AxisR(3); AxisR[0] = -1.; AxisR[1] = 0.; AxisR[2] = 0.;
  (dynamic_cast<TPZElasticityAxiMaterial*>(mat.operator->()))->SetOrigin(Orig, AxisZ, AxisR);

  //MohrCoulomb data
//   double phi = M_PI/6.;
//   double fc = 15./1.4; //fck = 15 MPa, fc = fck/1.4
//   double c = fc*(1. - sin(phi))/(2.*cos(phi));
//   (dynamic_cast<TPZElasticityAxiMaterial*>(mat.operator->()))->SetMohrCoulomb(c,phi);

  TPZFMatrix Base1(2,2,0.), Base2(2,1,0.);
  Base1(1,1) = 1.E20;
  TPZAutoPointer<TPZMaterial> bcmatBase = mat->CreateBC(mat, -8, 2, Base1, Base2);

  TPZFMatrix Force1(2,2,0.), Force2(2,1,0.);
  Force2(0,0) = -700.;
  TPZAutoPointer<TPZMaterial> bcmatForce = mat->CreateBC(mat, -7, 3, Force1, Force2);

//   TPZFMatrix Dot1(2,2,0.), Dot2(2,1,0.);
//   Dot1(1,1) = 0.01;
//   TPZAutoPointer<TPZMaterial> bcmatDot = mat->CreateBC(mat, -3, 2, Dot1, Dot2);
//**************************

///Computational Mesh
  const int dim = 2;
  TPZCompEl::SetgOrder(3);
  TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
  cmesh->InsertMaterialObject(mat);
  cmesh->InsertMaterialObject(bcmatBase);
  cmesh->InsertMaterialObject(bcmatForce);
//   cmesh->InsertMaterialObject(bcmatDot);
//**************************

  cmesh->AutoBuild();

  return cmesh; 
}
