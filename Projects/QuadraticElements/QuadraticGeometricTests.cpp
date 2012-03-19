#include <iostream>
#include "pzfmatrix.h"
#include "pzelasmat.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzbndcond.h"
#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "pzlog.h"
#include "tpzgeoblend.h"
#include "tpzellipse3d.h"
#include "tpzarc3d.h"
#include "TPZGeoLinear.h"
#include "pzgeotriangle.h"
#include "pzgeopoint.h"
#include "pzgeopyramid.h"
#include "pzgeoquad.h"
#include "TPZRefPatternDataBase.h"
#include "TPZRefPatternTools.h"
#include "TPZVTKGeoMesh.h"
#include "pzgmesh.h"
#include "tpzmathtools.h"

#include "tpzchangeel.h"
#include "tpzquadraticline.h"
#include "tpzquadratictrig.h"
#include "tpzquadraticquad.h"
#include "tpzquadratictetra.h"
#include "tpzquadraticpyramid.h"
#include "tpzquadraticprism.h"
#include "tpzquadraticcube.h"

#include "pzgengrid.h"

using namespace std;

//----------------------------------------------------------------------------------------------------------------------------------

const int matElId = 1;

//Line
/*
 int main(int argc, char * const argv[])
 {
 gRefDBase.InitializeUniformRefPattern(EOned);
 
 TPZGeoMesh * gmesh = new TPZGeoMesh;
 
 const int Qnodes = 2;
 TPZVec < TPZVec <REAL> > NodeCoord(Qnodes);
 for(int i = 0; i < Qnodes; i++) NodeCoord[i].Resize(3,0.);
 
 //setting nodes coords
 int nodeId = 0;
 double d = 0.2;
 
 NodeCoord[nodeId][0] = -1.+d;
 NodeCoord[nodeId][1] =  0.;
 NodeCoord[nodeId][2] =  0.-d;
 nodeId++;
 
 NodeCoord[nodeId][0] = +1.;
 NodeCoord[nodeId][1] =  0.+d;
 NodeCoord[nodeId][2] =  0.+d;
 
 //initializing gmesh->NodeVec()
 gmesh->NodeVec().Resize(Qnodes);
 TPZVec <TPZGeoNode> Node(Qnodes);
 for(int n = 0; n < Qnodes; n++)
 {
 Node[n].SetNodeId(n);
 Node[n].SetCoord(NodeCoord[n]);
 gmesh->NodeVec()[n] = Node[n]; 
 }
 
 int elId = 0;
 TPZVec <int> Topol(Qnodes);
 for(int n = 0; n < Qnodes; n++) Topol[n] = n;
 new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elId,Topol,matElId,*gmesh);
 elId++;
 
 gmesh->BuildConnectivity(); 
 
 std::ofstream outAntes("meshAntes.txt");
 gmesh->Print(outAntes);
 
 TPZChangeEl::ChangeToQuadratic(gmesh, 0);
 
 std::ofstream outDepois("meshDepois.txt");
 gmesh->Print(outDepois);
 
 TPZVec<REAL> nodeCH(3);
 for(int c = 0; c < 3; c++) nodeCH[c] = gmesh->NodeVec()[2].Coord(c) + 0.1;
 gmesh->NodeVec()[2].SetCoord(nodeCH);
 
 //Teste de convergência (tem que dar muito proximo de 2)
 std::cout << "Convergence Order:\n";
 TPZGeoEl * myGel = gmesh->ElementVec()[0];
 TPZVec<REAL> qsi(myGel->Dimension(),0.);
 TPZMathTools conv; 
 conv.JacobianConv(*myGel, qsi);
 
 //Teste somatoria das phi_s (tem que dar 1 em qq pto)
 std::cout << "Summ of phi's at many points:\n";
 int npts = 4;
 
 pzgeom::TPZQuadraticLine * myQuadraticGel = dynamic_cast<pzgeom::TPZQuadraticLine*>(myGel);
 TPZFMatrix<REAL> phi(myGel->NNodes(),1), dphi(myGel->Dimension(),myGel->NNodes());
 for(int xi = 0; xi <= npts; xi++)
 {
 qsi[0] = -1. + xi*2./npts;
 
 myQuadraticGel->Shape(qsi, phi, dphi);
 double sum = 0.;
 for(int s = 0; s < myGel->NNodes(); s++) sum += phi[s];
 std::cout << sum << std::endl;
 }    
 
 //Teste do refinamento
 std::cout << "Uniform refinement:\n";
 int nDiv = 5;
 for(int D = 0; D < nDiv; D++)
 {
 int nels = gmesh->NElements();
 for(int elem = 0; elem < nels; elem++)
 {    
 TPZVec< TPZGeoEl * > filhos;
 TPZGeoEl * gel = gmesh->ElementVec()[elem];
 if(!gel) continue;
 gel->Divide(filhos);
 }
 }
 std::ofstream out("QuadraticLine.vtk");
 TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, false);
 
 std::cout << "FINISHED!!!" << std::endl;
 
 return 0;
 }
 */

//Triangle
/*
 int main(int argc, char * const argv[])
 {
 gRefDBase.InitializeUniformRefPattern(ETriangle);
 
 TPZGeoMesh * gmesh = new TPZGeoMesh;
 
 const int Qnodes = 3;
 TPZVec < TPZVec <REAL> > NodeCoord(Qnodes);
 for(int i = 0; i < Qnodes; i++) NodeCoord[i].Resize(3,0.);
 
 //setting nodes coords
 int nodeId = 0;
 double d = 0.1;
 
 NodeCoord[nodeId][0] =  0.+d;
 NodeCoord[nodeId][1] =  0.;
 NodeCoord[nodeId][2] =  0.-d;
 nodeId++;
 
 NodeCoord[nodeId][0] = +1.;
 NodeCoord[nodeId][1] =  0.+d;
 NodeCoord[nodeId][2] =  0.+d;
 nodeId++;
 
 NodeCoord[nodeId][0] =  0.-d;
 NodeCoord[nodeId][1] = +1.+d;
 NodeCoord[nodeId][2] =  0.+d;
 nodeId++;
 
 //initializing gmesh->NodeVec()
 gmesh->NodeVec().Resize(Qnodes);
 TPZVec <TPZGeoNode> Node(Qnodes);
 for(int n = 0; n < Qnodes; n++)
 {
 Node[n].SetNodeId(n);
 Node[n].SetCoord(NodeCoord[n]);
 gmesh->NodeVec()[n] = Node[n]; 
 }
 
 //inserting quadrilaterals
 int elId = 0;
 TPZVec <int> Topol(Qnodes);
 for(int n = 0; n < Qnodes; n++) Topol[n] = n;
 new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (elId,Topol,matElId,*gmesh);
 elId++;
 
 gmesh->BuildConnectivity();
 
 std::ofstream outAntes("meshAntes.txt");
 gmesh->Print(outAntes);
 
 TPZChangeEl::ChangeToQuadratic(gmesh, 0);
 
 std::ofstream outDepois("meshDepois.txt");
 gmesh->Print(outDepois);
 
 TPZVec<REAL> nodeCH(3);
 for(int n = 0; n < 6; n++)
 {  
 for(int c = 0; c < 3; c++)
 {
 double num = (2.*double(rand()%11) - 10.)/100.;
 nodeCH[c] = gmesh->NodeVec()[n].Coord(c) + num;   
 }
 gmesh->NodeVec()[n].SetCoord(nodeCH);
 }
 
 
 //Teste de convergência (tem que dar muito proximo de 2)
 std::cout << "Convergence Order:\n";
 TPZGeoEl * myGel = gmesh->ElementVec()[0];
 TPZVec<REAL> qsi(myGel->Dimension(),0.);
 TPZMathTools conv; 
 conv.JacobianConv(*myGel, qsi);
 
 //Teste somatoria das phi_s (tem que dar 1 em qq pto)
 std::cout << "Summ of phi's at many points:\n";
 int npts = 4;
 
 pzgeom::TPZQuadraticTrig * myQuadraticGel = dynamic_cast<pzgeom::TPZQuadraticTrig*>(myGel);
 TPZFMatrix<REAL> phi(myGel->NNodes(),1), dphi(myGel->Dimension(),myGel->NNodes());
 for(int xi = 0; xi <= npts; xi++)
 {
 for(int et = 0; et <= npts; et++)
 {
 qsi[0] = 0. + xi*2./npts;
 qsi[1] = 0. + et*(1.- qsi[0])/npts;
 
 myQuadraticGel->Shape(qsi, phi, dphi);
 double sum = 0.;
 for(int s = 0; s < myGel->NNodes(); s++) sum += phi[s];
 std::cout << sum << std::endl;
 }
 }    
 
 //Teste do refinamento
 std::cout << "Uniform refinement:\n";
 int nDiv = 5;
 for(int D = 0; D < nDiv; D++)
 {
 int nels = gmesh->NElements();
 for(int elem = 0; elem < nels; elem++)
 {    
 TPZVec< TPZGeoEl * > filhos;
 TPZGeoEl * gel = gmesh->ElementVec()[elem];
 gel->Divide(filhos);
 }
 }
 std::ofstream out("QuadraticTriang.vtk");
 TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, false);
 
 std::cout << "FINISHED!!!" << std::endl;
 
 return 0;
 }
 */

//Quadrilateral

int main(int argc, char * const argv[])
{
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    
    const int Qnodes = 6;
    TPZVec < TPZVec <REAL> > NodeCoord(Qnodes);
    for(int i = 0; i < Qnodes; i++) NodeCoord[i].Resize(3,0.);
    
    //setting nodes coords
    int nodeId = 0;
    
    NodeCoord[nodeId][0] = 0.;
    NodeCoord[nodeId][1] = 0.;
    NodeCoord[nodeId][2] = 0.;
    nodeId++;
    
    NodeCoord[nodeId][0] = 1.;
    NodeCoord[nodeId][1] = 0.;
    NodeCoord[nodeId][2] = 0.;
    nodeId++;
    
    NodeCoord[nodeId][0] = 1.;
    NodeCoord[nodeId][1] = 1.;
    NodeCoord[nodeId][2] = 0.;
    nodeId++;
    
    NodeCoord[nodeId][0] = 0.;
    NodeCoord[nodeId][1] = 1.;
    NodeCoord[nodeId][2] = 0.;
    nodeId++;

    NodeCoord[nodeId][0] = 2.;
    NodeCoord[nodeId][1] = 0.;
    NodeCoord[nodeId][2] = 0.;
    nodeId++;
    
    NodeCoord[nodeId][0] = 2.;
    NodeCoord[nodeId][1] = 1.;
    NodeCoord[nodeId][2] = 0.;
    nodeId++;
    
    //initializing gmesh->NodeVec()
    gmesh->NodeVec().Resize(Qnodes);
    TPZVec <TPZGeoNode> Node(Qnodes);
    for(int n = 0; n < Qnodes; n++)
    {
        Node[n].SetNodeId(n);
        Node[n].SetCoord(NodeCoord[n]);
        gmesh->NodeVec()[n] = Node[n]; 
    }
    
    //inserting quadrilaterals
    int elId = 0;
    TPZVec <int> Topol(4);
    
    for(int n = 0; n < 4; n++) Topol[n] = n;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elId,Topol,matElId,*gmesh);
    elId++;
    
    Topol[0] = 1; Topol[1] = 4; Topol[2] = 5; Topol[3] = 2;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elId,Topol,matElId,*gmesh);
    elId++;    
    
    gmesh->BuildConnectivity();
    
    std::ofstream outAntes("meshAntes.txt");
    gmesh->Print(outAntes);
    
//    TPZChangeEl::ChangeToQuadratic(gmesh, 0);
//    TPZChangeEl::ChangeToQuadratic(gmesh, 1);
    
    TPZChangeEl::ChangeToQuarterPoint(gmesh, 0, 0);
    TPZChangeEl::ChangeToQuarterPoint(gmesh, 1, 1);
    
    std::ofstream outDepois("meshDepois.txt");
    gmesh->Print(outDepois);
    
//    TPZVec<REAL> nodeCH(3);
//    for(int n = 0; n < gmesh->NNodes(); n++)
//    {  
//        for(int c = 0; c < 3; c++)
//        {
//            double num = (2.*double(rand()%11) - 10.)/100.;
//            nodeCH[c] = gmesh->NodeVec()[n].Coord(c) + num;   
//        }
//        gmesh->NodeVec()[n].SetCoord(nodeCH);
//    }
    
    //Teste de convergência (tem que dar muito proximo de 2)
    std::cout << "Convergence Order:\n";
    TPZGeoEl * myGel = gmesh->ElementVec()[0];
    TPZVec<REAL> qsi(myGel->Dimension(),0.);
    TPZMathTools conv; 
    conv.JacobianConv(*myGel, qsi);
    
    //Teste somatoria das phi_s (tem que dar 1 em qq pto)
    std::cout << "Summ of phi's at many points:\n";
    int npts = 4;
    
    pzgeom::TPZQuadraticQuad * myQuadraticGel = dynamic_cast<pzgeom::TPZQuadraticQuad*>(myGel);
    TPZFMatrix<REAL> phi(myGel->NNodes(),1), dphi(myGel->Dimension(),myGel->NNodes());
    for(int xi = 0; xi <= npts; xi++)
    {
        for(int et = 0; et <= npts; et++)
        {
            qsi[0] = -1. + xi*2./npts;
            qsi[1] = -1. + et*2./npts;
            
            myQuadraticGel->Shape(qsi, phi, dphi);
            double sum = 0.;
            for(int s = 0; s < myGel->NNodes(); s++) sum += phi[s];
            std::cout << sum << std::endl;
        }
    }    
    
    //Teste do refinamento
    std::cout << "Uniform refinement:\n";
    int nDiv = 1;
    for(int D = 0; D < nDiv; D++)
    {
        int nels = gmesh->NElements();
        for(int elem = 0; elem < nels; elem++)
        {    
            TPZVec< TPZGeoEl * > filhos;
            TPZGeoEl * gel = gmesh->ElementVec()[elem];
            gel->Divide(filhos);
        }
    }
    std::ofstream out("QuadraticQuad.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, false);
    
    std::cout << "FINISHED!!!" << std::endl;
    
    return 0;
}


//Tetraedro
/*
 int main(int argc, char * const argv[])
 {
 //    gRefDBase.InitializeUniformRefPattern(ECube);
 //    gRefDBase.InitializeUniformRefPattern(EPiramide);
 
 TPZGeoMesh * gmesh = new TPZGeoMesh;
 
 const int Qnodes = 4;
 TPZVec < TPZVec <REAL> > NodeCoord(Qnodes);
 for(int i = 0; i < Qnodes; i++) NodeCoord[i].Resize(3,0.);
 
 //setting nodes coords
 int nodeId = 0;
 double d = 0.0;
 
 NodeCoord[nodeId][0] =  0.+d;
 NodeCoord[nodeId][1] =  0.;
 NodeCoord[nodeId][2] =  0.-d;
 nodeId++;
 
 NodeCoord[nodeId][0] = +1.;
 NodeCoord[nodeId][1] =  0.+d;
 NodeCoord[nodeId][2] =  0.+d;
 nodeId++;
 
 NodeCoord[nodeId][0] =  0.-d;
 NodeCoord[nodeId][1] = +1.+d;
 NodeCoord[nodeId][2] =  0.+d;
 nodeId++;
 
 NodeCoord[nodeId][0] =  0.-d;
 NodeCoord[nodeId][1] =  0.+d;
 NodeCoord[nodeId][2] = +1.+d;
 nodeId++;
 
 //initializing gmesh->NodeVec()
 gmesh->NodeVec().Resize(Qnodes);
 TPZVec <TPZGeoNode> Node(Qnodes);
 for(int n = 0; n < Qnodes; n++)
 {
 Node[n].SetNodeId(n);
 Node[n].SetCoord(NodeCoord[n]);
 gmesh->NodeVec()[n] = Node[n]; 
 }
 
 int elId = 0;
 TPZVec <int> Topol(Qnodes);
 for(int n = 0; n < Qnodes; n++) Topol[n] = n;
 new TPZGeoElRefPattern< pzgeom::TPZGeoTetrahedra > (elId,Topol,matElId,*gmesh);
 elId++;
 
 gmesh->BuildConnectivity(); 
 
 std::ofstream outAntes("meshAntes.txt");
 gmesh->Print(outAntes);
 
 TPZChangeEl::ChangeToQuadratic(gmesh, 0);
 
 std::ofstream outDepois("meshDepois.txt");
 gmesh->Print(outDepois);
 
 //    TPZVec<REAL> nodeCH(3);
 //    for(int n = 0; n < 10; n++)
 //    {  
 //        for(int c = 0; c < 3; c++)
 //        {
 //            double num = (2.*double(rand()%11) - 10.)/100.;
 //            nodeCH[c] = gmesh->NodeVec()[n].Coord(c) + num;   
 //        }
 //        gmesh->NodeVec()[n].SetCoord(nodeCH);
 //    }
 
 {
 int targetSide = 4;
 TPZGeoEl * quarterGel = TPZChangeEl::ChangeToQuarterPoint(gmesh, 0, targetSide);
 
 std::cout << "QuarterPoint:\n";
 int nDiv = 2;
 for(int D = 0; D < nDiv; D++)
 {
 int nels = gmesh->NElements();
 for(int elem = 0; elem < nels; elem++)
 {    
 TPZVec< TPZGeoEl * > filhos;
 quarterGel->Divide(filhos);
 }
 }
 std::ofstream outQPTetra("QuarterPointTetraedron.vtk");
 TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outQPTetra, false);
 }
 
 //Teste de convergência (tem que dar muito proximo de 2)
 std::cout << "Convergence Order:\n";
 TPZGeoEl * myGel = gmesh->ElementVec()[0];
 TPZVec<REAL> qsi(myGel->Dimension(),0.);
 TPZMathTools conv; 
 conv.JacobianConv(*myGel, qsi);
 
 //Teste somatoria das phi_s (tem que dar 1 em qq pto)
 std::cout << "Summ of phi's at many points:\n";
 int npts = 4;
 pzgeom::TPZQuadraticTetra * myQuadraticGel = dynamic_cast<pzgeom::TPZQuadraticTetra*>(myGel);
 TPZFMatrix<REAL> phi(myGel->NNodes(),1), dphi(myGel->Dimension(),myGel->NNodes());
 for(int xi = 0; xi <= npts; xi++)
 {
 for(int et = 0; et <= npts; et++)
 {
 for(int zet = 0; zet <= npts; zet++)
 {
 qsi[0] = -(1.-zet) + xi*2./npts;
 qsi[1] = -(1.-zet) + et*2./npts;
 qsi[2] = -(1.-zet) + zet*2./npts;
 myQuadraticGel->Shape(qsi, phi, dphi);
 double sum = 0.;
 for(int s = 0; s < myGel->NNodes(); s++) sum += phi[s];
 std::cout << sum << std::endl;
 }
 }
 }    
 
 //Teste do refinamento
 std::cout << "Uniform refinement:\n";
 int nDiv = 3;
 for(int D = 0; D < nDiv; D++)
 {
 int nels = gmesh->NElements();
 for(int elem = 0; elem < nels; elem++)
 {    
 TPZVec< TPZGeoEl * > filhos;
 TPZGeoEl * gel = gmesh->ElementVec()[elem];
 gel->Divide(filhos);
 }
 }
 std::ofstream outPyram("QuadraticTetra.vtk");
 TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outPyram, false);
 
 std::cout << "FINISHED!!!" << std::endl;
 
 return 0;
 }
 */

//Pyramid
/*
 int main(int argc, char * const argv[])
 {
 //    gRefDBase.InitializeUniformRefPattern(ECube);
 //    gRefDBase.InitializeUniformRefPattern(EPiramide);
 
 TPZGeoMesh * gmesh = new TPZGeoMesh;
 
 const int Qnodes = 5;
 TPZVec < TPZVec <REAL> > NodeCoord(Qnodes);
 for(int i = 0; i < Qnodes; i++) NodeCoord[i].Resize(3,0.);
 
 //setting nodes coords
 int nodeId = 0;
 double d = 0.0;
 
 NodeCoord[nodeId][0] = -1.;
 NodeCoord[nodeId][1] = -1.;
 NodeCoord[nodeId][2] =  0.;
 nodeId++;
 
 NodeCoord[nodeId][0] = +1.;
 NodeCoord[nodeId][1] = -1.;
 NodeCoord[nodeId][2] =  0.;
 nodeId++;
 
 NodeCoord[nodeId][0] = +1.;
 NodeCoord[nodeId][1] = +1.;
 NodeCoord[nodeId][2] =  0.;
 nodeId++;
 
 NodeCoord[nodeId][0] = -1.;
 NodeCoord[nodeId][1] = +1.;
 NodeCoord[nodeId][2] =  0.;
 nodeId++;
 
 NodeCoord[nodeId][0] =  0.;
 NodeCoord[nodeId][1] =  0.;
 NodeCoord[nodeId][2] =  1.;
 nodeId++;
 
 //initializing gmesh->NodeVec()
 gmesh->NodeVec().Resize(Qnodes);
 TPZVec <TPZGeoNode> Node(Qnodes);
 for(int n = 0; n < Qnodes; n++)
 {
 Node[n].SetNodeId(n);
 Node[n].SetCoord(NodeCoord[n]);
 gmesh->NodeVec()[n] = Node[n]; 
 }
 
 int elId = 0;
 TPZVec <int> Topol(Qnodes);
 for(int n = 0; n < Qnodes; n++) Topol[n] = n;
 new TPZGeoElRefPattern< pzgeom::TPZGeoPyramid > (elId,Topol,matElId,*gmesh);
 elId++;
 
 gmesh->BuildConnectivity(); 
 
 std::ofstream outAntes("meshAntes.txt");
 gmesh->Print(outAntes);
 
 TPZChangeEl::ChangeToQuadratic(gmesh, 0);
 
 std::ofstream outDepois("meshDepois.txt");
 gmesh->Print(outDepois);
 
 //    TPZVec<REAL> nodeCH(3);
 //    for(int n = 0; n < 13; n++)
 //    {  
 //        for(int c = 0; c < 3; c++)
 //        {
 //            double num = (2.*double(rand()%11) - 10.)/100.;
 //            nodeCH[c] = gmesh->NodeVec()[n].Coord(c) + num;   
 //        }
 //        gmesh->NodeVec()[n].SetCoord(nodeCH);
 //    }
 
 {
 int targetSide = 17;
 TPZGeoEl * quarterGel = TPZChangeEl::ChangeToQuarterPoint(gmesh, 0, targetSide);
 
 std::cout << "QuarterPoint:\n";
 int nDiv = 2;
 for(int D = 0; D < nDiv; D++)
 {
 int nels = gmesh->NElements();
 for(int elem = 0; elem < nels; elem++)
 {    
 TPZVec< TPZGeoEl * > filhos;
 quarterGel->Divide(filhos);
 }
 }
 std::ofstream outQPPir("QuarterPointPiramid.vtk");
 TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outQPPir, false);
 }
 
 //Teste de convergência (tem que dar muito proximo de 2)
 std::cout << "Convergence Order:\n";
 TPZGeoEl * myGel = gmesh->ElementVec()[0];
 TPZVec<REAL> qsi(myGel->Dimension(),0.);
 TPZMathTools conv; 
 conv.JacobianConv(*myGel, qsi);
 
 //Teste somatoria das phi_s (tem que dar 1 em qq pto)
 std::cout << "Summ of phi's at many points:\n";
 int npts = 4;
 pzgeom::TPZQuadraticPyramid * myPyramid = dynamic_cast<pzgeom::TPZQuadraticPyramid*>(myGel);
 TPZFMatrix<REAL> phi(myGel->NNodes(),1), dphi(myGel->Dimension(),myGel->NNodes());
 for(int zet = 0; zet <= npts; zet++)
 {
 for(int et = 0; et <= npts; et++)
 {
 for(int xi = 0; xi <= npts; xi++)
 {
 qsi[2] = zet*1./npts;
 qsi[0] = -(1.-qsi[2]) + xi*2.*(1.-qsi[2])/npts;
 qsi[1] = -(1.-qsi[2]) + et*2.*(1.-qsi[2])/npts;
 myPyramid->Shape(qsi, phi, dphi);
 double sum = 0.;
 for(int s = 0; s < myGel->NNodes(); s++)
 {
 sum += phi[s];       
 }
 std::cout << "sum = " << sum << std::endl;
 if(fabs(sum-1.) > 1.E-8)
 {
 std::cout << qsi[0] << " , " << qsi[1] << " , " << qsi[2] << std::endl;
 for(int s = 0; s < myGel->NNodes(); s++)
 {
 std::cout << phi[s] << std::endl;
 }
 std::cout << "===========================\n";
 }
 }
 }
 }    
 
 //Teste do refinamento
 std::cout << "Uniform refinement:\n";
 int nDiv = 3;
 for(int D = 0; D < nDiv; D++)
 {
 int nels = gmesh->NElements();
 for(int elem = 0; elem < nels; elem++)
 {    
 TPZVec< TPZGeoEl * > filhos;
 TPZGeoEl * gel = gmesh->ElementVec()[elem];
 gel->Divide(filhos);
 }
 }
 std::ofstream outPyram("QuadraticPyramid.vtk");
 TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outPyram, false);
 
 std::cout << "FINISHED!!!" << std::endl;
 
 
 
 return 0;
 }
 */

//Prism
/*
 int main(int argc, char * const argv[])
 {
 //    gRefDBase.InitializeUniformRefPattern(ECube);
 //    gRefDBase.InitializeUniformRefPattern(EPiramide);
 
 TPZGeoMesh * gmesh = new TPZGeoMesh;
 
 const int Qnodes = 6;
 TPZVec < TPZVec <REAL> > NodeCoord(Qnodes);
 for(int i = 0; i < Qnodes; i++) NodeCoord[i].Resize(3,0.);
 
 //setting nodes coords
 int nodeId = 0;
 double d = 0.0;
 
 NodeCoord[nodeId][0] =  0.+d;
 NodeCoord[nodeId][1] =  0.;
 NodeCoord[nodeId][2] = -1.-d;
 nodeId++;
 
 NodeCoord[nodeId][0] = +1.;
 NodeCoord[nodeId][1] =  0.+d;
 NodeCoord[nodeId][2] = -1.+d;
 nodeId++;
 
 NodeCoord[nodeId][0] =  0.-d;
 NodeCoord[nodeId][1] = +1.+d;
 NodeCoord[nodeId][2] = -1.+d;
 nodeId++;
 
 NodeCoord[nodeId][0] =  0.-d;
 NodeCoord[nodeId][1] =  0.+d;
 NodeCoord[nodeId][2] = +1.;
 nodeId++;
 
 NodeCoord[nodeId][0] = +1.-d;
 NodeCoord[nodeId][1] =  0.-d;
 NodeCoord[nodeId][2] = +1.+d;
 nodeId++;
 
 NodeCoord[nodeId][0] =  0.;
 NodeCoord[nodeId][1] = +1.-d;
 NodeCoord[nodeId][2] = +1.-d;
 nodeId++;
 
 //initializing gmesh->NodeVec()
 gmesh->NodeVec().Resize(Qnodes);
 TPZVec <TPZGeoNode> Node(Qnodes);
 for(int n = 0; n < Qnodes; n++)
 {
 Node[n].SetNodeId(n);
 Node[n].SetCoord(NodeCoord[n]);
 gmesh->NodeVec()[n] = Node[n]; 
 }
 
 //inserting quadrilaterals
 int elId = 0;
 TPZVec <int> Topol(Qnodes);
 for(int n = 0; n < Qnodes; n++) Topol[n] = n;
 new TPZGeoElRefPattern< pzgeom::TPZGeoPrism > (elId,Topol,matElId,*gmesh);
 elId++;
 
 gmesh->BuildConnectivity(); 
 
 std::ofstream outAntes("meshAntes.txt");
 gmesh->Print(outAntes);
 
 TPZChangeEl::ChangeToQuadratic(gmesh, 0);
 
 std::ofstream outDepois("meshDepois.txt");
 gmesh->Print(outDepois);
 
 TPZVec<REAL> nodeCH(3);
 for(int n = 0; n < 15; n++)
 {  
 for(int c = 0; c < 3; c++)
 {
 double num = (2.*double(rand()%11) - 10.)/100.;
 nodeCH[c] = gmesh->NodeVec()[n].Coord(c) + num;   
 }
 gmesh->NodeVec()[n].SetCoord(nodeCH);
 }
 
 //Teste de convergência (tem que dar muito proximo de 2)
 std::cout << "Convergence Order:\n";
 TPZGeoEl * myGel = gmesh->ElementVec()[0];
 TPZVec<REAL> qsi(myGel->Dimension(),0.);
 TPZMathTools conv; 
 conv.JacobianConv(*myGel, qsi);
 
 //Teste somatoria das phi_s (tem que dar 1 em qq pto)
 std::cout << "Summ of phi's at many points:\n";
 int npts = 4;
 pzgeom::TPZQuadraticPrism * myPrism = dynamic_cast<pzgeom::TPZQuadraticPrism*>(myGel);
 TPZFMatrix<REAL> phi(myGel->NNodes(),1), dphi(myGel->Dimension(),myGel->NNodes());
 for(int xi = 0; xi <= npts; xi++)
 {
 for(int et = 0; et <= npts; et++)
 {
 for(int zet = 0; zet <= npts; zet++)
 {
 qsi[0] = -(1.-zet) + xi*2./npts;
 qsi[1] = -(1.-zet) + et*2./npts;
 qsi[2] = -(1.-zet) + zet*2./npts;
 myPrism->Shape(qsi, phi, dphi);
 double sum = 0.;
 for(int s = 0; s < myGel->NNodes(); s++) sum += phi[s];
 std::cout << sum << std::endl;
 }
 }
 }    
 
 //Teste do refinamento
 std::cout << "Uniform refinement:\n";
 int nDiv = 3;
 for(int D = 0; D < nDiv; D++)
 {
 int nels = gmesh->NElements();
 for(int elem = 0; elem < nels; elem++)
 {    
 TPZVec< TPZGeoEl * > filhos;
 TPZGeoEl * gel = gmesh->ElementVec()[elem];
 gel->Divide(filhos);
 }
 }
 std::ofstream outPrism("QuadraticPrism.vtk");
 TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outPrism, false);
 
 std::cout << "FINISHED!!!" << std::endl;
 
 return 0;
 }
 */

//Cube
/*
 int main(int argc, char * const argv[])
 {
 gRefDBase.InitializeUniformRefPattern(ECube);
 
 TPZGeoMesh * gmesh = new TPZGeoMesh;
 
 const int Qnodes = 8;
 TPZVec < TPZVec <REAL> > NodeCoord(Qnodes);
 for(int i = 0; i < Qnodes; i++) NodeCoord[i].Resize(3,0.);
 
 //setting nodes coords
 int nodeId = 0;
 double d = 0.0;
 
 NodeCoord[nodeId][0] = -1.+d;
 NodeCoord[nodeId][1] = -1.;
 NodeCoord[nodeId][2] = -1.-d;
 nodeId++;
 
 NodeCoord[nodeId][0] = +1.;
 NodeCoord[nodeId][1] = -1.+d;
 NodeCoord[nodeId][2] = -1.+d;
 nodeId++;
 
 NodeCoord[nodeId][0] = +1.-d;
 NodeCoord[nodeId][1] = +1.+d;
 NodeCoord[nodeId][2] = -1.+d;
 nodeId++;
 
 NodeCoord[nodeId][0] = -1.;
 NodeCoord[nodeId][1] = +1.+d;
 NodeCoord[nodeId][2] = -1.+d;
 nodeId++;
 //
 //
 NodeCoord[nodeId][0] = -1.-d;
 NodeCoord[nodeId][1] = -1.;
 NodeCoord[nodeId][2] = +1.;
 nodeId++;
 
 NodeCoord[nodeId][0] = +1.;
 NodeCoord[nodeId][1] = -1.+d;
 NodeCoord[nodeId][2] = +1.;
 nodeId++;
 
 NodeCoord[nodeId][0] = +1.;
 NodeCoord[nodeId][1] = +1.;
 NodeCoord[nodeId][2] = +1.-d;
 nodeId++;
 
 NodeCoord[nodeId][0] = -1.;
 NodeCoord[nodeId][1] = +1.+d;
 NodeCoord[nodeId][2] = +1.;
 nodeId++;
 
 //initializing gmesh->NodeVec()
 gmesh->NodeVec().Resize(Qnodes);
 TPZVec <TPZGeoNode> Node(Qnodes);
 for(int n = 0; n < Qnodes; n++)
 {
 Node[n].SetNodeId(n);
 Node[n].SetCoord(NodeCoord[n]);
 gmesh->NodeVec()[n] = Node[n]; 
 }
 
 //inserting quadrilaterals
 int elId = 0;
 TPZVec <int> Topol(Qnodes);
 for(int n = 0; n < Qnodes; n++) Topol[n] = n;
 new TPZGeoElRefPattern< pzgeom::TPZGeoCube > (elId,Topol,matElId,*gmesh);
 elId++;
 
 gmesh->BuildConnectivity(); 
 
 std::ofstream outAntes("meshAntes.txt");
 gmesh->Print(outAntes);
 
 TPZChangeEl::ChangeToQuadratic(gmesh, 0);
 
 std::ofstream outDepois("meshDepois.txt");
 gmesh->Print(outDepois);
 
 //Dando um peteleco nos nohs
 //    TPZVec<REAL> nodeCH(3);
 //    for(int n = 0; n < 20; n++)
 //    {  
 //        for(int c = 0; c < 3; c++)
 //        {
 //            double num = (2.*double(rand()%11) - 10.)/100.;
 //            nodeCH[c] = gmesh->NodeVec()[n].Coord(c) + num;   
 //        }
 //        gmesh->NodeVec()[n].SetCoord(nodeCH);
 //    }
 
 {
 int targetSide = 20;
 TPZGeoEl * quarterGel = TPZChangeEl::ChangeToQuarterPoint(gmesh, 0, targetSide);
 
 std::cout << "QuarterPoint:\n";
 int nDiv = 2;
 for(int D = 0; D < nDiv; D++)
 {
 int nels = gmesh->NElements();
 for(int elem = 0; elem < nels; elem++)
 {    
 TPZVec< TPZGeoEl * > filhos;
 quarterGel->Divide(filhos);
 }
 }
 std::ofstream outQPCube("QuarterPointCube.vtk");
 TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outQPCube, false);
 }
 
 //Teste de convergência (tem que dar muito proximo de 2)
 std::cout << "Convergence Order:\n";
 TPZGeoEl * myGel = gmesh->ElementVec()[0];
 TPZVec<REAL> qsi(myGel->Dimension(),0.);
 TPZMathTools conv; 
 conv.JacobianConv(*myGel, qsi);
 
 //Teste somatoria das phi_s (tem que dar 1 em qq pto)
 std::cout << "Summ of phi's at many points:\n";
 int npts = 4;
 TPZQuadraticCube * myCube = dynamic_cast<TPZQuadraticCube*>(myGel);
 TPZFMatrix<REAL> phi(myGel->NNodes(),1), dphi(myGel->Dimension(),myGel->NNodes());
 for(int xi = 0; xi <= npts; xi++)
 {
 for(int et = 0; et <= npts; et++)
 {
 for(int zet = 0; zet <= npts; zet++)
 {
 qsi[0] = -1. + xi*2./npts;
 qsi[1] = -1. + et*2./npts;
 qsi[2] = -1. + zet*2./npts;
 myCube->Shape(qsi, phi, dphi);
 double sum = 0.;
 for(int s = 0; s < myGel->NNodes(); s++) sum += phi[s];
 std::cout << sum << std::endl;
 }
 }
 }    
 
 //Teste do refinamento
 std::cout << "Uniform refinement:\n";
 int nDiv = 2;
 for(int D = 0; D < nDiv; D++)
 {
 int nels = gmesh->NElements();
 for(int elem = 0; elem < nels; elem++)
 {    
 TPZVec< TPZGeoEl * > filhos;
 TPZGeoEl * gel = gmesh->ElementVec()[elem];
 gel->Divide(filhos);
 }
 }
 std::ofstream outCube("QuadraticCube.vtk");
 TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outCube, false);
 
 std::cout << "FINISHED!!!" << std::endl;
 
 return 0;
 }
 */
