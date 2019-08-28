/*
 *  GenerateSquare.cpp
 *  SubstructEigen
 *
 *  Created by Philippe Devloo on 29/11/08.
 *  Copyright 2008 UNICAMP. All rights reserved.
 *
 */

#include "GenerateSquare.h"
#include "pzgengrid.h"
#include "tpzquadraticquad.h"
#include "tpzgeoblend.h"
#include "tpzarc3d.h"
#include "tpzcurvedtriangle.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"

TPZGeoMesh *Square(TPZGeoMesh &model, TConfig &cf)
{
	TPZVec<REAL> x0(2,0.),x1(2,1.);
	x1[1]=0.5;
	TPZVec<int> nx(2,2);
	nx[0]=2;
	nx[1]=1;
	TPZGenGrid gen(nx,x0,x1,1,0.);
	TPZGeoMesh *geo = new TPZGeoMesh(model);
	gen.Read(*geo);
	gen.SetBC(geo,0,cf.materialcondense);
	gen.SetBC(geo,1,cf.materialborder1);
	gen.SetBC(geo,2,cf.materialborder2);
	gen.SetBC(geo,3,cf.materialborder1);
	/*	TPZVec<REAL> one(3,0.), two(3,0.), three(3,0.);
	 two[0] = 0.5;
	 three[0] = 1.;
	 gen.SetBC(geo, one, two, cf.materialborder);
	 gen.SetBC(geo, two, three, cf.materialcondense);
	 */
	return geo;

}

TPZGeoMesh *SquareSingular(TPZGeoMesh &model, TConfig &cf)
{
	TPZVec<REAL> x0(2,0.),x1(2,1.);
	x1[1] = 0.5;
	TPZVec<int> nx(2,1);
	nx[0] = 2;
	nx[1] = 1;
	TPZGenGrid gen(nx,x0,x1,1,0.);
	TPZGeoMesh *geo = new TPZGeoMesh(model);
	gen.Read(*geo);
	gen.SetBC(geo,1,cf.materialborder1);
	gen.SetBC(geo,2,cf.materialborder2);
	gen.SetBC(geo,3,cf.materialborder1);
	TPZVec<REAL> one(3,0.), two(3,0.), three(3,0.);
	two[0] = 0.5;
	three[0] = 1.;
	//	gen.SetBC(geo, one, two, cf.materialborder);
	gen.SetBC(geo, two, three, cf.materialcondense);
	//	gen.SetBC(geo,3,cf.materialborder);
	return geo;

}

TPZGeoMesh *SquareQuad(TPZGeoMesh &model, TConfig &cf)
{
	REAL co[15][3];
	int nodeindex1[8] = {0,2,12,10,1,7,11,5};
	int i,j;
	for(i=0; i<5; i++)
	{
		for(j=0; j<3; j++)
		{
			co[j*5+i][0] = (i%5)*1./4.;
			co[j*5+i][1] = (j)*1./4.;
			co[j*5+i][2] = 0.;
		}
	}
	co[1][0] = 0.375;
	co[3][0] = 0.625;
	co[7][1] = 0.125;
	TPZGeoMesh *gmesh = new TPZGeoMesh(model);
	gmesh->NodeVec().Resize(15);
	for(i=0; i<15; i++)
	{
		TPZManVector<REAL,3> coord(3,0.);
		coord[0] = co[i][0];
		coord[1] = co[i][1];
		gmesh->NodeVec()[i].Initialize(coord,*gmesh);
	}
	int index;
	TPZManVector<int> nodeindexes(8,0);
	TPZGeoEl *gel[2];
	for(i=0; i<8; i++) nodeindexes[i] = nodeindex1[i];
	gel[0] = new TPZGeoElRefPattern<TPZQuadraticQuad> (nodeindexes, cf.materialid, *gmesh, index);
	for(i=0; i<8; i++) nodeindexes[i] = nodeindex1[i]+2;
	gel[1] = new TPZGeoElRefPattern<TPZQuadraticQuad> (nodeindexes, cf.materialid, *gmesh, index);
	gmesh->BuildConnectivity();
	TPZGeoElBC(gel[1],4,cf.materialcondense,*gmesh);
	TPZGeoElBC(gel[1],5,cf.materialborder1,*gmesh);
	TPZGeoElBC(gel[1],6,cf.materialborder2,*gmesh);
	TPZGeoElBC(gel[0],6,cf.materialborder2,*gmesh);
	TPZGeoElBC(gel[0],7,cf.materialborder1,*gmesh);
	return gmesh;

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
        th = i*M_PI/8.;
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

    Mesh->BuildConnectivity();

    return Mesh;
}

/*
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
        th = i*M_PI/8.;
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
*/

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
        th = i*M_PI/8.;
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

    Mesh->BuildConnectivity();

    return Mesh;
}

/*
TPZGeoMesh * CurvedTriangleMesh(TPZGeoMesh &model, TConfig &cf)
{
	double R = 1.;
	TPZGeoMesh * Mesh = new TPZGeoMesh(model);
	int Qnodes = 7;

	TPZVec < TPZVec <REAL> > NodesCoords(Qnodes);
	for(int i = 0; i < Qnodes; i++) NodesCoords[i].Resize(3);

	NodesCoords[0][0] = R;
	NodesCoords[0][1] = 0.;
	NodesCoords[0][2] = 0.;

	NodesCoords[1][0] = R*sqrt(2.);
	NodesCoords[1][1] = 0.;
	NodesCoords[1][2] = 0.;

	NodesCoords[2][0] = R;
	NodesCoords[2][1] = R;
	NodesCoords[2][2] = 0.;

	NodesCoords[3][0] = 0.;
	NodesCoords[3][1] = R;
	NodesCoords[3][2] = 0.;

	NodesCoords[4][0] = 0.;
	NodesCoords[4][1] = R*sqrt(2.);
	NodesCoords[4][2] = 0.;

	NodesCoords[5][0] = 0.9238795325112867*R;
	NodesCoords[5][1] = 0.3826834323650898*R;
	NodesCoords[5][2] = 0.;

	NodesCoords[6][0] = 0.3826834323650898*R;
	NodesCoords[6][1] = 0.9238795325112867*R;
	NodesCoords[6][2] = 0.;

	Mesh->NodeVec().Resize(Qnodes);
	TPZVec <TPZGeoNode> Node(Qnodes);
	for(int n = 0; n < Qnodes; n++)
	{
		Node[n].SetNodeId(n);
		Node[n].SetCoord(&NodesCoords[n][0]);
		Mesh->NodeVec()[n] = Node[n];
	}

	// Materials IDs
	int triangleMat = cf.materialid;
	int ArcMat = -1;
	int id = 0;

	TPZVec <int> Topol(3);

	//el0
	Topol[0] = 1; Topol[1] = 2; Topol[2] = 0;
	new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoTriangle> > (id,Topol,triangleMat,*Mesh);
	id++;

	//el1
	Topol[0] = 3; Topol[1] = 2; Topol[2] = 4;
	new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoTriangle> > (id,Topol,triangleMat,*Mesh);
	id++;

	//el2
	Topol[0] = 2; Topol[1] = 3; Topol[2] = 0;
	new TPZGeoElRefPattern<TPZCurvedTriangle > (id,Topol,triangleMat,*Mesh);
	id++;

	//el3
	Topol[0] = 1; Topol[1] = 2; Topol[2] = 5;
	new TPZGeoElRefPattern<TPZArc3D > (id,Topol,ArcMat,*Mesh);
	id++;

	//el4
	Topol[0] = 2; Topol[1] = 4; Topol[2] = 6;
	new TPZGeoElRefPattern<TPZArc3D > (id,Topol,ArcMat,*Mesh);

	Mesh->BuildBlendConnectivity();

	return Mesh;
}
*/

TPZGeoMesh * CurvedTriangleMesh(TPZGeoMesh &model, TConfig &cf)
{
	double R = 1.;
	TPZGeoMesh * Mesh = new TPZGeoMesh(model);
	int Qnodes = 8;

	TPZVec < TPZVec <REAL> > NodesCoords(Qnodes);
	for(int i = 0; i < Qnodes; i++) NodesCoords[i].Resize(3);

	NodesCoords[0][0] = R;
	NodesCoords[0][1] = 0.;
	NodesCoords[0][2] = 0.;

	NodesCoords[1][0] = R*sqrt(2.);
	NodesCoords[1][1] = 0.;
	NodesCoords[1][2] = 0.;

	NodesCoords[2][0] = R;
	NodesCoords[2][1] = R;
	NodesCoords[2][2] = 0.;

	NodesCoords[3][0] = 0.;
	NodesCoords[3][1] = R;
	NodesCoords[3][2] = 0.;

	NodesCoords[4][0] = 0.;
	NodesCoords[4][1] = R*sqrt(2.);
	NodesCoords[4][2] = 0.;

	NodesCoords[5][0] = 0.9238795325112867*R*sqrt(2.);
	NodesCoords[5][1] = 0.3826834323650898*R*sqrt(2.);
	NodesCoords[5][2] = 0.;

	NodesCoords[6][0] = 0.3826834323650898*R*sqrt(2.);
	NodesCoords[6][1] = 0.9238795325112867*R*sqrt(2.);
	NodesCoords[6][2] = 0.;

	NodesCoords[7][0] = R*sqrt(2.)/2.;
	NodesCoords[7][1] = R*sqrt(2.)/2.;
	NodesCoords[7][2] = 0.;

	Mesh->NodeVec().Resize(Qnodes);
	TPZVec <TPZGeoNode> Node(Qnodes);
	for(int n = 0; n < Qnodes; n++)
	{
		Node[n].SetNodeId(n);
		Node[n].SetCoord(&NodesCoords[n][0]);
		Mesh->NodeVec()[n] = Node[n];
	}

	// Materials IDs
	int triangleMat = cf.materialid;
	int ArcMat1 = cf.materialborder3;
	int line1 = cf.materialborder1;
	int ArcMat2 = 0;
	int line2 = cf.materialborder1;
	int id = 0;

	TPZVec <int> Topol(3);

	//el0
	Topol[0] = 1; Topol[1] = 2; Topol[2] = 0;
	new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoTriangle> > (id,Topol,triangleMat,*Mesh);
	id++;

	//el1
	Topol[0] = 3; Topol[1] = 2; Topol[2] = 4;
	new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoTriangle> > (id,Topol,triangleMat,*Mesh);
	id++;

	//el2
	Topol[0] = 2; Topol[1] = 3; Topol[2] = 0;
//	new TPZGeoElRefPattern<TPZCurvedTriangle > (id,Topol,triangleMat,*Mesh);
	new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoTriangle>  > (id,Topol,triangleMat,*Mesh);
	id++;

	//el3
	Topol[0] = 1; Topol[1] = 2; Topol[2] = 5;
	new TPZGeoElRefPattern<TPZArc3D> (id,Topol,ArcMat1,*Mesh);
	id++;

	//el4
	Topol[0] = 2; Topol[1] = 4; Topol[2] = 6;
	new TPZGeoElRefPattern<TPZArc3D> (id,Topol,ArcMat1,*Mesh);
	id++;

	//el5
	Topol.Resize(2);
	Topol[0] = 4; Topol[1] = 3;
	new TPZGeoElRefPattern<pzgeom::TPZGeoLinear> (id,Topol,line1,*Mesh);
	id++;

	//el6
	Topol.Resize(3);
	Topol[0] = 3; Topol[1] = 0; Topol[2] = 7;
	new TPZGeoElRefPattern<TPZArc3D> (id,Topol,ArcMat2,*Mesh);
	id++;

	//el7
	Topol.Resize(2);
	Topol[0] = 0; Topol[1] = 1;
	new TPZGeoElRefPattern<TPZGeoLinear> (id,Topol,line2,*Mesh);

	Mesh->BuildConnectivity();

	return Mesh;
}

TPZGeoMesh * QuadrilateralMesh(TPZGeoMesh &model, TConfig &cf)
{
	double R = 1.;
	TPZGeoMesh * Mesh = new TPZGeoMesh(model);
	int Qnodes = 10;

	TPZVec < TPZVec <REAL> > NodesCoords(Qnodes);
	for(int i = 0; i < Qnodes; i++) NodesCoords[i].Resize(3);

	NodesCoords[0][0] = R;
	NodesCoords[0][1] = 0.;
	NodesCoords[0][2] = 0.;

	NodesCoords[1][0] = R*sqrt(2.);
	NodesCoords[1][1] = 0.;
	NodesCoords[1][2] = 0.;

	NodesCoords[2][0] = R*sqrt(2.)/2.;
	NodesCoords[2][1] = R*sqrt(2.)/2.;
	NodesCoords[2][2] = 0.;

	NodesCoords[3][0] = R;
	NodesCoords[3][1] = R;
	NodesCoords[3][2] = 0.;

	NodesCoords[4][0] = 0.;
	NodesCoords[4][1] = R;
	NodesCoords[4][2] = 0.;

	NodesCoords[5][0] = 0.;
	NodesCoords[5][1] = R*sqrt(2.);
	NodesCoords[5][2] = 0.;

	NodesCoords[6][0] = 0.9238795325112867*R;
	NodesCoords[6][1] = 0.3826834323650898*R;
	NodesCoords[6][2] = 0.;

	NodesCoords[7][0] = 1.3065629648763766*R;
	NodesCoords[7][1] = 0.541196100146197*R;
	NodesCoords[7][2] = 0.;

	NodesCoords[8][0] = 0.3826834323650898*R;
	NodesCoords[8][1] = 0.9238795325112867*R;
	NodesCoords[8][2] = 0.;

	NodesCoords[9][0] = 0.541196100146197*R;
	NodesCoords[9][1] = 1.3065629648763766*R;
	NodesCoords[9][2] = 0.;

	Mesh->NodeVec().Resize(Qnodes);
	TPZVec <TPZGeoNode> Node(Qnodes);
	for(int n = 0; n < Qnodes; n++)
	{
		Node[n].SetNodeId(n);
		Node[n].SetCoord(&NodesCoords[n][0]);
		Mesh->NodeVec()[n] = Node[n];
	}

	// Materials IDs
	int QuadrilMat = cf.materialid;
	int lineMat1 = cf.materialborder1;
	int ArcMat1 = cf.materialborder3;
	int lineMat2 = cf.materialborder1;
	int ArcMat2 = 0;
	int id = 0;

	TPZVec <int> Topol(4);

	//el0
	Topol[0] = 0; Topol[1] = 1; Topol[2] = 3;  Topol[3] = 2;
	new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,QuadrilMat,*Mesh);
	id++;

	//el1
	Topol[0] = 2; Topol[1] = 3; Topol[2] = 5;  Topol[3] = 4;
	new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (id,Topol,QuadrilMat,*Mesh);
	id++;

	//el2
	Topol.Resize(2);
	Topol[0] = 0; Topol[1] = 1;
	new TPZGeoElRefPattern<TPZGeoLinear > (id,Topol,lineMat1,*Mesh);
	id++;

	//el3
	Topol.Resize(3);
	Topol[0] = 1; Topol[1] = 3; Topol[2] = 7;
	new TPZGeoElRefPattern<TPZArc3D > (id,Topol,ArcMat1,*Mesh);
	id++;

	//el4
	Topol[0] = 3; Topol[1] = 5; Topol[2] = 9;
	new TPZGeoElRefPattern<TPZArc3D > (id,Topol,ArcMat1,*Mesh);
	id++;

	//el5
	Topol.Resize(2);
	Topol[0] = 5; Topol[1] = 4;
	new TPZGeoElRefPattern<TPZGeoLinear > (id,Topol,lineMat2,*Mesh);
	id++;

	//el6
	Topol.Resize(3);
	Topol[0] = 4; Topol[1] = 2; Topol[2] = 8;
	new TPZGeoElRefPattern<TPZArc3D > (id,Topol,ArcMat2,*Mesh);
	id++;

	//el7
	Topol[0] = 2; Topol[1] = 0; Topol[2] = 6;
	new TPZGeoElRefPattern<TPZArc3D > (id,Topol,ArcMat2,*Mesh);
	id++;

//	Topol.Resize(1);
//	Topol[0] = 0;
//	new TPZGeoElRefPattern<pzgeom::TPZGeoPoint> (id,Topol,cf.materialcorner,*Mesh);

	Mesh->BuildConnectivity();

	return Mesh;
}
