/*
 *  BrazilianTestGeoMesh.cpp
 *  ElastoPlasticModels
 *
 *  Created by Diogo Cecilio on 10/5/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
//#include "TPZTensor.h"
#include "BrazilianTestGeoMesh.h"
#include "pzgeotetrahedra.h"
#include "pzgeopyramid.h"
#include "TPZRefPatternTools.h"

using namespace pzgeom;


BrazilianTestGeoMesh::BrazilianTestGeoMesh()
{
}

BrazilianTestGeoMesh::~BrazilianTestGeoMesh()
{
}

TPZGeoMesh * BrazilianTestGeoMesh::GeoMesh(int h)
{
	int Qnodes = 8;
	
	REAL val= 5.3033;
	REAL H=30.;
	
	
	TPZGeoMesh * gMesh = new TPZGeoMesh;
	gMesh->SetMaxNodeId(Qnodes-1);
	gRefDBase.InitializeRefPatterns();
	gMesh->NodeVec().Resize(Qnodes);
	
	TPZVec <int> TopolCube(8);
	TPZVec <int> TopolLine(2);
	TPZVec <int> TopolQuad(4);
	TPZVec <int> TopolPoint(1);
	
	TPZVec<TPZGeoNode> Node(Qnodes);
	
	
	/////ZERO////////                      
	Node[0].SetNodeId(0);
	Node[0].SetCoord(0 , -val);//coord X
	Node[0].SetCoord(1 , -val);//coord Y 
	Node[0].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[0] = Node[0];
	
	///////UM////////
	Node[1].SetNodeId(1);
	Node[1].SetCoord(0 ,  val);//coord X
	Node[1].SetCoord(1 , -val);//coord Y
	Node[1].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[1] = Node[1];
	
	///////DOIS////////
	Node[2].SetNodeId(2);
	Node[2].SetCoord(0 ,  val);//coord X
	Node[2].SetCoord(1 ,  val);//coord Y
	Node[2].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[2] = Node[2];
	
	///////TRES////////
	Node[3].SetNodeId(3);
	Node[3].SetCoord(0 ,  -val);//coord X
	Node[3].SetCoord(1 ,   val);//coord Y
	Node[3].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[3] = Node[3];
	
	///////QUATRO////////
	Node[4].SetNodeId(4);
	Node[4].SetCoord(0 ,  -val);//coord X
	Node[4].SetCoord(1 ,  -val);//coord Y
	Node[4].SetCoord(2,    H );//coord Z
	gMesh->NodeVec()[4] = Node[4];
	
	///////CINCO////////
	Node[5].SetNodeId(5);
	Node[5].SetCoord(0 ,   val);//coord X
	Node[5].SetCoord(1 ,  -val);//coord Y
	Node[5].SetCoord(2,    H );//coord Z
	gMesh->NodeVec()[5] = Node[5];
	
	///////SEIS////////
	Node[6].SetNodeId(6);
	Node[6].SetCoord(0 ,  val);//coord X
	Node[6].SetCoord(1 ,  val);//coord Y
	Node[6].SetCoord(2,    H );//coord Z
	gMesh->NodeVec()[6] = Node[6];
	
	///////SETE////////
	Node[7].SetNodeId(7);
	Node[7].SetCoord(0 ,  -val);//coord X
	Node[7].SetCoord(1 ,   val);//coord Y
	Node[7].SetCoord(2,    H );//coord Z
	gMesh->NodeVec()[7] = Node[7];
	
	//GERA CUBO///
	TopolCube[0] = 0;	TopolCube[1] = 1;	TopolCube[4] = 4;	TopolCube[5] = 5;
	TopolCube[2] = 2; 	TopolCube[3] = 3;   TopolCube[6] = 6; 	TopolCube[7] = 7;
	new TPZGeoElRefPattern<  TPZGeoCube > (0,TopolCube,1,*gMesh);//id=0 // matind=1
	
	//GERA QUADRILATERO INFERIOR///
	TopolQuad[0] = 0;	TopolQuad[1] = 1;	TopolQuad[2] = 2;	TopolQuad[3] = 3;
	new TPZGeoElRefPattern< TPZGeoQuad >  (1,TopolQuad,-1,*gMesh);//id=1 // matind=-1
	
	//GERA QUADRILATERO SUPERIOR///
	TopolQuad[0] = 4;	TopolQuad[1] = 5;	TopolQuad[2] = 6;	TopolQuad[3] = 7;
	new TPZGeoElRefPattern< TPZGeoQuad >  (2,TopolQuad,-2,*gMesh);//id=2 // matind=-2
	
	
	///QUADRILATEROS PARA O TESTE DE CONFINAMENTO
	/////////
	
	//MATIDEX -11 //ID 7
	TopolQuad[0] = 0;	TopolQuad[1] = 1;	TopolQuad[2] = 5;	TopolQuad[3] = 4;
	new TPZGeoElRefPattern< TPZGeoQuad >  (7,TopolQuad,-11,*gMesh);
	
	//MATIDEX -12 //ID 8
	TopolQuad[0] = 1;	TopolQuad[1] = 2;	TopolQuad[2] = 6;	TopolQuad[3] = 5;
	new TPZGeoElRefPattern< TPZGeoQuad >  (8,TopolQuad,-12,*gMesh);
	
	//MATIDEX -13 //ID 9
	TopolQuad[0] = 3;	TopolQuad[1] = 2;	TopolQuad[2] = 6;	TopolQuad[3] = 7;
	new TPZGeoElRefPattern< TPZGeoQuad >  (9,TopolQuad,-13,*gMesh);
	
	//MATIDEX -14 //ID 10
	TopolQuad[0] = 0;	TopolQuad[1] = 3;	TopolQuad[2] = 7;	TopolQuad[3] = 4;
	new TPZGeoElRefPattern< TPZGeoQuad >  (10,TopolQuad,-14,*gMesh);
	
	
	
	//GERA OBJETO PONTO NO 0 ///
	TopolPoint[0] = 0;
	new TPZGeoElRefPattern< TPZGeoPoint >  (3,TopolPoint,-3,*gMesh);
	
	//GERA OBJETO PONTO NO 1 ///
	TopolPoint[0] = 1;
	new TPZGeoElRefPattern< TPZGeoPoint >  (4,TopolPoint,-4,*gMesh);
	
	//GERA OBJETO PONTO NO 2 ///
	TopolPoint[0] = 2;
	new TPZGeoElRefPattern< TPZGeoPoint >  (5,TopolPoint,-5,*gMesh);
	
	
	gMesh->BuildConnectivity();
	
	for(int stage = 0; stage < h; stage++)
	{
		int nel = gMesh->NElements();
		for ( int iref = 0; iref < nel; iref++ )
		{
			TPZVec<TPZGeoEl*> filhos;
			TPZGeoEl * gelP1 = gMesh->ElementVec()[iref];
			if(!gelP1) continue;
			gelP1->Divide(filhos);
		}
	}
	
	return gMesh;
}

TPZGeoMesh *BrazilianTestGeoMesh::GeoBlendMesh(int h)
{
	int Qnodes = 48;	
	REAL H=30.;
	
	
	
	
	//gMesh->SetMaxNodeId(Qnodes-1);
	TPZGeoMesh *gMesh = new TPZGeoMesh;
	gRefDBase.InitializeRefPatterns();
//	gMesh->InitializeRefPatterns();
	//	gMesh->ImportRefPattern ( );
	gMesh->NodeVec().Resize(Qnodes);
	
	TPZVec <int> TopolCube(8);
	TPZVec <int> TopolPrism(6);
	TPZVec <int> TopolArc(3);
	TPZVec <int> TopolLine(2);
	TPZVec <int> TopolQuad(4);
	TPZVec <int> TopolPoint(1);
	TPZVec <int> TopolTringle(3);
	
	TPZVec<TPZGeoNode> Node(Qnodes);
	
	
	/////ZERO////////
	Node[0].SetNodeId(0);
	Node[0].SetCoord(0 , -0.75);//coord X
	Node[0].SetCoord(1 , -0.75);//coord Y
	Node[0].SetCoord(2,    0.);//coord Z
	gMesh->NodeVec()[0] = Node[0];
	
	///////UM////////
	Node[1].SetNodeId(1);
	Node[1].SetCoord(0 ,  0.75);//coord X
	Node[1].SetCoord(1 , -0.75);//coord Y
	Node[1].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[1] = Node[1];
	
	
	///////DOIS////////
	Node[2].SetNodeId(2);
	Node[2].SetCoord(0 ,  0.75);//coord X
	Node[2].SetCoord(1 ,  0.75);//coord Y
	Node[2].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[2] = Node[2];
	
	///////TRES////////
	Node[3].SetNodeId(3);
	Node[3].SetCoord(0 ,  -0.75);//coord X
	Node[3].SetCoord(1 ,   0.75);//coord Y
	Node[3].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[3] = Node[3];
	
	///////QUATRO////////
	Node[4].SetNodeId(4);
	Node[4].SetCoord(0 ,  -5.3033);//coord X
	Node[4].SetCoord(1 ,  -5.3033);//coord Y
	Node[4].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[4] = Node[4];
	
	///////CINCO////////
	Node[5].SetNodeId(5);
	Node[5].SetCoord(0 ,  -2.87013);//coord X
	Node[5].SetCoord(1 , -6.92909);//coord Y
	Node[5].SetCoord(2,    0.);//coord Z
	gMesh->NodeVec()[5] = Node[5];
	
	///////SEIS////////
	Node[6].SetNodeId(6);
	Node[6].SetCoord(0 ,  -0.75);//coord X
	Node[6].SetCoord(1 ,  -7.46);//coord Y
	Node[6].SetCoord(2,    0.);//coord Z
	gMesh->NodeVec()[6] = Node[6];
	
	///////SETE////////
	Node[7].SetNodeId(7);
	Node[7].SetCoord(0 ,  0.);//coord X
	Node[7].SetCoord(1 ,   -7.5);//coord Y
	Node[7].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[7] = Node[7];
	
	///////OITO////////
	Node[8].SetNodeId(8);
	Node[8].SetCoord(0 ,  0.75);//coord X
	Node[8].SetCoord(1 ,  -7.46);//coord Y
	Node[8].SetCoord(2,   0. );//coord Z
	gMesh->NodeVec()[8] = Node[8];
	
	///////NOVE////////
	Node[9].SetNodeId(9);
	Node[9].SetCoord(0 ,  2.87013);//coord X
	Node[9].SetCoord(1 ,  -6.92909);//coord Y
	Node[9].SetCoord(2,   0. );//coord Z
	gMesh->NodeVec()[9] = Node[9];
	
	///////DEZ////////
	Node[10].SetNodeId(10);
	Node[10].SetCoord(0 ,  5.3033);//coord X
	Node[10].SetCoord(1 ,  -5.3033);//coord Y
	Node[10].SetCoord(2,   0.);//coord Z
	gMesh->NodeVec()[10] = Node[10];
	
	///////ONZE////////
	Node[11].SetNodeId(11);
	Node[11].SetCoord(0 ,  7.5);//coord X
	Node[11].SetCoord(1 ,  0.);//coord Y
	Node[11].SetCoord(2,   0. );//coord Z
	gMesh->NodeVec()[11] = Node[11];
	
	
	///////DOZE////////
	Node[12].SetNodeId(12);
	Node[12].SetCoord(0 ,  5.3033);//coord X
	Node[12].SetCoord(1 ,  5.3033);//coord Y
	Node[12].SetCoord(2,   0.);//coord Z
	gMesh->NodeVec()[12] = Node[12];
	
	///////TREZE////////
	Node[13].SetNodeId(13);
	Node[13].SetCoord(0 ,  2.87013);//coord X
	Node[13].SetCoord(1 ,  6.92909);//coord Y
	Node[13].SetCoord(2,   0. );//coord Z
	gMesh->NodeVec()[13] = Node[13];
	
	///////QUATORZE////////
	Node[14].SetNodeId(14);
	Node[14].SetCoord(0 ,  0.75);//coord X
	Node[14].SetCoord(1 ,  7.46);//coord Y
	Node[14].SetCoord(2,   0 );//coord Z
	gMesh->NodeVec()[14] = Node[14];
	
	///////QUINZE////////
	Node[15].SetNodeId(15);
	Node[15].SetCoord(0 ,  0.);//coord X
	Node[15].SetCoord(1 ,  7.5);//coord Y
	Node[15].SetCoord(2,   0. );//coord Z
	gMesh->NodeVec()[15] = Node[15];
	
	
	///////DEZESSEIS////////
	Node[16].SetNodeId(16);
	Node[16].SetCoord(0 ,  -0.75);//coord X
	Node[16].SetCoord(1 ,  7.46);//coord Y
	Node[16].SetCoord(2,   0. );//coord Z
	gMesh->NodeVec()[16] = Node[16];
	
	///////DEZESSETE////////
	Node[17].SetNodeId(17);
	Node[17].SetCoord(0 ,  -2.87013);//coord X
	Node[17].SetCoord(1 ,  6.92909);//coord Y
	Node[17].SetCoord(2,   0. );//coord Z
	gMesh->NodeVec()[17] = Node[17];
	
	
	///////DEZOITO////////
	Node[18].SetNodeId(18);
	Node[18].SetCoord(0 ,  -5.3033);//coord X
	Node[18].SetCoord(1 ,  5.3033);//coord Y
	Node[18].SetCoord(2,   0.);//coord Z
	gMesh->NodeVec()[18] = Node[18];
	
	///////DEZENOVE////////
	Node[19].SetNodeId(19);
	Node[19].SetCoord(0 ,  -7.5);//coord X
	Node[19].SetCoord(1 ,  0.);//coord Y
	Node[19].SetCoord(2,   0.);//coord Z
	gMesh->NodeVec()[19] = Node[19];
	
	///////VINTE////////
	Node[20].SetNodeId(20);
	Node[20].SetCoord(0 ,  -0.75);//coord X
	Node[20].SetCoord(1 ,  -8.46);//coord Y
	Node[20].SetCoord(2,   0.);//coord Z
	gMesh->NodeVec()[20] = Node[20];
	
	///////VINTE E UM////////
	Node[21].SetNodeId(21);
	Node[21].SetCoord(0 ,  0.75);//coord X
	Node[21].SetCoord(1 ,  -8.46);//coord Y
	Node[21].SetCoord(2,   0.);//coord Z
	gMesh->NodeVec()[21] = Node[21];
	
	///////VINTE E DOIS////////
	Node[22].SetNodeId(22);
	Node[22].SetCoord(0 ,  0.75);//coord X
	Node[22].SetCoord(1 ,  8.46);//coord Y
	Node[22].SetCoord(2,   0.);//coord Z
	gMesh->NodeVec()[22] = Node[22];
	
	///////VINTE E TRES////////
	Node[23].SetNodeId(23);
	Node[23].SetCoord(0 ,  -0.75);//coord X
	Node[23].SetCoord(1 ,  8.46);//coord Y
	Node[23].SetCoord(2,   0.);//coord Z
	gMesh->NodeVec()[23] = Node[23];
	
	
	/////VINTE E QUATRO////////
	Node[24].SetNodeId(24);
	Node[24].SetCoord(0 , -0.75);//coord X
	Node[24].SetCoord(1 , -0.75);//coord Y
	Node[24].SetCoord(2,    H);//coord Z
	gMesh->NodeVec()[24] = Node[24];
	
	///////VINTE E CINCO////////
	Node[25].SetNodeId(25);
	Node[25].SetCoord(0 ,  0.75);//coord X
	Node[25].SetCoord(1 , -0.75);//coord Y
	Node[25].SetCoord(2,    H );//coord Z
	gMesh->NodeVec()[25] = Node[25];
	
	
	///////VINTE E SEIS////////
	Node[26].SetNodeId(26);
	Node[26].SetCoord(0 ,  0.75);//coord X
	Node[26].SetCoord(1 ,  0.75);//coord Y
	Node[26].SetCoord(2,    H );//coord Z
	gMesh->NodeVec()[26] = Node[26];
	
	///////VINTE E SETE////////
	Node[27].SetNodeId(27);
	Node[27].SetCoord(0 ,  -0.75);//coord X
	Node[27].SetCoord(1 ,   0.75);//coord Y
	Node[27].SetCoord(2,    H );//coord Z
	gMesh->NodeVec()[27] = Node[27];
	
	///////VINTE E OITO////////
	Node[28].SetNodeId(28);
	Node[28].SetCoord(0 ,  -5.3033);//coord X
	Node[28].SetCoord(1 ,  -5.3033);//coord Y
	Node[28].SetCoord(2,    H );//coord Z
	gMesh->NodeVec()[28] = Node[28];
	
	///////VINTE E NOVE////////
	Node[29].SetNodeId(29);
	Node[29].SetCoord(0 ,  -2.87013);//coord X
	Node[29].SetCoord(1 , -6.92909);//coord Y
	Node[29].SetCoord(2,    H);//coord Z
	gMesh->NodeVec()[29] = Node[29];
	
	///////SEIS////////
	Node[30].SetNodeId(30);
	Node[30].SetCoord(0 ,  -0.75);//coord X
	Node[30].SetCoord(1 ,  -7.46);//coord Y
	Node[30].SetCoord(2,    H);//coord Z
	gMesh->NodeVec()[30] = Node[30];
	
	///////SETE////////
	Node[31].SetNodeId(31);
	Node[31].SetCoord(0 ,  0.);//coord X
	Node[31].SetCoord(1 ,   -7.5);//coord Y
	Node[31].SetCoord(2,    H );//coord Z
	gMesh->NodeVec()[31] = Node[31];
	
	///////OITO////////
	Node[32].SetNodeId(32);
	Node[32].SetCoord(0 ,  0.75);//coord X
	Node[32].SetCoord(1 ,  -7.46);//coord Y
	Node[32].SetCoord(2,   H );//coord Z
	gMesh->NodeVec()[32] = Node[32];
	
	///////NOVE////////
	Node[33].SetNodeId(33);
	Node[33].SetCoord(0 ,  2.87013);//coord X
	Node[33].SetCoord(1 ,  -6.92909);//coord Y
	Node[33].SetCoord(2,   H );//coord Z
	gMesh->NodeVec()[33] = Node[33];
	
	///////DEZ////////
	Node[34].SetNodeId(34);
	Node[34].SetCoord(0 ,  5.3033);//coord X
	Node[34].SetCoord(1 ,  -5.3033);//coord Y
	Node[34].SetCoord(2,   H);//coord Z
	gMesh->NodeVec()[34] = Node[34];
	
	///////ONZE////////
	Node[35].SetNodeId(35);
	Node[35].SetCoord(0 ,  7.5);//coord X
	Node[35].SetCoord(1 ,  0.);//coord Y
	Node[35].SetCoord(2,   H );//coord Z
	gMesh->NodeVec()[35] = Node[35];
	
	
	///////DOZE////////
	Node[36].SetNodeId(36);
	Node[36].SetCoord(0 ,  5.3033);//coord X
	Node[36].SetCoord(1 ,  5.3033);//coord Y
	Node[36].SetCoord(2,   H);//coord Z
	gMesh->NodeVec()[36] = Node[36];
	
	///////TREZE////////
	Node[37].SetNodeId(37);
	Node[37].SetCoord(0 ,  2.87013);//coord X
	Node[37].SetCoord(1 ,  6.92909);//coord Y
	Node[37].SetCoord(2,   H );//coord Z
	gMesh->NodeVec()[37] = Node[37];
	
	///////QUATORZE////////
	Node[38].SetNodeId(38);
	Node[38].SetCoord(0 ,  0.75);//coord X
	Node[38].SetCoord(1 ,  7.46);//coord Y
	Node[38].SetCoord(2,   H );//coord Z
	gMesh->NodeVec()[38] = Node[38];
	
	///////QUINZE////////
	Node[39].SetNodeId(39);
	Node[39].SetCoord(0 ,  0.);//coord X
	Node[39].SetCoord(1 ,  7.5);//coord Y
	Node[39].SetCoord(2,  H );//coord Z
	gMesh->NodeVec()[39] = Node[39];
	
	
	///////DEZESSEIS////////
	Node[40].SetNodeId(40);
	Node[40].SetCoord(0 ,  -0.75);//coord X
	Node[40].SetCoord(1 ,  7.46);//coord Y
	Node[40].SetCoord(2,   H );//coord Z
	gMesh->NodeVec()[40] = Node[40];
	
	///////DEZESSETE////////
	Node[41].SetNodeId(41);
	Node[41].SetCoord(0 ,  -2.87013);//coord X
	Node[41].SetCoord(1 ,  6.92909);//coord Y
	Node[41].SetCoord(2,   H );//coord Z
	gMesh->NodeVec()[41] = Node[41];
	
	
	///////DEZOITO////////
	Node[42].SetNodeId(42);
	Node[42].SetCoord(0 ,  -5.3033);//coord X
	Node[42].SetCoord(1 ,  5.3033);//coord Y
	Node[42].SetCoord(2,   H);//coord Z
	gMesh->NodeVec()[42] = Node[42];
	
	///////DEZENOVE////////
	Node[43].SetNodeId(43);
	Node[43].SetCoord(0 ,  -7.5);//coord X
	Node[43].SetCoord(1 ,  0.);//coord Y
	Node[43].SetCoord(2,  H);//coord Z
	gMesh->NodeVec()[43] = Node[43];
	
	///////VINTE////////
	Node[44].SetNodeId(44);
	Node[44].SetCoord(0 ,  -0.75);//coord X
	Node[44].SetCoord(1 ,  -8.46);//coord Y
	Node[44].SetCoord(2,   H);//coord Z
	gMesh->NodeVec()[44] = Node[44];
	
	///////VINTE E UM////////
	Node[45].SetNodeId(45);
	Node[45].SetCoord(0 ,  0.75);//coord X
	Node[45].SetCoord(1 ,  -8.46);//coord Y
	Node[45].SetCoord(2,   H);//coord Z
	gMesh->NodeVec()[45] = Node[45];
	
	///////VINTE E DOIS////////
	Node[46].SetNodeId(46);
	Node[46].SetCoord(0 ,  0.75);//coord X
	Node[46].SetCoord(1 ,  8.46);//coord Y
	Node[46].SetCoord(2,  H);//coord Z
	gMesh->NodeVec()[46] = Node[46];
	
	///////VINTE E TRES////////
	Node[47].SetNodeId(47);
	Node[47].SetCoord(0 ,  -0.75);//coord X
	Node[47].SetCoord(1 ,  8.46);//coord Y
	Node[47].SetCoord(2,   H);//coord Z
	gMesh->NodeVec()[47] = Node[47];
	
	
	//ELEMENTO 0///
	TopolCube[0] = 0;	TopolCube[1] = 1;	TopolCube[4] = 24;	TopolCube[5] = 25;
	TopolCube[2] = 2; 	TopolCube[3] = 3;   TopolCube[6] = 26; 	TopolCube[7] = 27;//////////////////0
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoCube> > (0,TopolCube,1,*gMesh);//id=0 // matind=1
	
	//ELEMENTO 1///
	TopolCube[0] = 1;	TopolCube[1] = 10;	TopolCube[4] = 25;	TopolCube[5] = 34;
	TopolCube[2] = 12; 	TopolCube[3] = 2;   TopolCube[6] = 36; 	TopolCube[7] = 26;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoCube> > (1,TopolCube,1,*gMesh);//id=0 // matind=1
	
	//ELEMENTO 2///
	TopolPrism[0] = 2;	TopolPrism[1] = 12;	  	TopolPrism[3] = 26;  TopolPrism[4] = 36;
	TopolPrism[2] = 14; 	                    TopolPrism[5] = 38;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoPrism> > (2,TopolPrism,1,*gMesh);//id=0 // matind=1
	
	//ELEMENTO 3///
	TopolCube[0] = 3;	TopolCube[1] = 2;	TopolCube[4] = 27;	TopolCube[5] = 26;
	TopolCube[2] = 14; 	TopolCube[3] = 16;   TopolCube[6] = 38; 	TopolCube[7] = 40;///////////////////3
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoCube> > (3,TopolCube,1,*gMesh);//id=0 // matind=1
	
	//ELEMENTO 4///
	TopolPrism[0] = 3;	TopolPrism[1] = 16;	  	TopolPrism[3] = 27; TopolPrism[4] = 40;
	TopolPrism[2] = 18; 	                    TopolPrism[5] = 42;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoPrism> > (4,TopolPrism,1,*gMesh);//id=0 // matind=1
	
	//ELEMENTO 5///
	TopolCube[0] = 0;	TopolCube[1] = 3;	TopolCube[4] = 24;	TopolCube[5] = 27;
	TopolCube[2] = 18; 	TopolCube[3] = 4;   TopolCube[6] = 42; 	TopolCube[7] = 28;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoCube> > (5,TopolCube,1,*gMesh);//id=0 // matind=1
	
	//ELEMENTO 6///
	TopolPrism[0] = 6;	TopolPrism[1] = 0;	  	TopolPrism[3] = 30; TopolPrism[4] = 24;
	TopolPrism[2] = 4; 	                    TopolPrism[5] = 28;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoPrism> > (6,TopolPrism,1,*gMesh);//id=0 // matind=1
	
	//ELEMENTO 7///
	TopolCube[0] = 6;	TopolCube[1] = 8;	TopolCube[4] = 30;	TopolCube[5] = 32;
	TopolCube[2] = 1; 	TopolCube[3] = 0;   TopolCube[6] = 25; 	TopolCube[7] = 24;/////////////////////////7
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoCube> > (7,TopolCube,1,*gMesh);//id=0 // matind=1
	
	//ELEMENTO 8///
	TopolPrism[0] = 8;	TopolPrism[1] = 10;	  	TopolPrism[3] = 32; TopolPrism[4] = 34;
	TopolPrism[2] = 1; 	                    TopolPrism[5] = 25;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoPrism> > (8,TopolPrism,1,*gMesh);//id=0 // matind=1
	
	
	///ARCOS BASE////////
    // ARC 1 //
	TopolArc[0] = 18;
	TopolArc[1] = 4;
	TopolArc[2] = 19;
	new TPZGeoElRefPattern< TPZArc3D >  (11,TopolArc,-1,*gMesh);//id=1 // matind=-1
	
	// ARC 2 //
	TopolArc[0] = 4;
	TopolArc[1] = 6;
	TopolArc[2] = 5;
	new TPZGeoElRefPattern< TPZArc3D > (12,TopolArc,-1,*gMesh);//id=2 // matind=-2
	
	// ARC 3 //
	TopolArc[0] = 6;
	TopolArc[1] = 8;
	TopolArc[2] = 7;
	new TPZGeoElRefPattern< TPZArc3D > (13,TopolArc,-1,*gMesh);//id=3 // matind=-3
	
	// ARC 4 //
	TopolArc[0] = 8;
	TopolArc[1] = 10;
	TopolArc[2] = 9;
	new TPZGeoElRefPattern< TPZArc3D > (14,TopolArc,-1,*gMesh);//id=4 // matind=-4
	
	// ARC 5 //
	TopolArc[0] = 10;
	TopolArc[1] = 12;
	TopolArc[2] = 11;
	new TPZGeoElRefPattern< TPZArc3D >  (15,TopolArc,-1,*gMesh);//id=1 // matind=-1
	
	// ARC 6 //
	TopolArc[0] = 12;
	TopolArc[1] = 14;
	TopolArc[2] = 13;
	new TPZGeoElRefPattern< TPZArc3D > (16,TopolArc,-1,*gMesh);//id=2 // matind=-2
	
	// ARC 7 //
	TopolArc[0] = 14;
	TopolArc[1] = 16;
	TopolArc[2] = 15;
	new TPZGeoElRefPattern< TPZArc3D > (17,TopolArc,-1,*gMesh);//id=3 // matind=-3
	
	// ARC 8 //
	TopolArc[0] = 16;
	TopolArc[1] = 18;
	TopolArc[2] = 17;
	new TPZGeoElRefPattern< TPZArc3D > (18,TopolArc,-1,*gMesh);//id=4 // matind=-4
	
	
	///ARCOS TOPO////////
    // ARC 1 //
	TopolArc[0] = 42;
	TopolArc[1] = 28;
	TopolArc[2] = 43;
	new TPZGeoElRefPattern< TPZArc3D >  (19,TopolArc,-1,*gMesh);//id=1 // matind=-1
	
	// ARC 2 //
	TopolArc[0] = 28;
	TopolArc[1] = 30;
	TopolArc[2] = 29;
	new TPZGeoElRefPattern< TPZArc3D > (20,TopolArc,-1,*gMesh);//id=2 // matind=-2
	
	// ARC 3 //
	TopolArc[0] = 30;
	TopolArc[1] = 32;
	TopolArc[2] = 31;
	new TPZGeoElRefPattern< TPZArc3D > (21,TopolArc,-1,*gMesh);//id=3 // matind=-3
	
	// ARC 4 //
	TopolArc[0] = 32;
	TopolArc[1] = 34;
	TopolArc[2] = 33;
	new TPZGeoElRefPattern< TPZArc3D > (22,TopolArc,-1,*gMesh);//id=4 // matind=-4
	
	// ARC 5 //
	TopolArc[0] = 34;
	TopolArc[1] = 36;
	TopolArc[2] = 35;
	new TPZGeoElRefPattern< TPZArc3D >  (23,TopolArc,-1,*gMesh);//id=1 // matind=-1
	
	// ARC 6 //
	TopolArc[0] = 36;
	TopolArc[1] = 38;
	TopolArc[2] = 37;
	new TPZGeoElRefPattern< TPZArc3D > (24,TopolArc,-1,*gMesh);//id=2 // matind=-2
	
	// ARC 7 //
	TopolArc[0] = 38;
	TopolArc[1] = 40;
	TopolArc[2] = 39;
	new TPZGeoElRefPattern< TPZArc3D > (25,TopolArc,-1,*gMesh);//id=3 // matind=-3
	
	// ARC 8 //
	TopolArc[0] = 40;
	TopolArc[1] = 42;
	TopolArc[2] = 41;
	new TPZGeoElRefPattern< TPZArc3D > (26,TopolArc,-1,*gMesh);//id=4 // matind=-4
	
	
	//GERA OBJETO PONTO NO 6 ///
	TopolPoint[0] = 6;
	new TPZGeoElRefPattern< TPZGeoPoint >  (29,TopolPoint,-2,*gMesh);
	
	//GERA OBJETO PONTO NO 8 ///
	TopolPoint[0] = 8;
	new TPZGeoElRefPattern< TPZGeoPoint >  (30,TopolPoint,-3,*gMesh);
	
	//	//GERA OBJETO PONTO NO 30 ///
	//	TopolPoint[0] = 30;
	//	new TPZGeoElRefPattern< TPZGeoPoint >  (31,TopolPoint,-4,*gMesh);
	
	
	
	//ELEMENTO 13///
	TopolQuad[0] = 6;	TopolQuad[1] = 8;
	TopolQuad[2] = 32; 	TopolQuad[3] = 30;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoQuad> > (27,TopolQuad,-5,*gMesh);
	
	//ELEMENTO 14///
	TopolQuad[0] = 16;	TopolQuad[1] = 14;
	TopolQuad[2] = 38; 	TopolQuad[3] = 40;
	new TPZGeoElRefPattern<TPZGeoBlend< TPZGeoQuad> >(28,TopolQuad,-6,*gMesh);
	
	///FACE Z = 30
	
	//ELEMENTO 32///
	TopolQuad[0] = 24;	TopolQuad[1] = 25;
	TopolQuad[2] = 26; 	TopolQuad[3] = 27;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoQuad> > (32,TopolQuad,-7,*gMesh);
	
	//ELEMENTO 33///
	TopolQuad[0] = 25;	TopolQuad[1] = 34;
	TopolQuad[2] = 36; 	TopolQuad[3] = 26;
	new TPZGeoElRefPattern<TPZGeoBlend< TPZGeoQuad> >(33,TopolQuad,-7,*gMesh);
	
	//ELEMENTO 34///
	TopolTringle[0] = 26;	TopolTringle[1] = 36; TopolTringle[2] = 38;
	new TPZGeoElRefPattern<TPZGeoBlend< TPZGeoTriangle> >(34,TopolTringle,-7,*gMesh);
	
	//ELEMENTO 35///
	TopolQuad[0] = 27;	TopolQuad[1] = 26;
	TopolQuad[2] = 38; 	TopolQuad[3] = 40;
	new TPZGeoElRefPattern<TPZGeoBlend< TPZGeoQuad> >(35,TopolQuad,-7,*gMesh);
	
	//ELEMENTO 36///
	TopolTringle[0] = 27;	TopolTringle[1] = 40; TopolTringle[2] = 42;
	new TPZGeoElRefPattern<TPZGeoBlend< TPZGeoTriangle> >(36,TopolTringle,-7,*gMesh);
	
	//ELEMENTO 37///
	TopolQuad[0] = 24;	TopolQuad[1] = 27;
	TopolQuad[2] = 42; 	TopolQuad[3] = 28;
	new TPZGeoElRefPattern<TPZGeoBlend< TPZGeoQuad> >(37,TopolQuad,-7,*gMesh);
	
	//ELEMENTO 38///
	TopolTringle[0] = 28;	TopolTringle[1] = 30; TopolTringle[2] = 24;
	new TPZGeoElRefPattern<TPZGeoBlend< TPZGeoTriangle> >(38,TopolTringle,-7,*gMesh);
	
	//ELEMENTO 39///
	TopolQuad[0] = 30;	TopolQuad[1] = 32;
	TopolQuad[2] = 25; 	TopolQuad[3] = 24;
	new TPZGeoElRefPattern<TPZGeoBlend< TPZGeoQuad> >(39,TopolQuad,-7,*gMesh);
	
	//ELEMENTO 40///
	TopolTringle[0] = 32;	TopolTringle[1] = 34; TopolTringle[2] = 25;
	new TPZGeoElRefPattern<TPZGeoBlend< TPZGeoTriangle> >(40,TopolTringle,-7,*gMesh);
	
	///FACE Z = 0
	
	//ELEMENTO 41///
	TopolQuad[0] = 0;	TopolQuad[1] = 1;
	TopolQuad[2] = 2; 	TopolQuad[3] = 3;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoQuad> > (41,TopolQuad,-8,*gMesh);
	
	//ELEMENTO 42///
	TopolQuad[0] = 1;	TopolQuad[1] = 10;
	TopolQuad[2] = 12; 	TopolQuad[3] = 2;
	new TPZGeoElRefPattern<TPZGeoBlend< TPZGeoQuad> >(42,TopolQuad,-8,*gMesh);
	
	//ELEMENTO 43///
	TopolTringle[0] = 2;	TopolTringle[1] = 12; TopolTringle[2] = 14;
	new TPZGeoElRefPattern<TPZGeoBlend< TPZGeoTriangle> >(43,TopolTringle,-7,*gMesh);
	
	//ELEMENTO 44///
	TopolQuad[0] = 3;	TopolQuad[1] = 2;
	TopolQuad[2] = 14; 	TopolQuad[3] = 16;
	new TPZGeoElRefPattern<TPZGeoBlend< TPZGeoQuad> >(44,TopolQuad,-8,*gMesh);
	
	//ELEMENTO 45///
	TopolTringle[0] = 3;	TopolTringle[1] = 16; TopolTringle[2] = 18;
	new TPZGeoElRefPattern<TPZGeoBlend< TPZGeoTriangle> >(45,TopolTringle,-8,*gMesh);
	
	//ELEMENTO 46///
	TopolQuad[0] = 0;	TopolQuad[1] = 3;
	TopolQuad[2] = 18; 	TopolQuad[3] = 4;
	new TPZGeoElRefPattern<TPZGeoBlend< TPZGeoQuad> >(46,TopolQuad,-8,*gMesh);
	
	//ELEMENTO 47///
	TopolTringle[0] = 4;	TopolTringle[1] = 6; TopolTringle[2] = 0;
	new TPZGeoElRefPattern<TPZGeoBlend< TPZGeoTriangle> >(47,TopolTringle,-8,*gMesh);
	
	//ELEMENTO 48///
	TopolQuad[0] = 6;	TopolQuad[1] = 8;
	TopolQuad[2] = 1; 	TopolQuad[3] = 0;
	new TPZGeoElRefPattern<TPZGeoBlend< TPZGeoQuad> >(48,TopolQuad,-8,*gMesh);
	
	//ELEMENTO 49///
	TopolTringle[0] = 8;	TopolTringle[1] = 10; TopolTringle[2] = 1;
	new TPZGeoElRefPattern<TPZGeoBlend< TPZGeoTriangle> >(49,TopolTringle,-8,*gMesh);
	
	//LATERAIS
	
	//ELEMENTO 50///
	TopolQuad[0] = 8;	TopolQuad[1] = 10;
	TopolQuad[2] = 34; 	TopolQuad[3] = 32;
	new TPZGeoElRefPattern<TPZGeoBlend< TPZGeoQuad> >(50,TopolQuad,-9,*gMesh);	//ELEMENTO 42///
	//ELEMENTO 51///
	TopolQuad[0] = 12;	TopolQuad[1] = 10;
	TopolQuad[2] = 34; 	TopolQuad[3] = 36;
	new TPZGeoElRefPattern<TPZGeoBlend< TPZGeoQuad> >(51,TopolQuad,-9,*gMesh);	//ELEMENTO 42///
	//ELEMENTO 52///
	TopolQuad[0] = 14;	TopolQuad[1] = 12;
	TopolQuad[2] = 36; 	TopolQuad[3] = 38;
	new TPZGeoElRefPattern<TPZGeoBlend< TPZGeoQuad> >(52,TopolQuad,-9,*gMesh);	//ELEMENTO 42///
	//ELEMENTO 53///
	TopolQuad[0] = 18;	TopolQuad[1] = 16;
	TopolQuad[2] = 40; 	TopolQuad[3] = 42;
	new TPZGeoElRefPattern<TPZGeoBlend< TPZGeoQuad> >(53,TopolQuad,-9,*gMesh);
	//ELEMENTO 54///
	TopolQuad[0] = 4;	TopolQuad[1] = 28;
	TopolQuad[2] = 42; 	TopolQuad[3] = 18;
	new TPZGeoElRefPattern<TPZGeoBlend< TPZGeoQuad> >(54,TopolQuad,-9,*gMesh);
	//ELEMENTO 55///
	TopolQuad[0] = 4;	TopolQuad[1] = 6;
	TopolQuad[2] = 30; 	TopolQuad[3] = 28;
	new TPZGeoElRefPattern<TPZGeoBlend< TPZGeoQuad> >(55,TopolQuad,-9,*gMesh);
	
	gMesh->BuildConnectivity();
	
	//	int nel = gMesh->NElements();
	//	for(int i=0;i<nel;i++)
	//	{
	//		TPZGeoEl * gelP1 = gMesh->ElementVec()[i];
	//		/**return 1 if the element has subelements along side*/
	//		//virtual int HasSubElement() = 0;
	//		if(gelP1->HasSubElement() == 0 /* dont have subelements*/)
	//		{
	//			
	//		}
	//			
	//	}
	
	std::set<int> matids;
	matids.insert(-17);
	matids.insert(-18);
	//matids.insert(-21);
	//matids.insert(-22);
	//matids.insert(-23);
	//matids.insert(-24);
	//	 for(int i = 0; i < h; i++)
	//	 {
	//			int nel = gMesh->NElements();
	//			for(int eel = 0; eel < nel; eel++)	
	//			{
	//				TPZGeoEl *elemento = gMesh->ElementVec()[eel];
	//				TPZRefPattern::RefineDirectional(elemento, matids);
	//			}
	//	}
	
	
	for(int stage = 0; stage < h; stage++)
	{
		int nel = gMesh->NElements();
		for ( int iref = 0; iref < nel; iref++ )
		{
			TPZVec<TPZGeoEl*> tatara;
			TPZGeoEl * gelP1 = gMesh->ElementVec()[iref];
			if(!gelP1) continue;
			gelP1->Divide(tatara);
			
		}
	}
	//		int nel1 = gMesh->NElements();
	//		for ( int iref = 0; iref < nel1; iref++ )
	//		{
	//			TPZVec<TPZGeoEl*> netos;
	//			TPZGeoEl * gelP1 = gMesh->ElementVec()[iref];
	//			if(gelP1->TypeName() == "Prism")
	//			{
	//				if(!gelP1) continue;
	//				gelP1->Divide(netos);
	//			}
	//			
	//		}
	//		int nel2 = gMesh->NElements();
	//		for ( int iref = 0; iref < nel2; iref++ )
	//		{
	//			TPZVec<TPZGeoEl*> bisnetos;
	//			TPZGeoEl * gelP1 = gMesh->ElementVec()[iref];
	//			if(gelP1->Id() == 0 || gelP1->Id() == 3 || gelP1->Id() == 7 )
	//			{
	//				if(!gelP1) continue;
	//				gelP1->Divide(bisnetos);
	//			}
	//			
	//			
	//		}
	
	
	
	
	//	std::ofstream out("MalhaBraziliamTestDEBUG.vtk");
	//	TPZVTKGeoMesh::PrintGMeshVTK(gMesh, out);
	
	return gMesh;
}


void BrazilianTestGeoMesh::TransformBlendToLinearMesh(TPZGeoMesh *newlinearmesh, int h)
{
	TPZGeoMesh *MESHblend = new TPZGeoMesh;
	MESHblend = BrazilianTestGeoMesh::GeoBlendMesh(h);
	
	newlinearmesh->NodeVec() = MESHblend->NodeVec();
	int nel = MESHblend->NElements();
	
	int gelMat, gelId = 0;
	
	for(int el = 0; el < nel; el++)
	{
		TPZGeoEl * gel = MESHblend->ElementVec()[el];
		if(gel->HasSubElement())
		{
			continue;
		}
		MElementType gelType = gel->Type();
		
		int nnodes = gel->NNodes();
		TPZVec<int> Topol(nnodes);
		for(int n = 0; n < nnodes; n++)
		{
			Topol[n] = gel->NodeIndex(n);
		}
		
		switch(gelType)
		{
			case(EPoint) :
			{
				gelMat = gel->MaterialId();
				new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (gelId,Topol,gelMat,*newlinearmesh);
				//
				gelId++;
				break;
			}
			case (EOned) :
			{
				gelMat = gel->MaterialId();
				new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (gelId,Topol,gelMat,*newlinearmesh);
				//
				gelId++;
				break;
			}
			case (ETriangle) :
			{
				gelMat = gel->MaterialId();
				new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (gelId,Topol,gelMat,*newlinearmesh);
				//
				gelId++;
				break;
			}
			case (EQuadrilateral) :
			{
				gelMat = gel->MaterialId();
				new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (gelId,Topol,gelMat,*newlinearmesh);
				//
				gelId++;
				break;
			}
			case (ETetraedro) :
			{
				gelMat = gel->MaterialId();
				new TPZGeoElRefPattern< pzgeom::TPZGeoTetrahedra > (gelId,Topol,gelMat,*newlinearmesh);
				//
				gelId++;
				break;
			}
			case (EPrisma) :
			{
				gelMat = gel->MaterialId();
				new TPZGeoElRefPattern< pzgeom::TPZGeoPrism > (gelId,Topol,gelMat,*newlinearmesh);
				//
				gelId++;
				break;
			}
			case (EPiramide) :
			{
				gelMat = gel->MaterialId();
				new TPZGeoElRefPattern< pzgeom::TPZGeoPyramid > (gelId,Topol,gelMat,*newlinearmesh);
				//
				gelId++;
				break;
			}
			case (ECube) :
			{
				gelMat = gel->MaterialId();
				new TPZGeoElRefPattern< pzgeom::TPZGeoCube > (gelId,Topol,gelMat,*newlinearmesh);
				//
				gelId++;
				break;
			}
			default :
			{
				std::cout << "Elemento nao achado no metodo " << __PRETTY_FUNCTION__ << std::endl;
				DebugStop();
			}
				
		}
	}
	
	newlinearmesh->BuildConnectivity();
	
	
}

void BrazilianTestGeoMesh::TransformBlendToLinearMesh2(TPZGeoMesh *newlinearmesh, int h)
{
	TPZGeoMesh *MESHblend = new TPZGeoMesh;
	//	MESHblend = BrazilianTestGeoMesh::GeoBlendMesh(h);
	int dir = 0;
	MESHblend = BrazilianTestGeoMesh::TwoDMesh(h,dir);
	newlinearmesh->NodeVec() = MESHblend->NodeVec();
	int nel = MESHblend->NElements();
	
	int gelMat, gelId = 0;
	
	for(int el = 0; el < nel; el++)
	{
		TPZGeoEl * gel = MESHblend->ElementVec()[el];
		if(gel->HasSubElement())
		{
			continue;
		}
		//bool isOnZ0 = false;
		//		for(int n = 0; n < gel->NNodes(); n++)
		//		{
		//			double coordZ = gel->NodePtr(n)->Coord(2);
		//			if(fabs(coordZ) < /*3.741*/1.E-5 )
		//			{
		//				isOnZ0 = true;
		//				break;
		//			}
		//		}
		//		if(isOnZ0 == false)
		//		{
		//			continue;
		//		}
		MElementType gelType = gel->Type();
		
		int nnodes = gel->NNodes();
		TPZVec<int> Topol(nnodes);
		for(int n = 0; n < nnodes; n++)
		{
			Topol[n] = gel->NodeIndex(n);
		}
		
		switch(gelType)
		{
			case(EPoint) :
			{
				gelMat = gel->MaterialId();
				new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (gelId,Topol,gelMat,*newlinearmesh);
				//
				gelId++;
				break;
			}
			case (EOned) :
			{
				gelMat = gel->MaterialId();
				new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (gelId,Topol,gelMat,*newlinearmesh);
				//
				gelId++;
				break;
			}
			case (ETriangle) :
			{
				gelMat = gel->MaterialId();
				new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (gelId,Topol,gelMat,*newlinearmesh);
				//
				gelId++;
				break;
			}
			case (EQuadrilateral) :
			{
				gelMat = gel->MaterialId();
				new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (gelId,Topol,gelMat,*newlinearmesh);
				//
				gelId++;
				break;
			}
			case (ETetraedro) :
			{
				gelMat = gel->MaterialId();
				new TPZGeoElRefPattern< pzgeom::TPZGeoTetrahedra > (gelId,Topol,gelMat,*newlinearmesh);
				//
				gelId++;
				break;
			}
			case (EPrisma) :
			{
				gelMat = gel->MaterialId();
				new TPZGeoElRefPattern< pzgeom::TPZGeoPrism > (gelId,Topol,gelMat,*newlinearmesh);
				//
				gelId++;
				break;
			}
			case (EPiramide) :
			{
				gelMat = gel->MaterialId();
				new TPZGeoElRefPattern< pzgeom::TPZGeoPyramid > (gelId,Topol,gelMat,*newlinearmesh);
				//
				gelId++;
				break;
			}
			case (ECube) :
			{
				gelMat = gel->MaterialId();
				new TPZGeoElRefPattern< pzgeom::TPZGeoCube > (gelId,Topol,gelMat,*newlinearmesh);
				//
				gelId++;
				break;
			}
			default :
			{
				std::cout << "Elemento nao achado no metodo " << __PRETTY_FUNCTION__ << std::endl;
				DebugStop();
			}
				
		}
	}
	
	
	//TPZVec<int> topopoint(1);
	//	topopoint[0]=4096;
	//	new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (200000,topopoint,-4,*newlinearmesh);
	newlinearmesh->BuildConnectivity();
	
	
}

void BrazilianTestGeoMesh::ReadMesh(TPZGeoMesh &mesh)
{
	int numnodes=-1;
	int numelements=-1;
	
	string FileName;
	FileName = "MalhaFinal.txt"; 
	
	{
		bool countnodes = false;
		bool countelements = false;
		
		ifstream read(FileName.c_str());
		
		while(read)
		{
			char buf[1024];
			read.getline(buf, 1024);
			std::string str(buf);
			if(str == "POINTS 33280 float")
			{
				countnodes = true;
			}
			
			
			if(str == "end coordinates") 
			{
				countnodes = false;
			}
			
			
			
			if(countnodes)
			{
				numnodes++;
			}
			
			
			if(str == "CELLS 4736 38016") 
			{
				countelements = true;
			}
			if(str == "end elements") 
			{
				countelements = false;
			}
			if(countelements)
			{
				numelements++;
			}
		}
		
	}
	
	
	mesh.NodeVec().Resize(numnodes);
	
	TPZVec <int> TopolTetra(4);
	
	const int Qnodes = numnodes;
	TPZVec <TPZGeoNode> Node(Qnodes);
	
	ifstream read;
	read.open(FileName.c_str());
	
	double nodecoordX , nodecoordY , nodecoordZ ;
	
	char buf[1024];
	read.getline(buf, 1024);
	read.getline(buf, 1024);
	std::string str(buf);
	int in;
	for(in=0; in<numnodes; in++)
	{ 
		read >> nodecoordX;
		read >> nodecoordY;
		read >> nodecoordZ;
		Node[in].SetNodeId(in);
		Node[in].SetCoord(0,nodecoordX);
		Node[in].SetCoord(1,nodecoordY);
		Node[in].SetCoord(2,nodecoordZ);
		mesh.NodeVec()[in] = Node[in];
	}
	
	{
		
		read.close();
		read.open(FileName.c_str());
		
		
		
		int l , m = numnodes+5;
		for(l=0; l<m; l++)
		{
			read.getline(buf, 1024);
		}
		
		TPZVec <int> TopolPoint(1);
		TopolPoint[0]=22527;
		new TPZGeoElRefPattern< TPZGeoPoint > (100000,TopolPoint,-3,mesh);
		TopolPoint[0]=26184;
		new TPZGeoElRefPattern< TPZGeoPoint > (100001,TopolPoint,-4,mesh);
		TopolPoint[0]=24282;
		new TPZGeoElRefPattern< TPZGeoPoint > (100002,TopolPoint,-5,mesh);
		
		
		int topol=0;
		TPZVec <int> TopolCube(8);
		TPZVec <int> TopolPrism(6);
		TPZVec <int> TopolQuad(4);
		
		int el;
		
		for(el=0; el<numelements; el++)
		{
			
			read >> topol;
			if(topol==4)
			{
				int topo1,topo2,topo3,topo4;
				read >> topo1;
				read >> topo2;
				read >> topo3;
				read >> topo4;
				TopolQuad[0] = topo1; 
				TopolQuad[1] = topo2; 
				TopolQuad[2] = topo3; 
				TopolQuad[3] = topo4;
				TopolQuad[0]--;
				TopolQuad[1]--;
				TopolQuad[2]--;
				TopolQuad[3]--;
				
				double Ycoord = Node[topo1].Coord(1);
				if(Ycoord < 0.)
				{
					new TPZGeoElRefPattern< TPZGeoQuad > (el,TopolQuad,-2,mesh);
				}
				if(Ycoord > 0.)
				{
					new TPZGeoElRefPattern< TPZGeoQuad > (el,TopolQuad,-1,mesh);
				}
			}
			
			if(topol==6)
			{
				read >> TopolPrism[0]; 
				read >> TopolPrism[1]; 
				read >> TopolPrism[2]; 
				read >> TopolPrism[3];
				read >> TopolPrism[4]; 
				read >> TopolPrism[5];
				TopolPrism[0]--;
				TopolPrism[1]--;
				TopolPrism[2]--;
				TopolPrism[3]--;
				TopolPrism[4]--;
				TopolPrism[5]--;
				new TPZGeoElRefPattern< TPZGeoPrism > (el,TopolPrism,1,mesh);
				
			}
			if(topol==8)
			{
				read >> TopolCube[0]; 
				read >> TopolCube[1]; 
				read >> TopolCube[2]; 
				read >> TopolCube[3];
				read >> TopolCube[4]; 
				read >> TopolCube[5];
				read >> TopolCube[6];
				read >> TopolCube[7];
				TopolCube[0]--; 
				TopolCube[1]--; 
				TopolCube[2]--; 
				TopolCube[3]--;
				TopolCube[4]--; 
				TopolCube[5]--;
				TopolCube[6]--;
				TopolCube[7]--; 
				new TPZGeoElRefPattern< TPZGeoCube > (el,TopolCube,1,mesh);
				
			}
			
		}
		mesh.BuildConnectivity();
	}
	
	ofstream arg("BLENDERA.txt");
	mesh.Print(arg);
	
	std::ofstream out("BLEND.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(&mesh, out, true);
}

TPZGeoMesh * BrazilianTestGeoMesh::TwoDMesh(int h, int dir)
{
	int Qnodes = 24;	
	
	
	TPZGeoMesh *gMesh = new TPZGeoMesh;
	
	//gMesh->InitializeRefPatterns();
	gRefDBase.InitializeRefPatterns();
	
	//gMesh->ImportRefPattern ( );
	gMesh->NodeVec().Resize(Qnodes);
	
	TPZVec <int> TopolArc(3);
	TPZVec <int> TopolLine(2);
	TPZVec <int> TopolQuad(4);
	TPZVec <int> TopolPoint(1);
	TPZVec <int> TopolTringle(3);
	
	TPZVec<TPZGeoNode> Node(Qnodes);
	
	
	/////ZERO////////
	Node[0].SetNodeId(0);
	Node[0].SetCoord(0 , -0.75);//coord X
	Node[0].SetCoord(1 , -0.75);//coord Y
	Node[0].SetCoord(2,    0.);//coord Z
	gMesh->NodeVec()[0] = Node[0];
	
	///////UM////////
	Node[1].SetNodeId(1);
	Node[1].SetCoord(0 ,  0.75);//coord X
	Node[1].SetCoord(1 , -0.75);//coord Y
	Node[1].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[1] = Node[1];
	
	
	///////DOIS////////
	Node[2].SetNodeId(2);
	Node[2].SetCoord(0 ,  0.75);//coord X
	Node[2].SetCoord(1 ,  0.75);//coord Y
	Node[2].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[2] = Node[2];
	
	///////TRES////////
	Node[3].SetNodeId(3);
	Node[3].SetCoord(0 ,  -0.75);//coord X
	Node[3].SetCoord(1 ,   0.75);//coord Y
	Node[3].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[3] = Node[3];
	
	///////QUATRO////////
	Node[4].SetNodeId(4);
	Node[4].SetCoord(0 ,  -5.3033);//coord X
	Node[4].SetCoord(1 ,  -5.3033);//coord Y
	Node[4].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[4] = Node[4];
	
	///////CINCO////////
	Node[5].SetNodeId(5);
	Node[5].SetCoord(0 ,  -2.87013);//coord X
	Node[5].SetCoord(1 , -6.92909);//coord Y
	Node[5].SetCoord(2,    0.);//coord Z
	gMesh->NodeVec()[5] = Node[5];
	
	///////SEIS////////
	Node[6].SetNodeId(6);
	Node[6].SetCoord(0 ,  -0.75);//coord X
	Node[6].SetCoord(1 ,  -7.46);//coord Y
	Node[6].SetCoord(2,    0.);//coord Z
	gMesh->NodeVec()[6] = Node[6];
	
	///////SETE////////
	Node[7].SetNodeId(7);
	Node[7].SetCoord(0 ,  0.);//coord X
	Node[7].SetCoord(1 ,   -7.5);//coord Y
	Node[7].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[7] = Node[7];
	
	///////OITO////////
	Node[8].SetNodeId(8);
	Node[8].SetCoord(0 ,  0.75);//coord X
	Node[8].SetCoord(1 ,  -7.46);//coord Y
	Node[8].SetCoord(2,   0. );//coord Z
	gMesh->NodeVec()[8] = Node[8];
	
	///////NOVE////////
	Node[9].SetNodeId(9);
	Node[9].SetCoord(0 ,  2.87013);//coord X
	Node[9].SetCoord(1 ,  -6.92909);//coord Y
	Node[9].SetCoord(2,   0. );//coord Z
	gMesh->NodeVec()[9] = Node[9];
	
	///////DEZ////////
	Node[10].SetNodeId(10);
	Node[10].SetCoord(0 ,  5.3033);//coord X
	Node[10].SetCoord(1 ,  -5.3033);//coord Y
	Node[10].SetCoord(2,   0.);//coord Z
	gMesh->NodeVec()[10] = Node[10];
	
	///////ONZE////////
	Node[11].SetNodeId(11);
	Node[11].SetCoord(0 ,  7.5);//coord X
	Node[11].SetCoord(1 ,  0.);//coord Y
	Node[11].SetCoord(2,   0. );//coord Z
	gMesh->NodeVec()[11] = Node[11];
	
	
	///////DOZE////////
	Node[12].SetNodeId(12);
	Node[12].SetCoord(0 ,  5.3033);//coord X
	Node[12].SetCoord(1 ,  5.3033);//coord Y
	Node[12].SetCoord(2,   0.);//coord Z
	gMesh->NodeVec()[12] = Node[12];
	
	///////TREZE////////
	Node[13].SetNodeId(13);
	Node[13].SetCoord(0 ,  2.87013);//coord X
	Node[13].SetCoord(1 ,  6.92909);//coord Y
	Node[13].SetCoord(2,   0. );//coord Z
	gMesh->NodeVec()[13] = Node[13];
	
	///////QUATORZE////////
	Node[14].SetNodeId(14);
	Node[14].SetCoord(0 ,  0.75);//coord X
	Node[14].SetCoord(1 ,  7.46);//coord Y
	Node[14].SetCoord(2,   0 );//coord Z
	gMesh->NodeVec()[14] = Node[14];
	
	///////QUINZE////////
	Node[15].SetNodeId(15);
	Node[15].SetCoord(0 ,  0.);//coord X
	Node[15].SetCoord(1 ,  7.5);//coord Y
	Node[15].SetCoord(2,   0. );//coord Z
	gMesh->NodeVec()[15] = Node[15];
	
	
	///////DEZESSEIS////////
	Node[16].SetNodeId(16);
	Node[16].SetCoord(0 ,  -0.75);//coord X
	Node[16].SetCoord(1 ,  7.46);//coord Y
	Node[16].SetCoord(2,   0. );//coord Z
	gMesh->NodeVec()[16] = Node[16];
	
	///////DEZESSETE////////
	Node[17].SetNodeId(17);
	Node[17].SetCoord(0 ,  -2.87013);//coord X
	Node[17].SetCoord(1 ,  6.92909);//coord Y
	Node[17].SetCoord(2,   0. );//coord Z
	gMesh->NodeVec()[17] = Node[17];
	
	
	///////DEZOITO////////
	Node[18].SetNodeId(18);
	Node[18].SetCoord(0 ,  -5.3033);//coord X
	Node[18].SetCoord(1 ,  5.3033);//coord Y
	Node[18].SetCoord(2,   0.);//coord Z
	gMesh->NodeVec()[18] = Node[18];
	
	///////DEZENOVE////////
	Node[19].SetNodeId(19);
	Node[19].SetCoord(0 ,  -7.5);//coord X
	Node[19].SetCoord(1 ,  0.);//coord Y
	Node[19].SetCoord(2,   0.);//coord Z
	gMesh->NodeVec()[19] = Node[19];
	
	///////VINTE////////
	Node[20].SetNodeId(20);
	Node[20].SetCoord(0 ,  -0.75);//coord X
	Node[20].SetCoord(1 ,  -8.46);//coord Y
	Node[20].SetCoord(2,   0.);//coord Z
	gMesh->NodeVec()[20] = Node[20];
	
	///////VINTE E UM////////
	Node[21].SetNodeId(21);
	Node[21].SetCoord(0 ,  0.75);//coord X
	Node[21].SetCoord(1 ,  -8.46);//coord Y
	Node[21].SetCoord(2,   0.);//coord Z
	gMesh->NodeVec()[21] = Node[21];
	
	///////VINTE E DOIS////////
	Node[22].SetNodeId(22);
	Node[22].SetCoord(0 ,  0.75);//coord X
	Node[22].SetCoord(1 ,  8.46);//coord Y
	Node[22].SetCoord(2,   0.);//coord Z
	gMesh->NodeVec()[22] = Node[22];
	
	///////VINTE E TRES////////
	Node[23].SetNodeId(23);
	Node[23].SetCoord(0 ,  -0.75);//coord X
	Node[23].SetCoord(1 ,  8.46);//coord Y
	Node[23].SetCoord(2,   0.);//coord Z
	gMesh->NodeVec()[23] = Node[23];
	
	//ELEMENTO 0///
	TopolQuad[0] = 0;	TopolQuad[1] = 1;
	TopolQuad[2] = 2; 	TopolQuad[3] = 3;  
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoQuad> > (0,TopolQuad,1,*gMesh);//id=0 // matind=1
	
	//ELEMENTO 1///
	TopolQuad[0] = 1;	TopolQuad[1] = 10;
	TopolQuad[2] = 12; 	TopolQuad[3] = 2; 
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoQuad> > (1,TopolQuad,1,*gMesh);//id=1 // matind=1
	
	//ELEMENTO 2///
	TopolTringle[0] = 2;	TopolTringle[1] = 12;
	TopolTringle[2] = 14; 	                    
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoTriangle> > (2,TopolTringle,1,*gMesh);//id=2 // matind=1
	
	//ELEMENTO 3///
	TopolQuad[0] = 3;	TopolQuad[1] = 2;
	TopolQuad[2] = 14; 	TopolQuad[3] = 16; 
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoQuad> > (3,TopolQuad,1,*gMesh);//id=3 // matind=1
	
	//ELEMENTO 4///
	TopolTringle[0] = 3;	TopolTringle[1] = 16;
	TopolTringle[2] = 18; 
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoTriangle> > (4,TopolTringle,1,*gMesh);//id=0 // matind=1
	
	//ELEMENTO 5///
	TopolQuad[0] = 0;	TopolQuad[1] = 3;
	TopolQuad[2] = 18; 	TopolQuad[3] = 4;   
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoQuad> > (5,TopolQuad,1,*gMesh);//id=0 // matind=1
	
	//ELEMENTO 6///
	TopolTringle[0] = 6;	TopolTringle[1] = 0;
	TopolTringle[2] = 4; 
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoTriangle> > (6,TopolTringle,1,*gMesh);//id=0 // matind=1
	
	//ELEMENTO 7///
	TopolQuad[0] = 6;	TopolQuad[1] = 8;
	TopolQuad[2] = 1; 	TopolQuad[3] = 0;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoQuad> > (7,TopolQuad,1,*gMesh);//id=0 // matind=1
	
	//ELEMENTO 8///
	TopolTringle[0] = 8;	TopolTringle[1] = 10;
	TopolTringle[2] = 1;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoTriangle> > (8,TopolTringle,1,*gMesh);//id=0 // matind=1
	
	// ARC 1 //
	TopolArc[0] = 18;
	TopolArc[1] = 4;
	TopolArc[2] = 19;
	new TPZGeoElRefPattern< TPZArc3D >  (9,TopolArc,-1,*gMesh);//id=1 // matind=-1
	
	// ARC 2 //
	TopolArc[0] = 4;//ESTE
	TopolArc[1] = 6;
	TopolArc[2] = 5;
	new TPZGeoElRefPattern< TPZArc3D > (10,TopolArc,-1,*gMesh);//id=2 // matind=-2
	
	// ARC 3 //
	TopolArc[0] = 6;
	TopolArc[1] = 8;
	TopolArc[2] = 7;
	new TPZGeoElRefPattern< TPZArc3D > (11,TopolArc,-1,*gMesh);//id=3 // matind=-3
	
	// ARC 4 //
	TopolArc[0] = 8;//ESTE
	TopolArc[1] = 10;
	TopolArc[2] = 9;
	new TPZGeoElRefPattern< TPZArc3D > (12,TopolArc,-1,*gMesh);//id=4 // matind=-4
	
	// ARC 5 //
	TopolArc[0] = 10;
	TopolArc[1] = 12;
	TopolArc[2] = 11;
	new TPZGeoElRefPattern< TPZArc3D >  (13,TopolArc,-1,*gMesh);//id=1 // matind=-1
	
	// ARC 6 //
	TopolArc[0] = 12;//ESTE
	TopolArc[1] = 14;
	TopolArc[2] = 13;
	new TPZGeoElRefPattern< TPZArc3D > (14,TopolArc,-1,*gMesh);//id=2 // matind=-2
	
	// ARC 7 //
	TopolArc[0] = 14;
	TopolArc[1] = 16;
	TopolArc[2] = 15;
	new TPZGeoElRefPattern< TPZArc3D > (15,TopolArc,-1,*gMesh);//id=3 // matind=-3
	
	// ARC 8 //
	TopolArc[0] = 16;//ESTE
	TopolArc[1] = 18;
	TopolArc[2] = 17;
	new TPZGeoElRefPattern< TPZArc3D > (16,TopolArc,-1,*gMesh);//id=4 // matind=-4
	
	//DownLine
	TopolLine[0]=6;
	TopolLine[1]=8;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (17, TopolLine,-2,*gMesh);
	
	//UperLine
	TopolLine[0]=14;
	TopolLine[1]=16;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (18, TopolLine,-3,*gMesh);
	
	
	
	//EDGE
	TopolLine[0]=4;
	TopolLine[1]=6;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (22, TopolLine,-7,*gMesh);
	
	//EDGE
	TopolLine[0]=8;
	TopolLine[1]=10;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (23, TopolLine,-8,*gMesh);
	//EDGE
	TopolLine[0]=12;
	TopolLine[1]=14;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (24, TopolLine,-9,*gMesh);
	
	//EDGE
	TopolLine[0]=16;
	TopolLine[1]=18;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (25, TopolLine,-10,*gMesh);
	
	
	
	//CENTER
	TopolLine[0]=0;
	TopolLine[1]=1;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (26, TopolLine,-11,*gMesh);
	
	//CENTER
	TopolLine[0]=1;
	TopolLine[1]=2;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (27, TopolLine,-12,*gMesh);
	////CENTER
	TopolLine[0]=2;
	TopolLine[1]=3;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (28, TopolLine,-13,*gMesh);
	
	////CENTER
	TopolLine[0]=3;
	TopolLine[1]=0;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (29, TopolLine,-14,*gMesh);
	
	
	//DIAGONALS
	TopolLine[0]=3;
	TopolLine[1]=18;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (30, TopolLine,-15,*gMesh);
	
	//DIAGONALS
	TopolLine[0]=2;
	TopolLine[1]=12;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (31, TopolLine,-16,*gMesh);
	////DIAGONALS
	TopolLine[0]=0;
	TopolLine[1]=4;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (32, TopolLine,-17,*gMesh);
	
	////DIAGONALS
	TopolLine[0]=1;
	TopolLine[1]=10;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (33, TopolLine,-18,*gMesh);
	
	//VERTICALS
	TopolLine[0]=3;
	TopolLine[1]=16;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (34, TopolLine,-19,*gMesh);
	
	//VERTICALS
	TopolLine[0]=2;
	TopolLine[1]=14;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (35, TopolLine,-20,*gMesh);
	////VERTICALS
	TopolLine[0]=1;
	TopolLine[1]=8;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (36, TopolLine,-21,*gMesh);
	
	////VERTICALS
	TopolLine[0]=0;
	TopolLine[1]=6;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (37, TopolLine,-22,*gMesh);
	
	
	
//	//left down point
//	TopolPoint[0]=6;
//	new TPZGeoElRefPattern< TPZGeoPoint > (19, TopolPoint,-4,*gMesh);
//	
//	//right down point
//	TopolPoint[0]=8;
//	new TPZGeoElRefPattern<  TPZGeoPoint  > (20, TopolPoint,-5,*gMesh);
//	
//	//upper right point
//	TopolPoint[0]=14;
//	new TPZGeoElRefPattern<  TPZGeoPoint  > (21, TopolPoint,-6,*gMesh);
    
  //CENTER
	TopolPoint[0]=0;
	new TPZGeoElRefPattern< TPZGeoPoint > (38, TopolPoint,-4,*gMesh);
	
	//CENTER
	TopolPoint[0]=1;
	new TPZGeoElRefPattern<  TPZGeoPoint  > (39, TopolPoint,-5,*gMesh);
	
	//CENTER
	TopolPoint[0]=2;
	new TPZGeoElRefPattern<  TPZGeoPoint  > (40, TopolPoint,-6,*gMesh);
    
    //CENTER
	TopolPoint[0]=3;
	new TPZGeoElRefPattern<  TPZGeoPoint  > (41, TopolPoint,-7,*gMesh);
	
	
	
	gMesh->BuildConnectivity();
	
	for(int stage = 0; stage < h; stage++)
	{
		int nel = gMesh->NElements();
		for ( int iref = 0; iref < nel; iref++ )
		{
			TPZVec<TPZGeoEl*> tatara;
			TPZGeoEl * gelP1 = gMesh->ElementVec()[iref];
			if(!gelP1) continue;
			gelP1->Divide(tatara);
			
		}
	}
	
	std::set<int> matids;
	//EDGES
	matids.insert(-2);
	matids.insert(-3);
	//	matids.insert(-7);
	//	matids.insert(-8);
	//	matids.insert(-9);
	//	matids.insert(-10);
	
	//	//CENTER
	//	matids.insert(-11);
	//	matids.insert(-12);
	//	matids.insert(-13);
	//	matids.insert(-14);
	
	//DIAGONALS
	
	//	matids.insert(-15);
	//	matids.insert(-16);
	//	matids.insert(-17);
	//	matids.insert(-18);
	
	
	//VERTICALS
	//	matids.insert(-19);
	//	matids.insert(-20);
	//	matids.insert(-21);
	//	matids.insert(-22);
	
	for(int i = 0; i < dir; i++)
	{
		int nel = gMesh->NElements();
		for(int eel = 0; eel < nel; eel++)	
		{
			TPZGeoEl *elemento = gMesh->ElementVec()[eel];
			TPZRefPatternTools::RefineDirectional(elemento,matids);
//			TPZRefPattern::RefineDirectional(elemento, matids);
		}
	}
	
	
	return gMesh;
	
	
}

TPZGeoMesh * BrazilianTestGeoMesh::TwoDMeshII(int h, int dir)
{
	int Qnodes = 25;	
	
	
	TPZGeoMesh *gMesh = new TPZGeoMesh;
	
	//gMesh->InitializeRefPatterns();
	gRefDBase.InitializeRefPatterns();
	
	//gMesh->ImportRefPattern ( );
	gMesh->NodeVec().Resize(Qnodes);
	
	TPZVec <int> TopolArc(3);
	TPZVec <int> TopolLine(2);
	TPZVec <int> TopolQuad(4);
	TPZVec <int> TopolPoint(1);
	TPZVec <int> TopolTringle(3);
	
	TPZVec<TPZGeoNode> Node(Qnodes);
	
	
    /////ZERO////////
	Node[24].SetNodeId(24);
	Node[24].SetCoord(0 ,0.*10);//coord X
	Node[24].SetCoord(1 ,0.*10);//coord Y
	Node[24].SetCoord(2, 0.);//coord Z
	gMesh->NodeVec()[24] = Node[24];
    
	/////ZERO////////
	Node[0].SetNodeId(0);
	Node[0].SetCoord(0 , -0.75*10);//coord X
	Node[0].SetCoord(1 , -0.75*10);//coord Y
	Node[0].SetCoord(2,    0.);//coord Z
	gMesh->NodeVec()[0] = Node[0];
	
	///////UM////////
	Node[1].SetNodeId(1);
	Node[1].SetCoord(0 ,  0.75*10);//coord X
	Node[1].SetCoord(1 , -0.75*10);//coord Y
	Node[1].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[1] = Node[1];
	
	
	///////DOIS////////
	Node[2].SetNodeId(2);
	Node[2].SetCoord(0 ,  0.75*10);//coord X
	Node[2].SetCoord(1 ,  0.75*10);//coord Y
	Node[2].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[2] = Node[2];
	
	///////TRES////////
	Node[3].SetNodeId(3);
	Node[3].SetCoord(0 ,  -0.75*10);//coord X
	Node[3].SetCoord(1 ,   0.75*10);//coord Y
	Node[3].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[3] = Node[3];
	
	///////QUATRO////////
	Node[4].SetNodeId(4);
	Node[4].SetCoord(0 ,  -5.3033*10);//coord X
	Node[4].SetCoord(1 ,  -5.3033*10);//coord Y
	Node[4].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[4] = Node[4];
	
	///////CINCO////////
	Node[5].SetNodeId(5);
	Node[5].SetCoord(0 ,  -2.87013*10);//coord X
	Node[5].SetCoord(1 , -6.92909*10);//coord Y
	Node[5].SetCoord(2,    0.);//coord Z
	gMesh->NodeVec()[5] = Node[5];
	
	///////SEIS////////
	Node[6].SetNodeId(6);
	Node[6].SetCoord(0 ,  -0.75*10);//coord X
	Node[6].SetCoord(1 ,  -7.46*10);//coord Y
	Node[6].SetCoord(2,    0.);//coord Z
	gMesh->NodeVec()[6] = Node[6];
	
	///////SETE////////
	Node[7].SetNodeId(7);
	Node[7].SetCoord(0 ,  0.*10);//coord X
	Node[7].SetCoord(1 ,   -7.5*10);//coord Y
	Node[7].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[7] = Node[7];
	
	///////OITO////////
	Node[8].SetNodeId(8);
	Node[8].SetCoord(0 ,  0.75*10);//coord X
	Node[8].SetCoord(1 ,  -7.46*10);//coord Y
	Node[8].SetCoord(2,   0. );//coord Z
	gMesh->NodeVec()[8] = Node[8];
	
	///////NOVE////////
	Node[9].SetNodeId(9);
	Node[9].SetCoord(0 ,  2.87013*10);//coord X
	Node[9].SetCoord(1 ,  -6.92909*10);//coord Y
	Node[9].SetCoord(2,   0. );//coord Z
	gMesh->NodeVec()[9] = Node[9];
	
	///////DEZ////////
	Node[10].SetNodeId(10);
	Node[10].SetCoord(0 ,  5.3033*10);//coord X
	Node[10].SetCoord(1 ,  -5.3033*10);//coord Y
	Node[10].SetCoord(2,   0.);//coord Z
	gMesh->NodeVec()[10] = Node[10];
	
	///////ONZE////////
	Node[11].SetNodeId(11);
	Node[11].SetCoord(0 ,  7.5*10);//coord X
	Node[11].SetCoord(1 ,  0.);//coord Y
	Node[11].SetCoord(2,   0. );//coord Z
	gMesh->NodeVec()[11] = Node[11];
	
	
	///////DOZE////////
	Node[12].SetNodeId(12);
	Node[12].SetCoord(0 ,  5.3033*10);//coord X
	Node[12].SetCoord(1 ,  5.3033*10);//coord Y
	Node[12].SetCoord(2,   0.);//coord Z
	gMesh->NodeVec()[12] = Node[12];
	
	///////TREZE////////
	Node[13].SetNodeId(13);
	Node[13].SetCoord(0 ,  2.87013*10);//coord X
	Node[13].SetCoord(1 ,  6.92909*10);//coord Y
	Node[13].SetCoord(2,   0. );//coord Z
	gMesh->NodeVec()[13] = Node[13];
	
	///////QUATORZE////////
	Node[14].SetNodeId(14);
	Node[14].SetCoord(0 ,  0.75*10);//coord X
	Node[14].SetCoord(1 ,  7.46*10);//coord Y
	Node[14].SetCoord(2,   0 );//coord Z
	gMesh->NodeVec()[14] = Node[14];
	
	///////QUINZE////////
	Node[15].SetNodeId(15);
	Node[15].SetCoord(0 ,  0.*10);//coord X
	Node[15].SetCoord(1 ,  7.5*10);//coord Y
	Node[15].SetCoord(2,   0. );//coord Z
	gMesh->NodeVec()[15] = Node[15];
	
	
	///////DEZESSEIS////////
	Node[16].SetNodeId(16);
	Node[16].SetCoord(0 ,  -0.75*10);//coord X
	Node[16].SetCoord(1 ,  7.46*10);//coord Y
	Node[16].SetCoord(2,   0. );//coord Z
	gMesh->NodeVec()[16] = Node[16];
	
	///////DEZESSETE////////
	Node[17].SetNodeId(17);
	Node[17].SetCoord(0 ,  -2.87013*10);//coord X
	Node[17].SetCoord(1 ,  6.92909*10);//coord Y
	Node[17].SetCoord(2,   0. );//coord Z
	gMesh->NodeVec()[17] = Node[17];
	
	
	///////DEZOITO////////
	Node[18].SetNodeId(18);
	Node[18].SetCoord(0 ,  -5.3033*10);//coord X
	Node[18].SetCoord(1 ,  5.3033*10);//coord Y
	Node[18].SetCoord(2,   0.);//coord Z
	gMesh->NodeVec()[18] = Node[18];
	
	///////DEZENOVE////////
	Node[19].SetNodeId(19);
	Node[19].SetCoord(0 ,  -7.5*10);//coord X
	Node[19].SetCoord(1 ,  0.*10);//coord Y
	Node[19].SetCoord(2,   0.);//coord Z
	gMesh->NodeVec()[19] = Node[19];
	
	///////VINTE////////
	Node[20].SetNodeId(20);
	Node[20].SetCoord(0 ,  -0.75*10);//coord X
	Node[20].SetCoord(1 ,  -8.46*10);//coord Y
	Node[20].SetCoord(2,   0.);//coord Z
	gMesh->NodeVec()[20] = Node[20];
	
	///////VINTE E UM////////
	Node[21].SetNodeId(21);
	Node[21].SetCoord(0 ,  0.75*10);//coord X
	Node[21].SetCoord(1 ,  -8.46*10);//coord Y
	Node[21].SetCoord(2,   0.);//coord Z
	gMesh->NodeVec()[21] = Node[21];
	
	///////VINTE E DOIS////////
	Node[22].SetNodeId(22);
	Node[22].SetCoord(0 ,  0.75*10);//coord X
	Node[22].SetCoord(1 ,  8.46*10);//coord Y
	Node[22].SetCoord(2,   0.);//coord Z
	gMesh->NodeVec()[22] = Node[22];
	
	///////VINTE E TRES////////
	Node[23].SetNodeId(23);
	Node[23].SetCoord(0 ,  -0.75*10);//coord X
	Node[23].SetCoord(1 ,  8.46*10);//coord Y
	Node[23].SetCoord(2,   0.);//coord Z
	gMesh->NodeVec()[23] = Node[23];
	
	//ELEMENTO 0///
	TopolQuad[0] = 0;	TopolQuad[1] = 1;
	TopolQuad[2] = 2; 	TopolQuad[3] = 3;  
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoQuad> > (0,TopolQuad,1,*gMesh);//id=0 // matind=1
	
	//ELEMENTO 1///
	TopolQuad[0] = 1;	TopolQuad[1] = 10;
	TopolQuad[2] = 12; 	TopolQuad[3] = 2; 
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoQuad> > (1,TopolQuad,1,*gMesh);//id=1 // matind=1
	
	//ELEMENTO 2///
	TopolTringle[0] = 2;	TopolTringle[1] = 12;
	TopolTringle[2] = 14; 	                    
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoTriangle> > (2,TopolTringle,1,*gMesh);//id=2 // matind=1
	
	//ELEMENTO 3///
	TopolQuad[0] = 3;	TopolQuad[1] = 2;
	TopolQuad[2] = 14; 	TopolQuad[3] = 16; 
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoQuad> > (3,TopolQuad,1,*gMesh);//id=3 // matind=1
	
	//ELEMENTO 4///
	TopolTringle[0] = 3;	TopolTringle[1] = 16;
	TopolTringle[2] = 18; 
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoTriangle> > (4,TopolTringle,1,*gMesh);//id=0 // matind=1
	
	//ELEMENTO 5///
	TopolQuad[0] = 0;	TopolQuad[1] = 3;
	TopolQuad[2] = 18; 	TopolQuad[3] = 4;   
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoQuad> > (5,TopolQuad,1,*gMesh);//id=0 // matind=1
	
	//ELEMENTO 6///
	TopolTringle[0] = 6;	TopolTringle[1] = 0;
	TopolTringle[2] = 4; 
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoTriangle> > (6,TopolTringle,1,*gMesh);//id=0 // matind=1
	
	//ELEMENTO 7///
	TopolQuad[0] = 6;	TopolQuad[1] = 8;
	TopolQuad[2] = 1; 	TopolQuad[3] = 0;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoQuad> > (7,TopolQuad,1,*gMesh);//id=0 // matind=1
	
	//ELEMENTO 8///
	TopolTringle[0] = 8;	TopolTringle[1] = 10;
	TopolTringle[2] = 1;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoTriangle> > (8,TopolTringle,1,*gMesh);//id=0 // matind=1
	
	// ARC 1 //
	TopolArc[0] = 18;
	TopolArc[1] = 4;
	TopolArc[2] = 19;
	new TPZGeoElRefPattern< TPZArc3D >  (9,TopolArc,-1,*gMesh);//id=1 // matind=-1
	
	// ARC 2 //
	TopolArc[0] = 4;//ESTE
	TopolArc[1] = 6;
	TopolArc[2] = 5;
	new TPZGeoElRefPattern< TPZArc3D > (10,TopolArc,-1,*gMesh);//id=2 // matind=-2
	
	// ARC 3 //
	TopolArc[0] = 6;
	TopolArc[1] = 8;
	TopolArc[2] = 7;
	new TPZGeoElRefPattern< TPZArc3D > (11,TopolArc,-1,*gMesh);//id=3 // matind=-3
	
	// ARC 4 //
	TopolArc[0] = 8;//ESTE
	TopolArc[1] = 10;
	TopolArc[2] = 9;
	new TPZGeoElRefPattern< TPZArc3D > (12,TopolArc,-1,*gMesh);//id=4 // matind=-4
	
	// ARC 5 //
	TopolArc[0] = 10;
	TopolArc[1] = 12;
	TopolArc[2] = 11;
	new TPZGeoElRefPattern< TPZArc3D >  (13,TopolArc,-1,*gMesh);//id=1 // matind=-1
	
	// ARC 6 //
	TopolArc[0] = 12;//ESTE
	TopolArc[1] = 14;
	TopolArc[2] = 13;
	new TPZGeoElRefPattern< TPZArc3D > (14,TopolArc,-1,*gMesh);//id=2 // matind=-2
	
	// ARC 7 //
	TopolArc[0] = 14;
	TopolArc[1] = 16;
	TopolArc[2] = 15;
	new TPZGeoElRefPattern< TPZArc3D > (15,TopolArc,-1,*gMesh);//id=3 // matind=-3
	
	// ARC 8 //
	TopolArc[0] = 16;//ESTE
	TopolArc[1] = 18;
	TopolArc[2] = 17;
	new TPZGeoElRefPattern< TPZArc3D > (16,TopolArc,-1,*gMesh);//id=4 // matind=-4
	
	//DownLine
	TopolLine[0]=6;
	TopolLine[1]=8;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (17, TopolLine,-2,*gMesh);
	
	//UperLine
	TopolLine[0]=14;
	TopolLine[1]=16;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (18, TopolLine,-3,*gMesh);
	
	
	
	//EDGE
	TopolLine[0]=4;
	TopolLine[1]=6;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (22, TopolLine,-7,*gMesh);
	
	//EDGE
	TopolLine[0]=8;
	TopolLine[1]=10;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (23, TopolLine,-8,*gMesh);
	//EDGE
	TopolLine[0]=12;
	TopolLine[1]=14;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (24, TopolLine,-9,*gMesh);
	
	//EDGE
	TopolLine[0]=16;
	TopolLine[1]=18;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (25, TopolLine,-10,*gMesh);
	
	
	
	//CENTER
	TopolLine[0]=0;
	TopolLine[1]=1;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (26, TopolLine,-11,*gMesh);
	
	//CENTER
	TopolLine[0]=1;
	TopolLine[1]=2;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (27, TopolLine,-12,*gMesh);
	////CENTER
	TopolLine[0]=2;
	TopolLine[1]=3;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (28, TopolLine,-13,*gMesh);
	
	////CENTER
	TopolLine[0]=3;
	TopolLine[1]=0;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (29, TopolLine,-14,*gMesh);
	
	
	//DIAGONALS
	TopolLine[0]=3;
	TopolLine[1]=18;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (30, TopolLine,-15,*gMesh);
	
	//DIAGONALS
	TopolLine[0]=2;
	TopolLine[1]=12;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (31, TopolLine,-16,*gMesh);
	////DIAGONALS
	TopolLine[0]=0;
	TopolLine[1]=4;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (32, TopolLine,-17,*gMesh);
	
	////DIAGONALS
	TopolLine[0]=1;
	TopolLine[1]=10;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (33, TopolLine,-18,*gMesh);
	
	//VERTICALS
	TopolLine[0]=3;
	TopolLine[1]=16;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (34, TopolLine,-19,*gMesh);
	
	//VERTICALS
	TopolLine[0]=2;
	TopolLine[1]=14;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (35, TopolLine,-20,*gMesh);
	////VERTICALS
	TopolLine[0]=1;
	TopolLine[1]=8;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (36, TopolLine,-21,*gMesh);
	
	////VERTICALS
	TopolLine[0]=0;
	TopolLine[1]=6;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (37, TopolLine,-22,*gMesh);
	
	
	
    //	//left down point
    //	TopolPoint[0]=6;
    //	new TPZGeoElRefPattern< TPZGeoPoint > (19, TopolPoint,-4,*gMesh);
    //	
    //	//right down point
    //	TopolPoint[0]=8;
    //	new TPZGeoElRefPattern<  TPZGeoPoint  > (20, TopolPoint,-5,*gMesh);
    //	
    //	//upper right point
    //	TopolPoint[0]=14;
    //	new TPZGeoElRefPattern<  TPZGeoPoint  > (21, TopolPoint,-6,*gMesh);
    
    //CENTER
	TopolPoint[0]=0;
	new TPZGeoElRefPattern< TPZGeoPoint > (38, TopolPoint,-4,*gMesh);
	
	//CENTER
	TopolPoint[0]=1;
	new TPZGeoElRefPattern<  TPZGeoPoint  > (39, TopolPoint,-5,*gMesh);
	
	//CENTER
	TopolPoint[0]=2;
	new TPZGeoElRefPattern<  TPZGeoPoint  > (40, TopolPoint,-6,*gMesh);
    
    //CENTER
	TopolPoint[0]=3;
	new TPZGeoElRefPattern<  TPZGeoPoint  > (41, TopolPoint,-7,*gMesh);
    
    //CENTER
	TopolPoint[0]=24;
	new TPZGeoElRefPattern<  TPZGeoPoint  > (42, TopolPoint,-8,*gMesh);
    
    //CENTER
	TopolPoint[0]=11;
	new TPZGeoElRefPattern<  TPZGeoPoint  > (43, TopolPoint,-9,*gMesh);
    
    //CENTER
	TopolPoint[0]=19;
	new TPZGeoElRefPattern<  TPZGeoPoint  > (44, TopolPoint,-10,*gMesh);
    
    //CENTER
	TopolPoint[0]=7;
	new TPZGeoElRefPattern<  TPZGeoPoint  > (45, TopolPoint,-11,*gMesh);
    
    //CENTER
	TopolPoint[0]=15;
	new TPZGeoElRefPattern<  TPZGeoPoint  > (46, TopolPoint,-12,*gMesh);
	
	
	
	gMesh->BuildConnectivity();
	
	for(int stage = 0; stage < h; stage++)
	{
		int nel = gMesh->NElements();
		for ( int iref = 0; iref < nel; iref++ )
		{
			TPZVec<TPZGeoEl*> tatara;
			TPZGeoEl * gelP1 = gMesh->ElementVec()[iref];
			if(!gelP1) continue;
			gelP1->Divide(tatara);
			
		}
	}
	
	std::set<int> matids;
	//EDGES
	matids.insert(-2);
	matids.insert(-3);
	//	matids.insert(-7);
	//	matids.insert(-8);
	//	matids.insert(-9);
	//	matids.insert(-10);
	
	//	//CENTER
	//	matids.insert(-11);
	//	matids.insert(-12);
	//	matids.insert(-13);
	//	matids.insert(-14);
	
	//DIAGONALS
	
	//	matids.insert(-15);
	//	matids.insert(-16);
	//	matids.insert(-17);
	//	matids.insert(-18);
	
	
	//VERTICALS
	//	matids.insert(-19);
	//	matids.insert(-20);
	//	matids.insert(-21);
	//	matids.insert(-22);
	
	for(int i = 0; i < dir; i++)
	{
		int nel = gMesh->NElements();
		for(int eel = 0; eel < nel; eel++)	
		{
			TPZGeoEl *elemento = gMesh->ElementVec()[eel];
			TPZRefPatternTools::RefineDirectional(elemento,matids);
            //			TPZRefPattern::RefineDirectional(elemento, matids);
		}
	}
	
	
	return gMesh;
	
	
}


TPZGeoMesh * BrazilianTestGeoMesh::TwoDMeshSlopeStability(int h, int dir)
{
	int Qnodes = 10;	
	
	
	TPZGeoMesh *gMesh = new TPZGeoMesh;
	
	gRefDBase.InitializeRefPatterns();
	//gMesh->ImportRefPattern ( );
	gMesh->NodeVec().Resize(Qnodes);
	
	TPZVec <int> TopolArc(3);
	TPZVec <int> TopolLine(2);
	TPZVec <int> TopolQuad(4);
	TPZVec <int> TopolPoint(1);
	TPZVec <int> TopolTringle(3);
	
	TPZVec<TPZGeoNode> Node(Qnodes);
	
	/*	
	 /////ZERO////////
	 Node[0].SetNodeId(0);
	 Node[0].SetCoord(0 , 0.);//coord X
	 Node[0].SetCoord(1 , 0.);//coord Y
	 Node[0].SetCoord(2,  0.);//coord Z
	 gMesh->NodeVec()[0] = Node[0];
	 
	 ///////UM////////
	 Node[1].SetNodeId(1);
	 Node[1].SetCoord(0 ,75.);//coord X
	 Node[1].SetCoord(1 ,0.);//coord Y
	 Node[1].SetCoord(2, 0.);//coord Z
	 gMesh->NodeVec()[1] = Node[1];
	 
	 
	 ///////DOIS////////
	 Node[2].SetNodeId(2);
	 Node[2].SetCoord(0 ,  75.);//coord X
	 Node[2].SetCoord(1 ,  30.);//coord Y
	 Node[2].SetCoord(2,    0. );//coord Z
	 gMesh->NodeVec()[2] = Node[2];
	 
	 ///////TRES////////
	 Node[3].SetNodeId(3);
	 Node[3].SetCoord(0 ,  45.);//coord X
	 Node[3].SetCoord(1 ,  30.);//coord Y
	 Node[3].SetCoord(2,    0. );//coord Z
	 gMesh->NodeVec()[3] = Node[3];
	 
	 ///////QUATRO////////
	 Node[4].SetNodeId(4);
	 Node[4].SetCoord(0 ,  25.);//coord X
	 Node[4].SetCoord(1 ,  30.);//coord Y
	 Node[4].SetCoord(2,    0. );//coord Z
	 gMesh->NodeVec()[4] = Node[4];
	 
	 ///////CINCO////////
	 Node[5].SetNodeId(5);
	 Node[5].SetCoord(0 ,  15.);//coord X
	 Node[5].SetCoord(1 ,  30.);//coord Y
	 Node[5].SetCoord(2,    0.);//coord Z
	 gMesh->NodeVec()[5] = Node[5];
	 
	 ///////SEIS////////
	 Node[6].SetNodeId(6);
	 Node[6].SetCoord(0 ,  0.);//coord X
	 Node[6].SetCoord(1 ,  30.);//coord Y
	 Node[6].SetCoord(2,    0.);//coord Z
	 gMesh->NodeVec()[6] = Node[6];
	 
	 ///////SETE////////
	 Node[7].SetNodeId(7);
	 Node[7].SetCoord(0 ,  35.);//coord X
	 Node[7].SetCoord(1 ,   40.);//coord Y
	 Node[7].SetCoord(2,    0. );//coord Z
	 gMesh->NodeVec()[7] = Node[7];
	 
	 ///////OITO////////
	 Node[8].SetNodeId(8);
	 Node[8].SetCoord(0 , 25.);//coord X
	 Node[8].SetCoord(1 , 40.);//coord Y
	 Node[8].SetCoord(2,   0. );//coord Z
	 gMesh->NodeVec()[8] = Node[8];
	 
	 ///////NOVE////////
	 Node[9].SetNodeId(9);
	 Node[9].SetCoord(0 ,  15.);//coord X
	 Node[9].SetCoord(1 ,  40.);//coord Y
	 Node[9].SetCoord(2,   0. );//coord Z
	 gMesh->NodeVec()[9] = Node[9];
	 
	 ///////DEZ////////
	 Node[10].SetNodeId(10);
	 Node[10].SetCoord(0 ,  0.);//coord X
	 Node[10].SetCoord(1 ,  40.);//coord Y
	 Node[10].SetCoord(2,   0.);//coord Z
	 gMesh->NodeVec()[10] = Node[10];
	 
	 
	 //ELEMENTO 0///
	 TopolQuad[0] = 0;	TopolQuad[1] = 1;
	 TopolQuad[2] = 2; 	TopolQuad[3] = 6;  
	 new TPZGeoElRefPattern< TPZGeoQuad > (0,TopolQuad,1,*gMesh);
	 
	 //ELEMENTO 1///
	 TopolQuad[0] = 4;	TopolQuad[1] = 3;
	 TopolQuad[2] = 7; 	TopolQuad[3] = 8; 
	 new TPZGeoElRefPattern< TPZGeoQuad > (1,TopolQuad,1,*gMesh);
	 
	 
	 //ELEMENTO 2///
	 TopolQuad[0] = 5;	TopolQuad[1] = 4;
	 TopolQuad[2] = 8; 	TopolQuad[3] = 9; 
	 new TPZGeoElRefPattern< TPZGeoQuad > (2,TopolQuad,1,*gMesh);
	 
	 
	 //ELEMENTO 3///
	 TopolQuad[0] = 6;	TopolQuad[1] = 5;
	 TopolQuad[2] = 9; 	TopolQuad[3] = 10;   
	 new TPZGeoElRefPattern<  TPZGeoQuad > (3,TopolQuad,1,*gMesh);
	 
	 
	 //DownLine
	 TopolLine[0]=0;
	 TopolLine[1]=1;
	 new TPZGeoElRefPattern< TPZGeoLinear  > (4, TopolLine,-1,*gMesh);
	 
	 //RightLine
	 TopolLine[0]=1;
	 TopolLine[1]=2;
	 new TPZGeoElRefPattern<  TPZGeoLinear  > (5, TopolLine,-2,*gMesh);
	 
	 //LeftLine
	 TopolLine[0]=6;
	 TopolLine[1]=0;
	 new TPZGeoElRefPattern<  TPZGeoLinear  > (6, TopolLine,-3,*gMesh);
	 
	 //LeftLine
	 TopolLine[0]=7;
	 TopolLine[1]=8;
	 new TPZGeoElRefPattern<  TPZGeoLinear  > (7, TopolLine,-4,*gMesh);
	 //LeftLine
	 TopolLine[0]=8;
	 TopolLine[1]=9;
	 new TPZGeoElRefPattern<  TPZGeoLinear  > (8, TopolLine,-5,*gMesh);
	 //LeftLine
	 TopolLine[0]=9;
	 TopolLine[1]=10;
	 new TPZGeoElRefPattern<  TPZGeoLinear  > (9, TopolLine,-6,*gMesh);
	 
	 */	
	/////ZERO////////
	Node[0].SetNodeId(0);
	Node[0].SetCoord(0 , 0.);//coord X
	Node[0].SetCoord(1 , 0.);//coord Y
	Node[0].SetCoord(2,  0.);//coord Z
	gMesh->NodeVec()[0] = Node[0];
	
	///////UM////////
	Node[1].SetNodeId(1);
	Node[1].SetCoord(0 ,60.);//coord X
	Node[1].SetCoord(1 ,0.);//coord Y
	Node[1].SetCoord(2, 0.);//coord Z
	gMesh->NodeVec()[1] = Node[1];
	
	
	///////DOIS////////
	Node[2].SetNodeId(2);
	Node[2].SetCoord(0 ,  60.);//coord X
	Node[2].SetCoord(1 ,  10.);//coord Y
	Node[2].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[2] = Node[2];
	
	///////TRES////////
	Node[3].SetNodeId(3);
	Node[3].SetCoord(0 ,  40.);//coord X
	Node[3].SetCoord(1 ,  10.);//coord Y
	Node[3].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[3] = Node[3];
	
	///////QUATRO////////
	Node[4].SetNodeId(4);
	Node[4].SetCoord(0 ,  15.);//coord X
	Node[4].SetCoord(1 ,  10.);//coord Y
	Node[4].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[4] = Node[4];
	
	///////CINCO////////
	Node[5].SetNodeId(5);
	Node[5].SetCoord(0 ,  20.);//coord X
	Node[5].SetCoord(1 ,  20.);//coord Y
	Node[5].SetCoord(2,    0.);//coord Z
	gMesh->NodeVec()[5] = Node[5];
	
	///////SEIS////////
	Node[6].SetNodeId(6);
	Node[6].SetCoord(0 ,  15.);//coord X
	Node[6].SetCoord(1 ,  20.);//coord Y
	Node[6].SetCoord(2,    0.);//coord Z
	gMesh->NodeVec()[6] = Node[6];
	
	///////SETE////////
	Node[7].SetNodeId(7);
	Node[7].SetCoord(0 ,    0.);//coord X
	Node[7].SetCoord(1 ,   20.);//coord Y
	Node[7].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[7] = Node[7];
	
	//	///////OITO////////
	//	Node[8].SetNodeId(8);
	//	Node[8].SetCoord(0 , 17.);//coord X
	//	Node[8].SetCoord(1 , 20.);//coord Y
	//	Node[8].SetCoord(2,   0. );//coord Z
	//	gMesh->NodeVec()[8] = Node[8];
	//	
	//	///////OITO////////
	//	Node[9].SetNodeId(9);
	//	Node[9].SetCoord(0 , 0.);//coord X
	//	Node[9].SetCoord(1 , 20.);//coord Y
	//	Node[9].SetCoord(2,   0. );//coord Z
	//	gMesh->NodeVec()[9] = Node[9];
	
	//ELEMENTO 0///
	TopolQuad[0] = 0;	TopolQuad[1] = 1;
	TopolQuad[2] = 3; 	TopolQuad[3] = 4;  
	new TPZGeoElRefPattern< TPZGeoQuad > (0,TopolQuad,1,*gMesh);
	
	//ELEMENTO 1///
	TopolTringle[0] = 1;	TopolTringle[1] = 2;
	TopolTringle[2] = 3;
	new TPZGeoElRefPattern< TPZGeoTriangle > (1,TopolTringle,1,*gMesh);
	
	
	//ELEMENTO 2///
	TopolQuad[0] = 3;	TopolQuad[1] = 5;
	TopolQuad[2] = 6; 	TopolQuad[3] = 4; 
	new TPZGeoElRefPattern< TPZGeoQuad > (2,TopolQuad,1,*gMesh);
	
	
	//ELEMENTO 3///
	TopolQuad[0] = 0;	TopolQuad[1] = 4;
	TopolQuad[2] = 6; 	TopolQuad[3] = 7;   
	new TPZGeoElRefPattern<  TPZGeoQuad > (3,TopolQuad,1,*gMesh);
	
	
	//DownLine
	TopolLine[0]=0;
	TopolLine[1]=1;
	new TPZGeoElRefPattern< TPZGeoLinear  > (4, TopolLine,-1,*gMesh);
	
	//RightLine
	TopolLine[0]=1;
	TopolLine[1]=2;
	new TPZGeoElRefPattern<  TPZGeoLinear  > (5, TopolLine,-2,*gMesh);
	
	//LeftLine
	TopolLine[0]=7;
	TopolLine[1]=0;
	new TPZGeoElRefPattern<  TPZGeoLinear  > (6, TopolLine,-3,*gMesh);
	
	
	
	gMesh->BuildConnectivity();
	
	for(int stage = 0; stage < h; stage++)
	{
		int nel = gMesh->NElements();
		for ( int iref = 0; iref < nel; iref++ )
		{
			TPZVec<TPZGeoEl*> tatara;
			TPZGeoEl * gelP1 = gMesh->ElementVec()[iref];
			if(!gelP1) continue;
			gelP1->Divide(tatara);
			
		}
	}
	
	//	std::set<int> matids;
	//	matids.insert(-2);
	//	matids.insert(-3);
	//
	//	
	//	for(int i = 0; i < dir; i++)
	//	{
	//		int nel = gMesh->NElements();
	//		for(int eel = 0; eel < nel; eel++)	
	//		{
	//			TPZGeoEl *elemento = gMesh->ElementVec()[eel];
	//			TPZRefPattern::RefineDirectional(elemento, matids);
	//		}
	//	}
	//	
	
	return gMesh;
	
}



TPZGeoMesh * BrazilianTestGeoMesh::TwoDMeshSlopeStability45(int h, int dir)
{
	int Qnodes = 16;	
	
	
	TPZGeoMesh *gMesh = new TPZGeoMesh;
	
	gRefDBase.InitializeRefPatterns();
	//	gMesh->ImportRefPattern ( );
	gMesh->NodeVec().Resize(Qnodes);
	
	TPZVec <int> TopolArc(3);
	TPZVec <int> TopolLine(2);
	TPZVec <int> TopolQuad(4);
	TPZVec <int> TopolPoint(1);
	TPZVec <int> TopolTringle(3);
	
	TPZVec<TPZGeoNode> Node(Qnodes);
	/////ZERO////////
	Node[0].SetNodeId(0);
	Node[0].SetCoord(0 , 0.);//coord X
	Node[0].SetCoord(1 , 0.);//coord Y
	Node[0].SetCoord(2,  0.);//coord Z
	gMesh->NodeVec()[0] = Node[0];
	
	///////UM////////
	Node[1].SetNodeId(1);
	Node[1].SetCoord(0 ,45.);//coord X
	Node[1].SetCoord(1 ,0.);//coord Y
	Node[1].SetCoord(2, 0.);//coord Z
	gMesh->NodeVec()[1] = Node[1];
	
	
	///////DOIS////////
	Node[2].SetNodeId(2);
	Node[2].SetCoord(0 ,  55.);//coord X
	Node[2].SetCoord(1 ,  0.);//coord Y
	Node[2].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[2] = Node[2];
	
	///////TRES////////
	Node[3].SetNodeId(3);
	Node[3].SetCoord(0 ,  75.);//coord X
	Node[3].SetCoord(1 ,  0.);//coord Y
	Node[3].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[3] = Node[3];
	
	///////QUATRO////////
	Node[4].SetNodeId(4);
	Node[4].SetCoord(0 ,  75.);//coord X
	Node[4].SetCoord(1 ,  25.);//coord Y
	Node[4].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[4] = Node[4];
	
	///////CINCO////////
	Node[5].SetNodeId(5);
	Node[5].SetCoord(0 ,  55.);//coord X
	Node[5].SetCoord(1 ,  25.);//coord Y
	Node[5].SetCoord(2,    0.);//coord Z
	gMesh->NodeVec()[5] = Node[5];
	
	///////SEIS////////
	Node[6].SetNodeId(6);
	Node[6].SetCoord(0 ,  45.);//coord X
	Node[6].SetCoord(1 ,  25.);//coord Y
	Node[6].SetCoord(2,    0.);//coord Z
	gMesh->NodeVec()[6] = Node[6];
	
	///////SETE////////
	Node[7].SetNodeId(7);
	Node[7].SetCoord(0 ,    20.);//coord X
	Node[7].SetCoord(1 ,    25.);//coord Y
	Node[7].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[7] = Node[7];
	
	//	///////OITO////////
	Node[8].SetNodeId(8);
	Node[8].SetCoord(0 , 75.);//coord X
	Node[8].SetCoord(1 , 30.);//coord Y
	Node[8].SetCoord(2,   0. );//coord Z
	gMesh->NodeVec()[8] = Node[8];
	//	
	//	///////NOVE////////
	Node[9].SetNodeId(9);
	Node[9].SetCoord(0 , 55.);//coord X
	Node[9].SetCoord(1 , 30.);//coord Y
	Node[9].SetCoord(2,   0. );//coord Z
	gMesh->NodeVec()[9] = Node[9];
	
	//	///////DEZ////////
	Node[10].SetNodeId(10);
	Node[10].SetCoord(0 , 45.);//coord X
	Node[10].SetCoord(1 , 30.);//coord Y
	Node[10].SetCoord(2,   0. );//coord Z
	gMesh->NodeVec()[10] = Node[10];
	
	//	///////ONZE////////
	Node[11].SetNodeId(11);
	Node[11].SetCoord(0 , 25.);//coord X
	Node[11].SetCoord(1 , 30.);//coord Y
	Node[11].SetCoord(2,   0. );//coord Z
	gMesh->NodeVec()[11] = Node[11];
	
	//	///////DOZE////////
	Node[12].SetNodeId(12);
	Node[12].SetCoord(0 , 35.);//coord X
	Node[12].SetCoord(1 , 40.);//coord Y
	Node[12].SetCoord(2,   0. );//coord Z
	gMesh->NodeVec()[12] = Node[12];
	
	//	///////TREZE////////
	Node[13].SetNodeId(13);
	Node[13].SetCoord(0 , 25.);//coord X
	Node[13].SetCoord(1 , 40.);//coord Y
	Node[13].SetCoord(2,   0. );//coord Z
	gMesh->NodeVec()[13] = Node[13];
	
	//	///////QUATORZE////////
	Node[14].SetNodeId(14);
	Node[14].SetCoord(0 , 20.);//coord X
	Node[14].SetCoord(1 , 40.);//coord Y
	Node[14].SetCoord(2,   0. );//coord Z
	gMesh->NodeVec()[14] = Node[14];
	
	//	///////QUINZE////////
	Node[15].SetNodeId(15);
	Node[15].SetCoord(0 , 0.);//coord X
	Node[15].SetCoord(1 , 40.);//coord Y
	Node[15].SetCoord(2,   0. );//coord Z
	gMesh->NodeVec()[15] = Node[15];
	
	
	TopolQuad[0] = 0;	TopolQuad[1] = 7;
	TopolQuad[2] = 14; 	TopolQuad[3] = 15;  
	new TPZGeoElRefPattern< TPZGeoQuad > (0,TopolQuad,1,*gMesh);
	
	TopolQuad[0] = 0;	TopolQuad[1] = 1;
	TopolQuad[2] = 6; 	TopolQuad[3] = 7;  
	new TPZGeoElRefPattern< TPZGeoQuad > (1,TopolQuad,1,*gMesh);
	
	TopolQuad[0] = 1;	TopolQuad[1] = 2;
	TopolQuad[2] = 5; 	TopolQuad[3] = 6; 
	new TPZGeoElRefPattern< TPZGeoQuad > (2,TopolQuad,1,*gMesh);
	
	TopolQuad[0] = 2;	TopolQuad[1] = 3;
	TopolQuad[2] = 4; 	TopolQuad[3] = 5;   
	new TPZGeoElRefPattern<  TPZGeoQuad > (3,TopolQuad,1,*gMesh);
	
	TopolQuad[0] = 5;	TopolQuad[1] = 4;
	TopolQuad[2] = 8; 	TopolQuad[3] = 9;   
	new TPZGeoElRefPattern<  TPZGeoQuad > (4,TopolQuad,1,*gMesh);
	
	TopolQuad[0] = 6;	TopolQuad[1] = 5;
	TopolQuad[2] = 9; 	TopolQuad[3] = 10;   
	new TPZGeoElRefPattern<  TPZGeoQuad > (5,TopolQuad,1,*gMesh);
	
	TopolQuad[0] = 7;	TopolQuad[1] = 6;
	TopolQuad[2] = 10; 	TopolQuad[3] = 11;   
	new TPZGeoElRefPattern<  TPZGeoQuad > (6,TopolQuad,1,*gMesh);
	
	TopolQuad[0] = 7;	TopolQuad[1] = 11;
	TopolQuad[2] = 13; 	TopolQuad[3] = 14;   
	new TPZGeoElRefPattern<  TPZGeoQuad > (7,TopolQuad,1,*gMesh);
	
	TopolQuad[0] = 11;	TopolQuad[1] = 10;
	TopolQuad[2] = 12; 	TopolQuad[3] = 13;   
	new TPZGeoElRefPattern<  TPZGeoQuad > (8,TopolQuad,1,*gMesh);
	
	//DownLine1
	TopolLine[0]=0;
	TopolLine[1]=1;
	new TPZGeoElRefPattern< TPZGeoLinear  > (9, TopolLine,-1,*gMesh);
	
	//DownLine2
	TopolLine[0]=1;
	TopolLine[1]=2;
	new TPZGeoElRefPattern<  TPZGeoLinear  > (10, TopolLine,-1,*gMesh);
	
	//DownLine3
	TopolLine[0]=2;
	TopolLine[1]=3;
	new TPZGeoElRefPattern<  TPZGeoLinear  > (11, TopolLine,-1,*gMesh);
	
	//LeftLine
	TopolLine[0]=15;
	TopolLine[1]=0;
	new TPZGeoElRefPattern<  TPZGeoLinear  > (12, TopolLine,-2,*gMesh);
	
	//RightLine
	TopolLine[0]=3;
	TopolLine[1]=4;
	new TPZGeoElRefPattern<  TPZGeoLinear  > (13, TopolLine,-3,*gMesh);
	
	//RightLine
	TopolLine[0]=4;
	TopolLine[1]=8;
	new TPZGeoElRefPattern<  TPZGeoLinear  > (14, TopolLine,-3,*gMesh);
	
	
	gMesh->BuildConnectivity();
	
	TPZVec<TPZGeoEl *> filhos;
	int n = gMesh->NElements();
	for(int i = 0; i < n; i++)
	{
		TPZGeoEl * gel = gMesh->ElementVec()[i];
		if(gel->HasSubElement()) continue;
		if (gel->Id() == 5 || gel->Id() == 6 || gel->Id() == 7 || gel->Id() == 8) 
		{
			gel->Divide(filhos);
		}
	}
	
	for(int stage = 0; stage < h; stage++)
	{
		int nel = gMesh->NElements();
		for ( int iref = 0; iref < nel; iref++ )
		{
			TPZVec<TPZGeoEl*> tatara;
			TPZGeoEl * gelP1 = gMesh->ElementVec()[iref];
			if(!gelP1) continue;
			gelP1->Divide(tatara);
			
		}
	}
	
	return gMesh;
	
}

TPZGeoMesh * BrazilianTestGeoMesh::TwoDMeshSlopeStability452(int h, int dir)
{
	int Qnodes = 9;	
	
	
	TPZGeoMesh *gMesh = new TPZGeoMesh;
	
	gRefDBase.InitializeRefPatterns();
	//gMesh->ImportRefPattern ( );
	gMesh->NodeVec().Resize(Qnodes);
	
	TPZVec <int> TopolArc(3);
	TPZVec <int> TopolLine(2);
	TPZVec <int> TopolQuad(4);
	TPZVec <int> TopolPoint(1);
	TPZVec <int> TopolTringle(3);
	
	TPZVec<TPZGeoNode> Node(Qnodes);
	/////ZERO////////
	Node[0].SetNodeId(0);
	Node[0].SetCoord(0 , 0.);//coord X
	Node[0].SetCoord(1 , 0.);//coord Y
	Node[0].SetCoord(2,  0.);//coord Z
	gMesh->NodeVec()[0] = Node[0];
	
	///////UM////////
	Node[1].SetNodeId(1);
	Node[1].SetCoord(0 ,45.);//coord X
	Node[1].SetCoord(1 ,0.);//coord Y
	Node[1].SetCoord(2, 0.);//coord Z
	gMesh->NodeVec()[1] = Node[1];
	
	
	///////DOIS////////
	Node[2].SetNodeId(2);
	Node[2].SetCoord(0 ,  75.);//coord X
	Node[2].SetCoord(1 ,  0.);//coord Y
	Node[2].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[2] = Node[2];
	
	///////TRES////////
	Node[3].SetNodeId(3);
	Node[3].SetCoord(0 ,  75.);//coord X
	Node[3].SetCoord(1 ,  30.);//coord Y
	Node[3].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[3] = Node[3];
	
	///////QUATRO////////
	Node[4].SetNodeId(4);
	Node[4].SetCoord(0 ,  45.);//coord X
	Node[4].SetCoord(1 ,  30.);//coord Y
	Node[4].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[4] = Node[4];
	
	///////CINCO////////
	Node[5].SetNodeId(5);
	Node[5].SetCoord(0 ,  28.);//coord X
	Node[5].SetCoord(1 ,  30.);//coord Y
	Node[5].SetCoord(2,    0.);//coord Z
	gMesh->NodeVec()[5] = Node[5];
	
	///////SEIS////////
	Node[6].SetNodeId(6);
	Node[6].SetCoord(0 ,  35.);//coord X
	Node[6].SetCoord(1 ,  40.);//coord Y
	Node[6].SetCoord(2,    0.);//coord Z
	gMesh->NodeVec()[6] = Node[6];
	
	///////SETE////////
	Node[7].SetNodeId(7);
	Node[7].SetCoord(0 ,    28.);//coord X
	Node[7].SetCoord(1 ,   40.);//coord Y
	Node[7].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[7] = Node[7];
	
	///////SETE////////
	Node[8].SetNodeId(8);
	Node[8].SetCoord(0 ,  0.);//coord X
	Node[8].SetCoord(1 ,  40.);//coord Y
	Node[8].SetCoord(2,   0. );//coord Z
	gMesh->NodeVec()[8] = Node[8];
	
	//ELEMENTO 0///
	TopolQuad[0] = 0;	TopolQuad[1] = 1;
	TopolQuad[2] = 4; 	TopolQuad[3] = 5;  
	new TPZGeoElRefPattern< TPZGeoQuad > (0,TopolQuad,1,*gMesh);
	
	//ELEMENTO 0///
	TopolQuad[0] = 1;	TopolQuad[1] = 2;
	TopolQuad[2] = 3; 	TopolQuad[3] = 4;  
	new TPZGeoElRefPattern< TPZGeoQuad > (1,TopolQuad,1,*gMesh);
	
	
	//ELEMENTO 2///
	TopolQuad[0] = 4;	TopolQuad[1] = 6;
	TopolQuad[2] = 7; 	TopolQuad[3] = 5; 
	new TPZGeoElRefPattern< TPZGeoQuad > (2,TopolQuad,1,*gMesh);
	
	
	//ELEMENTO 3///
	TopolQuad[0] = 5;	TopolQuad[1] = 7;
	TopolQuad[2] = 8; 	TopolQuad[3] = 0;   
	new TPZGeoElRefPattern<  TPZGeoQuad > (3,TopolQuad,1,*gMesh);
	
	//LeftLine
	TopolLine[0]=0;
	TopolLine[1]=1;
	new TPZGeoElRefPattern<  TPZGeoLinear  > (4, TopolLine,-1,*gMesh);
	//LeftLine
	TopolLine[0]=1;
	TopolLine[1]=2;
	new TPZGeoElRefPattern<  TPZGeoLinear  > (5, TopolLine,-2,*gMesh);
	//LeftLine
	TopolLine[0]=2;
	TopolLine[1]=3;
	new TPZGeoElRefPattern<  TPZGeoLinear  > (6, TopolLine,-3,*gMesh);
	//LeftLine
	TopolLine[0]=8;
	TopolLine[1]=0;
	new TPZGeoElRefPattern<  TPZGeoLinear  > (7, TopolLine,-4,*gMesh);
	
	TopolLine[0]=3;
	TopolLine[1]=4;
	new TPZGeoElRefPattern<  TPZGeoLinear  > (8, TopolLine,-5,*gMesh);
	
	TopolPoint[0] = 4;
	new TPZGeoElRefPattern<  TPZGeoPoint  > (9, TopolPoint,11,*gMesh);
	
	gMesh->BuildConnectivity();
	
	TPZVec<TPZGeoEl *> filhos;
	int n = gMesh->NElements();
	for(int i = 0; i < n; i++)
	{
		TPZGeoEl * gel = gMesh->ElementVec()[i];
		if(gel->HasSubElement()) continue;
		if (gel->Id() == 2) 
		{
			gel->Divide(filhos);
		}
	}
	/*
	 for(int ref = 0; ref < h; ref++)
	 {
	 TPZVec<TPZGeoEl *> tatara;
	 int n = gMesh->NElements();
	 for(int i = 0; i < n; i++)
	 {
	 TPZGeoEl * gel = gMesh->ElementVec()[i];
	 if(gel->HasSubElement()) continue;
	 gel->Divide(tatara);
	 }
	 }
	 
	 */
	for(int cc = 0; cc < h; cc++)
	{
		int nElel = gMesh->NElements();
		TPZVec<TPZGeoEl *> tatara2;
		for(int el = 0; el < nElel; el++)
		{
			TPZGeoEl * gel = gMesh->ElementVec()[el];
			TPZGeoElSide  sidegel;
			
			//			virtual int SideDimension(int side) = 0;
			//			 virtual TPZGeoElSide HigherDimensionSides(int side,int targetdimension);//S�PARA TESTAR CONTINUIDADE - APAGAR DEPOIS
			//			int WhichSubel();
			//			virtual int NSideSubElements2(int side) = 0;
			int sides=gel->NSides(); 
			for(int i=0;i<sides;i++)
			{
				//int nsubelemens = gel->NSideSubElements2(i);
				int nsubelements = gel->NSubElements();
				sidegel = gel->Neighbour(i);
				TPZGeoEl *vizinho = sidegel.Element();
				int sidesubelements = vizinho->NSubElements();
				//	int sidesubelements = sidegel.NSubElements2();
				if(nsubelements<sidesubelements)
				{
					gel->Divide(tatara2);
				}
			}
		}
	}
	
	//	std::set<int> matids;
	//	matids.insert(-5);
	//
	//	for(int i = 0; i < dir; i++)
	//	{
	//		int nel = gMesh->NElements();
	//		for(int eel = 0; eel < nel; eel++)	
	//		{
	//			TPZGeoEl *elemento = gMesh->ElementVec()[eel];
	//			TPZRefPattern::RefineDirectional(elemento, matids);
	//		}
	//	}
	
	
	return gMesh;
	
}

TPZGeoMesh * BrazilianTestGeoMesh::MisesPressure(int h, int dir)
{
	int Qnodes = 6;	
	
	
	TPZGeoMesh *gMesh = new TPZGeoMesh;
	
	gRefDBase.InitializeRefPatterns();
	//	gMesh->ImportRefPattern ( );
	gMesh->NodeVec().Resize(Qnodes);
	
	TPZVec <int> TopolArc(3);
	TPZVec <int> TopolLine(2);
	TPZVec <int> TopolQuad(4);
	TPZVec <int> TopolPoint(1);
	
	TPZVec<TPZGeoNode> Node(Qnodes);
	/////ZERO////////
	Node[0].SetNodeId(0);
	Node[0].SetCoord(0 , 100.);//coord X
	Node[0].SetCoord(1 , 0.);//coord Y
	Node[0].SetCoord(2,  0.);//coord Z
	gMesh->NodeVec()[0] = Node[0];
	
	///////UM////////
	Node[1].SetNodeId(1);
	Node[1].SetCoord(0 ,200.);//coord X
	Node[1].SetCoord(1 ,0.);//coord Y
	Node[1].SetCoord(2, 0.);//coord Z
	gMesh->NodeVec()[1] = Node[1];
	
	
	///////DOIS////////
	Node[2].SetNodeId(2);
	Node[2].SetCoord(0 ,  100*sqrt(2.));//coord X
	Node[2].SetCoord(1 ,  100*sqrt(2.));//coord Y
	Node[2].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[2] = Node[2];
	
	///////TRES////////
	Node[3].SetNodeId(3);
	Node[3].SetCoord(0 ,  0.);//coord X
	Node[3].SetCoord(1 ,  200.);//coord Y
	Node[3].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[3] = Node[3];
	
	
	///////TRES////////
	Node[4].SetNodeId(4);
	Node[4].SetCoord(0 ,  0.);//coord X
	Node[4].SetCoord(1 ,  100.);//coord Y
	Node[4].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[4] = Node[4];
	
	///////TRES////////
	Node[5].SetNodeId(5);
	Node[5].SetCoord(0 ,  50*sqrt(2.));//coord X
	Node[5].SetCoord(1 ,  50*sqrt(2.));//coord Y
	Node[5].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[5] = Node[5];
	
	//ELEMENTO 0///
	TopolQuad[0] = 0;	TopolQuad[1] = 1;
	TopolQuad[2] = 3; 	TopolQuad[3] = 4;  
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoQuad > > (0,TopolQuad,1,*gMesh);
	
	
	//ELEMENTO 2///
	TopolLine[0] = 0;	TopolLine[1] = 1;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (1,TopolLine,-1,*gMesh);
	
	//ELEMENTO 2///
	TopolLine[0] = 3;	TopolLine[1] = 4;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (2,TopolLine,-2,*gMesh);
	
	//ELEMENTO 2///
	TopolLine[0] = 4;	TopolLine[1] = 0;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (3,TopolLine,-3,*gMesh);
	
	//ELEMENTO 2///
	TopolArc[0] = 4;	TopolArc[1] =0;		TopolArc[2] =5;
	new TPZGeoElRefPattern<  TPZArc3D  > (4,TopolArc,-4,*gMesh);
	
	//ELEMENTO 2///
	TopolArc[0] = 1;	TopolArc[1] =3;		TopolArc[2] =2;
	new TPZGeoElRefPattern< TPZArc3D  > (5,TopolArc,-5,*gMesh);
	/*	
	 //ELEMENTO 2///
	 TopolPoint[0] = 0;
	 new TPZGeoElRefPattern< TPZGeoPoint  > (6,TopolPoint,-6,*gMesh);
	 
	 TopolPoint[0] = 1;
	 new TPZGeoElRefPattern< TPZGeoPoint  > (7,TopolPoint,-7,*gMesh);
	 
	 //ELEMENTO 2///
	 TopolPoint[0] = 3;
	 new TPZGeoElRefPattern< TPZGeoPoint  > (8,TopolPoint,-8,*gMesh);
	 
	 TopolPoint[0] = 4;
	 new TPZGeoElRefPattern< TPZGeoPoint  > (9,TopolPoint,-9,*gMesh);
	 */	
	gMesh->BuildConnectivity();
	
	for(int ref = 0; ref < h; ref++)
	{
		TPZVec<TPZGeoEl *> tatara;
		int n = gMesh->NElements();
		for(int i = 0; i < n; i++)
		{
			TPZGeoEl * gel = gMesh->ElementVec()[i];
			//if(gel->HasSubElement()) continue;
			gel->Divide(tatara);
		}
	}
	/*	
	 for(int cc = 0; cc < dir; cc++)
	 {
	 int nElel = gMesh->NElements();
	 std::set<int> matToRefine;
	 matToRefine.insert(-3);
	 for(int el = 0; el < nElel; el++)
	 {
	 TPZGeoEl * gel = gMesh->ElementVec()[el];
	 TPZRefPattern::RefineDirectional(gel, matToRefine);
	 }
	 }
	 */
	//	std::set<int> matids;
	//	matids.insert(-5);
	//
	//	for(int i = 0; i < dir; i++)
	//	{
	//		int nel = gMesh->NElements();
	//		for(int eel = 0; eel < nel; eel++)	
	//		{
	//			TPZGeoEl *elemento = gMesh->ElementVec()[eel];
	//			TPZRefPattern::RefineDirectional(elemento, matids);
	//		}
	//	}
	
	
	return gMesh;
	
}

TPZGeoMesh * BrazilianTestGeoMesh::MisesPressure2(int h, int dir)
{
	int Qnodes = 16;	
	
	
	TPZGeoMesh *gMesh = new TPZGeoMesh;
	
	gRefDBase.InitializeRefPatterns();
	//	gMesh->ImportRefPattern ( );
	gMesh->NodeVec().Resize(Qnodes);
	
	TPZVec <int> TopolArc(3);
	TPZVec <int> TopolLine(2);
	TPZVec <int> TopolQuad(4);
	TPZVec <int> TopolPoint(1);
	
	TPZVec<TPZGeoNode> Node(Qnodes);
	/////ZERO////////
	Node[0].SetNodeId(0);
	Node[0].SetCoord(0 , 100.);//coord X
	Node[0].SetCoord(1 , 0.);//coord Y
	Node[0].SetCoord(2,  0.);//coord Z
	gMesh->NodeVec()[0] = Node[0];
	
	///////UM////////
	Node[1].SetNodeId(1);
	Node[1].SetCoord(0 ,50*sqrt(2.));//coord X
	Node[1].SetCoord(1 ,50*sqrt(2.));//coord Y
	Node[1].SetCoord(2, 0.);//coord Z
	gMesh->NodeVec()[1] = Node[1];
	
	
	///////DOIS////////
	Node[2].SetNodeId(2);
	Node[2].SetCoord(0 ,  0.);//coord X
	Node[2].SetCoord(1 ,  100.);//coord Y
	Node[2].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[2] = Node[2];
	
	///////TRES////////
	Node[3].SetNodeId(3);
	Node[3].SetCoord(0 ,  -50*sqrt(2.));//coord X
	Node[3].SetCoord(1 ,  50*sqrt(2.));//coord Y
	Node[3].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[3] = Node[3];
	
	
	///////TRES////////
	Node[4].SetNodeId(4);
	Node[4].SetCoord(0 ,  -100.);//coord X
	Node[4].SetCoord(1 ,  0.);//coord Y
	Node[4].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[4] = Node[4];
	
	///////TRES////////
	Node[5].SetNodeId(5);
	Node[5].SetCoord(0 , -50*sqrt(2.));//coord X
	Node[5].SetCoord(1 ,  -50*sqrt(2.));//coord Y
	Node[5].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[5] = Node[5];
	
	
	///////TRES////////
	Node[6].SetNodeId(6);
	Node[6].SetCoord(0 , 0.);//coord X
	Node[6].SetCoord(1 ,  -100.);//coord Y
	Node[6].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[6] = Node[6];
	
	
	
	///////TRES////////
	Node[7].SetNodeId(7);
	Node[7].SetCoord(0 , 50*sqrt(2.));//coord X
	Node[7].SetCoord(1 ,  -50*sqrt(2.));//coord Y
	Node[7].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[7] = Node[7];
	
	///////TRES////////
	Node[8].SetNodeId(8);
	Node[8].SetCoord(0 , 200.);//coord X
	Node[8].SetCoord(1 ,  0.);//coord Y
	Node[8].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[8] = Node[8];
	
	
	Node[9].SetNodeId(9);
	Node[9].SetCoord(0 , 100*sqrt(2.));//coord X
	Node[9].SetCoord(1 , 100*sqrt(2.));//coord Y
	Node[9].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[9] = Node[9];
	
	Node[10].SetNodeId(10);
	Node[10].SetCoord(0 , 0.);//coord X
	Node[10].SetCoord(1 , 200.);//coord Y
	Node[10].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[10] = Node[10];
	
	Node[11].SetNodeId(11);
	Node[11].SetCoord(0 , -100*sqrt(2.));//coord X
	Node[11].SetCoord(1 , 100*sqrt(2.));//coord Y
	Node[11].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[11] = Node[11];
	
	Node[12].SetNodeId(12);
	Node[12].SetCoord(0 , -200.);//coord X
	Node[12].SetCoord(1 ,    0.);//coord Y
	Node[12].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[12] = Node[12];
	
	Node[13].SetNodeId(13);
	Node[13].SetCoord(0 , -100*sqrt(2.));//coord X
	Node[13].SetCoord(1 ,    -100*sqrt(2.));//coord Y
	Node[13].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[13] = Node[13];
	
	Node[14].SetNodeId(14);
	Node[14].SetCoord(0 , 0.);//coord X
	Node[14].SetCoord(1 ,   -200.);//coord Y
	Node[14].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[14] = Node[14];
	
	Node[15].SetNodeId(14);
	Node[15].SetCoord(0 , 100*sqrt(2.));//coord X
	Node[15].SetCoord(1 ,   -100*sqrt(2.));//coord Y
	Node[15].SetCoord(2,    0. );//coord Z
	gMesh->NodeVec()[15] = Node[15];
	
	//ELEMENTO 0///
	TopolQuad[0] = 0;	TopolQuad[1] = 8;
	TopolQuad[2] = 10; 	TopolQuad[3] = 2;  
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoQuad > > (0,TopolQuad,1,*gMesh);
	
	//ELEMENTO 0///
	TopolQuad[0] = 2;	TopolQuad[1] = 10;
	TopolQuad[2] = 12; 	TopolQuad[3] = 4;  
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoQuad > > (1,TopolQuad,1,*gMesh);
	
	//ELEMENTO 0///
	TopolQuad[0] = 4;	TopolQuad[1] = 12;
	TopolQuad[2] = 14; 	TopolQuad[3] = 6;  
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoQuad > > (2,TopolQuad,1,*gMesh);
	
	//ELEMENTO 0///
	TopolQuad[0] = 6;	TopolQuad[1] = 14;
	TopolQuad[2] = 8; 	TopolQuad[3] = 0;  
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoQuad > > (3,TopolQuad,1,*gMesh);
	
	
	//ELEMENTO 2///
	TopolLine[0] = 0;	TopolLine[1] = 2;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (4,TopolLine,-1,*gMesh);
	
	//ELEMENTO 2///
	TopolLine[0] = 2;	TopolLine[1] = 4;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (5,TopolLine,-1,*gMesh);
	
	//ELEMENTO 2///
	TopolLine[0] = 4;	TopolLine[1] = 6;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (6,TopolLine,-1,*gMesh);
	
	
	//ELEMENTO 2///
	TopolLine[0] = 6;	TopolLine[1] = 0;
	new TPZGeoElRefPattern< TPZGeoBlend< TPZGeoLinear > > (7,TopolLine,-1,*gMesh);
	
	//ELEMENTO 2///
	TopolArc[0] = 0;	TopolArc[1] =2;		TopolArc[2] =1;
	new TPZGeoElRefPattern<  TPZArc3D  > (8,TopolArc,-5,*gMesh);
	
	//ELEMENTO 2///
	TopolArc[0] = 2;	TopolArc[1] =4;		TopolArc[2] =3;
	new TPZGeoElRefPattern< TPZArc3D  > (9,TopolArc,-5,*gMesh);
	
	//ELEMENTO 2///
	TopolArc[0] = 4;	TopolArc[1] =6;		TopolArc[2] =5;
	new TPZGeoElRefPattern< TPZArc3D  > (10,TopolArc,-5,*gMesh);
	
	//ELEMENTO 2///
	TopolArc[0] = 6;	TopolArc[1] =0;		TopolArc[2] =7;
	new TPZGeoElRefPattern< TPZArc3D  > (11,TopolArc,-5,*gMesh);
	
	//ELEMENTO 2///
	TopolArc[0] = 8;	TopolArc[1] =10;		TopolArc[2] =9;
	new TPZGeoElRefPattern< TPZArc3D  > (12,TopolArc,-6,*gMesh);
	//ELEMENTO 2///
	TopolArc[0] = 10;	TopolArc[1] =12;		TopolArc[2] =11;
	new TPZGeoElRefPattern< TPZArc3D  > (13,TopolArc,-6,*gMesh);
	//ELEMENTO 2///
	TopolArc[0] = 12;	TopolArc[1] =14;		TopolArc[2] =13;
	new TPZGeoElRefPattern< TPZArc3D  > (14,TopolArc,-6,*gMesh);
	//ELEMENTO 2///
	TopolArc[0] = 14;	TopolArc[1] =8;		TopolArc[2] =15;
	new TPZGeoElRefPattern< TPZArc3D  > (15,TopolArc,-6,*gMesh);
	
	//ELEMENTO 2///
	TopolPoint[0] = 8;
	new TPZGeoElRefPattern< TPZGeoPoint  > (16,TopolPoint,-7,*gMesh);
	
	TopolPoint[0] = 12;
	new TPZGeoElRefPattern< TPZGeoPoint  > (17,TopolPoint,-8,*gMesh);
	
	//ELEMENTO 2///
	TopolPoint[0] = 10;
	new TPZGeoElRefPattern< TPZGeoPoint  > (18,TopolPoint,-9,*gMesh);
	
	TopolPoint[0] = 14;
	new TPZGeoElRefPattern< TPZGeoPoint  > (19,TopolPoint,-10,*gMesh);
	
	
	gMesh->BuildConnectivity();
	
	for(int ref = 0; ref < h; ref++)
	{
		TPZVec<TPZGeoEl *> tatara;
		int n = gMesh->NElements();
		for(int i = 0; i < n; i++)
		{
			TPZGeoEl * gel = gMesh->ElementVec()[i];
			//if(gel->HasSubElement()) continue;
			gel->Divide(tatara);
		}
	}
	/*	
	 for(int cc = 0; cc < dir; cc++)
	 {
	 int nElel = gMesh->NElements();
	 std::set<int> matToRefine;
	 matToRefine.insert(-3);
	 for(int el = 0; el < nElel; el++)
	 {
	 TPZGeoEl * gel = gMesh->ElementVec()[el];
	 TPZRefPattern::RefineDirectional(gel, matToRefine);
	 }
	 }
	 */
	//	std::set<int> matids;
	//	matids.insert(-5);
	//
	//	for(int i = 0; i < dir; i++)
	//	{
	//		int nel = gMesh->NElements();
	//		for(int eel = 0; eel < nel; eel++)	
	//		{
	//			TPZGeoEl *elemento = gMesh->ElementVec()[eel];
	//			TPZRefPattern::RefineDirectional(elemento, matids);
	//		}
	//	}
	
	
	return gMesh;
	
}


//TPZGeoMesh *BrazilianTestGeoMesh::MalhaPredio()
//{
//
//	int numnodes;//=-1;
//	int numelements;//=-1;
//	
//	string FileName;
//	FileName = "gidmeshnodes.txt";
//    ifstream read (FileName.c_str());
//
//    	
//    int nodeId = 0, elementId = 0, matElId = 1;
//    
//    double nodecoordX , nodecoordY , nodecoordZ ;
//    read >> numnodes;
//    
//    TPZGeoMesh * gMesh = new TPZGeoMesh;
//    
//    gMesh -> NodeVec().Resize(numnodes);
//    
//    const int Qnodes = numnodes;
//    TPZVec <TPZGeoNode> Node(Qnodes);
//    
//    for(int in=0; in<numnodes; in++)
//    {
//        read >> nodeId;
//        read >> nodecoordX;
//        read >> nodecoordY;
//        read >> nodecoordZ;
//        Node[nodeId-1].SetNodeId(nodeId);
//        Node[nodeId-1].SetCoord(0,nodecoordX);
//        Node[nodeId-1].SetCoord(1,nodecoordY);
//        Node[nodeId-1].SetCoord(2,nodecoordZ);
//        gMesh->NodeVec()[nodeId-1] = Node[nodeId-1];
//            
//    }
//
//    string FileName2;
//	FileName2 = "gidmeshelements.txt";
//    ifstream read2 (FileName2.c_str());
//    read2 >> numelements;
//    TPZVec <int> TopolTriangle(3);
//	
//    for(int el=0; el<numelements; el++)
//    {
//        int topol1,topol2,topol3;
//        read2 >> elementId;
//        read2 >> topol1; //node 1
//        read2 >> topol2; //node 2
//        read2 >> topol3; //node 3
//        
//        TopolTriangle[0]=topol1;
//        TopolTriangle[1]=topol2;
//        TopolTriangle[2]=topol3;
//        
//        TopolTriangle[0]--;
//        TopolTriangle[1]--;
//        TopolTriangle[2]--;
//        
//        int index = el;
//        new TPZGeoElRefPattern< TPZGeoTriangle  > (elementId,TopolTriangle,1,*gMesh);
//       // TPZGeoEl * tri = new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (index, TopolTriangle, matElId, *gMesh);
//        
//    }
//  //  gMesh->BuildConnectivity();
//
//            
//    	
// 
// //   TPZVec <int> TopolPoint(1);
//    
////    
////	TopolPoint[0] = 66;
////	new TPZGeoElRefPattern< TPZGeoPoint  > (520,TopolPoint,-7,*gMesh);
////	
////	TopolPoint[0] = 56;
////	new TPZGeoElRefPattern< TPZGeoPoint  > (521,TopolPoint,-8,*gMesh);
////	
////	//ELEMENTO 2///
////	TopolPoint[0] = 50;
////	new TPZGeoElRefPattern< TPZGeoPoint  > (522,TopolPoint,-9,*gMesh);
////	
////	TopolPoint[0] = 40;
////	new TPZGeoElRefPattern< TPZGeoPoint  > (523,TopolPoint,-10,*gMesh);
////    
////    
////    
////    
////    TopolPoint[0] = 292;
////	new TPZGeoElRefPattern< TPZGeoPoint  > (524,TopolPoint,-11,*gMesh);
////	
////	TopolPoint[0] = 286;
////	new TPZGeoElRefPattern< TPZGeoPoint  > (525,TopolPoint,-12,*gMesh);
////	
////	//ELEMENTO 2///
////	TopolPoint[0] = 278;
////	new TPZGeoElRefPattern< TPZGeoPoint  > (526,TopolPoint,-13,*gMesh);
////	
////	TopolPoint[0] = 274;
////	new TPZGeoElRefPattern< TPZGeoPoint  > (527,TopolPoint,-14,*gMesh);
//    
//    TPZVec <int> TopolPoint(1);
//    
//    
//    TopolPoint[0] = 138;
//	new TPZGeoElRefPattern< TPZGeoPoint  > (521,TopolPoint,-15,*gMesh);
//	
//	//ELEMENTO 2///
//	TopolPoint[0] = 125;
//	new TPZGeoElRefPattern< TPZGeoPoint  > (522,TopolPoint,-16,*gMesh);
//	
//	TopolPoint[0] = 148;
//	new TPZGeoElRefPattern< TPZGeoPoint  > (523,TopolPoint,-17,*gMesh);
//    
//    
//    
//    //ELEMENTO 2///
//    TPZVec <int> TopolLine(2);
//	TopolLine[0] =6;	TopolLine[1] = 4;
//	new TPZGeoElRefPattern<TPZGeoLinear> (524,TopolLine,-1,*gMesh);
//
//    TopolLine[0] = 4;	TopolLine[1] = 2;
//	new TPZGeoElRefPattern<TPZGeoLinear> (525,TopolLine,-1,*gMesh);
//    
//    TopolLine[0] = 2;	TopolLine[1] = 5;
//	new TPZGeoElRefPattern<TPZGeoLinear> (526,TopolLine,-1,*gMesh);
//    
//    TopolLine[0] =5;	TopolLine[1] = 7;
//	new TPZGeoElRefPattern<TPZGeoLinear> (527,TopolLine,-1,*gMesh);
//    
//    
//    TopolLine[0] = 312;	TopolLine[1] = 314;
//	new TPZGeoElRefPattern<TPZGeoLinear> (528,TopolLine,-2,*gMesh);
//    
//    TopolLine[0] = 314;	TopolLine[1] = 316;
//	new TPZGeoElRefPattern<TPZGeoLinear> (529,TopolLine,-2,*gMesh);
//    
//    TopolLine[0] = 316;	TopolLine[1] = 315;
//	new TPZGeoElRefPattern<TPZGeoLinear> (530,TopolLine,-2,*gMesh);
//    
//    TopolLine[0] = 315;	TopolLine[1] = 313;
//	new TPZGeoElRefPattern<TPZGeoLinear> (531,TopolLine,-2,*gMesh);
//    
//    
//    TopolPoint[0] = 176;
//	new TPZGeoElRefPattern< TPZGeoPoint  > (532,TopolPoint,-18,*gMesh);
//    TopolPoint[0] = 175;
//	new TPZGeoElRefPattern< TPZGeoPoint  > (533,TopolPoint,-19,*gMesh);
//    
//    
//    //-20 em cima
//    TopolPoint[0] = 316;
//	new TPZGeoElRefPattern< TPZGeoPoint  > (534,TopolPoint,-20,*gMesh);
//    
//
//    //-21 em baixo
//    TopolPoint[0] = 2;
//	new TPZGeoElRefPattern< TPZGeoPoint  > (535,TopolPoint,-21,*gMesh);
//
//    
//    gMesh->BuildConnectivity();
//    
//	
//	// identificando as superficies que terao cond de contorno. Coord z dos 3 nos = 0
////		for(int el = 0; el<numnodes-1; el++) 
////		{
////			Nodefind[el] = gMesh->NodeVec()[el];
////	
////		}
////		Nodefind.Print(std::cout);
////		std::cout.flush();
//	
////	TPZGeoElBC(TPZGeoEl *el,int side,int matid, TPZGeoMesh &mesh);
////	TPZGeoElBC(TPZGeoElSide &elside,int matid, TPZGeoMesh &mesh);
//    ofstream arg("malhaPZ.txt");
//    gMesh->Print(arg);
//    ofstream predio("GeoPredio.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(gMesh,predio,true); 
//    
////	{
////		bool countnodes = false;
////		bool countelements = false;
////		
////		ifstream read (FileName.c_str());
////		
////		while(read)
////		{
////			char buf[1024];
////			read.getline(buf, 1024);
////			std::string str(buf);
////			if(str == "Coordinates") countnodes = true;
////			if(str == "end coordinates") countnodes = false;
////			if(countnodes) numnodes++;
////			
////			if(str == "Elements") countelements = true;
////			if(str == "end elements") countelements = false;
////			if(countelements) numelements++;
////		}
////	}
////	
////	TPZGeoMesh * gMesh = new TPZGeoMesh;
////	
////	gMesh -> NodeVec().Resize(numnodes);
////	
////	TPZVec <int> TopolTriangle(3);
////	
////	const int Qnodes = numnodes;
////	TPZVec <TPZGeoNode> Node(Qnodes);
////	
////	//setting nodes coords
////	int nodeId = 0, elementId = 0, matElId = 1;
////	
////	ifstream read;
////	read.open(FileName.c_str());
////	
////	double nodecoordX , nodecoordY , nodecoordZ ;
////	
////	char buf[1024];
////	read.getline(buf, 1024);
////	read.getline(buf, 1024);
////	std::string str(buf);
////	int in;
////	for(in=0; in<numnodes; in++)
////	{ 
////		read >> nodeId;
////		read >> nodecoordX;
////		read >> nodecoordY;
////		read >> nodecoordZ;
////		Node[nodeId-1].SetNodeId(nodeId);
////		Node[nodeId-1].SetCoord(0,nodecoordX);
////		Node[nodeId-1].SetCoord(1,nodecoordY);
////		Node[nodeId-1].SetCoord(2,nodecoordZ);
////		gMesh->NodeVec()[nodeId-1] = Node[nodeId-1];
////		
////		
////	}
////	
////	{
////		
////		read.close();
////		read.open(FileName.c_str());
////		
////		
////		
////		int l , m = numnodes+5;
////		for(l=0; l<m; l++)
////		{
////			read.getline(buf, 1024);
////		}
////		
////		
////		int el;
////		int matBCid = -1;
////		//std::set<int> ncoordz; //jeitoCaju
////		for(el=0; el<numelements; el++)
////		{
////			read >> elementId;
////			read >> TopolTriangle[0]; //node 1
////			read >> TopolTriangle[1]; //node 2
////			read >> TopolTriangle[2]; //node 3
////			
////			// O GID comeca com 1 na contagem dos nodes, e nao zero como no PZ, assim o node 1 na verdade é o node 0
////			TopolTriangle[0]--;
////			TopolTriangle[1]--;
////			TopolTriangle[2]--;
////			
////			int index = el;
////			
////			TPZGeoEl * tetra = new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (index, TopolTriangle, matElId, *gMesh);
////		}
////		
////		gMesh->BuildConnectivity();
////		// Colocando as condicoes de contorno
////        
//////		for(el=0; el<numelements; el++)
//////		{
//////			TPZManVector <TPZGeoNode,4> Nodefinder(4);
//////			TPZManVector <REAL,3> nodecoord(3);
//////			TPZGeoEl *tetra = gMesh->ElementVec()[el];
//////			// na face z = 0
//////			TPZVec<int> ncoordzVec(0); int sizeOfVec = 0;
//////			for (int i = 0; i < 4; i++) 
//////			{
//////				int pos = tetra->NodeIndex(i);
//////				Nodefinder[i] = gMesh->NodeVec()[pos];
//////				Nodefinder[i].GetCoordinates(nodecoord);
//////				if (nodecoord[2] == 0.)
//////				{
//////					sizeOfVec++;
//////					ncoordzVec.Resize(sizeOfVec);
//////					ncoordzVec[sizeOfVec-1] = pos;
//////				}
//////			}
//////			if(ncoordzVec.NElements() == 3)
//////			{
//////				int lado = tetra->WhichSide(ncoordzVec);
//////				TPZGeoElSide tetraSide(tetra, lado);
//////				TPZGeoElBC(tetraSide,matBCid);		
//////			}
//////		}
////	}
////	
//	
//	
//	// identificando as superficies que terao cond de contorno. Coord z dos 3 nos = 0
//	//	for (int el = 0; el < numnodes-1; el++) 
//	//	{
//	//		Nodefind[el] = gMesh->NodeVec()[el];
//	//
//	//	}
//	//	Nodefind.Print(std::cout);
//	//	std::cout.flush();
//	
//	//TPZGeoElBC(TPZGeoEl *el,int side,int matid, TPZGeoMesh &mesh);
//	//TPZGeoElBC(TPZGeoElSide &elside,int matid, TPZGeoMesh &mesh);
//	
////	ofstream arg("malhaPZ.txt");
////	gMesh->Print(arg);
////	ofstream predio("GeoPredio.vtk");
////	TPZVTKGeoMesh::PrintGMeshVTK(gMesh,predio,true); 
////	
//	return gMesh;
//	
//}
//
//
//TPZGeoMesh *BrazilianTestGeoMesh::MalhaPredio()
//{
//    
//	int numnodes;//=-1;
//	int numelements;//=-1;
//	
//	string FileName;
//	FileName = "gidmeshnodes.txt";
//    ifstream read (FileName.c_str());
//    
//    
//    int nodeId = 0, elementId = 0, matElId = 1;
//    
//    double nodecoordX , nodecoordY , nodecoordZ ;
//    read >> numnodes;
//    
//    TPZGeoMesh * gMesh = new TPZGeoMesh;
//    
//    gMesh -> NodeVec().Resize(numnodes);
//    
//    const int Qnodes = numnodes;
//    TPZVec <TPZGeoNode> Node(Qnodes);
//    
//    for(int in=0; in<numnodes; in++)
//    {
//        read >> nodeId;
//        read >> nodecoordX;
//        read >> nodecoordY;
//        read >> nodecoordZ;
//        Node[nodeId-1].SetNodeId(nodeId);
//        Node[nodeId-1].SetCoord(0,nodecoordX);
//        Node[nodeId-1].SetCoord(1,nodecoordY);
//        Node[nodeId-1].SetCoord(2,nodecoordZ);
//        gMesh->NodeVec()[nodeId-1] = Node[nodeId-1];
//        
//    }
//    
//    string FileName2;
//	FileName2 = "gidmeshelements.txt";
//    ifstream read2 (FileName2.c_str());
//    read2 >> numelements;
//    TPZVec <int> TopolTriangle(3);
//	
//    for(int el=0; el<numelements; el++)
//    {
//        int topol1,topol2,topol3;
//        read2 >> elementId;
//        read2 >> topol1; //node 1
//        read2 >> topol2; //node 2
//        read2 >> topol3; //node 3
//        
//        TopolTriangle[0]=topol1;
//        TopolTriangle[1]=topol2;
//        TopolTriangle[2]=topol3;
//        
//        TopolTriangle[0]--;
//        TopolTriangle[1]--;
//        TopolTriangle[2]--;
//        
//        int index = el;
//        TPZVec<REAL> vec1(3,0.),vec2(3,0.),vec3(3,0.);
//        gMesh->NodeVec()[TopolTriangle[0]].GetCoordinates(vec1);
//        gMesh->NodeVec()[TopolTriangle[1]].GetCoordinates(vec2);
//        gMesh->NodeVec()[TopolTriangle[2]].GetCoordinates(vec3);
//        if(vec1[1] < -76.  || vec2[1]< -76. || vec3[1] < -76. || vec1[1] > 76.  || vec2[1] > 76. || vec3[1] >  76.)
//        {
//            new TPZGeoElRefPattern< TPZGeoTriangle  > (elementId,TopolTriangle,2,*gMesh);
//            continue;
//        }
//        new TPZGeoElRefPattern< TPZGeoTriangle  > (elementId,TopolTriangle,1,*gMesh);
//
//    }
//
//    
//    TPZVec <int> TopolPoint(1);
//    
//    
//    //NOS CENTRAIS
//    TopolPoint[0] = 224;
//	new TPZGeoElRefPattern< TPZGeoPoint  > (935,TopolPoint,-1,*gMesh);
//	
//    
//    //NOS CENTRAIS
//	//ELEMENTO 2///
//	TopolPoint[0] = 252;
//	new TPZGeoElRefPattern< TPZGeoPoint  > (936,TopolPoint,-2,*gMesh);
//	
//    
//    //NOS CENTRAIS
//	TopolPoint[0] = 234;
//	new TPZGeoElRefPattern< TPZGeoPoint  > (937,TopolPoint,-3,*gMesh);
//    
//    //NOS EXTREMIDADES -75 e 75 em X
//	TopolPoint[0] = 124;
//	new TPZGeoElRefPattern< TPZGeoPoint  > (938,TopolPoint,-4,*gMesh);
//	
//	TopolPoint[0] = 480;
//	new TPZGeoElRefPattern< TPZGeoPoint  > (939,TopolPoint,-5,*gMesh);
//    
//    
//    TopolPoint[0] = 490;
//	new TPZGeoElRefPattern< TPZGeoPoint  > (945,TopolPoint,-12,*gMesh);
//    
//    
//    //LINHA SUPERIOR
//    
//    TPZVec <int> TopolLine(2);
//	TopolLine[0] =410;	TopolLine[1] = 566;
//	new TPZGeoElRefPattern<TPZGeoLinear> (940,TopolLine,-6,*gMesh);
//    
//    TopolLine[0] = 1;	TopolLine[1] = 389;
//	new TPZGeoElRefPattern<TPZGeoLinear> (941,TopolLine,-7,*gMesh);
//    
//    //Linhas laterias para impedir as barras em x
//    //Superior esquerdo
//    TopolLine[0] =410;	TopolLine[1] = 433;
//	new TPZGeoElRefPattern<TPZGeoLinear> (941,TopolLine,-8,*gMesh);
//    
//    //superior direito
//    TopolLine[0] = 566;	TopolLine[1] = 567;
//	new TPZGeoElRefPattern<TPZGeoLinear> (942,TopolLine,-9,*gMesh);
//    
//    //Inferior esquerdo
//	TopolLine[0] =1;	TopolLine[1] = 0;
//	new TPZGeoElRefPattern<TPZGeoLinear> (943,TopolLine,-10,*gMesh);
//    
//    //Inferior direito
//    TopolLine[0] = 388;	TopolLine[1] = 389;
//	new TPZGeoElRefPattern<TPZGeoLinear> (944,TopolLine,-11,*gMesh);
//    
//
//    
//    gMesh->BuildConnectivity();
//    
//    ofstream arg("malhaPZ.txt");
//    gMesh->Print(arg);
//    ofstream predio("GeoPredio.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(gMesh,predio,true); 
//    
//       //	
//	return gMesh;
//	
//}

TPZGeoMesh *BrazilianTestGeoMesh::MalhaPredio()
{
    
	int numnodes;//=-1;
	int numelements;//=-1;
	
	string FileName;
	FileName = "gidmeshnodesrefinado.txt";
    ifstream read (FileName.c_str());
    
   // gRefDBase.InitializeRefPatterns();
    
    
    int nodeId = 0, elementId = 0, matElId = 1;
    
    double nodecoordX , nodecoordY , nodecoordZ ;
    read >> numnodes;
    
    TPZGeoMesh * gMesh = new TPZGeoMesh;
    
    gMesh -> NodeVec().Resize(numnodes);
    
    const int Qnodes = numnodes;
    TPZVec <TPZGeoNode> Node(Qnodes);
    
    for(int in=0; in<numnodes; in++)
    {
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
    
    string FileName2;
	FileName2 = "gidmeshelementsrefinado.txt";
    ifstream read2 (FileName2.c_str());
    read2 >> numelements;
    TPZVec <int> TopolTriangle(3);
    TPZVec <int> TopolLine(2);
	
    for(int el=0; el<numelements; el++)
    {
        
        //  if(el<=1124)
        //  {
        int topol1,topol2,topol3;
        read2 >> elementId;
        read2 >> topol1; //node 1
        read2 >> topol2; //node 2
        read2 >> topol3; //node 3
        
        TopolTriangle[0]=topol1;
        TopolTriangle[1]=topol2;
        TopolTriangle[2]=topol3;
        
        TopolTriangle[0]--;
        TopolTriangle[1]--;
        TopolTriangle[2]--;
        
        
       

        
        TPZVec<REAL> vec1(3,0.),vec2(3,0.),vec3(3,0.);
        gMesh->NodeVec()[TopolTriangle[0]].GetCoordinates(vec1);
        gMesh->NodeVec()[TopolTriangle[1]].GetCoordinates(vec2);
        gMesh->NodeVec()[TopolTriangle[2]].GetCoordinates(vec3);
        if(vec1[1] < -75.  || vec2[1]< -75. || vec3[1] < -75. || vec1[1] > 75.  || vec2[1] > 75. || vec3[1] >  75.)
        {
            new TPZGeoElRefPattern< TPZGeoTriangle  > (elementId,TopolTriangle,2,*gMesh);
            //continue;
            cout << "\n elements with mat  = 2 "<< elementId << endl;
            
        }
        else 
        {
             new TPZGeoElRefPattern< TPZGeoTriangle  > (elementId,TopolTriangle,1,*gMesh);
        }
        
    }
    
    TPZVec <int> TopolPoint(1);
    //NO CENTRAL
    TopolPoint[0] = 333;
	new TPZGeoElRefPattern< TPZGeoPoint  > (1063,TopolPoint,-1,*gMesh);
	
    //NOS Laterais direito
	//ELEMENTO 2///
	TopolPoint[0] = 374;
	new TPZGeoElRefPattern< TPZGeoPoint  > (1064,TopolPoint,-2,*gMesh);
    
    //NOS Laterais esquerdo
	TopolPoint[0] = 382;
	new TPZGeoElRefPattern< TPZGeoPoint  > (1065,TopolPoint,-3,*gMesh);
    
    //LINHA SUPERIOR
	TopolLine[0] =1;	TopolLine[1] = 22;
	new TPZGeoElRefPattern<TPZGeoLinear> (1066,TopolLine,-5,*gMesh);
    
  //  TopolLine[0] =599;	TopolLine[1] = 598;
//	new TPZGeoElRefPattern<TPZGeoLinear> (1067,TopolLine,-5,*gMesh);
    
    //LINHA Inferior
    TopolLine[0] = 720;	TopolLine[1] =716;
	new TPZGeoElRefPattern<TPZGeoLinear> (1067,TopolLine,-6,*gMesh);
    
  //  TopolLine[0] = 3;	TopolLine[1] =1;
//	new TPZGeoElRefPattern<TPZGeoLinear> (1069,TopolLine,-6,*gMesh);
    
   // std::set<int> matids;
//	matids.insert(2);
//	matids.insert(1);


    
    gMesh->BuildConnectivity();
    
    for(int ref = 0; ref <0; ref++)
	{
		TPZVec<TPZGeoEl *> tatara;
		int n = gMesh->NElements();
		for(int i = 0; i < n; i++)
		{
			TPZGeoEl * gel = gMesh->ElementVec()[i];
			//if(gel->HasSubElement()) continue;
			gel->Divide(tatara);
		}
	}
    
    ofstream arg("malhaPZ.txt");
    gMesh->Print(arg);
    ofstream predio("malha.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gMesh,predio,true); 
    
    //	
	return gMesh;
	
}

