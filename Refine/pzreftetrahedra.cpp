/**
 * @file
 * @brief Contains the implementation of the TPZRefTetrahedra methods. 
 */
#include "pzrefpyram.h"
#include "pzgeopyramid.h"
#include "pzshapepiram.h"
#include "pzreftetrahedra.h"
#include "pzgeotetrahedra.h"
#include "pzshapetetra.h"
#include "pzgeoelside.h"
#include "pzgeoel.h"
#include "pzgmesh.h"

using namespace pzshape;
using namespace std;

namespace pzrefine {
	
	
	static int nsubeldata[15] = {1,1,1,1,3,3,3,3,3,3,7,7,7,7,11};
	
	static int subeldata[15][11][2] = {//TAMANHO DISTINTO
		//lados do pai:{{conectividades dos filhos}}
		/*00*/{{0,0}},//os lados dos filhos formam
		/*01*/{{1,1}},//uma particao do lado do pai
		/*02*/{{2,2}},
		/*03*/{{3,3}},
		/*04*/{{0,4},{0,1},{1,4}},
		/*05*/{{1,5},{1,2},{2,5}},
		/*06*/{{2,6},{2,0},{0,6}},
		/*07*/{{0,7},{0,3},{3,7}},
		/*08*/{{1,8},{1,3},{3,8}},
		/*09*/{{2,9},{2,3},{3,9}},
		/*10*/{{0,5},{1,6},{2,4},{0,10},{1,10},{2,10},{5,15}},
		/*11*/{{0,8},{1,7},{3,4},{0,11},{1,11},{3,11},{4,14}},
		/*12*/{{1,9},{2,8},{3,5},{1,12},{2,12},{3,12},{5,17}},
		/*13*/{{0,9},{2,7},{3,6},{0,13},{2,13},{3,13},{4,16}},
		/*14*/{{0,14},{1,14},{2,14},{3,14},{4,18},{5,18},{4,13},{0,12},{1,13},{2,11},{3,10}}
	};
	
	
	static int MidSideNodes[6][2]  = { 
		{0,1},
		{1,2},
		{2,0},
		{3,0},
		{3,1},
		{3,2} 
	};
	
	static REAL MidCoord[6][3] = { 
		{.5,0.,0.},
		{.5,.5,0.},
		{0.,.5,0.},
		{0.,0.,.5},
		{.5,0.,.5},
		{0.,.5,.5} };
	
	/**
	 * define as conectividades entre sub-elementos
	 * linha i é filho i, {a,b,c} = {lado do filho atual,
	 * irmão vizinho,lado do vizinho}
	 */
	const int NumInNeigh = 16;
	static int InNeigh[6][NumInNeigh][3] = { 
		{{1,1,0},{2,4,3},{3,4,4},{5,5,6},{8,4,9},{9,4,12},{12,4,17},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
		{{0,5,1},{2,2,1},{3,3,1},{6,5,10},{7,5,5},{9,5,9},{13,5,14},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
		{{0,0,2},{1,5,4},{3,3,2},{4,5,11},{7,4,7},{8,5,12},{11,5,16},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
		{{0,0,3},{1,4,1},{2,4,2},{4,4,10},{5,4,6},{6,4,11},{10,4,15},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
		{{0,0,1},{1,5,0},{2,5,3},{3,5,2},{4,3,0},{5,1,7},{6,5,8},{7,5,7},{8,0,5},{9,0,8},{10,3,4},{11,3,6},{12,0,9},{13,5,13},{15,3,10},{17,0,12}},
		{{0,1,3},{1,4,0},{2,2,0},{3,2,3},{4,1,2},{5,4,5},{6,4,8},{7,2,7},{8,3,5},{9,1,9},{10,1,6},{11,2,4},{12,2,8},{13,4,13},{14,1,13},{16,2,11}} 
	};
	
	/**
	 * define os cantos locais dos fihos
	 */
	static int CornerSons[6][5] = { 
		{0,4,6,7,-1},
		{4,1,5,8,-1},
		{6,5,2,9,-1},
		{7,8,9,3,-1},
		{4,8,9,6,7},
		{8,4,6,9,5} 
	};
	
	static REAL buildt[6][19][4][3] = {//por colunas
		/*S0*/{
			/*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
			/*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
			/*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
			/*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
			/*04*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
			/*05*/{{-0.25,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
			/*06*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
			/*07*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
			/*08*/{{-0.25,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
			/*09*/{{-0.25,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
			/*10*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.,0.}},
			/*11*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.,0.}},
			/*12*/{{-0.5,0.5,0},{-0.5,0,0.5},{0.,0.,0.},{0.5,0,0}},
			/*13*/{{0.5,0,0.},{0,0.5,0.},{0.,0.,0.},{0,0,0.}},
			/*14*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,.5},{0.,0.,0.}},
			/*15*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
			/*16*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
			/*17*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
			/*18*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}},
		/*S1*/{
			/*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
			/*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{1,0.,0.}},
			/*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
			/*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
			/*04*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
			/*05*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
			/*06*/{{0.,-0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
			/*07*/{{0.,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
			/*08*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
			/*09*/{{-0.25,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
			/*10*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.5,0.,0.}},
			/*11*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.5,0.,0.}},
			/*12*/{{0.5,0,0.},{0,0.5,0.},{0.,0.,0.},{0,0,0.}},
			/*13*/{{0,0.5,0},{0,0,0.5},{0.,0.,0.},{0.5,0,0}},
			/*14*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,.5},{.5,0.,0.}},
			/*15*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
			/*16*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
			/*17*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
			/*18*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}},
		/*S2*/{
			/*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
			/*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
			/*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,1,0.}},
			/*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
			/*04*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
			/*05*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
			/*06*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
			/*07*/{{0,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
			/*08*/{{0,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
			/*09*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
			/*10*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.5,0.}},
			/*11*/{{0.5,0.,0.},{0.,0.,0.5},{0.,0.,0.},{0.,0.5,0.}},
			/*12*/{{0.5,0,0.},{0,0.5,0.},{0.,0.,0.},{0.5,0,0.}},
			/*13*/{{0.5,0,0.},{0,0.5,0.},{0.,0.,0.},{0.5,0,0.}},
			/*14*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,.5},{0.,.5,0.}},
			/*15*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
			/*16*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
			/*17*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
			/*18*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}},
		/*S3*/{
			/*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
			/*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
			/*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
			/*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,1}},
			/*04*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
			/*05*/{{0.25,0,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
			/*06*/{{-0.25,0,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
			/*07*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
			/*08*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
			/*09*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
			/*10*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.,0.5}},
			/*11*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.5,0.}},
			/*12*/{{0.5,0,0.},{0,0.5,0.},{0.,0.,0.},{0,0.5,0.}},
			/*13*/{{0.5,0,0.},{0,0.5,0.},{0.,0.,0.},{0,0.5,0.}},
			/*14*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,.5},{0.,0.,.5}},
			/*15*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
			/*16*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
			/*17*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
			/*18*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}},
		/*S4*/{
			/*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
			/*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
			/*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
			/*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
			/*04*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
			/*05*/{{0.,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
			/*06*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
			/*07*/{{0.,-0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
			/*08*/{{0.25,-0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
			/*09*/{{-0.25,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
			/*10*/{{-0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
			/*11*/{{-0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
			/*12*/{{-0.25,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
			/*13*/{{0.,0.,0.25},{-0.25,0.25,0.},{0.,0.,0.},{0.25,0.25,0.25}},
			/*14*/{{0.,0.5,0.},{-0.5,0.5,0.},{0.,0.,0.},{0.5,0.,0.}},
			/*15*/{{-0.5,0.5,0.},{-0.5,0.,0.},{0.,0.,0.},{0.5,0.,0.5}},
			/*16*/{{0.,0.5,0.},{-0.5,0.5,0.},{0.,0.,0.},{0.5,0.,0.}},
			/*17*/{{-0.5,0.5,0.},{-0.5,0.,0.5},{0.,0.,0.},{0.5,0.,0.}},
			/*18*/{{0.,0.,.25},{-.25,.25,0.},{-.25,-.25,.25},{.25,.25,.25}}},
		/*S5*/{
			/*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
			/*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
			/*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
			/*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
			/*04*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
			/*05*/{{0.,-0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
			/*06*/{{-0.25,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
			/*07*/{{0.,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
			/*08*/{{-0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
			/*09*/{{0.25,-0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
			/*10*/{{0.,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
			/*11*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
			/*12*/{{0.,-0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
			/*13*/{{0.,0.,-0.25},{-0.25,0.25,0.},{0.,0.,0.},{0.25,0.25,0.25}},
			/*14*/{{0.,0.,-0.5},{0.,0.5,-0.5},{0.,0.,0.},{0.5,0.,0.5}},
			/*15*/{{-0.5,0.5,0.},{0.,0.5,0.},{0.,0.,0.},{0.5,0.,0.}},
			/*16*/{{0.,0.,-0.5},{0.5,0.,-0.5},{0.,0.,0.},{0.,0.5,0.5}},
			/*17*/{{0.5,0.,0.},{0.5,-0.5,0.},{0.,0.,0.},{0.,0.5,0.}},
			/*18*/{{0.,0.,-.25},{-.25,.25,0.},{.25,.25,-.25},{.25,.25,.25}}}
	};
	
	//Into Divides is necesary to consider the connectivity with the all neighboards
	void TPZRefTetrahedra::Divide(TPZGeoEl *geo,TPZVec<TPZGeoEl *> &SubElVec) {
		int i;
		SubElVec.Resize(NSubEl);
		if(geo->HasSubElement()) {
			for(i=0;i<NSubEl;i++) SubElVec[i] = geo->SubElement(i);
			return;//If exist fSubEl return this sons
		}
		int j,sub,matid=geo->MaterialId(),index;
		int np[TPZShapeTetra::NSides];//guarda conectividades dos 8 subelementos
		for(j=0;j<TPZShapeTetra::NCornerNodes;j++) np[j] = geo->NodeIndex(j);
		for(sub=TPZShapeTetra::NCornerNodes;sub<10;sub++) {
			NewMidSideNode(geo,sub,index);
			np[sub] = index;
		}
		// creating new subelements
		for (i=0;i<4;i++){
			TPZManVector<int> cornerindexes(TPZShapeTetra::NCornerNodes);
			for(int j=0;j<TPZShapeTetra::NCornerNodes;j++) cornerindexes[j] = np[CornerSons[i][j]];
			TPZGeoEl *t3sub = geo->Mesh()->CreateGeoElement(ETetraedro,cornerindexes,matid,index);
			geo->SetSubElement(i,t3sub);
			t3sub->SetFather(geo);
			t3sub->SetFather(geo->Index());
			SubElVec[i] = t3sub;
		}
		for (;i<6;i++){
			TPZManVector<int> cornerindexes(TPZShapePiram::NCornerNodes);
			for(int j=0;j<TPZShapePiram::NCornerNodes;j++) 
				cornerindexes[j] = np[CornerSons[i][j]];
			TPZGeoEl *pi3sub = geo->Mesh()->CreateGeoElement(EPiramide,cornerindexes,matid,index);
			geo->SetSubElement(i,pi3sub);
			pi3sub->SetFather(geo);
			pi3sub->SetFather(geo->Index());
			SubElVec[i] = pi3sub;
		}
		for(i=0;i<NSubEl;i++) {//conectividades entre os filhos : viz interna
			for(j=0;j<NumInNeigh;j++) {        //lado do subel  numero do filho viz.             lado do viz.
				int elside = InNeigh[i][j][0];//lado do subel
				if(elside == -1) break;
				geo->SubElement(i)->SetNeighbour(elside,TPZGeoElSide(geo->SubElement(InNeigh[i][j][1]),InNeigh[i][j][2]));
			}
		}
		geo->SetSubElementConnectivities();
	}
	
	void TPZRefTetrahedra::NewMidSideNode(TPZGeoEl *gel,int side,int &index) {
		
		MidSideNodeIndex(gel,side,index);
		if(index < 0) {
			TPZGeoElSide gelside = gel->Neighbour(side);
			if(gelside.Element()) {
				while(gelside.Element() != gel) {
					gelside.Element()->MidSideNodeIndex(gelside.Side(),index);
					if(index!=-1) return;
					gelside = gelside.Neighbour();
				}	
			}
			TPZVec<REAL> par(3,0.);
			TPZVec<REAL> coord(3,0.);
			if(side < TPZShapeTetra::NCornerNodes) {
				index = gel->NodeIndex(side); 
				return;
			}
			//aqui side = 8 a 26
			side-=TPZShapeTetra::NCornerNodes;//0,1,..,18
			par[0] = MidCoord[side][0];
			par[1] = MidCoord[side][1];
			par[2] = MidCoord[side][2];
			gel->X(par,coord);
			index = gel->Mesh()->NodeVec().AllocateNewElement();
			gel->Mesh()->NodeVec()[index].Initialize(coord,*gel->Mesh());
		}
	}
	
	void TPZRefTetrahedra::MidSideNodeIndex(TPZGeoEl *gel,int side,int &index) {
		index = -1;
		if(side<0 || side>TPZShapeTetra::NSides-1) {
			PZError << "TPZRefTetrahedra::MidSideNodeIndex. Bad parameter side = " << side << endl;
			return;
		}
		//sides 0 a 7
		if(side<TPZShapeTetra::NCornerNodes) {//o nó medio do lado 0 é o 0 etc.
			index = (gel)->NodeIndex(side);
			return; 
		}
		//o nó medio da face é o centro e o nó medio do centro é o centro
		//como nó de algum filho se este existir
		//caso tenha filhos é o canto de algum filho, se não tiver filhos retorna -1
		if(gel->HasSubElement()) {
			side-=TPZShapeTetra::NCornerNodes;
			index=(gel->SubElement(MidSideNodes[side][0]))->NodeIndex(MidSideNodes[side][1]);
		}
	}
	
	void TPZRefTetrahedra::GetSubElements(TPZGeoEl *father,int side, TPZStack<TPZGeoElSide> &subel){
		
		subel.Resize(0);
		if(side<0 || side>TPZShapeTetra::NSides || !father->HasSubElement()){
			PZError << "TPZRefTetrahedra::GetSubelements2 called with error arguments\n";
			return;
		}
		int nsub = NSideSubElements(side);//nsubeldata[side];
		for(int i=0;i<nsub;i++)
			subel.Push(TPZGeoElSide(father->SubElement(subeldata[side][i][0]),
									subeldata[side][i][1]));
	}
	
	int TPZRefTetrahedra::NSideSubElements(int side) {  
		if(side<0 || side>TPZShapeTetra::NSides-1){
			PZError << "TPZRefTetrahedra::NSideSubelements2 called with error arguments\n";
			return -1;
		}
		return nsubeldata[side];
	}
	
	//int TPZRefTetrahedra::NSideSubElements(int side) {
	//  if(side < 0 || side > 26) {
	//    PZError << "TPZRefTetrahedra::NSideSubElements called for side " << side << endl;
	//    return 0;
	//  }
	//  if(side==26) return 8;//centro
	//  if(side>19 && side<26) return 4;//faces
	//  if(side>7) return 2;//lados
	//  return 1;//cantos
	//}
	
	
	TPZTransform TPZRefTetrahedra::GetTransform(int side,int whichsubel){
		if(side<0 || (whichsubel < 4 && side>TPZShapeTetra::NSides-1) ||
		   (side >TPZShapePiram::NSides-1)){
			PZError << "TPZRefTetrahedra::GetTransform side out of range or father null\n";
			return TPZTransform(0,0);
		}
		int smalldim;
		if(whichsubel <4) smalldim = TPZShapeTetra::SideDimension(side);
		else smalldim = TPZShapePiram::SideDimension(side);
		int fatherside = FatherSide(side,whichsubel);
		int largedim = TPZShapeTetra::SideDimension(fatherside);
		TPZTransform trans(largedim,smalldim);
		int i,j;
		for(i=0; i<largedim; i++) {
			for(j=0; j<smalldim; j++) {
				trans.Mult()(i,j) = buildt[whichsubel][side][j][i];
			}
			trans.Sum() (i,0) = buildt[whichsubel][side][3][i];
		}
		return trans;
	}
	
	static int fatherside[6][19] = {
		/*00*/{0,4,6,7,4,10,6,7,11,13,10,11,14,13,14,-1,-1,-1,-1},
		/*01*/{4,1,5,8,4,5,10,11,8,12,10,11,12,14,14,-1,-1,-1,-1},
		/*02*/{6,5,2,9,10,5,6,13,12,9,10,14,12,13,14,-1,-1,-1,-1},
		/*03*/{7,8,9,3,11,12,13,7,8,9,14,11,12,13,14,-1,-1,-1,-1},
		/*04*/{4,8,9,6,7,11,12,13,10,11,11,13,13,14,11,14,13,14,14},
		/*05*/{8,4,6,9,5,11,10,13,12,12,10,10,12,14,14,10,14,12,14},
	};
	
	// static int fatherside2[4][19] = {//tetraedro com pai pirâmide
	// /*06*/{9,5,13,10,14,13,18,14,14,18,18,14,18,18,18,-1,-1,-1,-1},
	// /*07*/{6,10,11,13,15,15,15,13,18,18,15,18,18,18,18,-1,-1,-1,-1},
	// /*08*/{12,13,7,11,18,13,16,16,18,16,18,18,18,16,18,-1,-1,-1,-1},
	// /*09*/{13,9,12,8,18,17,18,13,17,17,18,18,17,18,18,-1,-1,-1,-1} };
	
	int TPZRefTetrahedra::FatherSide(int side,int whichsubel){
		
		return fatherside[whichsubel][side];
	}
	
};
