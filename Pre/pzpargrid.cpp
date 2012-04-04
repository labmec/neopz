/**
 * @file
 * @brief Contains the implementation of the TPZGenPartialGrid methods. 
 */

#include "pzpargrid.h"
#include "pzcmesh.h"
#include "pzgmesh.h"

#include "pzgnode.h"
#include "pzgeoel.h"
#include "pzgeoelbc.h"

#include "pzvec.h"
#include "pzstack.h"

using namespace std;

TPZGenPartialGrid::TPZGenPartialGrid(TPZVec<int> &nx, TPZVec<int> &rangex, TPZVec<int> &rangey, TPZVec<REAL> &x0,TPZVec<REAL> &x1) : fNx(nx), fX0(x0), fX1(x1),
fDelx(2)
{
	fDelx[0] = (x1[0]-x0[0])/(nx[0]);
	fDelx[1] = (x1[1]-x0[1])/(nx[1]);
	fNumNodes= (nx[0]+1)*(nx[1]+1);
	fElementType = 0;
	fRangex = rangex;
	fRangey = rangey;
}

TPZGenPartialGrid::~TPZGenPartialGrid() {
	
}

int TPZGenPartialGrid::NodeIndex(int x, int y){
	return y*(fNx[0]+1)+x;
}

int TPZGenPartialGrid::ElementIndex(int x, int y) {
	return y*fNx[0]+x;
}

short
TPZGenPartialGrid::Read (TPZGeoMesh &grid) {
	
	TPZVec<REAL> coor(2,0.);
	int i,j;
	int index;
	
	if(grid.NodeVec().NElements() < fNumNodes) grid.NodeVec().Resize(fNumNodes);
	for(i=fRangex[0]; i<fRangex[1]+1; i++) {
		for(j = fRangey[0]; j<fRangey[1]+1; j++){
			index = NodeIndex(i,j);
			Coord(index,coor);
			grid.NodeVec()[index].Initialize(coor,grid);
		}
	}
	//
	//numbering the nodes
	//
	// create the geometric elements (retangular)
	
	TPZVec<int> nos(4);
	for(i=fRangex[0]; i<fRangex[1]; i++) {
		for(j=fRangey[0]; j<fRangey[1]; j++) {
			index = ElementIndex(i,j);
			ElementConnectivity(index,nos);
			if(fElementType == 0) {
				int index;   
				grid.CreateGeoElement(EQuadrilateral,nos,1,index);
			} else {
				int index;   
				grid.CreateGeoElement(ETriangle,nos,1,index);
				nos[1] = nos[2];
				nos[2] = nos[3];
				grid.CreateGeoElement(ETriangle,nos,1,index);
			}
		}
	}
	grid.BuildConnectivity();
	return 0;
	
}

void TPZGenPartialGrid::Coord(int i, TPZVec<REAL> &coor) {
	int ix = i%(fNx[0]+1);
	int iy = i/(fNx[0]+1);
	coor[0] = fX0[0]+fDelx[0]*ix;
	coor[1] = fX0[1]+fDelx[1]*iy;
}

void TPZGenPartialGrid::ElementConnectivity(int i, TPZVec<int> &rectangle_nodes){
	int xel = i%(fNx[0]);
	int yel = i/(fNx[0]);
	rectangle_nodes[0] = yel*(fNx[0]+1)+xel;
	rectangle_nodes[1] = yel*(fNx[0]+1)+xel+1;
	rectangle_nodes[2] = (yel+1)*(fNx[0]+1)+xel+1;
	rectangle_nodes[3] = (yel+1)*(fNx[0]+1)+xel;
}

void TPZGenPartialGrid::Print( char *name , ostream &out  )
{
	
	out<<"\n"<<name<<"\n";
	out << "Number of divisions " << fNx[0] << ' ' << fNx[1] << endl;
	out << "Corner Coordinates " << endl << fX0[0] << ' ' << fX0[1] << endl;
	out << fX1[0] << ' ' << fX1[1] << endl;
	
}

void TPZGenPartialGrid::SetBC(TPZGeoMesh &g, int side, int bc) {
	if(side == 0 && fRangey[0] != 0) return;
	if(side == 1 && fRangex[1] != fNx[0]) return;
	if(side == 2 && fRangey[1] != fNx[1]) return;
	if(side == 3 && fRangex[0] != 0) return;
	TPZStack<TPZGeoEl*> ElementVec;
	TPZStack<int> Sides;
	TPZVec<int> cornernodes(4);
	cornernodes[0] = NodeIndex(fRangex[0],fRangey[0]);
	cornernodes[1] = NodeIndex(fRangex[1],fRangey[0]);
	cornernodes[2] = NodeIndex(fRangex[1],fRangey[1]);
	cornernodes[3] = NodeIndex(fRangex[0],fRangey[1]);
	g.GetBoundaryElements(cornernodes[side],cornernodes[(side+1)%4],ElementVec, Sides);
	int numel = ElementVec.NElements();
	for(int el=0; el<numel; el++) {
		TPZGeoEl *gel = (TPZGeoEl *) ElementVec[el];
		if(gel) {
			TPZGeoElBC(gel,Sides[el],bc);
		}
	}
}

