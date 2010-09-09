//
// Author: MISAEL LUIS SANTANA MANDUJANO.
//
// File:   TPZGenGrid.cpp
//
// Class:  TPZGenGrid
//
// Obs.:   Gera uma malha retangular:
//
// Versao: 10 / 1996.
//
#include "pzgengrid.h"
#include "pzcmesh.h"
#include "pzgmesh.h"

#include "pzgnode.h"
#include "pzgeoel.h"
#include "pzgeoelbc.h"
#include "pzconnect.h"

#include "pzvec.h"
#include "pzstack.h"

#include <fstream>

using namespace std;
TPZGenGrid::TPZGenGrid(TPZVec<int> &nx, TPZVec<REAL> &x0,TPZVec<REAL> &x1, int numl, REAL rot) : fNx(nx), fX0(x0), fX1(x1),
		fDelx(2), fNumLayers(numl), fRotAngle(rot)
{
	fDelx[0] = (x1[0]-x0[0])/(nx[0]);
	fDelx[1] = (x1[1]-x0[1])/(nx[1]);
	fNumNodes= (nx[0]+1)*(nx[1]+1)+(fNumLayers-1)*(nx[0])*(nx[1]+1);
	fElementType = 0;
	fSuperElemlength=0;
	fNumBlocks[0]=0;
	fNumBlocks[1]=0;
}

TPZGenGrid::~TPZGenGrid() {

}


short
TPZGenGrid::Read (TPZGeoMesh &grid) {

	GenerateNodes(grid);
   GenerateElements(grid);
   return 0;
}

void TPZGenGrid::GenerateNodes(TPZGeoMesh &grid) {


// create the geometric nodes
//
	TPZVec<REAL> coor(3,0.);
	int i;
	grid.NodeVec().Resize(fNumNodes);
	for(i=0; i<fNumNodes; i++) {
		Coord(i,coor);
		grid.NodeVec()[i].Initialize(coor,grid);
	}
}
//
// create the geometric elements (retangular)

void TPZGenGrid::GenerateElements(TPZGeoMesh &grid) {

	int num_rectangles=fNx[0]*fNx[1]*fNumLayers;
	TPZVec<int> nos(9);
	if(fElementType == 0) nos.Resize(4);
   int i, index;
	for(i=0; i<num_rectangles; i++) {
		ElementConnectivity(i,nos);
		if(fElementType == 0) {
      grid.CreateGeoElement(EQuadrilateral,nos, 1, index,0);
		} else if(fElementType == 1) {
      grid.CreateGeoElement(ETriangle,nos, 1, index,0);  
			nos[1] = nos[2];
			nos[2] = nos[3];
			grid.CreateGeoElement(ETriangle,nos, 1, index,0);  
		} else if(fElementType == 2) {
      std::cout << __PRETTY_FUNCTION__ << " - Quadratic interpolation is not available";
      exit(-1);        
			grid.CreateGeoElement(EQuadrilateral,nos, 1, index,0);  
    }
	}
	grid.BuildConnectivity();
}

void TPZGenGrid::Coord(int i, TPZVec<REAL> &coor) {
	int ix;
	int iy;
	int ilayer;
	if(fElementType == 0 || fElementType == 1) {
		if(i < (fNx[0]+1)*(fNx[1]+1)) {
			ilayer = 0;
		} else {
			ilayer = (i-(fNx[0]+1)*(fNx[1]+1))/((fNx[0])*(fNx[1]+1))+1;
		}
	} else if(fElementType == 2) {
		if(i < (2*fNx[0]+1)*(2*fNx[1]+1) ) {
			ilayer = 0;
		} else {
			ilayer = (i-(2*fNx[0]+1)*(2*fNx[1]+1))/((2*fNx[0])*(2*fNx[1]+1))+1;
		}
	}
	REAL Rot = fRotAngle*ilayer;
	if(ilayer != 0 && (fElementType == 0 || fElementType == 1)) {
		i -= ((fNx[0]+1)*(fNx[1]+1))+(ilayer-1)*((fNx[0])*(fNx[1]+1));
	} else if (ilayer != 0 && fElementType == 2) {
		i -= ((2*fNx[0]+1)*(2*fNx[1]+1))+(ilayer-1)*((2*fNx[0])*(2*fNx[1]+1));
	}
	if(ilayer == 0) {
		if(fElementType == 0 || fElementType == 1) {
			ix = i%(fNx[0]+1);
			iy = i/(fNx[0]+1);
		} else if(fElementType == 2) {
			ix = i%(2*fNx[0]+1);
			iy = i/(2*fNx[0]+1);
		}
	} else {
		if(fElementType == 0 || fElementType == 1) {
			ix = i%(fNx[0])+1;
			iy = i/(fNx[0]);
		} else if(fElementType == 2) {
			ix = i%(2*fNx[0])+1;
			iy = i/(2*fNx[0]);
		}
	}
//cout << "Coord i = " << i << " ix = " << ix << " iy = " << iy << " layer = " << ilayer << endl;
//cout.flush();
   REAL coorold[2];
   coorold[0] = fX0[0]+fDelx[0]*ix;
   coorold[1] = fX0[1]+fDelx[1]*iy;
	coor[0] = fX0[0]+(coorold[0]-fX0[0])*cos(Rot);
	coor[2] = coorold[0]*sin(Rot);
//	coor[0] = coorold[0];
	coor[1] = coorold[1];
//   coor[2] = 0.;// 0.1*(coorold[0]-fX0[0])*(fX1[0]-coorold[0]);
}

void TPZGenGrid::ElementConnectivity(int i, TPZVec<int> &rectangle_nodes){
	int xel = i%(fNx[0]);
	int yel = (i/(fNx[0]))%(fNx[1]);
	int layer = i/(fNx[0]*fNx[1]);
//cout << "ElConnectivity : xel = " << xel << " yel = " << yel << " layer = " << layer << endl;
//cout.flush();
   if(fElementType == 0 || fElementType == 1) {
		rectangle_nodes[0] = GlobalI(xel,yel,layer);
		rectangle_nodes[1] = GlobalI(xel+1,yel,layer);
		rectangle_nodes[2] = GlobalI(xel+1,yel+1,layer);
		rectangle_nodes[3] = GlobalI(xel,yel+1,layer);
//cout << "ElConnectivity : " << rectangle_nodes[0] << ' '<< rectangle_nodes[1] << ' '<<rectangle_nodes[2] << ' '<<rectangle_nodes[3] << endl;
//cout.flush();
   } else if(fElementType == 2) {
		rectangle_nodes[0] = GlobalI(2*xel,2*yel,layer);
		rectangle_nodes[1] = GlobalI(2*xel+2,2*yel,layer);
		rectangle_nodes[2] = GlobalI(2*xel+2,2*yel+2,layer);
		rectangle_nodes[3] = GlobalI(2*xel,2*yel+2,layer);
		rectangle_nodes[4] = GlobalI(2*xel+1,2*yel,layer);
		rectangle_nodes[5] = GlobalI(2*xel+2,2*yel+1,layer);
		rectangle_nodes[6] = GlobalI(2*xel+1,2*yel+2,layer);
		rectangle_nodes[7] = GlobalI(2*xel,2*yel+1,layer);
		rectangle_nodes[8] = GlobalI(2*xel+1,2*yel+1,layer);
   }
}



void
TPZGenGrid::Print( char *name , ostream &out  )
{

	out<<"\n"<<name<<"\n";
   out << "element type = " << fElementType << endl;
	out << "Number of divisions " << fNx[0] << ' ' << fNx[1] << endl;
	out << "Corner Coordinates " << endl << fX0[0] << ' ' << fX0[1] << endl;
	out << fX1[0] << ' ' << fX1[1] << endl;

}

void TPZGenGrid::SetBC(TPZGeoMesh*g, int side, int bc) {
	int ielfirst = 0;
	int iellast = 0;
	int ielinc;
	int layer;
	int iel;
	TPZGeoEl *gel;
	int elementside = side + 4;
	for(layer=0; layer<fNumLayers; layer++) {
		switch(side) {
		case 0:
			ielfirst = layer*fNx[0]*fNx[1];
			iellast = ielfirst+fNx[0];
			ielinc = 1;
			break;
		case 1:
			ielfirst = layer*fNx[0]*fNx[1]+fNx[0]-1;
			iellast = (layer+1)*fNx[0]*fNx[1];
			ielinc = fNx[0];
			break;
		case 2:
			ielfirst = layer*fNx[0]*fNx[1]+fNx[0]*(fNx[1]-1);
			ielinc = 1;
			iellast = (layer+1)*fNx[0]*fNx[1];
			break;
		case 3:
			ielfirst = layer*fNx[0]*fNx[1];
			ielinc = fNx[0];
			iellast = ielfirst+fNx[1]*ielinc;
			break;
		}
		if(fElementType == 1) {
			elementside = side+3;
			ielfirst *= 2;
			iellast *= 2;
			ielinc *= 2;
			if(side > 1) {
				ielfirst += 1;
				iellast += 1;
				side -= 1;
				elementside--;
			}
		}
		for(iel = ielfirst; iel<iellast; iel += ielinc) {
			gel = g->ElementVec()[iel];
			TPZGeoElBC(gel,elementside,bc,*g);
//			gel->SetSide(side,bc);
		}
	}

/*	VoidPtrVec ElementVec;
	TPZVec<int> Sides;
	TPZVec<int> cornernodes(4);
	for(layer =0; layer< fNumLayers; layer++) {
		if(fElementType == 0 || fElementType == 1) {
			cornernodes[0] = GlobalI(0,0,layer);
			cornernodes[1] = GlobalI(fNx[0],0,layer);
			cornernodes[2] = GlobalI(fNx[0],fNx[1],layer);
			cornernodes[3] = GlobalI(0,fNx[1],layer);
		} else if (fElementType == 2) {
			cornernodes[0] = GlobalI(0,0,layer);
			cornernodes[1] = GlobalI(2*fNx[0],0,layer);
			cornernodes[2] = GlobalI(2*fNx[0],2*fNx[1],layer);
			cornernodes[3] = GlobalI(0,2*fNx[1],layer);
		}
cout << "SetBC cornernodes[0] = " << cornernodes[0] << " cornernodes[1] = " << cornernodes[1] << " cornernodes[2] = " << cornernodes[2] << " cornernodes[3] = " << cornernodes[3] << endl; 
cout.flush();
		g->GetBoundaryElements(cornernodes[side],cornernodes[(side+1)%4],ElementVec, Sides);
		int numel = ElementVec.capacity();
cout << "Boundary elements ";
		for(int el=0; el<numel; el++) {
			TPZGeoEl *gel = (TPZGeoEl *) ElementVec[el];
			if(gel) {
cout << gel->Id() << " ";
				gel->SetSide(Sides[el],bc);
			}
		}
cout << endl;
	}
*/
}

void TPZGenGrid::SetBC(TPZGeoMesh*g, TPZVec<REAL> &start, TPZVec<REAL> &end, int bc){

	TPZGeoNode *gn1 = g->FindNode(start);
	TPZGeoNode *gn2 = g->FindNode(end);

	TPZStack<TPZGeoEl *> ElementVec;
	TPZStack<int> Sides;
	g->GetBoundaryElements(gn1->Id(),gn2->Id(),ElementVec, Sides);
	int numel = ElementVec.NElements();
	for(int el=0; el<numel; el++) {
		TPZGeoEl *gel = (TPZGeoEl *) ElementVec[el];
		if(gel) {
			TPZGeoElBC(gel,Sides[el],bc,*g);
      }
	}
}

int TPZGenGrid::GlobalI(int ix, int iy, int layer) {
	if(layer == 0 || ix == 0) {
		if(fElementType == 0 || fElementType == 1) {
			return ix+iy*(fNx[0]+1);
		} else if(fElementType == 2) {
			return ix+iy*(2*fNx[0]+1);
		}
	} else {
		if(fElementType == 0 || fElementType == 1) {
			return (fNx[0]+1)*(fNx[1]+1)+(fNx[0])*(fNx[1]+1)*(layer-1)+ix-1+iy*(fNx[0]);
		} else if(fElementType == 2) {
			return (2*fNx[0]+1)*(2*fNx[1]+1)+(2*fNx[0])*(2*fNx[1]+1)*(layer-1)+ix-1+iy*(2*fNx[0]);
		}
	}
	return 0;
}


int TPZGenGrid::ElemId(int iel,int jel, int layer){
  return (fElementType == 0 || fElementType == 2)? ( iel*fNx[0]+jel+fNx[0]*fNx[1]*layer ):(iel*2*fNx[0]+jel+layer*2*fNx[0]*fNx[1]) ;
   }

int TPZGenGrid::SuperElemId(int iel,int jel){
  return ( iel*fNumBlocks[0]+jel );
}

REAL
TPZGenGrid::Distance(TPZVec<REAL> &x1,TPZVec<REAL> &x2){
	REAL l1,l2;
	l1=x1[0]-x2[0];
	l2=x1[1]-x2[1];
	return( sqrt(l1*l1+l2*l2) );
}

void TPZGenGrid::SetElementType(int type) {
	fElementType = type;
	if(fElementType == 0 || fElementType == 1) {
		fNumNodes = (fNx[0]+1)*(fNx[1]+1)+(fNumLayers-1)*(fNx[0])*(fNx[1]+1);
		fDelx[0] = (fX1[0]-fX0[0])/(fNx[0]);
		fDelx[1] = (fX1[1]-fX0[1])/(fNx[1]);
	} else if(fElementType == 2) {
		fNumNodes = (2*fNx[0]+1)*(2*fNx[1]+1)+(fNumLayers-1)*(2*fNx[0])*(2*fNx[1]+1);
		fDelx[0] = (fX1[0]-fX0[0])/(2*fNx[0]);
		fDelx[1] = (fX1[1]-fX0[1])/(2*fNx[1]);
	}
}

