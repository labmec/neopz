/**
 * @file
 * @brief Contains the implementation of the TPZGenGrid methods. 
 */
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
fDelx(2), fGeometricProgression(2,1.), fNumLayers(numl), fRotAngle(rot)
{
	fDelx[0] = (x1[0]-x0[0])/(nx[0]);
	fDelx[1] = (x1[1]-x0[1])/(nx[1]);
	fNumNodes= (nx[0]+1)*(nx[1]+1)+(fNumLayers-1)*(nx[0])*(nx[1]+1);
	fElementType = 0;
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
            DebugStop();
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
    REAL coorold[2] = {fX0[0],fX0[1]};
    REAL elsize[2] = {fDelx[0],fDelx[1]};
    for (int i=0; i<ix; i++) {
        coorold[0] += elsize[0];
        elsize[0] *= fGeometricProgression[0];
    }
    for (int j=0; j<iy; j++) {
        coorold[1] += elsize[1];
        elsize[1] *= fGeometricProgression[1];
    }
	//    coorold[1] = fX0[1]+fDelx[1]*iy;
    
    // rotate along the y axis
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

/// compute the geometric progression such that the first elements have this size
REAL TPZGenGrid::GeometricProgression(REAL minsize, REAL domainsize, int numdiv)
{
    REAL progression = 1.;
    REAL factor = domainsize/minsize;
    REAL func = 0.;//pow(progression[idim],fNx[idim])-1.+factor[idim]*(1.-progression[idim]);
    REAL nextsize = 1.;
    for (int i=0; i<numdiv; i++) {
        func += nextsize;
        nextsize *= progression;
    }
    func -= factor;
    int iter = 0;
    int maxiter = 200;
    while (fabs(func/factor) >= 1.e-10 && iter < 200) {
        REAL dfunc = 0.;// fNx[idim]*pow(progression[idim], fNx[idim]-1)-factor[idim];
        func = 0.;
        nextsize = 1.;
        for (int i=0; i<numdiv; i++) {
            func += nextsize;
            dfunc += i*nextsize/progression;
            nextsize *= progression;
        }
        func -= factor;
        //std::cout << "func = " << func << std::endl;
        progression -= func/dfunc;
        func = 0.;//pow(progression[idim],fNx[idim])-1.+factor[idim]*(1.-progression[idim]);
        nextsize = 1.;
        dfunc = 0.;
        for (int i=0; i<numdiv; i++) {
            func += nextsize;
            dfunc += i*nextsize/progression;
            nextsize *= progression;
        }
        func -= factor;
        
        iter++;
    }
    if (iter == maxiter)
    {
        DebugStop();
    }
    return progression;
	
}



/// compute the geometric progression such that the first elements have this size
void TPZGenGrid::ComputeGeometricProgression(TPZVec<REAL> &minsizes, TPZVec<REAL> &progression)
{
    int idim;
    progression.Resize(2, 1.);
    progression[0] = 1.;
    progression[1] = 1.;
    REAL factor[2];
    for (idim = 0; idim < 2; idim++) {
        // apply a simple newton method to find the geometric progression factors
        factor[idim] = (fX1[idim]-fX0[idim])/minsizes[idim];
        REAL func = 0.;//pow(progression[idim],fNx[idim])-1.+factor[idim]*(1.-progression[idim]);
        REAL nextsize = 1.;
        for (int i=0; i<fNx[idim]; i++) {
            func += nextsize;
            nextsize *= progression[idim];
        }
        func -= factor[idim];
        int iter = 0;
        int maxiter = 200;
        while (fabs(func/factor[idim]) >= 1.e-10 && iter < 200) {
            REAL dfunc = 0.;// fNx[idim]*pow(progression[idim], fNx[idim]-1)-factor[idim];
            func = 0.;
            nextsize = 1.;
            for (int i=0; i<fNx[idim]; i++) {
                func += nextsize;
                dfunc += i*nextsize/progression[idim];
                nextsize *= progression[idim];
            }
            func -= factor[idim];
            //std::cout << "func = " << func << std::endl;
            progression[idim] -= func/dfunc;
            func = 0.;//pow(progression[idim],fNx[idim])-1.+factor[idim]*(1.-progression[idim]);
            nextsize = 1.;
            for (int i=0; i<fNx[idim]; i++) {
                func += nextsize;
                dfunc += i*nextsize/progression[idim];
                nextsize *= progression[idim];
            }
            func -= factor[idim];
            
            iter++;
        }
        if (iter == maxiter)
        {
            DebugStop();
        }
    }
}

/// set the geometric progression of the mesh to be generated
void TPZGenGrid::SetGeometricProgression(TPZVec<REAL> &progression)
{
    for (int idim=0; idim<2; idim++) {
        REAL totalsize = 0;
        REAL nextsize = 1.;
        for (int i = 0; i<fNx[idim]; i++) {
            totalsize += nextsize;
            nextsize *= progression[idim];
        }
        fDelx[idim] = (fX1[idim]-fX0[idim])/totalsize;
    }
    fGeometricProgression = progression;
}

/// Generate a boundary geometric element at the indicated node
void TPZGenGrid::SetPointBC(TPZGeoMesh *gr, TPZVec<REAL> &x, int bc)
{
    // look for an element/corner node whose distance is close to start
    TPZGeoNode *gn1 = gr->FindNode(x);
    int iel;
    int nelem = gr->ElementVec().NElements();
    TPZGeoEl *gel;
    for (iel = 0; iel<nelem; iel++) {
        gel = gr->ElementVec()[iel];
        if(!gel) continue;
        int nc = gel->NCornerNodes();
        int c;
        for (c=0; c<nc; c++) {
            TPZGeoNode *gn = gel->NodePtr(c);
            if (gn == gn1) {
                break;
            }
        }
        if (c<nc) {
            TPZGeoElBC(gel, c, bc, *gr);
            return;
        }
    }
}

