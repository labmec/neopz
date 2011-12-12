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
// Obs.:   Gera uma malha sobre um dominio rectangular
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

#include "pzfmatrix.h"

#include "pzvec.h"
#include "pzstack.h"

#include <fstream>

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.gengrid.tpzgengrid"));
#endif

using namespace std;

TPZGenGrid::TPZGenGrid(TPZVec<int> &nx, TPZVec<REAL> &x0,TPZVec<REAL> &x1, int numl, REAL rot) : fNx(nx), fX0(x0), fX1(x1),
fDelx(2), fGeometricProgression(2,1.), fNumLayers(numl), fRotAngle(rot) {
	fDelx[0] = (x1[0]-x0[0])/(nx[0]);   // Delta x
	fDelx[1] = (x1[1]-x0[1])/(nx[1]);   // Delta y
	fNumNodes= (nx[0]+1)*(nx[1]+1)+(fNumLayers-1)*(nx[0])*(nx[1]+1);
	fElementType = 0;
}

TPZGenGrid::~TPZGenGrid() {    
}

short TPZGenGrid::Read(TPZGeoMesh &grid) {
	if(!GenerateNodes(grid))
		return 1;
    if(!GenerateElements(grid))
		return 1;
    return 0;
}
short TPZGenGrid::Read(TPZAutoPointer<TPZGeoMesh> grid) {
	if(!GenerateNodes(*(grid.operator->())))
		return 1;
    if(!GenerateElements(*(grid.operator->())))
		return 1;
    return 0;
}

/* This method merge gridtomerge into the gridinitial. The process is as follow:
 * 1. Each node in gridtomerge is checked wether it is in meshinitial. If it is then the node has its id modified as the id in gridinitial.
 * Elsewhere a new node is created into the gridinitial, then we get the id of the new node and change the id node (old) by the id of the new node in all the elements of the gridtomerge.
 * 2. For each node repeated in the gridtomerge (case 1. true), we get the id and is substitutived in the elements of the gridinitial.
 * 3. For each element in meshtomerge is constructed a new element copy in meshinitial.
 * 4. Constructed the connectivity again to meshinitial.
 */
bool TPZGenGrid::ReadAndMergeGeoMesh(TPZAutoPointer<TPZGeoMesh> gridinitial,TPZAutoPointer<TPZGeoMesh> tomerge) {
	// gridinitial is created by TPZGenGrid current
	if(Read(gridinitial))
		return false;
	if(!tomerge->NNodes() && !tomerge->NElements())
		return true;
	
	// Copy vectors for nodes and elements of the gridtomerge, then the original data is preserved
	TPZGeoMesh gtomerge(tomerge);
	TPZGeoMesh *gridtomerge = & gtomerge;
	
	int i,j,k,nnodestomerge = gridtomerge->NNodes();
	int nnodesinitial = gridinitial->NNodes();
	TPZVec<REAL> coordinitial(3,0.);
	TPZVec<REAL> coordtomerge(3,0.);
	TPZGeoNode *nodetomerge;
	TPZGeoEl *gel;
	// Verifing each node in gridtomerge if exist into the gridinitial (as same coordinates). It is inefficient.
	for(i=0;i<nnodestomerge;i++) {
		nodetomerge = &(gridtomerge->NodeVec()[i]);
		if(!nodetomerge) continue;
		nodetomerge->GetCoordinates(coordtomerge);
		for(j=0;j<nnodesinitial;j++) {
			gridinitial->NodeVec()[j].GetCoordinates(coordinitial);
			if(IsZero(Distance(coordtomerge,coordinitial))) {
				// In this case exists a node with same coordinates, then the id is update as id of the gridinitial with same coordinates
				// and the old id is stored in coord[0]
				nodetomerge->SetCoord(0,nodetomerge->Id());
				nodetomerge->SetNodeId(gridinitial->NodeVec()[j].Id());
				break;
			}
		}
		
		// If the node (i) not exists into the gridinitial is created a new node copy in this grid, and is substitutived in all the 
		// elements in gridinitial the id as old node with the id of the new node. At last the id of the node duplicated is put as -1
		if(j==nnodesinitial) {
			// resizing the vector of the nodes
			int index = gridinitial->NodeVec().AllocateNewElement();
			gridinitial->NodeVec()[index].Initialize(coordtomerge,gridinitial);
			index = gridinitial->NodeVec()[index].Id();
			int oldid = nodetomerge->Id();
			for(k=0;k<gridtomerge->NElements();k++) {
				gel = gridtomerge->ElementVec()[k];
				if(!gel) continue;
				for(int p=0;p<gel->NNodes();p++)
					if(gel->NodeIndex(p)==oldid)
						gel->SetNodeIndex(p,index);
			}
			nodetomerge->SetNodeId(-1);
		}
	}

	// changing the id of the repeated nodes into the geometric elements of the gridtomerge
	for(i=0;i<nnodestomerge;i++) {
		nodetomerge = &(gridtomerge->NodeVec()[i]);
		if(!nodetomerge || nodetomerge->Id()==-1) continue;
		int idnew = nodetomerge->Id(), idold = (int)(nodetomerge->Coord(0));
		for(k=0;k<gridtomerge->NElements();k++) {
			gel = gridtomerge->ElementVec()[k];
			if(!gel) continue;
			for(int p=0;p<gel->NNodes();p++)
				if(gel->NodeIndex(p)==idold)
					gel->SetNodeIndex(p,idnew);
		}
	}

	// creating new element into gridinitial corresponding for each element in gridtomerge
	for(i=0;i<gridtomerge->NElements();i++) {
		gel = gridtomerge->ElementVec()[i];
		if(!gel) continue;
		TPZVec<int> nos;
		int ngelnodes = gel->NNodes(), index;
		nos.Resize(gel->NNodes());
		for(j=0;j<ngelnodes;j++)
			nos[j] = gel->NodeIndex(j);
		if(!gridinitial->CreateGeoElement(gel->Type(),nos,gel->MaterialId(), index,0))
			DebugStop();
	}
	
	// computing the connectivity
	gridinitial->ResetConnectivities();
	gridinitial->BuildConnectivity();
	return true;
}

bool TPZGenGrid::GenerateNodes(TPZGeoMesh &grid) {
    // create the geometric nodes
	TPZVec<REAL> coor(3,0.);
	int i;
	// grid can not to contain other nodes and elements
	if(grid.NNodes()) {
#ifdef LOG4CXX
		LOGPZ_DEBUG(logger,"Mesh is not empty");
#endif
		return false;
	}

	// resizing the vector of the nodes
	grid.NodeVec().Resize(fNumNodes);
	for(i=0; i<fNumNodes; i++) {
		// computes the coordinates of the ith-node, depends on fElementType, layer and fRotAngle.
		Coord(i,coor);
		grid.NodeVec()[i].Initialize(coor,grid);
	}
	return true;
}

bool TPZGenGrid::GenerateElements(TPZGeoMesh &grid) {
	// create the geometric elements (retangular)    
	int num_rectangles=fNx[0]*fNx[1]*fNumLayers;
	TPZVec<int> nos(9);
	if(fElementType == 0) nos.Resize(4);
    int i, index;

	// grid can not to contain other elements
	if(grid.NElements()) {
#ifdef LOG4CXX
		LOGPZ_DEBUG(logger,"Mesh is not empty");
#endif
		return false;
	}
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
	return true;
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

void TPZGenGrid::Print( char *name , ostream &out  )
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
			TPZGeoElBC(gel,elementside,bc);
            //			gel->SetSide(side,bc);
		}
	}
}

void TPZGenGrid::SetBC(TPZGeoMesh *g, TPZVec<REAL> &start, TPZVec<REAL> &end, int bc) {
	TPZGeoNode *gn1 = g->FindNode(start);
	TPZGeoNode *gn2 = g->FindNode(end);
    
	TPZStack<TPZGeoEl *> ElementVec;
	TPZStack<int> Sides;
	g->GetBoundaryElements(gn1->Id(),gn2->Id(),ElementVec, Sides);
	int numel = ElementVec.NElements();
	for(int el=0; el<numel; el++) {
		TPZGeoEl *gel = (TPZGeoEl *) ElementVec[el];
		if(gel) {
			TPZGeoElBC(gel,Sides[el],bc);
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

REAL TPZGenGrid::Distance(TPZVec<REAL> &x1,TPZVec<REAL> &x2){
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
REAL TPZGenGrid::GeometricProgression(REAL minsize, REAL domainsize, int numdiv) {
    REAL progression = 1.;
    REAL factor = domainsize/minsize;
    REAL func = 0.;//pow(progression[idim],fNx[idim])-1.+factor[idim]*(1.-progression[idim]);
    REAL nextsize = 1.;
    for(int i=0; i<numdiv; i++) {
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
            TPZGeoElBC(gel, c, bc);
            return;
        }
    }
}

