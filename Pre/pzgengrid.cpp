/**
 * @file
 * @brief Contains the implementation of the TPZGenGrid methods. 
 */

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
	fElementType = EQuadrilateral;
}

TPZGenGrid::~TPZGenGrid() {    
}

void TPZGenGrid::SetData(TPZVec<int> &nx, TPZVec<REAL> &x0, TPZVec<REAL> &x1, MElementType eltype, int numl, REAL rot) {
	int i;
	fNx.Resize(nx.NElements());
	for(i=0;i<fNx.NElements();i++)
		fNx[i] = nx[i];
	fX0.Resize(x0.NElements());
	for(i=0;i<fX0.NElements();i++)
		fX0[i] = x0[i];
	fX1.Resize(x1.NElements());
	for(i=0;i<fX1.NElements();i++)
		fX1[i] = x1[i];
	fDelx.Resize(2);
	fDelx[0] = (x1[0]-x0[0])/(nx[0]);   // Delta x
	fDelx[1] = (x1[1]-x0[1])/(nx[1]);   // Delta y
	fNumNodes= (nx[0]+1)*(nx[1]+1)+(fNumLayers-1)*(nx[0])*(nx[1]+1);
	fElementType = eltype;
	fNumLayers = numl;
	fRotAngle = rot;
}

short TPZGenGrid::Read(TPZGeoMesh *grid,int matid) {
	if(!grid) return 1;
    grid->SetDimension(2);
	if(!GenerateNodes(grid))
		return 1;
    if(!GenerateElements(grid,matid))
		return 1;
	// computing the connectivity
	grid->ResetConnectivities();
	grid->BuildConnectivity();
    return 0;
}
short TPZGenGrid::Read(TPZAutoPointer<TPZGeoMesh> &grid) {
    if(!grid)
    {
        DebugStop();
    }
    grid->SetDimension(2);
	if(!GenerateNodes(grid.operator->()))
		return 1;
    if(!GenerateElements(grid.operator->()))
		return 1;
	// computing the connectivity
	grid->ResetConnectivities();
	grid->BuildConnectivity();
    return 0;
}

/* This method merge gridtomerge into the gridinitial. The process is as follow:
 * 1. Each node in gridtomerge is checked wether it is in meshinitial. If it is then the node has its id modified as the id in gridinitial.
 * Elsewhere a new node is created into the gridinitial, then we get the id of the new node and change the id node (old) by the id of the new node in all the elements of the gridtomerge.
 * 2. For each node repeated in the gridtomerge (case 1. true), we get the id and is substitutied in the elements of the gridinitial.
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
	
	long i,j,k,nnodestomerge = gridtomerge->NNodes();
	long nnodesinitial = gridinitial->NNodes();
//    int nneltomerge = gridtomerge->NElements();
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
			long index = gridinitial->NodeVec().AllocateNewElement();
			gridinitial->NodeVec()[index].Initialize(coordtomerge,gridinitial);
			index = gridinitial->NodeVec()[index].Id();
			long oldid = nodetomerge->Id();
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
		long idnew = nodetomerge->Id(), idold = (long)(nodetomerge->Coord(0));
		for(k=0;k<gridtomerge->NElements();k++) {
			gel = gridtomerge->ElementVec()[k];
			if(!gel) continue;
			for(int p=0;p<gel->NNodes();p++)
				if(gel->NodeIndex(p)==idold)
					gel->SetNodeIndex(p,idnew);
		}
	}
	
	// creating new element into gridinitial corresponding for each element in gridtomerge
    long nelmerge = gridtomerge->NElements();
	for(i=0; i< nelmerge; i++) {
		gel = gridtomerge->ElementVec()[i];
		if(!gel) continue;
		TPZVec<long> nos;
		long ngelnodes = gel->NNodes(), index;
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
bool TPZGenGrid::ReadAndMergeGeoMesh(TPZGeoMesh* gridinitial,TPZGeoMesh* tomerge,int matid) {
	// gridinitial is created by TPZGenGrid current
	if(Read(gridinitial,matid))
		return false;
	if(!tomerge->NNodes() && !tomerge->NElements())
		return true;
	
	// Copy vectors for nodes and elements of the gridtomerge, then the original data is preserved
	TPZGeoMesh gtomerge(*tomerge);
	TPZGeoMesh *gridtomerge = & gtomerge;
	
	long i,j,k,nnodestomerge = gridtomerge->NNodes();
	long nnodesinitial = gridinitial->NNodes();
	//    int nneltomerge = gridtomerge->NElements();
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
			long index = gridinitial->NodeVec().AllocateNewElement();
			gridinitial->NodeVec()[index].Initialize(coordtomerge,*gridinitial);
			index = gridinitial->NodeVec()[index].Id();
			long oldid = nodetomerge->Id();
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
		long idnew = nodetomerge->Id(), idold = (long)(nodetomerge->Coord(0));
		for(k=0;k<gridtomerge->NElements();k++) {
			gel = gridtomerge->ElementVec()[k];
			if(!gel) continue;
			for(int p=0;p<gel->NNodes();p++)
				if(gel->NodeIndex(p)==idold)
					gel->SetNodeIndex(p,idnew);
		}
	}
	
	// creating new element into gridinitial corresponding for each element in gridtomerge
    long nelmerge = gridtomerge->NElements();
	for(i=0; i< nelmerge; i++) {
		gel = gridtomerge->ElementVec()[i];
		if(!gel) continue;
		TPZVec<long> nos;
		long ngelnodes = gel->NNodes(), index;
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
bool TPZGenGrid::MergeGeoMesh(TPZGeoMesh* gridinitial,TPZGeoMesh* tomerge,int matid) {
	if(!tomerge->NNodes() && !tomerge->NElements())
		return true;
	
	// Copy vectors for nodes and elements of the gridtomerge, then the original data is preserved
	TPZGeoMesh gtomerge(*tomerge);
	TPZGeoMesh *gridtomerge = &gtomerge;

	long i,j,k,nnodestomerge = gridtomerge->NNodes();
	long nnodesinitial = gridinitial->NNodes();
	//    int nneltomerge = gridtomerge->NElements();
	TPZManVector<REAL> coordinitial(3,0.);
	TPZManVector<REAL> coordtomerge(3,0.);
	TPZGeoNode *nodetomerge;
	TPZGeoEl *gel;
	long newid = -1, oldid = -1;
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
				oldid = i;
//				nodetomerge->SetCoord(0,oldid);
				newid = j;
//				nodetomerge->SetNodeId(gridinitial->NodeVec()[j].Id());
				break;
			}
		}
		
		// If the node (i) not exists into the gridinitial is created a new node copy in this grid, and is substitutived in all the 
		// elements in gridinitial the id as old node with the id of the new node. At last the id of the node duplicated is put as -1
		if(j==nnodesinitial) {
			// resizing the vector of the nodes
			newid = gridinitial->NodeVec().AllocateNewElement();
			gridinitial->NodeVec()[newid].Initialize(coordtomerge,*gridinitial);
//			newid = gridinitial->NodeVec()[newid].Id();
			oldid = i;
//			nodetomerge->SetNodeId(newid);
		}
		for(k=0;k<gridtomerge->NElements();k++) {
			gel = tomerge->ElementVec()[k];
			if(!gel) continue;
			for(int p=0;p<gel->NNodes();p++) {
				if(gel->NodeIndex(p)==oldid) {
					gel = gridtomerge->ElementVec()[k];
					gel->SetNodeIndex(p,newid);
				}
			}
		}
	}
	
	// Must to check if geometric elements into gridtomerge are boundary elements and whether merged it is no boundary more
    long nelmerge = gridtomerge->NElements();
	for(i=0; i< nelmerge; i++) {
		TPZManVector<bool> founded(50);
		gel = gridtomerge->ElementVec()[i];
		if(!gel) continue;
		TPZManVector<long> nos;
		int ngelnodes = gel->NNodes();
		nos.Resize(gel->NNodes());
		for(j=0;j<ngelnodes;j++)
			nos[j] = gel->NodeIndex(j);
		long ii, jj, kk;
		// Comparating with all elements into gridinitial
		for(ii=0;ii<gridinitial->NElements();ii++) {
			founded.Fill(false);
			TPZGeoEl *gelinit = gridinitial->ElementVec()[ii];
			if(!gelinit) continue;
			for(kk=0;kk<nos.NElements();kk++) {
				for(jj=0;jj<gelinit->NNodes();jj++) {
					long nodegelinit = gelinit->NodeIndex(jj);
					if(nos[kk] == nodegelinit)
						break;
				}
				if(jj!=gelinit->NNodes()) {
					founded[kk] = true;
				}
			}
			for(kk=0;kk<nos.NElements();kk++) {
				if(!founded[kk]) {
					break;
				}
			}
			if(kk==nos.NElements()) {
				break;
			}
		}
		if(ii!=gridinitial->NElements()) {
			delete gel;
			gridtomerge->ElementVec()[i] = 0;
		}
	}
		
	// creating new element into gridinitial corresponding for each element in gridtomerge
    nelmerge = gridtomerge->NElements();
	for(i=0; i< nelmerge; i++) {
		gel = gridtomerge->ElementVec()[i];
		if(!gel) continue;
		TPZVec<long> nos;
		long ngelnodes = gel->NNodes(), index;
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

bool TPZGenGrid::GenerateNodes(TPZGeoMesh *grid) {
	if(!grid) return false;
    // create the geometric nodes
	TPZManVector<REAL> coor(3,0.);
	long i;
	// grid can not to contain other nodes and elements
	if(grid->NNodes()) {
#ifdef LOG4CXX
		LOGPZ_DEBUG(logger,"Mesh is not empty");
#endif
		return false;
	}
	
	// resizing the vector of the nodes
	grid->NodeVec().Resize(fNumNodes);
	for(i=0; i<fNumNodes; i++) {
		// computes the coordinates of the ith-node, depends on fElementType, layer and fRotAngle.
		Coord(i,coor);
		grid->NodeVec()[i].Initialize(coor,(*grid));
	}
	return true;
}

bool TPZGenGrid::GenerateElements(TPZGeoMesh *grid,int matid) {
	if(!grid) return false;
	// create the geometric elements (retangular)    
	int num_rectangles=fNx[0]*fNx[1]*fNumLayers;
	TPZVec<long> nos(9);
	if(fElementType == EQuadrilateral) nos.Resize(4);
    long i, index;
	
	// grid can not to contain other elements
	if(grid->NElements()) {
#ifdef LOG4CXX
		LOGPZ_DEBUG(logger,"Mesh is not empty");
#endif
		return false;
	}
	for(i=0; i<num_rectangles; i++) {
		ElementConnectivity(i,nos);
		if(fElementType == EQuadrilateral) {
		  grid->CreateGeoElement(EQuadrilateral,nos, matid, index,0);
		} else if(fElementType == ETriangle) {
            TPZManVector<long> nodloc(3);
            nodloc[0] = nos[0];
            nodloc[1] = nos[1];
            nodloc[2] = nos[2];
		  grid->CreateGeoElement(ETriangle,nodloc, matid, index,0);
            nodloc[0] = nos[0];
            nodloc[1] = nos[2];
            nodloc[2] = nos[3];
		  grid->CreateGeoElement(ETriangle,nodloc, matid, index,0);
		} else if(fElementType == ENoType) {
            std::cout << __PRETTY_FUNCTION__ << " - Quadratic interpolation is not available";
            DebugStop();
			grid->CreateGeoElement(EQuadrilateral,nos, matid, index,0);  
        }
	}
	grid->BuildConnectivity();
	return true;
}

void TPZGenGrid::Coord(int i, TPZVec<REAL> &coor) {
	int ix = 0;
	int iy = 0;
	int ilayer = 0;
	if(fElementType == EQuadrilateral || fElementType == ETriangle) {
		if(i < (fNx[0]+1)*(fNx[1]+1)) {
			ilayer = 0;
		} else {
			ilayer = (i-(fNx[0]+1)*(fNx[1]+1))/((fNx[0])*(fNx[1]+1))+1;
		}
	} else if(fElementType == ENoType) {
		if(i < (2*fNx[0]+1)*(2*fNx[1]+1) ) {
			ilayer = 0;
		} else {
			ilayer = (i-(2*fNx[0]+1)*(2*fNx[1]+1))/((2*fNx[0])*(2*fNx[1]+1))+1;
		}
	}
	REAL Rot = fRotAngle*ilayer;
	if(ilayer != 0 && (fElementType == EQuadrilateral || fElementType == ETriangle)) {
		i -= ((fNx[0]+1)*(fNx[1]+1))+(ilayer-1)*((fNx[0])*(fNx[1]+1));
	} else if (ilayer != 0 && fElementType == ENoType) {
		i -= ((2*fNx[0]+1)*(2*fNx[1]+1))+(ilayer-1)*((2*fNx[0])*(2*fNx[1]+1));
	}
	if(ilayer == 0) {
		if(fElementType == EQuadrilateral || fElementType == ETriangle) {
			ix = i%(fNx[0]+1);
			iy = i/(fNx[0]+1);
		} else if(fElementType == ENoType) {
			ix = i%(2*fNx[0]+1);
			iy = i/(2*fNx[0]+1);
		}
	} else {
		if(fElementType == EQuadrilateral || fElementType == ETriangle) {
			ix = i%(fNx[0])+1;
			iy = i/(fNx[0]);
		} else if(fElementType == ENoType) {
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
	coor[2] = fX0[2] + coorold[0]*sin(Rot);
    //	coor[0] = coorold[0];
	coor[1] = coorold[1];
    //   coor[2] = 0.;// 0.1*(coorold[0]-fX0[0])*(fX1[0]-coorold[0]);
}

void TPZGenGrid::ElementConnectivity(long i, TPZVec<long> &rectangle_nodes){
	int xel = i%(fNx[0]);
	int yel = (i/(fNx[0]))%(fNx[1]);
	int layer = i/(fNx[0]*fNx[1]);
    //cout << "ElConnectivity : xel = " << xel << " yel = " << yel << " layer = " << layer << endl;
    //cout.flush();
    if(fElementType == EQuadrilateral || fElementType == ETriangle) {
		rectangle_nodes[0] = GlobalI(xel,yel,layer);
		rectangle_nodes[1] = GlobalI(xel+1,yel,layer);
		rectangle_nodes[2] = GlobalI(xel+1,yel+1,layer);
		rectangle_nodes[3] = GlobalI(xel,yel+1,layer);
        //cout << "ElConnectivity : " << rectangle_nodes[0] << ' '<< rectangle_nodes[1] << ' '<<rectangle_nodes[2] << ' '<<rectangle_nodes[3] << endl;
        //cout.flush();
    } else if(fElementType == ENoType) {
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
	long ielfirst = 0;
	long iellast = 0;
	long ielinc;
	int layer;
	long iel;
	TPZGeoEl *gel;
	int elementside = side;
	for(layer=0; layer<fNumLayers; layer++) {
		switch(side) {
            case 4:
                ielfirst = layer*fNx[0]*fNx[1];
                iellast = ielfirst+fNx[0];
                ielinc = 1;
                break;
            case 5:
                ielfirst = layer*fNx[0]*fNx[1]+fNx[0]-1;
                iellast = (layer+1)*fNx[0]*fNx[1];
                ielinc = fNx[0];
                break;
            case 6:
                ielfirst = layer*fNx[0]*fNx[1]+fNx[0]*(fNx[1]-1);
                ielinc = 1;
                iellast = (layer+1)*fNx[0]*fNx[1];
                break;
            case 7:
                ielfirst = layer*fNx[0]*fNx[1];
                ielinc = fNx[0];
                iellast = ielfirst+fNx[1]*ielinc;
                break;
			default:
				cout << "It is not implemented for side = " << side << endl;
                DebugStop();
				return;
		}
		if(fElementType == ETriangle) {
			elementside -= 1;
			ielfirst *= 2;
			iellast *= 2;
			ielinc *= 2;
			if(side > 5) {
				ielfirst += 1;
				iellast += 1;
				side -= 1;
				elementside--;
			}
		}
		for(iel = ielfirst; iel<iellast; iel += ielinc) {
			gel = g->ElementVec()[iel];
			if(gel->HasSubElement()) continue;
			TPZGeoElBC(gel,elementside,bc);
		}
	}
}

void TPZGenGrid::SetBC(TPZGeoMesh *g, TPZVec<REAL> &start, TPZVec<REAL> &end, int bc) {
	TPZGeoNode *gn1 = g->FindNode(start);
	TPZGeoNode *gn2 = g->FindNode(end);
    
	TPZStack<TPZGeoEl *> ElementVec;
	TPZStack<int> Sides;
	g->GetBoundaryElements(gn1->Id(),gn2->Id(),ElementVec, Sides);
	long numel = ElementVec.NElements();
	for(long el=0; el<numel; el++) {
		TPZGeoEl *gel = (TPZGeoEl *) ElementVec[el];
		if(gel) {
			TPZGeoElBC(gel,Sides[el],bc);
        }
	}
}

int TPZGenGrid::GlobalI(int ix, int iy, int layer) {
	if(layer == 0 || ix == 0) {
		if(fElementType == EQuadrilateral || fElementType == ETriangle) {
			return ix+iy*(fNx[0]+1);
		} else if(fElementType == ENoType) {
			return ix+iy*(2*fNx[0]+1);
		}
	} else {
		if(fElementType == EQuadrilateral || fElementType == ETriangle) {
			return (fNx[0]+1)*(fNx[1]+1)+(fNx[0])*(fNx[1]+1)*(layer-1)+ix-1+iy*(fNx[0]);
		} else if(fElementType == ENoType) {
			return (2*fNx[0]+1)*(2*fNx[1]+1)+(2*fNx[0])*(2*fNx[1]+1)*(layer-1)+ix-1+iy*(2*fNx[0]);
		}
	}
	return 0;
}


long TPZGenGrid::ElemId(long iel,long jel, int layer){
    return (fElementType == EQuadrilateral || fElementType == ENoType)? ( iel*fNx[0]+jel+fNx[0]*fNx[1]*layer ):(iel*2*fNx[0]+jel+layer*2*fNx[0]*fNx[1]) ;
}

REAL TPZGenGrid::Distance(TPZVec<REAL> &x1,TPZVec<REAL> &x2){
	REAL l1,l2, l3;
	l1=x1[0]-x2[0];
	l2=x1[1]-x2[1];
	l3=x1[2]-x2[2];
	return( sqrt(l1*l1+l2*l2+l3*l3) );
}

void TPZGenGrid::SetElementType(MElementType type) {
    
    if (type != EQuadrilateral && type != ETriangle) {
        DebugStop();
    }
	fElementType = type;
	if(fElementType == EQuadrilateral || fElementType == ETriangle) {
		fNumNodes = (fNx[0]+1)*(fNx[1]+1)+(fNumLayers-1)*(fNx[0])*(fNx[1]+1);
		fDelx[0] = (fX1[0]-fX0[0])/(fNx[0]);
		fDelx[1] = (fX1[1]-fX0[1])/(fNx[1]);
	} else if(fElementType == ENoType) {
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
    long iel;
    long nelem = gr->ElementVec().NElements();
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

