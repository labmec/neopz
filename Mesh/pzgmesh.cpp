// -*- c++ -*-
/**File : pzgmesh.c

Method definition for class TPZGeoMesh.*/

#include "pzgmesh.h"
#include "pzvec.h"
//template class TPZVec<REAL>;
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzgnode.h"
#include "pzmaterial.h"
#include "pzerror.h"
#include "pzgeoel.h"
#include "pzcosys.h"
#include "pzmatrix.h"
#include "pzavlmap.h"


TPZGeoMesh::TPZGeoMesh() : fElementVec(0), fNodeVec(0), fCosysVec(0),
  fBCElementVec(0), fBCNodeVec(0) {

  fName[0] = '\0';
  fReference = 0;
  fNodeMaxId = -1;
  fElementMaxId = -1;
}

TPZGeoMesh::~TPZGeoMesh() {
  CleanUp();
}

/**Delete element, nodes, Cosys, boundary elements and boundary nodes in list*/
void TPZGeoMesh::CleanUp() {
  int i, nel = fElementVec.NElements();
  for(i=0; i<nel; i++) {
    TPZGeoEl *el = fElementVec[i];
    if(el) delete el;
    fElementVec[i] = 0;
  }
  fElementVec.Resize(0);
  fElementVec.CompactDataStructure(1);
  fNodeVec.Resize(0);
  fNodeVec.CompactDataStructure(1);
  fCosysVec.Resize(0);
  fCosysVec.CompactDataStructure(1);
  fBCElementVec.Resize(0);
  fBCElementVec.CompactDataStructure(1);
  fBCNodeVec.Resize(0);
  fBCNodeVec.CompactDataStructure(1);

}

void TPZGeoMesh::SetName (char *nm) {
  if(nm != NULL) {
    strncpy(fName,nm,62);
    fName[62] = '\0';
  }
}


void TPZGeoMesh::Print (ostream & out) {
  out << "\n\t\t GEOMETRIC TPZGeoMesh INFORMATIONS:\n\n";
  out << "TITLE-> " << fName << "\n\n";
  out << "number of nodes               = " << fNodeVec.NElements() << "\n";
  out << "number of elements            = " << fElementVec.NElements() << "\n";

  out << "\n\tGeometric Node Information:\n\n";
  int i;
  int nnodes = fNodeVec.NElements();
  for(i=0; i<nnodes; i++) {
    fNodeVec[i].Print(out);
    out << "\n";
  }
  out << "\n\tGeometric Element Information:\n\n";
  int nelem = fElementVec.NElements();
  for(i=0; i<nelem; i++) {
    if(fElementVec[i]) fElementVec[i]->Print(out);
    out << "\n";
  }

  out << "\nBoundary Element Information : \n\n";
  for(i=0; i<fBCElementVec.NElements () ; i++) {
	  fBCElementVec[i].Print (out);
  }
}

void TPZGeoMesh::GetNodePtr(TPZVec<int> &nos,TPZVec<TPZGeoNode *> &nodep) {

  int i,nnodes=nos.NElements();
  for(i=0;i<nnodes;i++) nodep[i]=&fNodeVec[nos[i]];
}

void  TPZGeoMesh::ResetReference() {

  TPZGeoEl *elp;
  int i,nelements=fElementVec.NElements();
  for(i=0;i<nelements;i++) {
    elp = fElementVec[i];
    if(elp) elp->ResetReference();
  }
  fReference = 0;
}

void TPZGeoMesh::RestoreReference(TPZCompMesh *cmesh) {

  ResetReference();
  fReference = cmesh;
  TPZGeoEl *gel;
  TPZCompEl *cel;
  int i,nelem = cmesh->ElementVec().NElements();
  for(i=0;i<nelem;i++) {
    cel = cmesh->ElementVec()[i];
    if(cel) {
      gel = cel->Reference();
      if(!gel) {
	PZError << "RestoreReference incomplete. Exist computational element with geometrical\n";
	PZError << "element not belongs to the current geometrical mesh.\n";
	return;
      }
      gel->SetReference(cel);
    }
  }
}

// GetBoundaryElements returns all elements beweeen NodFrom and NodTo counterclock wise
//		this method uses the connectivity of the elements
//		BuildConnectivity should be called to initialize the connectivity information
// 	this method will only work for grid with 2-D topology
//		the current version will only work for a grid with only one level
void TPZGeoMesh::GetBoundaryElements(int NodFrom, int NodTo,TPZStack<TPZGeoEl *> &ElementVec,TPZStack<int> &Sides) {
  // Find a first element whose first node on the side is NodFrom
  TPZGeoEl *def = 0;
  TPZAVLMap<int,TPZGeoEl *> elmap(def);
  int i,nelements=NElements();
  for(i=0;i<nelements;i++) {
    TPZGeoEl *el = fElementVec[i];
    if(el) elmap[el->Id()]=fElementVec[i];
  }

  int currentnode = NodFrom;
  TPZGeoEl *candidate = 0;
  int candidateside = 0;
  while(currentnode != NodTo) {
    // put all elements connected to currentnode in elmap, eliminate the elements
    //		from elmap which do not contain the node
    BuildElementsAroundNode(currentnode,elmap);
    // find, within elmap the element which has currentnode as its first boundary side
    //  	node
    FindElement(elmap, currentnode, candidate, candidateside);
    //	if the element found is already contained in the list, we have a circular list
    // if no element was found, the topology may not be two dimensional
    if(!candidate) break;
    int index = 0;
    int nelvec = ElementVec.NElements();
    while(index<nelvec && ElementVec[index] != candidate) index++;
    if(index <nelvec && Sides[index]==candidateside) break;
    ElementVec.Push(candidate);
    Sides.Push(candidateside);
    elmap.CleanUp();
    elmap[candidate->Id()] = candidate;
    // initialize the list in which to look for connected elements
    currentnode = candidate->SideNodeIndex(candidateside,1);
  }
}

// Find all elements in elmap or neighbour of elements in elmap which contain a node
void TPZGeoMesh::BuildElementsAroundNode(int currentnode,TPZAVLMap<int,TPZGeoEl*> &elmap){

  // first eliminate all elements which do not contain currentnode
  TPZPix iel = elmap.First();
  TPZGeoEl *el;
  int i;
  while(iel) {
    el = elmap.Contents(iel);
    elmap.Next(iel);
    int numnode = el->NNodes();
    for(i=0; i< numnode; i++) {
      if(el->NodeIndex(i) == currentnode) break;
    }
    if(i == numnode) elmap.Delete(el->Id());
  }
  iel = elmap.First();
  while(iel) {
    el = elmap.Contents(iel);
    elmap.Next(iel);
    int nside = el->NSides();
    for(int is=0; is<nside; is++) {
      TPZGeoElSide neigh = el->Neighbour(is);
      if(!neigh.Exists()) continue;
      int numnode = neigh.Element()->NNodes();
      for(i=0; i< numnode; i++) {
	if(neigh.Element()->NodeIndex(i) == currentnode){
	  if(!elmap.Contains(neigh.Element()->Id())) {
	    // this should be implemented as a stack, so that we dont have to
	    // 	go through the list again each time
	    elmap[neigh.Element()->Id()] = neigh.Element();
	    iel = elmap.First();
	  }
	  break;	// get out of the loop over the nodes
	}
      }
    }
  }
}

// find, within elmap the element which has currentnode as its first boundary side
//  	node
void TPZGeoMesh::FindElement(TPZAVLMap<int,TPZGeoEl *> &elmap,int currentnode,TPZGeoEl* &candidate,int &candidateside) {

  candidate = 0;
  TPZPix iel = elmap.First();
  while(iel) {
    TPZGeoEl *el = elmap.Contents(iel);
    elmap.Next(iel);
    int ns = el->NSides();
    int is = el->NCornerNodes();
    for(; is < ns; is++) {
      TPZGeoElSide neigh = el->Neighbour(is);
      TPZGeoElSide father = el->Father2(is);
      if(!neigh.Exists() && !father.Exists() && el->SideNodeIndex(is,0) == currentnode) {
	candidate = el;
	candidateside = is;
	return;
      }
    }
  }
}

TPZGeoNode *TPZGeoMesh::FindNode(TPZVec<REAL> &co) {

  int i=0, in, nnodes = fNodeVec.NElements();
  while(i<nnodes && fNodeVec[i].Id() == -1) i++;
  if(i == nnodes) return 0;
  TPZGeoNode *gnkeep = &fNodeVec[i];
  REAL distkeep = 0.;
  for(in=0;in<3;in++)
    distkeep += (co[in]-(gnkeep->Coord(in)))*(co[in]-(gnkeep->Coord(in)));
  while(i< nnodes) {
    TPZGeoNode *gn = &fNodeVec[i];
    REAL dist = 0.;
    for(in=0;in<3;in++)
      dist += (co[in]-gn->Coord(in))*(co[in]-gn->Coord(in));
    if(dist < distkeep) {
      gnkeep = gn;
      distkeep = dist;
    }
    i++;
    while(i<nnodes && fNodeVec[i].Id() == -1) i++;
  }
  return gnkeep;
}
void TPZGeoMesh::BuildConnectivity() {

  TPZVec<int> SideNum(NNodes(),-1);
  TPZVec<TPZGeoEl *> NeighNode(NNodes(),0);
  int nelem = NElements();
  int iel = 0;
  while(iel<nelem && fElementVec[iel] == 0) iel++;

  long numsearch =1;
  // if there are no elements, do nothing
  while(iel < nelem) {
    TPZGeoEl *el = fElementVec[iel];
    int numsides = el->NSides();
    int side;
    for(side = 0;side<numsides;side++) {

      // check whether all entries in NeighNode are equal

	int equalnode = 1;
    int numsidenodes = el->NSideNodes(side);
    int sidenode = el->SideNodeIndex(side,0);
    TPZGeoEl *neigh = NeighNode[sidenode];
    int sidenumber = SideNum[sidenode];
    for(int sn = 0;sn < numsidenodes; sn++) {
		sidenode = el->SideNodeIndex(side,sn);
		if (neigh != NeighNode[sidenode]){
			equalnode=0;
			break;
		}
	}

    if(equalnode && neigh == 0) {
		if(el->SideIsUndefined(side)) {
			int elloaded = 0;
			for(int in=0; in<el->NNodes(); in++) {
				if(NeighNode[el->NodeIndex(in)] == el) elloaded = 1;
			}
			// this element is not loaded and its side is undefined

			// load the element side in the NeighNode vector
			for(int sn=0;!elloaded && sn < numsidenodes; sn++) {
				sidenode = el->SideNodeIndex(side,sn);
				NeighNode[sidenode] = el;
				SideNum[sidenode] = side;
			}
			numsearch++;
		}
	} else if(equalnode && side == sidenumber && neigh == el) {
	// unload the element side
		for(int sn=0;sn < numsidenodes; sn++) {
			sidenode = el->SideNodeIndex(side,sn);
			NeighNode[sidenode] = 0;
			SideNum[sidenode] = -1;
		}
	// if no neighbouring element was detected during the loop
	//    define the element side as undefined
	TPZGeoElSide neighbour = el->Neighbour(side);
	if(!neighbour.Exists()) el->SetSideDefined(side);
	numsearch++;
      } else if(equalnode && neigh != el) {
	// we found a neigbour
	TPZManVector<int> SideNodes(numsidenodes);
	// detect which side of the neigbour is loaded witin NeighNode
	for(int sn=0;sn < numsidenodes; sn++) {
		sidenode = el->SideNodeIndex(side,sn);
		SideNodes[sn] = sidenode;
	}
	int neighside = neigh->WhichSide(SideNodes);
	TPZGeoElSide neighbour(neigh,neighside);
	// WhichSide will tell the side number which contains the vector
	//    of node ids SideNodes
	if(neighbour.Side() != -1 && !el->NeighbourExists(side,neighbour)){
	  TPZGeoElSide(el,side).SetConnectivity(neighbour);
	  numsearch++;
	}
      }
    } // loop over the sides
    iel++;
    while(iel<nelem && fElementVec[iel] == 0) iel++;
    if(iel==nelem && numsearch) {
      numsearch = 0;
      iel = 0;
      while(iel<nelem && fElementVec[iel] == 0) iel++;
    }
  }
}

//Cedric : 03/03/99
TPZGeoEl *TPZGeoMesh::FindElement(int elid) {

	int nel = fElementVec.NElements();
   TPZGeoEl *gel = 0;
	for(int i=0;i<nel;i++) {
    	gel = fElementVec[i];
      if(gel && gel->Id()==elid) break;
   }
   return gel;
}

int TPZGeoMesh::ElementIndex(TPZGeoEl *gel){
	int i=0;
	int numel = ElementVec().NElements();
	while ( i < numel ) {
		if (ElementVec()[i] == gel) break;
		i++;
	}
	return i;
}

int TPZGeoMesh::NodeIndex(TPZGeoNode *nod){
	int i=0;
	int numel = NodeVec().NElements();
	while ( i < numel ) {
		if (&NodeVec()[i] == nod) break;
		i++;
	}
	return i;
}
