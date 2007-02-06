//$Id: pzgmesh.cpp,v 1.33 2007-02-06 17:43:11 cesar Exp $

// -*- c++ -*-
/**File : pzgmesh.c

Method definition for class TPZGeoMesh.*/

#include "pzgmesh.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "pzvec.h"
//template class TPZVec<REAL>;
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzgnode.h"
#include "pzmaterial.h"
#include "pzerror.h"
#include "pzgeoel.h"
//#include "pzcosys.h"
#include "pzmatrix.h"
//#include "pzavlmap.h"

#include <TPZRefPattern.h>
#include <tpzgeoelrefpattern.h>

#include <sstream>

#include <string>
#ifdef BORLAND
#include <io.h>
#include <fcntl.h>
#endif


using namespace std;

TPZGeoMesh::TPZGeoMesh() : fElementVec(0), fNodeVec(0){

//  fName[0] = '\0';
  fReference = 0;
  fNodeMaxId = -1;
  fElementMaxId = -1;
//  InitializeRefPatterns();
}

TPZGeoMesh::TPZGeoMesh(const TPZGeoMesh &cp):
              TPZSaveable(cp){
  this->operator =(cp);
}

TPZGeoMesh & TPZGeoMesh::operator= (const TPZGeoMesh &cp ){
  this->CleanUp();

  this->fName = cp.fName;
  int i, n = cp.fNodeVec.NElements();
  this->fNodeVec.Resize( n );
  for(i = 0; i < n; i++){
    this->fNodeVec[i] = cp.fNodeVec[i];
  }//for
  n = cp.fElementVec.NElements();
  this->fElementVec.Resize( n );
  for(i = 0; i < n; i++){
    this->fElementVec[i] = cp.fElementVec[i]->Clone(*this);
  }

  this->fNodeMaxId = cp.fNodeMaxId;
  this->fElementMaxId = cp.fElementMaxId;
  this->fInterfaceMaterials = cp.fInterfaceMaterials;
  this->fRefPatterns = cp.fRefPatterns;
//  this->fCosysVec = cp.fCosysVec;
  this->fReference = NULL;
  return *this;
}

TPZGeoMesh::~TPZGeoMesh() {
  CleanUp();
}

/**Delete element, nodes, Cosys, boundary elements and boundary nodes in list*/
void TPZGeoMesh::CleanUp() {
  int i, nel = fElementVec.NElements();
  for(i=0; i<nel; i++) {
    TPZGeoEl *el = fElementVec[i];
    if(el) {
      delete el;
      fElementVec[i] = 0;
    }
  }
  fElementVec.Resize(0);
  fElementVec.CompactDataStructure(1);
  fNodeVec.Resize(0);
  fNodeVec.CompactDataStructure(1);
//   fCosysVec.Resize(0);
//   fCosysVec.CompactDataStructure(1);

  fRefPatterns.clear();
}

void TPZGeoMesh::SetName (char *nm) {
  fName = nm;
//  if(nm != NULL) {
//    strncpy(fName,nm,62);
//    fName[62] = '\0';
//  }
}


// void TPZGeoMesh::PatternSidesFile(std::ofstream &filename){
//
//   int count=0;
//   filename << std::endl;
//   filename << std::endl;
//   std::map< MElementType,list<TPZRefPattern *> >::iterator first = fRefPatterns.begin();
//   std::map< MElementType,list<TPZRefPattern *> >::iterator last = fRefPatterns.end();
//   std::map< MElementType,list<TPZRefPattern *> >::iterator iter;
//   for (iter = first; iter != last; iter++){
//     list<TPZRefPattern *> &map_el = (*iter).second;
//     list<TPZRefPattern *>::iterator name_first = map_el.begin();
//     list<TPZRefPattern *>::iterator name_last = map_el.end();
//     list<TPZRefPattern *>::iterator name_iter;
//     for(name_iter = name_first; name_iter != name_last; name_iter++){
//       TPZRefPattern *refpat = (*name_iter);
//       TPZGeoEl *elemento = refpat->Element(0);
//       int nsides = elemento->NSides();
//       TPZVec<int> indices;
//       //TPZVec<int> selected(nsides,0);
//       filename << "Novo Padrao" << std::endl;
//       filename << "Tipo do elemento " << elemento->Type() << std::endl;
//       filename << "Lados refinados " ;
//       for (int p=0 ; p<nsides; p++){
//         if (elemento->SideDimension(p)==1){
//           int referencia = refpat->NSideNodes(p);
//           //elemento->MidSideNodeIndices(p,indices);
//           if (referencia!=0){
//              //selected[p]=1;
//             filename << p << " ";
//           }
//         }
//       }
//       //std::ofstream teste("refpatterns_geral.txt");
//       //refpat->CreateFile(filename);
//       count++;
//       //delete refpat;
//       filename << std::endl ;
//     }
//   }
//   filename.seekp (0);
//   filename << count ;
//   /*filename << std::endl;*/
//
// }

void TPZGeoMesh::PatternSidesFile(std::ofstream &filename){

  int count=0;
  filename << std::endl;
  filename << std::endl;
  std::map< MElementType, std::map<int, TPZAutoPointer<TPZRefPattern> > >::iterator first = fRefPatterns.begin(),
                                                                     last = fRefPatterns.end(),
                                                                     iter;
  for (iter = first; iter != last; iter++){
    std::map<int, TPZAutoPointer<TPZRefPattern> > &map_el = (*iter).second;
    std::map<int, TPZAutoPointer<TPZRefPattern> >::iterator name_first = map_el.begin();
    std::map<int, TPZAutoPointer<TPZRefPattern> >::iterator name_last = map_el.end();
    std::map<int, TPZAutoPointer<TPZRefPattern> >::iterator name_iter;
    for(name_iter = name_first; name_iter != name_last; name_iter++){
      TPZAutoPointer<TPZRefPattern> refpat = name_iter->second;
      refpat->ShortPrint(filename);
      count++;
      filename << std::endl ;
    }
  }
  //filename.seekp (0);
  //filename << count ;
  /*filename << std::endl;*/

}

void TPZGeoMesh::PatternFileLoad(std::ifstream &file){
  int i,k;
  file >> k ;
  std::vector< TPZAutoPointer<TPZRefPattern> > collect(k);
  for(i=0; i<k; i++)
  {
    collect[i] = TPZAutoPointer<TPZRefPattern>(new TPZRefPattern(this));
  }
  for(i=0; i<k; i++){
//    std::cout << i+1 << endl;
    collect[i]->ReadPattern(file,collect);
    this->InsertRefPattern(collect[i]);
    //A->InsertPermuted(*this);
  }
}
int TPZGeoMesh::NRefPatterns (){
  //return fRefPatterns.size();
  int count = 0;
  std::map< MElementType, std::map<int, TPZAutoPointer<TPZRefPattern> > >::iterator it;
  for (it=fRefPatterns.begin();it!=fRefPatterns.end();it++)
  {
    count += (*it).second.size();
  }
  return count;
}


void TPZGeoMesh::RefPatternFile(std::ofstream &filename){

  int count=0;
  filename << NRefPatterns() << std::endl ;
  std::map< MElementType, std::map<int, TPZAutoPointer<TPZRefPattern> > >::iterator first = fRefPatterns.begin(),
                                                                     last = fRefPatterns.end(),
                                                                     iter;
  for (iter = first; iter != last; iter++){
    std::map<int, TPZAutoPointer<TPZRefPattern> > &map_el = (*iter).second;
    std::map<int, TPZAutoPointer<TPZRefPattern> >::iterator name_first = map_el.begin();
    std::map<int, TPZAutoPointer<TPZRefPattern> >::iterator name_last = map_el.end();
    std::map<int, TPZAutoPointer<TPZRefPattern> >::iterator name_iter;
    for(name_iter = name_first; name_iter != name_last; name_iter++){
      TPZAutoPointer<TPZRefPattern> refpat = name_iter->second;
      //std::ofstream teste("refpatterns_geral.txt");
      refpat->WritePattern(filename);
      count++;
      //delete refpat;
    }
  }
//  filename.seekp (0);
//  filename << count ;
  /*filename << std::endl*/;
}

void TPZGeoMesh::Print (ostream & out) {
  out << "\n\t\t GEOMETRIC TPZGeoMesh INFORMATIONS:\n\n";
  out << "TITLE-> " << fName << "\n\n";
  out << "number of nodes               = " << fNodeVec.NElements() << "\n";
  out << "number of elements            = " << fElementVec.NElements()-fElementVec.NFreeElements() << "\n";

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

  out << "\nInterface materials : \n\n";
  InterfaceMaterialsMap::iterator w, e = this->fInterfaceMaterials.end();
  const int n = this->fInterfaceMaterials.size();
  int l, r, m;
  out << "number = " << n << "\n";
  out << "left material / right material -> interface material\n";
  for(w = this->fInterfaceMaterials.begin(); w != e; w++){
    l = w->first.first;
    r = w->first.second;
    m = w->second;
    out << l << " / " << r << " -> " << m << "\n";
  }

  out << "\nPrinting refinement patterns:\n";
  std::map<MElementType,std::map< int, TPZAutoPointer<TPZRefPattern> > >::const_iterator itg, eg;
  eg = this->fRefPatterns.end();
  for(itg = this->fRefPatterns.begin(); itg != eg; itg++){
    const std::map<int, TPZAutoPointer<TPZRefPattern> > & mymap = itg->second;
    out << "Element type = " << itg->first << std::endl;
    out << "Number of refinement patterns for this element type: " << mymap.size() << std::endl;
    std::map<int, TPZAutoPointer<TPZRefPattern> >::const_iterator it, e;
    e = mymap.end();
    for(it = mymap.begin(); it != e; it++){
      it->second->/*ShortPrint*/Print1(out);
      out << "\n";
    }//for it
    out << "\n\n";
  }//for itg

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
//  TPZGeoEl *def = 0;
  //TPZAVLMap<int,TPZGeoEl *> elmap(def);
  map<int,TPZGeoEl *> elmap;
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
    elmap.erase(elmap.begin(), elmap.end());//CleanUp();
    elmap[candidate->Id()] = candidate;
    // initialize the list in which to look for connected elements
    currentnode = candidate->SideNodeIndex(candidateside,1);
  }
}

// Find all elements in elmap or neighbour of elements in elmap which contain a node
//void TPZGeoMesh::BuildElementsAroundNode(int currentnode,TPZAVLMap<int,TPZGeoEl*> &elmap){
void TPZGeoMesh::BuildElementsAroundNode(int currentnode,map<int,TPZGeoEl*> &elmap){
  // first eliminate all elements which do not contain currentnode
  //TPZPix iel = elmap.First();
	map<int, TPZGeoEl *>::iterator ielprev,iel=elmap.begin();
  TPZGeoEl *el;
  int i;
  while(iel!=elmap.end()) {
    el = iel->second;
	ielprev=iel;
    iel++;//elmap.Next(iel);
    int numnode = el->NNodes();
    for(i=0; i< numnode; i++) {
      if(el->NodeIndex(i) == currentnode) break;
    }
    if(i == numnode) {
		elmap.erase(ielprev);
	}
  }
  iel = elmap.begin();
  while(iel!=elmap.end()) {
    el = iel->second;//elmap.Contents(iel);
    iel++;//elmap.Next(iel);
    int nside = el->NSides();
    for(int is=0; is<nside; is++) {
      TPZGeoElSide neigh = el->Neighbour(is);
      if(!neigh.Exists()) continue;
      int numnode = neigh.Element()->NNodes();
      for(i=0; i< numnode; i++) {
	if(neigh.Element()->NodeIndex(i) == currentnode){
	  if(elmap.find(neigh.Element()->Id())!=elmap.end()) {
	    // this should be implemented as a stack, so that we dont have to
	    // 	go through the list again each time
	    elmap[neigh.Element()->Id()] = neigh.Element();
	    iel = elmap.begin();
	  }
	  break;	// get out of the loop over the nodes
	}
      }
    }
  }
}

// find, within elmap the element which has currentnode as its first boundary side
//  	node
//void TPZGeoMesh::FindElement(TPZAVLMap<int,TPZGeoEl *> &elmap,int currentnode,TPZGeoEl* &candidate,int &candidateside) {
void TPZGeoMesh::FindElement(map<int,TPZGeoEl *> &elmap,int currentnode,TPZGeoEl* &candidate,int &candidateside) {

  candidate = 0;
  //TPZPix iel = elmap.First();
  map<int , TPZGeoEl *>::iterator ielprev, iel = elmap.begin();
  while(iel!=elmap.end()) {
    TPZGeoEl *el = iel->second;//elmap.Contents(iel);
	ielprev=iel;
    iel++;//elmap.Next(iel);
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

void TPZGeoMesh::BuildConnectivity()
{
  TPZVec<int> SideNum(NNodes(),-1);
  TPZVec<TPZGeoEl *> NeighNode(NNodes(),0);
  int nelem = NElements();
  int iel = 0;
  for(iel=0; iel<nelem; iel++)
    {
      TPZGeoEl *gel = fElementVec[iel];
      if(!gel) continue;
      int ncor = gel->NCornerNodes();
      int in;
      for(in=0; in<ncor; in++) {
	int nod = gel->NodeIndex(in);
	if(SideNum[nod] == -1)
	  {
	    NeighNode[nod] = gel;
	    SideNum[nod] = in;
	    if(gel->SideIsUndefined(in)) gel->SetSideDefined(in);
	  } else
	    {
	      TPZGeoElSide neigh(NeighNode[nod],SideNum[nod]);
	      TPZGeoElSide gelside(gel,in);
              if(!neigh.NeighbourExists(gelside))
              {
                neigh.SetConnectivity(gelside);
              }
	    }
      }
    }
  for(iel=0; iel<nelem; iel++)
    {
      TPZGeoEl *gel = fElementVec[iel];
      if(!gel) continue;
      int ncor = gel->NCornerNodes();
      int nsides = gel->NSides();
      int is;
      for(is=ncor; is<nsides; is++)
	{
	  if( gel->SideIsUndefined(is))
	    {
	      gel->SetSideDefined(is);
	      TPZGeoElSide gelside(gel,is);
	      TPZStack<TPZGeoElSide> neighbours;
	      gelside.ComputeNeighbours(neighbours);
	      int nneigh = neighbours.NElements();
	      int in;
	      for(in=0; in<nneigh; in++) {
		if(neighbours[in].Side() == -1)
		  {
		    std::cout << "TPZGeoMesh::BuildConnectivity : Inconsistent mesh detected!\n";
		    continue;
		  }
		gelside.SetConnectivity(neighbours[in]);
	      }
	    }
	}
    }
}

void TPZGeoMesh::BuildConnectivity2() {

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
        int index = gel->Index();
        if (ElementVec()[index] == gel) return index;

	int numel = ElementVec().NElements();
	while ( i < numel ) {
		if (ElementVec()[i] == gel) break;
		i++;
	}
	if(i<numel) return i;
	else return -1;
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

#include "TPZGeoElement.h"
#include "TPZGeoCube.h"
#include "pzshapecube.h"
#include "TPZRefCube.h"
#include "pzshapelinear.h"
#include "TPZGeoLinear.h"
#include "TPZRefLinear.h"
#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "pzshapetriang.h"
#include "pzreftriangle.h"
#include "pzgeotriangle.h"
#include "pzshapeprism.h"
#include "pzrefprism.h"
#include "pzgeoprism.h"
#include "pzshapetetra.h"
#include "pzreftetrahedra.h"
#include "pzgeotetrahedra.h"
#include "pzshapepiram.h"
#include "pzrefpyram.h"
#include "pzgeopyramid.h"
#include "pzgeopoint.h"
#include "pzrefpoint.h"
#include "pzshapepoint.h"

using namespace pzgeom;
using namespace pzrefine;
using namespace pzshape;

TPZGeoEl *TPZGeoMesh::CreateGeoElement(MElementType type,
                                       TPZVec<int>& nodeindexes,
                                       int matid,
                                       int& index,
                                       int reftype){
  if (!reftype)  switch( type ){
    case 0://point
      return new TPZGeoElement< TPZShapePoint, TPZGeoPoint, TPZRefPoint>(
                              nodeindexes, matid, *this, index );
    case 1://line
      return new TPZGeoElement< TPZShapeLinear, TPZGeoLinear, TPZRefLinear>(
                              nodeindexes, matid, *this, index );
    case 2://triangle
      return new TPZGeoElement< TPZShapeTriang, TPZGeoTriangle, TPZRefTriangle >(
                              nodeindexes, matid, *this, index );
      //return new TPZGeoElT2d(nodeindexes,matid,*this);
    case 3://quadrilatera
      return  new TPZGeoElement< TPZShapeQuad, TPZGeoQuad, TPZRefQuad >(
                              nodeindexes, matid, *this, index );
      //return new TPZGeoElQ2d(nodeindexes,matid,*this);
    case 4://tetraedra
      //return new TPZGeoElT3d(nodeindexes,matid,*this);
      return new TPZGeoElement< TPZShapeTetra, TPZGeoTetrahedra, TPZRefTetrahedra >(
                              nodeindexes, matid, *this, index );
    case 5:
      //return new TPZGeoElPi3d(nodeindexes,matid,*this);
      return new TPZGeoElement< TPZShapePiram, TPZGeoPyramid, TPZRefPyramid >(
                              nodeindexes, matid, *this, index );
    case 6:
      return new TPZGeoElement< TPZShapePrism, TPZGeoPrism, TPZRefPrism >(
                              nodeindexes, matid, *this, index );
    case 7:
      return new TPZGeoElement< TPZShapeCube, TPZGeoCube, TPZRefCube >(
                              nodeindexes, matid, *this, index );
    default:
      PZError << "TPZGeoMesh::CreateGeoElement type element not exists:"
              << " type = " << type << std::endl;
    return NULL;
  } else {
    //TPZAutoPointer<TPZRefPattern> ref = GetUniformPattern(type);
    switch( type ){
      case 0://point
      {
        TPZGeoElRefPattern<TPZShapePoint, TPZGeoPoint> * gel =
            new TPZGeoElRefPattern<TPZShapePoint, TPZGeoPoint> (nodeindexes, matid, *this, index);
        return gel;
      }
      case 1://line
      {
        TPZGeoElRefPattern < TPZShapeLinear, TPZGeoLinear > *gel =
            new TPZGeoElRefPattern < TPZShapeLinear, TPZGeoLinear >
                (nodeindexes, matid, *this, index);
        //gel->SetRefPattern (ref);
        return gel;
      }
      case 2://triangle
      {
        TPZGeoElRefPattern < TPZShapeTriang, TPZGeoTriangle > *gel =
            new TPZGeoElRefPattern < TPZShapeTriang, TPZGeoTriangle >
                (nodeindexes, matid, *this, index);
        //gel->SetRefPattern (ref);
        return gel;
      }
      case 3://quadrilatera
      {
        TPZGeoElRefPattern < TPZShapeQuad, TPZGeoQuad > * gel =
            new TPZGeoElRefPattern < TPZShapeQuad, TPZGeoQuad >
                (nodeindexes, matid, *this, index);
        //gel->SetRefPattern (ref);
        return gel;
      }
      case 4://tetraedra
      {
        TPZGeoElRefPattern < TPZShapeTetra, TPZGeoTetrahedra > *gel =
            new TPZGeoElRefPattern < TPZShapeTetra, TPZGeoTetrahedra >
                (nodeindexes, matid, *this, index);
        //gel->SetRefPattern (ref);
        return gel;
      }
      case 5://pyramid
      {
        TPZGeoElRefPattern < TPZShapePiram, TPZGeoPyramid > *gel =
            new TPZGeoElRefPattern < TPZShapePiram, TPZGeoPyramid >
                (nodeindexes, matid, *this, index);
        //gel->SetRefPattern (ref);
        return gel;
      }
      case 6://prism
      {
        TPZGeoElRefPattern < TPZShapePrism, TPZGeoPrism > *gel =
            new TPZGeoElRefPattern < TPZShapePrism, TPZGeoPrism >
                (nodeindexes, matid, *this, index);
        //gel->SetRefPattern (ref);
        return gel;
      }
      case 7://cube
      {
        TPZGeoElRefPattern < TPZShapeCube, TPZGeoCube > *gel =
            new TPZGeoElRefPattern < TPZShapeCube, TPZGeoCube >
                (nodeindexes, matid, *this, index);
        //gel->SetRefPattern (ref);
        return gel;
      }
      default:
      {
        PZError << "TPZGeoMesh::CreateGeoElement type element not exists:"
                << " type = " << type << std::endl;
        return NULL;
      }
    }
  }
  //return NULL;
}

TPZAutoPointer<TPZRefPattern> TPZGeoMesh::GetUniformPattern(MElementType &type)
{
  if (!fRefPatterns[type].size())
  {
    TPZAutoPointer < TPZRefPattern > NULLrefpat;
    return NULLrefpat;
  }
  std::map<int, TPZAutoPointer<TPZRefPattern> >::iterator refIt;
  refIt = fRefPatterns[type].begin();
  return refIt->second;
}

/** check whether the refinement pattern already exists */
TPZAutoPointer<TPZRefPattern> TPZGeoMesh::FindRefPattern(TPZAutoPointer<TPZRefPattern> &refpat)
{
  TPZAutoPointer<TPZRefPattern> NULLRefPat;
  if(!(refpat) ) return NULLRefPat;
  //cout << "Looking for "; refpat->ShortPrint(cout); cout << endl;
  MElementType eltype = refpat->Element(0)->Type();
  std::map<int, TPZAutoPointer<TPZRefPattern> >::iterator it;
  for(it=fRefPatterns[eltype].begin(); it != fRefPatterns[eltype].end(); it++)
  {
    //std::cout << "Comparing with "; it->second.operator->()->ShortPrint(cout); cout << endl;
   // std::cout << it->second.operator->() << std::endl;
   // std::cout << refpat.operator->() << std::endl;
    if( *(it->second.operator->()) == refpat ) return it->second;
  }
  //cout << "Not found\n";
  return NULLRefPat;
}


void TPZGeoMesh::InsertRefPattern(TPZAutoPointer<TPZRefPattern> &refpat){
  if (!refpat) {
    PZError << "TPZGeoMesh::InsertRefPattern ERROR : NULL refinement pattern! " << std::endl;
    return;
  }
  MElementType eltype = refpat->Element(0)->Type();
  const int id = fRefPatterns[eltype].size();
  refpat->SetId(id);
  fRefPatterns[eltype][id] = refpat;
}

TPZAutoPointer<TPZRefPattern> TPZGeoMesh::GetRefPattern(MElementType eltype, const std::string &name)
{
  TPZAutoPointer < TPZRefPattern > NULLrefpat;
  std::map< MElementType, std::map<int, TPZAutoPointer<TPZRefPattern> > >::iterator eltype_iter = fRefPatterns.find(eltype);
  if (eltype_iter == fRefPatterns.end()) return NULLrefpat;
  std::map<int, TPZAutoPointer<TPZRefPattern> >::iterator name_iter = fRefPatterns[eltype].begin();
  NULLrefpat = name_iter->second;
  while(name_iter != fRefPatterns[eltype].end())
  {
    if(name_iter->second->GetName() == name) return name_iter->second;
    name_iter++;
  }
  return NULLrefpat;
}

/*
TPZGeoEl* TPZGeoMesh::CreateGeoElement( MElementType type, int* nodeindexes,
					int matid, int& index )
{
   switch( type )
   {
      case 0://point
	 return new TPZGeoElement< TPZShapeLinear, TPZGeoPoint, TPZRefPoint >(
	    nodeindexes, matid, *this, index );

      case 1://line
	 return new TPZGeoElement< TPZShapeLinear, TPZGeoLinear, TPZRefLinear >(
	    nodeindexes, matid, *this, index );

      case 2://triangle
	 return new TPZGeoElement<
	    TPZShapeTriang, TPZGeoTriangle, TPZRefTriangle >(
	       nodeindexes, matid, *this, index );

      case 3://quadrilatera
	 return  new TPZGeoElement< TPZShapeQuad, TPZGeoQuad, TPZRefQuad >(
	    nodeindexes, matid, *this, index );

      case 4://tetraedra
	 return new TPZGeoElement<
	    TPZShapeTetra, TPZGeoTetrahedra, TPZRefTetrahedra >(
	       nodeindexes, matid, *this, index );

      case 5:
	 return new TPZGeoElement<
	    TPZShapePiram, TPZGeoPyramid, TPZRefPyramid >(
	       nodeindexes, matid, *this, index );

      case 6:
	 return new TPZGeoElement< TPZShapePrism, TPZGeoPrism, TPZRefPrism >(
	    nodeindexes, matid, *this, index );

      case 7:
	 return new TPZGeoElement< TPZShapeCube, TPZGeoCube, TPZRefCube >(
	    nodeindexes, matid, *this, index );

      default:
	 PZError << "TPZGeoMesh::CreateGeoElement type element not exists:"
		 << " type = " << type << endl;
	 return NULL;
   }

   return NULL;
}
*/
void TPZGeoMesh::DeleteElement(TPZGeoEl *gel,int index){
  if(index < 0 || gel != fElementVec[index]){
    index = ElementIndex(gel);
    if(index < 0) {
      PZError << "TPZGeoMesh::DeleteElement index error\n";
      return;
    }
  }
  if (gel->HasSubElement())
  {
    for (int i=0;i<gel->NSubElements();i++)
    {
      TPZGeoEl* son = gel->SubElement(i);
      DeleteElement(son,son->Index());
    }
  }
  gel->RemoveConnectivities();
  if(gel) delete gel;
  fElementVec[index] = NULL;
  fElementVec.SetFree(index);
}

/** Verifies if the side based refinement pattern exists. If the refinement pattern doesn't exists return a Null refinement Pattern. */
TPZAutoPointer<TPZRefPattern> TPZGeoMesh::GetRefPattern (TPZGeoEl *gel, int side){
  MElementType type = gel->Type();
  //construct the named for the refinement pattern.
  string name = "SIDE_";
  TPZGeoElSide gelside (gel,side);
  int dimension = gelside.Dimension();
  switch (dimension){
    case (0) :{
      name += "NODE_";
      break;
    }
    case (1) :{
      name += "RIB_";
      break;
    }
    case (2) :{
      name += "FACE_";
      break;
    }
    case (3) :{
      name += "VOLUME_";
      break;
    }
    default :{
      name += "RIB_";
      break;
    }
  }
//  int size = 2;
  if (side/10 == 0){
    name += "0" ;
//    size = 1;
  }
  char aux[256];
  sprintf(aux,"%d",side);
  name += aux;

  switch (type){
    case (0) : {
      name += "_POINT";
      break;
    }
    case (1) : {
      name += "_LINE";
      break;
    }
    case (2) : {
      name += "_TRIANGLE";
      break;
    }
    case (3) : {
      name += "_QUAD";
      break;
    }
    case (4) : {
      name += "_TETRA";
      break;
    }
    case (5) : {
      name += "_PYRAMID";
      break;
    }
    case (6) : {
      name += "_PRISM";
      break;
    }
    case (7) : {
      name += "_HEXA";
      break;
    }
    default:{
      PZError << "TPZGeoMesh::GetRefPattern ERROR : Undefined element type " << type << std::endl;
      return 0;
    }
  }
  return GetRefPattern(type,name);
}

int TPZGeoMesh::ImportRefPattern(){
  std::string StartingPath;
#ifndef BORLAND
  StartingPath = REFPATTERNDIR;
#else
  StartingPath = "NeoPZ/Refine/RefPatterns";
#endif
  std::string FileTypes ("*.rpt");
  std::string Command = std::string("ls -1 ") + StartingPath + std::string("/") + FileTypes;
  std::cout << "Generated command: " << Command.c_str() << std::endl;
  FILE   *fp;
#ifndef BORLAND
  fp = popen(Command.c_str(), "r");
#else
  fp = (FILE *)open(Command.c_str(), O_RDONLY );
#endif

  if (!fp) return -1;
  int count = 0;
  char psBuffer[1024];
  while( !feof( fp ) )
  {
    if( fgets(psBuffer, sizeof(psBuffer), fp ) != NULL )
    {
      if (psBuffer[strlen(psBuffer)-1] == '\n') psBuffer[strlen(psBuffer)-1] = 0;
      std::cout << "Reading refinement patern file : " << psBuffer << std::endl;
      std::string filref (psBuffer);
      TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern( this,filref );
      this->InsertRefPattern(refpat);
      count++;
    }
  }
#ifndef BORLAND
  pclose(fp);
#else
  close((int)fp);
#endif
  return count;
}

int TPZGeoMesh::ClassId() const {
  return TPZGEOMESHID;
}

#ifndef BORLAND
template class TPZRestoreClass<TPZGeoMesh,TPZGEOMESHID>;
#endif

void TPZGeoMesh::Read(TPZStream &buf, void *context)
{
  TPZSaveable::Read(buf,context);
  buf.Read(&fName,1);
  ReadObjects(buf,fNodeVec,this);
  ReadObjectPointers(buf,fElementVec,this);
  buf.Read(&fNodeMaxId,1);
  buf.Read(&fElementMaxId,1);
  int ninterfacemaps;
  buf.Read(&ninterfacemaps,1);
  int c;
  for(c=0; c< ninterfacemaps; c++){
    int vals[3];
    buf.Read(vals,3);
    fInterfaceMaterials[pair<int,int>(vals[0],vals[1])]=vals[2];
  }

  //Reading TPZRefPattern's
  this->fRefPatterns.clear();
  int bigmapsize, iRef;
  buf.Read(&bigmapsize, 1);
  for(iRef = 0; iRef < bigmapsize; iRef++){
    int intElementType;
    buf.Read(&intElementType, 1);
    MElementType MElType = static_cast<MElementType>(intElementType);
    int smallmapsize, iMap;
    buf.Read(&smallmapsize, 1);
    for(iMap = 0; iMap < smallmapsize; iMap++){
      TPZRefPattern * refp = new TPZRefPattern(this);
      refp->Read(buf);
      this->fRefPatterns[MElType][refp->Id()] = refp;
    }//for
  }//for
  //Reading TPZRefPattern's
}

void TPZGeoMesh::Write(TPZStream &buf, int withclassid)
{
  TPZSaveable::Write(buf,withclassid);
  buf.Write(&fName,1);
  WriteObjects(buf,fNodeVec);
  WriteObjectPointers(buf,fElementVec);
  buf.Write(&fNodeMaxId,1);
  buf.Write(&fElementMaxId,1);
  int ninterfacemaps = fInterfaceMaterials.size();
  buf.Write(&ninterfacemaps,1);
  InterfaceMaterialsMap::iterator it = fInterfaceMaterials.begin();
  for(; it != fInterfaceMaterials.end(); it++){
    int vals[3];
    vals[0] = (it->first).first;
    vals[1] = (it->first).second;
    vals[2] = it->second;
    buf.Write(vals,3);
  }

  //Writing TPZRefPattern's
  std::map< MElementType, std::map<int, TPZAutoPointer<TPZRefPattern> > >::iterator eRef,itRef;
  int bigmapsize = fRefPatterns.size();
  buf.Write(&bigmapsize, 1);
  eRef = this->fRefPatterns.end();
  for(itRef = this->fRefPatterns.begin(); itRef != eRef; itRef++){
    int intElementType = static_cast<int>(itRef->first);
    buf.Write(&intElementType, 1);
    std::map<int, TPZAutoPointer<TPZRefPattern> >::iterator eMap, itMap;
    std::map<int, TPZAutoPointer<TPZRefPattern> > &SmallMap = itRef->second;
    int smallmapsize = SmallMap.size();
    buf.Write(&smallmapsize, 1);
    eMap = SmallMap.end();
    for(itMap = SmallMap.begin(); itMap != eMap; itMap++){
      TPZAutoPointer<TPZRefPattern> refp = itMap->second;
      refp->Write(buf);
    }//for
  }//for
  //Finishing writing TPZRefPattern's
}//method

int TPZGeoMesh::AddInterfaceMaterial(int leftmaterial, int rightmaterial, int interfacematerial){
  std::pair<int, int> leftright(leftmaterial, rightmaterial);
  InterfaceMaterialsMap::iterator w, e;
  e = fInterfaceMaterials.end();

  w = fInterfaceMaterials.find(leftright);
  if (w == e) { //std::pair leftright does not exist yet
      fInterfaceMaterials[leftright] = interfacematerial;
      return 1;
  }
  return 0;
}

int TPZGeoMesh::InterfaceMaterial(int leftmaterial, int rightmaterial){
  std::pair<int, int> leftright(leftmaterial, rightmaterial);
  InterfaceMaterialsMap::iterator w, e;
  e = fInterfaceMaterials.end();

  //trying to find an interface material associated to left and right materials
  w = fInterfaceMaterials.find(leftright);
  if (w != e) return w->second;

  //when left and right are equal and no exception case was inserted in interfacematerialmap return the same material of left and right
  if (leftmaterial == rightmaterial) return leftmaterial;

  //error message
  std::stringstream mess;
  mess << "\nTPZGeoMesh::InterfaceMaterial - Interface material not found ";
  PZError << mess.str()  << std::endl;
  return -9999;
}

void TPZGeoMesh::ClearInterfaceMaterialsMap(){
  InterfaceMaterialsMap::iterator b, e;
  b = fInterfaceMaterials.begin();
  e = fInterfaceMaterials.end();
  fInterfaceMaterials.erase(b, e);
}

void TPZGeoMesh::ResetConnectivities(){
  TPZGeoElSide side;
  int iel;
  const int nelem = this->NElements();
  for(iel = 0; iel < nelem; iel++){
    TPZGeoEl * gel = this->ElementVec()[iel];
    if (!gel) continue;
    const int nsides = gel->NSides();
    int is;
    for(is = 0; is < nsides; is++){
      gel->SetNeighbour(is, side);
    }
  }
}

const std::map<int, TPZAutoPointer<TPZRefPattern> > &TPZGeoMesh::RefPatternList(MElementType eltype){
  return fRefPatterns[eltype];
}

void TPZGeoMesh::InitializeRefPatterns()
{
  //line
  {
    char buf[] =
        "3	3\n"
        "202	UNIFORM_LINE\n"
        "-1.	0.	0. "
        "1.	0.	0. "
        "0.	0.	0. "
        "1	2		0	1 "
        "1	2		0	2 "
        "1	2		2	1 ";
    std::istringstream str(buf);
    TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(this,str);
    if(!FindRefPattern(refpat)) InsertRefPattern(refpat);
    refpat->InsertPermuted();
  }
  //triangle
  {
    char buf[] =
        "6  5\n"
        "399 UNIFORM_TRIANGLE\n"
        "0.  0.  0. "
        "1.  0.  0. "
        "0.  1.  0. "
        "0.5 0.  0. "
        "0.5 0.5 0. "
        "0.  0.5 0. "
        "2 3   0 1 2 "
        "2 3   0 3 5 "
        "2 3   3 1 4 "
        "2 3   5 4 2 "
        "2 3   4 5 3 ";
    std::istringstream str(buf);
    TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(this,str);
    if(!FindRefPattern(refpat)) InsertRefPattern(refpat);
    refpat->InsertPermuted();
  }
  //quadrilateral
  {
    char buf[] =
        "9 5\n"
        "499 UNIFORM_QUAD\n"
        "-1. -1. 0. "
        "1.  -1. 0. "
        "1.  1.  0. "
        "-1. 1.  0. "
        "0. -1.  0. "
        "1.  0.  0. "
        "0.  1.  0. "
        "-1. 0.  0. "
        "0.  0.  0. "
        "3 4   0 1 2 3 "
        "3 4   0 4 8 7 "
        "3 4   4 1 5 8 "
        "3 4   8 5 2 6 "
        "3 4   7 8 6 3 ";
    std::istringstream str(buf);
    TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(this,str);
    if(!FindRefPattern(refpat)) InsertRefPattern(refpat);
    refpat->InsertPermuted();
  }

    //hexahedre
  {
    std::cout << "\n\ninserting hexahedre\n";
    char buf[] =
        "27 9\n"
        "899 HALF_HEXA\n"
        "-1.0 -1.0 -1.0 "
        " 1.0 -1.0 -1.0 "
        " 1.0  1.0 -1.0 "
        "-1.0  1.0 -1.0 "
        "-1.0 -1.0  1.0 "
        " 1.0 -1.0  1.0 "
        " 1.0  1.0  1.0 "
        "-1.0  1.0  1.0 "
        " 0.0 -1.0 -1.0 "
        " 1.0  0.0 -1.0 "
        " 0.0  1.0 -1.0 "
        "-1.0  0.0 -1.0 "
        "-1.0 -1.0  0.0 "
        " 1.0 -1.0  0.0 "
        " 1.0  1.0  0.0 "
        "-1.0  1.0  0.0 "
        " 0.0 -1.0  1.0 "
        " 1.0  0.0  1.0 "
        " 0.0  1.0  1.0 "
        "-1.0  0.0  1.0 "
        " 0.0  0.0 -1.0 "
        " 0.0 -1.0  0.0 "
        " 1.0  1.0  0.0 "
        " 0.0  1.0  0.0 "
        "-1.0  0.0  0.0 "
        " 0.0  0.0  1.0 "
        " 0.0  0.0  0.0 "
        "7 8   0  1  2  3  4  5  6  7  "
        "7 8   0  8  20 11 12 21 26 24 "
        "7 8   8  1  9  20 21 13 22 26 "
        "7 8   20 9  2  10 26 22 14 23 "
        "7 8   11 20 10 3  24 26 23 15 "
        "7 8   12 21 26 24 4  16 25 19 "
        "7 8   21 13 22 26 16 5  17 25 "
        "7 8   26 22 14 23 25 17 6  18 "
        "7 8   24 26 23 15 19 25 18 7  ";
    std::istringstream str(buf);
    TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(this,str);
    if(!FindRefPattern(refpat)) InsertRefPattern(refpat);
    refpat->InsertPermuted();
  }


  //tetrahedre
  {
    std::cout << "\n\ninserting tetrahedre\n";
    char buf[] =
        "10  7\n"
        "799 UNIFORM_TETRA\n"
        "0.  0.  0. "
        "1.  0.  0. "
        "0.  1.  0. "
        "0.  0.  1. "
        "0.5 0.  0. "
        "0.5 0.5 0. "
        "0.  0.5 0. "
        "0.  0.  0.5 "
        "0.5 0.  0.5 "
        "0.  0.5 0.5 "
        "4 4   0 1 2 3 "
        "4 4   0 4 6 7 "
        "4 4   4 1 5 8 "
        "4 4   6 5 2 9 "
        "4 4   7 8 9 3 "
        "5 5   4 8 9 6 7 "
        "5 5   8 4 6 9 5 ";
        std::istringstream str(buf);
    TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(this,str);
    if(!FindRefPattern(refpat)) InsertRefPattern(refpat);
    refpat->InsertPermuted();
  }
  //pyramid
  {
    std::cout << "\n\ninserting pyramid\n";
    char buf[] =
        "14  11\n"
        "599 UNIFORM_PYRAMID\n"
        "-1.0 -1.0  0.0 "
         "1.0 -1.0  0.0 "
         "1.0  1.0  0.0 "
        "-1.0  1.0  0.0 "
         "0.0  0.0  1.0 "
         "0.5 -1.0  0.0 "
         "1.0  0.0  0.0 "
         "0.0  1.0  0.0 "
        "-1.0  0.0  0.0 "
        "-0.5 -0.5  0.5 "
         "0.5 -0.5  0.5 "
         "0.5  0.5  0.5 "
        "-0.5  0.5  0.5 "
         "0.0  0.0  0.0 "
        "5 5   0  1  2  3  4 "
        "5 5   0  5 13  8  9 "
        "5 5   5  1  6  13 10 "
        "5 5   13 6  2  7  11 "
        "5 5   8  13 7  3  12 "
        "5 5   9  10 11 12 4 "
        "5 5   9  12 11 10 13 "
        "4 4   9  5  13 10 "
        "4 4   10 6  13 11 "
        "4 4   11 7  12 13 "
        "4 4    8 13 12 9 ";
    std::istringstream str(buf);
    TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(this,str);
    if(!FindRefPattern(refpat)) InsertRefPattern(refpat);
    refpat->InsertPermuted();
  }
  //prism
  {
    std::cout << "\n\ninserting prism\n";
    char buf[] =
        "18  9\n"
        "699 UNIFORM_PRISM\n"
        "0.  0.  -1. "
        "1.  0.  -1. "
        "0.  1.  -1. "
        "0.  0.  1. "
        "1.  0.  1. "
        "0.  1.  1. "
        "0.5 0.  -1. "
        "0.5 0.5 -1. "
        "0.  0.5 -1. "
        "0.  0.  0. "
        "1.  0.  0. "
        "0.  1.  0. "
        "0.5 0.  1. "
        "0.5 0.5 1. "
        "0.  0.5 1. "
        "0.5 0.  0. "
        "0.5 0.5 0. "
        "0.  0.5 0. "
        "6 6 0 1 2 3 4 5 "
        "6 6 0 6 8 9 15  17 "
        "6 6 6 1 7 15  10  16 "
        "6 6 8 7 2 17  16  11 "
        "6 6 17  16  15  8 7 6 "
        "6 6 9 15  17  3 12  14 "
        "6 6 15  10  16  12  4 13 "
        "6 6 17  16  11  14  13  5 "
        "6 6 14  13  12  17  16  15";
    std::istringstream str(buf);
    TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(this,str);
    if(!FindRefPattern(refpat)) InsertRefPattern(refpat);
    refpat->InsertPermuted();
  }
}

