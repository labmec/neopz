// -*- c++ -*-
#include "pzgeoelside.h"
#include "pzgeoel.h"
#include "pztrnsform.h"
#include "pzstack.h"
#include "pzvec_extras.h"
#include "pzquad.h"
#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzcompel.h"
#include "pzintel.h"


void TPZGeoElSide::RemoveConnectivity(){

  if(!Exists()) return;
  if(fSide < 0 || fSide >= fGeoEl->NSides()) {
    PZError << "TPZGeoElSide::SetConnectivity Index out of bound\n";
  }
  //it removes the connectivity of the cycle where this inserted one: 
  //neighpre->this->neighpos => neighpre->neighpos
  TPZGeoElSide neighpre,neigh = Neighbour();
  if(neigh.Element() == NULL || neigh.Side() == -1){
    PZError << "TPZGeoElSide::SetConnectivity trying to remove null or inexistent connection";
  }
  TPZGeoElSide neighpos = neigh;
  while(neigh.Exists() && neigh != *this){
    neighpre = neigh;
    neigh = neigh.Neighbour();
  }
  if(neigh == *this){
    this->SetNeighbour(TPZGeoElSide());
    neighpre.SetNeighbour(neighpos);
  } else {
    PZError << "TPZGeoElSide::SetConnectivity neighbourhood cycle error";
  }
}

void TPZGeoElSide::SetConnectivity(const TPZGeoElSide &neighbour) const {
  
  if(!Exists()) return;
  if(fSide >= fGeoEl->NSides()) {
    PZError << "ERROR(TPZGeoEl::SetConnectivity)-> Index greater than number of sides.\n";
    PZError << " fNumSides = " << fGeoEl->NSides() << " side = " << fSide << "\n";
  }
  if(!neighbour.Exists()) {
    fGeoEl->SetSideDefined(fSide);
  }
  
  TPZGeoElSide neighneigh, currentneigh;
  neighneigh = neighbour.Neighbour();
  currentneigh = Neighbour();
  /**The neighbour has the same level as the current element
     the neighbour of the neighbour has a smaller level (i.e. larger element)
     the neighbour will now point to the current element and vice versa*/
  if (!neighneigh.Exists() && !currentneigh.Exists()) {
    SetNeighbour(neighbour);
    neighbour.SetNeighbour(*this);
  } else if (neighneigh.Exists() && currentneigh.Exists()) {
    // It would be convenient to check the consistency of both loops
    // insert the connectivity between two independent loops
    if(NeighbourExists(neighbour) || neighbour.NeighbourExists(*this)) {
      PZError << "TPZGeoEl::SetConnectivity Fourth untreated case, wrong data structure\n";
    } else {
      SetNeighbour(neighneigh);
      neighbour.SetNeighbour(currentneigh);
    }
    /**The neighbouring element has already a connectivity loop, insert in his loop*/
  } else if (neighneigh.Exists()) {
    // the neighbour already has a loop, insert this into the loop
    SetNeighbour(neighneigh);
    neighbour.SetNeighbour(*this);
    /**The current element has already a connectivity loop and the neighbour doesnt*/
  } else if (currentneigh.Exists()) {
    // this is already inserted in a loop, insert neighbour in the loop
    SetNeighbour(neighbour);
    neighbour.SetNeighbour(currentneigh);
  }
}

/**Null object*/
TPZGeoElSide::TPZGeoElSide() {
  fGeoEl = 0;
  fSide  = -1;
}
/**Copy object*/
TPZGeoElSide::TPZGeoElSide(const TPZGeoElSide &gelside) {
  fGeoEl = gelside.fGeoEl;
  fSide  = gelside.fSide;
if(fSide > 26)
  cout << "TPZGeoElSide bad side" << endl;
}
/**object construct with element gel and fSide side*/
TPZGeoElSide::TPZGeoElSide(TPZGeoEl *gel,int side) {
  fGeoEl = gel;
  fSide  = side;
if(fSide > 26)
   cout << "TPZGeoElSide bad parameter side" << endl;
}

TPZGeoEl * TPZGeoElSide::Element() {
  return fGeoEl;
}

int TPZGeoElSide::Side() {
  return fSide;
}

TPZGeoElSide TPZGeoElSide::Neighbour() const {
  if (!fGeoEl) return TPZGeoElSide();
  TPZGeoElSide neighbour = fGeoEl->Neighbour(fSide);
  return neighbour;
}

void TPZGeoElSide::AllNeighbours(TPZStack<TPZGeoElSide> &allneigh) {
	if(! Exists() || ! this->Neighbour().Exists()) 
	{
		return;
	}
	TPZGeoElSide neigh = Neighbour();
	while(neigh != *this)
	{
		allneigh.Push(neigh);
		neigh = neigh.Neighbour();
	}
}

void TPZGeoElSide::ComputeNeighbours(TPZStack<TPZGeoElSide> &compneigh) {
  if(fSide < fGeoEl->NCornerNodes()) 
    {
      AllNeighbours(compneigh);
      return;
    }
  int nsnodes = NSideNodes();
  TPZStack<TPZGeoElSide> GeoElSideSet;
  TPZStack<TPZGeoEl *> GeoElSet[27];
  int in;
  TPZManVector<int> nodeindexes(nsnodes);
  for(in=0; in<nsnodes; in++) 
    {
      nodeindexes[in] = SideNodeIndex(in);
      int locnod = fGeoEl->SideNodeLocIndex(fSide,in);
      GeoElSideSet.Resize(0);
      TPZGeoElSide locside(fGeoEl,locnod);
      locside.AllNeighbours(GeoElSideSet);
      int nel = GeoElSideSet.NElements();
      int el;
      for(el=0; el<nel; el++) {
	GeoElSet[in].Push(GeoElSideSet[el].Element());
      }
      Sort<TPZGeoEl *>(GeoElSet[in]);
    }
  TPZStack<TPZGeoEl *,100> result;
  switch(nsnodes) {
  case 1:
    {
      result = GeoElSet[0];
    }
    break;
  case 2:
    Intersect<TPZGeoEl *,100>(GeoElSet[0],GeoElSet[1],result);
    break;
  case 3:
    Intersect<TPZGeoEl *,100>(GeoElSet[0],GeoElSet[1],GeoElSet[2],result);
    break;
  case 4:
    {
      TPZStack<TPZGeoEl *> inter1, inter2;
      Intersect<TPZGeoEl *,100>(GeoElSet[0],GeoElSet[2],inter1);
      if(inter1.NElements()==0) break;
      Intersect<TPZGeoEl *,100>(GeoElSet[1],GeoElSet[3],inter2);
      if(inter2.NElements()==0) break;
      Intersect<TPZGeoEl *,100>(inter1,inter2,result);
    }
    break;
  default:
    {
      TPZStack<TPZGeoEl *> inter1, inter2;
      inter1 = GeoElSet[0];
      for(in=0; in<nsnodes-1; in++) {
	inter2.Resize(0);
	Intersect<TPZGeoEl *,100>(inter1,GeoElSet[in+1],inter2);
	if(inter2.NElements() == 0) break;
	inter1 = inter2;
      }
      result = inter2;
    }
  }
  int el,nel = result.NElements();
  for(el=0; el<nel; el++) {
    compneigh.Push(TPZGeoElSide(result[el],result[el]->WhichSide(nodeindexes)));
  }
}


/**Calcula a transformação entre o lado side do el. atual e o lado do vizinho que contem side*/
/*
void TPZGeoElSide::SideTransform(TPZGeoElSide neighbour,TPZTransform &t)	{
  
  TPZGeoElSide El;//vizinho do elemento El pelo seu lado side.
  if (!neighbour.Exists()) {
    PZError<<endl<<"Neighbour element does not exist"<<endl;
    return;
  }
  int sidedimension = fGeoEl->SideDimension(Side());//dimensao da transformacao
  while(!neighbour.NeighbourExists(El)) {
	  El = El.Father2();
	  if(!El.Exists()) {
		PZError << "TPZGeoEl::SideTransform called with wrong parameters, bye!\n";
		exit(-1);
	  }
  }
  TPZTransform tloc(sidedimension),tside(sidedimension);//transformação local

  Element()->BuildTransform2(Side(),El.Element(),tloc);

  t = El.NeighbourSideTransform(neighbour).Multiply(tloc);
}
*/

TPZTransform TPZGeoElSide::NeighbourSideTransform(TPZGeoElSide &neighbour) {
  int sidedimension = Dimension();
  TPZTransform tside(sidedimension);//transformação local
  switch (sidedimension) {
  case 0://canto para canto viz

    break;

  case 1://aresta para aresta viz
    if (SideNodeIndex(0) == neighbour.SideNodeIndex(0)) {
      tside.Mult()(0,0) =  1.;
    }
    else {
      tside.Mult()(0,0) = -1.;
    }
    break;

  case 2://transformacoes entre faces viz
  	 int i;
    //TPZCompEl *cel = Element()->Reference();
    TPZVec<int> idto(0),idfrom(0);
	 if(Element()->NSideNodes(Side()) == 4) {//faces quadrilaterais
     	idto.Resize(4);
      idfrom.Resize(4);
      tside.Mult()(0,0) = 0;
      tside.Mult()(1,1) = 0;
      for(i=0;i<4;i++) idto[i]=neighbour.Element()->SideNodeIndex(neighbour.Side(),i);
      //if(neighbour.Element()->NSides() > 9) {//NCornerNodes() == 8 , cubo
      //   neighbour.Element()->NodeFaceIds(idto,neighbour.Side());
      //}
      for(i=0;i<4;i++) idfrom[i]=Element()->SideNodeIndex(Side(),i);
      //	   if(Element()->NSides() > 9) {//NCornerNodes() == 8 , cubo
      //   Element()->NodeFaceIds(idfrom,Side());
      //}
      //if(neighbour.Element()->NCornerNodes() == 4) {//quadrilatero
      //	      for(i=0;i<4;i++) idto[i] = neighbour.Element()->NodeIndex(i);
      //}
      //if(Element()->NCornerNodes() == 4) {//quadrilatero
      //	      for(i=0;i<4;i++) idfrom[i] = Element()->NodeIndex(i);
      //}
      int transid = Element()->GetTransformId2dQ(idfrom,idto);
      //TPZCompEl *cel = Element()->Reference();
      tside.Mult()(0,0) = TPZShapeQuad::gTrans2dQ[transid][0][0];//cel->gTrans2dQ[transid][0][0];
      tside.Mult()(0,1) = TPZShapeQuad::gTrans2dQ[transid][0][1];
      tside.Mult()(1,0) = TPZShapeQuad::gTrans2dQ[transid][1][0];
      tside.Mult()(1,1) = TPZShapeQuad::gTrans2dQ[transid][1][1];
    } else if(Element()->NSideNodes(Side()) == 3) {//faces triangulares
     	idto.Resize(3);
      idfrom.Resize(3);
      tside.Mult()(0,0) = 0.;
      tside.Mult()(1,1) = 0.;
      for(i=0;i<3;i++) idfrom[i] = Element()->SideNodeIndex(Side(),i);
//	   if(Element()->NSides() > 9) {//elementos 3D
//         Element()->NodeFaceIds(idfrom,Side());
//      } else {//triângulos
//         for(i=0;i<3;i++) idfrom[i] = Element()->SideNodeIndex(Side(),i);
//      }
      for(i=0;i<3;i++) idto[i] = neighbour.Element()->SideNodeIndex(neighbour.Side(),i);
//      if(neighbour.Element()->NSides() > 9) {//elementos 3D
//         neighbour.Element()->NodeFaceIds(idto,neighbour.Side());
//      } else {//triângulos
//         for(i=0;i<3;i++) idto[i] = neighbour.Element()->NodeIndex(i);
//      }
  		int transid = Element()->GetTransformId2dT(idfrom,idto);
      //TPZCompEl *cel = Element()->Reference();
      tside.Mult()(0,0) = TPZShapeTriang::gTrans2dT[transid][0][0];
      tside.Mult()(0,1) = TPZShapeTriang::gTrans2dT[transid][0][1];
      tside.Mult()(1,0) = TPZShapeTriang::gTrans2dT[transid][1][0];
      tside.Mult()(1,1) = TPZShapeTriang::gTrans2dT[transid][1][1];
      tside.Sum()(0,0) = TPZShapeTriang::gVet2dT[transid][0];
      tside.Sum()(1,0) = TPZShapeTriang::gVet2dT[transid][1];
    } else {
     	PZError << "TPZGeoElSide::NeighbourSideTransform : elemento desconhecido" << endl;
    }
    break;
  }
  return tside;
}

int TPZGeoElSide::NeighbourExists(const TPZGeoElSide &gel) const {
  if(gel == *this) return 1;
  TPZGeoElSide neighbour = Neighbour();
  if(!neighbour.Exists()) return 0;
  while(neighbour != *this) {
    if(gel == neighbour) return 1;
    neighbour = neighbour.Neighbour();
  }
  return 0;
}

/** adds himself or his children to the stack if they have a reference element*/
/*
void TPZGeoElSide::SmallConnect(int level,TPZStack<TPZCompElSide> &elsidevec,
				int onlyinterpolated) {

  TPZCompElSide ref = Reference();
  TPZGeoElSide neighbour(*this);
  int locallevel = Level();
  if(locallevel >= level) {
    do {
      ref = neighbour.Reference();
      if(ref.Element() && (!onlyinterpolated || ref.Element()->IsInterpolated())) {
	elsidevec.Push(ref);
	return;
      }
      neighbour = neighbour.Neighbour();
    } while (neighbour != *this && neighbour.Element());
  }                              //Cedric
  neighbour = *this;
  do {
    if(!neighbour.HasSubElement())
      neighbour = neighbour.Neighbour();
    else {
      int is,numsidenodes = neighbour.NSideNodes();
      int refstore[8];
      TPZManVector<int> referencenode(numsidenodes,refstore,8);
      for(is=0; is<numsidenodes; is++)
        referencenode[is] = neighbour.SideNodeIndex(is);
      TPZGeoElSide subelvecstore[10];
      TPZManVector<TPZGeoElSide> subelvec(0,subelvecstore,10);
      neighbour.GetSubElement(referencenode,subelvec);
      int numsub = subelvec.NElements();
      for(is=0; is<numsub; is++)
        subelvec[is].SmallConnect(level,elsidevec,onlyinterpolated);
      return;
    }                                      //Cedric
  } while (neighbour.Element() != Element() && neighbour.Element());
}
*/
TPZCompElSide TPZGeoElSide::Reference() {
  TPZCompElSide compside(fGeoEl->Reference(),fSide);
  return compside;
}

/*
void TPZGeoElSide::GetSubElement(TPZVec<int> &referencenodes,
				 TPZVec<TPZGeoElSide> &subelements) {
  
  if(!fGeoEl || !fGeoEl->HasSubElement()) {   // Jorge 10/01/2000
    subelements.Resize(0);
  } else {
    fGeoEl->GetSubElement(fSide,referencenodes,subelements);
  }
}
*/
/*
void TPZGeoElSide::GetSubElements(TPZVec<TPZGeoElSide> &subelements) {
	if(!fGeoEl || !fGeoEl->HasSubElement()) {   // Jorge 10/01/2000
		subelements.Resize(0);
		return;
	}
	int nsub = fGeoEl->NSideSubElements(fSide);
	subelements.Resize(nsub);
	int is;
	for(is=0; is<nsub; is++) {
		subelements[is] = fGeoEl->SideSubElement(fSide,is);
	}
}
*/

int TPZGeoElSide::Dimension() {
  if (!fGeoEl) {
    PZError << "TPZGeoElSide::Dimension : null element\n";
    return -1;
  }
  return fGeoEl->SideDimension(fSide);
}

/*
TPZCompElSide TPZGeoElSide::HigherDimensionConnected(int targetdimension,int onlyinterpolated){

  if(!fGeoEl) return TPZCompElSide();
  int currentdim = Dimension();
  if(currentdim >= targetdimension) return TPZCompElSide();

  TPZGeoElSide father = *this;
  TPZGeoElSide nextfather = father.Father();
  while(nextfather.Element()) {
    father = nextfather;
    nextfather = father.Father();
  }
  //aqui o pai de menor nivel achado e' father
  TPZGeoElSide largeside;
  TPZGeoElSide neighbour = father;
  // check on all neighbours to see whether one of them has a
  // higher dimension side
  do {
    largeside = neighbour.Element()->HigherDimensionSides(neighbour.Side(),targetdimension);
    neighbour = neighbour.Neighbour();
  } while (!largeside.Exists() && neighbour.Exists() && neighbour != father);
  TPZGeoElSide nextneigh;
  neighbour = TPZGeoElSide();
  if (largeside.Exists())
    neighbour = largeside.Father();
  while(neighbour.Exists()) {
    // loop over the neighbours of neighbour
    TPZCompEl *ref = neighbour.Element()->Reference();
    if(ref && (!onlyinterpolated || ref->IsInterpolated())) {
      // the neighbour is the father of a geometric element with reference
      // therefore, the neighbour should not have a higher dimension connected
      // element (except for comparaisons between distinct grids)
      // PZError << "TPZGeoElSide::HigherDimensionConnected I don't understand\n";         //   !!!   ???
      return neighbour.Reference();
    }
    nextneigh = neighbour.Neighbour();
    if(nextneigh.Exists()) {
      while(nextneigh != neighbour) {
	TPZCompEl *ref = nextneigh.Element()->Reference();
	if(ref && (!onlyinterpolated || ref->IsInterpolated())) {
	  return nextneigh.Reference();
	}
	nextneigh = nextneigh.Neighbour();
      }
    }
    neighbour = neighbour.Father();
  }
  return TPZCompElSide();
}
*/
/*
void TPZGeoElSide::SideTransform2(TPZGeoElSide neighbour,TPZTransform &t)	{
  //t : atual -> neighbour
  TPZGeoElSide father = *this;
  if(!father.Exists()) return;
  if(father.NeighbourExists(neighbour)) {
    t = NeighbourSideTransform(neighbour).Multiply(t);
    return;
  }
  TPZGeoElSide nextfather = father.Father();
  while(nextfather.Exists()) {
    //father.Element()->BuildTransform(Side(),nextfather.Element(),t);
    father.Element()->BuildTransform2(father.Side(),nextfather.Element(),t);//Cedric 01/10/99 e  e 30/04/00
    father = nextfather;
    if(father.NeighbourExists(neighbour)) {
      t = father.NeighbourSideTransform(neighbour).Multiply(t);
      return;
    }
    nextfather = father.Father();
  }

  TPZGeoElSide neighfather = father;
  TPZGeoElSide largeside = father.HigherDimensionSide(neighbour.Dimension());
  if(!largeside.Exists() || largeside.Element()!=neighfather.Element()) {
    do {
      neighfather = neighfather.Neighbour();
      largeside = neighfather.HigherDimensionSide(neighbour.Dimension());
    } while(!largeside.Exists() && largeside.Element()!=neighfather.Element() && neighfather != father);
  }
  if(!largeside.Exists()) {
    PZError << "TPZGeoElSide:SideTranform2 did not find the neighbour\n";
    return;
  }

  t = father.NeighbourSideTransform(neighfather).Multiply(t);//transf. entre father e seu viz

  if(neighfather.Element() != largeside.Element()) {
     //os subelementos são de tipos diferentes mais são irmão e viz por ser obtido por HigherDim...
     //percorrer a viz de neighfather até achar largeside e o lado respectivo
     //o largeside contem a aresta ou face de chegada
     //então a transf. sera feita entre lados de dimensão diferente do mesmo elemento
     TPZGeoElSide largelowside = neighfather;
     while(largelowside.Exists() && largelowside.Element()!=largeside.Element()) {
         largelowside = largelowside.Neighbour();//viz de neighfather
     }//acha largeside como viz de neighfather
     if(!largelowside.Exists()) {
     	   PZError << "TPZGeoElSide::SideTransform2 did not find the neighbour\n";
         return;
     }
     t = neighfather.NeighbourSideTransform(largelowside).Multiply(t);//transf. entre father e seu viz

     neighfather = largelowside;
  }//Cedric 04/10/99 e 30/04/00

  TPZTransform temp = neighfather.SideToSideTransform(largeside);
  t = temp.Multiply(t);
  largeside.SideTransform2(neighbour,t);
}
*/

void TPZGeoElSide::SideTransform3(TPZGeoElSide neighbour,TPZTransform &t)	{
  //t : atual -> neighbour
  TPZGeoElSide father(*this);
  if(!father.Exists()) {
	  PZError << "TPZGeoElSide::SideTransform3 I dont understand\n";
	  return;
  }
  if(father.NeighbourExists(neighbour)) {
    t = NeighbourSideTransform(neighbour).Multiply(t);
    return;
  }
  TPZGeoElSide nextfather = father.Father2();
  while(nextfather.Exists()) {
    //father.Element()->BuildTransform(Side(),nextfather.Element(),t);
    t = father.Element()->BuildTransform2(father.Side(),nextfather.Element(),t);//Cedric 01/10/99 e  e 30/04/00
    father = nextfather;
    if(father.NeighbourExists(neighbour)) {
      t = father.NeighbourSideTransform(neighbour).Multiply(t);
      return;
    }
    nextfather = father.Father2();
  }

  PZError << "TPZGeoElSide:SideTranform3 did not find the neighbour\n";
  return;
}

void TPZGeoElSide::ConnectedCompElementList(TPZStack<TPZCompElSide> &ellist,
					 int onlyinterpolated, int removeduplicates) {
  if(!fGeoEl) return;
  TPZCompElSide father = LowerLevelCompElementList2(onlyinterpolated);
  if(father.Exists()) ellist.Push(father);
  EqualLevelCompElementList(ellist,onlyinterpolated,removeduplicates);
  HigherLevelCompElementList2(ellist,onlyinterpolated,removeduplicates);
}

/**retorna o elem. computacional de nivel menor (elemento grande) ao qual estou restrito*/
/*
TPZCompElSide TPZGeoElSide::LowerLevelCompElementList(int onlyinterpolated) {
  if(!Dimension()) return Dim0LowerLevelCompElement(onlyinterpolated);
  TPZGeoElSide gel,gelfather;
  TPZCompElSide ref;
  gel = Neighbour();
  gelfather = Father();
  // these statements will only become effective if the higher level
  // element will point to zero if it doesn't have an equal level neighbour
  while(gelfather.Exists()) {
    gel = gelfather.Neighbour();
    if(gel.Exists()) {
      while(gel != gelfather) {
	     ref = gel.Reference();
	     if(ref.Element() && (!onlyinterpolated || ref.Element()->IsInterpolated())) {             // (2)
          return ref;
        }
        gel = gel.Neighbour();
      }
    }
    gelfather = gelfather.Father();
  }
  int dimension = Dimension();
  for(int dim=dimension+1; dim < 3; dim++) {                       //   !!!   ???
    ref = HigherDimensionConnected(dim,onlyinterpolated);
    if(ref.Exists() && (!onlyinterpolated || ref.Element()->IsInterpolated())) {             // (2)
      return ref;
    }
  }
  return TPZCompElSide();
}
*/
/*
void TPZGeoElSide::HigherLevelCompElementList(TPZStack<TPZCompElSide> &elvec,
					   int onlyinterpolated, int removeduplicates) {

//Philippe 13/5/99 If the element/side combination has dimension 0, nothing to do
//  if(!fGeoEl->Dimension()) return;
  if(!Dimension()) return;
  /*
  int level = Level();
  int levelprocessed = level+1;

  SmallConnect(levelprocessed,elvec,onlyinterpolated);

  if(removeduplicates) Reference().RemoveDuplicates(elvec);
  
  TPZGeoElSide neighbour(*this);
  TPZStack<TPZGeoElSide> subel(20);
  do {
	  if(neighbour.HasSubElement()) {
		neighbour.GetSubElements(subel);
		int nsub = subel.NElements();
		int is;
		for(is=0; is<nsub; is++) {
			subel[is].EqualorHigherCompElementList(elvec,onlyinterpolated,removeduplicates);
		}
		neighbour = *this;
	  } else {
		  neighbour = neighbour.Neighbour();
	  }
	  if(!neighbour.Exists()) break;
  } while(neighbour != *this);
}
*/
void TPZGeoElSide::EqualLevelCompElementList(TPZStack<TPZCompElSide> &elsidevec,
					  int onlyinterpolated, int removeduplicates) {

//  if(!fGeoEl->Dimension())
// Philippe 26/4/2000
//  if(!Dimension())
//    Dim0EqualLevelCompElementList(elsidevec,onlyinterpolated,removeduplicates);

  TPZGeoElSide neighbour;
  TPZCompElSide ref;
  neighbour = Neighbour();
  if(!neighbour.Exists()) return;

  while(neighbour.Element() != this->Element()) {
    ref = neighbour.Reference();
    if(ref.Element() && ref.Element() != Reference().Element() && (!onlyinterpolated || ref.Element()->IsInterpolated())) {
      elsidevec.Push(ref);
      if(removeduplicates) return;
    }
    neighbour = neighbour.Neighbour();
  }
}

void TPZGeoElSide::HigherDimensionElementList(TPZStack<TPZCompElSide> &elsidevec, int onlyinterpolated) {

	TPZStack<TPZGeoElSide> gelsides;
	fGeoEl->AllHigherDimensionSides(fSide,2,gelsides);
	int il,cap = gelsides.NElements();
	for(il=0; il<cap; il++) {
		TPZCompElSide cels = gelsides[il].Reference();
		if(onlyinterpolated) {
			TPZInterpolatedElement *cel = dynamic_cast<TPZInterpolatedElement *> (cels.Element());
			if(!cel) continue;
			int conind = cel->ConnectIndex(cels.Side());
			if(conind < 0) continue;
		}
		elsidevec.Push(cels);
	}

}

/*
void TPZGeoElSide::Dim0EqualLevelCompElementList(TPZStack<TPZCompElSide> &elsidevec,int onlyinterpolated, int removeduplicates){

/**Look for the highest level element which contains this side
	TPZGeoElSide gel = *this, geltest = *this;
	TPZCompElSide ref;
   int nodeliststore[4];
   //int oldsize = elsidevec.NElements();
   TPZManVector<int> nodelist(1,nodeliststore,4);
   TPZGeoElSide subelstore[4];
   TPZManVector<TPZGeoElSide> subel(0,subelstore,4);
/**look for the highest level element which is linked to this side
   do {
   	geltest = gel;
/**go as far as possible along the tree structure
      while(gel.HasSubElement()) {
         nodelist[0] = SideNodeIndex(0);
         gel.GetSubElement(nodelist,subel);
         gel = subel[0];
      }
      TPZGeoElSide gelneigh = gel.Neighbour();
      if(!gelneigh.Exists()) continue;
/**at the end of the tree check is a neighbour has sub elements
      while(gelneigh != gel && !gelneigh.HasSubElement()) gelneigh = gelneigh.Neighbour();
      gel = gelneigh;
   } while (gel != geltest);
   while(gel.Exists()) {
/**put the equal level elements on the stack.
    if gel is a direct neighbour of this, call EqualLevelElementList for this so as to
    avoid putting this on the stack
		ref = gel.Reference();
		if(ref.Exists() && gel != *this && (!onlyinterpolated || ref.Element()->IsInterpolated())) {
			elsidevec.Push(ref);
/**if removeduplicates is true, only one element needs to be put on the stack
			if(removeduplicates) return;
		}
		TPZGeoElSide gelneigh = gel.Neighbour();
		while(gelneigh.Exists() && gelneigh != gel) {
			ref = gelneigh.Reference();
			if(ref.Exists() && gelneigh != *this && (!onlyinterpolated || ref.Element()->IsInterpolated())) {
				elsidevec.Push(ref);
/**if removeduplicates is true, only one element needs to be put on the stack
				if(removeduplicates) return;
			}
			gelneigh = gelneigh.Neighbour();
		}
		gel = gel.Father();
	}
}
*/
/**retorna o elem. computacional de nivel menor (elemento grande) ao qual estou restrito*/
/*
TPZCompElSide TPZGeoElSide::Dim0LowerLevelCompElement(int onlyinterpolated) {
  TPZGeoElSide gel,gelfather;
  TPZCompElSide ref;
  if(Dimension() != 0) {
  	PZError << "TPZGeoElSide::Dim0LowerLevelCompElement should be called for an\n"
   	"element/side combination of dimension 0\n";
   PZError << Element()->Id() << "/" << Side() << endl;
   return TPZCompElSide();
  }
  gel = Neighbour();
  gelfather = Father();
  while(gelfather.Exists()) {
    gel = gelfather;
    gelfather = gelfather.Father();
  }
  int dimension = Dimension();
  for(int dim=dimension+1; dim < 3; dim++) {                       //   !!!   ???
    ref = HigherDimensionConnected(dim,onlyinterpolated);
    if(ref.Exists() && (!onlyinterpolated || ref.Element()->IsInterpolated())) {             // (2)
      return ref;
    }
  }
  return TPZCompElSide();
}
*/


/*
TPZGeoElSide TPZGeoElSide::HigherDimensionSide(int targetdimension) {
    return fGeoEl->HigherDimensionSides(fSide,targetdimension);
}
*/
//EXCUSIVAMENTE PARA TESTAR CONTINUIDADE
//#include "pzelgt2d.h"
/*
void TPZGeoElSide::LargestNeighbour(TPZGeoElSide &largeel, int targetdimension) {

  TPZGeoElSide fatherkeep(Element(),Side());

  if(Dimension() > targetdimension) {
  		largeel = TPZGeoElSide();//fatherkeep;
      return;
  }

  TPZGeoElSide father = Father();

  while(father.Element()) {//father.Dimension() = 1,2
    fatherkeep = father;
    father = father.Father();
  }

  father = fatherkeep;
  if(Dimension() < targetdimension){
    TPZGeoEl *gel = Element();
    if(gel->NSides() == 6){
      father = ( Element())->HigherDimensionSides(Side(),targetdimension);
    }

  }
  if(!father.Element()) {
		largeel = *this;//fatherkeep;
  		return;
  }

  while(father.Element()) {//father.Dimension() = 1,2
    fatherkeep = father;
    father = father.Father();
  }

  largeel = fatherkeep.Neighbour();
  if(!largeel.Element()) largeel = fatherkeep;
}
*/
/**Jorge 13/01/2000
   Retorna todos os subelementos computacionais do atual elemento geometrico
   ao longo do lado fSide*/
/*
void TPZGeoElSide::GetCompSubElements(TPZStack<TPZCompElSide> &elvec,int onlyinterpolated,
     int removeduplicates) {

  TPZVec<TPZGeoElSide> subelements(10);
  TPZVec<TPZGeoElSide> subs;
  GetSubElements(subs);
  TPZCompEl *cel;
  int r, nallsubs = 0, nsubs = subs.NElements();
	
  while(nsubs) {
    for(r=0;r<nsubs-1;r++) {
      cel = subs[r].Element()->Reference();
      if(cel  && (!onlyinterpolated || cel->IsInterpolated())) elvec.Push(subs[r].Reference());
      else subelements[nallsubs++] = subs[r];
    }
    cel = subs[r].Element()->Reference();
    if(cel  && (!onlyinterpolated || cel->IsInterpolated())) {
      elvec.Push(subs[r].Reference());
      nsubs = 0;
    }
    else {
      subs[r].GetSubElements(subs);
      nsubs = subs.NElements();
    }
  }
  for(r=0;r<nallsubs;r++) {
    subelements[r].GetSubElements(subs);
    nsubs = subs.NElements();
//if(!nsubs) PZError << "TPZGeoElSide::GetCompSubElements. It has not sub-elements.\n";
    for(int p=0;p<nsubs;p++) {
      if(subs[p].Element()->Reference()) elvec.Push(subs[p].Reference());
      else subelements[nallsubs++] = subs[p];
    }
  }

  if(removeduplicates) Reference().RemoveDuplicates(elvec);
}
*/
int TPZGeoElSide::Id() {
	return fGeoEl->Id();
}

/**fill in the data structure for the neighbouring information*/
void TPZGeoElSide::SetNeighbour(const TPZGeoElSide &neighbour) const {
    if(!fGeoEl) {
      PZError << "TPZGeoElSide::SetNeighbour. Don't exist geometrical element.\n";
      return;
    }
    fGeoEl->SetNeighbour(fSide,neighbour);
}


TPZTransform TPZGeoElSide::SideToSideTransform(TPZGeoElSide &higherdimensionside){
    if(fGeoEl != higherdimensionside.fGeoEl) {
      PZError << "TPZGeoElSide::SideToSideTransform inconsistent id1 = " << fGeoEl->Id() << 
	" id2 = " << higherdimensionside.fGeoEl->Id() << endl;
    }
    return fGeoEl->SideToSideTransform(fSide,higherdimensionside.fSide);
}

//int TPZGeoElSide::Level() {
//		return (fGeoEl->Level());
//}
/*return 1 if the element has subelements along side*/
int TPZGeoElSide::HasSubElement() {
    if(!fGeoEl) return 0;
    return fGeoEl->HasSubElement();
}
/*return the number of nodes for a particular side*/
int TPZGeoElSide::NSideNodes() {
    if(!fGeoEl) return 0;
    return fGeoEl->NSideNodes(fSide);
}

/**returns the index of the nodenum node of side*/
int TPZGeoElSide::SideNodeIndex(int nodenum) {
    if(!fGeoEl) return -1;
    return ( fGeoEl->SideNodeIndex(fSide,nodenum) );
}


/*
void TPZGeoElSide::EqualorHigherCompElementList(TPZStack<TPZCompElSide> &celside, int onlyinterpolated, int removeduplicates)
{

	int ncelsides = celside.NElements();
	if(Reference().Exists()) {
		celside.Push(Reference());
		if(removeduplicates) {
			return;
		}
	}
	this->EqualLevelCompElementList(celside,onlyinterpolated,removeduplicates);
	if(ncelsides != celside.NElements()) return;
	TPZStack<TPZGeoElSide> gelsides(20);
	TPZGeoElSide neighbour(*this);
	do {
		if(neighbour.HasSubElement()) {
			neighbour.GetSubElements2(gelsides);
			int nsub = gelsides.NElements();
			int is;
			for(is=0; is<nsub; is++) {
				gelsides[is].EqualorHigherCompElementList(celside,onlyinterpolated,removeduplicates);
			}
			neighbour = *this;
		} else {
			neighbour = neighbour.Neighbour();
		}
		if(!neighbour.Exists()) break;
	} while(neighbour != *this);
}
*/
TPZCompElSide TPZGeoElSide::LowerLevelCompElementList2(int onlyinterpolated)
{
	TPZGeoElSide father = StrictFather();
	if(!father.Exists()) return TPZCompElSide();
	TPZStack<TPZCompElSide> equal;
	father.EqualLevelCompElementList(equal,onlyinterpolated,1);
	while(father.Exists() && equal.NElements() == 0) {
		father = father.StrictFather();
		father.EqualLevelCompElementList(equal,onlyinterpolated,1);
	}
	if(equal.NElements()) return equal[0];
	return TPZCompElSide();
}

TPZGeoElSide TPZGeoElSide::Father2()
{
	if(!fGeoEl) return TPZGeoElSide();
	return fGeoEl->Father2(fSide);
}

TPZGeoElSide TPZGeoElSide::StrictFather()
{
	TPZGeoElSide father = Father2();
	int nfathsub = 0;
	if(father.Exists()) nfathsub = father.fGeoEl->NSideSubElements2(father.fSide);
	while(father.Exists() && nfathsub == 1) {
		father = father.Father2();
		if(father.Exists()) nfathsub = father.fGeoEl->NSideSubElements2(father.fSide);
	}
	return father;
}

void TPZGeoElSide::GetSubElements2(TPZStack<TPZGeoElSide> &subelements)
{
	if(!fGeoEl || !fGeoEl->HasSubElement()) {   // Jorge 10/01/2000
		subelements.Resize(0);
		return;
	}
	fGeoEl->GetSubElements2(fSide,subelements);
}

void TPZGeoElSide::HigherLevelCompElementList2(TPZStack<TPZCompElSide> &elvec, int onlyinterpolated, int removeduplicates) {

  if(!Dimension()) return;
  TPZGeoElSide neighbour(*this);
  TPZStack<TPZGeoElSide> subel;
  do {
	  if(neighbour.HasSubElement() && neighbour.NSubElements2() > 1) {
		neighbour.GetSubElements2(subel);
		int nsub = subel.NElements();
		int is;
		for(is=0; is<nsub; is++) {
			subel[is].EqualorHigherCompElementList2(elvec,onlyinterpolated,removeduplicates);
		}
		neighbour = *this;
	  } else {
		  neighbour = neighbour.Neighbour();
	  }
	  if(!neighbour.Exists()) break;
  } while(neighbour != *this);
}

void TPZGeoElSide::EqualorHigherCompElementList2(TPZStack<TPZCompElSide> &celside, int onlyinterpolated, int removeduplicates){


	int ncelsides = celside.NElements();
	if(Reference().Exists()) {
		celside.Push(Reference());
		if(removeduplicates) {
			return;
		}
	}
	this->EqualLevelCompElementList(celside,onlyinterpolated,removeduplicates);
	if(ncelsides != celside.NElements()) return;
	TPZStack<TPZGeoElSide> gelsides;
	TPZGeoElSide neighbour(*this);
	do {
		if(neighbour.HasSubElement() && neighbour.NSubElements2() > 1) {
			neighbour.GetSubElements2(gelsides);
			int nsub = gelsides.NElements();
			int is;
			for(is=0; is<nsub; is++) {
				gelsides[is].EqualorHigherCompElementList2(celside,onlyinterpolated,removeduplicates);
			}
			neighbour = *this;
		} else {
			neighbour = neighbour.Neighbour();
		}
		if(!neighbour.Exists()) break;
	} while(neighbour != *this);

}


int TPZGeoElSide::NSubElements2()
{

	if(!Exists()) return -1;
	return fGeoEl->NSideSubElements2(fSide);
}

void TPZGeoElSide::BuildConnectivities(TPZVec<TPZGeoElSide> &sidevec,TPZVec<TPZGeoElSide> &neighvec){
/**
   os vetores trazem a particão do lado comum a  
   dois vizinhos segundo os seus proprios padrões de
   refinamento, a divisão é identica para este lado comum*///cout << "Sao iguais: acertar as vizinhancas!!!\n";
  int size = sidevec.NElements();
  int neighsize = neighvec.NElements();
  if(size!=neighsize || !size){
    PZError << "TPZGeoElSide::BuildConnectivities wrong vectors: abort!!!\n";
    exit(-1);
  }
  int iv,ivn,side,neighside,sidedim,neighsidedim;
  TPZGeoElSide subside,neighsubside;
  for(iv=0;iv<size;iv++){
    subside = sidevec[iv];
    side = subside.Side();
    sidedim = subside.Dimension();
	TPZGeoEl *elside = subside.Element();
    for(ivn=0;ivn<neighsize;ivn++){
      neighsubside = neighvec[ivn];
      neighside = neighsubside.Side();
      neighsidedim = neighsubside.Dimension();
	  TPZGeoEl *elneigh = neighsubside.Element();
      if(neighsidedim != sidedim) continue;
      //if(temp.Neighbour().Element()) continue;//?????? MELHORAR ESTA LINHA: pode ser viz. do irmão
      int in[4],face,i,j,num;
      int im[4],neighface;      
      switch(sidedim){
      case 0://canto	    
		if(elside->SideNodeIndex(side,0) == elneigh->SideNodeIndex(neighside,0)){
			if(subside.NeighbourExists(neighsubside)) {
				cout << "TPZGeoElSide::BuildConnectivities the neighbour already exists?";
			} else {
				subside.SetConnectivity(neighsubside);
			}
		}
		break;
      case 1://aresta
		in[0] = elside->SideNodeIndex(side,0);
		in[1] = elside->SideNodeIndex(side,1);
		im[0] = elneigh->SideNodeIndex(neighside,0);
		im[1] = elneigh->SideNodeIndex(neighside,1);
		if( (in[0] == im[0] && in[1] == im[1]) || (in[0] == im[1] && in[1] == im[0]) ){
			if(subside.NeighbourExists(neighsubside)) {
				cout << "TPZGeoElSide::BuildConnectivities the neighbour already exists?";
			} else {
				subside.SetConnectivity(neighsubside);
			}
		}
		break;
      case 2://face
		//face = NSideNodes(side);//original substituida pela seguinte
		face = elside->NSideNodes(side);
		neighface = elneigh->NSideNodes(neighside);
		if(face!=neighface) break;
		for(i=0;i<face;i++){
			in[i] = elside->SideNodeIndex(side,i);
			im[i] = elneigh->SideNodeIndex(neighside,i);
		}
		if(face==3) in[3] = im[3] = -1;
		num = 0;
		for(i=0;i<4;i++) for(j=0;j<4;j++) if(in[i]==im[j]) num++;
		if(num==4){
			if(subside.NeighbourExists(neighsubside)) {
				cout << "TPZGeoElSide::BuildConnectivities the neighbour already exists?";
			} else {
				subside.SetConnectivity(neighsubside);
			}
		}
		break;
      default:
		PZError << "TPZGeoElSide::BuildConnectivities error !!!\n";
	
      }
    }
  }
}
