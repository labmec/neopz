/**
 * @file
 * @brief Contains the implementation of the TPZGeoElSide methods.
 */

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

using namespace pzshape;
using namespace std;

#include "pzlog.h"

#include <sstream>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.tpzgeoelside"));
#endif

TPZGeoElSide::TPZGeoElSide(TPZGeoEl *gel, std::set<int> &sideCornerNodes)
{
	fGeoEl = 0; fSide = -1;
	
	std::set<int> nodes;
	for(int s = 0; s < gel->NSides(); s++)
	{
		nodes.clear();
		TPZGeoElSide actualSide(gel, s);
		if(actualSide.NSideNodes() != (int)sideCornerNodes.size())
		{
			continue;
		}
		for(int as = 0; as < actualSide.NSideNodes(); as++)
		{
			nodes.insert(actualSide.SideNodeIndex(as));
		}
		if(nodes == sideCornerNodes)
		{
			fGeoEl = gel;
			fSide = s;
			
			return;
		}
	}
}

bool TPZGeoElSide::IsAncestor(TPZGeoElSide other){
	if(*this == other) return true;
	TPZGeoElSide father = this->Father2();
	if(father.Element()){
		if(father.Element() == other.Element()){
			return true;
		}
		else{
			if(father.IsAncestor(other)){
				return true;
			}
		}
	}
	return false;
}

bool TPZGeoElSide::IsRelative(TPZGeoElSide other){
	if( this->IsAncestor(other) ) return true;
	if( other.IsAncestor(*this) ) return true;
	return false;
}

void TPZGeoElSide::X(TPZVec< REAL > &loc, TPZVec< REAL > &result) {
	
	TPZVec< REAL > locElement(fGeoEl->Dimension(), 0.);
	
	TPZTransform ElementDim = fGeoEl->SideToSideTransform(fSide, fGeoEl->NSides()-1);
	
	ElementDim.Apply(loc, locElement);
	
	fGeoEl->X(locElement, result);
}

void TPZGeoElSide::Jacobian(TPZVec<REAL> &param,TPZFMatrix &jacobian,TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv) {
	
	if(!fGeoEl) return;
	int DIM = fGeoEl->Dimension();
	
	TPZTransform ThisTransf = fGeoEl->SideToSideTransform(fSide, fGeoEl->NSides()-1);
	TPZManVector< REAL,3 > paramElement(DIM,0.);
	
	ThisTransf.Apply(param,paramElement);
	REAL DetElement;
	
	TPZFNMatrix<9> JCn(3,DIM,0.), JacElement(DIM,DIM,0.), AxesElement(DIM,3,0.), Temp(3,DIM,0.), InvElement(DIM,DIM,0.), axest(3,DIM);
	
	fGeoEl->Jacobian(paramElement,JacElement,AxesElement,DetElement,InvElement);
	
	AxesElement.Transpose();
	AxesElement.Multiply(JacElement,JCn);
	JCn.Multiply(ThisTransf.Mult(),Temp);
	
	Temp.GramSchmidt(axest,jacobian); 
	axest.Transpose(&axes);
	jacinv.Resize(jacobian.Rows(),jacobian.Cols());
	if(axes.Rows() == 1) 
	{ 
		detjac = jacobian(0,0); 
		if(detjac)
		{
			jacinv(0,0) = 1./jacobian(0,0); 
		}
		else
		{
			jacinv(0,0) = 0.;
		}
	}
	if(axes.Rows() == 2) {
		detjac = jacobian(0,0)*jacobian(1,1) - jacobian(1,0)*jacobian(0,1);
		if(detjac)
		{
			jacinv(0,0) =  jacobian(1,1)/detjac; jacinv(0,1) = -jacobian(0,1)/detjac;
			jacinv(1,0) = -jacobian(1,0)/detjac; jacinv(1,1) =  jacobian(0,0)/detjac;
		}
		else
		{
			jacinv.Zero();
		}
	}
	if(axes.Rows() == 3) {
		detjac = -(jacobian(0,2)*jacobian(1,1)*jacobian(2,0)) + jacobian(0,1)*jacobian(1,2)*jacobian(2,0) + jacobian(0,2)*jacobian(1,0)*jacobian(2,1) - jacobian(0,0)*jacobian(1,2)*jacobian(2,1) - jacobian(0,1)*jacobian(1,0)*jacobian(2,2) + jacobian(0,0)*jacobian(1,1)*jacobian(2,2);
		jacinv(0,0) = (-(jacobian(1,2)*jacobian(2,1)) + jacobian(1,1)*jacobian(2,2))/detjac;
		jacinv(0,1) = (jacobian(0,2)*jacobian(2,1) - jacobian(0,1)*jacobian(2,2))/detjac;
		jacinv(0,2) = (-(jacobian(0,2)*jacobian(1,1)) + jacobian(0,1)*jacobian(1,2))/detjac;
		jacinv(1,0) = (jacobian(1,2)*jacobian(2,0) - jacobian(1,0)*jacobian(2,2))/detjac;
		jacinv(1,1) = (-(jacobian(0,2)*jacobian(2,0)) + jacobian(0,0)*jacobian(2,2))/detjac;
		jacinv(1,2) = (jacobian(0,2)*jacobian(1,0) - jacobian(0,0)*jacobian(1,2))/detjac;
		jacinv(2,0) = (-(jacobian(1,1)*jacobian(2,0)) + jacobian(1,0)*jacobian(2,1))/detjac;
		jacinv(2,1) = (jacobian(0,1)*jacobian(2,0) - jacobian(0,0)*jacobian(2,1))/detjac;
		jacinv(2,2) = (-(jacobian(0,1)*jacobian(1,0)) + jacobian(0,0)*jacobian(1,1))/detjac;
	}
}

/// Returns the number of sides in which the current side can be decomposed
int TPZGeoElSide::NSides()
{
    TPZStack<int> lower;
    fGeoEl->LowerDimensionSides(fSide,lower);
    return lower.NElements()+1;
}


/// Area associated with the side
REAL TPZGeoElSide::Area()
{
	TPZManVector<REAL,3> elparam(fGeoEl->Dimension(),0.), sideparam(Dimension(),0);
    TPZTransform tr = fGeoEl->SideToSideTransform(fGeoEl->NSides()-1, fSide);
	REAL detjac;
	TPZFNMatrix<9> jacinv(3,3),jacobian(3,3),axes(3,3);
    //supondo jacobiano constante: X linear
	CenterPoint(elparam);
    tr.Apply(elparam, sideparam);
	Jacobian(sideparam,jacobian,axes,detjac,jacinv);
    TPZIntPoints *intrule = fGeoEl->CreateSideIntegrationRule(fSide, 0);
    REAL RefElVolume = 0.;
    int np = intrule->NPoints();
    TPZManVector<REAL,3> points(Dimension());
    REAL weight;
    for (int ip=0; ip<np; ip++) {
        intrule->Point(ip, points, weight);
        RefElVolume += weight;
    }
	return (RefElVolume*detjac);//RefElVolume(): volume do elemento de refer�ncia
	
}


int TPZGeoElSide::NNeighbours()
{
	int nneighbours = 0;
	TPZGeoElSide neigh = this->Neighbour();
	
	while(neigh.Element() != this->Element())
	{
		nneighbours++;
		neigh = neigh.Neighbour();
	}
	
	return nneighbours;
}


int TPZGeoElSide::NNeighboursButThisElem(TPZGeoEl *thisElem)
{
	int nneighbours = 0;
	TPZGeoElSide neigh = this->Neighbour();
	
	while(neigh.Element() != this->Element())
	{
		if(neigh.Element() != thisElem)
		{
			nneighbours++;	
		}
		neigh = neigh.Neighbour();
	}
	
	return nneighbours;
}

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
		if (neighpre.Exists()) neighpre.SetNeighbour(neighpos);
	} else {
		PZError << "TPZGeoElSide::SetConnectivity neighbourhood cycle error";
	}
}

using namespace std;

/**
 * This method inserts the element/side and all lowerdimension sides into the connectivity loop
 */
void TPZGeoElSide::InsertConnectivity(TPZGeoElSide &neighbour)
{
	if(!fGeoEl || !neighbour.fGeoEl) return;
	
	TPZStack<int> mylowerdimension, neighbourlowerdimension;
	fGeoEl->LowerDimensionSides(fSide,mylowerdimension);
	neighbour.fGeoEl->LowerDimensionSides(neighbour.fSide,neighbourlowerdimension);
	SetConnectivity(neighbour);
	int ns = mylowerdimension.NElements();
	int is;
	for(is=0; is<ns; is++)
	{
		TPZGeoElSide myl(fGeoEl,mylowerdimension[is]);
		TPZGeoElSide neighl(neighbour.fGeoEl,neighbourlowerdimension[is]);
		myl.SetConnectivity(neighl);
	}
	
}

void TPZGeoElSide::SetConnectivity(const TPZGeoElSide &neighbour) const{
	
	if(!Exists()) return;
	if(fSide >= fGeoEl->NSides()) {
		PZError << "ERROR(TPZGeoElSide::SetConnectivity)-> Index greater than number of sides.\n";
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
		int a = NeighbourExists(neighbour);
		int b = neighbour.NeighbourExists(*this);
		if(
		   (a && !b) ||
		   (!a && b)
		   )
		{
			std::cout << "This element side : " << fSide << std::endl;
			this->Element()->Print(std::cout);
			cout << "\nNeighbour side :"   << neighbour.Side() << std::endl;
			neighbour.Element()->Print(std::cout);
			PZError << "TPZGeoElSide::SetConnectivity Fourth untreated case, wrong data structure\n";
		} else if(!a) {
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

void TPZGeoElSide::CenterPoint(TPZVec<REAL> &center)
{
	if(!fGeoEl) return;
	fGeoEl->CenterPoint(fSide,center);
}

/*
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
 */

void TPZGeoElSide::ComputeNeighbours(TPZStack<TPZGeoElSide> &compneigh) {
	if(fSide < fGeoEl->NCornerNodes())
    {
		AllNeighbours(compneigh);
		return;
    }
	int nsnodes = NSideNodes();
	TPZStack<TPZGeoElSide> GeoElSideSet;
	TPZStack<int> GeoElSet[27];
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
			GeoElSet[in].Push(GeoElSideSet[el].Element()->Index());
		}
		Sort<int>(GeoElSet[in]);
	}
	TPZStack<int> result;
	switch(nsnodes) {
		case 1:
		{
			result = GeoElSet[0];
		}
			break;
		case 2:
			Intersect<int,DEFAULTVEC_ALLOC>(GeoElSet[0],GeoElSet[1],result);
			break;
		case 3:
			Intersect<int,DEFAULTVEC_ALLOC>(GeoElSet[0],GeoElSet[1],GeoElSet[2],result);
			break;
		case 4:
		{
			TPZStack<int> inter1, inter2;
			Intersect<int,DEFAULTVEC_ALLOC>(GeoElSet[0],GeoElSet[2],inter1);
			if(inter1.NElements()==0) break;
			Intersect<int,DEFAULTVEC_ALLOC>(GeoElSet[1],GeoElSet[3],inter2);
			if(inter2.NElements()==0) break;
			Intersect<int,DEFAULTVEC_ALLOC>(inter1,inter2,result);
		}
			break;
		default:
		{
			TPZStack<int> inter1, inter2;
			inter1 = GeoElSet[0];
			for(in=0; in<nsnodes-1; in++) {
				inter2.Resize(0);
				Intersect<int,DEFAULTVEC_ALLOC>(inter1,GeoElSet[in+1],inter2);
				if(inter2.NElements() == 0) break;
				inter1 = inter2;
			}
			result = inter2;
		}
	}
	int el,nel = result.NElements();
	TPZGeoMesh * geoMesh = fGeoEl->Mesh();
	for(el=0; el<nel; el++) {
		TPZGeoEl * gelResult = geoMesh->ElementVec()[result[el]];
		int whichSd = gelResult->WhichSide(nodeindexes);
		if(whichSd > 0)
		{
			compneigh.Push(TPZGeoElSide( gelResult, whichSd));
		}
	}
}



TPZTransform TPZGeoElSide::NeighbourSideTransform(TPZGeoElSide &neighbour) {
	
#ifdef DEBUG
	if(!NeighbourExists(neighbour))
	{
		stringstream sout;
		sout << __PRETTY_FUNCTION__ << "Neighbour does not exist : expect trouble";
		LOGPZ_ERROR(logger,sout.str());
		TPZTransform toto;
		return toto;
	}
#endif
	int sidedimension = Dimension();
	TPZTransform tside(sidedimension);//transforma�o local
	switch (sidedimension) {
		case 0://canto para canto viz
			
			break;
			
		case 1://aresta para aresta viz
			if (SideNodeIndex(0) == neighbour.SideNodeIndex(0))
			{
				tside.Mult()(0,0) =  1.;
			}
			else
			{
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
				//      } else {//tri�gulos
				//         for(i=0;i<3;i++) idfrom[i] = Element()->SideNodeIndex(Side(),i);
				//      }
				for(i=0;i<3;i++) idto[i] = neighbour.Element()->SideNodeIndex(neighbour.Side(),i);
				//      if(neighbour.Element()->NSides() > 9) {//elementos 3D
				//         neighbour.Element()->NodeFaceIds(idto,neighbour.Side());
				//      } else {//tri�gulos
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
				PZError << "TPZGeoElSide::NeighbourSideTransform : elemento desconhecido" << std::endl;
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

TPZCompElSide TPZGeoElSide::Reference() const {
	if (!fGeoEl) return TPZCompElSide();
	TPZCompElSide compside(fGeoEl->Reference(),fSide);
	return compside;
}

int TPZGeoElSide::Dimension() const {
	if (!fGeoEl) {
		PZError << "TPZGeoElSide::Dimension : null element\n";
		return -1;
	}
	return fGeoEl->SideDimension(fSide);
}

void TPZGeoElSide::SideTransform3(TPZGeoElSide neighbour,TPZTransform &t)	{
	//t : atual -> neighbour
	TPZGeoElSide father(*this);
	if(!father.Exists()) {
		PZError << "TPZGeoElSide::SideTransform3 I dont understand\n";
		return;
	}
	while(father.Exists())
	{
		if(father.NeighbourExists(neighbour)) {
			TPZTransform Temp =  father.NeighbourSideTransform(neighbour);
			//       t =  NeighbourSideTransform(neighbour).Multiply(t);
			t =  Temp.Multiply(t);
			return;
		} else {
			TPZGeoElSide nextfather = father.StrictFather();
			if(nextfather.Exists())
			{
				t = father.Element()->BuildTransform2(father.Side(),nextfather.Element(),t);
			}
			father = nextfather;
		}
	}  
	TPZGeoElSide start,neighbourwithfather,neighfather;
	start = *this;
	neighbourwithfather = *this;
	do {
		neighfather = neighbourwithfather.StrictFather();
		if(!neighfather.Exists()) neighbourwithfather = neighbourwithfather.Neighbour();
	} while(!neighfather.Exists() && neighbourwithfather.Exists() && neighbourwithfather != start);
	int secondcase = 0;
	
	//  TPZGeoElSide nextfather = father.Father2();
	if(neighfather.Exists()) {
		if(neighbourwithfather != start) 
		{
			secondcase++;
			t = start.NeighbourSideTransform(neighbourwithfather).Multiply(t);
		}
		neighbourwithfather.SideTransform3(neighbour,t);
		return;
		
		//     //father.Element()->BuildTransform(Side(),nextfather.Element(),t);
		//     t = neighbourwithfather.Element()->BuildTransform2(neighbourwithfather.Side(),neighfather.Element(),t);//Cedric 01/10/99 e  e 30/04/00
		//     start = neighfather;
		//     if(start.NeighbourExists(neighbour)) {
		//       t = start.NeighbourSideTransform(neighbour).Multiply(t);
		//       //      if(secondcase) {
		//       //	cout << "TPZGeoElSide:SideTranform3 secondcase = " << secondcase << "\n";
		//       //      }
		//       return;
		//     }
		//     neighbourwithfather = start;
		//     do {
		//       neighfather = neighbourwithfather.Father2();
		//       if(!neighfather.Exists()) neighbourwithfather = neighbourwithfather.Neighbour();
		//     } while(!neighfather.Exists() && neighbourwithfather.Exists() && neighbourwithfather != start);
		//     //    nextfather = father.Father2();
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
		if(ref.Element() && ref.Element() != Reference().Element() && (!onlyinterpolated || dynamic_cast<TPZInterpolatedElement*>(ref.Element()) )) {
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


TPZTransform TPZGeoElSide::SideToSideTransform(TPZGeoElSide &higherdimensionside) {
    if(fGeoEl != higherdimensionside.fGeoEl) {
		PZError << "TPZGeoElSide::SideToSideTransform inconsistent id1 = " << fGeoEl->Id() << 
		" id2 = " << higherdimensionside.fGeoEl->Id() << std::endl;
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
int TPZGeoElSide::NSideNodes() const {
    if(!fGeoEl) return 0;
    return fGeoEl->NSideNodes(fSide);
}

/**returns the index of the nodenum node of side*/
int TPZGeoElSide::SideNodeIndex(int nodenum) const {
    if(!fGeoEl) return -1;
    return ( fGeoEl->SideNodeIndex(fSide,nodenum) );
}

/**returns the index of the local nodenum  node of side*/
int TPZGeoElSide::SideNodeLocIndex(int nodenum) const {
    if(!fGeoEl) return -1;
    return ( fGeoEl->SideNodeLocIndex(fSide,nodenum) );
}

std::set<int> TPZGeoElSide::SideNodeIndexes()
{
	std::set<int> nodes;
	for(int n = 0; n < this->NSideNodes(); n++)
	{
		nodes.insert(this->SideNodeIndex(n));
	}
	
	return nodes;
}

TPZCompElSide TPZGeoElSide::LowerLevelCompElementList2(int onlyinterpolated)
{
	// This method was modified to look for the father of any neighbouring element
	// It is not suficient to look for the father of the current element only, because a neighbour
	// might have a father where the current element doesn t. This happens in the clone meshes. It probably
	// will happen when working with interface elements or any situation where an element is inserted in an already refined mesh
	TPZGeoElSide father,neighbour,start;
	start = *this;
	neighbour =  Neighbour();
	father = StrictFather();
	int secondcase = 0;
	while(!father.Exists() && neighbour.Exists() && neighbour != start) {
		father = neighbour.StrictFather();
		neighbour = neighbour.Neighbour();
		secondcase++;
	}
	if(!father.Exists()) return TPZCompElSide();
	TPZStack<TPZCompElSide> equal;
	//2003-12-03
	if (father.Reference().Exists())  equal.Push(father.Reference());

	father.EqualLevelCompElementList(equal,onlyinterpolated,1);
	
	
	while(father.Exists() && equal.NElements() == 0) {
		neighbour = father.Neighbour();
		start = father;
		father = father.StrictFather();
		while(!father.Exists() && neighbour.Exists() && neighbour != start) {
			secondcase++;
			father = neighbour.StrictFather();
			neighbour = neighbour.Neighbour();
		}
		father.EqualLevelCompElementList(equal,onlyinterpolated,1);
		//2003-12-03
		if (father.Reference().Exists())  equal.Push(father.Reference());
		
	}
	//  if(equal.NElements() && secondcase) {
	//    cout << "TPZGeoElSide::LowerLevelCompElementList second case " << secondcase << "\n";
	//  }
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
		if(neighbour.HasSubElement() && neighbour.NSubElements() > 1) {
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
		if(neighbour.HasSubElement() && neighbour.NSubElements() > 1) {
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


int TPZGeoElSide::NSubElements()
{
	
	if(!Exists()) return -1;
	return fGeoEl->NSideSubElements2(fSide);
}

void TPZGeoElSide::BuildConnectivities(TPZVec<TPZGeoElSide> &sidevec,TPZVec<TPZGeoElSide> &neighvec){
	/**
	 os vetores trazem a partic� do lado comum a  
	 dois vizinhos segundo os seus proprios padr�s de
	 refinamento, a divis� �identica para este lado comum*/ //cout << "Sao iguais: acertar as vizinhancas!!!\n";
	int size = sidevec.NElements();
	int neighsize = neighvec.NElements();
	if(size!=neighsize || !size){
		PZError << "TPZGeoElSide::BuildConnectivities wrong vectors: abort!!!\n";
		DebugStop();
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
			//if(temp.Neighbour().Element()) continue;//?????? MELHORAR ESTA LINHA: pode ser viz. do irm�
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

std::ostream &operator << (std::ostream & out,const TPZGeoElSide &geoside){
	out << "TPZGeoElSide : side = " << geoside.Side() << std::endl ;
	if (geoside.Element()) geoside.Element()->Print(out);
	return out;
}

bool TPZGeoElSide::IsLinearMapping() const
{
	if(!fGeoEl) return false;
	return fGeoEl->IsLinearMapping();
}
