/**
 * @file
 * @brief Contains the implementation of the TPZGeoEl methods.
 */

#include "pzgeoel.h"
#include "pzcompel.h"
#include "pzgmesh.h"
#include "pzgnode.h"
#include "pzerror.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pztrnsform.h"
#include "pzvec.h"
#include "pzstack.h"
#include "pzquad.h"
#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "TPZRefPattern.h"
#include "pzlog.h"

#include "TPZGeoLinear.h"
#include "pzgeotriangle.h"
#include "pzgeoquad.h"
#include "pzgeotetrahedra.h"
#include "pzgeopyramid.h"
#include "pzgeoprism.h"
#include "TPZGeoCube.h"
#include "pzcmesh.h"


#include <stdio.h>
#include <stdlib.h>

using namespace std;
using namespace pzgeom;

#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.tpzgeoel");
static TPZLogger loggerorient("pz.mesh.tpzgeoel.orient");
#endif


class TPZRefPattern;

TPZFMatrix<REAL> TPZGeoEl::gGlobalAxes;

// Destructor and Constructors
TPZGeoEl::~TPZGeoEl(){
    int64_t index = fIndex;
    if (this->fFatherIndex != -1) {
        if(!this->Father()){
            //Why did this element lose its father?
            DebugStop();
        } else {
            int subelindex = WhichSubel();
            if (subelindex == -1) {
                DebugStop();
            }
            Father()->SetSubElement(subelindex, 0);
        }
    }
    if (index != -1){
#ifdef PZDEBUG
        if(!fMesh) DebugStop();
#endif
        fMesh->ElementVec()[index] = NULL;
        fMesh->ElementVec().SetFree(index);  //the same line in TPZGeoMesh::DeleteElement was commented. Just call this once.
    }
    fIndex = -1;
    fId = -1;
    fMatId = 0;
    fNumInterfaces = 0;
    fReference = 0;
    fMesh = 0;
}


TPZGeoEl::TPZGeoEl(int64_t id,int materialid,TPZGeoMesh &mesh) {
	fMesh = &mesh;
	fId = id;
	mesh.SetElementIdUsed(id);
	fMatId = materialid;
	this->fReference = NULL;
	fFatherIndex = -1;
	int64_t index = fMesh->ElementVec().AllocateNewElement();
	fIndex = index;
	fMesh->ElementVec()[index] = this;
	this->fNumInterfaces = 0;
}

TPZGeoEl::TPZGeoEl(const TPZGeoEl &el):TPZSavable(el){
	this->fMesh = el.fMesh;
	this->fId = fMesh->CreateUniqueElementId();
	this->fMatId = el.fMatId;
	this->fReference = NULL;
	this->fFatherIndex = el.fFatherIndex;
	this->fIndex = fMesh->ElementVec().AllocateNewElement();
	this->fMesh->ElementVec()[fIndex] = this;
	this->fNumInterfaces = 0;
}

TPZGeoEl::TPZGeoEl(int materialid,TPZGeoMesh &mesh, int64_t &index) {
	this->fMesh = &mesh;
    this->fId = fMesh->CreateUniqueElementId();
	this->fMatId = materialid;
	this->fReference = NULL;
	this->fFatherIndex = -1;
	index = this->fMesh->ElementVec().AllocateNewElement();
	this->fIndex = index;
	this->fMesh->ElementVec()[fIndex] = this;
	this->fNumInterfaces = 0;
}

TPZGeoEl::TPZGeoEl(int materialid,TPZGeoMesh &mesh) {
    this->fMesh = &mesh;
    this->fId = fMesh->CreateUniqueElementId();
	this->fMatId = materialid;
	this->fReference = NULL;
	this->fFatherIndex = -1;
	const int64_t index = this->Mesh()->ElementVec().AllocateNewElement();
	this->fIndex = index;
	this->Mesh()->ElementVec()[index] = this;
	this->fNumInterfaces = 0;
}

void TPZGeoEl::Shape1d(REAL x,int num,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi){
	if(num != 2 && num != 3){
		PZError << "elcalc1d.shape, at this point only linear and quadratic elements\n";
		return;
	}
	
	if(num == 2) {
		phi(0,0) = (1-x)/2.;
		phi(1,0) = (1+x)/2.;
		dphi(0,0) = -0.5;
		dphi(0,1) = 0.5;
	} else {
		phi(0,0) = -x*(1.-x)*0.5;
		dphi(0,0) = x-0.5;
		phi(1,0) = (1.-x*x);
		dphi(0,1) = -2.*x;
		phi(2,0) = 0.5*x*(1.+x);
		dphi(0,2) = x+0.5;
	}
}

void TPZGeoEl::ShapePhi1d(REAL x,int num,TPZFMatrix<REAL> &phi) {
	if(num != 2 && num != 3){
		PZError << "TPZGeoEl ShapePhi1d, at this point only linear and quadratic elements\n";
		return;
	}
	
	if(num == 2) {
		phi(0,0) = (1-x)/2.;
		phi(1,0) = (1+x)/2.;
	} else {
		phi(0,0) = -x*(1.-x)*0.5;
		phi(1,0) = (1.-x*x);
		phi(2,0) = 0.5*x*(1.+x);
	}
}

int TPZGeoEl::WhichSide(TPZVec<int64_t> &SideNodeIds) {
	int64_t cap = SideNodeIds.NElements();
	int nums = NSides();
	for(int side=0; side<nums; side++) {
		if(NSideNodes(side)==2 && cap == 2) {
			int64_t isn1 = SideNodeIndex(side,0);
			int64_t isn2 = SideNodeIndex(side,1);//sao = para side<3
			if((isn1 == SideNodeIds[0] && isn2 == SideNodeIds[1]) ||
			   (isn2 == SideNodeIds[0] && isn1 == SideNodeIds[1]))    return side;
		} else if(NSideNodes(side)== 1 && cap ==1) {
			if(SideNodeIndex(side,0) == SideNodeIds[0]) return side;
			//completar
		} else if(NSideNodes(side) == 3 && cap==3) {
			int64_t sni[3],snx[3],k;
			for(k=0;k<3;k++) snx[k] = SideNodeIndex(side,k);//el atual
			for(k=0;k<3;k++) sni[k] = SideNodeIds[k];//el viz
			for(k=0;k<3;k++) {
				if(snx[0]==sni[k] && snx[1]==sni[(k+1)%3] && snx[2]==sni[(k+2)%3]) return side;
				if(snx[0]==sni[k] && snx[1]==sni[(k+2)%3] && snx[2]==sni[(k+1)%3]) return side;
			}//012 120 201 , 021 102 210
		} else if(NSideNodes(side) == 4 && cap == 4) {//face quadrilateral
			int64_t sni[4],snx[4],k;
			for(k=0;k<4;k++) snx[k] = SideNodeIndex(side,k);//el atual
			for(k=0;k<4;k++) sni[k] = SideNodeIds[k];//vizinho
			if(snx[0]==sni[0]) {
				for(k=1;k<4;k++) {
					if(snx[1]==sni[k] && snx[2]==sni[k%3+1]     && snx[3]==sni[(k+1)%3+1]) return side;
					if(snx[1]==sni[k] && snx[2]==sni[(k+1)%3+1] && snx[3]==sni[k%3+1])     return side;
				}// /* 0123 0231 0312 , 0132 0213 0321 */
			} else if(snx[1]==sni[0]) {
				for(k=1;k<4;k++) {
					if(snx[0]==sni[k] && snx[2]==sni[k%3+1]     && snx[3]==sni[(k+1)%3+1]) return side;
					if(snx[0]==sni[k] && snx[2]==sni[(k+1)%3+1] && snx[3]==sni[k%3+1])     return side;
				}// /* 0123 0231 0312 , 0132 0213 0321 */                               /* 1023 1230 1302 , 1032 1203 1320 */
			} else if(snx[2]==sni[0]) {
				for(k=0;k<4;k++) {
					if(snx[0]==sni[k] && snx[1]==sni[k%3+1]     && snx[3]==sni[(k+1)%3+1]) return side;
					if(snx[0]==sni[k] && snx[1]==sni[(k+1)%3+1] && snx[3]==sni[k%3+1])     return side;
				}// /* 0123 0231 0312 , 0132 0213 0321 */                               /* 2013 2130 2301 , 2031 2103 2310 */
			} else if(snx[3]==sni[0]) {
				for(k=0;k<4;k++) {
					if(snx[0]==sni[k] && snx[1]==sni[k%3+1]     && snx[2]==sni[(k+1)%3+1]) return side;
					if(snx[0]==sni[k] && snx[1]==sni[(k+1)%3+1] && snx[2]==sni[k%3+1])     return side;
				}// /* 0123 0231 0312 , 0132 0213 0321 */                               / * 3012 3120 3201 , 3021 3102 3210 * /
			}
		} else if(cap<1 || cap > 4) {
            int is;
            for (is=0; is<nums; is++) {
                if (NSideNodes(is) == cap) {
                    break;
                }
            }
            if (is != nums) {
                PZError << "TPZGeoEl::WhichSide must be extended\n";
                DebugStop();
            }
		}
	}
	return -1;
}

int TPZGeoEl::NeighbourExists(int side,const TPZGeoElSide &gel) {
	TPZGeoElSide thisside(this,side);
	if(gel == thisside) return 1;
	TPZGeoElSide neighbour = Neighbour(side);
	if(!neighbour.Exists()) return 0;
	while(neighbour != thisside) {
		if(gel == neighbour) return 1;
		neighbour = neighbour.Neighbour();
	}
	return 0;
}

void TPZGeoEl::Print(std::ostream & out) {
	
	out << "Element index      " << fIndex << endl;
	out << "Element id         " << fId << endl;
	out << "Is GeoBlend ? : ";
	if(this->IsGeoBlendEl())
	{
		out << "true" << endl;
	}
	else
	{
		out << "false" << endl;
	}
	out << "Is Linear Mapping ? : ";
	if(this->IsLinearMapping())
	{
		out << "true" << endl;
	}
	else
	{
		out << "false" << endl;
	}
	out << "Number of nodes    " << NNodes() << endl;
	out << "Corner nodes       " << NCornerNodes() << endl;
	out << "Nodes indexes          ";
	int i;
	for (i = 0;i < NNodes();i++) out << NodeIndex(i) << " ";
	out << "\nNumber of sides    " << NSides() << endl;
	if (fMatId < 0) out << "boundary condition " << fMatId << endl;
	else out << "Material id        " << fMatId << endl;
	if (!Father()) out << "no father\n";
	else out << "Father index          " << Father()->Index() << endl;
	if (!SubElement(0)) out << "no subelements";
	else {
		out << "Subelement indexes     ";
		for (i = 0;i < NSubElements();i++) {
			if (!SubElement(i) ) continue;
			out << SubElement(i)->Index() << ' ' ;
		}
	}
	out << endl;
	for (i = 0;i < NSides();i++) {
		out << "Neighbours for side   " << i << " : ";
		TPZGeoElSide neighbour = Neighbour(i);
		TPZGeoElSide thisside(this,i);
		if (!neighbour.Exists()) 
        {
            out << "No neighbour\n";
        }
		else {
			while (neighbour != thisside ) {
				out << neighbour.Element()->Index() << "/" << neighbour.Side() << ' ';
				neighbour = neighbour.Neighbour();
			}
			out << endl;
		}
	}
	out << "Reference element pointer address: " << fReference << endl;
	if (this->Reference()) out << "Reference element index : " << this->Reference()->Index() << endl;
	out << "fNumInterfaces = " << fNumInterfaces << endl;
}

void TPZGeoEl::PrintTopologicalInfo(std::ostream & out) {
    int64_t elIndex = this->Index();
    int nnodes = this->NNodes();
    out << "Element index: " << elIndex << "\n";
    
    for(int n = 0; n < nnodes; n++) {
        REAL nodeX = this->NodePtr(n)->Coord(0);
        REAL nodeY = this->NodePtr(n)->Coord(1);
        REAL nodeZ = this->NodePtr(n)->Coord(2);
        
        out << "Node " << n << " : " << nodeX << " , " << nodeY << " , " << nodeZ << "\n";
    }
}

std::ostream &operator<<(std::ostream &out,TPZGeoEl & el) {
	el.Print(out);
	return out;
}

int TPZGeoEl::Level() {
	int level = 0;
	TPZGeoEl *father = Father();
	while(father) {
		level++;
		father = father->Father();
	}
	return level;
}


int TPZGeoEl::GetTransformId2dQ(TPZVec<int> &idfrom,TPZVec<int> &idto) {
	
	if(idfrom[0]==idto[0] && idfrom[1]==idto[1]) return 0;//sentido horario     : 0123
	if(idfrom[0]==idto[0] && idfrom[1]==idto[3]) return 1;//sentido antihorario : 0321
	//0 esta na posicao 3 e 1 esta na posicao 0  i.e. 1,2,3,0
	if(idfrom[0]==idto[3] && idfrom[1]==idto[0]) return 2;//sentido horario     : 1230
	if(idfrom[0]==idto[1] && idfrom[1]==idto[0]) return 3;//sentido antihorario : 1032
	
	if(idfrom[0]==idto[2] && idfrom[1]==idto[3]) return 4;//sentido horario     : 2301
	if(idfrom[0]==idto[2] && idfrom[1]==idto[1]) return 5;//sentido antihorario : 2103
	
	if(idfrom[0]==idto[1] && idfrom[1]==idto[2]) return 6;//sentido horario     : 3012
	if(idfrom[0]==idto[3] && idfrom[1]==idto[2]) return 7;//sentido antihorario : 3210
	
	return 0;
}

int TPZGeoEl::GetTransformId2dT(TPZVec<int> &idfrom,TPZVec<int> &idto) {
	//REVISAR
	if(idto[0]==idfrom[0] && idto[1]==idfrom[1]) return 0;//sentido horario
	if(idto[0]==idfrom[0] && idto[1]==idfrom[2]) return 1;//sentido antihorario
	
	if(idto[0]==idfrom[1] && idto[1]==idfrom[2]) return 2;//sentido horario
	if(idto[0]==idfrom[1] && idto[1]==idfrom[0]) return 3;//sentido antihorario
	
	if(idto[0]==idfrom[2] && idto[1]==idfrom[0]) return 4;//sentido horario
	if(idto[0]==idfrom[2] && idto[1]==idfrom[1]) return 5;//sentido antihorario
	
	return 0;
}

int TPZGeoEl::ElementExists(TPZGeoEl *elem,int64_t id) {
	
	if(elem == 0 && id == 0) {
		PZError << "\nTPZGeoEl::ElementExists :  element id or element pointer will not be null\n";
		return -1;
	}
	int64_t nelg = fMesh->ElementVec().NElements(),index;
	for(index=0;index<nelg;index++) {
		TPZGeoEl *gel = fMesh->ElementVec()[index];
		if(!gel) continue;
		if(gel == elem || gel->Id() == id) return index;
	}
	return -1;
}

TPZGeoElSide TPZGeoEl::Father2(int /*side*/) const{//Augusto:09/01/01
	PZError << "TPZGeoEl::Father2 should never be called\n";
	return TPZGeoElSide();
}

int TPZGeoEl::FatherSide(int side, int son){
	PZError << "TPZGeoEl::FatherSide should never be called\n";
	return -1;
}

TPZTransform<> TPZGeoEl::BuildTransform2(int /*side*/, TPZGeoEl * /*father*/, TPZTransform<> & /* tr */){//Augusto:09/01/01
	PZError << "TPZGeoEl::BuildTransform2 should never be called\n";
	return TPZTransform<>(0,0);
}


void TPZGeoEl::GetSubElements2(int /*side*/, TPZStack<TPZGeoElSide> &/*subel*/) const
{//Augusto:09/01/01
	PZError << "TPZGeoEl::GetSubElements2 should never be called\n";
}

void TPZGeoEl::GetSubElements2(int side, TPZStack<TPZGeoElSide> &subel, int dimension) const
{
	TPZStack<TPZGeoElSide> subel2;
	GetSubElements2(side,subel2);
	int64_t cap = subel2.NElements();
	int64_t s;
	for(s=0; s<cap; s++) {
		if(subel2[s].Dimension() == dimension) {
			subel.Push(subel2[s]);
		}
	}
}

void TPZGeoEl::GetAllSiblings(TPZStack<TPZGeoEl*> &unrefinedSons)
{
#ifdef PZDEBUG
	PZError << __PRETTY_FUNCTION__ << "is deprecated. Use TPZGeoEl::YoungestChildren instead \n";
	DebugStop();
#endif // PZDEBUG
    int nsons = this->NSubElements();
    for(int s = 0; s < nsons; s++)
    {
        TPZGeoEl * son = this->SubElement(s);
        if(son->HasSubElement() == false)
        {
            unrefinedSons.Push(son);
        }
        else
        {
            son->GetAllSiblings(unrefinedSons);
        }
    }
}

void TPZGeoEl::YoungestChildren(TPZStack<TPZGeoEl*> &unrefinedSons)
{
    int nsons = this->NSubElements();
    for(int s = 0; s < nsons; s++)
    {
        TPZGeoEl * son = this->SubElement(s);
        if(son->HasSubElement() == false)
        {
            unrefinedSons.Push(son);
        }
        else
        {
            son->YoungestChildren(unrefinedSons);
        }
    }
}

int TPZGeoEl::WhichSubel() const{
	
	if(fFatherIndex == -1) {
		PZError << "TPZGeoEl::WhichSubel called with null element\n";
		return -1;
	}
	int son;
	TPZGeoEl *father = Father();
	int nsub = father->NSubElements();
	for(son=0;son<nsub;son++) if(father->SubElement(son) == this) break;
	if(son > (nsub-1)){
		PZError << "TPZGeoEl::WhichSubel son does not exist\n";
#ifdef PZ_LOG
		{
			std::stringstream sout;
			sout << "Father element\n";
			father->Print(sout);
			sout << "Son element\n";
            TPZGeoEl *foefel = (TPZGeoEl *) this;
			foefel->Print(sout);
            Mesh()->Print(sout);
			LOGPZ_ERROR(logger,sout.str())
		}
#endif
		DebugStop();
		return -1;
	}
	return son;
}

int TPZGeoEl::WhichSide(TPZVec<REAL> &pt){
	int dim = Dimension();
	int nums = NSides();//n�mero de lados
	REAL tol = 1.e-06;//toler�ncia do zero -> x � zero se -o<=x<=+o
	int is;
	for(is=0; is<nums; is++) {
		int sdim = SideDimension(is);
		TPZTransform<> t1 = SideToSideTransform(nums-1,is);
		TPZTransform<> t2 = SideToSideTransform(is,nums-1);
		TPZVec<REAL> pts(sdim),pt2(dim);
		t1.Apply(pt,pts);
		t2.Apply(pts,pt2);
		REAL dif=0;
		int d;
		for(d=0; d<dim; d++) {
			dif += (pt[d]-pt2[d])*(pt[d]-pt2[d]);
		}
		dif = sqrt(dif/dim);
		if(dif < tol) return is;
	}
	cout << "TPZGeoEl::WhichSide ERROR : side not found" << endl ;
	
	return -1;
}

void TPZGeoEl::CheckSubelDataStructure(){
	
	TPZVec<TPZGeoEl *> sub;
	Divide(sub);
	int side,nsubside,isub,ok;            //no m�ximo 9 subelementos/lados
	TPZStack<TPZGeoElSide> subside;//4 faces+4 arestas+1 canto
	ofstream out("Checksubdata.dat");
	for(side=0;side<NSides();side++){//percorre-se os lados do elemento atual = pai
		GetSubElements2(side,subside);
		nsubside = subside.NElements();
		for(isub=0;isub<nsubside;isub++){
			TPZGeoElSide neigh = subside[isub];//1o o subelemento depois os seus vizinhos
			while(neigh.Element()){//existe vizinho!
				TPZGeoElSide fatside = neigh.Element()->Father2(neigh.Side());
				if(!fatside.Element()){
					cout << "Inconsistencia de dados : \n";
					cout << "Pai atual/lado : " << Id() << "/" << side << endl;
					cout << "Pai do sub elemento e' nulo!\n";
					cout << "Sub elemento/lado : " << neigh.Element()->Id() << "/" << neigh.Side() << endl;
					cin >> ok;
				} else //o pai existe � igual ao atual e a dimens�o dos lados � a mesma OK!
					if(fatside.Element() == this){
						cout << "Pai atual/lado        : " << Id() << "/" << side << endl;
						cout << "Sub elemento/lado     : " << neigh.Element()->Id() << "/" << neigh.Side() << endl;
						cout << "Pai sub elemento/lado : " << fatside.Element()->Id() << "/" << fatside.Side() << endl;
						out << "Pai atual/lado        : " << Id() << "/" << side << endl;
						out << "Sub elemento/lado     : " << neigh.Element()->Id() << "/" << neigh.Side() << endl;
						out << "Pai sub elemento/lado : " << fatside.Element()->Id() << "/" << fatside.Side() << endl;
						if(fatside.Side()!=side){
							cout << "Dados acima inconsistentes : lados distintos\n";
							out << "Dados acima inconsistentes : lados distintos\n";
							cin >> ok;
						}
					} else {
						cout << "OK!: Vizinho com pai diferente\n";
					}
				neigh = neigh.Neighbour();
				if(!neigh.Element()){
					cout << "Vizinho nulo\n\n";
					out << "Vizinho nulo\n\n";
					break;
				}
				cout << "\n";
				out << "\n";
				if(neigh.Element()==subside[isub].Element()) break;
				//no pr�ximo while confia-se que o ciclo de conectividades � fechado e os viz sao filhos de alguem
				TPZGeoEl *gel = neigh.Element()->Father2(neigh.Element()->NSides()-1).Element();//pai pelo interior:sempre existe
				while(gel && !(gel== this)) neigh = neigh.Neighbour();
				if(neigh.Element()==subside[isub].Element()){
					cout << "Nao existe vizinho irmao!\n";
					cout << "Sub elemento/lado " << neigh.Element()->Id() << "/" << neigh.Side() << endl;
					cout << "OK!\n";
				}
			}
		}
	}
	out.flush();
	out.close();
}

void TPZGeoEl::SetSubElementConnectivities() {
	
	int side;
	for(side=0; side<NCornerNodes(); side++)
	{
		TPZGeoElSide thisside(this,side);
		TPZStack<TPZGeoElSide> subel;
		this->GetSubElements2(side,subel);
		if(!subel[0].NeighbourExists(thisside))
		{
			subel[0].SetConnectivity(thisside);
		}
	}
	for(side=NCornerNodes(); side<NSides(); side++)
	{
		TPZGeoElSide thisside(this,side);
		TPZGeoElSide neighbour = this->Neighbour(side);
		while(neighbour.Exists() && neighbour != thisside)
		{
			if(neighbour.HasSubElement() && neighbour.NSubElements() != 1) {
				TPZStack<TPZGeoElSide> elvec,neighvec;
				GetSubElements2(side,elvec);
				neighbour.GetSubElements2(neighvec);
				
				// The currently divided element may have only one element as a son. In this case, the son is already
				// neighbour of the father
				
				if(elvec.NElements() > 1) TPZGeoElSide::BuildConnectivities(elvec,neighvec);
				
				break;
			}
			neighbour = neighbour.Neighbour();
		}
	}
	int nsubel = NSubElements();
	for (int isub=0; isub<nsubel; isub++) 
	{
		this->SubElement(isub)->InitializeNeighbours();
	}
}

REAL TPZGeoEl::CharacteristicSize()
{
	int nn = NNodes();
	if(!nn)
    {
        return 0.;
    }
    TPZVec<REAL> xmin(3),xmax(3);
    TPZManVector<REAL,3> values(3);
    
	NodePtr(0)->GetCoordinates(values);
	for(int c = 0; c < 3; c++)
    {
        xmin[c] = values[c];
        xmax[c] = values[c];
    }
	for(int n = 1; n < nn; n++)
	{
		NodePtr(n)->GetCoordinates(values);
		for(int c = 0; c < 3; c++)
		{
			xmin[c] = Min(values[c],xmin[c]);
			xmax[c] = Max(values[c],xmax[c]);
		}
	}
    REAL diagVecNorm = 0.;
    for(int c = 0; c < 3; c++)
    {
        diagVecNorm += (xmax[c]-xmin[c])*(xmax[c]-xmin[c]);
    }
    diagVecNorm = sqrt(diagVecNorm);
    
	return diagVecNorm;
}

REAL TPZGeoEl::SmallerEdge()
{
    REAL norm = 0.;
    
    int firstEdge = this->NNodes();
    TPZGeoElSide edge0(this,firstEdge);
    
    TPZVec<REAL> coords0(3), coords1(3);
    int64_t node0id = edge0.SideNodeIndex(0);
    int64_t node1id = edge0.SideNodeIndex(1);
    this->Mesh()->NodeVec()[node0id].GetCoordinates(coords0);
    this->Mesh()->NodeVec()[node1id].GetCoordinates(coords1);
    for(int c = 0; c < 3; c++)
    {
        REAL delta = coords1[c] - coords0[c];
        norm += delta*delta;
    }
    norm = sqrt(norm);
    
    for(int s = firstEdge+1; s < this->NSides(); s++)
    {
        TPZGeoElSide edgeOther(this,s);
        if(edgeOther.Dimension() > 1)
        {
            break;
        }
        
        /////////////
        node0id = edgeOther.SideNodeIndex(0);
        node1id = edgeOther.SideNodeIndex(1);
        this->Mesh()->NodeVec()[node0id].GetCoordinates(coords0);
        this->Mesh()->NodeVec()[node1id].GetCoordinates(coords1);
        REAL normTemp = 0.;
        for(int c = 0; c < 3; c++)
        {
            REAL delta = coords1[c] - coords0[c];
            normTemp += delta*delta;
        }
        normTemp = sqrt(normTemp);
        
        norm = Min(norm,normTemp);
        /////////////
    }
    
    return norm;
}

bool TPZGeoEl::ComputeXInverse(TPZVec<REAL> &XD, TPZVec<REAL> &qsi, REAL Tol) {
	REAL error = 10.;
	int iter = 0;
	const int nMaxIter = 10000;
	REAL radius = CharacteristicSize();
	int dim = Dimension();
	TPZManVector<REAL,3> X0(3);
//	ZeroTolerance();
	// First verify if the entry qsi yields the right point
	if(qsi.NElements()!= dim)
	{
		PZError << "\nTPZGeoEl::ComputeXInverse vector dimension error\n";
		qsi.Resize(Dimension(),0.);//zero esta em todos os elementos mestres
	}
	X(qsi,X0);//qsi deve ter dimensao do elemento atual
	TPZFNMatrix<9> DelX(3,1);
	int i;
	for(i=0; i<3; i++)
	{
		DelX(i,0) = XD[i]-X0[i];
	}
	error = Norm(DelX)/radius;
	if(error <= Tol)
	{
		return ( this->IsInParametricDomain(qsi) );
	}
	
	// Verify if the point is not a corner node
	TPZManVector<REAL,3> values(3);
	int nn = NCornerNodes();
	if(nn <= 1)
	{
		return false;
	}
	int in;
	for(in=0; in<nn; in++)
	{
		NodePtr(in)->GetCoordinates(values);
		TPZFNMatrix<3> DelX(3,1);
		for(i=0; i<3; i++) DelX(i,0) = XD[i]-values[i];
		error = Norm(DelX)/radius;
		if(error <= Tol)
		{
			TPZVec<REAL> zero(0);
			TPZTransform<> tr = SideToSideTransform(in, NSides()-1);
			tr.Apply(zero, qsi);
			return true;
		}
	}
	
	while(error > Tol && iter < nMaxIter)
	{
		iter++;
		TPZFNMatrix<9> residual(dim,1),delqsi(dim,1);
		REAL detJ;
		TPZFNMatrix<9> J(dim,dim,0.),axes(dim,3,0.),Inv(dim,dim,0.);
		TPZFNMatrix<9> JXt(dim,3,0.),JX(3,dim,0.),JXtJX(dim,dim,0.);
        TPZFNMatrix<9,REAL> gradx(3,dim);
        GradX(qsi, gradx);
		Jacobian(gradx,J,axes,detJ,Inv);
		if(fabs(detJ) < 2.e-10)
		{
			TPZManVector<REAL,3> center(Dimension(),0.);
			CenterPoint(NSides()-1, center);
			cout << "ComputeXInverse found zero Jacobian Index " << this->fIndex << " qsi " << qsi << " detJ " << detJ << std::endl;
			for(int ik = 0; ik < qsi.NElements(); ik++)
            {
                qsi[ik] += Tol*(center[ik]-qsi[ik]);
            }
			residual(0,0) = 1.e12;
		}
		else
		{
			TPZFNMatrix<9> axest;
			axes.Transpose(&axest);
			axest.Resize(3,dim);//casos 1D e 2D onde JX espacial � 1x3 e 2x3 respectivamente
			if(dim==1)
            {
				JX(0,0) = axest(0,0)*J(0,0);
				JX(1,0) = axest(1,0)*J(0,0);
				JX(2,0) = axest(2,0)*J(0,0);
			}
            else
            {
				axest.Multiply(J,JX,0);
			}
			
			JX.Transpose(&JXt);
			JXt.Multiply(JX,JXtJX,0);//JXtJX = JXt*JX;
			JXt.Multiply(DelX,residual);//cout << "\nComputeXInverse: : \n";
			JXtJX.SolveDirect(residual,ELU);//cout << "Atual/dimensao : " << Id() << " / " << Dimension();
			for(i=0; i<dim; i++)
            {
                qsi[i] += residual(i,0);
            }
            
		}
		X(qsi,X0);
		for(i=0; i<3; i++)
        {
            DelX(i,0) = XD[i]-X0[i];
        }
		//A norma sobre coordenada parametrica eh mais objetiva pois os limites do
		//elemento mestre sao mais claros que os do elemento real
		//       error = Norm(DelX);
		error = Norm(DelX)/radius;
	}
	
#ifdef PZDEBUG
	if(iter == nMaxIter)
	{
		std::stringstream sout;
		sout << "Error at " << __PRETTY_FUNCTION__ << " - nMaxIter was reached before tolerance is achieved - ElementId" << this->Id() << std::endl;
        sout << "error = " << error << "\t | \t tol = " << Tol << std::endl;
		PZError << "\n" << sout.str() << "\n";
        Print(std::cout);
        int nnodes = NNodes();
        for (int i=0; i<NNodes(); i++) {
            NodePtr(i)->Print();
        }
        
#ifdef PZ_LOG
		LOGPZ_ERROR(logger,sout.str().c_str());
#endif
	}
#endif
	
	return ( this->IsInParametricDomain(qsi) );
}


void TPZGeoEl::TransformSonToFather(TPZGeoEl *ancestor, TPZVec<REAL> &qsiSon, TPZVec<REAL> &qsiAncestor){
	
	int dim = Dimension();
	if(qsiAncestor.NElements()!= dim)
	{
		PZError << "\nTPZGeoEl::TransformSonToFather vector dimension error\n";
		qsiAncestor.Resize(Dimension(),0.);//zero esta em todos os elementos mestres
	}
	
	TPZVec<REAL> xson;
	X(qsiSon,xson);
	
	TPZGeoEl *father = this;
	while(father != ancestor)
	{
		father = father->Father();
	}
	if(!father)
	{
		DebugStop();
	}
	TPZTransform<> tr(dim);
	tr = BuildTransform2(NSides()-1, father, tr);
	tr.Apply(qsiSon, qsiAncestor);
	REAL Tol;
	ZeroTolerance(Tol);
	father->ComputeXInverse(xson, qsiAncestor,Tol);
}

TPZTransform<> TPZGeoEl::ComputeParamTrans(TPZGeoEl *fat,int fatside, int sideson){
	
	//transformacao do lado de elemento pequeno para elemento grande que o contem
	int dimf = fat->Dimension();
	int dimsf = fat->SideDimension(fatside);
	int dim = Dimension();
	int dimss = SideDimension(sideson);
	if(dimsf < dimss){
		PZError << "\nTPZGeoEl::ComputeParamTrans called with sides error\n";
		DebugStop();
	}
	
	/**para o canto do pai n�o existe transformac�o definida*/
	if(!fat->SideDimension(fatside)) return TPZTransform<>(0,0);
	
	REAL weight;
	TPZFNMatrix<9> jac(dim,dim),axes(3,3,0.);
	TPZFNMatrix<9> jacinv(dim,dim);
	TPZManVector<REAL,3> x(3,0.);
	TPZManVector<REAL,3> intpoint(dimss,0.);
	int tam = (dimss+1);
	TPZFNMatrix<16> hess(tam,tam,0.),grad0(tam,1,0.);
	TPZIntPoints *intrule = CreateSideIntegrationRule(sideson,2);
	TPZManVector<int,3> order(dimss,2);
	intrule->SetOrder(order);
	//integra��o sobre o lado-filho contido no lado-pai
	int ij,ik,indp;
	REAL D2Edaikdaij,D2Edcidaij,D2Edci2;
	D2Edcidaij = 0.;
	D2Edci2 = 0.;
	for(ij=0;ij<dimss;ij++){
		for(ik=ij;ik<dimss;ik++){
			D2Edaikdaij = 0.;
			D2Edcidaij  = 0.;
			for(indp = 0; indp < intrule->NPoints(); ++indp){
				intrule->Point(indp,intpoint,weight);
				//Jacobian(intpoint,jac,axes,detjac,jacinv);/**estes passos n�o s�o precisos pois:*/
				//weight *= fabs(detjac);/**ambas as matrizes, A e b, multiplicam a mesma constante detjac*/
				D2Edaikdaij += intpoint[ik]*intpoint[ij]*weight;
				if(ik==ij) D2Edcidaij += intpoint[ij]*weight;
				if(ij==0 && ik==0) D2Edci2 += weight;
			}
			hess(ij,ik) = 2.*D2Edaikdaij;
			hess(ik,ij) = hess(ij,ik);/**basta repassar sendo ik>ij*/
			if(ik==ij) {
				hess(ij,dimss) = 2.*D2Edcidaij;
				hess(dimss,ij) = hess(ij,dimss);
			}
			if(ij==0 && ik==0) hess(dimss,dimss) = 2.*D2Edci2;
		}
	}// final do integral hess
	//do lado sideson para o elemento atual (filho)
	TPZTransform<> tsidetoson(Dimension());//identidade
	if(dimss<Dimension()) tsidetoson = SideToSideTransform(sideson,NSides()-1);
	TPZTransform<> fatelside = fat->SideToSideTransform(fat->NSides()-1,fatside);
	TPZManVector<REAL,3> sidepoint(Dimension());//dimensao do dominio da transformacao X do filho
	int j;//transf. para o lado do pai
	TPZFNMatrix<9> A(dimsf,dimss,0.),sol(dimsf,1,0.);
	REAL Tol;
	ZeroTolerance(Tol);
	for(int ifat=0;ifat<dimsf;ifat++){//numero de variaveis do pai
		REAL DEdci = 0.;
		for(j=0;j<(dimss+1);j++){
			REAL DEdaij = 0.;
			for(indp = 0; indp < intrule->NPoints(); ++indp){
				intrule->Point(indp,intpoint,weight);
				tsidetoson.Apply(intpoint,sidepoint);//lado do filho para o seu interior: mestre
				X(sidepoint,x);//ponto do mestre do filho para o filho deformado, 3 coordenadas
				TPZVec<REAL> csi(dimf,0.);/**o seguinte passo n�o � preciso dado que fat = elemento mestre*/
				fat->ComputeXInverse(x,csi,Tol);//csi ponto no pai mestre
				TPZVec<REAL> outcsi(dimsf,0.);
				fatelside.Apply(csi,outcsi);
				//        fat->ProjectPointElToSide(fatside,csi,outcsi);
				if(j<dimss) DEdaij += outcsi[ifat]*intpoint[j]*weight;
				if(j==0) DEdci += outcsi[ifat]*weight;
			}
			if(j<dimss) grad0(j,0) = 2.*DEdaij;
			if(j==0) grad0(dimss,0) = 2.*DEdci;
			if(!dimss) grad0(0,0) = DEdci;
		}//final integral gradiente
		//resolu��o do sistema para cada variavel ifat do pai
		if(dimss) hess.SolveDirect(grad0,ELU);
		for(int k=0;k<dimss;k++) A(ifat,k) = grad0(k,0);
		sol(ifat,0) = grad0(dimss,0);
	}//fim sistema ifat
	delete intrule;
	TPZTransform<> t(dimsf,dimss);
	t.SetMatrix(A,sol);
	return t;
}

REAL TPZGeoEl::Distance(TPZVec<REAL> &centel,TPZVec<REAL> &centface){
	
	if(centel.NElements() != 3 || centface.NElements() != 3)
		PZError << "TPZGeoEl::Distance point dimension error\n";
	
	int i;
	REAL distance = 0.;
	for(i=0;i<3;i++)
		distance += (centel[i]-centface[i])*(centel[i]-centface[i]);
	
	return sqrt(distance);
}

REAL TPZGeoEl::ElementRadius(){
	
	switch( this->Dimension() )
	{
		case 0:
			return 0.;
			//	 break;
			
		case 1:
		case 2:{
			TPZManVector<REAL, 3> centel  (this->Dimension(), 0.);
			TPZManVector<REAL, 3> centface(this->Dimension(), 0.);
			TPZManVector<REAL, 3> masscent(3,0.), massface(3, 0.);
			REAL mindist = 1000.;
			REAL dist;
			
			CenterPoint(NSides()-1,centel);
			X(centel,masscent);
			
			int nsides = NSides();
			for(int iside = 0; iside < nsides - 1; iside++){
				CenterPoint(iside, centface);
				X(centface,massface);
				dist = Distance(masscent,massface);
				if(mindist > dist) mindist = dist;
			}
			return mindist;
			//	 break;
		}
			
		case 3:{
			TPZManVector<REAL, 3> centel(3,0.),centface(3,0.);
			TPZManVector<REAL, 3> masscent(3,0.),massface(3,0.);
			REAL mindist = 1000.,dist;
			CenterPoint(NSides()-1,centel);
			X(centel,masscent);
			int nfaces,face,face0;
			if(NSides()==15) {nfaces = 14; face0 = 10;}//tetrahedron
			else if(NSides()==19) {nfaces = 18; face0 = 13;}//pyramid
			else if(NSides()==21) {nfaces = 20; face0 = 15;}//prism
			else if(NSides()==27) {nfaces = 26; face0 = 20;}//hexahedro
			else return 0.;//line, point, triangle, quadrilateral
			for(face=face0;face<nfaces;face++){
				CenterPoint(face,centface);
				X(centface,massface);
				dist = Distance(masscent,massface);
				if(mindist > dist) mindist = dist;
			}
			return mindist;
			//	 break;
		}
			
		default:
			PZError <<  "TPZGeoEl::ElementRadius - Dimension not implemented." << endl;
			return 0.;
			
	}//end of switch
	//   return 0.;
} //end of method

REAL TPZGeoEl::TriangleArea(TPZVec<TPZGeoNode *> &nodes){
	
	if(nodes.NElements() != 3 && nodes.NElements() != 4){
		PZError << "TPZGeoEl::AreaFromTheFaceT argument error size: nodes\n";
		nodes.Resize(0);
		return 0.;
	}
	
	REAL cb0 = nodes[2]->Coord(0) - nodes[1]->Coord(0);
	REAL cb1 = nodes[2]->Coord(1) - nodes[1]->Coord(1);
	REAL cb2 = nodes[2]->Coord(2) - nodes[1]->Coord(2);
	REAL ab0 = nodes[0]->Coord(0) - nodes[1]->Coord(0);
	REAL ab1 = nodes[0]->Coord(1) - nodes[1]->Coord(1);
	REAL ab2 = nodes[0]->Coord(2) - nodes[1]->Coord(2);
	
	//produto vetorial
	REAL coord0 = cb1*ab2-ab1*cb2;
	REAL coord1 = ab0*cb2-cb0*ab2;
	REAL coord2 = cb0*ab1-ab0*cb1;
	//norma da metade do vetor
	return ( 0.5*sqrt(coord0*coord0+coord1*coord1+coord2*coord2) );
}

REAL TPZGeoEl::QuadArea(TPZVec<TPZGeoNode *> &nodes){
	if(nodes.NElements() != 4){
		PZError << "TPZGeoEl::AreaFromTheFaceQ argument error size: nodes\n";
		nodes.Resize(0);
		return 0.;
	}
	
	REAL areat1 = TriangleArea(nodes);
	
	nodes[1] = nodes[0];
	
	nodes[0] = nodes[2];
	nodes[2] = nodes[1];//antigo node0
	nodes[1] = nodes[3];
	
	REAL areat2 = TriangleArea(nodes);
	
	return (areat1+areat2);
}

REAL TPZGeoEl::Volume(){

    int dim = Dimension();
	TPZManVector<REAL,3> param(dim,0.);
	REAL detjac;
	TPZFNMatrix<9> jacinv(dim,dim),jacobian(dim,dim),axes(dim,3), gradx(3,dim);
	//supondo jacobiano constante: X linear
	CenterPoint(NSides()-1,param);
    GradX(param, gradx);
	Jacobian(gradx,jacobian,axes,detjac,jacinv);
	return (RefElVolume()*detjac);//RefElVolume(): volume do elemento de refer�ncia
}

REAL TPZGeoEl::SideArea(int side){
	
	if(side < 0 || side > NSides()-1)
		PZError << "TPZGeoEl::AreaFromTheFace side error, side = " << side << endl;
	
	if(SideDimension(side) != 2)
    {
		TPZGeoElSide gelside(this,side);
        return gelside.Area();
    }
	
	if(SideDimension(side) == 2){
		
		int nsn = NSideNodes(side);
		
		TPZVec<TPZGeoNode *> nodes(nsn);
		int i;
		
		for(i=0;i</*3*/nsn;i++)
			nodes[i] = &Mesh()->NodeVec()[  SideNodeIndex(side,i) ];
		if (nsn==4)
			return ( QuadArea(nodes) );
		else
			return ( TriangleArea(nodes) );
	}
	return 0.;
}

TPZCompEl *TPZGeoEl::CreateBCCompEl(int side,int bc,TPZCompMesh &cmesh) {
	TPZGeoEl *gel = CreateBCGeoEl(side,bc);
	return cmesh.CreateCompEl(gel);
}

void TPZGeoEl::RemoveConnectivities(){
	
	int nsides = NSides(),side;
	for(side=0;side<nsides;side++){
		TPZGeoElSide thisside(this,side);
		TPZGeoElSide neighbour (thisside.Neighbour());
        if(!neighbour) DebugStop();
		thisside.RemoveConnectivity();
		if(neighbour != thisside){
		    TPZGeoElSide neighneigh = neighbour;
		    do{
                if(neighneigh.ResetBlendConnectivity(fIndex)){
                    neighneigh.Element()->SetNeighbourForBlending(neighneigh.Side());
                }
		        neighneigh = neighneigh.Neighbour();
            }   while (neighbour != neighneigh);
        }
	}
}

void TPZGeoEl::InitializeNeighbours(){
	int i,j;
	for (i=0;i<NSides();i++){		
		
		TPZStack <TPZGeoElSide> subel;
		if (HasSubElement()){
			GetSubElements2(i,subel);
			for (j=0;j<subel.NElements();j++){
				TPZGeoEl *el = subel[j].Element();
				el->InitializeNeighbours();
			}
		}
		TPZGeoElSide neighside = Neighbour(i);
		if (!neighside.Element() || neighside.Side() == -1){
			TPZGeoElSide thisside (this,i);
			SetNeighbour(i,thisside);
		}
	}
}

void TPZGeoEl::MidSideNodeIndices(int side,TPZVec<int64_t> &indices) const {
	indices.Resize(1);
	MidSideNodeIndex(side,indices[0]);
	if(indices[0] == -1) indices.Resize(0);
}

/** @brief Compute a decomposition of the gradient of the mapping function, as a rotation matrix (Jacobian) and orthonormal basis (axes)  */
void TPZGeoEl::Jacobian(TPZVec<REAL> &qsi,TPZFMatrix<REAL> &jac,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) const{
    TPZFNMatrix<9,REAL> gradx;
    GradX(qsi, gradx);
#ifdef PZ_LOG
    if(logger.isDebugEnabled())
    {
        std::stringstream sout;
        gradx.Print("gradx",sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    Jacobian(gradx, jac, axes, detjac, jacinv);
}

/** @brief Compute a decomposition of the gradient of the mapping function, as a rotation matrix (Jacobian) and orthonormal basis (axes)  */
void TPZGeoEl::JacobianXYZ(TPZVec<REAL> &qsi,TPZFMatrix<REAL> &jac,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) const{
    TPZFNMatrix<9,REAL> gradx;
    GradX(qsi, gradx);
    JacobianXYZ(gradx, jac, axes, detjac, jacinv);
}


void TPZGeoEl::Jacobian(const TPZFMatrix<REAL> &gradx, TPZFMatrix<REAL> &jac,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv){

    detjac = 0.0;
    int nrows = gradx.Rows();
    int ncols = gradx.Cols();
    int dim   = ncols;
    
    switch (dim) {
        case 0:
        {
            jac.Resize(dim,dim);
            axes.Resize(dim,3);
            jacinv.Resize(dim,dim);
            detjac=1.;
            
            break;
        }
            
        case 1:
        {
            jac.Resize(dim,dim);
            axes.Resize(dim,3);
            jacinv.Resize(dim,dim);
            jac.Zero();
            
            /**  Definitions: v1 -> is the xi_direction of the Gradient */
            TPZManVector<REAL,3> v_1(3,0.);
            
            for (int i = 0; i < nrows; i++) {
                v_1[i]  = gradx.GetVal(i,0);
            }
            
            REAL norm_v_1 = 0.;
            for(int i = 0; i < nrows; i++) {
                norm_v_1 += v_1[i]*v_1[i];
            }
            
            norm_v_1    = sqrt(norm_v_1);
            jac(0,0)    = norm_v_1;
            detjac      = norm_v_1;

            if(IsZero(detjac))
            {
                
#ifdef PZDEBUG
                std::stringstream sout;
                sout << "Singular Jacobian, 1 determinant of jacobian = " << detjac << std::endl;
                LOGPZ_ERROR(logger, sout.str())
                //DebugStop();
#endif
                detjac = ZeroTolerance();
            }
            
            jacinv(0,0) = 1.0/detjac;
            
            for(int i=0; i < 3; i++) {
                axes(0,i) = v_1[i]/norm_v_1;
            }
            
        }
            break;
        case 2:
        {
            
            jac.Resize(dim,dim);
            axes.Resize(dim,3);
            jacinv.Resize(dim,dim);
            jac.Zero();
            
            /**  Definitions: v1 -> is the xi_direction of the Gradient, v2 -> is the eta_direction of the Gradient */
            TPZManVector<REAL,3> v_1(3,0.), v_2(3,0.);
            
            /**  Definitions: v_1_til and v_2_til -> asscoiated orthonormal vectors to v_1 and v_2 */
            TPZManVector<REAL,3> v_1_til(3,0.), v_2_til(3,0.);

            for (int i = 0; i < nrows; i++) {
                v_1[i]  = gradx.GetVal(i,0);
                v_2[i]  = gradx.GetVal(i,1);
            }

            REAL norm_v_1_til = 0.0;
            REAL norm_v_2_til = 0.0;
            REAL v_1_dot_v_2  = 0.0;

            for(int i = 0; i < 3; i++) {
                norm_v_1_til    += v_1[i]*v_1[i];
                v_1_dot_v_2     += v_1[i]*v_2[i];
            }
            norm_v_1_til = sqrt(norm_v_1_til);
            
            for(int i=0 ; i < 3; i++) {
                v_1_til[i]          = v_1[i] / norm_v_1_til; // Normalizing
                v_2_til[i]          = v_2[i] - v_1_dot_v_2 * v_1_til[i] / norm_v_1_til;
                norm_v_2_til   += v_2_til[i]*v_2_til[i];
            }
            norm_v_2_til = sqrt(norm_v_2_til);
            
            
            jac(0,0) = norm_v_1_til;
            jac(0,1) = v_1_dot_v_2/norm_v_1_til;
            jac(1,1) = norm_v_2_til;
            
            detjac = jac(0,0)*jac(1,1)-jac(1,0)*jac(0,1);
            
            jacinv(0,0) = +jac(1,1)/detjac;
            jacinv(1,1) = +jac(0,0)/detjac;
            jacinv(0,1) = -jac(0,1)/detjac;
            jacinv(1,0) = -jac(1,0)/detjac;
            
            if(IsZero(detjac))
            {
                
#ifdef PZDEBUG
                std::stringstream sout;
                sout << "Singular Jacobian, 2 determinant of jacobian = " << detjac << std::endl;
                LOGPZ_ERROR(logger, sout.str())
                DebugStop();
#endif
                detjac = ZeroTolerance();
            }
            
            for(int i=0; i < 3; i++) {
                v_2_til[i] /= norm_v_2_til; // Normalizing
                axes(0,i)  = v_1_til[i];
                axes(1,i)  = v_2_til[i];
            }
            
        }
            break;
        case 3:
        {
            jac.Resize(dim,dim);
            axes.Resize(dim,3);
            jacinv.Resize(dim,dim);
            jac.Zero();
            
            for (int i = 0; i < nrows; i++) {
                jac(i,0)  = gradx.GetVal(i,0);
                jac(i,1)  = gradx.GetVal(i,1);
                jac(i,2)  = gradx.GetVal(i,2);
            }
            
            detjac -= jac(0,2)*jac(1,1)*jac(2,0);//- a02 a11 a20
            detjac += jac(0,1)*jac(1,2)*jac(2,0);//+ a01 a12 a20
            detjac += jac(0,2)*jac(1,0)*jac(2,1);//+ a02 a10 a21
            detjac -= jac(0,0)*jac(1,2)*jac(2,1);//- a00 a12 a21
            detjac -= jac(0,1)*jac(1,0)*jac(2,2);//- a01 a10 a22
            detjac += jac(0,0)*jac(1,1)*jac(2,2);//+ a00 a11 a22
            
            if(IsZero(detjac))
            {
                
#ifdef PZDEBUG
                std::stringstream sout;
                sout << "Singular Jacobian, 3 determinant of jacobian = " << detjac << std::endl;
                LOGPZ_ERROR(logger, sout.str())
                DebugStop();
#endif
                detjac = ZeroTolerance();
            }

            jacinv(0,0) = (-jac(1,2)*jac(2,1)+jac(1,1)*jac(2,2))/detjac;//-a12 a21 + a11 a22
            jacinv(0,1) = ( jac(0,2)*jac(2,1)-jac(0,1)*jac(2,2))/detjac;//a02 a21 - a01 a22
            jacinv(0,2) = (-jac(0,2)*jac(1,1)+jac(0,1)*jac(1,2))/detjac;//-a02 a11 + a01 a12
            jacinv(1,0) = ( jac(1,2)*jac(2,0)-jac(1,0)*jac(2,2))/detjac;//a12 a20 - a10 a22
            jacinv(1,1) = (-jac(0,2)*jac(2,0)+jac(0,0)*jac(2,2))/detjac;//-a02 a20 + a00 a22
            jacinv(1,2) = ( jac(0,2)*jac(1,0)-jac(0,0)*jac(1,2))/detjac;//a02 a10 - a00 a12
            jacinv(2,0) = (-jac(1,1)*jac(2,0)+jac(1,0)*jac(2,1))/detjac;//-a11 a20 + a10 a21
            jacinv(2,1) = ( jac(0,1)*jac(2,0)-jac(0,0)*jac(2,1))/detjac;//a01 a20 - a00 a21
            jacinv(2,2) = (-jac(0,1)*jac(1,0)+jac(0,0)*jac(1,1))/detjac;//-a01 a10 + a00 a11

            
            axes.Zero();
            axes(0,0) = 1.0;
            axes(1,1) = 1.0;
            axes(2,2) = 1.0;
            
        }
            break;
            
        default:
        {
            std::cout << " Object with wrong dimensions, unable to compute jacobian matrix. Dimension = " << dim << std::endl;
            DebugStop();
        }
            break;
    }
    
    
}

void TPZGeoEl::JacobianXYZ(const TPZFMatrix<REAL> &gradx, TPZFMatrix<REAL> &jac,TPZFMatrix<REAL> &axesXYZ,REAL &detjac,TPZFMatrix<REAL> &jacinv){

    int nrows = gradx.Rows();
    int ncols = gradx.Cols();
    int dim   = ncols;
    
    switch (dim) {
        case 1:
        {
            jac.Resize(dim,dim);
            axesXYZ.Resize(dim,3);
            jacinv.Resize(dim,dim);
            jac.Zero();
            

            if (gradx.GetVal(0,0) != 0.0) { // X oriented
#ifdef PZDEBUG
                if (! ( IsZero(gradx.GetVal(1,0)) && IsZero(gradx.GetVal(2,0)) ) ) {
                    std::cout << "TPZGeoEl::JacobianXYZ -> Geometric mesh is not oriented with x axis. Should be called TPZGeoEl::Jacobian " << std::endl;
                    DebugStop();
                }
#endif
                jac(0,0)    = gradx.GetVal(0,0);
                detjac      = gradx.GetVal(0,0);
                jacinv(0,0) = 1.0/gradx.GetVal(0,0);
                
                axesXYZ.Zero();
                axesXYZ(0,0) = 1.0;
                
                if(IsZero(detjac))
                {
                    
#ifdef PZDEBUG
                    std::stringstream sout;
                    sout << "Singular Jacobian, 4 determinant of jacobian = " << detjac << std::endl;
                    LOGPZ_ERROR(logger, sout.str())
                    DebugStop();
#endif
                    detjac = ZeroTolerance();
                }
                
                return;
            }
            

            if (gradx.GetVal(1,0) != 0.0) { // Y oriented
#ifdef PZDEBUG
                if (! ( IsZero(gradx.GetVal(0,0)) && IsZero(gradx.GetVal(2,0)) ) ) {
                    std::cout << "TPZGeoEl::JacobianXYZ -> Geometric mesh is not oriented with y axis. Should be called TPZGeoEl::Jacobian " << std::endl;
                    DebugStop();
                }
#endif
                jac(0,0)    = gradx.GetVal(1,0);
                detjac      = gradx.GetVal(1,0);
                jacinv(0,0) = 1.0/gradx.GetVal(1,0);
                
                axesXYZ.Zero();
                axesXYZ(0,1) = 1.0;
                
                if(IsZero(detjac))
                {
                    
#ifdef PZDEBUG
                    std::stringstream sout;
                    sout << "Singular Jacobian, 5 determinant of jacobian = " << detjac << std::endl;
                    LOGPZ_ERROR(logger, sout.str())
                    DebugStop();
#endif
                    detjac = ZeroTolerance();
                }
                
                return;
            }
            
            if (gradx.GetVal(2,0) != 0.0) { // Z oriented
#ifdef PZDEBUG
                if (! ( IsZero(gradx.GetVal(0,0)) && IsZero(gradx.GetVal(1,0)) ) ) {
                    std::cout << "TPZGeoEl::JacobianXYZ -> Geometric mesh is not oriented with z axis. Should be called TPZGeoEl::Jacobian " << std::endl;
                    DebugStop();
                }
#endif
                jac(0,0)    = gradx.GetVal(2,0);
                detjac      = gradx.GetVal(2,0);
                jacinv(0,0) = 1.0/gradx.GetVal(2,0);
                
                axesXYZ.Zero();
                axesXYZ(0,2) = 1.0;
                
                if(IsZero(detjac))
                {
                    
#ifdef PZDEBUG
                    std::stringstream sout;
                    sout << "Singular Jacobian, 6 determinant of jacobian = " << detjac << std::endl;
                    LOGPZ_ERROR(logger, sout.str())
                    DebugStop();
#endif
                    detjac = ZeroTolerance();
                }
                
                return;
            }
            
        }
            break;
        case 2:
        {
            jac.Resize(dim,dim);
            axesXYZ.Resize(dim,3);
            jacinv.Resize(dim,dim);
            jac.Zero();
            
            if (IsZero(gradx.GetVal(2,0)) && IsZero(gradx.GetVal(2,1))) { // XY oriented
#ifdef PZDEBUG
                if ( ( IsZero(gradx.GetVal(0,0)) && IsZero(gradx.GetVal(0,1)) ) ) {
                    std::cout << "TPZGeoEl::JacobianXYZ -> Geometric mesh is not oriented with XY plane. Should be called TPZGeoEl::Jacobian " << std::endl;
                    DebugStop();
                }
#endif
                
                jac(0,0) = gradx.GetVal(0,0);
                jac(1,0) = gradx.GetVal(1,0);
                jac(0,1) = gradx.GetVal(0,1);
                jac(1,1) = gradx.GetVal(1,1);
                
                detjac = jac(0,0)*jac(1,1)-jac(1,0)*jac(0,1);
                jacinv(0,0) = +jac(1,1)/detjac;
                jacinv(1,1) = +jac(0,0)/detjac;
                jacinv(0,1) = -jac(0,1)/detjac;
                jacinv(1,0) = -jac(1,0)/detjac;
                
                axesXYZ.Zero();
                axesXYZ(0,0) = 1.0;
                axesXYZ(1,1) = 1.0;
                
                if(IsZero(detjac))
                {
                    
#ifdef PZDEBUG
                    std::stringstream sout;
                    sout << "Singular Jacobian, 7 determinant of jacobian = " << detjac << std::endl;
                    LOGPZ_ERROR(logger, sout.str())
                    DebugStop();
#endif
                    detjac = ZeroTolerance();
                }
                
                return;
            }
          
            if (IsZero(gradx.GetVal(1,0)) && IsZero(gradx.GetVal(1,1))) { // XZ oriented
#ifdef PZDEBUG
                if ( ( IsZero(gradx.GetVal(0,0)) && IsZero(gradx.GetVal(0,1)) ) ) {
                    std::cout << "TPZGeoEl::JacobianXYZ -> Geometric mesh is not oriented with XZ plane. Should be called TPZGeoEl::Jacobian " << std::endl;
                    DebugStop();
                }
#endif
                
                jac(0,0) = gradx.GetVal(0,0);
                jac(1,0) = gradx.GetVal(2,0);
                jac(0,1) = gradx.GetVal(0,1);
                jac(1,1) = gradx.GetVal(2,1);
                
                detjac = jac(0,0)*jac(1,1)-jac(1,0)*jac(0,1);
                jacinv(0,0) = +jac(1,1)/detjac;
                jacinv(1,1) = +jac(0,0)/detjac;
                jacinv(0,1) = -jac(0,1)/detjac;
                jacinv(1,0) = -jac(1,0)/detjac;
                
                axesXYZ.Zero();
                axesXYZ(0,0) = 1.0;
                axesXYZ(1,2) = 1.0;
                
                if(IsZero(detjac))
                {
                    
#ifdef PZDEBUG
                    std::stringstream sout;
                    sout << "Singular Jacobian, 8 determinant of jacobian = " << detjac << std::endl;
                    LOGPZ_ERROR(logger, sout.str())
                    DebugStop();
#endif
                    detjac = ZeroTolerance();
                }
                
                return;
            }
            
            if (IsZero(gradx.GetVal(0,0)) && IsZero(gradx.GetVal(0,1))) { // YZ oriented
#ifdef PZDEBUG
                if ( ( IsZero(gradx.GetVal(1,0)) && IsZero(gradx.GetVal(1,1)) ) ) {
                    std::cout << "TPZGeoEl::JacobianXYZ -> Geometric mesh is not oriented with YZ plane. Should be called TPZGeoEl::Jacobian " << std::endl;
                    DebugStop();
                }
#endif
                
                jac(0,0) = gradx.GetVal(1,0);
                jac(1,0) = gradx.GetVal(2,0);
                jac(0,1) = gradx.GetVal(1,1);
                jac(1,1) = gradx.GetVal(2,1);
                
                detjac = jac(0,0)*jac(1,1)-jac(1,0)*jac(0,1);
                jacinv(0,0) = +jac(1,1)/detjac;
                jacinv(1,1) = +jac(0,0)/detjac;
                jacinv(0,1) = -jac(0,1)/detjac;
                jacinv(1,0) = -jac(1,0)/detjac;
                
                axesXYZ.Zero();
                axesXYZ(0,1) = 1.0;
                axesXYZ(1,2) = 1.0;
                
                if(IsZero(detjac))
                {
                    
#ifdef PZDEBUG
                    std::stringstream sout;
                    sout << "Singular Jacobian, 9 determinant of jacobian = " << detjac << std::endl;
                    LOGPZ_ERROR(logger, sout.str())
                    DebugStop();
#endif
                    detjac = ZeroTolerance();
                }
                
                return;
            }

            
        }
            break;
        case 3:
        {
            jac.Resize(dim,dim);
            axesXYZ.Resize(dim,3);
            jacinv.Resize(dim,dim);
            jac.Zero();
            
            for (int i = 0; i < nrows; i++) {
                jac(i,0)  = gradx.GetVal(i,0);
                jac(i,1)  = gradx.GetVal(i,1);
                jac(i,2)  = gradx.GetVal(i,2);
            }
            
            detjac = -jac(0,2)*jac(1,1)*jac(2,0);//-a02 a11 a20
            detjac += jac(0,1)*jac(1,2)*jac(2,0);//+ a01 a12 a20
            detjac += jac(0,2)*jac(1,0)*jac(2,1);//+ a02 a10 a21
            detjac -= jac(0,0)*jac(1,2)*jac(2,1);//- a00 a12 a21
            detjac -= jac(0,1)*jac(1,0)*jac(2,2);//- a01 a10 a22
            detjac += jac(0,0)*jac(1,1)*jac(2,2);//+ a00 a11 a22
            
            jacinv(0,0) = (-jac(1,2)*jac(2,1)+jac(1,1)*jac(2,2))/detjac;//-a12 a21 + a11 a22
            jacinv(0,1) = ( jac(0,2)*jac(2,1)-jac(0,1)*jac(2,2))/detjac;//a02 a21 - a01 a22
            jacinv(0,2) = (-jac(0,2)*jac(1,1)+jac(0,1)*jac(1,2))/detjac;//-a02 a11 + a01 a12
            jacinv(1,0) = ( jac(1,2)*jac(2,0)-jac(1,0)*jac(2,2))/detjac;//a12 a20 - a10 a22
            jacinv(1,1) = (-jac(0,2)*jac(2,0)+jac(0,0)*jac(2,2))/detjac;//-a02 a20 + a00 a22
            jacinv(1,2) = ( jac(0,2)*jac(1,0)-jac(0,0)*jac(1,2))/detjac;//a02 a10 - a00 a12
            jacinv(2,0) = (-jac(1,1)*jac(2,0)+jac(1,0)*jac(2,1))/detjac;//-a11 a20 + a10 a21
            jacinv(2,1) = ( jac(0,1)*jac(2,0)-jac(0,0)*jac(2,1))/detjac;//a01 a20 - a00 a21
            jacinv(2,2) = (-jac(0,1)*jac(1,0)+jac(0,0)*jac(1,1))/detjac;//-a01 a10 + a00 a11
            
            if(IsZero(detjac))
            {
                
#ifdef PZDEBUG
                std::stringstream sout;
                sout << "Singular Jacobian, 10 determinant of jacobian = " << detjac << std::endl;
                LOGPZ_ERROR(logger, sout.str())
                DebugStop();
#endif
                detjac = ZeroTolerance();
            }
            
            axesXYZ.Zero();
            axesXYZ(0,0) = 1.0;
            axesXYZ(1,1) = 1.0;
            axesXYZ(2,2) = 1.0;
            
        }
            break;
            
        default:
        {
            std::cout << " Object with wrong dimensions, unable to compute jacobian matrix. Dimension = " << dim << std::endl;
            DebugStop();
        }
            break;
    }
    
}

/** Defines the refinement pattern. It's used only in TPZGeoElRefPattern objects. */
void TPZGeoEl::SetRefPattern(TPZAutoPointer<TPZRefPattern> ){
	PZError << "TPZGeoEl::SetRefPattern ERROR : Should not be called in TPZGeoEl" << endl;
    DebugStop();
}

void TPZGeoEl::Read(TPZStream &buf, void *context) { //ok
    fMesh = dynamic_cast<TPZGeoMesh *>(TPZPersistenceManager::GetInstance(&buf));
    buf.Read(&fId,1);
    buf.Read(&fMatId,1);
    fReference = dynamic_cast<TPZCompEl *>(TPZPersistenceManager::GetInstance(&buf));
    buf.Read(&fFatherIndex,1);
    buf.Read(&fIndex,1);
    gGlobalAxes.Read(buf, 0);
    buf.Read(&fNumInterfaces,1);
}

void TPZGeoEl::Write(TPZStream &buf, int withclassid) const { //ok
    TPZPersistenceManager::WritePointer(fMesh, &buf);
    buf.Write(&fId, 1);
    buf.Write(&fMatId, 1);
    TPZPersistenceManager::WritePointer(fReference, &buf);
    buf.Write(&fFatherIndex, 1);
    buf.Write(&fIndex, 1);
    gGlobalAxes.Write(buf, 0);
    buf.Write(&fNumInterfaces, 1);
}

TPZGeoEl::TPZGeoEl(TPZGeoMesh & DestMesh, const TPZGeoEl &cp):TPZSavable(cp){
	this->fMesh = &DestMesh;
	this->fId = cp.fId;
	this->fMatId = cp.fMatId;
	this->fReference = cp.fReference;
	this->fFatherIndex = cp.fFatherIndex;
	this->fIndex = cp.fIndex;
	this->fMesh->ElementVec()[this->fIndex] = this;
	this->fNumInterfaces = 0;
}

TPZGeoEl::TPZGeoEl(TPZGeoMesh & DestMesh, const TPZGeoEl &cp, std::map<int64_t,int64_t> &org2clnMap):TPZSavable(cp){
	this->fMesh = &DestMesh;
	this->fId = cp.fId;
	this->fMatId = cp.fMatId;
	this->fReference = 0;
	if ( cp.fFatherIndex == -1) this->fFatherIndex = -1;
	else if (org2clnMap.find(cp.fFatherIndex) == org2clnMap.end())
	{
		std::stringstream sout;
		sout << "ERROR in - " << __PRETTY_FUNCTION__
        << " original father element index: " << cp.fIndex << " is not mapped!";
		LOGPZ_ERROR (logger,sout.str().c_str());
		DebugStop();
	}
	else this->fFatherIndex = org2clnMap[cp.fFatherIndex];
	
	if (org2clnMap.find(cp.fIndex) == org2clnMap.end())
	{
		std::stringstream sout;
		sout << "ERROR in - " << __PRETTY_FUNCTION__
		<< " original element index: " << cp.fIndex << " is not mapped!";
		LOGPZ_ERROR (logger,sout.str().c_str());
		DebugStop();
	}
	
	this->fIndex = org2clnMap[cp.fIndex];
	this->fMesh->ElementVec()[this->fIndex] = this;
	this->fNumInterfaces = 0;
}

// return the refinement pattern associated with the element
TPZAutoPointer<TPZRefPattern> TPZGeoEl::GetRefPattern() const
{
	TPZAutoPointer<TPZRefPattern> result;
	return result;
}

bool TPZGeoEl::VerifyNodeCoordinates(REAL tol){
	const int nnodes = this->NCornerNodes();
	TPZManVector<REAL,3> qsi(this->Dimension());
	TPZManVector<REAL,3> MappedX(3), NodeX(3);
	for(int inode = 0; inode < nnodes; inode++){
		for(int dim = 0; dim < 3; dim++){
			NodeX[dim] = this->NodePtr(inode)->Coord(dim);
		}//dim
		this->CenterPoint(inode,qsi);
		this->X(qsi,MappedX);
		REAL error = 0.;
		for(int dim = 0; dim < 3; dim++){
			error += (NodeX[dim]-MappedX[dim])*(NodeX[dim]-MappedX[dim]);
		}//dim
		error = sqrt(error);
		if(error > tol || !(error==error)){
			std::stringstream mess;
			mess << "FATAL ERROR AT " << __PRETTY_FUNCTION__ << " - Node coordinate differs from mapped node.\n";
			this->Print(mess);
			PZError << mess.str() << "\n";
#ifdef PZ_LOG
			LOGPZ_ERROR(logger,mess.str().c_str());
#endif
			DebugStop();
			return false;
		}
	}//for i
	return true;
}//method


#include "tpzgeoelrefpattern.h"
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


/** @brief Find a side which is not contained in allsides and whose dimension is one higher than the dimension of the side side at gel */
int ConjugateSide(TPZGeoEl *gel, int side, TPZStack<int> &allsides);
/** @brief Computes the normal vector goes from the center of the conjugate side (LC) to the center of the LS side (orthogonal this) */
void NormalVector(TPZGeoElSide &LC, TPZGeoElSide &LS, TPZVec<REAL> &normal);
/** @brief Normalize normlow vector (in self) */
void Normalize(TPZVec<REAL> &normlow, TPZVec<REAL> &normal);

//void TPZGeoEl::ComputeNormals(TPZMatrix<REAL> &normals)
//{
//    
//#ifdef PZ_LOG
//    {
//		std::stringstream sout;
//		sout<< "Metodo Compute normal \n";
//		
//		LOGPZ_DEBUG(logger,sout.str())
//    }
//#endif
//
//    
//	int numbernormals = 0;
//	int dimension = Dimension();
//	int is;
//	int nsides = NSides();
//	for(is=0; is<nsides; is++)
//	{
//        TPZStack<int> lowdim;
//		if(SideDimension(is) == dimension-1)
//		{
//			//TPZStack<int> lowdim;
//			LowerDimensionSides(is,lowdim);
//			numbernormals += lowdim.NElements()+1;
//		}
//        
//	}
//	normals.Redim(3, numbernormals);
//	int counter = 0;
//	for(is=0; is<nsides; is++)
//	{
//		if(SideDimension(is) == dimension-1)
//            
//		{
//            
//#ifdef PZ_LOG
//            {
//                std::stringstream sout;
//                sout<< "Side "<<is<<std::endl;
//            
//                LOGPZ_DEBUG(logger,sout.str())
//            }
//#endif
//
//            
//			TPZStack<int> lowdim;
//			LowerDimensionSides(is,lowdim);
//            lowdim.Push(is);
//#ifdef PZ_LOG
//            {
//                std::stringstream sout;
//                sout<< "LowerDimensionSides "<<lowdim<<std::endl;
//                
//                LOGPZ_DEBUG(logger,sout.str())
//            }
//#endif
//            
//			int nlowdim = lowdim.NElements();
//			int lowis;
//			for(lowis=0; lowis < nlowdim; lowis++)
//			{
//				int conj_side = ConjugateSide(this,lowdim[lowis],lowdim);
//				TPZGeoElSide LC(this,conj_side);
//				TPZGeoElSide LS(this,lowdim[lowis]);
//#ifdef PZ_LOG
//                {
//                    std::stringstream sout;
//                    sout<< "Side "<<is << " Conjugate Side "<< conj_side<< " Ls side "<< lowdim[lowis]<<std::endl;
//                    
//                    LOGPZ_DEBUG(logger,sout.str())
//                }
//#endif
//				TPZManVector<REAL> normal(3,0.);
//				NormalVector(LC,LS,normal);
//                
//#ifdef PZ_LOG
//                {
//                    std::stringstream sout;
//                    sout<< "Vetores do NormalVector "<<normal<<std::endl;
//                    
//                    LOGPZ_DEBUG(logger,sout.str())
//                }
//#endif
//                
//				int d;
//				for(d=0; d<3; d++) normals(d,counter) = normal[d];
//				counter++;
//			}
//            
//#ifdef PZ_LOG
//            {
//                std::stringstream sout;
//                sout<< "A partir daqui sera um processo de normalizacao dos vetores normals"<<std::endl;
//                normals.Print("vetor normals",sout);
//                
//                LOGPZ_DEBUG(logger,sout.str())
//            }
//#endif
//			TPZManVector<REAL,3> normal(3,0.);
//			int d;
//			for(d=0; d<3; d++) normal[d] = normals(d,counter-1);
//			for(lowis = counter-nlowdim; lowis < counter; lowis++)
//			{
//				TPZManVector<REAL,3> normlow(3,0.);
//				for(d=0; d<3; d++) normlow[d] = normals(d,lowis);
//				Normalize(normlow,normal);
//				for(d=0; d<3; d++) normals(d,lowis) = normlow[d];
//			}
//		}
//	}
//    
//#ifdef PZ_LOG
//    {
//        std::stringstream sout;
//        sout<< "Vetores Normais normalizados "<<normals<<std::endl;
//        
//        LOGPZ_DEBUG(logger,sout.str())
//    }
//#endif
//}

void TPZGeoEl::SetNeighbourForBlending(int side){
	if( !this->IsGeoBlendEl() ) return;
	
	TPZGeoElSide ElemSide(this,side);
	TPZGeoElSide NextSide(this,side);
	
	//trying to get a self non-linear neighbour
	TPZGeoElSide NextSideNeighbour = NextSide.Neighbour();
	while(NextSideNeighbour.Element() != ElemSide.Element())
	{
		if(NextSideNeighbour.Exists() && !NextSideNeighbour.Element()->IsLinearMapping() && !NextSideNeighbour.Element()->IsGeoBlendEl())
		{
			if(NextSideNeighbour.IsRelative(ElemSide) == false){
				if(NextSideNeighbour.Element()->IsGeoElMapped() == false){
					TPZGeoElSide NeighSide = NextSideNeighbour;
					TPZTransform<> NeighTransf(NeighSide.Dimension(),NeighSide.Dimension());
					ElemSide.SideTransform3(NeighSide,NeighTransf);
					this->SetNeighbourInfo(side,NeighSide,NeighTransf);
					return;
				}// !TPZGeoElMapped
			}//if IsRelative == false
		}
		NextSide = NextSide.Neighbour();
        NextSideNeighbour = NextSide.Neighbour();
	}
    NextSide = TPZGeoElSide(this,side);
    NextSideNeighbour = NextSide.Neighbour();
	//now TPZGeoElMapped are accepted
	while(NextSideNeighbour.Element() != ElemSide.Element())
	{
		if(NextSideNeighbour.Exists() && !NextSideNeighbour.Element()->IsLinearMapping() && !NextSideNeighbour.Element()->IsGeoBlendEl())
		{
			if(NextSideNeighbour.IsRelative(ElemSide) == false){
                if(NextSideNeighbour.Element()->IsGeoElMapped() == true) {
                    // only use this element as map if this element has a REAL father
                    // meaning the father is not a neighbour
                    // if the father is a neighbour, then we will get to the father looping over the neighbours
                    auto NextSideNeighbourFather = NextSideNeighbour.Father2();
                    if(!NextSideNeighbour.IsNeighbour(NextSideNeighbourFather))
                    {
                        TPZGeoElSide NeighSide = NextSideNeighbour;
                        TPZTransform<> NeighTransf(NeighSide.Dimension(),NeighSide.Dimension());
                        ElemSide.SideTransform3(NeighSide,NeighTransf);
                        this->SetNeighbourInfo(side,NeighSide,NeighTransf);
                        return;
                    }
                }
			}//if IsRelative == false
		}
		NextSide = NextSide.Neighbour();
        NextSideNeighbour = NextSide.Neighbour();
	}
	
}//method

void TPZGeoEl::BuildBlendConnectivity(){
	if( !this->IsGeoBlendEl() ) return;
	const int nsides = this->NSides();
	for(int byside = this->NNodes(); byside < nsides; byside++)
	{
		this->SetNeighbourForBlending(byside);
	}//for byside
}//method

// Projection of the point to the side
// Compute the projection of the point within the interior of the element to the side of the element
TPZTransform<> TPZGeoEl::Projection(int side)
{
	TPZTransform<> tr = SideToSideTransform(NSides()-1,side);
	TPZTransform<> tr2 = SideToSideTransform(side,NSides()-1);
	TPZTransform<> tr3 = tr2.Multiply(tr);
	//	std::cout << "interior to side " << side << " trans " << tr << std::endl;
	//	std::cout << "side to interior " << side << " trans " << tr2 << std::endl;
	//	std::cout << "side to side " << side << " trans " << tr3 << std::endl;
	return tr3;
}

TPZGeoEl * TPZGeoEl::EldestAncestor() const
{
    TPZGeoEl *father = this->Father();
    if(!father)
    {
        return const_cast<TPZGeoEl*>(this); // we have to cast away the constness given by the const qualifier of this function
    }
    TPZGeoEl *nextfather = NULL;
    if(father) nextfather = father->Father();
    while(nextfather)
    {
        father = nextfather;
        nextfather = father->Father();
    }
    
    return father;
}///method

void NormalVector(TPZGeoElSide &LC, TPZGeoElSide &LS, TPZVec<REAL> &normal)
{
	//TPZGeoEl *gel = LC.Element();
	// take the centerpoint of LC and the centerpoint of LS
	TPZManVector<REAL,3> LCCenter(3,0.), LSCenter(3,0.);
	TPZManVector<REAL,3> XLC(3,0.),XLS(3,0.);
	LC.CenterPoint(LCCenter);
	LC.X(LCCenter,XLC);
    
	LS.CenterPoint(LSCenter);
	LS.X(LSCenter,XLS);
	TPZManVector<REAL,3> dir(3,0.);
    
#ifdef PZ_LOG
    if(logger.isDebugEnabled())
    {
        std::stringstream sout;
        sout<< "Centro de LC "<<XLC<< " Centro de LS "<<XLS<<std::endl;
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    
    
	// The normal vector needs to be in the plane of LC and perpendicular to LS
	// A starting vector is the direction of one center to the next
	int i;
	for(i=0; i<3; i++) dir[i] = XLS[i]-XLC[i];
    
    
    
	TPZFNMatrix<10> jacobian,axes,jacinv;
	REAL detjac;
	LS.Jacobian(LSCenter,jacobian,axes,detjac,jacinv);
    
#ifdef PZ_LOG
    if(logger.isDebugEnabled())
    {
        std::stringstream sout;
        sout<< "vetor axes "<<axes<<std::endl;
         sout<< "vetor dir "<<dir<<std::endl;
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    
	TPZFNMatrix<20> axtrans(axes.Cols(),axes.Rows()+1,0.);
	int j;
	for(i=0; i<3; i++) for(j=0; j<axes.Rows(); j++)
	{
		axtrans(i,j) = axes(j,i);
	}
	int lastcol = j;
	for(i=0; i<3; i++)
	{
		axtrans(i,lastcol) = dir[i];
	}
    
    
#ifdef PZ_LOG
    if(logger.isDebugEnabled())
	{
        std::stringstream sout;
		sout << "axtrans = ";
        axtrans.Print("axtrans matrix",sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
    
    
	// Then orthogonalize the vectors of the LS side with this vector
	TPZFNMatrix<20> transf, ortho;
	axtrans.GramSchmidt(ortho,transf);
    
    
#ifdef PZ_LOG
    if(logger.isDebugEnabled())
	{
        std::stringstream sout;
		sout << "Vetr Apos GramSchmidt = ";
        ortho.Print("Orthog matrix",sout);
		transf.Print("TransfToOrthog matrix",sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
    
    
	normal.Resize(3);
	for(i=0; i<3; i++)
	{
		normal[i] = ortho(i,lastcol);
	}
    
    
}

void Normalize(TPZVec<REAL> &normlow, TPZVec<REAL> &normal)
{
#ifdef PZ_LOG
    if(logger.isDebugEnabled())
    {
        std::stringstream sout;
        sout<< "Normalize..u.n=1 "<<std::endl;
        sout<< "normlow "<<normlow<< " normal "<<normal<<std::endl;
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    
    
	REAL inner = 0.;
	int i;
	for(i=0; i<normlow.NElements(); i++)
	{
		inner += normlow[i]*normal[i];
	}
    
    if(inner==0)
	{
		LOGPZ_ERROR(logger,"Inner product zero to normalize vector")
		DebugStop();
	}
    
    
    
	for(i=0; i<normlow.NElements(); i++)
	{
        normlow[i] /= inner;
	}
}

//void TPZGeoEl::ComputeNormals(int side, TPZFMatrix<REAL> &normals, TPZVec<int> &vectorsides)
//{
//    int numbernormals = 0;
//    int dimension = Dimension();
//    int sidedimension = SideDimension(side);
//    TPZStack<int> lowdim;
//    LowerDimensionSides(side,lowdim);
//    lowdim.Push(side);
//    // the normals corresponding to the internal shape functions
//    // Compute the number of normals we need to compute
//    int nsides = NSides();
//    if(sidedimension == dimension-1)
//    {
//        numbernormals = lowdim.NElements();
//    }
//    else if(sidedimension == dimension)
//    {
//        numbernormals = nsides*2;
//    }
//    else
//    {
//        numbernormals = 0;
//    }
//    normals.Redim(3, numbernormals);
//    vectorsides.Resize(numbernormals);
//    vectorsides.Fill(0);
//    if(!numbernormals)
//    {
//        return;
//    }
//    int is = side;
//    if(sidedimension == dimension-1)
//    {
//        int nlowdim = lowdim.NElements();
//        int lowis;
//        // work from lowest dimension sides upward
//        for(lowis=0; lowis < nlowdim; lowis++)
//        {
//            // find a side which is not contained in the currently analysed side
//            // whose dimension is one higher than the dimension of lowdim[lowis]
//            int conj_side = ConjugateSide(this,lowdim[lowis],lowdim);
//            // the normal vector is in the alignment of the conjugate side
//            TPZGeoElSide LC(this,conj_side);
//            TPZGeoElSide LS(this,lowdim[lowis]);
//            TPZManVector<REAL> normal(3,0.);
//            // the normal vector goes from the center of the conjugate side to
//            // the center of the LS side
//            NormalVector(LC,LS,normal);
//            int d;
//            for(d=0; d<3; d++) normals(d,lowis) = normal[d];
//            vectorsides[lowis] = lowdim[lowis];
//        }
//        // Why mormalize a vector which has just been computed?
//        TPZManVector<REAL> normal(3,0.);
//        int d;
//        for(d=0; d<3; d++) normal[d] = normals(d,nlowdim-1);
//        for(lowis = 0; lowis < nlowdim; lowis++)
//        {
//            TPZManVector<REAL,3> normlow(3,0.);
//            for(d=0; d<3; d++) normlow[d] = normals(d,lowis);
//            Normalize(normlow,normal);
//            for(d=0; d<3; d++) normals(d,lowis) = normlow[d];
//        }
//        TPZManVector<int,9> sidepermutationgather(nlowdim);
////        HDivPermutation(is,sidepermutationgather);
//        for(int i=0; i< nlowdim; i++) sidepermutationgather[i] = i;
//        TPZFNMatrix<12> sidenormals(3,nlowdim);
//        TPZManVector<int> localvecsides(nlowdim);
//        // compute whether the side is from this element to the next or contrary
//        int sideorient = 1;//NormalOrientation(side);
//        int i;
//        for(i=0; i<nlowdim; i++)
//        {
//            for(d=0; d<3; d++)
//            {
//                sidenormals(d,i) = normals(d,sidepermutationgather[i]);
//            }
//            localvecsides[i] = vectorsides[sidepermutationgather[i]];
//        }
//        for(i=0; i<nlowdim; i++)
//        {
//            for(d=0; d<3; d++)
//            {
//                normals(d,i) = sidenormals(d,i)*sideorient;
//            }
//            vectorsides[i] = localvecsides[i];
//        }
//    }
//    else if(sidedimension == dimension)
//    {
//        int counter = 0;
//        for(is=0; is<nsides; is++)
//        {
//            if(SideDimension(is) > 0)
//            {
//                TPZManVector<REAL,3> Center(3,0.);
//                TPZManVector<REAL,3> X;
//                TPZFNMatrix<10> jacobian, jacinv, axes;
//                REAL detjac;
//                TPZGeoElSide gelside(this,is);
//                gelside.CenterPoint(Center);
//                gelside.Jacobian(Center,jacobian,axes,detjac,jacinv);
//                int d,s;
//                dimension = SideDimension(is);
//                for(d=0; d<dimension; d++)
//                {
//                    for(s=0; s<3; s++)
//                    {
//                        normals(s,counter) = axes(d,s);
//                    }
//                    vectorsides[counter] = is;
//                    counter++;
//                }
//            }
//        }
//        vectorsides.Resize(counter);
//        normals.Resize(3,counter);
//    }
//}


//void TPZGeoEl::ComputeNormalsDG(int side, TPZVec<REAL> &pt, TPZFMatrix<REAL> &normals, TPZVec<int> &vectorsides)
//{
//    if (SideDimension(side) >= Dimension()-1) {
//        Directions(side, pt, normals, vectorsides);
//    }
//    if (SideDimension(side) == Dimension()-1)
//    {
//        // we need to permute the normals and associated sides
//        TPZGeoElSide thisside(this,side);
//        int nlowdim = thisside.NSides();
//        
//        TPZManVector<int,9> sidepermutationgather(nlowdim);
//        HDivPermutation(side,sidepermutationgather);
//#ifdef PZ_LOG
//        if(logger.isDebugEnabled()){
//            std::stringstream sout;
//            sout << "Permutation for side " << side << " is " << sidepermutationgather;
//            LOGPZ_DEBUG(logger, sout.str())
//        }
//#endif
//        TPZFNMatrix<12> sidenormals(3,nlowdim);
//        TPZManVector<int> localvecsides(nlowdim);
//        // compute whether the side is from this element to the next or contrary
//        int sideorient = NormalOrientation(side);
//        int i;
//        for(i=0; i<nlowdim; i++)
//        {
//            for(int d=0; d<3; d++)
//            {
//                sidenormals(d,i) = normals(d,sidepermutationgather[i]);
//            }
//            localvecsides[i] = vectorsides[sidepermutationgather[i]];
//        }
//        for(i=0; i<nlowdim; i++)
//        {
//            for(int d=0; d<3; d++)
//            {
//                normals(d,i) = sidenormals(d,i)*sideorient;
//            }
//            vectorsides[i] = localvecsides[i];
//        }
//        
//    }
//
//}

//void TPZGeoEl::ComputeNormals(TPZFMatrix<REAL> &normals, TPZVec<int> &vectorsides)
//{
//    int numbernormals = 0;
//    int dimension = Dimension();
//    // the normals corresponding to the internal shape functions
//    int is;
//    // Compute the number of normals we need to compute
//    int nsides = NSides();
//    numbernormals = nsides*dimension; // @omar:: why two???
//    normals.Redim(3, numbernormals);
//    vectorsides.Resize(numbernormals);
//    vectorsides.Fill(0);
//    int counter = 0;
//    // effectively compute the normals
//    for(is=0; is<nsides; is++)
//    {
//        TPZFNMatrix<100> sidenormals;
//        TPZManVector<int> sidevectors;
//        ComputeNormals(is,sidenormals,sidevectors);
//        int numnormals = sidevectors.NElements();
//        int in;
//        for(in=0; in<numnormals; in++)
//        {
//            int d;
//            for(d=0; d<3; d++)
//            {
//                normals(d,counter) = sidenormals(d,in);
//            }
//            vectorsides[counter] = sidevectors[in];
//            counter++;
//        }
//    }
//#ifdef PZDEBUG
//    if(counter != numbernormals)
//    {
//        DebugStop();
//    }
//#endif
//}

//void TPZGeoEl::ComputeNormalsDG(TPZVec<REAL> &pt, TPZFMatrix<REAL> &normals, TPZVec<int> &vectorsides)
//{
//    int numbernormals = 0;
//    //    int dimension = Dimension();
//    // the normals corresponding to the internal shape functions
//    int is;
//    // Compute the number of normals we need to compute
//    int nsides = NSides();
//    numbernormals = nsides*Dimension();
//    normals.Redim(3, numbernormals);
//    vectorsides.Resize(numbernormals);
//    vectorsides.Fill(0);
//    int counter = 0;
//    // effectively compute the normals
//    for(is=0; is<nsides; is++)
//    {
//        TPZFNMatrix<100> sidenormals;
//        TPZManVector<int> sidevectors;
//        ComputeNormalsDG(is,pt, sidenormals,sidevectors);
//        int numnormals = sidevectors.NElements();
//        int in;
//        for(in=0; in<numnormals; in++)
//        {
//            int d;
//            for(d=0; d<3; d++)
//            {
//                normals(d,counter) = sidenormals(d,in);
//            }
//            vectorsides[counter] = sidevectors[in];
//            counter++;
//        }
//    }
//}


// Determine the orientation of the normal vector comparing the ids of the neighbouring elements
int TPZGeoEl::NormalOrientation(int side)
{
	int dimel = Dimension();
	int dimside = SideDimension(side);
	if(dimside != dimel-1)
	{
		LOGPZ_ERROR(logger,"NormalOrientation called with wrong side")
		return 0;
	}
#ifdef PZ_LOG
    if (loggerorient.isDebugEnabled())
    {
        std::stringstream sout;
        TPZCompEl *cel = Reference();
        if (cel) {
            sout << "Comp Element index " << cel->Index() << std::endl;
        }
        sout << "Element index " << this->Index() << " Id = " << this->Id() << " side = " << side;
        LOGPZ_DEBUG(loggerorient, sout.str())
    }
#endif
	TPZGeoElSide thisside(this,side);

	TPZGeoElSide fatherside = thisside.Father2();
	while(fatherside.Exists() && fatherside.Dimension() == dimside)//a segunda condicional eu inclui agora
	{
		thisside = fatherside;
		fatherside = fatherside.Father2();
	}

#ifdef PZ_LOG
    if (loggerorient.isDebugEnabled())
    {
        std::stringstream sout;
        sout << "largest father index " << thisside.Element()->Index() << " id = " << thisside.Element()->Id() << std::endl;
        LOGPZ_DEBUG(loggerorient, sout.str())
    }
#endif

    
    TPZGeoElSide neighbour(thisside.Neighbour());
    // case there is no neighbour (should put a debugstop here)
	if(!neighbour.Exists() ) 
	{
		return 1;
	}
    fatherside = neighbour.Father2();
    // look for a neighbour of equal dimension
    while ((neighbour.Element()->Dimension() != Dimension() && neighbour != thisside) || (fatherside && fatherside.Dimension() == dimside)) {
        neighbour = neighbour.Neighbour();
        fatherside = neighbour.Father2();
    }
    if (neighbour == thisside) {
        return 1;
    }
	
#ifdef PZDEBUG
	if(!thisside.NeighbourExists(neighbour))//inclui agora esta verificacao
	{
		std::stringstream sout;
#ifdef PZ_LOG
		LOGPZ_ERROR(logger,sout.str().c_str());
#endif
		DebugStop();
		
	}
#endif
#ifdef PZ_LOG
    if (loggerorient.isDebugEnabled())
    {
        std::stringstream sout;
        sout << "Element index " << Index() << std::endl;
        sout << "neighbour index " << neighbour.Element()->Index() << " id = " << neighbour.Element()->Id() << " side " << neighbour.Side() <<  std::endl;
        sout << "thisside index " << thisside.Element()->Index() << " id = " << thisside.Element()->Id() << " side " << thisside.Side() << std::endl;
        if(thisside.Element()->Id() < neighbour.Element()->Id())
        {
            sout << "returning 1\n";
        }
        else {
            sout << "returning -1\n";
        }
        LOGPZ_DEBUG(loggerorient, sout.str())
    }
#endif
	
	if(thisside.Element()->Id() < neighbour.Element()->Id())
	{
		return 1;
	}
	else {
		return -1;
	}
    
}

int ConjugateSide(TPZGeoEl *gel, int side, TPZStack<int> &allsides)
{
	std::set<int> allside;
	allside.insert(&allsides[0],&allsides[0]+allsides.NElements());
	TPZStack<TPZGeoElSide> highsides;
	// the dimension of the conjugate side
	int targetdimension = gel->SideDimension(side)+1;
	// find all sides connected to side which have target dimension
	gel->AllHigherDimensionSides(side,targetdimension,highsides);
	int nhigh = highsides.NElements();
	int is;
	// find the side which is not contained in allsides
	for(is=0; is<nhigh; is++)
	{
		int highside = highsides[is].Side();
		if(allside.find(highside) == allside.end())
		{
			return highside;
		}
	}
	return -1;
}

void TPZGeoEl::GetHigherSubElements(TPZVec<TPZGeoEl*> &unrefinedSons)
{
    int nsons = this->NSubElements();
    for(int s = 0; s < nsons; s++)
    {
        TPZGeoEl * son = this->SubElement(s);
        if(!son) continue;
    
        if(son->HasSubElement() == false)
        {
            int oldSize = unrefinedSons.NElements();
            unrefinedSons.Resize(oldSize+1, son);
        }
        else
        {
            son->GetHigherSubElements(unrefinedSons);
        }
    }
}

// Compute the permutation for an HDiv side
void TPZGeoEl::HDivPermutation(int side, TPZVec<int> &permutegather)
{
	int dimension = Dimension();
	int sidedimension = SideDimension(side);
    
	if(dimension != sidedimension+1)
	{
		std::stringstream sout;
		sout << "HDivPermutation called with wrong side parameter " << side;
#ifdef PZ_LOG
		LOGPZ_ERROR(logger,sout.str())
#endif
		cout << sout.str() << std::endl;
	}
    DebugStop();
}

void TPZGeoEl::GetNodeIndices( TPZVec<int64_t> &nodeindices ){
    const int nnodes = this->NNodes();
    nodeindices.Resize( nnodes );
    for(int i = 0; i < nnodes; i++){
        nodeindices[i] = this->NodeIndex(i);
    }
}///void

void TPZGeoEl::GetNodeIndices( std::set<int64_t> &nodeindices ){
    nodeindices.clear();
    const int nnodes = this->NNodes();
    for(int i = 0; i < nnodes; i++){
        nodeindices.insert(this->NodeIndex(i));
    }
}///void


void TPZGeoEl::NodesCoordinates(TPZFMatrix<REAL > &coordinates){
    coordinates.Resize(3,this->NNodes());
    coordinates.Zero();
    TPZManVector<REAL,3> co(3,0.0);
    
    for (int in = 0; in < this->NNodes(); in++) {
        TPZGeoNode inode = this->Node(in);
        inode.GetCoordinates(co);
        coordinates(0,in) = co[0];
        coordinates(1,in) = co[1];
        coordinates(2,in) = co[2];
    }
    
}

int TPZGeoEl::StaticClassId() {
    return Hash("TPZGeoEl");
}

int TPZGeoEl::ClassId() const {
    return StaticClassId();
}

TPZGeoNode* TPZGeoEl::NodePtr(int i) const {
    return &(fMesh->NodeVec()[NodeIndex(i)]);
}

TPZGeoNode& TPZGeoEl::Node(int i) const {
    return (fMesh->NodeVec()[NodeIndex(i)]);
}

TPZGeoEl* TPZGeoEl::Father() const {
    return (fFatherIndex == -1) ? 0 : Mesh()->ElementVec()[fFatherIndex];
}

TPZGeoNode* TPZGeoEl::SideNodePtr(int side, int nodenum) const {
    return &(fMesh->NodeVec()[SideNodeIndex(side, nodenum)]);
}
