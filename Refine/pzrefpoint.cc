#include "pzrefpoint.h"
#include "pzshapelinear.h"
#include "pzelgpoint.h"

static int NumInNeigh = 0;
static int InNeigh[1][1][1] = {
	{{0}}
};

static int CornerSons[1][1] = {
	{0}
};

static int MidSideNodes[1][1]  = {
	{0}
};

static REAL MidCoord[1][1] = { 
	{0.} 
};

static int subeldata[1][1][1] =
{
	{{0}}
};

static int nsubeldata[1] = {0};

static REAL buildt[1][1][1][1] = {//por colunas
/*S0*/
	{{{0.}}}
};

static int fatherside[1][1] = {
	{0}
};


void TPZRefPoint::Divide(TPZGeoEl *geo,TPZVec<TPZGeoEl *> &SubElVec) {
	int i;
	if(geo->HasSubElement()) {
		SubElVec.Resize(NSubEl);
		for(i=0;i<NSubEl;i++) 
			SubElVec[i] = geo->SubElement(i);
		return;//If exist fSubEl return this sons
	}
	//int j,index;
	int sub,matid=geo->MaterialId();
	int np[1/*TPZShapePoint::NSides*/];//guarda conectividades dos 8 subelementos

	np[0] = geo->NodeIndex(0);
//  for(j=0;j<1/*TPZShapeLinear::NNodes*/;j++) 
//	  np[j] = geo->NodeIndex(j);
//  for(sub=1/*TPZShapeLinear::NNodes*/;sub<1;sub++) {
//    NewMidSideNode(geo,sub,index);
//    np[sub] = index;
//  }
  // creating new subelements
	TPZGeoElPoint *pt0d = (TPZGeoElPoint *) geo;
	//for(i=0;i<TPZShapeLinear::NNodes;i++) {
		TPZManVector<int>  cornerindexes(1/*TPZShapeLinear::NNodes*/);
		cornerindexes[0]=np[0];
		//for(int j=0;j<TPZShapeLinear::NNodes;j++) 
		//	cornerindexes[j] = np[CornerSons[i][j]];
		TPZGeoElPoint *npt = new TPZGeoElPoint(cornerindexes,matid,*geo->Mesh());
		pt0d->SetSubElement(/*i*/0 , npt);
	//}

  SubElVec.Resize(NSubEl);
  for(sub=0;sub<NSubEl;sub++) {
    SubElVec[sub] = geo->SubElement(sub);
    SubElVec[sub]->SetFather(geo);
  }
  //for(i=0;i<NSubEl;i++) {//conectividades entre os filhos : viz interna
  //  for(j=0;j<1;j++) {        //lado do subel                                          numero do filho viz.             lado do viz.
  //    geo->SubElement(i)->SetNeighbour(InNeigh[i][j][0],TPZGeoElSide(geo->SubElement(InNeigh[i][j][1]),InNeigh[i][j][2]));
  //  }
  //}

  geo->SetSubElementConnectivities();

}
void TPZRefPoint::MidSideNodeIndex(TPZGeoEl *gel,int side,int &index){
  index = -1;
  if(side != 0) {
    PZError << "TPZRefCube::MidSideNodeIndex. Bad parameter side = " << side << endl;
    return;
  }
  //sides 0
  //if(side<TPZShapeLinear::NNodes) {//o nó medio do lado 0 é o 0 etc.
    index = (gel)->NodeIndex(side);
    return; 
  //}
  //o nó medio da face é o centro e o nó medio do centro é o centro
  //como nó de algum filho se este existir
  //caso tenha filhos é o canto de algum filho, se não tiver filhos retorna -1
//  if(gel->HasSubElement()) {
//	  side-=TPZShapeLinear::NNodes;
//    index=(gel->SubElement(MidSideNodes[side][0]))->NodeIndex(MidSideNodes[side][1]);
//  }
}
void TPZRefPoint::NewMidSideNode(TPZGeoEl *gel,int side,int &index) {

  MidSideNodeIndex(gel,side,index);
  if(index < 0) {
    TPZGeoElSide gelside = gel->Neighbour(side);
	if(gelside.Element()) {
	    while(gelside.Element() != gel) {
			gelside.Element()->MidSideNodeIndex(gelside.Side(),index);
			if(index!=-1) return;
			gelside = gelside.Neighbour();
	    }
	}
    TPZVec<REAL> par(3,0.);
    TPZVec<REAL> coord(3,0.);
    if(side == 0) {
		index = gel->NodeIndex(side); 
		return;
	}
    //aqui side = 8 a 26
//    side-=TPZShapeLinear::NNodes;//0,1,..,18
//    par[0] = MidCoord[side][0];
//    gel->X(par,coord);
//    index = gel->Mesh()->NodeVec().AllocateNewElement();
//    gel->Mesh()->NodeVec()[index].Initialize(coord,*gel->Mesh());
  }

}

void TPZRefPoint::GetSubElements(TPZGeoEl *father,int side, TPZStack<TPZGeoElSide> &subel) {
  subel.Resize(0);
  if(side<0 || side>1 || !father->HasSubElement()){
    PZError << "TPZRefPoint::GetSubElements called with bad arguments\n";
    return;
  }
  int nsub = NSideSubElements(side);//nsubeldata[side];
  for(int i=0;i<nsub;i++)
    subel.Push(TPZGeoElSide(father->SubElement(subeldata[side][i][0]),
												subeldata[side][i][1]));

}
int TPZRefPoint::NSideSubElements(int side) {
	if(side<0 || side>1){
		PZError << "TPZRefPoint::NSideSubelements called with error arguments\n";
		return -1;
	}
	return nsubeldata[side];
}

TPZTransform TPZRefPoint::GetTransform(int side,int whichsubel) {
	return TPZTransform (0,0);
	/*if(side<0 || side>1){
		PZError << "TPZRefPoint::GetTransform bad side\n";
		return TPZTransform(0,0);
	}
	int smalldim = 0;//TPZShapeLinear::SideDimension(side);
	int fatherside = FatherSide(side,whichsubel);
	int largedim = 0;//TPZShapeLinear::SideDimension(fatherside);
	TPZTransform trans(largedim,smalldim);
	int i,j;
	for(i=0; i<largedim; i++) {
		for(j=0; j<smalldim; j++) {
			trans.Mult()(i,j) = buildt[whichsubel][side][j][i];
		}
		trans.Sum() (i,0) = buildt[whichsubel][side][3][i];
	}
	return trans;
	*/
}

int TPZRefPoint::FatherSide(int side,int whichsubel) {
	if(side<0 || side>1){
		PZError << "TPZRefPoint::FatherSide called error" << endl;
		return -1;
	}
	return 0;//fatherside[whichsubel][side];
}
