#include "pzshapetriang.h"
#include "pzshapelinear.h"
#include "pzshapepoint.h"
#include "pzmanvector.h"
#include "pzerror.h"
#include "pzreal.h"


/**Transformation of the point within a triangular face */
REAL TPZShapeTriang::gTrans2dT[6][2][2] = {//s* , t*
  { { 1., 0.},{ 0., 1.} },
  { { 0., 1.},{ 1., 0.} },
  { { 0., 1.},{-1.,-1.} },//s* = t   t* = -s-t-1 ,  etc
  { {-1.,-1.},{ 0., 1.} },
  { {-1.,-1.},{ 1., 0.} },
  { { 1., 0.},{-1.,-1.} }
};

REAL TPZShapeTriang::gVet2dT[6][2] = {  {0.,0.},{0.,0.},{0.,1.},{1.,0.},{1.,0.},{0.,1.} };

REAL TPZShapeTriang::gRibTrans2dT1d[3][2] = { {2.,1.},{-1.,1.},{-1.,-2.} };//Cedric : 06/03/99

REAL TPZShapeTriang::gVet1dT[3] = {-1.,0.,1.};

static int sidedimension[7] = {0,0,0,1,1,1,2};

static int nhighdimsides[7] = {3,3,3,1,1,1,0};

static int highsides[7][3] = {
{3,5,6},
{3,4,6},
{4,5,6},
{6},
{6},
{6},
{-999}
};

static REAL sidetosidetransforms[7][3][4][3] = {
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,0,-99}}
},
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}}
},
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}}
},
{
{{0.5,0,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0,-99}}
},
{
{{-0.5,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0.5,-99}}
},
{
{{0,-0.5,-99},{-99,-99,-99},{-99,-99,-99},{0,0.5,-99}}
},
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-99,-99,-99}}
}
};

static REAL MidSideNode[7][3] = {
/*00*/{.0,0.},/*01*/{1.0,.0},/*02*/{0.,1.0},
/*03*/{.5,0.},/*04*/{0.5,.5},/*05*/{0.,0.5},
/*06*/{ 1./3.,1./3.} };

void TPZShapeTriang::ShapeCorner(TPZVec<REAL> &pt, TPZFMatrix &phi, TPZFMatrix &dphi) {

  phi(0,0) =  1.-pt[0]-pt[1];
  phi(1,0) =  pt[0];
  phi(2,0) =  pt[1];
  dphi(0,0) = -1.;
  dphi(1,0) = -1.;
  dphi(0,1) =  1.;
  dphi(1,1) =  0.;
  dphi(0,2) =  0.;
  dphi(1,2) =  1.;
}

void TPZShapeTriang::Shape(TPZVec<REAL> &pt, TPZVec<int> &id, TPZVec<int> &order,
			    TPZFMatrix &phi,TPZFMatrix &dphi) {
  ShapeCorner(pt,phi,dphi);
  if (order[0] < 2 && order[1] < 2 && order[2] < 2) return;
  REAL out;
  int shape = 3;
  for (int rib = 0; rib < 3; rib++) {

    ProjectPoint2dTriangToRib(rib,pt,out);
    TPZManVector<REAL,1> outvec(1,out);
    TPZVec<int> ids(2);
    ids[0] = id[rib%3];
    ids[1] = id[(rib+1)%3];
    REAL store1[20],store2[40];
    int ord2 = order[rib]-1;//numero de shapes por lado rib
    TPZFMatrix phin(ord2,1,store1,20),dphin(2,ord2,store2,40);
    TPZShapeLinear *shplin=0;
    shplin->ShapeInternal(outvec,order[rib],phin,dphin,shplin->GetTransformId1d(ids));
    TransformDerivativeFromRibToTriang(rib,ord2,dphin);
    for (int i = 0; i < ord2; i++) {
      phi(shape,0) = phi(rib,0)*phi((rib+1)%3,0)*phin(i,0);
      for(int xj=0;xj<2;xj++) {
           dphi(xj,shape) = dphi(xj,rib)* phi((rib+1)%3, 0 )* phin( i, 0)+
                            phi(rib, 0 )*dphi(xj,(rib+1)%3) * phin( i, 0)+
                            phi(rib, 0 )* phi((rib+1)%3, 0 )*dphin(xj,i);
      }
      shape++;
    }
  }
  if (order[3] < 3) return;//ordem na face
  REAL store1[20],store2[40];
  int ord =  order[3]-2;//num de shapes da face
  int nsh = (ord*(ord+1))/2;
  TPZFMatrix phin(nsh,1,store1,20),dphin(2,nsh,store2,40);
  ShapeInternal(pt,order[3]-2,phin,dphin,GetTransformId2dT(id));
  for(int i=0;i<nsh;i++)	{//number of internal shape equal maximal order
    phi(shape,0) = phi(0,0)*phi(1,0)*phi(2,0)*phin(i,0);
    for(int d=0;d<2;d++) {
      dphi(d,shape) = dphi(d,0)* phi(1,0)* phi(2,0)* phin(i,0) +
                       phi(0,0)*dphi(d,1)* phi(2,0)* phin(i,0) +
                       phi(0,0)* phi(1,0)*dphi(d,2)* phin(i,0) +
                       phi(0,0)* phi(1,0)* phi(2,0)*dphin(d,i);
    }
    shape++;
  }
}

void TPZShapeTriang::SideShape(int side, TPZVec<REAL> &pt, TPZVec<int> &id, TPZVec<int> &order,
			    TPZFMatrix &phi,TPZFMatrix &dphi) {
  if(side<0 || side>6) PZError << "TPZShapeTriang::SideShape. Bad paramenter side.\n";
  else if(side==6) Shape(pt,id,order,phi,dphi);
  else if(side<3) {
    TPZShapePoint::Shape(pt,id,order,phi,dphi);
  } else {
    TPZShapeLinear::Shape(pt,id,order,phi,dphi);
  }


}
void TPZShapeTriang::ShapeInternal(TPZVec<REAL> &x, int order,TPZFMatrix &phi,
					TPZFMatrix &dphi,int triangle_transformation_index) {

  if(order < 0) return;
  int ord1 = order;
  int numshape = (ord1*(ord1+1))/2;
  TPZManVector<REAL> out(2);
  TransformPoint2dT(triangle_transformation_index,x,out);

  if (phi.Rows() < numshape || dphi.Cols() < numshape) {
    PZError << "\nTPZCompEl::Shape2dTriangleInternal phi or dphi resized\n";
    phi.Resize(numshape,1);
    dphi.Resize(dphi.Rows(),numshape);
  }
  REAL store1[20],store2[20],store3[20],store4[20];
  TPZFMatrix phi0(numshape,1,store1,20),phi1(numshape,1,store2,20),dphi0(1,numshape,store3,20),dphi1(1,numshape,store4,20);

  TPZShapeLinear::fOrthogonal(2.*out[0]-1.,numshape,phi0,dphi0);
  TPZShapeLinear::fOrthogonal(2.*out[1]-1.,numshape,phi1,dphi1);
  int index = 0;
  int i;
  for (int iplusj=0;iplusj<ord1;iplusj++) {
    for (int j=0;j<=iplusj;j++) {
    	i = iplusj-j;
      phi(index,0) = phi0(i,0)*phi1(j,0);
      dphi(0,index) = 2.*dphi0(0,i)*phi1(j,0);
      dphi(1,index) = 2.*phi0(i,0)*dphi1(0,j);
      index++;
    }
  }
  TransformDerivative2dT(triangle_transformation_index,numshape,dphi);
}

void TPZShapeTriang::ProjectPoint2dTriangToRib(int rib, TPZVec<REAL> &in, REAL &out) {

  out = gRibTrans2dT1d[rib][0]*in[0]+gRibTrans2dT1d[rib][1]*in[1]+gVet1dT[rib];
}

void TPZShapeTriang::TransformDerivativeFromRibToTriang(int rib,int num,TPZFMatrix &dphi) {

  for (int j = 0;j<num;j++) {

    dphi(1,j) = gRibTrans2dT1d[rib][1]*dphi(0,j);
    dphi(0,j) = gRibTrans2dT1d[rib][0]*dphi(0,j);
  }
}

int TPZShapeTriang::GetTransformId2dT(TPZVec<int> &id) {

  int id0,id1,minid;
  id0 = (id[0] < id[1]) ? 0 : 1;
  minid = (id[2] < id[id0]) ? 2 : id0;
  id0 = (minid+1)%3;
  id1 = (minid+2)%3;

  if (id[id0] < id[id1]) {//antihorario

    if (minid == 0) return 0;
    if (minid == 1) return 2;
    if (minid == 2) return 4;

  } else {//horario

    if (minid == 0) return 1;
    if (minid == 1) return 3;
    if (minid == 2) return 5;
  }
  return 0;
}

//transf. o ponto dentro da face triangular
void TPZShapeTriang::TransformPoint2dT(int transid, TPZVec<REAL> &in, TPZVec<REAL> &out) {

  out[0] = gTrans2dT[transid][0][0]*in[0]+gTrans2dT[transid][0][1]*in[1]+gVet2dT[transid][0];
  out[1] = gTrans2dT[transid][1][0]*in[0]+gTrans2dT[transid][1][1]*in[1]+gVet2dT[transid][1];
}

void TPZShapeTriang::TransformDerivative2dT(int transid, int num, TPZFMatrix &in) {

  int i;
  for(i=0;i<num;i++) { //ds/dcsi
    REAL aux[2];
    aux[0] = in(0,i);
    aux[1] = in(1,i);
    in(0,i) = gTrans2dT[transid][0][0]*aux[0]+gTrans2dT[transid][1][0]*aux[1];
    in(1,i) = gTrans2dT[transid][0][1]*aux[0]+gTrans2dT[transid][1][1]*aux[1];
  }
}

int TPZShapeTriang::NConnectShapeF(int side, TPZVec<int> &order) {
  switch(side) {
  case 0:
  case 1:
  case 2:
    return 1;
  case 3:
  case 4:
  case 5:
    return order[side-3]-1;
  case 6:
    return (order[3]-2) < 0 ? 0 : ((order[3]-2)*(order[3]-1))/2;
  default:
    PZError << "TPZShapeTriang::NConnectShapeF, bad parameter iconnect " << side << endl;
    return 0;
  }
}

int TPZShapeTriang::NShapeF(TPZVec<int> &order) {
  int nn = NConnects();
  int in,res=0;
  for(in=0;in<nn;in++) res += NConnectShapeF(in,order);
  return res;
}

int TPZShapeTriang::NConnects() {
	return 7;
}

TPZTransform TPZShapeTriang::TransformElementToSide(int side){

	if(side<0 || side>6){
  	PZError << "TPZShapeTriang::TransformElementToSide called with side error\n";
    return TPZTransform(0,0);
  }

  TPZTransform t(sidedimension[side],2);//t(dimto,2)
  t.Mult().Zero();
  t.Sum().Zero();

  switch(side){
    case 0:
    case 1:
    case 2:
      return t;
    case 3:
    	t.Mult()(0,0) =  2.0;//par. var.
      t.Sum()(0,0)  = -1.0;
      return t;
    case 4:
      t.Mult()(0,0) = -1.0;
      t.Mult()(0,1) =  1.0;
      return t;
    case 5:
    	t.Mult()(0,1) = -2.0;
      t.Sum()(0,0)  =  1.0;
      return t;
    case 6:
      t.Mult()(0,0) =  1.0;
    	t.Mult()(1,1) =  1.0;
      return t;
  }
  return TPZTransform(0,0);
}

TPZTransform TPZShapeTriang::TransformSideToElement(int side){

	if(side<0 || side>6){
  	PZError << "TPZShapeTriang::TransformSideToElement side out range\n";
    return TPZTransform(0,0);
  }
  TPZTransform t(2,sidedimension[side]);
  t.Mult().Zero();
  t.Sum().Zero();

  switch(side){
  	case 0:
      return t;
    case 1:
    	t.Sum()(0,0) =  1.0;
      return t;
    case 2:
      t.Sum()(1,0) =  1.0;
      return t;
    case 3:
    	t.Mult()(0,0) =  0.5;
      t.Sum()(0,0)  =  0.5;
      return t;
    case 4:
      t.Mult()(0,0) = -0.5;
      t.Mult()(1,0) =  0.5;
      t.Sum() (0,0) =  0.5;
      t.Sum() (1,0) =  0.5;
      return t;
    case 5:
      t.Mult()(1,0) = -0.5;
      t.Sum() (1,0) =  0.5;
      return t;
    case 6:
      t.Mult()(0,0) =  1.0;
      t.Mult()(1,1) =  1.0;
      return t;
  }
	return TPZTransform(0,0);
}

int TPZShapeTriang::SideNodeLocId(int side, int node)
{
	if(side<3 && node == 0) return side;
	if(side>=3 && side<6 && node <2) return (side-3+node) %3;
	if(side==6 && node <3) return node;
	PZError << "TPZShapeTriang::SideNodeLocId inconsistent side or node " << side
		<< ' ' << node << endl;
	return -1;
}

static int nsidenodes[7] = {1,1,1,2,2,2,3};

int TPZShapeTriang::NSideNodes(int side)
{
	return nsidenodes[side];
}

void TPZShapeTriang::HigherDimensionSides(int side, TPZStack<int> &high)
{
	if(side <0 || side >= NSides) {
		PZError << "TPZShapeTriang::HigherDimensionSides side "<< side << endl;
	}
	int is;
	for(is=0; is<nhighdimsides[side]; is++) high.Push(highsides[side][is]);
	
}

TPZTransform TPZShapeTriang::SideToSideTransform(int sidefrom, int sideto)
{
	if(sidefrom <0 || sidefrom >= NSides || sideto <0 || sideto >= NSides) {
		PZError << "TPZShapeTriang::HigherDimensionSides sidefrom "<< sidefrom << 
			' ' << sideto << endl;
		return TPZTransform(0);
	}
	if(sidefrom == sideto) {
		return TPZTransform(sidedimension[sidefrom]);
	}
	if(sidefrom == NSides-1) {
		return TransformElementToSide(sideto);
	}
	int nhigh = nhighdimsides[sidefrom];
	int is;
	for(is=0; is<nhigh; is++) {
		if(highsides[sidefrom][is] == sideto) {
			int dfr = sidedimension[sidefrom];
			int dto = sidedimension[sideto];
			TPZTransform trans(dto,dfr);
			int i,j;
			for(i=0; i<dto; i++) {
				for(j=0; j<dfr; j++) {
					trans.Mult()(i,j) = sidetosidetransforms[sidefrom][is][j][i];
				}
				trans.Sum()(i,0) = sidetosidetransforms[sidefrom][is][3][i];
			}
			return trans;
		}
	}
	PZError << "TPZShapeTriang::SideToSideTransform highside not found sidefrom "
		<< sidefrom << ' ' << sideto << endl;
	return TPZTransform(0);
}

int TPZShapeTriang::SideDimension(int side) {
	if(side<0 || side >= NSides) {
		PZError << "TPZShapeTriang::SideDimension side " << side << endl;
		return -1;
	}
	return sidedimension[side];
}
int TPZShapeTriang::NSideConnects(int side) {
  if(side<0 || side>6) {
    PZError << "TPZShapeTriang::NSideConnects. Bad parameter i.\n";
    return 0;
  }
  if(side<3) return 1;
  if(side<6) return 3;
  return 7;
}

/**It do not verify the values of the c*/
int TPZShapeTriang::SideConnectLocId(int side, int c) {
  switch(side) {
  case 0:
  case 1:
  case 2:
    return side;
  case 3:
  case 4:
  case 5:
    if(!c) return side-3;
    if(c==1) return (side-2)%3;
    if(c==2) return side;
  case 6:
    return c;
  default:
    PZError << "TPZShapeTriang::SideConnectLocId, connect = " << c << endl;
    return -1;
  }
}

void TPZShapeTriang::CenterPoint(int side, TPZVec<REAL> &center) {
  //center.Resize(Dimension);
  int i;
  for(i=0; i<Dimension; i++) {
    center[i] = MidSideNode[side][i];
  }
}
