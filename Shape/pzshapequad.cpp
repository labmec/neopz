#include "pzshapequad.h"
#include "pzshapelinear.h"
//#include "pzelgq2d.h"
#include "pzmanvector.h"
#include "pzerror.h"
#include "pzreal.h"


/**Transformation of the point within a quadrilateral face */
REAL TPZShapeQuad::gTrans2dQ[8][2][2] = {//s* , t*
  { { 1., 0.},{ 0., 1.} },
  { { 0., 1.},{ 1., 0.} },
  { { 0., 1.},{-1., 0.} },
  { {-1., 0.},{ 0., 1.} },
  { {-1., 0.},{ 0.,-1.} },//s* = -s   t* = -t  , etc
  { { 0.,-1.},{-1., 0.} },
  { { 0.,-1.},{ 1., 0.} },
  { { 1., 0.},{ 0.,-1.} }
};

REAL TPZShapeQuad::gRibTrans2dQ1d[4][2] = { {1.,0.},{0.,1.},{-1.,0.},{0.,-1.} };

static int sidedimension[9] = {0,0,0,0,1,1,1,1,2};

static int nhighdimsides[9] = {3,3,3,3,1,1,1,1,0};

static int highsides[9][3] = {
{4,7,8},
{4,5,8},
{5,6,8},
{6,7,8},
{8},
{8},
{8},
{8},
{-999}
};

static REAL sidetosidetransforms[9][3][4][3] = {
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-1,-99}}
},
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-1,-99}}
},
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,1,-99}}
},
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,1,-99}}
},
{
{{1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,-1,-99}}
},
{
{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}}
},
{
{{-1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}}
},
{
{{0,-1,-99},{-99,-99,-99},{-99,-99,-99},{-1,0,-99}}
},
{
{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-99,-99,-99}}
}
};

static REAL MidSideNode[9][3] = {
/*00*/{-1.,-1.},/*01*/{ 1.,-1.},/*02*/{1.,1.},
/*03*/{-1., 1.},/*04*/{ 0.,-1.},/*05*/{1.,0.},
/*06*/{ 0., 1.},/*07*/{-1., 0.},/*08*/{0.,0.} };


void TPZShapeQuad::ShapeCornerQuad(TPZVec<REAL> &pt, TPZFMatrix &phi, TPZFMatrix &dphi) {

  REAL x[2],dx[2],y[2],dy[2];
  x[0]  =  (1.-pt[0])/2.;
  x[1]  =  (1.+pt[0])/2.;
  dx[0] = -0.5;
  dx[1] =  0.5;
  y[0]  =  (1.-pt[1])/2.;
  y[1]  =  (1.+pt[1])/2.;
  dy[0] = -0.5;
  dy[1] =  0.5;
  phi(0,0)  = x[0]*y[0];
  phi(1,0)  = x[1]*y[0];
  phi(2,0)  = x[1]*y[1];
  phi(3,0)  = x[0]*y[1];
  dphi(0,0) = dx[0]*y[0];
  dphi(1,0) = x[0]*dy[0];
  dphi(0,1) = dx[1]*y[0];
  dphi(1,1) = x[1]*dy[0];
  dphi(0,2) = dx[1]*y[1];
  dphi(1,2) = x[1]*dy[1];
  dphi(0,3) = dx[0]*y[1];
  dphi(1,3) = x[0]*dy[1];
}

void TPZShapeQuad::ShapeQuad(TPZVec<REAL> &pt, TPZVec<int> &id, TPZVec<int> &order,
			  TPZFMatrix &phi,TPZFMatrix &dphi) {
  ShapeCornerQuad(pt,phi,dphi);
  REAL out;
  int shape = 4;
  for (int rib = 0; rib < 4; rib++) {

    ProjectPoint2dQuadToRib(rib,pt,out);
    TPZVec<int> ids(2);
    ids[0] = id[rib%4];
    ids[1] = id[(rib+1)%4];
    REAL store1[20],store2[40];
    int ord2 = order[rib]-1;//two orders : order in x and order in y
    TPZFMatrix phin(ord2,1,store1,20),dphin(2,ord2,store2,40);
    TPZShapeLinear::Shape1dInternal(out,ord2,phin,dphin,TPZShapeLinear::GetTransformId1d(ids));
    TransformDerivativeFromRibToQuad(rib,ord2,dphin);
    for (int i = 0; i < ord2; i++) {
      phi(shape,0) = phi(rib,0)*phi((rib+1)%4,0)*phin(i,0);
      for(int xj=0;xj<2;xj++) {
        dphi(xj,shape) = dphi(xj,rib)*phi((rib+1)%4,0)*phin(i,0)+
                         phi(rib,0)*dphi(xj,(rib+1)%4)*phin(i,0)+
                         phi(rib,0)*phi((rib+1)%4,0)*dphin(xj,i);
      }
      shape++;
    }
  }
  REAL store1[20],store2[40];
  int ord = (order[4]-1)*(order[4]-1);
  TPZFMatrix phin(ord,1,store1,20),dphin(2,ord,store2,40);
  Shape2dQuadInternal(pt,order[4]-2,phin,dphin,GetTransformId2dQ(id));
  for(int i=0;i<ord;i++)	{//funcoes de interior são em numero ordem-1
    phi(shape,0) = phi(0,0)*phi(2,0)*phin(i,0);
    for(int xj=0;xj<2;xj++) {//x e y
      dphi(xj,shape) = dphi(xj,0)*phi(2,0)*phin(i,0) +
                       phi(0,0)*dphi(xj,2)*phin(i,0) +
                       phi(0,0)*phi(2,0)*dphin(xj,i);
    }
    shape++;
  }
}

void TPZShapeQuad::Shape2dQuadInternal(TPZVec<REAL> &x, int order,
					TPZFMatrix &phi, TPZFMatrix &dphi,int quad_transformation_index) {

  if(order < 0) return;
  int ord1 = order+1;
  int numshape = ord1*ord1;
  TPZManVector<REAL> out(2);
  TransformPoint2dQ(quad_transformation_index,x,out);

  if(numshape > phi.Rows() || phi.Cols() < 1) phi.Resize(numshape,1);
  if(dphi.Rows() < 2 || dphi.Cols() < numshape) dphi.Resize(2,numshape);
  REAL store1[20],store2[20],store3[20],store4[20];
  TPZFMatrix phi0(ord1,1,store1,20),phi1(ord1,1,store2,20),dphi0(1,ord1,store3,20),dphi1(1,ord1,store4,20);

  TPZShapeLinear::fOrthogonal(out[0],ord1,phi0,dphi0);
  TPZShapeLinear::fOrthogonal(out[1],ord1,phi1,dphi1);
  for (int i=0;i<ord1;i++) {
    for (int j=0;j<ord1;j++) {
      int index = i*ord1+j;
      phi(index,0) =  phi0(i,0)* phi1(j,0);
      dphi(0,index) = dphi0(0,i)* phi1(j,0);
      dphi(1,index) =  phi0(i,0)*dphi1(0,j);
    }
  }
  TransformDerivative2dQ(quad_transformation_index,numshape,dphi);
}

void TPZShapeQuad::TransformDerivative2dQ(int transid, int num, TPZFMatrix &in) {

  for(int i=0;i<num;i++) {
    double aux[2];
    aux[0] = in(0,i);
    aux[1] = in(1,i);
    in(0,i) = gTrans2dQ[transid][0][0]*aux[0]+gTrans2dQ[transid][1][0]*aux[1];
    in(1,i) = gTrans2dQ[transid][0][1]*aux[0]+gTrans2dQ[transid][1][1]*aux[1];
  }
}

//transf. o ponto dentro da face quadrilateral
void TPZShapeQuad::TransformPoint2dQ(int transid, TPZVec<REAL> &in, TPZVec<REAL> &out) {

  out[0] = gTrans2dQ[transid][0][0]*in[0]+gTrans2dQ[transid][0][1]*in[1];//Cedric 23/02/99
  out[1] = gTrans2dQ[transid][1][0]*in[0]+gTrans2dQ[transid][1][1]*in[1];//Cedric 23/02/99
}

void TPZShapeQuad::ProjectPoint2dQuadToRib(int rib, TPZVec<REAL> &in, REAL &out) {

  out = gRibTrans2dQ1d[rib][0]*in[0]+gRibTrans2dQ1d[rib][1]*in[1];
}

int TPZShapeQuad::GetTransformId2dQ(TPZVec<int> &id) {

  int id0,id1,minid;
  id0 = (id[0] < id[1]) ? 0 : 1;
  id1 = (id[2] < id[3]) ? 2 : 3;
  minid = (id[id0] < id[id1]) ? id0 : id1;//minid : menor id local
  id0 = (minid+1)%4;//id anterior local (sentido antihorario)
  id1 = (minid+3)%4;//id posterior local (sentido horario)
  minid = id[minid];//minid : menor id global

  if (id[id0] < id[id1]) {//antihorario

    if (minid == id[0]) return 0;
    if (minid == id[1]) return 2;
    if (minid == id[2]) return 4;
    if (minid == id[3]) return 6;

  } else {//horario

    if (minid == id[0]) return 1;
    if (minid == id[1]) return 3;
    if (minid == id[2]) return 5;
    if (minid == id[3]) return 7;
  }
  return 0;
}

void TPZShapeQuad::TransformDerivativeFromRibToQuad(int rib,int num,TPZFMatrix &dphi) {

  for (int j = 0;j<num;j++) {
    dphi(1,j) = gRibTrans2dQ1d[rib][1]*dphi(0,j);
    dphi(0,j) = gRibTrans2dQ1d[rib][0]*dphi(0,j);
  }
}

int TPZShapeQuad::NConnectShapeF(int side, TPZVec<int> &order) {
   if(side<4) return 1;//0 a 4
   int s = side-4;//s = 0 a 14 ou side = 6 a 20
   if(side<8) return (order[s]-1);//6 a 14
   if(side==8) {
      return ((order[s]-1)*(order[s]-1));
   }
   PZError << "TPZShapeQuad::NConnectShapeF, bad parameter side " << side << endl;
   return 0;
}

int TPZShapeQuad::NShapeF(TPZVec<int> &order) {
  int nn = NConnects();
  int in,res=0;
  for(in=0;in<nn;in++) res += NConnectShapeF(in,order);
  return res;
}
int TPZShapeQuad::NConnects() {
		return 9;
}

TPZTransform TPZShapeQuad::TransformElementToSide(int side){

	if(side<0 || side>8){
  	PZError << "TPZShapeQuad::TransformElementToSide called with side error\n";
    return TPZTransform(0,0);
  }

  TPZTransform t(sidedimension[side],2);//t(dimto,2)
  t.Mult().Zero();	//TPZGeoElQ2d *gq;
  t.Sum().Zero();//int dimto = gq->SideDimension(side);

  switch(side){
    case 0:
    case 1:
    case 2:
    case 3:
      return t;
    case 4:
    	t.Mult()(0,0) = 1.0;//par. var.
      return t;
    case 5 :
      t.Mult()(0,1) = 1.0;
      return t;
    case 6:
    	t.Mult()(0,0) = -1.0;
      return t;
    case 7:
    	t.Mult()(0,1) = -1.0;
      return t;
    case 8:
    	t.Mult()(0,0) = 1.0;
      t.Mult()(1,1) = 1.0;
      return t;
  }
  return TPZTransform(0,0);
}

TPZTransform TPZShapeQuad::TransformSideToElement(int side){

	if(side<0 || side>8){
  	PZError << "TPZShapeQuad::TransformSideToElement side out range\n";
    return TPZTransform(0,0);
  }
  TPZTransform t(2,sidedimension[side]);
  t.Mult().Zero();
  t.Sum().Zero();

  switch(side){
  	case 0:
    	t.Sum()(0,0) = -1.0;
      t.Sum()(1,0) = -1.0;
      return t;
    case 1:
    	t.Sum()(0,0) =  1.0;
      t.Sum()(1,0) = -1.0;
      return t;
    case 2:
    	t.Sum()(0,0) =  1.0;
      t.Sum()(1,0) =  1.0;
      return t;
    case 3:
    	t.Sum()(0,0) = -1.0;
      t.Sum()(1,0) =  1.0;
      return t;
    case 4:
      t.Mult()(0,0) =  1.0;
      t.Sum() (1,0) = -1.0;
      return t;
    case 5:
      t.Mult()(1,0) =  1.0;
      t.Sum() (0,0) =  1.0;
      return t;
    case 6:
      t.Mult()(0,0) = -1.0;
      t.Sum() (1,0) =  1.0;
      return t;
    case 7:
      t.Mult()(1,0) = -1.0;
      t.Sum() (0,0) = -1.0;
      return t;
    case 8:
      t.Mult()(0,0) =  1.0;
      t.Mult()(1,1) =  1.0;
      return t;
  }
	return TPZTransform(0,0);
}


int TPZShapeQuad::SideNodeLocId(int side, int node)
{
	if(side<4 && node==0) return side;
	if(side>=4 && side<8 && node <2) return (side+node)%4;
	if(side==8 && node <4) return node;
	PZError << "TPZShapeQuad::SideNodeLocId inconsistent side or node " << side
		<< ' ' << node << endl;
	return -1;
}

static int nsidenodes[9] = {
	1,1,1,1,2,2,2,2,4};

int TPZShapeQuad::NSideNodes(int side)
{
	return nsidenodes[side];
}

void TPZShapeQuad::HigherDimensionSides(int side, TPZStack<int> &high)
{
	if(side <0 || side >= NSides) {
		PZError << "TPZShapeQuad::HigherDimensionSides side "<< side << endl;
	}
	int is;
	for(is=0; is<nhighdimsides[side]; is++) high.Push(highsides[side][is]);
	
}

TPZTransform TPZShapeQuad::SideToSideTransform(int sidefrom, int sideto)
{
	if(sidefrom <0 || sidefrom >= NSides || sideto <0 || sideto >= NSides) {
		PZError << "TPZShapeQuad::HigherDimensionSides sidefrom "<< sidefrom << 
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
	PZError << "TPZShapeQuad::SideToSideTransform highside not found sidefrom "
		<< sidefrom << ' ' << sideto << endl;
	return TPZTransform(0);
}

int TPZShapeQuad::SideDimension(int side) {
	if(side<0 || side >= NSides) {
		PZError << "TPZShapeQuad::SideDimension side " << side << endl;
		return -1;
	}
	return sidedimension[side];
}

int TPZShapeQuad::NSideConnects(int side) {
  if(side<0 || side>8) {
    PZError << "TPZShapeQuad::NSideConnects. Bad parameter i.\n";
    return 0;
  }
  if(side<4) return 1;
  if(side<8) return 3;
  return 9;//Cedric
}

/**It do not verify the values of the c*/
int TPZShapeQuad::SideConnectLocId(int side,int c) {
  switch(side) {
  case 0:
  case 1:
  case 2:
  case 3:
    return side;
  case 4:
  case 5:
  case 6:
  case 7:
    if(!c) return side-4;
    if(c==1) return (side-3)%4;
    if(c==2) return side;
  case 8:
    return c;
  default:
    PZError << "TPZShapeQuad::SideConnectLocId, connect = " << c << endl;
    return -1;
  }
}

void TPZShapeQuad::CenterPoint(int side, TPZVec<REAL> &center) {
  //center.Resize(Dimension);
  int i;
  for(i=0; i<Dimension; i++) {
    center[i] = MidSideNode[side][i];
  }
}

#ifdef _AUTODIFF

void TPZShapeQuad::Shape2dQuadInternal(TPZVec<FADREAL> &x, int order,
					TPZVec<FADREAL> &phi,int quad_transformation_index) {

  const int ndim = 3;
  if(order < 0) return;
  int ord1 = order+1;
  int numshape = ord1*ord1;
  TPZVec<FADREAL> out(2);
  TransformPoint2dQ(quad_transformation_index,x,out);

  if(numshape > phi.NElements()/*Rows()*/ || phi[0].size()/*Cols()*/ < ndim) phi.Resize(numshape, FADREAL(ndim, 0.0));
  //if(dphi.Rows() < 2 || dphi.Cols() < numshape) dphi.Resize(2,numshape);
  //REAL store1[20],store2[20],store3[20],store4[20];
  //TPZFMatrix phi0(ord1,1,store1,20),phi1(ord1,1,store2,20),dphi0(1,ord1,store3,20),dphi1(1,ord1,store4,20);
  TPZVec<FADREAL> phi0(20, FADREAL(ndim, 0.0)),
                  phi1(20, FADREAL(ndim, 0.0));

  TPZShapeLinear::FADfOrthogonal(out[0],ord1,phi0);
  TPZShapeLinear::FADfOrthogonal(out[1],ord1,phi1);
  for (int i=0;i<ord1;i++) {
    for (int j=0;j<ord1;j++) {
      int index = i*ord1+j;
      //phi(index,0) =  phi0(i,0)* phi1(j,0);
      phi[index] =  phi0[i] * phi1[j];
      /*dphi(0,index) = dphi0(0,i)* phi1(j,0);
      dphi(1,index) =  phi0(i,0)*dphi1(0,j);*/
    }
  }
//  TransformDerivative2dQ(quad_transformation_index,numshape,phi);
}

void TPZShapeQuad::TransformPoint2dQ(int transid, TPZVec<FADREAL> &in, TPZVec<FADREAL> &out) {

  out[0] = gTrans2dQ[transid][0][0]*in[0]+gTrans2dQ[transid][0][1]*in[1];//Cedric 23/02/99
  out[1] = gTrans2dQ[transid][1][0]*in[0]+gTrans2dQ[transid][1][1]*in[1];//Cedric 23/02/99
}

/*
void TPZShapeQuad::TransformDerivative2dQ(int transid, int num, TPZVec<FADREAL> &in) {

  for(int i=0;i<num;i++) {
    double aux[2];
    aux[0] = in[i].d(0);
    aux[1] = in[i].d(1);
    in[i].fastAccessDx(0) = gTrans2dQ[transid][0][0]*aux[0]+gTrans2dQ[transid][1][0]*aux[1];
    in[i].fastAccessDx(1) = gTrans2dQ[transid][0][1]*aux[0]+gTrans2dQ[transid][1][1]*aux[1];
  }
}

void TPZShapeQuad::TransformDerivativeFromRibToQuad(int rib,int num,TPZVec<FADREAL> &phi) {

  for (int j = 0;j<num;j++) {
    //dphi(1,j) = gRibTrans2dQ1d[rib][1]*dphi(0,j);
    //dphi(0,j) = gRibTrans2dQ1d[rib][0]*dphi(0,j);
    phi[j].fastAccessDx(1) = gRibTrans2dQ1d[rib][1]*phi[j].d(0);
    phi[j].fastAccessDx(0) = gRibTrans2dQ1d[rib][0]*phi[j].d(0);
  }
}
*/

#endif
