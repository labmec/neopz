#include "pzshapelinear.h"
#include "pzerror.h"
#include "pzreal.h"

static int nhighdimsides[3] = {1,1,0};

static int highsides[3][1] = {
{2},
{2},
{0}
};

static int sidedimension[3] = {0,0,1};

static REAL sidetosidetransforms[3][1][4][3] = {
{
{{0,0,0},{0,0,0},{0,0,0},{-1,0,0}}
},
{
{{0,0,0},{0,0,0},{0,0,0},{1,0,0}}
},
{
{{0,0,0},{0,0,0},{0,0,0},{0,0,0}}}
};

static REAL MidSideNode[3][1] = {{-1.},{1.},{0.}};

void TPZShapeLinear::Chebyshev(REAL x,int num,TPZFMatrix &phi,TPZFMatrix &dphi){
  // Quadratic or higher shape functions
  if(num <= 0) return;
  phi.Put(0,0,1.0);
  dphi.Put(0,0, 0.0);
  if(num == 1) return;
  phi.Put(1,0, x);
  dphi.Put(0,1, 1.0);
  int ord;
  for(ord = 2;ord<num;ord++) {
    phi.Put(ord,0, 2.0*x*phi(ord-1,0) - phi(ord-2,0));
    dphi.Put(0,ord, 2.0*x*dphi(0,ord-1) + 2.0*phi(ord-1,0) - dphi(0,ord-2));
  }
}

void (*TPZShapeLinear::fOrthogonal)(REAL, int, TPZFMatrix &, TPZFMatrix &) = TPZShapeLinear::Chebyshev/*(REAL ,int ,TPZFMatrix& ,TPZFMatrix& )*/;//Chebyshev;

//  REAL TPZCompEl::gTrans1d[2] = {1.,-1.};

void TPZShapeLinear::Shape1d(REAL x,int order,TPZFMatrix &phi,TPZFMatrix &dphi,TPZVec<int> &id) {
  //	num = number of functions to compute
  if ( order < 0 ) {
    PZError << "Compelbas::shape --> Invalid dimension for arguments: order = " << order
	    << " phi.Rows = " << (int) phi.Rows() << " dphi.Cols = " << (int) dphi.Cols() << "\n";
    return;
  }
  if(phi.Rows() < order+1) phi.Resize(order, phi.Cols());
  if(dphi.Cols() < order+1) dphi.Resize(dphi.Rows(),order);

  if ( order == 0) {
    phi(0,0) = 1.;
    dphi(0,0) = 0.;
  } else if (order == 1) {		// Linear shape functions
    phi.Put(0,0, (1-x)/2.);
    phi.Put(1,0, (1+x)/2.);
    dphi.Put(0,0, -0.5);
    dphi.Put(0,1, 0.5);
    return;
  }

  // Quadratic or higher shape functions
  int num2 = order-1;
  int transformationindex = GetTransformId1d(id);
  if(num2 > 0) Shape1dInternal(x,num2,phi,dphi,transformationindex);
  int ord;
  for (ord = order; ord > 1; ord--) {
    phi(ord,0) = phi(ord-2,0);
    dphi(0,ord) = dphi(0,ord-2);
  }
  phi(0,0) = (1-x)/2.;
  phi(1,0) = (1+x)/2.;
  dphi(0,0) = -0.5;
  dphi(0,1) = 0.5;
  for (ord = 2; ord < order+1; ord++) {    // even functions
    dphi(0,ord) *= phi(0,0)*phi(1,0);
    dphi(0,ord) += 0.5*phi(0,0)*phi(ord,0)-0.5*phi(1,0)*phi(ord,0);
    phi(ord,0) *= phi(0,0)*phi(1,0);
  }
}

void TPZShapeLinear::Shape1dInternal(REAL x,int num,TPZFMatrix &phi,TPZFMatrix &dphi,int transformation_index){
  // Quadratic or higher shape functions
  if(num <= 0) return;
  REAL y;
  TransformPoint1d(transformation_index,x,y);
  fOrthogonal(y,num,phi,dphi);
  TransformDerivative1d(transformation_index,num,dphi);
}

void TPZShapeLinear::TransformPoint1d(int transid,REAL in,REAL &out) {
  if (!transid) out =  in;
  else          out = -in;
}

void TPZShapeLinear::TransformDerivative1d(int transid,int num,TPZFMatrix &in) {

  if(transid == 0) return;
  int i;
  for(i=0;i<num;i++) in(0,i) = -in(0,i);
}

int TPZShapeLinear::GetTransformId1d(TPZVec<int> &id) {
  if (id[1] < id[0]) return 1;
  else               return 0;
}

int TPZShapeLinear::NConnectShapeF(int side, TPZVec<int> &order) {
   if(side<2) return 1;//0 a 4
   if(side<3) return (order[0]-1);//6 a 14
   PZError << "TPZShapeLinear::NConnectShapeF, bad parameter side " << side << endl;
   return 0;
}

int TPZShapeLinear::NShapeF(TPZVec<int> &order) {
  int nn = NConnects();
  int in,res=0;
  for(in=0;in<nn;in++) res += NConnectShapeF(in,order);
  return res;
}

int TPZShapeLinear::NConnects() {
	return 3;
}

TPZTransform TPZShapeLinear::TransformElementToSide(int side){

	if(side<0 || side>2){
  	PZError << "TPZShapeLinear::TransformElementToSide called with side error\n";
    return TPZTransform(0,0);
  }

//  int sidedim = SideDimension(side);
  TPZTransform t(0,1);//t(dimto,2)
  t.Mult().Zero();	//TPZGeoElQ2d *gq;
  t.Sum().Zero();//int dimto = gq->SideDimension(side);

  switch(side){
    case 0:
    case 1:
      return t;
    case 2:
    	t.Mult()(0,0) = 1.0;//par. var.
      return t;
  }
  return TPZTransform(0,0);
}

TPZTransform TPZShapeLinear::TransformSideToElement(int side){

	if(side<0 || side>2){
  		PZError << "TPZShapeLinear::TransformSideToElement side out range\n";
		return TPZTransform(0,0);
	}
	int sidedim = 1;
	if(side <2) sidedim = 0;

  TPZTransform t(1,sidedim);
  t.Mult().Zero();
  t.Sum().Zero();

  switch(side){
  	case 0:
    	t.Sum()(0,0) = -1.0;
      return t;
    case 1:
      t.Sum()(0,0) =  1.0;
      return t;
    case 2:
    	t.Mult()(0,0) =  1.0;
      return t;
  }
	return TPZTransform(0,0);
}

int TPZShapeLinear::SideNodeLocId(int side, int node)
{
	if(side <2 && node == 0) return side;
	if(side == 2 && node <2) return node;
	PZError << "TPZShapeLinear::SideNodeLocId inconsistent side or node " << side
		<< ' ' << node << endl;
	return -1;
}

static int nsidenodes[3] = {1,1,2};
int TPZShapeLinear::NSideNodes(int side)
{
	return nsidenodes[side];
}

void TPZShapeLinear::HigherDimensionSides(int side, TPZStack<int> &high)
{
	if(side <0 || side >= NSides) {
		PZError << "TPZShapeLinear::HigherDimensionSides side "<< side << endl;
	}
	int is;
	for(is=0; is<nhighdimsides[side]; is++) high.Push(highsides[side][is]);
	
}

TPZTransform TPZShapeLinear::SideToSideTransform(int sidefrom, int sideto)
{
	if(sidefrom <0 || sidefrom >= NSides || sideto <0 || sideto >= NSides) {
		PZError << "TPZShapeLinear::HigherDimensionSides sidefrom "<< sidefrom << 
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
	PZError << "TPZShapeLinear::SideToSideTransform highside not found sidefrom "
		<< sidefrom << ' ' << sideto << endl;
	return TPZTransform(0);
}

int TPZShapeLinear::SideDimension(int side) {
	if(side<0 || side >= NSides) {
		PZError << "TPZShapeLinear::SideDimension side " << side << endl;
		return -1;
	}
	return sidedimension[side];
}

int TPZShapeLinear::NSideConnects(int i) {
  if(i==0 || i==1) return 1;
  else if(i==2) return 3;
  PZError << "TPZShapeLinear::NSideConnects. Bad parameter i = " << i << " .\n";
  return 0;
}

int TPZShapeLinear::SideConnectLocId(int side,int c) {
  switch(side) {
  case 0:
    if(c != 0)
      PZError << "TPZShapeLinear::SideConnectLocId, connect = " << c << endl;
    return 0;
  case 1:
    return 1;
  case 2:
    return c;
  default:
    PZError << "TPZShapeLinear::SideConnectLocId called with side = " << side << endl;
    return 0;
  }
}
void TPZShapeLinear::CenterPoint(int side, TPZVec<REAL> &center) {
  //center.Resize(Dimension);
  int i;
  for(i=0; i<Dimension; i++) {
    center[i] = MidSideNode[side][i];
  }
}

#ifdef _AUTODIFF
void TPZShapeLinear::Shape1dInternal(FADREAL & x,int num,TPZVec<FADREAL> & phi,int transformation_index){
  // Quadratic or higher shape functions
  if(num <= 0) return;
  FADREAL y;
  TransformPoint1d(transformation_index,x,y);
  FADfOrthogonal(y,num,phi);
//  TransformDerivative1d(transformation_index,num,phi);
}

void TPZShapeLinear::TransformPoint1d(int transid,FADREAL & in,FADREAL &out) {
  if (!transid) out =  in;
  else          out = -in;
}

void TPZShapeLinear::Chebyshev(FADREAL & x,int num,TPZVec<FADREAL> &phi){
  // Quadratic or higher shape functions
  if(num <= 0) return;
  //phi.Put(0,0,1.0);
  //dphi.Put(0,0, 0.0);
  phi[0] = 1.0; // <!> Remark: the derivatives other than the 0th are set to null
  if(num == 1) return;
  //phi.Put(1,0, x);
  //dphi.Put(0,1, 1.0);
  phi[1] = x;
  //phi[1].fastAccessDx(0)=1.0; // <!> Remark: the derivatives other than the 0th aren't set to null
  //just ensuring the derivatives are properly initialized and that FAD objects of more than
  // one derivative are used
  int ord;
  for(ord = 2;ord<num;ord++) {
    //phi.Put(ord,0, 2.0*x*phi(ord-1,0) - phi(ord-2,0));
    //dphi.Put(0,ord, 2.0*x*dphi(0,ord-1) + 2.0*phi(ord-1,0) - dphi(0,ord-2));
    phi[ord] = x * phi[ord-1] * 2.0 - phi[ord-2];
  }
}

void (*TPZShapeLinear::FADfOrthogonal)(FADREAL&,int ,TPZVec<FADREAL> &) =  TPZShapeLinear::Chebyshev/*(FADREAL&, int, TPZVec<FADREAL>&)*/;//Chebyshev;
/*
void TPZShapeLinear::TransformDerivative1d(int transid,int num,TPZVec<FADREAL> &in) {

  if(transid == 0) return;
  int i;
  for(i=0;i<num;i++) in[i].fastAccessDx(0) = -in[i].d(0);//(0,i) = -in(0,i);
}
*/
#endif
