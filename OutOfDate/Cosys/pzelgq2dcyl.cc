#include "pzelgq2dcyl.h"
#include "pzgmesh.h"

const REAL PI=3.1415926536;
 /**Constructors. Parameters: id - element id, nodeindexes - vector containing node indexes,
     matind - material index, refind - index of the node which indicates the reference direction*/

TPZGeoElQ2dCyl::TPZGeoElQ2dCyl(int id,TPZVec<int> &nodeindexes,int matind,TPZGeoMesh &mesh, int cosysindex):TPZGeoElQ2d(id,nodeindexes,matind,mesh){

	fCosys = mesh.CosysVec()[cosysindex];
	//  TPZCylinsys *cosys = mesh.CosysVec()[cosysindex];
	if (fCosys->Type() != cylindric && fCosys->Type() != esferic) {
		fCosysIndex = -1;
		return;
	}
	int i,j;
	for (i=0; i<4; i++){
		TPZGeoNode &gn = mesh.NodeVec()[nodeindexes[i]];
		TPZVec<REAL> co(3,0.); //,cotrans[3];
		for(j=0; j<3; j++) co[j]=gn.Coord(j);
		fCosys->FromReference(co);
		fX[i][0]=co[0] ;
		fX[i][1]=co[1];
		fX[i][2]=co[2];
		fCosysIndex=cosysindex;
	}
	TPZFMatrix wrap(4,3);
	for(i=0;i<4;i++)
		for(j=0;j<3;j++)
				wrap(i,j)=fX[i][j];

	fCosys->VerifyRange(wrap); 
	for(i=0;i<4;i++)
		for(j=0;j<3;j++)
				fX[i][j] = wrap(i,j);
}
     

TPZGeoElQ2dCyl::TPZGeoElQ2dCyl(TPZVec<int> &nodeindexes,int matind,TPZGeoMesh &mesh, int cosysindex):TPZGeoElQ2d(nodeindexes,matind,mesh){
	
	TPZCosys *cosys = mesh.CosysVec()[cosysindex];
	fCosys = mesh.CosysVec()[cosysindex];
	
	if (cosys->Type() != cylindric  && fCosys->Type() != esferic) {
		fCosysIndex = -1;
	}
	
	int i,j;
	
	for (i=0; i<4; i++){
		TPZGeoNode &gn = mesh.NodeVec()[nodeindexes[i]];
		TPZVec<REAL> co(3,0.);// ,cotrans[3];
		for(j=0; j<3; j++) co[j]=gn.Coord(j);
		cosys->FromReference(co);
		fX[i][0]=co[0] ;
		fX[i][1]=co[1];
		fX[i][2]=co[2];
		fCosysIndex=0;
	}
	TPZFMatrix wrap(4,3);
	for(i=0;i<4;i++)
		for(j=0;j<3;j++)
				wrap(i,j)=fX[i][j];

	fCosys->VerifyRange(wrap); 
	for(i=0;i<4;i++)
		for(j=0;j<3;j++)
				fX[i][j] = wrap(i,j);
}



TPZGeoElQ2dCyl::~TPZGeoElQ2dCyl() {
	if (fCosys) delete fCosys;
}


TPZGeoElQ2d *TPZGeoElQ2dCyl::CreateGeoEl(TPZVec<int> &nodeindexes,int matind,TPZGeoMesh &mesh) {
  
		return new TPZGeoElQ2dCyl(nodeindexes,matind,mesh,fCosysIndex);
}

/*
void TPZGeoElQ2dCyl::VerifyTheta(){
	int i,j,k=0;
	for(i=0;i<4;i++){
		for(j=i+1;j<4;j++){
			if ((fTheta[i]-fTheta[j])>PI || (fTheta[i]-fTheta[j])< -(PI) ){
			k=1;
			break;
			}
		}
		if(k) break;
	}
	if(k) {
		for (i=0;i<4;i++){
			if (fTheta[i]>PI) fTheta[i]=fTheta[i]-2.*PI;
		}
	}
}
*/

void TPZGeoElQ2dCyl::Jacobian(TPZVec<REAL> &param,TPZFMatrix &jacobian,TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv){

	int nnodes = NNodes();
#ifdef DEBUG
	if (nnodes != 4) {
		PZError << "TPZGeoElQ2dCyl.jacobian only implemented for"
			" 4, 8 or 9 nodes, NumberOfNodes = " << nnodes << "\n";
	}
	if(param.NElements() != 2 || param[0] < -1. || param[0] > 1. ||
		param[1] < -1. || param[1] > 1.) {
		PZError << "TPZGeoElQ2dCyl.jacobian. param out of range : "
			" param.NElements() = " << param.NElements() <<
			"\nparam[0] = " << param[0] << " param[1] = " << param[1] << "\n";
		return;
	}
#endif
  
	REAL spacephi[4];
	TPZFMatrix phi(4,1,spacephi,4);
	REAL spacedphi[8];
	TPZFMatrix dphi(2,4,spacedphi,8);
	Shape(param,phi,dphi);
	int i,j,k;
	for(i=0;i<2;i++)
		for(j=0;j<2;j++)
			jacobian(i,j)=0.;

	//TPZGeoNode *np;
	TPZVec<REAL> V1(3,0.),v1aux(3,0.),V2(3,0.),v2aux(3,0.),V2til(3,0.),V3(3,0.);
	REAL V1Norm=0.,V1V2=0.,V2tilNorm=0.;//,theta=0.,z=0.,r=0.;
	//R Theta e Z
	TPZVec<REAL> Xco(3,0.);
	for (i=0;i<4;i++){
		for(j=0; j<3;j++) {
			Xco[j] += fX[i][j] * phi(i,0);
		}
	}
  
	/**Parte a ser calculada pelo sistema de coordenada **/

// ===Cesar=== start grad rtz/gradxyzhat
	TPZVec<REAL> x(3,0.);
	TPZFMatrix gradX(3,2,0.);
	TPZFMatrix gradx(3,2,0.);
	
	//Cálculo de phi
	for (i=0;i<4;i++){
		for(j=0; j<3;j++) {
			for(k=0; k<2; k++) {
				gradX(j,k)+= fX[i][j] * dphi(k,i);
			}
		}
	}

	if (fCosys) fCosys->TransformGradient(Xco,gradX,x,gradx,NULL);
	else{
		x = Xco;
		gradx = gradX;
	}
	for (i=0;i<3;i++){
		V1[i]=gradx(i,0);
		V2[i]=gradx(i,1);
	}
	/** Fim do trecho a ser calculado pelo sistema de coordenadas **/
	for(i=0;i<3;i++) {
		V1Norm += V1[i]*V1[i];
		V1V2 += V1[i]*V2[i];
	}

	V1Norm = sqrt(V1Norm);
	for(i=0;i<3;i++) {
		V1[i] /= V1Norm;
		V2til[i] = V2[i] - V1V2*V1[i]/V1Norm;
		V2tilNorm += V2til[i]*V2til[i];
	}

	V2tilNorm = sqrt(V2tilNorm);
	jacobian(0,0) = V1Norm;
	jacobian(0,1) = V1V2/V1Norm;
	jacobian(1,1) = V2tilNorm;
  
	for(i=0;i<3;i++) {
		axes(0,i) = V1[i];
		axes(1,i) = V2til[i]/V2tilNorm;
	}
  
	detjac = jacobian(0,0)*jacobian(1,1)-jacobian(1,0)*jacobian(0,1);
	jacinv(0,0) = +jacobian(1,1)/detjac;
	jacinv(1,1) = +jacobian(0,0)/detjac;
	jacinv(0,1) = -jacobian(0,1)/detjac;
	jacinv(1,0) = -jacobian(1,0)/detjac;

	axes(2,0) = axes(0,1)*axes(1,2)-axes(0,2)*axes(1,1);
	axes(2,1) = -axes(0,0)*axes(1,2)+axes(0,2)*axes(1,0);
	axes(2,2) = axes(0,0)*axes(1,1)-axes(0,1)*axes(1,0);
}


/* return the cartesian coordinate in real element, given the master
   element point coordinate */ 
void TPZGeoElQ2dCyl::X(TPZVec<REAL> & loc,TPZVec<REAL> &result){
	REAL spacephi[4],spacedphi[8];
	int i,j;
//	REAL point[3] = {0.};
	TPZFMatrix phi(4,1,spacephi,4);
	TPZFMatrix dphi(2,4,spacedphi,8);
	Shape(loc,phi,dphi);

	
	/** O método X do TPZCosys calcula result dao phi e dphi **/

//	fCosys->X(phi,dphi,result);
	//X(phi,dphi,result);
	for(j=0; j<3;j++) {
		result[j]=0.;
		for (i=0;i<4;i++){
			result[j] += fX[i][j] * phi(i,0);
		}
	}
	
	if (fCosys) fCosys->ToGlobal(result);
}




/**It's necessary to define the normal vector to side 4, that is the orthogonal
   vector to the surface*/
void TPZGeoElQ2dCyl::NormalVector(int side,TPZVec<REAL> &param,TPZVec<REAL> &normal,TPZFMatrix &axes,TPZFMatrix &jac1d) {

	cout << "ERROR - old NormalVector method is not implemented. \n";
/*	int nnodes = NNodes();
#ifndef NODEBUG
  if (nnodes != 4) {
    PZError << "TPZGeoElQ2d.NormalVector, only implemented for"
      " 4 nodes, NumberOfNodes = " << nnodes << "\n";
  }
  if(param.NElements() != 2 || param[0] < -1. || param[0] > 1. ||
     param[1] < -1. || param[1] > 1.) {
    PZError << "TPZGeoElQ2d.NormalVector, fl out of range : "
      " point.NElements() = " << param.NElements() <<
      "\npoint[0] = " << param[0] << " point[1] = " << param[1] << "\n";
    return;
  }
  if(normal.NElements() != 3) {
    PZError << "elgq2d.NormalVector normal.capacity() = " << normal.NElements() <<
      "\n";
    return;
  }
  if(side < 0 || side >= 4) {
    PZError << "TPZGeoElQ2d.jacobian invalid side : "
      " side = " << side << "\n";
    return;
  }
#endif

  REAL spacephi[4],spacedphi[8];
//  TPZFMatrix phi(4,1,spacedphi,4);
// Philippe 31;3;99
  TPZFMatrix phi(4,1,spacephi,4);
  TPZFMatrix dphi(2,4,spacedphi,8);
  Shape(param,phi,dphi);
  //====Rem by Cesar 18/05/00====
  //  TPZGeoNode *np;
  TPZVec<REAL> t(3,0.);
  int i,j,ider = 0;
  if(side==1 || side==3) ider = 1;

  for(i=0;i<nnodes;i++) {
    //======Cesar=======   
    //    np = NodePtr(i);
    REAL point[3];
    point[0]=fX[i];
    point[1]=fTheta[i];
    point[2]=fZ[i];
    fCosys->ToCart(point);
    for(j=0;j<3;j++){
      //      t[j] += np->Coord(j)*dphi(ider,i);
      t[j] += point[j]*dphi(ider,i);
    }
  }
  //=====End Cesar

  //      note that ||t|| != 1 , ||t|| = |J|

  jac1d(0,0) = sqrt( t[0]*t[0] + t[1]*t[1] + t[2]*t[2] );

  // consistent axes computation Philippe 17/4/97
  TPZVec<REAL> V1(3,0.),V2(3,0.),V2til(3,0.),V3(3,0.),V1til(3,0.);
  REAL V1Norm=0.,V2Norm=0.,V1V2=0.,V2tilNorm=0.,V1tilNorm =0.;
  for(i=0;i<nnodes;i++) {
    //=====Cesar======
    //    np = NodePtr(i);
    REAL point[3];
    point[0]=fRadius[i];
    point[1]=fTheta[i];
    point[2]=fZ[i];
    fCosys->ToCart(point);
    for(j=0;j<3;j++) {
      //      V1[j] += np->Coord(j)*dphi(0,i);
      //      V2[j] += np->Coord(j)*dphi(1,i);
      V1[j] += point[j]*dphi(0,i);
      V2[j] += point[j]*dphi(1,i);
      //=====End Cesar=====
    }
  }
  for(j=0;j<3;j++) {
    V1Norm += V1[j]*V1[j];
    V2Norm += V2[j]*V2[j];
    V1V2 += V1[j]*V2[j];
  }
  V1Norm = sqrt(V1Norm);
  V2Norm = sqrt(V2Norm);
  for(j=0;j<3;j++) {
    V1[j] /= V1Norm;
    V2[j] /= V2Norm;
    V2til[j] = V2[j] - V1V2*V1[j]/V1Norm/V2Norm;
    V1til[j] = V1[j] - V1V2*V2[j]/V1Norm/V2Norm;
    V2tilNorm += V2til[j]*V2til[j];
    V1tilNorm += V1til[j]*V1til[j];
  }
  V2tilNorm = sqrt(V2tilNorm);
  V1tilNorm = sqrt(V1tilNorm);
  for(j=0;j<3;j++) {
    axes(0,j) = V1[j];
    axes(1,j) = V2til[j]/V2tilNorm;
  }
  switch(side) {
  case 0:
  case 2:
    for(i=0;i<3;i++)
      normal[i] = V2til[i]/V2tilNorm;
    break;
  case 1:
  case 3:
    for(i=0;i<3;i++)
      normal[i] = V1til[i]/V1tilNorm;
    break;
  }
  if(side == 0 || side == 3) for(i=0;i<3;i++) normal[i] *= -1.;

  axes(2,0) = axes(0,1)*axes(1,2)-axes(0,2)*axes(1,1);
  axes(2,1) = -axes(0,0)*axes(1,2)+axes(0,2)*axes(1,0);
  axes(2,2) = axes(0,0)*axes(1,1)-axes(0,1)*axes(1,0);
  return;
  */
}


void TPZGeoElQ2dCyl::Print(ostream & out) {
    TPZGeoElQ2d::Print(out);
  int i;
  //=====Cesar=====
  out << "Nodes coordinates  ";
  for (i = 0;i < NNodes();i++)
    out << "\n Coordinate: \t" << i << "\t Radius: \t" << fX[i][0] << "\t Theta: \t" << fX[i][1] << "\t Z: \t" << fX[i][2];
}
