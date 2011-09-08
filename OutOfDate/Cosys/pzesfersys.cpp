/**
 * \file
 * @brief DEPRECATED FILE. Contains implementation of the methods to Spherical Coordinate System.
 */
//METHODS DEFINITION FOR CLASS COSYS

#include "pzesfersys.h"
#include "pzcartsys.h"

//***************************************
//***************************************
TPZEsfersys::TPZEsfersys() : TPZCosys() {
	
}
//***************************************
//***************************************
TPZEsfersys::TPZEsfersys(int num, TPZCartsys* ref) : TPZCosys(num,ref) {
	
}

//***************************************
//***************************************
void TPZEsfersys::ToReference(TPZVec<REAL> &point) {
	
    REAL phi   = point[2];//*asin(1.)/90.;
    REAL theta = point[1];//*asin(1.)/90.;
	
    point[2] = point[0] * cos(phi);
    point[0] *= sin(phi);
    point[1] = point[0] * sin(theta);
    point[0] *= cos(theta);
}
//***************************************
//***************************************
void TPZEsfersys::FromReference(TPZVec<REAL> &point) {
	
    REAL theta;
    if (point[0] == 0.0)
        theta = (point[1] < 0.) ? asin(-1.0) : asin(1.0);
    else
        theta = atan(point[1]/point[0]);
	
    if (theta >= 0.0 && point[0] < 0.0) theta += 2.*asin(1.0);
    if (theta <= 0.0 && point[1] > 0.0) theta += 2.*asin(1.0);
	
    point[0] = sqrt(point[0]*point[0] + point[1]*point[1]);
    point[1] = theta;// * 90./asin(1.0);
	
    REAL phi;
    if (point[2] == 0.0)
        phi = asin(1.0);
    else
        phi = atan(point[0]/point[2]);
	
    if (phi <= 0.0) phi += 2.*asin(1.0);
	
    point[0] = sqrt(point[0]*point[0] + point[2]*point[2]);
    point[2] = phi;// * 90./asin(1.0);
	
	
}
//***************************************
//***************************************
void TPZEsfersys::TransformGradient(TPZVec<REAL> &X, TPZFMatrix &GradX, TPZVec<REAL> &x, TPZFMatrix &Gradx, TPZCosys *dest){
    REAL ct,st,cf,sf;
	
    ct = cos(X[1]);
    st = sin(X[1]);
    cf = cos(X[2]);
    sf = sin(X[2]);
	
	int col,ncol=GradX.Cols();
	
	x=X;
	for(col=0; col<ncol; col++) {
		Gradx(0,col)=(GradX(0,col)*sf*ct)+(X[0]*ct*cf*GradX(1,col))-(X[0]*st*sf*GradX(2,col));
		Gradx(1,col)=(GradX(0,col)*sf*st)+(X[0]*ct*sf*GradX(1,col))+(X[0]*st*cf*GradX(2,col));
		Gradx(2,col)=(GradX(0,col)*cf)-(X[0]*sf*GradX(2,col));
	}
	
    if (dest != fReference) {
        TPZFMatrix gradin(GradX);
		FromReference(x);
		TransformGradient(x,gradin,x,Gradx,dest);
	}
}
//***************************************
//***************************************
void TPZEsfersys::VerifyRange(TPZFMatrix &points){
	int i,j,k=0;
    for(i=0;i<points.Rows();i++){
        for(j=i+1;j<points.Rows();j++){
            if ((points(i,1)-points(j,1))>M_PI || (points(i,1)-points(j,1))< -(M_PI) ){
                k=1;
	            break;
            }
        }
        if(k) break;
    }
    if(k) {
        for (i=0;i<points.Rows();i++){
            if (points(i,1)>0) points(i,1)=points(i,1)-2.*M_PI;
        }
    }
	//verify phi
	/*	k=0;
	 for(i=0;i<points.Rows();i++){
	 for(j=i+1;j<points.Rows();j++){
	 if ((points(i,2)-points(j,1))>PI || (points(i,2)-points(j,2))< -(PI) ){
	 k=1;
	 break;
	 }
	 }
	 if(k) break;
	 }
	 if(k) {
	 for (i=0;i<points.Rows();i++){
	 if (points(i,2)>PI) points(i,2)=points(i,2)-2.*PI;
	 }
	 }*/
}
