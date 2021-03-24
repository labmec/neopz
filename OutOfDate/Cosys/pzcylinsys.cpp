/**
 * \file
 * @brief DEPRECATED FILE. Contains implementation of the methods to Cylindrical Coordinate System.
 */
#include "pzcylinsys.h"
#include "pzcartsys.h"

/* Default empty constructor
 */
TPZCylinsys::TPZCylinsys() : TPZCosys() {
}

/* Create one object from other changing the origin
 */
TPZCylinsys::TPZCylinsys(int num, TPZCartsys* ref) : TPZCosys(num,ref){
}

/* Return in the reference coordinate system one point given in current
 coordinate system */
void TPZCylinsys::ToReference(TPZVec<REAL> &point) {
    REAL theta = point[1];
    point[1] = point[0] * sin(theta);
    point[0] *= cos(theta);
}

/* Return the current coordinate system of one point given in reference
 coordinate system */
void TPZCylinsys::FromReference(TPZVec<REAL> &point) {
    REAL theta;
    if (point[1] == 0.0)
        theta = 0.0;
    else if (point[0] == 0.0)
        theta = (point[1] < 0.) ? asin(-1.0) : asin(1.0);
	else {
		REAL ratio = point[1]/point[0];
		theta = atan(ratio);
	}
	
    if (theta >= 0.0 && point[0] < 0.0) theta += 2.*asin(1.0);
    if (theta <= 0.0 && point[1] > 0.0) theta += 2.*asin(1.0);
	
    point[0] = sqrt(point[0]*point[0] + point[1]*point[1]);
    point[1] = theta;
    //  point[1] = theta * 90./asin(1.0);
}

/* Verify the element theta vector consitancy
 */
void TPZCylinsys::VerifyRange(TPZFMatrix &points){
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
}

/* Calculates the Changing coordinate system transformation gradient
 */
void TPZCylinsys::TransformGradient(TPZVec<REAL> &X, TPZFMatrix &GradX, TPZVec<REAL> &x, TPZFMatrix &Gradx, TPZCosys *dest){
	
	int col,ncol=GradX.Cols();
	for(col=0; col<ncol; col++) {
		x[col]=X[col];
		Gradx(0,col) = (GradX(0,col)*cos(X[1]))+(X[0]*sin(X[1]*GradX(1,col)));
		Gradx(1,col) = (GradX(0,col)*sin(X[1]))+(X[0]*cos(X[1]*GradX(1,col)));
		Gradx(2,col) = GradX(2,col);
	}
	
	if (dest != fReference) {
		TPZFMatrix gradin(GradX);
		FromReference(x);
		fReference->TransformGradient(x,gradin,x,Gradx,dest);
	}
	
}
