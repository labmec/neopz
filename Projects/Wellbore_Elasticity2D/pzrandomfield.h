#ifndef PZRANDOMFIELD_H
#define PZRANDOMFIELD_H

#include "pzfunction.h"
#include <iostream>
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzerror.h"
#include "tpzverysparsematrix.h"
#include "pzsfulmat.h"

#include <math.h>
#include <complex>
#include <string>


////#include <nrutil.h>

//////******** NRUTIL.H **************//
////
//////#include <nrutil.h>
////
//////#ifndef _NR_UTILS_H_
//////#define _NR_UTILS_H_
////
//static float sqrarg;
//#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
//
//static double dsqrarg;
//#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)
//
//static double dmaxarg1,dmaxarg2;
//#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
//(dmaxarg1) : (dmaxarg2))
//
//static double dminarg1,dminarg2;
//#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?\
//(dminarg1) : (dminarg2))
//
//static float maxarg1,maxarg2;
//#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
//(maxarg1) : (maxarg2))
//
//static float minarg1,minarg2;
//#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
//(minarg1) : (minarg2))
//
//static long lmaxarg1,lmaxarg2;
//#define LMAX(a,b) (lmaxarg1=(a),lmaxarg2=(b),(lmaxarg1) > (lmaxarg2) ?\
//(lmaxarg1) : (lmaxarg2))
//
//static long lminarg1,lminarg2;
//#define LMIN(a,b) (lminarg1=(a),lminarg2=(b),(lminarg1) < (lminarg2) ?\
//(lminarg1) : (lminarg2))
//
//static int imaxarg1,imaxarg2;
//#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
//(imaxarg1) : (imaxarg2))
//
//static int iminarg1,iminarg2;
//#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
//(iminarg1) : (iminarg2))
//
//#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
//
////#if defined(__STDC__) || defined(ANSI) || defined(NRANSI) /* ANSI */
//
//void nrerror(char error_text[]);
//float *vector(long nl, long nh);
//int *ivector(long nl, long nh);
//unsigned char *cvector(long nl, long nh);
//unsigned long *lvector(long nl, long nh);
//double *dvector(long nl, long nh);
//float **matrix(long nrl, long nrh, long ncl, long nch);
//double **dmatrix(long nrl, long nrh, long ncl, long nch);
//int **imatrix(long nrl, long nrh, long ncl, long nch);
//float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch,
//                  long newrl, long newcl);
//float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch);
//float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
//void free_vector(float *v, long nl, long nh);
//void free_ivector(int *v, long nl, long nh);
//void free_cvector(unsigned char *v, long nl, long nh);
//void free_lvector(unsigned long *v, long nl, long nh);
//void free_dvector(double *v, long nl, long nh);
//void free_matrix(float **m, long nrl, long nrh, long ncl, long nch);
//void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
//void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
//void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch);
//void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch);
//void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,
//                   long ndl, long ndh);
////
//#else /* ANSI */
///* traditional - K&R */
////
//void nrerror();
//float *vector();
//float **matrix();
//float **submatrix();
//float **convert_matrix();
//float ***f3tensor();
//double *dvector();
//double **dmatrix();
//int *ivector();
//int **imatrix();
//unsigned char *cvector();
//unsigned long *lvector();
//void free_vector();
//void free_dvector();
//void free_ivector();
//void free_cvector();
//void free_lvector();
//void free_matrix();
//void free_submatrix();
//void free_convert_matrix();
//void free_dmatrix();
//void free_imatrix();
//void free_f3tensor();
////
//////#endif /* ANSI */
//////
//////#endif /* _NR_UTILS_H_ */
////
//////****************************** FIm NRUTIL.H ***************************//


//#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
//
//void free_vector(float *v, long nl, long nh);
//void free_vector();
//
//static float maxarg1,maxarg2;
//#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
//(maxarg1) : (maxarg2))
//
//static int imaxarg1,imaxarg2;
//#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
//(imaxarg1) : (imaxarg2))
//
//static int iminarg1,iminarg2;
//#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
//(iminarg1) : (iminarg2))
//
//float *vector(long nl, long nh);
//float *vector();
//
//void nrerror(char error_text[]);
//void nrerror();


template<class TVar>
class TPZRandomField : public TPZFunction<TVar>
{
	
    void (*fFunc)(const TPZVec<REAL> &x, TPZVec<TVar> &f, TPZFMatrix<REAL> &K);
    void (*fFunc2)(const TPZVec<REAL> &x, TPZVec<TVar> &f, TPZFMatrix<TVar> &gradf);
    void (*fFunc3)(const TPZVec<REAL> &x, REAL ftime, TPZVec<TVar> &f, TPZFMatrix<TVar> &gradf);
    
    int fPorder;
public:
	
	/** @brief Class constructor */
	TPZRandomField() : TPZFunction<TVar>(), fPorder(-1)
    {
        fFunc = 0;
		fFunc2 = 0;
		fFunc3 = 0;		
    }
	
	/** @brief Class destructor */
	virtual ~TPZRandomField()
    {
        
    }
    
    TPZRandomField(void (*FuncPtr)(const TPZVec<REAL> &x, TPZVec<TVar> &val, TPZFMatrix<REAL> &K))
    {
        fFunc = FuncPtr;
		fFunc2 = 0;
		fFunc3 = 0;
		fPorder = -1;
    }
    
    TPZRandomField(void (*FuncPtr)(const TPZVec<REAL> &x, TPZVec<TVar> &val, TPZFMatrix<TVar> &gradf))
    {
		fFunc = 0;
        fFunc2 = FuncPtr;		
		fFunc3 = 0;
		fPorder = -1;
    }
	
    TPZRandomField(void (*FuncPtr)(const TPZVec<REAL> &x, REAL ftime, TPZVec<TVar> &val, TPZFMatrix<TVar> &gradf))
    {
		fFunc = 0;
        fFunc2 = 0;		
		fFunc3 = FuncPtr;
		fPorder = -1;
    }	
    
    TPZRandomField(const TPZRandomField &cp) : fFunc(cp.fFunc), fFunc2(cp.fFunc2), fFunc3(cp.fFunc3), fPorder(cp.fPorder)
    {
        
    }
    

    TPZRandomField &operator=(const TPZRandomField &cp)
    {
        fFunc = cp.fFunc;
		fFunc2 = cp.fFunc2;
		fFunc3 = cp.fFunc3;
        fPorder = cp.fPorder;
        return *this;
    }
	/**
	 * @brief Performs function computation
	 * @param x point coordinate which is suppose to be in real coordinate system but can be in master coordinate system in derived classes.
	 * @param f function values
	 * @param df function derivatives
	 */
	virtual void Execute(const TPZVec<REAL> &x, TPZVec<TVar> &f, TPZFMatrix<TVar> &df)
    {
        if (!fFunc2) {
			DebugStop();
		}
        fFunc2(x, f, df);
    }
	
	/**
	 * @brief Performs time dependent function computation
	 * @param x point coordinate which is suppose to be in real coordinate system but can be in master coordinate system in derived classes.
	 * @param ftime  time to evaluate	 
	 * @param f function values
	 * @param gradf function derivatives
	 */	
	virtual void Execute(const TPZVec<REAL> &x, REAL ftime, TPZVec<TVar> &f, TPZFMatrix<TVar> &gradf)
    {
        if (!fFunc3) {
			DebugStop();
		}
        fFunc3(x, ftime, f, gradf);
    }	
    
	/**
	 * @brief Execute method receiving axes. It is used in shape functions
	 * @note NOT IMPLEMENTED
	 */
	virtual void Execute(const TPZVec<REAL> &x, const TPZFMatrix<REAL> &axes, TPZVec<TVar> &f, TPZFMatrix<TVar> &df){
        DebugStop();
    }
    
    /**
	 * @brief Simpler version of Execute method which does not compute function derivatives 
	 * @param x point coordinate which is suppose to be in real coordinate system but can be in master coordinate system in derived classes.
	 * @param f function values
	 */
    virtual void Execute(const TPZVec<REAL> &x, TPZVec<TVar> &f, TPZFMatrix<REAL> &K){
		//if (!fFunc) {
		//	DebugStop();
		//}
        //fFunc(x,f);
        
		//AQUI NATHALIA
        //REAL E = rand() % 3000 + 15300+1; // not uniform
        
        REAL E = arc4random_uniform(3001) + 15300; //uniform
        
       /* The arc4random() function uses the key stream generator employed by the arc4 cipher, which uses 8*8 8 bit S-Boxes.  The S-Boxes can be in about (2**1700) states.  The arc4random() function returns pseudo-random numbers in the range of 0 to (2**32)-1, and therefore has twice the range of rand(3) and random(3). */
        /* arc4random_uniform() will return a uniformly distributed random number less than upper_bound.
         arc4random_uniform() is recommended over constructions like ``arc4random() % upper_bound'' as it avoids "modulo bias" when the upper bound is not a power of two. */
        
		f[0] = E;
        
        
    }
    
    
    
        //******** Cholesky ******//
    
//        //TPZSMatrix (const TPZSMatrix< TVar > &test);
//        TPZSFMatrix<STATE> test;
//
//        test.Decompose_Cholesky();
//        

        
        
//        //*******  SVD ********//
//        
//     virtual void svdcmp(float **a, int m, int n, float w[], float **v)
////        Given a matrix a[1..m][1..n], this routine computes its singular value decomposition, A =
////        U ·W ·V T . The matrix U replaces a on output. The diagonal matrix of singular values W is output
////        as a vector w[1..n]. The matrix V (not the transpose V T ) is output as v[1..n][1..n].
//        {
//            float pythag(float a, float b);
//            int flag,i,its,j,jj,k,l,nm;
//            float anorm,c,f,g,h,s,scale,x,y,z,*rv1;
//            rv1=vector(1,n);
//            g=scale=anorm=0.0;                      //Householder reduction to bidiagonal form.
//            for (i=1;i<=n;i++) {
//                l=i+1;
//                rv1[i]=scale*g;
//                g=s=scale=0.0;
//                if (i <= m) {
//                    for (k=i;k<=m;k++) scale += fabs(a[k][i]);
//                    if (scale) {
//                        for (k=i;k<=m;k++) {
//                            a[k][i] /= scale;
//                            s += a[k][i]*a[k][i];
//                        }
//                        f=a[i][i];
//                        g = -SIGN(sqrt(s),f);
//                        h=f*g-s;
//                        a[i][i]=f-g;
//                        for (j=l;j<=n;j++) {
//                            for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
//                            f=s/h;
//                            for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
//                        }
//                        for (k=i;k<=m;k++) a[k][i] *= scale;
//                    }
//                }
//                w[i]=scale *g;
//                g=s=scale=0.0;
//                if (i <= m && i != n) {
//                    for (k=l;k<=n;k++) scale += fabs(a[i][k]);
//                    if (scale) {
////        
//                        for (k=l;k<=n;k++) {
//                            a[i][k] /= scale;
//                            s += a[i][k]*a[i][k];
//                        }
//                        f=a[i][l];
//                        g = -SIGN(sqrt(s),f);
//                        h=f*g-s;
//                        a[i][l]=f-g;
//                        for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
//                        for (j=l;j<=m;j++) {
//                            for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
//                            for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
//                        }
//                        for (k=l;k<=n;k++) a[i][k] *= scale;
//                    }
//                }
//                anorm=FMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
//            }
//            for (i=n;i>=1;i--) {                    //Accumulation of right-hand transformations.
//                if (i < n) {
//                    if (g) {
//                        for (j=l;j<=n;j++)          //Double division to avoid possible underflow.
//                            v[j][i]=(a[i][j]/a[i][l])/g;
//                        for (j=l;j<=n;j++) {
//                            for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
//                            for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
//                        }
//                    }
//                    for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
//                }
//                v[i][i]=1.0;
//                g=rv1[i];
//                l=i;
//            }
//            for (i=IMIN(m,n);i>=1;i--) {            //Accumulation of left-hand transformations.
//                l=i+1;
//                g=w[i];
//                for (j=l;j<=n;j++) a[i][j]=0.0;
//                if (g) {
//                    g=1.0/g;
//                    for (j=l;j<=n;j++) {
//                        for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
//                        f=(s/a[i][i])*g;
//                        for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
//                    }
//                    for (j=i;j<=m;j++) a[j][i] *= g;
//                } else for (j=i;j<=m;j++) a[j][i]=0.0;
//                ++a[i][i];
//            }
//            for (k=n;k>=1;k--) {                    //Diagonalization of the bidiagonal form: Loop over
//                for (its=1;its<=30;its++) {         //singular values, and over allowed iterations.
//                    flag=1;
//                    for (l=k;l>=1;l--) {            //Test for splitting.
//                        nm=l-1;                     //Note that rv1[1] is always zero.
//                        if ((float)(fabs(rv1[l])+anorm) == anorm) {
//                            flag=0;
//                            break;
//                        }
//                        if ((float)(fabs(w[nm])+anorm) == anorm) break;
//                    }
//                    if (flag) {
//                        c=0.0;                      //Cancellation of rv1[l], if l > 1.
//                            s=1.0;
//                        for (i=l;i<=k;i++) {
//                            
//                            f=s*rv1[i];
//                            rv1[i]=c*rv1[i];
//                            if ((float)(fabs(f)+anorm) == anorm) break;
//                            g=w[i];
//                            h=pythag(f,g);
//                            w[i]=h;
//                            h=1.0/h;
//                            c=g*h;
//                            s = -f*h;
//                            for (j=1;j<=m;j++) {
//                                y=a[j][nm];
//                                z=a[j][i];
//                                a[j][nm]=y*c+z*s;
//                                a[j][i]=z*c-y*s;
//                            }
//                        }
//                    }
//                    z=w[k];
//                    if (l == k) {                       //Convergence.
//                        if (z < 0.0) {                  //Singular value is made nonnegative.
//                            w[k] = -z;
//                            for (j=1;j<=n;j++) v[j][k] = -v[j][k];
//                        }
//                        break;
//                    }
//                    if (its == 30) nrerror("no convergence in 30 svdcmp iterations");
//                    x=w[l];                             //Shift from bottom 2-by-2 minor.
//                    nm=k-1;
//                    y=w[nm];
//                    g=rv1[nm];
//                    h=rv1[k];
//                    f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
//                    g=pythag(f,1.0);
//                    f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
//                    c=s=1.0;                            //Next QR transformation:
//                    for (j=l;j<=nm;j++) {
//                        i=j+1;
//                        g=rv1[i];
//                        y=w[i];
//                        h=s*g;
//                        g=c*g;
//                        z=pythag(f,h);
//                        rv1[j]=z;
//                        c=f/z;
//                        s=h/z;
//                        f=x*c+g*s;
//                        g = g*c-x*s;
//                        h=y*s;
//                        y *= c;
//                        for (jj=1;jj<=n;jj++) {
//                            x=v[jj][j];
//                            z=v[jj][i];
//                            v[jj][j]=x*c+z*s;
//                            v[jj][i]=z*c-x*s;
//                        }
//                        z=pythag(f,h);
//                        w[j]=z;                         //Rotation can be arbitrary if z = 0.
//                            if (z) {
//                                z=1.0/z;
//                                c=f*z;
//                                s=h*z;
//                            }
//                        f=c*g+s*y;
//                        x=c*y-s*g;
//                        
//                        for (jj=1;jj<=m;jj++) {
//                            y=a[jj][j];
//                            z=a[jj][i];
//                            a[jj][j]=y*c+z*s;
//                            a[jj][i]=z*c-y*s;
//                        }
//                    }
//                    rv1[l]=0.0;
//                    rv1[k]=f;
//                    w[k]=x;
//                }
//            }
//            free_vector(rv1,1,n);
//        }
    



	/** @brief Returns number of functions. */
	virtual int NFunctions()
    {
        return 1;
    }
    
    void SetPolynomialOrder(int porder)
    {
        fPorder = porder;
    }
	
	/** @brief Polynomial order of this function. */
	/** In case of non-polynomial function it can be a reasonable approximation order. */
	virtual int PolynomialOrder() 
    {
#ifdef DEBUG
        if (fPorder == -1) {
            DebugStop();
        }
#endif
        return fPorder;
    }
	
	/** @brief Unique identifier for serialization purposes */
	virtual int ClassId() const
    {
        return -1;
    }
	
	/** @brief Saves the element data to a stream */
	virtual void Write(TPZStream &buf, int withclassid)
    {
//        DebugStop();
        TPZSaveable::Write(buf,withclassid);
    }
	
	/** @brief Reads the element data from a stream */
	virtual void Read(TPZStream &buf, void *context)
    {
        DebugStop();
    }
};

#endif
