// Emacs will be in -*- Mode: c++ -*-
//
// ********** DO NOT REMOVE THIS BANNER **********
//
//  Nicolas Di Cesare <Nicolas.Dicesare@ann.jussieu.fr>
//  http://www.ann.jussieu.fr/~dicesare
//
//  Benchmark between ADOL-C and Fad
//  in forward mode
// 
//********************************************************
//
//  KCC 3.3f (kai C++) doesn't compile ADOL-C 
//  egcs 1.1.2 (cygnus) is OK.
//
//********************************************************

// C++ includes
#include <iomanip>
              
// C include
#include <cmath>

// ADOL-C includes
#include <adouble.h>                   // use of active double and taping
#include <DRIVERS/drivers.h>           // use of "Easy To Use" drivers 
                                       // gradient(.) and hessian(.)
// Fad include
#include <Fad/fad.h>

// vector and timer includes
#include <utils/timer.h>
#include <utils/vectors.h>



double testADOLC(Timer& timer, const int n, const int nloop, double & yp, double& errg);
double testFAD(Timer& timer, const int n, const int nloop, double & yp, double& errg);

/****************************************************************************/
/*                                                             MAIN PROGRAM */
int main() 
{ 

  //int n,k;
  const int nloop = 10000;
  const int nbench = 3;

  double errg1 = 0., errg2 = 0.;
  double  yp1 = 0.0, yp2 = 0.;
  double time1 = 0., time2 = 0.;


  Timer timer(123123000.12);

  for (int n=1; n<100; ++n) {

    //clear cache
    testFAD(timer, n, nloop, yp1, errg1);
    testADOLC(timer, n, nloop, yp2, errg2);

    time1 = time2 = 0.;
    for (int j=0; j<nbench; ++j)
      {
	errg1 = errg2 = 0.;
	yp1 = yp2 = 0.;
	time1 += testFAD(timer, n, nloop, yp1, errg1);
	time2 += testADOLC(timer, n, nloop, yp2, errg2);

      }

    cout.setf(ios::fixed,ios::floatfield);
    cout << setw(4) << n << setw(12) << time1/(double)nbench << setw(12) << time2/(double)nbench << endl;
  
  }


  return 1;
} // end main


double testADOLC(Timer& timer, const int n, const int nloop, double & yp, double & errg)
{
  double mtime;
  int i,k;
  double *xp = new double[n];          
  adouble *x = new adouble[n];         // or: adoublev x(n);
  adouble  y = 1., tmp;
  double* g = new double[n];           
  
  for(i=0; i<n; i++)
    xp[i] = (i+1.0)/(2.0+i);           // some initialization

  timer.start();
  // repeat nloop time the loop to obtain
  // a mesurable time
  for (k=0; k<nloop; ++k) {
    
    trace_on(1);                         // tag = 1, keep = 0 by default

    y = 1.;
    for(i=0; i<n; i++){
      tmp <<= xp[i];
      y = y+tmp+tmp+tmp+tmp;// y = y*x[i]*x[i]+x[i]+x[i];
    }

    y >>= yp;
    trace_off();
  
    gradient(1,n,xp,g);                  // gradient evaluation
  }

  timer.stop();
  mtime = timer.elapsedSeconds();

  for(i=0; i<n; i++) 
    errg += g[i];       // vanishes analytically.

  delete[] x;                        // not needed if x adoublev
  delete [] g;
  delete [] xp;
  
  return mtime;
}


double testFAD(Timer& timer, const int n, const int nloop, double & yp, double& errg)
{
  int i,k;
  double mtime;
  double *xp = new double[n];          
  Fad<double> *x = new Fad<double>[n];         // or: adoublev x(n);
  Fad<double>  y(n,1.), tmp;
  
  for(i=0; i<n; i++)
    xp[i] = (i+1.0)/(2.0+i);           // some initialization
    
  timer.start();
  // repeat nloop time the loop to obtain
  // a mesurable time
  for (k=0; k<nloop; ++k) {
    
    y = 1.;
    for(i=0; i<n; i++)
      { 
	tmp = Fad<double>(n, i, xp[i]);
	y = y+tmp+tmp+tmp+tmp;// y = y*x[i]*x[i]+x[i]+x[i];
      } // end for

    yp = y.val();
  }
  timer.stop();
  mtime = timer.elapsedSeconds();

  
  for(i=0; i<n; i++) 
    errg += y.d(i);       // vanishes analytically.

  delete [] x;                        // not needed if x adoublev
  delete [] xp;

  return mtime;
}
