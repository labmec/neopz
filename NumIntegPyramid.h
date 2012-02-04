/*
 *  NumIntegPyramid.h
 *  PZ
 *
 *  Created by Jorge on 2/3/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef NUMINTEGPYRAMIDHH
#define NUMINTEGPYRAMIDHH

# include <cstdlib>
# include <iomanip>
# include <iostream>
# include <cmath>
# include <ctime>
# include <string>
# include <fstream>

using namespace std;

void jacobi_compute ( int order, double alpha, double beta, double x[], 
					 double w[] );
void jacobi_recur ( double *p2, double *dp2, double *p1, double x, int order, 
				   double alpha, double beta, double b[], double c[] );
void jacobi_root ( double *x, int order, double alpha, double beta, 
				  double *dp2, double *p1, double b[], double c[] );
void legendre_compute ( int order, double x[], double w[] );
void pyramid_handle ( int legendre_order, int jacobi_order, string filename );
double r8_abs ( double x );
double r8_epsilon ( );
double r8_gamma ( double x );
void r8mat_write ( string output_filename, int m, int n, double table[] );
void timestamp ( );

#endif