/** standard C I/O functions */
#include <cstdio>

/** pz full matrix class */
#include "pzfmatrix.h"

#include <stdlib.h>

#include <sys/time.h>
/** util function - return wall clock time in seconds */
double mysecond()
{
    struct timeval tp;
    struct timezone tzp;
    gettimeofday(&tp,&tzp);
    return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
} /* mysecond */

int main()
{    
    setenv("VECLIB_MAXIMUM_THREADS", "1", true);

    int n = 1500;
    TPZFMatrix<double> A(n, n, 0.0);
    TPZFMatrix<double> B(n, n, 0.0);
    
    /** putting values diferent from zero in the matrices */
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            A(i, j) = ( i * n + j ) + 1;
            if (i == j) B(i, j) = 1.0; //( i * n + j ) * 0.2;
        }
    }
    double begin = mysecond();
    TPZFMatrix<double> C = A * B;
    double end = mysecond();
    
    
    //C.Print(std::cout);
    
    printf("Multiply Time: %.3lf\n", end-begin);
    
} /* main */