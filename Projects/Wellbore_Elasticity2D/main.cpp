#include "GeometricMesh.hpp"
#include "ComputationalMesh.hpp"
#include "Problem2D.hpp"
#include "Problem3D.hpp"
#include "AproximationRates.hpp"


int main(int argc, char *argv[])
{
    int ncircle = 30;
    int nradial = 25;
    
    REAL rw = 0.1;
    REAL rext = 3.0;
    REAL drdcirc = 0.5;
    
    // Define Posicao do Poco
    REAL direction   = 0.; // Azimuth em graus
    REAL inclination = 0.; // Polar Inclination em graus
    
    bool isStochastic = true;
    
    Problem2D(rw, rext, ncircle, nradial, drdcirc, direction, inclination,
              isStochastic);
    
    //Problem3D();
    
    //ApproximationRates();
    
    return 0;
}
