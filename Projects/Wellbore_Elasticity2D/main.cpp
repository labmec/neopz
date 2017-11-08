#include "GeometricMesh.hpp"
#include "ComputationalMesh.hpp"
#include "Problem2D.hpp"
#include "Problem3D.hpp"
#include "AproximationRates.hpp"
//#include <clapack.h>


int main(int argc, char *argv[])
{
    int ncircle = 30;
    int nradial = 25;
    
    // Define se a sol. analitica sera utilizada ou nao
    // analytic=0 nao usa sol analitica como prestress e BC
    // analytic=1 usa sol analitica como prestress e BC (zerar BCond0 e BCond1)
    // analytic=2 nao usa sol analitica como prestress mas usa como BC (zerar BCond0 e BCond1)
    int analytic = 0;
    
    // define se havera projecao no plano horizontal
    int projection = 0;
    
    // define disposicao do poco - inclined == 1
    int inclinedwellbore = 0;
    
    // pressao da lama de perfuracao - MPa
    REAL Pwb = -19.5; //-10.5
    
    REAL rw = 0.10795;
    REAL rext = 3.0;
    REAL drdcirc = 0.5;
    
    // Define Posicao do Poco
    REAL direction   = 0.; // Azimuth em graus
    REAL inclination = 0.; // Polar Inclination em graus / wellbore inclination
    
    // Tensoes in Situ, horizontais e vertical em MPa
    REAL SigmaV = -48.2; // tensao vertical
    REAL Sigmah = -45.9; // tensao horizontal menor
    REAL SigmaH = -62.1; // tensao horizontal maior
    
    bool isStochastic = true;
    
    std::ofstream solutionfile("f1_solution.csv");
    solutionfile << "Case,Total plastified area" << std::endl;
    
    int ncases = 1;
    
    for(int i=0; i < ncases; i++){
        Problem2D(rw, rext, ncircle, nradial, projection, inclinedwellbore, analytic, SigmaV, Sigmah,
                  SigmaH, Pwb, drdcirc, direction, inclination, isStochastic, solutionfile, i);
    }
    solutionfile.close();
    //Problem3D();
    
    //ApproximationRates();
    
    return 0;
}
