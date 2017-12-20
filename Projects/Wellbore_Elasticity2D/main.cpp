#include "GeometricMesh.hpp"
#include "ComputationalMesh.hpp"
#include "Problem2D.hpp"
#include "Problem3D.hpp"
#include "AproximationRates.hpp"
//#include <clapack.h>
#include "pzrandomfield.h"

// Read Decomposed Matrix from File
template<typename TVar>
TPZFMatrix<TVar> readDecomposedMatrixFromFile(int nSquareElements, int matsize, int stochasticInclined) {

	TPZFMatrix<TVar> M(nSquareElements, nSquareElements, 0.);

	if (stochasticInclined == 1) {
		M.Resize(matsize, matsize);

		// Setar valores de M obtidos do Mathematica (Decomposed Matrix)
		std::ifstream DecMatFile("../decomposed_matrix/decomposed_matrix.tbl");

		if (!DecMatFile.good()) {
			std::cout << "Decomposed Matrix (.tbl) file does not exist!\n" << std::endl;
			DebugStop();
		}

		std::string line;

		int i = 0;

		while (std::getline(DecMatFile, line)) {
			REAL value;
			int j = 0;
			std::stringstream ss(line);

			while (ss >> value) {
				M(i, j) = value;
				j++;
			}
			i++;
		}

		M.Resize(nSquareElements, matsize);

		return M;

	}

	else {
		// Setar valores de M obtidos do Mathematica (Decomposed Matrix)
		std::ifstream DecMatFile("../decomposed_matrix/decomposed_matrix.tbl");

		if (!DecMatFile.good()) {
			std::cout << "Decomposed Matrix (.tbl) file does not exist!\n" << std::endl;
			DebugStop();
		}

		std::string line;

		int i = 0;

		while (std::getline(DecMatFile, line)) {
			REAL value;
			int j = 0;
			std::stringstream ss(line);

			while (ss >> value) {
				M(i, j) = value;
				j++;
			}
			i++;
		}

		return M;

	}
}

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
    int inclinedwellbore = 1;
    
    // pressao da lama de perfuracao - MPa
    REAL Pwb = -19.5; //-10.5
    
    REAL rw = 0.10795;
    REAL rext = 2.0;
    REAL drdcirc = 0.5;
    
    // Define Posicao do Poco
    REAL direction   = 30.; // Azimuth em graus
    REAL inclination = 45.; // Polar Inclination em graus / wellbore inclination
    
    // Tensoes in Situ, horizontais e vertical em MPa
    REAL SigmaV = -48.2; // tensao vertical
    REAL Sigmah = -45.9; // tensao horizontal menor
    REAL SigmaH = -62.1; // tensao horizontal maior
    
    bool isStochastic = true;
    
    std::ofstream solutionfile("f1_solution_inclined.csv");
    solutionfile << "Case,Total plastified area" << std::endl;
    
    int ncases = 10000;
	
	int nLayers = 8;
	REAL fH = 2 * rext; // altura total do cilindro em metros
	REAL fh = fH / nLayers; // altura de cada cubo (elemento) em metros 
	int nSquareElements = nradial * ncircle;
	int matsize = nSquareElements * (fH / fh) + nSquareElements;

	std::cout << "Read decomposed Matrix" << std::endl;
	TPZFMatrix<STATE> M = readDecomposedMatrixFromFile<STATE>(nSquareElements, matsize, inclinedwellbore);
    
    for(int i=0; i < ncases; i++){
        Problem2D(rw, rext, ncircle, nradial, projection, inclinedwellbore, analytic, SigmaV, Sigmah,
                  SigmaH, Pwb, drdcirc, direction, inclination, isStochastic, solutionfile, i, M);
    }
    solutionfile.close();
    //Problem3D();
    
    //ApproximationRates();
    
    return 0;
}
