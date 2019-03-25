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
		std::ifstream DecMatFile("../decomposed_matrix/decomposed_matrix_inclined.tbl");

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
        std::ifstream DecMatFile("../decomposed_matrix/decomposed_matrix_vertical.tbl");
        
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
    int nradial = 25; //25; //NANANAN teste malha progressao
    
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
    REAL Pwb = -19.5;//-19.5; //-10.52
    
    REAL rw = 0.10795;
    REAL rext = 2.0; //Inclinado: 2m
    REAL drdcirc = 0.5;
    
    // Define Posicao do Poco
    REAL direction   = 0.;//30.; // Azimuth em graus
    REAL inclination = 0.;//45.; // Polar Inclination em graus / wellbore inclination
    
    // Tensoes in Situ, horizontais e vertical em MPa
    REAL SigmaV = -48.2; //-48.053;  // tensao vertical
    REAL Sigmah = -45.9; //-48.0107; // tensao horizontal menor
    REAL SigmaH = -62.1; //-68.3251; // tensao horizontal maior
    
    bool isStochastic = true; // Stochastic?
    bool isInSituStoch = false; // In-situ stresses stochastic?
    
    REAL cv = 0.1; // variation coefficient of In-situ stresses
    
    std::ofstream solutionfile("Vertical_Stoch_Pw19_5_2m.csv");
    solutionfile << "Case,Total plastified area" << std::endl;
    
    int ncases = 1;
	
    int nLayers = 8;
	REAL fH = 2 * rext; // altura total do cilindro em metros
	REAL fh = fH / nLayers; // altura de cada cubo (elemento) em metros 
	int nSquareElements = nradial * ncircle;
	int matsize = nSquareElements * (fH / fh) + nSquareElements;

    std::cout << "Read decomposed Matrix" << std::endl;
    TPZFMatrix<STATE> M = readDecomposedMatrixFromFile<STATE>(nSquareElements, matsize, inclinedwellbore);
    
    if (isStochastic == true  && isInSituStoch==true) {
        // Random In-Situ Stresses - Normal Distribution
        TPZFMatrix<REAL> Sig_V (ncases, 1, 0.);
        TPZFMatrix<REAL> Sig_h (ncases, 1, 0.);
        TPZFMatrix<REAL> Sig_H (ncases, 1, 0.);
        
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::default_random_engine generator (seed);
        
        std::normal_distribution<double> distribution_SigV(SigmaV,SigmaV*cv);
        std::normal_distribution<double> distribution_Sigh(Sigmah,Sigmah*cv);
        std::normal_distribution<double> distribution_SigH(SigmaH,SigmaH*cv);
        
        std::cout << "Create Random In-Situ Stresses with condition" << std::endl;
        
        for (int i = 0; i < ncases; i++) {
            // SigV distribution
            Sig_V(i,0) = distribution_SigV(generator);
            // Sigh distribution
            Sig_h(i,0) = distribution_Sigh(generator);
            // SigH distribution
            Sig_H(i,0) = distribution_SigH(generator);
            
            // Check if SigH > SigV > Sigh
            
            while (fabs(Sig_H(i,0)) <= fabs(Sig_V(i,0))) {
                Sig_H(i,0) = distribution_SigH(generator);
                Sig_V(i,0) = distribution_SigV(generator);
            }
            
            while (fabs(Sig_V(i,0)) <= fabs(Sig_h(i,0))) {
                Sig_h(i,0) = distribution_Sigh(generator);
            }
            
            distribution_SigV.reset();
            distribution_Sigh.reset();
            distribution_SigH.reset();
        }
        
        std::ofstream out_VecSigh("Sig_h.txt");
        Sig_h.Print("Sig_h = ", out_VecSigh, EMathematicaInput);
        
        // Exporta In-SItu Stresses
        std::ofstream out_VecSigV("Sig_V.txt");
        Sig_V.Print("Sig_V = ", out_VecSigV, EMathematicaInput);
        
        std::ofstream out_VecSigH_("Sig_H_.txt");
        Sig_H.Print("Sig_H = ", out_VecSigH_, EMathematicaInput);
        
        //
        //    // In-Situ Stresses varying
        //    for(int i=0; i < ncases; i++){
        //
        //        std::cout << "SigV= " << Sig_V(i) << std::endl;
        //        std::cout << "Sigh= " << Sig_h(i) << std::endl;
        //        std::cout << "SigH= " << Sig_H(i) << std::endl;
        //
        //        TPZFMatrix<REAL> Sig_h_Sim (ncases, 1, 0.);
        //        TPZFMatrix<REAL> Sig_H_Sim (ncases, 1, 0.);
        //
        //
        //        if (fabs(Sig_h(i)) > fabs(Sig_H(i))) {
        //            Sig_h_Sim(i) = Sig_H(i);
        //            Sig_H_Sim(i) = Sig_h(i);
        //        }
        //
        //        else{
        //
        //            Sig_h_Sim(i) = Sig_h(i);
        //            Sig_H_Sim(i) = Sig_H(i);
        //        }
        //
        //        std::cout << "SigV= " << Sig_V(i) << std::endl;
        //        std::cout << "Sigh= " << Sig_h_Sim(i) << std::endl;
        //        std::cout << "SigH= " << Sig_H_Sim(i) << std::endl;
        
        
        // Stpchastic cases and In-Situ Stresses varying
        for(int i=0; i < ncases; i++){
            
            std::cout << "Stochastic case = " << i+1 << std::endl;
            
            std::cout << "SigV= " << Sig_V(i) << std::endl;
            std::cout << "Sigh= " << Sig_h(i) << std::endl;
            std::cout << "SigH= " << Sig_H(i) << std::endl;
            
            Problem2D(rw, rext, ncircle, nradial, projection, inclinedwellbore, analytic, Sig_V(i), Sig_h(i),
                      Sig_H(i), Pwb, drdcirc, direction, inclination, isStochastic, solutionfile, i, M);
            
        }

    }
    
    else {
        
        std::cout << "SigV= " << SigmaV << std::endl;
        std::cout << "Sigh= " << Sigmah << std::endl;
        std::cout << "SigH= " << SigmaH << std::endl;
        
        for(int i=0; i < ncases; i++){
            Problem2D(rw, rext, ncircle, nradial, projection, inclinedwellbore, analytic, SigmaV, Sigmah,
                      SigmaH, Pwb, drdcirc, direction, inclination, isStochastic, solutionfile, i, M);
        }
    }
    
    solutionfile.close();
    
    //Problem3D();
    
    //ApproximationRates();
    
    std::cout << "Verifique sempre: q da progressao, scale, e mais importante, se a matrix KCorr calculada se refere a decomposicao que esta sendo lida. Dica: rodar sempre uma vez, decompor no mathematica e rodar novamente" << std::endl;
    
    return 0;
}
