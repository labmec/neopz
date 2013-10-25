#include <iostream>
#include <cstdlib>


#include "TPZPlasticityTest.h"

using namespace std;
void MaterialPointTests();

 int main()
 {
     InitializePZLOG();
     TPZPlasticityTest test;
     std::string filename("../ensaio_all_columns.txt");
     test.ReadInputStrainStress(filename);
     
     TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> sandler;
     //sandler.McCormicRanchSand(sandler);
     REAL poisson, E, A, B, C, R, D, W;
     
     
     //cin >> E;
     E = 29269; //100; //ksi
     
     
     //  cin>> poisson;
     poisson = 0.203; //0.40;
     
     
     A = 616.67; //0.25;
     B = 0.0036895; //0.67;
     C = 111.48; //0.18;
     D = 0.018768; //0.67 / 2.; 
     R = 0.91969; //2.5;
     W = 0.006605; //0.066;
     
     
     sandler.SetUp(poisson, E,A, B, C, R, D, W);
     
     
     REAL inttol = 1.e-4;
     
     sandler.SetIntegrTol(inttol);
     
     test.SetSandlerDimaggio(sandler);
     
     test.SetSimulationInitialStep(0);

     test.PerformSimulation();
     
     
// 	cout << "QT GUI" << endl;
//     MaterialPointTests();
//     cout << "\n End runing ! "<< endl;

 return EXIT_SUCCESS;
 }

#include "TPZPlasticityTest.h"
//
void MaterialPointTests()
{
    
//    cout << "\nChoose Plasticity test:";
//    cout << "\n0 - Isotropic compression ";
//    cout << "\n1 - Biaxial Tests ";
//    cout << "\n2 - Uniaxial traction ";
//    cout << "\n";
//    int choice;
//  //  cin >> choice;
//    choice =0;
//    switch(choice)
//    {
//        case(0):
//            cout << "\n Choose the Plastic model tou need to run Isotropic compression: ";
//            cout << "\n0 - Lade - Kim ";
//            cout << "\n1 - Sandler Dimaggio ";
//            cout << "\n2 - Drucker Prager ";
//            cout << "\n";
//            int choice2;
//         //   cin >> choice2;
//            choice2 = 1;
//            switch(choice2)
//        {
//            case(0):
//                LKIsotropicCompression();
//                break;
//            case(1):
//                SandlerDimaggioIsotropicCompression();
//                break;
//            case(2):
//                DruckerIsotropicCompression();
//                break;
//        }
//            
//            
//            break;
//        case(1):
//            LKBiaxialTest();
//            break;
//        case(2):
//            cout << "NOT IMPLEMENTED YET";
//            //    LadeKim_ReversalTest();
//            break;
//        default:
//            cout << "Unknown Test Type. Exiting...";
//    }
    
    
    ///TESTE UCS
    
    ofstream outUCS("UCS.txt");
    ofstream outConfinamento("Confinamento03.txt");
    TPZTensor<REAL> deltaeps,eps,sigma,deltasigma;
    
    
    TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1> sandler;
    //sandler.McCormicRanchSand(sandler);
    REAL poisson, E, A, B, C, R, D, W;
    
    
    //cin >> E;
    E = 100; //ksi
    
    
  //  cin>> poisson;
    poisson = 0.40;

    
    A = 0.25;
    B = 0.67;
    C = 0.18;
    D = 0.67 / 2.; 
    R = 2.5;
    W = 0.066;
    

    sandler.SetUp(poisson, E,A, B, C, R, D, W);
    
    
    REAL inttol = 1.e-4;
    
    sandler.SetIntegrTol(inttol);
    //IntTolerance
    
    REAL epsx;
    epsx=-0.0013;
    
    deltaeps.XX()= epsx;
    deltaeps.YY()=0;
    deltaeps.ZZ()=0;
    eps=deltaeps;
    
    int numberofsteps=100;
    int unloadstep=49;
    
    for(int i=0;i<numberofsteps;i++)
    {
        
        sandler.ApplyStrainComputeSigma(eps, sigma);//UCS
        // sandler.ApplyLoad(sigma, eps);
        
        TPZPlasticState<REAL> state = sandler.GetState();
        //state.Print(cout);
        if(i==unloadstep)
        {
            deltaeps*=-1;
        }
        
        outUCS << fabs(eps.XX()) << " " << fabs(sigma.XX()) << "\n";
        eps+=deltaeps;
    }
    
    
    ///////////////////FIM UCS////////////////////////
    
    //////// ENSAIO TRIAXIAL /////////////////////////
    
    
    TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1> sandler2;
    
    //cin >> E;
    E = 100; //ksi
    
    //  cin>> poisson;
    poisson = 0.40;
    
    
    A = 0.25;
    B = 0.67;
    C = 0.18;
    D = 0.67 / 2.;
    R = 2.5;
    W = 0.066;
    
    
    sandler2.SetUp(poisson, E,A, B, C, R, D, W);
    
    
    inttol = 1.e-8;
    
    sandler2.SetIntegrTol(inttol);
    
    ///TESTE DE CONFINAMENTO
    
    deltasigma.XX()=-0.004;
    deltasigma.YY()=deltasigma.XX()*0.41;
    deltasigma.ZZ()=deltasigma.YY();
    sigma=deltasigma;

     numberofsteps=100;
    
    for(int i=0;i<numberofsteps;i++)
    {
        sandler2.ApplyLoad(sigma,eps);
        outConfinamento << fabs(eps.XX()) << " " << fabs(sigma.XX()) << "\n";
        sigma+=deltasigma;
    }
    //VisualizeSandlerDimaggio(fileName,pSD);

  
}

