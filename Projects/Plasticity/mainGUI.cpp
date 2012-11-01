#include <iostream>
#include <cstdlib>





#ifdef HAVE_CONFIG_H
#include <config.h>
#endif




 using namespace std;
void MaterialPointTests();

 int main()
 {
 	cout << "QT GUI" << endl;
     MaterialPointTests();
     cout << "\n End runing ! "<< endl;

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
    TPZSandlerDimaggio sandler;
    sandler.McCormicRanchSand(sandler);
    ofstream outfiletxty("SandlerDimaggioApplyStrain.txt");
    ofstream outfiletxtS("SandlerApplyLoad.txt");
    TPZTensor<REAL> deltaeps,eps,sigma,deltasigma;
    
    sandler.fIntegrTol = 1.e-12;
    sandler.SetIntegrTol(1.e-4);
    
    deltaeps.XX()= -0.0013;
    deltaeps.YY()=0;
    deltaeps.ZZ()=0;
    eps=deltaeps;
    
    
    for(int i=0;i<100;i++)
    {
        
        sandler.ApplyStrainComputeSigma(eps, sigma);//UCS
        // sandler.ApplyLoad(sigma, eps);
        
        TPZPlasticState<REAL> state = sandler.GetState();
        //state.Print(cout);
        if(i==49)
        {
            deltaeps*=-1;
        }
        
        outfiletxty << fabs(eps.XX()) << " " << fabs(sigma.XX()) << "\n";
        eps+=deltaeps;
    }
    
    TPZSandlerDimaggio sandler2;
    sandler2.McCormicRanchSand(sandler2);
    sandler2.SetIntegrTol(0.0000001);
    deltaeps.XX()= 1.5*-0.065;
    deltaeps.YY()=0;
    deltaeps.ZZ()=0;
    eps=deltaeps;
    sandler2.ApplyStrainComputeSigma(eps, sigma);
    //outfiletxty << fabs(eps.XX()) << " " << fabs(sigma.XX()) << "\n";
    
    ///TESTE DE COMPRESSAO SIMPLES
    deltasigma.XX()=-0.004;
    deltasigma.YY()=deltasigma.XX()*0.4;
    deltasigma.ZZ()=deltasigma.YY();
    sigma=deltasigma;
    
    TPZSandlerDimaggio sandler3;
    sandler3.McCormicRanchSand(sandler3);
    
    for(int i=0;i<100;i++)
    {
        sandler3.ApplyLoad(sigma,eps);
        outfiletxtS << fabs(eps.XX()) << " " << fabs(sigma.XX()) << "\n";
        sigma+=deltasigma;
    }
    //VisualizeSandlerDimaggio(fileName,pSD);

  
}

