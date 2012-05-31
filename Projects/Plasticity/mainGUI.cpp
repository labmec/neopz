#include <iostream>
#include <cstdlib>
 using namespace std;
void MaterialPointTests();

 int main(int argc, char *argv[])
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
    
    cout << "\nChoose Plasticity test:";
    cout << "\n0 - Isotropic compression ";
    cout << "\n1 - Biaxial Tests ";
    cout << "\n2 - Uniaxial traction ";
    cout << "\n";
    int choice;
  //  cin >> choice;
    choice =0;
    switch(choice)
    {
        case(0):
            cout << "\n Choose the Plastic model tou need to run Isotropic compression: ";
            cout << "\n0 - Lade - Kim ";
            cout << "\n1 - Sandler Dimaggio ";
            cout << "\n2 - Drucker Prager ";
            cout << "\n";
            int choice2;
         //   cin >> choice2;
            choice2 = 1;
            switch(choice2)
        {
            case(0):
                LKIsotropicCompression();
                break;
            case(1):
                SandlerDimaggioIsotropicCompression();
                break;
            case(2):
                DruckerIsotropicCompression();
                break;
        }
            
            
            break;
        case(1):
            LKBiaxialTest();
            break;
        case(2):
            cout << "NOT IMPLEMENTED YET";
            //    LadeKim_ReversalTest();
            break;
        default:
            cout << "Unknown Test Type. Exiting...";
    }
    
}

