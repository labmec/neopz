// $Id: TPZPlasticityTest.h,v 1.32 2010-08-20 20:56:45 diogo Exp $
// This file implements the interfaces for the
// most general plasticity load tests

#include <iostream>

using namespace std;

#include "pzlog.h"

#ifdef LOG4CXX // LOG4CXX may be defined alone or with LOG4CXX_PLASTICITY. The latter shall not be used alone.
#include <log4cxx/logger.h>
#include <log4cxx/basicconfigurator.h>
#include <log4cxx/propertyconfigurator.h>

static LoggerPtr plasticIntegrLogger(Logger::getLogger("plasticity.plasticIntegr"));
#endif


#ifdef LOG4CXX_PLASTICITY
static LoggerPtr testLogger(Logger::getLogger("plasticity.test"));
#endif



#ifdef LOG4CXX_PLASTICITY
static LoggerPtr MaterialPoint(Logger::getLogger("MaterialPointTest"));
#endif

#include "pzvec.h"
#include "pzfmatrix.h"
#include "TPZTensor.h" 

#include "TPZLadeKim.h"
#include "TPZSandlerDimaggio.h"

#include "TPZYCDruckerPrager.h"
#include "TPZThermoForceA.h"
#include "TPZElasticResponse.h"


#include "TPZYCMohrCoulomb.h"
#include "TPZYCModifiedMohrCoulomb.h"

#include "TPZYCWillamWarnke.h"

#include "TPZYCVonMises.h"
#include "TPZVonMises.h"
#include "TPZDruckerPrager.h"



class TPZPlasticTest
	{
	public:
		
		static void InitializeLOG();
		
		template <class T>
		static void ReciprocityTest(T & plasticModel, TPZTensor<REAL> strain1);
		
		template <class T>
		static void StressTest(T & plasticModel, const char * filename, REAL  stressMultiplier=1);
		
		template <class T>
		static void StrainTest(T & plasticModel, const char * filename, REAL  strainMultiplier=1);
		
		static void LoadTest(const char * filename);
		
		template <class T>
		static int CreatePlasticModel(T * ( & plasticModel), const char * line);
		
		template <class T>
		static void GlobalCheckConv(T & plasticModel, TPZTensor<REAL> & strain, REAL maxDeltaStrain = 0.01);
		
		static void DruckerPragerTest();
		
		static void MohrCoulombTest();
		
		static void ModifiedMohrCoulombTest();
		
		static void WillamWarnkeTest();
		
		static void VonMisesTest();
		
		static void UndocumentedTest2();
		
		static void UndocumentedTest3();
		
		static void UndocumentedTest4();
        
        static void SandlerDimaggioIsotropicCompression();
        
        static void LKFineSilicaLoadTest();
		
        static void LKIsotropicCompression();
        
        static void LKKoCompressionLoadTest();
        
        static void LKLoadingTest();
        
        static void DruckerIsotropicCompression();
        
        static void DruckerTest();
        
        static void LKBiaxialTest();
        
        //static void MaterialPointTests();
		//////////////////CheckConv related methods/////////////////////
		
		/**
		 number of types of residuals
		 */
		int NumCases() 
		{
			return 9;
		}
		
		TPZTensor<REAL> gRefTension;
		
	//	TPZMatElastoPlastic<TPZMaterial> mate;
		TPZVonMises gPlasticModel;
		
		/**
			LoadState will keep a given state as static variable of the class
		*/
		
		void LoadState(TPZFMatrix &state)
		{
			int i;
			for(i=0; i<6; i++) gRefTension.fData[i] = state(i,0);
		}
		
		void ComputeTangent(TPZFMatrix &tangent, TPZVec<REAL> &coefs, int icase)
		{
			switch(icase)
			{
				case 0:
				{
					TPZTensor<REAL> grad,EpsT,DiagonalStress;
					TPZFNMatrix<6*6> Dep(6,6,0.);
					gPlasticModel.ApplyStrainComputeDep(gRefTension,DiagonalStress,Dep);
					tangent.Redim(1,1);
					REAL norm = Norm(Dep);
					tangent(0,0) = norm;

					break;
				}
					
			}
		}
		
		void Residual(TPZFMatrix &res,int icase)
		{
			
			res.Redim(1,1);
			TPZTensor<REAL> grad,DiagonalStress;
			TPZFNMatrix<36> Dep(6,6,0.);
			gPlasticModel.ApplyStrainComputeDep(gRefTension,DiagonalStress,Dep);
			res.Redim(1,1);
			REAL norm = DiagonalStress.Norm();
			switch(icase)
			{
				case 0:
				{
					res(0,0) = norm;
					break;
				}

			}
			
		}
		
		
		
		
		static void RotateMatrix(TPZFMatrix &Mat, double thetaRad,int rotateaboutaxes);
		
		//static void RotationMatrix(TPZFMatrix &R, double thetaRad, int axis);
		
		template <class T>
		static void PlasticIntegratorCheck(int thetaintervals,T mat);
		
		static void VerifyIntegrationAtPoint(TPZVec< TPZTensor<REAL> > vectensor);
	//	static void DruckerTest();
        
        
	};


//inline void MaterialPointTests()
//{
//    
//    cout << "\nChoose Plasticity test:";
//    cout << "\n0 - Isotropic compression ";
//    cout << "\n1 - Biaxial Tests ";
//    cout << "\n2 - Uniaxial traction ";
//    cout << "\n";
//    int choice;
//    cin >> choice;
//    
//    switch(choice)
//    {
//        case(0):
//            cout << "\n Choose the Plastic model tou need to run Isotropic compression: ";
//            cout << "\n0 - Lade - Kim ";
//            cout << "\n1 - Sandler Dimaggio ";
//            cout << "\n2 - Drucker Prager ";
//            cout << "\n";
//            int choice2;
//            cin >> choice2;
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
//    
//}



inline void SandlerDimaggioIsotropicCompression()//
{
    ofstream outfiletxt5("SDMcCormicRanchSand.txt");
    TPZTensor<REAL> stress, strain, deltastress, deltastrain;
    TPZFNMatrix<6*6> Dep(6,6,0.);
    
    
    TPZSandlerDimaggio SD;
    
    cout << "\n Put the value of strain you want to add in each step of your loat test: (sugg.0.0001) ";
    REAL straininput;
    cin >> straininput;
    deltastrain.XX() = -straininput;
    deltastrain.XY() = 0.;
    deltastrain.XZ() = 0.;
    deltastrain.YY() = -straininput;
    deltastrain.YZ() = 0.;
    deltastrain.ZZ() = -straininput;
    strain=deltastrain;    
    cout << "Choose the material pareameters you want to set to SandlerDimaggio Test :";
    cout << "\n0 - McCormicRanchSandMod";
    cout << "\n1 - McCormicRanchSandMod2 ";
    cout << "\n2 - UncDeepSandRes ";
    cout << "\n3 - UncDeepSandResPSI";
    cout << "\n4 - UncDeepSandResMPa";
    cout << "\n5 - Put the material parameters you want ";
    int choice;
    cin >> choice;
    
    cout << "\n Put the numbers of steps you want: (sugg. 20)";
    int length;
    cin >> length;
    
    switch (choice) {
        case(0):
        {
            TPZSandlerDimaggio::McCormicRanchSandMod(SD);
            std::ofstream outfiletxt("TPZSandlerDimaggioMcCormicRanchSandMod(SD).txt");
            
            for(int step=0;step<length;step++)
            {
                cout << "\nstep "<< step;
                SD.ApplyStrainComputeDep(strain, stress,Dep);
                outfiletxt << fabs(strain.XX()) << " " << fabs(stress.XX()) << "\n";
                strain += deltastrain;
            }
            break;
        }
        case(1):
        {
            TPZSandlerDimaggio::McCormicRanchSandMod2(SD);
            std::ofstream outfiletxt("TPZSandlerDimaggioMcCormicRanchSandMod2.txt");
            for(int step=0;step<length;step++)
            {
                cout << "\nstep "<< step;
                SD.ApplyStrainComputeDep(strain, stress,Dep);
                outfiletxt << fabs(strain.XX()) << " " << fabs(stress.XX()) << "\n";
                strain += deltastrain;
                
            }
            break;
        }
        case(2):
        {
            TPZSandlerDimaggio::UncDeepSandRes(SD);
            std::ofstream outfiletxt("TPZLadeKimUncDeepSandRes.txt");
            for(int step=0;step<length;step++)
            {
                cout << "\nstep "<< step;
                SD.ApplyStrainComputeDep(strain, stress,Dep);
                outfiletxt << fabs(strain.XX()) << " " << fabs(stress.XX()) << "\n";
                strain += deltastrain;
                
            }
            break;
        }
        case(3):
        {
            TPZSandlerDimaggio::UncDeepSandResPSI(SD);
            std::ofstream outfiletxt("TPZLadeKimUncDeepSandResPSI.txt");
            for(int step=0;step<length;step++)
            {
                cout << "\nstep "<< step;
                SD.ApplyStrainComputeDep(strain, stress,Dep);
                outfiletxt << fabs(strain.XX()) << " " << fabs(stress.XX()) << "\n";
                strain += deltastrain;
                
            }
            break;
        }
        case(4):
        {
            TPZSandlerDimaggio::UncDeepSandResMPa(SD);
            std::ofstream outfiletxt("TPZLadeKimUncDeepSandResMPa.txt");
            for(int step=0;step<length;step++)
            {
                cout << "\nstep "<< step;
                SD.ApplyStrainComputeDep(strain, stress,Dep);
                outfiletxt << fabs(strain.XX()) << " " << fabs(stress.XX()) << "\n";
                strain += deltastrain;
                
            }
            break;
            
            break;
            
        }    
        case(5):
        {
            
//            REAL E = 9000, 
//            poisson = 0.25;
//            
//            material.fER.SetUp(E, poisson);
//            
//            REAL A = 18, 
//            B = 0.0245, 
//            C = 17.7, 
//            D = 0.00735, 
//            R = 1.5,
//            W = 0.0908;
            cout<< "\n Young Modulus 9000.";
            REAL E; 
            cin >> E;
            
            cout<< "\n poisson 0.25 ";
            REAL poisson;
            cin >> poisson;
            
            
            SD.fER.SetUp(E, poisson);
            
            cout<< "\n A (sugg.18)"; 
            REAL A;
            cin >> A;
            
            cout << "\n B (sugg. 0.0245) ";
            REAL B;
            cin >> B;
            
            cout << "\n C (sugg. 17.7) ";
            REAL C;
            cin >> C;
            
            cout << "\n D (sugg.0.00735) ";
            REAL D;
            cin >> D;
            
            cout << "\n R (sugg. 1.5) ";
            REAL R;
            cin >> R;
            
            cout << "\n W (sugg. 0.0908) ";
            REAL W;
            cin >> W;
            
            SD.fYC.SetUp(A, B, C, D, R, W);
            
            std::ofstream outfiletxt("SandlerDimaggioYOURMODEL.txt");
            for(int step=0;step<length;step++)
            {
                cout << "\nstep "<< step;
                SD.ApplyStrainComputeDep(strain, stress,Dep);
                outfiletxt << fabs(strain.XX()) << " " << fabs(stress.XX()) << "\n";
                strain += deltastrain;
                
            }
            break;
        }
        default:
        {
            cout << "Unknown Test Type. Exiting...";
            break;
        }
    }
    
    
    
}

inline void LKFineSilicaLoadTest()//
{
    ofstream outfiletxt1("FineSilica.txt");
    TPZLadeKim LK;
    TPZLadeKim::FineSilicaSand(LK);
    TPZTensor<REAL> stress, strain, deltastress, deltastrain;
    TPZFNMatrix<6*6> Dep(6,6,0.);
    deltastress.XX() = -4.;
    deltastress.XY() = 0.;
    deltastress.XZ() = 0.;
    deltastress.YY() = -4.;
    deltastress.YZ() = 0.;
    deltastress.ZZ() = -4.;
    stress = deltastress;
    
    
    int length =30;
    for(int step=0;step<length;step++)
    {
        cout << "\nstep "<< step;
        if(step == 69)deltastress *=-1.;
        LK.ApplyLoad(stress,strain);            
        outfiletxt1 << fabs(strain.XX()) << " " << fabs(stress.XX()) << "\n";
        stress +=  deltastress;
        cout << "strain = "<<strain <<"\n";
        cout << "sigma = "<< stress <<"\t "<< "I1 = "<< stress.I1() <<"\n";
        
    }
    
    
    
    
}

inline void LKIsotropicCompression()
{
    
    TPZTensor<REAL> stress, strain, deltastress, deltastrain;
    TPZFNMatrix<6*6> Dep(6,6,0.);
    
    cout << "\n Put the value of strain you want to add in each step of your loat test: (sugg. 0.0001) ";
    REAL straininput;
    cin >> straininput;
    
    deltastrain.XX() = -straininput;
    deltastrain.XY() = 0.;
    deltastrain.XZ() = 0.;
    deltastrain.YY() = -straininput;
    deltastrain.YZ() = 0.;
    deltastrain.ZZ() = -straininput;
    strain=deltastrain;
    
    cout << "Choose the material pareameters you want to set to Lade Kim Test :";
    cout << "\n0 - Plain Concrete ";
    cout << "\n1 - Loose Sacramento River Sand ";
    cout << "\n2 - Dense Sacramento River Sand ";
    cout << "\n3 - Fine Silica Sand";
    cout << "\n4 - Put the material parameters you want";
    int choice;
    cin >> choice;
    
    cout << "\n Put the numbers of steps you want:(sugg. 20)";
    int length;
    cin >> length;
    
    TPZLadeKim LK2;
    switch (choice) {
        case(0):
        {
            TPZLadeKim::PlainConcrete(LK2);
            std::ofstream outfiletxt("TPZLadeKim::PlainConcrete.txt");
            
            for(int step=0;step<length;step++)
            {
                cout << "\nstep "<< step;
                LK2.ApplyStrainComputeDep(strain, stress,Dep);
                outfiletxt << fabs(strain.XX()) << " " << fabs(stress.XX()) << "\n";
                strain += deltastrain;
            }
            break;
        }
        case(1):
        {
            TPZLadeKim::LooseSacrRiverSand(LK2);
            std::ofstream outfiletxt("TPZLadeKim::LooseSacrRiverSand.txt");
            for(int step=0;step<length;step++)
            {
                cout << "\nstep "<< step;
                LK2.ApplyStrainComputeDep(strain, stress,Dep);
                outfiletxt << fabs(strain.XX()) << " " << fabs(stress.XX()) << "\n";
                strain += deltastrain;
                
            }
            break;
        }
        case(2):
        {
            TPZLadeKim::DenseSacrRiverSand(LK2);
            std::ofstream outfiletxt("TPZLadeKim::DenseSacrRiverSand.txt");
            for(int step=0;step<length;step++)
            {
                cout << "\nstep "<< step;
                LK2.ApplyStrainComputeDep(strain, stress,Dep);
                outfiletxt << fabs(strain.XX()) << " " << fabs(stress.XX()) << "\n";
                strain += deltastrain;
                
            }
            break;
        }
        case(3):
        {
            TPZLadeKim::FineSilicaSand(LK2);
            std::ofstream outfiletxt("TPZLadeKim::FineSilicaSand.txt");
            for(int step=0;step<length;step++)
            {
                cout << "\nstep "<< step;
                LK2.ApplyStrainComputeDep(strain, stress,Dep);
                outfiletxt << fabs(strain.XX()) << " " << fabs(stress.XX()) << "\n";
                strain += deltastrain;
                
            }
            break;
        }
        case(4):
        {
            cout << "\n poisson (sugg. 0.18)";
            REAL poisson;// = 0.18;
            cin>>poisson; 
            
            cout << "\n M (sugg. 361800.)";
            REAL M;//       = 361800.;
            cin >> M;
            
            cout << "\nlambda (sugg. 0.)";
            REAL lambda;//  = 0.;
            cin >> lambda;
            
            cout << "\n a (sugg. 28.5)";
            REAL a;//       = 28.5;
            cin >> a;
            
            cout << "\n m (sugg. 1.113)";
            REAL m;//       = 1.113;
            cin >> m;
            
            cout << "\n neta1 (sugg. 159800.)";
            REAL neta1;//   = 159800.;
            cin >> neta1;
            
            cout << "\n ksi2 (sugg. -2.92)";
            REAL ksi2; //   = -2.92;
            cin >> ksi2;
            
            cout << "\n mu (sugg. 5.06)";
            REAL mu;//     = 5.06;
            cin >> mu;
            
            cout << "\n C (sugg. 0.712E-12)";
            REAL C;//       = 0.712E-12;
            cin >> C;
            
            cout << "\n p (sugg. 3.8)";
            REAL p;//       = 3.8;
            cin >> p;
            
            cout <<"\n h (sugg. 1.990) ";
            REAL h;//       = 1.990;
            cin >> h;
            
            cout << "\n alpha (sugg. 0.75) ";
            REAL alpha;//   = 0.75;
            cin >> alpha;
            
            cout << "\n pa (sugg. 14.7) ";
            REAL pa;//      = 14.7;
            cin >> pa;
            
            REAL restol;
            cout << "\n Tolerance (sugg. 0.0001) ";
            cin >> restol;
            
            LK2.fResTol = restol;
            
            LK2.SetUp(poisson, M, lambda,
                      a, m, neta1,
                      ksi2, mu,
                      C, p,
                      h, alpha,
                      pa);
            std::ofstream outfiletxt("TPZLadeKim::YOURMODEL.txt");
            for(int step=0;step<length;step++)
            {
                cout << "\nstep "<< step;
                LK2.ApplyStrainComputeDep(strain, stress,Dep);
                outfiletxt << fabs(strain.XX()) << " " << fabs(stress.XX()) << "\n";
                strain += deltastrain;
                
            }
            break;
        }
        default:
        {
            cout << "Unknown Test Type. Exiting...";
            break;
        }
    }
    
}


inline void LKKoCompressionLoadTest()
{
    ofstream outfiletxt1("LKKoCompressionLoadTest.txt");
    TPZTensor<REAL> stress, strain, deltastress, deltastrain;
    TPZFNMatrix<6*6> Dep(6,6,0.);
    
    deltastress.XX()=-1.;
    deltastress.YY()=-0.5;
    deltastress.ZZ()=-0.5;
    stress=deltastress;
    
    
    TPZLadeKim LK2;
    TPZLadeKim::FineSilicaSand(LK2);
    
    int length2 =20;
    for(int step=0;step<length2;step++)
    {
        cout << "\nstep "<< step;
        LK2.ApplyLoad(stress,strain);
        LK2.ApplyStrainComputeDep(strain, stress, Dep);
        outfiletxt1 << fabs(strain.I1()) << " " << fabs(stress.XX()/stress.ZZ()) << "\n";
        stress += deltastress;
        cout << "strain = "<<strain <<"\n";
        cout << "sigma = "<< stress <<"\t "<< "I1 = "<< stress.I1() <<"\n";
        
    }
    
}


inline void LKLoadingTest()
{
    
    ofstream outfiletxt1("FineSilica.txt");
    ofstream outfiletxt2("PlainConcretee1.txt");
    ofstream outfiletxt3("PlainConcretee2.txt");
    ofstream outfiletxt4("PlainConcretee3.txt");
    ofstream outfiletxt5("SDMcCormicRanchSand.txt");
    TPZTensor<REAL> stress, strain, deltastress, deltastrain;
    TPZFNMatrix<6*6> Dep(6,6,0.);
    
    
    
    
    deltastress.XX() = -0.001;
    deltastress.XY() = 0.;
    deltastress.XZ() = 0.;
    deltastress.YY() = -0.001;
    deltastress.YZ() = 0.;
    deltastress.ZZ() = -0.001;
    stress = deltastress;
    
    TPZLadeKim LK;
    
    TPZLadeKim::FineSilicaSand(LK);
    //    LK.ApplyLoad(stress,deltastrain);
    
    deltastress.XX() = -4.;
    deltastress.XY() = 0.;
    deltastress.XZ() = 0.;
    deltastress.YY() = -4.;
    deltastress.YZ() = 0.;
    deltastress.ZZ() = -4.;
    stress = deltastress;
    
    
    //    deltastrain.XX() = -0.0001;
    //    deltastrain.XY() = 0.;
    //    deltastrain.XZ() = 0.;
    //    deltastrain.YY() = -0.0001;
    //    deltastrain.YZ() = 0.;
    //    deltastrain.ZZ() = -0.0001;
    //    strain=deltastrain;
    
    int length =72;
    for(int step=0;step<length;step++)
    {
        cout << "\nstep "<< step;
        if(step == 69)deltastress *=-1.;
        LK.ApplyLoad(stress,strain);
        //        if(step == 16)deltastrain *=-1.;
        //        LK.ApplyStrainComputeDep(strain, stress,Dep);
        if(step==0)
        {
            
            //outfiletxt1 << 0. << " " << 0. << "\n";
            
        }
        
        else
        {
            
            outfiletxt1 << fabs(strain.XX()) << " " << fabs(stress.XX()) << "\n";
            
        }
        
        //        deltastress.Multiply(1.1, 1.);
        stress +=  deltastress;
        //        strain +=deltastrain;
        cout << "strain = "<<strain <<"\n";
        cout << "sigma = "<< stress <<"\t "<< "I1 = "<< stress.I1() <<"\n";
        
    }
    
    
    deltastrain.XX() = -0.0001;
    deltastrain.XY() = 0.;
    deltastrain.XZ() = 0.;
    deltastrain.YY() = -0.;
    deltastrain.YZ() = 0.;
    deltastrain.ZZ() = -0.;
    strain=deltastrain;
    TPZLadeKim LK2;
    TPZLadeKim::PlainConcrete(LK2);
    
    int length2 =30;
    for(int step=0;step<length2;step++)
    {
        cout << "\nstep "<< step;
        //    if(step == 10 || step == 16|| step == 40 || step==51 || step==80)deltastrain *=-1.;
        LK2.ApplyStrainComputeDep(strain, stress,Dep);
        outfiletxt2 << fabs(strain.XX()) << " " << fabs(stress.XX()) << "\n";
        outfiletxt3 << fabs(strain.YY()) << " " << fabs(stress.XX()) << "\n";
        outfiletxt4 << fabs(strain.ZZ()) << " " << fabs(stress.XX()) << "\n";
        strain += deltastrain;
        cout << "strain = "<<strain <<"\n";
        cout << "sigma = "<< stress <<"\t "<< "I1 = "<< stress.I1() <<"\n";
        
    }
    
    TPZSandlerDimaggio SD;
    TPZSandlerDimaggio::McCormicRanchSand(SD);
    
    deltastrain.XX() = -0.005;
    deltastrain.XY() = 0.;
    deltastrain.XZ() = 0.;
    deltastrain.YY() = -0.;
    deltastrain.YZ() = 0.;
    deltastrain.ZZ() = -0.;
    strain=deltastrain;
    
    int length3 =23;
    for(int step=0;step<length3;step++)
    {
        cout << "\nstep "<< step;
        if(step == 14 )deltastrain *=-1.;
        SD.ApplyStrainComputeDep(strain, stress,Dep);
        outfiletxt5 << fabs(strain.XX()) << " " << fabs(stress.XX()) << "\n";
        strain += deltastrain;
        cout << "strain = "<<strain <<"\n";
        cout << "sigma = "<< stress <<"\t "<< "I1 = "<< stress.I1() <<"\n";
        
    }    
    
    
}



//int main()
//{
//
//    
//    
//    
//    PressureCilinder();
//    return 0;
//}


inline void DruckerIsotropicCompression()
{
    TPZDruckerPrager DP;
    ofstream outfiletxt("DruckerPragerIsotropicCompression.txt");
    cout << "\n Put the value of strain you want to add in each step of your loat test: (sugg. 0.0001)";
    REAL straininput;
    cin >> straininput;
    TPZFNMatrix<6*6> Dep(6,6,0.);
    TPZTensor<REAL> deltastrain,strain,stress,deltastress;
    deltastrain.XX() = -straininput;
    deltastrain.XY() = 0.;
    deltastrain.XZ() = 0.;
    deltastrain.YY() = -straininput;
    deltastrain.YZ() = 0.;
    deltastrain.ZZ() = -straininput;
    strain=deltastrain;
    
    cout << "\n4 - Put the material parameters you want";
    
    cout << "\n Young modulus (sugg. 20000.)";
    REAL E;
    cin >> E;
    
    cout << "\n Poisson (sugg. 0.2)";
    REAL poisson;
    cin >> poisson;
    
    int mcfit;
    cout << "\n choose 0 for Iner Morh-Coulomb fit or 1 for outer Morh-Coulomb Fit (sugg. 0) ";
    cin >> mcfit;
    
    if(mcfit!= 0 || mcfit!= 1)
    {
        cout << "\n wrong choice in Morh-Coulomb fit tipe 0 or 1";
        return;
    }
    
    REAL phi;
    cout << "\n Type the internal frictional angle in degrees(sugg. 20.)";
    cin >> phi;
    
    REAL c;
    cout << "\n Type the material coesion (sugg. 9.)";
    cin >> c;
    
    REAL h;
    cout << "\n Type the material hardening modulus (sugg. 1000.)";
    cin >> h;
    
    DP.fYC.SetUp(phi/180. *M_PI ,mcfit);
    DP.fTFA.SetUp(c,h);
    DP.fER.SetUp(E,poisson);
    
    int length =30;
    for(int step=0;step<length;step++)
    {
        cout << "\nstep "<< step;    
        DP.ApplyStrainComputeDep(strain, stress, Dep);
        strain += deltastrain;
        outfiletxt << fabs(strain.XX()) << " " << fabs(stress.XX()) << "\n";
    }
    
}

inline void LKBiaxialTest()
{
    TPZTensor<REAL> stress, strain, deltastress, deltastrain;
    TPZFNMatrix<6*6> Dep(6,6,0.);
    
    cout << "\n Put the value of stress you want to put sx ";
    REAL stressinputx;
    cin >> stressinputx;
    
    cout << "\n Put the value of stress you want to put sy ";
    REAL stressinputy;
    cin >> stressinputy;
    
    cout << "\n Put the value of stress you want to put sz ";
    REAL stressinputz;
    cin >> stressinputz;
    
    deltastress.XX() = -stressinputx;
    deltastress.XY() = 0.;
    deltastress.XZ() = 0.;
    deltastress.YY() = -stressinputy;
    deltastress.YZ() = 0.;
    deltastress.ZZ() = -stressinputz;
    stress=deltastress;
    
    cout << "Choose the material pareameters you want to set to Lade Kim Test :";
    cout << "\n0 - Plain Concrete ";
    cout << "\n1 - Loose Sacramento River Sand ";
    cout << "\n2 - Dense Sacramento River Sand ";
    cout << "\n3 - Fine Silica Sand";
    cout << "\n4 - Put the material parameters you want";
    int choice;
    cin >> choice;
    
    cout << "\n Put the numbers of steps you want: (sugg. 20)";
    int length;
    cin >> length;
    
    TPZLadeKim LK2;
    switch (choice) {
        case(0):
        {
            TPZLadeKim::PlainConcrete(LK2);
            std::ofstream outfiletxt1("BiaxialXTPZLadeKimPlainConcrete.txt");
            std::ofstream outfiletxt2("BiaxialYTPZLadeKimPlainConcrete.txt");
            std::ofstream outfiletxt3("BiaxialZTPZLadeKimPlainConcrete.txt");
            
            for(int step=0;step<length;step++)
            {
                cout << "\nstep "<< step;
                LK2.ApplyLoad(stress, strain);
                outfiletxt1 << strain.XX() << " " << fabs(stress.XX()) << "\n";
                outfiletxt2 << strain.YY() << " " << fabs(stress.XX()) << "\n";
                outfiletxt3 << strain.ZZ() << " " << fabs(stress.XX()) << "\n";
                stress += deltastress;
            }
            break;
        }
        case(1):
        {
            TPZLadeKim::LooseSacrRiverSand(LK2);
            std::ofstream outfiletxt1("BiaxialXTPZLadeKimLooseSacrRiverSand.txt");
            std::ofstream outfiletxt2("BiaxialYTPZLadeKimLooseSacrRiverSand.txt");
            std::ofstream outfiletxt3("BiaxialZTPZLadeKimLooseSacrRiverSand.txt");
            
            for(int step=0;step<length;step++)
            {
                cout << "\nstep "<< step;
                LK2.ApplyLoad(stress, strain);
                outfiletxt1 << strain.XX() << " " << fabs(stress.XX()) << "\n";
                outfiletxt2 << strain.YY() << " " << fabs(stress.XX()) << "\n";
                outfiletxt3 << strain.ZZ() << " " << fabs(stress.XX()) << "\n";
                stress += deltastress;
            }
            break;
        }
        case(2):
        {
            TPZLadeKim::DenseSacrRiverSand(LK2);
            std::ofstream outfiletxt("TPZLadeKim::DenseSacrRiverSand.txt");
            for(int step=0;step<length;step++)
            {
                cout << "\nstep "<< step;
                LK2.ApplyStrainComputeDep(strain, stress,Dep);
                outfiletxt << fabs(strain.XX()) << " " << fabs(stress.XX()) << "\n";
                outfiletxt << fabs(strain.YY()) << " " << fabs(stress.XX()) << "\n";
                outfiletxt << fabs(strain.ZZ()) << " " << fabs(stress.XX()) << "\n";
                strain += deltastrain;
                
            }
            break;
        }
        case(3):
        {
            TPZLadeKim::FineSilicaSand(LK2);
            std::ofstream outfiletxt("TPZLadeKim::FineSilicaSand.txt");
            for(int step=0;step<length;step++)
            {
                cout << "\nstep "<< step;
                LK2.ApplyStrainComputeDep(strain, stress,Dep);
                outfiletxt << fabs(strain.XX()) << " " << fabs(stress.XX()) << "\n";
                outfiletxt << fabs(strain.YY()) << " " << fabs(stress.XX()) << "\n";
                outfiletxt << fabs(strain.ZZ()) << " " << fabs(stress.XX()) << "\n";
                strain += deltastrain;
                
            }
            break;
        }
        case(4):
        {
            cout << "\n poisson (sugg. 0.18)";
            REAL poisson;// = 0.18;
            cin>>poisson; 
            
            cout << "\n M (sugg. 361800.)";
            REAL M;//       = 361800.;
            cin >> M;
            
            cout << "\nlambda (sugg. 0.)";
            REAL lambda;//  = 0.;
            cin >> lambda;
            
            cout << "\n a (sugg. 28.5)";
            REAL a;//       = 28.5;
            cin >> a;
            
            cout << "\n m (sugg. 1.113)";
            REAL m;//       = 1.113;
            cin >> m;
            
            cout << "\n neta1 (sugg. 159800.)";
            REAL neta1;//   = 159800.;
            cin >> neta1;
            
            cout << "\n ksi2 (sugg. -2.92)";
            REAL ksi2; //   = -2.92;
            cin >> ksi2;
            
            cout << "\n mu (sugg. 5.06)";
            REAL mu;//     = 5.06;
            cin >> mu;
            
            cout << "\n C (sugg. 0.712E-12)";
            REAL C;//       = 0.712E-12;
            cin >> C;
            
            cout << "\n p (sugg. 3.8)";
            REAL p;//       = 3.8;
            cin >> p;
            
            cout <<"\n h (sugg. 1.990) ";
            REAL h;//       = 1.990;
            cin >> h;
            
            cout << "\n alpha (sugg. 0.75) ";
            REAL alpha;//   = 0.75;
            cin >> alpha;
            
            cout << "\n pa (sugg. 14.7) ";
            REAL pa;//      = 14.7;
            cin >> pa;
            
            REAL restol;
            cout << "\n Tolerance (sugg. 0.0001) ";
            cin >> restol;
            
            LK2.fResTol = restol;
            
            LK2.SetUp(poisson, M, lambda,
                      a, m, neta1,
                      ksi2, mu,
                      C, p,
                      h, alpha,
                      pa);
            std::ofstream outfiletxt("TPZLadeKim::YOURMODEL.txt");
            for(int step=0;step<length;step++)
            {
                cout << "\nstep "<< step;
                LK2.ApplyStrainComputeDep(strain, stress,Dep);
                outfiletxt << fabs(strain.XX()) << " " << fabs(stress.XX()) << "\n";
                outfiletxt << fabs(strain.YY()) << " " << fabs(stress.XX()) << "\n";
                outfiletxt << fabs(strain.ZZ()) << " " << fabs(stress.XX()) << "\n";
                strain += deltastrain;
                
            }
            break;
        }
        default:
        {
            cout << "Unknown Test Type. Exiting...";
            break;
        }
    }
    
    
    
    
}


//inline void MaterialPointTests()
//{
//    
//    cout << "\nChoose Plasticity test:";
//    cout << "\n0 - Isotropic compression ";
//    cout << "\n1 - Biaxial Tests ";
//    cout << "\n2 - Uniaxial traction ";
//    cout << "\n";
//    int choice;
//    cin >> choice;
//    
//    switch(choice)
//    {
//        case(0):
//            cout << "\n Choose the Plastic model tou need to run Isotropic compression: ";
//            cout << "\n0 - Lade - Kim ";
//            cout << "\n1 - Sandler Dimaggio ";
//            cout << "\n2 - Drucker Prager ";
//            cout << "\n";
//            int choice2;
//            cin >> choice2;
//            switch(choice2)
//        {
//            case(0):
//                TPZPlasticTest::LKIsotropicCompression();
//                break;
//            case(1):
//                TPZPlasticTest::SandlerDimaggioIsotropicCompression();
//                break;
//            case(2):
//                TPZPlasticTest::DruckerIsotropicCompression();
//                break;
//        }
//            
//            
//            break;
//        case(1):
//            TPZPlasticTest::LKBiaxialTest();
//            break;
//        case(2):
//            cout << "NOT IMPLEMENTED YET";
//            //    LadeKim_ReversalTest();
//            break;
//        default:
//            cout << "Unknown Test Type. Exiting...";
//    }
//    
//}



inline void TPZPlasticTest::InitializeLOG()
{
#ifdef LOG4CXX
	
	std::string path;
	std::string configfile;
#ifdef HAVE_CONFIG_H
	path = PLASTICITYSOURCEDIR;
	path += "/src/";
	cout << path.c_str() << endl;
	cout.flush();
#else
	path = "";
#endif
	configfile = path;
	
	configfile += "log4cxx.cfg";
	log4cxx::PropertyConfigurator::configure(configfile.c_str());
	
	std::stringstream sout;
	sout << __PRETTY_FUNCTION__ << "\nLOG4CXX configured.\n"
	<< "LOG4CXX config file:" << configfile;
	
	LOGPZ_INFO(plasticIntegrLogger,sout.str().c_str());
	
#ifdef LOG4CXX_PLASTICITY
	LOGPZ_INFO(testLogger,sout.str().c_str());
#endif
	
#endif
}


template <class T>
inline void TPZPlasticTest::ReciprocityTest(T & plasticModel, TPZTensor<REAL> strain1)
{
	
	TPZPlasticTest::InitializeLOG();
	
	T plasticModelCopy = plasticModel;	
	
	TPZTensor<REAL> stress1, stress2, strain2;
	TPZFNMatrix < 6*6 > tangent (6,6,0.);
	
	TPZVec<REAL> phi(2,0.), phi2(2,0.);
	
	plasticModel.Phi(strain1,phi);
	
#ifdef LOG4CXX_PLASTICITY
	{
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ 
		<< "\nImposing strain = \n" << strain1
		<< "\nPhi before plastic loop = \n" << phi ;
		LOGPZ_INFO(testLogger,sout.str().c_str());
		
	}
#endif
	
	// applying the strain to the element and processing plastic loop
	
	plasticModel.ApplyStrainComputeDep(strain1, stress1, tangent);
	
	plasticModel.Phi(strain1,phi);
	
#ifdef LOG4CXX_PLASTICITY
	{
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__
		<< "\nphi  after plastic loop =\n" << phi
		<< "\ncausing stress = " << stress1
		<< "\nNow applying the opposite operation...";
		LOGPZ_INFO(testLogger,sout.str().c_str());
		
	}
#endif
	
	stress1.CopyTo(stress2);
	
	plasticModelCopy.ApplyLoad(stress2, strain2);
	
	plasticModelCopy.Phi(strain2,phi2);
	
#ifdef LOG4CXX_PLASTICITY
	{
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__
		<< "\nphi  after plastic loop =\n" << phi2
		<< "\ncausing strain = " << strain2
		<< "\nEnd of test";
		LOGPZ_INFO(testLogger,sout.str().c_str());
		
	}
#endif
	
}

template <class T>
inline void TPZPlasticTest::StressTest(T & plasticModel, const char * filename, REAL  stressMultiplier)
{
	
	// Load test: in this test the material is submitted to a stress loading
	// path read from a file in the format below.
	
	TPZPlasticTest::InitializeLOG();
	
#ifdef LOG4CXX_PLASTICITY
	{
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__
		<< "\nApply Load Test\nPlease enter file name containing load path (sigmaxx xy xz yy yz zz) for all iterations:\n"
		<< "\nReading file named \"" << filename << "\"\n";	
		LOGPZ_INFO(testLogger,sout.str().c_str());
	}
#endif
	
	const int linelen = 1024;
	char outfilename[linelen], line[linelen];
	
	ifstream file;
	ofstream outFile;
	
	file.open(filename);
	
	if(!file.is_open())
	{
		cout << "\nFile not open.\nExiting...\n";
		return;
	}
	
	file.getline(line, linelen);
	strncpy(outfilename, filename, 120);
	strcpy(outfilename+strlen(outfilename), ".out");
	
#ifdef LOG4CXX_PLASTICITY
	{
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__
		<< "\nReading file named \"" << filename << "\"\n"	
		<< "\nTest data description: \"" << line << "\"\n"
		<< "\nOutput file named \"" << outfilename << "\"\n";	
		LOGPZ_INFO(testLogger,sout.str().c_str());
	}
#endif
	
	outFile.open(outfilename); 
	outFile << line << endl;
	
	TPZManVector<TPZTensor<REAL>, 2> loadPath, totalStrainPath;
	int i, j=0;
	
	while(file.good())
	{
		file.getline(line, linelen);
		stringstream strLine(line);
		loadPath.Resize(j+1);
		for(i = 0; i < 6; i++)
		{
			strLine >> loadPath[j].fData[i];
		}
		
		j++;
	}
	int nsteps = loadPath.NElements();
	totalStrainPath.Resize(nsteps);
	
	for(i = 0; i < nsteps; i++)
	{
		TPZTensor<REAL> tempLoad(loadPath[i]);
		for(j=0;j<6;j++)tempLoad.fData[j]*=stressMultiplier;
		if(i > 1)totalStrainPath[i] = totalStrainPath[i-1];
		
		plasticModel.ApplyLoad(tempLoad, totalStrainPath[i]);
		
		std::stringstream outputLine;
		outputLine << "step " << i << ", sigma: " << loadPath[i] ;
		outputLine << "\nfN:/n" << plasticModel.GetState() << endl; 
		
		outFile << outputLine.str();
#ifdef LOG4CXX_PLASTICITY
		LOGPZ_INFO(testLogger,outputLine.str().c_str());
#endif
		outFile.flush();
	}
	
	file.close();
	outFile.close();
	
}


template <class T>
inline void TPZPlasticTest::StrainTest(T & plasticModel, const char * filename, REAL  strainMultiplier)
{
	
	// Strain test: in this test the material is submitted to a strain
	// path according to that read from a file in the format below.
	
	TPZPlasticTest::InitializeLOG();
	
#ifdef LOG4CXX_PLASTICITY
	{
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__
		<< "\nApply Strain Test\nPlease enter file name containing strain path (sigmaxx xy xz yy yz zz) for all iterations:\n"
		<< "\nReading file named \"" << filename << "\"\n";	
		LOGPZ_INFO(testLogger,sout.str().c_str());
	}
#endif
	
	const int linelen = 1024;
	char outfilename[linelen], line[linelen];
	
	ifstream file;
	ofstream outFile;
	
	file.open(filename);
	
	if(!file.is_open())
	{
		cout << "\nFile not open.\nExiting...\n";
		return;
	}
	
	file.getline(line, linelen);
	strncpy(outfilename, filename, 120);
	strcpy(outfilename+strlen(outfilename), ".out");
	
#ifdef LOG4CXX_PLASTICITY
	{
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__
		<< "\nReading file named \"" << filename << "\"\n"	
		<< "\nTest data description: \"" << line << "\"\n"
		<< "\nOutput file named \"" << outfilename << "\"\n";	
		LOGPZ_INFO(testLogger,sout.str().c_str());
	}
#endif
	
	outFile.open(outfilename); 
	outFile << line << endl;
	
	TPZManVector<TPZTensor<REAL>, 2> strainPath, stressPath;
	TPZFNMatrix < 6*6 > tangent (6,6,0.); // satisfying interface only
	int i, j=0;
	
	while(file.good())
	{
		file.getline(line, linelen);
		stringstream strLine(line);
		strainPath.Resize(j+1);
		for(i = 0; i < 6; i++)
		{
			strLine >> strainPath[j].fData[i];
		}
		
		j++;
	}
	int nsteps = strainPath.NElements() - 1;// the last line loaded always return null path
	stressPath.Resize(nsteps);
	
	for(i = 0; i < nsteps; i++)
	{
		TPZTensor<REAL> tempStrain(strainPath[i]);
		for(j=0;j<6;j++)tempStrain.fData[j]*=strainMultiplier;
		if(i > 1)stressPath[i] = stressPath[i-1];
		
		plasticModel.ApplyStrain(tempStrain);
		plasticModel.Sigma(tempStrain, stressPath[i], tangent);
		
		TPZTensor<REAL> epsp;
		plasticModel.GetPlasticStrain(epsp);
		REAL alpha = plasticModel.GetAlpha();
		std::stringstream outputLine;
		outputLine << "step " << i << ", imposedEps: " << strainPath[i] 
		<< ", sigma: " << stressPath[i]
		<< ", epsP= " << epsp
		<< ", alpha= " << alpha << endl; 
		
		outFile << outputLine.str();
#ifdef LOG4CXX_PLASTICITY
		LOGPZ_INFO(testLogger,outputLine.str().c_str());
#endif
		outFile.flush();
	}
	
	file.close();
	outFile.close();
	
}

///////////////////

inline void TPZPlasticTest::LoadTest(const char * filename)
{
	
	// Load test: in this test the material is submitted to a loading path
	// read from a file in the format below. The keywords strain and stress
	// switches to the proper loading conditions.
	
	TPZPlasticTest::InitializeLOG();
	
#ifdef LOG4CXX_PLASTICITY
	{
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__
		<< "\nApply Load Test\nPlease enter file name containing load path (strain/stress sigmaxx xy xz yy yz zz) for all iterations:\n"
		<< "\nReading file named \"" << filename << "\"\n";	
		LOGPZ_INFO(testLogger,sout.str().c_str());
	}
#endif
	
	const int linelen = 1024;
	char outfilename[linelen], line[linelen];
	
	ifstream file;
	ofstream outFile;
	
	file.open(filename);
	
	if(!file.is_open())
	{
		cout << "\nFile not open.\nExiting...\n";
		return;
	}
	
	//reading test description line
	file.getline(line, linelen);
	strncpy(outfilename, filename, 120);
	strcpy(outfilename+strlen(outfilename), ".out");
	
#ifdef LOG4CXX_PLASTICITY
	{
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__
		<< "\nReading file named \"" << filename << "\"\n"	
		<< "\nTest data description: \"" << line << "\"\n"
		<< "\nOutput file named \"" << outfilename << "\"\n";	
		LOGPZ_INFO(testLogger,sout.str().c_str());
	}
#endif
	
	outFile.open(outfilename); 
	outFile << line << endl;
	
	//reading second line description (CVS label data for instance)
	file.getline(line, linelen);
	outFile << line << endl;
	
	//reading plasticity model
	file.getline(line, linelen);
	TPZPlasticBase * pPlasticModel;
	
	if(!CreatePlasticModel(pPlasticModel, line))
	{
		std::stringstream sout;
		sout << "Could not create plastic model named " << line;
#ifdef LOG4CXX_PLASTICITY
		{
			LOGPZ_INFO(testLogger,sout.str().c_str());
		}
#endif	
		
		outFile << endl << sout.str();
		outFile.close();
		file.close();
		return;
	}
	
	
	TPZManVector<TPZTensor<REAL>, 2> loadPath, totalStrainPath;
	TPZTensor<REAL> stress, tempStress, strain, tempStrain, epsp, strainP;
	//TPZFNMatrix < 6*6 > tangent (6,6,0.); // satisfying interface only
	REAL stressMultiplier = 1., strainMultiplier = 1., integrTol;
	int i, j = 0;
	
	while(file.good())
	{
		file.getline(line, linelen);
		if(!strncmp(line, "strainLoad",10))
		{
			stringstream strLine(line + 10);
			for(i = 0; i < 6; i++)
			{
				strLine >> strain.fData[i];
			}
			strain.CopyTo(tempStrain);
			for(i = 0; i < 6; i++)tempStrain.fData[i]*=strainMultiplier;
			
			pPlasticModel->ApplyStrainComputeSigma(tempStrain, tempStress);
			
			tempStress.CopyTo(stress);
			for(i = 0; i < 6; i++)stress.fData[i]/=stressMultiplier;
			
			pPlasticModel->GetState().fEpsP.CopyTo(strainP);
			for(i = 0; i < 6; i++)strainP.fData[i]/=strainMultiplier;
			
			std::stringstream outputLine;
			
			outputLine << "stress step " << j 
			<< ", strain: " << strain 
			<< ", stress: " << stress
			<< ", epsP = " << strainP
			<< ", alpha = " << pPlasticModel->GetState().fAlpha
			<< ", integrationSteps = " << pPlasticModel->IntegrationSteps() << endl; 		 
			
			outFile << outputLine.str();
			outFile.flush();
#ifdef LOG4CXX_PLASTICITY
			LOGPZ_INFO(testLogger,outputLine.str().c_str());
#endif
			j++;
		}
		if(!strncmp(line, "stressLoad",10))
		{
			stringstream strLine(line + 10);
			for(i = 0; i < 6; i++)
			{
				strLine >> stress.fData[i];
			}
			stress.CopyTo(tempStress);
			for(i = 0; i < 6; i++)tempStress.fData[i]*=stressMultiplier;
			
			pPlasticModel->ApplyLoad(tempStress, tempStrain);
			
			int intSteps = pPlasticModel->IntegrationSteps();
			
			//evaluating the converged sigma
			pPlasticModel->ApplyStrainComputeSigma(tempStrain, tempStress);
			
			tempStrain.CopyTo(strain);
			for(i = 0; i < 6; i++)strain.fData[i]/=strainMultiplier;
			
			pPlasticModel->GetState().fEpsP.CopyTo(strainP);
			for(i = 0; i < 6; i++)strainP.fData[i]/=strainMultiplier;
			
			tempStress.CopyTo(stress);
			for(i = 0; i < 6; i++)stress.fData[i]/=stressMultiplier;		 
			
			std::stringstream outputLine;
			
			outputLine << "stress step " << j 
			<< ", strain: " << strain 
			<< ", stress: " << stress
			<< ", epsP = " << strainP
			<< ", alpha = " << pPlasticModel->GetState().fAlpha
			<< ", integrationSteps = " << intSteps << endl; 	
			
			outFile << outputLine.str();
			outFile.flush();
#ifdef LOG4CXX_PLASTICITY
			LOGPZ_INFO(testLogger,outputLine.str().c_str());
#endif
			j++;
		}
		if(!strncmp(line, "stressMult",10))
		{
			stringstream strLine(line + 10);
			strLine >> stressMultiplier; 
			
			std::stringstream outputLine;
			outputLine << "stressMult = " << stressMultiplier << endl;
			
			outFile << outputLine.str();
			outFile.flush();
#ifdef LOG4CXX_PLASTICITY
			LOGPZ_INFO(testLogger,outputLine.str().c_str());
#endif
		}
		if(!strncmp(line, "strainMult",10))
		{
			stringstream strLine(line + 10);
			strLine >> strainMultiplier;
			
			std::stringstream outputLine;
			outputLine << "strainMult = " << strainMultiplier << endl;
			
			outFile << outputLine.str();
			outFile.flush();
#ifdef LOG4CXX_PLASTICITY
			LOGPZ_INFO(testLogger,outputLine.str().c_str());
#endif
		}
		if(!strncmp(line, "resetMater",10))
		{
			TPZPlasticState<REAL> nullState;
			
			pPlasticModel->SetState(nullState);
			
			std::stringstream outputLine;
			outputLine << "ResetMaterial" << endl;
			
			outFile << outputLine.str();
			outFile.flush();
#ifdef LOG4CXX_PLASTICITY
			LOGPZ_INFO(testLogger,outputLine.str().c_str());
#endif
		}
		if(!strncmp(line, "integrTol ",10))
		{
			stringstream strLine(line + 10);
			strLine >> integrTol;
			pPlasticModel->SetIntegrTol(integrTol);
			
			std::stringstream outputLine;
			outputLine << "integrTol " << integrTol << endl;
			
			outFile << outputLine.str();
			outFile.flush();
#ifdef LOG4CXX_PLASTICITY
			LOGPZ_INFO(testLogger,outputLine.str().c_str());
#endif
		}
	}
	
	file.close();
	outFile.close();
	
	delete pPlasticModel;
	
}

template <class T>
inline int TPZPlasticTest::CreatePlasticModel(T * ( & pPlasticModel), const char * line)
{
	if(!strncmp(line,"LadeKim.PlainConcrete",21))
	{
		TPZLadeKim * pLK = new TPZLadeKim();
		TPZLadeKim::PlainConcrete( *pLK);
		pPlasticModel = pLK;
		return 1;
	}
	if(!strncmp(line,"LadeKim.LooseSacrRiverSand",26))
	{
		TPZLadeKim * pLK = new TPZLadeKim();
		TPZLadeKim::LooseSacrRiverSand(*pLK);
		pPlasticModel = pLK;
		return 1;
	}
	if(!strncmp(line,"LadeKim.DenseSacrRiverSand",26))
	{
		TPZLadeKim * pLK = new TPZLadeKim();
		TPZLadeKim::DenseSacrRiverSand(*pLK);
		pPlasticModel = pLK;
		return 1;
	}
	if(!strncmp(line,"LadeKim.FineSilicaSand",22))
	{
		TPZLadeKim * pLK = new TPZLadeKim();
		TPZLadeKim::FineSilicaSand(*pLK);
		pPlasticModel = pLK;
		return 1;
	}
	if(!strncmp(line,"SandlerDimaggio.McCormicRanchSandMod",36))
	{
		TPZSandlerDimaggio * pSD = new TPZSandlerDimaggio();
		TPZSandlerDimaggio::McCormicRanchSandMod(*pSD);
		pPlasticModel = pSD;
		return 1;
	}
	if(!strncmp(line,"SandlerDimaggio.McCormicRanchSand",33))
	{
		TPZSandlerDimaggio * pSD = new TPZSandlerDimaggio();
		TPZSandlerDimaggio::McCormicRanchSand(*pSD);
		pPlasticModel = pSD;
		return 1;
	}
	if(!strncmp(line,"DruckerPrager.DummyConcrete",27))
	{
		typedef TPZPlasticStep<TPZYCDruckerPrager, TPZThermoForceA, TPZElasticResponse> TPZDruckerPrager;
		TPZDruckerPrager * pDP = new TPZDruckerPrager();
		pDP->fYC.SetUp(/*phi=30*/asin(0.5),/*innerMCFit*/0);
		pDP->fTFA.SetUp(/*yield- coesao inicial correspondeno a fck igual 32 Mpa */ 9.2376, /*k Modulo de hardening da coesao equivante 1 Mpa a cada 0.1% de deformacao */1000.);
		pDP->fER.SetUp(/*young*/ 20000, /*poisson*/ 0.2);
		pPlasticModel = pDP;
		return 1;
	}
	if(!strncmp(line,"MohrCoulomb.DummyConcrete",25))
	{
		typedef TPZPlasticStep<TPZYCDruckerPrager, TPZThermoForceA, TPZElasticResponse> TPZDruckerPrager;
		TPZDruckerPrager * pDP = new TPZDruckerPrager();
		pDP->fYC.SetUp(/*phi=30*/asin(0.5),/*innerMCFit*/0);
		pDP->fTFA.SetUp(/*yield- coesao inicial correspondeno a fck igual 32 Mpa */ 9.2376, /*k Modulo de hardening da coesao equivante 1 Mpa a cada 0.1% de deformacao */1000.);
		pDP->fER.SetUp(/*young*/ 20000, /*poisson*/ 0.2);
		pPlasticModel = pDP;
		return 1;
	}
	
	
	return 0;
}

template <class T>
inline void TPZPlasticTest::GlobalCheckConv(T & plasticModel, TPZTensor<REAL> & strain, REAL maxDeltaStrain)
{
	
	// GlobalCheckConv Test: tests if the elastoplastic evaluated stiffness matrix
	// is really the jacobian of the stress tensor with respect to the total
	// strain.
	// The imported plasticModel is submitted to the imposed strain. The evaluated stress
	// and jacobian matrix are kept in memory.
	// Sucessive copies of this plastic model are subjected to further strains and
	// the check conv method is evaluated.
	
	const int nRetries = 30; // must be at least 2
	const int nVars = 6;
	TPZFNMatrix<nVars*nVars> Dep(nVars,nVars), tempMatrix(nVars,nVars);
	TPZTensor<REAL> stress, tempStress;
	//The  below matrices hold the values of the stress vector (column wise) due to the
	//strain component changes.
	TPZManVector<TPZFNMatrix<nVars*nVars>, nRetries> FcnHistory(nRetries);
	// The below strain tensors hold the changes in strain components.
	TPZManVector<TPZTensor<REAL>, nRetries> StrainHistory(nRetries);
	
	plasticModel.ApplyStrainComputeDep(strain, stress, Dep);
	
	
	//TPZPlasticTest::InitializeLOG();
	
#ifdef LOG4CXX_PLASTICITY
	{
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__
		<< "\nGlobalCheckConv Test with " << nRetries << " samplings\n";
		sout << "\n Stress = " << stress
		<< "\n Dep=\n" << Dep;
		LOGPZ_INFO(testLogger,sout.str().c_str());
	}
#endif
	
	int i,j,k;
	
	
	//Evaluating the plastic integration for all cases.
	for(k = 0; k < nRetries; k++)
	{
		FcnHistory[k].Resize(nVars,nVars);
		for(i = 0; i < nVars; i++)
		{
			T plasticModelCopy(plasticModel);
			REAL random = (rand()&1000)/1000.;
			REAL deltaStrain = random * maxDeltaStrain;
			TPZTensor<REAL> tempStrain;
			tempStrain.fData[i] = deltaStrain;
			StrainHistory[k].fData[i] = deltaStrain;
			tempStrain.Add(strain,1.);
			plasticModelCopy.ApplyStrainComputeDep(tempStrain, tempStress, tempMatrix);
			for(j = 0; j < nVars; j++)
				FcnHistory[k](j,i) = tempStress.fData[j] 
				- stress.fData[j]
				- /*Dep(j,i)*/tempMatrix(j,i) * deltaStrain;
			
			/*
			 for(j = 0; j < nVars; j++)
			 if(i==j)
			 {
			 FcnHistory[k](j,i) = pow(-strain.fData[j],1.5) +
			 1.5 * pow(-strain.fData[j],0.5) * (-deltaStrain)
			 - pow(-tempStrain.fData[j],1.5);
			 }else{
			 FcnHistory[k](j,i) = pow(-strain.fData[j],1.5) +
			 1.5 * pow(-strain.fData[j],0.5) * 0.
			 - pow(-tempStrain.fData[j],1.5);
			 }
			 */
			
#ifdef LOG4CXX_PLASTICITY
			{
				std::stringstream sout;
				sout << "\nCase " << k
				<< "\nStrain= " << tempStrain
				<< "\nStress= " << tempStress
				<< "\ntempMatrix= " << tempMatrix;
				LOGPZ_INFO(testLogger,sout.str().c_str());
			}
#endif
		}
		
	}
	
	// Evaluating the checkConv routines
	
	std::stringstream output;
	for(k = 1; k < nRetries; k++)
	{
		output << "\n\n Result of checkconv with steps " << k-1 << " and " << k << "\n";
		output << "\nStrainHistory[" << k-1 << "]=" << StrainHistory[k-1];
		output << "\nFcnHistory[" << k-1 << "]=" << FcnHistory[k-1];
		output << "\nStrainHistory[" << k << "]=" << StrainHistory[k];
		output << "\nFcnHistory[" << k << "]=" << FcnHistory[k];
		output << "\nCheckConv:\n" ;
		for(j = 0; j < nVars; j++)
		{
			for(i = 0; i < nVars; i++)
			{
				REAL ykm1 = FcnHistory[k-1](j,i);
				REAL yk   = FcnHistory[k  ](j,i);
				REAL xkm1 = StrainHistory[k-1].fData[i];
				REAL xk   = StrainHistory[k  ].fData[i];
				
				output.width(12);
				if(ykm1 * yk < 1.e-24)
				{
					if(ykm1 * yk < -1.e-24)
					{
						output << "+/-y?"; // They should be of the same sign
					}else
					{
					    output << "Exact"; // at least one of the evaluations was exact
					}
					
				}else
				{
					if(xkm1 * xk < 1.e-24)
					{
						output << "+/-x?";
					}else
					{
						output << log(ykm1/yk) / log(xkm1 / xk);
					}
				}
			}
			output << "\n";
		}
		output << "\n" ;
	}
	
#ifdef LOG4CXX_PLASTICITY
	{
//		LOGPZ_INFO(testLogger,output.str().c_str());
      LOGPZ_INFO(MaterialPoint,output.str().c_str());
	}
#endif
	
}

inline void TPZPlasticTest::DruckerPragerTest()
{
	
	TPZPlasticTest::InitializeLOG();
	ofstream outfile("comparamathematica1.nb"); 
	
	int choice;
	choice =0;
    cout << "\nChoose Load Case:";
    cout << "\n0 - Strain Step";
    cout << "\n1 - Stress Step";
	cout << "\n";
	
	//cin >> choice;
	TPZTensor<REAL> stress, strain, deltastress, deltastrain;
	
	//deltastress.XX() = -1.;
	//	deltastress.YY() = 0.;
	//	deltastress.ZZ()=  0.;
	//	deltastress.XY()=  0.;
	//	deltastress.XZ()=  0.;
	//	deltastress.YZ()=  0.;
	
    deltastress.XX() = -10.;
	deltastress.ZZ() =  1.;
	stress = deltastress;
	
	/*
	 deltastrain.XX() = -0.0001;
	 //	deltastrain.XY() = 0.0001;
	 //	deltastrain.XZ() = -0.0001;
	 //	deltastrain.YZ() = -0.0001;
	 deltastrain.ZZ() = 0.00002;
	 strain=deltastrain;
	 */	
	
	typedef TPZPlasticStep<TPZYCDruckerPrager, TPZThermoForceA, TPZElasticResponse> TPZDruckerPrager;
	TPZDruckerPrager Pstep;
	REAL pi = M_PI;
	/*innerMCFit = 0*/
	/*OuterMCFit = 1*/
	Pstep.fYC.SetUp(/*phi=20*/ 20./180. * pi ,/*MCFit*/1);
	REAL coesao = 9.2376;
	Pstep.fTFA.SetUp(/*yield- coesao inicial correspondeno a fck igual 32 Mpa */ coesao, /*k Modulo de hardening da coesao equivante 1 Mpa a cada 0.1% de deformacao */1000.);
	Pstep.fER.SetUp(/*young*/ 20000., /*poisson*/ 0.2);
	
	int length = 300;
	for(int step=0;step<length;step++)
	{
		cout << "\nstep "<< step;
		
		
		if( step==100 || step==200)
		{
			deltastrain *= -1.;
			deltastress *= -1.;
		}
		
		stress += deltastress;
		Pstep.ApplyLoad(stress,strain);
		
		
		//	strain += deltastrain;
		//	Pstep.ApplyStrainComputeSigma(strain,stress);
		
		if(step==0)outfile <<"points={";
		
		if(step!=length && step!=length-1)
		{
			outfile <<"{"<<-strain.XX()<< ","<<-stress.XX()<<" }, ";
			//outfile <<"{"<<-stress.J2()<< ","<<-stress.I1()<<" }, ";
			//outfile <<"{"<<stress.ZZ()<< ","<<stress.XX()<<" }, ";
		}
		
		if(step==length-2)outfile <<"{"<<-strain.XX()<< ","<<-stress.XX()<<"}";
		//if(step==length-2)outfile <<"{"<<-stress.J2()<< ","<<-stress.I1()<<"}";
		//if(step==length-2)outfile <<"{"<<stress.ZZ()<< ","<<stress.XX()<<"}";
		
		if(step==length-1)
		{
			//do = SolAprox = ListPlot[points, PlotRange -> All, AxesLabel -> {epsilonx, sigmax}, PlotLabel -> DRUCKER*PRAGER*OUTERMCFIT, AxesOrigin -> {0, 0},PlotStyle -> {Thick, Blue, Dashed}, Joined -> True, FillingStyle -> Directive[Opacity[0.5], Orange]]
			outfile << " };\n";
			outfile<<"SolAprox = ListPlot[points, PlotRange-> All,AxesLabel->{epsilon,sigma}, PlotLabel->DRUCKER PRAGER,AxesOrigin->{0,0},Filling -> Axis , FillingStyle -> Directive[Opacity[0.5], Orange]]";
		}
		
	}
	
}

inline void TPZPlasticTest::MohrCoulombTest()
{
	int choice =1;
	//TPZPlasticTest::InitializeLOG();
	ofstream outfiletxt("mohrcoulomb.txt"); 
	
	TPZTensor<REAL> stress, strain, deltastress, deltastrain;
	
		deltastress.XX() = -13.5;
		deltastress.XY() = -0.01;
		deltastress.XZ() = 0.;
		deltastress.YY() = 0.;
		deltastress.YZ() = 0.;
		deltastress.ZZ() =  -0.01;
		stress = deltastress;
	
//	TPZFNMatrix<6*6> Dep(6,6,0.);
//    deltastrain.XX() = -0.0002;
//	deltastrain.XY() = 0.;
//	deltastrain.XZ() = 0.;
//	deltastrain.YY() = -0.00000001;
//	deltastrain.YZ() = 0.;
//	deltastrain.ZZ() = -0.00000001;
//	strain=deltastrain;
	
	typedef TPZPlasticStep<TPZYCModifiedMohrCoulomb, TPZThermoForceA, TPZElasticResponse> TPZMohrCoulomb;
	TPZMohrCoulomb Pstep;
	Pstep.fYC.SetUp(/*phi=20*/ 20./180. * M_PI,1);
	REAL coesao = 9.2376;
	Pstep.fTFA.SetUp(/*yield- coesao inicial correspondeno a fck igual 32 Mpa */ coesao, /*k Modulo de hardening da coesao equivante 1 Mpa a cada 0.1% de deformacao */1000.);
	Pstep.fER.SetUp(/*young*/ 20000., /*poisson*/ 0.);
	
	int length =1;
	for(int step=0;step<length;step++)
	{
		cout << "\nstep "<< step;
		
//		
//		if(step == 30 || step== 30 || step == 60)
//		{
//			deltastrain *= -1.;
//			deltastress *= -1.;
//		}
		
		stress += deltastress;
	//	strain+=deltastrain;
		Pstep.ApplyLoad(stress,strain);
	//	Pstep.ApplyStrainComputeDep(strain,stress, Dep);
		//TPZTensor<REAL> eigenval,dSigma1,dSigma2,dSigma3;
		//stress.Eigenvalue(eigenval,dSigma1,dSigma2,dSigma3);
		cout<<  "\nstress " << stress << endl;
		cout<<  "\nstrain " << strain << endl;
		Pstep.Print(cout);
		TPZVec<REAL> phis(1);
		Pstep.Phi(strain, phis);
		cout << "\nphis " << phis << endl;
 		//cout<<  "\nEigen " << eigenval << endl;

		outfiletxt << fabs(strain.XX()) << " " << fabs(stress.XX()) << "\n"; 
		
	}
	
	
	
}


inline void TPZPlasticTest::ModifiedMohrCoulombTest()
{
	
	TPZPlasticTest::InitializeLOG();
	ofstream outfile("mohr"); 
	
	int choice;
	choice =1;
    cout << "\nChoose Load Case:";
    cout << "\n0 - Strain Step";
    cout << "\n1 - Stress Step";
	cout << "\n";
	
	//cin >> choice;
	TPZTensor<REAL> stress, strain, deltastress, deltastrain;
	
	deltastrain.XX() = -0.000005;
	deltastrain.YY() = -0.000005;
	deltastress.XX()= -10.;
	deltastress.YY()= -20.;
	
	typedef TPZPlasticStep<TPZYCModifiedMohrCoulomb, TPZThermoForceA, TPZElasticResponse> TPZModifiedMohrCoulomb;
	TPZModifiedMohrCoulomb Pstep;
	REAL pi = M_PI;
	Pstep.fYC.SetUp(/*phi=20*/ 20./180. * pi ,/*innerMCFit*/0);
	REAL coesao = 9.2376;
	Pstep.fTFA.SetUp(/*yield- coesao inicial correspondeno a fck igual 32 Mpa */ coesao, /*k Modulo de hardening da coesao equivante 1 Mpa a cada 0.1% de deformacao */1000.);
	Pstep.fER.SetUp(/*young*/ 20000., /*poisson*/ 0.);
	
	
	int length = 1;
	for(int step=0;step<length;step++)
	{
		
		if(choice == 0)
		{
			strain += deltastrain;
			Pstep.ApplyStrainComputeSigma(strain,stress);
		}
		
		if(choice == 1)
		{
			stress += deltastress;
			Pstep.ApplyLoad(stress, strain);
		}
		
#ifdef LOG4CXX
		{
			std::stringstream sout;
			sout << "\nMOHRCOULOMB\n";
			sout<<"\n";	
			sout << "STRAIN" <<std::endl;
			sout<<strain<<std::endl;
			sout << "STRESS"<<std::endl;
			sout<<stress<<std::endl;
			
			sout<<"STRESS-I1"<<std::endl;
			sout<<stress.I1()<<std::endl;
			sout<<"STRESS-I2"<<std::endl;
			sout<<stress.I2()<<std::endl;
			sout<<"STRESS-I3"<<std::endl;
			sout<<stress.I3()<<std::endl;
			
			sout<<"STRESS-J2"<<std::endl;
			sout<<stress.J2()<<std::endl;
			sout<<"STRESS-J3"<<std::endl;
			sout<<stress.J3()<<std::endl;
			
			
			//		LOGPZ_INFO(MaterialPoint,sout.str());
		}
#endif
		
	}
}

inline void TPZPlasticTest::WillamWarnkeTest()
{
	int choice =1;
	//TPZPlasticTest::InitializeLOG();
	ofstream outfiletxt("TPZYCWillamWarnke.txt"); 
	
	TPZTensor<REAL> stress, strain, deltastress, deltastrain;
	TPZFNMatrix<6*6> Dep(6,6,0.);

	deltastress.XX() = -1.;
	deltastress.YY()= 0.;
	deltastress.ZZ()= 0.;
	deltastress.XY() = 0.;
	deltastress.XZ() =  0.;
	deltastress.YZ()=  0.;
	stress = deltastress;
	
//	 deltastrain.XX() = -0.0001;
//	 deltastrain.XY() = 0.0000001;
//	 deltastrain.XZ() = -0.0000001;
//	 deltastrain.YY() = -0.0000001;
//	 deltastrain.YZ() = -0.0000001;
//	 deltastrain.ZZ() = 0.0000002;
//	 strain=deltastrain;
	 
	typedef TPZPlasticStep<TPZYCWillamWarnke, TPZThermoForceA, TPZElasticResponse> TPZWillamWarnke;
	TPZWillamWarnke WW;
	WW.fYC.SetUp(1.,1.,20.); 
	REAL coesao = 9.2376;
	WW.fTFA.SetUp(/*yield- coesao inicial correspondeno a fck igual 32 Mpa */ coesao, /*k Modulo de hardening da coesao equivante 1 Mpa a cada 0.1% de deformacao */1000.);
	WW.fER.SetUp(/*young*/ 20000., /*poisson*/ 0.2);
	
	
	int length = 40;
	for(int step=0;step<length;step++)
	{
		cout << "\nstep "<< step;
		
		
		//		if(step == 50 || step==100 || step == 150)
		//		{
		//			deltastrain *= -1.;
		//			deltastress *= -1.;
		//		}
		
		stress += deltastress;
//		strain+=deltastrain;
		WW.ApplyLoad(stress,strain);
#ifdef LOG4CXX
		{
			std::stringstream sout;
			sout << "\n *********************************** STEP LOAD IN WILLANWERNAKE  = "<< step << "*********************************** \n";
			LOGPZ_INFO(plasticIntegrLogger,sout.str());
		}
#endif
//		WW.ApplyStrainComputeDep(strain,stress,Dep);
		TPZTensor<REAL> eigenval,dSigma1,dSigma2,dSigma3;
		stress.Eigenvalue(eigenval,dSigma1,dSigma2,dSigma3);
		cout<<  "\nstress " << stress << endl;
		cout<<  "\nEigen " << eigenval << endl;
		
		outfiletxt << strain.XX() << " " << stress.XX() << "\n"; 
		
	}

	

	
	
	
}



#include <math.h>
inline void RotationMatrix(TPZFMatrix &R, double thetaRad, int axis)
{
	R.Resize(3,3);
	
	
	switch (axis) 
	{
			
		case 0://ROTATE ABOUT X
			
			R.Put(0,0,1.);
			R.Put(1,1,cos(thetaRad));R.Put(1,2,sin(thetaRad));
			R.Put(2,1,-sin(thetaRad));R.Put(2,2,cos(thetaRad));
			
			break;
			
		case 1://ROTATE ABOUT Y
			
			R.Put(1,1,1.);
			R.Put(0,0,cos(thetaRad));R.Put(0,2,sin(thetaRad));
			R.Put(2,0,-sin(thetaRad));R.Put(2,2,cos(thetaRad));
			
			break;
			
		case 2://ROTATE ABOUT Z
			
			R.Put(0,0,cos(thetaRad)); R.Put(0, 1,sin(thetaRad));
			R.Put(1,0,-sin(thetaRad)); R.Put(1, 1,cos(thetaRad));
			R.Put(2,2,1.);
			
			break;
			
		case 3:
			
			
			R.Put(0,0,(1./3.)*(1+2*cos(thetaRad)));R.Put( 0,1, (1./3.) * (1-cos(thetaRad) - sqrt(3)*sin(thetaRad)) );R.Put( 0,2, (1./3.) * (1-cos(thetaRad) + sqrt(3)*sin(thetaRad)) );
			R.Put( 1,0, (1./3.) * (1-cos(thetaRad) + sqrt(3)*sin(thetaRad)) );R.Put(1,1,(1./3.)*(1+2*cos(thetaRad)));R.Put( 1,2, (1./3.) * (1-cos(thetaRad) - sqrt(3)*sin(thetaRad) ));
			R.Put( 2,0, (1./3.) * ( 1-cos(thetaRad) - sqrt(3)*sin(thetaRad) ) );R.Put( 2,1, (1./3.) * ( 1-cos(thetaRad) + sqrt(3)*sin(thetaRad) ) );R.Put(2,2,(1./3.)*(1+2*cos(thetaRad)));
			
			break;
			
		default:
			
			std::cout << " NON SPECIFIED AXIS "<<std::endl;
			break;
	}
	
}



//Metodo que aplica a tensao e calcula a deformacao de um estado de tensoes um pouco fora da zona elastica.
//Este metodo tambem rotaciona este estado de tensoes de forma a cobrir toda a superficie do plano em questao oriunda das constantes materiais escolhida pelo usuario.
//A variacao na direcao do eixo hidrostatico pode ser modificada alterando as constantes materiais no caso de materiais sensiveis a pressao hidrostatica.
//Neste metodo tambem apos a aplicacao da tensao ligeiramente fora da zona elastica, e obtida a correspondente deformacao total, e chamado o metodo GlobalCheckConv,
//que verificase a matriz de rigidez elastoplastica e realmente a jacobiana dSigma/dEpsT. Neste metodo tambem o modelo plastico importado e submetido a deformacao imposta. A tensao 
//e a matriz jacobiana calculadas sao mantidas na memoria. Sucessivas copias deste modelo de plastico sao submetidos a mais deformacao e a verificação de convergencia do método é avaliada.
template <class T>
inline void TPZPlasticTest::PlasticIntegratorCheck(int thetaintervals, T mat)
{
	
	std::ifstream input("SnubDodecahedron.txt");
	int sizedirs;
	input >>sizedirs; 
	TPZFMatrix directions(sizedirs,3,0.);
	
	TPZTensor<REAL> DiagonalStress,epst,epsp;
	
	epst.fData[_XX_] =  0.0000001;
	epst.fData[_XY_] =  0.0000001;
	epst.fData[_XZ_] =  0.0000001;
	epst.fData[_YY_] =  0.0000001;
	epst.fData[_YZ_] =  0.0000001;
	epst.fData[_ZZ_] =  0.0000001;
	
	epsp.fData[_XX_] =  0.;
	epsp.fData[_XY_] =  0.;
	epsp.fData[_XZ_] =  0.;
	epsp.fData[_YY_] =  0.;
	epsp.fData[_YZ_] =  0.;
	epsp.fData[_ZZ_] =  0.;
	
	int nyield = mat.NYield;
	
	TPZVec<REAL>  funcs(nyield);
	int checkForcedYield;
	
	
	
	
	for(int i=0;i<sizedirs;i++)
	{
/*
			REAL coordxx,coordyy,coordzz;
			input >> coordxx >> coordyy >> coordzz;
			DiagonalStress.XX()=coordxx;
			DiagonalStress.YY()=coordyy;
			DiagonalStress.ZZ()=coordzz;
			mat.ApplyLoad(DiagonalStress,epst);
			TPZPlasticState<REAL> stateone = mat.GetState();
		  //  stateone.EpsT() = epst;
 */
			REAL coordxx,coordyy,coordzz;
			input >> coordxx >> coordyy >> coordzz;
		
			DiagonalStress.XX()=coordxx;
			DiagonalStress.YY()=coordyy;
			DiagonalStress.ZZ()=coordzz;
			DiagonalStress*=0.5;
			cout << " i "<< i <<endl;
			T plasticModelCopy(mat);
#ifdef LOG4CXX
		{
			std::stringstream sout;
			sout << " \n\n i = " << i << endl;
			//LOGPZ_INFO(MaterialPoint,sout.str());
			
		}
#endif
		REAL func;
		do{
			
			//CONTINUA CONSTANTE COM 1.1 Comentado;
			//PHIS MAIORES DO QUE ZERO SEMPRE com 1.1 nao comentado
			plasticModelCopy.fYC.Compute(DiagonalStress,plasticModelCopy.fTFA.Compute( plasticModelCopy.GetState().Alpha() ),funcs,checkForcedYield);
			plasticModelCopy.ApplyLoad(DiagonalStress,epst);
			TPZFNMatrix<6*6> Dep(6,6,0.);
			plasticModelCopy.ApplyStrainComputeDep(plasticModelCopy.GetState().EpsT(),DiagonalStress,Dep);
			DiagonalStress*=1.1;

			
		
#ifdef LOG4CXX
			{
				
				std::stringstream sout;
				sout << " \n Dep = " << Dep << endl;
				sout << " \n funcs[0] = " << funcs[0] << endl;
				sout << " \n DiagonalStress = " << DiagonalStress <<endl;
				sout << " \n fTFA = " << plasticModelCopy.fTFA.Compute( plasticModelCopy.GetState().Alpha() ) <<endl;
				sout << " \n Alpha() = " << plasticModelCopy.GetState().Alpha() << endl;
				sout << " \n epst = " << plasticModelCopy.GetState().EpsT() <<endl;
				sout << " \n epsP = " << plasticModelCopy.GetState().EpsP() <<endl;
			//	sout << "\n Dep  = " << Dep << endl;
				//LOGPZ_INFO(MaterialPoint,sout.str());
				
			}
#endif
			
			
		/*	
			//CONTINUA CONSTANTE COM 1.1 Comentado;
			mat.ApplyLoad(DiagonalStress,epst);
			cout << "DiagonalStress = " << DiagonalStress <<endl;
			mat.fYC.Compute(DiagonalStress,mat.fTFA.Compute( mat.GetState().Alpha() ),funcs,checkForcedYield);
		    DiagonalStress*=1.1;
		    cout << "\nPHI = "<< funcs[0] <<endl;
		 
		 
		 
		 */	
			if(funcs[0] < -100.)break;

		   }while( 0.> funcs[0]);

		TPZFNMatrix<6> input(6,1), Range(6,1);
		input(_XX_) = epst.XX();
		input(_YY_) = epst.XY();
		input(_ZZ_) = epst.XZ();
		input(_XY_) = epst.YY();
		input(_XZ_) = epst.YZ();
		input(_YZ_) = epst.ZZ();
		Range = input * (1./19.);
		TPZVec< REAL > Coefs(1,1.);
		TPZPlasticTest test;
		CheckConvergence(test, input, Range, Coefs);
		
		epst = plasticModelCopy.GetState().EpsT();
	
	//	GlobalCheckConv(plasticModelCopy, epst, 0.0001);
		
		
		/*
			mat.fYC.Compute(DiagonalStress,mat.fTFA.Compute( mat.GetState().Alpha() ),funcs,checkForcedYield);
		    cout << funcs << endl; 
			REAL deltaeps = 0.0001;
			GlobalCheckConv(mat,epst, deltaeps);
			cout << DiagonalStress << endl; 
		 */
		
	}
	
/*	
	
	
	TPZTensor<REAL> stress, strain, deltastress, deltastrain,straincompress;
	
	TPZTensor<REAL> sigma;
	REAL yieldfunc;
	
	mat.SetUp();
	mat.SetIntegrTol(0.0001);
	mat.Print(cout);
		
	int nyield = mat.NYield;
	
	TPZVec<REAL>  funcs(nyield);
	int checkForcedYield;
	TPZFMatrix sigma1(3,1,0.);//Cria um vetor no espaco das tensoes principais que quase ultrapassa a superficie de plastificacao
	
	REAL thetatotal = 2.* M_PI;
	REAL theta = 0.;
	int rotateaboutaxes = 3;//Rotacao em torno do eixo deviatorico {1,1,1}
	TPZTensor<REAL> Tensor,epst,epsp,Tensor2;
	
	TPZPlasticState<REAL> state,state2;
	
	epst.fData[_XX_] =  0.0000001;
	epst.fData[_XY_] =  0.0000001;
	epst.fData[_XZ_] =  0.0000001;
	epst.fData[_YY_] =  0.0000001;
	epst.fData[_YZ_] =  0.0000001;
	epst.fData[_ZZ_] =  0.0000001;
	
	epsp.fData[_XX_] =  0.;
	epsp.fData[_XY_] =  0.;
	epsp.fData[_XZ_] =  0.;
	epsp.fData[_YY_] =  0.;
	epsp.fData[_YZ_] =  0.;
	epsp.fData[_ZZ_] =  0.;

	
	
	
	state.fEpsT = epst;
	state.fEpsP = epsp;
	state = mat.GetState();
	REAL yieldradius = mat.YieldRadius(state);
	sigma1(0,0) = yieldradius*1.2;

		for(int i=0;i<thetaintervals;i++)
		{
			
			TPZFMatrix R,resp;
			TPZFNMatrix<6> Dep(6,6,0.);
			RotationMatrix(R, theta,rotateaboutaxes);
			R.Multiply(sigma1,resp);
			Tensor.XX() = resp(0,0);
			Tensor.YY() = resp(1,0);
			Tensor.ZZ() = resp(2,0);
			theta += thetatotal/thetaintervals;
			mat.ApplyLoad(Tensor,epst);
			mat.ApplyStrainComputeDep(epst,Tensor2,Dep);
			state2 = mat.GetState();
			mat.fYC.Compute(Tensor,mat.fTFA.Compute( state2.Alpha() ),funcs,checkForcedYield);
			REAL deltaeps = 0.0001;
			GlobalCheckConv(mat,epst, deltaeps);
			
#ifdef LOG4CXX
			{
				
				std::stringstream sout;
				sout << " \n looptheta = " << i << endl;
				sout << " \n theta = " << theta <<endl;
				sout << " \n Tensor = " << Tensor <<endl;
				sout << "\n Tensor2 = " << Tensor2 << endl;
				sout << " \n epsT = " << state2.fEpsT <<endl;
				sout << " \n epsP = " << state2.fEpsP <<endl;
				sout << " \n phi = " << funcs <<endl;
				sout << " \n Alpha = " << state2.Alpha() <<endl;
				sout << "\n Dep  = " << Dep << endl;
				LOGPZ_INFO(MaterialPoint,sout.str());
				
			}
#endif
			
		}
*/
	
}

#include "TPZDruckerPrager.h"

inline void TPZPlasticTest::DruckerTest()
{
	
	ofstream outfiletxt1("e1LK.txt");
	ofstream outfiletxt2("e2LK.txt");
	ofstream outfiletxt3("e3LK.txt"); 
    ofstream outfiletxt4("VolLK.txt"); 	
	TPZTensor<REAL> stress, strain, deltastress, deltastrain;
//
//	deltastress.XX() = -0.5;
//	deltastress.XY() = -0.001;
//	deltastress.XZ() = -0.001;
//	deltastress.YY() = -0.001;
//	deltastress.YZ() = -0.001;
//	deltastress.ZZ() = -0.001;
	
	deltastress.XX() = -0.1;
	deltastress.XY() = -0.001;
	deltastress.XZ() = -0.003;
	deltastress.YY() = -0.15;
	deltastress.YZ() = -0.0015;
	deltastress.ZZ() = -0.17;
	
	
	
//	deltastrain.XX() = -0.00001;
//	deltastrain.XY() = -0.0000001;
//	deltastrain.XZ() = -0.0000003;
//	deltastrain.YY() = -0.0000015;
//	deltastrain.YZ() = -0.00000015;
//	deltastrain.ZZ() = -0.0000017;
//	strain=deltastrain;
	
//	TPZSandlerDimaggio SD;
//	TPZSandlerDimaggio::McCormicRanchSand(SD);
	
	TPZLadeKim LK;
	//TPZLadeKim::FineSilicaSand(LK);
	//TPZLadeKim::DenseSacrRiverSand(LK);
	TPZLadeKim::PlainConcrete(LK);
	
	typedef TPZPlasticStep<TPZYCDruckerPrager, TPZThermoForceA, TPZElasticResponse> TPZDruckerPrager;
	TPZDruckerPrager DP;
	REAL pi = M_PI;
	/*innerMCFit = 0*/
	/*OuterMCFit = 1*/
	DP.fYC.SetUp(/*phi=20*/ 20./180. * pi ,/*MCFit*/0);
	REAL coesao = 9.2376;
	DP.fTFA.SetUp(/*yield- coesao inicial correspondeno a fck igual 32 Mpa */ coesao, /*k Modulo de hardening da coesao equivante 1 Mpa a cada 0.1% de deformacao */1000.);
	DP.fER.SetUp(/*young*/ 20000., /*poisson*/ 0.2);
	
//	LK.Print(cout);
	LK.ApplyLoad(stress, strain);
	
//	deltastress.XX() = -0.8;
//	deltastress.XY() = -0.001;
//	deltastress.XZ() = -0.001;
//	deltastress.YY() = -0.001;
//	deltastress.YZ() = -0.001;
//	deltastress.ZZ() = -0.001;
	deltastress.XX() = -1000.;
	deltastress.XY() = -0.001;
	deltastress.XZ() = -0.003;
	deltastress.YY() = -0.15;
	deltastress.YZ() = -0.0015;
	deltastress.ZZ() = -0.17;
	stress = deltastress;
	
	int length =30;
	for(int step=0;step<length;step++)
	{
		cout << "\nstep "<< step;	
		LK.ApplyLoad(stress,strain);
		REAL pa = stress.I1()/3.;
		outfiletxt1 << strain.XX() << " " << fabs(stress.XX()) << "\n";
		outfiletxt2 << strain.YY() << " " << fabs(stress.XX()) << "\n";
		outfiletxt3 << strain.ZZ() << " " << fabs(stress.XX()) << "\n"; 
		outfiletxt4 << strain.I1() << " " << fabs(stress.XX()) << "\n"; 
		stress += deltastress;
		cout << "strain = "<<strain <<"\n";
		cout << "sigma = "<< stress <<"\t "<< "I1 = "<< stress.I1() <<"\n";
 
	}

	
}


inline void TPZPlasticTest::VonMisesTest()
{	
	
	ofstream outfiletxt1("Mises25.txt");
	ofstream outfiletxt2("Mises50.txt");
	TPZTensor<REAL> stress, strain, deltastress, deltastrain,stress2,strain2;
	
		deltastress.XX() = 250.;
		deltastress.XY() = 0.0001;
		deltastress.XY() = 0.0001;
		deltastress.YY() =  0.0001;
		deltastress.YZ() =  0.0001;
		deltastress.ZZ() =  0.0001;
		stress = deltastress;
	    stress2= deltastress;

	
	TPZVonMises VM1;/*CA25*/
	VM1.fTFA.SetUp(250.,2500./*(300/0.12)*/);
	VM1.fER.SetUp(/*young*/ 210000., /*poisson*/ 0.3);
	
	TPZVonMises VM2;/*CA50*/
	VM2.fTFA.SetUp(500.,7812.5/*(625/0.08)*/);
	VM2.fER.SetUp(/*young*/ 210000., /*poisson*/ 0.3);
	
	int length = 2;
	for(int step=0;step<length;step++)
	{
		cout << "\nstep "<< step;
	
		if(step==350|| step ==700)
		{
			deltastrain *= -1.;
			deltastress *= -1.;	
		}
		
		stress += deltastress;
		VM1.ApplyLoad(stress,strain);

		outfiletxt1 << strain.XX() << " " << stress.XX() << "\n"; 

		
	}
/*	
	int length2 = 2000;
	for(int step=0;step<length2;step++)
	{
		cout << "\nstep "<< step;
		
		if(step==600|| step ==1200)
		{
			deltastrain *= -1.;
			deltastress *= -1.;	
		}
		
		stress2 += deltastress;
		VM2.ApplyLoad(stress2,strain2);
		outfiletxt2 << strain2.XX() << " " << stress2.XX() << "\n";
		
		
	}
	*/
}



inline void TPZPlasticTest::UndocumentedTest2()
{
	TPZTensor<REAL> epsT, epsTNp1, deltaEpsT;
	TPZTensor<REAL> epsP;
	REAL alpha;
	
	epsT.fData[_XX_] = -0.000279312;
	epsT.fData[_XY_] = -0.0033376;
	epsT.fData[_XZ_] =  3.43561e-10;
	epsT.fData[_YY_] =  0.0180517;
	epsT.fData[_YZ_] =  3.96378e-10;
	epsT.fData[_ZZ_] =  0.0318554;
	
	epsP.fData[_XX_] = -0.00170072;
	epsP.fData[_XY_] = -0.00130612;
	epsP.fData[_XZ_] =  1.4072e-10;
	epsP.fData[_YY_] =  0.0140323;
	epsP.fData[_YZ_] =  1.48982e-10;
	epsP.fData[_ZZ_] =  0.028329;
	
	alpha = 50.4337;
	
	epsTNp1.fData[_XX_] = -0.00142113;
	epsTNp1.fData[_XY_] = -0.00528525;
	epsTNp1.fData[_XZ_] =  9.33585e-10;
	epsTNp1.fData[_YY_] =  0.0184717;
	epsTNp1.fData[_YZ_] =  9.64833e-10;
	epsTNp1.fData[_ZZ_] =  0.0318554;
	
	TPZPlasticState<REAL> state, state2;
	
	state.fEpsT = epsT;
	state.fEpsP = epsP;
	state.fAlpha = alpha;
	
	TPZLadeKim LK;
	TPZLadeKim::FineSilicaSand(LK);
	TPZLadeKim::DenseSacrRiverSand(LK);
	LK.SetState(state);
	
	LK.SetIntegrTol(0.0001);
	
	//LK.Print(cout);
	
	TPZTensor<REAL> sigma;
	TPZFNMatrix<6*6> Dep(6,6,0.);
	
	LK.ApplyStrainComputeDep(epsTNp1, sigma,Dep);
	//cout << " \n DEP \n";
	//cout << Dep << endl;
	
	//	LK.ApplyStrain(epsTNp1);
	//	LK.Print(cout);
}

//Verifica a se a derivada N = dphi/dsigma e consistente
inline void TPZPlasticTest::UndocumentedTest4()
{
	typedef TFad<6, REAL> TFAD;
	
	REAL A;
	TFAD A_FAD;
	
	TPZTensor<REAL> sigma;
	sigma.XX() = -0.1;
	sigma.XY() = -0.001;
	sigma.XZ() = -0.003;
	sigma.YY() = -0.15;
	sigma.YZ() = -0.0015;
	sigma.ZZ() = -0.17;
	
	TPZTensor<TFAD> sigma_FAD;
	sigma.CopyTo(sigma_FAD);
	
	int i;
	for(i = 0; i < 6; i++)sigma_FAD.fData[i].diff(i,6);
	
	TPZYCSandlerDimaggio YCSD;
	TPZYCSandlerDimaggio::McCormicRanchSand(YCSD);
	
	TPZVec<TFAD> res_FAD(2);
	TPZVec<TPZTensor<REAL> > NDir(2);
	
	YCSD.N(sigma, A, NDir,0);
	YCSD.Compute(sigma_FAD, A_FAD, res_FAD,0);
	
	cout << "\nNDir = " << NDir;
	cout << "\nres_FAD = " << res_FAD;
	
	
}

inline void TPZPlasticTest::UndocumentedTest3()
{
	typedef TFad<6, REAL> TFAD;
	
	//REAL A;
	TFAD A_FAD;
	
	TPZTensor<REAL> sigma, strain;
	strain.XX() = -0.1;
	strain.XY() = -0.001;
	strain.XZ() = -0.003;
	strain.YY() = -0.15;
	strain.YZ() = -0.0015;
	strain.ZZ() = -0.17;
	
	strain.XX() = -1;
	strain.XY() =  0;
	strain.XZ() =  0;
	strain.YY() =  0.5;
	strain.YZ() =  0;
	strain.ZZ() =  0.5;
	
	REAL multipl;
	
	cin >> multipl;
	
	strain *= multipl;
	
	TPZFNMatrix<6*6> Dep(6,6,0.);
	
	TPZSandlerDimaggio SD;
	TPZSandlerDimaggio::McCormicRanchSand(SD);
	SD.SetIntegrTol(1.e-0);
	
	SD.ApplyStrainComputeDep(strain,sigma,Dep);
	
	cout << "\n sigma: " << sigma;
	
	cout << "\n Dep: " << Dep;
	
	cout << "\n PlasticState: " << SD.GetState();
	
	strain *= 2.;
	
	SD.ApplyStrainComputeDep(strain,sigma,Dep);
	
	cout << "\n sigma: " << sigma;
	
	cout << "\n Dep: " << Dep;
	
	cout << "\n PlasticState: " << SD.GetState();
	
	strain.XX() -= 0.00001;
	
	SD.ApplyStrainComputeDep(strain,sigma,Dep);
	
	cout << "\n sigma: " << sigma;
	
	cout << "\n Dep: " << Dep;
	
	cout << "\n PlasticState: " << SD.GetState();
	
}
