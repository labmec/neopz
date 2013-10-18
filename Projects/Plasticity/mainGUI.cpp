#include <iostream>
#include <cstdlib>


#include "mainGUI.h"

/// Default constructor
TPZPlasticityTest::TPZPlasticityTest() : fZStressKnown(false), fRStressKnown(false), fPoreStressRZ(2,0.), fPoreStrainRZ(2,0.), fStressRZInput(0,2,0.),
                fStrainRZInput(0,2,0.), fStressRZSimulated(0,2,0.), fStrainRZSimulated(0,2,0.)
{
    fNumSteps = 200;
}

/// Get the stress to the pore stress
void TPZPlasticityTest::ApplyInitialStress()
{
    // load hydrostatic till I1 is in the middle of the ellips
    STATE I1 = fPoreStressRZ[0]*2.+fPoreStressRZ[1];
    STATE I1Initial = I1-fSandler.fYC.fA*fSandler.fYC.fR;
    
cout << "----------------------------> " << I1Initial << " " << I1 << " " << fSandler.fYC.fA << " " << fSandler.fYC.fR << " " << fPoreStressRZ[0] << " " << fPoreStressRZ[1] << endl;
    
    // should improve!!defLateral->value(i)
    TPZTensor<STATE> hidro, epstotal;
    for (STATE i1 = I1Initial/10.; i1>=I1Initial; i1 += I1Initial/10.) {
        hidro.XX() = i1/3.;
        hidro.YY() = i1/3.;
        hidro.ZZ() = i1/3.;
cout << "----------------------------> i1 " << i1  << " " << hidro.XX() << endl; 
        fSandler.ApplyLoad(hidro, epstotal);
    }
    
    // load elastic to the initial stress
    TPZTensor<STATE> porestress;
    porestress.XX() = fPoreStressRZ[0];
    porestress.YY() = fPoreStressRZ[0];
    porestress.ZZ() = fPoreStressRZ[1];
    fSandler.ApplyLoad(porestress, epstotal);
    fPoreStrainRZ[0] = epstotal.XX();
    if (fabs(epstotal.XX()-epstotal.YY()) > 1.e-8) {
        DebugStop();
    }
    fPoreStrainRZ[1] = epstotal.ZZ();
    
}

/// Evoluate the stress and strain to step ist
void TPZPlasticityTest::EvoluateToStep(TPZVec<STATE> &strainRZ, TPZVec<STATE> &stressRZ)
{
    TPZTensor<STATE> strain,stressgoal;
    TPZPlasticState<STATE> locstate;
    locstate = fSandler.GetState();
    STATE residual;
    
    strain.XX() = strainRZ[0];
    strain.YY() = strain.XX();
    strain.ZZ() = strainRZ[1];
    // if the stress is imposed, then this is its target value
    stressgoal.XX() = stressRZ[0];
    stressgoal.YY() = stressgoal.XX();
    stressgoal.ZZ() = stressRZ[1];
    do {
        
        fSandler.SetState(locstate);
        fSandler.fYC.fIsonCap = false;
        TPZFNMatrix<36,STATE> dep(6,6,0.);
        TPZTensor<STATE> stress;
        if(fRStressKnown || fZStressKnown)
        {
            fSandler.ApplyStrainComputeDep(strain, stress, dep);
        }
        else {
            fSandler.ApplyStrainComputeSigma(strain, stress);
        }
        residual = 0.;
        if (fRStressKnown) {
            residual += fabs(stress.XX()-stressgoal.XX());
        }
        if (fZStressKnown) {
            residual += fabs(stress.ZZ()-stressgoal.ZZ());
        }
        if (fabs(stress.XX()-stress.YY()) > 1.e-8) {
            DebugStop();
        }
        TPZFNMatrix<4,STATE> depsmall(2,2,0.);
        depsmall(0,0) = dep(_XX_,_XX_)+dep(_XX_,_YY_);
        depsmall(0,1) = dep(_XX_,_ZZ_);
        depsmall(1,0) = dep(_ZZ_,_XX_)+dep(_ZZ_,_YY_);
        depsmall(1,1) = dep(_ZZ_,_ZZ_);
        
        if (!fRStressKnown) {
            depsmall(0,0) = 1.;
            depsmall(0,1) = 0.;
            depsmall(1,0) = 0.;
            stress.XX() = stressgoal.XX();
        }
        if (!fZStressKnown) {
            depsmall(1,1) = 1.;
            depsmall(0,1) = 0.;
            depsmall(1,0) = 0.;
            stress.ZZ() = stressgoal.ZZ();
        }
        
        STATE det = depsmall(0,0)*depsmall(1,1)-depsmall(0,1)*depsmall(1,0);
        TPZFNMatrix<4,STATE> depsmallinv(2,2);
        depsmallinv(0,0) = depsmall(1,1)/det;
        depsmallinv(1,1) = depsmall(0,0)/det;
        depsmallinv(1,0) = -depsmall(1,0)/det;
        depsmallinv(0,1) = -depsmall(0,1)/det;
        
        TPZManVector<STATE,2> stressres(2,0.);
        stressres[0] = -stressgoal.XX()+stress.XX();
        stressres[1] = -stressgoal.ZZ()+stress.ZZ();
        
        strain.XX() -= depsmallinv(0,0)*stressres[0]+depsmallinv(0,1)*stressres[1];
        strain.YY() = strain.XX();
        strain.ZZ() -= depsmallinv(1,0)*stressres[0]+depsmallinv(1,1)*stressres[1];
        
        // maybe a line search is in order??
        
        
        
    } while (residual > 1.e-8);
    
    fSandler.SetState(locstate);
    TPZTensor<STATE> sigma;
    fSandler.ApplyStrainComputeSigma(strain, sigma);
    
    strainRZ[0] = strain.XX();
    strainRZ[1] = strain.ZZ();
    stressRZ[0] = sigma.XX();
    stressRZ[1] = sigma.ZZ();
    
}

/// compute the stress strain curve
void TPZPlasticityTest::PerformSimulation()
{
    ApplyInitialStress();
    int nloadsteps = fStressRZInput.Rows();
    fStrainRZSimulated.Resize(nloadsteps, 2);
    fStressRZSimulated.Resize(nloadsteps, 2);
    for (int istep = 0; istep < nloadsteps; istep++) {
        TPZManVector<STATE,2> strainRZ(2),stressRZ(2);
        strainRZ[0] = fStrainRZInput(istep,0);
        stressRZ[0] = fStressRZInput(istep,0);
        strainRZ[1] = fStrainRZInput(istep,1);
        stressRZ[1] = fStressRZInput(istep,1);
        if(istep == 1494)
        {
            std::cout << __FUNCTION__ << std::endl;
        }
        EvoluateToStep(strainRZ, stressRZ);
        fStrainRZSimulated(istep,0) = strainRZ[0];
        fStrainRZSimulated(istep,1) = strainRZ[1];
        fStressRZSimulated(istep,0) = stressRZ[0];
        fStressRZSimulated(istep,1) = stressRZ[1];
        
    }
    
    
}

/// read the input strain and stress from the laboratory file
void TPZPlasticityTest::ReadInputStrainStress(const std::string &filename)
{
    std::ifstream input(filename.c_str());
    if (!input) {
        DebugStop();
    }
    int numlines = 0;
    char buf[1024];
    input.getline(buf , 1024);
    STATE x, sig_ax_t, tempo, sig_ax_dev, sig_r, eps_ax, eps_r, eps_v;
    while (input) {
        input >> x >> sig_ax_t >> tempo >> sig_ax_dev >> tempo >> sig_r >> tempo >> eps_ax >> tempo >> eps_r >> tempo >> eps_v;
        if (!input) {
            break;
        }
        if(numlines >= fStressRZInput.Rows())
        {
            fStressRZInput.Resize(numlines+100, 2);
            fStrainRZInput.Resize(numlines+100, 2);
        }
        fStressRZInput(numlines,1) = -sig_ax_t;
        fStressRZInput(numlines,0) = -sig_r;
        fStrainRZInput(numlines,1) = -eps_ax/100.;
        fStrainRZInput(numlines,0) = -eps_r/100.;
	
	cout << "i= " << numlines << " " <<  fStressRZInput(numlines,1) << " " <<  fStressRZInput(numlines,0)<< " " << fStrainRZInput(numlines,1)<< " "<< fStrainRZInput(numlines,0) << endl;
	
        numlines++;
    }
    fStressRZInput.Resize(numlines, 2);
    fStrainRZInput.Resize(numlines, 2);
    
}

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
     
     
     A = 116.67; //0.25;
     B = 0.0036895; //0.67;
     C = 111.48; //0.18;
     D = 0.018768; //0.67 / 2.; 
     R = 0.91969; //2.5;
     W = 0.006605; //0.066;
     
     
     sandler.SetUp(poisson, E,A, B, C, R, D, W);
     
     
     REAL inttol = 1.e-4;
     
     sandler.SetIntegrTol(inttol);
     
     test.SetSandlerDimaggio(sandler);

     test.PerformSimulation();
     
     
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

