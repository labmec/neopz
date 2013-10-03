#include <iostream>
#include <cstdlib>


#include "simulation.h"

/// Default constructor
TPZPlasticityTest::TPZPlasticityTest() : fZStressKnown(false), fRStressKnown(false), fPoreStressRZ(2,0.), fPoreStrainRZ(2,0.), fStressRZInput(0,2,0.),
                fStrainRZInput(0,2,0.), fStressRZSimulated(0,2,0.), fStrainRZSimulated(0,2,0.), fNumSteps(200)
{

}

/// Get the stress to the pore stress
void TPZPlasticityTest::ApplyInitialStress()
{
    // load hydrostatic till I1 is in the middle of the ellips
    STATE I1 = fPoreStressRZ[0]*2.+fPoreStressRZ[1];
    STATE I1Initial = I1-fSandler.fYC.fA*fSandler.fYC.fR;
    // should improve!!
    TPZTensor<STATE> hidro, epstotal;
    for (STATE i1 = I1Initial/10.; i1>=I1Initial; i1 += I1Initial/10.) {
        hidro.XX() = i1/3.;
        hidro.YY() = i1/3.;
        hidro.ZZ() = i1/3.;
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
        fSandler.ApplyStrainComputeDep(strain, stress, dep);

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
        if(istep == 2)
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
        numlines++;
    }
    fStressRZInput.Resize(numlines, 2);
    fStrainRZInput.Resize(numlines, 2);

}

//#include "TPZPlasticityTest.h"
//
void MaterialPointTests()
{
    ///TESTE UCS

    ofstream outUCS("UCS.txt");
    ofstream outConfinamento("Confinamento03.txt");
    TPZTensor<REAL> deltaeps,eps,sigma,deltasigma;


    TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1> sandler;
    REAL poisson, E, A, B, C, R, D, W;
    E = 100; //ksi
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



using namespace std;
void MaterialPointTests();

 int main2()
 {
     //InitializePZLOG();
     TPZPlasticityTest test;
     std::string filename("../ensaio_all_columns.txt");
     test.ReadInputStrainStress(filename);

     TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> sandler;
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

     test.SetSandlerDimaggio(sandler);

     test.PerformSimulation();


     cout << "QT GUI" << endl;
     MaterialPointTests();
     cout << "\n End runing ! "<< endl;

 return EXIT_SUCCESS;
 }


































//    Simulation s;
//    s.setupTriaxial();

//    qDebug() << "SIZEEEEEEEEE " <<s.getXsize();

//    for (int iii=0; iii<= s.getXsize(); iii++){
//        qDebug() << "value X[] " << s.getXvalues()[iii] << "value Y[] " << s.getYvalues()[iii];
//    }

//    TXT arquivo_3;
//    arquivo_3.name = "Triaxial";
//    arquivo_3.X = s.getXvalues();
//    arquivo_3.Y = s.getYvalues();
//    arquivo_3.size = s.getXsize();
//    position = this->FilesList->size();
//    this->FilesList->insert(position, arquivo_3 );
//    //Criando entrada na listWidget
//    QListWidgetItem *item3;
//    item3 = new QListWidgetItem();
//    item3->setFlags(item3->flags() | Qt::ItemIsUserCheckable);
//    item3->setCheckState(Qt::Unchecked);
//    item3->setText(this->FilesList->value(position).name);
//    ui->listWidget->addItem(item3);


//    s.setupUCS();


//    qDebug() << "SIZEEEEEEEEE " <<s.getXsize();

//    for (int iii=0; iii<= s.getXsize(); iii++){
//        qDebug() << "value X[] " << s.getXvalues()[iii] << "value Y[] " << s.getYvalues()[iii];
//    }

//    TXT arquivo_4;
//    arquivo_4.name = "UCS";
//    arquivo_4.X = s.getXvalues();
//    arquivo_4.Y = s.getYvalues();
//    arquivo_4.size = s.getXsize();
//    position = this->FilesList->size();
//    this->FilesList->insert(position, arquivo_4 );
//    //Criando entrada na listWidget
//    QListWidgetItem *item4;
//    item4 = new QListWidgetItem();
//    item4->setFlags(item4->flags() | Qt::ItemIsUserCheckable);
//    item4->setCheckState(Qt::Unchecked);
//    item4->setText(this->FilesList->value(position).name);
//    ui->listWidget->addItem(item4);
































//#include "simulation.h"
//#include <QDebug>

//using namespace std;

//Simulation::Simulation()
//{
//    _E = 100; //ksi
//    _poisson = 0.40;
//    _A = 0.25;
//    _B = 0.67;
//    _C = 0.18;
//    _D = 0.67 / 2.; //(sugg. 0.335)
//    _R = 2.5;
//    _W = 0.066;

//    _inttol = 1.e-4;

//    _nofsteps = 100;
//    _unloadstep = 49;

//    // P/ UCS
//    _epsx = -0.0013;

//    // P/ Triaxial
//    _deltasigmaXX = -0.004;
//    _percent = 41;


//    sandler = new TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1>;
//}

//double * Simulation::getXvalues() {
//    return x.data();
//}

//double * Simulation::getYvalues() {
//    return y.data();
//}

//int Simulation::getXsize() {
//    return x.size();
//}

//void Simulation::setupUCS() {

//    // x y solution storage
//    x.clear();
//    y.clear();

//    TPZTensor<REAL> deltaeps;
//    TPZTensor<REAL> eps;
//    TPZTensor<REAL> sigma;
//    TPZTensor<REAL> deltasigma;

//    sandler->SetUp(_poisson, _E, _A, _B, _C, _R, _D, _W);
//    sandler->SetIntegrTol(_inttol);

//    deltaeps.XX() = _epsx;
//    deltaeps.YY() = 0;
//    deltaeps.ZZ() = 0;

//    eps=deltaeps;

//    qDebug() << "SAIDA UCS: eps.XX()  sigma.XX()" << endl;

//    for(int i=0;i<_nofsteps;i++)
//    {

//        sandler->ApplyStrainComputeSigma(eps, sigma);//UCS
//        // sandler.ApplyLoad(sigma, eps);

//        TPZPlasticState<REAL> state = sandler->GetState();
//        //state.Print(cout);
//        if(i==_unloadstep)
//        {
//            deltaeps*=-1;
//        }

//        //qDebug() << eps.XX() << " " << sigma.XX() << endl; // SAIDA UCS!!!

//        x.push_back(eps.XX());
//        y.push_back(sigma.XX());

//        eps+=deltaeps;
//    }
//}


//void Simulation::setupTriaxial() {

//    // x y solution storage
//    x.clear();
//    y.clear();

//    TPZTensor<REAL> deltaeps;
//    TPZTensor<REAL> eps;
//    TPZTensor<REAL> sigma;
//    TPZTensor<REAL> deltasigma;

//    sandler->SetUp(_poisson, _E, _A, _B, _C, _R, _D, _W);
//    sandler->SetIntegrTol(_inttol);

//    deltasigma.XX()=_deltasigmaXX; //deltasigma.XX()=-0.004;

//    deltasigma.YY()=deltasigma.XX()*(_percent/100);
//    deltasigma.ZZ()=deltasigma.YY();
//    sigma=deltasigma;

//    qDebug() << "SAIDA CONFINAMENTO: fabs(eps.XX())  fabs(sigma.XX())" << endl;
//    for(int i=0;i<_nofsteps;i++)
//    {
//        sandler->ApplyLoad(sigma,eps);

//        //qDebug() << fabs(eps.XX()) << " " << fabs(sigma.XX()) << endl; // SAIDA CONFINAMENTO!!!

//        x.push_back(fabs(eps.XX()));
//        y.push_back(fabs(sigma.XX()));

//        sigma+=deltasigma;
//    }
//    //VisualizeSandlerDimaggio(fileName,pSD);
//}
