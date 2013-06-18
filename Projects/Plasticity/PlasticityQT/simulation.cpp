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
































#include "simulation.h"
#include <QDebug>

using namespace std;

Simulation::Simulation()
{
    _E = 100; //ksi
    _poisson = 0.40;
    _A = 0.25;
    _B = 0.67;
    _C = 0.18;
    _D = 0.67 / 2.; //(sugg. 0.335)
    _R = 2.5;
    _W = 0.066;

    _inttol = 1.e-4;

    _nofsteps = 100;
    _unloadstep = 49;

    // P/ UCS
    _epsx = -0.0013;

    // P/ Triaxial
    _deltasigmaXX = -0.004;
    _percent = 41;


    sandler = new TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1>;
}

double * Simulation::getXvalues() {
    return x.data();
}

double * Simulation::getYvalues() {
    return y.data();
}

int Simulation::getXsize() {
    return x.size();
}

void Simulation::setupUCS() {

    // x y solution storage
    x.clear();
    y.clear();

    TPZTensor<REAL> deltaeps;
    TPZTensor<REAL> eps;
    TPZTensor<REAL> sigma;
    TPZTensor<REAL> deltasigma;

    sandler->SetUp(_poisson, _E, _A, _B, _C, _R, _D, _W);
    sandler->SetIntegrTol(_inttol);

    deltaeps.XX() = _epsx;
    deltaeps.YY() = 0;
    deltaeps.ZZ() = 0;

    eps=deltaeps;

    qDebug() << "SAIDA UCS: eps.XX()  sigma.XX()" << endl;

    for(int i=0;i<_nofsteps;i++)
    {

        sandler->ApplyStrainComputeSigma(eps, sigma);//UCS
        // sandler.ApplyLoad(sigma, eps);

        TPZPlasticState<REAL> state = sandler->GetState();
        //state.Print(cout);
        if(i==_unloadstep)
        {
            deltaeps*=-1;
        }

        //qDebug() << eps.XX() << " " << sigma.XX() << endl; // SAIDA UCS!!!

        x.push_back(eps.XX());
        y.push_back(sigma.XX());

        eps+=deltaeps;
    }
}


void Simulation::setupTriaxial() {

    // x y solution storage
    x.clear();
    y.clear();

    TPZTensor<REAL> deltaeps;
    TPZTensor<REAL> eps;
    TPZTensor<REAL> sigma;
    TPZTensor<REAL> deltasigma;

    sandler->SetUp(_poisson, _E, _A, _B, _C, _R, _D, _W);
    sandler->SetIntegrTol(_inttol);

    deltasigma.XX()=_deltasigmaXX; //deltasigma.XX()=-0.004;

    deltasigma.YY()=deltasigma.XX()*(_percent/100);
    deltasigma.ZZ()=deltasigma.YY();
    sigma=deltasigma;

    qDebug() << "SAIDA CONFINAMENTO: fabs(eps.XX())  fabs(sigma.XX())" << endl;
    for(int i=0;i<_nofsteps;i++)
    {
        sandler->ApplyLoad(sigma,eps);

        //qDebug() << fabs(eps.XX()) << " " << fabs(sigma.XX()) << endl; // SAIDA CONFINAMENTO!!!

        x.push_back(fabs(eps.XX()));
        y.push_back(fabs(sigma.XX()));

        sigma+=deltasigma;
    }
    //VisualizeSandlerDimaggio(fileName,pSD);
}
