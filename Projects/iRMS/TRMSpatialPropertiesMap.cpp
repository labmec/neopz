//
//  TRMSpatialPropertiesMap.cpp
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#include "TRMSpatialPropertiesMap.h"
#include <boost/math/distributions/normal.hpp>

/** @brief default constructor */
TRMSpatialPropertiesMap::TRMSpatialPropertiesMap(){
    
}

/** @brief default destructor */
TRMSpatialPropertiesMap::~TRMSpatialPropertiesMap(){
    
}

/** @brief Absolute Permeability m2  $\kappa$ */
void TRMSpatialPropertiesMap::Kappa(TPZManVector<STATE,3> &x, TPZFMatrix<STATE> &kappa, TPZFMatrix<STATE> &inv_kappa, TPZManVector<STATE,10> &state_vars){
    
    switch (this->MapModel()) {
        case 0:
        {
            this->Kappa_c(x, kappa, inv_kappa, state_vars);
        }
            break;
        case 1:
        {
            DebugStop();
        }
            break;
        case 2:
        {
            DebugStop();
        }
            break;
        default:
        {
            std::cout << "Error: Model not implemented." << std::endl;
            DebugStop();
        }
            break;
    }
    
}

/** @brief Absolute Permeability m2  $\kappa$ */
void TRMSpatialPropertiesMap::phi(TPZManVector<STATE,3> &x, TPZManVector<STATE,10> &phi, TPZManVector<STATE,10> &state_vars){
    
    switch (this->MapModel()) {
        case 0:
        {
            this->phi_c(x, phi, state_vars);
        }
            break;
        case 1:
        {
            DebugStop();
        }
            break;
        case 2:
        {
            DebugStop();
        }
            break;
        default:
        {
            std::cout << "Error: Model not implemented." << std::endl;
            DebugStop();
        }
            break;
    }
    
}

/** @brief Absolute Permeability m2  $\kappa$ */
void TRMSpatialPropertiesMap::Kappa_c(TPZManVector<STATE,3> &x, TPZFMatrix<STATE> &kappa, TPZFMatrix<STATE> &inv_kappa, TPZManVector<STATE,10> &state_vars){

//    kappa.Resize(3,3);
//    kappa.Zero();
//    STATE val = 1.0e-13;
//    REAL kx = sin(100.0*((x[0]*x[1])/1000.0)) + 1.1;
//    REAL ky = cos(5.0*(x[1]/100.0)) + 1.1;
//    kappa(0,0) = val*fabs(kx*phi_cdf(x[0]/1000.0));
//    kappa(1,1) = val*fabs(ky*phi_cdf(x[1]/100.0));
//    kappa(2,2) = val;
//    
//    inv_kappa.Resize(3,3);
//    inv_kappa.Zero();
//    inv_kappa(0,0) = 1.0/kappa(0,0);
//    inv_kappa(1,1) = 1.0/kappa(1,1);
//    inv_kappa(2,2) = 1.0/kappa(2,2);
    
    kappa.Resize(3,3);
    kappa.Zero();
    STATE val = 1.0e-13;
    REAL epsilon = 500.0;
    REAL kx = (2.0 + 1.8*sin(2.0*M_PI*x[0]*x[1]/epsilon))/(2.0 + 1.8*sin(2.0*M_PI*x[1]/epsilon));
    REAL ky = (2.0 + 1.8*sin(2.0*M_PI*x[0]*x[1]/epsilon))/(2.0 + 1.8*sin(2.0*M_PI*x[0]/epsilon));
    kappa(0,0) = val*fabs(kx+ky);
    kappa(1,1) = val*fabs(kx+ky);
    kappa(2,2) = val;
    
    inv_kappa.Resize(3,3);
    inv_kappa.Zero();
    inv_kappa(0,0) = 1.0/kappa(0,0);
    inv_kappa(1,1) = 1.0/kappa(1,1);
    inv_kappa(2,2) = 1.0/kappa(2,2);
    
}




/** @brief Porosity fraction  $\phi$ */
void TRMSpatialPropertiesMap::phi_c(TPZManVector<STATE,3> &x, TPZManVector<STATE,10> &phi, TPZManVector<STATE,10> &state_vars){
    
//    phi.Resize(10, 0.0);
//    STATE val = 0.25;
//    REAL kx = sin(100.0*((x[0]*x[1])/1000.0)) + 1.1;
//    REAL ky = cos(5.0*(x[1]/100.0)) + 1.1;
//    val *= fabs(kx*ky)*0.75;
//    phi[0] = val;
    
    phi.Resize(10, 0.0);
    STATE val = 0.25;
    REAL epsilon = 500.0;
    REAL kx = (2.0 + 1.8*sin(2.0*M_PI*x[0]*x[1]/epsilon))/(2.0 + 1.8*sin(2.0*M_PI*x[1]/epsilon));
    REAL ky = (2.0 + 1.8*sin(2.0*M_PI*x[0]*x[1]/epsilon))/(2.0 + 1.8*sin(2.0*M_PI*x[0]/epsilon));
    val *= fabs(kx + ky)*0.1;
    phi[0] = val;
    
}

/** @brief Absolute Permeability m2  $\kappa$ */
void TRMSpatialPropertiesMap::Kappa_f(TPZManVector<STATE,3> &x, TPZFMatrix<STATE> &kappa, TPZFMatrix<STATE> &inv_kappa, TPZManVector<STATE,10> &state_vars){

    DebugStop();
}

/** @brief Porosity fraction  $\phi$ */
void TRMSpatialPropertiesMap::phi_f(TPZManVector<STATE,3> &x, TPZManVector<STATE,10> &phi, TPZManVector<STATE,10> &state_vars){
    
    DebugStop();
}

double TRMSpatialPropertiesMap::phi_cdf(double x)
{
    static const double RT2PI = sqrt(4.0*acos(0.0));
    
    static const double SPLIT = 7.07106781186547;
    
    static const double N0 = 220.206867912376;
    static const double N1 = 221.213596169931;
    static const double N2 = 112.079291497871;
    static const double N3 = 33.912866078383;
    static const double N4 = 6.37396220353165;
    static const double N5 = 0.700383064443688;
    static const double N6 = 3.52624965998911e-02;
    static const double M0 = 440.413735824752;
    static const double M1 = 793.826512519948;
    static const double M2 = 637.333633378831;
    static const double M3 = 296.564248779674;
    static const double M4 = 86.7807322029461;
    static const double M5 = 16.064177579207;
    static const double M6 = 1.75566716318264;
    static const double M7 = 8.83883476483184e-02;
    
    const double z = fabs(x);
    double c = 0.0;
    
    if(z<=37.0)
    {
        const double e = exp(-z*z/2.0);
        if(z<SPLIT)
        {
            const double n = (((((N6*z + N5)*z + N4)*z + N3)*z + N2)*z + N1)*z + N0;
            const double d = ((((((M7*z + M6)*z + M5)*z + M4)*z + M3)*z + M2)*z + M1)*z + M0;
            c = e*n/d;
        }
        else
        {
            const double f = z + 1.0/(z + 2.0/(z + 3.0/(z + 4.0/(z + 13.0/20.0))));
            c = e/(RT2PI*f);
        }
    }
    return x<=0.0 ? c : 1-c;
}