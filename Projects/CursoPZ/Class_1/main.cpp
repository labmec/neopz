/**
 * @file
 * @brief Implements the use of the vectors and integration rules routines as tutorial example of the util and integral NeoPZ modules
 */

#include "pzvec.h"
#include "pzquad.h"
#include "tpzintrulep3d.h"

#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

// First example: Polinomial functions
REAL Funcao1D (TPZVec<REAL> &pt,int degree);
REAL Funcao2D (TPZVec<REAL> &pt,int degree);
REAL Funcao3D (TPZVec<REAL> &pt,int degree);

void PrintMaxOrder();

int main() {
    PrintMaxOrder();
    // variables
    TPZVec<REAL> point(3,0.);
    REAL weight = 0.;
    REAL integral;
    
    // output file
    ofstream integrate("integratedPolinomial.txt");
    //	integrate << setprecision(14);
    
    // degree of the numeric integration
    int p = 2;
    
    // degree d -> The function is f(x) = x^d - 2x^(d-1)
    int degree = 2;   // f(x) = x^2 - 2x
    int it, npoints;
    
    //=====1D Rule===Computing $f \int_{-1} f(x) $f==================================
    while(degree < 26) {
        TPZInt1d ordem1d (p);
        //	for(it2=0;it2<3;it2++) {
        integral = 0.0;
        //ordem1d.SetType(it2,p);
        ordem1d.SetType(0,p);
        npoints = ordem1d.NPoints();
        for (it=0;it<npoints;it++){
            ordem1d.Point(it,point,weight);
            integral += weight * Funcao1D(point,degree);
        }
        if(!npoints) {
            std::string nome;
            ordem1d.Name(nome);
            integrate << "Type = " << nome << endl;
        }
        //		ordem1d.Print(integrate);
        integrate << "1D : p = " << p << " - Integral = " << integral << endl;
        //}
        degree++;
        p=degree;
    }
    // Integral = 1.16666666666666667 (d = 2)
    // Integral = -1.2500000000000000 (d = 3) (minimum p = 2)
    // Integral = -1.8958333333333333 (d = 5)
    // Integral = -3.2031250000000000 (d = 7)
    // Integral = -5.7664062500000000 (d = 9)
    // Integral = 14.9707406850961538 (d = 12)
    // Integral = -41.052551269531250 (d = 15)
    // Integral = 116.675674840023643 (d = 18)
    // Integral = 237.518337885538766 (d = 20)
    // Integral = -1456.7981708095624 (d = 25)
    //=====End of 1D Rule================================
    
    //=====2D Rule===Computing $f \int f(x,y) $f=====
    // The function is f(x) = x^d + y^d - 2 x y
    //=====Triangle Rule==============================
    //point.Resize(2);
    degree = p = 2;
    while(degree < 12) {
        TPZIntTriang ordem2dt (p);
        npoints = ordem2dt.NPoints();
        integral = 0.0;
        for (it=0;it<npoints;it++){
            ordem2dt.Point(it,point,weight);
            integral += weight * Funcao2D(point,degree);
        }
        integrate << "2D - Triangle Integral = " << integral << endl;
        degree++;
        p= degree;
    }
    
    //=====Quad Rule==================================
    degree = p = 2;
    while(degree < 12) {
        TPZIntQuad ordem2dq (p,p);
        npoints = ordem2dq.NPoints();
        integral = 0.0;
        for (it=0;it<npoints;it++){
            ordem2dq.Point(it,point,weight);
            integral += weight * Funcao2D(point,degree);
        }
        integrate << "2D - Quad Integral = " << integral << endl;
        //=====End 2D integration============================
        degree++;
        p= degree;
    }
    
    //=====3D Rule===Computing $f \int f(x,y,z) $f=====
    // The function is f(x) = x^d + y^d + z^d - 2 x y z
    //=====Tetrahedra Rule==================================
    //point.Resize(3);
    TPZIntTetra3D ordem3dt (p);
    npoints = ordem3dt.NPoints();
    integral = 0.0;
    for (it=0;it<npoints;it++){
        ordem3dt.Point(it,point,weight);
        integral += weight * Funcao3D(point,degree);
    }
    integrate << "3D - Tetra Integral = " << integral << endl;
    
    //=====Pyramid Rule==================================
    TPZIntPyram3D ordem3dpy (p);
    npoints = ordem3dpy.NPoints();
    integral = 0.0;
    for (it=0;it<npoints;it++){
        ordem3dpy.Point(it,point,weight);
        integral += weight * Funcao3D(point,degree);
    }
    integrate << "3D - Pyramid Integral = " << integral << endl;
    
    //=====Prism Rule===================================
    TPZIntPrism3D ordem3dpr (p);
    npoints = ordem3dpr.NPoints();
    integral = 0.0;
    for (it=0;it<npoints;it++){
        ordem3dpr.Point(it,point,weight);
        integral += weight * Funcao3D(point,degree);
    }
    integrate << "3D - Prism Integral = " << integral << endl;
    
    //=====Hexahedra Rule====================================
    TPZIntCube3D ordem3dc (p);
    npoints = ordem3dc.NPoints();
    integral = 0.0;
    for (it=0;it<npoints;it++){
        ordem3dc.Point(it,point,weight);
        integral += weight * Funcao3D(point,degree);
    }
    integrate << "3D - Hexa Integral = " << integral << endl;
    
    // close the output file
    integrate.close();
    
    return 0;
}

REAL Funcao1D (TPZVec<REAL> &pt,int degree) {
    if(pt.NElements() < 0)
        DebugStop();
    REAL a=1.;
    
    for(int i=0;i<degree;i++) {
        a *= (pt[0]);
    }
    return (a + 3.*pt[0]);
}
REAL Funcao2D (TPZVec<REAL> &pt,int degree) {
    if(pt.NElements() < 0)
        DebugStop();
    REAL a=1., b=1.;
    
    for(int i=0;i<degree;i++) {
        a *= (pt[0] - 0.5);
        b *= pt[1];
    }
    return (a+b + 3.*pt[0]*pt[1]);
}
REAL Funcao3D (TPZVec<REAL> &pt,int degree) {
    if(pt.NElements() < 0)
        DebugStop();
    REAL a=1., b=1., c=1.;
    
    for(int i=0;i<degree;i++) {
        a *= (pt[0] - 0.5);
        b *= pt[1];
        c *= pt[2];
    }
    return (a+b+c + 3.*pt[0]*pt[1]*pt[2]);
}

//Functions declaration
REAL Funcao (TPZVec<REAL> &pt, int i, int j, int k, int p);
REAL Alfa (int i, int j, int k, int p);

int main_2() {
    
    int i = 2;
    int j = 2;
    int k = 2;
    int p = 6;
    TPZVec<REAL> point(1,0.);
    REAL weight = 0.;
    REAL integral;
    
    //=====1D Rule=====================================
    TPZInt1d ordem1d (p);
    int it, it2, npoints;
    for(it2=0;it2<2;it2++) {
        integral = 0.0;
        ordem1d.SetType(it2,p);
        npoints = ordem1d.NPoints();
        for (it=0;it<npoints;it++){
            ordem1d.Point(it,point,weight);
            integral += weight * Funcao(point,i,j,k,p);
        }
        ordem1d.Print();
        cout << "1D - Integral = " << integral << endl;
    }
    //=====End of 1D Rule================================
    
    //=====2D Triangle Rule==============================
    point.Resize(2);
    TPZIntTriang ordem2dt (p);
    npoints = ordem2dt.NPoints();
    integral = 0.0;
    for (it=0;it<npoints;it++){
        ordem2dt.Point(it,point,weight);
        integral += weight * Funcao(point,i,j,k,p);
    }
    //=====End of 2D Triangle Rule=======================
    cout << "2D - Triangle Integral = " << integral << endl;
    
    //=====2D Quad Rule==================================
    TPZIntQuad ordem2dq (p,p);
    npoints = ordem2dq.NPoints();
    integral = 0.0;
    for (it=0;it<npoints;it++){
        ordem2dq.Point(it,point,weight);
        integral += weight * Funcao(point,i,j,k,p);
    }
    //=====End of 2D Quad Rule===========================
    cout << "2D - Quad Integral = " << integral << endl;
    
    //=====3D Tetra Rule==================================
    point.Resize(3);
    TPZIntTetra3D ordem3dt (p);
    npoints = ordem3dt.NPoints();
    integral = 0.0;
    for (it=0;it<npoints;it++){
        ordem3dt.Point(it,point,weight);
        integral += weight * Funcao(point,i,j,k,p);
    }
    //=====End of 3D Tetra Rule==========================
    cout << "3D - Tetra Integral = " << integral << endl;
    
    //=====3D Pyramid Rule==================================
    TPZIntPyram3D ordem3dpy (p);
    npoints = ordem3dpy.NPoints();
    integral = 0.0;
    for (it=0;it<npoints;it++){
        ordem3dpy.Point(it,point,weight);
        integral += weight * Funcao(point,i,j,k,p);
    }
    //=====End of 3D Pyramid Rule==========================
    cout << "3D - Pyramid Integral = " << integral << endl;
    
    //=====3D Prism Rule===================================
    TPZIntPrism3D ordem3dpr (p);
    npoints = ordem3dpr.NPoints();
    integral = 0.0;
    for (it=0;it<npoints;it++){
        ordem3dpr.Point(it,point,weight);
        integral += weight * Funcao(point,i,j,k,p);
    }
    //=====End of 3D Prism Rule============================
    cout << "3D - Prism Integral = " << integral << endl;
    
    //=====3D Hexa Rule====================================
    TPZIntCube3D ordem3dc (p);
    npoints = ordem3dc.NPoints();
    integral = 0.0;
    for (it=0;it<npoints;it++){
        ordem3dc.Point(it,point,weight);
        integral += weight * Funcao(point,i,j,k,p);
    }
    //=====End of 3D Hexa Rule=============================
    cout << "3D - Hexa Integral = " << integral << endl;
    
    /** Knowing the points and weight to some Gaussian rules */
    int np = 1;
    std::ofstream nome("PyramidQuad.txt",ios::app);
    while(np < 17) {
        //		cout << "\nOrder of the polinomials integrated exactly by gaussian rule (0 - quit) : ";
        //		cin >> np;
        TPZGaussRule *grule;
        grule = TPZIntRuleList::gIntRuleList.GetRule(np);
        grule->Print(nome);
        np++;
    }
    np = 1;
    while(np < 7) {
        TPZIntRuleP3D *rule;
        rule = TPZIntRuleList::gIntRuleList.GetRuleP3D(np);
        rule->Print(nome);
        //		rule->ComputePyramidPointsAndWeights(np,np+1);
        //		rule->Print(nome);
        np++;
    }
    nome.close();
    
    return 0;
}


//Functions definition
REAL Funcao(TPZVec<REAL> &pt, int i, int j, int k, int p){
    int r,s,t;
    REAL result=0.;
    REAL x = 0.;
    REAL y = 0.;
    REAL z = 0.;
    
    //To treat 1D and 2D
    int dim = pt.NElements();
    if (dim > 2) z = pt[2];
    else z = 1.;
    if (dim > 1) y = pt[1];
    else y = 1.;
    x = pt[0];
    //cout << "Point coordinate:\n";
    //cout << "x = " << x << "  y = " << y << "  z = " << z << endl;
    
    for (r=0;r<=i;r++){
        for(s=0;s<=j;s++){
            for (t=0;t<=k;t++){
                result +=  Alfa(r,s,t,p) * (pow (x,r) * pow (y,s) * sin(z));
            }
        }
    }
    return result;
};

REAL Alfa(int i, int j, int k, int p){
    int n = i*p*p + j*p + k;
    REAL alfa = sin((REAL)n);
    //cout << "alfa = " << alfa << endl;
    return alfa;
}

void PrintMaxOrder()
{
    int64_t p = 500;
    TPZInt1d ordem1d (p);
    TPZVec<int> order(1);
    ordem1d.GetOrder(order);
    std::cout << "Numero de pontos 1D " << ordem1d.NPoints() << " Ordem obtida " << order << " Max Order " << ordem1d.GetMaxOrder() << std::endl;
    
    order.resize(2);
    TPZIntTriang ordem2dt (p);
    ordem2dt.GetOrder(order);
    std::cout << "Numero de pontos Triangle " << ordem2dt.NPoints() << " Ordem obtida " << order << " Max Order " << ordem2dt.GetMaxOrder() << std::endl;
    
    TPZIntQuad ordem2dq (p,p);
    ordem2dq.GetOrder(order);
    std::cout << "Numero de pontos Quadrilateral " << ordem2dq.NPoints() << " Ordem obtida " << order << " Max Order " << ordem2dq.GetMaxOrder() << std::endl;
    
    order.resize(3);
    
    TPZIntTetra3D ordem3dt (p);
    ordem3dt.GetOrder(order);
    std::cout << "Numero de pontos Tetraedron " << ordem3dt.NPoints() << " Ordem obtida " << order << " Max Order " << ordem3dt.GetMaxOrder() << std::endl;
    
    TPZIntPyram3D ordem3dpy (p);
    ordem3dpy.GetOrder(order);
    std::cout << "Numero de pontos Pyramid " << ordem3dpy.NPoints() << " Ordem obtida " << order << " Max Order " << ordem3dpy.GetMaxOrder() << std::endl;
    
    
    TPZIntPrism3D ordem3dpr (p,p);
    ordem3dpr.GetOrder(order);
    std::cout << "Numero de pontos Prism " << ordem3dpr.NPoints() << " Ordem obtida " << order << " Max Order " << ordem3dpr.GetMaxOrder() << std::endl;
    
    TPZIntCube3D ordem3dc (p,p,p);
    ordem3dc.GetOrder(order);
    std::cout << "Numero de pontos Hexahedron " << ordem3dc.NPoints() << " Ordem obtida " << order << " Max Order " << ordem3dc.GetMaxOrder() << std::endl;
    
}