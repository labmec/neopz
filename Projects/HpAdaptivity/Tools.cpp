/*
 *  Tools.cpp
 *  PZ
 *
 *  Created by Denise de Siqueira on 9/5/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include "Tools.h"
#include "pzpoisson3d.h"
#include "pzgeoelbc.h"
#include "tpzgeoelrefpattern.h"
#include "TPZGeoLinear.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "pzgengrid.h"
#include "TPZSkylineNSymStructMatrix.h"

#include <cmath>

int const bc0=-1; //em y=0
int const bc1=-2; //em x=1
int const bc2=-3; //em y=1
int const bc3=-4; //em x=0

int const dirichlet =0;
int const neumann = 1;
int matId =1;


#ifdef LOG4CXX

static LoggerPtr logger(Logger::getLogger("HpAdaptivity.main"));

#endif

void PrintGMeshVTK(TPZGeoMesh * gmesh, std::ofstream &file)
{
    file.clear();
    int nelements = gmesh->NElements();
    
    std::stringstream node, connectivity, type;
    
    //Header
    file << "# vtk DataFile Version 3.0" << std::endl;
    file << "TPZGeoMesh VTK Visualization" << std::endl;
    file << "ASCII" << std::endl << std::endl;
    
    file << "DATASET UNSTRUCTURED_GRID" << std::endl;
    file << "POINTS ";
    
    int actualNode = -1, size = 0, nVALIDelements = 0;
    
    for(int el = 0; el < nelements; el++)
    {
        if(gmesh->ElementVec()[el]->Type() == EPoint)//Exclude Lines and Arc3D
        {
            continue;
        }
        if(gmesh->ElementVec()[el]->Type() == EOned)//Exclude Lines and Arc3D
        {
            continue;
        }
        if(gmesh->ElementVec()[el]->HasSubElement())
        {
            continue;
        }
        
        int elNnodes = gmesh->ElementVec()[el]->NNodes();
        size += (1+elNnodes);
        connectivity << elNnodes;
        
        for(int t = 0; t < elNnodes; t++)
        {
            for(int c = 0; c < 3; c++)
            {
                double coord = gmesh->NodeVec()[gmesh->ElementVec()[el]->NodeIndex(t)].Coord(c);
                node << coord << " ";
            }
            node << std::endl;
            
            actualNode++;
            connectivity << " " << actualNode;
        }
        connectivity << std::endl;
        
        int elType = -1;
        switch (gmesh->ElementVec()[el]->Type())
        {
            case (ETriangle):
            {
                elType = 5;
                break;
            }
            case (EQuadrilateral ):
            {
                elType = 9;
                break;
            }
            case (ETetraedro):
            {
                elType = 10;
                break;
            }
            case (EPiramide):
            {
                elType = 14;
                break;
            }
            case (EPrisma):
            {
                elType = 13;
                break;
            }
            case (ECube):
            {
                elType = 12;
                break;
            }
            default:
            {
                //ElementType NOT Found!!!
                DebugStop();
                break;
            }
        }
        
        type << elType << std::endl;
        nVALIDelements++;
    }
    node << std::endl;
    actualNode++;
    file << actualNode << " float" << std::endl << node.str();
    
    file << "CELLS " << nVALIDelements << " ";
    
    file << size << std::endl;
    file << connectivity.str() << std::endl;
    
    file << "CELL_TYPES " << nVALIDelements << std::endl;
    file << type.str();
    
    file.close();
}
const REAL MyPi=4.*atan(1.);
//const REAL epsilon=1000.;
void Forcing1(const TPZVec<REAL> &pt, TPZVec<STATE> &disp) {
//    REAL x = pt[0];
//    REAL y = pt[1];
//    disp[0]= 2.*pow(MyPi,2)*sin(MyPi*x)*sin(MyPi*y);
  
    return;
}
void Forcing2(const TPZVec<REAL> &pt, TPZVec<STATE> &disp) {
    double x = pt[0];
    double y = pt[1];
  	disp[0]=(-1)*(-10*(pow(-1 + y,2)*pow(y,2) - 6*x*pow(-1 + y,2)*pow(y,2) + 5*exp(10*pow(y,2))*(1 + exp(10*pow(x,2)))*(-1 + x)*pow(x,2)*(-5 + x*(9 + 20*(-1 + x)*x))*pow(-1 + y,2)*pow(y,2) -
                       2*pow(x,3)*(1 + 6*(-1 + y)*y) + pow(x,4)*(1 + 6*(-1 + y)*y) + pow(x,2)*(1 + 6*(-1 + y)*y*(1 + (-1 + y)*y)) + exp(5*pow(x,2))*(-(exp(5*pow(x,2))*	(pow(-1 + x,2)*pow(x,2) - 6*pow(-1 + x,2)*pow(x,2)*y +
                                                                                                                                                                         (1 + 2*(-1 + x)*x*(3 + 4*x*(-7 + x*(12 + 25*(-1 + x)*x))))*pow(y,2) - 2*(1 + 2*(-1 + x)*x*(3 + 5*x*(-5 + x*(9 + 20*(-1 + x)*x))))*pow(y,3) +
                                                                                                                                                                         (1 + 2*(-1 + x)*x*(3 + 5*x*(-5 + x*(9 + 20*(-1 + x)*x))))*pow(y,4))) +	2*exp(10*pow(y,2))*(pow(-1 + y,2)*pow(y,2) - 6*x*pow(-1 + y,2)*pow(y,2) - 200*pow(x,5)*pow(-1 + y,2)*pow(y,2) + 100*pow(x,6)*pow(-1 + y,2)*pow(y,2) -
                                                                                                                                                                                                                                                                    2*pow(x,3)*(1 + (-1 + y)*y*(6 + 5*y*(-17 + 5*y*(5 + 8*(-1 + y)*y)))) + pow(x,4)*(1 + (-1 + y)*y*(6 + 5*y*(-39 + y*(47 + 40*(-1 + y)*y)))) + pow(x,2)*(1 + (-1 + y)*y*(6 + y*(-81 + y*(121 + 200*(-1 + y)*y)))))*
                                                                                                                                                     sinh(5*pow(x,2)))));
    return;
}

void SolExata(const TPZVec<REAL> &pt, TPZVec<STATE> &p, TPZFMatrix<STATE> &flux ) {
    double x = pt[0];
    double y = pt[1];
    TPZVec<REAL> disp;
    //double fator=(-1.)*(x*x+y*y);
    p[0]= sin(MyPi*x)*sin(MyPi*y);//x*exp(fator);//Solucao
    flux(0,0)= (-1.)*MyPi*cos(MyPi*x)*sin(MyPi*y);//(-1.)*exp(fator)+2.*x*x*exp(fator);//dx
    flux(1,0)= (-1.)*MyPi*cos(MyPi*y)*sin(MyPi*x);//2.*x*y*exp(fator);// dy
    flux(2,0)= 2*MyPi*MyPi*sin(MyPi*x)*sin(MyPi*y);//(-4.)*x*exp(fator)*(-2.+x*x+y*y);//coloco o divergetne aq para testar
    
    
    
    
}
void SolExata2(const TPZVec<REAL> &pt, TPZVec<STATE> &p, TPZFMatrix<STATE> &flux ) {
    double x=pt[0];
    double y=pt[1];
    p[0]=5*(-1 + exp(10*pow(x,2)))*(-1 + exp(10*pow(y,2)))*pow(1 - x,2)*pow(x,2)*pow(1 - y,2)*pow(y,2);
    flux(0,0)= -10*(-1 + exp(10*pow(y,2)))*(-1 + x)*x*(1 - 2*x + exp(10*pow(x,2))*(-1 + 2*x*(1 + 5*(-1 + x)*x)))*pow(-1 + y,2)*pow(y,2);
    flux(1,0)=-10*(-1 + exp(10*pow(x,2)))*pow(-1 + x,2)*pow(x,2)*(-1 + y)*y*(1 - 2*y + exp(10*pow(y,2))*(-1 + 2*y*(1 + 5*(-1 + y)*y)));
    flux(2,0)=(-1)*(-10*(pow(-1 + y,2)*pow(y,2) - 6*x*pow(-1 + y,2)*pow(y,2) + 5*exp(10*pow(y,2))*(1 + exp(10*pow(x,2)))*(-1 + x)*pow(x,2)*(-5 + x*(9 + 20*(-1 + x)*x))*pow(-1 + y,2)*pow(y,2) -
                         2*pow(x,3)*(1 + 6*(-1 + y)*y) + pow(x,4)*(1 + 6*(-1 + y)*y) + pow(x,2)*(1 + 6*(-1 + y)*y*(1 + (-1 + y)*y)) + exp(5*pow(x,2))*(-(exp(5*pow(x,2))*	(pow(-1 + x,2)*pow(x,2) - 6*pow(-1 + x,2)*pow(x,2)*y +
                                                                                                                                                                             (1 + 2*(-1 + x)*x*(3 + 4*x*(-7 + x*(12 + 25*(-1 + x)*x))))*pow(y,2) - 2*(1 + 2*(-1 + x)*x*(3 + 5*x*(-5 + x*(9 + 20*(-1 + x)*x))))*pow(y,3) +
                                                                                                                                                                             (1 + 2*(-1 + x)*x*(3 + 5*x*(-5 + x*(9 + 20*(-1 + x)*x))))*pow(y,4))) +	2*exp(10*pow(y,2))*(pow(-1 + y,2)*pow(y,2) - 6*x*pow(-1 + y,2)*pow(y,2) -
                                                                                                                                                                                                                                                                        200*pow(x,5)*pow(-1 + y,2)*pow(y,2) + 100*pow(x,6)*pow(-1 + y,2)*pow(y,2) - 	2*pow(x,3)*(1 + (-1 + y)*y*(6 + 5*y*(-17 + 5*y*(5 + 8*(-1 + y)*y)))) + pow(x,4)*(1 + (-1 + y)*y*(6 + 5*y*(-39 + y*(47 + 40*(-1 + y)*y)))) + pow(x,2)*(1 + (-1 + y)*y*(6 + y*(-81 + y*(121 + 200*(-1 + y)*y)))))*
                                                                                                                                                       sinh(5*pow(x,2)))));
}
void SolExata3(const TPZVec<REAL> &pt, TPZVec<STATE> &p, TPZFMatrix<STATE> &flux ) {
    double x = pt[0];
    double y= pt[1];
    p[0]=5.+3.*x+2.*y+4.*x*y;
    flux(0,0)=(-1.)*(3.+4.*y);
    flux(1,0)=(-1.)*(2.+4.*x);
    flux(2,0)=0;
    
}

#define Power pow
#define ArcTan atan
#define Sqrt sqrt
#define eps epsilon
#define Pi MyPi

void SolArcTan(const TPZVec<REAL> &pt, TPZVec<STATE> &p, TPZFMatrix<STATE> &flux){
//    REAL x = pt[0];
//    REAL y = pt[1];
    p[0]=0;
    flux(0,0)=0;
    flux(1,0)=0;
    flux(2,0)=0;

    
//    p[0]=5*(-1 + x)*x*(-1 + y)*y*(Pi/2. + ArcTan(Sqrt(eps)*(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2))));//8*(1-x)*x*(1-y)*y*atan(0.0625+2*sqrt(epsilon)*(0.0625-pow(x-0.5, 2)-pow(y-0.5, 2)));
//    
//    //px
//    flux(0,0)=(10*Sqrt(eps)*(-1 + x)*(-0.5 + x)*x*(-1 + y)*y)/(1 + eps*Power(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2),2)) -5*(-1 + x)*(-1 + y)*y*(Pi/2. + ArcTan(Sqrt(eps)*                                        (0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2)))) -  5*x*(-1 + y)*y*(Pi/2. + ArcTan(Sqrt(eps)*(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2))));
//    //(-1.)*((-32*(1. - x)*(-0.5 + x)*x*(1. - y)*y*sqrt(epsilon))/(1 + pow(0.0625 + 2*(0.0625 - pow(-0.5 + x,2) - pow(-0.5 + y,2))*sqrt(epsilon),2)) + 8*(1. - x)*(1. - y)*y*atan(0.0625 + 2*(0.0625 - pow(-0.5 + x,2) - pow(-0.5 + y,2))*sqrt(epsilon)) -                           8*x*(1. - y)*y*atan(0.0625 + 2*(0.0625 - pow(-0.5 + x,2) - pow(-0.5 + y,2))*sqrt(epsilon)));
//    //py
//    
//    flux(1,0)=(10*Sqrt(eps)*(-1 + x)*x*(-1 + y)*(-0.5 + y)*y)/
//    (1 + eps*Power(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2),2)) -
//    5*(-1 + x)*x*(-1 + y)*(Pi/2. + ArcTan(Sqrt(eps)*
//                                          (0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2)))) -
//    5*(-1 + x)*x*y*(Pi/2. + ArcTan(Sqrt(eps)*
//                                   (0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2))));
//    //(-1.)*((-32*(1. - x)*x*(1. - y)*(-0.5 + y)*y*sqrt(epsilon))/(1 + pow(0.0625 + 2*(0.0625 - pow(-0.5 + x,2) - pow(-0.5 + y,2))*sqrt(epsilon),2)) +8*(1. - x)*x*(1. - y)*atan(0.0625 + 2*(0.0625 - pow(-0.5 + x,2) - pow(-0.5 + y,2))*sqrt(epsilon)) - 8*(1. - x)*x*y*atan(0.0625 + 2*(0.0625 - pow(-0.5 + x,2) - pow(-0.5 + y,2))*sqrt(epsilon)));
//    //div
//    flux(2,0)=-5*(-1 + x)*x*((-2*Sqrt(eps))/
//                             (1 + eps*Power(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2),2)) -
//                             (8*Power(eps,1.5)*Power(-0.5 + x,2)*
//                              (0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2)))/
//                             Power(1 + eps*Power(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2),2),2))*
//    (-1 + y)*y + (4*Sqrt(eps)*(-0.5 + x)*(5*(-1 + x) + 5*x)*(-1 + y)*y)/
//    (1 + eps*Power(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2),2)) -
//    10*(-1 + y)*y*(Pi/2. + ArcTan(Sqrt(eps)*
//                                  (0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2)))) -
//    5*(-1 + x)*x*(((-2*Sqrt(eps))/
//                   (1 + eps*Power(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2),2)) -
//                   (8*Power(eps,1.5)*(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2))*
//                    Power(-0.5 + y,2))/
//                   Power(1 + eps*Power(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2),2),2))*
//                  (-1 + y)*y - (4*Sqrt(eps)*(-0.5 + y)*(-1 + 2*y))/
//                  (1 + eps*Power(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2),2)) +
//                  2*(Pi/2. + ArcTan(Sqrt(eps)*(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2)))));
    
    
    //(-1.)*((-8*(8*(1. - x) - 8*x)*(-0.5 + x)*(1. - y)*y*sqrt(epsilon))/ (1 + pow(0.0625 + 2*(0.0625 - pow(-0.5 + x,2) - pow(-0.5 + y,2))*sqrt(epsilon),2)) +8*(1. - x)*x*(1. - y)*y*((-4*sqrt(epsilon))/(1 + pow(0.0625 + 2*(0.0625 - pow(-0.5 + x,2) - pow(-0.5 + y,2))*sqrt(epsilon),2)) -(32*pow(-0.5 + x,2)*(0.0625 + 2*(0.0625 - pow(-0.5 + x,2) - pow(-0.5 + y,2))*sqrt(epsilon))*epsilon)/ pow(1 + pow(0.0625 + 2*(0.0625 - pow(-0.5 + x,2) - pow(-0.5 + y,2))*sqrt(epsilon),2),2)) + 8*(1. - x)*x*((-8*(1. - 2*y)*(-0.5 + y)*sqrt(epsilon))/ (1 + pow(0.0625 + 2*(0.0625 - pow(-0.5 + x,2) - pow(-0.5 + y,2))*sqrt(epsilon),2)) + (1. - y)*y*((-4*sqrt(epsilon))/(1 + pow(0.0625 + 2*(0.0625 - pow(-0.5 + x,2) - pow(-0.5 + y,2))*sqrt(epsilon),2)) - (32*pow(-0.5 + y,2)*(0.0625 + 2*(0.0625 - pow(-0.5 + x,2) - pow(-0.5 + y,2))*sqrt(epsilon))*epsilon)/pow(1 + pow(0.0625 + 2*(0.0625 - pow(-0.5 + x,2) - pow(-0.5 + y,2))*sqrt(epsilon),2),2)) -  2*atan(0.0625 + 2*(0.0625 - pow(-0.5 + x,2) - pow(-0.5 + y,2))*sqrt(epsilon))) -  16*(1. - y)*y*atan(0.0625 + 2*(0.0625 - pow(-0.5 + x,2) - pow(-0.5 + y,2))*sqrt(epsilon)));
    
    
    
    
    /*
     REAL F = 2.*sqrt(epsilon);
     REAL arc = F*(0.25*0.25 - (x - 0.5)*(x - 0.5) - (y - 0.5)*(y - 0.5));
     REAL prodx = x*(x-1.);
     REAL prody = y*(y-1.);
     REAL prod = prodx*prody;
     p[0] = 8.*prod*(1.+(2./MyPi)*atan(arc));//solucao p
     //gradiente
     
     REAL temp = (MyPi + 2.*atan(arc))*(2.*x-1);
     REAL frac = 0;
     REAL aux2x= (0.25*0.25)-(x-0.5)*(x-0.5)-(y-0.5)*(y-0.5);
     frac= (8.*prodx*(x-0.5)*sqrt(epsilon))/(1.+(4.*aux2x*aux2x)*epsilon);
     //componente ux
     flux(0,0) = (-1.)*(8./MyPi)*prody*(temp - frac);
     //componente uy
     
     REAL temp2 = (MyPi + 2.*atan(arc))*(2.*y-1.);
     frac= (8.*prody*(y-0.5)*sqrt(epsilon))/(1.+(4.*aux2x*aux2x*epsilon));
     //componente ux
     flux(1,0) = (-1.)*(8.*prodx)/MyPi*(temp2 - frac);
     
     
     //divergente
     //uxx
     
     REAL fracaux=(1.+(4.*aux2x*aux2x)*epsilon);
     REAL termo2=fracaux*fracaux;
     REAL fator1=prody*(MyPi+2.*atan(arc));
     REAL fator2=(4./fracaux)*prody*sqrt(epsilon)*((2*x-1)*(2*x-1)+x*(x-1));
     REAL fator3=(64./termo2)*(prod*(x-0.5)*(x-0.5)*aux2x*(pow(epsilon, 1.5)));
     REAL aux2xx= (16./MyPi)*(fator1-fator2-fator3);
     
     //uyy
     
     REAL fator1y=prodx*(MyPi+2.*atan(arc));
     REAL fator2y=(4./fracaux)*prodx*sqrt(epsilon)*((2*y-1)*(2*y-1)+y*(y-1));
     REAL fator3y=(64./termo2)*(prod*(y-0.5)*(y-0.5)*aux2x*(pow(epsilon,1.5)));
     REAL aux2yy= (16./MyPi)*(fator1y-fator2y-fator3y);
     
     //uxx+uyy
     
     flux(2,0)=(-1.)*(aux2xx+aux2yy);
     
     
     */
}
#define Power pow
#define ArcTan atan
#define Sqrt sqrt
void ForcingTang(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    //void ForcingTang(const TPZVec<REAL> &pt, TPZVec<REAL> &res,TPZFMatrix<STATE> &disp){
    //   disp.Redim(2,1);
//    double x = pt[0];
//    double y = pt[1];
    
    disp[0]=0;
    
//    disp[0]=-5*(-1 + x)*x*((-2*Sqrt(eps))/(1 + eps*Power(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2),2)) -(8*Power(eps,1.5)*Power(-0.5 + x,2)*(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2)))/
//                           Power(1 + eps*Power(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2),2),2))*(-1 + y)*y + (4*Sqrt(eps)*(-0.5 + x)*(5*(-1 + x) + 5*x)*(-1 + y)*y)/(1 + eps*Power(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2),2)) -10*(-1 + y)*y*(Pi/2. + ArcTan(Sqrt(eps)*(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2)))) -5*(-1 + x)*x*(((-2*Sqrt(eps))/(1 + eps*Power(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2),2)) - (8*Power(eps,1.5)*(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2))*Power(-0.5 + y,2))/Power(1 + eps*Power(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2),2),2))*(-1 + y)*y - (4*Sqrt(eps)*(-0.5 + y)*(-1 + 2*y))/(1 + eps*Power(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2),2)) +
//                                                                                                                                                                                                                                                                                                                                                      2*(Pi/2. + ArcTan(Sqrt(eps)*(0.0625 - Power(-0.5 + x,2) - Power(-0.5 + y,2)))));//(-1.)*((-8*(8*(1. - x) - 8*x)*(-0.5 + x)*(1. - y)*y*sqrt(epsilon))/ (1 + pow(0.0625 + 2*(0.0625 - pow(-0.5 + x,2) - pow(-0.5 + y,2))*sqrt(epsilon),2)) +8*(1. - x)*x*(1. - y)*y*((-4*sqrt(epsilon))/(1 + pow(0.0625 + 2*(0.0625 - pow(-0.5 + x,2) - pow(-0.5 + y,2))*sqrt(epsilon),2)) -(32*pow(-0.5 + x,2)*(0.0625 + 2*(0.0625 - pow(-0.5 + x,2) - pow(-0.5 + y,2))*sqrt(epsilon))*epsilon)/ pow(1 + pow(0.0625 + 2*(0.0625 - pow(-0.5 + x,2) - pow(-0.5 + y,2))*sqrt(epsilon),2),2)) + 8*(1. - x)*x*((-8*(1. - 2*y)*(-0.5 + y)*sqrt(epsilon))/ (1 + pow(0.0625 + 2*(0.0625 - pow(-0.5 + x,2) - pow(-0.5 + y,2))*sqrt(epsilon),2)) + (1. - y)*y*((-4*sqrt(epsilon))/(1 + pow(0.0625 + 2*(0.0625 - pow(-0.5 + x,2) - pow(-0.5 + y,2))*sqrt(epsilon),2)) - (32*pow(-0.5 + y,2)*(0.0625 + 2*(0.0625 - pow(-0.5 + x,2) - pow(-0.5 + y,2))*sqrt(epsilon))*epsilon)/pow(1 + pow(0.0625 + 2*(0.0625 - pow(-0.5 + x,2) - pow(-0.5 + y,2))*sqrt(epsilon),2),2)) -  2*atan(0.0625 + 2*(0.0625 - pow(-0.5 + x,2) - pow(-0.5 + y,2))*sqrt(epsilon))) -  16*(1. - y)*y*atan(0.0625 + 2*(0.0625 - pow(-0.5 + x,2) - pow(-0.5 + y,2))*sqrt(epsilon)));
    /*
     REAL F = 2*sqrt(epsilon);
     REAL arc = F*((0.25*0.25) - (x - 0.5)*(x - 0.5) - (y - 0.5)*(y - 0.5));
     REAL prodx = x*(x-1.);
     REAL prody = y*(y-1.);
     REAL prod = prodx*prody;
     REAL aux=8*prod*sqrt(epsilon)*(2*x-1.)*(x-0.5);
     REAL frac = aux/(1+arc*arc);
     
     
     
     //uxx
     REAL aux2x= (0.25*0.25)-(x-0.5)*(x-0.5)-(y-0.5)*(y-0.5);
     REAL fracaux=(1.+(4.*aux2x*aux2x)*epsilon);
     REAL termo2=fracaux*fracaux;
     REAL fator1=prody*(MyPi+2.*atan(arc));
     REAL fator2=(4./fracaux)*prody*sqrt(epsilon)*((2*x-1)*(2*x-1)+x*(x-1));
     REAL fator3=(64./termo2)*(prod*(x-0.5)*(x-0.5)*aux2x*(pow(epsilon, 1.5)));
     REAL aux2xx= (16./MyPi)*(fator1-fator2-fator3);
     
     //uyy
     
     REAL fator1y=prodx*(MyPi+2.*atan(arc));
     REAL fator2y=(4./fracaux)*prodx*sqrt(epsilon)*((2*y-1)*(2*y-1)+y*(y-1));
     REAL fator3y=(64./termo2)*(prod*(y-0.5)*(y-0.5)*aux2x*(pow(epsilon,1.5)));
     REAL aux2yy= (16./MyPi)*(fator1y-fator2y-fator3y);
     
     //divergente uxx+uyy
     
     disp[0]=(-1.)*(aux2xx+aux2yy);
     */
    
    
    
}

//Dados para o problema do juan
void ForcingMista(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
	double x = pt[0];
    double y = pt[1];
    disp[0]= (1.- 4.*MyPi*MyPi)/(2.*MyPi)*(cos(MyPi*x)*exp(y/2.));
}
void NeumannAbaixo(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
    double x = pt[0];
    //double y = pt[1];
    disp[0] = -cos(MyPi*x)/MyPi;
}

void NeumannAcima(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
    double x = pt[0];
    //double y = pt[1];
    disp[0] = 0.4087179842429694*cos(MyPi*x);
}

void DirichletEsquerda(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
//    TPZFMatrix<STATE> flux;
//    SolExataMista(pt, disp, flux);
    double y=pt[1];
    disp[0]=-(2./MyPi)*exp(y/2);
}

void SolExataMista(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux){
    
    solp.Resize(1, 0.);
    flux.Resize(3, 1.);
    flux(0,0)=flux(1,0)=flux(2,0)=0.;
    double x = pt[0];
    double y = pt[1];
    
    solp[0] = (-2./MyPi)*cos(MyPi*x)*exp(y/2.);
    flux(0,0)= -2*sin(MyPi*x)*exp(y/2.);
    flux(1,0)= (1./MyPi)*cos(MyPi*x)*exp(y/2.);
    flux(2,0)= (1.- 4.*MyPi*MyPi)/(2.*MyPi)*(cos(MyPi*x)*exp(y/2.));
}





//



void CC1(const TPZVec<REAL> &pt, TPZVec<STATE> &f) {
//    REAL x=pt[0];
    //    REAL prodx=x*(x-1.);
    //    REAL arc=2.*sqrt(epsilon)*((-0.1875)-(x-0.5)*(x-0.5));
    //    f[0]= 8.*prodx*(1. + (2./MyPi)*atan(arc));
    
    f[0]=0.;//8.*(1. - x)*x*atan(0.0625 + 2*sqrt(epsilon)*(-0.1875 - pow(-0.5 + x,2)));
    
    
    
}

void CC2(const TPZVec<REAL> &pt, TPZVec<STATE> &f) {
    
//    REAL y=pt[1];
    f[0]=0.;//8.*(1. - y)*y*atan(0.0625 + 2*sqrt(epsilon)*(-0.1875 - pow(-0.5 + y,2)));
    //    REAL prody=y*(y-1.);
    //    REAL arc=2.*sqrt(epsilon)*((-0.1875)-(y-0.5)*(y-0.5));
    //    f[0]= 8.*prody*(1 + (2./MyPi)*atan(arc));
}
void CC3(const TPZVec<REAL> &pt, TPZVec<STATE> &f) {
    double x=pt[0];
    
    //double fator=-4.-x*x;
    f[0]=2.+4.*x;//x*exp(fator);//0.;//2.*exp(x)*(1. - pow(x,2.));	//0.;//
    //f[0]=7.*x+7.;
}
void CC4(const TPZVec<REAL> &pt, TPZVec<STATE> &f) {
    
    double y=pt[0];
    //double fator=-4.-y*y;
    f[0]=(-1.)*(3.+4.*y);//0.;//-2.*exp(fator);//-MyPi*cos(MyPi*y);//2.*exp(x)*(1. - pow(x,2.));	//0.;
    //f[0]=2.-2.*y;
}


TPZCompMesh *CompMeshPAdap(TPZGeoMesh &gmesh,int porder,bool prefine){
    
    
    TPZCompMesh *comp = new TPZCompMesh(&gmesh);
    
    comp->SetDefaultOrder(porder);
    comp->SetDimModel(2);
    // Criar e inserir os materiais na malha
    TPZMatPoisson3d *mat = new TPZMatPoisson3d(matId,2);
    mat-> SetTrueFShapeHdiv();
    TPZMaterial * automat(mat);
    comp->InsertMaterialObject(automat);
    
    
    //funcao do lado direito da equacao do problema
    TPZAutoPointer<TPZFunction<STATE> > forcef;
    TPZDummyFunction<STATE> *dum = new TPZDummyFunction<STATE>(Forcing1, 5);
    dum->SetPolynomialOrder(20);
    forcef = dum;
    mat->SetForcingFunction(forcef);
    
    
    TPZAutoPointer<TPZFunction<STATE> > exata1 = new TPZDummyFunction<STATE>(SolExata, 5);
    mat->SetForcingFunctionExact(exata1);
    
    
    
    
    
    
    ///Criar condicoes de contorno
    
    //
    //		TPZAutoPointer<TPZFunction<STATE> > fCC1 = new TPZDummyFunction<STATE>(CC1);
    //        TPZAutoPointer<TPZFunction<STATE> > fCC2 = new TPZDummyFunction<STATE>(CC2);
    //		TPZAutoPointer<TPZFunction<STATE> > fCC3 = new TPZDummyFunction<STATE>(CC1);
    //		TPZAutoPointer<TPZFunction<STATE> > fCC4 = new TPZDummyFunction<STATE>(CC2);
    
    
    
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    TPZMaterial *bnd1 = mat->CreateBC (automat,-1,dirichlet,val1,val2);
    TPZMaterial *bnd2 = mat->CreateBC (automat,-2,dirichlet,val1,val2);
    TPZMaterial *bnd3 = mat->CreateBC (automat,-3,dirichlet,val1,val2);
    TPZMaterial *bnd4 = mat->CreateBC (automat,-4,dirichlet,val1,val2);
    TPZMaterial *bnd5 = mat->CreateBC (automat,-5,0,val1,val2);
    TPZMaterial *bnd6 = mat->CreateBC (automat,-6,0,val1,val2);
    TPZMaterial *bnd7 = mat->CreateBC (automat,-7,0,val1,val2);
    TPZMaterial *bnd8 = mat->CreateBC (automat,-8,0,val1,val2);
    
    //
    //            bnd->SetForcingFunction(fCC1);
    //            bnd2->SetForcingFunction(fCC2);
    //			bnd3->SetForcingFunction(fCC3);
    //            bnd4->SetForcingFunction(fCC4);
    
	
    ///Inserir condicoes de contorno
    comp->InsertMaterialObject(bnd1);
    comp->InsertMaterialObject(bnd2);
    comp->InsertMaterialObject(bnd3);
    comp->InsertMaterialObject(bnd4);
    comp->InsertMaterialObject(bnd5);
    comp->InsertMaterialObject(bnd6);
    comp->InsertMaterialObject(bnd7);
    comp->InsertMaterialObject(bnd8);
    
    
    // comp->SetAllCreateFunctionsHDivFull();
    comp->SetAllCreateFunctionsHDivPressure();
	//	comp->SetAllCreateFunctionsContinuous();
    
    
    
    // Ajuste da estrutura de dados computacional
    comp->AutoBuild();
    comp->AdjustBoundaryElements();//ajusta as condicoes de contorno
    comp->CleanUpUnconnectedNodes();//deleta os nos que nao tem elemntos conectados
    
    comp->SetName("Malha Computacional com ordem inicializada");
#ifdef LOG4CXX
    {
        std::stringstream sout;
        comp->Print(sout);
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
	/*
     
     int nel = comp->NElements();
     int iel;
     for(iel=0; iel<nel; iel++){
     
     TPZInterpolatedElement *intel;
     TPZCompEl *cel = comp->ElementVec()[iel];
     
     if (prefine) {
     
     
     intel = dynamic_cast<TPZInterpolatedElement *>(cel);
     if(intel){
     
     int fator=iel%2;
     
     if (cel->Dimension()==2 && fator==0) {
     
     intel->PRefine(porder+1);
     
     }
     
     
     if(cel->Dimension()==2 && fator!=0) {
     
     intel->PRefine(porder);
     
     
     }
     
     
     }
     }
     
     
     }
     
     
     comp->LoadReferences();
     comp->ExpandSolution();
     comp->AdjustBoundaryElements();
     */
    comp->SetName("Malha Computacional");
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        gmesh.Print(sout);
        comp->Print(sout);
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    
    
    
    return comp;
    
}

TPZCompMesh *CompMeshPAdapJuan(TPZGeoMesh &gmesh,int porder,bool prefine){
    
    
    TPZCompMesh *comp = new TPZCompMesh(&gmesh);
    
    comp->SetDefaultOrder(porder);
    comp->SetDimModel(2);
    // Criar e inserir os materiais na malha
    TPZMatPoisson3d *mat = new TPZMatPoisson3d(1,2);
    mat-> SetTrueFShapeHdiv();
    TPZMaterial * automat(mat);
    comp->InsertMaterialObject(automat);
    
    
    //funcao do lado direito da equacao do problema
    TPZAutoPointer<TPZFunction<STATE> > forcef;
    TPZDummyFunction<STATE> *dum = new TPZDummyFunction<STATE>(ForcingMista, 5);
    dum->SetPolynomialOrder(20);
    forcef = dum;
    mat->SetForcingFunction(forcef);
    
    
    TPZAutoPointer<TPZFunction<STATE> > exata1 = new TPZDummyFunction<STATE>(SolExataMista, 5);
    mat->SetForcingFunctionExact(exata1);
    
    
    
    
    
    
    ///Criar condicoes de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    TPZMaterial * BCond0 = automat->CreateBC(mat, bc0,neumann, val1, val2);
    TPZMaterial * BCond1 = automat->CreateBC(mat, bc1,neumann, val1, val2);
    TPZMaterial * BCond2 = automat->CreateBC(mat, bc2,neumann, val1, val2);
    TPZMaterial * BCond3 = automat->CreateBC(mat, bc3,dirichlet, val1, val2);
    
    TPZAutoPointer<TPZFunction<STATE> > bcmatNeumannAcima;
    bcmatNeumannAcima = new TPZDummyFunction<STATE>(NeumannAcima, 5);
    BCond2->SetForcingFunction(bcmatNeumannAcima);
    
    TPZAutoPointer<TPZFunction<STATE> > bcmatNeumannAbaixo;
    bcmatNeumannAbaixo = new TPZDummyFunction<STATE>(NeumannAbaixo, 5);
    BCond0->SetForcingFunction(bcmatNeumannAbaixo);
    
    TPZAutoPointer<TPZFunction<STATE> > bcmatDirichletEsquerda;
    bcmatDirichletEsquerda = new TPZDummyFunction<STATE>(DirichletEsquerda, 5);
    BCond3->SetForcingFunction(bcmatDirichletEsquerda);
    
    comp->InsertMaterialObject(BCond0);
    comp->InsertMaterialObject(BCond1);
    comp->InsertMaterialObject(BCond2);
    comp->InsertMaterialObject(BCond3);

    
    //
    
        comp->SetAllCreateFunctionsHDivPressure();

    
    
    
    // Ajuste da estrutura de dados computacional
    comp->AutoBuild();
    comp->AdjustBoundaryElements();//ajusta as condicoes de contorno
    comp->CleanUpUnconnectedNodes();//deleta os nos que nao tem elemntos conectados
    
    comp->SetName("Malha Computacional com ordem inicializada");
#ifdef LOG4CXX
    {
        std::stringstream sout;
        comp->Print(sout);
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
	/*
     
     int nel = comp->NElements();
     int iel;
     for(iel=0; iel<nel; iel++){
     
     TPZInterpolatedElement *intel;
     TPZCompEl *cel = comp->ElementVec()[iel];
     
     if (prefine) {
     
     
     intel = dynamic_cast<TPZInterpolatedElement *>(cel);
     if(intel){
     
     int fator=iel%2;
     
     if (cel->Dimension()==2 && fator==0) {
     
     intel->PRefine(porder+1);
     
     }
     
     
     if(cel->Dimension()==2 && fator!=0) {
     
     intel->PRefine(porder);
     
     
     }
     
     
     }
     }
     
     
     }
     
     
     comp->LoadReferences();
     comp->ExpandSolution();
     comp->AdjustBoundaryElements();
     */
    comp->SetName("Malha Computacional com ordens diferentes");
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        gmesh.Print(sout);
        comp->Print(sout);
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    
    
    
    return comp;
    
}
TPZCompMeshReferred *CreateCompMesh2d(TPZGeoMesh &gmesh,int porder){
    TPZCompEl::SetgOrder(porder);
    TPZCompMeshReferred *comp = new TPZCompMeshReferred(&gmesh);
    
    // Criar e inserir os materiais na malha
    TPZMatPoisson3d *mat = new TPZMatPoisson3d(1,2);
    TPZMaterial * automat(mat);
    comp->InsertMaterialObject(automat);
    
    TPZAutoPointer<TPZFunction<STATE> > force1 = new TPZDummyFunction<STATE>(Forcing1, 5);
    mat->SetForcingFunction(force1);
    TPZAutoPointer<TPZFunction<STATE> > exata1 = new TPZDummyFunction<STATE>(SolExata3, 5);
    mat->SetForcingFunctionExact(exata1);
    ///Inserir condicoes de contorno
    
    TPZFMatrix<STATE> val1(1,1,0.),val2(1,1,0.);
    TPZFMatrix<STATE> val11(1,1,0.), val22(1,1,0.);
    TPZMaterial *bnd = automat->CreateBC (automat,-1,1,val1,val2);//1
    TPZMaterial *bnd2 = automat->CreateBC (automat,-2,1,val1,val2);
    TPZMaterial *bnd3 = automat->CreateBC (automat,-3,1,val1,val2);//1
    TPZMaterial *bnd4 = automat->CreateBC (automat,-4,1,val1,val2);
    
    TPZAutoPointer<TPZFunction<STATE> > fCC1 = new TPZDummyFunction<STATE>(CC1, 5);
    //TPZAutoPointer<TPZFunction<STATE> > fCC2 = new TPZDummyFunction<STATE>(CC2, 5);
    bnd->SetForcingFunction(fCC1);
    bnd2->SetForcingFunction(fCC1);
    bnd3->SetForcingFunction(fCC1);
    bnd4->SetForcingFunction(fCC1);
    
    
    
    comp->InsertMaterialObject(bnd);
    comp->InsertMaterialObject(bnd2);
    comp->InsertMaterialObject(bnd3);
    comp->InsertMaterialObject(bnd4);
    //espaco de aproximacao
    
    comp->SetAllCreateFunctionsHDivPressure();
    
    
    // Ajuste da estrutura de dados computacional
    comp->AutoBuild();
    comp->AdjustBoundaryElements();//ajusta as condicoes de contorno
    comp->CleanUpUnconnectedNodes();//deleta os nos que nao tem elemntos conectados
    comp->SetName("Malha Computacional Original");
    
#ifdef LOG4CXX
    {
        std::stringstream sout;
        comp->Print(sout);
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    
    
    
    return comp;
    
}


TPZGeoMesh * MalhaGeoT(const int h,bool hrefine){//malha triangulo
    
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    //Criar ns
    const int nnode = 4;//AQUI
    const int nelem = 2;
    TPZGeoEl *elvec[nelem];
    const int dim = 2;//AQUI
    
    REAL co[nnode][dim] ={{0.,0.},{1.,0.},{1.,1.},{0.,1.}};//{{-1.,-1},{1.,-1},{1.,1.},{-1.,1.}};// {{0.,0.},{2.,0},{2.,2.},{0.,2.}};//{{-2.,-2},{2.,-2},{2.,2.},{-2.,2.}};//
    int64_t indices[2][nnode];//como serao enumerados os nos
    
    
    //el 1
    indices[0][0] = 0;
    indices[0][1] = 1;
    indices[0][2] = 3;
    //el2
    indices[1][0] = 2;
    indices[1][1] = 3;
    indices[1][2] = 1;
    
    
    int nod;
    TPZVec<REAL> coord(dim);
    for(nod=0; nod<nnode; nod++) {
        int64_t nodind = gmesh->NodeVec().AllocateNewElement();
        
        for(int d = 0; d < dim; d++)
        {
            coord[d] = co[nod][d];
        }
        gmesh->NodeVec()[nod].Initialize(nodind,coord,*gmesh);
    }
    //Criacao de elementos
    
    
    TPZVec<int64_t> nodind1(3);
    TPZVec<int64_t> nodind2(3);
    for(int i=0; i<3; i++){
        nodind1[i] = indices[0][i];
        nodind2[i] = indices[1][i];
    }
    
    int64_t index;
    elvec[0] = gmesh->CreateGeoElement(ETriangle,nodind1,1,index); //AQUI
    elvec[1] = gmesh->CreateGeoElement(ETriangle,nodind2,1,index); //AQUI
    
    
    gmesh->BuildConnectivity();
    
    
    //Cria as condicoes de contorno
    TPZGeoElBC gbc1(elvec[0],3,bc0);// condicao de fronteira tipo -1:
    TPZGeoElBC gbc2(elvec[0],5,bc3);// condicao de fronteira tipo -2:
    
    TPZGeoElBC gbc3(elvec[1],3,bc1);// condicao de fronteira tipo -3:
    TPZGeoElBC gbc4(elvec[1],5,bc2);// condicao de fronteira tipo -4:
    
    const std::string nameref;
    
    TPZAutoPointer<TPZRefPattern> ref;
    
    /*refinamento nao uniforme
    if (hrefine) {
        NoUniformRefine(gmesh, h);
    }
    
    
    else UniformRefine(gmesh, h);
    
    {
     ///Refinamento uniforme
     for ( int ref = 0; ref < h; ref++ )
     {// h indica o numero de refinamentos
     TPZVec<TPZGeoEl *> filhos;
     int64_t n = gmesh->NElements();
     for ( int64_t i = 0; i < n; i++ )
     {
     TPZGeoEl * gel = gmesh->ElementVec() [i];
     //if ( gel->Dimension() == 2 ) gel->Divide ( filhos );
     if(!gel->HasSubElement())
     {
     gel->Divide(filhos);
     }
     }//for i
     }
     
     
     
     }*/
    
    
    
#ifdef LOG4CXX
    {
        std::stringstream sout;
        gmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str().c_str());
    }
#endif
    return gmesh;
}
TPZGeoMesh *GMesh(bool ftriang, REAL Lx, REAL Ly){
    
    int Qnodes = 4;
	
	TPZGeoMesh * gmesh = new TPZGeoMesh;
	gmesh->SetMaxNodeId(Qnodes-1);
	gmesh->NodeVec().Resize(Qnodes);
	TPZVec<TPZGeoNode> Node(Qnodes);
	
	TPZVec <int64_t> TopolQuad(4);
    TPZVec <int64_t> TopolTriang(3);
	TPZVec <int64_t> TopolLine(2);
    TPZVec <int64_t> TopolPoint(1);
	
	//indice dos nos
	int64_t id = 0;
	REAL valx;
	for(int xi = 0; xi < Qnodes/2; xi++)
	{
		valx = xi*Lx;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx );//coord X
		Node[id].SetCoord(1 ,0. );//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
	
	for(int xi = 0; xi < Qnodes/2; xi++)
	{
		valx = Lx - xi*Lx;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx );//coord X
		Node[id].SetCoord(1 ,Ly);//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
	
	//indice dos elementos
	id = 0;
    
    if(ftriang==true)
    {
        TopolTriang[0] = 0;
        TopolTriang[1] = 1;
        TopolTriang[2] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (id,TopolTriang,matId,*gmesh);
        id++;
        
        TopolTriang[0] = 2;
        TopolTriang[1] = 1;
        TopolTriang[2] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,matId,*gmesh);
        id++;
        
        TopolLine[0] = 0;
        TopolLine[1] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc0,*gmesh);
        id++;
        
        TopolLine[0] = 2;
        TopolLine[1] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc1,*gmesh);
        id++;
        
        TopolLine[0] = 3;
        TopolLine[1] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
        id++;
        
        TopolLine[0] = 3;
        TopolLine[1] = 0;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc3,*gmesh);
    }
    else{
        TopolQuad[0] = 0;
        TopolQuad[1] = 1;
        TopolQuad[2] = 2;
        TopolQuad[3] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matId,*gmesh);
        id++;
        
        TopolLine[0] = 0;
        TopolLine[1] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc0,*gmesh);
        id++;
        
        TopolLine[0] = 1;
        TopolLine[1] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc1,*gmesh);
        id++;
        
        TopolLine[0] = 2;
        TopolLine[1] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
        id++;
        
        TopolLine[0] = 3;
        TopolLine[1] = 0;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc3,*gmesh);
    }
    
	gmesh->BuildConnectivity();
    
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout<<"\n\n Malha Geometrica Inicial\n ";
        gmesh->Print(sout);
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    
	return gmesh;
}






TPZGeoMesh * MalhaGeo/*QUADRILATEROS*/ ( const int h, bool nouniform)
{
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    REAL co[4][2] = {{0.,0.},{1.,0.},{1.,1.},{0.,1.}};//{{-1.,-1},{1.,-1},{1.,1.},{-1.,1.}};//
    int indices[1][4] = {{0,1,2,3}};
    
    int nnode = 4;
    const int nelem = 1;
    TPZGeoEl *elvec[nelem];
    int nod;
    for ( nod=0; nod<nnode; nod++ )
    {
        int64_t nodind = gmesh->NodeVec().AllocateNewElement();
        TPZVec<REAL> coord ( 2 );
        coord[0] = co[nod][0];
        coord[1] = co[nod][1];
        gmesh->NodeVec() [nodind].Initialize ( nod,coord,*gmesh );
    }
    
    int el;
    for ( el=0; el<nelem; el++ )
    {
        TPZVec<int64_t> nodind ( 4 );
        for ( nod=0; nod<4; nod++ ) nodind[nod]=indices[el][nod];
        int64_t index;
        elvec[el] = gmesh->CreateGeoElement ( EQuadrilateral,nodind,1,index );
    }
    
    gmesh->BuildConnectivity();
    
    TPZGeoElBC gbc1 ( elvec[0],4,-1 );// condicao de fronteira tipo -1: (x,y=0)
    TPZGeoElBC gbc2 ( elvec[0],5,-2 );// condicao de fronteira tipo -2: (x=1,y)
    TPZGeoElBC gbc3 ( elvec[0],6,-3 );// condicao de fronteira tipo -3: (x,y=1)
    TPZGeoElBC gbc4 ( elvec[0],7,-4 );// condicao de fronteira tipo -4: (x=0,y)
    
  /*  if (nouniform) {
        NoUniformRefine(gmesh, h);
    }
    else UniformRefine(gmesh, h);*/
    
/////Refinamento uniforme
//for ( int ref = 0; ref < h; ref++ )
//{// h indica o numero de refinamentos
//TPZVec<TPZGeoEl *> filhos;
//int64_t n = gmesh->NElements();
//for ( int64_t i = 0; i < n; i++ )
//{
//TPZGeoEl * gel = gmesh->ElementVec() [i];
////if ( gel->Dimension() == 2 ) gel->Divide ( filhos );
//if(!gel->HasSubElement())
//{
//gel->Divide(filhos);
//}
//}//for i
//}
//
////        //refinamento nao uniforme
////		if (hrefine) {
////            NoUniformRefine(gmesh, h);
////        }
////
////    else UniformRefine(gmesh, h);
    
    
#ifdef LOG4CXX
    {
        std::stringstream sout;
        gmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str().c_str());
    }
#endif
    
    return gmesh;
}


TPZGeoMesh * MalhaGeo2(const int h){//malha quadrilatero com 4 elementos
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    const int nelem=4;
    TPZGeoEl *elvec[nelem];
    //Criar ns
    const int nnode = 9;//AQUI
    const int dim = 2;//AQUI
    
    REAL co[nnode][dim] = {{0.,0.},{1.,0.},{1.,1.},{0.,1.},{0.5,0},{1.,0.5},{0.5,1.},{0.,0.5},{0.5,0.5}};//{{-1.,-1},{1.,-1},{1.,1.},{-1.,1.},{0.,-1.},{0.,1.}};
    
    
    int nodindAll[4][4]={{0,4,8,7},{4,1,5,8},{8,5,2,6},{7,8,6,3}};//como serao enumerados os nos
    
    
    int nod;
    TPZVec<REAL> coord(dim);
    for(nod=0; nod<nnode; nod++) {
        int64_t nodind = gmesh->NodeVec().AllocateNewElement();
        
        for(int d = 0; d < dim; d++)
        {
            coord[d] = co[nod][d];
        }
        gmesh->NodeVec()[nod].Initialize(nodind,coord,*gmesh);
    }
    /*Criacao de elementos
     int matId=40;
     int64_t id=0;
     TPZVec <int64_t> TopolQuad(4);
     TPZVec <int64_t> TopolLine(2);
     //-----
     
     TopolQuad[0] = 0;
     TopolQuad[1] = 4;
     TopolQuad[2] = 8;
     TopolQuad[3] = 7;
     new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matId,*gmesh);
     id++;
     matId++;
     
     TopolQuad[0] = 4;
     TopolQuad[1] = 1;
     TopolQuad[2] = 5;
     TopolQuad[3] = 8;
     new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matId,*gmesh);
     id++;
     matId++;
     
     TopolQuad[0] = 8;
     TopolQuad[1] = 5;
     TopolQuad[2] = 2;
     TopolQuad[3] = 6;
     new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matId,*gmesh);
     id++;
     matId++;
     
     TopolQuad[0] = 7;
     TopolQuad[1] = 8;
     TopolQuad[2] = 6;
     TopolQuad[3] = 3;
     new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matId,*gmesh);
     id++;
     matId++;
     
     
     TopolLine[0] = 0;
     TopolLine[1] = 4;
     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,-1,*gmesh);
     id++;
     
     TopolLine[0] = 4;
     TopolLine[1] = 1;
     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,-2,*gmesh);
     id++;
     
     TopolLine[0] = 1;
     TopolLine[1] = 5;
     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,-3,*gmesh);
     id++;
     
     TopolLine[0] = 5;
     TopolLine[1] = 2;
     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,-4,*gmesh);
     id++;
     
     TopolLine[0] = 2;
     TopolLine[1] = 6;
     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,-5,*gmesh);
     id++;
     
     TopolLine[0] = 6;
     TopolLine[1] = 3;
     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,-6,*gmesh);
     id++;
     
     TopolLine[0] = 3;
     TopolLine[1] = 7;
     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,-7,*gmesh);
     id++;
     
     TopolLine[0] = 7;
     TopolLine[1] = 0;
     new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,-8,*gmesh);
     
     */
    //Criacao de elementos
    //int matId=10;
    int64_t index;
	for ( int el=0; el<nelem; el++ )
	{
		TPZVec<int64_t> nodind(4);
		nodind[0]=nodindAll[el][0];
		nodind[1]=nodindAll[el][1];
		nodind[2]=nodindAll[el][2];
		nodind[3]=nodindAll[el][3];
		
		elvec[el] = gmesh->CreateGeoElement (EQuadrilateral,nodind,matId,index );
        
        //matId++;
        index++;
        
        
	}
    
    
    
    gmesh->BuildConnectivity();
    
    //Cria as condicoes de contorno
	TPZGeoElBC gbc1(elvec[0],4,-1);
    TPZGeoElBC gbc2(elvec[0],7,-2);
	TPZGeoElBC gbc3(elvec[1],4,-3);
	TPZGeoElBC gbc4(elvec[1],5,-4);
	TPZGeoElBC gbc5(elvec[2],5,-5);
	TPZGeoElBC gbc6(elvec[2],6,-6);
	TPZGeoElBC gbc7(elvec[3],6,-7);
    TPZGeoElBC gbc8(elvec[3],7,-8);
    
    
    
	/*	Refinamento uniforme
    for(int ref = 0; ref < h; ref++){// h indica o numero de refinamentos
        TPZVec<TPZGeoEl *> filhos;
        int64_t n = gmesh->NElements();
        for(int64_t i = 0; i < n; i++){
            TPZGeoEl * gel = gmesh->ElementVec()[i];
            if(!gel->HasSubElement())
            {
                gel->Divide(filhos);
            }
        }}
    
    if (hrefine) {
        //refinamento diferencialvel
        {
            
            TPZVec<TPZGeoEl *> filhos;
            int64_t n = gmesh->NElements();
            
            
            for(int64_t i = 0; i < n; i++){
                TPZGeoEl * gel = gmesh->ElementVec()[i];
                if(!gel->HasSubElement() && gel->Dimension()==2 && i%2==0)
                    
                {
                    gel->Divide(filhos);
                }
            }
        }
        
        //refinamento 1D--irei refinar tambem os elementos 1D
        
        {
            
            TPZVec<TPZGeoEl *> filhos;
            int64_t n = gmesh->NElements();
            
            
            for(int64_t i = 0; i < n; i++){
                TPZGeoEl * gel = gmesh->ElementVec()[i];
                if (gel->Dimension()!=1) {
                    continue;
                }
                TPZGeoElSide Elside=gel->Neighbour(2);
                TPZGeoEl *NeighEl=Elside.Element();
                if (NeighEl->HasSubElement()) {
                    gel->Divide(filhos);
                }
                
                
                
            }
        }
    }
    
    
    */
    
#ifdef LOG4CXX
    {
        std::stringstream sout;
        gmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str().c_str());
    }
#endif
    return gmesh;
    
}
void SolGraf(TPZCompMesh *malha, std::ofstream &GraficoSol){
    const int64_t nelem = malha->NElements();
    
    //TPZFMatrix<REAL> sol(4,1,0.);
    TPZManVector<REAL,4> solP(1);
    TPZManVector<REAL,4> solF(3,0.);
    
    TPZFMatrix<REAL> axes(3,1,0.);
    TPZFMatrix<REAL> phi;
    TPZFMatrix<REAL> dphix;
    
    //TPZManVector<REAL> sol(1);
    TPZSolVec sol;
    sol[0].Resize(1);
    //TPZFMatrix<REAL> dsol(3,1,0.);
    TPZGradSolVec dsol;
    dsol[0].Redim(3,1);
    
    ///Percorrer todos elementos
    for(int64_t el=0; el < nelem; el++){
        TPZCompEl * Cel = malha->ElementVec()[el];
        
        if(!Cel) continue;
        
        
        TPZGeoEl *gel = Cel->Reference();
        if(	malha->ElementVec()[el]->Dimension()== 2) {
            TPZIntPoints *intr = gel-> CreateSideIntegrationRule(gel->NSides() - 1,2);//side and order
            
            for(int in=0; in < intr->NPoints(); in++)
            {
                REAL peso;
                TPZVec<REAL> pto(2),pto2(1);
                intr->Point(in, pto, peso);
                
                TPZFMatrix<REAL> jac(2,2),invjac(2,2),axes(3,3);
                REAL jacdet;
                gel->Jacobian(pto,jac,axes,jacdet,invjac);
                
                //	malha->LoadSolution(sol);
                
                
                TPZManVector< REAL,3 > xco(3);
                TPZManVector<STATE> p(1);
                TPZFMatrix<STATE> fluxo(3,0);
                TPZManVector<REAL,4> solF(3,0.);
                
                gel->X(pto,xco);
                SolExata(xco, p, fluxo);
                Cel->ComputeSolution(xco,sol,dsol,axes);
                GraficoSol<<"{ "<<xco[0]<< " , "<< xco[1]<< ","<<sol[0]<<"},"<<std::endl;
                
                
            }
            delete intr;
            
        }
        else continue;
        
        
    }
    
    
}

void SolveLU ( TPZAnalysis &an ){
    TPZCompMesh *malha = an.Mesh();
    time_t tempoinicial;
    time_t tempofinal;
    time(&tempoinicial);
    //TPZFrontStructMatrix<TPZFrontNonSym> mat ( malha );// não funciona com método iterativo
    TPZSkylineStructMatrix mat(malha);
    //	TPZFStructMatrix mat( malha );
    //	TPZSpStructMatrix mat( malha );
    TPZStepSolver<STATE> solv;
    
    //solv.SetDirect ( ELU );
    //solv.SetDirect(ECholesky);
    solv.SetDirect(ELDLt);
    
    std::cout << "ELU " << std::endl;
    an.SetSolver ( solv );
    an.SetStructuralMatrix ( mat );
    std::cout << std::endl;
    an.Solution().Redim ( 0,0 );
    std::cout << "Assemble " << std::endl;
    time(&tempofinal);
    std::cout << "Antes do assemble " << tempofinal - tempoinicial << std::endl;
    
    time(&tempoinicial);
    an.Assemble();
    time(&tempofinal);
    std::cout << "  Tempo do assemble " << tempofinal - tempoinicial << std::endl;
    
    time(&tempoinicial);
    an.Solve();
    time(&tempofinal);
    std::cout << "  Tempo do assemble " << tempofinal - tempoinicial << std::endl;
    
    std::cout << std::endl;
    std::cout << "  No equacoes = " << malha->NEquations() << std::endl;
    // cout << "Banda = " << malha->BandWidth() << endl;
}

void ChangeP(TPZCompMesh * cmesh, TPZCompEl * cel, int newP)
{
    //cmesh->RemoveCompElfromPorderContainer(cel);
    //alterar P
    TPZInterpolationSpace * sp = dynamic_cast<TPZInterpolatedElement*>(cel);
    if(sp)
    {
        sp->PRefine(newP);
    }
    else
    {
        //FUDEU!!!
        DebugStop();
    }
    //cmesh->AddCompElfromPorderContainer(cel);
}

TPZGeoMesh * GeoMeshGrid( int h){
    // Rectangular geometric mesh using TPZGenGrid
    
    TPZVec < int > refin(2);
    TPZVec < REAL > corx0(3);
    TPZVec < REAL > corx1(3);
    //int  	numlayer = 1; // Layers Numbers
    //REAL  	rotation = 0.5; // For testing purpose
    int numlayer=1;
    // refinement level
    refin[0] = 1;
    refin[1] = 1;
    //	x0	lower left coordinate
    corx0[0] = 0.0;
    corx0[1] = 0.0;
    corx0[2] = 0.0;
    //	x1	upper right coordinate
    corx1[0] =1.0;// 2.0;//1.0;
    corx1[1] = 1.0;//2.0;//1.0;
    corx1[2] = 0.0;
    TPZGenGrid geomesh(refin,corx0,corx1,numlayer);
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    geomesh.Read(gmesh);
    
    // Setting BC conditions
    geomesh.SetBC(gmesh,4,0);
    geomesh.SetBC(gmesh,5,0);
    geomesh.SetBC(gmesh,6,0);
    geomesh.SetBC(gmesh,7,0);
    
    // refinamento
    for ( int ref = 0; ref < h; ref++ )
    {// h indica o numero de refinamentos
        TPZVec<TPZGeoEl *> filhos;
        int64_t n = gmesh->NElements();
        for ( int64_t i = 0; i < n; i++ )
        {
            TPZGeoEl * gel = gmesh->ElementVec() [i];
            //if ( gel->Dimension() == 2 ) gel->Divide ( filhos );
            if(!gel->HasSubElement())
            {
                gel->Divide(filhos);
            }
        }//for i
    }//ref
    
    
    
    
#ifdef LOG4CXX
    {
        std::stringstream sout;
        gmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str().c_str());
    }
#endif
    return gmesh;
}

void RefineGeoElements(int dim,TPZGeoMesh *gmesh,TPZVec<REAL> &point,REAL r,REAL &distance,bool &isdefined) {
	TPZManVector<REAL> centerpsi(3), center(3);
	// Refinamento de elementos selecionados
	TPZGeoEl *gel;
	TPZVec<TPZGeoEl *> sub;
	
	int nelem = 0;
	int ngelem=gmesh->NElements();
	// na esquina inferior esquerda Nó = (0,-1,0)
	while(nelem<ngelem) {
		gel = gmesh->ElementVec()[nelem++];
		if(gel->Dimension()!=dim || gel->HasSubElement()) continue;
		gel->CenterPoint(gel->NSides()-1,centerpsi);
		gel->X(centerpsi,center);
		if(!isdefined) {
			TPZVec<REAL> FirstNode(3,0.);
			gel->CenterPoint(0,centerpsi);
			gel->X(centerpsi,FirstNode);
			REAL distancia = TPZGeoEl::Distance(center,FirstNode);
			if(distancia > distance) distance = distancia;
			isdefined = true;
		}
		REAL centerdist = TPZGeoEl::Distance(center,point);
		if(fabs(r-centerdist) < distance) {
			gel->Divide(sub);
		}
	}
}
void GetPointsOnCircunference(int npoints,TPZVec<REAL> &center,REAL radius,TPZVec<TPZManVector<REAL,3> > &Points) {
	Points.Resize(npoints);
	TPZManVector<REAL,3> point(3,0.);
	REAL angle = (2*M_PI)/npoints;
	for(int i=0;i<npoints;i++) {
		point[0] = center[0]+radius*cos(i*angle);
		point[1] = center[1]+radius*sin(i*angle);
		Points[i] = point;
	}
}

void UniformRefine(TPZGeoMesh* gmesh, int nDiv)
{
    /*   TPZVec<TPZGeoEl *> filhos;
     int64_t n = gmesh->NElements();
     for(int64_t i = 0; i < n; i++){
     TPZGeoEl * gel = gmesh->ElementVec()[i];
     
     if(!gel->HasSubElement() )
     {
     gel->Divide(filhos);
     }
     }
     */
    
    
    ///Refinamento uniforme
    for ( int ref = 0; ref < nDiv; ref++ )
    {// h indica o numero de refinamentos
        TPZVec<TPZGeoEl *> filhos;
        int64_t n = gmesh->NElements();
        for ( int64_t i = 0; i < n; i++ )
        {
            TPZGeoEl * gel = gmesh->ElementVec() [i];
            //if ( gel->Dimension() == 2 ) gel->Divide ( filhos );
            if(!gel->HasSubElement())
            {
                gel->Divide(filhos);
            }
        }//for i
    }
    
    gmesh->ResetConnectivities();
    gmesh->BuildConnectivity();
    
    
}


void NoUniformRefine(TPZGeoMesh* gmesh, int nDiv)
{
    
    //refinamento nao uniforme
    
    
    TPZVec<TPZGeoEl *> filhos;
    int64_t n = gmesh->NElements();
    
    
    for(int64_t i = 0; i < n; i++){
        TPZGeoEl * gel = gmesh->ElementVec()[i];
        if(!gel->HasSubElement() && gel->Dimension()==2 && i%2==0)
            
        {
            gel->Divide(filhos);
        }
    }
    
    //refinamento 1D--irei refinar tambem os elementos 1D
    
    
    for(int64_t i = 0; i < n; i++){
        TPZGeoEl * gel = gmesh->ElementVec()[i];
        if (gel->Dimension()!=1) {
            continue;
        }
        TPZGeoElSide Elside=gel->Neighbour(2);
        TPZGeoEl *NeighEl=Elside.Element();
        if (NeighEl->HasSubElement()) {
            gel->Divide(filhos);
        }
        
        
        
    }
    
    //    gmesh->ResetConnectivities();
    //	gmesh->BuildConnectivity();
    //
    //    return gmesh;
}

void RefiningNearCircunference(int dim,TPZGeoMesh *gmesh,int nref,int ntyperefs) {
    
	int i;
	bool isdefined = false;
	
	// Refinando no local desejado
	int npoints = 1000;
	TPZVec<REAL> point(3);
	point[0] = point[1] = 0.5; point[2] = 0.0;
	REAL r = 0.25;
	TPZVec<TPZManVector<REAL,3> > Points(npoints);
	GetPointsOnCircunference(npoints,point,r,Points);
	
	if(ntyperefs==2) {
		REAL radius = 0.19;
		for(i=0;i<nref;i+=2) {
			// To refine elements with center near to points than radius
			RefineGeoElements(dim,gmesh,point,r,radius,isdefined);
			RefineGeoElements(dim,gmesh,point,r,radius,isdefined);
			if(i < 5) radius *= 0.35;
			else if(i < 7) radius *= 0.2;
			else radius *= 0.1;
		}
		if(i==nref) {
			RefineGeoElements(dim,gmesh,point,r,radius,isdefined);
		}
	}
	else {
		REAL radius = 0.2;
		for(i=0;i<nref+1;i++) {
			// To refine elements with center near to points than radius
			RefineGeoElements(dim,gmesh,point,r,radius,isdefined);
			if(i < 8) radius *= 0.6;
			else if(i < 7) radius *= 0.3;
			else radius *= 0.15;
		}
	}
	// Constructing connectivities
	//gmesh->ResetConnectivities();
    RegularizeMesh(gmesh);
	gmesh->BuildConnectivity();
}

void RegularizeMesh(TPZGeoMesh *gmesh)
{
    bool changed = true;
    while (changed)
    {
        changed = false;
        int nel = gmesh->NElements();
        for (int64_t el=0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->ElementVec()[el];
            if (gel->HasSubElement()) {
                continue;
            }
            int dim = gel->Dimension();
            if (dim != 2) {
                continue;
            }
            int nsides = gel->NSides();
            int nrefined = 0;
            int nsidedim = 0;
            for (int is=0; is<nsides; is++) {
                int sidedim = gel->SideDimension(is);
                if (sidedim != dim-1) {
                    continue;
                }
                nsidedim++;
            }
            for (int is=0; is<nsides; is++) {
                int sidedim = gel->SideDimension(is);
                if (sidedim != dim-1) {
                    continue;
                }
                TPZGeoElSide thisside(gel,is);
                TPZGeoElSide neighbour = thisside.Neighbour();
                if (neighbour != thisside) {
                    TPZStack<TPZGeoElSide> subelements;
                    neighbour.GetSubElements2(subelements);
                    int nsub = subelements.size();
                    if (nsub > 0) {
                        nrefined++;
                    }
                    for (int isub=0; isub<nsub; isub++) {
                        TPZGeoElSide sub = subelements[isub];
                        if (sub.Dimension() != dim-1) {
                            continue;
                        }
                        if (sub.HasSubElement()) {
                            TPZManVector<TPZGeoEl *> newsub;
                            gel->Divide(newsub);
                            changed = true;
                            break;
                        }
                    }
                }
                if (gel->HasSubElement()) {
                    break;
                }
            }
            if (nrefined >= nsidedim-1) {
                TPZManVector<TPZGeoEl *> newsub;
                gel->Divide(newsub);
                changed = true;
            }
        }
    }
}


void ErrorHDiv(TPZCompMesh *hdivmesh, std::ostream &out)
{
    int64_t nel = hdivmesh->NElements();
    int dim = hdivmesh->Dimension();
    TPZStack<REAL> vech;
    TPZManVector<REAL,10> globerrors(10,0.);
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = hdivmesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        //TPZMaterial *material= hdivmesh-> MaterialVec()[el];
        TPZMaterialData data;
        data.fNeedsHSize=true;
        TPZInterpolationSpace *sp = dynamic_cast<TPZInterpolationSpace * >(cel);
        TPZVec<REAL> qsi(2,0.);
        sp->ComputeRequiredData(data, qsi);
        REAL &hsize = data.HSize;
        vech.push_back(hsize);
        
        TPZGeoEl *gel = cel->Reference();
        
        if (!gel || gel->Dimension() != dim) {
            continue;
        }
        TPZManVector<REAL,10> elerror(10,0.);
        cel->EvaluateError(SolExata, elerror, false);
        
        int nerr = elerror.size();
        for (int i=0; i<nerr; i++) {
            globerrors[i] += elerror[i]*elerror[i];
        }
        
    }
    
    int nh = vech.size();
    REAL hmax=0;
    for(int i=0; i<nh; i++){
        if(vech[i]>hmax){
            hmax=vech[i];
        }
    }
    // out << "Errors associated with HDiv space\n";
    out << "Hmax = "    << hmax << std::endl;
    out << "L2 Norm for pressure = "    << sqrt(globerrors[0]) << std::endl;
    out << "L2 Norm for flux = "    << sqrt(globerrors[1]) << std::endl;
    out << "L2 Norm for divergence = "    << sqrt(globerrors[2])  <<std::endl;
    out << "Hdiv Norm for flux = "    << sqrt(globerrors[3])  <<std::endl;
    
}

void ErrorL2(TPZCompMesh *l2mesh, std::ostream &out)
{
    int64_t nel = l2mesh->NElements();
    int dim = l2mesh->Dimension();
    TPZManVector<REAL,10> globerrors(10,0.);
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = l2mesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel || gel->Dimension() != dim) {
            continue;
        }
        TPZManVector<REAL,10> elerror(10,0.);
        cel->EvaluateError(SolExataMista, elerror, false);
        int nerr = elerror.size();
        globerrors.resize(nerr);
#ifdef LOG4CXX
        if (logger->isDebugEnabled()) {
            std::stringstream sout;
            sout << "L2 Error sq of element " << el << elerror[0]*elerror[0];
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        for (int i=0; i<nerr; i++) {
            globerrors[i] += elerror[i]*elerror[i];
        }
        
    }
    out << "Errors associated with L2 space\n";
    out << "L2 Norm for pressure = "    << sqrt(globerrors[0]) << std::endl;
}

void SetDifferentOrderP(TPZCompMesh *comp,int porder){
    int nel = comp->NElements();
    int iel;
    for(iel=0; iel<nel; iel++){
        
        TPZInterpolatedElement *intel;
        TPZCompEl *cel = comp->ElementVec()[iel];
        
               
            
            intel = dynamic_cast<TPZInterpolatedElement *>(cel);
            if(intel){
                
                int fator=iel%2;
                
                if (cel->Dimension()==2 && fator==0) {
                    
                    intel->PRefine(porder+1);
                    
                }
                
                
                if(cel->Dimension()==2 && fator!=0) {
                    
                    intel->PRefine(porder);
                    
                    
                }
                
                
            }
        }
        
    
    
    comp->LoadReferences();
    comp->ExpandSolution();
    comp->AdjustBoundaryElements();
    
    comp->SetName("Malha Computacional cOm diferentes Ordens");
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        comp->Print(sout);
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif

}

void SetDifferentOrderPMesh4Elem(TPZCompMesh *comp,int porder){
    int nel = comp->NElements();
    int iel;
    for(iel=0; iel<nel; iel++)
    {
        
        TPZInterpolatedElement *intel;
        TPZCompEl *cel = comp->ElementVec()[iel];
        if(!cel) continue;
        
        TPZGeoEl *gel = cel->Reference();
        
        if(gel->Dimension()==2)
        {
            //pegar o elemento geometrico mais grosso que contenha gel
            TPZGeoEl *fgel = gel->LowestFather();
            
            if(fgel->Index()==0)
            {
                intel = dynamic_cast<TPZInterpolatedElement *>(cel);
                 intel->PRefine(porder);
            }
            
            if(fgel->Index()==1)
            {
                intel = dynamic_cast<TPZInterpolatedElement *>(cel);
                intel->PRefine(porder+1);
            }
            
            if(fgel->Index()==2)
            {
                intel = dynamic_cast<TPZInterpolatedElement *>(cel);
                intel->PRefine(porder+2);
            }
            
            if(fgel->Index()==3)
            {
                intel = dynamic_cast<TPZInterpolatedElement *>(cel);
                intel->PRefine(porder+3);
            }
        }
        
    }
    
    comp->LoadReferences();
    comp->ExpandSolution();
    comp->AdjustBoundaryElements();
    
    comp->SetName("Malha Computacional cOm diferentes Ordens");
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        comp->Print(sout);
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    
}



