/*******************************************
 *   Made by P. Cesar de A. Lucci          *
 *   LabMeC - 2007                         *
 *******************************************/
#include "tpzcurvedtriangle.h"

#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzgeotriangle.h"
#include "pznoderep.h"
#include "pzshapetriang.h"

#include <iostream>
#include <cmath>

using namespace std;
using namespace pzshape;
using namespace pzgeom;
using namespace pztopology;

TPZCurvedTriangle::~TPZCurvedTriangle(){
}

void TPZCurvedTriangle::X(TPZFMatrix &coord, TPZVec<REAL>& par, TPZVec< REAL >& result)
{
     double qsi = par[0]; double eta = par[1];
     double qsiTil, etaTil;
     if(fabs(qsi) < 0.001) qsi = 0.;
     if(fabs(eta) < 0.001) eta = 0.;
     if(qsi > 0.999 && qsi < 1.001) qsi = 1.;
     if(eta > 0.999 && eta < 1.001) eta = 1.;
     if(eta > qsi)
     {
          double tempMap = qsi;
          qsi = eta;
          eta = tempMap;
     }

     if(qsi == 1. || eta == 1.)
     {
          qsiTil = qsi;
          etaTil = eta;
     }

     else
     {
          qsiTil = (-2.*eta - qsi*eta*(qsi + eta) + eta*(3*qsi + eta) + (qsi + eta)*sqrt(pow(1. - qsi,2)*(2. + pow(qsi,2) + pow(eta,2) - 2.*(qsi + eta))))/sqrt(pow(1. - qsi,2)*(2 + pow(qsi,2) + pow(eta,2) - 2.*(qsi + eta)));

          etaTil = ((2. - qsi - eta)*eta*(-1. - qsi*eta + (qsi + eta) + sqrt(pow(1. - qsi,2)*(2. + pow(qsi,2) + pow(eta,2) - 2.*(qsi + eta)))))/((1. - qsi)*sqrt(pow(1. - qsi,2)*(2. + pow(qsi,2) + pow(eta,2) - 2.*(qsi + eta))));
     }

     qsi = par[0]; eta = par[1];
     if(eta > qsi)
     {
          double tempMap = qsiTil;
          qsiTil = etaTil;
          etaTil = tempMap;
     }
     TPZVec<REAL> parMap(2);
     parMap[0] = qsiTil; parMap[1] = etaTil;
     TPZGeoTriangle::X(coord,parMap,result);
}

void TPZCurvedTriangle::Jacobian(TPZFMatrix & coord, TPZVec<REAL>& par, TPZFMatrix &jacobian, TPZFMatrix &axes, REAL &detjac, TPZFMatrix &jacinv)
{
     jacobian.Resize(2,2); axes.Resize(2,3); jacinv.Resize(2,2);
     double qsiTil = par[0];
     double etaTil = par[1];
     double J00, J01, J10, J11;

     if(etaTil > qsiTil)
     {
          double tempEtaMap = qsiTil;
          qsiTil = etaTil;
          etaTil = tempEtaMap;
     }

     if(qsiTil >= .9999 && qsiTil <= 1.0001 && etaTil >= -0.0001 && etaTil <= 0.0001)
     {
          J00 = 1.; J01 = 0.; J10 = 0.; J11 = 0.;
     }

     else
     {
          J00 = 1. + (pow(1. - qsiTil,3)*(1. - etaTil)*(qsiTil - etaTil)*etaTil)/pow(pow(1. - qsiTil,2)*(2. + pow(qsiTil,2) + pow(etaTil,2) - 2.*(qsiTil + etaTil)),1.5);
          J01 = 1. - (pow(1. - qsiTil,3)*(4. - pow(qsiTil,3) - 2.*pow(qsiTil,2)*etaTil - pow(etaTil,3) - 6.*(qsiTil + etaTil) + (4.*pow(qsiTil,2) + 5.*qsiTil*etaTil +         3.*pow(etaTil,2))))/pow(pow(1. - qsiTil,2)*(2. + pow(qsiTil,2) + pow(etaTil,2) - 2.*(qsiTil + etaTil)),1.5);
          J10 = ((1. - etaTil)*etaTil*(1. - qsiTil + (sqrt(pow(1. - qsiTil,2)*(2. + pow(qsiTil,2) + pow(etaTil,2) - 2.*(qsiTil + etaTil)))*(-4. + pow(qsiTil,3) + 2.*pow(qsiTil,2)*etaTil + pow(etaTil,3) + (7.*qsiTil + 5.*etaTil) - (5.*pow(qsiTil,2) + 4.*qsiTil*etaTil + 3.*pow(etaTil,2))))/(pow(2. + pow(qsiTil,2) + pow(etaTil,2) - 2.*(qsiTil + etaTil),2))))/pow(1. - qsiTil,3);
          J11 = ((1. - qsiTil)*(2. - qsiTil - 2.*etaTil) - (sqrt(pow(1. - qsiTil,2)*(2. + pow(qsiTil,2) + pow(etaTil,2) - 2.*(qsiTil + etaTil)))*(4. + (qsiTil + 3.*etaTil)*(4.*qsiTil + 5.*etaTil) - 2.*(3.*qsiTil + 7.*etaTil) + etaTil*(2.*pow(qsiTil,3) + 3.*pow(qsiTil,2)*etaTil + qsiTil*pow(etaTil,2) + 2.*pow(etaTil,3)) - (pow(qsiTil,3) + 10.*pow(qsiTil,2)*etaTil + 9.*qsiTil*pow(etaTil,2) + 8.*pow(etaTil,3))))/(pow(2. + pow(qsiTil,2) + pow(etaTil,2) - 2.*(qsiTil + etaTil),2)))/pow(1. - qsiTil,2);
     }

     qsiTil = par[0]; etaTil = par[1];
     if(etaTil > qsiTil)
     {
          double tempGrad = J00;
          J00 = J11; J11 = tempGrad;
          tempGrad = J01;
          J01 = J10; J10 = tempGrad;
     }
     jacobian(0,0) = J00; jacobian(0,1) = J01;
     jacobian(1,0) = J10; jacobian(1,1) = J11;

     TPZFMatrix jacobianPZ(2,2);
     TPZVec< REAL > out(4);
     X(coord,par,out);

     TPZGeoTriangle::Jacobian(coord, out, jacobianPZ, axes, detjac, jacinv);
     double jac1, jac2, jac3, jac4;
     jac1 = jacobianPZ(0,0)*jacobian(0,0) + jacobianPZ(0,1)*jacobian(1,0);
     jac2 = jacobianPZ(0,0)*jacobian(0,1) + jacobianPZ(0,1)*jacobian(1,1);
     jac3 = jacobianPZ(1,0)*jacobian(0,0) + jacobianPZ(1,1)*jacobian(1,0);
     jac4 = jacobianPZ(1,0)*jacobian(0,1) + jacobianPZ(1,1)*jacobian(1,1);
     jacobian(0,0) = jac1; jacobian(0,1) = jac2;
     jacobian(1,0) = jac3; jacobian(1,1) = jac4;

     detjac = (jacobian(0,0)*jacobian(1,1)) - (jacobian(0,1)*jacobian(1,0));
     if( fabs(detjac) < 1.e-8  )
     {
          jacinv(0,0) =  1.; jacinv(0,1) =  0.;
          jacinv(1,0) =  0.; jacinv(1,1) =  1.;
     }

     else
     {
          jacinv(0,0) =  jacobian(1,1) / detjac;
          jacinv(0,1) = -jacobian(0,1) / detjac;
          jacinv(1,0) = -jacobian(1,0) / detjac;
          jacinv(1,1) =  jacobian(0,0) / detjac;
     }
}

TPZGeoEl *TPZCurvedTriangle::CreateBCGeoEl(TPZGeoEl *orig,int side,int bc)
{
     if(side==6)
     {
          TPZManVector<int> nodes(3);
          int i;
          for (i=0;i<3;i++) nodes[i] = orig->SideNodeIndex(side,i);
          int index;
          TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(ETriangle,nodes,bc,index); int iside;
          for (iside = 0; iside <3; iside++)
          {
               TPZGeoElSide(gel,iside).SetConnectivity(TPZGeoElSide(orig,TPZShapeTriang::SideConnectLocId(side,iside)));
          }
          TPZGeoElSide(gel,6).SetConnectivity(TPZGeoElSide(orig,side));
          return gel;
     }

     else if(side>-1 && side<3)
     {
          TPZManVector<int> nodeindexes(1);
          nodeindexes[0] = orig->SideNodeIndex(side,0); int index;
          TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EPoint,nodeindexes,bc,index);
          TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,side));
          return gel;
     }

     else if(side > 2 && side < 3)
     {
          TPZManVector<int> nodes(2);
          nodes[0] = orig->SideNodeIndex(side,0);
          nodes[1] = orig->SideNodeIndex(side,1); int index;
          TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EOned,nodes,bc,index);
          TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,TPZShapeTriang::SideConnectLocId(side,0)));
          TPZGeoElSide(gel,1).SetConnectivity(TPZGeoElSide(orig,TPZShapeTriang::SideConnectLocId(side,1)));
          TPZGeoElSide(gel,2).SetConnectivity(TPZGeoElSide(orig,side));
          return gel;
     }
     else PZError << "TPZGeoTriangle::CreateBCGeoEl has no bc.\n";
     return 0;
}
