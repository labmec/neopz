#include "includes.h"
#include "pzbndcond.h"

using namespace std;
static TPZCompMesh *Create3DExpMesh();
static void Exact3DExp(TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix &dsol);
static void ForcingFunction3DExp (TPZVec<REAL> &x,TPZVec<REAL> &sol);

//*************************************
//************Option 15****************
//*****Discontinuous 3D Cube Mesh******
//*************************************
TPZCompMesh *Create3DExpMesh() {
  REAL co[8][3] = {
    {0.,0.,0.},
    {1.,0.,0.},
    {1.,1.,0.},
    {0.,1.,0.},
    {0.,0.,1.},
    {1.,0.,1.},
    {1.,1.,1.},
    {0.,1.,1.}
  };
  int indices[1][8] = {{0,1,2,3,4,5,6,7}};
  const int nelem = 1;
  int nnode = 8;
  TPZGeoEl *elvec[nelem];
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  int nod;
  for(nod=0; nod<nnode; nod++) {
    int nodind = gmesh->NodeVec().AllocateNewElement();
    TPZVec<REAL> coord(3);
    coord[0] = co[nod][0];
    coord[1] = co[nod][1];
    coord[2] = co[nod][2];
    TPZGeoNode pznode (nod,coord,*gmesh);
    gmesh->NodeVec()[nodind] = pznode;
  }
  int el;
  for(el=0; el<nelem; el++) {
    TPZVec<int> nodind(8);
    for(nod=0; nod<8; nod++) nodind[nod]=indices[el][nod];
    int index;
    elvec[el] = gmesh->CreateGeoElement(ECube,nodind,1,index);
  }
  gmesh->BuildConnectivity2();

  // bc -1 -> Dirichlet at all faces
  TPZGeoElBC gbc1(elvec[0],20,-1,*gmesh);
  TPZGeoElBC gbc3(elvec[0],21,-2,*gmesh);
  TPZGeoElBC gbc4(elvec[0],22,-3,*gmesh);
  TPZGeoElBC gbc5(elvec[0],23,-4,*gmesh);
  TPZGeoElBC gbc6(elvec[0],24,-5,*gmesh);
  TPZGeoElBC gbc2(elvec[0],25,-6,*gmesh);

  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);

  TPZAutoPointer<TPZMaterial> mat;
  TPZMaterialTest3D *mat3 = new TPZMaterialTest3D(1);
  TPZFMatrix mp (3,1,1.);
  TPZMaterialTest3D::geq3=1;
  mat3->SetMaterial(mp);
  mat3->SetForcingFunction(ForcingFunction3DExp);
  mat = mat3;
  cmesh->InsertMaterialObject(mat);

  TPZFMatrix val1(1,1,0.),val2(1,1,0.);
  TPZAutoPointer<TPZMaterial> bc[6];
  int i;
  for (i=0;i<6;i++) {
    bc[i] = mat->CreateBC(mat,-(i+1),0,val1,val2);
    cmesh->InsertMaterialObject(bc[i]);
  }
  cmesh->AutoBuild();
  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();
  return cmesh;
}

void ForcingFunction3DExp(TPZVec<REAL> &pt, TPZVec<REAL> &sol){
if (pt.NElements()<3 || sol.NElements()<1){
	cout << "Chabu Aqui...";
 }
  REAL x = pt[0];
  REAL y = pt[1];
  REAL z = pt[2];
 sol[0] = -6*(1 - x)*x*(1 - y)*y*(1 - z) + 6*(1 - x)*x*(1 - y)*y*z - 6*(1 - x)*x*(1 - y)*(1 - z)*z +
   6*(1 - x)*x*y*(1 - z)*z - 6*(1 - x)*(1 - y)*y*(1 - z)*z + 6*x*(1 - y)*y*(1 - z)*z +
   2*(1 - x)*x*(1 - y)*y*(x + y + z) - 2*(1 - x)*x*(1 - y)*(1 - z)*(x + y + z) +
   2*(1 - x)*x*y*(1 - z)*(x + y + z) - 2*(1 - x)*(1 - y)*y*(1 - z)*(x + y + z) +
   2*x*(1 - y)*y*(1 - z)*(x + y + z) + 2*(1 - x)*x*(1 - y)*z*(x + y + z) -
   2*(1 - x)*x*y*z*(x + y + z) + 2*(1 - x)*(1 - y)*y*z*(x + y + z) -
   2*x*(1 - y)*y*z*(x + y + z) + 2*(1 - x)*x*(1 - z)*z*(x + y + z) -
   2*(1 - x)*(1 - y)*(1 - z)*z*(x + y + z) + 2*x*(1 - y)*(1 - z)*z*(x + y + z) +
   2*(1 - x)*y*(1 - z)*z*(x + y + z) - 2*x*y*(1 - z)*z*(x + y + z) +
   2*(1 - y)*y*(1 - z)*z*(x + y + z) + (1 - x)*x*(1 - y)*y*(1 - z)*z*(x + y + z);


  /*
sol[0] = -((-1 + pow(x,2))*(-1 + pow(y,2))*z*(-1 + pow(z,2))*
      pow((4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4),
       1/(1. + 2*pow(z,2)))) - (pow((4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))/
        pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4),1/(1. + 2*pow(z,2)))*
      (-((pow(x,2)*(-1 + pow(y,2))*z*(-1 + pow(z,2))*
             (51.2 - 12.799999999999997*pow(z,2) + 64.*pow(z,4) + pow(x,2)*(-64. + 64.*pow(z,2)) +
               pow(y,2)*(-16. + 64.*pow(z,2))))/
           ((0.5 + pow(z,2))*(0.1 + pow(x,2) + pow(y,2) + pow(z,2))*
             (4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2)))) -
        (x*(-1 + pow(x,2))*y*z*(-1 + pow(z,2))*
           (49.6 + 35.2*pow(z,2) + 64.*pow(z,4) + pow(y,2)*(-16. + 64.*pow(z,2)) +
             pow(x,2)*(32. + 64.*pow(z,2))))/
         ((0.5 + pow(z,2))*(0.1 + pow(x,2) + pow(y,2) + pow(z,2))*
           (4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))) +
        (8.*pow(x,2)*(-1 + pow(y,2))*z*(-1 + pow(z,2))*
           (-0.8 + 24.8*pow(z,2) + 9.6*pow(z,4) + 16.*pow(z,6) + pow(x,4)*(-16. + 16.*pow(z,2)) +
             pow(y,4)*(8. + 16.*pow(z,2)) + pow(y,2)*(24.8 + 17.6*pow(z,2) + 32.*pow(z,4)) +
             pow(x,2)*(25.6 - 6.399999999999999*pow(z,2) + 32.*pow(z,4) +
                pow(y,2)*(-8. + 32.*pow(z,2)))))/
         ((0.5 + pow(z,2))*(0.1 + pow(x,2) + pow(y,2) + pow(z,2))*
           pow(4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2),2)) +
        (2*pow(x,2)*(-1 + pow(y,2))*z*(-1 + pow(z,2))*
           (-0.8 + 24.8*pow(z,2) + 9.6*pow(z,4) + 16.*pow(z,6) + pow(x,4)*(-16. + 16.*pow(z,2)) +
             pow(y,4)*(8. + 16.*pow(z,2)) + pow(y,2)*(24.8 + 17.6*pow(z,2) + 32.*pow(z,4)) +
             pow(x,2)*(25.6 - 6.399999999999999*pow(z,2) + 32.*pow(z,4) +
                pow(y,2)*(-8. + 32.*pow(z,2)))))/
         ((0.5 + pow(z,2))*pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),2)*
           (4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))) +
        ((-1 + pow(y,2))*z*(1 - pow(z,2))*
           (-0.8 + 24.8*pow(z,2) + 9.6*pow(z,4) + 16.*pow(z,6) + pow(x,4)*(-16. + 16.*pow(z,2)) +
             pow(y,4)*(8. + 16.*pow(z,2)) + pow(y,2)*(24.8 + 17.6*pow(z,2) + 32.*pow(z,4)) +
             pow(x,2)*(25.6 - 6.399999999999999*pow(z,2) + 32.*pow(z,4) +
                pow(y,2)*(-8. + 32.*pow(z,2)))))/
         ((0.5 + pow(z,2))*(0.1 + pow(x,2) + pow(y,2) + pow(z,2))*
           (4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))) +
        (8.*x*(-1 + pow(x,2))*y*z*(-1 + pow(z,2))*
           (-0.8 + 24.8*pow(z,2) + 9.6*pow(z,4) + 16.*pow(z,6) + pow(y,4)*(-16. + 16.*pow(z,2)) +
             pow(x,4)*(8. + 16.*pow(z,2)) +
             pow(y,2)*(25.6 - 6.399999999999999*pow(z,2) + 32.*pow(z,4)) +
             pow(x,2)*(24.8 + 17.6*pow(z,2) + 32.*pow(z,4) + pow(y,2)*(-8. + 32.*pow(z,2)))))/
         ((0.5 + pow(z,2))*(0.1 + pow(x,2) + pow(y,2) + pow(z,2))*
           pow(4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2),2)) +
        (2*x*(-1 + pow(x,2))*y*z*(-1 + pow(z,2))*
           (-0.8 + 24.8*pow(z,2) + 9.6*pow(z,4) + 16.*pow(z,6) + pow(y,4)*(-16. + 16.*pow(z,2)) +
             pow(x,4)*(8. + 16.*pow(z,2)) +
             pow(y,2)*(25.6 - 6.399999999999999*pow(z,2) + 32.*pow(z,4)) +
             pow(x,2)*(24.8 + 17.6*pow(z,2) + 32.*pow(z,4) + pow(y,2)*(-8. + 32.*pow(z,2)))))/
         ((0.5 + pow(z,2))*pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),2)*
           (4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))) -
        (2*x*y*z*(-1 + pow(z,2))*(-0.8 + 24.8*pow(z,2) + 9.6*pow(z,4) + 16.*pow(z,6) +
             pow(y,4)*(-16. + 16.*pow(z,2)) + pow(x,4)*(8. + 16.*pow(z,2)) +
             pow(y,2)*(25.6 - 6.399999999999999*pow(z,2) + 32.*pow(z,4)) +
             pow(x,2)*(24.8 + 17.6*pow(z,2) + 32.*pow(z,4) + pow(y,2)*(-8. + 32.*pow(z,2)))))/
         ((0.5 + pow(z,2))*(0.1 + pow(x,2) + pow(y,2) + pow(z,2))*
           (4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))) -
        (x*(-1 + pow(x,2))*(-1 + pow(y,2))*pow(z,2)*(-1 + pow(z,2))*
           (-0.32000000000000006 + 192.*pow(x,6) + 192.*pow(y,6) - 7.680000000000003*pow(z,2) +
             95.99999999999999*pow(z,4) + 384.*pow(z,6) + pow(y,4)*(108.8 + 768.*pow(z,2)) +
             pow(x,4)*(108.8 + 576.*pow(y,2) + 768.*pow(z,2)) +
             pow(y,2)*(-7.040000000000001 + 204.79999999999995*pow(z,2) + 960.*pow(z,4)) +
             pow(x,2)*(-7.040000000000001 + 576.*pow(y,4) + 204.79999999999995*pow(z,2) +
                960.*pow(z,4) + pow(y,2)*(217.6 + 1536.*pow(z,2)))))/
         (pow(0.5 + pow(z,2),2)*pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),2)*
           pow(4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2),2)) +
        4*x*(-1 + pow(y,2))*(1 - 3*pow(z,2) +
           (z*(1 - pow(z,2))*((z*(0.8 - 24.*pow(x,2) - 24.*pow(y,2) - 24.*pow(z,2))*(0.5 + pow(z,2)))/
                 ((0.1 + pow(x,2) + pow(y,2) + pow(z,2))*(4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2)))\
                 - 2*z*log((4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))/
                   pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4))))/(2.*pow(0.5 + pow(z,2),2))) +
        (x*(0.8 - 24.*pow(x,2) - 24.*pow(y,2) - 24.*pow(z,2))*
           ((x*(-1 + pow(y,2))*z*(1 - pow(z,2))*
                (-0.8 + 24.8*pow(z,2) + 9.6*pow(z,4) + 16.*pow(z,6) +
                  pow(x,4)*(-16. + 16.*pow(z,2)) + pow(y,4)*(8. + 16.*pow(z,2)) +
                  pow(y,2)*(24.8 + 17.6*pow(z,2) + 32.*pow(z,4)) +
                  pow(x,2)*(25.6 - 6.399999999999999*pow(z,2) + 32.*pow(z,4) +
                     pow(y,2)*(-8. + 32.*pow(z,2)))))/
              ((0.5 + pow(z,2))*(0.1 + pow(x,2) + pow(y,2) + pow(z,2))*
                (4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))) +
             ((-1 + pow(x,2))*y*z*(1 - pow(z,2))*
                (-0.8 + 24.8*pow(z,2) + 9.6*pow(z,4) + 16.*pow(z,6) +
                  pow(y,4)*(-16. + 16.*pow(z,2)) + pow(x,4)*(8. + 16.*pow(z,2)) +
                  pow(y,2)*(25.6 - 6.399999999999999*pow(z,2) + 32.*pow(z,4)) +
                  pow(x,2)*(24.8 + 17.6*pow(z,2) + 32.*pow(z,4) + pow(y,2)*(-8. + 32.*pow(z,2)))))/
              ((0.5 + pow(z,2))*(0.1 + pow(x,2) + pow(y,2) + pow(z,2))*
                (4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))) +
             2*(-1 + pow(x,2))*(-1 + pow(y,2))*
              (1 - 3*pow(z,2) + (z*(1 - pow(z,2))*
                   ((z*(0.8 - 24.*pow(x,2) - 24.*pow(y,2) - 24.*pow(z,2))*(0.5 + pow(z,2)))/
                      ((0.1 + pow(x,2) + pow(y,2) + pow(z,2))*
                        (4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))) -
                     2*z*log((4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))/
                        pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4))))/(2.*pow(0.5 + pow(z,2),2)))))
          /((0.1 + pow(x,2) + pow(y,2) + pow(z,2))*(1. + 2*pow(z,2))*
           (4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2)))))/2. -
   (pow((4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4),
       1/(1. + 2*pow(z,2)))*(-(((-1 + pow(x,2))*pow(y,2)*z*(-1 + pow(z,2))*
             (51.2 - 12.799999999999997*pow(z,2) + 64.*pow(z,4) + pow(y,2)*(-64. + 64.*pow(z,2)) +
               pow(x,2)*(-16. + 64.*pow(z,2))))/
           ((0.5 + pow(z,2))*(0.1 + pow(x,2) + pow(y,2) + pow(z,2))*
             (4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2)))) -
        (x*y*(-1 + pow(y,2))*z*(-1 + pow(z,2))*
           (49.6 + 35.2*pow(z,2) + 64.*pow(z,4) + pow(x,2)*(-16. + 64.*pow(z,2)) +
             pow(y,2)*(32. + 64.*pow(z,2))))/
         ((0.5 + pow(z,2))*(0.1 + pow(x,2) + pow(y,2) + pow(z,2))*
           (4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))) +
        (8.*x*y*(-1 + pow(y,2))*z*(-1 + pow(z,2))*
           (-0.8 + 24.8*pow(z,2) + 9.6*pow(z,4) + 16.*pow(z,6) + pow(x,4)*(-16. + 16.*pow(z,2)) +
             pow(y,4)*(8. + 16.*pow(z,2)) + pow(y,2)*(24.8 + 17.6*pow(z,2) + 32.*pow(z,4)) +
             pow(x,2)*(25.6 - 6.399999999999999*pow(z,2) + 32.*pow(z,4) +
                pow(y,2)*(-8. + 32.*pow(z,2)))))/
         ((0.5 + pow(z,2))*(0.1 + pow(x,2) + pow(y,2) + pow(z,2))*
           pow(4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2),2)) +
        (2*x*y*(-1 + pow(y,2))*z*(-1 + pow(z,2))*
           (-0.8 + 24.8*pow(z,2) + 9.6*pow(z,4) + 16.*pow(z,6) + pow(x,4)*(-16. + 16.*pow(z,2)) +
             pow(y,4)*(8. + 16.*pow(z,2)) + pow(y,2)*(24.8 + 17.6*pow(z,2) + 32.*pow(z,4)) +
             pow(x,2)*(25.6 - 6.399999999999999*pow(z,2) + 32.*pow(z,4) +
                pow(y,2)*(-8. + 32.*pow(z,2)))))/
         ((0.5 + pow(z,2))*pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),2)*
           (4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))) -
        (2*x*y*z*(-1 + pow(z,2))*(-0.8 + 24.8*pow(z,2) + 9.6*pow(z,4) + 16.*pow(z,6) +
             pow(x,4)*(-16. + 16.*pow(z,2)) + pow(y,4)*(8. + 16.*pow(z,2)) +
             pow(y,2)*(24.8 + 17.6*pow(z,2) + 32.*pow(z,4)) +
             pow(x,2)*(25.6 - 6.399999999999999*pow(z,2) + 32.*pow(z,4) +
                pow(y,2)*(-8. + 32.*pow(z,2)))))/
         ((0.5 + pow(z,2))*(0.1 + pow(x,2) + pow(y,2) + pow(z,2))*
           (4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))) +
        (8.*(-1 + pow(x,2))*pow(y,2)*z*(-1 + pow(z,2))*
           (-0.8 + 24.8*pow(z,2) + 9.6*pow(z,4) + 16.*pow(z,6) + pow(y,4)*(-16. + 16.*pow(z,2)) +
             pow(x,4)*(8. + 16.*pow(z,2)) +
             pow(y,2)*(25.6 - 6.399999999999999*pow(z,2) + 32.*pow(z,4)) +
             pow(x,2)*(24.8 + 17.6*pow(z,2) + 32.*pow(z,4) + pow(y,2)*(-8. + 32.*pow(z,2)))))/
         ((0.5 + pow(z,2))*(0.1 + pow(x,2) + pow(y,2) + pow(z,2))*
           pow(4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2),2)) +
        (2*(-1 + pow(x,2))*pow(y,2)*z*(-1 + pow(z,2))*
           (-0.8 + 24.8*pow(z,2) + 9.6*pow(z,4) + 16.*pow(z,6) + pow(y,4)*(-16. + 16.*pow(z,2)) +
             pow(x,4)*(8. + 16.*pow(z,2)) +
             pow(y,2)*(25.6 - 6.399999999999999*pow(z,2) + 32.*pow(z,4)) +
             pow(x,2)*(24.8 + 17.6*pow(z,2) + 32.*pow(z,4) + pow(y,2)*(-8. + 32.*pow(z,2)))))/
         ((0.5 + pow(z,2))*pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),2)*
           (4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))) +
        ((-1 + pow(x,2))*z*(1 - pow(z,2))*
           (-0.8 + 24.8*pow(z,2) + 9.6*pow(z,4) + 16.*pow(z,6) + pow(y,4)*(-16. + 16.*pow(z,2)) +
             pow(x,4)*(8. + 16.*pow(z,2)) +
             pow(y,2)*(25.6 - 6.399999999999999*pow(z,2) + 32.*pow(z,4)) +
             pow(x,2)*(24.8 + 17.6*pow(z,2) + 32.*pow(z,4) + pow(y,2)*(-8. + 32.*pow(z,2)))))/
         ((0.5 + pow(z,2))*(0.1 + pow(x,2) + pow(y,2) + pow(z,2))*
           (4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))) -
        ((-1 + pow(x,2))*y*(-1 + pow(y,2))*pow(z,2)*(-1 + pow(z,2))*
           (-0.32000000000000006 + 192.*pow(x,6) + 192.*pow(y,6) - 7.680000000000003*pow(z,2) +
             95.99999999999999*pow(z,4) + 384.*pow(z,6) + pow(y,4)*(108.8 + 768.*pow(z,2)) +
             pow(x,4)*(108.8 + 576.*pow(y,2) + 768.*pow(z,2)) +
             pow(y,2)*(-7.040000000000001 + 204.79999999999995*pow(z,2) + 960.*pow(z,4)) +
             pow(x,2)*(-7.040000000000001 + 576.*pow(y,4) + 204.79999999999995*pow(z,2) +
                960.*pow(z,4) + pow(y,2)*(217.6 + 1536.*pow(z,2)))))/
         (pow(0.5 + pow(z,2),2)*pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),2)*
           pow(4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2),2)) +
        4*(-1 + pow(x,2))*y*(1 - 3*pow(z,2) +
           (z*(1 - pow(z,2))*((z*(0.8 - 24.*pow(x,2) - 24.*pow(y,2) - 24.*pow(z,2))*(0.5 + pow(z,2)))/
                 ((0.1 + pow(x,2) + pow(y,2) + pow(z,2))*(4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2)))\
                 - 2*z*log((4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))/
                   pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4))))/(2.*pow(0.5 + pow(z,2),2))) +
        (y*(0.8 - 24.*pow(x,2) - 24.*pow(y,2) - 24.*pow(z,2))*
           ((x*(-1 + pow(y,2))*z*(1 - pow(z,2))*
                (-0.8 + 24.8*pow(z,2) + 9.6*pow(z,4) + 16.*pow(z,6) +
                  pow(x,4)*(-16. + 16.*pow(z,2)) + pow(y,4)*(8. + 16.*pow(z,2)) +
                  pow(y,2)*(24.8 + 17.6*pow(z,2) + 32.*pow(z,4)) +
                  pow(x,2)*(25.6 - 6.399999999999999*pow(z,2) + 32.*pow(z,4) +
                     pow(y,2)*(-8. + 32.*pow(z,2)))))/
              ((0.5 + pow(z,2))*(0.1 + pow(x,2) + pow(y,2) + pow(z,2))*
                (4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))) +
             ((-1 + pow(x,2))*y*z*(1 - pow(z,2))*
                (-0.8 + 24.8*pow(z,2) + 9.6*pow(z,4) + 16.*pow(z,6) +
                  pow(y,4)*(-16. + 16.*pow(z,2)) + pow(x,4)*(8. + 16.*pow(z,2)) +
                  pow(y,2)*(25.6 - 6.399999999999999*pow(z,2) + 32.*pow(z,4)) +
                  pow(x,2)*(24.8 + 17.6*pow(z,2) + 32.*pow(z,4) + pow(y,2)*(-8. + 32.*pow(z,2)))))/
              ((0.5 + pow(z,2))*(0.1 + pow(x,2) + pow(y,2) + pow(z,2))*
                (4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))) +
             2*(-1 + pow(x,2))*(-1 + pow(y,2))*
              (1 - 3*pow(z,2) + (z*(1 - pow(z,2))*
                   ((z*(0.8 - 24.*pow(x,2) - 24.*pow(y,2) - 24.*pow(z,2))*(0.5 + pow(z,2)))/
                      ((0.1 + pow(x,2) + pow(y,2) + pow(z,2))*
                        (4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))) -
                     2*z*log((4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))/
                        pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4))))/(2.*pow(0.5 + pow(z,2),2)))))
          /((0.1 + pow(x,2) + pow(y,2) + pow(z,2))*(1. + 2*pow(z,2))*
           (4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2)))))/2. -
   (pow((4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4),
       1/(1. + 2*pow(z,2)))*(-((x*(-1 + pow(y,2))*pow(z,2)*(-1 + pow(z,2))*
             (49.6 + 32.*pow(x,4) + 32.*pow(y,4) + 38.4*pow(z,2) + 96.*pow(z,4) +
               pow(y,2)*(35.2 + 128.*pow(z,2)) +
               pow(x,2)*(-12.799999999999997 + 64.*pow(y,2) + 128.*pow(z,2))))/
           ((0.5 + pow(z,2))*(0.1 + pow(x,2) + pow(y,2) + pow(z,2))*
             (4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2)))) -
        ((-1 + pow(x,2))*y*pow(z,2)*(-1 + pow(z,2))*
           (49.6 + 32.*pow(x,4) + 32.*pow(y,4) + 38.4*pow(z,2) + 96.*pow(z,4) +
             pow(y,2)*(-12.799999999999997 + 128.*pow(z,2)) +
             pow(x,2)*(35.2 + 64.*pow(y,2) + 128.*pow(z,2))))/
         ((0.5 + pow(z,2))*(0.1 + pow(x,2) + pow(y,2) + pow(z,2))*
           (4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))) +
        (8.*x*(-1 + pow(y,2))*pow(z,2)*(-1 + pow(z,2))*
           (-0.8 + 24.8*pow(z,2) + 9.6*pow(z,4) + 16.*pow(z,6) + pow(x,4)*(-16. + 16.*pow(z,2)) +
             pow(y,4)*(8. + 16.*pow(z,2)) + pow(y,2)*(24.8 + 17.6*pow(z,2) + 32.*pow(z,4)) +
             pow(x,2)*(25.6 - 6.399999999999999*pow(z,2) + 32.*pow(z,4) +
                pow(y,2)*(-8. + 32.*pow(z,2)))))/
         ((0.5 + pow(z,2))*(0.1 + pow(x,2) + pow(y,2) + pow(z,2))*
           pow(4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2),2)) +
        (2*x*(-1 + pow(y,2))*pow(z,2)*(-1 + pow(z,2))*
           (-0.8 + 24.8*pow(z,2) + 9.6*pow(z,4) + 16.*pow(z,6) + pow(x,4)*(-16. + 16.*pow(z,2)) +
             pow(y,4)*(8. + 16.*pow(z,2)) + pow(y,2)*(24.8 + 17.6*pow(z,2) + 32.*pow(z,4)) +
             pow(x,2)*(25.6 - 6.399999999999999*pow(z,2) + 32.*pow(z,4) +
                pow(y,2)*(-8. + 32.*pow(z,2)))))/
         ((0.5 + pow(z,2))*pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),2)*
           (4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))) +
        (2*x*(-1 + pow(y,2))*pow(z,2)*(-1 + pow(z,2))*
           (-0.8 + 24.8*pow(z,2) + 9.6*pow(z,4) + 16.*pow(z,6) + pow(x,4)*(-16. + 16.*pow(z,2)) +
             pow(y,4)*(8. + 16.*pow(z,2)) + pow(y,2)*(24.8 + 17.6*pow(z,2) + 32.*pow(z,4)) +
             pow(x,2)*(25.6 - 6.399999999999999*pow(z,2) + 32.*pow(z,4) +
                pow(y,2)*(-8. + 32.*pow(z,2)))))/
         (pow(0.5 + pow(z,2),2)*(0.1 + pow(x,2) + pow(y,2) + pow(z,2))*
           (4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))) -
        (2*x*(-1 + pow(y,2))*pow(z,2)*(-0.8 + 24.8*pow(z,2) + 9.6*pow(z,4) + 16.*pow(z,6) +
             pow(x,4)*(-16. + 16.*pow(z,2)) + pow(y,4)*(8. + 16.*pow(z,2)) +
             pow(y,2)*(24.8 + 17.6*pow(z,2) + 32.*pow(z,4)) +
             pow(x,2)*(25.6 - 6.399999999999999*pow(z,2) + 32.*pow(z,4) +
                pow(y,2)*(-8. + 32.*pow(z,2)))))/
         ((0.5 + pow(z,2))*(0.1 + pow(x,2) + pow(y,2) + pow(z,2))*
           (4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))) +
        (x*(-1 + pow(y,2))*(1 - pow(z,2))*
           (-0.8 + 24.8*pow(z,2) + 9.6*pow(z,4) + 16.*pow(z,6) + pow(x,4)*(-16. + 16.*pow(z,2)) +
             pow(y,4)*(8. + 16.*pow(z,2)) + pow(y,2)*(24.8 + 17.6*pow(z,2) + 32.*pow(z,4)) +
             pow(x,2)*(25.6 - 6.399999999999999*pow(z,2) + 32.*pow(z,4) +
                pow(y,2)*(-8. + 32.*pow(z,2)))))/
         ((0.5 + pow(z,2))*(0.1 + pow(x,2) + pow(y,2) + pow(z,2))*
           (4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))) +
        (8.*(-1 + pow(x,2))*y*pow(z,2)*(-1 + pow(z,2))*
           (-0.8 + 24.8*pow(z,2) + 9.6*pow(z,4) + 16.*pow(z,6) + pow(y,4)*(-16. + 16.*pow(z,2)) +
             pow(x,4)*(8. + 16.*pow(z,2)) +
             pow(y,2)*(25.6 - 6.399999999999999*pow(z,2) + 32.*pow(z,4)) +
             pow(x,2)*(24.8 + 17.6*pow(z,2) + 32.*pow(z,4) + pow(y,2)*(-8. + 32.*pow(z,2)))))/
         ((0.5 + pow(z,2))*(0.1 + pow(x,2) + pow(y,2) + pow(z,2))*
           pow(4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2),2)) +
        (2*(-1 + pow(x,2))*y*pow(z,2)*(-1 + pow(z,2))*
           (-0.8 + 24.8*pow(z,2) + 9.6*pow(z,4) + 16.*pow(z,6) + pow(y,4)*(-16. + 16.*pow(z,2)) +
             pow(x,4)*(8. + 16.*pow(z,2)) +
             pow(y,2)*(25.6 - 6.399999999999999*pow(z,2) + 32.*pow(z,4)) +
             pow(x,2)*(24.8 + 17.6*pow(z,2) + 32.*pow(z,4) + pow(y,2)*(-8. + 32.*pow(z,2)))))/
         ((0.5 + pow(z,2))*pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),2)*
           (4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))) +
        (2*(-1 + pow(x,2))*y*pow(z,2)*(-1 + pow(z,2))*
           (-0.8 + 24.8*pow(z,2) + 9.6*pow(z,4) + 16.*pow(z,6) + pow(y,4)*(-16. + 16.*pow(z,2)) +
             pow(x,4)*(8. + 16.*pow(z,2)) +
             pow(y,2)*(25.6 - 6.399999999999999*pow(z,2) + 32.*pow(z,4)) +
             pow(x,2)*(24.8 + 17.6*pow(z,2) + 32.*pow(z,4) + pow(y,2)*(-8. + 32.*pow(z,2)))))/
         (pow(0.5 + pow(z,2),2)*(0.1 + pow(x,2) + pow(y,2) + pow(z,2))*
           (4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))) -
        (2*(-1 + pow(x,2))*y*pow(z,2)*(-0.8 + 24.8*pow(z,2) + 9.6*pow(z,4) + 16.*pow(z,6) +
             pow(y,4)*(-16. + 16.*pow(z,2)) + pow(x,4)*(8. + 16.*pow(z,2)) +
             pow(y,2)*(25.6 - 6.399999999999999*pow(z,2) + 32.*pow(z,4)) +
             pow(x,2)*(24.8 + 17.6*pow(z,2) + 32.*pow(z,4) + pow(y,2)*(-8. + 32.*pow(z,2)))))/
         ((0.5 + pow(z,2))*(0.1 + pow(x,2) + pow(y,2) + pow(z,2))*
           (4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))) +
        ((-1 + pow(x,2))*y*(1 - pow(z,2))*
           (-0.8 + 24.8*pow(z,2) + 9.6*pow(z,4) + 16.*pow(z,6) + pow(y,4)*(-16. + 16.*pow(z,2)) +
             pow(x,4)*(8. + 16.*pow(z,2)) +
             pow(y,2)*(25.6 - 6.399999999999999*pow(z,2) + 32.*pow(z,4)) +
             pow(x,2)*(24.8 + 17.6*pow(z,2) + 32.*pow(z,4) + pow(y,2)*(-8. + 32.*pow(z,2)))))/
         ((0.5 + pow(z,2))*(0.1 + pow(x,2) + pow(y,2) + pow(z,2))*
           (4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))) +
        (-1 + pow(x,2))*(-1 + pow(y,2))*(-12*z +
           (z*(1 - pow(z,2))*((-8.*pow(z,2)*(0.8 - 24.*pow(x,2) - 24.*pow(y,2) - 24.*pow(z,2))*
                   (0.5 + pow(z,2)))/
                 ((0.1 + pow(x,2) + pow(y,2) + pow(z,2))*
                   pow(4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2),2)) -
                (2*pow(z,2)*(0.8 - 24.*pow(x,2) - 24.*pow(y,2) - 24.*pow(z,2))*(0.5 + pow(z,2)))/
                 (pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),2)*
                   (4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))) -
                (48.*pow(z,2)*(0.5 + pow(z,2)))/
                 ((0.1 + pow(x,2) + pow(y,2) + pow(z,2))*(4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2)))\
                 + ((0.8 - 24.*pow(x,2) - 24.*pow(y,2) - 24.*pow(z,2))*(0.5 + pow(z,2)))/
                 ((0.1 + pow(x,2) + pow(y,2) + pow(z,2))*(4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2)))\
                 - 2*log((4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))/
                   pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4))))/pow(0.5 + pow(z,2),2) -
           (4*pow(z,2)*(1 - pow(z,2))*((z*(0.8 - 24.*pow(x,2) - 24.*pow(y,2) - 24.*pow(z,2))*
                   (0.5 + pow(z,2)))/
                 ((0.1 + pow(x,2) + pow(y,2) + pow(z,2))*(4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2)))\
                 - 2*z*log((4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))/
                   pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4))))/pow(0.5 + pow(z,2),3) -
           (2*pow(z,2)*((z*(0.8 - 24.*pow(x,2) - 24.*pow(y,2) - 24.*pow(z,2))*(0.5 + pow(z,2)))/
                 ((0.1 + pow(x,2) + pow(y,2) + pow(z,2))*(4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2)))\
                 - 2*z*log((4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))/
                   pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4))))/pow(0.5 + pow(z,2),2) +
           ((1 - pow(z,2))*((z*(0.8 - 24.*pow(x,2) - 24.*pow(y,2) - 24.*pow(z,2))*(0.5 + pow(z,2)))/
                 ((0.1 + pow(x,2) + pow(y,2) + pow(z,2))*(4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2)))\
                 - 2*z*log((4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))/
                   pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4))))/pow(0.5 + pow(z,2),2)) +
        (((z*(0.8 - 24.*pow(x,2) - 24.*pow(y,2) - 24.*pow(z,2))*(1. + 2*pow(z,2)))/
              ((0.1 + pow(x,2) + pow(y,2) + pow(z,2))*(4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))) -
             4*z*log((4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))/
                pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4)))*
           ((x*(-1 + pow(y,2))*z*(1 - pow(z,2))*
                (-0.8 + 24.8*pow(z,2) + 9.6*pow(z,4) + 16.*pow(z,6) +
                  pow(x,4)*(-16. + 16.*pow(z,2)) + pow(y,4)*(8. + 16.*pow(z,2)) +
                  pow(y,2)*(24.8 + 17.6*pow(z,2) + 32.*pow(z,4)) +
                  pow(x,2)*(25.6 - 6.399999999999999*pow(z,2) + 32.*pow(z,4) +
                     pow(y,2)*(-8. + 32.*pow(z,2)))))/
              ((0.5 + pow(z,2))*(0.1 + pow(x,2) + pow(y,2) + pow(z,2))*
                (4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))) +
             ((-1 + pow(x,2))*y*z*(1 - pow(z,2))*
                (-0.8 + 24.8*pow(z,2) + 9.6*pow(z,4) + 16.*pow(z,6) +
                  pow(y,4)*(-16. + 16.*pow(z,2)) + pow(x,4)*(8. + 16.*pow(z,2)) +
                  pow(y,2)*(25.6 - 6.399999999999999*pow(z,2) + 32.*pow(z,4)) +
                  pow(x,2)*(24.8 + 17.6*pow(z,2) + 32.*pow(z,4) + pow(y,2)*(-8. + 32.*pow(z,2)))))/
              ((0.5 + pow(z,2))*(0.1 + pow(x,2) + pow(y,2) + pow(z,2))*
                (4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))) +
             2*(-1 + pow(x,2))*(-1 + pow(y,2))*
              (1 - 3*pow(z,2) + (z*(1 - pow(z,2))*
                   ((z*(0.8 - 24.*pow(x,2) - 24.*pow(y,2) - 24.*pow(z,2))*(0.5 + pow(z,2)))/
                      ((0.1 + pow(x,2) + pow(y,2) + pow(z,2))*
                        (4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))) -
                     2*z*log((4.*pow(x,2) + 4.*pow(y,2) + 4.*pow(z,2))/
                        pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4))))/(2.*pow(0.5 + pow(z,2),2)))))
          /pow(1. + 2*pow(z,2),2)))/2.;

*/
/*
  sol[0] = (-1 + pow(x,2))*(-1 + pow(y,2))*z*(1 - pow(z,2))*
    pow((4.*pow(x,2))/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4) +
     (4.*pow(y,2))/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4) +
	  (4.*pow(z,2))/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4),1/(2.*(0.5 + pow(z,2))));
          */
}

void Exact3DExp(TPZVec<REAL> &pt, TPZVec<REAL> &sol, TPZFMatrix &dsol) {

  REAL x = pt[0];
  REAL y = pt[1];
  REAL z = pt[2];

  sol[0] = (1 - x) * (1 - y) * (1 - z) * (x + y + z) * x * y * z;

  dsol(0,0) = (1 - x)*x*(1 - y)*y*(1 - z)*z + (1 - x)*(1 - y)*y*(1 - z)*z*(x + y + z) -
   x*(1 - y)*y*(1 - z)*z*(x + y + z);

  dsol(1,0) = (1 - x)*x*(1 - y)*y*(1 - z)*z + (1 - x)*x*(1 - y)*(1 - z)*z*(x + y + z) -
   (1 - x)*x*y*(1 - z)*z*(x + y + z);

  dsol(2,0) = (1 - x)*x*(1 - y)*y*(1 - z)*z + (1 - x)*x*(1 - y)*y*(1 - z)*(x + y + z) -
   (1 - x)*x*(1 - y)*y*z*(x + y + z);

/*  sol[0] = (-1 + pow(x,2))*(-1 + pow(y,2))*z*(1 - pow(z,2))*
   pow((4.*pow(x,2))/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4) +
     (4.*pow(y,2))/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4) +
     (4.*pow(z,2))/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4),1/(2.*(0.5 + pow(z,2))));

  dsol(0,0) = 2*x*(-1 + pow(y,2))*z*(1 - pow(z,2))*pow((4.*pow(x,2))/
       pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4) +
      (4.*pow(y,2))/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4) +
      (4.*pow(z,2))/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4),1/(2.*(0.5 + pow(z,2)))) +
   ((-1 + pow(x,2))*(-1 + pow(y,2))*z*(1 - pow(z,2))*
      ((-32.*pow(x,3))/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),5) -
        (32.*x*pow(y,2))/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),5) -
        (32.*x*pow(z,2))/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),5) +
        (8.*x)/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4))*
      pow((4.*pow(x,2))/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4) +
        (4.*pow(y,2))/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4) +
        (4.*pow(z,2))/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4),-1 + 1/(2.*(0.5 + pow(z,2)))))/
    (2.*(0.5 + pow(z,2)));

  dsol(1,0) = 2*(-1 + pow(x,2))*y*z*(1 - pow(z,2))*pow((4.*pow(x,2))/
       pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4) +
      (4.*pow(y,2))/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4) +
      (4.*pow(z,2))/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4),1/(2.*(0.5 + pow(z,2)))) +
   ((-1 + pow(x,2))*(-1 + pow(y,2))*z*(1 - pow(z,2))*
      ((-32.*pow(x,2)*y)/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),5) -
        (32.*pow(y,3))/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),5) -
        (32.*y*pow(z,2))/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),5) +
        (8.*y)/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4))*
      pow((4.*pow(x,2))/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4) +
        (4.*pow(y,2))/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4) +
        (4.*pow(z,2))/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4),-1 + 1/(2.*(0.5 + pow(z,2)))))/
    (2.*(0.5 + pow(z,2)));

  dsol(2,0) = -2*(-1 + pow(x,2))*(-1 + pow(y,2))*pow(z,2)*
    pow((4.*pow(x,2))/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4) +
      (4.*pow(y,2))/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4) +
      (4.*pow(z,2))/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4),1/(2.*(0.5 + pow(z,2)))) +
   (-1 + pow(x,2))*(-1 + pow(y,2))*(1 - pow(z,2))*
    pow((4.*pow(x,2))/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4) +
      (4.*pow(y,2))/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4) +
      (4.*pow(z,2))/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4),1/(2.*(0.5 + pow(z,2)))) +
   (-1 + pow(x,2))*(-1 + pow(y,2))*z*(1 - pow(z,2))*
    pow((4.*pow(x,2))/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4) +
      (4.*pow(y,2))/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4) +
      (4.*pow(z,2))/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4),1/(2.*(0.5 + pow(z,2))))*
    (((-32.*pow(x,2)*z)/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),5) -
         (32.*pow(y,2)*z)/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),5) -
         (32.*pow(z,3))/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),5) +
         (8.*z)/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4))/
       (2.*(0.5 + pow(z,2))*((4.*pow(x,2))/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4) +
           (4.*pow(y,2))/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4) +
           (4.*pow(z,2))/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4))) -
      (z*log((4.*pow(x,2))/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4) +
           (4.*pow(y,2))/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4) +
           (4.*pow(z,2))/pow(0.1 + pow(x,2) + pow(y,2) + pow(z,2),4)))/pow(0.5 + pow(z,2),2));
           */
}
