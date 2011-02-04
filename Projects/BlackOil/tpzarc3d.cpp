#include "tpzarc3d.h"
#include "pzshapelinear.h"

#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzgmesh.h"
#include "pzgeoel.h"
#include "pznoderep.h"
#include "pzgnode.h"
#include "pzreal.h"

using namespace std;
using namespace pzgeom;
using namespace pzshape;
using namespace pztopology;

//////////////////
void TPZArc3D::ComputeAtributes(TPZFMatrix &coord)
{
     /** fIBaseCn -> Basis Change Matrix: from Base(R2) to Canonic(R3) | fIBaseCn 1st column = BaseX | fIBaseCn 2nd column = BaseY */
     TPZFMatrix IBaseCnCP(3,3,0.), NotUsedHere(3,3,0.);
     fIBaseCn.Resize(3,3); fICnBase.Resize(3,3); fCenter3D.Resize(3);
     for(int i = 0; i < 3; i++)
     {
          IBaseCnCP(i,0) = coord(i,0) - coord(i,2);
          IBaseCnCP(i,1) = coord(i,1) - coord(i,2);
     }

     /** fIBaseCn 3rd column = BaseZ = Cross[BaseX,BaseY] |  */
     IBaseCnCP(0,2) = -coord(1,1)*coord(2,0) + coord(1,2)*coord(2,0) + coord(1,0)*coord(2,1) - coord(1,2)*coord(2,1) - coord(1,0)*coord(2,2) + coord(1,1)*coord(2,2);
     IBaseCnCP(1,2) =  coord(0,1)*coord(2,0) - coord(0,2)*coord(2,0) - coord(0,0)*coord(2,1) + coord(0,2)*coord(2,1) + coord(0,0)*coord(2,2) - coord(0,1)*coord(2,2);
     IBaseCnCP(2,2) = -coord(0,1)*coord(1,0) + coord(0,2)*coord(1,0) + coord(0,0)*coord(1,1) - coord(0,2)*coord(1,1) - coord(0,0)*coord(1,2) + coord(0,1)*coord(1,2);

     TPZFMatrix axest;
     IBaseCnCP.GramSchmidt(axest,NotUsedHere);
     axest.Transpose(&fICnBase);

     /** fICnBase -> Basis Change Matrix: from Canonic(R3) to Base(R2) | fICnBase(i,0) = Inverse[fIBaseCn] */
     fIBaseCn = fICnBase; fIBaseCn.Transpose();

     double Xa, Ya, Xb, Yb, Angle;
     ComputeR2Points(coord,Xa,Ya,Xb,Yb,Angle);

     /** Computing the Center Coordinates in R2 base */
     double Xcenter = Xa/2.;
     double Ycenter = (-Xa*Xb + Xb*Xb + Yb*Yb) / (2.*Yb);

     /** Computing the Center Coordinates in R3 base */
     TPZVec<REAL> Temp(3,0.);
     Temp[0] = Xcenter; Temp[1] = Ycenter; Temp[2] = 0.; double temp = 0.;
     for(int i = 0; i < 3; i++)
     {
          for(int j = 0; j < 3; j++) temp += fIBaseCn(i,j)*Temp[j];
          fCenter3D[i] = temp + coord(i,2);
          temp = 0.;
     }

     /** Computing Radius */
     fRadius = sqrt( Xcenter*Xcenter + Ycenter*Ycenter );
}

/////////////
/** This method compute the 3 given points with respect to R2 Basis */
void TPZArc3D::ComputeR2Points(TPZFMatrix &coord, double &xa, double &ya, double &xb, double &yb, double &angle)
{
     /** vector (ini - middle) written in new R2 base */
     TPZVec<REAL> Axe(3,0.), Temp(3,0.);
     for(int i = 0; i < 3; i++) Axe[i] = coord(i,0) - coord(i,2);
     for(int i = 0; i < 3; i++)
     {
          for(int j = 0; j < 3; j++) Temp[i] += fICnBase(i,j)*Axe[j];
     }
     xa = Temp[0]; ya = Temp[1];
     Temp.Fill(0.);

     /** vector (final - middle) written in new R2 base */
     for(int i = 0; i < 3; i++) Axe[i] = coord(i,1) - coord(i,2);
     for(int i = 0; i < 3; i++)
     {
          for(int j = 0; j < 3; j++) Temp[i] += fICnBase(i,j)*Axe[j];
     }
     xb = Temp[0]; yb = Temp[1];

     angle = ArcAngle(coord,xa, ya, xb, yb);
}

///////////////
/** This method return the absolute angle with respect of the arc formed between (ini - center) and (fin - center), passing by midnode
    Note: (xm,ym) don't appear because this coordinates are always (0,0) - it's the origin of R2 basis */
double TPZArc3D::ArcAngle(TPZFMatrix &coord, double xa, double ya, double xb, double yb)
{
     double cos, Angle1, Angle2, Angle3, Xcenter, Ycenter;

     /** Computing the Center Coordinates in this R2 base */
     Xcenter = xa/2.;
     Ycenter = (-xa*xb + xb*xb + yb*yb) / (2.*yb);

     /** angle between (ini-center) 'n' (mid-center) */
     cos = (-((xa - Xcenter)*Xcenter) - (ya - Ycenter)*Ycenter) / (sqrt(pow(xa - Xcenter,2) + pow(ya - Ycenter,2))*sqrt(pow(Xcenter,2) + pow(Ycenter,2)));
     if(cos < -0.999) cos = -1.; if(cos > 0.999) cos = 1.;
     Angle1 = acos(cos);

     /** angle between (fin-center) 'n' (mid-center) */
     cos = (-((xb - Xcenter)*Xcenter) - (yb - Ycenter)*Ycenter) / (sqrt(pow(xb - Xcenter,2) + pow(yb - Ycenter,2))*sqrt(pow(Xcenter,2) + pow(Ycenter,2)));
     if(cos < -0.999) cos = -1.; if(cos > 0.999) cos = 1.;
     Angle2 = acos(cos);

     /** angle between (ini-center) 'n' (fin-center) */
     cos = ((xa - Xcenter)*(xb - Xcenter) + (ya - Ycenter)*(yb - Ycenter)) / (sqrt(pow(xa - Xcenter,2) + pow(ya - Ycenter,2))*sqrt(pow(xb - Xcenter,2) + pow(yb - Ycenter,2)));
     if(cos < -0.999) cos = -1.; if(cos > 0.999) cos = 1.;
     Angle3 = acos(cos);

     /** verification if midpoint is in smaller arc angle (<= pi) or in the bigger arc angle (>= pi)
         Note: smaller and bigger arc angles reffers to the angle formed between
         (ini-center) and (fin-center) vectors, where [smaller + bigger = 2PI] */
     if( fabs(Angle3/(Angle1 + Angle2)) > 0.99 ) return Angle3; /** Smaller Arc Angle = Angle3 */
     else return (2.*acos(-1) - Angle3); /** Bigger Arc Angle = 2Pi - Angle3 */
}

///////////////
/** Mapping -> result = f(NodesCoord,qsi) */
void TPZArc3D::X(TPZFMatrix &nodes,TPZVec<REAL> &loc,TPZVec<REAL> &result)
{

     /** Cross[(mid-ini),(fin-ini)] */
     double CrossX = fabs(-(nodes(1,1)*nodes(2,0)) + nodes(1,2)*nodes(2,0) + nodes(1,0)*nodes(2,1) - nodes(1,2)*nodes(2,1) - nodes(1,0)*nodes(2,2) + nodes(1,1)*nodes(2,2));
     double CrossY = fabs(nodes(0,1)*nodes(2,0) - nodes(0,2)*nodes(2,0) - nodes(0,0)*nodes(2,1) + nodes(0,2)*nodes(2,1) + nodes(0,0)*nodes(2,2) - nodes(0,1)*nodes(2,2));
     double CrossZ = fabs(-(nodes(0,1)*nodes(1,0)) + nodes(0,2)*nodes(1,0) + nodes(0,0)*nodes(1,1) - nodes(0,2)*nodes(1,1) - nodes(0,0)*nodes(1,2) + nodes(0,1)*nodes(1,2));
     /** If Cross[(mid-ini),(fin-ini)] == 0, than the 3 given points are co-linear */
     if(CrossX <= 1.E-3 && CrossY <= 1.E-3 && CrossZ <= 1.E-3)
     {
          /** (fin - ini) */
          double dx1 = nodes(0,1) - nodes(0,0); double dy1 = nodes(1,1) - nodes(1,0); double dz1 = nodes(2,1) - nodes(2,0);
          double norm1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);

          /** (mid - ini) */
          double dx2 = nodes(0,2) - nodes(0,0); double dy2 = nodes(1,2) - nodes(1,0); double dz2 = nodes(2,2) - nodes(2,0);
          double norm2 = sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);

          /** (mid - fin) */
          double dx3 = nodes(0,2) - nodes(0,1); double dy3 = nodes(1,2) - nodes(1,1); double dz3 = nodes(2,2) - nodes(2,1);
          double norm3 = sqrt(dx3*dx3 + dy3*dy3 + dz3*dz3);

          if(fabs(norm1) < 0.001)
          {
               cout << "\nInitial and Final Nodes Are Coincidents! There is no ARC!\nMethod Aborted!\n";
               exit(-1);
          }
          if(norm1 < (norm2 + norm3))
          {
               cout << "\n*** The 3 Points are co-linear and mid point\nis outside initial and final range! ***\nThis Results an Infinite Arc... \nMethod Aborted!\n";
               exit(-1);
          }
          double desloc = (loc[0] + 1.)*norm1/2.;
          result[0] = nodes(0,0) + desloc*dx1/norm1;
          result[1] = nodes(1,0) + desloc*dy1/norm1;
          result[2] = nodes(2,0) + desloc*dz1/norm1;
          return;
     }

     /** Computing Atributes (Center, Radius), since the 3 poinst are NOT co-linear! */
     ComputeAtributes(nodes);

     double Xa, Ya, Xb, Yb, Angle;
     ComputeR2Points(nodes,Xa,Ya,Xb,Yb,Angle);

     /** Computing the Center Coordinates in R2 base and the Angle between Va and Vb */
     double Xcenter = Xa*Xa*Yb / (2.*Xa*Yb);
     double Ycenter = (-Xa*Xa*Xb + Xa*(Xb*Xb + Yb*Yb)) / (2.*Xa*Yb);

     TPZFMatrix RotMatrix(2,2); double deflection = Angle * fRadius * (loc[0] + 1.) / (2.*fRadius);
     RotMatrix(0,0) =  cos(deflection); RotMatrix(0,1) = sin(deflection);
     RotMatrix(1,0) = -sin(deflection); RotMatrix(1,1) = cos(deflection);

     /** Computing initialVector = (iniR2 - CenterR2) */
     TPZVec<REAL> initialVector(2,0.);
     TPZVec<REAL> MappedBASE2D(3,0.);
     initialVector[0] = Xa - Xcenter; initialVector[1] = Ya - Ycenter; double vectRotated = 0.;

     /** MappedPoint_R2 = centerCoord + vectorRotated , where Vx = RotationMatrix . Va */
     double centerCoord;
     for(int i = 0; i < 2; i++)
     {
          for(int j = 0; j < 2; j++) vectRotated += RotMatrix(i,j)*initialVector[j];
          centerCoord = (1-i)*Xcenter + i*Ycenter;
          MappedBASE2D[i] = centerCoord + vectRotated;
          vectRotated = 0.;
     }

     /** Changing Basis of Obtained MappedPoint from R2 to R3 */
     MappedBASE2D[2] = 0.;
     for(int i = 0; i < 3; i++)
     {
          for(int j = 0; j < 3; j++) vectRotated += fIBaseCn(i,j)*MappedBASE2D[j];
          result[i] = vectRotated + nodes(i,2);
          if( fabs(result[i]) < 0.001) result[i] = 0.;
          vectRotated = 0.;
     }
}

///////////////
void TPZArc3D::Jacobian(TPZFMatrix &coord, TPZVec<REAL> &par, TPZFMatrix &jacobian, TPZFMatrix &axes, REAL &detjac, TPZFMatrix &jacinv)
{
     jacobian.Resize(1,1); axes.Resize(1,3); jacinv.Resize(1,1);
     /** Cross[(mid-ini),(fin-ini)] */
     double CrossX = fabs(-(coord(1,1)*coord(2,0)) + coord(1,2)*coord(2,0) + coord(1,0)*coord(2,1) - coord(1,2)*coord(2,1) - coord(1,0)*coord(2,2) + coord(1,1)*coord(2,2));
     double CrossY = fabs(coord(0,1)*coord(2,0) - coord(0,2)*coord(2,0) - coord(0,0)*coord(2,1) + coord(0,2)*coord(2,1) + coord(0,0)*coord(2,2) - coord(0,1)*coord(2,2));
     double CrossZ = fabs(-(coord(0,1)*coord(1,0)) + coord(0,2)*coord(1,0) + coord(0,0)*coord(1,1) - coord(0,2)*coord(1,1) - coord(0,0)*coord(1,2) + coord(0,1)*coord(1,2));

     /** Verifying if the 3 points are co-linear, i.e. Cross[(mid-ini),(fin-ini)] == 0 */
     if(CrossX <= 1.E-3 && CrossY <= 1.E-3 && CrossZ <= 1.E-3)
     {
          /** (fin - ini) */
          double dx1 = coord(0,1) - coord(0,0); double dy1 = coord(1,1) - coord(1,0); double dz1 = coord(2,1) - coord(2,0);
          double norm1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);

          /** (mid - ini) */
          double dx2 = coord(0,2) - coord(0,0); double dy2 = coord(1,2) - coord(1,0); double dz2 = coord(2,2) - coord(2,0);
          double norm2 = sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);

          /** (mid - fin) */
          double dx3 = coord(0,2) - coord(0,1); double dy3 = coord(1,2) - coord(1,1); double dz3 = coord(2,2) - coord(2,1);
          double norm3 = sqrt(dx3*dx3 + dy3*dy3 + dz3*dz3);

          if(norm1 < 1.E-12)
          {
               cout << "\nInitial and Final Nodes Are Coincidents! There is no ARC!\nMethod Aborted!\n";
               exit(-1);
          }
          if(norm1 < (norm2 + norm3))
          {
               cout << "\n*** The 3 Points are co-linear and mid point\nis outside initial and final range! ***\nThis Results an Infinite Arc... \nMethod Aborted!\n"; exit(-1);
          }
          jacobian(0,0) = norm1/2.; jacinv(0,0) = 1./jacobian(0,0); detjac = jacobian(0,0);
          axes(0,0) = dx1/norm1; axes(0,1) = dy1/norm1; axes(0,2) = dz1/norm1;
          return;
     }

     /** Computing Atributes (Center, Radius, ArcLengh, ArcAngle), since the 3 poinst are NOT co-linear! */
     ComputeAtributes(coord);

     double Xa, Ya, Xb, Yb, Angle; ComputeR2Points(coord,Xa,Ya,Xb,Yb,Angle);
     jacobian(0,0) = Angle * fRadius/2.; jacinv(0,0) = 1./jacobian(0,0); detjac = jacobian(0,0);

     /** Computing Axes */
     TPZVec< REAL > Vpc(3), Vpa(3), Vpb(3), Vt(3), OUT(3);

     TPZVec< REAL > middle(1); middle[0] = 0.;
     X(coord,middle,OUT);

     /** Vector From MappedPoint to Ini */
     Vpa[0] = coord(0,0) - OUT[0]; Vpa[1] = coord(1,0) - OUT[1]; Vpa[2] = coord(2,0) - OUT[2];

     /** Vector From MappedPoint to Fin */
     Vpb[0] = coord(0,1) - OUT[0]; Vpb[1] = coord(1,1) - OUT[1]; Vpb[2] = coord(2,1) - OUT[2];

     X(coord,par,OUT);

     /** Vector From MappedPoint to Center */
     Vpc[0] = fCenter3D[0] - OUT[0]; Vpc[1] = fCenter3D[1] - OUT[1]; Vpc[2] = fCenter3D[2] - OUT[2];

     /** Tangent Vector From Point in the Arc */
     Vt[0] =  Vpa[1]*Vpb[0]*Vpc[1] - Vpa[0]*Vpb[1]*Vpc[1] + Vpa[2]*Vpb[0]*Vpc[2] - Vpa[0]*Vpb[2]*Vpc[2];
     Vt[1] = -Vpa[1]*Vpb[0]*Vpc[0] + Vpa[0]*Vpb[1]*Vpc[0] + Vpa[2]*Vpb[1]*Vpc[2] - Vpa[1]*Vpb[2]*Vpc[2];
     Vt[2] = -Vpa[2]*Vpb[0]*Vpc[0] + Vpa[0]*Vpb[2]*Vpc[0] - Vpa[2]*Vpb[1]*Vpc[1] + Vpa[1]*Vpb[2]*Vpc[1];

     double Vtnorm = 0.;
     for(int i = 0; i < 3; i++)
     {
          if( fabs(Vt[i]) < 1.E-12 ) Vt[i] = 0.;
          Vtnorm += Vt[i]*Vt[i];
     }
     for(int j = 0; j < 3; j++) axes(0,j) = Vt[j]/sqrt(Vtnorm);
}

///////////////
TPZGeoEl *TPZArc3D::CreateBCGeoEl(TPZGeoEl *orig, int side,int bc)
{
     if(side==2)
     {
          TPZManVector<int> nodes(3);
          nodes[0] = orig->SideNodeIndex(side,0); nodes[1] = orig->SideNodeIndex(side,1); nodes[2] = orig->SideNodeIndex(side,2); int index;
          TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EOned,nodes,bc,index);
          TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,TPZShapeLinear::SideConnectLocId(side,0)));
          TPZGeoElSide(gel,1).SetConnectivity(TPZGeoElSide(orig,TPZShapeLinear::SideConnectLocId(side,1)));
          TPZGeoElSide(gel,2).SetConnectivity(TPZGeoElSide(orig,side));
          return gel;
     }

     else if(side==0 || side==1)
     {
          TPZManVector<int> nodeindexes(1);
          nodeindexes[0] = orig->SideNodeIndex(side,0); int index;
          TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EPoint,nodeindexes,bc,index);
          TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,side));
          return gel;
     }

     else PZError << "\nTPZGeoLinear::CreateBCGeoEl. Side = " << side << endl;
     return 0;
}
