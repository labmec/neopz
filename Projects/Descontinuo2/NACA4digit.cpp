#include "NACA4digit.h"
#include "ratio.h"
#include "pzeuleranalysis.h"
#include "pzconslaw.h"
#include "pzmaterial.h"
#include "pzbndcond.h"
#include "pzeulerconslaw.h"
#include "pzartdiff.h"
#include "pzreal.h"
#include "pzvec.h"
#include "pzflowcmesh.h"
#include "pzgmesh.h"
#include "pzgeoelbc.h"
#include <iostream>
#include <fstream>
#include "TPZGeoElement.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "pzrefquad.h"
#include "TPZGeoLinear.h"
#include "TPZRefLinear.h"
#include "pzbstrmatrix.h"
#include "pzstepsolver.h"
#include "pzblock.h"

using namespace std;

 using namespace pzgeom;
 using namespace pzshape;
 using namespace pzrefine;
REAL TPZNACAXXXX::NearestParameter(TPZVec<REAL> &pt, int &uplow, int maxPt) {

  REAL distminlow = 20.*fCord;
  REAL distminup = 20.*fCord;
  REAL distlow,distup;
  REAL ptl[2],ptu[2];
  int ip,maxp=maxPt;
  REAL par,parlow,parup;
  for(ip=0; ip<=maxp; ip++) {
    par = ip*fCord/maxp;
    ptu[0] = xua(par);
    ptu[1] = yua(par);
    ptl[0] = xla(par);
    ptl[1] = yla(par);
    distlow = (ptl[0]-pt[0])*(ptl[0]-pt[0])+(ptl[1]-pt[1])*(ptl[1]-pt[1]);
    distup = (ptu[0]-pt[0])*(ptu[0]-pt[0])+(ptu[1]-pt[1])*(ptu[1]-pt[1]);
    if(distlow < distminlow) parlow = par;
    if(distup < distminup) parup = par;
    distminlow = distminlow < distlow ? distminlow : distlow;
    distminup = distminup < distup ? distminup : distup;
  }
  REAL delpar = 0.1/maxp;
  if(distminlow < distminup) {
    uplow = 0;
    REAL distprev = distminlow;
    par = parlow;
    while(fabs(delpar) > 0.00001/maxp) {
      ptl[0] = xla(par+delpar);
      ptl[1] = yla(par+delpar);
      distlow = (ptl[0]-pt[0])*(ptl[0]-pt[0])+(ptl[1]-pt[1])*(ptl[1]-pt[1]);
      if(distlow < distprev) {
	par += delpar;
	distprev = distlow;
      } else if (delpar > 0.) {
	delpar *= -1.;
      } else {
	delpar *= -0.1;
      }
    }
  } else {
    uplow = 1;
    REAL distprev = distminup;
    par = parup;
    while(fabs(delpar) > 0.001/maxp) {
      ptu[0] = xua(par+delpar);
      ptu[1] = yua(par+delpar);
      distup = (ptu[0]-pt[0])*(ptu[0]-pt[0])+(ptu[1]-pt[1])*(ptu[1]-pt[1]);
      if(distup < distprev) {
	par += delpar;
	distprev = distup;
      } else if (delpar > 0.) {
	delpar *= -1.;
      } else {
	delpar *= -0.1;
      }
    }
  }
  return par;
}
// creates an one-quadrilateral element mesh

// This file generates a mesh for the NACA airfoils

/*
double scale = 4.;
double entrance = 5. * scale,
             exitlength = 4. * scale,
	     cord = 5.,
	     height = 10. * scale,
	     extraheight = scale,
	     q = 1.5,
	     qn = 1.5;
	     */
/*
double scale = 10.; //4;
double entrance = 5. * scale,
             exitlength = 4. * scale,
	     cord = 1.,//5.,
	     height = 10. * scale,
	     extraheight = scale,
	     q = 1.5,
	     qn = 1.5;
*/

double
entrance = 25.,
exitlength = 25.,
cord = 1.,
extraheight = 5.,
height = 40.,
qn, ql, q,
scale = 5.;

int l, m, n, p, k;
void NACAPoints(TPZNACAXXXX &profile, TPZVec< TPZVec<REAL> > & pt, TPZVec< TPZVec< int> > &elms, int nSubdiv)
{

    
 m = nSubdiv;

 n = 5 * m / 9 * (int) pow(scale, .6);
 l = n / 2;
 p = n / 4;

 q  = pow(7., 1./(double)m);
 qn = pow(196., 1./(double)n);
 ql = pow(7., 1./(double)l);




 n = 3 * m / 2;
 l = m / 4;
 p = n / 2;


 q  = pow(7., 1./(double)m);
 qn = pow(98., 1./(double)n);
 ql = pow(28., 1./(double)l);


   int index, indexPt, indexPt2;

   TPZVec<REAL> coord(3), coordBC(3), coordNACA(3);
   elms.Resize(m*n*2);
   pt.Resize(2 * m + n * (2*m+1) + /*exit elements*/ p*(2*l + 1) + (2*(l+p) - 1)*(n-l) );

   cout << "\nNumber of Points: "<< pt.NElements();

   // defining points on the surface of the airfoil
   // i: airfoil index;
   // j: radial index
   int i, j;
   for(i = 1; i < m; i++)
      {
         REAL x = xtrig(i, m + 1);//xpg(q, i, m);
	 coord[0] = profile.xua(x * cord);
	 coord[1] = profile.yua(x * cord);
	 coord[2] = 0.;
	 pt[i] = coord;

	 coord[0] = profile.xla(x * cord);
	 coord[1] = profile.yla(x * cord);
	 coord[2] = 0.;
	 pt[2*m -i] = coord;

      }
   coord[0] = profile.xua(0.);
   coord[1] = profile.yua(0. * cord);
   coord[2] = 0.;
   pt[0] = coord;

   coord[0] = profile.xla(1. * cord);
   coord[1] = ( profile.yla(1. * cord) +
                profile.yua(1. * cord) ) / 2.;//height/2.;
   coord[2] = 0.;
   pt[m] = coord;

   //defining points at the boundary
   // index of the leftmost centered point
   index = 2*m + (n-1) * ((2*m)+1);
   coord[0] = 0.;
   coord[1] = height/2.;
   coord[2] = 0.;
   pt[index] = coord;

   for(i = 1; i < l; i++)
   {
      coord[0] = 0.;
      coord[1] =  xpg(q, i, l) * height/2. + height/2.;
      coord[2] = 0.;
      pt[index + i] = coord;

      coord[0] = 0.;
      coord[1] = -xpg(q, i, l) * height/2. + height/2.;
      coord[2] = 0.;
      pt[index + (2*m+1) - i] = coord;
   }

   for(i = l; i <= m; i++)
   {
      coord[0] = xpg(1./q, i-l, m-l) * (entrance + cord);
      coord[1] = height + xpg(1./q, i-l, m-l) * extraheight;
      coord[2] = 0.;
      pt[index + i] = coord;

      coord[0] = xpg(1./q, i-l, m-l) * (entrance + cord);
      coord[1] = 0. - xpg(1./q, i-l, m-l) * extraheight;
      coord[2] = 0.;
      pt[index + (2*m+1) - i] = coord;
   }

   // defining intermediate points
   for(j = 1; j < n; j++)
   {
      //creating division rule
      REAL ratio = xpg(qn, j, n);

      // resolving entrance centered point
      index = n * ((2*m)+1) - 1; // index of BC point
      indexPt = j * ((2*m)+1) - 1; // index of leftmost layer point
      coordNACA = pt[0];
      coordBC   = pt[index];
      coord[0] = ratio * coordBC[0] + (1.-ratio) * coordNACA[0];
      coord[1] = ratio * coordBC[1] + (1.-ratio) * coordNACA[1];
      coord[2] = ratio * coordBC[2] + (1.-ratio) * coordNACA[2];
      pt[indexPt] = coord;

      // resolving exit points
      coordNACA = pt[m];

      indexPt2 = indexPt + m;
      coordBC   = pt[index + m];// upper point
      coord[0] = ratio * coordBC[0] + (1.-ratio) * coordNACA[0];
      coord[1] = ratio * coordBC[1] + (1.-ratio) * coordNACA[1];
      coord[2] = ratio * coordBC[2] + (1.-ratio) * coordNACA[2];
      pt[indexPt2] = coord;

      indexPt2 = indexPt + m + 1;
      coordBC   = pt[index + m + 1];// bottom point
      coord[0] = ratio * coordBC[0] + (1.-ratio) * coordNACA[0];
      coord[1] = ratio * coordBC[1] + (1.-ratio) * coordNACA[1];
      coord[2] = ratio * coordBC[2] + (1.-ratio) * coordNACA[2];
      pt[indexPt2] = coord;

      // resolving other points
      for(i = 1; i < m; i++)
      {
         // upper points
         coordNACA = pt[i];
         indexPt2 = indexPt + i;
         coordBC   = pt[index + i];// BC upper point
         coord[0] = ratio * coordBC[0] + (1.-ratio) * coordNACA[0];
         coord[1] = ratio * coordBC[1] + (1.-ratio) * coordNACA[1];
         coord[2] = ratio * coordBC[2] + (1.-ratio) * coordNACA[2];
         pt[indexPt2] = coord;

         // bottom points
         coordNACA = pt[2*m-i];
         indexPt2 = indexPt+2*m+1-i;
         coordBC   = pt[index+2*m+1-i];// BC bottom point
         coord[0] = ratio * coordBC[0] + (1.-ratio) * coordNACA[0];
         coord[1] = ratio * coordBC[1] + (1.-ratio) * coordNACA[1];
         coord[2] = ratio * coordBC[2] + (1.-ratio) * coordNACA[2];
         pt[indexPt2] = coord;

      }

   }

   //Creating Exit points with fewer elements

   // Creating exit points
   int firstExitPt = 2*m + (2 * m + 1) * n; // first exit point index
   // resolving ponts closer to NACA
   for(i = 1; i <= p; i++)
   {

      REAL ratio = xpg(sqrt(qn), i, p);

      // center point
      index = m; // index of existent point
      indexPt = firstExitPt + (2*l+1) * (i-1);
      coord[0] = /*pt[index][0] */ entrance + cord + 1.5 * cord * ratio;
      coord[1] = pt[index][1];
      coord[2] = 0.;
      pt[indexPt] = coord;

      for(j = 1; j <= l; j++)
      {
         //upper points
         index = j * (2 * m + 1) + m - 1; // index of existent point
         indexPt2 = indexPt + j;
         coord[0] = entrance + cord + 1.5 * cord * ratio;
	 coord[1] = pt[index][1];// - extraheight * ratio * xpg(sqrt(qn), j, n);
         coord[2] = 0.;
         pt[indexPt2] = coord;

         // bottom points
         index = j * (2 * m + 1) + m; // index of existent point
         indexPt2 = indexPt + l + j;
         coord[0] = entrance + cord + 1.5 * cord * ratio;
         coord[1] = pt[index][1];// + extraheight * ratio * xpg(sqrt(qn), j, n);
         coord[2] = 0.;
         pt[indexPt2] = coord;
      }
   }

   // creating outer upper and lower exit points
   int firstExitPt2 = 2*m + (2 * m + 1) * n + p * (l*2+1); // first exit point index of this series
   int outerCenterExitPt = firstExitPt2 + (2*(l+p)-1) * (n-l - 1);

   // center point
   indexPt2 = firstExitPt2 + (2*(l + p)-1)*(n-l-1);
   coord[0] = entrance + cord + exitlength;
   coord[1] = height / 2.;
   coord[2] = 0.;
   pt[indexPt2] = coord;
   // vertical exit face
   for(j = 1; j <= l; j++)
   {
      REAL ratio = xpg(pow(4., 1./l), j, l);
      //upper points
      indexPt2 = outerCenterExitPt + j;
      coord[0] = entrance + cord + exitlength;
      coord[1] = height / 2. * (1. + ratio);
      coord[2] = 0.;
      pt[indexPt2] = coord;

      // bottom points
      indexPt2 = outerCenterExitPt + l + p - 1 + j;
      coord[0] = entrance + cord + exitlength;
      coord[1] = height / 2. * (1. - ratio);
      coord[2] = 0.;
      pt[indexPt2] = coord;
   }
   // upper/lower angled exit faces
   for(j = 1; j < p; j++)
   {
      REAL ratio = xpg(qn, p-j, p);
      //upper points
      indexPt = m + (2 * m + 1) * n - 1; // first exit point index//2*m + (2 * m + 1) * n + (2*l+1) * (p-j);
      indexPt2 = outerCenterExitPt + l;

      coord[0] = pt[indexPt2][0] * ratio + pt[indexPt][0] * (1.-ratio);
      coord[1] = pt[indexPt2][1] * ratio + pt[indexPt][1] * (1.-ratio);
      coord[2] = pt[indexPt2][2] * ratio + pt[indexPt][2] * (1.-ratio);
      pt[indexPt2 + j] = coord;

      // bottom points
      indexPt = m + (2 * m + 1) * n; // first exit point index;
      indexPt2 = outerCenterExitPt + 2*l + p - 1;

      coord[0] = pt[indexPt2][0] * ratio + pt[indexPt][0] * (1.-ratio);
      coord[1] = pt[indexPt2][1] * ratio + pt[indexPt][1] * (1.-ratio);
      coord[2] = pt[indexPt2][2] * ratio + pt[indexPt][2] * (1.-ratio);
      pt[indexPt2 + j] = coord;
   }

   // creating intermediate points
   for(i = 1; i < n - l; i++)
   {
      REAL ratio = xpg(/*sqrt(qn)*/qn, i, n - l);
      //centered points
      index = firstExitPt2 + (i - 1) * (2 * (l+p) - 1);
      indexPt = firstExitPt2 - (2 * l + 1);
      indexPt2 = outerCenterExitPt;
      coord[0] = pt[indexPt2][0] * ratio + pt[indexPt][0] * (1-ratio);
      coord[1] = pt[indexPt2][1] * ratio + pt[indexPt][1] * (1-ratio);
      coord[2] = pt[indexPt2][2] * ratio + pt[indexPt][2] * (1-ratio);
      pt[index] = coord;

      // upper and lower points
      for(j = 1; j <= l; j++)
      {
         // upper points
         indexPt = firstExitPt2 - (2 * l + 1) + j;
	 indexPt2 = outerCenterExitPt + j;
         coord[0] = pt[indexPt2][0] * ratio + pt[indexPt][0] * (1-ratio);
         coord[1] = pt[indexPt2][1] * ratio + pt[indexPt][1] * (1-ratio);
         coord[2] = pt[indexPt2][2] * ratio + pt[indexPt][2] * (1-ratio);
         pt[index + j] = coord;

	 // lower points
	 indexPt = firstExitPt2 - l + j - 1;
	 indexPt2 = outerCenterExitPt + l + p - 1 + j;
         coord[0] = pt[indexPt2][0] * ratio + pt[indexPt][0] * (1-ratio);
         coord[1] = pt[indexPt2][1] * ratio + pt[indexPt][1] * (1-ratio);
         coord[2] = pt[indexPt2][2] * ratio + pt[indexPt][2] * (1-ratio);
	 pt[index + l + p - 1 + j] = coord;
      }
      for(j = 1; j < p; j++)
      {
         // upper points
         indexPt = 2*m + (m*2 + 1) * n + (p-j-1) * (2* l + 1) + l;
	 indexPt2 = outerCenterExitPt + j + l;
         coord[0] = pt[indexPt2][0] * ratio + pt[indexPt][0] * (1-ratio);
         coord[1] = pt[indexPt2][1] * ratio + pt[indexPt][1] * (1-ratio);
         coord[2] = pt[indexPt2][2] * ratio + pt[indexPt][2] * (1-ratio);
         pt[index + j + l] = coord;

	 // lower points
	 indexPt = 2*m + (m*2 + 1) * n + (p-j) * (2* l + 1) - 1;
	 indexPt2 = outerCenterExitPt + 2*l + p - 1 + j;
         coord[0] = pt[indexPt2][0] * ratio + pt[indexPt][0] * (1-ratio);
         coord[1] = pt[indexPt2][1] * ratio + pt[indexPt][1] * (1-ratio);
         coord[2] = pt[indexPt2][2] * ratio + pt[indexPt][2] * (1-ratio);
	 pt[index + 2*l + p - 1 + j] = coord;
      }
   }

   //////////

   // definig elements
   // quadrilateral data
   TPZVec< int > nodes(4);
   elms.Resize(m * n * 2 + p * l * 2 + (l+p) * (n-l) * 2);

   cout << "\nNumber of Elements: "<< elms.NElements() << endl;

   // the first row (closest to NACA is special)
   //upper elements
   nodes[0] = 0;
   nodes[1] = 1;
   nodes[2] = 1+2*m;
   nodes[3] = 2*m;
   elms[0] = nodes;
   // bottom elements
   nodes[0] = 2*m-1;
   nodes[1] = 0;
   nodes[2] = 2*m;
   nodes[3] = 4*m;
   elms[2*m-1] = nodes;
   for(i = 1; i < m; i++)
   {
   //upper elements
      index = i;
      nodes[0] = index;
      nodes[1] = nodes[0]+1;
      nodes[2] = nodes[1]+2*m;
      nodes[3] = nodes[2]-1;
      elms[index] = nodes;
   // bottom elements
      index = 2*m-i-1;
      nodes[0] = index;
      nodes[1] = nodes[0] + 1;
      nodes[2] = nodes[1] + 2 * m + 1;
      nodes[3] = nodes[2] - 1;
      elms[index] = nodes;
   }

   // elements at rows farther from NACA profile
   for(j = 1; j < n; j++)
   {
         // leftmost upper elements
         index = j * 2*m;
	 nodes[0] = index+ j - 1;
	 nodes[1] = nodes[0] + 1;
	 nodes[2] = nodes[1] + 2*m + 1;
	 nodes[3] = nodes[2] - 1;
	 elms[index] = nodes;

	 // leftmost lower elements
         index = (j+1) * 2*m - 1;
	 nodes[0] = index + j;
	 nodes[1] = nodes[0] - 2*m;//j * 2*m;//index + 2;
	 nodes[2] = nodes[1] + 2 * m + 1;
	 nodes[3] = nodes[0] + 2 * m + 1;
	 elms[index] = nodes;
      for(i = 1; i < m; i++)
      {
         // upper elements
         index = i + j * 2*m;
	 nodes[0] = index + j - 1;
	 nodes[1] = nodes[0] + 1;
	 nodes[2] = nodes[1] + 2*m + 1;
	 nodes[3] = nodes[2] - 1;
	 elms[index] = nodes;

	 // lower elements
         index = (j+1) * 2*m - i - 1;
	 nodes[0] = index + j;
	 nodes[1] = nodes[0] + 1;
	 nodes[2] = nodes[1] + 2 * m + 1;
	 nodes[3] = nodes[2] - 1;
	 elms[index] = nodes;
      }
   }

   index = 2 * m * n;
   // exit elements
   // indexing exit elements closer to NACA
   // upper element
   nodes[0] = m;
   nodes[1] = firstExitPt;
   nodes[2] = firstExitPt + 1;
   nodes[3] = 3 * m;
   elms[index] = nodes;
   index++;
   //lower element
   nodes[0] = 3*m+1;
   nodes[1] = firstExitPt + l + 1;
   nodes[2] = firstExitPt;
   nodes[3] = m;
   elms[index] = nodes;
   index++;
   for(i = 1; i < l; i++)
   {
      // upper element
      nodes[0] = m + (m*2+1)*i -1;
      nodes[1] = firstExitPt + i;
      nodes[2] = nodes[1] + 1;
      nodes[3] = nodes[0] + m*2+1;
      elms[index] = nodes;
      index++;
      // lower element
      nodes[3] = nodes[0] + 1;
      nodes[0] = nodes[0] + 2*m + 2;
      nodes[1] = firstExitPt + l + 1 + i;
      nodes[2] = nodes[1] - 1;
      elms[index] = nodes;
      index++;
   }
   // elements at right of the above
   for(i = 1; i < p; i++)
   {
      // centered elements
      indexPt = firstExitPt + (2*l+1)*(i-1); // centered point of such row
      // upper elements
      nodes[0] = indexPt;
      nodes[1] = nodes[0] + (2*l+1);
      nodes[2] = nodes[1] + 1;
      nodes[3] = nodes[0] + 1;
      elms[index] = nodes;
      index++;
      // lower elements
      nodes[0] = indexPt + l + 1;
      nodes[1] = nodes[0] + (2*l+1);
      nodes[2] = indexPt + (2*l+1);
      nodes[3] = indexPt;
      elms[index] = nodes;
      index++;
      // other elements
      for(j = 1; j < l; j++)
      {
         // upper elements
         nodes[0] = indexPt + j;
	 nodes[1] = nodes[0] + (2*l+1);
	 nodes[2] = nodes[1] + 1;
	 nodes[3] = nodes[0] + 1;
         elms[index] = nodes;
         index++;
	 // lower elements
         nodes[0] = indexPt + l + j + 1;
	 nodes[1] = nodes[0] + (2*l+1);
	 nodes[2] = nodes[1] - 1;
	 nodes[3] = nodes[0] - 1;
         elms[index] = nodes;
         index++;
      }
   }

   int indexEl = 2*m*n + 2*l*p; // index should be such ->test purposes
   // indexig rightmost elements
   // resolving exception points (first row)
   //firstExitPt2 = firstExitPt + p * (2*l+1);
   // centered points
   // upper point
   index = indexEl + l + p - 1;
   nodes[0] = firstExitPt2 - (2*l+1);
   nodes[1] = firstExitPt2;
   nodes[2] = nodes[1] + 1;
   nodes[3] = nodes[0] + 1;
   elms[index] = nodes;
   // lower point
   index ++;
   nodes[3] = nodes[0];
   nodes[0] += l + 1;
   nodes[1] = firstExitPt2 + l+p;
   nodes[2] = firstExitPt2;
   elms[index] = nodes;
   // other elements (yet vertical)
   for(j = 1; j < l; j++)
   {
      // upper elements
      index = indexEl + l + p - 1 - j;
      nodes[0] = firstExitPt2 - (2*l+1) + j;
      nodes[1] = firstExitPt2 + j;
      nodes[2] = nodes[1] + 1;
      nodes[3] = nodes[0] + 1;
      elms[index] = nodes;
      // lower elements
      index = indexEl + l + p + j;
      nodes[0] = firstExitPt2 - l + j;
      nodes[1] = firstExitPt2 + l + p + j;
      nodes[2] = nodes[1] - 1;
      nodes[3] = nodes[0] - 1;
      elms[index] = nodes;
   }
   // indexing angled elements
   //upper point
   index = indexEl;
   nodes[0] = m + (m*2+1)*l - 1;
   nodes[1] = firstExitPt + l;
   nodes[2] = firstExitPt2 + l + p -1;
   nodes[3] = nodes[0] + 2 * m + 1;
   elms[index] = nodes;
   //lower point
   index = indexEl + 2*(l+p) - 1;
   nodes[0] = m + (m*2+1)*(l+1);
   nodes[1] = firstExitPt2 + 2*(l+p-1);
   nodes[2] = firstExitPt + 2*l;
   nodes[3] = nodes[0] - (2 * m + 1);
   elms[index] = nodes;
   //other elements
   for(j = 1; j < p; j++)
   {
      // upper points
      index = indexEl + j;
      nodes[0] = firstExitPt + (2*l+1) * (j-1) + l;
      nodes[1] = nodes[0] + (2*l+1);
      nodes[2] = firstExitPt2 + l + p - j - 1;
      nodes[3] = nodes[2] + 1;
      elms[index] = nodes;
      // lower points
      index = indexEl + 2*(l + p) - j - 1;
      nodes[0] = firstExitPt + (2*l+1)*j + 2*l;
      nodes[1] = nodes[0] - (2*l+1);
      nodes[2] = firstExitPt2 + 2*(l+p) - j - 1;
      nodes[3] = nodes[2] - 1;
      elms[index] = nodes;
   }

   // creating elements of right parenthesis -> )
   // at the rightmost boundary.
   for(i = 1; i < n-l; i++)
   {
      indexEl = m*n*2 + 2*l*p + 2*(p+l)*i;
      // centered elements
      // upper elements
      index = indexEl + l + p -1;
      nodes[0] = firstExitPt2 + (2*(l+p) -1)*(i-1) + 1;
      nodes[1] = nodes[0] - 1;
      nodes[2] = nodes[1] + (l+p)*2 - 1;;
      nodes[3] = nodes[2] + 1;
      elms[index] = nodes;
      // lower elements
      index = indexEl + l + p;
      nodes[0] = firstExitPt2 + (2*(l+p) -1)*(i-1);
      nodes[1] = nodes[0] + (l+p);
      nodes[2] = nodes[1] + (l+p)*2 - 1;
      nodes[3] = nodes[2] - (l+p);
      elms[index] = nodes;

      //leftmost elements
      // upper elements
      index = indexEl;
      nodes[0] = m + (m*2 +1)*(i+l) - 1;
      nodes[1] = firstExitPt2 + (l+p) - 1 + (2*(l+p) -1)*(i-1);
      nodes[2] = nodes[1] + (2*(l+p) -1);
      nodes[3] = nodes[0] + (m*2 +1);
      elms[index] = nodes;
      // lower elements
      index = indexEl + (2*(l+p) -1);
      nodes[0] = firstExitPt2 + (l+p-1)*2 + (2*(l+p) -1)*(i-1);
      nodes[1] = m + (m*2 +1)*(i+l);
      nodes[2] = nodes[1] + (m*2 +1);
      nodes[3] = nodes[0] + (2*(l+p) -1);
      elms[index] = nodes;

      //other elements
      for(j = 0; j < l + p - 2; j++)
      {
         // upper elements
         index = indexEl + j + 1;
	 nodes[0] = firstExitPt2 + (l+p) - 1 - j + (2*(l+p) -1)*(i-1);
	 nodes[1] = nodes[0] - 1;
	 nodes[2] = nodes[1] + (2*(l+p) -1);
	 nodes[3] = nodes[0] + (2*(l+p) -1);
         elms[index] = nodes;
	 // lower elements
         index = indexEl + (l+p) * 2 - 2 - j;
	 nodes[0] = firstExitPt2 + (l+p)*2 - 3 - j + (2*(l+p) -1)*(i-1);
	 nodes[1] = nodes[0] + 1;
	 nodes[2] = nodes[1] + (2*(l+p) -1);
	 nodes[3] = nodes[0] + (2*(l+p) -1);
         elms[index] = nodes;
      }

   }
}

TPZGeoMesh * CreateNACAGeoMesh(TPZGeoMesh *gmesh, TPZNACAXXXX &profile, TPZVec< TPZVec< REAL > > & nodes,
                           TPZVec< TPZVec< int > > & elms,
			   MElementType ElType, int matId,
			   TPZVec<TPZGeoEl *> & gEls,
			   int nSubdiv)
{
   //TPZGeoMesh * gmesh = new TPZGeoMesh;

   gEls.Resize(elms.NElements());
   gmesh->NodeVec().Resize(nodes.NElements());
   int i,j;

   for(i = 0; i < nodes.NElements(); i++)
   {
      gmesh->NodeVec()[i].Initialize(nodes[i],*gmesh);
   }

   for( i = 0; i < elms.NElements(); i++)
   {
      gEls[i] = gmesh->CreateGeoElement(ElType, elms[i], matId, i);
   }

// Constructing neighborhood

   gmesh->BuildConnectivity();



   {// Dividing elements to create a mesh of 4 elems around NACA surface

     //first row near naca
     TPZVec<TPZGeoEl * > firstDiv, secondDiv, thirdDiv;
     TPZManVector<REAL, 3> pt(3,0.);
     int ii;//, jj;
     for(i = 0; i < 2*m; i++)
       {
	 gEls[i]->Divide(firstDiv);
	 pt[0] = firstDiv[0]->NodePtr(1)->Coord(0);
	 pt[1] = firstDiv[0]->NodePtr(1)->Coord(1);
	 profile.ProjectPoint(pt, 32*m + 1);
	 firstDiv[0]->NodePtr(1)->SetCoord(0,pt[0]);
	 firstDiv[0]->NodePtr(1)->SetCoord(1,pt[1]);

	 if(fabs(i-m+.5) > 3*m/4 || fabs(i-m+.5) < m/5)
	 {
            for(j = 0; j < 2; j++) {
	      firstDiv[j]->Divide(secondDiv);
	      pt[0] = secondDiv[0]->NodePtr(1)->Coord(0);
	      pt[1] = secondDiv[0]->NodePtr(1)->Coord(1);
	      profile.ProjectPoint(pt, 64*m + 1);
	      secondDiv[0]->NodePtr(1)->SetCoord(0,pt[0]);
	      secondDiv[0]->NodePtr(1)->SetCoord(1,pt[1]);
	      for(ii = 0; ii < 2; ii++)
	      {
	        secondDiv[ii]->Divide(thirdDiv);
	        pt[0] = thirdDiv[0]->NodePtr(1)->Coord(0);
	        pt[1] = thirdDiv[0]->NodePtr(1)->Coord(1);
	        profile.ProjectPoint(pt, 256*m + 1);
	        thirdDiv[0]->NodePtr(1)->SetCoord(0,pt[0]);
	        thirdDiv[0]->NodePtr(1)->SetCoord(1,pt[1]);
	      }
	      if(i == 0 || i == 2*m-1)
	         for(ii = 2; ii < 4; ii++)
	         {
	           secondDiv[ii]->Divide(thirdDiv);
	         }
	    }
	    for(j = 2; j < 4; j++)
	    {
              firstDiv[j]->Divide(secondDiv);
	    }
	 }else{
            for(j = 0; j < 2; j++) {
	      firstDiv[j]->Divide(secondDiv);
	      pt[0] = secondDiv[0]->NodePtr(1)->Coord(0);
	      pt[1] = secondDiv[0]->NodePtr(1)->Coord(1);
	      profile.ProjectPoint(pt, 8*m + 1);
	      secondDiv[0]->NodePtr(1)->SetCoord(0,pt[0]);
	      secondDiv[0]->NodePtr(1)->SetCoord(1,pt[1]);
	    }
         }
       }

     // first n/3 rows of elements near naca
     for(i = 2*m; i < 2*(n/3)*m; i++)
     {
        gEls[i]->Divide(firstDiv);
	if( ((i==2*m)||(i==4*m-1)) )
	    for(j = 0; j < 4; j++)
	    {
              firstDiv[j]->Divide(secondDiv);
	    }
     }

     // exit mesh
     for(j = 0; j < p; j++)
	{
	   for(i = 0; i < min(2*(p-j+1),2*l); i++)
           {
	     gEls[2*m*n+j*(l*2)+i]->Divide(firstDiv);
             if(i < 2 && j < p )
	       if(j < p/4)
	       {
	          for(ii = 0; ii < 2; ii++)
	          {
                     firstDiv[ii + 2*i]->Divide(secondDiv);
		     if(j < p/4)
		     {
                        secondDiv[2*i]->Divide(thirdDiv);
		        secondDiv[2*i+1]->Divide(thirdDiv);
		     }
		     firstDiv[ii + 2 -2*i]->Divide(secondDiv);
                  }
	       }else{
	          for(ii = 0; ii < 2; ii++)
	          {
                     firstDiv[ii + 2 * i]->Divide(secondDiv);
		  }

               }
	   }
	}
     // conic exit mesh
     for(i = 0; i < n-l; i++)
     {
        int index = 2*m*n + 2*l*p + l + p + 2*(l+p)*i;
        gEls[index-1]->Divide(firstDiv);
	{
	   if(i == 0)
	   {
	      firstDiv[0]->Divide(secondDiv);
	      firstDiv[1]->Divide(secondDiv);
	   }else{
	      firstDiv[1]->Divide(secondDiv);
	      firstDiv[2]->Divide(secondDiv);
	   }
	}
	gEls[index]->Divide(firstDiv);
	{
	   if(i == 0)
	   {
	      firstDiv[2]->Divide(secondDiv);
	      firstDiv[3]->Divide(secondDiv);
	   }else{
	      firstDiv[0]->Divide(secondDiv);
	      firstDiv[3]->Divide(secondDiv);
	   }
	}
     }

     // shock wings
     for(i = n/3; i < n/2-1; i++)
        for(j = l + i - n/3; j < m; j++)
        {
	   gEls[i*2*m+j]->Divide(firstDiv);
	   gEls[(i+1)*2*m-j-1]->Divide(firstDiv);
        }

   }

  // if(nSubdiv > 1)PZError << "CreateOneElGeoMesh unsupported number of subdivisions";

   return gmesh;
}


// Creating all the geometric and computational meshes
// for the reflected shock problem.

TPZFlowCompMesh *
   NACACompMesh(TPZFlowCompMesh *cmesh, REAL CFL, REAL delta,
                 int degree, int nSubdiv,
		 TPZArtDiffType DiffType,
		 TPZTimeDiscr Diff_TD,
		 TPZTimeDiscr ConvVol_TD,
   TPZTimeDiscr ConvFace_TD, std::ostream &options)
{
   TPZCompEl::SetgOrder(degree);
   REAL gamma = 1.4;
   int i;
   REAL Mach;

// Configuring the PZ to generate discontinuous elements
//    TPZGeoElement<TPZShapeQuad,TPZGeoQuad,TPZRefQuad>
//                 ::SetCreateFunction(TPZCompElDisc::CreateDisc);
//
//    TPZGeoElement<TPZShapeLinear,TPZGeoLinear,TPZRefLinear>
//                 ::SetCreateFunction(TPZCompElDisc::CreateDisc);

   const int dim = 2;
//   int interfdim = dim -1;
//   TPZCompElDisc::gInterfaceDimension = interfdim;


// Retrieving the point coordinates and element references
   TPZVec< TPZVec< REAL > > nodes;
   TPZVec< TPZVec< int  > > elms;
   TPZVec< TPZGeoEl *> gElem;
   const int digits = 12;
   TPZManVector<REAL> x0(3,0.);
   x0[0] = entrance;
   x0[1] = height/2.;
   REAL profangle = 0.;

   cout << "\nAirfoil angle [degrees]\n";
   cin >> profangle;
   options << profangle << std::endl;
   profangle *= M_PI/180.;

   TPZNACAXXXX profile(cord,digits,profangle,x0);

   NACAPoints(profile, nodes, elms, nSubdiv);

// Creating the geometric mesh
   TPZGeoMesh * gmesh = CreateNACAGeoMesh(cmesh->Reference(), profile, nodes, elms, EQuadrilateral, 1, gElem, nSubdiv);

   //TPZFlowCompMesh * cmesh = new TPZFlowCompMesh(gmesh);

// Creating the materials
   TPZEulerConsLaw * matp = new TPZEulerConsLaw(1/*nummat*/,
                                            0/*timeStep*/,
					    gamma /*gamma*/,
					    dim /* dim*/,
					    DiffType);
// Setting initial solution
   matp->SetForcingFunction(NULL);
   // Setting the time discretization method
   matp->SetTimeDiscr(Diff_TD,
                     ConvVol_TD,
		     ConvFace_TD);
   //mat->SetDelta(0.1); // Not necessary, since the artDiff
   // object computes the delta when it equals null.

   matp->SetCFL(CFL);
   matp->SetDelta(delta);

   TPZAutoPointer<TPZMaterial> mat(matp);
   cmesh -> InsertMaterialObject(mat);

// Boundary conditions

   TPZAutoPointer<TPZMaterial>  bc;
   TPZFMatrix val1(4,4), val2(4,1);

   //aresta interna NACA: Wall
   val1.Zero();
   val2.Zero();
   for( i = 0; i < 2*m; i++)
   {
      TPZGeoElBC((TPZGeoEl *)gElem[i],4,-1);
   }
   bc = mat->CreateBC(mat,-1,/*11*/5,val1,val2);
   cmesh->InsertMaterialObject(bc);

   REAL angle = 0.;

   cout << "\nMach number\n";
   cin >> Mach;
   options << Mach << std::endl;

   // leftmost bc face: Inlet
   val1.Zero();
   val2.Zero();
   val2(0,0) = 1.;// rho
   val2(1,0) = Mach * cos(angle);// Machx
   val2(2,0) = Mach * sin(angle);// Machy
   val2(3,0) = 1.;// pressure
   for( i = 0; i < l; i++)
   {
      TPZGeoElBC((TPZGeoEl *)gElem[(n-1)*2*m+i],6,-2);
      TPZGeoElBC((TPZGeoEl *)gElem[n*2*m-i-1]  ,6,-2);
   }
   bc = mat->CreateBC(mat,-2,10,val1,val2); // inflow
   cmesh->InsertMaterialObject(bc);


   // Material was -2 - directional
   // Apparently nothing changed...
   // upper and lower extern NACA BC faces: inflow/outflow
   for( i = (n-1)*2*m + l; i < n*2*m - l; i++)
   {
      TPZGeoElBC((TPZGeoEl *)gElem[i],6,-3);
   }

   int lastElement = 2*(m*n + p*n + (n-l) * l);
   // exit upper and bottom faces:  inflow/outflow
   for(i = 0; i < 2*(p+l); i++)
   {
      TPZGeoElBC((TPZGeoEl *)gElem[lastElement-i-1],6,-3);
   }

   bc = mat->CreateBC(mat,-3,9,val1,val2); // inflow/outflow
   cmesh->InsertMaterialObject(bc);

   cmesh->AutoBuild();
   //cmesh->AdjustBoundaryElements();

// printing meshes

   ofstream geoOut("geomesh.out");
   gmesh->Print(geoOut);
   geoOut.close();

   ofstream compOut("compmesh.out");
   cmesh->Print(compOut);
   compOut.close();

// generating initial guess for the mesh solution
   TPZFMatrix Solution = cmesh->Solution();
   Solution.Zero();

   int nVars = Solution.Rows();
   for(i = 0; i < nVars; i++)Solution(i,0) = 0;//.05;
   int j, NSolutionBlocks;
   //TPZBlock * pBlock = cmesh->Block();
   NSolutionBlocks = cmesh->Block().NBlocks();
   int nShape = Solution.Rows() / NSolutionBlocks / (dim + 2);
   int lastShapeFun = (nShape - 1)*(dim+2);
   for(j = 0; j < NSolutionBlocks; j++)
   {
      int blockOffset = cmesh->Block().Position(j) + lastShapeFun;

      REAL rho = 1.0,
           p = 1.,
	   u = sqrt(1.4 * p / rho) * Mach,
	   v = 0.,
	   vel2 = u*u + v*v;
      Solution(blockOffset  ,0) = rho;
      Solution(blockOffset+1,0) = rho * u;
      Solution(blockOffset+2,0) = rho * v;
      Solution(blockOffset+3,0) = p/(gamma-1.0) + 0.5 * rho * vel2;
   }

   cmesh->LoadSolution(Solution);

   return cmesh;
}
