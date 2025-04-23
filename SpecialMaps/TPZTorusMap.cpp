//
//  TPZTorusMap.cpp
//  PZ
//
//  Created by Philippe Devloo on 3/21/14.
//
//

#include "TPZTorusMap.h"
#include "tpzgeoelmapped.h"
#include "tpzgeoelrefpattern.h"

namespace pzgeom {

  template<class TGeo>
  int TPZTorusMap<TGeo>::ClassId() const{
    return Hash("TPZTorusMap") ^ TGeo::ClassId() << 1;
  }

  /// compute the corner coordinates of the corner nodes
    template<class TGeo>
    void TPZTorusMap<TGeo>::ComputeCornerCoordinates(TPZGeoMesh &gmesh)
    {

      auto TestX = [this](const REAL phi, const REAL theta,
                          const TPZVec<REAL> &codiff){
        TPZManVector<REAL,3> xcalc(3);
        xcalc[0] = (fR + fr*cos(theta))*cos(phi);
        xcalc[1] = (fR + fr*cos(theta))*sin(phi);
        xcalc[2] = fr*sin(theta);

        TPZManVector<REAL,3> diff = codiff;
        REAL normdiff{0};
        for(int ix = 0; ix < 3; ix++){
          diff[ix] -= xcalc[ix];
          normdiff += diff[ix]*diff[ix];
        }
        normdiff = sqrt(normdiff);
        return normdiff;
      };
      constexpr int nnodes = TGeo::NNodes;
      for (int in=0; in<nnodes; in++) {
        const int64_t nodeindex = TGeo::fNodeIndexes[in];
        TPZManVector<REAL,3> co(3), codiff(3);
        gmesh.NodeVec()[nodeindex].GetCoordinates(co);
        //now we need to test whether co lies in the torus surface
        codiff = co;
        //first we subtract the origin, the center of the torus
        for(int ix = 0; ix < 3; ix++){codiff[ix] -= fOrigin[ix];}
        //first guess for phi, between -pi/2 and pi/2
        REAL phi = atan2(codiff[1],codiff[0]);
        REAL theta = asin(codiff[2]/fr);
        //now we compute the point based on phi and theta and check if it
        //matches codiff
        auto normdiff = TestX(phi,theta,codiff);

        constexpr REAL eps = 1e-8;
        if(normdiff < eps){
          fPhiTheta(0,in) = phi;
          fPhiTheta(1,in) = theta;
          //next point
          continue;
        }

        //now we test for theta over pi_2, 3pi_2
        theta += M_PI;
        normdiff = TestX(phi,theta,codiff);

        if(normdiff < eps){
          fPhiTheta(0,in) = phi;
          fPhiTheta(1,in) = theta;
          //next point
          continue;
        }

        PZError<<__PRETTY_FUNCTION__
               <<"\nNode "<<in<<" with coords \n";
        for(int ix = 0; ix < 3; ix++){PZError<<co[ix]<<' ';}
        PZError<<"\nDoes not lie in the torus surface for radii "
               <<fr<<" and "<<fR
               <<"\nand origin: ";
          for(int ix = 0; ix < 3; ix++){PZError<<fOrigin[ix]<<' ';}
        PZError<<std::endl;
        DebugStop();
            
      }
    }

  template<class TGeo>
  void TPZTorusMap<TGeo>::InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size)
  {
    
    REAL R = 0.5, r = 0.2;
    TPZManVector<int64_t,4> indexes(TGeo::NNodes);
    TPZFNMatrix<TGeo::NNodes*3,REAL> phitheta(TGeo::NNodes,2);
    if constexpr (std::is_same_v<TGeo, pzgeom::TPZGeoTriangle>){

      phitheta = {
        {0,0},{M_PI/3,0},{0,M_PI}
      };
    }else{
      phitheta = {
        {0., 0},
        {M_PI/3, 0},
        {M_PI/3, M_PI},    
        {0., M_PI},
      };
    }
        
    for (int i=0; i<TGeo::NNodes; i++) {
      indexes[i] = gmesh.NodeVec().AllocateNewElement();
      const REAL phi = phitheta.Get(i,0);
      const REAL theta = phitheta.Get(i,1);
      TPZManVector<REAL,3> xco(3,0);
      xco[0] = (R + r*cos(theta))*cos(phi);
      xco[1] = (R + r*cos(theta))*sin(phi);
      xco[2] = r*sin(theta);

      for(int ix = 0; ix < 3; ix++){
        xco[ix]+=lowercorner[ix];
      }
      gmesh.NodeVec()[indexes[i]].Initialize(xco, gmesh);
    }
    TPZGeoElRefPattern<pzgeom::TPZTorusMap<TGeo>>
      *gel = new TPZGeoElRefPattern<pzgeom::TPZTorusMap<TGeo>>(indexes,matid,gmesh);
    gel->Geom().SetOrigin(lowercorner);
    gel->Geom().SetRadii(R, r);
    gel->Geom().ComputeCornerCoordinates(gmesh);
  }
};//namespace geom

#define IMPLEMENTTORUSMAP(TGEO)                                   \
  template class pzgeom::TPZTorusMap<TGEO>;                       \
  template class TPZGeoElRefLess<pzgeom::TPZTorusMap<TGEO> >;     \
  template class TPZGeoElRefPattern<pzgeom::TPZTorusMap<TGEO> >;

  IMPLEMENTTORUSMAP(pzgeom::TPZGeoTriangle)
  IMPLEMENTTORUSMAP(pzgeom::TPZGeoQuad)

#undef IMPLEMENTTORUSMAP