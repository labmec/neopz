//
// C++ Interface: tpzgeoelmapped
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef TPZGEOELMAPPED_H
#define TPZGEOELMAPPED_H

#include "pzmanvector.h"
#include "pzfmatrix.h"
#include "pzgeoel.h"
#include "pzgeoelside.h"
class TPZGeoMesh;

/**
This class implements a geometric element which uses its ancestral to compute its jacobian. Its main intent is to make the division of specially mapped elements easier: if the coarse grid map is consistent, then so will all refined meshes

@author Philippe R. B. Devloo
*/
template<class TGeo, class TShape, class TFather>
class TPZGeoElMapped : public TFather {
public:
  TPZGeoElMapped() : TFather(), fCornerCo(TShape::Dimension,TGeo::NNodes,0.)
  {
  }
  TPZGeoElMapped(int id,TPZVec<int> &nodeindexes,int matind,TPZGeoMesh &mesh) :
      TFather(id,nodeindexes,matind,mesh), fCornerCo(TShape::Dimension,TGeo::NNodes,0.)
      {
      }
  TPZGeoElMapped(TPZVec<int> &nodeindices,int matind,TPZGeoMesh &mesh) :
      TFather(nodeindices,matind,mesh), fCornerCo(TShape::Dimension,TGeo::NNodes,0.)
      {
      }
  TPZGeoElMapped(TPZVec<int> &nodeindices,int matind,TPZGeoMesh &mesh,int &index) :
      TFather(nodeindices,matind,mesh,index), fCornerCo(TShape::Dimension,TGeo::NNodes,0.)
      {
      }


  ~TPZGeoElMapped()
  {
  }

  /**Sets the father element index*/
  virtual void SetFather(int fatherindex)
  {
    TFather::SetFather(fatherindex);
    TPZGeoEl *father = TFather::Father();
    if(!father) return;
    TPZGeoEl *nextfather = 0;
    if(father) nextfather = father->Father();
    while(nextfather)
    {
      father = nextfather;
      nextfather = father->Father();
    }
    int in, nnodes = TGeo::NNodes;
    for(in=0; in<nnodes; in++)
    {
      TPZTransform tr = TShape::SideToSideTransform(in,TGeo::NSides-1);
      TPZManVector<REAL,TShape::Dimension> ptin(0),ptout(TShape::Dimension,0.);
      tr.Apply(ptin,ptout);
      int nfs = father->NSides();
      TPZGeoElSide thisside(this,TGeo::NSides-1);
      TPZGeoElSide ancestor(father,nfs-1);
      TPZTransform trfat(TShape::Dimension);
      TPZManVector<REAL,3> ptancestor(father->Dimension(),0.);
      thisside.SideTransform3(ancestor,trfat);
      trfat.Apply(ptout,ptancestor);
      int id;
      for(id=0; id<TShape::Dimension; id++)
      {
        fCornerCo(id,in) = ptancestor[id];
      }
    }
    
  }

  /**return the Jacobian matrix at the point*/
  virtual void Jacobian(TPZVec<REAL> &coordinate,TPZFMatrix &jac,TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv)
  {
    TPZGeoEl *father = TFather::Father();
    if(!father) return;
    TPZGeoEl *nextfather = 0;
    if(father) nextfather = father->Father();
    while(nextfather)
    {
      father = nextfather;
      nextfather = father->Father();
    }
    const int dim = TShape::Dimension;
    TPZManVector<REAL,3> ksibar(father->Dimension());
    TPZFNMatrix<dim*dim> jaclocal(dim,dim,0.),jacinvlocal(dim,dim,0.),jacfather(dim,dim,0.),
    jacinvfather(dim,dim,0.);
    TPZFNMatrix<9> axeslocal(3,3,0.),axesfinal(3,3,0.),axesfather(3,3,0.);
    REAL detjaclocal, detjacfather;
    TGeo::Jacobian(fCornerCo,coordinate,jaclocal,axeslocal,detjaclocal,jacinvlocal);
    axeslocal.Transpose();
    TGeo::X(fCornerCo,coordinate,ksibar);
    father->Jacobian(ksibar,jacfather,axesfather,detjacfather,jacinvfather);
    axeslocal.Multiply(axes,axesfinal);
    jacfather.Multiply(jaclocal,jac);
    jacinvlocal.Multiply(jacinvfather,jacinv);
    detjac = detjaclocal*detjacfather;
  }

  /**return the coordinate in real space of the point coordinate in the master element space*/
  virtual void X(TPZVec<REAL> &coordinate,TPZVec<REAL> &result)
  {
    TPZGeoEl *father = TFather::Father();
    if(!father) return;
    TPZGeoEl *nextfather = 0;
    if(father) nextfather = father->Father();
    while(nextfather)
    {
      father = nextfather;
      nextfather = father->Father();
    }

    TPZManVector<REAL,3> ksibar(father->Dimension());
    KsiBar(coordinate,ksibar);
    father->X(ksibar,result);
  }

private:
    TPZFNMatrix<TShape::Dimension*TGeo::NNodes> fCornerCo;
    
    /// compute the map of the point ksi to the ancestor ksibar and the gradient of the ancestor ksibar with respect to ksi
    void KsiBar(TPZVec<REAL> &ksi, TPZVec<REAL> &ksibar, TPZFMatrix &jac)
    {
      const int dim = TShape::Dimension;
      TPZFNMatrix<TGeo::NNodes> phi(TGeo::NNodes,1,0.);
      TPZFNMatrix<dim*TGeo::NNodes> dphi(dim,TGeo::NNodes,0.);
      TGeo::Shape(ksi,phi,dphi);
      jac.Redim(dim,dim);
      ksibar.Fill(0.);
      int in,id,jd;
      for(in=0; in<TGeo::NNodes; in++)
      {
        for(id=0; id<dim; id++)
        {
          ksibar[id] += phi(in,0)*fCornerCo(id,in);
          for(jd=0; jd<dim; jd++)
          {
            jac(id,jd) += dphi(jd,in)*fCornerCo(id,in);
          }
        }
      }
    }
    
    /// compute the map of the point ksi to the ancestor ksibar
    void KsiBar(TPZVec<REAL> &ksi, TPZVec<REAL> &ksibar)
    {
      const int dim = TShape::Dimension;
      TPZFNMatrix<TGeo::NNodes> phi(TGeo::NNodes,1,0.);
      TPZFNMatrix<dim*dim> jac(dim,dim,0.);
      TPZFNMatrix<dim*TGeo::NNodes> dphi(dim,TGeo::NNodes,0.);
      TGeo::Shape(ksi,phi,dphi);
      ksibar.Fill(0.);
      int in,id;
      for(in=0; in<TGeo::NNodes; in++)
      {
        for(id=0; id<dim; id++)
        {
          ksibar[id] += phi(in,0)*fCornerCo(id,in);
        }
      }
    }

};

#endif
