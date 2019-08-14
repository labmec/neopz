#include "tpzgeoblend.h"

#include "pzgeoelside.h"
#include "tpzgeoelmapped.h"

#include "pzlog.h"
#include <sstream>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.caju.geoblend"));
#endif

template <class TGeo>
void TPZGeoBlend<TGeo>::MapToNeighSide(int side, TPZVec<REAL> &InternalPar, TPZVec<REAL> &NeighPar, TPZFMatrix &JacNeighSide)
{
     int SideDim    = fNeighbours[side-TGeo::NNodes].Dimension();
     TPZFMatrix JacSide;

     TPZVec< REAL > SidePar(SideDim);
     this->MapToSide(side, InternalPar, SidePar,JacSide);

     NeighPar.Resize(SideDim);
     TransfBetweenNeigh(side).Apply(SidePar,NeighPar);

     JacNeighSide.Resize(TransfBetweenNeigh(side).Mult().Rows(),JacSide.Cols());
     TransfBetweenNeigh(side).Mult().Multiply(JacSide,JacNeighSide);
}

template <class TGeo>
void TPZGeoBlend<TGeo>::X(TPZFMatrix & coord, TPZVec<REAL>& par, TPZVec<REAL> &result)
{
    ///Criando uma cópia das coordenadas dos nós
    TPZVec< TPZVec<REAL> > NodeCoord(TGeo::NNodes);
    for(int i = 0; i < TGeo::NNodes; i++)
    {
        NodeCoord[i].Resize(3);
        for(int j = 0; j < 3; j++)
        {
            NodeCoord[i][j] = coord(j,i);
        }
    }
    ///
    result.Fill(0.);
    TPZVec<REAL> NeighPar, SidePar, Xside(3,0.);

    int majorSide = TGeo::NSides - 1;

    TPZVec<REAL> SidesCounter(majorSide,0);
    TPZStack<int> LowNodeSides, LowAllSides;

    TPZFMatrix blend(TGeo::NNodes,1), Dblend(TGeo::Dimension,TGeo::NNodes), NotUsedHere;
    TGeo::Shape(par,blend,Dblend);

    for(int byside = majorSide; byside >= TGeo::NNodes; byside--)
    {
        if(fNeighbours[byside-TGeo::NNodes].Exists())
        {
            TGeo::LowerDimensionSides(byside,LowNodeSides,0);
            TGeo::LowerDimensionSides(byside,LowAllSides);
            REAL blendTemp = 0.;
            for(int a = 0; a < LowNodeSides.NElements(); a++)
            {
                ///Mapeando os Nós ( para contemplar casos em que CoordNó != X(CoordNó) )
                TPZVec<REAL> coordTemp(3);
                TGeo::CenterPoint(LowNodeSides[a],coordTemp);
                MapToNeighSide(byside,coordTemp,NeighPar,NotUsedHere);
                Neighbour(byside).X(NeighPar,Xside);
                NodeCoord[LowNodeSides[a]] = Xside;
                ///
                blendTemp += blend(LowNodeSides[a],0);
            }

            MapToNeighSide(byside,par,NeighPar,NotUsedHere);
            Neighbour(byside).X(NeighPar,Xside);

            for(int c = 0; c < 3; c++)
            {
                result[c] += (1 - SidesCounter[byside]) * Xside[c]*blendTemp;
            }

            for(int b = 0; b < LowAllSides.NElements(); b++)
            {
                SidesCounter[LowAllSides[b]] += (1 - SidesCounter[byside]);
            }
        }
    }

    for(int a = 0; a < TGeo::NNodes; a++)
    {
        for(int b = 0; b < 3; b++)
        {
                /// Aqui usa-se os Nós mapeados --------.
                /// (ao inves de coord(a,b))            |
                ///                                     v
            result[b] += (1 - SidesCounter[a]) * NodeCoord[a][b]*blend(a,0);
        }
    }
}
/////////////////////////////////////////////////////////////////////////////////////////////

template <class TGeo>
void TPZGeoBlend<TGeo>::Jacobian(TPZFMatrix & coord, TPZVec<REAL>& par, TPZFMatrix &jacobian, TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv)
{
    TPZManVector<REAL,3> NeighPar, SidePar, Xside(3,0.), XNode(3,0.);
    int majorSide = TGeo::NSides - 1;

    TPZManVector<REAL> SidesCounter(majorSide,0);
    TPZStack<int> LowNodeSides, LowAllSides;

    TPZFNMatrix<24> blend(TGeo::NNodes,1), Dblend(TGeo::Dimension,TGeo::NNodes), NotUsedHere;
    TGeo::Shape(par,blend,Dblend);

    ///Criando uma cópia das coordenadas dos nós
    TPZVec< TPZVec<REAL> > NodeCoord(TGeo::NNodes);
    for(int i = 0; i < TGeo::NNodes; i++)
    {
        NodeCoord[i].Resize(3);
        for(int j = 0; j < 3; j++)
        {
            NodeCoord[i][j] = coord(j,i);
        }
    }
    ///

    TPZFNMatrix<9> J1, J2, Ax, JacTemp(3,TGeo::Dimension, 0.), Jneighbourhood; 
    REAL Det;
    for(int byside = majorSide; byside >= TGeo::NNodes; byside--)
    {
        if(Neighbour(byside).Exists())
        {
            TGeo::LowerDimensionSides(byside,LowNodeSides,0);
            TGeo::LowerDimensionSides(byside,LowAllSides);
            TPZFNMatrix<9> Inv(fNeighbours[byside-TGeo::NNodes].Dimension(),fNeighbours[byside-TGeo::NNodes].Dimension());
            MapToNeighSide(byside,par,NeighPar, Jneighbourhood);
            Neighbour(byside).X(NeighPar,Xside);
            Neighbour(byside).Jacobian(NeighPar,J1,Ax,Det,Inv);
            Ax.Transpose(); 
            Ax.Multiply(J1,J2);
            J2.Multiply(Jneighbourhood,J1);
            REAL blendTemp = 0.; 
            TPZManVector<REAL,3> DblendTemp(TGeo::Dimension,0.);
            for(int a = 0; a < LowNodeSides.NElements(); a++)
            {
                    ///Mapeando os Nós ( para contemplar casos em que CoordNó != X(CoordNó) )
                    TPZVec<REAL> coordTemp(3);
                    TGeo::CenterPoint(LowNodeSides[a],coordTemp);
                    MapToNeighSide(byside,coordTemp,NeighPar,NotUsedHere);
                    Neighbour(byside).X(NeighPar,XNode);
                    NodeCoord[LowNodeSides[a]] = XNode;
                    ///
                    TPZVec<REAL> parChanged(par.NElements());
                    TGeo::FixSingularity(byside,par,parChanged);
                    TGeo::Shape(parChanged,blend,Dblend);

                    blendTemp += blend(LowNodeSides[a],0);
                    for(int b = 0; b < TGeo::Dimension; b++) 
                    {
                        DblendTemp[b] += Dblend(b,LowNodeSides[a]);
                    }
            }

            for(int a = 0; a < 3; a++)
            {
                    for(int b = 0; b < TGeo::Dimension; b++) 
                    {
                        JacTemp(a,b) += (1 - SidesCounter[byside]) * (J1(a,b)*blendTemp + Xside[a]*DblendTemp[b]);
                    }
            }
            for(int a = 0; a < LowAllSides.NElements(); a++) 
            {
                SidesCounter[LowAllSides[a]] += (1 - SidesCounter[byside]);
            }
        }
    }

    for(int a = 0; a < TGeo::NNodes; a++)
    {
        for(int b = 0; b < 3; b++) 
        { 
            for(int c = 0; c < TGeo::Dimension; c++)
            { 
                /// Aqui usa-se os Nós mapeados --------------.
                /// (ao inves de coord(a,b))                  |
                ///                                           v
                JacTemp(b,c) += (1 - SidesCounter[a]) * NodeCoord[a][b]*Dblend(c,a);
            }
        }
    }
    TPZFMatrix axest;
    JacTemp.GramSchmidt(axest,jacobian);
    axest.Transpose(&axes);

    if(TGeo::Dimension == 2)
    {
        detjac = jacobian(0,0)*jacobian(1,1) - jacobian(1,0)*jacobian(0,1);
        jacinv(0,0) =  jacobian(1,1) / detjac;
        jacinv(1,1) =  jacobian(0,0) / detjac;
        jacinv(0,1) = -jacobian(0,1) / detjac;
        jacinv(1,0) = -jacobian(1,0) / detjac;
    }
    if(TGeo::Dimension == 3)
    {
        detjac = -jacobian(0,2)*jacobian(1,1)*jacobian(2,0);//- a02 a11 a20
        detjac += jacobian(0,1)*jacobian(1,2)*jacobian(2,0);//+ a01 a12 a20
        detjac += jacobian(0,2)*jacobian(1,0)*jacobian(2,1);//+ a02 a10 a21
        detjac -= jacobian(0,0)*jacobian(1,2)*jacobian(2,1);//- a00 a12 a21
        detjac -= jacobian(0,1)*jacobian(1,0)*jacobian(2,2);//- a01 a10 a22
        detjac += jacobian(0,0)*jacobian(1,1)*jacobian(2,2);//+ a00 a11 a22

        jacinv(0,0) = (-jacobian(1,2)*jacobian(2,1)+jacobian(1,1)*jacobian(2,2)) / detjac;//-a12 a21 + a11 a22
        jacinv(0,1) = ( jacobian(0,2)*jacobian(2,1)-jacobian(0,1)*jacobian(2,2)) / detjac;// a02 a21 - a01 a22
        jacinv(0,2) = (-jacobian(0,2)*jacobian(1,1)+jacobian(0,1)*jacobian(1,2)) / detjac;//-a02 a11 + a01 a12
        jacinv(1,0) = ( jacobian(1,2)*jacobian(2,0)-jacobian(1,0)*jacobian(2,2)) / detjac;// a12 a20 - a10 a22
        jacinv(1,1) = (-jacobian(0,2)*jacobian(2,0)+jacobian(0,0)*jacobian(2,2)) / detjac;//-a02 a20 + a00 a22
        jacinv(1,2) = ( jacobian(0,2)*jacobian(1,0)-jacobian(0,0)*jacobian(1,2)) / detjac;// a02 a10 - a00 a12
        jacinv(2,0) = (-jacobian(1,1)*jacobian(2,0)+jacobian(1,0)*jacobian(2,1)) / detjac;//-a11 a20 + a10 a21
        jacinv(2,1) = ( jacobian(0,1)*jacobian(2,0)-jacobian(0,0)*jacobian(2,1)) / detjac;// a01 a20 - a00 a21
        jacinv(2,2) = (-jacobian(0,1)*jacobian(1,0)+jacobian(0,0)*jacobian(1,1)) / detjac;//-a01 a10 + a00 a11
    }
}

template <class TGeo>
void TPZGeoBlend<TGeo>::Print(std::ostream &out)
{
  TGeo::Print(out);
  out << "Neighbours/transformations used for mapping the sides :\n";
  int is;
  for(is=TGeo::NNodes; is<TGeo::NSides-1; is++)
  {
    if(fNeighbours[is-TGeo::NNodes].Element())
    {
          out << "Side: " << is << " El/side: " << fNeighbours[is-TGeo::NNodes].Element()->Index() << ":" <<
          fNeighbours[is-TGeo::NNodes].Side() << '\n';
    }
  }
}

template <class TGeo>
void TPZGeoBlend<TGeo>::Initialize(TPZGeoEl *refel)
{
  for(int byside = TGeo::NNodes; byside < (TGeo::NSides - 1); byside++)
  {
    TPZGeoElSide ElemSide(refel,byside);
    TPZGeoElSide NextSide(refel,byside);
    while(NextSide.Neighbour().Element() != ElemSide.Element())
    {
      if(NextSide.Neighbour().Exists() && !NextSide.Neighbour().Element()->IsLinearMapping() && !NextSide.Neighbour().Element()->IsGeoBlendEl())
      {
        TPZGeoElSide NeighSide = NextSide.Neighbour();
        TPZTransform NeighTransf(NeighSide.Dimension(),NeighSide.Dimension());
        ElemSide.SideTransform3(NeighSide,NeighTransf);
        SetNeighbourInfo(byside,NeighSide,NeighTransf);
        break;
      }
      NextSide = NextSide.Neighbour();
    }
  }
}

///inseri este método pois o compilador não encontrou
template <class TGeo>
void TPZGeoBlend<TGeo>::Initialize(TPZVec<int> &nodeindexes, TPZGeoMesh &mesh)
{
  TGeo::Initialize(nodeindexes, mesh);
}

#include "TPZGeoCube.h"
#include "TPZGeoLinear.h"
#include "pzgeoquad.h"
#include "pzgeotriangle.h"
#include "pzgeoprism.h"
#include "pzgeotetrahedra.h"
#include "pzgeopyramid.h"
#include "pzgeopoint.h"
#include "pzcmesh.h"


using namespace pzgeom;

template <class TGeo>
TPZGeoEl *TPZGeoBlend<TGeo>::CreateBCGeoEl(TPZGeoEl *orig, int side,int bc)
{
  TPZStack<int> LowAllSides;
  TGeo::LowerDimensionSides(side,LowAllSides);
  if(side < 0 || side > TGeo::NSides-1)
  {
    return 0;
  }
  bool straight = true;
  for(int lowside = 0; lowside < LowAllSides.NElements(); lowside++)
  {
    if(LowAllSides[lowside] >= TGeo::NNodes && fNeighbours[LowAllSides[lowside]-TGeo::NNodes].Element())
    {
      straight = false;
    }
  }
  if(straight)
  {
      return TGeo::CreateBCGeoEl(orig,side,bc);
  }
  else
  {
    TPZGeoEl *newel = CreateBCGeoBlendEl(orig,side,bc);
    return newel;
  }

}

template <class TGeo>
TPZGeoEl *TPZGeoBlend<TGeo>::CreateBCGeoBlendEl(TPZGeoEl *orig,int side,int bc)
{
  int ns = orig->NSideNodes(side);
  TPZManVector<int> nodeindices(ns);
  int in;
  for(in=0; in<ns; in++)
  {
    nodeindices[in] = orig->SideNodeIndex(side,in);
  }
  int index;

  TPZGeoMesh *mesh = orig->Mesh();
  MElementType type = orig->Type(side);

  TPZGeoEl *newel = CreateGeoBlendElement(*mesh, type, nodeindices, bc, index);
  TPZGeoElSide me(orig,side);
  TPZGeoElSide newelside(newel,newel->NSides()-1);

  newelside.InsertConnectivity(me);
  newel->Initialize();

  return newel;
}

TPZGeoEl *CreateGeoBlendElement(TPZGeoMesh &mesh, MElementType type, TPZVec<int>& nodeindexes, int matid, int& index)
{
  if(!&mesh) return 0;

  
  switch( type ){
    case 0://point
    {
      TPZGeoEl *gel = new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoPoint> > (nodeindexes,matid,mesh);
      return gel;
    }
    case 1://line
    {
      TPZGeoEl *gel = new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoLinear> > (nodeindexes,matid,mesh);
      return gel;
    }
    case 2://triangle
    {
      TPZGeoEl *gel = new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoTriangle> > (nodeindexes,matid,mesh);
      return gel;
    }
    case 3://quadrilateral
    {
      TPZGeoEl *gel = new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoQuad> > (nodeindexes,matid,mesh);
      return gel;
    }
    case 4://tetraedron
    {
      TPZGeoEl *gel = new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoTetrahedra> > (nodeindexes,matid,mesh);
      return gel;
    }
    case 5://pyramid
    {
      TPZGeoEl *gel = new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoPyramid> > (nodeindexes,matid,mesh);
      return gel;
    }
    case 6://prism
    {
      TPZGeoEl *gel = new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoPrism> > (nodeindexes,matid,mesh);
      return gel;
    }
    case 7://cube
    {
      TPZGeoEl *gel = new TPZGeoElRefPattern<TPZGeoBlend<pzgeom::TPZGeoCube> > (nodeindexes,matid,mesh);
      return gel;
    }
    default:
    {
      PZError << "TPZGeoMesh::CreateGeoElementRefPattern type element not exists:"
          << " type = " << type << std::endl;
      return NULL;
    }
  }
}

template class TPZGeoBlend<TPZGeoCube>;
template class TPZGeoBlend<TPZGeoTriangle>;
template class TPZGeoBlend<TPZGeoPrism>;
template class TPZGeoBlend<TPZGeoPyramid>;
template class TPZGeoBlend<TPZGeoTetrahedra>;
template class TPZGeoBlend<TPZGeoQuad>;
template class TPZGeoBlend<TPZGeoLinear>;
template class TPZGeoBlend<TPZGeoPoint>;
