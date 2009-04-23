#include "tpzgeoblend.h"

#include "pzgeoelside.h"
#include "tpzgeoelmapped.h"

#include "pzlog.h"
#include <sstream>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.geoblend"));
#endif

template <class TGeo>
bool TPZGeoBlend<TGeo>::MapToNeighSide(int side, TPZVec<REAL> &InternalPar, TPZVec<REAL> &NeighPar, TPZFMatrix &JacNeighSide)
{
     int SideDim    = fNeighbours[side-TGeo::NNodes].Dimension();
     TPZFMatrix JacSide;

     TPZVec< REAL > SidePar(SideDim);
     const bool check = this->MapToSide(side, InternalPar, SidePar,JacSide);

#ifdef LOG4CXX
  if(logger->isDebugEnabled())
  {
    std::stringstream sout;
    sout << "side " << side << std::endl;
    sout << "InternalPar: ";
    for(int i = 0; i < InternalPar.NElements(); i++) sout << InternalPar[i] << "\t";
    sout << "\n";
    
    sout << "SidePar: ";
    for(int i = 0; i < SidePar.NElements(); i++) sout << SidePar[i] << "\t";
    sout << "\n";

    JacSide.Print("JacSide = ",sout);
    LOGPZ_DEBUG(logger,sout.str())
  }
#endif

     if(!check){
		  return false;
	   }

     NeighPar.Resize(SideDim);
     TransfBetweenNeigh(side).Apply(SidePar,NeighPar);

     JacNeighSide.Resize(TransfBetweenNeigh(side).Mult().Rows(),JacSide.Cols());
     TransfBetweenNeigh(side).Mult().Multiply(JacSide,JacNeighSide);

#ifdef LOG4CXX
  if(logger->isDebugEnabled())
  {
    std::stringstream sout;

    sout << "NeighPar: ";
    for(int i = 0; i < NeighPar.NElements(); i++) sout << NeighPar[i] << "\t";
    sout << "\n";

    JacNeighSide.Print("JacNeighSide = ",sout);
    LOGPZ_DEBUG(logger,sout.str())
  }
#endif

	return true;
}

template <class TGeo>
void TPZGeoBlend<TGeo>::X(TPZFMatrix & coord, TPZVec<REAL>& par, TPZVec<REAL> &result)
{
    ///Criando uma cópia das coordenadas dos nós
    TPZManVector< TPZManVector<REAL,3>,27 > NodeCoord(TGeo::NNodes);
    for(int i = 0; i < TGeo::NNodes; i++)
    {
        NodeCoord[i].Resize(3);
        for(int j = 0; j < 3; j++)
        {
            NodeCoord[i][j] = coord(j,i);
        }
    }
#ifdef LOG4CXX
	if(logger->isDebugEnabled())
	{
		std::stringstream sout;
		sout << "input parameter par " << par << std::endl;
		sout << "NodeCoord " << NodeCoord;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	
    ///
    result.Fill(0.);
    TPZManVector<REAL,3> NeighPar, SidePar, Xside(3,0.);

    int majorSide = TGeo::NSides - 1;

    TPZManVector<REAL,27> SidesCounter(TGeo::NSides,0);
    TPZStack<int> LowNodeSides, LowAllSides;

    TPZFNMatrix<9> blend(TGeo::NNodes,1), Dblend(TGeo::Dimension,TGeo::NNodes), NotUsedHere;
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

            if(!MapToNeighSide(byside,par,NeighPar,NotUsedHere))
			{
#ifdef LOG4CXX
				if(logger->isDebugEnabled())
				{
					std::stringstream sout;
					sout << "MapToNeighSide is singular for par " << par << " and side " << byside << " skipping the side ";
					LOGPZ_DEBUG(logger,sout.str())
				}
#endif
				continue;
			}
            Neighbour(byside).X(NeighPar,Xside);
#ifdef LOG4CXX
			if(logger->isDebugEnabled())
			{
				std::stringstream sout;
				sout << "NeighPar " << NeighPar << ' ';
				sout << "Xside " << Xside << ' ';
				sout << "blendTemp " << blendTemp;
				LOGPZ_DEBUG(logger,sout.str())
			}
#endif
			
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
#ifdef LOG4CXX
	if(logger->isDebugEnabled())
	{
		std::stringstream sout;
		sout << "sidescounter before contributing linear map " << SidesCounter;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	
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
#ifdef LOG4CXX
	if(logger->isDebugEnabled())
	{
		std::stringstream sout;
		sout << "result " << result;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	
}
/////////////////////////////////////////////////////////////////////////////////////////////

template <class TGeo>
void TPZGeoBlend<TGeo>::Jacobian(TPZFMatrix & coord, TPZVec<REAL>& par, TPZFMatrix &jacobian, TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv)
{
    TPZManVector<REAL,3> NeighPar, SidePar, Xside(3,0.), XNode(3,0.);
    int majorSide = TGeo::NSides - 1;

    TPZManVector<REAL> SidesCounter(TGeo::NSides,0);
    TPZStack<int> LowNodeSides, LowAllSides;
#ifdef LOG4CXX
	if(logger->isDebugEnabled())
	{
		std::stringstream sout;
		sout << "input parameter par " << par;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
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
            if(!MapToNeighSide(byside,par,NeighPar, Jneighbourhood))
			{
				continue;
			}
            Neighbour(byside).X(NeighPar,Xside);
            Neighbour(byside).Jacobian(NeighPar,J1,Ax,Det,Inv);
            Ax.Transpose(); 
            Ax.Multiply(J1,J2);
#ifdef LOG4CXX
			if(logger->isDebugEnabled())
			{
				std::stringstream sout;
				sout << "byside " << byside << std::endl;
        sout << "side of the neighbour " << Neighbour(byside).Side() << std::endl;
        Neighbour(byside).Element()->Print(sout);
				sout << "neighbour parameter(NeighPar) " << NeighPar << std::endl;
				sout << "Jacobian of the map(Jneighborhood) " << Jneighbourhood << std::endl;
				sout << "Xside " << Xside << std::endl;
				sout << "jacobian neighbour(J1) " << J1 << std::endl;
        Ax.Print("Ax of the neighbour (Ax) ",sout);        
				sout << "jacobian of the neighbour multiplied by the axes(J2) " << J2 << std::endl;
				LOGPZ_DEBUG(logger,sout.str())
			}
#endif
            J2.Multiply(Jneighbourhood,J1);
#ifdef LOG4CXX
			if(logger->isDebugEnabled())
			{
				std::stringstream sout;
				sout << "acumulated jacobian(J1) " << J1 << std::endl;
				LOGPZ_DEBUG(logger,sout.str())
			}
#endif
            REAL blendTemp = 0.; 
            TPZManVector<REAL,3> DblendTemp(TGeo::Dimension,0.);
            for(int a = 0; a < LowNodeSides.NElements(); a++)
            {
                    ///Mapeando os Nós ( para contemplar casos em que CoordNó != X(CoordNó) )
                    TPZVec<REAL> coordTemp(3);
                    TGeo::CenterPoint(LowNodeSides[a],coordTemp);
                    if(!(MapToNeighSide(byside,coordTemp,NeighPar,NotUsedHere)))
					{
						cout << __PRETTY_FUNCTION__ << " This should never happen\n";
						DebugStop();
					}
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
#ifdef LOG4CXX
	{
		std::stringstream sout;
		JacTemp.Print("Jabobian before contributing the nodes",sout);
		sout << "SidesCounter " << SidesCounter << std::endl;
		sout << "DBlend " << Dblend << std::endl;
		sout << "NodeCoord " << NodeCoord << std::endl;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
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

    if(TGeo::Dimension == 1)
    {
      detjac = jacobian(0,0);
      jacinv(0,0) = 1./detjac;
    }
    else if(TGeo::Dimension == 2)
    {
        detjac = jacobian(0,0)*jacobian(1,1) - jacobian(1,0)*jacobian(0,1);
        jacinv(0,0) =  jacobian(1,1) / detjac;
        jacinv(1,1) =  jacobian(0,0) / detjac;
        jacinv(0,1) = -jacobian(0,1) / detjac;
        jacinv(1,0) = -jacobian(1,0) / detjac;
    }
    else if(TGeo::Dimension == 3)
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
  for(int byside = TGeo::NNodes; byside < (TGeo::NSides); byside++)
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
  if(side < 0 || side >= TGeo::NSides-1)
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

  TPZGeoEl *newel = mesh->CreateGeoBlendElement(type, nodeindices, bc, index);
  TPZGeoElSide me(orig,side);
  TPZGeoElSide newelside(newel,newel->NSides()-1);

  newelside.InsertConnectivity(me);
  newel->Initialize();

  return newel;
}

template class TPZGeoBlend<TPZGeoCube>;
template class TPZGeoBlend<TPZGeoTriangle>;
template class TPZGeoBlend<TPZGeoPrism>;
template class TPZGeoBlend<TPZGeoPyramid>;
template class TPZGeoBlend<TPZGeoTetrahedra>;
template class TPZGeoBlend<TPZGeoQuad>;
template class TPZGeoBlend<TPZGeoLinear>;
template class TPZGeoBlend<TPZGeoPoint>;

#ifndef WIN32
#include "pzgeoelrefless.h.h"
#include "tpzgeoelrefpattern.h.h"
#endif

///CreateGeoElement -> TPZGeoBlend
#define IMPLEMENTBLEND(TGEO,CLASSID,CREATEFUNCTION) \
template< > \
TPZGeoEl *TPZGeoElRefLess<TPZGeoBlend<TGEO> >::CreateGeoElement(MElementType type, TPZVec<int>& nodeindexes, int matid, int& index) \
{ \
TPZGeoMesh &mesh = *(this->Mesh()); \
if(!&mesh) return 0; \
return CreateGeoElementMapped(mesh,type,nodeindexes,matid,index); \
} \
\
template<> \
int TPZGeoElRefPattern<TPZGeoBlend<TGEO>  >::ClassId() const { \
return CLASSID; \
} \
template class \
TPZRestoreClass< TPZGeoElRefPattern<TPZGeoBlend<TGEO> >, CLASSID>; \
\
template<> \
TPZCompEl *(*TPZGeoElRefLess<TPZGeoBlend<TGEO> >::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CREATEFUNCTION; \
\
template class TPZGeoElRefLess<TPZGeoBlend<TGEO> >;


IMPLEMENTBLEND(pzgeom::TPZGeoPoint,TPZGEOBLENDPOINTID,CreatePointEl)
IMPLEMENTBLEND(pzgeom::TPZGeoLinear,TPZGEOBLENDLINEARID,CreateLinearEl)
IMPLEMENTBLEND(pzgeom::TPZGeoQuad,TPZGEOBLENDQUADID,CreateQuadEl)
IMPLEMENTBLEND(pzgeom::TPZGeoTriangle,TPZGEOBLENDTRIANGLEID,CreateTriangleEl)
IMPLEMENTBLEND(pzgeom::TPZGeoCube,TPZGEOBLENDCUBEID,CreateCubeEl)
IMPLEMENTBLEND(pzgeom::TPZGeoPrism,TPZGEOBLENDPRISMID,CreatePrismEl)
IMPLEMENTBLEND(pzgeom::TPZGeoPyramid,TPZGEOBLENDPYRAMIDID,CreatePyramEl)
IMPLEMENTBLEND(pzgeom::TPZGeoTetrahedra,TPZGEOBLENDTETRAHEDRAID,CreateTetraEl)

#include "pznoderep.h.h"
template class pzgeom::TPZNodeRep<8,TPZGeoBlend<TPZGeoCube> >;
template class pzgeom::TPZNodeRep<6,TPZGeoBlend<TPZGeoPrism> >;


