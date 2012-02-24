/**
 * @file
 * @brief Contains the implementation of the TPZGeoBlend methods. 
 */
#include "tpzgeoblend.h"

#include "pzgeoelside.h"
#include "tpzgeoelmapped.h"
#include "pzgeoelrefless.h"
#include "tpzgeoelrefpattern.h"

#include "pzlog.h"
#include <sstream>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.geoblend"));
#endif


template <class TGeo>
void pzgeom::TPZGeoBlend<TGeo>::SetNeighbourInfo(int side, TPZGeoElSide &neigh, TPZTransform &trans)
{
	if(!fNeighbours[side-TGeo::NNodes].ElementIndex() != -1)
	{
		fNeighbours[side-TGeo::NNodes] = neigh;
		fTrans[side - TGeo::NNodes] = trans;
	}
	else
	{
#ifdef LOG4CXX
		std::stringstream mess;
		mess << "Trying to SetNeighbourInfo for an already set element\n";
		mess << "* this * = " << __PRETTY_FUNCTION__ << "\n";
		this->Print(mess);
		mess << "* neigh * = \n";
		neigh.Element()->Print(mess);
		LOGPZ_DEBUG(logger,mess.str());
#endif
	}
}

template <class TGeo>
bool pzgeom::TPZGeoBlend<TGeo>::MapToNeighSide(int side, int SideDim, TPZVec<REAL> &InternalPar, TPZVec<REAL> &NeighPar, TPZFMatrix &JacNeighSide)
{
	TPZFNMatrix<9> JacSide;
	
	TPZManVector< REAL, 3 > SidePar(SideDim);
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
	
	if(!check)
	{
	 	return false;
	}
	
	NeighPar.Resize(SideDim);
	TransfBetweenNeigh(side).Apply(SidePar,NeighPar);
	
	JacNeighSide.Resize(0, 0);
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
void pzgeom::TPZGeoBlend<TGeo>::X(const TPZGeoEl &gel, TPZVec<REAL>& par, TPZVec<REAL> &result)
{
    TPZFNMatrix<45> coord(3,TGeo::NNodes);
    this->CornerCoordinates(gel,coord);
    
#ifdef LOG4CXX
	if(logger->isDebugEnabled())
	{
		std::stringstream sout;
		sout << "input parameter par " << par << std::endl;
		sout << "NodeCoord " << coord;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	
    result.Fill(0.);
    TPZManVector<REAL,3> NeighPar, SidePar, Xside(3,0.);
	
    int majorSide = TGeo::NSides - 1;
	
    TPZManVector<REAL,27> SidesCounter(TGeo::NSides,0);
    TPZStack<int> LowNodeSides, LowAllSides;
	
    TPZFNMatrix<9> blend(TGeo::NNodes,1), Dblend(TGeo::Dimension,TGeo::NNodes), NotUsedHere;
    TGeo::Shape(par,blend,Dblend);
    TPZGeoMesh *gmesh = gel.Mesh();
	
    for(int byside = majorSide; byside >= TGeo::NNodes; byside--)
    {
        TPZGeoElSide gelside(fNeighbours[byside-TGeo::NNodes],gmesh);
        if(gelside.Exists())
        {
            TGeo::LowerDimensionSides(byside,LowNodeSides,0);
            TGeo::LowerDimensionSides(byside,LowAllSides);
            REAL blendTemp = 0.;
            for(int a = 0; a < LowNodeSides.NElements(); a++)
            {
                blendTemp += blend(LowNodeSides[a],0);
            }
            int sidedim = gelside.Dimension();
			
            if(!MapToNeighSide(byside,sidedim,par,NeighPar,NotUsedHere))
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
            
            Neighbour(byside,gmesh).X(NeighPar,Xside);
			
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
            result[b] += (1 - SidesCounter[a]) * coord(b,a)*blend(a,0);
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


template <class TGeo>
void pzgeom::TPZGeoBlend<TGeo>::Jacobian(const TPZGeoEl &gel, TPZVec<REAL>& par, TPZFMatrix &jacobian, TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv)
{
    TPZFNMatrix<45> coord(3,TGeo::NNodes);
    this->CornerCoordinates(gel,coord);

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
	
    TPZFNMatrix<9> J1, J2, Ax, JacTemp(3,TGeo::Dimension, 0.), Jneighbourhood;
    REAL Det;
    TPZGeoMesh *gmesh = gel.Mesh();
    for(int byside = majorSide; byside >= TGeo::NNodes; byside--)
    {
        TPZGeoElSide neighbyside = Neighbour(byside, gmesh);
        if(neighbyside.Exists())
        {
            TGeo::LowerDimensionSides(byside,LowNodeSides,0);
            TGeo::LowerDimensionSides(byside,LowAllSides);
            int dim = Neighbour(byside, gmesh).Dimension();
            TPZFNMatrix<9> Inv(dim,dim);
            int sidedim = neighbyside.Dimension();
            if(!MapToNeighSide(byside,sidedim,par,NeighPar, Jneighbourhood))
			{
				continue;
			}
            Neighbour(byside,gmesh).X(NeighPar,Xside);
            Neighbour(byside,gmesh).Jacobian(NeighPar,J1,Ax,Det,Inv);
            Ax.Transpose(); 
            Ax.Multiply(J1,J2);
			
#ifdef LOG4CXX
			if(logger->isDebugEnabled())
			{
				std::stringstream sout;
				sout << "byside " << byside << std::endl;
        		sout << "side of the neighbour " << Neighbour(byside,gmesh).Side() << std::endl;
        		Neighbour(byside,gmesh).Element()->Print(sout);
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
                TPZManVector<REAL> parChanged(par.NElements());
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
		sout << "NodeCoord " << coord << std::endl;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	
    for(int a = 0; a < TGeo::NNodes; a++)
    {
        for(int b = 0; b < 3; b++) 
        { 
            for(int c = 0; c < TGeo::Dimension; c++)
            {
            	JacTemp(b,c) += (1 - SidesCounter[a]) * coord(b,a)*Dblend(c,a);
            }
        }
    }
	
    if(TGeo::Dimension == 1)
    {
        TPZFNMatrix<9> axest;
        JacTemp.GramSchmidt(axest,jacobian);
        axest.Transpose(&axes);
        
        detjac = jacobian(0,0);
        if(IsZero(detjac)){
            detjac = ZeroTolerance();
        }
        jacinv(0,0) = 1./detjac;
    } else if(TGeo::Dimension == 2)
    {
        TPZFNMatrix<9> axest;
        JacTemp.GramSchmidt(axest,jacobian);
        axest.Transpose(&axes);
        
        detjac = jacobian(0,0)*jacobian(1,1) - jacobian(1,0)*jacobian(0,1);
        if(IsZero(detjac)){
            detjac = ZeroTolerance();
        }
        jacinv(0,0) =  jacobian(1,1) / detjac;
        jacinv(1,1) =  jacobian(0,0) / detjac;
        jacinv(0,1) = -jacobian(0,1) / detjac;
        jacinv(1,0) = -jacobian(1,0) / detjac;
    }
    else
    {
        jacobian = JacTemp;
        axes.Resize(3,3); axes.Zero();
        axes(0,0) = 1.; axes(1,1) = 1.; axes(2,2) = 1.;
        detjac = -jacobian(0,2)*jacobian(1,1)*jacobian(2,0);//- a02 a11 a20
        detjac += jacobian(0,1)*jacobian(1,2)*jacobian(2,0);//+ a01 a12 a20
        detjac += jacobian(0,2)*jacobian(1,0)*jacobian(2,1);//+ a02 a10 a21
        detjac -= jacobian(0,0)*jacobian(1,2)*jacobian(2,1);//- a00 a12 a21
        detjac -= jacobian(0,1)*jacobian(1,0)*jacobian(2,2);//- a01 a10 a22
        detjac += jacobian(0,0)*jacobian(1,1)*jacobian(2,2);//+ a00 a11 a22
        
        if(IsZero(detjac)){
#ifdef LOG4CXX
            {
                std::stringstream sout;
                Print(sout);
                sout << "Parameter " << par << std::endl;
                coord.Print("Corner coordinates ", sout);
                /*					JacTemp.Print("Gradient of the coordinates",sout);
                 axes.Print("axes matrix", sout);
                 jacobian.Print("Jacobian", sout);
                 sout << "detjac " << detjac << std::endl;
                 int is;
                 for(is=TGeo::NNodes; is<TGeo::NSides-1; is++)
                 {
                 if(fNeighbours[is-TGeo::NNodes].Element())
                 {
                 sout << "Side: " << is << " El/side: " << fNeighbours[is-TGeo::NNodes].Element()->Index() << ":" <<
                 fNeighbours[is-TGeo::NNodes].Side() << '\n';
                 TPZGeoEl *gel = fNeighbours[is-TGeo::NNodes].Element();
                 gel->Print(sout);
                 }
                 }
                 */
                sout << "Singular jacobian " << detjac;
                LOGPZ_ERROR(logger,sout.str())
            }
#endif
            detjac = ZeroTolerance();
        }
        
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
void pzgeom::TPZGeoBlend<TGeo>::Print(std::ostream &out)
{
	TGeo::Print(out);
	out << "Neighbours/transformations used for mapping the sides :\n";
	int is;
	for(is=TGeo::NNodes; is<TGeo::NSides-1; is++)
	{
        out << "Side: " << is << " El/side: " << fNeighbours[is-TGeo::NNodes].ElementIndex() << ":" <<
        fNeighbours[is-TGeo::NNodes].Side() << '\n';
	}
}

template <class TGeo>
void pzgeom::TPZGeoBlend<TGeo>::Initialize(TPZGeoEl *refel)
{
    for(int byside = TGeo::NNodes; byside < (TGeo::NSides); byside++)
    {
        TPZGeoElSide ElemSide(refel,byside);
        TPZGeoElSide NextSide(refel,byside);
        if(!NextSide.Neighbour().Element()) continue;
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

/*
 //inseri este método pois o compilador não encontrou
 template <class TGeo>
 void TPZGeoBlend<TGeo>::Initialize(TPZVec<int> &nodeindexes)
 {
 TGeo::Initialize(nodeindexes);
 }
 */

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
TPZGeoEl *pzgeom::TPZGeoBlend<TGeo>::CreateBCGeoEl(TPZGeoEl *orig, int side,int bc)
{
	TPZStack<int> LowAllSides;
	TGeo::LowerDimensionSides(side,LowAllSides);
	if(side < 0 || side > TGeo::NSides-1)
	{
		DebugStop();
		return 0;
	}
	bool straight = true;
    TPZGeoMesh *gmesh = orig->Mesh();
	for(int lowside = 0; lowside < LowAllSides.NElements(); lowside++)
	{
		if(LowAllSides[lowside] >= TGeo::NNodes && Neighbour(LowAllSides[lowside],gmesh).Element())
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
TPZGeoEl *pzgeom::TPZGeoBlend<TGeo>::CreateBCGeoBlendEl(TPZGeoEl *orig,int side,int bc)
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


/**
 * Creates a geometric element according to the type of the father element
 */
template <class TGeo>
TPZGeoEl *pzgeom::TPZGeoBlend<TGeo>::CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
													  TPZVec<int>& nodeindexes,
													  int matid,
													  int& index)
{
	return CreateGeoElementMapped(mesh,type,nodeindexes,matid,index);
}





template class pzgeom::TPZGeoBlend<TPZGeoCube>;
template class pzgeom::TPZGeoBlend<TPZGeoTriangle>;
template class pzgeom::TPZGeoBlend<TPZGeoPrism>;
template class pzgeom::TPZGeoBlend<TPZGeoPyramid>;
template class pzgeom::TPZGeoBlend<TPZGeoTetrahedra>;
template class pzgeom::TPZGeoBlend<TPZGeoQuad>;
template class pzgeom::TPZGeoBlend<TPZGeoLinear>;
template class pzgeom::TPZGeoBlend<TPZGeoPoint>;

///CreateGeoElement -> TPZGeoBlend
#define IMPLEMENTBLEND(TGEO,CLASSID,CREATEFUNCTION) \
\
template<> \
int TPZGeoElRefPattern<TPZGeoBlend<TGEO>  >::ClassId() const { \
return CLASSID; \
} \
template class \
TPZRestoreClass< TPZGeoElRefPattern<TPZGeoBlend<TGEO> >, CLASSID>; \
\
template class TPZGeoElRefLess<TPZGeoBlend<TGEO> >;\
template class TPZGeoElRefPattern<TPZGeoBlend<TGEO> >;


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

