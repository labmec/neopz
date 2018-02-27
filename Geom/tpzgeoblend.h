/**
 * @file
 * @brief Contains the TPZGeoBlend class which implements a blending map from curved boundaries to the interior of the element.
 */

#ifndef TPZGEOBLEND_H
#define TPZGEOBLEND_H

#include "pzgeotriangle.h"
#include "TPZGeoCube.h"
#include "pzgeoprism.h"
#include "pzgeoel.h"
#include "pznoderep.h"
#include "pzgeoelside.h"

#include <iostream>

namespace pzgeom 
{
	
	/**
	 * @ingroup geometry
	 * @brief Implements a blending map from curved boundaries to the interior of the element. \ref geometry "Geometry"
	 * @author Paulo Cesar de Alvarenga Lucci
	 * @since 2007
	 */
	template <class TGeo>
	class TPZGeoBlend : public TGeo {
		
	public:
            
            virtual int ClassId() const;
		
        bool IsLinearMapping(int side) const;
        
		bool IsGeoBlendEl() const { 
			return true; 
		}
		
		/** @brief Constructor with list of nodes */
		TPZGeoBlend(TPZVec<int64_t> &nodeindexes) : TPZRegisterClassId(&TPZGeoBlend::ClassId),
        TGeo(nodeindexes) {
		}
		/** @brief Empty constructor */
		TPZGeoBlend() : TPZRegisterClassId(&TPZGeoBlend::ClassId),
        TGeo() {
		}
		
		/** @brief Constructor with node map */
		TPZGeoBlend(const TPZGeoBlend &cp,std::map<int64_t,int64_t> & gl2lcNdMap) : TPZRegisterClassId(&TPZGeoBlend::ClassId),
        TGeo(cp,gl2lcNdMap) {
		}
		/** @brief Copy constructor */
		TPZGeoBlend(const TPZGeoBlend &cp) : TPZRegisterClassId(&TPZGeoBlend::ClassId),
        TGeo(cp) {
		}
		/** @brief Copy constructor */
		TPZGeoBlend(const TPZGeoBlend &cp, TPZGeoMesh &) : TPZRegisterClassId(&TPZGeoBlend::ClassId),
        TGeo(cp) {
		}
        
        void Read(TPZStream &buf, void *context) {
            TGeo::Read(buf, context);
            for (int is = 0; is < 1 + TGeo::NSides - TGeo::NNodes; is++) {
                fNeighbours[is].Read(buf, context);
            }
            for (int is = 0; is < 1 + TGeo::NSides - TGeo::NNodes; is++) {
                fTrans[is].Read(buf, context);
            }
        }

        void Write(TPZStream &buf, int withclassid) const {
            TGeo::Write(buf, withclassid);
            for (int is = 0; is < 1 + TGeo::NSides - TGeo::NNodes; is++) {
                fNeighbours[is].Write(buf, withclassid);
            }
            for (int is = 0; is < 1 + TGeo::NSides - TGeo::NNodes; is++) {
                fTrans[is].Write(buf, withclassid);
            }
        }

		void SetNeighbourInfo(int side, TPZGeoElSide &neigh, TPZTransform<> &trans);
		
		TPZGeoElSide Neighbour(int side,TPZGeoMesh *gmesh) const {
            if (side < TGeo::NNodes) {
                DebugStop();
            }
			return TPZGeoElSide(fNeighbours[side-TGeo::NNodes],gmesh);
		}
		
        template<class T>
		void TransfBetweenNeigh(int side, TPZTransform<T> &tr) const {
			tr.CopyFrom(fTrans[side - TGeo::NNodes]);
		}
		
		/** @brief Returns the type name of the element */
		static std::string TypeName() { return TGeo::TypeName();} 
		
        /** @brief Get the coordinates of the point at geometric elements from coordinates of the parametric point at the master element */
        template<class T>
		void X(const TPZGeoEl &gel, TPZVec<T>& par, TPZVec<T> &result) const;
        
		/** @brief Computes the Jacobian for parametric point at master element */
		void Jacobian(const TPZGeoEl &gel, TPZVec<REAL>& par, TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) const;
		
		/** @brief Computes the gradient of the transformation for parametric point at master element */
        template<class T>
        void GradX(const TPZGeoEl &gel, TPZVec<T> &par, TPZFMatrix<T> &gradx) const;
        
		/** @brief Print all relevant data of the element to cout*/
		void Print(std::ostream & out = std::cout) const;
		
		/**
		 * @brief Initialize the element checking connectivity on all sides.
		 *
		 / Este método percorre todos os lados do elemento fornecido (el) e, para cada lado,
		 / percorre todos os vizinhos até retornar a si próprio. Se neste processo encontrar algum elemento que
		 / seja quadrático e não seja do tipo TPZGeoBlend, recorre ao método SetNeighbourInfo() para estabelecer que este
		 / elemento encontrado será seu vizinho pelo respectivo lado.
		 */
		void Initialize(TPZGeoEl *refel);
		
		//void Initialize(TPZVec<int> &nodeindexes, TPZGeoMesh &mesh);
		
		/**
		 * @brief Method which creates a geometric boundary condition 
		 * element based on the current geometric element, 
		 * a side and a boundary condition number
		 */
		TPZGeoEl *CreateBCGeoEl(TPZGeoEl *orig, int side,int bc);
		
		TPZGeoEl *CreateBCGeoBlendEl(TPZGeoEl *orig,int side,int bc);
		
//		TPZGeoEl *CreateGeoBlend(TPZGeoMesh &mesh, MElementType type, TPZVec<int64_t>& nodeindexes, int matid, int64_t& index);
		
        
		
	public:
        
        /// create an example element based on the topology
        /* @param gmesh mesh in which the element should be inserted
         @param matid material id of the element
         @param lowercorner (in/out) on input lower corner o the cube where the element should be created, on exit position of the next cube
         @param size (in) size of space where the element should be created
         */
        static void InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size);

		/**
		 * @brief Creates a geometric element according to the type of the father element
		 */
		static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
										  TPZVec<int64_t>& nodeindexes,
										  int matid,
										  int64_t& index);
		
	protected:
		
		/// Project the InternalPar parameter to the parameter of the neighbour along side. Return true if the map is nonsingular
        template<class T>
		bool MapToNeighSide(int side, int sidedim, TPZVec<T> &InternalPar, TPZVec<T> &NeighPar, TPZFMatrix<T> &JacNeighSide) const;
		/** @brief Vector of indexes of the neighbours */
		//TPZGeoElSideIndex fNeighbours[1+TGeo::NSides - TGeo::NNodes];
		TPZGeoElSideIndex fNeighbours[1+TGeo::NSides - TGeo::NNodes];
		TPZTransform<> fTrans[1+TGeo::NSides - TGeo::NNodes];
	};
    
    
    template <class TGeo>
    template<class T>
    inline void pzgeom::TPZGeoBlend<TGeo>::X(const TPZGeoEl &gel, TPZVec<T>& par, TPZVec<T> &result) const
    {
        TPZFNMatrix<45> coord(3,TGeo::NNodes);
        this->CornerCoordinates(gel,coord);
        
        result.Resize(3);
        result.Fill(0);
        
        TPZManVector<T,3> NeighPar, SidePar, Xside(3,0.);
        
        int majorSide = TGeo::NSides - 1;
        
        TPZManVector<REAL,27> SidesCounter(TGeo::NSides,0);
        TPZStack<int> LowNodeSides, LowAllSides;
        
        TPZFNMatrix<9,T> blend(TGeo::NNodes,1), Dblend(TGeo::Dimension,TGeo::NNodes), NotUsedHere;
        TGeo::TShape(par,blend,Dblend);
        TPZGeoMesh *gmesh = gel.Mesh();
        
        for(int byside = majorSide; byside >= TGeo::NNodes; byside--)
        {
            TPZGeoElSide gelside(fNeighbours[byside-TGeo::NNodes],gmesh);
            if(gelside.Exists())
            {
                TGeo::LowerDimensionSides(byside,LowNodeSides,0);
                TGeo::LowerDimensionSides(byside,LowAllSides);
                T blendTemp = 0.;
                for(int a = 0; a < LowNodeSides.NElements(); a++)
                {
                    blendTemp += blend(LowNodeSides[a],0);
                }
                int sidedim = gelside.Dimension();
                
                if(!MapToNeighSide(byside,sidedim,par,NeighPar,NotUsedHere))
                {
#ifdef LOG4CXX2
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
                
#ifdef LOG4CXX2
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
        
#ifdef LOG4CXX2
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
        
#ifdef LOG4CXX2
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            sout << "result " << result;
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif
        
    }
    
    template <class TGeo>
    template<class T>
    inline bool pzgeom::TPZGeoBlend<TGeo>::MapToNeighSide(int side, int SideDim, TPZVec<T> &InternalPar, TPZVec<T> &NeighPar, TPZFMatrix<T> &JacNeighSide) const
    {
        TPZFNMatrix<9,T> JacSide;
        
        TPZManVector< T, 3 > SidePar(SideDim);
        const bool check = this->MapToSide(side, InternalPar, SidePar,JacSide);
        
#ifdef LOG4CXX2
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
        TPZTransform<T> tr;
        TransfBetweenNeigh(side, tr);
        tr.Apply(SidePar,NeighPar);
        
        JacNeighSide.Resize(0, 0);
        TransfBetweenNeigh(side, tr);
        JacNeighSide.Resize(tr.Mult().Rows(),JacSide.Cols());
        tr.Mult().Multiply(JacSide,JacNeighSide);
        
#ifdef LOG4CXX2
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
    
//    /** @brief Computes the gradient of the transformation for parametric point at master element */
//    template <class TGeo>
//    template<class T>
//    inline void pzgeom::TPZGeoBlend<TGeo>::GradX(const TPZGeoEl &gel, TPZVec<T> &par, TPZFMatrix<T> &gradx) const
//    {
//        
//        TPZVec<REAL> parloc;
//        
//        TPZFNMatrix<45> coord(3,TGeo::NNodes);
//        this->CornerCoordinates(gel,coord);
//        
//        TPZManVector<REAL,3> NeighPar, SidePar, Xside(3,0.), XNode(3,0.);
//        int majorSide = TGeo::NSides - 1;
//        
//        TPZManVector<REAL> SidesCounter(TGeo::NSides,0);
//        TPZStack<int> LowNodeSides, LowAllSides;
//        
//#ifdef LOG4CXX
//        if(logger->isDebugEnabled())
//        {
//            std::stringstream sout;
//            sout << "input parameter par " << par;
//            LOGPZ_DEBUG(logger,sout.str())
//        }
//#endif
//        
//        TPZFNMatrix<24> blend(TGeo::NNodes,1), Dblend(TGeo::Dimension,TGeo::NNodes), NotUsedHere;
//        TGeo::Shape(parloc,blend,Dblend);
//        
//        TPZFNMatrix<9> J1, J2, Ax, JacTemp(3,TGeo::Dimension, 0.), Jneighbourhood;
//        REAL Det;
//        TPZGeoMesh *gmesh = gel.Mesh();
//        for(int byside = majorSide; byside >= TGeo::NNodes; byside--)
//        {
//            TPZGeoElSide neighbyside = Neighbour(byside, gmesh);
//            if(neighbyside.Exists())
//            {
//                TGeo::LowerDimensionSides(byside,LowNodeSides,0);
//                TGeo::LowerDimensionSides(byside,LowAllSides);
//                int dim = Neighbour(byside, gmesh).Dimension();
//                TPZFNMatrix<9> Inv(dim,dim);
//                int sidedim = neighbyside.Dimension();
//                if(!MapToNeighSide(byside,sidedim,par,NeighPar, Jneighbourhood))
//                {
//                    continue;
//                }
//                Neighbour(byside,gmesh).X(NeighPar,Xside);
//                Neighbour(byside,gmesh).Jacobian(NeighPar,J1,Ax,Det,Inv);
//                Ax.Transpose();
//                Ax.Multiply(J1,J2);
//                
//#ifdef LOG4CXX
//                if(logger->isDebugEnabled())
//                {
//                    std::stringstream sout;
//                    sout << "byside " << byside << std::endl;
//                    sout << "side of the neighbour " << Neighbour(byside,gmesh).Side() << std::endl;
//                    Neighbour(byside,gmesh).Element()->Print(sout);
//                    sout << "neighbour parameter(NeighPar) " << NeighPar << std::endl;
//                    sout << "Jacobian of the map(Jneighborhood) " << Jneighbourhood << std::endl;
//                    sout << "Xside " << Xside << std::endl;
//                    sout << "jacobian neighbour(J1) " << J1 << std::endl;
//                    Ax.Print("Ax of the neighbour (Ax) ",sout);
//                    sout << "jacobian of the neighbour multiplied by the axes(J2) " << J2 << std::endl;
//                    LOGPZ_DEBUG(logger,sout.str())
//                }
//#endif
//                
//                J2.Multiply(Jneighbourhood,J1);
//                
//#ifdef LOG4CXX
//                if(logger->isDebugEnabled())
//                {
//                    std::stringstream sout;
//                    sout << "acumulated jacobian(J1) " << J1 << std::endl;
//                    LOGPZ_DEBUG(logger,sout.str())
//                }
//#endif
//                
//                REAL blendTemp = 0.;
//                TPZManVector<REAL,3> DblendTemp(TGeo::Dimension,0.);
//                for(int a = 0; a < LowNodeSides.NElements(); a++)
//                {
//                    TPZManVector<REAL> parChanged(par.NElements());
//                    TGeo::FixSingularity(byside,parloc,parChanged);
//                    TGeo::Shape(parChanged,blend,Dblend);
//                    
//                    blendTemp += blend(LowNodeSides[a],0);
//                    for(int b = 0; b < TGeo::Dimension; b++)
//                    {
//                        DblendTemp[b] += Dblend(b,LowNodeSides[a]);
//                    }
//                }
//                
//                for(int a = 0; a < 3; a++)
//                {
//                    for(int b = 0; b < TGeo::Dimension; b++)
//                    {
//                        JacTemp(a,b) += (1 - SidesCounter[byside]) * (J1(a,b)*blendTemp + Xside[a]*DblendTemp[b]);
//                    }
//                }
//                for(int a = 0; a < LowAllSides.NElements(); a++)
//                {
//                    SidesCounter[LowAllSides[a]] += (1 - SidesCounter[byside]);
//                }
//            }
//        }
//        
//#ifdef LOG4CXX
//        if(logger->isDebugEnabled())
//        {
//            std::stringstream sout;
//            JacTemp.Print("Jabobian before contributing the nodes",sout);
//            sout << "SidesCounter " << SidesCounter << std::endl;
//            sout << "DBlend " << Dblend << std::endl;
//            sout << "NodeCoord " << coord << std::endl;
//            LOGPZ_DEBUG(logger,sout.str())
//        }
//#endif
//        
//        for(int a = 0; a < TGeo::NNodes; a++)
//        {
//            for(int b = 0; b < 3; b++) 
//            { 
//                for(int c = 0; c < TGeo::Dimension; c++)
//                {
//                    JacTemp(b,c) += (1 - SidesCounter[a]) * coord(b,a)*Dblend(c,a);
//                }
//            }
//        }
//
//    }
    
    template <class TGeo>
    int TPZGeoBlend<TGeo>::ClassId() const{
        return Hash("TPZGeoBlend") ^ TGeo::ClassId() << 1;
    }
	
};
#endif
