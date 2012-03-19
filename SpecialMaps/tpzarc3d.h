/**
 * @file
 * @brief Contains the TPZArc3D class which implements three dimensional arc.
 */
#ifndef TPZARC3D_H
#define TPZARC3D_H

#include "pzgeoel.h"
#include "pznoderep.h"
// #include "pzgnode.h"
#include "tpzline.h"
#include "tpzgeoelrefpattern.h"

#include <iostream>

namespace pzgeom
{
	
	/** 
	 * @author Paulo Cesar de Alvarenga Lucci (Caju)
	 * @ingroup geometry
	 * @brief Implements three dimensional arc. \ref geometry "Geometry"
	 * @since 2007
	 */
	class TPZArc3D : public pzgeom::TPZNodeRep<3,pztopology::TPZLine> {
		
	public:
		
		enum {NNodes = 3};
		
		bool IsLinearMapping() const { return false; }
		
		TPZArc3D(const TPZArc3D &cp,std::map<int,int> & gl2lcNdMap) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(cp,gl2lcNdMap){
			this->fICnBase = cp.fICnBase;
			this->fIBaseCn = cp.fIBaseCn;
			this->fCenter3D = cp.fCenter3D;
			this->fRadius = cp.fRadius;		
		}
		
		TPZArc3D() : pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(),fICnBase(3,3),fIBaseCn(3,3) {
		}
		
		TPZArc3D(const TPZArc3D &cp) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(cp){
			this->fICnBase = cp.fICnBase;
			this->fIBaseCn = cp.fIBaseCn;
			this->fCenter3D = cp.fCenter3D;
			this->fRadius = cp.fRadius;
		}
		
		TPZArc3D(const TPZArc3D &cp, TPZGeoMesh &) : pzgeom::TPZNodeRep<NNodes, pztopology::TPZLine>(cp){
			this->fICnBase  = cp.fICnBase;
			this->fIBaseCn  = cp.fIBaseCn;
			this->fCenter3D = cp.fCenter3D;
			this->fRadius   = cp.fRadius;
		}
		
		TPZArc3D(TPZVec<int> &nodeindexes) : pzgeom::TPZNodeRep<NNodes,pztopology::TPZLine>(nodeindexes), fICnBase(3,3), fIBaseCn(3,3) {
			int nnod = nodeindexes.NElements();
			if(nnod != 3)
			{
				std::cout << "Arc geometry created with " << nnod << " nodes, bailing out\n";
				DebugStop();
			}
		}
		
		TPZArc3D(TPZFMatrix<REAL> &coord){
			ComputeAtributes(coord);
		}
		
		/** @brief Initialize the internal data structure of the arc using the coordinates of the nodes */
		void Initialize(TPZGeoEl *refel)
		{
			int nnod = 3;
			TPZFMatrix<REAL> coord(3,nnod);
			int nod, co;
			for(nod=0; nod<3; nod++)
			{
				for(co=0; co<3; co++)
				{
					coord(co,nod) = refel->NodePtr(nod)->Coord(co);
				}
			}
			ComputeAtributes(coord);
			
		}
		
		
		void X(TPZFMatrix<REAL> &coord,TPZVec<REAL> &loc,TPZVec<REAL> &result) const;
		void Jacobian(TPZFMatrix<REAL> &coord, TPZVec<REAL> &par, TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes, REAL &detjac, TPZFMatrix<REAL> &jacinv) const;
		
        void X(const TPZGeoEl &gel,TPZVec<REAL> &loc,TPZVec<REAL> &result) const
        {
            TPZFNMatrix<3*NNodes> coord(3,NNodes);
            CornerCoordinates(gel, coord);
            X(coord,loc,result);
        }
		
		void Jacobian(const TPZGeoEl &gel,TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) const
        {
            TPZFNMatrix<3*NNodes> coord(3,NNodes);
            CornerCoordinates(gel, coord);
            Jacobian(coord, param, jacobian, axes, detjac, jacinv);
        }

		
        static std::string TypeName() { return "Linear";}
		static TPZGeoEl * CreateBCGeoEl(TPZGeoEl *orig, int side,int bc);
		
	public:
		/**
		 * @brief Creates a geometric element according to the type of the father element
		 */
		static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
										  TPZVec<int>& nodeindexes,
										  int matid,
										  int& index);
		void Print(std::ostream &out)
		{
			pzgeom::TPZNodeRep<3,pztopology::TPZLine>::Print(out);
			out << "fCenter3D " << fCenter3D << " finitialVector " << finitialVector << std::endl;
			out << "fAngle " << fAngle << " fRadius " << fRadius << " fXcenter " << fXcenter << " fYcenter " << fYcenter << std::endl;
			fICnBase.Print("fICnBase", out);
			fIBaseCn.Print("fIBaseCn", out);
		}
        
        void Read(TPZStream &buf,void *context)
        {
            pzgeom::TPZNodeRep<3,pztopology::TPZLine>::Read(buf,0);
            fICnBase.Read(buf,0);
            fIBaseCn.Read(buf,0);
            TPZSaveable::ReadObjects<3>(buf, fCenter3D);
            TPZSaveable::ReadObjects<3>(buf,finitialVector);
            buf.Read(&fAngle,1);
            buf.Read(&fRadius,1);
            buf.Read(&fXcenter,1);
            buf.Read(&fYcenter,1);
        }
        
        void Write(TPZStream &buf)
        {
            pzgeom::TPZNodeRep<3,pztopology::TPZLine>::Write(buf);
            fICnBase.Write(buf,0);
            fIBaseCn.Write(buf,0);
            TPZSaveable::WriteObjects(buf, fCenter3D);
            TPZSaveable::WriteObjects(buf,finitialVector);
            buf.Write(&fAngle,1);
            buf.Write(&fRadius,1);
            buf.Write(&fXcenter,1);
            buf.Write(&fYcenter,1);
		}
        

		
	protected:
		
		void ComputeAtributes(TPZFMatrix<REAL> &coord);
		void ComputeR2Points(TPZFMatrix<REAL> &coord, double &xa, double &ya, double &xb, double &yb);
		double ArcAngle(TPZFMatrix<REAL> &coord, double xa, double ya, double xb, double yb) const;
		
		/** @name Atributes */
		// @{
		TPZFNMatrix<9> fICnBase, fIBaseCn;
		TPZManVector< REAL,3 > fCenter3D, finitialVector;
#ifdef contar
		double fAngle, fRadius, fXcenter, fYcenter;
#else
		REAL fAngle, fRadius, fXcenter, fYcenter;
#endif
		// @}
	};
	
};

/// Id for three dimensional arc element
#define TPZGEOELEMENTARC3DID 350
template<>
inline int TPZGeoElRefPattern<pzgeom::TPZArc3D>::ClassId() const {
	return TPZGEOELEMENTARC3DID;
}

#endif
