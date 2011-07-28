#ifndef TPZARC3D_H
#define TPZARC3D_H

// #include "pzfmatrix.h"
// #include "pzvec.h"
// #include "pzgmesh.h"


#include "pzgeoel.h"
#include "pznoderep.h"
// #include "pzgnode.h"
#include "tpzline.h"
#include "tpzgeoelrefpattern.h"


#include <iostream>

/**
 / Class made by Paulo Cesar de Alvarenga Lucci (Caju)
 / LabMeC - FEC - UNICAMP
 / 2007
 */
namespace pzgeom
{
	
	/** 
	 @ingroup geometry
	 @brief Implements three dimensional arc (Jorge?)
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
		
		TPZArc3D(TPZFMatrix &coord){
			ComputeAtributes(coord);
		}
		
		/** @brief Initialize the internal data structure of the arc using the coordinates of the nodes */
		void Initialize(TPZGeoEl *refel)
		{
			int nnod = 3;
			TPZFMatrix coord(3,nnod);
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
		
		
		void X(TPZFMatrix &coord,TPZVec<REAL> &loc,TPZVec<REAL> &result);
		void Jacobian(TPZFMatrix &coord, TPZVec<REAL> &par, TPZFMatrix &jacobian, TPZFMatrix &axes, REAL &detjac, TPZFMatrix &jacinv);
		
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
		
	protected:
		
		void ComputeAtributes(TPZFMatrix &coord);
		void ComputeR2Points(TPZFMatrix &coord, double &xa, double &ya, double &xb, double &yb);
		double ArcAngle(TPZFMatrix &coord, double xa, double ya, double xb, double yb) const;
		
		/** @name Atributes */
		// @{
		TPZFNMatrix<9> fICnBase, fIBaseCn;
		TPZManVector< REAL,3 > fCenter3D, finitialVector;
		double fAngle, fRadius, fXcenter, fYcenter;
		// @}
	};
	
};

#define TPZGEOELEMENTARC3DID 350
template<>
inline int TPZGeoElRefPattern<pzgeom::TPZArc3D>::ClassId() const {
	return TPZGEOELEMENTARC3DID;
}


//         template class pzgeom::TPZNodeRep<3,TPZArc3D>;
//template class TPZGeoElRefLess<TPZArc3D>;


#endif
