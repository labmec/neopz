/**
 * @file
 * @brief Contains the implementation of the TPZLine methods. 
 */

#include "tpzline.h"

#include "pzerror.h"
#include "pzreal.h"
#include "pztrnsform.h"
#include "pzquad.h"
#include "pzeltype.h"

#include "fad.h"


#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.topology.pzline");
#endif

using namespace std;

namespace pztopology {

	static constexpr int nhighdimsides[3] = {1,1,0};
	
	static constexpr int sidedimension[3] = {0,0,1};
	
	static constexpr int highsides[3][1] = {
		{2},
		{2},
		{0}
	};
	
	static constexpr REAL sidetosidetransforms[3][1][4][3] = {
		{
			{{0,0,0},{0,0,0},{0,0,0},{-1,0,0}}
		},
		{
			{{0,0,0},{0,0,0},{0,0,0},{1,0,0}}
		},
		{
			{{0,0,0},{0,0,0},{0,0,0},{0,0,0}}}
	};
	
	static constexpr REAL MidSideNode[3][1] = {{-1.},{1.},{0.}};
	
	static constexpr int nsidenodes[3] = {1,1,2};
    
    int TPZLine::NBilinearSides()
    {return 0;}
	
    static constexpr int vectorsideorder [3] = {0,1,2};
    
    static constexpr int bilinearounao [3] =   {0,0,1};
    static constexpr int direcaoksioueta [3] = {0,0,0};

    template<class T>
    inline void TPZLine::TShape(const TPZVec<T> &loc,TPZFMatrix<T> &phi,TPZFMatrix<T> &dphi) {
        T x = loc[0];
        phi(0,0) = (1.0-x)/2.;
        phi(1,0) = (1.0+x)/2.;
        dphi(0,0) = -0.5;
        dphi(0,1) = 0.5;
    }

    template<class T>
    void TPZLine::BlendFactorForSide(const int &side, const TPZVec<T> &xi, T &blendFactor,
                                        TPZVec<T> &blendFactorDxi){
        const REAL tol = pztopology::GetTolerance();
#ifdef PZDEBUG
        std::ostringstream sout;
        if(side < NCornerNodes || side >= NSides){
            sout<<"The side\t"<<side<<"is invalid. Aborting..."<<std::endl;
            PZError<<std::endl<<sout.str()<<std::endl;
            DebugStop();
        }

        if(!TPZLine::IsInParametricDomain(xi,tol)){
            sout<<"The method BlendFactorForSide expects the point xi to correspond to coordinates of a point";
            sout<<" inside the parametric domain. Aborting...";
            PZError<<std::endl<<sout.str()<<std::endl;
            #ifdef PZ_LOG
            LOGPZ_FATAL(logger,sout.str().c_str());
            #endif
            DebugStop();
        }
#endif
        TPZFNMatrix<4,T> phi(NCornerNodes,1);
        TPZFNMatrix<8,T> dphi(Dimension,NCornerNodes);
        TPZLine::TShape(xi,phi,dphi);
        blendFactorDxi.Resize(TPZLine::Dimension, (T) 0);
        switch(side){
            case 0:
                blendFactor = phi(0,0);
                blendFactorDxi[0] = dphi(0,0);
                break;
            case 1:
                blendFactor = phi(1,0);
                blendFactorDxi[0] = dphi(0,1);
                break;
            case 2:
                blendFactor = 1;
                blendFactorDxi[0] = 0;
                break;
        }
        return;
    }

	int TPZLine::NSideNodes(int side)
	{
		return nsidenodes[side];
	}
	
	int TPZLine::NumSides(int dimension) {
		if(dimension<0 || dimension> 1) {
			PZError << "TPZLine::NumSides. Bad parameter i.\n";
			return 0;
		}
		if(dimension==0) return 2;
		if(dimension==1) return 1;
		return -1;
	}
	
	
	void TPZLine::LowerDimensionSides(int side,TPZStack<int> &smallsides)
	{
		smallsides.Resize(0);
		int nsidecon = NContainedSides(side);
		int is;
		for(is=0; is<nsidecon-1; is++)
			smallsides.Push(ContainedSideLocId(side,is));
	}
	
	void TPZLine::LowerDimensionSides(int side,TPZStack<int> &smallsides, int DimTarget)
	{
		smallsides.Resize(0);
		int nsidecon = NContainedSides(side);
		for(int is = 0; is < nsidecon - 1; is++) {
			if (SideDimension(ContainedSideLocId(side,is)) == DimTarget) smallsides.Push(ContainedSideLocId(side,is));
		}
	}
	
	bool TPZLine::IsInParametricDomain(const TPZVec<REAL> &pt, REAL tol) {
		const REAL qsi = pt[0];
		if( fabs(qsi) <= 1. + tol){
			return true;
		}
		else{
			return false;
		}  
	}//method

    template<class T>
    bool TPZLine::CheckProjectionForSingularity(const int &side, const TPZVec<T> &xiInterior) {
        return true;
    }

    template<class T>
    void TPZLine::MapToSide(int side, TPZVec<T> &InternalPar, TPZVec<T> &SidePar, TPZFMatrix<T> &JacToSide) {
		TPZTransform<> TransfR = SideToSideTransform(NSides - 1, side);
        TPZTransform<T> Transf;
        Transf.CopyFrom(TransfR);
		SidePar.Resize(SideDimension(side));
		Transf.Apply(InternalPar,SidePar);
		
		JacToSide = Transf.Mult();
	}
    
    void TPZLine::ParametricDomainNodeCoord(int node, TPZVec<REAL> &nodeCoord)
    {
        if(node > NCornerNodes)
        {
            DebugStop();
        }
        nodeCoord.Resize(Dimension, 0.);
        switch (node) {
            case (0):
            {
                nodeCoord[0] = -1.;
                break;
            }
            case (1):
            {
                nodeCoord[0] = 1.;
                break;
            }
            default:
            {
                DebugStop();
                break;
            }
        }
    }
	
	void TPZLine::HigherDimensionSides(int side, TPZStack<int> &high)
	{
		if(side <0 || side >= NSides) {
			PZError << "TPZLine::HigherDimensionSides side "<< side << endl;
		}
		int is;
		for(is=0; is<nhighdimsides[side]; is++) high.Push(highsides[side][is]);
		
	}
	
	int TPZLine::SideNodeLocId(int side, int node)
	{
		if(side <2 && node == 0) return side;
		if(side == 2 && node <2) return node;
		PZError << "TPZLine::SideNodeLocId inconsistent side or node " << side
		<< ' ' << node << endl;
		return -1;
	}
	void TPZLine::CenterPoint(int side, TPZVec<REAL> &center) {
        if (center.size()!=Dimension) {
            DebugStop();
        }
		int i;
		for(i=0; i<Dimension; i++) {
			center[i] = MidSideNode[side][i];
		}
	}
	
    
    /** @brief Generates a random point in the master domain */
    void TPZLine::RandomPoint(TPZVec<REAL> &pt)
    {
        for(int i=0; i<1; i++)
        {
            REAL val = -1. + 2.*(REAL) rand() / (RAND_MAX);
            pt[i] = val;
        }
    }
    
	int TPZLine::SideDimension(int side) {
		if(side<0 || side >= NSides) {
			PZError << "TPZLine::SideDimension side " << side << endl;
			return -1;
		}
		return sidedimension[side];
	}
	
	TPZTransform<> TPZLine::SideToSideTransform(int sidefrom, int sideto)
	{
		if(sidefrom <0 || sidefrom >= NSides || sideto <0 || sideto >= NSides) {
			PZError << "TPZLine::HigherDimensionSides sidefrom "<< sidefrom << 
			' ' << sideto << endl;
			return TPZTransform<>(0);
		}
		if(sidefrom == sideto) {
			return TPZTransform<>(sidedimension[sidefrom]);
		}
		if(sidefrom == NSides-1) {
			return TransformElementToSide(sideto);
		}
        
        if (sideto == NSides -1) {
            return TransformSideToElement(sidefrom);
        }
        
		int nhigh = nhighdimsides[sidefrom];
		int is;
		for(is=0; is<nhigh; is++) {
			if(highsides[sidefrom][is] == sideto) {
				int dfr = sidedimension[sidefrom];
				int dto = sidedimension[sideto];
				TPZTransform<> trans(dto,dfr);
				int i,j;
				for(i=0; i<dto; i++) {
					for(j=0; j<dfr; j++) {
						trans.Mult()(i,j) = sidetosidetransforms[sidefrom][is][j][i];
					}
					trans.Sum()(i,0) = sidetosidetransforms[sidefrom][is][3][i];
				}
				return trans;
			}
		}
		PZError << "TPZLine::SideToSideTransform highside not found sidefrom "
		<< sidefrom << ' ' << sideto << endl;
		return TPZTransform<>(0);
	}
	
	TPZTransform<> TPZLine::TransformElementToSide(int side){
		
		if(side<0 || side>2){
			PZError << "TPZLine::TransformElementToSide called with side error\n";
			return TPZTransform<>(0,0);
		}
		
		//  int sidedim = SideDimension(side);
		TPZTransform<> t(SideDimension(side),1);//t(dimto,2)
		t.Mult().Zero();	//TPZGeoElQ2d *gq;
		t.Sum().Zero();//int dimto = gq->SideDimension(side);
		
		switch(side){
			case 0:
			case 1:
				return t;
			case 2:
				t.Mult()(0,0) = 1.0;//par. var.
				return t;
		}
		return TPZTransform<>(0,0);
	}
	
	TPZTransform<> TPZLine::TransformSideToElement(int side){
		
		if(side<0 || side>2){
			PZError << "TPZLine::TransformSideToElement side out range\n";
			return TPZTransform<>(0,0);
		}
		int sidedim = 1;
		if(side <2) sidedim = 0;
		
		TPZTransform<> t(1,sidedim);
		t.Mult().Zero();
		t.Sum().Zero();
		
		switch(side){
			case 0:
				t.Sum()(0,0) = -1.0;
				return t;
			case 1:
				t.Sum()(0,0) =  1.0;
				return t;
			case 2:
				t.Mult()(0,0) =  1.0;
				return t;
		}
		return TPZTransform<>(0,0);
	}
	
	
	TPZIntPoints *TPZLine::CreateSideIntegrationRule(int side, int order) {
		
		if(side<0 || side>2) {
			PZError << "TPZLine::CreateSideIntegrationRule wrong side " << side << endl;
			return 0;
		}
		if(side != 2) return new TPZInt1Point(order);   // sides 0 and 1 are vertices (corners)
		return new IntruleType(order);
		
		
	}
	MElementType TPZLine::Type(int side)
	{
		switch(side) {
			case 0:
			case 1:
				return EPoint;
			case 2:
				return EOned;
			default:
				return ENoType;
		}
	}
	
	
	int TPZLine::NumSides() {
		return 3;
	}
	
	
	int TPZLine::NContainedSides(int side) {
		if(side==0 || side==1) return 1;
		else if(side==2) return 3;
		PZError << "TPZLine::NContainedSides. Bad parameter side = " << side << " .\n";
        DebugStop();
		return 0;
	}
	
	int TPZLine::ContainedSideLocId(int side,int c) {
		switch(side) {
			case 0:
			case 1:
				if(c != 0)
				{
					PZError << "TPZLine::ContainedSideLocId, connect = " << c << endl;
				}
				return side;
			case 2:
				return c;
			default:
				PZError << "TPZLine::ContainedSideLocId called with side = " << side << endl;
				return 0;
		}
	}
	
	/**
	 * Method which identifies the transformation based on the IDs
	 * of the corner nodes
	 * @param id indexes of the corner nodes
	 * @return index of the transformation of the point corresponding to the topology
	 */
	int TPZLine::GetTransformId(const TPZVec<int64_t> &id)
	{
		return id[0] < id[1] ? 0 : 1;
	}
	
	/**
	 * Method which identifies the transformation of a side based on the IDs
	 * of the corner nodes
	 * @param id indexes of the corner nodes
	 * @return index of the transformation of the point corresponding to the topology
	 */	
	int TPZLine::GetTransformId(const int side, const TPZVec<int64_t> &id)
	{
		switch (side) {
			case 0:
			case 1:
				return 0;
				break;
			case 2:
				return id[0] < id[1] ? 0 : 1;
			default:
				break;
		}
		LOGPZ_ERROR(logger,"Wrong input parameter")
		return -1;
	}
    
    
	/**
	 * Identifies the permutation of the nodes needed to make neighbouring elements compatible 
	 * in terms of order of shape functions
	 * @param side : side for which the permutation is needed
	 * @param id : ids of the corner nodes of the elements
	 * @param permgather : permutation vector in a gather order
	 */
	void TPZLine::GetSideHDivPermutation(int transformationid, TPZVec<int> &permgather)
	{
        permgather.Resize(3);
        for (int i=0; i<3; i++)
        {
            permgather[i] = fPermutations[transformationid][i];
        }
		/*
         switch (side) {
			case 0:
			case 1:
				permgather[0] = 0;
				break;
			case 2:
				if(id[0]<id[1])
				{
					permgather[0] = 0;
					permgather[1] = 1;
					permgather[2] = 2;
				}
				else 
				{
					permgather[0] = 1;
					permgather[1] = 0;
					permgather[2] = 2;
				}
				
			default:
				break;
		}
		LOGPZ_ERROR(logger,"Wrong input parameter")*/
	}
    
    
    void TPZLine::ComputeDirections(int side, TPZFMatrix<REAL> &gradx, TPZFMatrix<REAL> &directions, TPZVec<int> &sidevectors)
    {
        
        if(gradx.Cols()!=1)
        {
            DebugStop();
        }
        TPZFMatrix<REAL> bvec = gradx;
        REAL norma = 0.0;
        for (int i = 0; i<3; i++)
        {
            norma += bvec(i,0)*bvec(i,0);
        }
        norma = sqrt(norma);

        directions.Redim(3, 1);
        for (int i=0; i<3; i++)
        {
            directions(i,0) = bvec(i,0)/norma;
        }
        sidevectors.Resize(1);
        sidevectors[0] = -1; // ??????????????????????

    }
    
    void TPZLine::GetSideHDivDirections(TPZVec<int> &sides, TPZVec<int> &dir, TPZVec<int> &bilounao)
    {
        int nsides = NumSides();
        
        sides.Resize(nsides);
        dir.Resize(nsides);
        bilounao.Resize(nsides);
        
        
        for (int is = 0; is<nsides; is++)
        {
            sides[is] = vectorsideorder[is];
            dir[is] = direcaoksioueta[is];
            bilounao[is] = bilinearounao[is];
        }
    }

    void TPZLine::GetSideHDivDirections(TPZVec<int> &sides, TPZVec<int> &dir, TPZVec<int> &bilounao, TPZVec<int> &sidevectors)
    {
        int nsides = NumSides();
        
        sides.Resize(nsides);
        dir.Resize(nsides);
        bilounao.Resize(nsides);
        
        for (int is = 0; is<nsides; is++)
        {
            sides[is] = vectorsideorder[is];
            dir[is] = direcaoksioueta[is];
            bilounao[is] = bilinearounao[is];
            sidevectors[is] = is;
        }
    }
    
    int TPZLine::ClassId() const{
        return Hash("TPZLine");
    }

    void TPZLine::Read(TPZStream& buf, void* context) {

    }
    
    void TPZLine::Write(TPZStream& buf, int withclassid) const {

    }
    
}

/**********************************************************************************************************************
 * The following are explicit instantiation of member function template of this class, both with class T=REAL and its
 * respective FAD<REAL> version. In other to avoid potential errors, always declare the instantiation in the same order
 * in BOTH cases.    @orlandini
 **********************************************************************************************************************/
template bool pztopology::TPZLine::CheckProjectionForSingularity<REAL>(const int &side, const TPZVec<REAL> &xiInterior);

template void pztopology::TPZLine::MapToSide<REAL>(int side, TPZVec<REAL> &InternalPar, TPZVec<REAL> &SidePar, TPZFMatrix<REAL> &JacToSide);

template void pztopology::TPZLine::BlendFactorForSide<REAL>(const int &, const TPZVec<REAL> &, REAL &, TPZVec<REAL> &);

template void pztopology::TPZLine::TShape<REAL>(const TPZVec<REAL> &loc,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);

template void pztopology::TPZLine::ComputeHDivDirections<REAL>(TPZFMatrix<REAL> &gradx, TPZFMatrix<REAL> &directions);

template bool pztopology::TPZLine::CheckProjectionForSingularity<Fad<REAL>>(const int &side, const TPZVec<Fad<REAL>> &xiInterior);

template void pztopology::TPZLine::MapToSide<Fad<REAL> >(int side, TPZVec<Fad<REAL> > &InternalPar, TPZVec<Fad<REAL> > &SidePar, TPZFMatrix<Fad<REAL> > &JacToSide);

template void pztopology::TPZLine::BlendFactorForSide<Fad<REAL>>(const int &, const TPZVec<Fad<REAL>> &, Fad<REAL> &,
                                                                   TPZVec<Fad<REAL>> &);
template void pztopology::TPZLine::TShape<Fad<REAL>>(const TPZVec<Fad<REAL>> &loc,TPZFMatrix<Fad<REAL>> &phi,TPZFMatrix<Fad<REAL>> &dphi);

template void pztopology::TPZLine::ComputeHDivDirections<Fad<REAL>>(TPZFMatrix<Fad<REAL>> &gradx, TPZFMatrix<Fad<REAL>> &directions);
