/**
 * @file
 * @brief Contains the implementation of the TPZQuadrilateral methods. 
 */

#include "tpzquadrilateral.h"

#include "pzmanvector.h"
#include "pzerror.h"
#include "pzreal.h"
#include "pzquad.h"
#include "pzeltype.h"

#include "pzcreateapproxspace.h"
#include "pzshapequad.h"

#include "pznumeric.h"

#ifdef _AUTODIFF
#include "fad.h"
#endif

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.topology.pzquadrilateral"));
#endif

using namespace std;

namespace pztopology {

	static int nsidenodes[9] = {
		1,1,1,1,2,2,2,2,4};
	
	int TPZQuadrilateral::NSideNodes(int side)
	{
		return nsidenodes[side];
	}
	
	static int nhighdimsides[9] = {3,3,3,3,1,1,1,1,0};
	
	
	static int sidedimension[9] = {0,0,0,0,1,1,1,1,2};
	
	static REAL sidetosidetransforms[9][3][4][3] = {
		{
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-1,-99}}
		},
		{
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-1,-99}}
		},
		{
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,1,-99}}
		},
		{
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,1,-99}}
		},
		{
			{{1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,-1,-99}}
		},
		{
			{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}}
		},
		{
			{{-1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}}
		},
		{
			{{0,-1,-99},{-99,-99,-99},{-99,-99,-99},{-1,0,-99}}
		},
		{
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-99,-99,-99}}
		}
	};
	
	
	static int highsides[9][3] = {
		{4,7,8},
		{4,5,8},
		{5,6,8},
		{6,7,8},
		{8},
		{8},
		{8},
		{8},
		{-999}
	};
	
	static REAL MidSideNode[9][3] = {
		/*00*/{-1.,-1.},/*01*/{ 1.,-1.},/*02*/{1.,1.},
		/*03*/{-1., 1.},/*04*/{ 0.,-1.},/*05*/{1.,0.},
		/*06*/{ 0., 1.},/*07*/{-1., 0.},/*08*/{0.,0.} };
    
    static REAL bQuad[18][2] = 
    {
        {0,-1},
        {0,-1},
        {0,-1},
        {1,0},
        {1,0},
        {1,0},
        {0,1},
        {0,1},
        {0,1},
        {-1,0},
        {-1,0},
        {-1,0},
        {1,0},
        {0,1},
        {-1,0},
        {0,-1},
        {1,0},
        {0,1}
        
    };
    
    static REAL tQuad[18][2] = 
    {
        {-1,0},
        {-1,0},
        {-1,0},
        {0,-1},
        {0,-1},
        {0,-1},
        {1,0},
        {1,0},
        {1,0},
        {0,1},
        {0,1},
        {0,1},
        {0,-1},
        {1,0},
        {0,1},
        {-1,0},
        {0,-1},
        {1,0}
    };
    
    static int vectorsideorder [18] = {0,1,4,1,2,5,2,3,6,3,0,7,4,5,6,7,8,8};
    
    static int bilinearounao [18] =   {0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1};
    static int direcaoksioueta [18] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1};
    
    static int permutationsQ [8][9] =
    {
        {0,1,2,3,4,5,6,7,8}, // id 0
        {0,3,2,1,7,6,5,4,8}, // id 1
        {1,2,3,0,5,6,7,4,8}, // id 2
        {1,0,3,2,4,7,6,5,8}, // id 3
        {2,3,0,1,6,7,4,5,8}, // id 4
        {2,1,0,3,5,4,7,6,8}, // id 5
        {3,0,1,2,7,4,5,6,8}, // id 6
        {3,2,1,0,6,5,4,7,8}  // id 7
    };

    int TPZQuadrilateral::NBilinearSides()
    {
        return 6;
    }
	
	int TPZQuadrilateral::SideNodeLocId(int side, int node)
	{
		if(side<4 && node==0) return side;
		if(side>=4 && side<8 && node <2) return (side+node)%4;
		if(side==8 && node <4) return node;
		PZError << "TPZQuadrilateral::SideNodeLocId inconsistent side or node " << side
		<< ' ' << node << endl;
		return -1;
	}
	
	void TPZQuadrilateral::LowerDimensionSides(int side,TPZStack<int> &smallsides)
	{
		smallsides.Resize(0);
		int nsidecon = NContainedSides(side);
		int is;
		for(is=0; is<nsidecon-1; is++)
			smallsides.Push(ContainedSideLocId(side,is));
	}
	
	void TPZQuadrilateral::LowerDimensionSides(int side,TPZStack<int> &smallsides, int DimTarget)
	{
		smallsides.Resize(0);
		int nsidecon = NContainedSides(side);
		for(int is = 0; is < nsidecon - 1; is++) {
			if (SideDimension(ContainedSideLocId(side,is)) == DimTarget) smallsides.Push(ContainedSideLocId(side,is));
		}
	}
	
	void TPZQuadrilateral::HigherDimensionSides(int side, TPZStack<int> &high)
	{
		if(side <0 || side >= NSides) {
			PZError << "TPZQuadrilateral::HigherDimensionSides side "<< side << endl;
		}
		int is;
		for(is=0; is<nhighdimsides[side]; is++) high.Push(highsides[side][is]);
		
	}
	
	void TPZQuadrilateral::CenterPoint(int side, TPZVec<REAL> &center) {
		//center.Resize(Dimension);
		int i;
		for(i=0; i<Dimension; i++) {
			center[i] = MidSideNode[side][i];
		}
	}
	
	TPZIntPoints * TPZQuadrilateral::CreateSideIntegrationRule(int side, int order){
		if(side<0 || side>8) {
			PZError << "TPZQuadrilateral::CreateSideIntegrationRule wrong side " << side << endl;
			return 0;
		}
		if(side<4) return new TPZInt1Point(order);              // sides 0 to 3 are vertices (corners)
		if(side<8) return new TPZInt1d(order);             // sides 4 to 7 are lines
		if(side==8) return new IntruleType(order,order);   // integration of the element
		return 0;
	}
	
	
	TPZTransform<> TPZQuadrilateral::TransformElementToSide(int side){
		
		if(side<0 || side>8){
			PZError << "TPZShapeQuad::TransformElementToSide called with side error\n";
			return TPZTransform<>(0,0);
		}
		
		TPZTransform<> t(sidedimension[side],2);//t(dimto,2)
		t.Mult().Zero();
		t.Sum().Zero();
		
		switch(side){
			case 0:
			case 1:
			case 2:
			case 3:
				return t;
			case 4:
				t.Mult()(0,0) = 1.0;
				return t;
			case 5 :
				t.Mult()(0,1) = 1.0;
				return t;
			case 6:
				t.Mult()(0,0) = -1.0;
				return t;
			case 7:
				t.Mult()(0,1) = -1.0;
				return t;
			case 8:
				t.Mult()(0,0) = 1.0;
				t.Mult()(1,1) = 1.0;
				return t;
		}
		return TPZTransform<>(0,0);
	}
	
	bool TPZQuadrilateral::IsInParametricDomain(TPZVec<REAL> &pt, REAL tol){
		const REAL qsi = pt[0];
		const REAL eta = pt[1];
		if( ( fabs(qsi) <= 1. + tol ) && ( fabs(eta) <= 1. + tol ) ){
			return true;
		}
		else{
			return false;
		}  
	}//method

    /** @brief Generates a random point in the master domain */
    void TPZQuadrilateral::RandomPoint(TPZVec<REAL> &pt)
    {
        for(int i=0; i<2; i++)
        {
            REAL val = -1. + 2.*(REAL) rand() / (RAND_MAX);
            pt[i] = val;
        }
    }
    

    template<class T>
    bool TPZQuadrilateral::MapToSide(int side, TPZVec<T> &InternalPar, TPZVec<T> &SidePar, TPZFMatrix<T> &JacToSide) {
		bool regularmap = true;
        TPZTransform<T> Transf;
        Transf.CopyFrom(pztopology::TPZQuadrilateral::SideToSideTransform(NSides - 1, side));
		SidePar.Resize(SideDimension(side));
		Transf.Apply(InternalPar,SidePar);
		
		int R = Transf.Mult().Rows();
		int C = Transf.Mult().Cols();
		
		JacToSide.Resize(R,C);
		for(int i = 0; i < R; i++)
		{
			for(int j = 0; j < C; j++) JacToSide(i,j) = Transf.Mult()(i,j);
		}
		return regularmap;
	}
    
    void TPZQuadrilateral::ParametricDomainNodeCoord(int node, TPZVec<REAL> &nodeCoord)
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
                nodeCoord[1] = -1.;
                break;
            }
            case (1):
            {
                nodeCoord[0] =  1.;
                nodeCoord[1] = -1.;
                break;
            }
            case (2):
            {
                nodeCoord[0] = 1.;
                nodeCoord[1] = 1.;
                break;
            }
            case (3):
            {
                nodeCoord[0] = -1.;
                nodeCoord[1] =  1.;
                break;
            }
            default:
            {
                DebugStop();
                break;
            }
        }
    }
	
	MElementType TPZQuadrilateral::Type()
	{
		return EQuadrilateral;
	}
	
	MElementType TPZQuadrilateral::Type(int side)
	{
		switch(side) {
			case 0:
			case 1:
			case 2:
			case 3:
				return EPoint;
			case 4:
			case 5:
			case 6:
			case 7:
				return EOned;
			case 8:
				return EQuadrilateral;
			default:
				return ENoType;
		}
	}
	
	TPZTransform<> TPZQuadrilateral::SideToSideTransform(int sidefrom, int sideto)
	{
		if(sidefrom <0 || sidefrom >= NSides || sideto <0 || sideto >= NSides) {
			PZError << "TPZShapeQuad::SideToSideTransform sidefrom "<< sidefrom << 
			' ' << sideto << endl;
			return TPZTransform<>(0);
		}
		if(sidefrom == sideto) {
			return TPZTransform<>(sidedimension[sidefrom]);
		}
		if(sidefrom == NSides-1) {
			return TransformElementToSide(sideto);
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
		PZError << "TPZShapeQuad::SideToSideTransform highside not found sidefrom "
		<< sidefrom << ' ' << sideto << ".\n";
		return TPZTransform<>(0);
	}
	
	
	int TPZQuadrilateral::SideDimension(int side) {
		if(side<0 || side >= NSides) {
			PZError << "TPZShapeQuad::SideDimension. Out of scope side " << side << ".\n";
			return -1;
		}
		return sidedimension[side];
	}
	
	int TPZQuadrilateral::NContainedSides(int side) {
		if(side<0 || side>8) {
			PZError << "TPZShapeQuad::NContainedSides. Bad parameter side = " << side << ".\n";
			return 0;
		}
		if(side<4) return 1;
		if(side<8) return 3;
		return 9;//Cedric
	}
	
	int TPZQuadrilateral::NumSides(int dimension) {
		if(dimension<0 || dimension> 2) {
			PZError << "TPZShapeQuad::NumSides. Bad parameter dimension = " << dimension << ".\n";
			return 0;
		}
		if(dimension==0) return 4;
		if(dimension==1) return 4;
		if(dimension==2) return 1;
		return -1;
		
	}
	
	/**It do not verify the values of the c*/
	// side é o lado do elemento, c é o noh do lado
	int TPZQuadrilateral::ContainedSideLocId(int side,int c) {
		switch(side) {
			case 0:
			case 1:
			case 2:
			case 3:
				return side;
			case 4:
			case 5:
			case 6:
			case 7:
				if(!c) return side-4;
				if(c==1) return (side-3)%4;
				if(c==2) return side;
			case 8:
				return c;
			default:
				PZError << "TPZShapeQuad::ContainedSideLocId, connect = " << c << endl;
				return -1;
		}
	}
	
	
	TPZTransform<> TPZQuadrilateral::TransformSideToElement(int side){
		
		if(side<0 || side>8){
			PZError << "TPZShapeQuad::TransformSideToElement side out range\n";
			return TPZTransform<>(0,0);
		}
		TPZTransform<> t(2,sidedimension[side]);
		t.Mult().Zero();
		t.Sum().Zero();
		
		switch(side){
			case 0:
				t.Sum()(0,0) = -1.0;
				t.Sum()(1,0) = -1.0;
				return t;
			case 1:
				t.Sum()(0,0) =  1.0;
				t.Sum()(1,0) = -1.0;
				return t;
			case 2:
				t.Sum()(0,0) =  1.0;
				t.Sum()(1,0) =  1.0;
				return t;
			case 3:
				t.Sum()(0,0) = -1.0;
				t.Sum()(1,0) =  1.0;
				return t;
			case 4:
				t.Mult()(0,0) =  1.0;
				t.Sum() (1,0) = -1.0;
				return t;
			case 5:
				t.Mult()(1,0) =  1.0;
				t.Sum() (0,0) =  1.0;
				return t;
			case 6:
				t.Mult()(0,0) = -1.0;
				t.Sum() (1,0) =  1.0;
				return t;
			case 7:
				t.Mult()(1,0) = -1.0;
				t.Sum() (0,0) = -1.0;
				return t;
			case 8:
				t.Mult()(0,0) =  1.0;
				t.Mult()(1,1) =  1.0;
				return t;
		}
		return TPZTransform<>(0,0);
	}
	
	
	/**
	 * Method which identifies the transformation based on the IDs
	 * of the corner nodes
	 * @param id indexes of the corner nodes
	 * @return index of the transformation of the point corresponding to the topology
	 */
	int TPZQuadrilateral::GetTransformId(TPZVec<int64_t> &id)
	{
		return pzshape::TPZShapeQuad::GetTransformId2dQ(id);
	}
    
    /**
     * @brief return the vector which permutes the connects according to the transformation id
     */
    void TPZQuadrilateral::GetGatherPermute(int transformid, TPZVec<int> &permute)
    {
#ifdef PZDEBUG
        if (permute.size() != 9) {
            DebugStop();
        }
#endif
        int dir = 1;
        if (transformid%2 ==1) dir = -1;
        int runsmall = transformid/4;
        int runlarge = runsmall+4;
        permute[8] = 8;
        for (int is=0; is<4; is++) {
            permute[is] = runsmall;
            permute[is+4] = runlarge;
            runsmall += dir;
            runlarge += dir;
            if (dir == 1 && runsmall > 3) {
                runsmall -= 4;
                runlarge -= 4;
            }
            else if(dir == -1 && runsmall < 0)
            {
                runsmall += 4;
                runlarge += 4;
            }
        }
    }

	/**
	 * Method which identifies the transformation of a side based on the IDs
	 * of the corner nodes
	 * @param id indexes of the corner nodes
	 * @return index of the transformation of the point corresponding to the topology
	 */	
	int TPZQuadrilateral::GetTransformId(int side, TPZVec<int64_t> &id)
	{
		switch (side) {
			case 0:
			case 1:
			case 2:
			case 3:
				return 0;
				break;
			case 4:
			case 5:
			case 6:
			case 7:
			{
				int in1 = ContainedSideLocId(side,0);
				int in2 = ContainedSideLocId(side,1);
				return id[in1]<id[in2] ? 0 : 1;
			}
				break;
			case 8:
			{
				return pzshape::TPZShapeQuad::GetTransformId2dQ(id);
			}
				break;
			default:
				break;
		}
		LOGPZ_ERROR(logger,"Wrong side parameter")
		return -1;
	}
	
	/**
	 * Identifies the permutation of the nodes needed to make neighbouring elements compatible 
	 * in terms of order of shape functions
	 * @param side : side for which the permutation is needed
	 * @param id : ids of the corner nodes of the elements
	 * @param permgather : permutation vector in a gather order
	 */
	void TPZQuadrilateral::GetSideHDivPermutation(int transformationid, TPZVec<int> &permgather)
	{
//        permgather.Resize(9);
#ifdef PZDEBUG
        if (transformationid < 0 || transformationid > 8 || permgather.size() != 9) {
            DebugStop();
        }
#endif
        
         for (int i=0; i<9; i++)
         {
             permgather[i] = permutationsQ[transformationid][i];
         }
        return;
        int i;
        
        if(transformationid%2 == 0)
        {
            switch (transformationid)
            {
                case 0:
                    for(i=0; i<9; i++) permgather[i] = i;
                    break;
                case 2:
                    for(i=0; i<4; i++) permgather[i] = (i+1)%4;
                    for(i=4; i<8; i++) permgather[i] = 4+(i+1)%4;
                    permgather[8] = 8;
                    break;
                case 4:
                    for(i=0; i<4; i++) permgather[i] = (i+2)%4;
                    for(i=4; i<8; i++) permgather[i] = 4+(i+2)%4;
                    permgather[8] = 8;
                    break;
                case 6:
                    for(i=0; i<4; i++) permgather[i] = (i+3)%4;
                    for(i=4; i<8; i++) permgather[i] = 4+(i+3)%4;
                    permgather[8] = 8;
                    break;
            }
        }
        else
        {
            TPZManVector<int,4> invid(4);
            invid[0] = 0;
            invid[1] = 3;
            invid[2] = 2;
            invid[3] = 1;
            switch (transformationid) {
                case 1:
                    for(i=0; i<4; i++) permgather[i] = invid[i];
                    for(i=4; i<8; i++) permgather[i] = 4+invid[(i+1)%4];
                    permgather[8] = 8;
                    break;
                case 3:
                    for(i=0; i<4; i++) permgather[i] = invid[(i+1)%4];
                    for(i=4; i<8; i++) permgather[i] = 4+invid[(i+2)%4];
                    permgather[8] = 8;
                    break;
                case 5:
                    for(i=0; i<4; i++) permgather[i] = invid[(i+2)%4];
                    for(i=4; i<8; i++) permgather[i] = 4+invid[(i+3)%4];
                    permgather[8] = 8;
                    break;
                case 7:
                    for(i=0; i<4; i++) permgather[i] = invid[(i+3)%4];
                    for(i=4; i<8; i++) permgather[i] = 4+invid[(i+0)%4];
                    permgather[8] = 8;
                    break;
                default:
                    break;
            }
        }

        

        
        /*
        switch (transformationid)
        {
            case 0:
            {
                for (int i=0; i<9; i++)
                {
                    permgather[i] = permutationsQ[0][i];
                }
            }
                break;
            case 1:
            {
                for (int i=0; i<9; i++)
                {
                    permgather[i] = permutationsQ[1][i];
                }
            }
                break;
            case 2:
            {
                for (int i=0; i<9; i++)
                {
                    permgather[i] = permutationsQ[2][i];
                }
            }
                break;
            case 3:
            {
                for (int i=0; i<9; i++)
                {
                    permgather[i] = permutationsQ[3][i];
                }
            }
                break;
            case 4:
            {
                for (int i=0; i<9; i++)
                {
                    permgather[i] = permutationsQ[4][i];
                }
            }
                break;
            case 5:
            {
                for (int i=0; i<9; i++)
                {
                    permgather[i] = permutationsQ[5][i];
                }
            }
                break;
            case 6:
            {
                for (int i=0; i<9; i++)
                {
                    permgather[i] = permutationsQ[6][i];
                }
            }
                break;
            case 7:
            {
                for (int i=0; i<9; i++)
                {
                    permgather[i] = permutationsQ[7][i];
                }
            }
                break;
            default:
                DebugStop();
                break;
        }*/

        
        /*
		switch (side) {
			case 0:
			case 1:
			case 2:
			case 3:
				permgather[0] = 0;
				break;
			case 4:
			case 5:
			case 6:
			case 7:
			{
				int in1 = ContainedSideLocId(side,0);
				int in2 = ContainedSideLocId(side,1);
				if(in1<in2)
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
			}
				break;
			case 8:
			{
				int i;
				int tid = pzshape::TPZShapeQuad::GetTransformId2dQ(id);
				if(tid%2 == 0)
				{
					switch (tid)
					{
						case 0:
							for(i=0; i<9; i++) permgather[i] = i;
							break;
						case 2:
							for(i=0; i<4; i++) permgather[i] = (i+1)%4;
							for(i=4; i<8; i++) permgather[i] = 4+(i+1)%4;
							permgather[8] = 8;
							break;
						case 4:
							for(i=0; i<4; i++) permgather[i] = (i+2)%4;
							for(i=4; i<8; i++) permgather[i] = 4+(i+2)%4;
							permgather[8] = 8;
							break;
						case 6:
							for(i=0; i<4; i++) permgather[i] = (i+3)%4;
							for(i=4; i<8; i++) permgather[i] = 4+(i+3)%4;
							permgather[8] = 8;
							break;
					}
				}
				else
				{
					TPZManVector<int,4> invid(4);
					invid[0] = 0;
					invid[1] = 3;
					invid[2] = 2;
					invid[3] = 1;
					switch (tid) {
						case 1:
							for(i=0; i<4; i++) permgather[i] = invid[i];
							for(i=4; i<8; i++) permgather[i] = 4+invid[(i+0)%4];
							permgather[8] = 8;
							break;
						case 3:
							for(i=0; i<4; i++) permgather[i] = invid[(i+1)%4];
							for(i=4; i<8; i++) permgather[i] = 4+invid[(i+1)%4];
							permgather[8] = 8;
							break;
						case 5:
							for(i=0; i<4; i++) permgather[i] = invid[(i+2)%4];
							for(i=4; i<8; i++) permgather[i] = 4+invid[(i+2)%4];
							permgather[8] = 8;
							break;
						case 7:
							for(i=0; i<4; i++) permgather[i] = invid[(i+3)%4];
							for(i=4; i<8; i++) permgather[i] = 4+invid[(i+3)%4];
							permgather[8] = 8;
							break;
						default:
							break;
					}
				}
			}
				break;
			default:
				break;
		}
		LOGPZ_ERROR(logger,"Wrong side parameter")*/
	}
    
    void computedirectionsq(int inicio, int fim, TPZFMatrix<REAL> &bvec, TPZFMatrix<REAL> &t1vec,
                           TPZFMatrix<REAL> &gradx, TPZFMatrix<REAL> &directions);
    void computedirectionsq(int inicio, int fim, TPZFMatrix<REAL> &bvec, TPZFMatrix<REAL> &t1vec,
                           TPZFMatrix<REAL> &gradx, TPZFMatrix<REAL> &directions)
    {
        REAL detgrad = 0.0;
        TPZVec<REAL> u(3);
        TPZVec<REAL> v(3);
        TPZVec<REAL> uxv(3);// result
        
        for (int ilin=0; ilin<3; ilin++)
        {
            u[ilin] = gradx(ilin, 0);
            v[ilin] = gradx(ilin, 1);
        }
        
        //TPZNumeric::ProdVetorial(u,v,uxv);
        uxv[0] = u[1]*v[2]-u[2]*v[1];
        uxv[1] = -(u[0]*v[2]-v[0]*u[2]);
        uxv[2] = u[0]*v[1]-v[0]*u[1];
        
        for (int pos=0; pos<3; pos++)
        {
            detgrad += uxv[pos]*uxv[pos];
        }
        detgrad = sqrt(fabs(detgrad));
        
        int cont = 0;
        
        for (int ivet=inicio; ivet<=fim; ivet++)
        {
            TPZFMatrix<REAL> Wvec(3,1);
            TPZVec<REAL> uxvtmp(3);
            REAL acumng = 0.0;
            // calc do g gradx*t
            TPZManVector<REAL,3> gvec(3,0.),Vvec(3,0.);
            REAL gvecnorm;
            for (int il=0; il<3; il++)
            {
                for (int i = 0 ; i<2; i++)
                {
                    gvec[il] += gradx(il,i) * t1vec(i,ivet);
                    Vvec[il] += gradx(il,i) * bvec(i,ivet);
                }
                u[il] = gvec[il];
                acumng += gvec[il]*gvec[il];
            }
            gvecnorm = sqrt(acumng);
            
            for (int il=0; il<3; il++)
            {
                Wvec(il,0) = Vvec[il]*gvecnorm/detgrad;
                directions(il,cont) = Wvec(il,0);
            }
            cont++;
        }
        
        
    }

    void TPZQuadrilateral::ComputeDirections(int side, TPZFMatrix<REAL> &gradx, TPZFMatrix<REAL> &directions, TPZVec<int> &sidevectors)
    {
        if(gradx.Cols()!=2)
        {
            DebugStop();
        }
        TPZFMatrix<REAL> bvec(2,18);
        TPZFMatrix<REAL> t1vec(2,18);
       
        bvec.Redim(2, 18);
        t1vec.Redim(2, 18);
        directions.Redim(3, 18);
        for (int lin = 0; lin<18; lin++)
        {
            for(int col = 0;col<2;col++)
            {
                bvec.PutVal(col, lin, bQuad[lin][col]);
                t1vec.PutVal(col, lin, tQuad[lin][col]);
            }
        }
        // calcula os vetores
        
        switch (side) {
            case 0:
            {
            }
                break;
            case 1:
            {
            }
                break;
            case 2:
            {
            }
                break;
            case 3:
            {
            }
                break;
            case 4:
            {
                directions.Redim(3, 3);
                sidevectors.Resize(3);
                int inumvec = 0, fnumvec = 2;
                computedirectionsq(inumvec, fnumvec, bvec, t1vec, gradx, directions);
                for (int ip = 0; ip < 3; ip++) {
                    sidevectors[ip] = vectorsideorder[ip];
                }
                
            }
                break;
            case 5:
            {
                directions.Redim(3, 3);
                sidevectors.Resize(3);
                int inumvec = 3, fnumvec = 5;
                computedirectionsq(inumvec, fnumvec, bvec, t1vec, gradx, directions);
                for (int ip = 0; ip < 3; ip++) {
                    sidevectors[ip] = vectorsideorder[ip+inumvec];
                }
            }
                break;
            case 6:
            {
                directions.Redim(3, 3);
                sidevectors.Resize(3);
                int inumvec = 6, fnumvec = 8;
                computedirectionsq(inumvec, fnumvec, bvec, t1vec, gradx, directions);
                for (int ip = 0; ip < 3; ip++) {
                    sidevectors[ip] = vectorsideorder[ip+inumvec];
                }
            }
                break;
            case 7:
            {
                directions.Redim(3, 3);
                sidevectors.Resize(3);
                int inumvec = 9, fnumvec = 11;
                computedirectionsq(inumvec, fnumvec, bvec, t1vec, gradx, directions);
                for (int ip = 0; ip < 3; ip++) {
                    sidevectors[ip] = vectorsideorder[ip+inumvec];
                }
            }
                break;
            case 8:
            {
                directions.Redim(3, 6);
                sidevectors.Resize(6);
                int inumvec = 12, fnumvec = 17;
                computedirectionsq(inumvec, fnumvec, bvec, t1vec, gradx, directions);
                for (int ip = 0; ip < 6; ip++) {
                    sidevectors[ip] = vectorsideorder[ip+inumvec];
                }
            }
                break;
                
            default:
                DebugStop();
                break;
        }

        
	}
    
    void TPZQuadrilateral::ComputeDirections(TPZFMatrix<REAL> &gradx, REAL detjac, TPZFMatrix<REAL> &directions)
    {
        TPZManVector<REAL, 3> v1(3),v2(3);
        for (int i=0; i<3; i++) {
            v1[i] = gradx(i,0);
            v2[i] = gradx(i,1);
        }

        
        REAL Nv1 = TPZNumeric::Norma(v1);
        REAL Nv2 = TPZNumeric::Norma(v2);
        
        /**
         * @file
         * @brief Computing mapped vector with scaling factor equal 1.0.
         * using contravariant piola mapping.
         */
        TPZManVector<REAL,3> NormalScales(2,1.);
        
        if (HDivPiola == 1)
        {
            NormalScales[0] = 1./Nv1;
            NormalScales[1] = 1./Nv2;
        }
        
        
        for (int i=0; i<3; i++) {
            v1[i] *= Nv2/detjac;
            v2[i] *= Nv1/detjac;
        }
        
        for (int i=0; i<3; i++)
        {
            for (int v=0; v<3; v++)
            {
                directions(i,v)     = -v2[i]*NormalScales[0];
                directions(i,v+3)   = v1[i]*NormalScales[1];
                directions(i,v+6)   = v2[i]*NormalScales[0];
                directions(i,v+9)   = -v1[i]*NormalScales[1];
            }
            
            directions(i,12)        =  v1[i]*NormalScales[1];
            directions(i,13)        =  v2[i]*NormalScales[0];
            directions(i,14)        = -v1[i]*NormalScales[1];
            directions(i,15)        = -v2[i]*NormalScales[0];
            
            directions(i,16)        = v1[i]*NormalScales[1];
            directions(i,17)        = v2[i]*NormalScales[0];
        }
    }
    
    void TPZQuadrilateral::GetSideDirections(TPZVec<int> &sides, TPZVec<int> &dir, TPZVec<int> &bilounao)
    {
        int nsides = NumSides()*2;
        
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
    void TPZQuadrilateral::GetSideDirections(TPZVec<int> &sides, TPZVec<int> &dir, TPZVec<int> &bilounao, TPZVec<int> &sidevectors)
    {
        int nsides = NumSides()*2;
        
        sides.Resize(nsides);
        dir.Resize(nsides);
        bilounao.Resize(nsides);
        
        for (int is = 0; is<nsides; is++)
        {
            sides[is] = vectorsideorder[is];
            dir[is] = direcaoksioueta[is];
            bilounao[is] = bilinearounao[is];
        }
        for (int i=0; i<Dimension*NumSides(); i++) {
            sidevectors[i] = vectorsideorder[i];
        }
    }
    
    int TPZQuadrilateral::ClassId() const{
        return Hash("TPZQuadrilateral");
    }

    void TPZQuadrilateral::Read(TPZStream& buf, void* context) {

    }
    
    void TPZQuadrilateral::Write(TPZStream& buf, int withclassid) const {

    }

}

template
bool pztopology::TPZQuadrilateral::MapToSide<REAL>(int side, TPZVec<REAL> &InternalPar, TPZVec<REAL> &SidePar, TPZFMatrix<REAL> &JacToSide);

#ifdef _AUTODIFF
template
bool pztopology::TPZQuadrilateral::MapToSide<Fad<REAL> >(int side, TPZVec<Fad<REAL> > &InternalPar, TPZVec<Fad<REAL> > &SidePar, TPZFMatrix<Fad<REAL> > &JacToSide);
#endif
