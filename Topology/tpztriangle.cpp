/**
 * @file
 * @brief Contains the implementation of the TPZTriangle methods. 
 */

#include "tpztriangle.h"
#include "pzquad.h"

#include "pzshapetriang.h"
#include "pzcreateapproxspace.h"

#ifdef _AUTODIFF
#include "fad.h"
#endif

#include "pzlog.h"
#include <cmath>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.topology.pztriangle"));
#endif

using namespace std;

namespace pztopology {
	
	static int sidedimension[7] = {0,0,0,1,1,1,2};
	
	static int nhighdimsides[7] = {3,3,3,1,1,1,0};
	
	static int highsides[7][3] = {
		{3,5,6},
		{3,4,6},
		{4,5,6},
		{6},
		{6},
		{6},
		{-999}
	};
	
	static REAL sidetosidetransforms[7][3][4][3] = {
		{
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,0,-99}}
		},
		{
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}}
		},
		{
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}}
		},
		{
			{{0.5,0,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0,-99}}
		},
		{
			{{-0.5,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0.5,-99}}
		},
		{
			{{0,-0.5,-99},{-99,-99,-99},{-99,-99,-99},{0,0.5,-99}}
		},
		{
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-99,-99,-99}}
		}
	};
	
	static REAL MidSideNode[7][3] = {
		/*00*/{.0,0.},/*01*/{1.0,.0},/*02*/{0.,1.0},
		/*03*/{.5,0.},/*04*/{0.5,.5},/*05*/{0.,0.5},
		/*06*/{ 1./3.,1./3.} };
	
	static int nsidenodes[7] = {1,1,1,2,2,2,3};
    
    static REAL bTriang[14][2] = 
    {
        {0,-1},//0
        {1,-1},
        {0,-1},
        {1,0},
        {0,1},
        {0.5,0.5},//5
        {-1,1},
        {-1,0},
        {-1,0},//8
        {1,0},
        {-0.5,0.5},
        {0,-1},
        {1,0},
        {0,1}
        
    };
    
    static REAL tTriang[14][2] = 
    {
        {-1,0},
        {-1,0},
        {-1,0},
        {1,-1},
        {1,-1},
        {1,-1},
        {0,1},
        {0,1},
        {0,1},//
        {0,-1},
        {1,1},
        {-1,0},
        {0,-1},
        {1,0}
    };

    static int vectorsideorder [14] = {0,1,3,1,2,4,2,0,5,3,4,5,6,6};
    
    //static int bilinearounao [14] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0};//Pk Pk-1
    static int bilinearounao [14] = {0,0,0,0,0,0,0,0,0,1,1,1,1,1};//P*k Pk
    
    static int direcaoksioueta [14] = {0,0,0,0,0,0,0,0,0,0,0,0,0,1};
    
    static int permutationsT [6][7] =
    {
        {0,1,2,3,4,5,6}, // id 0
        {0,2,1,5,4,3,6}, // id 1
        {1,2,0,4,5,3,6}, // id 2
        {1,0,2,3,5,4,6}, // id 3
        {2,0,1,5,3,4,6}, // id 4
        {2,1,0,4,3,5,6}  // id 5
    };

	int TPZTriangle::NBilinearSides()
    {
        DebugStop();
        return 5;
    }
    
    
	int TPZTriangle::NSideNodes(int side)
	{
		return nsidenodes[side];
	}
	
	bool TPZTriangle::IsInParametricDomain(TPZVec<REAL> &pt, REAL tol){
		const REAL qsi = pt[0];
		const REAL eta = pt[1];
		if( ( qsi <= 1. + tol ) && ( qsi >= 0. - tol ) &&
		   ( eta <= 1. + tol ) && ( eta >= 0. - tol ) &&
		   ( eta <= 1. - qsi + tol ) ){
			return true;
		}
		else{
			return false;
		}  
	}//method
    
    template<class T>
    bool TPZTriangle::MapToSide(int side, TPZVec<T> &InternalPar, TPZVec<T> &SidePar, TPZFMatrix<T> &JacToSide) {
		
		double zero = 1.E-5;
		T qsi = InternalPar[0]; T eta = InternalPar[1];
		SidePar.Resize(1); JacToSide.Resize(1,2);

		bool regularmap = true;
		
		switch(side)
		{
            case 0:
            case 1:
            case 2:
            {
                SidePar.Resize(0); JacToSide.Resize(0,0);
                break;
            }
			case 3:
				if(fabs((T)(eta - 1.)) < zero)
				{
                    SidePar[0] = 0.;
                    JacToSide(0,0) = 0.; JacToSide(0,1) = 0.;
					regularmap = false;
				}
				else
				{
                    SidePar[0] = 2.*qsi/(1.-eta) - 1.;
                    JacToSide(0,0) = 2./(1.-eta); JacToSide(0,1) = 2.*qsi/((1.-eta)*(1.-eta));
				}
				break;
				
			case 4:
				if((T)(qsi+eta) < (T)zero)
				{
                    SidePar[0] = 0.;
                    JacToSide(0,0) = 0.; JacToSide(0,1) = 0.;
					regularmap = false;
				}
				else
				{
                    SidePar[0] = 1. - 2.*qsi/(qsi + eta);
                    JacToSide(0,0) = -2.*eta/((qsi+eta)*(qsi+eta)); JacToSide(0,1) = 2.*qsi/((qsi+eta)*(qsi+eta));
				}
				break;
				
			case 5:
				if(fabs((T)(qsi - 1.)) < zero)
				{
                    SidePar[0] = 0.;
                    JacToSide(0,0) = 0.; JacToSide(0,1) = 0.;
					regularmap = false;
				}
				else
				{
                    SidePar[0] = 1. - 2.*eta/(1.-qsi);
                    JacToSide(0,0) = -2.*eta/((1.-qsi)*(1.-qsi)); JacToSide(0,1) = -2./(1.-qsi);
				}
				break;
            case 6:
                SidePar = InternalPar;
                JacToSide.Resize(2, 2);
                JacToSide.Identity();
                regularmap = true;
		}
		if(side > 6)
		{
			cout << "Cant compute MapToSide method in TPZGeoTriangle class!\nParameter (SIDE) must be 3, 4 or 5!\nMethod Aborted!\n";
			DebugStop();
		}
		return regularmap;
	}
    
    void TPZTriangle::ParametricDomainNodeCoord(int node, TPZVec<REAL> &nodeCoord)
    {
        if(node > NCornerNodes)
        {
            DebugStop();
        }
        nodeCoord.Resize(Dimension, 0.);
        switch (node) {
            case (0):
            {
                nodeCoord[0] = 0.;
                nodeCoord[1] = 0.;
                break;
            }
            case (1):
            {
                nodeCoord[0] = 1.;
                nodeCoord[1] = 0.;
                break;
            }
            case (2):
            {
                nodeCoord[0] = 0.;
                nodeCoord[1] = 1.;
                break;
            }
            default:
            {
                DebugStop();
                break;
            }
        }
    }
	
    /** @brief Generates a random point in the master domain */
    void TPZTriangle::RandomPoint(TPZVec<REAL> &pt)
    {
        REAL val = (REAL) rand() / (RAND_MAX);
        pt[0] = val;
        val = (1.-pt[0]) * (REAL) rand() / (RAND_MAX);
        pt[1] = val;
    }
    

	TPZTransform<> TPZTriangle::SideToSideTransform(int sidefrom, int sideto)
	{
		if(sidefrom <0 || sidefrom >= NSides || sideto <0 || sideto >= NSides) {
			PZError << "TPZTriangle::HigherDimensionSides sidefrom "<< sidefrom << 
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
		PZError << "TPZTriangle::SideToSideTransform highside not found sidefrom "
		<< sidefrom << ' ' << sideto << endl;
		return TPZTransform<>(0);
	}
	
	int TPZTriangle::SideNodeLocId(int side, int node)
	{
		if(side<3 && node == 0) return side;
		if(side>=3 && side<6 && node <2) return (side-3+node) %3;
		if(side==6 && node <3) return node;
		PZError << "TPZTriangle::SideNodeLocId inconsistent side or node " << side
		<< ' ' << node << endl;
		return -1;
	}
	
	void TPZTriangle::LowerDimensionSides(int side,TPZStack<int> &smallsides)
	{
		smallsides.Resize(0);
		int nsidecon = NContainedSides(side);
		int is;
		for(is=0; is<nsidecon-1; is++)
			smallsides.Push(ContainedSideLocId(side,is));
	}
	
	void TPZTriangle::LowerDimensionSides(int side,TPZStack<int> &smallsides, int DimTarget)
	{
		smallsides.Resize(0);
		int nsidecon = NContainedSides(side);
		for(int is = 0; is < nsidecon - 1; is++) {
			if (SideDimension(ContainedSideLocId(side,is)) == DimTarget) smallsides.Push(ContainedSideLocId(side,is));
		}
	}
	
	void TPZTriangle::HigherDimensionSides(int side, TPZStack<int> &high)
	{
		if(side <0 || side >= NSides) {
			PZError << "TPZTriangle::HigherDimensionSides side "<< side << endl;
		}
		int is;
		for(is=0; is<nhighdimsides[side]; is++) high.Push(highsides[side][is]);
		
	}
	
	//Tentando criar o metodo
	int TPZTriangle::NumSides(int dimension) {	
		if(dimension<0 || dimension> 2) {
			PZError << "TPZTriangle::NumSides. Bad parameter i.\n";
			return 0;
		}
		if(dimension==0) return 3;
		if(dimension==1) return 3;
		if(dimension==2) return 1;
		return -1;
	}
	int TPZTriangle::SideDimension(int side) {
		if(side<0 || side >= NSides) {
			PZError << "TPZTriangle::SideDimension side " << side << endl;
			return -1;
		}
		return sidedimension[side];
	}
	
	void TPZTriangle::CenterPoint(int side, TPZVec<REAL> &center) {
		//center.Resize(Dimension);
		int i;
		for(i=0; i<Dimension; i++) {
			center[i] = MidSideNode[side][i];
		}
	}
	
	TPZTransform<> TPZTriangle::TransformElementToSide(int side) {
		if(side<0 || side>6){
			PZError << "TPZTriangle::TransformElementToSide called with side error\n";
			return TPZTransform<>(0,0);
		}
		
		TPZTransform<> t(sidedimension[side],2);//t(dimto,2)
		t.Mult().Zero();
		t.Sum().Zero();
		
		switch(side){
			case 0:
			case 1:
			case 2:
				return t;
			case 3:
				t.Mult()(0,0) =  2.0;//par. var.
				t.Sum()(0,0)  = -1.0;
				return t;
			case 4:
				t.Mult()(0,0) = -1.0;
				t.Mult()(0,1) =  1.0;
				return t;
			case 5:
				t.Mult()(0,1) = -2.0;
				t.Sum()(0,0)  =  1.0;
				return t;
			case 6:
				t.Mult()(0,0) =  1.0;
				t.Mult()(1,1) =  1.0;
				return t;
		}
		return TPZTransform<>(0,0);
	}
	
	TPZTransform<> TPZTriangle::TransformSideToElement(int side){
		
		if(side<0 || side>6){
			PZError << "TPZTriangle::TransformSideToElement side out range\n";
			return TPZTransform<>(0,0);
		}
		TPZTransform<> t(2,sidedimension[side]);
		t.Mult().Zero();
		t.Sum().Zero();
		
		switch(side){
			case 0:
				return t;
			case 1:
				t.Sum()(0,0) =  1.0;
				return t;
			case 2:
				t.Sum()(1,0) =  1.0;
				return t;
			case 3:
				t.Mult()(0,0) =  0.5;
				t.Sum()(0,0)  =  0.5;
				return t;
			case 4:
				t.Mult()(0,0) = -0.5;
				t.Mult()(1,0) =  0.5;
				t.Sum() (0,0) =  0.5;
				t.Sum() (1,0) =  0.5;
				return t;
			case 5:
				t.Mult()(1,0) = -0.5;
				t.Sum() (1,0) =  0.5;
				return t;
			case 6:
				t.Mult()(0,0) =  1.0;
				t.Mult()(1,1) =  1.0;
				return t;
		}
		return TPZTransform<>(0,0);
	}
	
	TPZIntPoints * TPZTriangle::CreateSideIntegrationRule(int side, int order){
		if(side < 0 || side>6) {
			PZError << "TPZTriangle::CreateSideIntegrationRule wrong side = " << side << ".\n";
			return 0;
		}
		if(side<3) return new TPZInt1Point(order);     // sides 0 to 2 are vertices (corners)
		if(side<6) return new TPZInt1d(order);    // sides 3 to 5 are lines
		if(side==6)return new IntruleType(order); // integration of the element
		return 0;
	}
	
	
	MElementType TPZTriangle::Type()
	{
		return ETriangle;
	}
	
	MElementType TPZTriangle::Type(int side)
	{
		switch(side) {
			case 0:
			case 1:
			case 2:
				return EPoint;
			case 3:
			case 4:
			case 5:
				return EOned;
			case 6:
				return ETriangle;
			default:
				return ENoType;
		}
	}
	
	
	int TPZTriangle::NumSides() {
		return NSides;
	}
	
	
	
	int TPZTriangle::NContainedSides(int side) {
		if(side<0 || side>6) {
			PZError << "TPZShapeTriang::NContainedSides. Bad parameter side = " << side << ".\n";
			return 0;
		}
		if(side<3) return 1;
		if(side<6) return 3;
		return 7;
	}
	
	/**It do not verify the values of the c*/
	// side eh o lado do elemento, c eh o noh do lado
	int TPZTriangle::ContainedSideLocId(int side, int c) {
		switch(side) {
			case 0:
			case 1:
			case 2:
				return side;
			case 3:
			case 4:
			case 5:
				if(!c) return side-3;
				if(c==1) return (side-2)%3;
				if(c==2) return side;
			case 6:
				return c;
			default:
				PZError << "TPZShapeTriang::ContainedSideLocId, connect = " << c << ".\n";
				return -1;
		}
	}
	
	/**
	 * Method which identifies the transformation based on the IDs
	 * of the corner nodes
	 * @param id indexes of the corner nodes
	 * @return index of the transformation of the point corresponding to the topology
	 */
	int TPZTriangle::GetTransformId(TPZVec<int64_t> &id)
	{
		return pzshape::TPZShapeTriang::GetTransformId2dT(id);
	}
	/**
	 * Method which identifies the transformation of a side based on the IDs
	 * of the corner nodes
	 * @param id indexes of the corner nodes
	 * @return index of the transformation of the point corresponding to the topology
	 */	
	int TPZTriangle::GetTransformId(int side, TPZVec<int64_t> &id)
	{
		switch (side) {
			case 0:
			case 1:
			case 2:
				return 0;
				break;
			case 3:
			case 4:
			case 5:
			{
				int in1 = ContainedSideLocId(side,0);
				int in2 = ContainedSideLocId(side,1);
				return id[in1]<id[in2] ? 0 : 1;
			}
				break;
			case 6:
			{
				return pzshape::TPZShapeTriang::GetTransformId2dT(id);
			}
				break;
			default:
				break;
		}
		LOGPZ_ERROR(logger,"Wrong side parameter")
		return -1;
	}
    
/**
 * @brief return the vector which permutes the connects according to the transformation id
 */
void TPZTriangle::GetHDivGatherPermute(int transformid, TPZVec<int> &permute)
{
#ifdef PZDEBUG
    if (permute.size() != 7 || transformid >= 6 || transformid < 0) {
        DebugStop();
    }
#endif
    for (int i=0; i<7; i++) {
        permute[i] = permutationsT[transformid][i];
    }
    return;
    int dir = 1;
    if (transformid%2 ==1) dir = -1;
    int runsmall = transformid/2;
    int runlarge = runsmall+3;
    if (dir == -1) {
        runlarge += dir;
        if (runlarge < 3) {
            runlarge += 3;
        }
    }
    permute[6] = 6;
    for (int is=0; is<3; is++) {
        permute[is] = runsmall;
        permute[is+3] = runlarge;
        runsmall += dir;
        runlarge += dir;
        if (dir == 1 && runsmall > 2) {
            runsmall -= 3;
            runlarge -= 3;
        }
        else if(dir == -1 && runsmall < 0)
        {
            runsmall += 3;
        }
        if (dir == -1 && runlarge < 3) {
            runlarge += 3;
        }
    }
}
    
	
	/**
	 * Identifies the permutation of the nodes needed to make neighbouring elements compatible 
	 * in terms of order of shape functions
	 * @param side : side for which the permutation is needed
	 * @param id : ids of the corner nodes of the elements
	 * @param permgather : permutation vector in a gather order
	 */
	void TPZTriangle::GetSideHDivPermutation(int transformationid, TPZVec<int> &permgather)
	{
        permgather.Resize(7);
        
        for (int i=0; i<7; i++)
        {
            permgather[i] = permutationsT[transformationid][i];
        }
        
        /*
        switch (transformationid)
        {
            case 0:
            {
                for (int i=0; i<7; i++)
                {
                    permgather[i] = permutationsT[0][i];
                }
            }
                break;
            case 1:
            {
                for (int i=0; i<7; i++)
                {
                    permgather[i] = permutationsT[1][i];
                }
            }
                break;
            case 2:
            {
                for (int i=0; i<7; i++)
                {
                    permgather[i] = permutationsT[2][i];
                }
            }
                break;
            case 3:
            {
                for (int i=0; i<7; i++)
                {
                    permgather[i] = permutationsT[3][i];
                }
            }
                break;
            case 4:
            {
                for (int i=0; i<7; i++)
                {
                    permgather[i] = permutationsT[4][i];
                }
            }
                break;
            case 5:
            {
                for (int i=0; i<7; i++)
                {
                    permgather[i] = permutationsT[5][i];
                }
            }
                break;
            default:
                DebugStop();
                break;
        }
        */
        
		/*switch (side) {
			case 0:
			case 1:
			case 2:
				permgather[0] = 0;
				break;
			case 3:
			case 4:
			case 5:
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
			case 6:
			{
				int i;
				int tid = pzshape::TPZShapeTriang::GetTransformId2dT(id);
				if(tid%2 == 0)
				{
					switch (tid)
					{
						case 0:
							for(i=0; i<7; i++) permgather[i] = i;
							break;
						case 2:
							for(i=0; i<3; i++) permgather[i] = (i+1)%3;
							for(i=4; i<6; i++) permgather[i] = 3+(i+1)%3;
							permgather[6] = 6;
							break;
						case 4:
							for(i=0; i<3; i++) permgather[i] = (i+2)%3;
							for(i=4; i<6; i++) permgather[i] = 3+(i+2)%3;
							permgather[6] = 6;
							break;
					}
				}
				else
				{
					TPZManVector<int,3> invid(3);
					invid[0] = 0;
					invid[1] = 2;
					invid[2] = 1;
					switch (tid) {
						case 1:
							for(i=0; i<3; i++) permgather[i] = invid[i];
							for(i=4; i<6; i++) permgather[i] = 3+invid[(i+0)%3];
							permgather[6] = 6;
							break;
						case 3:
							for(i=0; i<3; i++) permgather[i] = invid[(i+1)%3];
							for(i=4; i<6; i++) permgather[i] = 3+invid[(i+1)%3];
							permgather[6] = 6;
							break;
						case 5:
							for(i=0; i<3; i++) permgather[i] = invid[(i+2)%3];
							for(i=4; i<6; i++) permgather[i] = 3+invid[(i+2)%3];
							permgather[6] = 6;
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
    
    
    
    void computedirectionst(int inicio, int fim, TPZFMatrix<REAL> &bvec, TPZFMatrix<REAL> &t1vec,
                            TPZFMatrix<REAL> &gradx, TPZFMatrix<REAL> &directions);
    void computedirectionst(int inicio, int fim, TPZFMatrix<REAL> &bvec, TPZFMatrix<REAL> &t1vec,
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
    
    void TPZTriangle::ComputeDirections(int side, TPZFMatrix<REAL> &gradx, TPZFMatrix<REAL> &directions, TPZVec<int> &sidevectors)
    {
    
        if(gradx.Cols()!=2)
        {
            DebugStop();
        }
        
        TPZFMatrix<REAL> bvec(2,14);
        TPZFMatrix<REAL> t1vec(2,14);

        bvec.Redim(2, 14);
        t1vec.Redim(2, 14);
        
        for (int lin = 0; lin<14; lin++)
        {
            for(int col = 0;col<2;col++)
            {
                bvec.PutVal(col, lin, bTriang[lin][col]);
                t1vec.PutVal(col, lin, tTriang[lin][col]);
            }
        }
        // calcula os vetores
        //int numvec = bvec.Cols();
        
        switch (side)
        {
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
                directions.Redim(3, 3);
                sidevectors.Resize(3);
                int inumvec = 0, fnumvec = 2;
                computedirectionst(inumvec, fnumvec, bvec, t1vec, gradx, directions);
                for (int ip = 0; ip < 3; ip++) {
                    sidevectors[ip] = vectorsideorder[ip+inumvec];
                }
                
            }
                break;
            case 4:
            {
                directions.Redim(3, 3);
                sidevectors.Resize(3);
                int inumvec = 3, fnumvec = 5;
                computedirectionst(inumvec, fnumvec, bvec, t1vec, gradx, directions);
                for (int ip = 0; ip < 3; ip++) {
                    sidevectors[ip] = vectorsideorder[ip+inumvec];
                }
            }
                break;
            case 5:
            {
                directions.Redim(3, 3);
                sidevectors.Resize(3);
                int inumvec = 6, fnumvec = 8;
                computedirectionst(inumvec, fnumvec, bvec, t1vec, gradx, directions);
                for (int ip = 0; ip < 3; ip++) {
                    sidevectors[ip] = vectorsideorder[ip+inumvec];
                }
            }
                break;
            case 6:
            {
                directions.Redim(3, 5);
                sidevectors.Resize(5);
                int inumvec = 9, fnumvec = 13;
                computedirectionst(inumvec, fnumvec, bvec, t1vec, gradx, directions);
                for (int ip = 0; ip < 5; ip++) {
                    sidevectors[ip] = vectorsideorder[ip+inumvec];
                }
            }
                break;
            default:
                DebugStop();
                break;
        }
        
    }
    void TPZTriangle::ComputeDirections(TPZFMatrix<REAL> &gradx, REAL detjac, TPZFMatrix<REAL> &directions)
    {
        TPZManVector<REAL, 3> v1(3),v2(3), vdiag(3);
        for (int i=0; i<3; i++) {
            v1[i] = gradx(i,0);
            v2[i] = gradx(i,1);
            vdiag[i] = (gradx(i,0)-gradx(i,1));
       }
        
        REAL Nv1 = TPZNumeric::Norma(v1);
        REAL Nv2 = TPZNumeric::Norma(v2);
        REAL Nvdiag = TPZNumeric::Norma(vdiag);
        
        
        /**
         * @file
         * @brief Computing mapped vector with scaling factor equal 1.0.
         * using contravariant piola mapping.
         */
        TPZManVector<REAL,3> NormalScales(3,1.);
        
        if (HDivPiola == 1)
        {
            NormalScales[0] = 2./Nv1;
            NormalScales[1] = 2./Nvdiag;
            NormalScales[2] = 2./Nv2;
        }
        
        for (int i=0; i<3; i++) {
            v1[i] /= detjac;
            v2[i] /= detjac;
            
        }
        
        for (int i=0; i<3; i++)
        {
            directions(i,0) = -v2[i]*Nv1*NormalScales[0];
            directions(i,1) = (v1[i]-v2[i])*Nv1*NormalScales[0];
            directions(i,2) = (directions(i,0)+directions(i,1))/2.;
            directions(i,3) = v1[i]*Nvdiag*NormalScales[1];
            directions(i,4) = v2[i]*Nvdiag*NormalScales[1];
            directions(i,5) = (directions(i,3)+directions(i,4))/2.;
            directions(i,6) = (v2[i]-v1[i])*Nv2*NormalScales[2];
            directions(i,7) = -v1[i]*Nv2*NormalScales[2];
            directions(i,8) = (directions(i,6)+directions(i,7))/2.;

            directions(i,9) = v1[i]*Nv1*NormalScales[0];
            directions(i,10) = (v2[i]-v1[i])*Nvdiag*NormalScales[1];
            directions(i,11) = -v2[i]*Nv2*NormalScales[2];

            directions(i,12) = v1[i]*Nv2*NormalScales[0];
            directions(i,13) = v2[i]*Nv1*NormalScales[1];
        }
    }
    
    void TPZTriangle::GetSideDirections(TPZVec<int> &sides, TPZVec<int> &dir, TPZVec<int> &bilounao)
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

    void TPZTriangle::GetSideDirections(TPZVec<int> &sides, TPZVec<int> &dir, TPZVec<int> &bilounao, TPZVec<int> &sidevectors)
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

    int TPZTriangle::ClassId() const{
        return Hash("TPZTriangle");
    }

    void TPZTriangle::Read(TPZStream& buf, void* context) {

    }

    void TPZTriangle::Write(TPZStream& buf, int withclassid) const {

    }
   
}

template
bool pztopology::TPZTriangle::MapToSide<REAL>(int side, TPZVec<REAL> &InternalPar, TPZVec<REAL> &SidePar, TPZFMatrix<REAL> &JacToSide);

#ifdef _AUTODIFF
template
bool pztopology::TPZTriangle::MapToSide<Fad<REAL> >(int side, TPZVec<Fad<REAL> > &InternalPar, TPZVec<Fad<REAL> > &SidePar, TPZFMatrix<Fad<REAL> > &JacToSide);
#endif
