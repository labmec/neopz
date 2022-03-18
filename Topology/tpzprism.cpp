/**
 * @file
 * @brief Contains the implementation of the TPZPrism methods. 
 */

#include "tpzprism.h"

#include "pzmanvector.h"
#include "pznumeric.h"
#include "pzerror.h"
#include "pzreal.h"
#include "pzquad.h"
#include "pzeltype.h"
#include "tpzquadrilateral.h"
#include "tpztriangle.h"
#include "pzlog.h"

#include "fad.h"


#ifdef PZ_LOG
static TPZLogger logger("pz.topology.pzprism");
#endif

using namespace std;

namespace pztopology {

	static constexpr int FaceConnectLocId[5][9] = { {0,1,2,6,7,8,15,-1,-1},{0,1,4,3,6,10,12,9,16},
		{1,2,5,4,7,11,13,10,17},{0,2,5,3,8,11,14,9,18},{3,4,5,12,13,14,19,-1,-1} };

	static constexpr int sidedimension[21] = {0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,2,2,2,2,2,3};

    static constexpr int fSideOrient[5] = {-1,1,1,-1,1};
	
	static constexpr int nhighdimsides[21] = {7,7,7,7,7,7,3,3,3,3,3,3,3,3,3,1,1,1,1,1,0};
	
	static constexpr int highsides[21][7] = {
		{6,8,9,15,16,18,20},
		{6,7,10,15,16,17,20},
		{7,8,11,15,17,18,20},
		{9,12,14,16,18,19,20},
		{10,12,13,16,17,19,20},
		{11,13,14,17,18,19,20},
		{15,16,20},
		{15,17,20},
		{15,18,20},
		{16,18,20},
		{16,17,20},
		{17,18,20},
		{16,19,20},
		{17,19,20},
		{18,19,20},
		{20},
		{20},
		{20},
		{20},
		{20},
		{-999}
	};
	
	static constexpr REAL sidetosidetransforms[21][7][4][3] = {
		{
			//0
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,0,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,0,-1}}
		},
		{
			//1
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-1,-99}},
			//{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,0,-99}},//est� errado deve ser {-1,-1,-99}
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-1,-99}},//CEDRIC
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-1}}
		},
		{
			//2
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-1}}
		},
		{
			//3
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,0,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,0,1}}
		},
		{
			//4
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,0,1}}
		},
		{
			//5
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{0,1,1}}
		},
		{
			//6
			{{0.5,0,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0,-99}},
			{{1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,-1,-99}},
			{{0.5,0,0},{-99,-99,-99},{-99,-99,-99},{0.5,0,-1}}
		},
		{
			//7
			{{-0.5,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0.5,-99}},
			{{1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,-1,-99}},
			{{-0.5,0.5,0},{-99,-99,-99},{-99,-99,-99},{0.5,0.5,-1}}
		},
		{
			//8
			{{0,-0.5,-99},{-99,-99,-99},{-99,-99,-99},{0,0.5,-99}},
			{{-1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,-1,-99}},
			{{0,-0.5,0},{-99,-99,-99},{-99,-99,-99},{0,0.5,-1}}
		},
		{
			//9
			{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{-1,0,-99}},
			{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{-1,0,-99}},
			{{0,0,1},{-99,-99,-99},{-99,-99,-99},{0,0,0}}
		},
		{
			//10
			{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
			{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{-1,0,-99}},
			{{0,0,1},{-99,-99,-99},{-99,-99,-99},{1,0,0}}
		},
		{
			//11
			{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
			{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
			{{0,0,1},{-99,-99,-99},{-99,-99,-99},{0,1,0}}
		},
		{
			//12
			{{1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
			{{0.5,0,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0,-99}},
			{{0.5,0,0},{-99,-99,-99},{-99,-99,-99},{0.5,0,1}}
		},
		{
			//13
			{{1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
			{{-0.5,0.5,-99},{-99,-99,-99},{-99,-99,-99},{0.5,0.5,-99}},
			{{-0.5,0.5,0},{-99,-99,-99},{-99,-99,-99},{0.5,0.5,1}}
		},
		{
			//14
			{{-1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
			{{0,-0.5,-99},{-99,-99,-99},{-99,-99,-99},{0,0.5,-99}},
			{{0,-0.5,0},{-99,-99,-99},{-99,-99,-99},{0,0.5,1}}
		},
		{
			//15
			{{1,0,0},{0,1,0},{-99,-99,-99},{0,0,-1}}
		},
		{
			//16
			{{0.5,0,0},{0,0,1},{-99,-99,-99},{0.5,0,0}}
		},
		{
			//17
			{{-0.5,0.5,0},{0,0,1},{-99,-99,-99},{0.5,0.5,0}}
		},
		{
			//18
			{{0,0.5,0},{0,0,1},{-99,-99,-99},{0,0.5,0}}
		},
		{
			//19
			{{1,0,0},{0,1,0},{-99,-99,-99},{0,0,1}}
		},
		{
			//20
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-99,-99,-99}}
		}
	};
	
	static constexpr REAL MidSideNode[21][3] = {
		/*00*/{0.,.0,-1.},/*01*/{1.,0.,-1.},/*02*/{.0,1.,-1.},/*03*/{.0,.0, 1.},
		/*04*/{1.,.0, 1.},/*05*/{0.,1., 1.},/*06*/{.5,.0,-1.},/*07*/{.5,.5,-1.},
		/*08*/{.0,.5,-1.},/*09*/{0.,.0, 0.},/*10*/{1.,.0, 0.},/*11*/{.0,1., 0.},
		/*12*/{.5,.0, 1.},/*13*/{.5,.5, 1.},/*14*/{.0,.5, 1.},/*15*/{1./3.,1./3.,-1.},
		/*16*/{.5,.0, 0.},/*17*/{.5,.5, 0.},/*18*/{0.,.5, 0.},/*19*/{1./3.,1./3., 1.},
		/*20*/{1./3.,1./3.,0.} };
    
    static constexpr REAL bPrism[63][3] =
    {
        {0,0,-1}, {0,0,-1}, {0,0,-1}, {0,0,-1}, {0,0,-1}, {0,0,-1}, {0,0,-1},   // face 0
        {0,-1,0}, {1,-1,0}, {1,-1,0}, {0,-1,0}, {0,-1,0},  {1,-1,0}, {0,-1,0}, {0,-1,0}, {0,-1,0},// face 1
        {1,0,0}, {0,1,0}, {0,1,0}, {1,0,0}, {1,1,0}, {0,1,0}, {1,1,0}, {1,0,0}, {1,1,0},// face 2
        {-1,0,0}, {-1,1,0}, {-1,1,0}, {-1,0,0}, {-1,0,0}, {-1,1,0}, {-1,0,0}, {-1,0,0}, {-1,0,0},// face 3
        {0,0,1}, {0,0,1}, {0,0,1}, {0,0,1}, {0,0,1}, {0,0,1}, {0,0,1}, // face 4
        //internos
        // arestas
        {1,0,0},{-1,1,0},{0,-1,0},   {0,0,1},{0,0,1},{0,0,1},   {1,0,0},{-1,1,0},{0,-1,0},
        //faces
        {-1,0,0}, {0,1,0}, // tang da face 0	229	        {-1,0,0}, {0,1,0}, // tang da face 0
        {1,0,0 }, {0,0,1}, //face 1	230	        {1,0,0 }, {0,0,1}, //face 1
        {-1,1,0}, {0,0,1}, //face 2	231	        {-1,1,0}, {0,0,1}, //face 2
        {0,0,1},  {0,1,0}, //face 3	232	        {0,0,1},  {0,1,0}, //face 3
        {1,0,0},  {0,1,0}, //face 4
                //interior
        {1,0,0} ,
        {0,1,0} ,
        {0,0,1}
    };
    static constexpr REAL t1Prism[63][3] =
    {
        {1,0,0},{1,0,0},{1,0,0},{1,0,0},{1,0,0},{1,0,0},{1,0,0},  // face 0
        {1,0,0}, {1,0,0},{1,0,0},{1,0,0}, {1,0,0}, {1,0,0}, {1,0,0}, {1,0,0},{1,0,0},// face 1
        {0,0,-1}, {0,0,-1}, {0,0,-1},{0,0,-1},{0,0,-1},{0,0,-1}, {0,0,-1},{0,0,-1},{0,0,-1}, // face 2
        {0,0,1}, {0,0,1}, {0,0,1}, {0,0,1}, {0,0,1}, {0,0,1}, {0,0,1},{0,0,1},{0,0,1}, // fsce 3
        {1,0,0}, {1,0,0}, {1,0,0}, {1,0,0}, {1,0,0},{1,0,0},{1,0,0}, // face 4
        //internos
        // arestas
        {0,1,0},{-1,-1,0},{1,0,0},   {-1,0,0},{0,-1,0},{1,1,0},  {0,1,0},{-1,-1,0},{1,0,0},
        //faces
        {0,1,0}, {0,0,1},  //face 0
        {0,1,0}, {-1,1,0}, //face 1
        {0,0,1}, {1,-1,0}, //face 2
        {0,1,0}, {0,0,-1}, //face 3
        {0,1,0}, {-1,0,0}, //face 4
        
        //interior
        {0,1,0} ,
        {0,0,1} ,
        {1,0,0}
    };
    
    static constexpr REAL t2Prism[63][3] =
    {
        {0,-1,0}, {0,-1,0}, {0,-1,0}, {0,-1,0}, {0,-1,0}, {0,-1,0}, {0,-1,0},// face 0
        {0,0,1}, {0,0,1}, {0,0,1}, {0,0,1}, {0,0,1}, {0,0,1}, {0,0,1},{0,0,1},{0,0,1}, // face 1
        {-1,1,0}, {-1,1,0}, {-1,1,0}, {-1,1,0}, {-1,1,0}, {-1,1,0}, {-1,1,0},{-1,1,0},{-1,1,0}, // face 2
        {0,1,0}, {0,1,0}, {0,1,0}, {0,1,0}, {0,1,0},{0,1,0}, {0,1,0},{0,1,0},{0,1,0}, // fsce 3
        {0,1,0}, {0,1,0}, {0,1,0}, {0,1,0}, {0,1,0},{0,1,0},{0,1,0}, // face 4
        //internos
        // arestas
        {0,0,1},{0,0,1},{0,0,1},  {0,-1,0},{1,1,0},{-1,0,0},  {0,0,1},{0,0,1},{0,0,1},
        //faces
        {0,0,1}, {1,0,0}, //face 0
        {0,0,1}, {0,1,0}, //face 1
        {1,1,0},  {1,1,0}, //face 2
        {-1,0,0}, {-1,0,0}, //face 3
        {0,0,1},  {0,0,1}, //face 4
        //interior
        {0,0,1} ,
        {1,0,0} ,
        {0,1,0}
    };
    
    static constexpr int vectorsideorderPr[63] =
    {
        0,1,2,6,7,8,15,
        0,1,4,3,6,10,12,9,16,
        1,2,5,4,7,11,13,10,17,
        0,2,5,3,8,11,14,9,18,
        3,4,5,12,13,14,19,
        6,7,8,9,10,11,12,13,14,
        15,15,
        16,16,
        17,17,
        18,18,
        19,19,
        20,20,20
    };
	
//    static int bilinearounao [63] =   {
//        0,0,0,0,0,0,0,0,0,0,
//        0,0,0,0,0,0,0,0,0,0,
//        0,0,0,0,0,0,0,0,0,0,
//        0,0,0,0,0,0,0,0,0,0,
//        0,0,0,0,1,1,1,0,0,0,
//        0,0,1,1,1,1,1,1,0,0,
//        1,1,1//1,1,1
//    };
    static constexpr int bilinearounao [63] =   {
        0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,
        0,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,
        0,0,1
    };
    
    static constexpr int direcaoksioueta [63] = {
        0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,
        0,1,0,1,0,1,0,1,0,1,
        0,1,2};

    template<class T>
    inline void TPZPrism::TShape(const TPZVec<T> &loc,TPZFMatrix<T> &phi,TPZFMatrix<T> &dphi) {
        T qsi = loc[0], eta = loc[1] , zeta  = loc[2];

        phi(0,0)  = .5*(1.-qsi-eta)*(1.-zeta);
        phi(1,0)  = .5*qsi*(1.-zeta);
        phi(2,0)  = .5*eta*(1.-zeta);
        phi(3,0)  = .5*(1.-qsi-eta)*(1.+zeta);
        phi(4,0)  = .5*qsi*(1.+zeta);
        phi(5,0)  = .5*eta*(1.+zeta);

        dphi(0,0) = -.5*(1.-zeta);
        dphi(1,0) = -.5*(1.-zeta);
        dphi(2,0) = -.5*(1.-qsi-eta);
        dphi(0,1) =  .5*(1.-zeta);
        dphi(1,1) =  .0;
        dphi(2,1) = -.5*qsi;
        dphi(0,2) =  .0;
        dphi(1,2) =  .5*(1.-zeta);
        dphi(2,2) = -.5*eta;
        dphi(0,3) = -.5*(1.+zeta);
        dphi(1,3) = -.5*(1.+zeta);
        dphi(2,3) =  .5*(1.-qsi-eta);
        dphi(0,4) =  .5*(1.+zeta);
        dphi(1,4) =  .0;
        dphi(2,4) =  .5*qsi;
        dphi(0,5) =  .0;
        dphi(1,5) =  .5*(1.+zeta);
        dphi(2,5) =  .5*eta;

    }

    template<class T>
    void TPZPrism::BlendFactorForSide(const int &side, const TPZVec<T> &xiVec, T &blendFactor,
                                        TPZVec<T> &blendFactorDxi){
        blendFactorDxi.Resize(TPZPrism::Dimension, (T) 0);
        blendFactor = 0;
        const REAL tol = pztopology::GetTolerance();
        #ifdef PZDEBUG
        std::ostringstream sout;
        if(side < NCornerNodes || side >= NSides){
            sout<<"The side\t"<<side<<"is invalid. Aborting..."<<std::endl;
        }

        if(!pztopology::TPZPrism::IsInParametricDomain(xiVec,tol)){
            sout<<"The method BlendFactorForSide expects the point xi to correspond to coordinates of a point";
            sout<<" inside the parametric domain. Aborting...";
        }

        if(!sout.str().empty()){
            PZError<<std::endl<<sout.str()<<std::endl;
#ifdef PZ_LOG
            LOGPZ_FATAL(logger,sout.str().c_str());
#endif
            DebugStop();
        }
        #endif
        //if the point is singular, the blend factor and its derivatives should be zero
        if(!CheckProjectionForSingularity(side,xiVec)){
            std::cout<<"Side projection is not regular and it should have been checked earlier. Aborting.."<<std::endl;
            DebugStop();
            blendFactor = 0;
            for(int i = 0; i < blendFactorDxi.size(); i++) blendFactorDxi[i] = 0;
            return;
        }
        blendFactorDxi.Resize(TPZPrism::Dimension, (T) 0);
        const T &xi = xiVec[0];
        const T &eta = xiVec[1];
        const T &zeta = xiVec[2];
        switch(side){
            case  0:
            case  1:
            case  2:
            case  3:
            case  4:
            case  5:
                blendFactor = 0;
                return;
            case  6:
                blendFactor = 0.5*(1.-zeta)*(1.-eta)*(1.-eta);
                blendFactorDxi[0] = 0;
                blendFactorDxi[1] = -1.*(1 - eta)*(1 - zeta);
                blendFactorDxi[2] = -0.5*((1 - eta)*(1 - eta));
                return;
            case  7:
                blendFactor = 0.5*(1.-zeta)*(xi+eta)*(xi+eta);
                blendFactorDxi[0] =1.*(eta + xi)*(1 - zeta);
                blendFactorDxi[1] =1.*(eta + xi)*(1 - zeta);
                blendFactorDxi[2] =-0.5*((eta + xi)*(eta + xi));
                return;
            case  8:
                blendFactor = 0.5*(1.-zeta)*(1.-xi)*(1.-xi);
                blendFactorDxi[0] =-1.*(1 - xi)*(1 - zeta);
                blendFactorDxi[1] =0;
                blendFactorDxi[2] =-0.5*((1 - xi)*(1 - xi));
                return;
            case  9:
                blendFactor = 1. - xi - eta;
                blendFactorDxi[0] = -1;
                blendFactorDxi[1] = -1;
                blendFactorDxi[2] = 0;
                return;
            case 10:
                blendFactor = xi;
                blendFactorDxi[0] = 1;
                blendFactorDxi[1] = 0;
                blendFactorDxi[2] = 0;
                return;
            case 11:
                blendFactor = eta;
                blendFactorDxi[0] = 0;
                blendFactorDxi[1] = 1;
                blendFactorDxi[2] = 0;
                return;
            case 12:
                blendFactor = 0.5*(1.+zeta)*(1.-eta)*(1.-eta);
                blendFactorDxi[0] = 0;
                blendFactorDxi[1] = -1.*(1 - eta)*(1 + zeta);
                blendFactorDxi[2] = 0.5*((1 - eta)*(1 - eta));
                return;
            case 13:
                blendFactor = 0.5*(1.+zeta)*(xi+eta)*(xi+eta);
                blendFactorDxi[0] = 1.*(eta + xi)*(1 + zeta);
                blendFactorDxi[1] = 1.*(eta + xi)*(1 + zeta);
                blendFactorDxi[2] = 0.5*((eta + xi)*(eta + xi));
                return;
            case 14:
                blendFactor = 0.5*(1.+zeta)*(1.-xi)*(1.-xi);
                blendFactorDxi[0] = -1.*(1 - xi)*(1 + zeta);
                blendFactorDxi[1] = 0;
                blendFactorDxi[2] = 0.5*((1 - xi)*(1 - xi));
                return;
            case 15:
                blendFactor = 0.5*(1.-zeta);
                blendFactorDxi[0] = 0;
                blendFactorDxi[1] = 0;
                blendFactorDxi[2] = -0.5;
                return;
            case 16:
                blendFactor = 1.-eta;
                blendFactorDxi[0] = 0;
                blendFactorDxi[1] = -1;
                blendFactorDxi[2] = 0;
                return;
            case 17:
                blendFactor = xi+eta;
                blendFactorDxi[0] = 1;
                blendFactorDxi[1] = 1;
                blendFactorDxi[2] = 0;
                return;
            case 18:
                blendFactor = 1.-xi;
                blendFactorDxi[0] = -1;
                blendFactorDxi[1] = 0;
                blendFactorDxi[2] = 0;
                return;
            case 19:
                blendFactor = 0.5*(1.+zeta);
                blendFactorDxi[0] = 0;
                blendFactorDxi[1] = 0;
                blendFactorDxi[2] = 0.5;
                return;
        }
    }

    int TPZPrism:: NBilinearSides()
    {
        DebugStop();
        return 9;
    }
    
	void TPZPrism::LowerDimensionSides(int side,TPZStack<int> &smallsides)
	{
		smallsides.Resize(0);
		int nsidecon = NContainedSides(side);
		int is;
		for(is=0; is<nsidecon-1; is++)
			smallsides.Push(ContainedSideLocId(side,is));
	}
	
	void TPZPrism::LowerDimensionSides(int side,TPZStack<int> &smallsides, int DimTarget)
	{
		smallsides.Resize(0);
		int nsidecon = NContainedSides(side);
		for(int is = 0; is < nsidecon - 1; is++) {
			if (SideDimension(ContainedSideLocId(side,is)) == DimTarget) smallsides.Push(ContainedSideLocId(side,is));
		}
	}
	
	void TPZPrism::HigherDimensionSides(int side, TPZStack<int> &high)
	{
		if(side <0 || side >= NSides) {
			PZError << "TPZPrism::HigherDimensionSides side "<< side << endl;
		}
		int is;
		for(is=0; is<nhighdimsides[side]; is++) high.Push(highsides[side][is]);
		
	}
	
	static constexpr int nsidenodes[21] = {
		1,1,1,1,1,1,
		2,2,2,2,2,2,2,2,2,
		3,4,4,4,3,
		6};
	
	int TPZPrism::NSideNodes(int side)
	{
		return nsidenodes[side];
	}
	//Tentando criar o metodo
	int TPZPrism::NumSides(int dimension) {
		if(dimension<0 || dimension> 3) {
			PZError << "TPZPrism::NumSides. Bad parameter i.\n";
			return 0;
		}
		if(dimension==0) return 6;
		if(dimension==1) return 9;
		if(dimension==2) return 5;
		if(dimension==3) return 1;
		return -1;
	}
	
	int TPZPrism::SideNodeLocId(int side, int node)
	{
		if(side<6 && node == 0) return side;
		if(side>= 6 && side < 15 && node<2) return SideNodes[side-6][node];
		if(side==15) {
			if(node < 3) return FaceNodes[side-15][node];
			else if(node == 3) return -1; //previsto para faces triangulares
		}
		if(side>15 && side <19 && node <4) return FaceNodes[side-15][node];
		if(side==19) {
			if(node<3) return FaceNodes[side-15][node];
			else if(node == 3) return -1; // Previsto p/ faces triangulares
        }
		
		if(side==20 && node<6) return node;
		PZError << "TPZPrism::SideNodeLocId inconsistent side or node " << side << ' ' << node << endl;
		return -1;
	}
	
	void TPZPrism::CenterPoint(int side, TPZVec<REAL> &center) {
        if (center.size()!=Dimension) {
            DebugStop();
        }
		int i;
		for(i=0; i<Dimension; i++) {
			center[i] = MidSideNode[side][i];
		}
	}
	
	int TPZPrism::SideDimension(int side) {
		if(side<0 || side >= NSides) {
			PZError << "TPZPrism::SideDimension side " << side << endl;
			return -1;
		}
		return sidedimension[side];
	}
	
	TPZTransform<> TPZPrism::SideToSideTransform(int sidefrom, int sideto)
	{
		if(sidefrom <0 || sidefrom >= NSides || sideto <0 || sideto >= NSides) {
			PZError << "TPZPrism::HigherDimensionSides sidefrom "<< sidefrom << 
			' ' << sideto << endl;
			return TPZTransform<>(0);
		}
		if(sidefrom == sideto) {
			return TPZTransform<>(sidedimension[sidefrom]);
		}
		if(sidefrom == NSides-1) {
			return TransformElementToSide(sideto);
		}
        if (sideto==NSides-1) {
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
		PZError << "TPZPrism::SideToSideTransform highside not found sidefrom "
		<< sidefrom << ' ' << sideto << endl;
		return TPZTransform<>(0);
	}
	
	
	TPZTransform<> TPZPrism::TransformElementToSide(int side){
		
		if(side<0 || side>20){
			PZError << "TPZPrism::TransformElementToSide called with side error\n";
			return TPZTransform<>(0,0);
		}
		
        if(side<0 || side>20){
            PZError << "TPZPrism::TransformElementToSide called with side error\n";
            return TPZTransform<>(0,0);
        }
        
        TPZTransform<> t(sidedimension[side],3);
        t.Mult().Zero();
        t.Sum().Zero();
   
        switch(side){
            case 0:
            case 1:
            case 2:
            case 3:
            case 4:
            case 5:
                return t;

            case  6:
            case 12:
                t.Mult()(0,0) =  2.0;
                t.Mult()(0,1) =  1.0;

                t.Sum()(0,0)  = -1.0;
                return t;
            case  7:
            case 13:
                t.Mult()(0,0) = -1.0;
                t.Mult()(0,1) =  1.0;
                return t;

            case  8:
            case 14:
                t.Mult()(0,0) = -1.0;
                t.Mult()(0,1) = -2.0;
                t.Sum()(0,0)  =  1.0;
                return t;

            case 9:
            case 10:
            case 11:
                t.Mult()(0,2) = 1.0;
                return t;

            case  15:

                t.Mult()(0,0) =  1.0;
                t.Mult()(1,1) =  1.0;

//                t.Sum()(0,0) = -1.0;
//                t.Sum()(1,0) = -1.0;
                return t;

            case 16:

                t.Mult()(0,0) =  2.0;
                t.Mult()(1,2) =  1.0;

                t.Sum()(0,0) = -1;
                return t;

            case 17:

                t.Mult()(0,0) =  -1.0;
                t.Mult()(0,1) =  1.0;
                t.Mult()(1,2) =  1.0;

                return t;

            case 18:

                t.Mult()(0,1) =  2.0;
                t.Mult()(1,2) =  1.0;

                t.Sum()(0,0)  = -1.0;
                return t;

            case 19:
                
                
                t.Mult()(0,0) =  1.0;
                t.Mult()(1,1) =  1.0;
                
//                t.Sum()(0,0) = -1.0;
//                t.Sum()(1,0) = -1.0;
                

                return t;
            case 20:
                t.Mult()(0,0) = 1.0;
                t.Mult()(1,1) = 1.0;
                t.Mult()(2,2) = 1.0; //rev
                return t;
        }
        return TPZTransform<>(0,0);
	}
	
	TPZTransform<> TPZPrism::TransformSideToElement(int side){
		
        if(side<0 || side>20){
            PZError << "TPZPrism::TransformSideToElement side out range\n";
            return TPZTransform<>(0,0);
        }
        TPZTransform<> t(3,sidedimension[side]);
        t.Mult().Zero();
        t.Sum().Zero();
        
        switch(side){
            case 0:
                t.Sum()(0,0) =  0.0;
                t.Sum()(1,0) =  0.0;
                t.Sum()(2,0) = -1.0;
                return t;
            case 1:
                t.Sum()(0,0) =  1.0;
                t.Sum()(1,0) =  0.0;
                t.Sum()(2,0) = -1.0;
                return t;
            case 2:
                t.Sum()(0,0) =  0.0;
                t.Sum()(1,0) =  1.0;
                t.Sum()(2,0) = -1.0;
                return t;
            case 3:
                t.Sum()(0,0) =  0.0;
                t.Sum()(1,0) =  0.0;
                t.Sum()(2,0) =  1.0;
                return t;
            case 4:
                t.Sum()(0,0) =  1.0;
                t.Sum()(1,0) =  0.0;
                t.Sum()(2,0) =  1.0;
                return t;
            case 5:
                t.Sum()(0,0) =  0.0;
                t.Sum()(1,0) =  1.0;
                t.Sum()(2,0) =  1.0;
                return t;
            case 6:
                t.Mult()(0,0) =  0.5;
                t.Sum() (0,0) =  0.5;
                t.Sum() (2,0) = -1.0;
                return t;
            case 7:
                t.Mult()(0,0) = -0.5;
                t.Mult()(1,0) =  0.5;
                t.Sum() (0,0) =  0.5;
                t.Sum() (1,0) =  0.5;
                t.Sum() (2,0) = -1.0;
                return t;
            case 8:
                t.Mult()(1,0) = -0.5;
                t.Sum() (1,0) =  0.5;
                t.Sum() (2,0) = -1.0;
                return t;
            case 9:
                t.Mult()(2,0) =  1.0;
                return t;
            case 10:
                t.Mult()(2,0) =  1.0;
                t.Sum() (0,0) =  1.0;
                return t;
            case 11:
                t.Mult()(2,0) =  1.0;
                t.Sum() (1,0) =  1.0;
                return t;
            case 12:
                t.Mult()(0,0) =  0.5;
                t.Sum() (0,0) =  0.5;
                t.Sum() (2,0) =  1.0;
                return t;
            case 13:
                t.Mult()(0,0) = -0.5;
                t.Mult()(1,0) =  0.5;
                t.Sum() (0,0) =  0.5;
                t.Sum() (1,0) =  0.5;
                t.Sum() (2,0) =  1.0;
                return t;
            case 14:
                t.Mult()(1,0) = -0.5;
                t.Sum() (1,0) =  0.5;
                t.Sum() (2,0) =  1.0;
                return t;
            case 15:
                t.Mult()(0,0) =  1.0;
                t.Mult()(1,1) =  1.0;
//                t.Sum() (0,0) =  0.5;
//                t.Sum() (1,0) =  0.5;
                
                t.Sum() (2,0) = -1.0;
                return t;
            case 16:
                t.Mult()(0,0) =  0.5;
                t.Mult()(2,1) =  1.0;
                t.Sum() (0,0) =  0.5;
                return t;
            case 17:
                t.Mult()(0,0) = -0.5;
                t.Mult()(1,0) =  0.5;
                t.Mult()(2,1) =  1.0;
                t.Sum() (0,0) =  0.5;
                t.Sum() (1,0) =  0.5;
                return t;
            case 18:
                t.Mult()(1,0) =  0.5;
                t.Mult()(2,1) =  1.0;
                t.Sum() (1,0) =  0.5;
                return t;
            case 19:
                t.Mult()(0,0) =  1.0;
                t.Mult()(1,1) =  1.0;
//                t.Sum() (0,0) =  0.5;
//                t.Sum() (1,0) =  0.5;
                
                t.Sum() (2,0) =  1.0;
                return t;
            case 20:
                t.Mult()(0,0) =  1.0;
                t.Mult()(1,1) =  1.0;
                t.Mult()(2,2) =  1.0;
                return t;
        }
		return TPZTransform<>(0,0);
	}
	
	
	TPZIntPoints * TPZPrism::CreateSideIntegrationRule(int side, int order){
		if(side<0 || side>20) {
			PZError << "TPZPrism::CreateSideIntegrationRule. bad side number.\n";
			return 0;
		}
		if(side<6)   return new TPZInt1Point(order);                   // sides 0 to 7 are vertices
		if(side<15)  return new TPZInt1d(order);                  // sides 7 to 14 are lines
		if(side==15 || side==19) return new TPZIntTriang(order);  // sides 15 and 19 are triangles
		if(side<20)  {
			return new TPZIntQuad(order,order);                   // sides 16 to 18 are quadrilaterals
		}
		if(side==20) {
			return new IntruleType(order,order);                  // integration of the element
		}
		return 0;
	}
	
	
	MElementType TPZPrism::Type(int side)
	{
		switch(side) {
			case 0:
			case 1:
			case 2:
			case 3:
			case 4:
			case 5:
				return EPoint;
			case 6:
			case 7:
			case 8:
			case 9:
			case 10:
			case 11:
			case 12:
			case 13:
			case 14:
				return EOned;
			case 15:
				return ETriangle;
			case 16:
			case 17:
			case 18:
				return EQuadrilateral;
			case 19:
				return ETriangle;
			case 20:
				return EPrisma;
			default:
				return ENoType;
		}
	}
	
	int TPZPrism::NumSides() {
		return NSides;
	}
	
	
	int TPZPrism::NContainedSides(int side) {
		if(side<0)   return -1;
		if(side<6)   return 1;//cantos : 0 a 5
		if(side<15)  return 3;//arestas
		if(side==15 || side==19)  return 7;//faces : 15,19 , triangulares
		if(side<19) return 9;//faces : 16 a 18  quadrilaterais
		if(side==20) return 21;//centro : 20
		return -1;
	}
	
	int TPZPrism::ContainedSideLocId(int side, int node) {
		if(side<0 || side>20 || node < 0) return -1;
		if(side<6) {
			if(node==0) return side;
		} else
			if(side<9) {//6,7,8
				int s = side-6;//0,1,2
				if(!node) return s;//0,1,2
				if(node==1) return (s+1)%3;//1,2,0
				if(node==2) return side;//6,7,8
			} else
				if(side<12) {//9,10,11
					int s = side-9;//0,1,2
					if(!node) return s;//0,1,2,4
					if(node==1) return s+3;//3,4,5
					if(node==2) return side;//5,6,7
				} else
					if(side<15) {//12,13,14
						int s = side-9;//3,4,5
						if(!node) return s;//3,4,5
						if(node==1) return (s+1)%3+3;//4,5,3
						if(node==2) return side;//12,13,14
					} else
						if(side==15 || side==19) {
							int s = side-15;
							if(side==15 && node<7) return FaceConnectLocId[s][node];
							if(side==19 && node<7) return FaceConnectLocId[s][node];
						} else
							if(side<20) {//16,17,18
								int s = side-15;
								if(node<9) return FaceConnectLocId[s][node];
							} else
								if(side==20 && node<21){
									return node;
								}
		PZError << "TPZShapePrism::ContainedSideLocId called for node = "
		<< node << " and side = " << side << "\n";
		return -1;
	}
	
	bool TPZPrism::IsInParametricDomain(const TPZVec<REAL> &pt, REAL tol){
		const REAL qsi = pt[0];
		const REAL eta = pt[1];
		const REAL zeta = pt[2];
		if( ( qsi <= 1. + tol ) && ( qsi >= 0. - tol ) &&
		   ( eta <= 1. + tol ) && ( eta >= 0. - tol ) &&
		   ( eta <= 1. - qsi + tol ) && 
		   ( zeta <= 1. + tol ) && ( zeta >= -1. - tol) ){
			return true;
		}
		else{
			return false;
		}  
		
	}//method
    
    /** @brief Generates a random point in the master domain */
    void TPZPrism::RandomPoint(TPZVec<REAL> &pt)
    {
        REAL val = (REAL) rand() / (RAND_MAX);
        pt[0] = val;
        val = (1.-pt[0]) * (REAL) rand() / (RAND_MAX);
        pt[1] = val;
        val = -1. + 2. * (REAL) rand() / (RAND_MAX);
        pt[2] = val;
    }

    template<class T>
    bool TPZPrism::CheckProjectionForSingularity(const int &side, const TPZVec<T> &xiInterior) {
        double zero = pztopology::gTolerance;

        T qsi = xiInterior[0];
        T eta = xiInterior[1];
        T zeta = xiInterior[2];

        bool regularmap = true;
        switch(side)
        {
            case 0:
            case 1:
            case 2:
            case 3:
            case 4:
            case 5:
                regularmap = true;
                break;
            case 6://1D
                if(fabs((T)(eta-1.)) < zero) regularmap = false;
                break;
            case 7://1D
                if(fabs((T)(qsi+eta)) < zero) regularmap = false;
                break;
            case 8://1D
                if(fabs((T)(qsi-1.)) < zero) regularmap = false;
                break;
            case 9://1D
            case 10://1D
            case 11://1D
                regularmap = true;
                break;
            case 12://1D
                if(fabs((T)(eta-1.)) < zero) regularmap = false;
                break;
            case 13://1D
                if(fabs((T)(qsi+eta)) < zero) regularmap = false;
                break;
            case 14://1D
                if(fabs((T)(qsi-1.)) < zero) regularmap = false;
                break;
            case 15://2D - triangle
                regularmap = true;
                break;
            case 16://2D - quadrilateral
                if(fabs((T)(eta-1.)) < zero) regularmap = false;
                break;
            case 17://2D - quadrilateral
                if(fabs((T)(qsi+eta)) < zero) regularmap = false;
                break;
            case 18://2D - quadrilateral
                if(fabs((T)(qsi-1.)) < zero) regularmap = false;
                break;
            case 19:
            case 20:
                regularmap = true;
        }
        if(side > 20)
        {
            cout << "Cant compute CheckProjectionForSingularity method in TPZPrism class!\nParameter (SIDE) must be between 6 and 19!\nMethod Aborted!\n";
            DebugStop();
        }
        return regularmap;
    }

    template<class T>
    void TPZPrism::MapToSide(int side, TPZVec<T> &InternalPar, TPZVec<T> &SidePar, TPZFMatrix<T> &JacToSide) {
        double zero = pztopology::GetTolerance();

        T qsi = InternalPar[0];
        T eta = InternalPar[1];
        T zeta = InternalPar[2];

        if(!CheckProjectionForSingularity(side,InternalPar)){
            std::cout<<"Side projection is not regular and it should have been checked earlier. Aborting.."<<std::endl;
            DebugStop();
        }
        switch(side)
        {
            case 0:
            case 1:
            case 2:
            case 3:
            case 4:
            case 5:
            {
                SidePar.Resize(0); JacToSide.Resize(0,0);
                break;
            }
            case 6://1D
                SidePar.Resize(1); JacToSide.Resize(1,3);
                {
                    T den = -1 + eta;
                    SidePar[0] = -1 - (2*qsi)/den;
                    JacToSide(0,0) = -2/den;
                    JacToSide(0,1) = (2*qsi)/(den*den);
                    JacToSide(0,2) = 0;
                }
                break;

            case 7://1D
                SidePar.Resize(1); JacToSide.Resize(1,3);
                {
                    T den = eta + qsi;
                    SidePar[0] = -1 + (2*eta)/den;
                    JacToSide(0,0) = (-2*eta)/(den*den);
                    JacToSide(0,1) = (2*qsi)/(den*den);
                    JacToSide(0,2) = 0;
                }
                break;

            case 8://1D
                SidePar.Resize(1); JacToSide.Resize(1,3);
                {
                    T den = -1 + qsi;
                    SidePar[0] = 1 + (2*eta)/den;
                    JacToSide(0,0) = (-2*eta)/(den*den);
                    JacToSide(0,1) = 2/den;
                    JacToSide(0,2) = 0;
                }
                break;

            case 9://1D
                SidePar.Resize(1); JacToSide.Resize(1,3);
                SidePar[0] = zeta;
                JacToSide(0,0) = 0;
                JacToSide(0,1) = 0;
                JacToSide(0,2) = 1;
                break;

            case 10://1D
                SidePar.Resize(1); JacToSide.Resize(1,3);
                SidePar[0] = zeta;
                JacToSide(0,0) = 0;
                JacToSide(0,1) = 0;
                JacToSide(0,2) = 1;
                break;

            case 11://1D
                SidePar.Resize(1); JacToSide.Resize(1,3);
                SidePar[0] = zeta;
                JacToSide(0,0) = 0;
                JacToSide(0,1) = 0;
                JacToSide(0,2) = 1;
                break;

            case 12://1D
                SidePar.Resize(1); JacToSide.Resize(1,3);
                {
                    T den = 1 - eta;
                    SidePar[0] = -1 +(2*qsi)/den;
                    JacToSide(0,0) = 2/den;
                    JacToSide(0,1) = (2*qsi)/(den*den);
                    JacToSide(0,2) = 0;
                }
                break;

            case 13://1D
                SidePar.Resize(1); JacToSide.Resize(1,3);
                {
                    SidePar[0] = -1 + (2*eta)/(eta + qsi);
                    JacToSide(0,0) = (-2*eta)/((eta + qsi)*(eta + qsi));
                    JacToSide(0,1) = (2*qsi)/((eta + qsi)*(eta + qsi));
                    JacToSide(0,2) = 0;
                }
                break;

            case 14://1D
                SidePar.Resize(1); JacToSide.Resize(1,3);
                {
                    SidePar[0] = 1 + (2*eta)/(-1 + qsi);
                    JacToSide(0,0) = (-2*eta)/((-1 + qsi)*(-1 + qsi));
                    JacToSide(0,1) = 2/(-1 + qsi);
                    JacToSide(0,2) = 0;
                }
                break;

            case 15://2D - triangle
                SidePar.Resize(2); JacToSide.Resize(2,3);
                SidePar[0] = qsi;
                SidePar[1] = eta;
                JacToSide(0,0) = 1;
                JacToSide(0,1) = 0;
                JacToSide(0,2) = 0;
                JacToSide(1,0) = 0;
                JacToSide(1,1) = 1;
                JacToSide(1,2) = 0;
                break;

            case 16://2D - quadrilateral
                SidePar.Resize(2); JacToSide.Resize(2,3);
                {
                    SidePar[0] = -1 - (2*qsi)/(-1 + eta);
                    SidePar[1] = zeta;
                    JacToSide(0,0) = -2/(-1 + eta);
                    JacToSide(0,1) = (2*qsi)/((-1 + eta)*(-1 + eta));
                    JacToSide(0,2) = 0;
                    JacToSide(1,0) = 0;
                    JacToSide(1,1) = 0;
                    JacToSide(1,2) = 1;
                }
                break;

            case 17://2D - quadrilateral
                SidePar.Resize(2); JacToSide.Resize(2,3);
                {
                    SidePar[0] = -1 + (2*eta)/(eta + qsi);
                    SidePar[1] = zeta;
                    JacToSide(0,0) = (-2*eta)/((eta + qsi)*(eta + qsi));
                    JacToSide(0,1) = (2*qsi)/((eta + qsi)*(eta + qsi));
                    JacToSide(0,2) = 0;
                    JacToSide(1,0) = 0;
                    JacToSide(1,1) = 0;
                    JacToSide(1,2) = 1;
                }
                break;

            case 18://2D - quadrilateral
                SidePar.Resize(2); JacToSide.Resize(2,3);
                {
                    SidePar[0] = -1 - (2*eta)/(-1 + qsi);
                    SidePar[1] = zeta;
                    JacToSide(0,0) = (2*eta)/((-1 + qsi)*(-1 + qsi));
                    JacToSide(0,1) = -2/(-1 + qsi);
                    JacToSide(0,2) = 0;
                    JacToSide(1,0) = 0;
                    JacToSide(1,1) = 0;
                    JacToSide(1,2) = 1;
                }
                break;

            case 19://2D - triangle
                SidePar.Resize(2); JacToSide.Resize(2,3);
                SidePar[0] = qsi;
                SidePar[1] = eta;
                JacToSide(0,0) = 1;
                JacToSide(0,1) = 0;
                JacToSide(0,2) = 0;
                JacToSide(1,0) = 0;
                JacToSide(1,1) = 1;
                JacToSide(1,2) = 0;
                break;
            case 20:
                SidePar = InternalPar;
                JacToSide.Resize(3, 3);
                JacToSide.Identity();
                break;
        }
        if(side > 20)
        {
            cout << "Cant compute MapToSide method in TPZPrism class!\nParameter (SIDE) must be between 6 and 19!\nMethod Aborted!\n";
            cout << "This should have been caught earlier in the execution, there is something wrong.\n";
            cout << "Check method TPZTetrahedron::CheckProjectionForSingularity<T>\n";
            DebugStop();
        }
    }
    
    void TPZPrism::ParametricDomainNodeCoord(int node, TPZVec<REAL> &nodeCoord)
    {
        if(node > NCornerNodes)
        {
            DebugStop();
        }
        nodeCoord.Resize(Dimension, 0.);
        switch (node) {
            case (0):
            {
                nodeCoord[0] =  0.;
                nodeCoord[1] =  0.;
                nodeCoord[2] = -1.;
                break;
            }
            case (1):
            {
                nodeCoord[0] =  1.;
                nodeCoord[1] =  0.;
                nodeCoord[2] = -1.;
                break;
            }
            case (2):
            {
                nodeCoord[0] =  0.;
                nodeCoord[1] =  1.;
                nodeCoord[2] = -1.;
                break;
            }
            case (3):
            {
                nodeCoord[0] = 0.;
                nodeCoord[1] = 0.;
                nodeCoord[2] = 1.;
                break;
            }
            case (4):
            {
                nodeCoord[0] = 1.;
                nodeCoord[1] = 0.;
                nodeCoord[2] = 1.;
                break;
            }
            case (5):
            {
                nodeCoord[0] = 0.;
                nodeCoord[1] = 1.;
                nodeCoord[2] = 1.;
                break;
            }
            default:
            {
                DebugStop();
                break;
            }
        }
    }
	
	/**
	 * Method which identifies the transformation based on the IDs
	 * of the corner nodes
	 * @param id indexes of the corner nodes
	 * @return index of the transformation of the point corresponding to the topology
	 */
	int TPZPrism::GetTransformId(const TPZVec<int64_t> &id)
	{
		LOGPZ_ERROR(logger,"Please implement me")
		return -1;
	}
	
	/**
	 * Method which identifies the transformation of a side based on the IDs
	 * of the corner nodes
	 * @param id indexes of the corner nodes
	 * @return index of the transformation of the point corresponding to the topology
	 */	
	int TPZPrism::GetTransformId(const int side, const TPZVec<int64_t> &id)
	{
        switch (side) {
            case 0:
            case 1:
            case 2:
            case 3:
            case 4:
            case 5:
            
                return 0;
                break;
            case 6:
            case 7:
            case 8:
            case 9:
            case 10:
            case 11:
            case 12:
            case 13:
            case 14:
            {
                int in1 = ContainedSideLocId(side,0);
                int in2 = ContainedSideLocId(side,1);
                return id[in1]<id[in2] ? 0 : 1;
            }
                break;
                
            case 15:
            case 19:
            {
                TPZManVector<int64_t,3> locid(3);
                int i;
                for(i=0; i<3; i++) locid[i] = id[ContainedSideLocId(side,i)];
                return pztopology::TPZTriangle::GetTransformId(locid);
            }
                break;
               
            case 16:
            case 17:
            case 18:
            {
                TPZManVector<int64_t,4> locid(4);
                int i;
                for(i=0; i<4; i++) locid[i] = id[ContainedSideLocId(side,i)];
                return pztopology::TPZQuadrilateral::GetTransformId(locid);
               
            }
                break;
                
            case 20:
                return 0;
            default:
                break;
        }
        return -1;
        
//        LOGPZ_ERROR(logger,"Please implement me")
//        return -1;
	}
	
	/**
	 * Identifies the permutation of the nodes needed to make neighbouring elements compatible 
	 * in terms of order of shape functions
	 * @param side : side for which the permutation is needed
	 * @param id : ids of the corner nodes of the elements
	 * @param permgather : permutation vector in a gather order
	 */
	void TPZPrism::GetSideHDivPermutation(int transformationid, TPZVec<int> &permgather)
	{
	// Not complete
        DebugStop();
		LOGPZ_ERROR(logger,"Please implement me")
		int nel = permgather.NElements();
		int iel;
		for(iel=0; iel<nel; iel++)
			permgather[iel]=iel;
	}
    
    void computedirectionsPr(int inicio, int fim, TPZFMatrix<REAL> &bvec, TPZFMatrix<REAL> &t1vec,
                           TPZFMatrix<REAL> &t2vec, TPZFMatrix<REAL> &gradx, TPZFMatrix<REAL> &directions);
    
    void computedirectionsPr(int inicio, int fim, TPZFMatrix<REAL> &bvec, TPZFMatrix<REAL> &t1vec,
                           TPZFMatrix<REAL> &t2vec, TPZFMatrix<REAL> &gradx, TPZFMatrix<REAL> &directions)
    {
        REAL detgrad = 0.0;
        TPZVec<REAL> u(3);
        TPZVec<REAL> v(3);
        TPZVec<REAL> uxv(3);// result
        int cont = 0;
        
        for (int ivet=inicio; ivet<=fim; ivet++)
        {
            if (inicio<41)
            {
                for (int ilin=0; ilin<3; ilin++)
                {
                    u[ilin] = t1vec(ilin,ivet);
                    v[ilin] = t2vec(ilin,ivet);
                }
                TPZVec<REAL> e2(3);
                detgrad = 0.0;
                REAL normaX0xX1 = 0.0;
                TPZNumeric::ProdVetorial(u,v,e2);
                //            e2[0] = u[1]*v[2]-u[2]*v[1];
                //            e2[1] = -(u[0]*v[2]-v[0]*u[2]);
                //            e2[2] = u[0]*v[1]-v[0]*u[1];
                
                // calc do v gradx*b
                TPZManVector<REAL,3> dxt1(3,0.),dxt2(3,0.),dxt3(3,0.),Vvec(3,0.);
                //            REAL be2 = 0.0, ne2 = 0.0;
                //            for(int i=0;i<3;i++)
                //            {
                //                ne2 += e2[i]*e2[i];
                //            }
                //            ne2 = sqrt(fabs(ne2));
                for (int il=0; il<3; il++)
                {
                    for (int i = 0 ; i<3; i++)
                    {
                        dxt1[il] += gradx(il,i) * t1vec(i,ivet);
                        dxt2[il] += gradx(il,i) * t2vec(i,ivet);
                        //dxt3[il] += gradx(il,i) * e2[i]/ne2;
                        Vvec[il] += gradx(il,i) * bvec(i,ivet);
                    }
                    //be2 += bvec(il,ivet)*e2[il]/ne2;
                }
                TPZManVector<REAL,3> normal(3,0.);
                TPZNumeric::ProdVetorial(dxt1,dxt2,normal);
                //            normal[0] = dxt1[1]*dxt2[2]-dxt1[2]*dxt2[1];
                //            normal[1] = -(dxt1[0]*dxt2[2]-dxt2[0]*dxt1[2]);
                //            normal[2] = dxt1[0]*dxt2[1]-dxt2[0]*dxt1[1];
                
                for (int pos=0; pos<3; pos++)
                {
                    //                detgrad += normal[pos]*dxt3[pos];//uxv[pos]*gradx.GetVal(pos, 2);
                    normaX0xX1 += normal[pos]*normal[pos]; //uxv[pos]*uxv[pos];
                }
                TPZFMatrix<REAL> Wvec(3,1);
                
                REAL detgrad = gradx(0,0)*gradx(1,1)*gradx(2,2) + gradx(0,1)*gradx(1,2)*gradx(2,0) + gradx(0,2)*gradx(1,0)*gradx(2,1) - gradx(0,2)*gradx(1,1)*gradx(2,0) - gradx(0,0)*gradx(1,2)*gradx(2,1) - gradx(0,1)*gradx(1,0)*gradx(2,2);
                //detgrad = fabs(detgrad);
                
                normaX0xX1 = sqrt(normaX0xX1);
                if (detgrad<0)
                {
                    DebugStop();
                }
                
                for (int il=0; il<3; il++)
                {
                    Wvec(il,0) = Vvec[il]*normaX0xX1/(detgrad/**be2*/);
                    directions(il,cont) = Wvec(il,0);
                }
                cont++;
            }
            else
            {
                // calc do v gradx*b
                TPZManVector<REAL,3> Vvec(3,0.);
                for (int il=0; il<3; il++)
                {
                    for (int i = 0 ; i<3; i++)
                    {
                        Vvec[il] += gradx(il,i) * bvec(i,ivet);
                    }
                }
                for (int il=0; il<3; il++)
                {
                    directions(il,cont) = Vvec[il];
                }
                cont++;
            }

        }
            
        
    }
    
    void TPZPrism::ComputeDirections(int side, TPZFMatrix<REAL> &gradx, TPZFMatrix<REAL> &directions, TPZVec<int> &sidevectors)
    {
        // this method is out of date
        std::cout << __PRETTY_FUNCTION__ << "Deprecated method, not compatible with Piola transform\n";
        DebugStop();

        if(gradx.Cols()!=3)
        { std::cout << "Gradient dimensions are not compatible with this topology" << std::endl;
            DebugStop();
        }
        TPZFMatrix<REAL> bvec(3,63);
        int numvec = bvec.Cols();
        TPZFMatrix<REAL> t1vec(3,numvec);
        TPZFMatrix<REAL> t2vec(3,numvec);
        
        directions.Redim(3, numvec);
        for (int lin = 0; lin<numvec; lin++)
        {
            for(int col = 0;col<3;col++)
            {
                bvec.PutVal(col, lin, bPrism[lin][col]);
                t1vec.PutVal(col, lin, t1Prism[lin][col]);
                t2vec.PutVal(col, lin, t2Prism[lin][col]);
            }
        }
        
        // calcula os vetores
        switch (side) {
            case 15:
            {
                directions.Resize(3, 7);
                sidevectors.Resize(7);
                int inicio = 0, fim = 6;
                computedirectionsPr( inicio, fim, bvec, t1vec, t2vec, gradx, directions);
                int diff = fim-inicio+1;
                for (int ip = 0; ip < diff; ip++) {
                    sidevectors[ip] = vectorsideorderPr[ip+inicio];
                }
            }
                break;
            case 16:
            {
                directions.Resize(3, 9);
                sidevectors.Resize(9);
                int inicio = 7, fim = 15;
                computedirectionsPr( inicio, fim, bvec, t1vec, t2vec, gradx, directions);
                int diff = fim-inicio+1;
                for (int ip = 0; ip < diff; ip++) {
                    sidevectors[ip] = vectorsideorderPr[ip+inicio];
                }
            }
                break;
            case 17:
            {
                directions.Resize(3, 9);
                sidevectors.Resize(9);
                int inicio = 16, fim = 24;
                computedirectionsPr( inicio, fim, bvec, t1vec, t2vec, gradx, directions);
                int diff = fim-inicio+1;
                for (int ip = 0; ip < diff; ip++) {
                    sidevectors[ip] = vectorsideorderPr[ip+inicio];
                }
            }
                break;
            case 18:
            {
                directions.Resize(3, 9);
                sidevectors.Resize(9);
                int inicio = 25, fim = 33;
                computedirectionsPr( inicio, fim, bvec, t1vec, t2vec, gradx, directions);
                int diff = fim-inicio+1;
                for (int ip = 0; ip < diff; ip++) {
                    sidevectors[ip] = vectorsideorderPr[ip+inicio];
                }
            }
                break;
            case 19:
            {
                directions.Resize(3, 7);
                sidevectors.Resize(7);
                int inicio = 34, fim = 40;
                computedirectionsPr( inicio, fim, bvec, t1vec, t2vec, gradx, directions);
                int diff = fim-inicio+1;
                for (int ip = 0; ip < diff; ip++) {
                    sidevectors[ip] = vectorsideorderPr[ip+inicio];
                }
            }
                break;
            case 20:
            {
                directions.Resize(3, 22);
                sidevectors.Resize(22);
                int inicio = 41, fim = 62;
                computedirectionsPr( inicio, fim, bvec, t1vec, t2vec, gradx, directions);
                int diff = fim-inicio+1;
                for (int ip = 0; ip < diff; ip++) {
                    sidevectors[ip] = vectorsideorderPr[ip+inicio];
                }
            }
                break;
                
            default:
                break;
        }

                
	}
    
    template <class TVar>
    void TPZPrism::ComputeHDivDirections(TPZFMatrix<TVar> &gradx, TPZFMatrix<TVar> &directions)
    {
        TVar detjac = TPZAxesTools<TVar>::ComputeDetjac(gradx);
        
        TPZManVector<TVar,3> v1(3),v2(3),v3(3),v1v2(3),v3v1(3),v2v3(3),vec1(3),vec2(3),vec3(3), vdiag(3), v3vdiag(3);
        for (int i=0; i<3; i++) {
            v1[i] = gradx(i,0);
            v2[i] = gradx(i,1);
            v3[i] = gradx(i,2);
            vdiag[i] = (gradx(i,0)-gradx(i,1));
        }

        
        TPZNumeric::ProdVetorial(v1,v2,v1v2);
        TPZNumeric::ProdVetorial(v2,v3,v2v3);
        TPZNumeric::ProdVetorial(v3,v1,v3v1);
        TPZNumeric::ProdVetorial(v3,vdiag,v3vdiag);
        
        TVar Nv1v2 = TPZNumeric::Norm(v1v2);
        TVar Nv2v3 = TPZNumeric::Norm(v2v3);
        TVar Nv3v1 = TPZNumeric::Norm(v3v1);
        TVar Nv3vdiag = TPZNumeric::Norm(v3vdiag);
        
        /**
         * @file
         * @brief Computing mapped vector with scaling factor equal 1.0.
         * using contravariant piola mapping.
         */
        TPZManVector<TVar,3> NormalScales(4,1.);
        
        
        {
            Nv1v2 = 1.;
            Nv2v3 = 2.;
            Nv3v1 = 2.;
            Nv3vdiag = 2.;
        }
        
        
        for (int i=0; i<3; i++) {
            v1[i] /= detjac;
            v2[i] /= detjac;
            v3[i] /= detjac;
        }
        for (int i=0; i<3; i++)
        {
            for (int iv=0; iv<7; iv++)
            {
                directions(i,iv) = -v3[i]*NormalScales[0]*6.;
                directions(i,iv+34) = v3[i]*NormalScales[0]*6.;
            }
            //face 1
            directions(i,7) = -v2[i]*Nv3v1*NormalScales[2];
            directions(i,8) = (v1[i]-v2[i])*Nv3v1*NormalScales[2];
            directions(i,9) = (v1[i]-v2[i])*Nv3v1*NormalScales[2];
            directions(i,10) = -v2[i]*Nv3v1*NormalScales[2];
            directions(i,11) = ( directions(i,7)+directions(i,8) )/2.;
            directions(i,12) = (v1[i]-v2[i])*Nv3v1*NormalScales[2];
            directions(i,13) = ( directions(i,9)+directions(i,10) )/2.;
            directions(i,14) = -v2[i]*Nv3v1*NormalScales[2];
            directions(i,15) = ( directions(i,12)+directions(i,14) )/2.;
            //face 2
            directions(i,16) = v1[i]*Nv3vdiag*NormalScales[3];
            directions(i,17) = v2[i]*Nv3vdiag*NormalScales[3];
            directions(i,18) = v2[i]*Nv3vdiag*NormalScales[3];
            directions(i,19) = v1[i]*Nv3vdiag*NormalScales[3];
            directions(i,20) = (directions(i,16) + directions(i,17))/2.;
            directions(i,21) = v2[i]*Nv3vdiag*NormalScales[3];
            directions(i,22) = (directions(i,18) + directions(i,19))/2.;
            directions(i,23) = v1[i]*Nv3vdiag*NormalScales[3];
            directions(i,24) = (directions(i,21) + directions(i,23))/2.;
            //face 3
            directions(i,25) = -v1[i]*Nv2v3*NormalScales[1];
            directions(i,26) = (v2[i]-v1[i])*Nv2v3*NormalScales[1];
            directions(i,27) = (v2[i]-v1[i])*Nv2v3*NormalScales[1];
            directions(i,28) = -v1[i]*Nv2v3*NormalScales[1];
            directions(i,29) = ( directions(i,25)+directions(i,26) )/2.;
            directions(i,30) = (v2[i]-v1[i])*Nv2v3*NormalScales[1];
            directions(i,31) = ( directions(i,27)+directions(i,28) )/2.;
            directions(i,32) = -v1[i]*Nv2v3*NormalScales[1];
            directions(i,33) = ( directions(i,30)+directions(i,32) )/2.;
            
            //arestas
            directions(i,41) = v1[i]*Nv2v3*NormalScales[1];//
            directions(i,42) = (v2[i]-v1[i])/2.;//*Nvdiag
            directions(i,43) = -v2[i]*Nv3v1*NormalScales[2];//
            
            directions(i,44) = v3[i]*Nv1v2*NormalScales[0];
            directions(i,45) = v3[i]*Nv1v2*NormalScales[0];
            directions(i,46) = v3[i]*Nv1v2*NormalScales[0];
            
            directions(i,47) = v1[i]*Nv2v3*NormalScales[1];//
            directions(i,48) = (v2[i]-v1[i])/2.;//*Nvdiag
            directions(i,49) = -v2[i]*Nv3v1*NormalScales[2];//
            
            // internal in faces
            directions(i,50) = v1[i]*Nv2v3*NormalScales[1];
            directions(i,51) = v2[i]*Nv3v1*NormalScales[2];
            directions(i,52) = v1[i]*Nv2v3*NormalScales[1]/2.;//
            directions(i,53) = v3[i]*Nv1v2*NormalScales[0];
            directions(i,54) = (v2[i]-v1[i])/2.;//*Nvdiag
            directions(i,55) = v3[i]*Nv1v2*NormalScales[0];
            directions(i,56) = v2[i]*Nv3v1*NormalScales[2]/2.;
            directions(i,57) = v3[i]*Nv1v2*NormalScales[0];
            directions(i,58) = v1[i]*Nv2v3*NormalScales[1];
            directions(i,59) = v2[i]*Nv3v1*NormalScales[2];
            
            
            directions(i,60) = v1[i]*Nv2v3*NormalScales[1];
            directions(i,61) = v2[i]*Nv3v1*NormalScales[2];
            directions(i,62) = v3[i]*Nv1v2*NormalScales[0];
            
        }
        
    }

    /// Compute the directions of the HDiv vectors
    // template <class TVar>
    void TPZPrism::ComputeConstantHDiv(TPZVec<REAL> &point, TPZFMatrix<REAL> &RT0function, TPZVec<REAL> &div)
    {
        REAL scale;        
        
        REAL qsi = point[0];
        REAL eta = point[1];
        REAL zeta = point[2];

        //Face functions
        //For each face function: compute div = \nabla \cdot RT0function = d_RT0/d_qsi + d_RT0/d_eta 

        // Top and bottom is the same as cube
        scale = 0.5;
        RT0function(2,0) = -0.5 * (1. - zeta) / scale;
        RT0function(2,4) = 0.5 * (1. + zeta) / scale;
        div[0] = 0.5 / scale;
        div[4] = 0.5 / scale;

        //Faces are the same as triangles
        scale = 2.;
        RT0function(0,1) = qsi / scale;
        RT0function(1,1) = (eta - 1.) / scale;   
        div[1] = 2./scale;

        scale = M_SQRT2 * 2.;
        RT0function(0,2) = M_SQRT2 * qsi / scale;
        RT0function(1,2) = M_SQRT2 * eta / scale;
        div[2] = 2.* M_SQRT2 / scale;

        scale = 2.;
        RT0function(0,3) = (qsi - 1.) / scale;
        RT0function(1,3) = eta / scale;
        div[3] = 2. / scale;
    }

    /// Compute the directions of the HCurl vectors
    // template <class TVar>
    void TPZPrism::ComputeConstantHCurl(TPZVec<REAL> &point, TPZFMatrix<REAL> &N0function, TPZFMatrix<REAL> &curl, const TPZVec<int> &transformationIds)
    {
        REAL scale = 1.;    
        REAL qsi = point[0];
        REAL eta = point[1];
        REAL zeta = point[2];

        constexpr auto nEdges{9};
        TPZManVector<REAL,nEdges> edgeSign(nEdges,0);
        for(auto iEdge = 0; iEdge < nEdges; iEdge++){
            edgeSign[iEdge] = transformationIds[iEdge] == 0 ? 1 : -1;
        }

        //First type Nedelec functions
        //The three first and three last functions are the same as triangle, multiplied by z direction.
        N0function(0,0) = 0.5 * (1. - eta) * (1. - zeta) / scale * edgeSign[0];
        N0function(1,0) = 0.5 * qsi * (1. - zeta) / scale * edgeSign[0];
        curl(0,0) = 0.5 * qsi / scale * edgeSign[0];
        curl(1,0) = -0.5 * (1. - eta) / scale * edgeSign[0];
        curl(2,0) = (1. - zeta) / scale * edgeSign[0];

        N0function(0,1) = 0.5 * (-eta) * (1. - zeta) / scale * edgeSign[1];
        N0function(1,1) = 0.5 * qsi * (1. - zeta) / scale * edgeSign[1];
        curl(0,1) = 0.5 * qsi / scale * edgeSign[1];
        curl(1,1) = 0.5 * eta / scale * edgeSign[1];
        curl(2,1) = (1. - zeta) / scale * edgeSign[1];

        N0function(0,2) = -0.5 * eta * (1. - zeta) / scale * edgeSign[2];
        N0function(1,2) = 0.5 * (qsi - 1.) * (1. - zeta) / scale * edgeSign[2];
        curl(0,2) = 0.5 * (qsi - 1.) / scale * edgeSign[2];
        curl(1,2) = 0.5 * eta / scale * edgeSign[2];
        curl(2,2) = -(zeta - 1.) / scale * edgeSign[2];

        N0function(0,6) = 0.5 * (1. - eta) * (1. + zeta) / scale * edgeSign[6];
        N0function(1,6) = 0.5 * qsi * (1. + zeta) / scale * edgeSign[6];
        curl(0,6) = -0.5 * qsi / scale * edgeSign[6];
        curl(1,6) = 0.5 * (1. - eta) / scale * edgeSign[6];
        curl(2,6) = (1. + zeta) / scale * edgeSign[6];

        N0function(0,7) = 0.5 * (-eta) * (1. + zeta) / scale * edgeSign[7];
        N0function(1,7) = 0.5 * qsi * (1. + zeta) / scale * edgeSign[7];
        curl(0,7) = -0.5 * qsi / scale * edgeSign[7];
        curl(1,7) = -0.5 * eta / scale * edgeSign[7];
        curl(2,7) = (1. + zeta) / scale * edgeSign[7];

        N0function(0,8) = -0.5 * eta * (1. + zeta) / scale * edgeSign[8];
        N0function(1,8) = 0.5 * (qsi - 1.) * (1. + zeta) / scale * edgeSign[8];
        curl(0,8) = -0.5 * (qsi - 1.) / scale * edgeSign[8];
        curl(1,8) = -0.5 * eta / scale * edgeSign[8];
        curl(2,8) = -(-zeta - 1.) / scale * edgeSign[8];

        //The three vertical edges (only the z component is != 0)
        scale = 2.;
        N0function(2,3) = (1. - qsi - eta) / scale * edgeSign[3];
        curl(0,3) = -1. / scale * edgeSign[3];
        curl(1,3) =  1. / scale * edgeSign[3];
        curl(2,3) = 0.;

        N0function(2,4) = qsi / scale * edgeSign[4];
        curl(0,4) = 0.;
        curl(1,4) = -1. / scale * edgeSign[4];
        curl(2,4) = 0.;

        N0function(2,5) = eta / scale * edgeSign[5];
        curl(0,5) = 1. / scale * edgeSign[5];
        curl(1,5) = 0.;
        curl(2,5) = 0.;

    }

    // Get face orientation
    int TPZPrism::GetSideOrient(const int &face){
        return fSideOrient[face];
    }

    void TPZPrism::GetSideHDivDirections(TPZVec<int> &sides, TPZVec<int> &dir, TPZVec<int> &bilounao)
    {
        int nsides = NumSides()*3;
        
        sides.Resize(nsides);
        dir.Resize(nsides);
        bilounao.Resize(nsides);
        
        for (int is = 0; is<nsides; is++)
        {
            sides[is] = vectorsideorderPr[is];
            dir[is] = direcaoksioueta[is];
            bilounao[is] = bilinearounao[is];
        }
    }


    void TPZPrism::GetSideHDivDirections(TPZVec<int> &sides, TPZVec<int> &dir, TPZVec<int> &bilounao, TPZVec<int> &sidevectors) {
        int nsides = NumSides()*3;

        sides.Resize(nsides);
        dir.Resize(nsides);
        bilounao.Resize(nsides);

        for (int is = 0; is<nsides; is++)
        {
            sides[is] = vectorsideorderPr[is];
            dir[is] = direcaoksioueta[is];
            bilounao[is] = bilinearounao[is];
        }
        for (int i=0; i<Dimension*NumSides(); i++) {
            sidevectors[i] = vectorsideorderPr[i];
        }
    }

    template <class TVar>
    void TPZPrism::ComputeHCurlDirections(TPZFMatrix<TVar> &gradx, TPZFMatrix<TVar> &directions, const TPZVec<int> &transformationIds)
    {
        TPZManVector<TVar,3> v1(3),v2(3),v3(3);

        for (int i=0; i<3; i++) {
            v1[i] = gradx(i,0);
            v2[i] = gradx(i,1);
            v3[i] = gradx(i,2);
        }
        constexpr auto nEdges = 9;
        constexpr auto nFaces = 5;
        //edges                      6,  7    ,8,9,10,11,12,   13,  14
        constexpr REAL edgeLength[nEdges]{1,M_SQRT2,1,2, 2, 2, 1,M_SQRT2,1};
        //faces                    15,16,    17,   18,19
        constexpr REAL faceArea[nFaces]{0.5,2,2*M_SQRT2,2,0.5};
        TPZManVector<REAL,nEdges> edgeSign(nEdges,0);
        TPZManVector<REAL,nFaces>  faceOrient(nFaces,0);

        for(auto iEdge = 0; iEdge < nEdges; iEdge++){
            edgeSign[iEdge] = transformationIds[iEdge] == 0 ? 1 : -1;
        }
        for(auto iFace = 0; iFace < nFaces; iFace++){
            faceOrient[iFace] = transformationIds[nEdges + iFace] % 2 == 0 ? 1 : -1;
        }
        for (int i=0; i<3; i++) {
            //v^{e,a} constant vector fields associated with edge e and vertex a
            //they are defined in such a way that v^{e,a} is normal to the edge \hat{e}
            //adjacent to edge e by the vertex a. the tangential component is set to be 1 /edgeLength[e]

            directions(i, 0) = (v1[i]) * edgeSign[0] / edgeLength[0];//edge 6  vertex 0
            directions(i, 1) = (v1[i] + v2[i]) * edgeSign[0] / edgeLength[0];//edge 6 vertex 1
            directions(i, 2) = (v2[i] * M_SQRT2) * edgeSign[1] / edgeLength[1];//edge 7 vertex 1
            directions(i, 3) = (-v1[i] * M_SQRT2) * edgeSign[1] / edgeLength[1];//edge 7 vertex 2
            directions(i, 4) = (v1[i] + v2[i]) * -1 * edgeSign[2] / edgeLength[2]; //edge 8 vertex 2
            directions(i, 5) = (-v2[i]) * edgeSign[2] / edgeLength[2]; //edge 8 vertex 0

            directions(i, 6) = v3[i] * edgeSign[3] / edgeLength[3];//edge 9 vertex 0
            directions(i, 7) = v3[i] * edgeSign[3] / edgeLength[3];//edge 9 vertex 3
            directions(i, 8) = v3[i] * edgeSign[4] / edgeLength[4];//edge 10 vertex 1
            directions(i, 9) = v3[i] * edgeSign[4] / edgeLength[4];//edge 10 vertex 4
            directions(i, 10) = v3[i] * edgeSign[5] / edgeLength[5];//edge 11 vertex 2
            directions(i, 11) = v3[i] * edgeSign[5] / edgeLength[5];//edge 11 vertex 5

            directions(i, 12) = (v1[i]) * edgeSign[6] / edgeLength[6];//edge 12 vertex 3
            directions(i, 13) = (v1[i] + v2[i]) * edgeSign[6] / edgeLength[6];//edge 12 vertex 4
            directions(i, 14) = (v2[i] * M_SQRT2) * edgeSign[7] / edgeLength[7];//edge 13 vertex 4
            directions(i, 15) = (-v1[i] * M_SQRT2) * edgeSign[7] / edgeLength[7];//edge 13 vertex 5
            directions(i, 16) = (v1[i] + v2[i]) * -1 * edgeSign[8] / edgeLength[8];//edge 14 vertex 5
            directions(i, 17) = (-v2[i]) * edgeSign[8] / edgeLength[8];//edge 14 vertex 3

            //v^{e,T} constant vector fields associated with edge e and aligned with it
            directions(i, 18) = (v1[i]+0.5*v2[i]) * edgeSign[0] / edgeLength[0];//edge 6
            directions(i, 19) = ((v2[i] - v1[i]) * M_SQRT1_2) * edgeSign[1] / edgeLength[1];//edge 7
            directions(i, 20) = (-0.5*v1[i]-v2[i]) * edgeSign[2] / edgeLength[2];//edge 8

            directions(i, 21) = v3[i] * edgeSign[3] / edgeLength[3];//edge 9
            directions(i, 22) = v3[i] * edgeSign[4] / edgeLength[4];//edge 10
            directions(i, 23) = v3[i] * edgeSign[5] / edgeLength[5];//edge 11

            directions(i, 24) = (v1[i]+0.5*v2[i]) * edgeSign[6] / edgeLength[6];//edge 12
            directions(i, 25) = ((v2[i] - v1[i]) * M_SQRT1_2) * edgeSign[7] / edgeLength[7];//edge 13
            directions(i, 26) = (-0.5*v1[i]-v2[i]) * edgeSign[8] / edgeLength[8];//edge 14

            //v^{F,e} constant vector fields associated with face F and edge e
            //they are defined in such a way that v^{F,e} is normal to the face \hat{F}
            //adjacent to face F by edge e
            directions(i, 27) = faceOrient[0] * edgeSign[6 - NCornerNodes] * v2[i]  / faceArea[0];//face 15 edge 6
            directions(i, 28) = faceOrient[0] * edgeSign[7 - NCornerNodes] * -1 * (v1[i] + v2[i])  / faceArea[0];//face 15 edge 7
            directions(i, 29) = faceOrient[0] * edgeSign[8 - NCornerNodes] * v1[i]  / faceArea[0];//face 15 edge 8

            directions(i, 30) = faceOrient[1] * edgeSign[6 - NCornerNodes] *  0.5 * v3[i] / faceArea[1];//face 16 edge 6
            directions(i, 31) = faceOrient[1] * edgeSign[10 - NCornerNodes] *  -1 * (v1[i]+v2[i]) / faceArea[1];//face 16 edge 10
            directions(i, 32) = faceOrient[1] * edgeSign[12 - NCornerNodes] *  0.5 * v3[i] / faceArea[1];;//face 16 edge 12
            directions(i, 33) = faceOrient[1] * edgeSign[9 - NCornerNodes] *  -v1[i] / faceArea[1];//face 16 edge 9

            directions(i, 34) = faceOrient[2] * edgeSign[7 - NCornerNodes] *  M_SQRT1_2 * v3[i] / faceArea[2];//face 17 edge 7
            directions(i, 35) = faceOrient[2] * edgeSign[11 - NCornerNodes] *  M_SQRT2 * v1[i] / faceArea[2];//face 17 edge 11
            directions(i, 36) = faceOrient[2] * edgeSign[13 - NCornerNodes] *  M_SQRT1_2 * v3[i] / faceArea[2];//face 17 edge 13
            directions(i, 37) = faceOrient[2] * edgeSign[10 - NCornerNodes] *  -M_SQRT2 * v2[i] / faceArea[2];//face 17 edge 10

            directions(i, 38) = faceOrient[3] * edgeSign[8 - NCornerNodes] * -0.5 * v3[i] / faceArea[3];//face 18 edge 8
            directions(i, 39) = faceOrient[3] * edgeSign[11 - NCornerNodes] *  -1 * (v1[i] + v2[i]) / faceArea[3];//face 18 edge 11
            directions(i, 40) = faceOrient[3] * edgeSign[14 - NCornerNodes] * -0.5 * v3[i] / faceArea[3];//face 18 edge 14
            directions(i, 41) = faceOrient[3] * edgeSign[9 - NCornerNodes] * -v2[i] / faceArea[3];//face 18 edge 9

            directions(i, 42) = faceOrient[4] * edgeSign[12 - NCornerNodes] * v2[i]  / faceArea[4];//face 19 edge 12
            directions(i, 43) = faceOrient[4] * edgeSign[13 - NCornerNodes] * -1 * (v1[i] + v2[i] )  / faceArea[4];//face 19 edge 13
            directions(i, 44) = faceOrient[4] * edgeSign[14 - NCornerNodes] * v1[i]  / faceArea[4];//face 19 edge 14

            //v^{F,T} are calculated afterwards

            //v^{F,orth} vector associated with face F and normal to it
            directions(i, 55) = -v3[i];//face 15
            directions(i, 56) = -v2[i];//face 16
            directions(i, 57) = (v1[i] + v2[i]) * M_SQRT1_2;//face 17
            directions(i, 58) = -v1[i];//face 18
            directions(i, 59) = v3[i];//face 19

            //v^{K,3}
            directions(i, 60) = v1[i];
            directions(i, 61) = v2[i];
            directions(i, 62) = v3[i];
        }

        TPZManVector<REAL,2> vft1(2,0), vft2(2,0);
        constexpr auto firstVftVec = 45;
        //v^{F,T} orthonormal vectors associated with face F and tangent to it.
        for(auto iFace = 0; iFace < nFaces; iFace ++){
            switch(iFace){
                case 0:
                case 4:
                    TPZTriangle::ComputeHCurlFaceDirections(vft1,vft2,transformationIds[nEdges + iFace]);
                    break;
                case 1:
                case 2:
                case 3:
                    TPZQuadrilateral::ComputeHCurlFaceDirections(vft1,vft2,transformationIds[nEdges + iFace]);
                    break;
            }
            directions(0,firstVftVec+2*iFace) = 0;directions(1,firstVftVec+2*iFace) = 0;directions(2,firstVftVec+2*iFace) = 0;
            directions(0,firstVftVec+2*iFace+1) = 0;directions(1,firstVftVec+2*iFace+1) = 0;directions(2,firstVftVec+2*iFace+1) = 0;
            auto axes = TPZPrism::TransformElementToSide(NCornerNodes+nEdges+iFace).Mult();
            axes.Transpose();
            for(auto x = 0; x < Dimension; x++){
                for(auto i = 0; i < 2; i++) {
                    directions(x, firstVftVec + 2 * iFace) += axes(x,i) * vft1[i];
                    directions(x, firstVftVec + 2 * iFace + 1) += axes(x,i) * vft2[i];
                }
            }
        }
    }

    int TPZPrism::ClassId() const{
        return Hash("TPZPrism");
    }

    void TPZPrism::Read(TPZStream& buf, void* context) {

    }

    void TPZPrism::Write(TPZStream& buf, int withclassid) const {

    }

}

/**********************************************************************************************************************
 * The following are explicit instantiation of member function template of this class, both with class T=REAL and its
 * respective FAD<REAL> version. In other to avoid potential errors, always declare the instantiation in the same order
 * in BOTH cases.    @orlandini
 **********************************************************************************************************************/
template bool pztopology::TPZPrism::CheckProjectionForSingularity<REAL>(const int &side, const TPZVec<REAL> &xiInterior);

template void pztopology::TPZPrism::MapToSide<REAL>(int side, TPZVec<REAL> &InternalPar, TPZVec<REAL> &SidePar, TPZFMatrix<REAL> &JacToSide);

template void pztopology::TPZPrism::BlendFactorForSide<REAL>(const int &, const TPZVec<REAL> &, REAL &, TPZVec<REAL> &);

template void pztopology::TPZPrism::TShape<REAL>(const TPZVec<REAL> &loc,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);

template void pztopology::TPZPrism::ComputeHDivDirections<REAL>(TPZFMatrix<REAL> &gradx, TPZFMatrix<REAL> &directions);

template void pztopology::TPZPrism::ComputeHCurlDirections<REAL>(TPZFMatrix<REAL> &gradx, TPZFMatrix<REAL> &directions, const TPZVec<int> &transformationIds);

template bool pztopology::TPZPrism::CheckProjectionForSingularity<Fad<REAL>>(const int &side, const TPZVec<Fad<REAL>> &xiInterior);

template void pztopology::TPZPrism::MapToSide<Fad<REAL> >(int side, TPZVec<Fad<REAL> > &InternalPar, TPZVec<Fad<REAL> > &SidePar, TPZFMatrix<Fad<REAL> > &JacToSide);

template void pztopology::TPZPrism::BlendFactorForSide<Fad<REAL>>(const int &, const TPZVec<Fad<REAL>> &, Fad<REAL> &,
                                                                   TPZVec<Fad<REAL>> &);
template void pztopology::TPZPrism::TShape<Fad<REAL>>(const TPZVec<Fad<REAL>> &loc,TPZFMatrix<Fad<REAL>> &phi,TPZFMatrix<Fad<REAL>> &dphi);

template void pztopology::TPZPrism::ComputeHDivDirections<Fad<REAL>>(TPZFMatrix<Fad<REAL>> &gradx, TPZFMatrix<Fad<REAL>> &directions);

template void pztopology::TPZPrism::ComputeHCurlDirections<Fad<REAL>>(TPZFMatrix<Fad<REAL>> &gradx, TPZFMatrix<Fad<REAL>> &directions, const TPZVec<int> &transformationIds);
