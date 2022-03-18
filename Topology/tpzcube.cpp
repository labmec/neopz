/**
 * @file
 * @brief Contains the implementation of the TPZCube methods. 
 */

#include "tpzcube.h"
#include "pzmanvector.h"
#include "pzerror.h"
#include "pzreal.h"
#include "pzquad.h"
#include "pzeltype.h"
#include "tpzquadrilateral.h"
#include "pznumeric.h"


#include "fad.h"

#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.topology.pzcube");
#endif
using namespace std;

namespace pztopology {

	/**
	 * @brief For each face (quadrilateral sides) was enumerated in sequence the sides contained in the closure of them.
	 * For example: First face (side 20) in its closure contains side 0, 1, 2, 3 (vertices), 8, 9, 10, 11, (edges) and it self 
	 */
	static constexpr int FaceConnectLocId[6][9] = { {0,1,2,3,8,9,10,11,20},{0,1,5,4,8,13,16,12,21},
		{1,2,6,5,9,14,17,13,22},{3,2,6,7,10,14,18,15,23},//{2,3,7,6,10,15,18,14,23}
		{0,3,7,4,11,15,19,12,24},{4,5,6,7,16,17,18,19,25} };


	/** @brief Vector of the dimension for each side */
	static constexpr int sidedimension[27] = {0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,3};
	
    static constexpr int fSideOrient[6] = {-1,1,1,-1,-1,1};

	/** @brief Vector with the number of vertices contained in the closure of the side */
	static constexpr int nsidenodes[27] = {1,1,1,1,1,1,1,1,
		2,2,2,2,2,2,2,2,2,2,2,2,
		4,4,4,4,4,4,
		8};
	

	/** @brief For each side was stored the sides connected with it but of the higher dimension */ 
	static constexpr int highsides[27][7] = {
		{8,11,12,20,21,24,26},
		{8,9,13,20,21,22,26},
		{9,10,14,20,22,23,26},
		{10,11,15,20,23,24,26},
		{12,16,19,21,24,25,26},
		{13,16,17,21,22,25,26},
		{14,17,18,22,23,25,26},
		{15,18,19,23,24,25,26},
		{20,21,26},
		{20,22,26},
		{20,23,26},
		{20,24,26},
		{21,24,26},
		{21,22,26},
		{22,23,26},
		{23,24,26},
		{21,25,26},
		{22,25,26},
		{23,25,26},
		{24,25,26},
		{26},
		{26},
		{26},
		{26},
		{26},
		{26},
		{-999}
	};

	/**
	 * @brief For each side was stored the number of sides connected of the higher dimension than it self 
	 * For example: First side (side 0 - vertice) was connected with the sides 8, 11, 12 (edges) 20, 21, 24 (faces) and 26 (the hexahedra). At total 7 sides */
	static constexpr int nhighdimsides[27] = {7,7,7,7,7,7,7,7,3,3,3,3,3,3,3,3,3,3,3,3,1,1,1,1,1,1,0};

	/** @brief The transformations for each side over neighboard side with higher dimension */
	static constexpr REAL sidetosidetransforms[27][7][4][3] = {
		{
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-1,-1}}
		},
		{
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-1,-1}}
		},
		{
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,1,-1}}
		},
		{
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,1,-1}}
		},
		{
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-1,1}}
		},
		{
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-1,1}}
		},
		{
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,1,1}}
		},
		{
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,-99,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{1,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,1,-99}},
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-1,1,1}}
		},
		{
			{{1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,-1,-99}},
			{{1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,-1,-99}},
			{{1,0,0},{-99,-99,-99},{-99,-99,-99},{0,-1,-1}}
		},
		{
			{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
			{{1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,-1,-99}},
			{{0,1,0},{-99,-99,-99},{-99,-99,-99},{1,0,-1}}
		},
		{
			{{-1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
			{{-1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,-1,-99}},
			{{-1,0,0},{-99,-99,-99},{-99,-99,-99},{0,1,-1}}
		},
		{
			{{0,-1,-99},{-99,-99,-99},{-99,-99,-99},{-1,0,-99}},
			{{-1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,-1,-99}},
			{{0,-1,0},{-99,-99,-99},{-99,-99,-99},{-1,0,-1}}
		},
		{
			{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{-1,0,-99}},
			{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{-1,0,-99}},
			{{0,0,1},{-99,-99,-99},{-99,-99,-99},{-1,-1,0}}
		},
		{
			{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
			{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{-1,0,-99}},
			{{0,0,1},{-99,-99,-99},{-99,-99,-99},{1,-1,0}}
		},
		{
			{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
			{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
			{{0,0,1},{-99,-99,-99},{-99,-99,-99},{1,1,0}}
		},
		{
			{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{-1,0,-99}},
			{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
			{{0,0,1},{-99,-99,-99},{-99,-99,-99},{-1,1,0}}
		},
		{
			{{1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
			{{1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,-1,-99}},
			{{1,0,0},{-99,-99,-99},{-99,-99,-99},{0,-1,1}}
		},
		{
			{{1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
			{{0,1,-99},{-99,-99,-99},{-99,-99,-99},{1,0,-99}},
			{{0,1,0},{-99,-99,-99},{-99,-99,-99},{1,0,1}}
		},
		{
			{{-1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
			{{-1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
			{{-1,0,0},{-99,-99,-99},{-99,-99,-99},{0,1,1}}
		},
		{
			{{-1,0,-99},{-99,-99,-99},{-99,-99,-99},{0,1,-99}},
			{{0,-1,-99},{-99,-99,-99},{-99,-99,-99},{-1,0,-99}},
			{{0,-1,0},{-99,-99,-99},{-99,-99,-99},{-1,0,1}}
		},
		{
			{{1,0,0},{0,1,0},{-99,-99,-99},{0,0,-1}}
		},
		{
			{{1,0,0},{0,0,1},{-99,-99,-99},{0,-1,0}}
		},
		{
			{{0,1,0},{0,0,1},{-99,-99,-99},{1,0,0}}
		},
		{
			{{1,0,0},{0,0,1},{-99,-99,-99},{0,1,0}}
		},
		{
			{{0,1,0},{0,0,1},{-99,-99,-99},{-1,0,0}}
		},
		{
			{{1,0,0},{0,1,0},{-99,-99,-99},{0,0,1}}
		},
		{
			{{-99,-99,-99},{-99,-99,-99},{-99,-99,-99},{-99,-99,-99}}
		}
	};
	
	/** @brief For each side the vector related the coordinates of the point considered middle of the side */
	static constexpr REAL MidSideNode[27][3] = {
		/*00*/{-1.,-1.,-1.},/*01*/{1.,-1.,-1.},/*02*/{1.,1.,-1.},/*03*/{-1.,1.,-1.},
		/*04*/{-1.,-1., 1.},/*05*/{1.,-1., 1.},/*06*/{1.,1., 1.},/*07*/{-1.,1., 1.},
		/*08*/{ 0.,-1.,-1.},/*09*/{1., 0.,-1.},/*10*/{0.,1.,-1.},/*11*/{-1.,0.,-1.},
		/*12*/{-1.,-1., 0.},/*13*/{1.,-1., 0.},/*14*/{1.,1., 0.},/*15*/{-1.,1., 0.},
		/*16*/{ 0.,-1., 1.},/*17*/{1., 0., 1.},/*18*/{0.,1., 1.},/*19*/{-1.,0., 1.},
		/*20*/{ 0., 0.,-1.},/*21*/{0.,-1., 0.},/*22*/{1.,0., 0.},/*23*/{ 0.,1., 0.},
		/*24*/{-1., 0., 0.},/*25*/{0., 0., 1.},/*26*/{0.,0., 0.} };
    
    static constexpr REAL bCubo[81][3] = // direcao perpendicular ao lado
    {
        {0,0,-1}, {0,0,-1}, {0,0,-1}, {0,0,-1}, {0,0,-1}, {0,0,-1}, {0,0,-1}, {0,0,-1}, {0,0,-1},// face 0
        {0,-1,0}, {0,-1,0}, {0,-1,0}, {0,-1,0}, {0,-1,0}, {0,-1,0}, {0,-1,0}, {0,-1,0}, {0,-1,0},// face 1
        {1,0,0} , {1,0,0} , {1,0,0} , {1,0,0} , {1,0,0} , {1,0,0} , {1,0,0} , {1,0,0} , {1,0,0} ,// face 2
        {0,1,0} , {0,1,0} , {0,1,0} , {0,1,0} , {0,1,0} , {0,1,0} , {0,1,0} , {0,1,0} , {0,1,0} ,// face 3
        {-1,0,0}, {-1,0,0}, {-1,0,0}, {-1,0,0}, {-1,0,0}, {-1,0,0}, {-1,0,0}, {-1,0,0}, {-1,0,0},// face 4
        {0,0,1} , {0,0,1} , {0,0,1} , {0,0,1} , {0,0,1} , {0,0,1} , {0,0,1} , {0,0,1} , {0,0,1}, // face 5
        //interiores
        //arestas
        //{1,0,0},{0,1,0},{-1,0,0},{0,-1,0},  {0,0,1},{0,0,1},{0,0,1},{0,0,1},  {1,0,0},{0,1,0},{-1,0,0},{0,-1,0},
        {1,0,0},{0,1,0},{-1,0,0},{0,-1,0},  {0,0,1},{0,0,1},{0,0,1},{0,0,1},  {1,0,0},{0,1,0},{-1,0,0},{0,-1,0},
        //faces
        {1,0,0}, {0,1,0}, // tang da face 0
        {1,0,0}, {0,0,1},  // tang da face 1
        {0,1,0}, {0,0,1}, // tang da face 2
        {1,0,0}, {0,0,1}, // tang da face 3
        {0,1,0}, {0,0,1}, // tang da face 4
        {1,0,0}, {0,1,0},  // tang da face 5
        {1,0,0}, // volume
        {0,1,0}, // volume
        {0,0,1}  // volume
    };
    
    static constexpr REAL t1Cubo[81][3] = // diretor da aresta (escolhido para formar uma base positivamente orientada)
    {
        {-1,0,0},{-1,0,0},{-1,0,0},{-1,0,0},{-1,0,0},{-1,0,0},{-1,0,0},{-1,0,0},{-1,0,0},//face 0
        {1,0,0}, {1,0,0}, {1,0,0}, {1,0,0}, {1,0,0}, {1,0,0}, {1,0,0}, {1,0,0}, {1,0,0}, //face 1
        {0,0,-1},{0,0,-1},{0,0,-1},{0,0,-1},{0,0,-1},{0,0,-1},{0,0,-1},{0,0,-1},{0,0,-1},//face 2
        {1,0,0} ,{1,0,0} ,{1,0,0} ,{1,0,0} ,{1,0,0} ,{1,0,0} ,{1,0,0} ,{1,0,0} ,{1,0,0} ,//face 3
        {0,0,1} ,{0,0,1} ,{0,0,1} ,{0,0,1} ,{0,0,1} ,{0,0,1} ,{0,0,1} ,{0,0,1} ,{0,0,1} ,//face 4
        {1,0,0} ,{1,0,0} ,{1,0,0} ,{1,0,0} ,{1,0,0} ,{1,0,0} ,{1,0,0} ,{1,0,0} ,{1,0,0}, //face 5
        //interiores
        //arestas
        {0,-1,0},{1,0,0},{0,1,0},{1,0,0},  {-1,0,0},{0,-1,0},{1,0,0},{0,1,0},  {0,0,1},{0,0,1},{0,0,1},{0,0,1},
        //faces
        {0,-1,0}, {1,0,0},   //  complementar da face 0
        {0,0,1},  {-1,0,0},  //  complementar da face 1
        {1,0,0},  {1,0,0},   //  complementar da face 2
        {0,0,-1}, {1,0,0},   //  complementar da face 3
        {0,0,-1}, {0,1,0},   //  complementar da face 4
        {0,1,0},  {1,0,0},   //  complementar da face 5

        
        {0,1,0},  // volume
        {0,0,1},  // volume
        {1,0,0}   // volume
        
    };
    static constexpr REAL t2Cubo[81][3] = // diretor da aresta (escolhido para formar uma base positivamente orientada)
    {
        {0,1,0}, {0,1,0}, {0,1,0}, {0,1,0}, {0,1,0}, {0,1,0}, {0,1,0}, {0,1,0}, {0,1,0},// face 0
        {0,0,1}, {0,0,1}, {0,0,1}, {0,0,1}, {0,0,1}, {0,0,1}, {0,0,1}, {0,0,1}, {0,0,1},// face 1
        {0,1,0}, {0,1,0}, {0,1,0}, {0,1,0}, {0,1,0}, {0,1,0}, {0,1,0}, {0,1,0}, {0,1,0},// face 2
        {0,0,-1},{0,0,-1},{0,0,-1},{0,0,-1},{0,0,-1},{0,0,-1},{0,0,-1},{0,0,-1},{0,0,-1},// face 3
        {0,1,0}, {0,1,0}, {0,1,0}, {0,1,0}, {0,1,0}, {0,1,0}, {0,1,0}, {0,1,0}, {0,1,0},// face 4
        {0,1,0}, {0,1,0}, {0,1,0}, {0,1,0}, {0,1,0}, {0,1,0}, {0,1,0}, {0,1,0}, {0,1,0},// face 5
        //interiores
        //arestas
        {0,0,-1},{0,0,-1},{0,0,-1},{0,0,1},  {0,-1,0},{1,0,0},{0,1,0},{-1,0,0},  {0,-1,0},{1,0,0},{0,1,0},{-1,0,0},
        //faces
        {0,0,-1}, {0,0,-1}, // 2 complementar da face 0
        {0,-1,0}, {0,-1,0}, // 2 complementar da face 1
        {0,0,-1},  {0,1,0},  // 2 complementar da face 2
        {0,1,0},  {0,1,0},  // 2 complementar da face 3
        {-1,0,0}, {-1,0,0}, // 2 complementar da face 4
        {0,0,1},  {0,0,-1},  // 2 complementar da face 5
        {0,0,1},  // volume
        {1,0,0},  // volume
        {0,1,0}   // volume
    };
    
	static constexpr int vectorsideorderC [81] =
    {
        0,1,2,3,8,9,10,11,20, //face 0
        0,1,5,4,8,13,16,12,21,//face 1
        1,2,6,5,9,14,17,13,22,//face 2
        3,2,6,7,10,14,18,15,23,//face 3
        //2,3,7,6,10,15,18,14,23,//face 3
        0,3,7,4,11,15,19,12,24,//face 4
        4,5,6,7,16,17,18,19,25,//face 5
        8,9,10,11,
        12,13,14,15,
        16,17,18,19,
        20,20,//tg face 0
        21,21,//tg face 1
        22,22,//tg face 2
        23,23,//tg face 3
        24,24,//tg face 4
        25,25,//tg face 5
        26,26,26
    };

    static constexpr int bilinearounao [81] =   {
        0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,
        1
    };
    
    static constexpr int direcaoksioueta [81] = {
        0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,
        0,1,0,1,0,1,0,1,0,1,0,1,//0,1,0,2,1,2,0,2,1,2,0,1,
        0,1,2};

    template<class T>
    inline void TPZCube::TShape(const TPZVec<T> &loc,TPZFMatrix<T> &phi,TPZFMatrix<T> &dphi) {
        T qsi = loc[0], eta = loc[1] , zeta  = loc[2];

        T x[2],dx[2],y[2],dy[2],z[2],dz[2];
        x[0]    = (1.-qsi)/2.;
        x[1]    = (1.+qsi)/2.;
        dx[0]   = -0.5;
        dx[1]   = +0.5;
        y[0]    = (1.-eta)/2.;
        y[1]    = (1.+eta)/2.;
        dy[0]   = -0.5;
        dy[1]   = +0.5;
        z[0]    = (1.-zeta)/2.;
        z[1]    = (1.+zeta)/2.;
        dz[0]   = -0.5;
        dz[1]   = +0.5;

        phi(0,0) = x[0]*y[0]*z[0];
        phi(1,0) = x[1]*y[0]*z[0];
        phi(2,0) = x[1]*y[1]*z[0];
        phi(3,0) = x[0]*y[1]*z[0];
        phi(4,0) = x[0]*y[0]*z[1];
        phi(5,0) = x[1]*y[0]*z[1];
        phi(6,0) = x[1]*y[1]*z[1];
        phi(7,0) = x[0]*y[1]*z[1];
        dphi(0,0) = dx[0]*y[0]*z[0];
        dphi(1,0) = x[0]*dy[0]*z[0];
        dphi(2,0) = x[0]*y[0]*dz[0];
        dphi(0,1) = dx[1]*y[0]*z[0];
        dphi(1,1) = x[1]*dy[0]*z[0];
        dphi(2,1) = x[1]*y[0]*dz[0];
        dphi(0,2) = dx[1]*y[1]*z[0];
        dphi(1,2) = x[1]*dy[1]*z[0];
        dphi(2,2) = x[1]*y[1]*dz[0];
        dphi(0,3) = dx[0]*y[1]*z[0];
        dphi(1,3) = x[0]*dy[1]*z[0];
        dphi(2,3) = x[0]*y[1]*dz[0];
        dphi(0,4) = dx[0]*y[0]*z[1];
        dphi(1,4) = x[0]*dy[0]*z[1];
        dphi(2,4) = x[0]*y[0]*dz[1];
        dphi(0,5) = dx[1]*y[0]*z[1];
        dphi(1,5) = x[1]*dy[0]*z[1];
        dphi(2,5) = x[1]*y[0]*dz[1];
        dphi(0,6) = dx[1]*y[1]*z[1];
        dphi(1,6) = x[1]*dy[1]*z[1];
        dphi(2,6) = x[1]*y[1]*dz[1];
        dphi(0,7) = dx[0]*y[1]*z[1];
        dphi(1,7) = x[0]*dy[1]*z[1];
        dphi(2,7) = x[0]*y[1]*dz[1];

    }

    template<class T>
    void TPZCube::BlendFactorForSide(const int &side, const TPZVec<T> &xi, T &blendFactor,
                                       TPZVec<T> &corrFactorDxi){
        const REAL tol = pztopology::GetTolerance();
        std::ostringstream sout;
        if(side < NCornerNodes || side >= NSides){
            sout<<"The side\t"<<side<<"is invalid. Aborting..."<<std::endl;

            PZError<<std::endl<<sout.str()<<std::endl;
            DebugStop();
        }
        #ifdef PZDEBUG

        if(!IsInParametricDomain(xi,tol)){
            sout<<"The method BlendFactorForSide expects the point xi to correspond to coordinates of a point";
            sout<<" inside the parametric domain. Aborting...";
            PZError<<std::endl<<sout.str()<<std::endl;
            #ifdef PZ_LOG
            LOGPZ_FATAL(logger,sout.str().c_str());
            #endif
            DebugStop();
        }
        #endif
        corrFactorDxi.Resize(TPZCube::Dimension,(T)0);
        if(side < NSides - 1){
            TPZFNMatrix<4,T> phi(NCornerNodes,1);
            TPZFNMatrix<8,T> dphi(Dimension,NCornerNodes);
            TPZCube::TShape(xi,phi,dphi);
            blendFactor = 0;
            for(int i = 0; i < TPZCube::NSideNodes(side);i++){
                const int currentNode = TPZCube::SideNodeLocId(side, i);
                blendFactor += phi(currentNode,0);
                corrFactorDxi[0] +=  dphi(0,currentNode);
                corrFactorDxi[1] +=  dphi(1,currentNode);
                corrFactorDxi[2] +=  dphi(2,currentNode);
            }

        }else{
            blendFactor = 1;
        }
    }

    int TPZCube::NBilinearSides()
    {
        return 27;
    }
    
//    static int permutacoesC [48][27] =
//    {
//        {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26},
//        {1,2,3,0,5,6,7,4,9,10,11,8,13,14,15,12,17,18,19,16,20,22,23,24,21,25,26}
//        
//    };
    
	void TPZCube::LowerDimensionSides(int side,TPZStack<int> &smallsides)
	{
		smallsides.Resize(0);
		int nsidecon = NContainedSides(side);
		for(int is = 0; is < nsidecon - 1; is++)
			smallsides.Push(ContainedSideLocId(side,is));
	}
	
	void TPZCube::LowerDimensionSides(int side,TPZStack<int> &smallsides, int DimTarget)
	{
		smallsides.Resize(0);
		int nsidecon = NContainedSides(side);
		for(int is = 0; is < nsidecon - 1; is++) {
			if (SideDimension(ContainedSideLocId(side,is)) == DimTarget) smallsides.Push(ContainedSideLocId(side,is));
		}
	}
	
	void TPZCube::HigherDimensionSides(int side, TPZStack<int> &high)
	{
		if(side <0 || side >= NSides) {
			PZError << "TPZCube::HigherDimensionSides side "<< side << endl;
		}
		int is;
		for(is=0; is<nhighdimsides[side]; is++) high.Push(highsides[side][is]);
		
	}
	
	int TPZCube::NSideNodes(int side)
	{
		return nsidenodes[side];
	}
	
	int TPZCube::SideNodeLocId(int side, int node)
	{
		if(side<8 && node == 0) return side;
		if(side>=8 && side < 20 && node < 2) return SideNodes[side-8][node];
		if(side>=20 && side < 26 && node < 4) return FaceNodes[side-20][node];
		if(side == 26 && node < 8) return node;
		PZError << "TPZCube::SideNodeLocId inconsistent side or node " << side
		<< ' ' << node << endl;
		return -1;
	}
	
	void TPZCube::CenterPoint(int side, TPZVec<REAL> &center) {
        if (center.size()!=Dimension) {
            DebugStop();
        }
		int i;
		for(i=0; i<Dimension; i++) {
			center[i] = MidSideNode[side][i];
		}
	}
	
	int TPZCube::SideDimension(int side) {
		if(side<0 || side >= NSides) {
			PZError << "TPZCube::SideDimension side " << side << endl;
			return -1;
		}
		return sidedimension[side];
	}
	
	TPZTransform<> TPZCube::SideToSideTransform(int sidefrom, int sideto)
	{
		if(sidefrom <0 || sidefrom >= NSides || sideto <0 || sideto >= NSides) {
			PZError << "TPZCube::HigherDimensionSides sidefrom "<< sidefrom << 
			' ' << sideto << endl;
			return TPZTransform<>(0);
		}
		if(sidefrom == sideto) {
			return TPZTransform<>(sidedimension[sidefrom]);
		}
		if(sidefrom == NSides-1) {
			return TransformElementToSide(sideto);
		}
        if (sideto== NSides -1) {
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
		PZError << "TPZCube::SideToSideTransform highside not found sidefrom "
		<< sidefrom << ' ' << sideto << endl;
		return TPZTransform<>(0);
	}
	
    /** @brief Generates a random point in the master domain */
    void TPZCube::RandomPoint(TPZVec<REAL> &pt)
    {
        for(int i=0; i<3; i++)
        {
            REAL val = -1. + 2.*(REAL) rand() / (RAND_MAX);
            pt[i] = val;
        }
    }

	TPZTransform<> TPZCube::TransformElementToSide(int side){
		
		if(side<0 || side>26){
			PZError << "TPZCube::TransformElementToSide called with side error\n";
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
			case 5 :
			case 6:
			case 7:
				return t;
                
                

			case  8:
			case 16:
				t.Mult()(0,0) =  1.0;
				return t;
                
			case  9:
			case 17:
				t.Mult()(0,1) =  1.0;
				return t;
			case 10:
			case 18:
                
				t.Mult()(0,0) = -1.0;
				return t;
			case 11:
			case 19:
                
				t.Mult()(0,1) = -1.0;
				return t;
			case 12:
			case 13:
			case 14:
			case 15:
                
				t.Mult()(0,2) = 1.0;
				return t;
                
			case 20:
			case 25:
				t.Mult()(0,0) =  1.0;
				t.Mult()(1,1) =  1.0;
				return t;
                
			case 21:
			case 23:
				t.Mult()(0,0) =  1.0;
				t.Mult()(1,2) =  1.0;
				return t;
                
			case 22:
			case 24:
				t.Mult()(0,1) =  1.0;
				t.Mult()(1,2) =  1.0;
				return t;
			case 26:
				t.Mult()(0,0) =  1.0;
				t.Mult()(1,1) =  1.0;
				t.Mult()(2,2) =  1.0;
				return t;
		}
		return TPZTransform<>(0,0);
	}
	
	TPZTransform<> TPZCube::TransformSideToElement(int side){
		
		if(side<0 || side>26){
			PZError << "TPZCube::TransformSideToElement side out range\n";
			return TPZTransform<>(0,0);
		}
		TPZTransform<> t(3,sidedimension[side]);
		t.Mult().Zero();
		t.Sum().Zero();
		
		switch(side){
			case 0:
				t.Sum()(0,0) = -1.0;
				t.Sum()(1,0) = -1.0;
				t.Sum()(2,0) = -1.0;
				return t;
			case 1:
				t.Sum()(0,0) =  1.0;
				t.Sum()(1,0) = -1.0;
				t.Sum()(2,0) = -1.0;
				return t;
			case 2:
				t.Sum()(0,0) =  1.0;
				t.Sum()(1,0) =  1.0;
				t.Sum()(2,0) = -1.0;
				return t;
			case 3:
				t.Sum()(0,0) = -1.0;
				t.Sum()(1,0) =  1.0;
				t.Sum()(2,0) = -1.0;
				return t;
			case 4:
				t.Sum()(0,0) = -1.0;
				t.Sum()(1,0) = -1.0;
				t.Sum()(2,0) =  1.0;
				return t;
			case 5:
				t.Sum()(0,0) =  1.0;
				t.Sum()(1,0) = -1.0;
				t.Sum()(2,0) =  1.0;
				return t;
			case 6:
				t.Sum()(0,0) =  1.0;
				t.Sum()(1,0) =  1.0;
				t.Sum()(2,0) =  1.0;
				return t;
			case 7:
				t.Sum()(0,0) = -1.0;
				t.Sum()(1,0) =  1.0;
				t.Sum()(2,0) =  1.0;
				return t;
			case 8:
				t.Mult()(0,0) =  1.0;
				t.Sum()(1,0)  = -1.0;
				t.Sum()(2,0)  = -1.0;
				return t;
			case 9:
				t.Mult()(1,0) =  1.0;
				t.Sum()(0,0)  =  1.0;
				t.Sum()(2,0)  = -1.0;
				return t;
			case 10:
				t.Mult()(0,0) = -1.0;
				t.Sum()(1,0)  =  1.0;
				t.Sum()(2,0)  = -1.0;
				return t;
			case 11:
				t.Mult()(1,0) = -1.0;
				t.Sum()(0,0)  = -1.0;
				t.Sum()(2,0)  = -1.0;
				return t;
			case 12:
				t.Mult()(2,0) =  1.0;
				t.Sum()(0,0)  = -1.0;
				t.Sum()(1,0)  = -1.0;
				return t;
			case 13:
				t.Mult()(2,0) =  1.0;
				t.Sum()(0,0)  =  1.0;
				t.Sum()(1,0)  = -1.0;
				return t;
			case 14:
				t.Mult()(2,0) =  1.0;
				t.Sum()(0,0)  =  1.0;
				t.Sum()(1,0)  =  1.0;
				return t;
			case 15:
				t.Mult()(2,0) =  1.0;
				t.Sum()(0,0)  = -1.0;
				t.Sum()(1,0)  =  1.0;
				return t;
			case 16:
				t.Mult()(0,0) =  1.0;
				t.Sum()(1,0)  = -1.0;
				t.Sum()(2,0)  =  1.0;
				return t;
			case 17:
				t.Mult()(1,0) =  1.0;
				t.Sum()(0,0)  =  1.0;
				t.Sum()(2,0)  =  1.0;
				return t;
			case 18:
				t.Mult()(0,0) = -1.0;
				t.Sum()(1,0)  =  1.0;
				t.Sum()(2,0)  =  1.0;
				return t;
			case 19:
				t.Mult()(1,0) = -1.0;
				t.Sum()(0,0)  = -1.0;
				t.Sum()(2,0)  =  1.0;
				return t;
			case 20:
				t.Mult()(0,0) =  1.0;
				t.Mult()(1,1) =  1.0;
				t.Sum()(2,0)  = -1.0;
				return t;
			case 21:
				t.Mult()(0,0) =  1.0;
				t.Mult()(2,1) =  1.0;
				t.Sum()(1,0)  = -1.0;
				return t;
			case 22:
				t.Mult()(1,0) =  1.0;
				t.Mult()(2,1) =  1.0;
				t.Sum()(0,0)  =  1.0;
				return t;
			case 23:
				t.Mult()(0,0) =  1.0;
				t.Mult()(2,1) =  1.0;
				t.Sum()(1,0)  =  1.0;
				return t;
			case 24:
				t.Mult()(1,0) =  1.0;
				t.Mult()(2,1) =  1.0;
				t.Sum()(0,0)  = -1.0;
				return t;
			case 25:
				t.Mult()(0,0) =  1.0;
				t.Mult()(1,1) =  1.0;
				t.Sum()(2,0)  =  1.0;
				return t;
			case 26:
				t.Mult()(0,0) =  1.0;
				t.Mult()(1,1) =  1.0;
				t.Mult()(2,2) =  1.0;
				return t;
		}
		return TPZTransform<>(0,0);
	}
	
	
	TPZIntPoints *TPZCube::CreateSideIntegrationRule(int side, int order){
		
		if(side<0 || side>26) {
			PZError << "TPZCube::CreateSideIntegrationRule. bad side number.\n";
			return 0;
		}
		if(side<8)   return new TPZInt1Point(order);            // sides 0 to 7 are vertices (corners)
		if(side<20)  return new TPZInt1d(order);           // sides 8 to 19 are lines
		if(side<26)  {                                     // sides 20 to 25 are quadrilaterals
			return new TPZIntQuad(order,order);
		}
		if(side==26) {
			return new IntruleType(order,order,order);     // integration of the element
		}
		return 0;
		
	}
	
	
	MElementType TPZCube::Type(int side)
	{
		switch(side) {
			case 0:
			case 1:
			case 2:
			case 3:
			case 4:
			case 5:
			case 6:
			case 7:
				return EPoint;
			case 8:
			case 9:
			case 10:
			case 11:
			case 12:
			case 13:
			case 14:
			case 15:
			case 16:
			case 17:
			case 18:
			case 19:        
				return EOned;
			case 20:
			case 21:
			case 22:
			case 23:
			case 24:
			case 25:
				return EQuadrilateral;
			case 26:
				return ECube;
			default:
				return ENoType;
		}
	}
	
	
	int TPZCube::NumSides() {
		return 27;
	}
	
	
	int TPZCube::NContainedSides(int side) {
		if(side<0) return -1;
		if(side<8) return 1;//cantos : 0 a 7
		if(side<20)	return 3;//lados : 8 a 19
		if(side<26)	return 9;//faces : 20 a 25
		if(side==26)	return 27;//centro : 26
		return -1;
	}
	/**
	 return number of sides of dimension dimension
	 **/
	int TPZCube::NumSides(int dimension) {
		if(dimension<0 || dimension> 3) {
			PZError << "TPZCube::NumSides. Bad parameter i.\n";
			return 0;
		}
		if(dimension==0) return 8;
		if(dimension==1) return 12;
		if(dimension==2) return 6;
		if(dimension==3) return 1;
		return -1;
	}
	
	// Pronto 23/04/98
	int TPZCube::ContainedSideLocId(int side, int node) {
		if(side<0 || side>26) return -1;
		if(side<8) {
			if(node==0) return side;
		} 
		else if(side<12) {//8,9,10,11
			int s = side-8;//0,1,2,3
			if(!node) return s;
			if(node==1) return (s+1)%4;
			if(node==2) return side;
		} 
		else if(side<16) {//12,13,14,15
			int s = side-12;//0,1,2,3
			if(!node) return s;
			if(node==1) return s+4;
			if(node==2) return side;
		} 
		else if(side<20) {//16,17,18,19
			int s = side-16;//0,1,2,3
			if(!node) return s+4;
			if(node==1) return (s+1)%4+4;
			if(node==2) return side;
		} 
		else if(side<26) {//20 a 25
			int s = side-20;
			if(node<9) return FaceConnectLocId[s][node];
		} 
		else if(side==26){
			return node;
		}
		PZError << "TPZShapeCube::ContainedSideLocId called for node = "
		<< node << " and side = " << side << "\n";
		return -1;
	}
	
	bool TPZCube::IsInParametricDomain(const TPZVec<REAL> &pt, REAL tol){
		const REAL qsi = pt[0];
		const REAL eta = pt[1];
		const REAL zeta = pt[2];
		if( ( fabs(qsi) <= 1. + tol ) && ( fabs(eta) <= 1. + tol ) && ( fabs(zeta) <= 1. + tol ) ){
			return true;
		}
		else{
			return false;
		}  
	}//method

    template<class T>
    bool TPZCube::CheckProjectionForSingularity(const int &side, const TPZVec<T> &xiInterior) {
        return true;
    }

    template<class T>
    void TPZCube::MapToSide(int side, TPZVec<T> &InternalPar, TPZVec<T> &SidePar, TPZFMatrix<T> &JacToSide) {
		TPZTransform<> TransfR = pztopology::TPZCube::SideToSideTransform(NSides - 1, side);
        TPZTransform<T> Transf;
        Transf.CopyFrom(TransfR);
		SidePar.Resize(SideDimension(side));
		Transf.Apply(InternalPar,SidePar);
		
		int R = Transf.Mult().Rows();
		int C = Transf.Mult().Cols();
		
		JacToSide.Resize(R,C);
		for(int i = 0; i < R; i++)
		{
			for(int j = 0; j < C; j++) JacToSide(i,j) = Transf.Mult()(i,j);
		}
	}
    
    void TPZCube::ParametricDomainNodeCoord(int node, TPZVec<REAL> &nodeCoord)
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
                nodeCoord[2] = -1.;
                break;
            }
            case (1):
            {
                nodeCoord[0] =  1.;
                nodeCoord[1] = -1.;
                nodeCoord[2] = -1.;
                break;
            }
            case (2):
            {
                nodeCoord[0] =  1.;
                nodeCoord[1] =  1.;
                nodeCoord[2] = -1.;
                break;
            }
            case (3):
            {
                nodeCoord[0] = -1.;
                nodeCoord[1] =  1.;
                nodeCoord[2] = -1.;
                break;
            }
            case (4):
            {
                nodeCoord[0] = -1.;
                nodeCoord[1] = -1.;
                nodeCoord[2] =  1.;
                break;
            }
            case (5):
            {
                nodeCoord[0] =  1.;
                nodeCoord[1] = -1.;
                nodeCoord[2] =  1.;
                break;
            }
            case (6):
            {
                nodeCoord[0] = 1.;
                nodeCoord[1] = 1.;
                nodeCoord[2] = 1.;
                break;
            }
            case (7):
            {
                nodeCoord[0] = -1.;
                nodeCoord[1] =  1.;
                nodeCoord[2] =  1.;
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
	int TPZCube::GetTransformId(const TPZVec<int64_t> &id)
	{
		LOGPZ_ERROR(logger,"GetTransformId not implemented")
		return -1;
	}
	/**
	 * Method which identifies the transformation of a side based on the IDs
	 * of the corner nodes
	 * @param id indexes of the corner nodes
	 * @return index of the transformation of the point corresponding to the topology
	 */	
	int TPZCube::GetTransformId(const int side, const TPZVec<int64_t> &id)
	{
		switch (side) {
			case 0:
			case 1:
			case 2:
			case 3:
			case 4:
			case 5:
			case 6:
			case 7:
				return 0;
				break;
			case 8:
			case 9:
			case 10:
			case 11:
			case 12:
			case 13:
			case 14:
			case 15:
			case 16:
			case 17:
			case 18:
			case 19:
			{
				int in1 = ContainedSideLocId(side,0);
				int in2 = ContainedSideLocId(side,1);
				return id[in1]<id[in2] ? 0 : 1;
			}
				break;
			case 20:
			case 21:
			case 22:
			case 23:
			case 24:
			case 25:
			{
				TPZManVector<int64_t,4> locid(4);
				int i;
				for(i=0; i<4; i++) locid[i] = id[ContainedSideLocId(side,i)];
				return pztopology::TPZQuadrilateral::GetTransformId(locid);
				//			return pzshape::TPZShapeQuad::GetTransformId2dQ(locid);
			}
				break;			
			case 26:
				return 0;//that is not really true
			default:
				break;
		}
		return -1;
	}
	
	/**
	 * Identifies the permutation of the nodes needed to make neighbouring elements compatible 
	 * in terms of order of shape functions
	 * @param side : side for which the permutation is needed
	 * @param id : ids of the corner nodes of the elements
	 * @param permgather : permutation vector in a gather order
	 */
	void TPZCube::GetSideHDivPermutation(int transformationid, TPZVec<int> &permgather)
	{
        DebugStop();
	}
    
    void computedirectionsC(int inicio, int fim, TPZFMatrix<REAL> &bvec, TPZFMatrix<REAL> &t1vec,
                           TPZFMatrix<REAL> &t2vec, TPZFMatrix<REAL> &gradx, TPZFMatrix<REAL> &directions);
    
    void computedirectionsC(int inicio, int fim, TPZFMatrix<REAL> &bvec, TPZFMatrix<REAL> &t1vec,
                           TPZFMatrix<REAL> &t2vec, TPZFMatrix<REAL> &gradx, TPZFMatrix<REAL> &directions)
    {
        TPZVec<REAL> u(3);
        TPZVec<REAL> v(3);
        TPZVec<REAL> uxv(3);// result
        int cont = 0;
        
        for (int ivet=inicio; ivet<=fim; ivet++)
        { 
            if(inicio < 54)
            {
                for (int ilin=0; ilin<3; ilin++)
                {
                    u[ilin] = t1vec(ilin,ivet);
                    v[ilin] = t2vec(ilin,ivet);
                }
                TPZVec<REAL> e2(3);
                TPZNumeric::ProdVetorial(u,v,e2);
                //            e2[0] = u[1]*v[2]-u[2]*v[1];
                //            e2[1] = -(u[0]*v[2]-v[0]*u[2]);
                //            e2[2] = u[0]*v[1]-v[0]*u[1];
                
                // calc do v gradx*b
                TPZManVector<REAL,3> dxt1(3,0.),dxt2(3,0.),dxt3(3,0.),Vvec(3,0.);
                for (int il=0; il<3; il++)
                {
                    for (int i = 0 ; i<3; i++)
                    {
                        dxt1[il] += gradx(il,i) * t1vec(i,ivet);
                        dxt2[il] += gradx(il,i) * t2vec(i,ivet);
                        //                    dxt3[il] += gradx(il,i) * e2[i];
                        Vvec[il] += gradx(il,i) * bvec(i,ivet);
                    }
                    //be2 += bvec(il,ivet)*e2[il];
                }
                REAL normaX0xX1 = 0.0;
                TPZManVector<REAL,3> normal(3,0.);
                TPZNumeric::ProdVetorial(dxt1,dxt2,normal);
                for (int pos=0; pos<3; pos++)
                {
                    normaX0xX1 += normal[pos]*normal[pos]; //uxv[pos]*uxv[pos];
                }
                
                TPZFMatrix<REAL> Wvec(3,1);
                
                REAL detgrad = gradx(0,0)*gradx(1,1)*gradx(2,2) + gradx(0,1)*gradx(1,2)*gradx(2,0) + gradx(0,2)*gradx(1,0)*gradx(2,1) - gradx(0,2)*gradx(1,1)*gradx(2,0) - gradx(0,0)*gradx(1,2)*gradx(2,1) - gradx(0,1)*gradx(1,0)*gradx(2,2);
                
                normaX0xX1 = sqrt(normaX0xX1);
                
                if (detgrad<0)
                {
                    DebugStop();
                }
                
                for (int il=0; il<3; il++)
                {
                    Wvec(il,0) = Vvec[il]*normaX0xX1/(detgrad);
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
    
//    void computedirectionsC(int inicio, int fim, TPZFMatrix<REAL> &bvec, TPZFMatrix<REAL> &t1vec,
//                            TPZFMatrix<REAL> &t2vec, TPZFMatrix<REAL> &gradx, TPZFMatrix<REAL> &directions)
//    {
//        REAL detgrad = 0.0;
//        TPZVec<REAL> u(3);
//        TPZVec<REAL> v(3);
//        TPZVec<REAL> uxv(3);// result
//        int cont = 0;
//        
//        for (int ivet=inicio; ivet<=fim; ivet++)
//        {
//            for (int ilin=0; ilin<3; ilin++)
//            {
//                u[ilin] = t1vec(ilin,ivet);
//                v[ilin] = t2vec(ilin,ivet);
//            }
//            TPZVec<REAL> e2(3);
//            detgrad = 0.0;
//            REAL normaX0xX1 = 0.0;
//            //TPZNumeric::ProdVetorial(u,v,e2);
//            e2[0] = u[1]*v[2]-u[2]*v[1];
//            e2[1] = -(u[0]*v[2]-v[0]*u[2]);
//            e2[2] = u[0]*v[1]-v[0]*u[1];
//            
//            // calc do v gradx*b
//            TPZManVector<REAL,3> dxt1(3,0.),dxt2(3,0.),dxt3(3,0.),Vvec(3,0.);
//            REAL be2 = 0.0, ne2 = 0.0;
//            for(int i=0;i<3;i++)
//            {
//                ne2 += e2[i]*e2[i];
//            }
//            ne2 = sqrt(fabs(ne2));
//            for (int il=0; il<3; il++)
//            {
//                for (int i = 0 ; i<3; i++)
//                {
//                    dxt1[il] += gradx(il,i) * t1vec(i,ivet);
//                    dxt2[il] += gradx(il,i) * t2vec(i,ivet);
//                    dxt3[il] += gradx(il,i) * e2[i]/ne2;
//                    Vvec[il] += gradx(il,i) * bvec(i,ivet);
//                }
//                be2 += bvec(il,ivet)*e2[il]/ne2;
//            }
//            TPZManVector<REAL,3> normal(3,0.);
//            //TPZNumeric::ProdVetorial(dxt1,dxt2,normal);
//            normal[0] = dxt1[1]*dxt2[2]-dxt1[2]*dxt2[1];
//            normal[1] = -(dxt1[0]*dxt2[2]-dxt2[0]*dxt1[2]);
//            normal[2] = dxt1[0]*dxt2[1]-dxt2[0]*dxt1[1];
//            
//            for (int pos=0; pos<3; pos++)
//            {
//                detgrad += normal[pos]*dxt3[pos];//uxv[pos]*gradx.GetVal(pos, 2);
//                normaX0xX1 += normal[pos]*normal[pos]; //uxv[pos]*uxv[pos];
//            }
//            TPZFMatrix<REAL> Wvec(3,1);
//            detgrad = fabs(detgrad);
//            normaX0xX1 = sqrt(normaX0xX1);
//            
//            for (int il=0; il<3; il++)
//            {
//                Wvec(il,0) = Vvec[il]*normaX0xX1/(detgrad*be2);
//                directions(il,cont) = Wvec(il,0);
//            }
//            cont++;
//        }
//        
//    }

    
//    static REAL bCubo[81][3] = // direcao perpendicular ao lado
//    {
//        {0,0,-1}, {0,0,-1}, {0,0,-1}, {0,0,-1}, {0,0,-1}, {0,0,-1}, {0,0,-1}, {0,0,-1}, {0,0,-1},// face 0
//        {0,-1,0}, {0,-1,0}, {0,-1,0}, {0,-1,0}, {0,-1,0}, {0,-1,0}, {0,-1,0}, {0,-1,0}, {0,-1,0},// face 1
//        {1,0,0} , {1,0,0} , {1,0,0} , {1,0,0} , {1,0,0} , {1,0,0} , {1,0,0} , {1,0,0} , {1,0,0} ,// face 2
//        {0,1,0} , {0,1,0} , {0,1,0} , {0,1,0} , {0,1,0} , {0,1,0} , {0,1,0} , {0,1,0} , {0,1,0} ,// face 3
//        {-1,0,0}, {-1,0,0}, {-1,0,0}, {-1,0,0}, {-1,0,0}, {-1,0,0}, {-1,0,0}, {-1,0,0}, {-1,0,0},// face 4
//        {0,0,1} , {0,0,1} , {0,0,1} , {0,0,1} , {0,0,1} , {0,0,1} , {0,0,1} , {0,0,1} , {0,0,1}, // face 5
//        //interiores
//        //arestas
//        //{1,0,0},{0,1,0},{-1,0,0},{0,-1,0},  {0,0,1},{0,0,1},{0,0,1},{0,0,1},  {1,0,0},{0,1,0},{-1,0,0},{0,-1,0},
//        {1,0,0},{0,1,0},{-1,0,0},{0,-1,0},  {0,0,1},{0,0,1},{0,0,1},{0,0,1},  {1,0,0},{0,1,0},{-1,0,0},{0,-1,0},
//        //faces
//        {1,0,0}, {0,1,0}, // tang da face 0
//        {1,0,0}, {0,0,1},  // tang da face 1
//        {0,1,0}, {0,0,1}, // tang da face 2
//        {1,0,0}, {0,0,1}, // tang da face 3
//        {0,1,0}, {0,0,1}, // tang da face 4
//        {1,0,0}, {0,1,0},  // tang da face 5
//        {1,0,0}, // volume
//        {0,1,0}, // volume
//        {0,0,1}  // volume
//    };

    template <class TVar>
    void TPZCube::ComputeHDivDirections(TPZFMatrix<TVar> &gradx, TPZFMatrix<TVar> &directions)
    {
        TVar detjac = TPZAxesTools<TVar>::ComputeDetjac(gradx);
        
        TPZManVector<TVar,3> v1(3),v2(3),v3(3),v1v2(3),v3v1(3),v2v3(3),vec1(3),vec2(3),vec3(3);
        for (int i=0; i<3; i++) {
            v1[i] = gradx(i,0);
            v2[i] = gradx(i,1);
            v3[i] = gradx(i,2);
        }
        
        TPZNumeric::ProdVetorial(v1,v2,v1v2);
        TPZNumeric::ProdVetorial(v2,v3,v2v3);
        TPZNumeric::ProdVetorial(v3,v1,v3v1);
        
        TVar Nv1v2 = TPZNumeric::Norm(v1v2);
        TVar Nv2v3 = TPZNumeric::Norm(v2v3);
        TVar Nv3v1 = TPZNumeric::Norm(v3v1);
        
        /**
         * @file
         * @brief Computing mapped vector with scaling factor equal 1.0.
         * using contravariant piola mapping.
         */
        TPZManVector<TVar,3> NormalScales(3,1.);
        
        
        {
            for (int i=0; i<3; i++) {
                v1[i] *= 1./detjac;
                v2[i] *= 1./detjac;
                v3[i] *= 1./detjac;
            }
            
        }
        
        for (int i=0; i<3; i++) {
            for (int iv=0; iv<9; iv++) {
                directions(i,iv)        = -v3[i]*NormalScales[1];
                directions(i,iv+9)      = -v2[i]*NormalScales[2];
                directions(i,iv+18)     = v1[i]*NormalScales[1];
                directions(i,iv+27)     = v2[i]*NormalScales[2];
                directions(i,iv+36)     = -v1[i]*NormalScales[1];
                directions(i,iv+45)     = v3[i]*NormalScales[1];
            }
            directions(i,54)    = v1[i]*NormalScales[1];
            directions(i,55)    = v2[i]*NormalScales[2];
            directions(i,56)    = -v1[i]*NormalScales[1];
            directions(i,57)    = -v2[i]*NormalScales[2];

            directions(i,58)    = v3[i]*NormalScales[1];
            directions(i,59)    = v3[i]*NormalScales[1];
            directions(i,60)    = v3[i]*NormalScales[1];
            directions(i,61)    = v3[i]*NormalScales[1];

            directions(i,62)    = v1[i]*NormalScales[1];
            directions(i,63)    = v2[i]*NormalScales[2];
            directions(i,64)    = -v1[i]*NormalScales[1];
            directions(i,65)    = -v2[i]*NormalScales[2];
            
            directions(i,66)    = v1[i]*NormalScales[1];
            directions(i,67)    = v2[i]*NormalScales[2];
            
            directions(i,68)    = v1[i]*NormalScales[1];
            directions(i,69)    = v3[i]*NormalScales[1];
            directions(i,70)    = v2[i]*NormalScales[2];
            directions(i,71)    = v3[i]*NormalScales[1];
            directions(i,72)    = v1[i]*NormalScales[1];
            directions(i,73)    = v3[i]*NormalScales[1];
            directions(i,74)    = v2[i]*NormalScales[2];
            directions(i,75)    = v3[i]*NormalScales[1];
            
            directions(i,76)    = v1[i]*NormalScales[1];
            directions(i,77)    = v2[i]*NormalScales[2];

            directions(i,78)    = v1[i]*NormalScales[1];
            directions(i,79)    = v2[i]*NormalScales[2];
            directions(i,80)    = v3[i]*NormalScales[1];
            
        }
        

    }

    /// Compute the directions of the HDiv vectors
    // template <class TVar>
    void TPZCube::ComputeConstantHDiv(TPZVec<REAL> &point, TPZFMatrix<REAL> &RT0function, TPZVec<REAL> &div)
    {

        REAL scale = 4.;    
        REAL qsi = point[0];
        REAL eta = point[1];
        REAL zeta = point[2];

        //Face functions
        //For each face function: compute div = \nabla \cdot RT0function = d_RT0/d_qsi + d_RT0/d_eta 
        RT0function(2,5) = 0.5 * (1. + zeta) / scale;
        div[5] = 0.5 / scale;
        RT0function(1,3) = 0.5 * (1. + eta) / scale;
        div[3] = 0.5 / scale;
        RT0function(0,2) = 0.5 * (1. + qsi) / scale;
        div[2] = 0.5 / scale;

        RT0function(0,4) = -0.5 * (1. - qsi) / scale;
        div[4] = 0.5 / scale;
        RT0function(1,1) = -0.5 * (1. - eta) / scale;
        div[1] = 0.5 / scale;
        RT0function(2,0) = -0.5 * (1. - zeta) / scale;
        div[0] = 0.5 / scale;


    }

    /// Compute the directions of the HDiv vectors
    // template <class TVar>
    void TPZCube::ComputeConstantHCurl(TPZVec<REAL> &point, TPZFMatrix<REAL> &N0function, TPZFMatrix<REAL> &curl, const TPZVec<int> &transformationIds)
    {
        REAL scale = 2.;    
        REAL qsi = point[0];
        REAL eta = point[1];
        REAL zeta = point[2];

        constexpr auto nEdges{12};
        TPZManVector<REAL,nEdges> edgeSign(nEdges,0);
        for(auto iEdge = 0; iEdge < nEdges; iEdge++){
            edgeSign[iEdge] = transformationIds[iEdge] == 0 ? 1 : -1;
        }

        //First type Nedelec functions
        //X direction
        //Edge 16
        N0function(0,8) = 0.25 * (1. - eta) * (1. + zeta) / scale * edgeSign[8];
        curl(0,8) = 0.;
        curl(1,8) = 0.25 * (1. - eta) / scale * edgeSign[8];
        curl(2,8) = 0.25 * (1. + zeta) / scale * edgeSign[8];
        //Edge 18
        N0function(0,10) = -0.25 * (1. + eta) * (1. + zeta) / scale * edgeSign[10];
        curl(0,10) = 0.;
        curl(1,10) = -0.25 * (1. + eta) / scale * edgeSign[10];
        curl(2,10) =  0.25 * (1. + zeta) / scale * edgeSign[10];
        //Edge 8
        N0function(0,0) = 0.25 * (1. - eta) * (1. - zeta) / scale * edgeSign[0];
        curl(0,0) = 0.;
        curl(1,0) = -0.25 * (1. - eta) / scale * edgeSign[0];
        curl(2,0) =  0.25 * (1. - zeta) / scale * edgeSign[0];
        //Edge 10
        N0function(0,2) = -0.25 * (1. + eta) * (1. - zeta) / scale * edgeSign[2];
        curl(0,2) = 0.;
        curl(1,2) = 0.25 * (1. + eta) / scale * edgeSign[2];
        curl(2,2) = 0.25 * (1. - zeta) / scale * edgeSign[2];

        //Y direction
        //Edge 17
        N0function(1,9) = 0.25 * (1. + qsi) * (1. + zeta) / scale * edgeSign[9];
        curl(0,9) = -0.25 * (1. + qsi) / scale * edgeSign[9];
        curl(1,9) = 0.;
        curl(2,9) =  0.25 * (1. + zeta) / scale * edgeSign[9];
        //Edge 19
        N0function(1,11) = -0.25 * (1. - qsi) * (1. + zeta) / scale * edgeSign[11];
        curl(0,11) = 0.25 * (1. - qsi) / scale * edgeSign[11];
        curl(1,11) = 0.;
        curl(2,11) = 0.25 * (1. + zeta) / scale * edgeSign[11];
        //Edge 9
        N0function(1,1) = 0.25 * (1. + qsi) * (1. - zeta) / scale * edgeSign[1];
        curl(0,1) = 0.25 * (1. + qsi) / scale * edgeSign[1];
        curl(1,1) = 0.;
        curl(2,1) = 0.25 * (1. - zeta) / scale * edgeSign[1];
        //Edge 11
        N0function(1,3) = -0.25 * (1. - qsi) * (1. - zeta) / scale * edgeSign[3];
        curl(0,3) = -0.25 * (1. - qsi) / scale * edgeSign[3];
        curl(1,3) = 0.;
        curl(2,3) = 0.25 * (1. - zeta) / scale * edgeSign[3];
               
        //Z direction
        //Edge 12
        N0function(2,4) = 0.25 * (1. - qsi) * (1. - eta) / scale * edgeSign[4];
        curl(0,4) = -0.25 * (1. - qsi) / scale * edgeSign[4];
        curl(1,4) =  0.25 * (1. - eta) / scale * edgeSign[4];
        curl(2,4) = 0.;
        //Edge 13
        N0function(2,5) = 0.25 * (1. + qsi) * (1. - eta) / scale * edgeSign[5];
        curl(0,5) = -0.25 * (1. + qsi) / scale * edgeSign[5];
        curl(1,5) = -0.25 * (1. - eta) / scale * edgeSign[5];
        curl(2,5) = 0.;
        //Edge 14
        N0function(2,6) = 0.25 * (1. + qsi) * (1. + eta) / scale * edgeSign[6];
        curl(0,6) =  0.25 * (1. + qsi) / scale * edgeSign[6];
        curl(1,6) = -0.25 * (1. + eta) / scale * edgeSign[6];
        curl(2,6) = 0.;
        //Edge 15
        N0function(2,7) = 0.25 * (1. - qsi) * (1. + eta) / scale * edgeSign[7];
        curl(0,7) = 0.25 * (1. - qsi) / scale * edgeSign[7];
        curl(1,7) = 0.25 * (1. + eta) / scale * edgeSign[7];
        curl(2,7) = 0.;

    }

    // Get face orientation
    int TPZCube::GetSideOrient(const int &face){
        return fSideOrient[face];
    }

    void TPZCube::ComputeDirections(int side, TPZFMatrix<REAL> &gradx, TPZFMatrix<REAL> &directions, TPZVec<int> &sidevectors)
    {
        if(gradx.Cols()!=3)
        { std::cout << "Gradient dimensions are not compatible with this topology" << std::endl;
            DebugStop();
        }
        TPZFMatrix<REAL> bvec(3,81);
        int numvec = bvec.Cols();
        TPZFMatrix<REAL> t1vec(3,numvec);
        TPZFMatrix<REAL> t2vec(3,numvec);
        
        directions.Redim(3, numvec);
        for (int lin = 0; lin<numvec; lin++)
        {
            for(int col = 0;col<3;col++)
            {
                bvec.PutVal(col, lin, bCubo[lin][col]);
                t1vec.PutVal(col, lin, t1Cubo[lin][col]);
                t2vec.PutVal(col, lin, t2Cubo[lin][col]);
            }
        }
        
        // calcula os vetores
       
        switch (side) {
            case 20:
            {
                directions.Resize(3, 9);
                sidevectors.Resize(9);
                int inicio = 0, fim = 8;
                computedirectionsC( inicio, fim, bvec, t1vec, t2vec, gradx, directions);
                for (int ip = 0; ip < 9; ip++) {
                    sidevectors[ip] = vectorsideorderC[ip+inicio];
                }
            }
                break;
            case 21:
            {
                directions.Resize(3, 9);
                sidevectors.Resize(9);
                int inicio = 9, fim = 17;
                computedirectionsC( inicio, fim, bvec, t1vec, t2vec, gradx, directions);
                for (int ip = 0; ip < 9; ip++) {
                    sidevectors[ip] = vectorsideorderC[ip+inicio];
                }
            }
                break;
            case 22:
            {
                directions.Resize(3, 9);
                sidevectors.Resize(9);
                int inicio = 18, fim = 26;
                computedirectionsC( inicio, fim, bvec, t1vec, t2vec, gradx, directions);
                for (int ip = 0; ip < 9; ip++) {
                    sidevectors[ip] = vectorsideorderC[ip+inicio];
                }
            }
                break;
            case 23:
            {
                directions.Resize(3, 9);
                sidevectors.Resize(9);
                int inicio = 27, fim = 35;
                computedirectionsC( inicio, fim, bvec, t1vec, t2vec, gradx, directions);
                for (int ip = 0; ip < 9; ip++) {
                    sidevectors[ip] = vectorsideorderC[ip+inicio];
                }
            }
                break;
            case 24:
            {
                directions.Resize(3, 9);
                sidevectors.Resize(9);
                int inicio = 36, fim = 44;
                computedirectionsC( inicio, fim, bvec, t1vec, t2vec, gradx, directions);
                for (int ip = 0; ip < 9; ip++) {
                    sidevectors[ip] = vectorsideorderC[ip+inicio];
                }
            }
                break;
            case 25:
            {
                directions.Resize(3, 9);
                sidevectors.Resize(9);
                int inicio = 45, fim = 53;
                computedirectionsC( inicio, fim, bvec, t1vec, t2vec, gradx, directions);
                for (int ip = 0; ip < 9; ip++) {
                    sidevectors[ip] = vectorsideorderC[ip+inicio];
                }
            }
                break;
            case 26:
            {
                directions.Resize(3, 27);
                sidevectors.Resize(27);
                int inicio = 54, fim = 80;
                computedirectionsC( inicio, fim, bvec, t1vec, t2vec, gradx, directions);
                for (int ip = 0; ip < 27; ip++) {
                    sidevectors[ip] = vectorsideorderC[ip+inicio];
                }
            }
                break;

                
            default:
                break;
        }
#ifdef PZDEBUG
        if (SideDimension(side) == 2) {
            TPZStack<int> lowerdim;
            LowerDimensionSides(side, lowerdim);
            lowerdim.Push(side);
            if (sidevectors.size() != lowerdim.size()) {
                DebugStop();
            }
            int nwrong = 0;
            for (int i=0; i<lowerdim.size(); i++) {
                if (lowerdim[i] != sidevectors[i]) {
                    nwrong++;
                }
            }
            if (nwrong)
            {
                std::cout << "sidevectors = " << sidevectors << " lowerdim = " << lowerdim << std::endl;
                DebugStop();
            }
        }
#endif
        
	}
    
    void TPZCube::GetSideHDivDirections(TPZVec<int> &sides, TPZVec<int> &dir, TPZVec<int> &bilounao)
    {
        int nsides = NumSides()*3;
        
        sides.Resize(nsides);
        dir.Resize(nsides);
        bilounao.Resize(nsides);
        
        for (int is = 0; is<nsides; is++)
        {
            sides[is] = vectorsideorderC[is];
            dir[is] = direcaoksioueta[is];
            bilounao[is] = bilinearounao[is];
        }
    }

    void TPZCube::GetSideHDivDirections(TPZVec<int> &sides, TPZVec<int> &dir, TPZVec<int> &bilounao, TPZVec<int> &sidevectors)
    {
        int nsides = NumSides()*3;
        
        sides.Resize(nsides);
        dir.Resize(nsides);
        bilounao.Resize(nsides);
        
        for (int is = 0; is<nsides; is++)
        {
            sides[is] = vectorsideorderC[is];
            dir[is] = direcaoksioueta[is];
            bilounao[is] = bilinearounao[is];
        }
        
        for (int i=0; i<Dimension*NumSides(); i++) {
            sidevectors[i] = vectorsideorderC[i];
        }
    }

    template <class TVar>
    void TPZCube::ComputeHCurlDirections(TPZFMatrix<TVar> &gradx, TPZFMatrix<TVar> &directions, const TPZVec<int> &transformationIds)
    {
        TPZManVector<TVar,3> v1(3),v2(3),v3(3);

        for (int i=0; i<3; i++) {
            v1[i] = gradx(i,0);
            v2[i] = gradx(i,1);
            v3[i] = gradx(i,2);
        }
        constexpr int nEdges = 12;
        constexpr REAL edgeLength[nEdges]{2,2,2,2,2,2,2,2,2,2,2,2};
        constexpr int nFaces = 6;
        constexpr REAL faceArea[nFaces]{4,4,4,4,4,4};
        TPZManVector<REAL,nEdges> edgeSign(nEdges,0);
        TPZManVector<REAL,nFaces>  faceOrient(nFaces,0);

        for(auto iEdge = 0; iEdge < nEdges; iEdge++){
            edgeSign[iEdge] = transformationIds[iEdge] == 0 ? 1 : -1;
        }
        for(auto iFace = 0; iFace < nFaces; iFace++){
            faceOrient[iFace] = transformationIds[nEdges + iFace] % 2 == 0 ? 1 : -1;
        }

        for(int iSide = 0; iSide < nEdges; iSide ++){
            int sign = (iSide < 4 || iSide > 7) ?
                    ( (iSide%4) /2 ? -1 : 1)
                    :
                    1;// sign will be : 1 1 -1 -1     1 1 1 1     1 1 -1 -1
            sign *= edgeSign[iSide];
            TPZVec<TVar>& vec1 = (iSide < 4 || iSide > 7) ?
                    ( (iSide%4) % 2 ? v2 : v1)
                    :
                    v3;// vec1 will be : v1 v2 v1 v2     v3 v3 v3 v3     v1 v2 v1 v2
            for (int i=0; i<3; i++){
                //v^{e,a} constant vector fields associated with edge e and vertex a
                //they are defined in such a way that v^{e,a} is normal to the edge \hat{e}
                //adjacent to edge e by the vertex a. the tangential component is set to be 1 /edgeLength[e] = 0.5
                directions(i,iSide * 2) =
                directions(i,iSide * 2 + 1) =
                        //v^{e,T} constant vector fields associated with edge e and aligned with it
                directions(i,nEdges*2 + iSide) = sign * vec1[i] / edgeLength[iSide];
            }
        }
        for (int i=0; i<3; i++) {
            //v^{F,e} constant vector fields associated with face F and edge e
            //they are defined in such a way that v^{F,e} is normal to the face \hat{F}
            //adjacent to face F by edge e
            directions(i, 36) =  v2[i] * faceOrient[0] * edgeSign[ 8-NCornerNodes] / faceArea[0];//face 20 edge 8
            directions(i, 37) = -v1[i] * faceOrient[0] * edgeSign[ 9-NCornerNodes] / faceArea[0];//face 20 edge 9
            directions(i, 38) = -v2[i] * faceOrient[0] * edgeSign[10-NCornerNodes] / faceArea[0];//face 20 edge 10
            directions(i, 39) =  v1[i] * faceOrient[0] * edgeSign[11-NCornerNodes] / faceArea[0];//face 20 edge 11

            directions(i, 40) =  v3[i] * faceOrient[1] * edgeSign[ 8-NCornerNodes] / faceArea[1];//face 21 edge 8
            directions(i, 41) = -v1[i] * faceOrient[1] * edgeSign[13-NCornerNodes] / faceArea[1];//face 21 edge 13
            directions(i, 42) =  v3[i] * faceOrient[1] * edgeSign[16-NCornerNodes] / faceArea[1];//face 21 edge 16
            directions(i, 43) = -v1[i] * faceOrient[1] * edgeSign[12-NCornerNodes] / faceArea[1];//face 21 edge 12

            directions(i, 44) =  v3[i] * faceOrient[2] * edgeSign[ 9-NCornerNodes] / faceArea[2];//face 22 edge 9
            directions(i, 45) = -v2[i] * faceOrient[2] * edgeSign[14-NCornerNodes] / faceArea[2];//face 22 edge 14
            directions(i, 46) =  v3[i] * faceOrient[2] * edgeSign[17-NCornerNodes] / faceArea[2];//face 22 edge 17
            directions(i, 47) = -v2[i] * faceOrient[2] * edgeSign[13-NCornerNodes] / faceArea[2];//face 22 edge 13

            directions(i, 48) = -v3[i] * faceOrient[3] * edgeSign[10-NCornerNodes] / faceArea[3];//face 23 edge 10
            directions(i, 49) = -v1[i] * faceOrient[3] * edgeSign[14-NCornerNodes] / faceArea[3];//face 23 edge 14
            directions(i, 50) = -v3[i] * faceOrient[3] * edgeSign[18-NCornerNodes] / faceArea[3];//face 23 edge 18
            directions(i, 51) = -v1[i] * faceOrient[3] * edgeSign[15-NCornerNodes] / faceArea[3];//face 23 edge 15

            directions(i, 52) = -v3[i] * faceOrient[4] * edgeSign[11-NCornerNodes] / faceArea[4];//face 24 edge 11
            directions(i, 53) = -v2[i] * faceOrient[4] * edgeSign[15-NCornerNodes] / faceArea[4];//face 24 edge 15
            directions(i, 54) = -v3[i] * faceOrient[4] * edgeSign[19-NCornerNodes] / faceArea[4];//face 24 edge 19
            directions(i, 55) = -v2[i] * faceOrient[4] * edgeSign[12-NCornerNodes] / faceArea[4];//face 24 edge 12

            directions(i, 56) =  v2[i] * faceOrient[5] * edgeSign[16-NCornerNodes] / faceArea[5];//face 25 edge 16
            directions(i, 57) = -v1[i] * faceOrient[5] * edgeSign[17-NCornerNodes] / faceArea[5];//face 25 edge 17
            directions(i, 58) = -v2[i] * faceOrient[5] * edgeSign[18-NCornerNodes] / faceArea[5];//face 25 edge 18
            directions(i, 59) =  v1[i] * faceOrient[5] * edgeSign[19-NCornerNodes] / faceArea[5];//face 25 edge 19

            //v^{F,T} are calculated afterwards

            //v^{F,orth} vector associated with face F and normal to it
            directions(i, 72) = -v3[i];//face 20
            directions(i, 73) = -v2[i];//face 21
            directions(i, 74) = v1[i];//face 22
            directions(i, 75) = v2[i];//face 23
            directions(i, 76) = -v1[i];//face 24
            directions(i, 77) = v3[i];//face 25

            //v^{K,3}
            directions(i, 78) = v1[i];
            directions(i, 79) = v2[i];
            directions(i, 80) = v3[i];
        }
        TPZManVector<REAL,2> vft1(2,0), vft2(2,0);
        constexpr auto firstVftVec = 60;
        //v^{F,T} orthonormal vectors associated with face F and tangent to it.
        for(auto iFace = 0; iFace < nFaces; iFace ++){
            TPZQuadrilateral::ComputeHCurlFaceDirections(vft1,vft2,transformationIds[nEdges + iFace]);
            directions(0,firstVftVec+2*iFace) = 0;directions(1,firstVftVec+2*iFace) = 0;directions(2,firstVftVec+2*iFace) = 0;
            directions(0,firstVftVec+2*iFace+1) = 0;directions(1,firstVftVec+2*iFace+1) = 0;directions(2,firstVftVec+2*iFace+1) = 0;
            auto axes = TPZCube::TransformElementToSide(NCornerNodes+nEdges+iFace).Mult();
            axes.Transpose();
            for(auto x = 0; x < Dimension; x++){
                for(auto i = 0; i < 2; i++) {
                    directions(x, firstVftVec + 2 * iFace) += axes(x,i) * vft1[i];
                    directions(x, firstVftVec + 2 * iFace + 1) += axes(x,i) * vft2[i];
                }
            }
        }
    }

    int TPZCube::ClassId() const{
        return Hash("TPZCube");
    }

    void TPZCube::Read(TPZStream& buf, void* context) {

    }

    void TPZCube::Write(TPZStream& buf, int withclassid) const {

    }

}

/**********************************************************************************************************************
 * The following are explicit instantiation of member function template of this class, both with class T=REAL and its
 * respective FAD<REAL> version. In other to avoid potential errors, always declare the instantiation in the same order
 * in BOTH cases.    @orlandini
 **********************************************************************************************************************/
template bool pztopology::TPZCube::CheckProjectionForSingularity<REAL>(const int &side, const TPZVec<REAL> &xiInterior);

template void pztopology::TPZCube::MapToSide<REAL>(int side, TPZVec<REAL> &InternalPar, TPZVec<REAL> &SidePar, TPZFMatrix<REAL> &JacToSide);

template void pztopology::TPZCube::BlendFactorForSide<REAL>(const int &, const TPZVec<REAL> &, REAL &, TPZVec<REAL> &);

template void pztopology::TPZCube::TShape<REAL>(const TPZVec<REAL> &loc,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);

template void pztopology::TPZCube::ComputeHDivDirections<REAL>(TPZFMatrix<REAL> &gradx, TPZFMatrix<REAL> &directions);

template void pztopology::TPZCube::ComputeHCurlDirections<REAL>(TPZFMatrix<REAL> &gradx, TPZFMatrix<REAL> &directions, const TPZVec<int> &transformationIds);

template bool pztopology::TPZCube::CheckProjectionForSingularity<Fad<REAL>>(const int &side, const TPZVec<Fad<REAL>> &xiInterior);

template void pztopology::TPZCube::MapToSide<Fad<REAL> >(int side, TPZVec<Fad<REAL> > &InternalPar, TPZVec<Fad<REAL> > &SidePar, TPZFMatrix<Fad<REAL> > &JacToSide);

template void pztopology::TPZCube::BlendFactorForSide<Fad<REAL>>(const int &, const TPZVec<Fad<REAL>> &, Fad<REAL> &,
                                                                   TPZVec<Fad<REAL>> &);
template void pztopology::TPZCube::TShape<Fad<REAL>>(const TPZVec<Fad<REAL>> &loc,TPZFMatrix<Fad<REAL>> &phi,TPZFMatrix<Fad<REAL>> &dphi);

template void pztopology::TPZCube::ComputeHDivDirections<Fad<REAL>>(TPZFMatrix<Fad<REAL>> &gradx, TPZFMatrix<Fad<REAL>> &directions);

template void pztopology::TPZCube::ComputeHCurlDirections<Fad<REAL>>(TPZFMatrix<Fad<REAL>> &gradx, TPZFMatrix<Fad<REAL>> &directions, const TPZVec<int> &transformationIds);
