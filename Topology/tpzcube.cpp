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


#ifdef _AUTODIFF
#include "fad.h"
#endif

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.topology.pzcube"));
#endif
using namespace std;

namespace pztopology {

	/**
	 * @brief For each face (quadrilateral sides) was enumerated in sequence the sides contained in the closure of them.
	 * For example: First face (side 20) in its closure contains side 0, 1, 2, 3 (vertices), 8, 9, 10, 11, (edges) and it self 
	 */
	static int FaceConnectLocId[6][9] = { {0,1,2,3,8,9,10,11,20},{0,1,5,4,8,13,16,12,21},
		{1,2,6,5,9,14,17,13,22},{3,2,6,7,10,14,18,15,23},//{2,3,7,6,10,15,18,14,23}
		{0,3,7,4,11,15,19,12,24},{4,5,6,7,16,17,18,19,25} };
	
	/** @brief For each face was enumerated the pontoal sides (vertices) */
	int TPZCube::FaceNodes[6][4]  = { {0,1,2,3},{0,1,5,4},{1,2,6,5},{3,2,6,7},{0,3,7,4},{4,5,6,7} };
	
	/** @brief For each edge was enumerated the pontoal sides (vertices) */
	int TPZCube::SideNodes[12][2]  = { {0,1},{1,2},{2,3},{3,0},{0,4},{1,5},
		{2,6},{3,7},{4,5},{5,6},{6,7},{7,4} };
	
	/** @brief For each face was enumerated the vertice sides on its main diagonal */
	int TPZCube::ShapeFaceId[6][2] = { {0,2},{0,5},{1,6},{3,6},{0,7},{4,6} };


    int TPZCube::fPermutations[48][27] =
            {
                    {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26},/*000*/
                    {0,1,5,4,3,2,6,7,8,13,16,12,11,9,17,19,10,14,18,15,21,20,22,25,24,23,26},/*001*/
                    {0,3,2,1,4,7,6,5,11,10,9,8,12,15,14,13,19,18,17,16,20,24,23,22,21,25,26},/*002*/
                    {0,3,7,4,1,2,6,5,11,15,19,12,8,10,18,16,9,14,17,13,24,20,23,25,21,22,26},/*003*/
                    {0,4,5,1,3,7,6,2,12,16,13,8,11,19,17,9,15,18,14,10,21,24,25,22,20,23,26},/*004*/
                    {0,4,7,3,1,5,6,2,12,19,15,11,8,16,18,10,13,17,14,9,24,21,25,23,20,22,26},/*005*/
                    {1,0,3,2,5,4,7,6,8,11,10,9,13,12,15,14,16,19,18,17,20,21,24,23,22,25,26},/*006*/
                    {1,0,4,5,2,3,7,6,8,12,16,13,9,11,19,17,10,15,18,14,21,20,24,25,22,23,26},/*007*/
                    {1,2,3,0,5,6,7,4,9,10,11,8,13,14,15,12,17,18,19,16,20,22,23,24,21,25,26},/*008*/
                    {1,2,6,5,0,3,7,4,9,14,17,13,8,10,18,16,11,15,19,12,22,20,23,25,21,24,26},/*009*/
                    {1,5,4,0,2,6,7,3,13,16,12,8,9,17,19,11,14,18,15,10,21,22,25,24,20,23,26},/*010*/
                    {1,5,6,2,0,4,7,3,13,17,14,9,8,16,18,10,12,19,15,11,22,21,25,23,20,24,26},/*011*/
                    {2,1,0,3,6,5,4,7,9,8,11,10,14,13,12,15,17,16,19,18,20,22,21,24,23,25,26},/*012*/
                    {2,1,5,6,3,0,4,7,9,13,17,14,10,8,16,18,11,12,19,15,22,20,21,25,23,24,26},/*013*/
                    {2,3,0,1,6,7,4,5,10,11,8,9,14,15,12,13,18,19,16,17,20,23,24,21,22,25,26},/*014*/
                    {2,3,7,6,1,0,4,5,10,15,18,14,9,11,19,17,8,12,16,13,23,20,24,25,22,21,26},/*015*/
                    {2,6,5,1,3,7,4,0,14,17,13,9,10,18,16,8,15,19,12,11,22,23,25,21,20,24,26},/*016*/
                    {2,6,7,3,1,5,4,0,14,18,15,10,9,17,19,11,13,16,12,8,23,22,25,24,20,21,26},/*017*/
                    {3,0,1,2,7,4,5,6,11,8,9,10,15,12,13,14,19,16,17,18,20,24,21,22,23,25,26},/*018*/
                    {3,0,4,7,2,1,5,6,11,12,19,15,10,8,16,18,9,13,17,14,24,20,21,25,23,22,26},/*019*/
                    {3,2,1,0,7,6,5,4,10,9,8,11,15,14,13,12,18,17,16,19,20,23,22,21,24,25,26},/*020*/
                    {3,2,6,7,0,1,5,4,10,14,18,15,11,9,17,19,8,13,16,12,23,20,22,25,24,21,26},/*021*/
                    {3,7,4,0,2,6,5,1,15,19,12,11,10,18,16,8,14,17,13,9,24,23,25,21,20,22,26},/*022*/
                    {3,7,6,2,0,4,5,1,15,18,14,10,11,19,17,9,12,16,13,8,23,24,25,22,20,21,26},/*023*/
                    {4,0,1,5,7,3,2,6,12,8,13,16,19,11,9,17,15,10,14,18,21,24,20,22,25,23,26},/*024*/
                    {4,0,3,7,5,1,2,6,12,11,15,19,16,8,10,18,13,9,14,17,24,21,20,23,25,22,26},/*025*/
                    {4,5,1,0,7,6,2,3,16,13,8,12,19,17,9,11,18,14,10,15,21,25,22,20,24,23,26},/*026*/
                    {4,5,6,7,0,1,2,3,16,17,18,19,12,13,14,15,8,9,10,11,25,21,22,23,24,20,26},/*027*/
                    {4,7,3,0,5,6,2,1,19,15,11,12,16,18,10,8,17,14,9,13,24,25,23,20,21,22,26},/*028*/
                    {4,7,6,5,0,3,2,1,19,18,17,16,12,15,14,13,11,10,9,8,25,24,23,22,21,20,26},/*029*/
                    {5,1,0,4,6,2,3,7,13,8,12,16,17,9,11,19,14,10,15,18,21,22,20,24,25,23,26},/*030*/
                    {5,1,2,6,4,0,3,7,13,9,14,17,16,8,10,18,12,11,15,19,22,21,20,23,25,24,26},/*031*/
                    {5,4,0,1,6,7,3,2,16,12,8,13,17,19,11,9,18,15,10,14,21,25,24,20,22,23,26},/*032*/
                    {5,4,7,6,1,0,3,2,16,19,18,17,13,12,15,14,8,11,10,9,25,21,24,23,22,20,26},/*033*/
                    {5,6,2,1,4,7,3,0,17,14,9,13,16,18,10,8,19,15,11,12,22,25,23,20,21,24,26},/*034*/
                    {5,6,7,4,1,2,3,0,17,18,19,16,13,14,15,12,9,10,11,8,25,22,23,24,21,20,26},/*035*/
                    {6,2,1,5,7,3,0,4,14,9,13,17,18,10,8,16,15,11,12,19,22,23,20,21,25,24,26},/*036*/
                    {6,2,3,7,5,1,0,4,14,10,15,18,17,9,11,19,13,8,12,16,23,22,20,24,25,21,26},/*037*/
                    {6,5,1,2,7,4,0,3,17,13,9,14,18,16,8,10,19,12,11,15,22,25,21,20,23,24,26},/*038*/
                    {6,5,4,7,2,1,0,3,17,16,19,18,14,13,12,15,9,8,11,10,25,22,21,24,23,20,26},/*039*/
                    {6,7,3,2,5,4,0,1,18,15,10,14,17,19,11,9,16,12,8,13,23,25,24,20,22,21,26},/*040*/
                    {6,7,4,5,2,3,0,1,18,19,16,17,14,15,12,13,10,11,8,9,25,23,24,21,22,20,26},/*041*/
                    {7,3,0,4,6,2,1,5,15,11,12,19,18,10,8,16,14,9,13,17,24,23,20,21,25,22,26},/*042*/
                    {7,3,2,6,4,0,1,5,15,10,14,18,19,11,9,17,12,8,13,16,23,24,20,22,25,21,26},/*043*/
                    {7,4,0,3,6,5,1,2,19,12,11,15,18,16,8,10,17,13,9,14,24,25,21,20,23,22,26},/*044*/
                    {7,4,5,6,3,0,1,2,19,16,17,18,15,12,13,14,11,8,9,10,25,24,21,22,23,20,26},/*045*/
                    {7,6,2,3,4,5,1,0,18,14,10,15,19,17,9,11,16,13,8,12,23,25,22,20,24,21,26},/*046*/
                    {7,6,5,4,3,2,1,0,18,17,16,19,15,14,13,12,10,9,8,11,25,23,22,21,24,20,26} /*047*/
            };

	/** @brief Vector of the dimension for each side */
	static int sidedimension[27] = {0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,3};
	
	/** @brief Vector with the number of vertices contained in the closure of the side */
	static int nsidenodes[27] = {1,1,1,1,1,1,1,1,
		2,2,2,2,2,2,2,2,2,2,2,2,
		4,4,4,4,4,4,
		8};
	

	/** @brief For each side was stored the sides connected with it but of the higher dimension */ 
	static int highsides[27][7] = {
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
	static int nhighdimsides[27] = {7,7,7,7,7,7,7,7,3,3,3,3,3,3,3,3,3,3,3,3,1,1,1,1,1,1,0};

	/** @brief The transformations for each side over neighboard side with higher dimension */
	static REAL sidetosidetransforms[27][7][4][3] = {
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
	static REAL MidSideNode[27][3] = {
		/*00*/{-1.,-1.,-1.},/*01*/{1.,-1.,-1.},/*02*/{1.,1.,-1.},/*03*/{-1.,1.,-1.},
		/*04*/{-1.,-1., 1.},/*05*/{1.,-1., 1.},/*06*/{1.,1., 1.},/*07*/{-1.,1., 1.},
		/*08*/{ 0.,-1.,-1.},/*09*/{1., 0.,-1.},/*10*/{0.,1.,-1.},/*11*/{-1.,0.,-1.},
		/*12*/{-1.,-1., 0.},/*13*/{1.,-1., 0.},/*14*/{1.,1., 0.},/*15*/{-1.,1., 0.},
		/*16*/{ 0.,-1., 1.},/*17*/{1., 0., 1.},/*18*/{0.,1., 1.},/*19*/{-1.,0., 1.},
		/*20*/{ 0., 0.,-1.},/*21*/{0.,-1., 0.},/*22*/{1.,0., 0.},/*23*/{ 0.,1., 0.},
		/*24*/{-1., 0., 0.},/*25*/{0., 0., 1.},/*26*/{0.,0., 0.} };
    
    static REAL bCubo[81][3] = // direcao perpendicular ao lado
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
    
    static REAL t1Cubo[81][3] = // diretor da aresta (escolhido para formar uma base positivamente orientada)
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
    static REAL t2Cubo[81][3] = // diretor da aresta (escolhido para formar uma base positivamente orientada)
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
    
	static int vectorsideorderC [81] =
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

    static int bilinearounao [81] =   {
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
    
    static int direcaoksioueta [81] = {
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
            #ifdef LOG4CXX
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
	
	
	MElementType TPZCube::Type()
	{
		return ECube;
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
	int TPZCube::GetTransformId(TPZVec<int64_t> &id)
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
	int TPZCube::GetTransformId(int side, TPZVec<int64_t> &id)
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
				LOGPZ_ERROR(logger,"Please Implement me")
				return -1;
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
#ifdef _AUTODIFF

template bool pztopology::TPZCube::CheckProjectionForSingularity<Fad<REAL>>(const int &side, const TPZVec<Fad<REAL>> &xiInterior);

template void pztopology::TPZCube::MapToSide<Fad<REAL> >(int side, TPZVec<Fad<REAL> > &InternalPar, TPZVec<Fad<REAL> > &SidePar, TPZFMatrix<Fad<REAL> > &JacToSide);

template void pztopology::TPZCube::BlendFactorForSide<Fad<REAL>>(const int &, const TPZVec<Fad<REAL>> &, Fad<REAL> &,
                                                                   TPZVec<Fad<REAL>> &);
template void pztopology::TPZCube::TShape<Fad<REAL>>(const TPZVec<Fad<REAL>> &loc,TPZFMatrix<Fad<REAL>> &phi,TPZFMatrix<Fad<REAL>> &dphi);

template void pztopology::TPZCube::ComputeHDivDirections<Fad<REAL>>(TPZFMatrix<Fad<REAL>> &gradx, TPZFMatrix<Fad<REAL>> &directions);
#endif
