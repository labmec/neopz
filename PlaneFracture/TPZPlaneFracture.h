/**
 * @file
 * @brief Contains the TPZPlaneFracture class which defines a plane of the fracture.
 */
#ifndef TPZPLANEFRACTUREH
#define TPZPLANEFRACTUREH

/*
 *  TPZPlaneFracture.h
 *  Crack
 *
 *  Created by Cesar Lucci on 09/08/10.
 *  Copyright 2010 LabMeC. All rights reserved.
 *
 */

//using namespace std;

#include <set>
#include "pzfmatrix.h"
#include "pzgmesh.h"
#include "pzvec.h"
#include "pzgeoelside.h"

#include "TPZPlaneFractureMaterials.h"



/// Real as tolerance 
const double __smallNum = 1.E-10;

/** @brief maximum element edge length */
const double __maxLength = 4.;

/** @brief plane fracture mesh height(Z) multiplier to set width(X), i.e.: (width = __lengthFactor x height)  */
const double __lengthFactor = 1.5;

/** @brief RefPatterns will be modulated to reduce the amount of options in the library */
/** @note Quantity of stretches for coarse edge intersection modulation */
const int __minTrimQTD = 4; //will be used for refpatterns

/** @brief Multiplier of __minTrimQTD for fine edge intersection modulation */
const int __TrimQTDmultiplier = 5; //will be used for searching direction between dots from poligonal chain



/** 
 * @brief Plane of the fracture. 
 * @author Cesar Lucci
 * @since 09/08/2010
 */

struct TPZFracture2DEl
{
public:
    
    TPZFracture2DEl()
    {
        fEdge.clear();
        fElem2D = NULL;
    }
    
    TPZFracture2DEl(TPZGeoEl * gel)
    {
        #ifdef DEBUG
        if(gel->HasSubElement())
        {
            std::cout << "The special kind of element TPZFracture2DEl cant be already refined!\n";
            std::cout << "See " << __PRETTY_FUNCTION__ << std::endl;
            DebugStop();
        }
        #endif
        
        fEdge.clear();
        fElem2D = gel;
        if(gel->Type() == ETriangle)
        {
            fEdge.insert(3);
            fEdge.insert(4);
            fEdge.insert(5);
        }
        else if(gel->Type() == EQuadrilateral)
        {
            fEdge.insert(4);
            fEdge.insert(5);
            fEdge.insert(6);
            fEdge.insert(7);
        }
        #ifdef DEBUG
        else
        {
            std::cout << "The special kind of element TPZFracture2DEl can be just 2D!\n";
            std::cout << "See " << __PRETTY_FUNCTION__ << std::endl;
            DebugStop();
        }
        #endif
    }
    
    ~TPZFracture2DEl()
    {
        fEdge.clear();
        fElem2D = NULL;
    }
    
    void RemoveThisEdge(int edg)
    {
        std::set<int>::iterator edgIt = fEdge.find(edg);
        if(edgIt != fEdge.end())
        {
            fEdge.erase(edgIt);
        }
    }
    
    int Id()
    {
        int elId = fElem2D->Id();
        return elId;
    }
    
    bool IsOver()
    {
        int nEdges = fEdge.size();
        if(nEdges == 0)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    
    std::set<int> fEdge;
    TPZGeoEl * fElem2D;
};

class TPZPlaneFracture
{
	public:
	
	/**
	 * @brief Constructor
	 * @param lw [in] : (positive) well length (positive)
	 * @param bulletDepthIni [in] : bullets perforation initial (positive) depth
	 * @param bulletDepthFin [in] : bullets perforation final (positive) depth
	 * @param pos_stress [in] : stress profile described by stretches
     *              Obs.: Stress profile in each stretch is linear
	 */
    TPZPlaneFracture(double lw, double bulletDepthIni, double bulletDepthFin, TPZVec< std::map<double,double> > & pos_stress);
    
	~TPZPlaneFracture();
    
    void GeneratePlaneMesh(std::set<double> & espacamento, double lengthFactor = __lengthFactor);
    void GenerateFullMesh(std::set<double> & espacamento, double lengthFactor = __lengthFactor);
    
    void GenerateNodesAtPlaneY(std::set<double> & espacamento, double lengthFactor, 
                               TPZVec< TPZVec<REAL> > & NodeCoord, int & nrows, int & ncols,
                               double Y);
	
	/**
	 * @brief Returns an GeoMesh based on original planeMesh, contemplating the poligonalChains geometry by refined elements
	 * @param poligonalChain [in] vector of boundary points coordinates
	 *
	 * @note Each vector position store x, y and z coordinates IN SEQUENCE of an poligonalChain geometry.
	 *
	 * Example:
	 *
	 *		x coordinate of first point of crack boundary: poligonalChain[0]\n
	 *		y coordinate of first point of crack boundary: poligonalChain[1]\n
	 *		z coordinate of first point of crack boundary: poligonalChain[2]\n
	 *		//
	 *		x coordinate of second point of crack boundary: poligonalChain[3]\n
	 *		y coordinate of second point of crack boundary: poligonalChain[4]\n
	 *		z coordinate of second point of crack boundary: poligonalChain[5]
	 */
	TPZGeoMesh * GetFractureMesh(const TPZVec<REAL> &poligonalChain);
		
//-----------------------------------------------------------------------------------------------------------------------------------------------------
	
	private:
	
	/*
	 * @brief Computes the edges of elements of fractMesh that are intercepted by the crack tip defined by poligonalChain points (defined by a vector coordinates)
	 * @param poligonalChain [in] vector of boundary points coordinates
	 *
	 * @note Each vector position store x, y and z coordinates IN SEQUENCE of an poligonalChain geometry.
	 *
	 * Example:
	 *
	 *		x coordinate of first point of crack boundary: poligonalChain[0]\n
	 *		y coordinate of first point of crack boundary: poligonalChain[1]\n
	 *		z coordinate of first point of crack boundary: poligonalChain[2]\n
	 *		//
	 *		x coordinate of second point of crack boundary: poligonalChain[3]\n
	 *		y coordinate of second point of crack boundary: poligonalChain[4]\n
	 *		z coordinate of second point of crack boundary: poligonalChain[5]
	 *
	 * @param fractMesh [in] geomesh of considered elements
	 * @param elId_TrimCoords [out] map that contains 1D element Id and a set of it trim 1D coordinates
	 * @param elIdSequence [out] the same of elId_TrimCoords, but keeps the trim 1D coordinates in generation sequence order
	 */
	void DetectEdgesCrossed(const TPZVec<REAL> &poligonalChain, TPZGeoMesh * fractMesh,
							std::map< int, std::set<double> > &elId_TrimCoords, std::list< std::pair<int,double> > &elIdSequence);
	
	/**
	 * @brief Returns (gel->Neighbour) relative to the side intersepted by the x+alpha.dx line
	 * @param gel [in] gel crossed by the line
	 * @param x [input and output data]  
	 *				x [as input] start point of line \n
	 *				x [as output] end point of line in gel and gel->Neighbour interface
	 * @param dx [in] direction of line from point x (input)
	 * @param alphaMin [in] if an start point (x) is in already in one edge of gel, it might be included or not in the intersections \n
	 *				        so, using alphaMin=0, in this case the first intersection (the x itself) is included...
	 *							   using alphaMin=1.E-10 (for example), in this case the first intersection (the x itself) is NOT included.
	 * @param elId_TrimCoords [out] map that contains the trim coordinates of 1D element, indexed by its Id (once 1D element was inserted in gel->Mesh)\n
	 *							    obs.: elId_TrimCoords was idealized to work in accumulative conception, i.e.,
	 *								    each time this method is called, this map grows!
	 * @param elIdSequence [out] the same of elId_TrimCoords, but keeps the trim 1D coordinates in generation sequence order
	 * @param pushback [in] set if dots on element edges will be inserted at the end, or not (i.e.: at beggining), of fCrackBoundary list
	 */
	TPZGeoEl * CrossToNextNeighbour(TPZGeoEl * gel, TPZVec<REAL> &x, TPZVec<REAL> dx, double alphaMin, std::map< int,
									std::set<double> > &elId_TrimCoords, std::list< std::pair<int,double> > &elIdSequence, bool pushback);
	
	/**
	 * @brief For a given element and internal point and an direction, computes the intersection
	 * coordinates with respect to it edges, and the respective intersected edge.
	 *
	 * @param gel [in] 2D geometric element whose edge will be intersected by (x+alphaX.dx) line
	 * @param x [in] element internal point coordinates
	 * @param dx [in] direction from point p
	 * @param edge [out] side Id of element edge that will be intersected by (x+alphaX.dx) line
	 * @param ExactIntersect [out] exact intersection coordinates with respect to edges parameter
	 * @param ModulatedIntersect [out] exact intersection coordinates, dragged to the nearest module (defined by fTrimQTD atribute)
	 * @param alphaMin [in] if an start point (x) is in already in one edge of gel, it might be included or not in the intersections\n
	 *					    so, using alphaMin=0, in this case the first intersection (the x itself) is included...\n
	 *						using alphaMin=1.E-10 (for example), in this case the first intersection (the x itself) is NOT included.
	 */
	bool EdgeIntersection(TPZGeoEl * gel, TPZVec<REAL> &x, TPZVec<REAL> &dx, TPZVec<int> &edge,
						  TPZVec< TPZVec<REAL> > &ExactIntersect, TPZVec< TPZVec<REAL> > &ModulatedIntersect, double alphaMin);
	
	/*
	 * @brief Returns an pointer to element of gMesh that contains the given point p (i.e.: point "p" belongs to element domain)
	 * @param p [in] point whose elements is going to be localized
	 * @param fractMesh [in] geomesh of elements candidates
	 * @param poligonalChain [in] vector of boundary points coordinates
	 *
	 * @note Each vector position store x, y and z coordinates IN SEQUENCE of an poligonalChain geometry.
	 *
	 * Example:
	 *
	 *		x coordinate of first point of crack boundary: poligonalChain[0]\n
	 *		y coordinate of first point of crack boundary: poligonalChain[1]\n
	 *		z coordinate of first point of crack boundary: poligonalChain[2]\n
	 *		//
	 *		x coordinate of second point of crack boundary: poligonalChain[3]\n
	 *		y coordinate of second point of crack boundary: poligonalChain[4]\n
	 *		z coordinate of second point of crack boundary: poligonalChain[5]
	 */
	static TPZGeoEl * PointElement(int p, TPZGeoMesh * fractMesh, const TPZVec<REAL> &poligonalChain);
	
	// alphaNode eh uma das solucoes do sistema: {x + alphaX.dx == node + alphaNode.dnode}, ou seja,
	// a norma que multiplica o vetor dnode e cruza a reta (x+alphaX.dx)
	/**
	 * @brief Given two vectorial lines \f$ x + alphaX.dx\f$ and \f$ node + alphaNode.dnode\f$, \n
	 * this method returns the alphaNode (norm that multiplies the unit vector dnode to intersect the line (x + alphax.dx) )
	 * @param x given point
	 * @param dx direction of line from point x
	 * @param node connect
	 * @param dnode
	 * @param norm [in] norm of edge that belongs to (node + alphaNode.dnode) line
	 * @param modulate [in] set if alphaNode will be modulated by stretches
	 * @param smooth [in] if alphaNode will be modulated, set if the stretches will be (norm/fTrimQTD) or smallest stretches \n
	 *							defined by (norm/(fTrimQTD*__TrimQTDmultiplier)).
	 * @note Obs.: alphaNode modulation is useful to reduce the possibilities of non regular refpatterns.
 	 * @note OBS.: dx and dnode MUST BE UNIT VECTORS!!!
	 */
	double ComputeAlphaNode(TPZVec<REAL> &x, TPZVec<REAL> &dx, TPZVec<REAL> &node, TPZVec<REAL> &dnode, double norm, bool modulate, bool smooth);
	
	// alphaX eh uma das solucoes do sistema: {x + alphaX.dx == node + alphaNode.dnode}, ou seja,
	// a norma que multiplica o vetor dx e cruza a reta (node+alphaNode.dnode)
	/**
	 * @brief Given two vectorial lines \f$ x + alphaX.dx\f$ and \f$ node + alphaNode.dnode\f$, \n
	 * this method returns the alphaX (norm that multiplies the unit vector dx to intersect the line (none + alphaNode.dnode) )
	 * @note dx and dnode MUST BE UNIT VECTORS!!!
	 */
	double ComputeAlphaX(TPZVec<REAL> &x, TPZVec<REAL> &dx, TPZVec<REAL> &node, TPZVec<REAL> &dnode);
	
	/**
	 * @brief Given 2 nodes (n0 and n1) and one point (x) in \f$ n0->n1 \f$ line, returns the point x in the line parametric space \f$ [-1,+1]\f$
	 */
	static double LinearComputeXInverse(TPZVec<REAL> x, TPZVec<REAL> n0, TPZVec<REAL> n1);
	
	/**
	 * @brief This method return a reffpattern of an unidimentional element that matches with the trim coordinates.
	 * @param TrimCoord Set of 1D element trimmed coordinates ( \f$[ -1 , +1 ]\f$ domain )
	 */
	static TPZAutoPointer<TPZRefPattern> Generate1DRefPatt(std::set<double> &TrimCoord);
	
	/**
	 * @brief Updates poligonal chain.
	 * @note The original Poligonal Chain (input data on GetFractureMesh method) is dots coordinates in the 2D mesh. \n
	 * This points normally are inside elements domain.\n
	 * The edges intersections of the original Poligonal Chain originate a new Poligonal Chain named poligonalChainUpdated 
	 */
	static void UpdatePoligonalChain(TPZGeoMesh * gmesh, std::list< std::pair<int,double> > &elIdSequence,
							  TPZVec<REAL> &poligonalChainUpdated);
	
	/**
	 * @param gmesh2D geometric mesh bi-dimensional
	 * @param gmesh3D geometric mesh three-dimensional
	 * @param elIdSequence - output data: list that contains 1D element Id and it trim 1D coordinates in generation sequence order
	 */
	void GenerateCrackBoundary(TPZGeoMesh * gmesh2D,
                               TPZGeoMesh * gmesh3D,
                               std::list< std::pair<int,double> > &elIdSequence,
                               TPZVec<int> &crackBoundaryElementsIds);
    
    void ChangeMaterialsOfFractureInterior(TPZGeoMesh * fullMesh, TPZVec<int> &crackBoundaryElementsIds);

	
//--------------------------------------------------------------------------------------------------------------------------------------------------
	
	/** @brief Original plane mesh (keeped intact for any poligonalChain configuration) */
	TPZGeoMesh * fPlaneMesh;
    
    /** @brief Original 3D mesh (keeped intact for any poligonalChain configuration) */
	TPZGeoMesh * fFullMesh;
    
    /** @brief Map that holds stress profile (position,<stressUp,stressDown>) */
    TPZVec< std::map<double,double> > fpos_stress;
    
	/** @brief It limits the amount of possible points in the edge of the elements */
	int fTrimQTD;
    
//--------------------------------------------------------------------------------------------------------------------------------------------------

public:
    /**
     * Metodo solicitado pelo Philippe para acidificacao
     */
    void GetSidesCrossedByPoligonalChain(const TPZVec<REAL> &poligonalChain, std::list< std::pair<TPZGeoElSide,double> > &sidesCrossed);
};

#endif