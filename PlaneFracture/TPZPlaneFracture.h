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
#include "pzgeoel.h"


#include "TPZPlaneFractureMaterials.h"



/// Real as tolerance 
const double __smallNum = 1.E-5;

/** @brief maximum element edge length */
const double __maxLength = 4.;

/** @brief plane fracture mesh height(Z) multiplier to set width(X), i.e.: (width = __lengthFactor x height)  */
const double __lengthFactor = 1.;

/** @brief RefPatterns will be modulated to reduce the amount of options in the library */
/** @note Quantity of stretches for coarse edge intersection modulation */
const int __EdgeStretchesQTD = 4; //will be used for refpatterns

/** @brief Multiplier of __EdgeStretchesQTD for fine edge intersection modulation */
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
    
    TPZPlaneFracture();
	
	/**
	 * @brief Constructor
	 * @param lw [in] : (positive) well length (positive)
	 * @param bulletDepthIni [in] : bullets perforation initial (TVD) depth
	 * @param bulletDepthFin [in] : bullets perforation final (TVD) depth
	 * @param pos_stress [in] : stress profile described by stretches (TVD)
     *              Obs.: Stress profile in each stretch is linear
     *
     * TVD: Total vertical depth (positive positions)
	 */
    TPZPlaneFracture(double lw, double bulletDepthIni, double bulletDepthFin, TPZVec< std::map<double,double> > & pos_stress);
    
	~TPZPlaneFracture();
    
    void RunThisFractureGeometry(const TPZVec<REAL> &poligonalChain, std::string vtkFile);
		
//-----------------------------------------------------------------------------------------------------------------------------------------------------
	
	private:
    
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
	TPZGeoMesh * GetFractureGeoMesh(const TPZVec<REAL> &poligonalChain);
    
    TPZCompMesh * GetFractureCompMesh(const TPZVec<REAL> &poligonalChain, int porder);
    
    /** @brief Generation of the persistent 2D mesh that contains the fracture
     *  @note This method set the fPlaneMesh atribute that will not be changed for every fracture time step
     */
    void GeneratePlaneMesh(std::list<double> & espacamento, double lengthFactor = __lengthFactor);
    
    /** @brief Generation of the persistent full mesh (2D and 3D) that contains the fracture and its porous media
     *  @note This method set the fFullMesh atribute that will not be changed for every fracture time step
     */
    void GenerateFullMesh(std::list<double> & espacamento, double lengthFactor = __lengthFactor);
    
    /** @brief Method used for the mesh generator methods GeneratePlaneMesh and GenerateFullMesh
     *  @note For a given xz plane (defined by Y coordinate), generate the node grid coordinates
     */
    void GenerateNodesAtPlaneY(std::list<double> & espacamento, double lengthFactor,
                               TPZVec< TPZVec<REAL> > & NodeCoord, int & nrows, int & ncols,
                               double Y);
	
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
	 *				x [as output] end point of line in gel interface
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
	 * @param ModulatedIntersect [out] exact intersection coordinates, dragged to the nearest module (defined by __EdgeStretchesQTD constant)
	 * @param alphaMin [in] if an start point (x) is in already in one edge of gel, it might be included or not in the intersections\n
	 *					    so, using alphaMin=0, in this case the first intersection (the x itself) is included...\n
	 *						using alphaMin=1.E-10 (for example), in this case the first intersection (the x itself) is NOT included.
	 */
	bool EdgeIntersection(TPZGeoEl * gel, TPZVec<REAL> &x, TPZVec<REAL> &dx, TPZVec<int> &edge,
						  TPZVec< TPZVec<REAL> > &ExactIntersect, TPZVec< TPZVec<REAL> > &ModulatedIntersect, double alphaMin);
	
	/*
	 * @brief Returns the Id of the element of fPlaneMesh (atribute) that contains the given coordinates (x).
	 */
	int PointElementOnPlaneMesh(const TPZVec<REAL> & x);
    
    /**
     * @brief Returns an pointer to element of given mesh (fullMesh) that contains the given coordinates (x).
     * @param x [in] coordinates whose elements is going to be localized.
     * @param fullMesh [in] geomesh of elements candidates.
     */
    TPZGeoEl * PointElementOnFullMesh(const TPZVec<REAL> & x, TPZGeoMesh * fullMesh);
    
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
	 * @param smooth [in] if alphaNode will be modulated, set if the stretches will be (norm/__EdgeStretchesQTD) or smallest stretches \n
	 *							defined by (norm/(fEdgeStretchesQTD*__TrimQTDmultiplier)).
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
                               std::list< std::pair<int,double> > &elIdSequence);
    
    /**
     * @brief Fill fcrackQpointsElementsIds atribute with the elements (and its sides) that toutch cracktip
     */
    void SeparateElementsInMaterialSets(TPZGeoMesh * fullMesh);

    /**
     * @brief Once the fcrackQpointsElementsIds atribute is filled (by HuntElementsSurroundingCrackTip method),\n
     *          the respective elements must be turned into quarterpoints,\n
     *          ruled by involved sides (that could be more than one by element)
     */
    void TurnIntoQuarterPoint(TPZGeoMesh * fullMesh);
    
    /**
     * @brief Returns if a given element touch cracktip and respective sides ids (in case of return true)
     */
    bool TouchCrackTip(TPZGeoEl * gel, std::set<int> &bySides);
    
    /**
     * @brief Returns if a given element is from boundary of domain
     */
    static bool IsBoundaryMaterial(TPZGeoEl * gel);

	
//--------------------------------------------------------------------------------------------------------------------------------------------------
	
	/** @brief Original plane mesh (keeped intact for any poligonalChain configuration) */
	TPZGeoMesh * fPlaneMesh;
    
    /** @brief Original 3D mesh (keeped intact for any poligonalChain configuration) */
	TPZGeoMesh * fFullMesh;
    
    /** @brief Map that holds stress profile (position,<stressUp,stressDown>) */
    TPZVec< std::map<double,double> > fpos_stress;
    
    /** @brief 1D elements Ids that compose crack boundary */
    TPZVec<int> fcrackBoundaryElementsIds;
    
    /** @brief Quarter points 3D elements Ids that surround crack boundary */
    //map< elementId , set< sides of this element that touch 1d cracktip > >
    std::map< int,std::set<int> > fcrackQpointsElementsIds;
    
    /** @brief smaller radius that defines the cylinder that involves the
     crack boundary inside quarter points elements (J integral) */
    double fMinimumRadius;
    
//--------------------------------------------------------------------------------------------------------------------------------------------------

    //** just for visualize given dots in vtk */
    static void InsertDots4VTK(TPZGeoMesh * gmesh, const TPZVec<REAL> &fractureDots);
    
public:
    /**
     * Metodo solicitado pelo Philippe para acidificacao
     */
    void GetSidesCrossedByPoligonalChain(const TPZVec<REAL> &poligonalChain, std::list< std::pair<TPZGeoElSide,double> > &sidesCrossed);
};

#endif