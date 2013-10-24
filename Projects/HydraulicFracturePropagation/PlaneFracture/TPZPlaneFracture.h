/**
 * @file
 * @brief Contains the TPZPlaneFracture class which defines a plane of the fracture.
 */
#ifndef TPZPLANEFRACTUREH
#define TPZPLANEFRACTUREH

/*
 *  TPZPlaneFracture.h
 *
 *  Created by Cesar Lucci on 09/08/10.
 *  Copyright 2010 LabMeC. All rights reserved.
 */

//using namespace std;

#include <set>
#include "pzgeoel.h"


#include "TPZPlaneFractureMaterials.h"
#include "tpzcompmeshreferred.h"


/// Real as tolerance 
const REAL __smallNum = 1.E-5;

/** @brief RefPatterns will be modulated to reduce the amount of options in the library */
/** @note Quantity of stretches for coarse edge intersection modulation */
const int __EdgeStretchesQTD = 5; //will be used for refpatterns

/** @brief Multiplier of __EdgeStretchesQTD for fine edge intersection modulation */
const int __TrimQTDmultiplier = 3; //will be used for searching direction between dots from poligonal chain



/** 
 * @brief Plane of the fracture. 
 * @author Cesar Lucci (Caju)
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
	 * @param bulletDepthTVDIni [in] : bullets perforation initial (TVD) depth
	 * @param bulletDepthTVDFin [in] : bullets perforation final (TVD) depth
	 * @param posTVD_stress [in] : stress profile described by stretches (TVD)
     *              Obs.: Stress profile in each stretch is linear, defined by (initialTVD,initialStress)~(finalTVD,finalStress)
     * @param xLength [in] : Reservoir length in x direction (crack propagation direction)
     * @param yLength [in] : Reservoir length in y direction (tickness of reservoir that couple fracture plane)
     * @param Lmax    [in] : Maximum element edge length
     *
     * TVD: True Vertical Depth (positive positions)
	 */
    TPZPlaneFracture(REAL lw, REAL bulletDepthTVDIni, REAL bulletDepthTVDFin,
                     TPZVec< std::map<REAL,REAL> > & posTVD_stress, REAL xLength, REAL yLength, REAL Lmax);
    
	~TPZPlaneFracture();
    
    /**
     * @brief Method that will run a FEM simmulation of a classical vertical plane fracture
     * @param poligonalChain [in] : Poligonal chain that represents the crack boundary
     * @param pressureInsideCrack [in] : uniform pressure inside crack
     * @param sigmaTraction [in] : uniform traction on the farfield surface (farfield surface have normal {0,1,0})
     * @param vtkFile [in] : VTK file name for post processing
     * @param printVTKfile [in] : flag that will determine if vtkFile will be generated (true) or not (false)
     */
    void RunThisFractureGeometry(const TPZVec<std::pair<REAL,REAL> > &poligonalChain,
                                 REAL pressureInsideCrack,
                                 std::string vtkFile,
                                 bool printVTKfile = false);
		
    //---------------------------------------------------------------------

protected:
    
    TPZCompMesh * GetFractureCompMesh(const TPZVec<std::pair<REAL,REAL> > &poligonalChain,
                                              int porder, REAL pressureInsideCrack);
    
    TPZCompMesh * GetFractureCompMeshNLinear(const TPZVec<std::pair<REAL,REAL> > &poligonalChain,
                                             int porder, REAL pressureInsideCrack);
    
    TPZCompMeshReferred * GetFractureCompMeshReferred(const TPZVec<std::pair<REAL,REAL> > &poligonalChain,
                                                      int porder, REAL pressureInsideCrack,
                                                      TPZCompMesh * cmeshRef);
    
    TPZCompMeshReferred * GetFractureCompMeshReferredNLinear(const TPZVec<std::pair<REAL,REAL> > &poligonalChain,
                                                             int porder, REAL pressureInsideCrack,
                                                             TPZCompMesh * cmeshRef);
    
    /**
	 * @brief Returns an GeoMesh based on original planeMesh, contemplating the poligonalChains geometry by refined elements
	 * @param poligonalChain [in] : vector of boundary points coordinates
	 *
	 * @note Each vector position store x and z coordinates IN SEQUENCE of an poligonalChain geometry.
	 *
	 * Example:
	 *
	 *		x,z coordinates of first point of crack boundary: poligonalChain[0]\n
	 *		x,z coordinates of second point of crack boundary: poligonalChain[1]\n
     *      x,z coordinates of third point of crack boundary: poligonalChain[2]\n
     *      etc...
	 */
	TPZGeoMesh * GetFractureGeoMesh(const TPZVec<std::pair<REAL,REAL> > &poligonalChain);

    /** @brief Generation of the persistent full mesh (2D and 3D) that contains the fracture and its porous media
     *  @note This method set the fPreservedMesh atribute that will not be changed for every fracture time step
     *  @param espacamento [in] : espacamento vertical que define interfaces entre camadas horizontais
     *  @param xLength [in] : Reservoir length in x direction (crack propagation direction)
     *  @param yLength [in] : Reservoir length in y direction (tickness that couple fracture plane)
     */
    void GeneratePreservedMesh(std::list<REAL> & espacamento, REAL xLength, REAL yLength);

    /** @brief Method used for the mesh generator methods GeneratePlaneMesh and GeneratePreservedMesh
     *  @note For a given xz plane (defined by Y coordinate), generate the node grid coordinates
     *  @param espacamento [in] : espacamento vertical que define interfaces entre camadas horizontais
     *  @param xLength [in] : Reservoir length in x direction (crack propagation direction)
     */
    void GenerateNodesAtPlaneY(std::list<REAL> & espacamento, REAL xLength,
                               TPZVec< TPZVec<REAL> > & NodeCoord, long & nrows, long & ncols,
                               REAL Y);

	/*
	 * @brief Computes the edges of elements of fractMesh that are intercepted by the crack tip defined by poligonalChain points (defined by a vector coordinates)
	 * @param poligonalChain [in] : vector of boundary points coordinates
	 *
	 * @note Each vector position store x, y and z coordinates IN SEQUENCE of an poligonalChain geometry.
	 *
	 * Example:
	 *
	 *		x coordinate of first point of crack boundary: poligonalChain[0].first\n
	 *		z coordinate of first point of crack boundary: poligonalChain[0].second\n
	 *		//
	 *		x coordinate of second point of crack boundary: poligonalChain[1].first\n
	 *		z coordinate of second point of crack boundary: poligonalChain[1].second
	 *
	 * @param fractMesh [in] : geomesh of considered elements
	 * @param elId_TrimCoords [out] : map that contains 1D element Id and a set of it trim 1D coordinates
	 * @param elIdSequence [out] : the same of elId_TrimCoords, but keeps the trim 1D coordinates in generation sequence order
	 */
	void DetectEdgesCrossed(const TPZVec<std::pair<REAL,REAL> > &poligonalChain, TPZGeoMesh * fractMesh,
							std::map< long, std::set<REAL> > &elId_TrimCoords, std::list< std::pair<long,REAL> > &elIdSequence);

	/**
	 * @brief Returns the next geoel (from given geoel) relative to the given direction by the x+alpha.dx
	 * @param gel [in] : initial geoel
	 * @param x [input and output data]  
	 *				x [as input] : start point of line \n
	 *				x [as output] : end point of line in gel interface
	 * @param dx [in] : direction of line from point x (input)
	 * @param alphaMin [in] : if an start point (x) is in already in one edge of gel, it might be included or not in the intersections \n
	 *				        so, using alphaMin=0, in this case the first intersection (the x itself) is included...
	 *							   using alphaMin=1.E-10 (for example), in this case the first intersection (the x itself) is NOT included.
	 * @param elId_TrimCoords [out] : map that contains the trim coordinates of 1D element, indexed by its Id (once 1D element was inserted in gel->Mesh)\n
	 *							    obs.: elId_TrimCoords was idealized to work in accumulative conception, i.e.,
	 *								    each time this method is called, this map grows!
	 * @param elIdSequence [out] : the same of elId_TrimCoords, but keeps the trim 1D coordinates in generation sequence order
	 * @param pushback [in] : set if dots on element edges will be inserted at the end, or not (i.e.: at beggining), of fCrackBoundary list
	 */
	static TPZGeoEl * CrossToNextNeighbour(TPZGeoEl * gel, TPZVec<REAL> &x, TPZVec<REAL> dx, REAL alphaMin, std::map< long,
									std::set<REAL> > &elId_TrimCoords, std::list< std::pair<long,REAL> > &elIdSequence, bool pushback,
                                    int planeAxe0, int planeAxe1, int planeNormal,
                                    bool closingFracture);

	static bool EdgeIntersection(TPZGeoEl * gel, TPZVec<REAL> &x, TPZVec<REAL> &dx, TPZVec<int> &edge, 
                                 TPZVec< TPZVec<REAL> > &ExactIntersect, REAL alphaMin,
                                 int planeAxe0, int planeAxe1, int planeNormal);

    /**
	 * @brief For a given element and internal point and an direction, computes the intersection
	 * coordinates with respect to it edges, and the respective intersected edge.
	 *
	 * @param gel [in] : 2D geometric element whose edge will be intersected by (x+alphaX.dx) line
	 * @param x [in] : element internal point coordinates
	 * @param dx [in] : direction from point p
	 * @param edge [out] : side Id of element edge that will be intersected by (x+alphaX.dx) line
     *          @brief :    this vector normally assumes size=1, but when alphaMin=0, assumes size=2
     *                      (the edge where the given coordinate lies and the next edge where the given direction points)
	 * @param ExactIntersect [out] : exact intersection coordinates with respect to edges parameter
	 * @param ModulatedIntersect [out] : exact intersection coordinates, dragged to the nearest module (defined by __EdgeStretchesQTD constant)
	 * @param alphaMin [in] : if an start point (x) is in already in one edge of gel, it might be included or not in the intersections\n
	 *					    so, using alphaMin=0, in this case the first intersection (the x itself) is included...\n
	 *						using alphaMin=1.E-10 (for example), in this case the first intersection (the x itself) is NOT included.
	 */
	static bool EdgeIntersection(TPZGeoEl * gel, TPZVec<REAL> &x, TPZVec<REAL> &dx, TPZVec<int> &edge,
                                 TPZVec< TPZVec<REAL> > &ExactIntersect, TPZVec< TPZVec<REAL> > &ModulatedIntersect, REAL alphaMin,
                                 int planeAxe0, int planeAxe1, int planeNormal);

	// alphaNode eh uma das solucoes do sistema: {x + alphaX.dx == node + alphaNode.dnode}, ou seja,
	// a norma que multiplica o vetor dnode e cruza a reta (x+alphaX.dx)
	/**
	 * @brief Given two vectorial lines \f$ x + alphaX.dx\f$ and \f$ node + alphaNode.dnode\f$, \n
	 * this method returns the alphaNode (norm that multiplies the unit vector dnode to intersect the line (x + alphax.dx) )
	 * @param x given point
	 * @param dx direction of line from point x
	 * @param node connect
	 * @param dnode
	 * @param norm [in] : norm of edge that belongs to (node + alphaNode.dnode) line
	 * @param modulate [in] : set if alphaNode will be modulated by stretches
	 * @param smooth [in] : if alphaNode will be modulated, set if the stretches will be (norm/__EdgeStretchesQTD) or smallest stretches \n
	 *							defined by (norm/(fEdgeStretchesQTD*__TrimQTDmultiplier)).
	 * @note Obs.: alphaNode modulation is useful to reduce the possibilities of non regular refpatterns.
 	 * @note OBS.: dx and dnode MUST BE UNIT VECTORS!!!
	 */
	static REAL ComputeAlphaNode(TPZVec<REAL> &x, TPZVec<REAL> &dx, TPZVec<REAL> &node, TPZVec<REAL> &dnode, REAL norm, bool modulate, bool smooth,
                                   int planeAxe0, int planeAxe1, int planeNormal);
	
	// alphaX eh uma das solucoes do sistema: {x + alphaX.dx == node + alphaNode.dnode}, ou seja,
	// a norma que multiplica o vetor dx e cruza a reta (node+alphaNode.dnode)
	/**
	 * @brief Given two vectorial lines \f$ x + alphaX.dx\f$ and \f$ node + alphaNode.dnode\f$, \n
	 * this method returns the alphaX (norm that multiplies the unit vector dx to intersect the line (none + alphaNode.dnode) )
	 * @note dx and dnode MUST BE UNIT VECTORS!!!
	 */
	static REAL ComputeAlphaX(TPZVec<REAL> &x, TPZVec<REAL> &dx, TPZVec<REAL> &node, TPZVec<REAL> &dnode,
                                int planeAxe0, int planeAxe1, int planeNormal);
	
	/**
	 * @brief Given 2 nodes (n0 and n1) and one point (x) in \f$ n0->n1 \f$ line, returns the point x in the line parametric space \f$ [-1,+1]\f$
	 */
	static REAL LinearComputeXInverse(TPZVec<REAL> x, TPZVec<REAL> n0, TPZVec<REAL> n1);
	
	/**
	 * @brief This method return a reffpattern of an unidimentional element that matches with the trim coordinates.
	 * @param TrimCoord Set of 1D element trimmed coordinates ( \f$[ -1 , +1 ]\f$ domain )
	 */
	static TPZAutoPointer<TPZRefPattern> Generate1DRefPatt(std::set<REAL> &TrimCoord);
	
	/**
	 * @brief Updates poligonal chain.
	 * @note The original Poligonal Chain (input data on GetFractureMesh method) is dots coordinates in the 2D mesh. \n
	 * This points normally are inside elements domain.\n
	 * The edges intersections of the original Poligonal Chain originate a new Poligonal Chain named poligonalChainUpdated 
	 */
	static void UpdatePoligonalChain(TPZGeoMesh * gmesh, std::list< std::pair<long,REAL> > &elIdSequence,
							  TPZVec<std::pair<REAL,REAL> > &poligonalChainUpdated);
	
	/**
	 * @param gmesh2D geometric mesh bi-dimensional
	 * @param gmesh3D geometric mesh three-dimensional
	 * @param elIdSequence - output data: list that contains 1D element Id and it trim 1D coordinates in generation sequence order
	 */
	void GenerateCrackBoundary(TPZGeoMesh * gmesh3D,
                               std::list< std::pair<long,REAL> > &elIdSequence);
    
    /**
     * @brief Fill fcrackQpointsElementsIds atribute with the elements (and its sides) that toutch cracktip
     */
    void SeparateElementsInMaterialSets(TPZGeoMesh * refinedMesh);

    /**
     * @brief Once the fcrackQpointsElementsIds atribute is filled (by HuntElementsSurroundingCrackTip method),\n
     *          the respective elements must be turned into quarterpoints,\n
     *          ruled by involved sides (that could be more than one by element)
     */
    void TurnIntoQuarterPoint(TPZGeoMesh * refinedMesh);
    
    /**
     * @brief Refinement proceedings to increase solution quality
     */
    void RefinementProceedings(TPZGeoMesh * refinedMesh);
    
    /**
     * @brief Returns if a given element is from boundary of domain
     */
    static bool IsBoundaryMaterial(TPZGeoEl * gel);

    //** just for visualize given dots in vtk */
    static void InsertDots4VTK(TPZGeoMesh * gmesh, const TPZVec<REAL> &fractureDots);
	
    //-------------------------------------------------------------------------------
	
    /** @brief Original 3D mesh (keeped intact for any poligonalChain configuration) */
	TPZGeoMesh * fPreservedMesh;
    
    /** @brief Map that holds stress profile (TVD position,<stressUp,stressDown>) */
    TPZVec< std::map<REAL,REAL> > fposTVD_stress;
    
    /** @brief 1D elements Ids that compose crack boundary */
    TPZVec<long> fcrackBoundaryElementsIds;
    
    /** @brief initial element index in point search (just for optimization) */
    long fInitialElIndex;
    
    /** @brief maximum element edge length */
    REAL fLmax;
};

#endif