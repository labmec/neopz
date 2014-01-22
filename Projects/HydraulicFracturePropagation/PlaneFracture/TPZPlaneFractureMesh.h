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


#include "TPZPlaneFractureData.h"
#include "tpzcompmeshreferred.h"
#include "TPZPlaneFractCouplingMat.h"


/// Real as tolerance 
const REAL __smallNum = 1.E-5;

/** @brief RefPatterns will be modulated to reduce the amount of options in the library */
/** @note Quantity of stretches for coarse edge intersection modulation */
const int __EdgeStretchesQTD = 10; //will be used for refpatterns

/** @brief Multiplier of __EdgeStretchesQTD for fine edge intersection modulation */
const int __TrimQTDmultiplier = 2; //will be used for searching direction between dots from poligonal chain


struct TPZLayerProperties
{
public:
    TPZLayerProperties()
    {
        fYoung = 0.;
        fPoisson = 0.;
        fSigmaMax = 0.;
        fSigmaMin = 0.;
        fSigmaConf = 0.;
        fTVDini = 0.;
        fTVDfin = 0.;
        fKIc = 0.;
        fCl = 0.;
        fPe = 0.;
        fgradPref = 0.;
        fvsp = 0.;
    }
    
    TPZLayerProperties(REAL Young, REAL Poisson, REAL SigMax, REAL SigMin, REAL SigConf, REAL TVDi, REAL TVDf,
                       REAL KIc, REAL Cl, REAL Pe, REAL gradPref, REAL vsp)
    {
#ifdef DEBUG
        REAL tol = 1.E-9;
        if(SigMax > 0. + tol || SigMin > 0. + tol || SigConf > 0. + tol)
        {
            std::cout << "\n\nPre-stress must be null or compressive!\n";
            std::cout << "\n\nSee " << __PRETTY_FUNCTION__ << "\n\n";
            DebugStop();
        }
        if(Pe < 0. - tol)
        {
            std::cout << "\n\nStatic pressure (Pe) must be null or positive!\n";
            std::cout << "\n\nSee " << __PRETTY_FUNCTION__ << "\n\n";
            DebugStop();
        }
        if(gradPref < 0. + tol)
        {
            std::cout << "\n\nReference pressure (gradPref) must be positive!\n";
            std::cout << "\n\nSee " << __PRETTY_FUNCTION__ << "\n\n";
            DebugStop();
        }
#endif
        
        fYoung = Young;
        fPoisson = Poisson;
        fSigmaMax = SigMax;
        fSigmaMin = SigMin;
        fSigmaConf = SigConf;
        fTVDini = TVDi;
        fTVDfin = TVDf;
        fKIc = KIc;
        fCl = Cl;
        fPe = Pe;
        fgradPref = gradPref;
        fvsp = vsp;
    }
    TPZLayerProperties(const TPZLayerProperties & cp)
    {
        fYoung = cp.fYoung;
        fPoisson = cp.fPoisson;
        fSigmaMax = cp.fSigmaMax;
        fSigmaMin = cp.fSigmaMin;
        fSigmaConf = cp.fSigmaConf;
        fTVDini = cp.fTVDini;
        fTVDfin = cp.fTVDfin;
        fKIc = cp.fKIc;
        fCl = cp.fCl;
        fPe = cp.fPe;
        fgradPref = 0.;
        fvsp = cp.fvsp;
    }
    ~TPZLayerProperties()
    {
        fYoung = 0.;
        fPoisson = 0.;
        fSigmaMax = 0.;
        fSigmaMin = 0.;
        fSigmaConf = 0.;
        fTVDini = 0.;
        fTVDfin = 0.;
        fKIc = 0.;
        fCl = 0.;
        fPe = 0.;
        fgradPref = 0.;
        fvsp = 0.;
    }
    
    //Elastic 3D
    REAL fYoung;
    REAL fPoisson;
    REAL fSigmaMax;
    REAL fSigmaMin;
    REAL fSigmaConf;
    
    //TVD limits
    REAL fTVDini;
    REAL fTVDfin;
    
    //SIF
    REAL fKIc;
    
    //leafoff
    REAL fCl;
    REAL fPe;
    REAL fgradPref;
    REAL fvsp;
};

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
    
    int Index()
    {
        int elIndex = fElem2D->Index();
        return elIndex;
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




class TPZPlaneFractureMesh
{
	public:
    
    TPZPlaneFractureMesh();
	
	/**
	 * @brief Constructor
	 * @param layerVec [in] : vector of layers
	 * @param bulletTVDIni [in] : bullets perforation initial (TVD) depth
	 * @param bulletTVDFin [in] : bullets perforation final (TVD) depth
     * @param xLength [in] : Reservoir length in x direction (crack propagation direction)
     * @param yLength [in] : Reservoir length in y direction (tickness of reservoir that couple fracture plane)
     * @param Lmax    [in] : Maximum element edge length
     * @param nstripes [in] : Amount of stripes in Y direction for applied pressure for reduced space
     *
     * TVD: True Vertical Depth (positive positions)
	 */
    TPZPlaneFractureMesh(TPZVec<TPZLayerProperties> & layerVec, REAL bulletTVDIni, REAL bulletTVDFin,
                         REAL xLength, REAL yLength, REAL Lmax, int nstripes);
    
	~TPZPlaneFractureMesh();
    
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
	void InitializeFractureGeoMesh(TPZVec<std::pair<REAL,REAL> > &poligonalChain);
    
    void WriteRefinedGeoMesh();
    void ReadRefinedGeomesh();
    
    TPZCompMesh * GetFractureCompMesh(int porder);
    
    TPZCompMeshReferred * GetFractureCompMeshReferred(TPZCompMesh * cmeshRef, int porder);
    
    TPZCompMesh * GetPressureCompMesh(REAL Qinj, int porder);
    
    TPZCompMesh * GetMultiPhysicsCompMesh(TPZVec<TPZCompMesh *> & meshvec, REAL Qinj, REAL visc, int porder);
    
    void SetSigmaNStripeNum(TPZCompMesh * cmeshref, int actStripe);
    
    int NStripes();
    
    void SetActualState();
    
    void SetPastState();
    
    int NCrackTipElements();
    
    TPZGeoEl * GetCrackTipGeoElement(int pos);
    
    REAL GetKIcFromLayerOfThisZcoord(REAL zCoord);
    
    void GetYoung_and_PoissonFromLayerOfThisZcoord(REAL zCoord, REAL & young, REAL & poisson);
    
    REAL GetPreStressYYOfThisLayer(int lay);
    
    bool GeoElementIsOnPreservedMesh(TPZGeoEl * gel);
    
protected:
    
    //** Returns the layer index based on zMed (will be used fabs(zMed), i.e., TVD) */
    int GetLayer(REAL zMed);
    
    /** @brief Generation of the persistent full mesh (2D and 3D) that contains the fracture and its porous media
     *  @note This method set the fPreservedMesh atribute that will not be changed for every fracture time step
     *  @param espacamentoVerticalDEPTH [in] : espacamento vertical que define interfaces entre camadas horizontais
     *  @param bulletTVDIni [in] : bullets perforation initial (TVD) depth
	 *  @param bulletTVDFin [in] : bullets perforation final (TVD) depth
     *  @param xLength [in] : Reservoir length in x direction (crack propagation direction)
     *  @param yLength [in] : Reservoir length in y direction (tickness that couple fracture plane)
     */
    void GeneratePreservedMesh(std::list<REAL> & espacamentoVerticalDEPTH,
                               REAL bulletTVDIni, REAL bulletTVDFin,
                               REAL xLength, REAL yLength);
    
    void RefineUniformAllFracturePlane(int ndiv);
    
    void RefineDirectionalToCrackTip(int ndiv);

    /** @brief Method used for the mesh generator methods GeneratePlaneMesh and GeneratePreservedMesh
     *  @note For a given xz plane (defined by Y coordinate), generate the node grid coordinates
     *  @param espacamentoVerticalDEPTH [in] : espacamento vertical que define interfaces entre camadas horizontais
     *  @param xLength [in] : Reservoir length in x direction (crack propagation direction)
     */
    void GenerateNodesAtPlaneY(std::list<REAL> & espacamentoVerticalDEPTH, REAL xLength,
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
	 * @param auxElIndex_TrimCoords [out] : map that contains 1D element Index and a set of it trim 1D coordinates
	 * @param auxElIndexSequence [out] : the same of auxElIndex_TrimCoords, but keeps the trim 1D coordinates in generation sequence order
	 */
	void DetectEdgesCrossed(TPZVec<std::pair<REAL,REAL> > &poligonalChain, TPZGeoMesh * fractMesh,
							std::map< long, std::set<REAL> > &auxElIndex_TrimCoords, std::list< std::pair<long,REAL> > &auxElIndexSequence);

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
	 * @param auxElIndex_TrimCoords [out] : map that contains the trim coordinates of 1D element, indexed by its Id (once 1D element was inserted in gel->Mesh)\n
	 *							    obs.: auxElIndex_TrimCoords was idealized to work in accumulative conception, i.e.,
	 *								    each time this method is called, this map grows!
	 * @param auxElIndexSequence [out] : the same of auxElIndex_TrimCoords, but keeps the trim 1D coordinates in generation sequence order
	 * @param pushback [in] : set if dots on element edges will be inserted at the end, or not (i.e.: at beggining), of fCrackBoundary list
	 */
	static TPZGeoEl * CrossToNextNeighbour(TPZGeoEl * gel, TPZVec<REAL> &x, TPZVec<REAL> dx, REAL alphaMin, std::map< long,
									std::set<REAL> > &auxElIndex_TrimCoords, std::list< std::pair<long,REAL> > &auxElIndexSequence, bool pushback,
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
	static void UpdatePoligonalChain(TPZGeoMesh * gmesh, std::list< std::pair<long,REAL> > &auxElIndexSequence,
							  TPZVec<std::pair<REAL,REAL> > &poligonalChainUpdated);
	
	/**
	 * @param gmesh2D geometric mesh bi-dimensional
	 * @param gmesh3D geometric mesh three-dimensional
	 * @param auxElIndexSequence - output data: list that contains 1D element index and it trim 1D coordinates in generation sequence order
	 */
	void GenerateCrackBoundary(TPZGeoMesh * gmesh3D,
                               std::list< std::pair<long,REAL> > &auxElIndexSequence);
    
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

    //** just for visualize given dots in vtk */
    static void InsertDots4VTK(TPZGeoMesh * gmesh, const TPZVec<REAL> &fractureDots);
	
    //-------------------------------------------------------------------------------
	
    /** @brief Original 3D mesh (keeped intact for any poligonalChain configuration) */
	TPZGeoMesh * fPreservedMesh;

    /** @brief Refined 3D mesh (preserved mesh after refinement provided from poligonalChain configuration) */
    TPZGeoMesh * fRefinedMesh;
    
    /** @brief Vector of layers struct */
    TPZVec<TPZLayerProperties> fLayerVec;
    
    /** @brief 1D elements Indexes that compose crack boundary */
    TPZVec<long> fcrackBoundaryElementsIndexes;
    
    /** @brief initial element index in point search (just for optimization) */
    long fInitialElIndex;
    
    /** @brief maximum element edge length */
    REAL fLmax;
    
    /** @brief maximum X coordinate of poligonal chain (is the max Lfracture) */
    REAL fLfrac;
    
    /** @brief Amount of stripes of pressure */
    int fnstripes;
    
    /** @brief Maximum compressive stress in fracture plane */
    REAL fmaxCompressiveStress;
    
    /** @brief Vector of coupling material */
    TPZVec<TPZPlaneFractCouplingMat *> fCouplingMatVec;
};

#endif