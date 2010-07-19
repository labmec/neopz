#ifndef PZREFPATTERNH
#define PZREFPATTERNH

#include "pzvec.h"
#include "pzgmesh.h"
#include <iostream>
#include <string>
#include <map>
#include <list>
#include <set>
#include "tpzpermutation.h"
#include "pztrnsform.h"
#include "TPZRefPatternDataBase.h"
#include "TPZRefPatternTools.h"

class TPZGeoNode;
class TPZCompMesh;
class TPZGeoElBC;
class TPZGeoEl;
class TPZGeoElSideIndex;
class TPZGeoElSide;

const int nonInitializedId = -50;
const std::string nonInitializedName = "noname";

/**
 * ***************************************** !!! @@@ READ ME @@@ !!! *****************************************
 *
 * The archive to be imported (of an Refinement Pattern) reffers to an geoMesh that contains the father element and its sub-elements.
 * Its format data is the following one:
 *
 * begin of file ....................................
 *
 * //block 1
 * int(nNodes)	int(nElements)
 * int(Id)	string(refName)
 *
 * //block 2 - [the line above repeat for each node]
 * double(nodeCoordinates)
 * 
 * //block 3 - [the line above repeat for each element]
 * int(elType)		int(matId)		int(topologySequence)
 * 
 * ................................................................. end of file
 *
 * Obs.:
 *
 * - The refName of refinement pattern must be "noname" if the refpattern is NOT regular (nodes in the middle of edges or c.g. of faces/volume)
 *
 * - The nodeCoordinates of block 2 MUST be in R3 dimension (x,y,z);
 *
 * - The sequence described for the nodeCoordinates will assumption that the first will be the node of Id=0, second node of Id=1 etc.
 *   This convention will be used to describe the topologySequence of elements;
 *
 * - The first element described on the block 3 will be considered as the father element of RefPatternMesh.
 *   The following ones will be its sons (sub-elements);
 *
 * - The values of elType, can be:  line(1), triangle(2), quadrilateral(3), tetrahedron(4), pyramid(5), prism(6), hexaedron(7).
 *
 * ***************************************** !!! @@@ THANKS @@@ !!! *****************************************
 */

class TPZRefPattern
{
	
public:
	
	TPZRefPattern();
	
    /**
     * Constructor whose argument is the name of the file with the definition
     * of the refinement standard
     */
	TPZRefPattern(std::istream &file);
	
	TPZRefPattern(const std::string &file);
	
	/**
     * Creates an TPZRefPattern from a given mesh
     */
    TPZRefPattern(TPZGeoMesh &RefPatternMesh);
    
    /**
	 * Copy constructor
	 */
    TPZRefPattern(const TPZRefPattern &copy);
    
    /**
     * Create a copy of the TPZRefPattern applying the permutation on the first element
     */
    TPZRefPattern(const TPZRefPattern &copy, const TPZPermutation &permute);
	
    /**
     * Destructor of the object
     */
    ~TPZRefPattern()
	{
	}
	
	int operator==(const TPZAutoPointer<TPZRefPattern> compare) const;
	
	void BuildName();
	
	bool NameInitialized()
	{
		return (fName != nonInitializedName);
	}
	
	void SetName(std::string name)
	{
		fName = name;
	}

    void Read(TPZStream &buf);
	
    void Write(TPZStream &buf);
	
	/**
     * Sides associates to the element father
     */
    int FatherSide(int side, int sub);
	
    /**
     * It returns the sub-element and side from the
	 * sub-element number position of side side 
     */
    void SideSubElement(int sidein, int position, int & sub, int & sideout);
		
    /**
     * It returns the hashing enters the side of an sub-element and the side of
     * the father who contains it.
     */
    TPZTransform Transform(int side, int sub);
	
    /**
     * It returns the stack from referring indices
	 * to the internal nodes to the side
     */
    void SideNodes(int side, TPZVec<int> &vecnodes);
	
    /**
     * It returns the number from internal nodes of side 
     */
    int NSideNodes(int side);
	
    /**
     * It returns the number of nodes of the mesh
     */
    int NNodes();
	
    /**
     * It returns the number from sub-elements of the division.
     */
    int NSubElements();
	
    /**
     * It fills the structure of data that determines the number of elements
     *  associates to the side of the element father
     */
    void NSideSubElements();
	
    /**
     * It returns the number from sub-elements associates to side.
     */
    int NSideSubElements(int side);
	
    /**
     * It prints the features of the standard of geometric refinement.
     */
    void MeshPrint();
    void PrintMore(std::ostream &out = std::cout);
	
	/**
	 * Print the useful information of a Refinement Pattern in a ostream file
     */
    void Print(std::ostream &out = std::cout);
    
    void ShortPrint(std::ostream &out = std::cout);
	
    /**
     * It accumulates the number of sides of son 0 until son ison.
     */
    int SizeOfSubsSides(int ison);
	
    /**
     * It defines the size of the partition.
     */
    void DefinitionOfSizePartition();
	
    /**
     * It is verified the son is neighboring of the father.
     * The side that corresponds to the son is returned.
     */
    int IsFatherNeighbour(TPZGeoElSide fathside,TPZGeoEl *son);
	
    /**
     * It returns the element number iel from the stack of elements of the
     * geometric mesh
     */
    TPZGeoEl *Element(int iel);
	
    /** 
     * It compares two hashings: in case that are equal returns 0,
     * case is distinct returns 1
     */
    int IsNotEqual(TPZTransform &Told, TPZTransform &Tnew);
	
	TPZAutoPointer<TPZRefPattern> SideRefPattern(int side);
    
    /**
     * Find the side refinement pattern corresponding to the parameter transformation
     */
    TPZAutoPointer<TPZRefPattern> SideRefPattern(int side, TPZTransform &trans);
	
	/**
	 * This method is used to create / identify the nodes of the refined elements.
	 * The method verify if the nodes are already created by the self element or by some neighbour.
	 * @param gel - pointer to the element which are being divided
	 * @param newnodeindexes - return all midside node indexes for the element division.
	 */
	void CreateNewNodes(TPZGeoEl *gel, TPZVec<int> &newnodeindexes);
	
	/**
	 * This method is used to create / identify the midside nodes for element elindex in its division process.
	 * The method verify if the nodes are already created by the self element or by some neighbour.
	 * @param gel - pointer to the element which are being divided
	 * @param side - Side along which the nodes will be identified/created
	 * @param newnodeindexes - return all midside node indexes for the element division.
	 */
	void CreateMidSideNodes(TPZGeoEl *gel, int side, TPZVec<int> &newnodeindexes);
	
	/**
	 * Returns the refinement pattern identifier
	 */
	std::string Name();
	
    /**
     *  Generate all permuted partitions and insert them in the mesh
     */
	void InsertPermuted();
    
    /**
	 * Generate the refinement patterns associated with the sides of the father element
	 */
    void GenerateSideRefPatterns();
	
    /**
     * Find the refinement pattern corresponding to the give transformation
     */
	TPZAutoPointer<TPZRefPattern> FindRefPattern(TPZTransform &trans);
	
	/**
	 * Return the id of the refinement pattern
	 */
	int Id()
	{
		return fId;
	}
	
	bool IdInitialized()
	{
		return (fId != nonInitializedId);
	}
	
	/**
	 * Set the id of the refinement pattern
	 */
	void SetId(int id)
	{
		fId = id;
		
		TPZGeoEl * el = fRefPatternMesh.ElementVec()[0];
		int nSides = el->NSides();
		fSideRefPattern[nSides-1] = id;
	}
	
	/**
	 * Fill a vector with TPZGeoElSideIndex_s  with respect to internal
	 * subelements of zero dimension of a given side
	 * @param side - side of father element that will be searched for internal nodes (generated from subelements)
	 * @param nodeIndexes - vector of subelements index/Sides (stored as TPZGeoElSideIndex objects that have dimension = 0) that
	 * belongs to the interior of given side of father element.
	 */
	void InternalNodesIndexes(int side, TPZVec<TPZGeoElSideIndex> &nodeIndexes);
	
	/**
	 * Fill a vector with TPZGeoElSideIndex_s  with respect to internal subelements of a given side
	 * @param side - side of father element that will be searched for internal sides (generated from subelements)
	 * @param nodeIndexes - vector of subelement index/sides (stored as TPZGeoElSideIndex objects) that belongs to the interior of given side of father element.
	 */
	void InternalSidesIndexes(int side, TPZVec<TPZGeoElSideIndex> &sideIndexes);
	
	TPZAutoPointer<TPZRefPattern> GetPermutation(int pos);
	
	TPZGeoMesh &RefPatternMesh()
	{
		return fRefPatternMesh;
	}
	
	void PrintVTK(std::ofstream &file);
	
private:
	
	/**
     * This method import a refpattern from a given *.rpt file
     */
    void ImportPattern(std::istream &in);
	
	/**
     * This method export a refpattern to a given *.rpt file
     */
	void ExportPattern(std::ostream &out);
	
	void ReadPattern(std::istream &in);
	
	/**
	 * Write out the Refinement Pattern to a file with the necessary  
     * information to read it with ReadPattern.
     */
    void WritePattern(std::ofstream &out);
	
	/**
     * It returns a stack from sub-elements and its sides. This stack partitions side of the element father.
	 * It returns the number from elements of the partition of the side
	 * @param side - side of father element that will be searched for internal sides (generated from subelements)
	 * @param gelvec - vector of sides (stored as TPZGeoElSide objects) that belongs to given side of father element.
     */
    int SidePartition(TPZVec<TPZGeoElSide> &gelvec, int side);
	
	/**
	 * Set the RefPatternMesh in (x,y,z)_coordinates to (qsi,eta,zeta)_coordinates, always respecting the R3 dimension.
	 */
	void SetRefPatternMeshToMasterDomain();
	
    /**
     * It calculates the hashings between the sides of the son and the father 
     */
    void ComputeTransforms();
	
    /**
     * It effects the partition of the sides of the
	 * element father using the sides of the children
     */
    void ComputePartition();
	
    /**
     * Geometric mesh which defines the topology of the current refinement pattern
     */
    TPZGeoMesh fRefPatternMesh;
	
    /**
	 * Refinamento para elementos de contorno associados ao refinamento atual
	 */
    TPZVec<int> fSideRefPattern;
	
	/**
	 * vector of refinement patterns for each permutation of the master element.
	 * The vector stores the id correspondent to the refinement pattern vector in fOwnerMesh
	 */
	std::vector<int> fPermutedRefPatterns;
    
    /**
	 * This should be available before the mesh initialization
	 */
    int fNSubEl;
    
    /**
     * Unique id given to the refpattern in order to read and write from disk
     */
	int fId;
	
	/**
     * Each side of the element father is gotten as a partition enters the
     * sides of its sub-elements.  The vector fTransformSides keeps
     * the respective hashing enters the side of the son and the side of 
     * the father who contains it.
     */
    struct TPZPartitionFatherSides
	{
        /**
         * Vector of position in fPartitionElSide of the side of the element to
         * be partitioned father
         */
        TPZVec<int> fInitSide;
		
        /**
         * Vector that contains the partition of each side of the element 
         * fFatherSides  means of sub-elements and its sides. As extension 
         * the partition associated with a vertex corresponds to the on
         * elements to this node
         */
        TPZVec<TPZGeoElSideIndex> fPartitionSubSide;
		
        /**
         * Number of asociados distinct sub-elements to the side of the father
         */
        TPZVec<int> fNSubSideFather;
		
        /**
         * It prints the properties of the structure
         */
        void Print(TPZGeoMesh &gmesh,std::ostream &out = std::cout);  
		
        void Read(TPZStream &buf);
        void Write(TPZStream &buf);
    };
	
    /**
     * This structure is defined with the intention to know which is the side 
     * of the element father who contains the side of the sub-element. A filled
     * time this information calculates it hashing enters the side of the sub-
     * element and the side of the respective element father
     */ 
    struct TPZSideTransform
	{
        /**
         * Vector of position of fSideFather
         */
        TPZVec<int> fInitSonSides;
		
        /**
         * Side of the element father associated with the side of the
         * sub-element
         */
        TPZVec<int> fFatherSide;
		
        /**
         * Hashing enters the side of the sub-element and the side of the
         * corresponding father
         */
        TPZVec<TPZTransform> fSideTransform;
		
        /**
         * It prints the properties of the structure
         */
        void Print(TPZGeoMesh &gmesh,std::ostream &out = std::cout);
        
        void Read(TPZStream &buf);
        
        void Write(TPZStream &buf);
    };

	/**
     * Partition of the element for objects sub-element-sides 
     */
	TPZPartitionFatherSides fFatherSides;
	
	/**
     * Partition of the element father in sides of sub-elements and
     * the respective hashing between these.
     */    
    TPZSideTransform fTransforms;
	
public:   
	
    /**
	 *  Automatically generate all permuted numberings for the father element of RefPatternMesh
	 */
    static void GeneratePermutations(TPZGeoEl *gel);
    
	protected: 
	
    /**
     * Copy the mesh structure applying the permutation on the nodenumbers of the first element
     */
	void PermuteMesh(const TPZPermutation &permute);
	
	/**
	 *  Build a geometric mesh associated with the side of the refinement pattern
	 */
	void BuildSideMesh(int side, TPZGeoMesh &SideRefPatternMesh);
	
	MElementType Type();
	
public:
	struct TPZRefPatternPermute
	{
		// permutation of the nodes
		TPZPermutation fPermute;
		
		// transformation to the nodes
		TPZTransform fTransform;
		
		TPZRefPatternPermute(): fPermute(0)
		{}

		TPZRefPatternPermute(const TPZRefPatternPermute &copy) : fPermute(copy.fPermute), fTransform(copy.fTransform)
		{
		}

		TPZRefPatternPermute &operator=(const TPZRefPatternPermute &copy)
		{
			fPermute = copy.fPermute;
			fTransform = copy.fTransform;
			return *this;
		}

		void Read(TPZStream &buf)
		{
			this->fPermute.Read(buf);
			this->fTransform.Read(buf);
		}

		void Write(TPZStream &buf)
		{
			this->fPermute.Write(buf);
			this->fTransform.Write(buf);
		}
	};
	
		//protected:
	
	/**
	 * Map of all valid permutations
	 */
	static std::map<MElementType, std::list<TPZRefPatternPermute> > fPermutations;
	
	std::string fName;
	
	friend class 
	TPZRefPatternDataBase;
	
	friend class
	TPZRefPatternTools;
};

#endif
