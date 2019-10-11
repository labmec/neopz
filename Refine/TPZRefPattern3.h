/**
 * @file
 * @brief Contains the TPZRefPattern class which defines the topology of the current refinement pattern to a mesh.
 */
#ifndef PZREFPATTERN3H //@TODOFran: remove the 3
#define PZREFPATTERN3H

#include "pzvec.h"
#include "pztrnsform.h"
#include "pzgeoelside.h"
#include "pzgmesh.h"
#include "tpzpermutation.h"

class TPZRefPattern;
// #include <iostream>
// #include <string>
// #include <map>
// #include <list>
// #include <set>
// #include "TPZRefPattern3DataBase.h"
// #include "TPZRefPattern3Tools.h"

// class TPZGeoNode;
// class TPZCompMesh;
// //class TPZGeoElBC;
// class TPZGeoEl;
// class TPZGeoElSideIndex;
// class TPZGeoElSide;

// /** \addtogroup refine
//  @{
//  */
// /// Id for non initialized pattern
// const int nonInitializedId = -50;
// /// Name for non initialized pattern
// const std::string nonInitializedName = "noname";

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

/**
 * @brief Defines the topology of the current refinement pattern to a mesh. \ref refine "Refine"
 */
class TPZRefPattern3 : public TPZSavable {
protected:
    /// Id for identifying a non initialized pattern
    static const int fNonInitializedId;
    /// Name for identifying a non initialized pattern
    static const std::string fNonInitializedName;
	/**
	* @brief Data structure with sub element info
	* @details fSubElSideInfo[i][j] contains a pair (int,TPZTransform) related to the j-th side of the i-th son.
	* the pair consists of
	 * a) the side of the father in which the son is contained
	 * b) the parametric transformation to the father's coordinates
	*/
	TPZVec<TPZVec<std::pair<int,TPZTransform<REAL>>>> fSubElSideInfo;

    /**
    * @brief Data structure for storage of refinement information for a father's side
    */
	struct SPZFatherSideInfo{
	     /// a vector of nodes contained in the i-th side
        TPZVec<int> fSideNodes;
        /// a vector of TPZGeoElSide relative to its sons
        TPZVec<TPZGeoElSide> fSideSons;
        ///if there is a midsidenode in the side, it's index. otherwise is -1
        int64_t fMidSideIndex;
        SPZFatherSideInfo(const TPZVec<int> &sideNodes,
                const TPZVec<TPZGeoElSide> &sideSons,
                const int64_t &midSideIndex) :
        fSideNodes(sideNodes), fSideSons(sideSons), fMidSideIndex(midSideIndex) {}
        SPZFatherSideInfo() = default;
        SPZFatherSideInfo(const SPZFatherSideInfo &) = default;
        ~SPZFatherSideInfo() = default;
	};
	/**
	* @brief Data structure for storage of refinement information for each father's side
	* @detail fFatherSideInfo[i] will contain a tuple made of three entities, stored in SPZFatherSideInfo.
	*/
	TPZVec<SPZFatherSideInfo> fFatherSideInfo;



	/////////////other stuff///////////////
	/// name of the refinement pattern
	std::string fName;

	/// id of the refinement pattern
	int fId;

    /**
    * @brief Refinement for elements associated with a certain side of the current element.
    * Usually, this is used for boundary/interface elements.
    */
    TPZVec<int> fSideRefPattern;

    /**
    * @brief Vector containing refinement patterns for each permutation of the master element.
    *
    * The vector stores the id correspondent to the refinement pattern vector in fRefPatternMesh
    */
    TPZVec<int> fPermutedRefPatterns;

    /**
    * @brief Geometric mesh which defines the topology of the current refinement pattern
    */
    TPZGeoMesh fRefPatternMesh;

    /**
     * @brief Number of subelements. Must be available before the mesh is generated.
     */
    int fNSubEl;//@TODOFran:: is it really necessary?

    /**
     * Given a certain geometric element corresponding to a sub-element, this method will return its index
     * @return index of the sub-element
     */
    int FindSubEl(TPZGeoEl *) const;

    /** @brief Auxiliar structure to permute nodes */
    class TPZRefPattern3Permute : public TPZSavable {
    public:
        /** @brief permutation of the nodes */
        TPZPermutation fPermute;

        /** @brief Transformation to the nodes */
        TPZTransform<> fTransform;
        /** @brief Default constructor */
        TPZRefPattern3Permute(): fPermute(0)
        {}
        /** @brief Copy constructor */
        TPZRefPattern3Permute(const TPZRefPattern3Permute &copy) : fPermute(copy.fPermute), fTransform(copy.fTransform)
        {
        }
        /** Assignment operator */
        TPZRefPattern3Permute &operator=(const TPZRefPattern3Permute &copy)
        {
            fPermute = copy.fPermute;
            fTransform = copy.fTransform;
            return *this;
        }
        int ClassId() const override;
        void Read(TPZStream &buf, void *context) override;
        void Write(TPZStream &buf, int withclassid) const override;
    };

    /**
     * @brief Map of all valid permutations
     */
    static std::map<MElementType, std::list<TPZRefPattern3Permute> > fPermutations;


    friend class
    TPZRefPatternDataBase;

    friend class
    TPZRefPatternTools;

public:

    TPZRefPattern3();
//    /**
//    * @brief Constructor whose argument is the name of the file with the definition
//    * of the refinement pattern
//    */
//    TPZRefPattern3(std::istream &file);//@TODOFran: Implement me!
//
//    /**
//    * @brief Constructor whose argument is a string containing the definition
//    * of the refinement pattern
//    */
//    TPZRefPattern3(const std::string &file);//@TODOFran: Implement me!
//
//    /**
//    * @brief Creates an TPZRefPattern3 from a given mesh
//    */
//    TPZRefPattern3(TPZGeoMesh &RefPatternMesh);//@TODOFran: Implement me!
//
//    /**
//    * @brief Create a copy of the TPZRefPattern3 applying the permutation on the first element
//    */
//    TPZRefPattern3(const TPZRefPattern3 &copy, const TPZPermutation &permute);//@TODOFran: Implement me!

// 	int operator==(const TPZAutoPointer<TPZRefPattern3> compare) const; //@TODOFran: Implement me!
    /**
 	* @brief Copy constructor
 	*/
    TPZRefPattern3(const TPZRefPattern3 &copy) = default;

	explicit TPZRefPattern3(TPZRefPattern &oldRef);

	virtual ~TPZRefPattern3() = default;

    //@TODOFran: Document me!
    void PrintMore(std::ostream &out = std::cout) const;

    /**
     * @brief Sets the name associated with the refinement pattern
     * @param name a string containing the name to be set
     */
 	void SetName(std::string name)
 	{
 		fName = name;
 	}
    int ClassId() const override;

//    void Read(TPZStream &buf, void *context) override;//@TODOFran: Implement me!
//
//    void Write(TPZStream &buf, int withclassid) const override;//@TODOFran: Implement me!

    /**
    * Gives the number of sub elements of this refinement pattern
    * @return number of sub elements
    */
    int NSubElements() const;
    /**
    * @brief Returns the father side associated with a given side of a certain sub element
    * @param side side of the sub element
    * @param sub index of the sub element
    * @return father's side
    */
    int FatherSide(int side, int sub) const;


    /**
    * @brief Returns the number of TPZGeoElSides associated with a father's side.
    */
    int NSideSubGeoElSides(int fatherSide) const;
    /**
     * @brief Gives information about the ith subelement contained in a given father's side
     * @param fatherSide father's side from which we want to know information about a sub element
     * @param subElPos the information regards the subElPos-th sub element of the fatherSide
     * @param subGeoEl TPZGeoElSide of the (sub-element,side) associated with the father's side
     */
    void SideSubGeoElSide(int fatherSide, int subElPos, TPZGeoElSide & subGeoEl) const;
	
    /**
    * @brief It returns the TPZTransform associated with a certain sub-element of a father's side
    * @param fatherSide father's side from which we want to know information about a sub element
    * @param sub the information regards the sub-th sub element of the fatherSide
    */
    TPZTransform<> Transform(int fatherSide, int sub);
	
     /**
      * @brief It returns a TPZVec containing the nodes contained in a given father's side
      * @param fatherSide father's side from which we want to know information about a sub element
      * @param vecNodes vector containing the nodes of side fatherSide
      */
     void SideNodes(int fatherSide, TPZVec<int> &vecNodes);
	
     /**
      * @brief Returns the number of internal nodes of side
      */
     int NSideNodes(int fatherSide);
	
//     /**
//      * @brief Returns the number of nodes of the mesh
//      */
//     int NNodes();
//     /**
//      * @brief It prints the features of the standard of geometric refinement.
//      */
//     void MeshPrint();
	
// 	/**
// 	 * @brief Prints the useful information of a Refinement Pattern in a ostream file
//      */
//     void Print(std::ostream &out = std::cout);
    
//     void ShortPrint(std::ostream &out = std::cout);
	
//     /**
//      * @brief It accumulates the number of sides of son 0 until son ison.
//      */
//     int SizeOfSubsSides(int ison);
	
//     /**
//      * @brief It defines the size of the partition.
//      */
//     void DefinitionOfSizePartition();
	
//     /**
//      * @brief It is verified the son is neighboring of the father.
//      * 
// 	 * The side that corresponds to the son is returned.
//      */
//     int IsFatherNeighbour(TPZGeoElSide fathside,TPZGeoEl *son);
	
//     /**
//      * @brief It returns the element number iel from the stack of elements of the
//      * geometric mesh
//      */
//     TPZGeoEl *Element(int iel);
	
//     /** 
//      * @brief It compares two hashings: in case that are equal returns 0,
//      * case is distinct returns 1
//      */
//     int IsNotEqual(TPZTransform<> &Told, TPZTransform<> &Tnew);
	
// 	TPZAutoPointer<TPZRefPattern3> SideRefPattern(int side);
    
//     /**
//      * @brief Find the side refinement pattern corresponding to the parameter transformation
//      */
//     TPZAutoPointer<TPZRefPattern3> SideRefPattern(int side, TPZTransform<> &trans);
	
// 	/**
// 	 * @brief This method is used to create / identify the nodes of the refined elements.
// 	 * 
// 	 * The method verify if the nodes are already created by the self element or by some neighbour.
// 	 * @param gel - pointer to the element which are being divided
// 	 * @param newnodeindexes - return all midside node indexes for the element division.
// 	 */
// 	void CreateNewNodes(TPZGeoEl *gel, TPZVec<int64_t> &newnodeindexes);
	
// 	/**
// 	 * @brief This method is used to create / identify the midside nodes for element elindex in its division process.
// 	 * @param gel - pointer to the element which are being divided
// 	 * @param side - Side along which the nodes will be identified/created
// 	 * @param newnodeindexes - return all midside node indexes for the element division.
// 	 */
// 	/**
// 	 * The method verify if the nodes are already created by the self element or by some neighbour.
// 	 */
// 	void CreateMidSideNodes(TPZGeoEl *gel, int side, TPZVec<int64_t> &newnodeindexes);
	
// 	/**
// 	 * @brief Returns the refinement pattern identifier
// 	 */
// 	std::string Name();
	
//     /**
//      * @brief Generate all permuted partitions and insert them in the mesh
//      */
// 	void InsertPermuted();
    
//     /**
// 	 * @brief Generate the refinement patterns associated with the sides of the father element
// 	 */
//     void GenerateSideRefPatterns();
	
//     /**
//      * @brief Find the refinement pattern corresponding to the give transformation
//      */
// 	TPZAutoPointer<TPZRefPattern3> FindRefPattern(TPZTransform<> &trans);
	
// 	/**
// 	 * @brief Return the id of the refinement pattern
// 	 */
// 	const int Id() const
// 	{
// 		return fId;
// 	}
	
// 	bool IdInitialized()
// 	{
// 		return (fId != nonInitializedId);
// 	}
	
// 	/**
// 	 * @brief Set the id of the refinement pattern
// 	 */
// 	void SetId(int id)
// 	{
// 		fId = id;
		
// 		TPZGeoEl * el = fRefPatternMesh.ElementVec()[0];
// 		int nSides = el->NSides();
// 		fSideRefPattern[nSides-1] = id;
// 	}
	
// 	*
// 	 * @brief Fill a vector with TPZGeoElSideIndex_s  with respect to internal
// 	 * subelements of zero dimension of a given side
// 	 * @param side - side of father element that will be searched for internal nodes (generated from subelements)
// 	 * @param nodeIndexes - vector of subelements index/Sides (stored as TPZGeoElSideIndex objects that have dimension = 0) that
// 	 * belongs to the interior of given side of father element.
	 
// 	void InternalNodesIndexes(int side, TPZVec<TPZGeoElSideIndex> &nodeIndexes);
	
// 	/**
// 	 * @brief Fill a vector with TPZGeoElSideIndex_s  with respect to internal subelements of a given side
// 	 * @param side side of father element that will be searched for internal sides (generated from subelements)
// 	 * @param sideIndexes vector of subelement index/sides (stored as TPZGeoElSideIndex objects) that belongs to the interior of given side of father element.
// 	 */
// 	void InternalSidesIndexes(int side, TPZVec<TPZGeoElSideIndex> &sideIndexes);
	
// 	TPZAutoPointer<TPZRefPattern3> GetPermutation(int pos);
	
// 	TPZGeoMesh &RefPatternMesh()
// 	{
// 		return fRefPatternMesh;
// 	}
	
// 	void PrintVTK(std::ofstream &file, bool matColor = false);
	
// private:
	
// 	/**
//      * @brief This method import a refpattern from a given *.rpt file
//      */
//     void ImportPattern(std::istream &in);
	
// 	/**
//      * @brief This method export a refpattern to a given *.rpt file
//      */
// 	void ExportPattern(std::ostream &out);
	
// 	void ReadPattern(std::istream &in);
	
// 	/**
// 	 * @brief Write out the Refinement Pattern to a file with the necessary  
//      * information to read it with ReadPattern.
//      */
//     void WritePattern(std::ofstream &out);
	
// 	/**
//      * @brief It returns a stack from sub-elements and its sides. This stack partitions side of the element father.
// 	 * 
// 	 * It returns the number from elements of the partition of the side
// 	 * @param side - side of father element that will be searched for internal sides (generated from subelements)
// 	 * @param gelvec - vector of sides (stored as TPZGeoElSide objects) that belongs to given side of father element.
//      */
//     int SidePartition(TPZVec<TPZGeoElSide> &gelvec, int side);
	
// 	/**
// 	 * @brief Sets the RefPatternMesh in (x,y,z)_coordinates to (qsi,eta,zeta)_coordinates, always respecting the R3 dimension.
// 	 */
// 	void SetRefPatternMeshToMasterDomain();
	
//     /**
//      * @brief Calculates the hashings between the sides of the son and the father 
//      */
//     void ComputeTransforms();
	
//     /**
//      * @brief It effects the partition of the sides of the
// 	 * element father using the sides of the children
//      */
//     void ComputePartition();
	
	
// 	/**
//      * @brief Each side of the element father is gotten as a partition enters the
//      * sides of its sub-elements.
// 	 */
// 	/** 
// 	 * The vector fTransformSides keeps
//      * the respective hashing enters the side of the son and the side of \n
//      * the father who contains it.
//      */
//     struct TPZPartitionFatherSides : public TPZSavable {
//         /**
//          * @brief Vector of position in fPartitionElSide of the side of the element to
//          * be partitioned father
//          */
//         TPZManVector<int,27> fInitSide;
		
//         /**
//          * @brief Vector that contains the partition of each side of the element 
//          * 
// 		 * fFatherSides  means of sub-elements and its sides. As extension 
//          * the partition associated with a vertex corresponds to the on
//          * elements to this node
//          */
//         TPZManVector<TPZGeoElSideIndex,10> fPartitionSubSide;
		
//         /**
//          * @brief Number of asociados distinct sub-elements to the side of the father
//          */
//         TPZManVector<int,27> fNSubSideFather;
		
//         /**
//          * @brief It prints the properties of the structure
//          */
//         void Print(TPZGeoMesh &gmesh,std::ostream &out = std::cout);  
// 		        int ClassId() const override;
//         void Read(TPZStream &buf, void *context) override;
//         void Write(TPZStream &buf, int withclassid) const override;
//     };
	
//     /**
//      * @brief This structure is defined with the intention to know which is the side 
//      * of the element father who contains the side of the sub-element.
// 	 */
// 	/** 
// 	 * A filled time this information calculates it hashing enters the side of the sub-element \n
// 	 * and the side of the respective element father
//      */ 
//     class TPZSideTransform : public TPZSavable {
//     public:
//         /**
//          * @brief Vector of position of fSideFather
//          */
//         TPZManVector<int> fInitSonSides;
		
//         /**
//          * @brief Side of the element father associated with the side of the
//          * sub-element
//          */
//         TPZManVector<int> fFatherSide;
		
//         /**
//          * @brief Hashing enters the side of the sub-element and the side of the
//          * corresponding father
//          */
//         TPZManVector<TPZTransform<> > fSideTransform;
		
//         /**
//          * @brief It prints the properties of the structure
//          */
//         void Print(TPZGeoMesh &gmesh,std::ostream &out = std::cout);
//                 int ClassId() const override;
//         void Read(TPZStream &buf, void *context) override;
//         void Write(TPZStream &buf, int withclassid) const override;
//     };
	
// 	/**
//      * @brief Partition of the element for objects sub-element-sides 
//      */
// 	TPZPartitionFatherSides fFatherSides;
	
// 	/**
//      * @brief Partition of the element father in sides of sub-elements and
//      * the respective hashing between these.
//      */    
//     TPZSideTransform fTransforms;
	
// public:   
	
//     /**
// 	 * @brief Automatically generate all permuted numberings for the father element of RefPatternMesh
// 	 */
//     static void GeneratePermutations(TPZGeoEl *gel);
    
// protected: 
	
//     /**
//      * @brief Copy the mesh structure applying the permutation on the nodenumbers of the first element
//      */
// 	void PermuteMesh(const TPZPermutation &permute);
	
// 	/**
// 	 * @brief Build a geometric mesh associated with the side of the refinement pattern
// 	 */
// 	void BuildSideMesh(int side, TPZGeoMesh &SideRefPatternMesh);
	
// 	MElementType Type();
	
// public:
// 	/** @brief Auxiliar structure to permute nodes */
// 	class TPZRefPattern3Permute : public TPZSavable {
//             public:
// 		/** @brief permutation of the nodes */
// 		TPZPermutation fPermute;
		
// 		/** @brief Transformation to the nodes */
// 		TPZTransform<> fTransform;
// 		/** @brief Default constructor */
// 		TPZRefPattern3Permute(): fPermute(0)
// 		{}
// 		/** @brief Copy constructor */
// 		TPZRefPattern3Permute(const TPZRefPattern3Permute &copy) : fPermute(copy.fPermute), fTransform(copy.fTransform)
// 		{
// 		}
// 		/** Assignment operator */
// 		TPZRefPattern3Permute &operator=(const TPZRefPattern3Permute &copy)
// 		{
// 			fPermute = copy.fPermute;
// 			fTransform = copy.fTransform;
// 			return *this;
// 		}
// 		                int ClassId() const override;
//                 void Read(TPZStream &buf, void *context) override;
//                 void Write(TPZStream &buf, int withclassid) const override;
// 	};
};

/** @} */

#endif
