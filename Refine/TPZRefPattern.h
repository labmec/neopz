#ifndef PZREFPATTERNH
#define PZREFPATTERNH

#include "pzvec.h"
#include "pzgmesh.h"
#include <iostream>

class TPZGeoNode;
class TPZCompMesh;
class TPZGeoElBC;
class TPZGeoEl;
class TPZGeoElSide;
class TPZTransform;

/* template <class T> */
/* class TPZVec; */


/**
 * The format of the archive of to be read data is the following one:
 *  
 * nodes_number  (n) 1+sub-elements_number  (m=father+childrens)
 * refinement_type (integer) 
 * co-ordinated_0 co-ordinated_1 co-ordinated_2 (first node)
 * co-ordinated_0 co-ordinated_1 co-ordinatvoid NSideSubElements();ed_2 (second node)
 * ... ... ...
 * co-ordinated_0 co-ordinated_1 co-ordinated_2 (node n)
 * typeofel idmat node_i1 node_j1 node_k1 node-l1 (nodal sequence of the first sub-element)
 * typeofel idmat node_i2 node_j2 node_k2 node_l2 (nodal sequence of the second sub-element)
 * ... ... ... ...
 * typeofel idmat node_im node_jm node_km node_lm (nodal sequence of the sub-element m)
 *
 * The first element, element 0, must be the element father.
 *
 * The values of typeofel, type of element, can be:  line (1), triangle (3), 
 * quadrilateral (4), tetrahedron (7), pyramid (5), prism (6), hexaedro (8).
 * 
 * The variable idmat is an entire number that identifies the material.
 *
 * The coordinates of we who define the conectividades of the sub-elements 
 * must be given in the element master of the element Father. The words 
 * between parenteses correspond the mere comments.
 */

class TPZRefPattern {

public:

    /**Empty constructor*/
    TPZRefPattern();
    /**
     * Constructor whose argument is the name of the archive with the definition
     * of the refinement standard
     */
    TPZRefPattern(char * file);

    /**
     * Destructor of the object
     */
    ~TPZRefPattern(){if(fMesh) delete fMesh;}

    /**
     * It returns the mesh of refinement pattern
     */
    TPZGeoMesh *Mesh() {return fMesh;}

    /**
     * It returns the refinement type
     */
    int Type() {return fRefineType;}

    /**
     * It effects the reading of the archive that defines the refinement standard
     */
    void ReadPattern();

    /**
     * It create geometric elements
     */
     TPZGeoEl *CreateGeoEl(int ntype,int mat,TPZVec<int> &nodes,TPZGeoMesh *gmesh);
    
    /**
     * It calculates the hashings between the sides of the son and the father 
     */
    void ComputeTransforms();

    /**
     * It effects the partition of the sides of the element father using the
     * sides of the children
     */
    void ComputePartition();

    /**
     * It returns a stack from sub-elements and its sides. This stack 
     * partitions  side of the element father. It returns the number 
     * from elements of the partition of the side
     */
    int SidePartition(TPZVec<TPZGeoElSide> &gelvec,int side);

    /**
     * It returns the sub-element and side from the sub-element number position
     * of side side 
     */
    void SideSubElement(int sidein, int position, int & sub, int & sideout);

    /**
     * Sides associates to the element father
     */
    int FatherSide(int side, int sub );

    /**
     * It returns the hashing enters the side of an sub-element and the side of
     * the father who contains it.
     */
    TPZTransform Transform(int side, int sub);

    /**
     * It returns the stack from referring indices to the internal nodes to the
     * side
     */
    void SideNodes(int side, TPZVec<int> & vecnodes);

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
     * It returns the number from sub-elements associates to side.
     */
    int NSideSubElements(int side);

    /**
     * It fills the structure of data that determines the number of elements
     *  associates to the side of the element father
     */
    void NSideSubElements();

    /**
     * It prints the features of the standard of geometric refinement.
     */
    void MeshPrint();
    void Print1(TPZGeoMesh &gmesh,ostream &out = cout);/////////////////////////???????????????????

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

   /** Algorithm that evaluates the veracity of the hashings between sides
     * of the elements children and corresponding sides of the father. A
     * point p in the parametric space of the side of the sub-element is
     * overcome and is calculated it mentioned hashing getting pf point in
     * the element father. One calculates for p and pf the corresponding
     * deformed point. Itself the hashing is consistent the deformed point
     * must the same be.
     */
    void TransformationTest();

    /** 
     * It compares two hashings: in case that are equal returns 0,
     * case is distinct returns 1
     */
    int IsNotEqual(TPZTransform &Told, TPZTransform &Tnew);

    /**
     * Each side of the element father is gotten as a partition enters the
     * sides of its sub-elements.  The vector fTransformSides keeps
     * the respective hashing enters the side of the son and the side of 
     * the father who contains it.
     */
    struct TPZPartitionFatherSides{
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
        TPZVec<TPZGeoElSide> fPartitionSubSide;

        /**
         * Number of asociados distinct sub-elements to the side of the father
         */
        TPZVec<int> fNSubSideFather;

        /**
         * It prints the properties of the structure
         */
        void Print(TPZGeoMesh &gmesh,ostream &out = cout);        
    };

    /**
     * This structure is defined with the intention to know which is the side 
     * of the element father who contains the side of the sub-element. A filled
     * time this information calculates it hashing enters the side of the sub-
     * element and the side of the respective element father
     */ 
    struct TPZSideTransform{

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
        void Print(TPZGeoMesh &gmesh,ostream &out = cout);
    };


private:

    /**
     * Associated geometric mesh to the standard of current refinement 
     */
    TPZGeoMesh *fMesh;

    /**
     * Name of the file that defines the refinement standard
     */
    char * fFileRefPatt;

    /**
     * Indice or type of the current refinement 
     */
    int fRefineType;

    /**
     * Partition of the element for objects sub-element-sides 
     */
     TPZPartitionFatherSides fFatherSides;

    /**
     * Partition of the element father in sides of sub-elements and
     * the respective hashing between these.
     */    
    TPZSideTransform fTransforms;

    /**refinamento para os sub-elementos associados ao refinamento atual*/
    TPZVec<TPZRefPattern *> fSubElRefPattern;

    /**refinamento para elementos de contorno associados ao refinamento atual*/
    TPZVec<TPZRefPattern *> fBCRefPattern;

 public:
    /**refinamento do contorno do element*/
    void SetBCRefPattern(TPZRefPattern *patt,int side){fBCRefPattern[side] = patt;}
    void ResizeBCPattern(int size){fBCRefPattern.Resize(size);}
    TPZRefPattern *BCRefPattern(int side){return fBCRefPattern[side];}
    /**refinamento dos sub-elementos*/
    void SetSubElRefPattern(TPZRefPattern *patt,int sub){fSubElRefPattern[sub] = patt;}
    void ResizeSubElPattern(int size){fSubElRefPattern.Resize(size);}
    TPZRefPattern *SubElRefPattern(int sub){return fSubElRefPattern[sub];}
    int NSubElPatt(){return fSubElRefPattern.NElements();}
  /** This method is used to create / identify the nodes of the refined elements.
   * The method verify if the nodes are already created by the self element or by some neighbour.
@param gel - pointer to the element which are being divided
@param newnodeindexes - return all midside node indexes for the element division. */
  void CreateNewNodes (TPZGeoEl *gel, TPZVec<int> &newnodeindexes);
  /** This method is used to create / identify the midside nodes for element elindex in its division process. The method verify if the nodes are already created by the self element or by some neighbour.
@param gel - pointer to the element which are being divided
@param side - Side along which the nodes will be identified/created
@param newnodeindexes - return all midside node indexes for the element division. */
  void CreateMidSideNodes (TPZGeoEl *gel, int side, TPZVec<int> &newnodeindexes);
};

#endif
