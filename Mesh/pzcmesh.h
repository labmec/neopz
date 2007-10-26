// -*- c++ -*-
//$Id: pzcmesh.h,v 1.37 2007-10-26 13:18:58 tiago Exp $
//HEADER FILE FOR CLASS MESH

#ifndef PZCMESHHPP
#define PZCMESHHPP


#include "pzadmchunk.h"
#include "pzfmatrix.h"
#include "pzblock.h"
#include "pzconnect.h"
#include "tpzautopointer.h"
#include "pzreal.h" // Added by ClassView
#include "pzsave.h"

//#include "pzanalysis.h"
//#include <iostream>
#include <map>
#include <iostream>
#include <set>
#include <string>

class TPZCompEl;
class TPZGeoEl;
struct TPZCompElBC;
class TPZConnect;
struct TPZConnectBC;
class TPZBndCond;
class TPZMaterial;
class TPZGeoMesh;
class TPZTransfer;
class TPZCoSys;
class TPZGeoEl;
class TPZStream;
class TPZCompElDisc;
class TPZInterpolatedElement;
template<class T> class TPZReferredCompEl;
template<class T> class TPZIntelGen;

/**
 * @brief Class TPZCompMesh implemments computational mesh
 The computational mesh is a repository for computational elements, nodes and
 material objects
 The computational mesh also contains the current solution of the mesh and an
 elementwise solution vector
 The data structure of this object is rather simple
 * @ingroup CompMesh
 */
class TPZCompMesh : public virtual TPZSaveable {

protected:
  /**
   * Geometric grid to which this grid refers
   */
  TPZGeoMesh	*fReference;

  /**
   * Grid name for model identification
   */
  std::string fName;


  /**
   * List of pointers to elements
   */
  TPZAdmChunkVector<TPZCompEl *>		fElementVec;

  /**
   * List of pointers to nodes
   */
  TPZAdmChunkVector<TPZConnect>			fConnectVec;

  /**
   * Map of pointers to materials
   */
  std::map<int, TPZAutoPointer<TPZMaterial> >	fMaterialVec;

  /**
   * List of nodes with associated boundary conditions
   */
  //  TPZAdmChunkVector<TPZConnectBC>		fBCConnectVec;

  /**
   * Block structure of the solution vector ????
   */
  TPZBlock		fSolutionBlock;

  /**Solution vector*/
  TPZFMatrix	fSolution;

  /**
   * Block structure to right construction of the
   * stiffness matrix and load vector
   */
  TPZBlock		fBlock;

  /**
   * Solution vectors organized by element
   */
  TPZFMatrix fElementSolution;

  /*set the dimension of the simulation or the model*/
  int fDimModel;

  /**
   * Default order for all elements of this mesh
   */
  int fDefaultOrder;

public:

  /**
   * Constructor from geometrical mesh
   * @param gr pointer to geometrical reference mesh
   */
  TPZCompMesh(TPZGeoMesh* gr=0);

  /**
   * Copy constructor
   */
  TPZCompMesh(const TPZCompMesh &copy);

  /**
   * Simple Destructor
   */
  virtual ~TPZCompMesh();

  /**
   * This method will initiate the comparison between the current computational
   * mesh and the mesh which is referenced by the geometric mesh
   * @param var state variable number
   * @param matname pointer to material name
   */
  virtual REAL CompareMesh(int var, char *matname);
  /**
   * Gives all patches of the mesh
   * @param grpatch - stack where will be inserted the geometric elements
   */
   void GetRefPatches(TPZStack<TPZGeoEl *> &grpatch);
  /**
   * Gives all patches of the mesh
   * @param grpatch - set with all reference geometric elements
   */
   void GetRefPatches(std::set<TPZGeoEl *>  &grpatch);

   /**
    * Gives the conects graphs
    * @param nodegraph
    * @param nodegraphindex
    */
   void GetNodeToElGraph(TPZVec<int> &nodtoelgraph, TPZVec<int> &nodtoelgraphinde,TPZStack<int> &elgraph, TPZVec<int> &elgraphindexx);

    /**
     * Gives the element patch
     * @param gel: geometric reference element for the patch
     * @param patch: patch of elements
     */
   void GetElementPatch(TPZVec<int> nodtoelgraph, TPZVec<int> nodtoelgraphindex, TPZStack<int> &elgraph, TPZVec<int> &elgraphindex,int elind ,TPZStack<int> &patch);

  /**
   * Set the mesh name
   */
  void SetName(const std::string &nm);

  /**set de dimension of the domain of the problem*/
  void SetDimModel(int dim){fDimModel = dim;}

  /** Sets the geometric reference mesh */
  void SetReference(TPZGeoMesh * gmesh);

  /**return the dimension of the simulation*/
  int Dimension() const {return fDimModel;}

  /**
   * Return the mesh name
   */
  std::string &Name() {return fName;}

  /**
   * Delete all the dynamically allocated data structures
   */
  void CleanUp();

  /**
   * @name Access_Private_Data
   * Routines for access to private data
   */
  //@{
  /**
   *Number of connects allocated including free nodes
   */
  int NConnects() {return fConnectVec.NElements();}

  /**
   * Number of independent connect objects
   */
  int NIndependentConnects();

  /**
   * Number of computational elements allocated
   */
  int NElements() {return fElementVec.NElements();}

  /**
   * Number of materials
   */
  int NMaterials() {return fMaterialVec.size();}

  /**
   * Number of connects with boundary condition
   */
  //  int NBCConnects() {return fBCConnectVec.NElements();}
  //@}

  /**
   * @name Access_Internal_Data_Structures
   * Methods for accessing the internal data structures
   */
  //@{
  /**
   * Return a reference to the element pointers vector
   */
  TPZAdmChunkVector<TPZCompEl *>   &ElementVec() { return fElementVec; }

  /**
   * Return a reference to the connect pointers vector
   */
  TPZAdmChunkVector<TPZConnect>    &ConnectVec() { return fConnectVec; }

  /**
   * Return a reference to the material pointers vector
   */
  std::map<int ,TPZAutoPointer<TPZMaterial> >	&MaterialVec() { return fMaterialVec; }

  /**
   * Return a reference to the vector of nodal boundary condition pointers
   */
  //  TPZAdmChunkVector<TPZConnectBC>  &BCConnectVec() { return fBCConnectVec; }

  /**
   * Return a pointer to the geometrical mesh associated
   */
  TPZGeoMesh *Reference() const { return fReference; }

  /**
   * Access the block structure of the solution vector
   */
  TPZBlock &Block() { return fBlock;}

  /**
   * Access the solution vector
   */
  TPZFMatrix &Solution(){ return fSolution;}

  /**
   * Access method for the element solution vectors
   */
  TPZFMatrix &ElementSolution() { return fElementSolution;}
  //@}

  /**
   * @name Modify_Mesh_Structure
   * Methods for modify mesh structure, materials, elements, node etc
   */
  //@{
  /**
   * Return an index to a new connect
   * @param blocksize new connect block size - default value = 0
   */
  virtual  int AllocateNewConnect(int blocksize=0);

  /**
   * Insert a material object in the datastructure
   @ @param mat pointer to the material
   */
  int InsertMaterialObject(TPZAutoPointer<TPZMaterial> & mat);

  /**
   * Resequence the block object, remove unconnected connect objects
   * and reset the dimension of the solution vector
   */
  void InitializeBlock();

  /**
   * Adapt the solution vector to new block dimensions
   */
  virtual void ExpandSolution();

  /**
   * Create degree of freedom boundary conditions
   */
  //  void CreateConnectBC();
  //@}

  /**
   * @name Access_Solution
   * Methods for access and manipulates solution
   */
  //@{
  /**
   * Set a ith element solution, expanding the element-solution matrix if necessary
   */
  void SetElementSolution(int i, TPZVec<REAL> &sol);
  //@}

  /**
   * @name Print
   * Print methods
   */
  //@{
  /**
   * Prints mesh data
   * @param out indicates the device where the data will be printed
   */
  virtual void Print(std::ostream & out = std::cout);

  /**
   * Print the solution by connect index
   * @param out indicates the device where the data will be printed
   */
   void ConnectSolution(std::ostream & out);
  //@}

  /**
   * Find the material with identity id
   * @param id material id to be found
   */
  TPZAutoPointer<TPZMaterial> FindMaterial(int id);

  /**
   * Map this grid in the geometric grid
   */
  void LoadReferences();

  /**
   * Compute the number of elements connected to each connect object
   */
  virtual void ComputeNodElCon();

  /**
   * Delete the nodes which have no elements connected to them
   */
  virtual void CleanUpUnconnectedNodes();


  /**
   * Computes the connectivity graph of the elements, as appropriate for the
   * TPZRenumbering class
   * @param elgraph stack of elements to create the grapho????
   * @param elgraphindex graphos indexes vector
   */
  void ComputeElGraph(TPZStack<int> &elgraph, TPZVec<int> &elgraphindex);

  /**
   * @name Submesh_methods
   * Methods required by submeshes
   */
  //@{
  /**
   * Get the father meshes stack
   */
  virtual TPZCompMesh *FatherMesh() {return NULL;}

  /**
   * Makes a specified connection a internal mesh connection.
   * @param local connection local number to be processed
   */
  virtual void MakeInternal(int local) {;}

  /**
   * Make all mesh connections internal mesh connections.
   * connects to an internal connection
   */
  virtual void MakeAllInternal(){;}

  /**
   * Returns the rootmesh who have the specified connection.
   */
  virtual TPZCompMesh* RootMesh (int local) {return this;}

  /**
   * Transfer one element from a specified mesh to the current submesh.
   * @param mesh pointer to the mesh whose the element from
   * @param elindex element index to transfer
   */
  virtual int TransferElementFrom(TPZCompMesh *mesh, int elindex) {return elindex;}

  /**
   * Transfer one element from a submesh to another mesh.
   * @param mesh mesh pointer who will receive the element
   * @param elindex element index to transfer
   */
  virtual int TransferElementTo(TPZCompMesh *mesh, int elindex) {return elindex;}

  /**
   * Transfer one element form a submesh to another mesh.
   * @param mesh mesh pointer who will receive the element
   * @param elindex element index to transfer
   */
  virtual int TransferElement(TPZCompMesh *mesh, int elindex) {return elindex;}

   /**
    * Put an local connection in the supermesh - Supermesh is one
    * mesh who contains the analised submesh.
    * @param local local index of the element to be trasnfered
    * @param super pointer to the destination mesh
    */
  virtual int PutinSuperMesh (int local, TPZCompMesh *super);

   /**
    * Get an external connection from the supermesh - Supermesh is one
    * mesh who contains the analised submesh.
    * @param superind index of the element to be trasnfered
    * @param super pointer to the destination mesh
    */
  virtual int GetFromSuperMesh (int superind, TPZCompMesh *super);
  //@}


  int GetDefaultOrder()
  {
    return fDefaultOrder;
  }

  void SetDefaultOrder( int order )
  {
    fDefaultOrder = order;
  }


  /**
   * @name SCIENTIFIC_ROUTINES
   * Scientific manipulating routines
   */

  //@{
  /**
   *This computes the number of equations associated with non-restrained nodes
   */
  int NEquations();

  /**
   *This method computes the bandwidth of the system of equations
   */
  int BandWidth();

  /**
   * This method computes the skyline of the system of equations
   * @param skyline vector where the skyline will be computed
   */
  void Skyline(TPZVec<int> &skyline);

  /**
   * Assemble the vector with errors estimators
   * @param estimator vector where will be assembled the errors
   * @param errorid index for dual or wheeler estimator
   */
  void AssembleError(TPZFMatrix &estimator, int errorid);


  /**
   * Builds the transfer matrix from the current grid to the coarse grid
   * @param coarsegrid grid for where the matrix will be transfered
   * @param transfer transfer matrix between the current mesh and the coarse mesh
   */
  void BuildTransferMatrix(TPZCompMesh &coarsemesh, TPZTransfer &transfer);

  /**
   * Para elementos decont�uos
   */
  void BuildTransferMatrixDesc(TPZCompMesh &transfermesh,TPZTransfer &transfer);
  void ProjectSolution(TPZFMatrix &projectsol);

  /**
   * Creates the computational elements, and the degree of freedom nodes
   */
  virtual void AutoBuild();

//   /**
//    * Enumerate to help AutoBuildContDisc() to be nice.
//    */
//   enum MCreationType{ ENone = 0, EContinuousEl = 1, EDiscontinuousEl = 2};

  /**
   * Creates the computational elements, and the degree of freedom nodes
   * Elements created may be TPZInterpolatedElement or TPZCompElDisc.
   * indices contains the type of the element. Element type are given by the enumerate MCreationType.
   */
  virtual void AutoBuildContDisc(const TPZVec<TPZGeoEl*> &continuous, const TPZVec<TPZGeoEl*> &discontinuous);

  void AutoBuild(std::set<int> &MaterialIDs);

static  void SetAllCreateFunctionsDiscontinuous();
static  void SetAllCreateFunctionsContinuous();
static  void SetAllCreateFunctionsDiscontinuousReferred();
static  void SetAllCreateFunctionsContinuousReferred();
static  void SetAllCreateFunctions(TPZCompEl &cel);

  /**
   * Will build the list of element boundary conditions build the list of connect boundary conditions.
   * Put material pointers into the elements. Check on the number of dof of the connects
   */
  int Consolidate();

  /**
   * ??
   */
  int Check();

  /**
   * ??
   */
  //  char IsChecked() { return fChecked; }

  /**
   * Given the solution of the global system of equations, computes and stores the
   * solution for the restricted nodes
   * @param sol given solution matrix
   */
  void LoadSolution(TPZFMatrix &sol);

  /**
   * Divide the element corresponding to index
   * @param index element to divide index
   * @param subindex divided elements indexes vector
   * @param interpolate ????
   */
  void Divide(int index,TPZVec<int> &subindex,int interpolate = 0);

  /**
   * Create a computational element father for the comp. elem. into elements
   * @param elements indices of subelements to be agrouped
   * @param index of new coarse element
   * @param CreateDiscontinuous = false indicates a TPZInterpolatedElement must be created. True indicates a TPZCompElDisc
   */
  void Coarsen(TPZVec<int> &elements, int &index, bool CreateDiscontinuous = false);

  /** Deletes all interfaces and rebuild them all */
  void RemakeAllInterfaceElements();

  /**
   * Will refine the elements associated with a boundary condition till there are
   * no elements constrained by boundary condition elements
   */
  void AdjustBoundaryElements();

  /**
   * Permute the sequence number of the connect objects
   * @param permute vector of elements to permute
   */
  void Permute(TPZVec<int> &permute);

  /**
   *Evaluates the error given the two vectors of the analised parameters
   * @param loc local vector of the analised parameter
   * @param val given vector to compare
   * @param deriv ????
   * @param true_error return the true error between the given vectors
   * @param L2_error return the L2 norm of the error between the given vectors
   * @param estimate ????
   */
  void EvaluateError(void (*fp)(TPZVec<REAL> &loc,TPZVec<REAL> &val,TPZFMatrix &deriv),
                                 TPZVec<REAL> &errorSum);

  /** This method compute the jump solution of interface and convert discontinuous elements with
   * jump less than eps in continuous elements.
   * It may be compared the following values to eps:
   * int val = 1: (leftsol - rightsol)^2
   * int val = 2: (DSolLeft - DSolRight)^2
   * int val = 0: (leftsol - rightsol)^2 + (DSolLeft - DSolRight)^2
   */
  void ConvertDiscontinuous2Continuous(REAL eps, int val, bool InterfaceBetweenContinuous = false);

  void Discontinuous2Continuous(int disc_index, int &new_index, bool InterfaceBetweenContinuous = false);

  //@}


  /**
   * Clone this mesh
   */
  TPZCompMesh * Clone() const;

  /**
   * Copies the materials of this mesh to the given mesh
   */
  void CopyMaterials(TPZCompMesh &mesh) const ;

  REAL DeltaX();

  REAL MaximumRadiusOfMesh();

  REAL LesserEdgeOfMesh();

  /** cria uma malha obtida por aglomera� de elementos,
   * accumlist relaciona a lista de elementos da malha fina que
   * ser� acumulados
   */
  TPZCompMesh *ComputeMesh(TPZVec<int> &accumlist,int numaggl);
  /** This method will fill the matrix passed as parameter with a representation of the fillin of the global stiffness matrix, based on the sequence number of the connects
@param resolution Number of rows and columns of the matrix
@param fillin Matrix which is mapped onto the global system of equations and represents the fillin be assigning a value between 0. and 1. in each element */
  void ComputeFillIn(int resolution, TPZFMatrix &fillin);

  /**
  * returns the unique identifier for reading/writing objects to streams
  */
  virtual int ClassId() const;
  /**
  Save the element data to a stream
  */
  virtual void Write(TPZStream &buf, int withclassid);

  /**
  Read the element data from a stream
  */
  virtual void Read(TPZStream &buf, void *context);


};


/**
 * Allocates new connecto in the mesh
 */
inline int TPZCompMesh::AllocateNewConnect(int blocksize) {
  int connectindex = fConnectVec.AllocateNewElement();
  int blocknum = fBlock.NBlocks();
  fBlock.SetNBlocks(blocknum+1);
  fBlock.Set(blocknum,blocksize);
  fConnectVec[connectindex].SetSequenceNumber(blocknum);
  return connectindex;
}

inline void TPZCompMesh::SetReference(TPZGeoMesh * gmesh){
  this->fReference = gmesh;
}

#endif

