// -*- c++ -*-
//$Id: pzcmesh.h,v 1.19 2004-05-25 12:58:48 erick Exp $
//HEADER FILE FOR CLASS MESH

#ifndef PZCMESHHPP
#define PZCMESHHPP


#include "pzadmchunk.h"
#include "pzfmatrix.h"
#include "pzblock.h"
#include "pzconnect.h"
#include "pzanalysis.h"
//#include <iostream>
//using namespace std;

#include <string.h>
#include "pzreal.h"	// Added by ClassView
#include "pzsave.h"

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
	/**
	 * Class TPZCompMesh implemments computational mesh
	 * @ingroup CompMesh
	 */
//-----------------------------------------------------------------------------
// class structure
//-----------------------------------------------------------------------------
class TPZCompMesh : public virtual TPZSaveable {

protected:
  /**
   * Geometric grid to which this grid refers
   */
  TPZGeoMesh	*fReference;
  
  /**
   * Grid name for model identification
   */
  string fName;
  

  /**
   * List of pointers to elements
   */
  TPZAdmChunkVector<TPZCompEl *>		fElementVec;

  /**
   * List of pointers to nodes
   */
  TPZAdmChunkVector<TPZConnect>			fConnectVec;

  /**
   * List of pointers to materials
   */
  TPZAdmChunkVector<TPZMaterial *>		fMaterialVec;
  
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
   * @param grefpatch - stack where will be inserted the geometric elements
   */
   void GetRefPatches(TPZStack<TPZGeoEl *> &grpatch);

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
  void SetName(const string &nm);

  /**set de dimension of the domain of the problem*/
  void SetDimModel(int dim){fDimModel = dim;}

  /**return the dimension of the simulation*/
  int Dimension() const {return fDimModel;}

  /**
   * Return the mesh name
   */
  string &Name() {return fName;}

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
  int NMaterials() {return fMaterialVec.NElements();}

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
  TPZAdmChunkVector<TPZMaterial *>	&MaterialVec() { return fMaterialVec; }

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
  int InsertMaterialObject(TPZMaterial *mat);

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
  void Print(ostream & out = cout);

  /**
   * Print the solution by connect index
   * @param out indicates the device where the data will be printed
   */
   void ConnectSolution(ostream & out);
  //@}

  /**
   * Find the material with identity id
   * @param id material id to be found
   */
  TPZMaterial* FindMaterial(int id);

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
   * Assemble the global stiffness matrix and right hand side
   * @param stiffness matrix where will be assembled the stiffness matrix
   * @param rhs vector where the load vector will be assembled
   */
  void Assemble(TPZMatrix &stiffness,TPZFMatrix &rhs);

  /**
   * Assemble only the right hand side
   * @param rhs vector where the load vector will be assembled
   */
  void Assemble(TPZFMatrix &rhs);

  /**
   * Builds the transfer matrix from the current grid to the coarse grid
   * @param coarsegrid grid for where the matrix will be transfered
   * @param transfer transfer matrix between the current mesh and the coarse mesh
   */
  void BuildTransferMatrix(TPZCompMesh &coarsemesh, TPZTransfer &transfer);

  /**
   * Para elementos decontínuos
   */
  void BuildTransferMatrixDesc(TPZCompMesh &transfermesh,TPZTransfer &transfer);
  void ProjectSolution(TPZFMatrix &projectsol);

  /**
   * Creates the computational elements, and the degree of freedom nodes
   */
  virtual void AutoBuild();

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
   * @param elements vector of elements whose father will be created????
   * @param index index for the father element????
   */
  void Coarsen(TPZVec<int> &elements, int &index);

  void CoarsenDisc(TPZVec<int> &elements, int &index);

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
  //@}


  /**
   * Clone this mesh
   */
  TPZCompMesh * Clone() const;

  /**
   * Copies the materials of this mesh to the given mesh
   */
  void CopyMaterials(TPZCompMesh *mesh) const ;

  REAL DeltaX();

  REAL MaximumRadiusOfMesh();

  REAL LesserEdgeOfMesh();

  /** cria uma malha obtida por aglomera¢ão de elementos,
   * accumlist relaciona a lista de elementos da malha fina que 
   * serão acumulados
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

#endif

