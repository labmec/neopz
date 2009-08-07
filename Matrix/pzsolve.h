#ifndef TPREH
#define TPREH

#include "pzfmatrix.h"

class TPZMatrixSolver;

/**
   @ingroup solver
   Defines a abstract class of solvers  which will be used by matrix classes
*/
class TPZSolver {

 public:

  // virtual void Solve(TPZFMatrix &F, TPZFMatrix &result) = 0;
  /**
     Solves the system of linear equations stored in current matrix
     @param F contains Force vector
     @param result contains the solution
     @param residual contains the residual for that linear system
  */
  virtual void Solve(const TPZFMatrix &F, TPZFMatrix &result, TPZFMatrix *residual = 0) = 0;

  /**
     Clones the current object returning a pointer of type TPZSolver
  */
  virtual TPZSolver *Clone() const = 0;
  /**
     Destructor
  */
  virtual ~TPZSolver();
  
  /**
  This method will reset the matrix associated with the solver
  This is useful when the matrix needs to be recomputed in a non linear problem
  */
  virtual void ResetMatrix() {}
  /**
  This method gives a preconditioner to share a matrix with the referring solver object
  */
//  virtual void SetMatrix(TPZMatrixSolver *solver);
  /**
  * Updates the values of the current matrix based on the values of the matrix
  */
  virtual void UpdateFrom(TPZMatrix *matrix)
  {
    std::cout << __PRETTY_FUNCTION__ << " called\n";
  }

};


/**
   @ingroup solver
   Defines a class of matrix solvers
*/
class TPZMatrixSolver : public TPZSolver {

 public:

  /**
     Constructor with initialization parameter
     @param RefMat Sets reference matrix to 0
  */
  TPZMatrixSolver(TPZAutoPointer<TPZMatrix> Refmat);
  
  TPZMatrixSolver();

  /**
     Copy constructor
     @param copy Model object to be copied from
  */
  TPZMatrixSolver(const TPZMatrixSolver &Source);

  /**
     Destructor
  */
  virtual ~TPZMatrixSolver();

  /**
     Sets a matrix to the current object
     @param RefMat Sets reference matrix to RefMat     
  */
virtual  void SetMatrix(TPZAutoPointer<TPZMatrix> Refmat);

  /**
  * Updates the values of the current matrix based on the values of the matrix
  */
  virtual void UpdateFrom(TPZAutoPointer<TPZMatrix> matrix)
  {
    if(fReferenceMatrix == matrix && matrix)
    {
      if(this->fContainer) this->fContainer->UpdateFrom(matrix.operator->());
    }
  }
  /**
     Resets current object
  */
  void ResetMatrix();

  /**
  This method gives a preconditioner to share a matrix with the referring solver object
  */
  virtual void SetReferenceMatrix(TPZAutoPointer<TPZMatrix> matrix)
  {
    fReferenceMatrix = matrix;
  }

  /**
     Returns a pointer to TPZMatrix
  */
  TPZAutoPointer<TPZMatrix> Matrix() const { return fContainer;}

  /**
     Shares the current matrix with another object of same type
     @param other Object that will share current matrix
  */
  void ShareMatrix(TPZMatrixSolver & other);

 protected:
  /**
     @enum MSolver
     Defines a series of solvers available in PZ
     @param ENoSolver No solver selected
     @param EJacobi Jacobi solver selected
     @param ESOR Successive Over Relaxation solver selected
     @param ESSOR Symmetric Successive Over Relaxation solver selected
     @param ECG Conjugate Gradiente solver selected
     @param EDirect Jacobi solver selected
  */
  enum MSolver {ENoSolver, EJacobi, ESOR, ESSOR, ECG, EGMRES, EBICGSTAB, EDirect, EMultiply};

 private:
  /**
     Defines a class of containers
     @ingroup solver
  */
  class TPZContainer {
    static int gnumcreated;
    static int gnumdeleted;
    int fRefCount;
    TPZMatrix *fRefMat;
  public:
    TPZContainer(TPZMatrix *mat);
    ~TPZContainer();
    void IncreaseRefCount();
    void DecreaseRefCount();
    TPZMatrix *Matrix();
    void SetMatrix(TPZMatrix *mat);
    static void Diagnose(std::ostream &out = std::cout);
  };
  /**
     Container classes
  */
  TPZAutoPointer<TPZMatrix> fContainer;
  /**
  * Reference matrix used to update the current matrix
  */
  TPZAutoPointer<TPZMatrix> fReferenceMatrix;
  //	TPZSolver *fPrecond;

  //misael

 protected:
  /**
     Manipulation matrix
  */
  TPZFMatrix fScratch;
};


;

#endif  // TPREH

