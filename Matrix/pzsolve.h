#ifndef TPREH
#define TPREH

#include "pzfmatrix.h"

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
  virtual ~TPZSolver() {}

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
  TPZMatrixSolver(TPZMatrix *Refmat = 0);

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
  void SetMatrix(TPZMatrix *Refmat);

  /**
     Resets current object
  */
  void ResetMatrix();

  /**
     Returns a pointer to TPZMatrix
  */
  TPZMatrix *Matrix() { return fContainer->Matrix();}

  /**
     Shares the current matrix with another object of same type
     @param other Object that will share current matrix
  */
  void ShareMatrix(TPZMatrixSolver & other);
  /**
     Produces some diagnoses of current object
     @param out Output device
  */
  static void Diagnose(ostream &out = cout) {
    out << "TPZMatrixSolver::Diagnose\nnumber of objects created " << gnumcreated << 
      "\nnumber of objects deleted " << gnumdeleted << endl;
    TPZContainer::Diagnose(out);
  }

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
  enum MSolver {ENoSolver, EJacobi, ESOR, ESSOR, ECG, EGMRES, EDirect, EMultiply};

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
    static void Diagnose(ostream &out = cout);
  };
  /**
     Container classes
  */
  TPZContainer *fContainer;
  static int gnumcreated;
  static int gnumdeleted;
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

