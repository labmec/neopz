/* Generated by Together */

#include "pzstrmatrix.h"
#include "pzcmesh.h" 

#include "TPZFrontNonSym.h"
#include "TPZFrontSym.h"

#include "pzelmat.h"


#ifndef TPZFRONTSTRUCTMATRIX_H
#define TPZFRONTSTRUCTMATRIX_H

struct TPZElementMatrix;


class TPZMatrix;
class TPZFMatrix;
class TPZCompMesh;


/**
 *Class responsible for a interface among Finite Element Package and Matrices package \n
 *Prevents users from all the necessary information to work with all matrices classes \n
 *It facilitates considerably the use of TPZAnalysis
 * @ingroup structural frontal
 */
template<class front> //! Type parameter for TPZFrontStructMatrix frontal matrix. \n It can assume values TPZFrontSym and TPZFrontNonSym for symmetric and non symmetric matrices
class TPZFrontStructMatrix : public TPZStructMatrix {

protected:
     /**
      * This vector contains an ordered list. \n
      * The elements must be asssembled in that order so the frontal works on its best \n
      * performance
      */
     TPZVec<int> fElementOrder;
    int f_quiet;
     /**
      * Returns a vector containing all elements connected to a degree of freedom.
      */
     void GetNumElConnected(
          TPZVec <int> &numelconnected //! Vector containing the number of connections for every ith dof.
          );
     /**
      * It is applied over fElementOrder putting it in the correct order.
      */
     void OrderElement();//TPZVec <int> &elorder);
public:    

    /**
     * Class constructor
     * @url http://www.fec.unicamp.br/~longhin
     * @url http://www.fec.unicamp.br/~phil 
     */ 
     TPZFrontStructMatrix(TPZCompMesh *);

     static int main();

    /**
     * Class destructor
     * @url http://www.fec.unicamp.br/~longhin
     * @url http://www.fec.unicamp.br/~phil 
     */ 

     ~TPZFrontStructMatrix();
  
     
     /**
      *Returns a pointer to TPZMatrix
      */
     TPZMatrix * Create();
     
     
     /**
      *It clones a TPZFrontStructMatrix
      */
     TPZStructMatrix * Clone();


     /**
      *Assemble a stiffness matrix according to rhs
      */ 	
     void AssembleNew(
     	TPZMatrix & stiffness //! Stiffness matrix to assembled
	     , TPZFMatrix & rhs //! Matrix contaning ???
	);
     
     /**
      *Assemble a stiffness matrix.
      */ 	
     void Assemble(
     	TPZMatrix & stiffness //! Stiffness matrix to assembled
     	, TPZFMatrix & rhs //! Vector contaning loads
          );
     
     /**
      * Computes element matrices. \n
      * Each computed element matrices would then be added to Stiffness matrix
      */
     void AssembleElement(
          TPZCompEl *el //! Actual element being computed
          , TPZElementMatrix & ek //! Formed element matrix
          , TPZElementMatrix & ef //! Formed element load matrix
          , TPZMatrix & stiffness //! Global Stiffness matrix
          , TPZFMatrix & rhs //!Global load matrix
          ); 
     
     /**
      * Returns a pointer to TPZMatrix. \n
      * This is a mandatory function, it is neded by all StructMatrix. \n
      * Except in frontal matrices, the returned matrix is not in its decomposed form.
      */
     TPZMatrix * CreateAssemble(
          TPZFMatrix &rhs //!Load matrix
          );
    void SetQuiet(int quiet);

   
};
#endif //TPZFRONTSTRUCTMATRIX_H
