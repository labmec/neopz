/**
 * doxsave.h
 */

/**
 * @defgroup save Classes and methods which support persistency
 *
 * Persistency within the PZ environment are implemented by deriving class from the TPZSave
 * class and implementing the Read and Write method
 *
 * The association of the Class Id and a unique function is implemented in the
 * TPZRestoreClass . The mere "instantiation" of the class will create a global object
 * which will create the association
 */
