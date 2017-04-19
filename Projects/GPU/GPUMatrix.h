/**
 * @file
 * @brief Contains GPUMatrix which implements  matrix (using column major representation) to be used in GPU devices.
 */

#include "pzmatrix.h"

#ifndef _TGPUMATRIX_
#define _TGPUMATRIX_

/**
 * @brief GPU matrix class. \ref matrix "Matrix"
 * @note The GPU matrix class is special in that the data is stored column wise and operations are done in GPU devices.
 * @author Thiago Dias dos Santos
 * @since 04/2017
 */

template<class TVar=REAL>
class TGPUMatrix: public TPZMatrix<TVar> {
    
protected:
    /**
     @brief Element stiffness matrix. All elements in the mesh have the same stiffness matrix.
     */
    TVar *fKe;
    
public:
    
    /**
     @brief Constructor with initialization parameters
     @param rows Initial number of rows
     @param columns Number of columns
     @param buf Preallocated memory area which can be used by the matrix object
     @param size Size of the area pointed to by buf
     */
    inline TGPUMatrix(const long rows ,const long columns, TVar* buf) : TPZMatrix<TVar>(rows,columns), fKe(0) {
        if(rows*columns && buf) fKe=new TVar[rows*columns];
        else DebugStop();
        fKe=buf;
    }
    /**
     @brief Constructor with initialization parameters
     @param rows Initial number of rows
     @param columns Number of columns
     @param val Inital value fill all elements
     */
    inline TGPUMatrix(const long rows ,const long columns,const TVar & val ) : TPZMatrix<TVar>(rows,columns), fKe(0){
        if(rows*columns) fKe=new TVar[rows*columns];
        else DebugStop();
        for(long i=0;i<rows*columns;i++) fKe[i]=val;
    }
    /**
     @brief Constructor with initialization parameters, the default value is 0.
     @param rows Initial number of rows
     @param columns Number of columns
     */
    inline  TGPUMatrix(const long rows ,const long columns = 1) : TPZMatrix<TVar>(rows,columns), fKe(0) {
        if(rows*columns) fKe=new TVar[rows*columns];
        else DebugStop();
    }
    /**
     * @brief Copy constructor
     * @param refmat Used as a model for current object
     */
    inline TGPUMatrix(const TGPUMatrix<TVar> & refmat) : TPZMatrix<TVar>(refmat.fRow,refmat.fCol), fKe(0){
        if(refmat.fRow*refmat.fCol) fKe=new TVar[refmat.fRow*refmat.fCol];
        else DebugStop();
        memcpy(fKe,refmat,(size_t)(refmat.fRow*refmat.fCol)*sizeof(TVar));
    }
    /** @brief Simple destructor */
    inline virtual  ~TGPUMatrix(){
        if(fKe) delete[] fKe;
    }
    
};

#endif