//
//  TPZMeshSolution.hpp
//  PZ
//
//  Created by Philippe Devloo on 30/10/16.
//
//

#ifndef TPZMeshSolution_hpp
#define TPZMeshSolution_hpp

#include <stdio.h>

#include "pzfunction.h"
#include "pzinterpolationspace.h"
#include "pzmultiphysicselement.h"


/// class which represents the solution and its derivative as a function
class TPZMeshSolution : public TPZFunction<STATE>
{

    /// Mesh for which the solution applies
    TPZCompMesh *fMesh;

    /// Geometric element index where the last point was found
    int64_t fGeoElIndex;
    
    /// Parametric coordinate where the last point was found
    TPZManVector<REAL,3> fLastLoc;
    
    /// dimension of the problem
    int fDimension;
    
    /// variable index of the solution variable
    int fSolutionVarindex;
    
    /// variable index of the gradient of the solution
    int fGradVarindex;
    
    /// number of state variables of the differential equation
    int fNumSolutions;
    
    /// the polynomial order of the function
    int fPolynomialOrder;
    
    /// index of the material object
    int fMaterialIndex;
    
public:
    
    /** @brief Class constructor */
    TPZMeshSolution(TPZCompMesh *cmesh, int materialid);
    
    /** @brief Class destructor */
    ~TPZMeshSolution();
    
    /**
     * @brief Performs function computation
     * @param x point coordinate which is suppose to be in real coordinate system but can be in master coordinate system in derived classes.
     * @param f function values
     * @param df function derivatives
     */
    virtual void Execute(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE> &df) override;
    
    /** @brief Returns number of functions. */
    virtual int NFunctions() const override
    {
        return fNumSolutions;
    }
    
    /** @brief Polynomial order of this function. */
    /** In case of non-polynomial function it can be a reasonable approximation order. */
    virtual int PolynomialOrder() const override
    {
        return fPolynomialOrder;
    }
    
    /** @brief Print a brief statement */
    virtual void Print(std::ostream &out) override;
    public:
int ClassId() const override;

    
};



#endif /* TPZMeshSolution_hpp */
