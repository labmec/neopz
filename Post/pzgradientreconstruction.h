//
//  pzgradientreconstruction.h
//  PZ
//
//  Created by Agnaldo Farias on 4/10/13.
//
//

#ifndef __PZ__pzgradientreconstruction__
#define __PZ__pzgradientreconstruction__

#include <iostream>

#include "pzcompel.h"
#include "pzvec.h"
#include "pzcmesh.h"

class TPZGradientReconstruction
{
    
protected:
    
    /** @param fvar: Index of the Solution Variable*/
    int fvar;
    
    /** @param fDistortedMesh: Parameter indicates if the mesh is distorted*/
    bool fDistortedMesh;
    
    /** @param fSlopeLimiter: Parameter indicates if is  used the slope limiter*/
    bool fSlopeLimiter;
    
public:
    TPZGradientReconstruction(int var);
    
    ~TPZGradientReconstruction();
    
    TPZGradientReconstruction(const TPZGradientReconstruction &cp);
    
    TPZGradientReconstruction &operator=(const TPZGradientReconstruction &copy);
    
    /*
     *@brief Method to reconstruction gradient by Least Squares
     *@param cel [in]:  Computational element
     *@param center [out]:  Center point of the element
     *@param Grad [out]: Value of the gradient reconstructed
     */
    void GradientReconstructionByLeastSquares(TPZCompEl *cel,TPZManVector<REAL,3> &center,TPZVec<REAL> &Grad);
    
    /*
     *@brief Method to return data gradients
     *@param cmesh [in]: Computational emesh
     *@param datagradients [out]: Matrix of data about gradient in each element: Center x0 of the element; gradient value in x0; index and volume of the element
     */
    void GetDataGradient(TPZCompMesh *cmesh,TPZFMatrix<REAL> &datagradients);
    
    /*
     *@brief Method to calculate the cell mean of solution
     *@param cel [in]: computational element
     *@param IntOrder [in]:Integration order 
     */

    REAL MeanCell(TPZCompEl *cel,int IntOrder);
    
    /**@brief Method indicates the choice of using weights in the distorted mesh*/
    void UseWeightCoefficients();
    
     /**@brief Method indicates the choice of using slope limiter*/
    void UseSlopeLimiter();
    
    /**@brief Method to choose the node of the cell closer of your center point
     *@param cel [in]: computational element of the cell
     *@param nodecelX [out]: node closer of the center point in coord X
     */
    void NodeCloserCenterX(TPZCompEl *cel, TPZManVector<REAL> &nodecelX);
    
    /**@brief Method to calculate the weights that we will use in distorted meshes
     *@param cel [in]: computational element of the cell
     *@param neighscel [in]: neighbors set of the  cell
     *@param paramk [in]: parameter used:  between 1 and 2
     *@param weights [out]: vector with weights to each neighbor
     */
    void CalcWeights(TPZCompEl *cel, std::set<TPZCompEl *> neighscel, int paramk, TPZManVector<REAL> &weights);
    
    /**@brief Method to insert the weights in the matrices of the system by least squares
     *@param weights [out]: vector with weights 
     *@param DeltaH [out]: matrix of distances between the center points 
     *@param DifSol [out]: matrix of differences  between the solutions 
    */
    void InsertWeights(TPZVec<REAL> weights, TPZFMatrix<REAL> &DeltaH, TPZFMatrix<REAL> &DifSol);
    
    /**@brief Method to calculate the slope limiter
     *@param cel [in]: computational element of the cellt
     *@param solcel [in]: solution of the cell
     *@param neighscell [in]: vector of the neighbors of the  cell
     *@param solneighscell [in] vector of the sol of the neighbors of the  cell
     *@param alpha [out]: parameter of the slope limiter
     */
    void CalcSlopeLimiter(TPZCompEl *cel, REAL solcel, TPZStack<TPZCompEl *> neighscell, TPZStack<REAL> solneighscell, REAL &alpha);
};

#endif /* defined(__PZ__pzgradientreconstruction__) */
