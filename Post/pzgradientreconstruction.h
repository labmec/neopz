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
    
    /*
     *@brief Class to calculate the gradient and slope limiter
     */
    class TPZGradientData
    {
    public:
        /*
         *@param cel: computational element to reconstruct the gradient
         *@param useweight [in]: if true indicates the choice of using weights in the distorted mesh
         *@param paramk[in]: value equal 1 or 2
         */
        TPZGradientData();
        
         ~TPZGradientData();
        
        TPZGradientData(const TPZGradientData &cp);
        
        TPZGradientData &operator=(const TPZGradientData &copy);
        
        void SetCel(TPZCompEl *cel, bool useweight, REAL paramK);
        
        
//        struct TNeighcell
//        {
//            int fIsBoundary;
//            STATE fNeighSol;
//            TPZManVector<REAL,3> fCenterNeighbour;
//            TPZManVector<REAL,3> fCenterInterface;
//        };
//        
//        TPZStack<TNeighcell> fData;
        
        void Print(std::ostream &out) const;
        
        void GetCenterPointAndCellAveraged(TPZCompEl *cel, TPZManVector<REAL,3> &xcenter, STATE &solcel);
        
        void InitializeGradData(TPZCompEl *cel);
        
        /*
         *@brief Method to reconstruction gradient by Least Squares and apply the slope limiter
         */
        void ComputeGradient();
        
        void QRFactorization(TPZFMatrix<REAL> &matA,TPZFMatrix<REAL> &vecb);
        
        /**
         *@brief Method to calculate the slope limiter (alphaK)
         */
        void ComputeSlopeLimiter();
        
        void ComputeSlopeLimiter2();
        
        void ComputeSlopeLimiter3();
        
        /**@brief Method to calculate the weights that we will use in distorted meshes
         *@param paramk [in]: parameter used in the harmonic mean, equal 1 or 2
         */
        void ComputeWeights(int paramk);
        
        /**@brief Method to choose the node of the cell closer of your center point
         *@param cel [in]: computational element of the cell
         *@param nodecelX [out]: node closer of the center point in coord X
         */
        void NodeCloserCenterX(TPZManVector<REAL,3> &nodecelX);
        
        /**@brief Method to insert the weights in the matrices of the system by least squares
         *@param DeltaH [out]: matrix of distances between the center points
         *@param DifSol [out]: matrix of differences  between the solutions
         */
        void InsertWeights(TPZFMatrix<REAL> &DeltaH, TPZFMatrix<REAL> &DifSol);
        
        
        /*
         *@brief Returns data from the computational element K (K cell)
         */
        void GetData(TPZManVector<REAL,3> &centerPoint, TPZManVector<STATE,3> &grad, STATE &cellAverage, STATE &slopeLimiter)
        {
            centerPoint = fCenterPointCellAndNeighbors[0];
            grad = fGradient;
            cellAverage = fSolCellAndNeighbors[0];
            slopeLimiter = fSlopeLimiter;
        }
      
        
    protected:
        
        /** @param SolCellAndNeighbors: Value of the cells averaged*/
        TPZStack<STATE> fSolCellAndNeighbors;
        
        /** @param CenterPointCellAndNeighbors: Center point of the element and neighbors*/
        TPZStack<TPZManVector<REAL,3> > fCenterPointCellAndNeighbors;
        
        /** @param CenterPointInterface: Center point of the interface between element and neighbors*/
        TPZStack<TPZManVector<REAL,3> > fCenterPointInterface;
        
        /** @param Gradient: gradient reconstructed*/
        TPZManVector<STATE,3> fGradient;
        
        /** @param SlopeLimiter: Value of the slope limiter*/
        STATE fSlopeLimiter;
        
        /** @param fdim: Dimension of the  element*/
        int fdim;
        
        /** @param fNeighborsCel: Neighbors of the cell*/
        TPZStack<TPZCompEl *> fCelAndNeighbors;
        
         /** @param fWeightsGrad: vector with weights to each neighbor*/
        TPZManVector<REAL,3> fWeightsGrad;
        
        /** @param fUseWeight: Parameter indicates if the mesh is distorted*/
        bool fUseWeight;
        
        /** @param fparamK: parameter used in the harmonic mean, equal 1 or 2*/
        REAL fparamK;
    };
    
public:
    
    /*
     *@param distmesh: Parameter indicates if the mesh is distorted
     *@param param: parameter used in the harmonic, equal 1 or 2, to calculate the weights
     */
    TPZGradientReconstruction(bool distmesh, REAL paramK);
    
    ~TPZGradientReconstruction();
    
    TPZGradientReconstruction(const TPZGradientReconstruction &cp);
    
    TPZGradientReconstruction &operator=(const TPZGradientReconstruction &copy);

    
    /*
     *@brief Method to replace the solution by finite element method from the reconstruction of the gradient.
     *Using the L2 projection
     *@param cmesh [in]: Computational mesh
     *@param datagradients [in]: Matrix of data about gradient in each element.
     *@param matidl2proj [in]: Id of the l2 projection material
     */
    void ProjectionL2GradientReconstructed(TPZCompMesh *cmesh, int matidl2proj);
    
    //change material id current to material id of the L2Projection
    void ChangeMaterialIdIntoCompElement(TPZCompEl *cel, int oldmatid, int newmatid);
    
    //assemble pos l2 projection
    void AssembleGlobalMatrix(TPZCompEl *el, TPZElementMatrix &ek, TPZElementMatrix &ef,TPZMatrix<STATE> & stiffmatrix, TPZFMatrix<STATE> &rhs);
    
    void GetInfoDistortedMesh(bool &useweight, REAL &paramK){
        useweight = fDistortedMesh;
        paramK = fparam;
    }
    
    TPZGradientData * fGradData;
    
protected:
    
    /** @param fDistortedMesh: Parameter indicates if the mesh is distorted*/
    bool fDistortedMesh;
    
    /** @param fparamK: parameter used in the harmonic mean, equal 1 or 2*/
    REAL fparam;
};

#endif /* defined(__PZ__pzgradientreconstruction__) */
