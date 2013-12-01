//
//  TPZPlaneFractureKernel.h
//  PZ
//
//  Created by Cesar Lucci on 18/11/13.
//
//

#ifndef __PZ__TPZPlaneFractureKernel__
#define __PZ__TPZPlaneFractureKernel__

#include <iostream>

#include "TPZPlaneFractureMesh.h"

class TPZPlaneFractureKernel
{
public:
    
    TPZPlaneFractureKernel();
    /**
	 * @brief Constructor
	 * @param layerVec [in] : vector of layers
	 * @param bulletTVDIni [in] : bullets perforation initial (TVD) depth
	 * @param bulletTVDFin [in] : bullets perforation final (TVD) depth
     * @param xLength [in] : Reservoir length in x direction (crack propagation direction)
     * @param yLength [in] : Reservoir length in y direction (tickness of reservoir that couple fracture plane)
     * @param Lmax    [in] : Maximum element edge length
     * @param nstripes [in] : Amounth of stripes in Y direction for applied pressure for reduced space
     *
     * TVD: True Vertical Depth (positive positions)
	 */
    TPZPlaneFractureKernel(TPZVec<TPZLayerProperties> & layerVec, REAL bulletTVDIni, REAL bulletTVDFin,
                           REAL xLength, REAL yLength, REAL Lmax, int nstripes, int porder);
    
    ~TPZPlaneFractureKernel();
    
    /**
     * @brief Method that will run a FEM simmulation of a classical vertical plane fracture
     * @param poligonalChain [in] : Poligonal chain that represents the crack boundary
     * @param sigmaTraction [in] : uniform traction on the farfield surface (farfield surface have normal {0,1,0})
     * @param vtkFile [in] : VTK file name for post processing
     * @param printVTKfile [in] : flag that will determine if vtkFile will be generated (true) or not (false)
     */
    bool RunThisFractureGeometry(const TPZVec<std::pair<REAL,REAL> > &poligonalChain,
                                 std::string vtkFile,
                                 bool printVTKfile = false);
    
    void ProcessElasticCMeshByStripes(TPZCompMesh * cmesh);
    
    void SolveInitialElasticity(TPZAnalysis &an, TPZCompMesh *cmesh);
    
    void StiffMatrixLoadVec(TPZAnalysis *an, TPZAutoPointer< TPZMatrix<REAL> > & matK1, TPZFMatrix<REAL> &fvec);
    
    void MassMatrix(TPZFMatrix<REAL> & Un);
    
    
protected:
    
    void CheckConv();
    
    TPZPlaneFractureMesh * fPlaneFractureMesh;
    
    TPZVec<TPZCompMesh *> fmeshVec;
    
    TPZCompMesh * fmphysics;
    
    int fpOrder;
};

#endif /* defined(__PZ__TPZPlaneFractureKernel__) */
