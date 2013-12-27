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
#include "TPZJIntegral.h"

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
     * @param nstripes [in] : Amount of stripes in Y direction for applied pressure for reduced space
     * @param Qinj_well [in] : Well injection flow rate
     * @param visc [in] : Injected fluid viscosity
     * @param Jradius [in] : J-Integral radius
     * @param porder [in] : polinomial order of simulation
     *
     * TVD: True Vertical Depth (positive positions)
	 */
    TPZPlaneFractureKernel(TPZVec<TPZLayerProperties> & layerVec, REAL bulletTVDIni, REAL bulletTVDFin,
                           REAL xLength, REAL yLength, REAL Lmax, int nstripes, REAL Qinj_well, REAL visc,
                           REAL Jradius,
                           int porder);
    
    ~TPZPlaneFractureKernel();
    
    void Run();
    
protected:
    
    void FillFractureDotsCircle(TPZVec< std::pair<REAL,REAL> > &fractureDots);
    
    /**
     * @brief Method that will run a FEM simmulation of a classical vertical plane fracture
     * @param poligonalChain [in] : Poligonal chain that represents the crack boundary
     * @param initialElasticKickIsNeeded [in] : When is the first fracture geometry that is running, a initial solution is needed
     * @param justNewFractGeometry [in] : Flag for transferring elastic solution (no leakoff and no propagation criterion)
     */
    void RunThisFractureGeometry(TPZVec<std::pair<REAL,REAL> > &poligonalChain,
                                 bool initialElasticKickIsNeeded,
                                 bool justNewFractGeometry,
                                 int &step);
    /** Restore Qinj to the simulation default (Qinj gived by user for this simulation) */
    void RestoreQinj1wing_Hbullet();
    
    /**
     * @brief Method that will initializate the JPath3D vector structure, based on 1D cracktip elements (available on fPlaneFractureMesh atribute)
     */
    void InitializePath3DVector();
    
    /**
     * @brief Method that will run a FEM simmulation of linear elasticity in many stripes of newman condition
     * @param cmesh [in] : cmesh that the FEM will run
     */
    void ProcessElasticCMeshByStripes(TPZCompMesh * cmesh);
    
    /**
     * @brief Method that will run a FEM simmulation of linear elasticity for the given cmesh
     * @param an [in] : Given TPZAnalysis, initializated already
     * @param cmesh [in] : cmesh that the FEM will run
     */
    void SolveElasticity(TPZAnalysis &an, TPZCompMesh *cmesh);
    
    /**
     * @brief Method that will compute the stiff matrix for actual time step
     * @param an [in] : Given TPZAnalysis, initializated already
     * @param matK1 [out] : stiff matrix
     * @param fvec [out] : load vector
     */
    void StiffMatrixLoadVec(TPZAnalysis *an, TPZAutoPointer< TPZMatrix<REAL> > & matK1, TPZFMatrix<REAL> &fvec);
    
    /**
     * @brief Method that will compute the mass matrix for last time step
     * @param Un [out] : Mass matrix
     */
    void MassMatrix(TPZFMatrix<REAL> & Un);
    
    /** During development, this is used to check the convergence order of the non linear system */
    void CheckConv();
    
    /** Compute the volume of the interior of the fracture (by w=2*uy integral) */
    void PostProcessAcumVolW();

    /** Compute the volume of the leakoof */
    void PostProcessVolLeakoff(int step);
    
    /** Generate vtk for displacement post process */
    void PostProcessElasticity(int step);

    /** Generate vtk for pressure post process */
    void PostProcessPressure(int step);
    
    /** Auxiliar method for the PostProcessAcumVolW() method*/
    REAL IntegrateW(TPZCompMesh * elasticCMesh);
    
    /** Will check if fracture propagate and, if true, return new geometry of poligonal chain */
    bool CheckPropagationCriteria(TPZVec<std::pair<REAL,REAL> > &newPoligonalChain);
    
    /** 
     * Auxiliar method for CheckPropagationCriteria(...) method that computes the new geometry of the poligonal chain
     * @param whoPropagate_KI [in] : map that holds the KI and KIc, indexed by poligonal chain index
     */
    void DefinePropagatedPoligonalChain(std::map< int, std::pair<REAL,REAL> > &whoPropagate_KI,
                                        TPZVec< std::pair< REAL,REAL > > &newPoligonalChain);
    
    /**
     *  After the fracture propagation, this method will transfer the leakoff from the old data structure (based on given cmesh)
     *  to the new data structure (based on the new cmesh, keeped in fmeshVec atribute in position 0)
     */
    void TransferLeakoff(TPZCompMesh * cmeshFrom);
    
    //Atributes:
    TPZPlaneFractureMesh * fPlaneFractureMesh;
    
    TPZVec<TPZCompMesh *> fmeshVec;
    
    TPZCompMesh * fmphysics;
    
    REAL fHbullet;
    REAL fQinj1wing_Hbullet;
    REAL fCenterTVD;
    REAL fPoligonalChainInitialRadius;
    REAL fJIntegralRadius;
    REAL fvisc;
    int fpOrder;
    
    JIntegral3D fPath3D;
};


#endif /* defined(__PZ__TPZPlaneFractureKernel__) */
