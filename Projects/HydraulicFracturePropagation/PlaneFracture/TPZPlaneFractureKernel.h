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
     * @param MaxDispl [in] : Maximum displacement when fracture propagate
     * @param MinDispl [in] : Minimum displacement when fracture propagate
     *
     * TVD: True Vertical Depth (positive positions)
	 */
    TPZPlaneFractureKernel(TPZVec<TPZLayerProperties> & layerVec, REAL bulletTVDIni, REAL bulletTVDFin,
                           REAL xLength, REAL yLength, REAL Lmax, int nstripes, REAL Qinj_well, REAL visc,
                           REAL Jradius,
                           int porder,
                           REAL MaxDispl,
                           REAL MinDispl);
    
    ~TPZPlaneFractureKernel();
    
    void Run();
    
protected:
    
    void InitializePoligonalChain();
    
    void InitializeMeshes();
    
    /**
     * @brief Method that will run a FEM simmulation of a classical vertical plane fracture
     * @param poligonalChain [in] : Poligonal chain that represents the crack boundary
     * @param step [in] : time step
     */
    void RunThisFractureGeometry(REAL & volAcum, TPZCompMesh * &lastPressureCMesh);
    
    void CloseActualTimeStep();
    
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
    
    void UpdateLeakoff();
    
    void PostProcessAll();
    
    void PostProcessSolutions();
    
    /** Compute the volume of the interior of the fracture (by w=2*uy integral) */
    void PostProcessAcumVolW();

    /** Compute the volume of the leakoof */
    void PostProcessVolLeakoff();
    
    /** Generate vtk for displacement post process */
    void PostProcessElasticity();

    /** Generate vtk for pressure post process */
    void PostProcessPressure();
    
    /** Insert on globFractOutput3DData the actual Lfrac, Hsup and Hinf */
    void PostProcessFractGeometry();
    
    /** Auxiliar method for the PostProcessAcumVolW() method*/
    REAL IntegrateW(TPZCompMesh * elasticCMesh);//<<<<<< precisa passar elasticCMesh??? Nao podia ser sempre fmeshVec[0]????
    
    REAL ComputeVlAcumLeakoff(TPZCompMesh * fluidCMesh);
    
    REAL ComputeVolInjected();
    
    /** Will check if fracture propagate and, if true, return new geometry of poligonal chain */
    bool CheckPropagationCriteria();
    
    /** 
     * Auxiliar method for CheckPropagationCriteria(...) method that computes the new geometry of the poligonal chain
     * @param whoPropagate_KI [in] : map that holds the KI and KIc, indexed by poligonal chain index
     */
    void DefinePropagatedPoligonalChain(REAL maxKI, REAL respectiveKIc,
                                        std::map< int, std::pair<REAL,REAL> > &whoPropagate_KI);
    
    /**
     * Remove zig-zag from given poligonal chain.
     * Note: Zig-zag is when (v1.v2 < 0).
     * If zig-zag was not removed, we will have trouble on PerfercMatchRefPattern on fracture geomesh generation.
     */
    bool RemoveZigZag(TPZVec< std::pair<REAL,REAL> > &newPoligonalChain);
    
    void TransferElasticSolution(REAL volAcum);
    
    /**
     *  After the fracture propagation, this method will transfer the leakoff from the old data structure (based on given cmesh)
     *  to the new data structure (based on the new cmesh, keeped in fmeshVec atribute in position 0)
     */
    void TransferLastLeakoff(TPZCompMesh * cmeshFrom);
    
    //Atributes:
    
    TPZVec< std::pair<REAL,REAL> > fpoligonalChain;
    int fstep;
    
    TPZPlaneFractureMesh * fPlaneFractureMesh;
    
    TPZVec<TPZCompMesh *> fmeshVec;
    TPZVec<TPZCompMesh *> fmeshVecElastic;
    
    TPZCompMesh * fmphysics;
    
    REAL fLmax;
    REAL fHbullet;
    REAL fQinj1wing_Hbullet;
    REAL fCenterTVD;
    REAL fPoligonalChainInitialRadius;
    REAL fJIntegralRadius;
    REAL fvisc;
    int fpOrder;
    
    REAL fMaxDispl;
    REAL fMinDispl;
    
    JIntegral3D fPath3D;
    
    int actColor = 0;
    std::string color[12] = {"Red","Green","Blue","Black","Gray","Cyan","Magenta","Yellow","Brown","Orange","Pink","Purple"};
};

class BezierCurve
{
public:
    BezierCurve();
    BezierCurve(TPZVec< std::pair< REAL,REAL > > &poligonalChain);
    ~BezierCurve();
    
    void F(REAL t, std::pair< REAL,REAL > & ft);
    
protected:
    
    REAL Bernstein(REAL t, REAL i);
    REAL Coef(int num1, int num2, int num3);
    
    int forder;
    TPZVec< std::pair< REAL,REAL > > fPoligonalChain;
};


#endif /* defined(__PZ__TPZPlaneFractureKernel__) */
