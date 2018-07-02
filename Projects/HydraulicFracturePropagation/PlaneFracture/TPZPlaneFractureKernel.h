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
#include "pzanalysis.h"


enum EWhoBlock {ENoBlock, EAllBlock, EBlockEntireFracure, EBlockStripes};

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
     * @param just1Stripe [in] : Requirement of just 1 stripe in Y direction for applied pressure for reduced space
     * @param Qinj_well [in] : Well injection flow rate
     * @param visc [in] : Injected fluid viscosity
     * @param Jradius [in] : J-Integral radius
     * @param porder [in] : polinomial order of simulation
     * @param MaxDispl [in] : Maximum displacement when fracture propagate
     * @param pressureIndependent [in] : Flag that mean if leakoff is pressure independent
     * @param just1stripe [in] : Flag that mean if reduced space will be a unique stripe
     *
     * TVD: True Vertical Depth (positive positions)
     */
    TPZPlaneFractureKernel(TPZVec<LayerProperties> & layerVec, REAL bulletTVDIni, REAL bulletTVDFin,
                           REAL xLength, REAL yLength, REAL Lmax, REAL Qinj_well, REAL visc,
                           REAL Jradius,
                           REAL MaxDispl,
                           bool pressureIndependent,
                           bool just1stripe,
                           bool layersStripesToo);
    
    ~TPZPlaneFractureKernel();
    
    void Run();
    
protected:
    
    void InitializePoligonalChain();
    
    void InitializeMeshes(TPZCompMesh * lastElastCMesh);
    
    /**
     * @brief Method that will run a FEM simmulation of linear elasticity in NULL newman condition (for reduced space)
     * @param cmesh [in] : cmesh that the FEM will run
     */
    void ProcessLinearElasticCMesh(TPZCompMesh * cmesh);
    
    /**
     * @brief Method that will run a FEM simmulation of a classical vertical plane fracture
     * @param volAcum [in] : Poligonal chain that represents the crack boundary
     * @param justTransferingElasticSolution [in] : time step
     */
    void RunThisFractureGeometry();
    
    void CloseActualTimeStep();
    
    /**
     * @brief Method that will initializate the JPath3D vector structure, based on 1D cracktip elements (available on fPlaneFractureMesh atribute)
     */
    void InitializePath3DVector();
    
    /**
     * @brief Method that will compute the stiff matrix for actual time step
     * @param an [in] : Given TPZAnalysis, initializated already
     * @param matK [out] : stiff matrix
     * @param matRes [out] : load vector (in matrix form)
     */
    void AssembleStiffMatrixLoadVec(TPZAnalysis * an,
                                    TPZAutoPointer< TPZMatrix<STATE> > & matK,
                                    TPZFMatrix<STATE> & matRes,
                                    EWhoBlock whoBlock);
    
    void ApplyInitialCondition(STATE val);
    
    void PutConstantPressureOnFluidSolution(STATE val);
    
    void ApplyEquationFilter(TPZAnalysis * an, EWhoBlock whoBlock);
    
    /**
     * @brief Method that will compute the mass matrix for last time step
     * @param massMat [out] : Mass matrix
     */
    void MassMatrix(TPZFMatrix<STATE> & massMat);
    
    /** During development, this is used to check the convergence order of the non linear system */
    void CheckConv();
    
    void UpdateLeakoffKernel();
    
    void PostProcessAll();
    
    void PostProcessSolutions(int num = -1);
    
    /** Compute the volume of the interior of the fracture (by w=2*uy integral) */
    void PostProcessAcumVolW();
    
    /** Compute the volume of the leakoff */
    void PostProcessLeakoff(int num = -1);
    
    /** Generate vtk for displacement post process */
    void PostProcessElasticity(int num = -1);
    
    /** Generate vtk for stress tensor invariants post process */
    void PostProcessInvariants(int num = -1);
    
    /** Generate vtk for pressure post process */
    void PostProcessPressure(int num = -1);
    
    /** Insert on globFractOutput3DData the actual Lfrac, Hsup and Hinf */
    void PostProcessOutputData(int num = -1);
    
    bool ThereIsNegativeNetPressure();
    
    /** Auxiliar method for the PostProcessAcumVolW() method*/
    // = 2. * Integral(uy_insideFracture)
    REAL IntegrateW(bool & thereIsNegW, REAL & negVol);
    
    /** 1 face from 1 wing fracture area */
    //Be aware that this area is of all fracture, not just permeable portion!!!
    REAL Fracture1wing_Area();
    
    // VlAcumLeakoff = Summ( elementArea * (2. * elementLeakoffPenetration) )
    REAL ComputeVlAcumLeakoff(TPZCompMesh * fluidCMesh);
    
    //Qinj * actTime
    REAL ComputeVolInjected();
    
    /** Will check if fracture propagate and, if true, return new geometry of poligonal chain */
    bool CheckPropagationCriteria(REAL & maxKI_KIc, std::set<int> & whoPropagate,
                                  int & maxKI_KIcPosition);
    
    /**
     * Auxiliar method for CheckPropagationCriteria(...) method that computes the new geometry of the poligonal chain
     * @param whoPropagate_KI [in] : map that holds propagated crackTip indexes
     */
    void DefinePropagatedPoligonalChain(REAL & maxKI_KIc, std::set<int> & whoPropagate);
    
    /**
     * Remove zig-zag from given poligonal chain.
     * Note: Zig-zag is when (v1.v2 < 0).
     * If zig-zag was not removed, we will have trouble on PerfercMatchRefPattern on fracture geomesh generation.
     */
    bool RemoveZigZag(TPZVec< std::pair<REAL,REAL> > &newPoligonalChain);
    
    void RemoveLayerInvasion(TPZVec< std::pair<REAL,REAL> > &newPoligonalChain);
    
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
    
    TPZCompMesh * fmphysics;
    
    REAL fLmax;
    REAL fHbullet;
    REAL fQinj1wing_Hbullet;
    REAL fCenterTVD;
    REAL fJIntegralRadius;
    REAL fvisc;
    
    REAL fResTop, fResBottom, fResRight;//Limits of reservoir domain in z and x
    REAL fZbulletTop, fZbulletBottom;
    
    REAL fMaxDispl;
    
    bool fjust1stripe;
    bool flayersStripesToo;
    
    JIntegral3D fPath3D;
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
    
    REAL falphaL;
};


#endif /* defined(__PZ__TPZPlaneFractureKernel__) */
