//
//  toolstransienttime.h
//  PZ
//
//  Created by Agnaldo Farias on 9/5/12.
//  Copyright (c) 2012 LabMec-Unicamp. All rights reserved.
//

#ifndef PZ_toolstransienttime_h
#define PZ_toolstransienttime_h

#include "pzcmesh.h"
#include "pzanalysis.h"
#include "tpzcompmeshreferred.h"
#include "pznlfluidstructure2d.h"




struct InputDataStruct
{
public:
    
    InputDataStruct()
    {
        
    }
    ~InputDataStruct()
    {
        
    }
    
    void SetData(REAL Lx, REAL Ly, REAL Lf, REAL Hf, REAL E, REAL Poiss, REAL Fx, REAL Fy, REAL Visc, REAL SigN,
                 REAL QinjTot, REAL Ttot, REAL deltaT, REAL Cl, REAL Pe, REAL SigmaConf, REAL Pref, REAL vsp, REAL KIc)
    {
        fLx = Lx;
        fLy = Ly;
        fLf = Lf;
        fHf = Hf;
        
        fE = E;
        fPoiss = Poiss;
        fFx = Fx;
        fFy = Fy;
        
        fVisc = Visc;
        
        fSigN = SigN;
        
        REAL Qinj1asa = QinjTot / 2.;
        REAL QinjSecao = Qinj1asa / Hf;
        fQinj = QinjSecao;
        
        fTtot = Ttot;
        fdeltaT = deltaT;
        
        fCl = Cl;
        fPe = Pe;
        fSigmaConf = SigmaConf;
        fPref = Pref;
        fvsp = vsp;
        
        fKIc = KIc;
    }
    
    void SetLf(REAL Lf)
    {
        fLf = Lf;
    }
    
    REAL Lx() { return fLx; }
    REAL Ly() { return fLy; }
    REAL Lf() { return fLf; }
    REAL Hf() { return fHf; }
    REAL E() { return fE; }
    REAL Poiss() { return fPoiss; }
    REAL Fx() { return fFx; }
    REAL Fy() { return fFy; }
    REAL Visc() { return fVisc; }
    REAL SigN() { return fSigN; }
    REAL Qinj() { return fQinj; }
    REAL Ttot() { return fTtot; }
    REAL deltaT() { return fdeltaT; }
    REAL Cl() { return fCl; }
    REAL Pe() { return fPe; }
    REAL SigmaConf() { return fSigmaConf; }
    REAL Pref() { return fPref; }
    REAL vsp() { return fvsp; }
    REAL KIc() { return fKIc; }
    
private:
    
    //Dimensions:
    REAL fLx;//Dimensao em x do domínio da malha do MEF
    REAL fLy;//Dimensao em y do domínio da malha do MEF
    REAL fLf;//Comprimento de 1/2 asa da fratura
    REAL fHf;//Altura da fratura
    
    //Elastic properties:
    REAL fE;//Modulo de elasticidade
    REAL fPoiss;//Poisson
    REAL fFx;//Bodyforces in x
    REAL fFy;//Bodyforces in y
    
    //Fluid property:
    REAL fVisc;//viscosidade do fluido de injecao
    
    //BCs:
    REAL fSigN;//Sigma.n no problema elastico que servira de espaco de aproximacao para o elastico multifisico
    REAL fQinj;//vazao de 1 asa de fratura dividido pela altura da fratura
    
    //time:
    REAL fTtot;//Tempo total da simulacao
    REAL fdeltaT;//deltaT
    
    //Leakoff:
    REAL fCl;//Carter
    REAL fPe;//Pressao estatica
    REAL fSigmaConf;//Tensao de confinamento
    REAL fPref;//Pressao de referencia da medicao do Cl
    REAL fvsp;//spurt loss
    
    //Propagation criterion
    REAL fKIc;
};







class ToolsTransient {
    
    public:
    
    ToolsTransient();
    ToolsTransient(int pOrder,
                   REAL Lx, REAL Ly, REAL Lf, REAL Hf, REAL E, REAL Poiss, REAL Fx, REAL Fy, REAL Visc, REAL SigN,
                   REAL Qinj, REAL Ttot, REAL Nsteps, REAL Cl, REAL Pe, REAL SigmaConf, REAL Pref, REAL vsp, REAL KIc);
    
    ~ToolsTransient();
    
    //---------------------------------------------------------------
    
    void Run();
    
        TPZGeoMesh * Mesh2D(REAL lmax);
        TPZCompMesh * CMeshElastic(TPZGeoMesh *gmesh);
        void SolveInitialElasticity(TPZAnalysis &an, TPZCompMesh *Cmesh);
        TPZCompMeshReferred * CMeshReduced(TPZGeoMesh *gmesh, TPZCompMesh *cmesh);
        TPZCompMesh * CMeshPressure(TPZGeoMesh *gmesh);
        TPZCompMesh * MalhaCompMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, TPZNLFluidStructure2d * &mymaterial);
        
        void SolveSistTransient(REAL deltaT, REAL maxTime, TPZFMatrix<REAL> InitialSolution, TPZAnalysis *an,
                                TPZNLFluidStructure2d * &mymaterial, TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics);
    
    //---------------------------------------------------------------
    
    void MassMatrix(TPZNLFluidStructure2d *mymaterial, TPZCompMesh *mphysics, TPZFMatrix<REAL> & Un);
    
    void StiffMatrixLoadVec(TPZNLFluidStructure2d *mymaterial, TPZCompMesh* mphysics, TPZAnalysis *an,
                                   TPZAutoPointer< TPZMatrix<REAL> > & matK1, TPZFMatrix<REAL> &fvec);

    void SaidaPressao(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics);
    
    void SaidaMathPressao(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZNLFluidStructure2d * &mymaterial, std::stringstream & outP);

    void FillPositionPressure(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, std::map<REAL,REAL> & pos_pressure);
    
    void PosProcessMult(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis *an, std::string plotfile);
    
    void PlotWIntegral(TPZCompMesh *cmesh, std::stringstream & outW, int solNum);
    
    REAL ComputeKIPlaneStrain(TPZCompMesh * elastMesh, REAL young, REAL poisson,
                              REAL radius, std::stringstream & outFile, int cent = -1,
                              REAL TimeValue = -1, bool firstCall = true);
    
    bool PropagateOneStep(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZNLFluidStructure2d * &mymaterial);
    
    void CheckConv(TPZFMatrix<REAL> InitialSolution, TPZAnalysis *an, TPZNLFluidStructure2d * &mymaterial,
                          TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics);
    
    //---------------------------------------------------------------
    
    int fpOrder;
    
    InputDataStruct fInputData;
};

#endif
