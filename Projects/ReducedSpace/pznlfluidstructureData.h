//
//  pznlfluidstructureData.h
//  PZ
//
//  Created by Cesar Lucci on 05/06/13.
//
//

#ifndef PZ_pznlfluidstructureData_h
#define PZ_pznlfluidstructureData_h

#include "pznlfluidstructureMaterials.h"
#include "pzreal.h"

class InputDataStruct
{
public:
    
    InputDataStruct();
    ~InputDataStruct();
    
    void SetData(REAL Lx, REAL Ly, REAL Lf, REAL Hf, REAL E, REAL Poisson, REAL Fx, REAL Fy, int NStripes, REAL Visc, REAL SigN,
                 REAL QinjTot, REAL Ttot, REAL maxDeltaT, int nTimes, REAL Cl, REAL Pe, REAL SigmaConf, REAL Pref, REAL vsp, REAL KIc);
    
    void SetLf(REAL Lf);
    
    REAL Lx();
    REAL Ly();
    REAL Lf();
    REAL Hf();
    REAL E();
    REAL Poisson();
    REAL Fx();
    REAL Fy();
    int NStripes();
    REAL Visc();
    REAL SigN();
    REAL Qinj();
    REAL Ttot();
    REAL actTime();
    REAL actDeltaT();
    REAL Cl();
    REAL Pe();
    REAL SigmaConf();
    REAL Pref();
    REAL vsp();
    REAL KIc();
    void SetMinDeltaT();
    void NextDeltaT();
    void NextActTime();
    
private:
    
    //Dimensions:
    REAL fLx;//Dimensao em x do domínio da malha do MEF
    REAL fLy;//Dimensao em y do domínio da malha do MEF
    REAL fLf;//Comprimento de 1/2 asa da fratura
    REAL fHf;//Altura da fratura
    
    //Elastic properties:
    REAL fE;//Modulo de elasticidade
    REAL fPoisson;//Poisson
    REAL fFx;//Bodyforces in x
    REAL fFy;//Bodyforces in y
    int fNStripes;//Amounth of pressure stripes for reduced space elastic references
    
    //Fluid property:
    REAL fVisc;//viscosidade do fluido de injecao
    
    //BCs:
    REAL fSigN;//Sigma.n no problema elastico que servira de espaco de aproximacao para o elastico multifisico
    REAL fQinj;//vazao de 1 asa de fratura dividido pela altura da fratura
    
    //time:
    REAL fTtot;//Tempo total da simulacao
    REAL factTime;//tempo atual (em segundos)
    REAL fmaxDeltaT;//delta T maximo
    REAL fminDeltaT;//delta T minimo
    REAL factDeltaT;//delta T atual
    int fNDeltaTsteps;//quantidade de incrementos do deltaT para definir o deltaT minimo
    
    //Leakoff:
    REAL fCl;//Carter
    REAL fPe;//Pressao estatica
    REAL fSigmaConf;//Tensao de confinamento
    REAL fPref;//Pressao de referencia da medicao do Cl
    REAL fvsp;//spurt loss
    
    //Propagation criterion
    REAL fKIc;
};


extern InputDataStruct globFractInputData;

#endif
