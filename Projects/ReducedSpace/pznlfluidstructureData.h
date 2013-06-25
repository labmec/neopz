//
//  pznlfluidstructureData.h
//  PZ
//
//  Created by Cesar Lucci on 05/06/13.
//
//

#ifndef PZ_pznlfluidstructureData_h
#define PZ_pznlfluidstructureData_h

class InputDataStruct
{
public:
    
    InputDataStruct()
    {
        
    }
    ~InputDataStruct()
    {
        
    }
    
    void SetData(REAL Lx, REAL Ly, REAL Lf, REAL Hf, REAL E, REAL Poisson, REAL Fx, REAL Fy, REAL Visc, TPZVec<REAL> & SigN,
                 REAL QinjTot, REAL Ttot, REAL deltaT, REAL Cl, REAL Pe, REAL SigmaConf, REAL Pref, REAL vsp, REAL KIc)
    {
        fLx = Lx;
        fLy = Ly;
        fLf = Lf;
        fHf = Hf;
        
        fE = E;
        fPoisson = Poisson;
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
    void DeltaT(REAL deltaT)
    {
        fdeltaT = deltaT;
    }
    
    REAL Lx() { return fLx; }
    REAL Ly() { return fLy; }
    REAL Lf() { return fLf; }
    REAL Hf() { return fHf; }
    REAL E() { return fE; }
    REAL Poisson() { return fPoisson; }
    REAL Fx() { return fFx; }
    REAL Fy() { return fFy; }
    REAL Visc() { return fVisc; }
    REAL SigN(int pos) { return fSigN[pos]; }
    REAL Qinj() { return fQinj; }
    REAL Ttot() { return fTtot; }
    REAL deltaT() { return fdeltaT; }
    REAL Cl() { return fCl; }
    REAL Pe() { return fPe; }
    REAL SigmaConf() { return fSigmaConf; }
    REAL Pref() { return fPref; }
    REAL vsp() { return fvsp; }
    REAL KIc() { return fKIc; }
    
//private:
    
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
    
    //Fluid property:
    REAL fVisc;//viscosidade do fluido de injecao
    
    //BCs:
    TPZVec<REAL> fSigN;//Sigma.n no problema elastico que servira de espaco de aproximacao para o elastico multifisico
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

#endif
