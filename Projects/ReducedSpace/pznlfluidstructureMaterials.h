//
//  pznlfluidstructureMaterials.h
//  PZ
//
//  Created by Cesar Lucci on 28/05/13.
//
//

#ifndef PZ_pznlfluidstructureMaterials_h
#define PZ_pznlfluidstructureMaterials_h

int const globReservMatId   = 1; //elastic
int const globPressureMatId = 2; //pressure
int const globMultiFisicMatId = 1;//multiphisics

int const globDirichletElastMatId = -1;
int const globBlockedXElastMatId  = -3;

int const globBCfluxIn  = -10; //bc pressure
int const globBCfluxOut = -20; //bc pressure

int const dirichlet = 0;
int const neumann   = 1;
int const mixed     = 2;

int const globDir_elast     = 11;
int const globMix_elast     = 20;
int const globNeum_pressure = 21;

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
                 REAL Qinj, REAL Ttot, REAL Nsteps, REAL Cl, REAL Pe, REAL SigmaConf, REAL Pref, REAL vsp, REAL KIc)
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
        fQinj = Qinj;

        fTtot = Ttot;
        fdeltaT = Ttot/Nsteps;

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

#endif
