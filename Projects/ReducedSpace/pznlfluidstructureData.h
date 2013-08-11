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
#include "pzanalysis.h"
#include "pzreal.h"
#include <map.h>
#include <fstream.h>

class InputDataStruct
{
public:
    
    InputDataStruct();
    ~InputDataStruct();
    
    void SetData(REAL Lx, REAL Ly, REAL Lf, REAL Hf, REAL E, REAL Poisson, REAL Fx, REAL Fy,
                 int NStripes, REAL Visc, REAL SigN, REAL QinjTot, REAL Ttot, REAL maxDeltaT, int nTimes,
                 REAL Cl, REAL Pe, REAL SigmaConf, REAL Pref, REAL vsp, REAL KIc, REAL Jradius);
    
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
    REAL Jradius();
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
    REAL fJradius;
    REAL fKIc;
};




class OutputDataStruct
{
public:
    
    OutputDataStruct();
    ~OutputDataStruct();
    
    int NTimes();
    void InsertTposP(int time, std::map<REAL,REAL> & posPmap);
    void InsertTposVolLeakoff(int time, REAL pos, REAL Ql);
    void InsertTAcumVolW(int time, REAL vol);
    void InsertTAcumVolLeakoff(int time, REAL vol);
    void InsertTKI(int time, REAL KI);
    void SetQinj1WingAndLfracmax(REAL Qinj1wing, REAL Lfracmax);
    
    void PlotElasticVTK(TPZAnalysis * an, int anCount = -1);
    void PrintMathematica(std::ofstream & outf);
    
    struct posP
    {
    public:
        
        posP()
        {
            fposP.clear();
        }
        ~posP()
        {
            fposP.clear();
        }
        void PrintMathematica(std::ofstream & outf)
        {
#ifdef DEBUG
            if(fposP.size() == 0)
            {
                DebugStop();
            }
#endif
            std::map<REAL,REAL>::iterator itposP;
            std::map<REAL,REAL>::iterator itposPLast = fposP.end();
            itposPLast--;
            
            outf << "{";
            for(itposP = fposP.begin(); itposP != fposP.end(); itposP++)
            {
                outf << "{" << itposP->first << "," << itposP->second << "}";
                if(itposP != itposPLast)
                {
                    outf << ",";
                }
            }
            outf << "}";
        }
        
        std::map<REAL,REAL> fposP;
    };
    
    struct posVolLeakoff
    {
    public:
        posVolLeakoff()
        {
            fposVolLeakoff.clear();
        }
        ~posVolLeakoff()
        {
            fposVolLeakoff.clear();
        }
        void InsertPoint(REAL pos, REAL Ql)
        {
            fposVolLeakoff[pos] = Ql;
        }
        void PrintMathematica(std::ofstream & outf)
        {
#ifdef DEBUG
            if(fposVolLeakoff.size() == 0)
            {
                DebugStop();
            }
#endif
            std::map<REAL,REAL>::iterator itposVolLeakoff;
            std::map<REAL,REAL>::iterator itposVolLeakoffLast = fposVolLeakoff.end();
            itposVolLeakoffLast--;
            
            outf << "{";
            for(itposVolLeakoff = fposVolLeakoff.begin(); itposVolLeakoff != fposVolLeakoff.end(); itposVolLeakoff++)
            {
                outf << "{" << itposVolLeakoff->first << "," << itposVolLeakoff->second << "}";
                if(itposVolLeakoff != itposVolLeakoffLast)
                {
                    outf << ",";
                }
            }
            outf << "}";
        }
        
        std::map<REAL,REAL> fposVolLeakoff;
    };
    
    //maps indexed by time
    std::map<int,posP> fTposP;
    std::map<int,posVolLeakoff> fTposVolLeakoff;
    std::map<int,REAL> fTAcumVolW;
    std::map<int,REAL> fTAcumVolLeakoff;
    std::map<int,REAL> fTKI;
    REAL fQinj1wing;
    REAL fLfracMax;
};


extern InputDataStruct globFractInputData;

extern OutputDataStruct globFractOutputData;

#endif
