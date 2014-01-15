//
//  pznlfluidstructureData.h
//  PZ
//
//  Created by Cesar Lucci on 05/06/13.
//
//

#ifndef PZ_TPZPlaneFractureData_h
#define PZ_TPZPlaneFractureData_h

#include <set>
#include <math.h>
#include "pzreal.h"
#include "pzerror.h"
#include "pzvec.h"
#include "pzcmesh.h"
#include "pzgeoel.h"
#include "pzcompel.h"
#include "pzintel.h"
#include <map>



struct TimeControl
{
public:
    TimeControl()
    {
        fTtot = 0.;
        factTime = 0.;
        fmaxDeltaT = 0.;
        fNDeltaTsteps = 0;
        fminDeltaT = 0.;
        factDeltaT = 0.;
    }
    
    ~TimeControl()
    {
        
    }
    
    void SetTimeControl(REAL Ttot, REAL maxDeltaT, int nTimes)
    {
        fTtot = Ttot;
        factTime = fminDeltaT;
        fmaxDeltaT = maxDeltaT;
        fNDeltaTsteps = nTimes;
        fminDeltaT = fmaxDeltaT/fNDeltaTsteps;
        factDeltaT = fminDeltaT;
    }
    
    void SetMinDeltaT()
    {
        factDeltaT = MIN(fminDeltaT,fTtot - factTime);
    }
    
    void SetNextDeltaT()
    {
        factDeltaT = MIN(fmaxDeltaT,(factDeltaT+fmaxDeltaT/fNDeltaTsteps));
        factDeltaT = MIN(factDeltaT,fTtot - factTime);
    }
    
    void UpdateActTime()
    {
        factTime += factDeltaT;
        std::cout << "\n\n=============== ActTime = " << factTime/60. << " min ===============\n\n";
    }
    
    REAL Ttot()
    {
        return fTtot;
    }
    
    REAL actTime()
    {
        return factTime;
    }
    
    REAL actDeltaT()
    {
        return factDeltaT;
    }
    
private:
    REAL fTtot;//Tempo total da simulacao
    REAL factTime;//tempo atual (em segundos)
    REAL fmaxDeltaT;//delta T maximo
    REAL fminDeltaT;//delta T minimo
    REAL factDeltaT;//delta T atual
    int fNDeltaTsteps;//quantidade de incrementos do deltaT para definir o deltaT minimo
};



struct LeakoffStorage
{
public:
    LeakoffStorage()
    {
        this->fGelId_Penetration.clear();
    }
    ~LeakoffStorage()
    {
        fGelId_Penetration.clear();
    }
    std::map<int,REAL> & GetLeakoffMap()
    {
        return fGelId_Penetration;
    }
    void SetLeakoffMap(std::map<int,REAL> & GelId_Penetration)
    {
        fGelId_Penetration = GelId_Penetration;
    }
    
    void UpdateLeakoff(TPZCompMesh * cmesh, int deltaT);
    
    REAL VlFtau(REAL pfrac, REAL tau, REAL Cl, REAL Pe, REAL gradPref, REAL vsp);
    
    REAL FictitiousTime(REAL VlAcum, REAL pfrac, REAL Cl, REAL Pe, REAL gradPref, REAL vsp);
    
    REAL QlFVl(int gelId, REAL pfrac, REAL deltaT, REAL Cl, REAL Pe, REAL gradPref, REAL vsp);
    
    REAL dQlFVl(int gelId, REAL pfrac, REAL deltaT, REAL Cl, REAL Pe, REAL gradPref, REAL vsp);
    
    std::map<int,REAL> fGelId_Penetration;
};




/**
 * To understand what is implemented here, see table (MaterialIds Table.xls) on folder (/NeoPZ/Projects/HydraulicFracturePropagation/PlaneFracture)
 */
struct MaterialIdGen
{
public:
    
    MaterialIdGen()
    {
        
    }
    
    ~MaterialIdGen()
    {
        
    }
    
    int Aux1DMatId()
    {
        return -1;
    }
    
    int CrackTipMatId()
    {
        return -2;
    }
    
    int RockMatId(int layer)
    {
#ifdef DEBUG
        if(layer < 0 || layer > 99)
        {
            //Soh pode ter 100 camadas (de 0 a 99)
            DebugStop();
        }
#endif
        
        return (layer+1)*10;
    }
    
    int BulletMatId(int layer)
    {
        return -RockMatId(layer);
    }
    
    int InsideFractMatId(int layer, int stripe)
    {
#ifdef DEBUG
        if(stripe < 0 || stripe > 9)
        {
            //Soh pode ter 10 faixas de pressao (de 0 a 9)
            DebugStop();
        }
#endif
        
        return -(RockMatId(layer) + 1000 + stripe);
    }
    
    int OutSideFractMatId(int layer)
    {
        return -(RockMatId(layer) + 2000);
    }
    
    int FarfieldMatId(int layer)
    {
        return -(RockMatId(layer) + 3000);
    }
    
    int LeftMatId(int layer)
    {
        return -(RockMatId(layer) + 4000);
    }
    
    int RightMatId(int layer)
    {
        return -(RockMatId(layer) + 5000);
    }
    
    int TopMatId()
    {
        return -6010;
    }
    
    int BottomMatId()
    {
        return -7010;
    }
    
    bool IsInsideFractMat(int matId)
    {
        if(matId <= -1010 && matId >= -2009)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    
    bool IsOutsideFractMat(int matId)
    {
        if(matId <= -2010 && matId >= -3000)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    
    bool IsBoundaryMaterial(int matId)
    {
        if(matId <= -3010 && matId >= -7010)
        {//is 2D BC between: left, right, farfield, top or bottom
            return true;
        }
        else
        {
            return false;
        }
    }
    
    bool IsBulletMaterial(int matId)
    {
        if(matId <= -10 && matId >= -1000)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    
    bool IsRockMaterial(int matId)
    {
        if(matId >= 10 && matId <= 1000)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    
    int WhatLayerFromInsideFracture(int insideMatId)
    {
#ifdef DEBUG
        if(IsInsideFractMat(insideMatId) == false)
        {//The given materialId IS NOT inside fracture
            DebugStop();
        }
#endif
        
        int inside = fabs(insideMatId);
        int stripe = inside - (inside/10)*10;
        
        int lay = (-1010-insideMatId-stripe)/10;
        
        return lay;
    }
    //------------------------------------------------------------------------------------------------------------
};

#include "pzanalysis.h"
#include "pzreal.h"
#include <map>
#include <fstream>



class Output3DDataStruct
{
public:
    
    Output3DDataStruct();
    ~Output3DDataStruct();

    void SetQinj1wing(REAL Qinj1wing);
    
    int NTimes();
    void InsertTAcumVolW(int time, REAL vol);
    void InsertTAcumVolLeakoff(int time, REAL vol);

    void PrintMathematica(std::ofstream & outf);
    
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
    
    REAL fQinj1wing;
    
    //maps indexed by time
    std::map<int,REAL> fTAcumVolW;
    std::map<int,REAL> fTAcumVolLeakoff;
};

extern TimeControl globTimeControl;

extern LeakoffStorage globLeakoffStorage;

extern MaterialIdGen globMaterialIdGen;

extern Output3DDataStruct globFractOutput3DData;

#endif
