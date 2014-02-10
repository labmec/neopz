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



const REAL globStressScale = 1.E-7;


class TimeControl
{
public:
    TimeControl()
    {
        fTtot = 0.;
        factTime = 0.;
        
        fDeltaT_left = 0.;
        fDeltaT_right = 0.;
        
        factDeltaT = 0.;
        
        fwasLeftLastTime = true;
        freachTime_right = false;
    }
    
    ~TimeControl()
    {
        
    }
    
    void SetTimeControl(REAL Ttot)
    {
        fTtot = Ttot;
        factTime = 0.;
        
        fDeltaT_left = 5.;
        fDeltaT_right = Ttot + 10.;
        
        fwasLeftLastTime = true;
        freachTime_right = false;
        
        ComputeActDeltaT();
    }
    
    void TimeisOnLeft()
    {
        fwasLeftLastTime = true;
        
        if(freachTime_right == false)
        {
            fDeltaT_left *= 2.;
        }
        else
        {
            fDeltaT_left = factDeltaT;
        }
    }
    
    void TimeisOnRight()
    {
        fwasLeftLastTime = false;
        
        if(freachTime_right == false)
        {//1st time reach time on right from KI=KIc
            fDeltaT_left = fDeltaT_left/2.;
            freachTime_right = true;
        }
        
        fDeltaT_right = factDeltaT;
    }
    
    void SetDeltaT(REAL deltaT)
    {
        factDeltaT = deltaT;
    }
    
    void ComputeActDeltaT()
    {
        if(freachTime_right == false)
        {
            factDeltaT = fDeltaT_left;
        }
        else
        {
            factDeltaT = (fDeltaT_left + fDeltaT_right)/2.;
        }
    }
    
    void UpdateActTime()
    {
        factTime += factDeltaT;
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
    
    void RestartBissection()
    {
        freachTime_right = false;
        factDeltaT = MAX(5.,factDeltaT/4.);
        fDeltaT_left = factDeltaT;
        fDeltaT_right = fTtot + 10.;
    }
    
    bool ReachEndOftime()
    {
        return (factTime >= fTtot - 0.01);
    }
    
    bool TimeLimitsIsCloseEnough()
    {
        if(fDeltaT_right < fDeltaT_left)
        {
            std::cout << "\n\n\nMetodo da bisseccao no tempo inverteu os limites!!!\n\n\n";
            DebugStop();
        }
        
        bool isCloseEnough = (fDeltaT_right - fDeltaT_left) < 1.;
        
        if(isCloseEnough && fwasLeftLastTime)
        {
            fDeltaT_left = fDeltaT_right;
            isCloseEnough = false;
        }
        
        return isCloseEnough;
    }
    
    REAL LeftDeltaT()
    {
        return fDeltaT_left;
    }
 
    REAL RightDeltaT()
    {
        return fDeltaT_right;
    }
    
private:
    REAL fTtot;//Tempo total da simulacao
    REAL factTime;//tempo atual (em segundos)
    REAL factDeltaT;//delta T atual
    
    REAL fDeltaT_left;//deltaT cujo factTime+dt nao propagou a fratura (serah utilizado no metodo da bisseccao).
    REAL fDeltaT_right;//deltaT cujo factTime+dt propagou a fratura (serah utilizado no metodo da bisseccao).
    
    bool fwasLeftLastTime;
    bool freachTime_right;
};



class LeakoffStorage
{
public:
    LeakoffStorage()
    {
        this->fGelId_Penetration.clear();
        this->fPressureIndependent = true;
        
        this->fDefaultLeakoffEnabled = true;
        this->fLeakoffEnabled = this->fDefaultLeakoffEnabled;
    }
    ~LeakoffStorage()
    {
        fGelId_Penetration.clear();
    }
    
    void SetPressureIndependent()
    {
        this->fPressureIndependent = true;
    }
    
    void SetPressureDependent()
    {
        this->fPressureIndependent = false;
    }
    
    bool IsPressureIndependent()
    {
        return this->fPressureIndependent;
    }
    
    std::map<int,REAL> & GetLeakoffMap()
    {
        return fGelId_Penetration;
    }
    
    void SetLeakoffMap(std::map<int,REAL> & GelId_Penetration)
    {
        fGelId_Penetration = GelId_Penetration;
    }
    
    void UpdateLeakoff(TPZCompMesh * cmesh, REAL deltaT);
    
    REAL VlFtau(REAL pfrac, REAL tau, REAL Cl, REAL Pe, REAL gradPref, REAL vsp);
    
    REAL FictitiousTime(REAL VlAcum, REAL pfrac, REAL Cl, REAL Pe, REAL gradPref, REAL vsp);
    
    REAL QlFVl(int gelId, REAL pfrac, REAL deltaT, REAL Cl, REAL Pe, REAL gradPref, REAL vsp);
    
    REAL dQlFVl(int gelId, REAL pfrac, REAL deltaT, REAL Cl, REAL Pe, REAL gradPref, REAL vsp);
    
    void Printleakoff(std::ofstream & outf);
    
    void SetDefaultLeakoffEnabled(bool isEnabled)
    {
        this->fDefaultLeakoffEnabled = isEnabled;
    }
    
    void DisableLeakoff()
    {
        fLeakoffEnabled = false;
    }
    
    void RestoreDefaultLeakoff()
    {
        fLeakoffEnabled = fDefaultLeakoffEnabled;
    }
    
protected:
    std::map<int,REAL> fGelId_Penetration;
    
    bool fPressureIndependent;
    
    bool fLeakoffEnabled;
    bool fDefaultLeakoffEnabled;
};




/**
 * To understand what is implemented here, see table (MaterialIds Table.xls) on folder (/NeoPZ/Projects/HydraulicFracturePropagation/PlaneFracture)
 */
class MaterialIdGen
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
    void InsertTAcumVolW(REAL time, REAL vol);
    void InsertTAcumVolLeakoff(REAL time, REAL vol);
    void InsertTL(REAL time, REAL L);
    void InsertTHsup(REAL time, REAL Hsup);
    void InsertTHinf(REAL time, REAL Hinf);

    void PrintMathematica(std::ofstream & outf);
        
    REAL fQinj1wing;
    
    //maps indexed by time
    std::map<REAL,REAL> fTAcumVolW;
    std::map<REAL,REAL> fTAcumVolLeakoff;
    std::map<REAL,REAL> fTL;
    std::map<REAL,REAL> fTHsup;
    std::map<REAL,REAL> fTHinf;
};

extern TimeControl globTimeControl;

extern LeakoffStorage globLeakoffStorage;

extern MaterialIdGen globMaterialIdGen;

extern Output3DDataStruct globFractOutput3DData;

#endif
