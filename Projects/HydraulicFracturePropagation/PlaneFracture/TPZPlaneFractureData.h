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



const REAL globStressScale = 1.E-7;//DO NOT TOUCH!!!

class TimeControl
{
public:
    TimeControl()
    {
        fTtot = 0.;
        factTime = 0.;
        factDeltaT = 0.;
    }
    
    ~TimeControl()
    {
        
    }
    
    void SetTimeControl(REAL Ttot)
    {
        fTtot = Ttot;
        factTime = 0.;
    }
    
    void SetDeltaT(REAL deltaT)
    {
        factDeltaT = deltaT;
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
    
    bool ReachEndOftime()
    {
        return (factTime >= fTtot - 0.01);
    }
    
private:
    REAL fTtot;//Tempo total da simulacao
    REAL factTime;//tempo atual (em segundos)
    REAL factDeltaT;//delta T atual
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
        this->fLeakoffEnabled = isEnabled;
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
    
    int WhatStripe(int matId)
    {
        if(IsInsideFractMat(matId) == false)
        {
            std::cout << "\n\nGiven materialId (" << matId << ") does NOT belong to InsideFracture!!!\n\n";
            DebugStop();
        }
        
        int val = matId/10.;
        int lastDigit = matId - val*10;
        
        return fabs(lastDigit);//Last digit corresponds to the stripe
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


class LayerProperties
{
public:
    LayerProperties()
    {
        this->fYoung = 0.;
        this->fPoisson = 0.;
        this->fSigYY = 0.;
        this->fTVDini = 0.;
        this->fTVDfin = 0.;
        this->fKIc = 0.;
        this->fCl = 0.;
        this->fPe = 0.;
        this->fgradPref = 0.;
        this->fvsp = 0.;
    }
    
    LayerProperties(REAL Young, REAL Poisson, REAL SigYY, REAL TVDi, REAL TVDf,
                    REAL KIc, REAL Cl, REAL Pe, REAL gradPref, REAL vsp)
    {
#ifdef DEBUG
        REAL tol = 1.E-9;
        if(SigYY > 0.)
        {
            std::cout << "\n\nPre-stress must be null or compressive!\n";
            std::cout << "\n\nSee " << __PRETTY_FUNCTION__ << "\n\n";
            DebugStop();
        }
        if(Pe < 0. - tol)
        {
            std::cout << "\n\nStatic pressure (Pe) must be null or positive!\n";
            std::cout << "\n\nSee " << __PRETTY_FUNCTION__ << "\n\n";
            DebugStop();
        }
        if(gradPref < 0. + tol)
        {
            std::cout << "\n\nReference pressure (gradPref) must be positive!\n";
            std::cout << "\n\nSee " << __PRETTY_FUNCTION__ << "\n\n";
            DebugStop();
        }
#endif
        
        this->fYoung = Young;
        this->fPoisson = Poisson;
        this->fSigYY = SigYY;
        this->fTVDini = TVDi;
        this->fTVDfin = TVDf;
        this->fKIc = KIc;
        this->fCl = Cl;
        this->fPe = Pe;
        this->fgradPref = gradPref;
        this->fvsp = vsp;
    }
    LayerProperties(const LayerProperties & cp)
    {
        this->fYoung = cp.fYoung;
        this->fPoisson = cp.fPoisson;
        this->fSigYY = cp.fSigYY;
        this->fTVDini = cp.fTVDini;
        this->fTVDfin = cp.fTVDfin;
        this->fKIc = cp.fKIc;
        this->fCl = cp.fCl;
        this->fPe = cp.fPe;
        this->fgradPref = 0.;
        this->fvsp = cp.fvsp;
    }
    ~LayerProperties()
    {

    }
    
    //Elastic 3D
    REAL fYoung;
    REAL fPoisson;
    REAL fSigYY;
    
    //TVD limits
    REAL fTVDini;
    REAL fTVDfin;
    
    //SIF
    REAL fKIc;
    
    //leafoff
    REAL fCl;
    REAL fPe;
    REAL fgradPref;
    REAL fvsp;
};


class LayerStruct
{
public:
    LayerStruct()
    {
        this->fLayerVec.Resize(0);
        this->fLayer_Stripe_SolutionRow.clear();
        this->fElastReducedSolution.Resize(0,0);
        this->fminPrestress = +1.E15;
        this->fmaxPrestress = -1.E15;
    }
    
    ~LayerStruct()
    {
        this->fLayerVec.Resize(0);
    }
    
    const REAL StressAppliedOnEntireFracture()
    {
        return 1.E7 * globStressScale;//DO NOT TOUCH!!!
    }
    
    void ResetData()
    {
        this->fElastReducedSolution.Resize(0,0);
        this->fLayer_Stripe_SolutionRow.clear();
    }
    
    void SetLayerVec(TPZVec<LayerProperties> & LayerVec)
    {
        this->fLayerVec = LayerVec;
        for(int lay = 0; lay < LayerVec.NElements(); lay++)
        {
            this->fminPrestress = Min(this->fminPrestress,-LayerVec[lay].fSigYY);
            this->fmaxPrestress = MAX(this->fmaxPrestress,-LayerVec[lay].fSigYY);
        }
    }
    
    void SetLayerStripe_ContactSolutionRow(int layer, int stripe, int row)
    {
        std::map< int, std::map<int,int> >::iterator itLay = this->fLayer_Stripe_SolutionRow.find(layer);
        if(itLay == this->fLayer_Stripe_SolutionRow.end())
        {
            this->fLayer_Stripe_SolutionRow[layer][stripe] = row;
        }
        else
        {
            std::map<int,int>::iterator itStripe = itLay->second.find(stripe);
            if(itStripe == itLay->second.end())
            {
                (itLay->second)[stripe] = row;
            }
            else
            {
                std::cout << "\n\n\nIncluindo 2 vezes a mesma layer e stripe???\n\n\n";
                DebugStop();
            }
        }
    }
    
    void SetElastSolutionMatrix(TPZFMatrix<REAL> & solution)
    {
        this->fElastReducedSolution = solution;
    }
    
    REAL GetStressAppliedForJIntegral(int layer, int stripe)
    {
        int pressAppliedRow = 1;
        //Como agora eh aplicado newman unitario, o alpha corresponde aa pressao aplicada!
        REAL pressApplied = this->fElastReducedSolution(pressAppliedRow,0);
        
        int contactRow = this->GetContactSolutionRow(layer,stripe);
        //Como agora eh aplicado newman unitario, o alpha corresponde aa pressao aplicada!
        REAL contactStress = this->fElastReducedSolution(contactRow,0);
        //Onde nao ha contato, contactStress=0.
        
        REAL preStress = -this->fLayerVec[layer].fSigYY;
        
        REAL cellStressApllied = pressApplied + contactStress - preStress;
        
        return MAX(0.,cellStressApllied);
    }
    
    REAL GetLowerPreStress()
    {
        return this->fminPrestress;
    }
    
    REAL GetHigherPreStress()
    {
        return this->fmaxPrestress;
    }
    
    int NLayers()
    {
        return this->fLayerVec.NElements();
    }
    
    const LayerProperties & GetLayer(int index)
    {
        return this->fLayerVec[index];
    }
    
    const LayerProperties & GetLayerFromZcoord(REAL zCoord)
    {
        int whatLayer = this->WhatLayer(zCoord);
        return this->fLayerVec[whatLayer];
    }
    
    int WhatLayer(int zCoord)
    {
        for(int lay = 0; lay < this->fLayerVec.NElements(); lay++)
        {
            if(fabs(zCoord) > this->fLayerVec[lay].fTVDini && fabs(zCoord) < this->fLayerVec[lay].fTVDfin)
            {
                return lay;
            }
        }
        DebugStop();//nao achou o layer
        return -1;
    }
    
    int GetContactSolutionRow(int layer, int stripe)
    {
        int contactSolutionRow = -1;
        
        std::map< int, std::map<int,int> >::iterator itLay = this->fLayer_Stripe_SolutionRow.find(layer);
        if(itLay == this->fLayer_Stripe_SolutionRow.end())
        {
            std::cout << "\n\n\nLayer nao encontrado no metodo " << __PRETTY_FUNCTION__ << ".\n";
            DebugStop();
        }
        else
        {
            std::map<int,int>::iterator itStripe = itLay->second.find(stripe);
            if(itStripe != itLay->second.end())
            {
                contactSolutionRow = itStripe->second;
            }
        }
        
        return contactSolutionRow;
    }
    
protected:

    /** Vetor de camadas. Posicao 0: camada mais acima. Ultima posicao: camada mais abaixo */
    TPZVec<LayerProperties> fLayerVec;
    
    /** Mapa que guarda a linha da solucao do contato baseado no layer e na stripe */
    std::map< int, std::map<int,int> > fLayer_Stripe_SolutionRow;
    
    /** Matriz_vetor solucao da elastica de espacos reduzidos */
    TPZFMatrix<REAL> fElastReducedSolution;
    
    REAL fminPrestress;
    REAL fmaxPrestress;
};


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
    std::stringstream fFractContour;
    
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

extern LayerStruct globLayerStruct;

#endif
