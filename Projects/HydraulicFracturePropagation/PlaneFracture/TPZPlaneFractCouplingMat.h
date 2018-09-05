//
//  TPZPlaneFractCouplingMat.h
//  PZ
//
//  Created by Cesar Lucci on 11/01/13.
//
//

#ifndef __PZ__TPZPlaneFractCouplingMat__
#define __PZ__TPZPlaneFractCouplingMat__

#include <iostream>


#include <iostream>
#include "TPZElast3Dnlinear.h"
#include "pzdiscgal.h"
#include "pzfunction.h"
#include "pzgmesh.h"
#include "pzcmesh.h"

class TPZLastElastFunction
{
public:
	
	/** @brief Class constructor */
	TPZLastElastFunction()
    {
        this->flastElastCMesh = NULL;
        this->finiElIndex = 0;
    }
    
    void SetLastElastCMesh(TPZCompMesh * LastElastCMesh)
    {
        this->flastElastCMesh = LastElastCMesh;
    }
	
	/** @brief Class destructor */
	~TPZLastElastFunction()
    {
        
    }
    
	/**
	 * @brief Performs function computation
	 * @param x point coordinate which is suppose to be in real coordinate system but can be in master coordinate system in derived classes.
	 * @param f function values
	 * @param df function derivatives
	 */
	virtual void Execute(TPZVec<REAL> &x, REAL &uy)
    {
        if(!this->flastElastCMesh)
        {
            uy = 0.;
            return;
        }
        
        this->flastElastCMesh->LoadReferences();
        
        TPZGeoMesh * gmesh = this->flastElastCMesh->Reference();
        TPZVec<REAL> qsi(3,0.);
        TPZGeoEl * gel = gmesh->FindElementCaju(x, qsi, this->finiElIndex, 3);
        
        if(!gel)
        {
            std::cout << "\n\n\ngeoEl not found on " << __PRETTY_FUNCTION__ << "\n\n\n";
            DebugStop();
        }
        if(!gel->Reference())
        {
            std::cout << "\n\n\nNULL geoEl->Reference() not found on " << __PRETTY_FUNCTION__ << "\n\n\n";
            DebugStop();
        }
        
        TPZVec<STATE> sol(3,0.);
        gel->Reference()->Solution(qsi,0,sol);
        uy = MAX(0.,sol[1]);
    }
    
    TPZCompMesh * flastElastCMesh;
    int64_t finiElIndex;
};


//===================================================================================


class TPZPlaneFractCouplingMat : public TPZElast3Dnlinear
{
protected:
    enum EState { EPastState = 0, EActualState = 1 };
	static EState gState;
    
public:
    TPZPlaneFractCouplingMat();
    TPZPlaneFractCouplingMat(int nummat, STATE E, STATE poisson, TPZVec<STATE> &force,
                             STATE preStressXX, STATE preStressYY, STATE preStressZZ,
                             STATE visc,
                             STATE Cl,
                             STATE Pe,
                             STATE gradPref,
                             STATE vsp);
    
    ~TPZPlaneFractCouplingMat();
    
    void SetLastElastCMesh(TPZCompMesh * LastElastCMesh);
    
    //Soh para nao ficar mostrando warnings!
    using TPZElast3Dnlinear::TPZMaterial::Contribute;
    using TPZElast3Dnlinear::TPZMaterial::ContributeBC;
    using TPZElast3Dnlinear::TPZMaterial::Solution;
    using TPZElast3Dnlinear::FillDataRequirements;
    
    virtual int NSolutionVariables(int var);
    virtual int VariableIndex(const std::string &name);
    
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec,
                            STATE weight,
                            TPZFMatrix<STATE> &ek,
                            TPZFMatrix<STATE> &ef);
    
    virtual void ContributePressure(TPZVec<TPZMaterialData> &datavec,
                                    REAL weight,
                                    TPZFMatrix<STATE> &ek,
                                    TPZFMatrix<STATE> &ef);
	
	virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec,
                              STATE weight,
                              TPZFMatrix<STATE> &ek,
                              TPZFMatrix<STATE> &ef,
                              TPZBndCond &bc);
    
    virtual void Solution(TPZVec<TPZMaterialData> &datavec,
                          int var,
                          TPZVec<STATE> &Solout);
    
    virtual void ApplyDirichlet_U(TPZVec<TPZMaterialData> &datavec,
                                  STATE weight,
                                  TPZFMatrix<STATE> &ek,
                                  TPZFMatrix<STATE> &ef,
                                  TPZBndCond &bc);
    
    virtual void ApplyNeumann_U(TPZVec<TPZMaterialData> &datavec,
                                STATE weight,
                                TPZFMatrix<STATE> &ek,
                                TPZFMatrix<STATE> &ef,
                                TPZBndCond &bc);
    
    virtual void ApplyBlockedDir_U(TPZVec<TPZMaterialData> &datavec,
                                   STATE weight,
                                   TPZFMatrix<STATE> &ek,
                                   TPZFMatrix<STATE> &ef,
                                   TPZBndCond &bc);
    
    virtual void ApplyDirichlet_P(TPZVec<TPZMaterialData> &datavec,
                                  STATE weight,
                                  TPZFMatrix<STATE> &ek,
                                  TPZFMatrix<STATE> &ef,
                                  TPZBndCond &bc);
    
    virtual void ApplyNeumann_P(TPZVec<TPZMaterialData> &datavec,
                                STATE weight,
                                TPZFMatrix<STATE> &ek,
                                TPZFMatrix<STATE> &ef,
                                TPZBndCond &bc);
    
    virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec);
    
    virtual void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec);
    
    static void SetPastState(){ gState = EPastState; }
	static void SetActualState(){ gState = EActualState; }
    static bool IsPastState(){ return (gState == EPastState); }
	static bool IsActualState(){ return (gState == EActualState); }

    
    REAL Cl()
    {
        return fCl;
    }
    REAL Pe()
    {
        return fPe;
    }
    REAL gradPref()
    {
        return fgradPref;
    }
    REAL vsp()
    {
        return fvsp;
    }

private:
    
    TPZLastElastFunction * fLastElastFunction;
    
    //Fluid
    STATE fVisc;
    
    //Leakoff
    REAL fCl;//Carter
    REAL fPe;//Pressao estatica
    REAL fgradPref;//Pressao de referencia da medicao do Cl
    REAL fvsp;//spurt loss
};




/*
class TPZPlaneFractBulletMat : public TPZDiscontinuousGalerkin
{
public:
    TPZPlaneFractBulletMat();
    
    TPZPlaneFractBulletMat(int nummat,
                           REAL Diameter,
                           REAL visc,
                           REAL Qinj_hbullet);
    
    ~TPZPlaneFractBulletMat();
    
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec,
                            STATE weight,
                            TPZFMatrix<STATE> &ek,
                            TPZFMatrix<STATE> &ef);
    
    virtual void ContributeInterface(TPZVec<TPZMaterialData> &datavec,
                                     TPZVec<TPZMaterialData> &dataleftvec,
                                     TPZVec<TPZMaterialData> &datarightvec,
                                     REAL weight,
                                     TPZFMatrix<STATE> &ek,
                                     TPZFMatrix<STATE> &ef);
    
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
    {
        DebugStop();
    }
    
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
    {
        DebugStop();
    }
    
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
                                       REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
    {
        DebugStop();
    }
    
    using TPZDiscontinuousGalerkin::FillDataRequirementsInterface;
    
    virtual void FillDataRequirementsInterface(TPZVec<TPZMaterialData> &datavec);
    
    virtual int Dimension()
    {
        return 1;
    }
    
    virtual int NStateVariables()
    {
        return 1;
    }
    
protected:
    
    REAL fDiameter;
    REAL fvisc;
    REAL fQinj_hbullet;
};
*/
#endif /* defined(__PZ__TPZPlaneFractCouplingMat__) */


