/*
    <this class implements a dynamic linear wave equation.>
    
    $ \frac{\delta^{2} u}{\delta t^{2}} - c^{2} \Delta u = f \text{ in } \Omega $

    With boundary conditions.
    
    $ u = 0 \text{ on } \delta \Omega \text{ or} $
    
    and initial conditions.
    
    $ u(0,x) = u_{0}(x) \text{ in } \Omega $
    
    
    Copyright (C) 2014  <copyright Omar> <omaryesiduran@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef TPZLINEARWAVE_H
#define TPZLINEARWAVE_H

/**
 * @file
 * @brief Contains the TPZLinearWave class.
 */

#include <iostream>
#include "pzdiscgal.h"
#include "pzfmatrix.h"

/**
 * @ingroup material
 * @brief this class implements a dynamic linear wave equation.
 */
/**
 * \f$ \frac{\delta^{2} u}{\delta t^{2}} - c^{2} \Delta u = f \text{ in } \Omega \f$
 */

/**
 * \f$ Bcs\;\; u = 0 \text{ on } \delta \Omega \text{ or} \f$
 */

/**
 * \f$ InitialCon\;\; u(0,x) = u_{0}(x) \text{ in } \Omega \f$
 */

/**
 * \f$ Weak formulation\;\; \int_{\Omega} \frac{\delta^{2} u} {\delta t^{2}}v dV + c^{2} \Big[ \int_{\Omega} \nabla u \cdot \nabla v dV 
- \int_{\delta \Omega} \nabla u \cdot \hat{n} \cdot v dA \Big] = \int fv \cdot dx. \f$
 */


class TPZLinearWave : public TPZDiscontinuousGalerkin 
{
    
protected:
    
    /** @brief Forcing function value */
    STATE fXf;
    
    /** @brief Problem dimension */
    int fDim;
    
    /** @brief Wave speed coefficient which multiplies the Laplace operator. */
    STATE fC;
    
    /** @brief Simulation time step */
    STATE fDeltaT;      
    
public:
    
     /** @brief Enumerate for timestep definitions */
    enum EState { EMinusOneState = 0, ENState = 1, ENPlusOneState = 2 };
    
    static EState gState;
    
    /** @brief Defines simulation time step. */
    void SetTimeStep(REAL timestep){ this->fDeltaT = timestep;}     
    
    /** @brief Defines simulation states. */    
    void SetMinusOneState(){ gState = EMinusOneState; }
    void SetNState(){ gState = ENState; }
    void SetPlusOneState(){ gState = ENPlusOneState; }    
    
    TPZLinearWave(int nummat, int dim);
    
    TPZLinearWave(int matid) : TPZDiscontinuousGalerkin(matid), fXf(0.), fDim(-1), fC(0.), fDeltaT(1.0)
    {
    
    }
    
    TPZLinearWave();
    
    TPZLinearWave(const TPZLinearWave &copy);
    
    virtual ~TPZLinearWave();
    
    TPZLinearWave &operator=(const TPZLinearWave &copy);
    
    virtual TPZMaterial * NewMaterial(){
        return new TPZLinearWave(*this);
    }
    
    /** 
     * @brief Fill material data parameter with necessary requirements for the
     * @since April 10, 2007
     */
    /** 
     * Contribute method. Here, in base class, all requirements are considered as necessary. 
     * Each derived class may optimize performance by selecting only the necessary data.
     */
    virtual void FillDataRequirements(TPZMaterialData &data)
    {
        data.SetAllRequirements(false);
        data.fNeedsSol=true;
    }
        
    /** @brief This method defines which parameters need to be initialized in order to compute the contribution of the boundary condition */
    virtual void FillBoundaryConditionDataRequirement(int type,TPZMaterialData &data)
    {
        data.SetAllRequirements(false);
        data.fNeedsSol=true;
    }
    
    
    int Dimension() const { return fDim;}
    
    int NStateVariables();
    
    void SetParameters(STATE &Cspeed);
    
    void GetParameters(STATE &Cspeed);
    
    void SetDimension(int dim)
    {
        if(dim<0 || dim >3)
        {
            DebugStop();
        }
        fDim = dim;
    }
    
    void SetInternalFlux(STATE flux)
    {
        fXf = flux;
    }
    
    
    virtual void Print(std::ostream & out);
    
    virtual std::string Name() { return "TPZLinearWave"; }

    /**
     * @name Contribute methods (weak formulation)
     * @{
     */
     
    virtual void Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ef);
    virtual void Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);
    virtual void ContributeBC(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);    
        
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) 
    {
        DebugStop();
    }
    
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
    {
        DebugStop();
    }

    
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) 
    {
        DebugStop();
    }
    
    virtual int VariableIndex(const std::string &name);
    
    virtual int NSolutionVariables(int var);
    
    virtual int NFluxes(){ return 3;}
    
protected:
    virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout) {
        TPZDiscontinuousGalerkin::Solution(datavec,var,Solout);
    }
    
    virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec) {
        TPZDiscontinuousGalerkin::FillDataRequirements(datavec);
    }
    
public:
    
    virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout);
    
    
//     virtual int ClassId() const{
//         return TPZLINEARWAVE;
//     }
    
    virtual void Write(TPZStream &buf, int withclassid) const;
    
    virtual void Read(TPZStream &buf, void *context);

};


#endif // TPZLINEARWAVE_H
