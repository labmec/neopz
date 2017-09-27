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

/**
 * @file
 * @brief Contains implementations of the TPZLinearWave methods.
 */

#include "tpzlinearwave.h"
#include "pzelmat.h"
#include "pzbndcond.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include "pzmaterialdata.h"
#include <math.h>
#include "pzlog.h"
#include "pzaxestools.h"

#include <cmath>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.material.linearwave"));
#endif

TPZLinearWave::EState TPZLinearWave::gState;


TPZLinearWave::TPZLinearWave(int nummat, int dim) : TPZDiscontinuousGalerkin(nummat), fXf(0.), fDim(dim) {
    fC = 1.;
    fDeltaT=1.0;
}

TPZLinearWave::TPZLinearWave():TPZDiscontinuousGalerkin(), fXf(0.), fDim(1){
    fC = 1.;
    fDeltaT=1.0;    
}

TPZLinearWave::TPZLinearWave(const TPZLinearWave &copy):TPZDiscontinuousGalerkin(copy){
    this->operator =(copy);
}

TPZLinearWave & TPZLinearWave::operator=(const TPZLinearWave &copy){
    TPZDiscontinuousGalerkin::operator = (copy);
    fXf  = copy.fXf;
    fDim = copy.fDim;
    fC   = copy.fC;
    fDeltaT = copy.fDeltaT;
    return *this;
}

void TPZLinearWave::SetParameters(STATE &Cspeed) {
    fC = Cspeed;
}

void TPZLinearWave::GetParameters(STATE &Cspeed) {
    Cspeed = fC;
}

TPZLinearWave::~TPZLinearWave() {
}

int TPZLinearWave::NStateVariables() {
    return 1;
}

void TPZLinearWave::Print(std::ostream &out) {
    out << "name of material : " << Name() << "\n";
    out << "Wave speed, laplace operator multiplier fC "<< fC << std::endl;
    out << "Forcing vector fXf " << fXf << std::endl;
    out << "Base Class properties :";
    TPZMaterial::Print(out);
    out << "\n";
}

// Tangent Matrix
void TPZLinearWave::Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) {
    
    TPZFMatrix<REAL>  &phi = data.phi;
    TPZFMatrix<REAL> &dphi = data.dphix;
    TPZVec<REAL>  &x = data.x;
    TPZFMatrix<REAL> &axes = data.axes;
    TPZFMatrix<REAL> &jacinv = data.jacinv;
    int phr = phi.Rows();

    // Linear Wave equation tangent matrix contribution
    for( int in = 0; in < phr; in++ ) 
    {
        for( int jn = 0; jn < phr; jn++ ) 
        {
            for(int kd=0; kd<fDim; kd++) 
            {
                ek(in,jn) += (STATE) weight *( (fDeltaT*fDeltaT) * (fC*fC)*(STATE)( dphi(kd,in) * dphi(kd,jn) ) 
                                       + (STATE)(phi(in) * phi(jn)) );
            }
        }
    }
    
    this->Contribute(data,weight,ef);    
    
}

// Residual Vector
void TPZLinearWave::Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ef) {
    
    TPZFMatrix<STATE>  &phi = data.phi;
    TPZFMatrix<STATE> &dphi = data.dphix;
    TPZVec<STATE>  &x = data.x;
    STATE sol_p =data.sol[0][0];
    TPZFMatrix<STATE> dsol_p =data.dsol[0];

    //  Compute grad(P)
    TPZManVector<STATE> GradientP(2,0);     
    GradientP[0] = dsol_p(0,0)*data.axes(0,0)+dsol_p(1,0)*data.axes(1,0);
    GradientP[1] = dsol_p(0,0)*data.axes(0,1)+dsol_p(1,0)*data.axes(1,1);        
    
    int phr = phi.Rows();

    
    if(gState == ENPlusOneState)
    {
        STATE fXfLoc = fXf;
        if(fForcingFunction) 
        {
            TPZManVector<STATE,1> res(1);
            TPZFMatrix<STATE> dres(Dimension(),1);
            fForcingFunction->Execute(x,res,dres);
            fXfLoc = res[0];
        }
        // Linear Wave equation Residual vector contribution at tn+1
        for( int in = 0; in < phr; in++ ) 
        {
            //  Compute grad(W)
            TPZManVector<STATE> GradientPhi(2,0);     
            GradientPhi[0] = dphi(0,in)*data.axes(0,0)+dphi(1,in)*data.axes(1,0);
            GradientPhi[1] = dphi(0,in)*data.axes(0,1)+dphi(1,in)*data.axes(1,1);
            
            STATE dot = (STATE) (GradientP[0] * GradientPhi[0] + GradientP[1] * GradientPhi[1]); 
            ef(in, 0) += (1.0) * (fDeltaT*fDeltaT) * (STATE) weight * ((fC*fC) * (STATE) dot - fXfLoc *((STATE) phi(in,0)));
            ef(in, 0) += 1.0 * (STATE) weight * sol_p *((STATE) phi(in,0));
        }
    }
 
     if(gState == ENState)
    {
        // Linear Wave equation Residual vector contribution at tn
        for( int in = 0; in < phr; in++ ) 
        {
            ef(in, 0) -= 2.0 * (STATE) weight * sol_p *((STATE) phi(in,0));
        }
    }
    
     if(gState == EMinusOneState)
    {
        // Linear Wave equation Residual vector contribution at tn-1
        for( int in = 0; in < phr; in++ ) 
        {
            ef(in, 0) += 1.0 * (STATE) weight * sol_p *((STATE) phi(in,0));
        }
    }    


}

void TPZLinearWave::ContributeBC(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
    TPZFMatrix<REAL>  &phi = data.phi;
    TPZVec<STATE>  &x = data.x;
    STATE sol_p =data.sol[0][0];
    TPZFMatrix<STATE> dsol_p =data.dsol[0];
    int phr = phi.Rows();
    short in,jn;
    STATE v2[1];
    v2[0] = bc.Val2()(0,0);
    
    if(bc.HasForcingFunction()) {            // phi(in, 0) = phi_in                          // JORGE 2013 01 26
        TPZManVector<STATE> res(1);
        bc.ForcingFunction()->Execute(data.x,res);       // dphi(i,j) = dphi_j/dxi
        v2[0] = res[0];
    }
    if(gState==ENPlusOneState)
    {
        switch (bc.Type()) {
        case 0 :            // Dirichlet condition
            for(in = 0 ; in < phr; in++) {
                ef(in,0) += (STATE)(gBigNumber* phi(in,0) * weight) *(sol_p -  v2[0]);
                for (jn = 0 ; jn < phr; jn++) 
                {
                    ek(in,jn) += gBigNumber * phi(in,0) * phi(jn,0) * weight;
                }
            }
            break;
        case 1 :            // Neumann condition
            for(in = 0 ; in < phi.Rows(); in++) {
                ef(in,0) -= v2[0] * (STATE)(phi(in,0) * weight);
            }
            break;
        default:
            {
                std::cout << "This BC doesn't exist." << std::endl;            
                DebugStop();            
            }
            break;
        }   
    }
}

/** Returns the variable index associated with the name */
int TPZLinearWave::VariableIndex(const std::string &name){
    if(!strcmp("Pressure",name.c_str()))        return  1;
    if(!strcmp("Velocity",name.c_str()))        return  2;
    if(!strcmp("WaveSpeed",name.c_str()))       return  3;
    return TPZMaterial::VariableIndex(name);
}

int TPZLinearWave::NSolutionVariables(int var){
    if(var == 1) return 1;
    if(var == 2) return fDim;
    if (var == 3) return 1;
    return TPZMaterial::NSolutionVariables(var);
}

void TPZLinearWave::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout){
    
    TPZVec<STATE> pressure(1);
    TPZVec<REAL> pto(3);
//     TPZFMatrix<STATE> flux(3,1);
      TPZFMatrix<STATE>  DSol = data.dsol[0], Axes = data.axes;
      STATE solP = data.sol[0][0];
    
    int numbersol = data.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }
    
    if(var == 1)
    {
        Solout[0] = solP;//Pressure
        return;
    }//var == 1
    
    if(var == 2) 
    {
        int id;
        for(id=0 ; id<fDim; id++) 
        {
            TPZFNMatrix<9,STATE> dsoldx;
            TPZAxesTools<STATE>::Axes2XYZ(DSol, dsoldx, Axes);
            Solout[id] = dsoldx(id,0);//Velocity
        }
        return;
    }//var == 2
    
    if (var == 3){
        Solout[0] =  this->fC;
        return;
    }//var ==3
    
}

void TPZLinearWave::Write(TPZStream &buf, int withclassid) const{
    TPZDiscontinuousGalerkin::Write(buf, withclassid);
    buf.Write(&fXf, 1);
    buf.Write(&fDim, 1);
    buf.Write(&fC, 1);
}

void TPZLinearWave::Read(TPZStream &buf, void *context){
    TPZDiscontinuousGalerkin::Read(buf, context);
    buf.Read(&fXf, 1);
    buf.Read(&fDim, 1);
    buf.Read(&fC, 1);
}