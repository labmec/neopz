
#include "TPZMatDualHybridPoisson.h"
#include "pzbndcond.h"
#include "pzaxestools.h"

//int TPZMatDualHybridPoisson::mydim = 2;

TPZMatDualHybridPoisson::TPZMatDualHybridPoisson(int nummat, REAL f, REAL betaZero)
: TPZRegisterClassId(&TPZMatDualHybridPoisson::ClassId), 
TPZMaterial(nummat),fXf(f), fBetaZero(betaZero){
   mydim = 0;
}

TPZMatDualHybridPoisson::TPZMatDualHybridPoisson(int matid) :
TPZRegisterClassId(&TPZMatDualHybridPoisson::ClassId), TPZMaterial(matid),
fXf(0.), fBetaZero(0.)
{
    mydim = 0;
    
}

TPZMatDualHybridPoisson::TPZMatDualHybridPoisson():
TPZRegisterClassId(&TPZMatDualHybridPoisson::ClassId), TPZMaterial(){
    mydim = 0;
}

TPZMatDualHybridPoisson::TPZMatDualHybridPoisson(const TPZMatDualHybridPoisson &copy)
: TPZRegisterClassId(&TPZMatDualHybridPoisson::ClassId), TPZMaterial(copy){
    mydim = copy.mydim;
    fXf = copy.fXf;
    fBetaZero = copy.fBetaZero;
}

TPZMatDualHybridPoisson::~TPZMatDualHybridPoisson(){
    
}

void TPZMatDualHybridPoisson::Print(std::ostream & out){
    out << "\n" << this->Name() << "\n";
    out << "fXf = " << fXf << "\n";
    out << "fBetaZero = " << fBetaZero << "\n";
    TPZMaterial::Print(out);
}

void TPZMatDualHybridPoisson::Contribute(TPZMaterialData &data,
                                         REAL weight,
                                         TPZFMatrix<STATE> &ek,
                                         TPZFMatrix<STATE> &ef){
    if(mydim==0) DebugStop();
    
    TPZFMatrix<REAL>  &phi = data.phi;
    TPZFMatrix<REAL> &dphi = data.dphix;
    TPZVec<REAL>  &x = data.x;
    const int nshape = phi.Rows();
    //    const REAL beta = this->Beta(data.p,data.HSize);
    
    if(dphi.Rows() == mydim-1){///estou no elemento 1D (ou 2D) do multiplicador de Lagrange
        
        return;
        
    }///1D
    else if(dphi.Rows() == mydim){///estou no elemento 2D (ou 3D)
        
        STATE Fval = fXf;
        if(this->HasForcingFunction()){
            TPZManVector<STATE,1> res(1);
            TPZFMatrix<STATE> dres(Dimension(),1);
            fForcingFunction->Execute(x,res,dres);
            Fval = res[0];
        }
        
        ///Equacao de Poisson
        for( int in = 0; in < nshape; in++ ) {
            ef(in, 0) +=  weight * Fval * phi(in,0);
            for( int jn = 0; jn < nshape; jn++ ) {
                //ek(in,jn) += weight * ( dphi(0,in) * dphi(0,jn) + dphi(1,in) * dphi(1,jn)  );
                
                REAL temp=0.;
                for(int k =0; k<mydim; k++){
                    temp +=  dphi(k,in) * dphi(k,jn);
                }
                ek(in,jn) += weight*temp;
            }
        }
        
        return;
        
    }///2D
    else DebugStop();
    
}///void

void TPZMatDualHybridPoisson::ContributeBC(TPZMaterialData &data,
                                           REAL weight,
                                           TPZFMatrix<STATE> &ek,
                                           TPZFMatrix<STATE> &ef,
                                           TPZBndCond &bc){
    
    TPZFMatrix<REAL>  &phi = data.phi;
    const int nshape = phi.Rows();
    STATE valBC;
    valBC = bc.Val2()(0,0);
    
    if(bc.HasForcingFunction()) {
        TPZManVector<STATE> res(1);
        bc.ForcingFunction()->Execute(data.x,res);
        valBC = res[0];
    }
    
    if(bc.Type() == 0)
    { // Dirichlet condition
        for(int i = 0; i < nshape; i++){
            ef(i,0) += weight * gBigNumber * phi(i,0) * valBC;
            for (int j = 0; j < nshape; j++){
                ek(i,j) += weight * gBigNumber * phi(i,0) * phi(j,0);
            }
        }
        return;
    }
    //        else{
    //            PZError << "TPZMatDualHybridPoisson::ContributeBC error. BC type not implemented\n";
    //            DebugStop();
    //        }
    
    
    if(bc.Type()==1){//Neumann
        for(int i = 0; i < nshape; i++) {
            ef(i,0) += weight*phi(i,0)*valBC;
        }
        return;
    }
    
    
}///void

void TPZMatDualHybridPoisson::ContributeInterface(TPZMaterialData &data,
                                                  TPZMaterialData &dataleft,
                                                  TPZMaterialData &dataright,
                                                  REAL weight,
                                                  TPZFMatrix<STATE> &ek,
                                                  TPZFMatrix<STATE> &ef){
    
    TPZFMatrix<REAL> &dphiLdAxes = dataleft.dphix;
    TPZFMatrix<REAL> &dphiRdAxes = dataright.dphix;
    TPZFMatrix<REAL> &phiL = dataleft.phi;
    TPZFMatrix<REAL> &phiR = dataright.phi;
    TPZManVector<REAL,3> &normal = data.normal;
    const REAL faceSize = data.HSize;
    const REAL beta = this->Beta(data.p,faceSize);
    
    const int nshapeL = phiL.Rows();
    const int nshapeR = phiR.Rows();
    
    TPZFNMatrix<660> dphiL, dphiR;
    TPZAxesTools<REAL>::Axes2XYZ(dphiLdAxes, dphiL, dataleft.axes);
    TPZAxesTools<REAL>::Axes2XYZ(dphiRdAxes, dphiR, dataright.axes);
    
    const REAL theta = +1;///para nao simetrico, colocar theta = -1
    
    if(dphiRdAxes.Rows() != mydim-1 || dphiLdAxes.Rows() != mydim) DebugStop(); //o multiplicador de Lagrange foi assumido como direito
    ///u, v : left
    ///lambda, mu : right
    
    ///-gradU.n v
    for(int il = 0; il < nshapeL; il++){
        for(int jl = 0; jl < nshapeL; jl++){
            //ek(il,jl) += weight * (-1.)* (dphiL(0,jl)*normal[0]+dphiL(1,jl)*normal[1]) * phiL(il,0);
            
            REAL temp = 0.;
            for(int k=0; k<mydim; k++){
                temp += dphiL(k,jl)*normal[k];
            }
            ek(il,jl) += weight * (-1.)* temp * phiL(il,0);
        }
    }
    
    ///-gradV.n u
    for(int il = 0; il < nshapeL; il++){
        for(int jl = 0; jl < nshapeL; jl++){
            //ek(il,jl) += weight * theta * (-1.)* (dphiL(0,il)*normal[0]+dphiL(1,il)*normal[1]) * phiL(jl,0);
            
            REAL temp = 0.;
            for(int k=0; k<mydim; k++){
                temp += dphiL(k,il)*normal[k];
            }
            ek(il,jl) += weight * theta * (-1.)* temp * phiL(jl,0);
        }
    }
    
    ///+gradV.n lambda
    for(int il = 0; il < nshapeL; il++){
        for(int jr = 0; jr < nshapeR; jr++){
            //ek(il,nshapeL+jr) += weight * theta * (+1.)* (dphiL(0,il)*normal[0]+dphiL(1,il)*normal[1]) * phiR(jr,0);
            
            REAL temp = 0.;
            for(int k=0; k<mydim; k++){
                temp += dphiL(k,il)*normal[k];
            }
            ek(il,nshapeL+jr) += weight * theta * (+1.)* temp * phiR(jr,0);
        }
    }
    
    ///beta u v
    for(int il = 0; il < nshapeL; il++){
        for(int jl = 0; jl < nshapeL; jl++){
            ek(il,jl) += weight * beta * phiL(il,0) * phiL(jl,0);
        }
    }
    
    ///-beta lambda v
    for(int il = 0; il < nshapeL; il++){
        for(int jr = 0; jr < nshapeR; jr++){
            ek(il,nshapeL+jr) += weight * (-1.) * beta * phiL(il,0) * phiR(jr,0);
        }
    }
    
    ///gradU.n mu - beta u mu
    for(int ir = 0; ir < nshapeR; ir++){
        for(int jl = 0; jl < nshapeL; jl++){
            //            ek(nshapeL+ir,jl) += weight * ( dphiL(0,jl)*normal[0] + dphiL(1,jl)*normal[1] ) * phiR(ir,0);
            //            ek(nshapeL+ir,jl) += weight * (-beta) * phiL(jl,0) * phiR(ir,0);
            
            REAL temp = 0.;
            for(int k=0; k<mydim; k++){
                temp += dphiL(k,jl)*normal[k];
            }
            ek(nshapeL+ir,jl) += weight * temp * phiR(ir,0);
            ek(nshapeL+ir,jl) += weight * (-beta) * phiL(jl,0) * phiR(ir,0);
            
        }
    }
    
    /// beta mu lambda
    for( int i = 0; i < nshapeR; i++ ) {
        for( int j = 0; j < nshapeR; j++ ) {
            ek(i+nshapeL,j+nshapeL) += weight * (beta) * phiR(i,0) * phiR(j,0);
        }
    }
    
    
}///void

void TPZMatDualHybridPoisson::ContributeBCInterface(TPZMaterialData &data,
                                                    TPZMaterialData &dataleft,
                                                    REAL weight,
                                                    TPZFMatrix<STATE> &ek,
                                                    TPZFMatrix<STATE> &ef,
                                                    TPZBndCond &bc){
    //PZError << "TPZMatDualHybridPoisson::ContributeBCInterface should never be called in this formulation\n";
    //DebugStop();
}

int TPZMatDualHybridPoisson::VariableIndex(const std::string &name){
    if(!strcmp("Solution",name.c_str())) return  1;
    if(!strcmp("ExactSolution",name.c_str())) return  2;
    if(!strcmp("Grad",name.c_str())) return  3;
    if(!strcmp("ExactGrad",name.c_str())) return  4;
    return TPZMaterial::VariableIndex(name);
}

int TPZMatDualHybridPoisson::NSolutionVariables(int var){
    if(var == 1 || var==99 || var==2) return 1;
    if(var == 3 || var==4) return Dimension();
    return TPZMaterial::NSolutionVariables(var);
}

void TPZMatDualHybridPoisson::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout){
    Solout.Resize( this->NSolutionVariables( var ) );
    
    if(var == 1){
        Solout[0] = data.sol[0][0];//solution - escalar
        return;
    }
    
    if(var == 99){
        Solout[0] = data.p;//solution - escalar
        return;
    }
    
    TPZVec<REAL> ptx(3);
    TPZVec<STATE> solExata(1);
    TPZFMatrix<STATE> flux(mydim,1);
    
    //Exact soluion
    if(var == 2){
        fExactSol->Execute(data.x, solExata,flux);
        Solout[0] = solExata[0];
        return;
    }//var6
    
    if(var == 3) {
        int id;
        for(id=0 ; id<Dimension(); id++) {
            TPZFNMatrix<9,STATE> dsoldx;
            TPZAxesTools<STATE>::Axes2XYZ(data.dsol[0], dsoldx, data.axes);
            Solout[id] = dsoldx(id,0);//derivate
        }
        return;
    }
    
    if(var == 4) {
        int id;
        fExactSol->Execute(data.x, solExata,flux);
        for(id=0 ; id<mydim; id++) {
            Solout[id] = flux(id,0);
        }
        return;
    }
    
    TPZMaterial::Solution(data,var,Solout);
}

void TPZMatDualHybridPoisson::Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,
                                     TPZFMatrix<STATE> &dudaxes, TPZFMatrix<REAL> &axes,
                                     TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values){
    values.Resize(3);
    values.Fill(0.);
    TPZFNMatrix<20,STATE> dudx(dudaxes.Rows(),dudaxes.Cols());
    TPZAxesTools<STATE>::Axes2XYZ(dudaxes, dudx, axes);
    if (dudx.Rows() < mydim) {
        // this is a lagrange multiplier element
        return;
    }
    ///L2 norm
    values[1] = (u[0] - u_exact[0])*(u[0] - u_exact[0]);
    ///semi norma de H1
    values[2] = 0.;
    for(int i = 0; i < mydim; i++){
        values[2] += (dudx(i,0) - du_exact(i,0))*(dudx(i,0) - du_exact(i,0));
    }
    ///H1 norm
    values[0] = values[1]+values[2];
}

int TPZMatDualHybridPoisson::ClassId() const{
    return Hash("TPZMatDualHybridPoisson") ^ TPZMaterial::ClassId() << 1;
}
