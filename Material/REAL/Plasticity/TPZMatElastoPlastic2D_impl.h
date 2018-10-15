//
//  TPZMatElastoPlastic2D_impl.h
//  pz
//
//  Created by Omar Durán on 9/8/18.
//

#include "TPZMatElastoPlastic2D.h"
#include "TPZMatElastoPlastic_impl.h"

#include "pzbndcond.h"

#ifndef WIN32
#include <fenv.h>//NAN DETECTOR
#endif

#ifdef LOG4CXX
#include "pzlog.h"
static LoggerPtr elastoplastic2dLogger(Logger::getLogger("material.pzElastoPlastic2D"));
#endif


template <class T, class TMEM>
TPZMatElastoPlastic2D<T,TMEM>::TPZMatElastoPlastic2D() : TPZMatElastoPlastic<T,TMEM>()
{
    fPlaneStrain = true;
}

template <class T, class TMEM>
TPZMatElastoPlastic2D<T,TMEM>::TPZMatElastoPlastic2D(int id , int PlaneStrainOrPlaneStress) : TPZMatElastoPlastic<T,TMEM>(id)
{
    fPlaneStrain = PlaneStrainOrPlaneStress;
}

template <class T, class TMEM>
TPZMatElastoPlastic2D<T,TMEM>::TPZMatElastoPlastic2D(const TPZMatElastoPlastic2D<T,TMEM> &mat) : TPZMatElastoPlastic<T,TMEM>(mat), fPlaneStrain(mat.fPlaneStrain)
{
#ifdef LOG4CXX
    if(elastoplastic2dLogger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << ">>> TPZMatElastoPlastic2D<T,TMEM>() copy constructor called ***";
        LOGPZ_INFO(elastoplastic2dLogger,sout.str().c_str());
    }
#endif
}

template <class T, class TMEM>
TPZMatElastoPlastic2D<T,TMEM>::~TPZMatElastoPlastic2D()
{
    
}

template <class T, class TMEM>
void TPZMatElastoPlastic2D<T,TMEM>::ApplyDeltaStrain(TPZMaterialData & data, TPZFMatrix<REAL> & DeltaStrain,TPZFMatrix<REAL> & Stress)
{
    
    
#ifdef PZDEBUG
    if (DeltaStrain.Rows() != 6) {
        DebugStop();
    }
#endif
    
    if (!fPlaneStrain) //
    {//
        DebugStop();//PlaneStress
    }
    
    TPZMatElastoPlastic<T,TMEM>::ApplyDeltaStrain(data,DeltaStrain,Stress);//
    
}


template <class T, class TMEM>
void TPZMatElastoPlastic2D<T,TMEM>::ApplyDeltaStrainComputeDep(TPZMaterialData & data, TPZFMatrix<REAL> & DeltaStrain,TPZFMatrix<REAL> & Stress, TPZFMatrix<REAL> & Dep)
{
#ifdef PZDEBUG
    if (DeltaStrain.Rows() != 6) {
        DebugStop();
    }
#endif
    if (!fPlaneStrain) //
    {//
        DebugStop();//PlaneStress
    }
    TPZMatElastoPlastic<T,TMEM>::ApplyDeltaStrainComputeDep(data,DeltaStrain,Stress,Dep);//
}

template <class T, class TMEM>
void TPZMatElastoPlastic2D<T, TMEM>::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef) {
    
    TPZFMatrix<REAL> &dphi = data.dphix, dphiXY;
    TPZFMatrix<REAL> &phi = data.phi;
    TPZFMatrix<REAL> &axes = data.axes, axesT;
    
    axes.Transpose(&axesT);
    axesT.Multiply(dphi, dphiXY);
    
    const int phr = phi.Rows();
    
    TPZFNMatrix<4> Deriv(2, 2);
    TPZFNMatrix<36> Dep(6, 6,0.0);
    TPZFNMatrix<6> DeltaStrain(6, 1);
    TPZFNMatrix<6> Stress(6, 1);
    int ptindex = data.intGlobPtIndex;
    
    if (TPZMatWithMem<TMEM>::fUpdateMem && data.sol.size() > 0) {
        // Loop over the solutions if update memory is true
        TPZSolVec locsol(data.sol);
        TPZGradSolVec locdsol(data.dsol);
        int numsol = locsol.size();
        
        for (int is = 0; is < numsol; is++) {
            data.sol[0] = locsol[is];
            data.dsol[0] = locdsol[is];
            
            this->ComputeDeltaStrainVector(data, DeltaStrain);
            this->ApplyDeltaStrainComputeDep(data, DeltaStrain, Stress, Dep);
        }
    } else {
        this->ComputeDeltaStrainVector(data, DeltaStrain);
        this->ApplyDeltaStrainComputeDep(data, DeltaStrain, Stress, Dep);
    }
    
#ifdef MACOS
    feclearexcept(FE_ALL_EXCEPT);
    if (fetestexcept(/*FE_DIVBYZERO*/ FE_ALL_EXCEPT)) {
        std::cout << "division by zero reported\n";
        DebugStop();
    }
#endif
    
#ifdef LOG4CXX
    if (elastoplastic2dLogger->isDebugEnabled()) {
        std::stringstream sout;
        sout << ">>> TPZMatElastoPlastic<T,TMEM>::Contribute ***";
        sout << "\nIntegration Local Point index = " << data.intGlobPtIndex;
        sout << "\nIntegration Global Point index = " << data.intGlobPtIndex;
        sout << "\ndata.axes = " << data.axes;
        sout << "\nDep " << endl;
        sout << Dep(_XX_, _XX_) << "\t" << Dep(_XX_, _YY_) << "\t" << Dep(_XX_, _XY_) << "\n";
        sout << Dep(_YY_, _XX_) << "\t" << Dep(_YY_, _YY_) << "\t" << Dep(_YY_, _XY_) << "\n";
        sout << Dep(_XY_, _XX_) << "\t" << Dep(_XY_, _YY_) << "\t" << Dep(_XY_, _XY_) << "\n";
        
        sout << "\nStress " << endl;
        sout << Stress(_XX_, 0) << "\t" << Stress(_YY_, 0) << "\t" << Stress(_XY_, 0) << "\n";
        
        sout << "\nDELTA STRAIN " << endl;
        sout << DeltaStrain(0, 0) << "\t" << DeltaStrain(1, 0) << "\t" << DeltaStrain(2, 0) << "\n";
        sout << "data.phi" << data.phi;
        
        LOGPZ_DEBUG(elastoplastic2dLogger, sout.str().c_str());
    }
#endif
    ptindex = 0;
    int nstate = NStateVariables();
    REAL val;
    
    TPZManVector<STATE, 5> ForceLoc(this->fForce);
    if (this->fForcingFunction) {
        this->fForcingFunction->Execute(data.x, ForceLoc);
    }
    
    int in;
    for (in = 0; in < phr; in++) {
        
        val = ForceLoc[0] * phi(in, 0);
        val -= Stress(_XX_, 0) * dphiXY(0, in);
        val -= Stress(_XY_, 0) * dphiXY(1, in);
        ef(in * nstate + 0, 0) += weight * val;
        
        val = ForceLoc[1] * phi(in, 0);
        val -= Stress(_XY_, 0) * dphiXY(0, in);
        val -= Stress(_YY_, 0) * dphiXY(1, in);
        ef(in * nstate + 1, 0) += weight * val;
        
        for (int jn = 0; jn < phr; jn++) {
            for (int ud = 0; ud < 2; ud++) {
                for (int vd = 0; vd < 2; vd++) {
                    Deriv(vd, ud) = dphiXY(vd, in) * dphiXY(ud, jn);
                }
            }
            
            val = 2. * Dep(_XX_, _XX_) * Deriv(0, 0); //dvdx*dudx
            val += Dep(_XX_, _XY_) * Deriv(0, 1); //dvdx*dudy
            val += 2. * Dep(_XY_, _XX_) * Deriv(1, 0); //dvdy*dudx
            val += Dep(_XY_, _XY_) * Deriv(1, 1); //dvdy*dudy
            val *= 0.5;
            ek(in * nstate + 0, jn * nstate + 0) += weight * val;
            
            val = Dep(_XX_, _XY_) * Deriv(0, 0);
            val += 2. * Dep(_XX_, _YY_) * Deriv(0, 1);
            val += Dep(_XY_, _XY_) * Deriv(1, 0);
            val += 2. * Dep(_XY_, _YY_) * Deriv(1, 1);
            val *= 0.5;
            ek(in * nstate + 0, jn * nstate + 1) += weight * val;
            
            val = 2. * Dep(_XY_, _XX_) * Deriv(0, 0);
            val += Dep(_XY_, _XY_) * Deriv(0, 1);
            val += 2. * Dep(_YY_, _XX_) * Deriv(1, 0);
            val += Dep(_YY_, _XY_) * Deriv(1, 1);
            val *= 0.5;
            ek(in * nstate + 1, jn * nstate + 0) += weight * val;
            
            val = Dep(_XY_, _XY_) * Deriv(0, 0);
            val += 2. * Dep(_XY_, _YY_) * Deriv(0, 1);
            val += Dep(_YY_, _XY_) * Deriv(1, 0);
            val += 2. * Dep(_YY_, _YY_) * Deriv(1, 1);
            val *= 0.5;
            ek(in * nstate + 1, jn * nstate + 1) += weight * val;
        }
    }
    
#ifdef LOG4CXX
    if (elastoplastic2dLogger->isDebugEnabled()) {
        std::stringstream sout;
        sout << "<<< TPZMatElastoPlastic2D<T,TMEM>::Contribute ***";
        sout << " Resultant rhs vector:\n" << ef;
        sout << " Resultant stiff vector:\n" << ek;
        LOGPZ_DEBUG(elastoplastic2dLogger, sout.str().c_str());
    }
#endif
}

template <class T, class TMEM>
void TPZMatElastoPlastic2D<T, TMEM>::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ef) {
    TPZFMatrix<REAL> &dphi = data.dphix;
    TPZFMatrix<REAL> &phi = data.phi;
    TPZFMatrix<REAL> &axes = data.axes;
    TPZFNMatrix<9, REAL> axesT;
    TPZFNMatrix<50, REAL> dphiXY;
    
    axes.Transpose(&axesT);
    axesT.Multiply(dphi, dphiXY);
    
    const int phr = phi.Rows();
    
    //TPZFNMatrix<36> Deriv(6, 6);
    TPZFNMatrix<6> DeltaStrain(6, 1);
    TPZFNMatrix<6> Stress(6, 1);
    int ptindex = data.intGlobPtIndex;
    
    if (TPZMatWithMem<TMEM>::fUpdateMem && data.sol.size() > 0) {
        // Loop over the solutions if update memory is true
        //TPZFNMatrix<9> Dep(3, 3);
        
        TPZSolVec locsol(data.sol);
        TPZGradSolVec locdsol(data.dsol);
        int numsol = locsol.size();
        
        for (int is = 0; is < numsol; is++) {
            data.sol[0] = locsol[is];
            data.dsol[0] = locdsol[is];
            
            this->ComputeDeltaStrainVector(data, DeltaStrain);
            this->ApplyDeltaStrain(data, DeltaStrain, Stress);
            //this->ApplyDeltaStrainComputeDep(data, DeltaStrain, Stress, Dep);
        }
    } else {
        this->ComputeDeltaStrainVector(data, DeltaStrain);
        this->ApplyDeltaStrain(data, DeltaStrain, Stress);
        //        this->ApplyDeltaStrainComputeDep(data, DeltaStrain, Stress, Dep);
    }
#ifdef MACOS
    feclearexcept(FE_ALL_EXCEPT);
    if (fetestexcept(/*FE_DIVBYZERO*/ FE_ALL_EXCEPT)) {
        std::cout << "division by zero reported\n";
        DebugStop();
    }
#endif
    
#ifdef LOG4CXX
    if (elastoplastic2dLogger->isDebugEnabled()) {
        std::stringstream sout;
        sout << ">>> TPZMatElastoPlastic<T,TMEM>::Contribute ***";
        sout << "\nIntegration Local Point index = " << data.intGlobPtIndex;
        sout << "\nIntegration Global Point index = " << data.intGlobPtIndex;
        sout << "\ndata.axes = " << data.axes;
        sout << "\nStress " << endl;
        sout << Stress(_XX_, 0) << "\t" << Stress(_YY_, 0) << "\t" << Stress(_XY_, 0) << "\n";
        sout << "\nDELTA STRAIN " << endl;
        sout << DeltaStrain(0, 0) << "\t" << DeltaStrain(1, 0) << "\t" << DeltaStrain(2, 0) << "\n";
        sout << "data.phi" << data.phi;
        
        LOGPZ_DEBUG(elastoplastic2dLogger, sout.str().c_str());
    }
#endif
    ptindex = 0;
    int nstate = NStateVariables();
    REAL val;
    
    TPZManVector<STATE, 3> ForceLoc(this->fForce);
    if (this->fForcingFunction) {
        this->fForcingFunction->Execute(data.x, ForceLoc);
    }
    
    int in;
    for (in = 0; in < phr; in++) {
        val = ForceLoc[0] * phi(in, 0);
        val -= Stress(_XX_, 0) * dphiXY(0, in);
        val -= Stress(_XY_, 0) * dphiXY(1, in);
        ef(in * nstate + 0, 0) += weight * val;
        
        val = ForceLoc[1] * phi(in, 0);
        val -= Stress(_XY_, 0) * dphiXY(0, in);
        val -= Stress(_YY_, 0) * dphiXY(1, in);
        ef(in * nstate + 1, 0) += weight * val;
    }
    
    
#ifdef LOG4CXX
    if (elastoplastic2dLogger->isDebugEnabled()) {
        std::stringstream sout;
        sout << "<<< TPZMatElastoPlastic2D<T,TMEM>::Contribute ***";
        sout << " Resultant rhs vector:\n" << ef;
        LOGPZ_DEBUG(elastoplastic2dLogger, sout.str().c_str());
    }
#endif
    
}


template <class T, class TMEM>
void TPZMatElastoPlastic2D<T,TMEM>::FillBoundaryConditionDataRequirement(int type,TPZMaterialData &data)
{
    
    
    //TPZMatWithMem<TMEM>::FillBoundaryConditionDataRequirement(type,data);
    data.fNeedsSol = true;
    if (type == 4 || type ==5 || type == 6) {
        data.fNeedsNormal = true;
    }
    else {
        data.fNeedsNormal = false;
    }
}

template <class T, class TMEM>
void TPZMatElastoPlastic2D<T, TMEM>::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef, TPZBndCond &bc) {
    TPZFMatrix<REAL> &phi = data.phi;
    const REAL BIGNUMBER = TPZMaterial::gBigNumber;
    int nstate = NStateVariables();
    
    const int phr = phi.Rows();
    int in, jn, idf, jdf;
    REAL v2[2];
    v2[0] = bc.Val2()(0, 0);
    v2[1] = bc.Val2()(1, 0);
    
    TPZFMatrix<REAL> &v1 = bc.Val1();
    switch (bc.Type()) {
        case 0: // Dirichlet condition
            for (in = 0; in < phr; in++) {
                ef(nstate * in + 0, 0) += BIGNUMBER * (v2[0] - data.sol[0][0]) * phi(in, 0) * weight;
                ef(nstate * in + 1, 0) += BIGNUMBER * (v2[1] - data.sol[0][1]) * phi(in, 0) * weight;
                
                for (jn = 0; jn < phr; jn++) {
                    ek(nstate * in + 0, nstate * jn + 0) += BIGNUMBER * phi(in, 0) * phi(jn, 0) * weight;
                    ek(nstate * in + 1, nstate * jn + 1) += BIGNUMBER * phi(in, 0) * phi(jn, 0) * weight;
                    
                }//jn
            }//in
            break;
            
        case 1: // Neumann condition
            for (in = 0; in < phi.Rows(); in++) {
                ef(nstate * in + 0, 0) += v2[0] * phi(in, 0) * weight;
                ef(nstate * in + 1, 0) += v2[1] * phi(in, 0) * weight;
            }
            break;
            
        case 2: // Mixed condition
        {
            TPZFNMatrix<2, STATE> res(2, 1, 0.);
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    res(i, 0) += bc.Val1()(i, j) * data.sol[0][j];
                }
            }
            
            for (in = 0; in < phi.Rows(); in++) {
                ef(nstate * in + 0, 0) += (v2[0] - res(0, 0)) * phi(in, 0) * weight;
                ef(nstate * in + 1, 0) += (v2[1] - res(1, 0)) * phi(in, 0) * weight;
                for (jn = 0; jn < phi.Rows(); jn++) {
                    for (idf = 0; idf < 2; idf++) {
                        for (jdf = 0; jdf < 2; jdf++) {
                            ek(nstate * in + idf, nstate * jn + jdf) += bc.Val1()(idf, jdf) * phi(in, 0) * phi(jn, 0) * weight;
                            //BUG FALTA COLOCAR VAL2
                            //DebugStop();
                        }
                    }
                }
            }//in
        }
            break;
            
        case 3: // Directional Null Dirichlet - displacement is set to null in the non-null vector component direction
            for (in = 0; in < phr; in++) {
                ef(nstate * in + 0, 0) += BIGNUMBER * (0. - data.sol[0][0]) * v2[0] * phi(in, 0) * weight;
                ef(nstate * in + 1, 0) += BIGNUMBER * (0. - data.sol[0][1]) * v2[1] * phi(in, 0) * weight;
                for (jn = 0; jn < phr; jn++) {
                    ek(nstate * in + 0, nstate * jn + 0) += BIGNUMBER * phi(in, 0) * phi(jn, 0) * weight * v2[0];
                    ek(nstate * in + 1, nstate * jn + 1) += BIGNUMBER * phi(in, 0) * phi(jn, 0) * weight * v2[1];
                }//jn
            }//in
            break;
            
        case 4: // stressField Neumann condition
            v2[0] = v1(0, 0) * data.normal[0] + v1(0, 1) * data.normal[1];
            v2[1] = v1(1, 0) * data.normal[0] + v1(1, 1) * data.normal[1];
            // The normal vector points towards the neighbor. The negative sign is there to
            // reflect the outward normal vector.
            for (in = 0; in < phi.Rows(); in++) {
                ef(nstate * in + 0, 0) += v2[0] * phi(in, 0) * weight;
                ef(nstate * in + 1, 0) += v2[1] * phi(in, 0) * weight;
            }
            break;
        case 5://PRESSAO DEVE SER POSTA NA POSICAO 0 DO VETOR v2
        {
            TPZFNMatrix<2, STATE> res(2, 1, 0.);
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    res(i, 0) += data.normal[i] * bc.Val1()(i, j) * data.sol[0][j] * data.normal[j];
                }
            }
            for (in = 0; in < phi.Rows(); in++) {
                ef(nstate * in + 0, 0) += (v2[0] * data.normal[0] - res(0, 0)) * phi(in, 0) * weight;
                ef(nstate * in + 1, 0) += (v2[0] * data.normal[1] - res(1, 0)) * phi(in, 0) * weight;
                for (jn = 0; jn < phi.Rows(); jn++) {
                    for (idf = 0; idf < 2; idf++) {
                        for (jdf = 0; jdf < 2; jdf++) {
                            ek(nstate * in + idf, nstate * jn + jdf) += bc.Val1()(idf, jdf) * data.normal[idf] * data.normal[jdf] * phi(in, 0) * phi(jn, 0) * weight;
                            // BUG FALTA COLOCAR VAL2
                            // DebugStop();
                        }
                    }
                }
                
            }
        }
            break;
            
        case 6://PRESSAO DEVE SER POSTA NA POSICAO 0 DO VETOR v2
        {
            TPZFNMatrix<2, STATE> res(2, 1, 0.);
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    res(i, 0) += bc.Val1()(i, j) * data.sol[0][j];
                }
            }
            for (in = 0; in < phi.Rows(); in++) {
                ef(nstate * in + 0, 0) += (v2[0] * data.normal[0] - res(0, 0)) * phi(in, 0) * weight;
                ef(nstate * in + 1, 0) += (v2[0] * data.normal[1] - res(1, 0)) * phi(in, 0) * weight;
                for (jn = 0; jn < phi.Rows(); jn++) {
                    for (idf = 0; idf < 2; idf++) {
                        for (jdf = 0; jdf < 2; jdf++) {
                            ek(nstate * in + idf, nstate * jn + jdf) += bc.Val1()(idf, jdf) * phi(in, 0) * phi(jn, 0) * weight;
                            // BUG FALTA COLOCAR VAL2
                            // DebugStop();
                        }
                    }
                }
                
            }
            
        }
            break;
            
        case 7 : {
            
            REAL v_null[2];
            v_null[0] = bc.Val1()(0, 0);
            v_null[1] = bc.Val1()(1, 1);
            
            for (in = 0; in < phr; in++) {
                ef(nstate * in + 0, 0) += BIGNUMBER * (v2[0] - data.sol[0][0]) * v_null[0] * phi(in, 0) * weight;
                ef(nstate * in + 1, 0) += BIGNUMBER * (v2[1] - data.sol[0][1]) * v_null[1] * phi(in, 0) * weight;
                for (jn = 0; jn < phr; jn++) {
                    ek(nstate * in + 0, nstate * jn + 0) += BIGNUMBER * phi(in, 0) * phi(jn, 0) * weight * v_null[0];
                    ek(nstate * in + 1, nstate * jn + 1) += BIGNUMBER * phi(in, 0) * phi(jn, 0) * weight * v_null[1];
                }//jn
            }//in
            break;
            
        }
            
        default:
#ifdef LOG4CXX
            if (elastoplastic2dLogger->isDebugEnabled()) {
                std::stringstream sout;
                sout << "<<< TPZMatElastoPlastic2D<T,TMEM>::ContributeBC *** WRONG BOUNDARY CONDITION TYPE = " << bc.Type();
                LOGPZ_ERROR(elastoplastic2dLogger, sout.str().c_str());
            }
#endif
            PZError << "TPZMatElastoPlastic2D::ContributeBC error - Wrong boundary condition type" << std::endl;
    }
}

template <class T, class TMEM>
void TPZMatElastoPlastic2D<T, TMEM>::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc) {
    TPZFMatrix<REAL> &phi = data.phi;
    const REAL BIGNUMBER = TPZMaterial::gBigNumber;
    int nstate = NStateVariables();
    
    const int phr = phi.Rows();
    int in;
    REAL v2[2];
    v2[0] = bc.Val2()(0, 0);
    v2[1] = bc.Val2()(1, 0);
    
    TPZFMatrix<REAL> &v1 = bc.Val1();
    switch (bc.Type()) {
        case 0: // Dirichlet condition
            for (in = 0; in < phr; in++) {
                ef(nstate * in + 0, 0) += BIGNUMBER * (v2[0] - data.sol[0][0]) * phi(in, 0) * weight;
                ef(nstate * in + 1, 0) += BIGNUMBER * (v2[1] - data.sol[0][1]) * phi(in, 0) * weight;
            }//in
            break;
            
        case 1: // Neumann condition
            for (in = 0; in < phi.Rows(); in++) {
                ef(nstate * in + 0, 0) += v2[0] * phi(in, 0) * weight;
                ef(nstate * in + 1, 0) += v2[1] * phi(in, 0) * weight;
            }
            break;
            
        case 2: // Mixed condition
        {
            TPZFNMatrix<2, STATE> res(2, 1, 0.);
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    res(i, 0) += bc.Val1()(i, j) * data.sol[0][j];
                }
            }
            
            for (in = 0; in < phi.Rows(); in++) {
                ef(nstate * in + 0, 0) += (v2[0] - res(0, 0)) * phi(in, 0) * weight;
                ef(nstate * in + 1, 0) += (v2[1] - res(1, 0)) * phi(in, 0) * weight;
            }//in
        }
            break;
            
        case 3: // Directional Null Dirichlet - displacement is set to null in the non-null vector component direction
            for (in = 0; in < phr; in++) {
                ef(nstate * in + 0, 0) += BIGNUMBER * (0. - data.sol[0][0]) * v2[0] * phi(in, 0) * weight;
                ef(nstate * in + 1, 0) += BIGNUMBER * (0. - data.sol[0][1]) * v2[1] * phi(in, 0) * weight;
            }//in
            break;
            
        case 4: // stressField Neumann condition
            v2[0] = v1(0, 0) * data.normal[0] + v1(0, 1) * data.normal[1];
            v2[1] = v1(1, 0) * data.normal[0] + v1(1, 1) * data.normal[1];
            // The normal vector points towards the neighbour. The negative sign is there to
            // reflect the outward normal vector.
            for (in = 0; in < phi.Rows(); in++) {
                ef(nstate * in + 0, 0) += v2[0] * phi(in, 0) * weight;
                ef(nstate * in + 1, 0) += v2[1] * phi(in, 0) * weight;
            }
            break;
            
        case 5://PRESSAO DEVE SER POSTA NA POSICAO 0 DO VETOR v2
        {
            TPZFNMatrix<2, STATE> res(2, 1, 0.);
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    res(i, 0) += data.normal[i] * bc.Val1()(i, j) * data.sol[0][j] * data.normal[j];
                }
            }
            for (in = 0; in < phi.Rows(); in++) {
                ef(nstate * in + 0, 0) += (v2[0] * data.normal[0] - res(0, 0)) * phi(in, 0) * weight;
                ef(nstate * in + 1, 0) += (v2[0] * data.normal[1] - res(1, 0)) * phi(in, 0) * weight;
            }
        }
            break;
            
        case 6://PRESSAO DEVE SER POSTA NA POSICAO 0 DO VETOR v2
        {
            TPZFNMatrix<2, STATE> res(2, 1, 0.);
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    res(i, 0) += bc.Val1()(i, j) * data.sol[0][j];
                }
            }
            for (in = 0; in < phi.Rows(); in++) {
                ef(nstate * in + 0, 0) += (v2[0] * data.normal[0] - res(0, 0)) * phi(in, 0) * weight;
                ef(nstate * in + 1, 0) += (v2[0] * data.normal[1] - res(1, 0)) * phi(in, 0) * weight;
            }
            
        }
            break;
            
        case 7 : {
            
            REAL v_null[2];
            v_null[0] = bc.Val1()(0, 0);
            v_null[1] = bc.Val1()(1, 1);
            
            for (in = 0; in < phr; in++) {
                ef(nstate * in + 0, 0) += BIGNUMBER * (v2[0] - data.sol[0][0]) * v_null[0] * phi(in, 0) * weight;
                ef(nstate * in + 1, 0) += BIGNUMBER * (v2[1] - data.sol[0][1]) * v_null[1] * phi(in, 0) * weight;
            }//in
            break;
            
        }
            
        default:
#ifdef LOG4CXX
            if (elastoplastic2dLogger->isDebugEnabled()) {
                std::stringstream sout;
                sout << "<<< TPZMatElastoPlastic2D<T,TMEM>::ContributeBC *** WRONG BOUNDARY CONDITION TYPE = " << bc.Type();
                LOGPZ_ERROR(elastoplastic2dLogger, sout.str().c_str());
            }
#endif
            PZError << "TPZMatElastoPlastic2D::ContributeBC error - Wrong boundary condition type" << std::endl;
    }
}



template <class T, class TMEM>
void TPZMatElastoPlastic2D<T,TMEM>::Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout)
{
    
    TPZMaterialData datalocal(data);
    datalocal.sol[0].Resize(3,0.);
    datalocal.dsol[0].Resize(3,3);
    datalocal.dsol[0](2,0) = 0.;
    datalocal.dsol[0](2,1) = 0.;
    datalocal.dsol[0](2,2) = 0.;
    datalocal.dsol[0](0,2) = 0.;
    datalocal.dsol[0](1,2) = 0.;
    TPZMatElastoPlastic<T,TMEM>::Solution(datalocal,var,Solout);
    
}

template <class T, class TMEM>
void TPZMatElastoPlastic2D<T, TMEM>::ComputeDeltaStrainVector(TPZMaterialData & data, TPZFMatrix<REAL> &DeltaStrain) {
    TPZFNMatrix<9> DSolXYZ(3, 3, 0.);
    data.axes.Multiply(data.dsol[0], DSolXYZ, 1/*transpose*/);
    if (DeltaStrain.Rows() != 6) {
        DebugStop();
    }
    //  DeltaStrain.Redim(3,1);
    DeltaStrain(_XX_, 0) = DSolXYZ(0, 0);
    DeltaStrain(_YY_, 0) = DSolXYZ(1, 1);
    DeltaStrain(_XY_, 0) = 0.5 * (DSolXYZ(1, 0) + DSolXYZ(0, 1));
    DeltaStrain(_XZ_, 0) = 0.;
    DeltaStrain(_YZ_, 0) = 0.;
    DeltaStrain(_ZZ_, 0) = 0.;
}


template <class T, class TMEM>
TPZMaterial * TPZMatElastoPlastic2D<T,TMEM>::NewMaterial()
{
    return new TPZMatElastoPlastic2D<T,TMEM>(*this);
}

#include "TPZSandlerExtended.h"
#include "TPZPlasticStepPV.h"
#include "TPZYCMohrCoulombPV.h"

template <class T, class TMEM>
std::string TPZMatElastoPlastic2D<T,TMEM>::Name()
{
    return "TPZMatElastoPlastic<T,TMEM>";
}

template <class T, class TMEM>
void TPZMatElastoPlastic2D<T,TMEM>::Write(TPZStream &buf, int withclassid) const{
    TPZMatElastoPlastic<T,TMEM>::Write(buf,withclassid);
    int classid = ClassId();
    buf.Write(&classid);
}

template <class T, class TMEM>
void TPZMatElastoPlastic2D<T,TMEM>::Read(TPZStream &buf, void *context)
{
    TPZMatElastoPlastic<T,TMEM>::Read(buf,context);
    int classid;
    buf.Read(&classid);
    if (classid != ClassId()) {
        DebugStop();
    }
}


template <class T, class TMEM>
void TPZMatElastoPlastic2D<T,TMEM>::Print(std::ostream &out, const int memory)
{
    out << __PRETTY_FUNCTION__ << std::endl;
    TPZMatElastoPlastic<T,TMEM>::Print(out,memory);
}

template <class T, class TMEM>
void TPZMatElastoPlastic2D<T,TMEM>::Print(std::ostream &out)
{
    out << __PRETTY_FUNCTION__ << std::endl;
    out << "Plane strain " << fPlaneStrain << std::endl;
    TPZMatElastoPlastic<T,TMEM>::Print(out);
}
