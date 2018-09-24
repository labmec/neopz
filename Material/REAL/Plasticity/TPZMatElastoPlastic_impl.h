//
//  TPZMatElastoPlastic_impl.h
//  pz
//
//  Created by Omar Durán on 9/8/18.
//

#include "TPZMatElastoPlastic.h"
#include "pzbndcond.h"
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr elastoplasticLogger(Logger::getLogger("pz.material.pzElastoPlastic"));
static LoggerPtr updatelogger(Logger::getLogger("pz.material.pzElastoPlastic.update"));
static LoggerPtr ceckconvlogger(Logger::getLogger("checkconvmaterial"));
#endif


template <class T, class TMEM>
TPZMatElastoPlastic<T,TMEM>::TPZMatElastoPlastic() : TPZMatWithMem<TMEM>(), fForce(), fRhoB(0), fPostProcessDirection(), fTol(1.e-6)
{
    fForce.Resize(3,0);
    fForce[1] = -9.8; // proper gravity acceleration in m/s^2
    fPostProcessDirection.Resize(3,0);
    fPostProcessDirection[0] = 1.;
    
#ifdef LOG4CXX
    if(elastoplasticLogger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << ">>> TPZMatElastoPlastic<T,TMEM>() constructor called ***";
        LOGPZ_DEBUG(elastoplasticLogger,sout.str().c_str());
    }
#endif
    
}

template <class T, class TMEM>
TPZMatElastoPlastic<T,TMEM>::TPZMatElastoPlastic(int id) : TPZMatWithMem<TMEM>(id), fForce(), fRhoB(0), fPostProcessDirection(), fTol(1.e-6)
{
    fForce.Resize(3,0);
    fForce[1] = -9.8; // proper gravity acceleration in m/s^2 -> 1=y 0=x 2=z
    fPostProcessDirection.Resize(3,0);
    fPostProcessDirection[0] = 1.;
    
    TPZPlasticState<STATE> def;
    
    
#ifdef LOG4CXX
    if (elastoplasticLogger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << ">>> TPZMatElastoPlastic<T,TMEM>(int id) constructor called with id = " << id << " ***";
        LOGPZ_DEBUG(elastoplasticLogger,sout.str().c_str());
    }
#endif
    
}

template <class T, class TMEM>
TPZMatElastoPlastic<T,TMEM>::TPZMatElastoPlastic(const TPZMatElastoPlastic &mat) : TPZMatWithMem<TMEM>(mat),
fForce(mat.fForce), fRhoB(mat.fRhoB), fPostProcessDirection(mat.fPostProcessDirection),
fPlasticity(mat.fPlasticity), fTol(mat.fTol)
{
#ifdef LOG4CXX
    if(elastoplasticLogger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << ">>> TPZMatElastoPlastic<T,TMEM>() copy constructor called ***";
        LOGPZ_DEBUG(elastoplasticLogger,sout.str().c_str());
    }
#endif
}


template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::SetPlasticity(T & plasticity)
{
#ifdef LOG4CXX
    if(elastoplasticLogger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << ">>> TPZMatElastoPlastic<T,TMEM>::SetUpPlasticity ***";
        sout << "\n with plasticity argument:\n";
        plasticity.Print(sout);
        LOGPZ_DEBUG(elastoplasticLogger,sout.str().c_str());
    }
#endif
    
    fPlasticity = plasticity;
    
    //fPlasticity.SetTensionSign(1);
    
    T plastloc(fPlasticity);
    
    TMEM memory;
    
    memory.fPlasticState = plastloc.GetState();
    
    plastloc.ApplyStrainComputeSigma(memory.fPlasticState.fEpsT, memory.fSigma);
    
    this->SetDefaultMem(memory);
    
#ifdef LOG4CXX
    if(elastoplasticLogger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "<< TPZMatElastoPlastic<T,TMEM>::SetUpPlasticity ***";
        sout << "\n with computed stresses:\n";
        sout << memory.fSigma;
        LOGPZ_DEBUG(elastoplasticLogger,sout.str().c_str());
    }
#endif
    
}


template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::SetBulkDensity(REAL & RhoB)
{
    fRhoB = RhoB;
}

template <class T, class TMEM>
TPZMatElastoPlastic<T,TMEM>::~TPZMatElastoPlastic()
{
    
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::Print(std::ostream &out, const int memory)
{
    out << this->Name();
    out << "\n with template argurment T = " << fPlasticity.Name();
    out << "\n Base material Data:\n";
    TPZMatWithMem<TMEM>::PrintMem(out, memory);
    out << "\n Localy defined members:";
    out << "\n Body Forces: " << fForce;
    out << "\n Post process direction: " << fPostProcessDirection;
    out << "\n Tolerance for internal post processing iterations: " << fTol;
    out << "\n Internal plasticity <T> member:\n";
    fPlasticity.Print(out);
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::Print(std::ostream &out)
{
    out << __PRETTY_FUNCTION__ << std::endl;
    out << this->Name();
    TPZMatWithMem<TMEM>::Print(out);
    out << "\nBody Forces: " << fForce;
    out << "\nfRhoB = " << fRhoB;
    out << "\nPost process direction: " << fPostProcessDirection;
    out << "\nTolerance for internal post processing iterations: " << fTol;
    out << "\nInternal plasticity <T> member:\n";
    fPlasticity.Print(out);
}

template <class T, class TMEM>
int TPZMatElastoPlastic<T,TMEM>::VariableIndex(const std::string &name)
{
    if(!strcmp("Displacement",             name.c_str()))  return TPZMatElastoPlastic<T,TMEM>::EDisplacement;
    if(!strcmp("DisplacementX",             name.c_str()))  return TPZMatElastoPlastic<T,TMEM>::EDisplacementX;
    if(!strcmp("DisplacementY",             name.c_str()))  return TPZMatElastoPlastic<T,TMEM>::EDisplacementY;
    if(!strcmp("DisplacementZ",             name.c_str()))  return TPZMatElastoPlastic<T,TMEM>::EDisplacementZ;
    if(!strcmp("NormalStress",             name.c_str()))  return TPZMatElastoPlastic<T,TMEM>::ENormalStress;
    if(!strcmp("ShearStress",              name.c_str()))  return TPZMatElastoPlastic<T,TMEM>::EShearStress;
    if(!strcmp("NormalStrain",             name.c_str()))  return TPZMatElastoPlastic<T,TMEM>::ENormalStrain;
    if(!strcmp("ShearStrain",              name.c_str()))  return TPZMatElastoPlastic<T,TMEM>::EShearStrain;
    if(!strcmp("PrincipalStress",          name.c_str()))  return TPZMatElastoPlastic<T,TMEM>::EPrincipalStress;
    if(!strcmp("PrincipalStrain",          name.c_str()))  return TPZMatElastoPlastic<T,TMEM>::EPrincipalStrain;
    
    if(!strcmp("I1Stress",         name.c_str()))  return TPZMatElastoPlastic<T,TMEM>::EI1Stress;
    if(!strcmp("J2Stress",         name.c_str()))  return TPZMatElastoPlastic<T,TMEM>::EJ2Stress;
    if(!strcmp("VolElasticStrain",         name.c_str()))  return TPZMatElastoPlastic<T,TMEM>::EVolElasticStrain;
    
    if(!strcmp("VolPlasticStrain",         name.c_str()))  return TPZMatElastoPlastic<T,TMEM>::EVolPlasticStrain;
    if(!strcmp("VolTotalStrain",           name.c_str()))  return TPZMatElastoPlastic<T,TMEM>::EVolTotalStrain;
    if(!strcmp("Alpha",                    name.c_str()))  return TPZMatElastoPlastic<T,TMEM>::EAlpha;
    if(!strcmp("PlasticSteps",             name.c_str()))  return TPZMatElastoPlastic<T,TMEM>::EPlasticSteps;
    if(!strcmp("PlasticSqJ2",             name.c_str()))  return TPZMatElastoPlastic<T,TMEM>::EPlasticSqJ2;
    if(!strcmp("YieldSurface",             name.c_str()))  return TPZMatElastoPlastic<T,TMEM>::EYield;
    if(!strcmp("TotalPlasticStrain",     name.c_str()))  return TPZMatElastoPlastic<T,TMEM>::ENormalPlasticStrain;
    if(!strcmp("EMisesStress",     name.c_str()))  return TPZMatElastoPlastic<T,TMEM>::EMisesStress;
    if(!strcmp("DisplacementMem",     name.c_str()))  return TPZMatElastoPlastic<T,TMEM>::EDisplacementMem;
    if(!strcmp("DisplacementMemX",     name.c_str()))  return TPZMatElastoPlastic<T,TMEM>::EDisplacementMemX;
    if(!strcmp("DisplacementMemY",     name.c_str()))  return TPZMatElastoPlastic<T,TMEM>::EDisplacementMemY;
    if(!strcmp("DisplacementMemZ",     name.c_str()))  return TPZMatElastoPlastic<T,TMEM>::EDisplacementMemZ;
    if(!strcmp("XStress",     name.c_str()))  return TPZMatElastoPlastic<T,TMEM>::EXStress;
    if(!strcmp("YStress",     name.c_str()))  return TPZMatElastoPlastic<T,TMEM>::EYStress;
    if(!strcmp("ZStress",     name.c_str()))  return TPZMatElastoPlastic<T,TMEM>::EZStress;
    if(!strcmp("PlasticSqJ2El",     name.c_str()))  return 100;
    return TPZMatWithMem<TMEM>::VariableIndex(name);
    PZError << "TPZMatElastoPlastic::VariableIndex Error\n";
    return -1;
}

template <class T, class TMEM>
int TPZMatElastoPlastic<T,TMEM>::NSolutionVariables(int var)
{
    if(var == TPZMatElastoPlastic<T,TMEM>::EDisplacement)              return 3;
    if(var == TPZMatElastoPlastic<T,TMEM>::EDisplacementX)              return 1;
    if(var == TPZMatElastoPlastic<T,TMEM>::EDisplacementY)              return 1;
    if(var == TPZMatElastoPlastic<T,TMEM>::EDisplacementZ)              return 1;
    if(var == TPZMatElastoPlastic<T,TMEM>::EDisplacementMem)           return 3;
    if(var == TPZMatElastoPlastic<T,TMEM>::EDisplacementMemX)           return 1;
    if(var == TPZMatElastoPlastic<T,TMEM>::EDisplacementMemY)           return 1;
    if(var == TPZMatElastoPlastic<T,TMEM>::EDisplacementMemZ)           return 1;
    if(var == TPZMatElastoPlastic<T,TMEM>::EPrincipalStress)           return 3;
    if(var == TPZMatElastoPlastic<T,TMEM>::ENormalStress)              return 3;
    if(var == TPZMatElastoPlastic<T,TMEM>::EShearStress)               return 3;
    if(var == TPZMatElastoPlastic<T,TMEM>::ENormalStrain)              return 3;
    if(var == TPZMatElastoPlastic<T,TMEM>::EShearStrain)               return 3;
    if(var == TPZMatElastoPlastic<T,TMEM>::ENormalStrain)              return 3;
    if(var == TPZMatElastoPlastic<T,TMEM>::ENormalPlasticStrain)       return 3;
    if(var == TPZMatElastoPlastic<T,TMEM>::EPrincipalStrain)           return 3;
    
    if(var == TPZMatElastoPlastic<T,TMEM>::EI1Stress)                  return 1;
    if(var == TPZMatElastoPlastic<T,TMEM>::EJ2Stress)                  return 1;
    if(var == TPZMatElastoPlastic<T,TMEM>::EVolElasticStrain)          return 1;
    if(var == TPZMatElastoPlastic<T,TMEM>::EVolPlasticStrain)          return 1;
    if(var == TPZMatElastoPlastic<T,TMEM>::EVolTotalStrain)            return 1;
    
    if(var == TPZMatElastoPlastic<T,TMEM>::EAlpha)                     return 1;
    if(var == TPZMatElastoPlastic<T,TMEM>::EPlasticSteps)              return 1;
    if(var == TPZMatElastoPlastic<T,TMEM>::EPlasticSqJ2)               return 1;
    if(var == TPZMatElastoPlastic<T,TMEM>::EYield)                     return T::fNYields::NYield;//Numero de funcoes falha
    if(var == TPZMatElastoPlastic<T,TMEM>::EMisesStress)               return 1;
    if(var == TPZMatElastoPlastic<T,TMEM>::EXStress)                   return 1;
    if(var == TPZMatElastoPlastic<T,TMEM>::EYStress)                   return 1;
    if(var == TPZMatElastoPlastic<T,TMEM>::EZStress)                   return 1;
    
    if(var == 100) return 1;
    return TPZMatWithMem<TMEM>::NSolutionVariables(var);
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::ApplyDirection(TPZFMatrix<REAL> &vectorTensor, TPZVec<REAL> &Out)
{
    Out.Resize(3);
    TPZVec<REAL> &Dir = this->fPostProcessDirection;
    Out[0] = Dir[0] * vectorTensor(_XX_,0) + Dir[1] * vectorTensor(_XY_,0) + Dir[2] * vectorTensor(_XZ_,0);
    Out[1] = Dir[0] * vectorTensor(_XY_,0) + Dir[1] * vectorTensor(_YY_,0) + Dir[2] * vectorTensor(_YZ_,0);
    Out[2] = Dir[0] * vectorTensor(_XZ_,0) + Dir[1] * vectorTensor(_YZ_,0) + Dir[2] * vectorTensor(_ZZ_,0);
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T, TMEM>::Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout) {
    
    if (var == TPZMatElastoPlastic<T, TMEM>::EDisplacement)
    {
        for (int i = 0; i < 3; ++i) {
            Solout[i] = data.sol[0][i];
        }//for
    }//EDisplacement from DoF
    
    int intPt = data.intGlobPtIndex;
    if (intPt == -1) {
        return;
    }
    
    TMEM &Memory = this->MemItem(intPt);
    T plasticloc(fPlasticity);
    plasticloc.SetState(Memory.fPlasticState);
    switch (var) {
        case TPZMatElastoPlastic<T, TMEM>::EDisplacementX:
        {
            Solout[0] = data.sol[0][0];
        }
            break;
        case TPZMatElastoPlastic<T, TMEM>::EDisplacementY:
        {
            Solout[0] = data.sol[0][1];
        }
            break;
        case TPZMatElastoPlastic<T, TMEM>::EDisplacementZ:
        {
            Solout[0] = data.sol[0][2];
        }
            break;
        case EDisplacementMem:
        {
            for (int i = 0; i < 3; i++) {
                Solout[i] = Memory.fDisplacement[i];
            }
        }
            break;
        case TPZMatElastoPlastic<T, TMEM>::EDisplacementMemX:
        {
            Solout[0] = Memory.fDisplacement[0];
        }
            break;
        case TPZMatElastoPlastic<T, TMEM>::EDisplacementMemY:
        {
            Solout[0] = Memory.fDisplacement[1];
        }
            break;
        case TPZMatElastoPlastic<T, TMEM>::EDisplacementMemZ:
        {
            Solout[0] = Memory.fDisplacement[2];
        }
            break;
        case TPZMatElastoPlastic<T, TMEM>::ENormalStrain:
        {
            TPZTensor<REAL> & totalStrain = Memory.fPlasticState.fEpsT;
            Solout[0] = totalStrain.XX();
            Solout[1] = totalStrain.YY();
            Solout[2] = totalStrain.ZZ();
        }//ENormalStrain
            break;
        case TPZMatElastoPlastic<T, TMEM>::EShearStrain:
        {
            TPZTensor<REAL> & totalStrain = Memory.fPlasticState.fEpsT;
            Solout[0] = totalStrain.XY();
            Solout[1] = totalStrain.XZ();
            Solout[2] = totalStrain.YZ();
        }//EShearStrain
            break;
        case TPZMatElastoPlastic<T, TMEM>::ENormalStress:
        {
            TPZTensor<REAL> & Sigma = Memory.fSigma;
            Solout[0] = Sigma.XX();
            Solout[1] = Sigma.YY();
            Solout[2] = Sigma.ZZ();
        }//ENormalStress
            break;
        case TPZMatElastoPlastic<T, TMEM>::EYStress:
        {
            TPZTensor<REAL> & Sigma = Memory.fSigma;
            Solout[0] = Sigma.YY();
        }
            break;
        case TPZMatElastoPlastic<T, TMEM>::EZStress:
        {
            TPZTensor<REAL> & Sigma = Memory.fSigma;
            Solout[0] = Sigma.ZZ();
        }
            break;
        case TPZMatElastoPlastic<T, TMEM>::EXStress:
        {
            TPZTensor<REAL> & Sigma = Memory.fSigma;
            Solout[0] = Sigma.XX();
        }
            break;
        case TPZMatElastoPlastic<T, TMEM>::EShearStress:
        {
            TPZTensor<REAL> & Sigma = Memory.fSigma;
            Solout[0] = Sigma.XY();
            Solout[1] = Sigma.XZ();
            Solout[2] = Sigma.YZ();
        }//EShearStress
            break;
        case TPZMatElastoPlastic<T, TMEM>::EPrincipalStress:
        {
            TPZTensor<REAL> & Sigma = Memory.fSigma;
            TPZTensor<REAL>::TPZDecomposed eigensystem;
            Sigma.EigenSystem(eigensystem);
            for (int i = 0; i < 3; i++)Solout[i] = eigensystem.fEigenvalues[i];
        }//EPrincipalStress - makes sense only if the evaluated point refers to an identified integration point
            break;
        case TPZMatElastoPlastic<T, TMEM>::EPrincipalStrain:
        {
            TPZTensor<REAL> & eps = Memory.fPlasticState.fEpsT;
            TPZTensor<REAL>::TPZDecomposed eigensystem;
            eps.EigenSystem(eigensystem);
            for (int i = 0; i < 3; i++)Solout[i] = eigensystem.fEigenvalues[i];
        }//EPrincipalStrain
            break;
        case TPZMatElastoPlastic<T, TMEM>::EI1Stress:
        {
            TPZTensor<REAL> Sigma = Memory.fSigma;
            Solout[0] = Sigma.I1();
        }//EI1Stress - makes sense only if the evaluated point refers to an identified integration point
            break;
        case TPZMatElastoPlastic<T, TMEM>::EJ2Stress:
        {
            TPZTensor<REAL> Sigma = Memory.fSigma;
            Solout[0] = Sigma.J2();
        }//EJ2Stress - makes sense only if the evaluated point refers to an identified integration point
            break;
        case TPZMatElastoPlastic<T, TMEM>::EVolPlasticStrain:
        {
            TPZTensor<REAL> & plasticStrain = Memory.fPlasticState.fEpsP;
            Solout[0] = plasticStrain.I1();
        }//EVolPlasticStrain - makes sense only if the evaluated point refers to an identified integration point
            break;
        case TPZMatElastoPlastic<T, TMEM>::EVolElasticStrain:
        {
            TPZTensor<REAL> & plasticStrain = Memory.fPlasticState.fEpsP;
            TPZTensor<REAL> & totalStrain = Memory.fPlasticState.fEpsT;
            Solout[0] = totalStrain.I1() - plasticStrain.I1();
        }//EVolElasticStrain - makes sense only if the evaluated point refers to an identified integration point
            break;
        case TPZMatElastoPlastic<T, TMEM>::EVolTotalStrain:
        {
            TPZTensor<REAL> & totalStrain = Memory.fPlasticState.fEpsT;
            Solout[0] = totalStrain.I1();
        }//EVolElasticStrain - makes sense only if the evaluated point refers to an identified integration point
            break;
        case TPZMatElastoPlastic<T, TMEM>::EPlasticSqJ2:
        {
            TPZTensor<REAL> & plasticStrain = Memory.fPlasticState.fEpsP;
            Solout[0] = sqrt(plasticStrain.J2());
        }//EVolTEPStrain - makes sense only if the evaluated point refers to an identified integration point
            break;
        case TPZMatElastoPlastic<T, TMEM>::EAlpha:
        {
            Solout[0] = Memory.fPlasticState.fAlpha;
        }//EAlpha - makes sense only if the evaluated point refers to an identified integration point
            break;
        case TPZMatElastoPlastic<T, TMEM>::EPlasticSteps:
        {
            Solout[0] = Memory.fPlasticSteps;
        }//EVolPlasticSteps - makes sense only if the evaluated point refers to an identified integration point
            break;
        case TPZMatElastoPlastic<T, TMEM>::EYield:
        {
            TPZTensor<REAL> & EpsT = Memory.fPlasticState.fEpsT;
            TPZTensor<STATE> epsElastic(EpsT);
            epsElastic -= Memory.fPlasticState.fEpsP;
            plasticloc.Phi(epsElastic, Solout);
        }//EVolPlasticSteps - makes sense only if the evaluated point refers to an identified integration point
            break;
        case TPZMatElastoPlastic<T, TMEM>::ENormalPlasticStrain:
        {
            TPZTensor<REAL> & plasticStrain = Memory.fPlasticState.fEpsP;
            Solout[0] = plasticStrain.XX();
            Solout[1] = plasticStrain.YY();
            Solout[2] = plasticStrain.ZZ();
        }
            break;
        case TPZMatElastoPlastic<T, TMEM>::EMisesStress:
        {
            TPZTensor<REAL> Sigma = Memory.fSigma;
            REAL J2 = Sigma.J2();
            REAL temp = sqrt(3. * J2);
            Solout[0] = temp;
        }//VonMisesStress
            break;
        default:
        {
            TPZMatWithMem<TMEM>::Solution(data, var, Solout);
        }
    }
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef)
{
    
#ifdef LOG4CXX
    if(elastoplasticLogger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << ">>> TPZMatElastoPlastic<T,TMEM>::Contribute ***";
        sout << "\nIntegration Point index = " << data.intGlobPtIndex;
        LOGPZ_DEBUG(elastoplasticLogger,sout.str().c_str());
    }
#endif
    
    TPZFMatrix<REAL> &dphi = data.dphix, dphiXYZ;
    TPZFMatrix<REAL> &phi  = data.phi;
    TPZFMatrix<REAL> &axes = data.axes, axesT;
    TPZManVector<REAL,3> &x = data.x;
    
    // rotating the shape functions to the XYZ coordinates
    axes.Transpose(&axesT);
    axesT.Multiply(dphi,dphiXYZ);
    /*
     cout << "\n phi(" << data.intPtIndex << ") =";
     for(int i = 0; i < data.phi.Rows(); i++)cout << " " << data.phi(i,0);
     cout << endl << dphiXYZ;
     cout << endl << axes;
     */
    const int phr = phi.Rows();
    if(this->fForcingFunction)
        this->fForcingFunction->Execute(x,this->fForce);
    
    //this matrix will store {{dvdx*dudx, dvdx*dudy, dvdx*dudz},
    //{dvdy*dudx, dvdy*dudy, dvdy*dudz},
    //{dvdz*dudx, dvdz*dudy, dvdz*dudz}}
    TPZFNMatrix<9>  Deriv(3,3);
    TPZFNMatrix<36> Dep(6,6);
    TPZFNMatrix<6>  DeltaStrain(6,1);
    TPZFNMatrix<6>  Stress(6,1);//, StressN(6,1);
    
    this->ComputeDeltaStrainVector(data, DeltaStrain);
    this->ApplyDeltaStrainComputeDep(data, DeltaStrain, Stress, Dep);
    
    //int dim = Dimension();
    int nstate = NStateVariables();
    REAL val,val2,val3,val4,val5,val6,val7,val8,val9,val10;
    
    TPZVec<STATE> ForceLoc(this->fForce);
    if(this->fForcingFunction)
    {
        this->fForcingFunction->Execute(data.x,ForceLoc);
    }
    
    int in;
    for(in = 0; in < phr; in++) { //in: test function index
        
        // fForce represents the gravity acceleration
        //First equation: fb and fk
        val  = fRhoB * ForceLoc[0] * phi(in,0); // fb
        val -= Stress(_XX_,0) * dphiXYZ(0,in); // |
        val -= Stress(_XY_,0) * dphiXYZ(1,in); // fk
        val -= Stress(_XZ_,0) * dphiXYZ(2,in); // |
        ef(in*nstate+0,0) += weight * val;
        
        //Second equation: fb and fk
        val  = fRhoB * ForceLoc[1] * phi(in,0); // fb
        val -= Stress(_XY_,0) * dphiXYZ(0,in); // |
        val -= Stress(_YY_,0) * dphiXYZ(1,in); // fk
        val -= Stress(_YZ_,0) * dphiXYZ(2,in); // |
        ef(in*nstate+1,0) += weight * val;
        
        //third equation: fb and fk
        val  = fRhoB * ForceLoc[2] * phi(in,0); // fb
        val -= Stress(_XZ_,0) * dphiXYZ(0,in); // |
        val -= Stress(_YZ_,0) * dphiXYZ(1,in); // fk
        val -= Stress(_ZZ_,0) * dphiXYZ(2,in); // |
        ef(in*nstate+2,0) += weight * val;
        
        for( int jn = 0; jn < phr; jn++ ) {
            //jn: trial function index
            //this matrix will store
            //{{dvdx*dudx, dvdx*dudy, dvdx*dudz},
            //{dvdy*dudx, dvdy*dudy, dvdy*dudz},
            //{dvdz*dudx, dvdz*dudy, dvdz*dudz}}
            //Compute Deriv matrix
            for(int ud = 0; ud < 3; ud++){
                for(int vd = 0; vd < 3; vd++){
                    Deriv(vd,ud) = dphiXYZ(vd,in)*dphiXYZ(ud,jn);
                }//ud
            }//vd
            
            
            //#define _XX_ 0
            //#define _XY_ 1
            //#define _XZ_ 2
            //#define _YY_ 3
            //#define _YZ_ 4
            //#define _ZZ_ 5
            //First equation Dot[Sigma1, gradV1]
            val2  = 2. * Dep(_XX_,_XX_) * Deriv(0,0);//dvdx*dudx
            val2 +=      Dep(_XX_,_XY_) * Deriv(0,1);//dvdx*dudy
            val2 +=       Dep(_XX_,_XZ_) * Deriv(0,2);//dvdx*dudz
            val2 += 2. * Dep(_XY_,_XX_) * Deriv(1,0);//dvdy*dudx
            val2 +=      Dep(_XY_,_XY_) * Deriv(1,1);//dvdy*dudy
            val2 +=      Dep(_XY_,_XZ_) * Deriv(1,2);//dvdy*dudz
            val2 += 2. * Dep(_XZ_,_XX_) * Deriv(2,0);//dvdz*dudx
            val2 +=      Dep(_XZ_,_XY_) * Deriv(2,1);//dvdz*dudy
            val2 +=      Dep(_XZ_,_XZ_) * Deriv(2,2);//dvdz*dudz
            val2 *= 0.5;
            ek(in*nstate+0,jn*nstate+0) += weight * val2;
            
            val3  =      Dep(_XX_,_XY_) * Deriv(0,0);
            val3 += 2. * Dep(_XX_,_YY_) * Deriv(0,1);
            val3 +=      Dep(_XX_,_YZ_) * Deriv(0,2);
            val3 +=      Dep(_XY_,_XY_) * Deriv(1,0);
            val3 += 2. * Dep(_XY_,_YY_) * Deriv(1,1);
            val3 +=      Dep(_XY_,_YZ_) * Deriv(1,2);
            val3 +=      Dep(_XZ_,_XY_) * Deriv(2,0);
            val3 += 2. * Dep(_XZ_,_YY_) * Deriv(2,1);
            val3 +=      Dep(_XZ_,_YZ_) * Deriv(2,2);
            val3 *= 0.5;
            ek(in*nstate+0,jn*nstate+1) += weight * val3;
            
            val4  =      Dep(_XX_,_XZ_) * Deriv(0,0);
            val4 +=      Dep(_XX_,_YZ_) * Deriv(0,1);
            val4 += 2. * Dep(_XX_,_ZZ_) * Deriv(0,2);//
            val4 +=      Dep(_XY_,_XZ_) * Deriv(1,0);
            val4 +=      Dep(_XY_,_YZ_) * Deriv(1,1);
            val4 += 2. * Dep(_XY_,_ZZ_) * Deriv(1,2);//
            val4 +=      Dep(_XZ_,_XZ_) * Deriv(2,0);
            val4 +=      Dep(_XZ_,_YZ_) * Deriv(2,1);
            val4 += 2. * Dep(_XZ_,_ZZ_) * Deriv(2,2);
            val4 *= 0.5;
            ek(in*nstate+0,jn*nstate+2) += weight * val4;
            
            //Second equation Dot[Sigma2, gradV2]
            val5  = 2. * Dep(_XY_,_XX_) * Deriv(0,0);
            val5 +=      Dep(_XY_,_XY_) * Deriv(0,1);
            val5 +=      Dep(_XY_,_XZ_) * Deriv(0,2);
            val5 += 2. * Dep(_YY_,_XX_) * Deriv(1,0);
            val5 +=      Dep(_YY_,_XY_) * Deriv(1,1);
            val5 +=      Dep(_YY_,_XZ_) * Deriv(1,2);
            val5 += 2. * Dep(_YZ_,_XX_) * Deriv(2,0);
            val5 +=      Dep(_YZ_,_XY_) * Deriv(2,1);
            val5 +=      Dep(_YZ_,_XZ_) * Deriv(2,2);
            val5 *= 0.5;
            ek(in*nstate+1,jn*nstate+0) += weight * val5;
            
            val6  =      Dep(_XY_,_XY_) * Deriv(0,0);
            val6 += 2. * Dep(_XY_,_YY_) * Deriv(0,1);
            val6 +=      Dep(_XY_,_YZ_) * Deriv(0,2);
            val6 +=      Dep(_YY_,_XY_) * Deriv(1,0);
            val6 += 2. * Dep(_YY_,_YY_) * Deriv(1,1);
            val6 +=      Dep(_YY_,_YZ_) * Deriv(1,2);
            val6 +=      Dep(_YZ_,_XY_) * Deriv(2,0);
            val6 += 2. * Dep(_YZ_,_YY_) * Deriv(2,1);
            val6 +=      Dep(_YZ_,_YZ_) * Deriv(2,2);
            val6 *= 0.5;
            ek(in*nstate+1,jn*nstate+1) += weight * val6;
            
            val7  =      Dep(_XY_,_XZ_) * Deriv(0,0);
            val7 +=      Dep(_XY_,_YZ_) * Deriv(0,1);
            val7 += 2. * Dep(_XY_,_ZZ_) * Deriv(0,2);//
            val7 +=      Dep(_YY_,_XZ_) * Deriv(1,0);
            val7 +=      Dep(_YY_,_YZ_) * Deriv(1,1);
            val7 += 2. * Dep(_YY_,_ZZ_) * Deriv(1,2);//
            val7 +=      Dep(_YZ_,_XZ_) * Deriv(2,0);
            val7 +=      Dep(_YZ_,_YZ_) * Deriv(2,1);
            val7 += 2. * Dep(_YZ_,_ZZ_) * Deriv(2,2);
            val7 *= 0.5;
            ek(in*nstate+1,jn*nstate+2) += weight * val7;
            
            //Third equation Dot[Sigma3, gradV3]
            val8  = 2. * Dep(_XZ_,_XX_) * Deriv(0,0);
            val8 +=      Dep(_XZ_,_XY_) * Deriv(0,1);
            val8 +=      Dep(_XZ_,_XZ_) * Deriv(0,2);
            val8 += 2. * Dep(_YZ_,_XX_) * Deriv(1,0);
            val8 +=      Dep(_YZ_,_XY_) * Deriv(1,1);
            val8 +=      Dep(_YZ_,_XZ_) * Deriv(1,2);
            val8 += 2. * Dep(_ZZ_,_XX_) * Deriv(2,0);//
            val8 +=      Dep(_ZZ_,_XY_) * Deriv(2,1);//
            val8 +=      Dep(_ZZ_,_XZ_) * Deriv(2,2);
            val8 *= 0.5;
            ek(in*nstate+2,jn*nstate+0) += weight * val8;
            
            val9  =      Dep(_XZ_,_XY_) * Deriv(0,0);
            val9 += 2. * Dep(_XZ_,_YY_) * Deriv(0,1);
            val9 +=      Dep(_XZ_,_YZ_) * Deriv(0,2);
            val9 +=      Dep(_YZ_,_XY_) * Deriv(1,0);
            val9 += 2. * Dep(_YZ_,_YY_) * Deriv(1,1);
            val9 +=      Dep(_YZ_,_YZ_) * Deriv(1,2);
            val9 +=      Dep(_ZZ_,_XY_) * Deriv(2,0);//
            val9 += 2. * Dep(_ZZ_,_YY_) * Deriv(2,1);//
            val9 +=      Dep(_ZZ_,_YZ_) * Deriv(2,2);
            val9 *= 0.5;
            ek(in*nstate+2,jn*nstate+1) += weight * val9;
            
            val10  =      Dep(_XZ_,_XZ_) * Deriv(0,0);
            val10 +=      Dep(_XZ_,_YZ_) * Deriv(0,1);
            val10 += 2. * Dep(_XZ_,_ZZ_) * Deriv(0,2);
            val10 +=      Dep(_YZ_,_XZ_) * Deriv(1,0);
            val10 +=      Dep(_YZ_,_YZ_) * Deriv(1,1);
            val10 += 2. * Dep(_YZ_,_ZZ_) * Deriv(1,2);
            val10 +=      Dep(_ZZ_,_XZ_) * Deriv(2,0);
            val10 +=      Dep(_ZZ_,_YZ_) * Deriv(2,1);
            val10 += 2. * Dep(_ZZ_,_ZZ_) * Deriv(2,2);//
            val10 *= 0.5;
            ek(in*nstate+2,jn*nstate+2) += weight * val10;
            
        }//jn
    }//in
    
#ifdef LOG4CXX
    if(elastoplasticLogger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "<<< TPZMatElastoPlastic<T,TMEM>::Contribute ***";
        //sout << " Resultant rhs vector:\n" << ef;
        LOGPZ_DEBUG(elastoplasticLogger,sout.str().c_str());
    }
    //#ifdef PZDEBUG
    //   if ( !ek.VerifySymmetry( 1.e-8 ) )
    //    {
    //        std::stringstream sout;
    //        sout << "<<< TPZMatElastoPlastic<T,TMEM>::Contribute *** NON SYMMETRIC CONTRIBUTE SUBMATRIX";
    //        LOGPZ_WARN(elastoplasticLogger,sout.str().c_str());
    //    }
    //#endif
#endif
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::ContributeBC(TPZMaterialData &data,
                                               REAL weight,
                                               TPZFMatrix<REAL> &ek,
                                               TPZFMatrix<REAL> &ef,
                                               TPZBndCond &bc)
{
#ifdef LOG4CXX
    if(elastoplasticLogger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << ">>> TPZMatElastoPlastic<T,TMEM>::ContributeBC *** with bc.Type()=" << bc.Type();
        LOGPZ_DEBUG(elastoplasticLogger,sout.str().c_str());
    }
#endif
    TPZFMatrix<REAL> &phi = data.phi;
    
    const REAL BIGNUMBER  = 1.e16;
    
    int dim = Dimension();
    int nstate = NStateVariables();
    
    const int phr = phi.Rows();
    int in,jn,idf,jdf;
    REAL v2[3];
    v2[0] = bc.Val2()(0,0);
    v2[1] = bc.Val2()(1,0);
    v2[2] = bc.Val2()(2,0);
    
    
    TPZFMatrix<REAL> &v1 = bc.Val1();
    //bc.Print(cout);
    //cout << "val2:  " << v2[0]          << ' ' << v2[1]          << ' ' << v2[2]          << endl;
    switch (bc.Type()) {
        case 0: // Dirichlet condition
            for(in = 0 ; in < phr; in++) {
                ef(nstate*in+0,0) += BIGNUMBER * (v2[0] - data.sol[0][0]) * phi(in,0) * weight;
                ef(nstate*in+1,0) += BIGNUMBER * (v2[1] - data.sol[0][1]) * phi(in,0) * weight;
                ef(nstate*in+2,0) += BIGNUMBER * (v2[2] - data.sol[0][2]) * phi(in,0) * weight;
                for (jn = 0 ; jn < phr; jn++) {
                    ek(nstate*in+0,nstate*jn+0) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight;
                    ek(nstate*in+1,nstate*jn+1) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight;
                    ek(nstate*in+2,nstate*jn+2) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight;
                }//jn
            }//in
            break;
            
        case 1: // Neumann condition
            for(in = 0 ; in < phi.Rows(); in++) {
                ef(nstate*in+0,0) += v2[0] * phi(in,0) * weight;
                ef(nstate*in+1,0) += v2[1] * phi(in,0) * weight;
                ef(nstate*in+2,0) += v2[2] * phi(in,0) * weight;
            }
            break;
            
        case 2: // Mixed condition
            for(in = 0 ; in < phi.Rows(); in++) {
                ef(nstate*in+0,0) += v2[0] * phi(in,0) * weight;
                ef(nstate*in+1,0) += v2[1] * phi(in,0) * weight;
                ef(nstate*in+2,0) += v2[2] * phi(in,0) * weight;
                for(jn=0; jn<phi.Rows(); jn++)
                {
                    for(idf=0; idf<3; idf++) for(jdf=0; jdf<3; jdf++)
                    {
                        ek(nstate*in+idf,nstate*jn+jdf) += bc.Val1()(idf,jdf);
                    }
                }
            }//in
            break;
            
        case 3: // Directional Null Dirichlet - displacement is set to null in the non-null vector component direction
            for(in = 0 ; in < phr; in++) {
                ef(nstate*in+0,0) += BIGNUMBER * (0. - data.sol[0][0]) * v2[0] * phi(in,0) * weight;
                ef(nstate*in+1,0) += BIGNUMBER * (0. - data.sol[0][1]) * v2[1] * phi(in,0) * weight;
                ef(nstate*in+2,0) += BIGNUMBER * (0. - data.sol[0][2]) * v2[2] * phi(in,0) * weight;
                for (jn = 0 ; jn < phr; jn++) {
                    ek(nstate*in+0,nstate*jn+0) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight * v2[0];
                    ek(nstate*in+1,nstate*jn+1) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight * v2[1];
                    ek(nstate*in+2,nstate*jn+2) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight * v2[2];
                }//jn
            }//in
            break;
            
        case 4: // stressField Neumann condition
            for(in = 0; in < dim; in ++)
                v2[in] = - ( v1(in,0) * data.normal[0] +
                            v1(in,1) * data.normal[1] +
                            v1(in,2) * data.normal[2] );
            // The normal vector points towards the neighbour. The negative sign is there to
            // reflect the outward normal vector.
            for(in = 0 ; in < phi.Rows(); in++) {
                ef(nstate*in+0,0) += v2[0] * phi(in,0) * weight;
                ef(nstate*in+1,0) += v2[1] * phi(in,0) * weight;
                ef(nstate*in+2,0) += v2[2] * phi(in,0) * weight;
                //    cout << "normal:" << data.normal[0] << ' ' << data.normal[1] << ' ' << data.normal[2] << endl;
                //    cout << "val2:  " << v2[0]  << endl;
            }
            break;
            
        case 5://PRESSAO
            for(in = 0 ; in < phi.Rows(); in++)
            {
                ef(nstate*in+0,0) += v2[0] * phi(in,0) * weight * (data.normal[0]);
                ef(nstate*in+1,0) += v2[0] * phi(in,0) * weight * (data.normal[1]);
                ef(nstate*in+2,0) += v2[0] * phi(in,0) * weight * (data.normal[2]);
            }
            break;
            
        case 6: // Directional Dirichlet - displacement is set to u_D in the non-null vector component direction
            
            REAL v_null[3];
            v_null[0] = bc.Val1()(0,0);
            v_null[1] = bc.Val1()(1,1);
            v_null[2] = bc.Val1()(2,2);
            
            for(in = 0 ; in < phr; in++) {
                ef(nstate*in+0,0) += BIGNUMBER * (v2[0] - data.sol[0][0]) * v_null[0] * phi(in,0) * weight;
                ef(nstate*in+1,0) += BIGNUMBER * (v2[1] - data.sol[0][1]) * v_null[1] * phi(in,0) * weight;
                ef(nstate*in+2,0) += BIGNUMBER * (v2[2] - data.sol[0][2]) * v_null[2] * phi(in,0) * weight;
                for (jn = 0 ; jn < phr; jn++) {
                    ek(nstate*in+0,nstate*jn+0) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight * v_null[0];
                    ek(nstate*in+1,nstate*jn+1) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight * v_null[1];
                    ek(nstate*in+2,nstate*jn+2) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight * v_null[2];
                }//jn
            }//in
            break;
            
        default:
#ifdef LOG4CXX
        {
            std::stringstream sout;
            sout << "<<< TPZMatElastoPlastic<T,TMEM>::ContributeBC *** WRONG BOUNDARY CONDITION TYPE = " << bc.Type();
            LOGPZ_ERROR(elastoplasticLogger,sout.str().c_str());
        }
#endif
            PZError << "TPZMatElastoPlastic::ContributeBC error - Wrong boundary condition type" << std::endl;
    }//switch
    
    //    cout << "normal:" << data.normal[0] << ' ' << data.normal[1] << ' ' << data.normal[2] << endl;
    //    cout << "val2:  " << v2[0] << endl;
    
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ef)
{
#ifdef LOG4CXX
    if(elastoplasticLogger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << ">>> TPZMatElastoPlastic<T,TMEM>::Contribute ***";
        sout << "\nIntegration Point index = " << data.intGlobPtIndex;
        LOGPZ_DEBUG(elastoplasticLogger,sout.str().c_str());
    }
#endif
    
    TPZFMatrix<REAL> &dphi = data.dphix, dphiXYZ;
    TPZFMatrix<REAL> &phi  = data.phi;
    TPZFMatrix<REAL> &axes = data.axes, axesT;
    TPZManVector<REAL,3> &x = data.x;
    
    // rotating the shape functions to the XYZ coordinates
    axes.Transpose(&axesT);
    axesT.Multiply(dphi,dphiXYZ);
    /*
     cout << "\n phi(" << data.intPtIndex << ") =";
     for(int i = 0; i < data.phi.Rows(); i++)cout << " " << data.phi(i,0);
     cout << endl << dphiXYZ;
     cout << endl << axes;
     */
    const int phr = phi.Rows();
    if(this->fForcingFunction)
        this->fForcingFunction->Execute(x,this->fForce);
    
    //this matrix will store {{dvdx*dudx, dvdx*dudy, dvdx*dudz},
    //{dvdy*dudx, dvdy*dudy, dvdy*dudz},
    //{dvdz*dudx, dvdz*dudy, dvdz*dudz}}
    TPZFNMatrix<9>  Deriv(3,3);
    //    TPZFNMatrix<36> Dep(6,6);
    TPZFNMatrix<6>  DeltaStrain(6,1);
    TPZFNMatrix<6>  Stress(6,1);//, StressN(6,1);
    
    this->ComputeDeltaStrainVector(data, DeltaStrain);
    this->ApplyDeltaStrain(data,DeltaStrain,Stress);
    //    this->ApplyDeltaStrainComputeDep(data, DeltaStrain, Stress, Dep);
    
    //int dim = Dimension();
    int nstate = NStateVariables();
    REAL val;
    
    int in;
    for(in = 0; in < phr; in++) { //in: test function index
        
        // fForce represents the gravity acceleration
        //First equation: fb and fk
        val  = fRhoB * fForce[0] * phi(in,0); // fb
        val -= Stress(_XX_,0) * dphiXYZ(0,in); // |
        val -= Stress(_XY_,0) * dphiXYZ(1,in); // fk
        val -= Stress(_XZ_,0) * dphiXYZ(2,in); // |
        ef(in*nstate+0,0) += weight * val;
        
        //Second equation: fb and fk
        val  = fRhoB * fForce[1] * phi(in,0); // fb
        val -= Stress(_XY_,0) * dphiXYZ(0,in); // |
        val -= Stress(_YY_,0) * dphiXYZ(1,in); // fk
        val -= Stress(_YZ_,0) * dphiXYZ(2,in); // |
        ef(in*nstate+1,0) += weight * val;
        
        //third equation: fb and fk
        val  = fRhoB * fForce[2] * phi(in,0); // fb
        val -= Stress(_XZ_,0) * dphiXYZ(0,in); // |
        val -= Stress(_YZ_,0) * dphiXYZ(1,in); // fk
        val -= Stress(_ZZ_,0) * dphiXYZ(2,in); // |
        ef(in*nstate+2,0) += weight * val;
        
    }//in
    
#ifdef LOG4CXX
    if(elastoplasticLogger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "<<< TPZMatElastoPlastic<T,TMEM>::Contribute ***";
        //sout << " Resultant rhs vector:\n" << ef;
        LOGPZ_DEBUG(elastoplasticLogger,sout.str().c_str());
    }
    //#ifdef PZDEBUG
    //   if ( !ek.VerifySymmetry( 1.e-8 ) )
    //    {
    //        std::stringstream sout;
    //        sout << "<<< TPZMatElastoPlastic<T,TMEM>::Contribute *** NON SYMMETRIC CONTRIBUTE SUBMATRIX";
    //        LOGPZ_WARN(elastoplasticLogger,sout.str().c_str());
    //    }
    //#endif
#endif
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::ContributeBC(TPZMaterialData &data,
                                               REAL weight,
                                               TPZFMatrix<REAL> &ef,
                                               TPZBndCond &bc)
{
    TPZMaterial::ContributeBC(data, weight, ef, bc);//not efficient but here to remember reimplementing it when ContributeBC becomes robust
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::Errors(TPZVec<REAL> &x,TPZVec<REAL> &u, TPZFMatrix<REAL> &dudx,
                                         TPZFMatrix<REAL> &axes, TPZVec<REAL> &flux,
                                         TPZVec<REAL> &u_exact,TPZFMatrix<REAL> &du_exact,TPZVec<REAL> &values)
{
    int i, j;
    
    /** L2 norm */
    REAL L2 = 0.;
    for(i = 0; i < 3; i++) L2 += (u[i] - u_exact[i]) * (u[i] - u_exact[i]);
    
    /** H1 semi-norm */
    REAL SemiH1 = 0.;
    for(i = 0; i < 3; i++) for(j = 0; j < 3; j++) SemiH1 += (dudx(i,j) - du_exact(i,j)) * (dudx(i,j) - du_exact(i,j));
    
    /** H1 norm */
    REAL H1 = L2 + SemiH1;
    
    //values[1] : eror em norma L2
    values[1]  = L2;
    
    //values[2] : erro em semi norma H1
    values[2] = SemiH1;
    
    //values[0] : erro em norma H1 <=> norma Energia
    values[0]  = H1;
    
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::ComputeStrainVector(TPZMaterialData & data, TPZFMatrix<REAL> &Strain)
{
    ComputeDeltaStrainVector(data, Strain);
    
    TPZTensor<REAL> & EpsT = this->MemItem(data.intGlobPtIndex).fPlasticState.fEpsT;
    
    int i;
    for( i = 0; i < 6; i++ )Strain(i,0) = Strain(i,0) + EpsT.fData[i];
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::ComputeDeltaStrainVector(TPZMaterialData & data, TPZFMatrix<REAL> &DeltaStrain)
{
    TPZFNMatrix<9> DSolXYZ(3,3,0.);
    data.axes.Multiply(data.dsol[0],DSolXYZ,1/*transpose*/);
    cout << "\n data dsol \n";
    data.dsol[0].Print("data.sol");
    DeltaStrain.Redim(6,1);
    DeltaStrain(_XX_,0) = DSolXYZ(0,0);
    DeltaStrain(_YY_,0) = DSolXYZ(1,1);
    DeltaStrain(_ZZ_,0) = DSolXYZ(2,2);
    DeltaStrain(_XY_,0) = 0.5 * ( DSolXYZ(1,0) + DSolXYZ(0,1) );
    DeltaStrain(_XZ_,0) = 0.5 * ( DSolXYZ(2,0) + DSolXYZ(0,2) );
    DeltaStrain(_YZ_,0) = 0.5 * ( DSolXYZ(2,1) + DSolXYZ(1,2) );
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::ComputeStressVector(TPZMaterialData & data, TPZFMatrix<REAL> &Stress)
{
    
    TPZFNMatrix<6> DeltaStrain;
    ComputeDeltaStrainVector(data, DeltaStrain);
    Stress.Redim(6,1);
    ApplyDeltaStrain(data, DeltaStrain, Stress);
    
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::CheckConvergence(TPZMaterialData & data, TPZFMatrix<REAL> & DeltaStrain)
{
    int intPt = data.intGlobPtIndex;//, plasticSteps;
    T plasticloc(fPlasticity);
    plasticloc.SetState(this->MemItem(intPt).fPlasticState);
    TPZTensor<REAL> deps;
    deps.CopyFrom(DeltaStrain);
    
    REAL alfa =1.e-6;
    REAL alfa2 = 2.e-6;
    TPZTensor<REAL> part1,part2,part3,temp;
    TPZTensor<REAL> Alfa1DeltaEps, Alfa2DeltaEps,Eps(plasticloc.GetState().fEpsT);
    Alfa1DeltaEps.CopyFrom(DeltaStrain);
    Alfa2DeltaEps.CopyFrom(DeltaStrain);
    TPZFNMatrix<36,REAL> DEP(6,6);
    
    Alfa1DeltaEps*=alfa;
    Alfa2DeltaEps*=alfa2;
    temp=Eps;
    temp+=Alfa1DeltaEps;
    plasticloc.ApplyStrainComputeSigma(temp,part1);
    plasticloc.ApplyStrainComputeDep(Eps,part2,DEP);
    TPZFNMatrix<6,REAL> part3temp(6,1),tempAlfa1DeltaEps(6,1);
    for(int i=0;i<6;i++)
    {
        tempAlfa1DeltaEps(i,0)=Alfa1DeltaEps.fData[i];
    }
    DEP.Multiply(tempAlfa1DeltaEps, part3temp);
    part3.CopyFrom(part3temp);
    TPZTensor<REAL> e1(part1);
    e1-=part2;
    e1-=part3;
    
    part1*=0.;
    part3*=0.;
    part3temp*=0.;
    temp*=0.;
    
    temp=Eps;
    temp+=Alfa2DeltaEps;
    plasticloc.ApplyStrainComputeSigma(temp,part1);
    for(int i=0;i<6;i++)
    {
        tempAlfa1DeltaEps(i,0)=Alfa2DeltaEps.fData[i];
    }
    DEP.Multiply(tempAlfa1DeltaEps, part3temp);
    part3.CopyFrom(part3temp);
    TPZTensor<REAL> e2(part1);
    e2-=part2;
    e2-=part3;
    REAL n = (log10(Norm(e1))-log10(Norm(e2)))/(log10(alfa)-log10(alfa2));
    
#ifdef LOG4CXX
    if(ceckconvlogger->isDebugEnabled())
    {
        std::stringstream sout;
        TPZManVector<REAL,3> phi(3,1);
        plasticloc.Phi(Eps,phi);
        sout << "DEP "<< DEP << std::endl;
        sout << "tempAlfa1DeltaEps "<< tempAlfa1DeltaEps << std::endl;
        sout << "Phi "<< phi << std::endl;
        sout << "Integration Point "<< intPt << std::endl;
        sout << "n = " << n << std::endl;
        LOGPZ_DEBUG(ceckconvlogger, sout.str())
    }
#endif
    
    
    
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::ApplyDeltaStrainComputeDep(TPZMaterialData & data, TPZFMatrix<REAL> & DeltaStrain,
                                                             TPZFMatrix<REAL> & Stress, TPZFMatrix<REAL> & Dep)
{
    
    int intPt = data.intGlobPtIndex;
    T plasticloc(fPlasticity);
    
    plasticloc.SetState(this->MemItem(intPt).fPlasticState);
    
    UpdateMaterialCoeficients(data.x,plasticloc);
    TPZTensor<REAL> EpsT, Sigma;
    
    EpsT.CopyFrom(DeltaStrain);
    EpsT.Add(plasticloc.GetState().fEpsT, 1.);
    
    plasticloc.ApplyStrainComputeSigma(EpsT, Sigma, &Dep);
    Sigma.CopyTo(Stress);
    
    if(TPZMatWithMem<TMEM>::fUpdateMem)
    {
        this->MemItem(intPt).fSigma        = Sigma;
        this->MemItem(intPt).fPlasticState = plasticloc.GetState();
        this->MemItem(intPt).fPlasticSteps = plasticloc.IntegrationSteps();
        int solsize = data.sol[0].size();
        for(int i=0; i<solsize; i++)
        {
            this->MemItem(intPt).fDisplacement[i] += data.sol[0][i];
        }
    }
    
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::UpdateMaterialCoeficients(TPZVec<REAL> &x,T & plasticity)
{
    /// Nothing to do
    return;
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::ApplyDeltaStrain(TPZMaterialData & data, TPZFMatrix<REAL> & Strain,
                                                   TPZFMatrix<REAL> & Stress)
{
    
    int intPt = data.intGlobPtIndex;
    T plasticloc(fPlasticity);
    
    plasticloc.SetState(this->MemItem(intPt).fPlasticState);
    
    UpdateMaterialCoeficients(data.x,plasticloc);
    TPZTensor<REAL> EpsT, Sigma;
    
    EpsT.CopyFrom(Strain);
    EpsT.Add(plasticloc.GetState().fEpsT, 1.);
    
    plasticloc.ApplyStrainComputeSigma(EpsT, Sigma);
    Sigma.CopyTo(Stress);
    
    if(TPZMatWithMem<TMEM>::fUpdateMem)
    {
        this->MemItem(intPt).fSigma        = Sigma;
        this->MemItem(intPt).fPlasticState = plasticloc.GetState();
        this->MemItem(intPt).fPlasticSteps = plasticloc.IntegrationSteps();
        int solsize = data.sol[0].size();
        for(int i=0; i<solsize; i++)
        {
            this->MemItem(intPt).fDisplacement[i] += data.sol[0][i];
        }
    }
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::EigenValues(TPZFMatrix<REAL> & vectorTensor, TPZVec<REAL> & ev)
{
    TPZFNMatrix<9> Tensor(3,3);
    ev.Resize(3);
    this->vectorToTensor(vectorTensor, Tensor);
    int64_t numiterations = 1000;
    
#ifdef PZDEBUG
    bool result = Tensor.SolveEigenvaluesJacobi(numiterations, fTol, &ev);
    if (result == false){
        PZError << __PRETTY_FUNCTION__ << " - ERROR! - result = false - numiterations = " << numiterations << " - tol = " << fTol << std::endl;
#ifdef LOG4CXX
        {
            std::stringstream sout;
            sout << "<<< TPZMatElastoPlastic<T,TMEM>::EigenValues *** not solved within " << numiterations << " iterations";
            sout << "\n vectorTensor = " << vectorTensor;
            LOGPZ_ERROR(elastoplasticLogger,sout.str().c_str());
        }
#endif
    }
#else
    Tensor.SolveEigenvaluesJacobi(numiterations, fTol, &ev);
#endif
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::EigenVectors(TPZFMatrix<REAL> &vectorTensor, TPZVec< REAL > &Solout, int direction)
{
    TPZFNMatrix<9> Tensor(3,3);
    this->vectorToTensor(vectorTensor, Tensor);
    
    TPZManVector<REAL,3> Eigenvalues(3);
    TPZFNMatrix<9> Eigenvectors(3,3);
    
    int64_t numiterations = 1000;
#ifdef PZDEBUG
    bool result = Tensor.SolveEigensystemJacobi(numiterations, fTol, Eigenvalues, Eigenvectors);
    if (result == false){
        PZError << __PRETTY_FUNCTION__ << " - ERROR! - result = false - numiterations = " << numiterations << " - tol = " << fTol << std::endl;
#ifdef LOG4CXX
        {
            std::stringstream sout;
            sout << "<<< TPZMatElastoPlastic<T,TMEM>::EigenVectors *** not solved within " << numiterations << " iterations";
            sout << "\n vectorTensor = " << vectorTensor;
            LOGPZ_ERROR(elastoplasticLogger,sout.str().c_str());
        }
#endif
    }
#else
    Tensor.SolveEigensystemJacobi(numiterations, fTol, Eigenvalues, Eigenvectors);
#endif
    Solout.Resize(3);
    for(int i = 0; i < 3; i++) Solout[i] = Eigenvectors(direction,i);
}

template <class T, class TMEM>
TPZMaterial * TPZMatElastoPlastic<T,TMEM>::NewMaterial()
{
    return new TPZMatElastoPlastic<T,TMEM>(*this);
}
/*
 void TPZMatElastoPlastic::SetData(std::istream &data)
 {
 TPZMaterial::SetData(data);
 data >> fDeltaT; // to be removed in the elastoplastic material and readded to the poroelastoplastic material
 }*/

#include "TPZSandlerExtended.h"
#include "TPZPlasticStepPV.h"
#include "TPZYCMohrCoulombPV.h"

template <class T, class TMEM>
std::string TPZMatElastoPlastic<T,TMEM>::Name() {
    return "TPZMatElastoPlastic<T,TMEM>";
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T, TMEM>::Write(TPZStream &buf, int withclassid) const {
    TPZMatWithMem<TMEM>::Write(buf, withclassid);
    
    buf.Write(&fForce[0], 3);
    buf.Write(&fPostProcessDirection[0], 3);
    fPlasticity.Write(buf, withclassid);
    buf.Write(&fTol, 1);
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T, TMEM>::Read(TPZStream &buf, void *context) {
    //    TPZSavable::Read(buf, context);
    
    TPZMatWithMem<TMEM>::Read(buf, context);
    
    buf.Read(&fForce[0], 3);
    buf.Read(&fPostProcessDirection[0], 3);
    fPlasticity.Read(buf, context);
    buf.Read(&fTol, 1);
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::SetTol(const REAL & tol)
{
    fTol = tol;
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::SetBulkDensity(const REAL & bulk)
{
    fRhoB = bulk;
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::vectorToTensor(const TPZFMatrix<REAL> & vectorTensor, TPZFMatrix<REAL> & Tensor)
{
    TPZTensor<REAL> vecT;
    vecT.CopyFrom(vectorTensor);
    vecT.CopyToTensor(Tensor);
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::FillDataRequirements(TPZMaterialData &data){
    
    TPZMatWithMem<TMEM>::FillDataRequirements(data);
    
    data.fNeedsSol = true;
    data.fNeedsNormal = false;
    data.fNeedsHSize = false;
    data.fNeedsNeighborCenter = false;
}

template <class T, class TMEM>
void TPZMatElastoPlastic<T,TMEM>::FillBoundaryConditionDataRequirement(int type,TPZMaterialData &data)
{
    
    
    //TPZMatWithMem<TMEM>::FillBoundaryConditionDataRequirement(type,data);
    data.fNeedsSol = true;
    data.fNeedsNormal = true;
}
