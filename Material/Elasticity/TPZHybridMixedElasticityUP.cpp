//
// Created by Giovane Avancini on 24/08/23.
//

#include <pzfmatrix.h>
#include <TPZBndCondT.h>
#include <pzaxestools.h>
#include <pzlog.h>

#include "TPZHybridMixedElasticityUP.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.elasticmaterial");
#endif

TPZHybridMixedElasticityUP::TPZHybridMixedElasticityUP(): TBase() {}

TPZHybridMixedElasticityUP::TPZHybridMixedElasticityUP(int matID,
                                                       int dimension,
                                                       REAL young_modulus,
                                                       REAL poisson,
                                                       AnalysisType analysisType,
                                                       REAL thickness) : 
                                                       TBase(matID),
                                                       fdimension(dimension),
                                                       fyoung(young_modulus),
                                                       fpoisson(poisson),
                                                       fAnalysisType(analysisType),
                                                       fthickness(thickness)
{
    flambda = fyoung * fpoisson / ((1.0 + fpoisson) * (1.0 - 2.0 * fpoisson));
    fmu = 0.5 * fyoung / (1.0 + fpoisson);
    feta = (fabs(fpoisson - 0.5) <= 1.0e-9) ? 1.e9 : 0.0;

    switch (fAnalysisType)
    {
        case AnalysisType::EGeneral:
        {
            fbulk = flambda + 2.0 * fmu / 3.0;;
            break;
        }
        case AnalysisType::EPlaneStrain:
        {
            fbulk = flambda + fmu;
            break;
        }
        case AnalysisType::EPlaneStress:
        {
            fbulk = (fabs(fpoisson - 0.5) <= 1.0e-9) ? 3.0*fmu : fmu * (2.0*fmu + 3.0*flambda) / (flambda + 2.0*fmu); //For full incompressible case. Replacing lambda with E and poisson didnt yield correct results
            break;
        }
        default:
        {
            std::cout << "Wrong analysis type." << std::endl;
            DebugStop();
            break;
        }
    }
}

TPZHybridMixedElasticityUP::~TPZHybridMixedElasticityUP() {}

void TPZHybridMixedElasticityUP::Contribute(const TPZVec<TPZMaterialDataT<STATE>>& datavec, REAL weight, TPZFMatrix<STATE>& ek, TPZFMatrix<STATE>& ef)
{
    int64_t nShapeU = datavec[EUindex].fVecShapeIndex.NElements(); // number of displacements Hdiv shape functions
    TPZFNMatrix<150, REAL> PhiU(fdimension, nShapeU, 0.0);
    TPZFNMatrix<20, REAL>& divPhiU = datavec[EUindex].divphi;

    TPZFMatrix<REAL>& PhiP = datavec[EPindex].phi;
    int64_t nShapeP = PhiP.Rows(); // number of pressure L2 shape functions

    const int n = fdimension * (fdimension + 1) / 2; //number of independent components of stress tensor in voight notation
    TPZFNMatrix<150, REAL> Strain(n, nShapeU, 0.0); //Using voight notation
    
    if (datavec[EUindex].fNeedsDeformedDirectionsFad)
    {
        for (int64_t j = 0; j < nShapeU; j++)
        {
            int cont = fdimension-1;
            for (int i = 0; i < fdimension; i++)
            {
                PhiU(i, j) = datavec[EUindex].fDeformedDirectionsFad(i, j).val();
                Strain(i,j) = datavec[EUindex].fDeformedDirectionsFad(i, j).fastAccessDx(i); //diagonal part of infinitesimal strain tensor
                for (int64_t k = i+1; k < fdimension && k != i; k++)
                {
                    Strain(++cont,j) = 0.5*(datavec[EUindex].fDeformedDirectionsFad(i, j).fastAccessDx(k) + datavec[EUindex].fDeformedDirectionsFad(k, j).fastAccessDx(i)); //off diagonal part of infinitesimal strain tensor
                }
            }
        }
    }

    TPZFNMatrix<3,REAL> SourceTerm(fdimension, 1.0, 0.0);
    TPZVec<REAL> sourceAux(3);
    if (this->HasForcingFunction())
    {
        this->ForcingFunction()(datavec[EUindex].x, sourceAux);
        for (int64_t i = 0; i < fdimension; i++)
        {
            SourceTerm(i,0) = sourceAux[i];
        }
    }

    TPZFNMatrix<36, REAL> D(n, n, 0.0); 
    DeviatoricElasticityTensor(D);
    for (int i = fdimension; i < n; i++) // The terms related to the off diagonal part of strain tensor is multiplied by 2 to account for its symmetry
        D(i,i) *= 2.0;
    
    //Body Forces contribution
    ef.AddContribution(0, 0, PhiU, true, SourceTerm, false, weight);

    //Stiffness Matrix A isochoric contribution
    TPZFNMatrix<150, REAL> aux;
    D.Multiply(Strain, aux);
    
    REAL factor = weight;
    ek.AddContribution(0, 0, Strain, true, aux, false, factor);

    //Stiffness Matrix A volumetric contribution (To get rid off zero-energy volumetric deformation modes under fully incompressibility)
    factor = feta * weight;
    ek.AddContribution(0 , 0, divPhiU, false, divPhiU, true, factor);

    //Divergence Matrix B contribution
    factor = -1.0 * weight;
    ek.AddContribution(0, nShapeU, divPhiU, false, PhiP, true, factor);

    //Divergence Matrix BT contribution
    ek.AddContribution(nShapeU, 0, PhiP, false, divPhiU, true, factor);

    //Bulk Matrix C contribution
    factor = -(1.0 / fbulk) * weight;
    ek.AddContribution(nShapeU, nShapeU, PhiP, false, PhiP, true, factor);
    
    if(datavec.size() > 2) //Static condensation in incompressibility regime
    {
        TPZFMatrix<REAL>& PhiUM = datavec[EVMindex].phi;
        TPZFMatrix<REAL>& phipM = datavec[EPMindex].phi;
        
        // Pressure and distributed displacement
        for(int j = 0; j < nShapeP; j++)
        {
            ek(nShapeU + nShapeP, nShapeU + j      ) += PhiP(j,0) * PhiUM(0,0) * weight;
            ek(nShapeU + j      , nShapeU + nShapeP) += PhiP(j,0) * PhiUM(0,0) * weight;
        }
        
        // Injection and average-pressure
        ek(nShapeU + nShapeP + 1, nShapeU + nShapeP    ) += PhiUM(0,0) * phipM(0,0) * weight;
        ek(nShapeU + nShapeP    , nShapeU + nShapeP + 1) += PhiUM(0,0) * phipM(0,0) * weight;
    }

#ifdef PZ_LOG
    if(logger.isDebugEnabled())
    {
        std::stringstream sout;
        D.Print("tensor D", sout, EMathematicaInput);
        ek.Print("ek", sout, EMathematicaInput);
        ef.Print("ef", sout, EMathematicaInput);
        sout << std::endl << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
}

void TPZHybridMixedElasticityUP::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc)
{
    int index = (bc.Type() == 0 || bc.Type() == 2)? EUindex : EPindex;
    index = EUindex;

    TPZFNMatrix<150, REAL> PhiU = datavec[EUindex].phi;
    TPZFMatrix<REAL>& PhiP = datavec[EPindex].phi;
    
    int64_t nShapeU = PhiU.Rows();
    int64_t nShapeP = PhiP.Rows();
            
    TPZFNMatrix<20,STATE> val1(3,3,0.0);
    TPZManVector<STATE, 3> val2(3, 0.0);

    if (bc.HasForcingFunctionBC())
    {
        TPZVec<STATE> uVal;
        TPZFMatrix<STATE> gradVal;
        bc.ForcingFunctionBC()(datavec[index].x, val2, val1);
    }
    else
    {
        val1 = bc.Val1();
        val2 = bc.Val2();
    }

    switch (bc.Type())
    {
        case 0: // Normal Displacement
        {
            REAL u_n = val2[0]; //if bc was set in .json file, the normal value was already prescribed
            if (bc.HasForcingFunctionBC()) //if the bc is set through an analytic solution, we need to compute its normal component
            {
                u_n = 0.0;
                for (int i = 0; i < fdimension; i++)
                    u_n += val2[i] * datavec[index].normal[i];
            }

            REAL factor = fBigNumber * weight;

            for (int64_t j = 0; j < nShapeU; j++)
            {
                ef(j) += u_n * PhiU(j,0) * factor;

                for (int64_t i = 0; i < nShapeU; i++)
                {
                    ek(i, j) += PhiU(i,0) * PhiU(j,0) * factor;
                }
            }
            break;
        }

        case 1: // Tangential displacement
        {
            TPZManVector<REAL,3> u_t = {0.0, 0.0, 0.0}; //for tangential bc, a vector is prescribed, so we take the inner product with the local tangential axe    
            for (int i = 0; i < fdimension-1; i++) //number of tangential components
                for (int j = 0 ; j < fdimension; j++)
                    u_t[i] += val2[j] * datavec[index].axes(i,j);
            
            // for (int64_t i = 0; i < nShapeP; i++)
            // {
            //     for (int j = 0; j < fdimension-1; j++)
            //     {
            //         int64_t index = (fdimension-1)*i+j;
            //         ef(index) += PhiP(i,0) * u_t[j] * weight;
            //     }
            // }
            REAL factor = fBigNumber * weight;

            for (int64_t j = 0; j < nShapeU; j++)
            {
                for (int k = 0; k < fdimension-1; k++)
                {
                    int64_t index1 = (fdimension-1)*j+k;
                    ef(index1) += -u_t[k] * PhiU(j,0) * factor;

                    for (int64_t i = 0; i < nShapeU; i++)
                    {
                        for (int l = 0; l < fdimension-1; l++)
                        {
                            int64_t index2 = (fdimension-1)*i+l;
                            if (k != l) continue;
                            ek(index1, index2) += PhiU(i,0) * PhiU(j,0) * factor;
                        }
                    }
                }
            }
            break;
        }

        case 2: // Normal Stress
        {
            REAL sigma_nn = val2[0]; //if bc was set in .json file, the normal value was already prescribed
            if (bc.HasForcingFunctionBC()) //if the bc is set through an analytic solution, we need to compute its normal component from the displacement gradient
            {
                const int n = fdimension * (fdimension + 1) / 2;

                TPZFNMatrix<6,REAL> sigmavoight(n,1,0.0), sigmavoight2(n,1,0.0);
                DeviatoricStressTensor(val1, sigmavoight);

                REAL p_exact = -datavec[index].x[1]*datavec[index].x[2] / 3.;

                TPZFNMatrix<9, STATE> sigma(3, 3, 0.0);
                int cont = fdimension-1;
                for (int i = 0; i < fdimension; i++)
                {
                    sigma(i,i) = sigmavoight(i,0);
                    for (int j = i+1; j < fdimension; j++)
                    {
                        sigma(i,j) = sigmavoight(++cont,0);
                        sigma(j,i) = sigmavoight(cont,0);
                    }
                }
                
                TPZFNMatrix<3,REAL> sigma_n(fdimension,1,0.0);
                
                for (int i = 0; i < fdimension; i++)
                    for (int j = 0; j < fdimension; j++)
                        sigma_n(i,0) += sigma(i,j) * datavec[index].normal[j];

                sigma_nn = 0.0;
                for (int i = 0; i < fdimension; i++)
                    sigma_nn += sigma_n[i] * datavec[index].normal[i];
            }

            for (int64_t i = 0; i < nShapeU; i++)
            {
                REAL phi = PhiU(i,0);
                ef(i) += sigma_nn * PhiU(i,0) * weight;
            }
            break;
        }

        case 3: // Tangential Stress
        {
            TPZManVector<REAL,3> sigma_nt(fdimension-1,0.);
            if (bc.HasForcingFunctionBC()) //if the bc is set through an analytic solution, we need to compute its tangential component from the displacement gradient
            {
                const int n = fdimension * (fdimension + 1) / 2;

                TPZFNMatrix<6,REAL> sigmavoight(n,1,0.0);
                DeviatoricStressTensor(val1, sigmavoight);

                REAL p_exact = -datavec[index].x[1]*datavec[index].x[2] / 3.;

                TPZFNMatrix<9, STATE> sigma(3, 3, 0.0);
                int cont = fdimension-1;
                for (int i = 0; i < fdimension; i++)
                {
                    sigma(i,i) = sigmavoight(i,0);
                    for (int j = i+1; j < fdimension; j++)
                    {
                        sigma(i,j) = sigmavoight(++cont,0);
                        sigma(j,i) = sigmavoight(cont,0);
                    }
                }
                
                TPZManVector<REAL,3> sigma_n(fdimension, 0.0);
                
                for (int i = 0; i < fdimension; i++)
                    for (int j = 0; j < fdimension; j++)
                        sigma_n[i] += sigma(i,j) * datavec[index].normal[j];

                for (int i = 0; i < fdimension-1; i++)
                {
                    for (int j = 0; j < fdimension; j++)
                    {
                        sigma_nt[i] += sigma_n[j] * datavec[index].axes(i,j);
                    }
                }
            }
            else
            {
                for (int i = 0; i < fdimension-1; i++)
                {
                    for (int j = 0; j < fdimension; j++)
                    {
                        sigma_nt[i] += val2[j] * datavec[index].axes(i,j);
                    }
                }
            }

            // REAL factor = fBigNumber * weight;

            // for (int64_t j = 0; j < nShapeP; j++)
            // {
            //     for (int k = 0; k < fdimension-1; k++)
            //     {
            //         int64_t index1 = (fdimension-1)*j+k;
            //         ef(index1) += -sigma_nt[k] * PhiP(j,0) * factor;

            //         for (int64_t i = 0; i < nShapeP; i++)
            //         {
            //             for (int l = 0; l < fdimension-1; l++)
            //             {
            //                 int64_t index2 = (fdimension-1)*i+l;
            //                 if (k != l) continue;
            //                 ek(index1, index2) += PhiP(i,0) * PhiP(j,0) * factor;
            //             }
            //         }
            //     }
            // }
            for (int64_t i = 0; i < nShapeU; i++)
            {
                for (int j = 0; j < fdimension-1; j++)
                {
                    int64_t index = (fdimension-1)*i+j;
                    ef(index) += PhiU(i,0) * -sigma_nt[j] * weight;
                }
            }

            break;
        }

        default:
        {
            std::cout << "ERROR: BOUNDARY NOT IMPLEMENTED" << std::endl;
            DebugStop();
            break;
        }
    }
}

int TPZHybridMixedElasticityUP::VariableIndex(const std::string& name) const {
    
    if(!strcmp("Pressure", name.c_str())) return EPressure;
    if(!strcmp("Displacement", name.c_str())) return EDisplacement;
    if(!strcmp("Force", name.c_str())) return EForce;
    if(!strcmp("Stress", name.c_str())) return EStress;
    if(!strcmp("Strain", name.c_str())) return EStrain;
    if(!strcmp("VonMises", name.c_str())) return EVonMises;
    
    std::cout << "\n\nVar index not implemented\n\n";
    DebugStop();
    
    return 0;
}

int TPZHybridMixedElasticityUP::NSolutionVariables(int var) const{
    
    int aux;
    switch (var) {
        case EPressure: // pressure  [scalar]
        case EVonMises: //VonMises
            aux = 1;
            break;
        case EDisplacement: // displacement [vector]
        case EForce: // external force [vector]
            aux = 3;
            break;
        case EStress: // stress tensor
        case EStrain: // strain tensor
            aux = 9;
            break;
        default:
            std::cout << "\n\nVar index not implemented!!!\n\n";
            DebugStop();
            break;
    }
    return aux;
}

void TPZHybridMixedElasticityUP::Solution(const TPZVec<TPZMaterialDataT<STATE>>& datavec, int var, TPZVec<STATE>& Solout) {
    
    TPZManVector<STATE, 3> u_h = datavec[EUindex].sol[0];
    TPZManVector<STATE, 3> p_h = datavec[EPindex].sol[0];

    TPZFNMatrix<10, STATE> gradU = datavec[EUindex].dsol[0];
    TPZFNMatrix<9, STATE> strain(3, 3, 0.0);

    const int n = fdimension * (fdimension + 1) / 2;

    for (int i = 0; i < fdimension; i++)
    {
        for (int j = 0; j < fdimension; j++)
        {
            strain(i, j) = 0.5 * (gradU(i, j) + gradU(j, i));
        }
    }

    if (fAnalysisType == AnalysisType::EPlaneStress)
        strain(2, 2) = (fabs(fpoisson - 0.5) <= 1.0e-9) ? -(strain(0, 0) + strain(1, 1)) : -(flambda / (flambda + 2.0 * fmu)) * (strain(0, 0) + strain(1, 1));

    Solout.Resize(NSolutionVariables(var));
    
    switch(var) {
            
        case EPressure:
        {
            Solout[0] = p_h[0];
            break;
        }
            
        case EDisplacement:
        {
            Solout[0] = u_h[0]; //Vx
            Solout[1] = u_h[1]; //Vy
            Solout[2] = u_h[2]; //Vz

            if (fAnalysisType==AnalysisType::EPlaneStress)
                Solout[2] = strain(2, 2) * fthickness;
            break;
        }
            
        case EForce:
        {
            TPZVec<STATE> f(3,0.);
            
            if(this->HasForcingFunction()){
                this->ForcingFunction()(datavec[EUindex].x, f);
            }
            
            Solout[0] = f[0];
            Solout[1] = f[1];
            Solout[2] = f[2];
            break;
        }
        case EStress:
        {
            TPZFNMatrix<6, STATE> sigmavoight(n, 1, 0.0);
            DeviatoricStressTensor(gradU, sigmavoight);

            TPZFNMatrix<9, STATE> sigma(3, 3, 0.0);
            int cont = fdimension-1;
            for (int i = 0; i < fdimension; i++)
            {
                sigma(i,i) = sigmavoight(i,0) - p_h[0];
                for (int j = i+1; j < fdimension; j++)
                {
                    sigma(i,j) = sigmavoight(++cont,0);
                    sigma(j,i) = sigmavoight(cont,0);
                }
            }
            
            if (fAnalysisType == AnalysisType::EPlaneStrain)
                sigma(2,2) = fthickness * fpoisson * (sigma(0,0) + sigma(1,1));

            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    Solout[i * 3 + j] = sigma(i, j);
                }
            }
            break;
        }
        case EStrain:
        {
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    Solout[i * 3 + j] = strain(i, j);
                }
            }
            break;
        }
        case EVonMises:
        {
            TPZManVector<STATE,3> PrincipalStress(3);
            TPZFNMatrix<6, STATE> sigmavoight(n, 1, 0.0);
            StressTensor(gradU, sigmavoight, p_h[0]);
            TPZFNMatrix<9, STATE> sigma(3, 3, 0.0);
            int cont = fdimension-1;
            for (int i = 0; i < fdimension; i++)
            {
                sigma(i,i) = sigmavoight(i,0);
                for (int j = i+1; j < fdimension; j++)
                {
                    sigma(i,j) = sigmavoight(++cont,0);
                    sigma(j,i) = sigmavoight(cont,0);
                }
            }
            if (fAnalysisType == AnalysisType::EPlaneStrain)
                sigma(2,2) = fthickness * fpoisson * (sigma(0,0) + sigma(1,1));
            
            REAL tol = 1.0e-9;
            int64_t numiterations = 1000;
            bool result;
            result = sigma.SolveEigenvaluesJacobi(numiterations, tol, &PrincipalStress);
    #ifdef PZDEBUG        
            if (result == false)
            {
                std::cout << "Error while computing the principal stresses: result == false." << std::endl;
                DebugStop();
            }
    #endif    
            Solout[0] = (PrincipalStress[0] - PrincipalStress[1]) * (PrincipalStress[0] - PrincipalStress[1]) 
                      + (PrincipalStress[1] - PrincipalStress[2]) * (PrincipalStress[1] - PrincipalStress[2])
                      + (PrincipalStress[2] - PrincipalStress[0]) * (PrincipalStress[2] - PrincipalStress[0]);
            Solout[0] = sqrt(0.5 * Solout[0]);
            break;
        }
 
        default:{
            std::cout << "\n\nVar index not implemented\n\n";
            DebugStop();
        }
    }
}

void TPZHybridMixedElasticityUP::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const
{
    int64_t ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        datavec[idata].fNeedsHSize = true;
        datavec[idata].fNeedsNormal = true;
    }
    datavec[0].fNeedsDeformedDirectionsFad = true;
}

void TPZHybridMixedElasticityUP::FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const
{
    datavec[EUindex].fNeedsSol = false;
    datavec[EPindex].fNeedsSol = false;
    datavec[EUindex].fNeedsNormal = true;
    datavec[EPindex].fNeedsNormal = true;
}

void TPZHybridMixedElasticityUP::Errors(const TPZVec<TPZMaterialDataT<STATE>>& data, TPZVec<REAL>& errors){
    
    // 0: L2 p, 1: L2 p_ex, 2: L2 u, 3: L2 u_ex, 4: L2 divu, 5: L2 divu_ex, 6: L2 sigma, 7: L2 sigma_Ex
    
    if(!HasExactSol()) DebugStop();

    errors.Resize(NEvalErrors());
    
    TPZManVector<STATE, 4> sol_exact(4);
    TPZFNMatrix<9,STATE> gradsol_exact(3,3);
    
    //Getting the exact solution for velocity, pressure and velocity gradient
    fExactSol(data[EUindex].x, sol_exact, gradsol_exact);
    REAL p_exact = -data[EUindex].x[1]*data[EUindex].x[2] / 3.; //Just for computing Bishop beam when poisson is 0.5. Remember to delete later.
    
    //Getting the numeric solution for velocity, pressure and velocity gradient
    TPZManVector<STATE> u_h(3, 0.0);
    TPZManVector<STATE> p_h(1, 0.0);
    TPZFNMatrix<10,STATE> gradv_h = data[EUindex].dsol[0];
    
    this->Solution(data, VariableIndex("Displacement"), u_h);
    this->Solution(data, VariableIndex("Pressure"), p_h);
    
    STATE diffv, diffp, diffdiv;

    diffp = p_h[0] - p_exact;
    errors[0] = diffp * diffp;
    errors[1] = p_exact * p_exact;
    
    errors[2] = 0.0;
    errors[3] = 0.0;
    for(int i = 0; i < fdimension; i++)
    {
        diffv = u_h[i] - sol_exact[i];
        errors[2] += diffv * diffv;
        errors[3] += sol_exact[i] * sol_exact[i];
    }
    
    STATE div_exact = 0.0, div_h = 0.0;
    for(int i = 0; i < fdimension; i++)
    {
        div_exact += gradsol_exact(i,i);
        div_h += gradv_h(i, i);
    }
    
    diffdiv = div_h - div_exact;
    errors[4] = diffdiv * diffdiv;
    errors[5] = div_exact * div_exact;

    const int n = fdimension * (fdimension + 1) / 2;
    TPZFNMatrix<6, REAL> sigma_exact(n,1), sigma_h(n,1);
    StressTensor(gradsol_exact, sigma_exact); //Just for Bishop beam. remember to delete later
    for (int i = 0; i < fdimension; i++)
        //sigma_exact(i,0) -= p_exact;
    StressTensor(gradv_h, sigma_h, p_h[0]);
    
    errors[6] = 0.0;
    errors[7] = 0.0;
    for (int i = 0; i < fdimension; i++)
    {
        const STATE diffsig = sigma_h(i,0) - sigma_exact(i,0);
        errors[6] += diffsig * diffsig;
        errors[7] += sigma_exact(i,0) * sigma_exact(i,0);
    }
    for(int i = fdimension; i < n; i++)
    {
        const STATE diffsig = sigma_h(i,0) - sigma_exact(i,0);
        errors[6] += 2. * diffsig * diffsig;
        errors[7] += sigma_exact(i,0) * sigma_exact(i,0);
    }
}

void TPZHybridMixedElasticityUP::DeviatoricElasticityTensor(TPZFNMatrix<36,REAL>& D)
{
    const int n = fdimension * (fdimension + 1) / 2; //number of independent components of Cauchy stress tensor

    switch (fAnalysisType)
    {
        case AnalysisType::EGeneral:
        {
            REAL c1 = 2.0 * fmu;
            REAL c2 = c1 * 2.0 / 3.0;
            REAL c3 = -c1 * 1.0 / 3.0;
            for (int64_t i = 0; i < fdimension; i++)
            {
                D(i, i) = c2;
                for (int64_t j = i+1; j < fdimension; j++)
                {
                    D(i, j) = c3; //volumetric part
                    D(j, i) = c3;
                }
            }
            for (int64_t k = fdimension; k < n; k++)
            {
                D(k, k) = c1; // off diagonal part of infinitesimal strain tensor
            }
            break;
        }
        case AnalysisType::EPlaneStrain:
        case AnalysisType::EPlaneStress:
        {
            REAL c1 = 2.0 * fmu;
            REAL c2 = c1 / 2.0;
            REAL c3 = -c2;
            for (int64_t i = 0; i < fdimension; i++)
            {
                D(i, i) = c2;
                for (int64_t j = i+1; j < fdimension; j++)
                {
                    D(i, j) = c3; //volumetric part
                    D(j, i) = c3;
                }
            }
            for (int64_t k = fdimension; k < n; k++)
            {
                D(k, k) = c1; // off diagonal part of infinitesimal strain tensor
            }
            break;
        }
        default:
        {
            std::cout << "Wrong analysis type." << std::endl;
            break;
        }
    }
}

void TPZHybridMixedElasticityUP::ElasticityTensor(TPZFNMatrix<36,REAL>& D)
{
    const int n = fdimension * (fdimension + 1) / 2; //number of independent components of Cauchy stress tensor

    switch (fAnalysisType)
    {
        case AnalysisType::EGeneral:
        case AnalysisType::EPlaneStrain:
        {
            REAL c1 = 2.0 * fmu;
            REAL c2 = c1 + flambda;
            for (int64_t i = 0; i < fdimension; i++)
            {
                D(i, i) = c2;
                for (int64_t j = i+1; j < fdimension; j++)
                {
                    D(i, j) = flambda; //volumetric part
                    D(j, i) = flambda;
                }
            }
            for (int64_t k = fdimension; k < n; k++)
            {
                D(k, k) = c1; // off diagonal part of infinitesimal strain tensor
            }
            break;
        }
        case AnalysisType::EPlaneStress:
        {
            const REAL c1 = fyoung / (1.0 - (fpoisson * fpoisson));
            const REAL c2 = c1 * fpoisson;
            const REAL c3 = c1 * (1.0 - fpoisson);
            
            for (int64_t i = 0; i < fdimension; i++)
            {
                D(i, i) = c1;
                for (int64_t j = i+1; j < fdimension; j++)
                {
                    D(i, j) = c2; //volumetric part
                    D(j, i) = c2;
                }
            }
            for (int64_t k = fdimension; k < n; k++)
            {
                D(k, k) = c3; // off diagonal part of infinitesimal strain tensor
            }
            break;
        }
    }
}

void TPZHybridMixedElasticityUP::StrainTensor(const TPZFNMatrix<10, STATE>& gradU, TPZFNMatrix<6,REAL>& epsilon)
{
    const int n = fdimension * (fdimension + 1) / 2;

    int cont = fdimension - 1;
    for (int i = 0; i < fdimension; i++)
    {
        epsilon(i, 0) = gradU(i, i); // diagonal part of infinitesimal strain tensor
        for (int64_t j = i + 1; j < fdimension; j++)
        {
            epsilon(++cont, 0) = 0.5 * (gradU(i, j) + gradU(j, i)); // off diagonal part of infinitesimal strain tensor
        }
    }
}

void TPZHybridMixedElasticityUP::DeviatoricStressTensor(const TPZFNMatrix<10, STATE>& gradU, TPZFNMatrix<6,REAL>& sigma)
{
    const int n = fdimension * (fdimension + 1) / 2;

    TPZFNMatrix<36, REAL> D(n, n, 0.0);
    DeviatoricElasticityTensor(D);
    
    TPZFNMatrix<6,REAL> strain(n,1,0.0);
    StrainTensor(gradU, strain);

    D.Multiply(strain,sigma);
}

void TPZHybridMixedElasticityUP::StressTensor(const TPZFNMatrix<10, STATE>& gradU, TPZFNMatrix<6,REAL>& sigma)
{
    const int n = fdimension * (fdimension + 1) / 2;

    TPZFNMatrix<36, REAL> D(n, n, 0.0);
    ElasticityTensor(D);
    
    TPZFNMatrix<6,REAL> strain(n,1,0.0);
    StrainTensor(gradU, strain);

    D.Multiply(strain,sigma);
}

void TPZHybridMixedElasticityUP::StressTensor(const TPZFNMatrix<10, STATE>& gradU, TPZFNMatrix<6,REAL>& sigma, REAL pressure)
{
    const int n = fdimension * (fdimension + 1) / 2;

    TPZFNMatrix<36, REAL> D(n, n, 0.0);
    DeviatoricElasticityTensor(D);
    
    TPZFNMatrix<6,REAL> strain(n,1,0.0);
    StrainTensor(gradU, strain);

    D.Multiply(strain,sigma);

    for (int i = 0; i < fdimension; i++)
        sigma(i,0) -= pressure;
}
