//
// Created by Giovane Avancini and Nathan Shauer on 18/12/23.
//

#include <pzfmatrix.h>
#include <TPZBndCondT.h>
#include <pzaxestools.h>
#include <pzlog.h>
#include <pzaxestools.h>

#include "TPZElasticityTH.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.elasticmaterial");
#endif

TPZElasticityTH::TPZElasticityTH() : TBase() {}

TPZElasticityTH::TPZElasticityTH(int matID,
                                 int dimension,
                                 REAL young_modulus,
                                 REAL poisson,
                                 AnalysisType analysisType,
                                 REAL thickness) : TBase(matID),
                                                   fdimension(dimension),
                                                   fyoung(young_modulus),
                                                   fpoisson(poisson),
                                                   fAnalysisType(analysisType),
                                                   fthickness(thickness)
{
    flambda = fyoung * fpoisson / ((1.0 + fpoisson) * (1.0 - 2.0 * fpoisson));
    fmu = 0.5 * fyoung / (1.0 + fpoisson);

    switch (fAnalysisType)
    {
    case AnalysisType::EGeneral:
    {
        fbulk = flambda + 2.0 * fmu / 3.0;
        break;
    }
    case AnalysisType::EPlaneStrain:
    {
        fbulk = flambda + fmu;
        break;
    }
    case AnalysisType::EPlaneStress:
    {
        fbulk = ((fpoisson - 0.5) <= 1.0e-9) ? 3.0 * fmu : fmu * (2.0 * fmu + 3.0 * flambda) / (flambda + 2.0 * fmu); // For full incompressible case. Replacing lambda with E and poisson didnt yield correct results
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

TPZElasticityTH::~TPZElasticityTH() {}

void TPZElasticityTH::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    /*
    This function computes the matrix contribution of each element, that has the following structure:
        |K   G|
        |GT  S|,
        where K is the stiffness matrix, G the gradient operator, GT the divergence operator and S the bulk matrix
    */

    TPZFMatrix<REAL> &PhiU = datavec[EUindex].phi;
    int64_t nShapeU = PhiU.Rows();

    TPZFMatrix<REAL> &PhiP = datavec[EPindex].phi;
    int64_t nShapeP = PhiP.Rows();

    TPZFNMatrix<60, REAL> dphi = datavec[EUindex].dphix;
    auto axes = datavec[EUindex].axes;

    TPZFNMatrix<3, REAL> dphiU(fdimension, nShapeU, 0.0);
    TPZAxesTools<REAL>::Axes2XYZ(dphi, dphiU, axes);

    TPZFNMatrix<3, REAL> SourceTerm(fdimension, 1, 0.0);
    TPZVec<REAL> sourceAux(3);

    if (this->HasForcingFunction())
    {
        this->ForcingFunction()(datavec[EUindex].x, sourceAux);
        for (int64_t i = 0; i < fdimension; i++)
            SourceTerm(i, 0) = sourceAux[i];
    }

    const int n = fdimension * (fdimension + 1) / 2; //number of independent variables using Voight notation

    //We shall divide by 2 the off diagonal part to account for the symmetry when doing the inner product of two tensors
    //and the 1/2 comming from each matrix B that represents the strain tensor, so 1/2*1/2*2 = 1/2
    TPZFNMatrix<36, REAL> D(n, n, 0.0); //Elasticity tensor D
    DeviatoricElasticityTensor(D);
    for (int i = fdimension; i < n; i++)
        D(i,i) *= 0.5;

    // divergence matrix
    TPZFNMatrix<150, STATE> divPhiU(fdimension * nShapeU, 1, 0.0);
    for (int j = 0; j < nShapeU; j++)
        for (int i = 0; i < fdimension; i++)
            divPhiU(j * fdimension + i, 0) = dphiU(i, j);

    // strain matrix B
    TPZFNMatrix<150, STATE> matrixB(n, nShapeU * fdimension, 0.0);

    for (int j = 0; j < nShapeU; j++)
    {
        int cont = fdimension;
        for (int i = 0; i < fdimension; i++)
        {
            matrixB(i, j * fdimension + i) = dphiU(i, j);

            for (int k = i + 1; k < fdimension; k++)
            {
                matrixB(cont, fdimension * j + i) = dphiU(k, j);
                matrixB(cont, fdimension * j + k) = dphiU(i, j);
                cont++;
            }
        }
    }

    TPZFMatrix<REAL> phiU_force(fdimension * nShapeU, fdimension, 0.0);
    for (int j = 0; j < nShapeU; j++)
    {
        for (int i = 0; i < fdimension; i++)
        {
            phiU_force(fdimension * j + i, i) = PhiU(j);
        }
    }

    // body forces contribution
    ef.AddContribution(0, 0, phiU_force, false, SourceTerm, false, weight);

    //Stiffness - Matrix K contribution
    TPZFNMatrix<150, REAL> aux;
    D.Multiply(matrixB, aux);
    
    REAL factor = weight;
    ek.AddContribution(0, 0, matrixB, true, aux, false, factor);

    // Gradient - matrix G contribution
    factor = -weight;
    ek.AddContribution(0, fdimension*nShapeU, divPhiU, false, PhiP, true, factor);

    // Divergence - Matrix GT contribution
    ek.AddContribution(fdimension*nShapeU, 0, PhiP, false, divPhiU, true, factor);

    //Bulk Matrix S contribution
    factor = -(1.0 / fbulk) * weight;
    ek.AddContribution(fdimension*nShapeU, fdimension*nShapeU, PhiP, false, PhiP, true, factor);

#ifdef PZ_LOG
    if(logger.isDebugEnabled()){
        std::stringstream sout;
        ek.Print("ek", sout, EMathematicaInput);
        ef.Print("ef", sout, EMathematicaInput);
        sout << std::endl << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
}

void TPZElasticityTH::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc)
{
    int index = (bc.Type() == 0 || bc.Type() == 2) ? EUindex : EPindex;

    TPZFNMatrix<150, REAL> PhiU = datavec[EUindex].phi;
    TPZFMatrix<REAL> &PhiP = datavec[EPindex].phi;

    int64_t nShapeU = PhiU.Rows();
    int64_t nShapeP = PhiP.Rows();

    TPZFNMatrix<20, STATE> val1(3, 3, 0.0);
    TPZManVector<STATE, 3> val2(3, 0.0);
    const auto &BIGNUMBER  = TPZMaterial::fBigNumber;

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
        case 0 : // Dirichlet condition at x and y direction
        {
            for(int i = 0 ; i < nShapeU; i++)
            {
                for (int j = 0; j < fdimension; j++)
                {
                    ef(fdimension*i+j, 0) += BIGNUMBER * val2[j] * PhiU(i,0) * weight;
                    for (int k = 0; k < nShapeU; k++)
                    {
                        ek(fdimension*i+j, fdimension*k+j) += BIGNUMBER * PhiU(i,0) *PhiU(k,0) * weight;
                    }
                }
            }
            break;
        }
        case 1 : // Neumann condition
        {
            for (int i = 0; i < nShapeU; i++) 
            {
                for (int j = 0; j < fdimension; j++)
                {
                    ef(fdimension*i, 0) += val2[j] * PhiU(i, 0) * weight;
                }
            }
            break;
        }
            
        case 2 : // Mixed Condition
        {
            DebugStop(); // Implement me            
            break;
        }
        
        case 3 : // Dirichlet X condition
        {
            for(int i = 0 ; i < nShapeU; i++)
            {
                ef(fdimension*i, 0) += BIGNUMBER * val2[0] * PhiU(i,0) * weight; // forced x displacement
                for (int j = 0 ; j < nShapeU; j++)
                {
                    ek(fdimension*i,fdimension*j) += BIGNUMBER * PhiU(i,0) *PhiU(j,0) * weight;
                }
            }
            break;
        }

        case 4 : // Dirichlet Y condition
        {
            for(int i = 0 ; i < nShapeU; i++)
            {
                ef(fdimension*i+1, 0) += BIGNUMBER * val2[1] * PhiU(i,0) * weight; // forced y displacement
                for (int j = 0 ; j < nShapeU; j++)
                {
                    ek(fdimension*i+1,fdimension*j+1) += BIGNUMBER * PhiU(i,0) *PhiU(j,0) * weight;
                }
            }
            break;
        }

        case 5 : // Dirichlet Z condition
        {   
            if (fdimension != 3)
                DebugStop();

            for(int i = 0 ; i < nShapeU; i++)
            {
                ef(fdimension*i+2, 0) += BIGNUMBER * val2[2] * PhiU(i,0) * weight; // forced z displacement
                for (int j = 0 ; j < nShapeU; j++)
                {
                    ek(fdimension*i+2,fdimension*j+2) += BIGNUMBER * PhiU(i,0) *PhiU(j,0) * weight;
                }
            }
            break;
        }          
		case 6: // stressField Neumann condition
			for(int in = 0; in < fdimension; in++)
                for(int jn = 0; jn < fdimension; jn++)
                    val2[in] += - val1(in,jn) * datavec[EUindex].normal[jn];
			// The normal vector points towards the neighbour. The negative sign is there to 
			// reflect the outward normal vector.
			for(int in = 0 ; in < nShapeU; in++) {
                for(int idim = 0; idim < fdimension; idim++){
                    ef(fdimension*in+idim,0) += val2[idim] * PhiU(in,0) * weight;    
                }				
			}
			break;        

        default:
        {
            std::cout << "ERROR: BOUNDARY NOT IMPLEMENTED" << std::endl;
            DebugStop();
            break;
        }
    }
}

int TPZElasticityTH::VariableIndex(const std::string &name) const
{

    if (!strcmp("Pressure", name.c_str()))
        return EPressure;
    if (!strcmp("Displacement", name.c_str()))
        return EDisplacement;
    if (!strcmp("Force", name.c_str()))
        return EForce;
    if (!strcmp("Stress", name.c_str()))
        return EStress;
    if (!strcmp("Strain", name.c_str()))
        return EStrain;
    if (!strcmp("VonMises", name.c_str()))
        return EVonMises;

    std::cout << "\n\nVar index not implemented\n\n";
    DebugStop();

    return 0;
}

int TPZElasticityTH::NSolutionVariables(int var) const
{

    int aux;
    switch (var)
    {
    case EPressure: // pressure  [scalar]
    case EVonMises: // VonMises
        aux = 1;
        break;
    case EDisplacement: // displacement [vector]
    case EForce:        // external force [vector]
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

void TPZElasticityTH::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &Solout)
{
    Solout.Fill(0.0);
    TPZManVector<STATE, 3> u_h = datavec[EUindex].sol[0];
    TPZManVector<STATE, 3> p_h = datavec[EPindex].sol[0];

    TPZFNMatrix<10, STATE> gradU_xsi = datavec[EUindex].dsol[0];
    auto axes = datavec[EUindex].axes;
    TPZFNMatrix<9, REAL> gradU(fdimension, fdimension, 0.0);
    TPZAxesTools<REAL>::Axes2XYZ(gradU_xsi, gradU, axes);

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
        strain(2, 2) = (fpoisson == 0.5) ? -(strain(0, 0) + strain(1, 1)) : -(flambda / (flambda + 2.0 * fmu)) * (strain(0, 0) + strain(1, 1));

    Solout.Resize(NSolutionVariables(var));

    switch (var)
    {

    case EPressure:
    {
        Solout[0] = p_h[0];
        break;
    }

    case EDisplacement:
    {
        Solout[0] = u_h[0]; // Vx
        Solout[1] = u_h[1]; // Vy
        if (Dimension() == 3)
            Solout[2] = u_h[2]; // Vz

        if (fAnalysisType == AnalysisType::EPlaneStress)
            Solout[2] = strain(2, 2) * fthickness;
        break;
    }

    case EForce:
    {
        TPZVec<STATE> f(3, 0.);

        if (this->HasForcingFunction())
        {
            this->ForcingFunction()(datavec[EUindex].x, f);
        }

        Solout[0] = f[0];
        Solout[1] = f[1];
        if (Dimension() == 3)
            Solout[2] = f[2];
        break;
    }
    case EStress:
    {
        TPZFNMatrix<6, STATE> sigmavoight(n, 1, 0.0);
        DeviatoricStressTensor(gradU, sigmavoight);

        TPZFNMatrix<9, STATE> sigma(3, 3, 0.0);
        int cont = fdimension - 1;
        for (int i = 0; i < fdimension; i++)
        {
            sigma(i, i) = sigmavoight(i, 0) - p_h[0];
            for (int j = i + 1; j < fdimension; j++)
            {
                sigma(i, j) = sigmavoight(++cont, 0);
                sigma(j, i) = sigmavoight(cont, 0);
            }
        }

        if (fAnalysisType == AnalysisType::EPlaneStrain)
            sigma(2, 2) = fthickness * fpoisson * (sigma(0, 0) + sigma(1, 1));

        for (int i = 0; i < Dimension(); i++)
        {
            for (int j = 0; j < Dimension(); j++)
            {
                Solout[i * 3 + j] = sigma(i, j);
            }
        }
        break;
    }
    case EStrain:
    {
        for (int i = 0; i < Dimension(); i++)
        {
            for (int j = 0; j < Dimension(); j++)
            {
                Solout[i * Dimension() + j] = strain(i, j);
            }
        }
        break;
    }
    case EVonMises:
    {
        TPZManVector<STATE, 3> PrincipalStress(3);
        TPZFNMatrix<6, STATE> sigmavoight(n, 1, 0.0);
        StressTensor(gradU, sigmavoight, p_h[0]);
        TPZFNMatrix<9, STATE> sigma(3, 3, 0.0);
        int cont = fdimension - 1;
        for (int i = 0; i < fdimension; i++)
        {
            sigma(i, i) = sigmavoight(i, 0);
            for (int j = i + 1; j < fdimension; j++)
            {
                sigma(i, j) = sigmavoight(++cont, 0);
                sigma(j, i) = sigmavoight(cont, 0);
            }
        }
        if (fAnalysisType == AnalysisType::EPlaneStrain)
            sigma(2, 2) = fthickness * fpoisson * (sigma(0, 0) + sigma(1, 1));

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
        Solout[0] = (PrincipalStress[0] - PrincipalStress[1]) * (PrincipalStress[0] - PrincipalStress[1]) + (PrincipalStress[1] - PrincipalStress[2]) * (PrincipalStress[1] - PrincipalStress[2]) + (PrincipalStress[2] - PrincipalStress[0]) * (PrincipalStress[2] - PrincipalStress[0]);
        Solout[0] = sqrt(0.5 * Solout[0]);
        break;
    }

    default:
    {
        std::cout << "\n\nVar index not implemented\n\n";
        DebugStop();
    }
    }
}

void TPZElasticityTH::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const
{
    int64_t ndata = datavec.size();
    for (int idata = 0; idata < ndata; idata++)
    {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        datavec[idata].fNeedsHSize = true;
        datavec[idata].fNeedsNormal = true;
    }
}

void TPZElasticityTH::FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const
{
    datavec[EUindex].fNeedsSol = false;
    datavec[EPindex].fNeedsSol = false;
    datavec[EUindex].fNeedsNormal = true;
    datavec[EPindex].fNeedsNormal = true;
}

void TPZElasticityTH::Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors)
{
    // 0: L2 p, 1: L2 p_ex, 2: L2 u, 3: L2 u_ex, 4: L2 divu, 5: L2 divu_ex, 6: L2 sigma, 7: L2 sigma_Ex
    if (!HasExactSol())
        DebugStop();

    errors.Resize(NEvalErrors());

    TPZManVector<STATE, 4> sol_exact(4);
    TPZFNMatrix<9, STATE> gradsol_exact(3, 3);

    // Getting the exact solution for velocity, pressure and velocity gradient
    fExactSol(data[EUindex].x, sol_exact, gradsol_exact);

    // Getting the numeric solution for velocity, pressure and velocity gradient
    TPZManVector<STATE> u_h(3, 0.0);
    TPZManVector<STATE> p_h(1, 0.0);
    TPZFNMatrix<10, STATE> gradU_xsi = data[EUindex].dsol[0];
    auto axes = data[EUindex].axes;
    TPZFNMatrix<9, REAL> gradu_h(fdimension, fdimension, 0.0);
    TPZAxesTools<REAL>::Axes2XYZ(gradU_xsi, gradu_h, axes);

    this->Solution(data, VariableIndex("Displacement"), u_h);
    this->Solution(data, VariableIndex("Pressure"), p_h);

    const int n = fdimension * (fdimension + 1) / 2;
    TPZFNMatrix<6, REAL> sigma_exact(n, 1), sigma_h(n, 1);
    StressTensor(gradsol_exact, sigma_exact);
    StressTensor(gradu_h, sigma_h, p_h[0]);

    STATE p_exact = 0.0; //If poisson == 0.5, this will break, maybe create a mor general Elasticity Analytic Solution class to deal with this issue
    for (int i = 0; i < fdimension; i++)
    {
        p_exact -= sigma_exact(i,i);
    }
    p_exact *= 1./fdimension;

    STATE div_exact = 0.0, div_h = 0.0;
    for (int i = 0; i < fdimension; i++)
    {
        div_exact += gradsol_exact(i, i);
        div_h += gradu_h(i, i);
    }

    STATE diffu, diffp, diffdiv;

    diffp = p_h[0] - p_exact;
    errors[0] = diffp * diffp;
    errors[1] = p_exact * p_exact;

    errors[2] = 0.0;
    errors[3] = 0.0;
    for (int i = 0; i < fdimension; i++)
    {
        diffu = u_h[i] - sol_exact[i];
        errors[2] += diffu * diffu;
        errors[3] += sol_exact[i] * sol_exact[i];
    }

    diffdiv = div_h - div_exact;
    errors[4] = diffdiv * diffdiv;
    errors[5] = div_exact * div_exact;

    errors[6] = 0.0;
    errors[7] = 0.0;
    for (int i = 0; i < fdimension; i++)
    {
        const STATE diffsig = sigma_h(i, 0) - sigma_exact(i, 0);
        errors[6] += diffsig * diffsig;
        errors[7] += sigma_exact(i, 0) * sigma_exact(i, 0);
    }
    for (int i = fdimension; i < n; i++)
    {
        const STATE diffsig = sigma_h(i, 0) - sigma_exact(i, 0);
        errors[6] += 2. * diffsig * diffsig;
        errors[7] += 2. * sigma_exact(i, 0) * sigma_exact(i, 0);
    }
}

void TPZElasticityTH::DeviatoricElasticityTensor(TPZFNMatrix<36, REAL> &D)
{
    const int n = fdimension * (fdimension + 1) / 2; // number of independent components of Cauchy stress tensor

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
            for (int64_t j = i + 1; j < fdimension; j++)
            {
                D(i, j) = c3; // volumetric part
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
            for (int64_t j = i + 1; j < fdimension; j++)
            {
                D(i, j) = c3; // volumetric part
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

void TPZElasticityTH::ElasticityTensor(TPZFNMatrix<36, REAL> &D)
{
    const int n = fdimension * (fdimension + 1) / 2; // number of independent components of Cauchy stress tensor

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
            for (int64_t j = i + 1; j < fdimension; j++)
            {
                D(i, j) = flambda; // volumetric part
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
            for (int64_t j = i + 1; j < fdimension; j++)
            {
                D(i, j) = c2; // volumetric part
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

void TPZElasticityTH::StrainTensor(const TPZFNMatrix<10, STATE> &gradU, TPZFNMatrix<6, REAL> &epsilon)
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

void TPZElasticityTH::DeviatoricStressTensor(const TPZFNMatrix<10, STATE> &gradU, TPZFNMatrix<6, REAL> &sigma)
{
    const int n = fdimension * (fdimension + 1) / 2;

    TPZFNMatrix<36, REAL> D(n, n, 0.0);
    DeviatoricElasticityTensor(D);

    TPZFNMatrix<6, REAL> strain(n, 1, 0.0);
    StrainTensor(gradU, strain);

    D.Multiply(strain, sigma);
}

void TPZElasticityTH::StressTensor(const TPZFNMatrix<10, STATE> &gradU, TPZFNMatrix<6, REAL> &sigma)
{
    const int n = fdimension * (fdimension + 1) / 2;

    TPZFNMatrix<36, REAL> D(n, n, 0.0);
    ElasticityTensor(D);

    TPZFNMatrix<6, REAL> strain(n, 1, 0.0);
    StrainTensor(gradU, strain);

    D.Multiply(strain, sigma);
}

void TPZElasticityTH::StressTensor(const TPZFNMatrix<10, STATE> &gradU, TPZFNMatrix<6, REAL> &sigma, REAL pressure)
{
    const int n = fdimension * (fdimension + 1) / 2;

    TPZFNMatrix<36, REAL> D(n, n, 0.0);
    DeviatoricElasticityTensor(D);

    TPZFNMatrix<6, REAL> strain(n, 1, 0.0);
    StrainTensor(gradU, strain);

    D.Multiply(strain, sigma);

    for (int i = 0; i < fdimension; i++)
        sigma(i, 0) -= pressure;
}
