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
    TPZFMatrix<REAL>& PhiU = datavec[EUindex].phi;
    int64_t nShapeU = PhiU.Rows(); // for p=2 should be 9?

    TPZFMatrix<REAL>& PhiP = datavec[EPindex].phi;
    int64_t nShapeP = PhiP.Rows(); // number of pressure H1 shape functions

    const int n = fdimension * (fdimension + 1) / 2; // number of independent components of stress tensor in voight notation
    // TPZFNMatrix<150, REAL> Strain(n, nShapeU, 0.0);  // Using voight notation

    TPZFNMatrix<60, REAL> dphi = datavec[EUindex].dphix;
    auto axes = datavec[EUindex].axes;

    TPZFNMatrix<60, REAL> dPhiU(fdimension, nShapeU, 0.);
    TPZAxesTools<REAL>::Axes2XYZ(dphi, dPhiU, axes);

    TPZFNMatrix<3, REAL> SourceTerm(fdimension, 1, 0.0);
    TPZVec<REAL> sourceAux(3);
    if (this->HasForcingFunction())
    {
        this->ForcingFunction()(datavec[EUindex].x, sourceAux);
        for (int64_t i = 0; i < fdimension; i++)
        {
            SourceTerm(i, 0) = sourceAux[i];
        }
    }

    // TPZFNMatrix<36, REAL> D(n, n, 0.0);
    // DeviatoricElasticityTensor(D);
    // for (int i = fdimension; i < n; i++) // The terms related to the off diagonal part of strain tensor is multiplied by 2 to account for its symmetry
    //     D(i, i) *= 2.0;

    // Body Forces contribution
    // ef.AddContribution(0, 0, PhiU, true, SourceTerm, false, weight);

    // if (datavec.size() > 2) // Static condensation in incompressibility regime
    // {
    //     TPZFMatrix<REAL> &PhiUM = datavec[EVMindex].phi;
    //     TPZFMatrix<REAL> &phipM = datavec[EPMindex].phi;

    // 2mu e(u) x e(v) diagonal term
    const int dim = Dimension();
    for (int i = 0; i < nShapeU; i++) {
        const STATE dvxdx = dPhiU(0, i), dvxdy = dPhiU(1, i);
        const STATE dvydx = dPhiU(0, i), dvydy = dPhiU(1, i);
        // Load vector test function times source term        
        ef(i*dim,0) += SourceTerm(0, 0) * PhiU(i, 0) * weight;
        ef(i*dim+1,0) += SourceTerm(1, 0) * PhiU(i, 0) * weight;

        for (int j = 0; j < nShapeU; j++) {
            // ek(i*dim,j*dim) += 2.0 * fmu * (dPhiU(0, i) * dPhiU(0, j) + 0.5 * dPhiU(1, i) * dPhiU(1, j)) * weight;
            // ek(i*dim,j*dim+1) += 2.0 * fmu * (dPhiU(1, i) * dPhiU(0, j)) * weight;
            // ek(j*dim,i*dim+1) += 2.0 * fmu * (dPhiU(1, i) * dPhiU(0, j)) * weight;
            // ek(i*dim+1,j*dim+1) += 2.0 * fmu * (dPhiU(1, i) * dPhiU(1, j) + 0.5 * dPhiU(0, i) * dPhiU(0, j)) * weight;
            const STATE duxdx = dPhiU(0, j), duxdy = dPhiU(1, j);
            const STATE duydx = dPhiU(0, j), duydy = dPhiU(1, j);           
            
            ek(i*dim,j*dim) += fmu*(duxdy*dvxdy + duxdx*dvxdx) * weight;
            ek(i*dim,j*dim+1) += fmu*(dvxdy*duydx - dvxdx*duydy) * weight;
            ek(i*dim+1,j*dim) += fmu*(dvydx*duxdy - dvydy*duxdx) * weight;
            ek(i*dim+1,j*dim+1) += fmu*(duydx*dvydx + duydy*dvydy) * weight;            
            // ek(i*dim,j*dim) += fmu*(duxdy*dvxdy + duxdx*dvxdx) * weight;
            // ek(i*dim,j*dim+1) += fmu*(dvxdy*duydx - dvxdx*duydy) * weight;
            // ek(j*dim,i*dim+1) += fmu*(dvxdy*duydx - dvxdx*duydy) * weight;
            // ek(i*dim+1,j*dim+1) += fmu*(duydx*dvydx + duydy*dvydy) * weight;            

        }
        for (int j = 0; j < nShapeP; j++) {
            const STATE contrib1 = -PhiP(j, 0) * dPhiU(0, i) * weight,
                contrib2 = -PhiP(j, 0) * dPhiU(1, i) * weight;
            ek(i*dim,j+dim*nShapeU) += contrib1;
            ek(i*dim+1,j+dim*nShapeU) += contrib2;
            
            // Complete the symmetric part of ek
            ek(j+dim*nShapeU,i*dim) += contrib1;
            ek(j+dim*nShapeU,i*dim+1) += contrib2;
        }
    }
    

    // q x p / kappa diagonal term 
    for (int i = 0; i < nShapeP; i++) {
      for (int j = 0; j < nShapeP; j++) {
        ek(dim*nShapeU + i, dim*nShapeU + j) += -PhiP(i, 0) * PhiP(j, 0) / fbulk * weight;
      }
    }

    //     // Pressure and distributed displacement
    //     for (int j = 0; j < nShapeP; j++)
    //     {
    //         ek(nShapeU + nShapeP, nShapeU + j) += PhiP(j, 0) * PhiUM(0, 0) * weight;
    //         ek(nShapeU + j, nShapeU + nShapeP) += PhiP(j, 0) * PhiUM(0, 0) * weight;
    //     }

    //     // Injection and average-pressure
    //     ek(nShapeU + nShapeP + 1, nShapeU + nShapeP) += PhiUM(0, 0) * phipM(0, 0) * weight;
    //     ek(nShapeU + nShapeP, nShapeU + nShapeP + 1) += PhiUM(0, 0) * phipM(0, 0) * weight;
    // }

#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        ek.Print("ek", sout, EMathematicaInput);
        ef.Print("ef", sout, EMathematicaInput);
        sout << std::endl
             << std::endl;
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
    switch (bc.Type()) {
        case 0 :			// Dirichlet condition
        {
            for(int in = 0 ; in < nShapeU; in++) {
                ef(2*in,0)   += BIGNUMBER * val2[0] * PhiU(in,0) * weight;        // forced v2 displacement
                ef(2*in+1,0) += BIGNUMBER * val2[1] * PhiU(in,0) * weight;        // forced v2 displacement
                for (int jn = 0 ; jn < nShapeU; jn++) {
                    ek(2*in,2*jn)     += BIGNUMBER * PhiU(in,0) *PhiU(jn,0) * weight;
                    ek(2*in+1,2*jn+1) += BIGNUMBER * PhiU(in,0) *PhiU(jn,0) * weight;
                }
            }
        }
            break;
            
        case 1 :		// Neumann condition
        {
            for (int in = 0; in < nShapeU; in++) 
            {
              ef(2 * in, 0) += val2[0] * PhiU(in, 0) * weight;      // force in x direction
              ef(2 * in + 1, 0) += val2[1] * PhiU(in, 0) * weight;  // force in y direction
            }
        }
            break;
            
        case 2 :		// Mixed Condition
        {
            DebugStop(); // Implement me            
            break;
        }
        
        case 3 :			// Dirichlet X condition
        {
            for(int in = 0 ; in < nShapeU; in++) {
                ef(2*in,0)   += BIGNUMBER * val2[0] * PhiU(in,0) * weight;        // forced x displacement
                for (int jn = 0 ; jn < nShapeU; jn++) {
                    ek(2*in,2*jn)     += BIGNUMBER * PhiU(in,0) *PhiU(jn,0) * weight;
                }
            }
        }
            break;

        case 4 :			// Dirichlet Y condition
        {
            for(int in = 0 ; in < nShapeU; in++) {
                ef(2*in+1,0) += BIGNUMBER * val2[1] * PhiU(in,0) * weight;        // forced y displacement
                for (int jn = 0 ; jn < nShapeU; jn++) {
                    ek(2*in+1,2*jn+1) += BIGNUMBER * PhiU(in,0) *PhiU(jn,0) * weight;
                }
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
    datavec[0].fNeedsDeformedDirectionsFad = true;
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
    REAL p_exact = -data[EUindex].x[1] * data[EUindex].x[2] / 3.; // Just for computing Bishop beam when poisson is 0.5. Remember to delete later.

    // Getting the numeric solution for velocity, pressure and velocity gradient
    TPZManVector<STATE> u_h(3, 0.0);
    TPZManVector<STATE> p_h(1, 0.0);
    TPZFNMatrix<10, STATE> gradv_h = data[EUindex].dsol[0];

    this->Solution(data, VariableIndex("Displacement"), u_h);
    this->Solution(data, VariableIndex("Pressure"), p_h);

    STATE diffv, diffp, diffdiv;

    diffp = p_h[0] - p_exact;
    errors[0] = diffp * diffp;
    errors[1] = p_exact * p_exact;

    errors[2] = 0.0;
    errors[3] = 0.0;
    for (int i = 0; i < fdimension; i++)
    {
        diffv = u_h[i] - sol_exact[i];
        errors[2] += diffv * diffv;
        errors[3] += sol_exact[i] * sol_exact[i];
    }

    STATE div_exact = 0.0, div_h = 0.0;
    for (int i = 0; i < fdimension; i++)
    {
        div_exact += gradsol_exact(i, i);
        div_h += gradv_h(i, i);
    }

    diffdiv = div_h - div_exact;
    errors[4] = diffdiv * diffdiv;
    errors[5] = div_exact * div_exact;

    const int n = fdimension * (fdimension + 1) / 2;
    TPZFNMatrix<6, REAL> sigma_exact(n, 1), sigma_h(n, 1);
    DeviatoricStressTensor(gradsol_exact, sigma_exact); // Just for Bishop beam. remember to delete later
    for (int i = 0; i < fdimension; i++)
        sigma_exact(i, 0) -= p_exact;
    StressTensor(gradv_h, sigma_h, p_h[0]);

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
        errors[7] += sigma_exact(i, 0) * sigma_exact(i, 0);
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
