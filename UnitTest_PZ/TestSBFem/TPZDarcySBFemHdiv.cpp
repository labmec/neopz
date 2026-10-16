/**
 * @file
 * @brief Contains the methods of the TPZHybridPoissonCollapsed class (multiphysics environment)
 * @author Karolinne Coelho
 * @date 2020/11/17
 */

#include "TPZDarcySBFemHdiv.h"
#include "pzlog.h"
#include "TPZBndCondT.h"
#include "pzfmatrix.h"
#include "pzaxestools.h"
#include "pzinterpolationspace.h"


#include <iostream>

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.mixedpoisson.data"));
static LoggerPtr logerror(Logger::getLogger("pz.mixedpoisson.error"));
#endif

TPZHybridPoissonCollapsed::TPZHybridPoissonCollapsed(): TPZMixedDarcyFlow() {}

TPZHybridPoissonCollapsed::TPZHybridPoissonCollapsed(int matid, int dim): TPZMixedDarcyFlow(matid,dim){}

TPZHybridPoissonCollapsed::TPZHybridPoissonCollapsed(const TPZHybridPoissonCollapsed &cp)
{
    new TPZHybridPoissonCollapsed(cp.Id(),cp.Dimension());
}

TPZHybridPoissonCollapsed::~TPZHybridPoissonCollapsed() = default;

TPZHybridPoissonCollapsed &TPZHybridPoissonCollapsed::operator=(const TPZHybridPoissonCollapsed &copy) {
    TPZHybridPoissonCollapsed::operator=(copy);
    return *this;
}

void TPZHybridPoissonCollapsed::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                                   TPZFMatrix<STATE> &ef)
{
    
    STATE force = 0;
    if (fForcingFunction) {
        TPZManVector<STATE> res(1);
        fForcingFunction(datavec[1].x, res);
        force = res[0];
    }

    const REAL perm = fPermeability;
    const REAL inv_perm = 1 / perm;

    // Setting the phis
    TPZFMatrix<REAL> &phiQ = datavec[0].fH1.fPhi;
    TPZFMatrix<REAL> &phip = datavec[1].fH1.fPhi;
    TPZFMatrix<REAL> &dphiQ = datavec[0].dphix;
    TPZFMatrix<REAL> &dphiP = datavec[1].dphix;
    TPZFNMatrix<9, REAL> dphiPXY(3, dphiP.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP, dphiPXY, datavec[1].axes);

    REAL &faceSize = datavec[0].HSize;

    int phrq, phrp;
    phrp = phip.Rows();
    phrq = datavec[0].fVecShapeIndex.NElements();

    //Calculate the matrix contribution for flux. Matrix A
    for (int iq = 0; iq < phrq; iq++) {
        //ef(iq, 0) += 0.;
        int ivecind = datavec[0].fVecShapeIndex[iq].first;
        int ishapeind = datavec[0].fVecShapeIndex[iq].second;
        TPZFNMatrix<3, REAL> ivec(3, 1, 0.);
        for (int id = 0; id < 3; id++) {
            ivec(id, 0) = datavec[0].fDeformedDirections(id, ivecind);
        }

        TPZFNMatrix<3, REAL> ivecZ(3, 1, 0.);
        TPZFNMatrix<3, REAL> jvecZ(3, 1, 0.);
        for (int jq = 0; jq < phrq; jq++) {
            TPZFNMatrix<3, REAL> jvec(3, 1, 0.);
            int jvecind = datavec[0].fVecShapeIndex[jq].first;
            int jshapeind = datavec[0].fVecShapeIndex[jq].second;

            //dot product between Kinv[u]v
            jvecZ.Zero();
            for (int id = 0; id < 3; id++) {
                jvecZ(id, 0) += inv_perm * datavec[0].fDeformedDirections(id, jvecind);
            }
            REAL prod1 = ivec(0, 0) * jvecZ(0, 0) + ivec(1, 0) * jvecZ(1, 0) + ivec(2, 0) * jvecZ(2, 0);
            ek(iq, jq) += weight * phiQ(ishapeind, 0) * phiQ(jshapeind, 0) * prod1;
        }
    }


    // Coupling terms between flux and pressure. Matrix B
    for (int iq = 0; iq < phrq; iq++)
    {
        for (int jp = 0; jp < phrp; jp++) {

            REAL fact = (-1.) * weight * phip(jp, 0) * datavec[0].divphi(iq,0);
            // Matrix B
            ek(iq, phrq + jp) += fact;

            // Matrix B^T
            ek(phrq + jp, iq) += fact;
        }
    }

    // source term related to the pressure equation
    for (int ip = 0; ip < phrp; ip++) {
        ef(phrq + ip, 0) += (-1.) * weight * force * phip(ip, 0);
    }

    std::ofstream ssout("localmatrices.txt");
    ek.Print("ek = ", ssout, EMathematicaInput);
}

void TPZHybridPoissonCollapsed::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                      TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc)
{
    int dim = Dimension();

    TPZFMatrix<REAL> &phiP = datavec[1].fH1.fPhi;
    int phrp = phiP.Rows();

    REAL v2 = bc.Val2()[0];
    REAL v1 = bc.Val1()(0, 0);
    REAL u_D = 0;
    REAL normflux = 0.;

    if (bc.HasForcingFunctionBC()) {
        TPZManVector<STATE> res(3);
        TPZFNMatrix<9, STATE> gradu(3, 1);
        bc.ForcingFunctionBC()(datavec[1].x, res, gradu);
        v1 = gradu(0,0);
        v2 = res[0];
    }

    TPZMaterial::fBigNumber = 1.e12;

	switch (bc.Type()) {
		case 0 :        // Dirichlet condition
            for (int iq = 0; iq < phrp; iq++) {
                //the contribution of the Dirichlet boundary condition appears in the flow equation
                ef(iq, 0) += (1.) * v2 * phiP(iq, 0) * weight * TPZMaterial::fBigNumber;
                for (int jq = 0; jq < phrp; jq++){
                    ek(iq,jq) += weight * TPZMaterial::fBigNumber * phiP(iq,0) * phiP(jq,0);
                }
            }
            break;

        case 1 :            // Neumann condition
            for (int iq = 0; iq < phrp; iq++) {
                ef(iq, 0) += v2 * phiP(iq, 0) * weight;
            }
            break;
        case 2 :			// mixed condition
            DebugStop();
            
            break;
	}
}

void TPZHybridPoissonCollapsed::Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors)
{
    
    /**
     * datavec[0]= Flux
     * datavec[1]= Pressure
     *
     * Errors:
     * [0] L2 for pressure
     * [1] L2 for flux
     * [2] L2 for div(flux)
     * [3] Grad pressure (Semi H1)
     * [4] Hdiv norm
    **/
    errors.Resize(NEvalErrors());
    errors.Fill(0.0);

    TPZManVector<STATE, 3> fluxfem(3), pressurefem(1);
    fluxfem = data[0].sol[0];
    STATE divsigmafem = data[0].divsol[0][0];

    auto dsol = data[1].dsol;

    TPZVec<STATE> u_exact(1, 0);
    TPZFMatrix<STATE> du_exact(3, 1, 0);
    TPZVec<STATE> divsigma(1);
    if (this->fExactSol) {
        this->fExactSol(data[1].x, u_exact, du_exact);
    }
    if (this->fForcingFunction) {
        this->fForcingFunction(data[1].x, divsigma);
    }

    REAL residual = (divsigma[0] - divsigmafem) * (divsigma[0] - divsigmafem);
    pressurefem[0] = data[1].sol[0][0];
    
    REAL perm = 1.;
    REAL inv_perm = 1 / perm;

    TPZManVector<STATE, 3> gradpressurefem(3, 0.);
    this->Solution(data, VariableIndex("GradPressure"), gradpressurefem);

    TPZManVector<STATE, 3> fluxexactneg(3, 0);
    TPZManVector<STATE, 3> gradpressure(3, 0);
    for (int i = 0; i < 3; i++) {
        gradpressure[i] = du_exact[i];
        fluxexactneg[i] = -perm * gradpressure[i];
    }

    REAL L2flux = 0., L2grad = 0.;
    for (int i = 0; i < Dimension(); i++) {
        L2flux += (fluxfem[i] + fluxexactneg[i]) * inv_perm * (fluxfem[i] + fluxexactneg[i]);
        L2grad += (du_exact[i] - gradpressurefem[i]) * (du_exact[i] - gradpressurefem[i]);
    }
    errors[0] = (pressurefem[0] - u_exact[0]) * (pressurefem[0] - u_exact[0]);//L2 error for pressure
    errors[1] = L2flux;//L2 error for flux
    errors[2] = residual;//L2 for div
    errors[3] = L2grad;
    errors[4] = L2flux + residual;
    if (isnan(errors[3]))
    {
        DebugStop();
    }
}
