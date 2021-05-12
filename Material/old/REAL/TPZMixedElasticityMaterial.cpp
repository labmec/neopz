/**
 * @file
 * @brief Contains implementations of the TPZMixedElasticityMaterial methods.
 */

#include "TPZMixedElasticityMaterial.h"
#include "pzelmat.h"
#include "pzbndcond.h"
#include "pzaxestools.h"
#include "TPZMatWithMem.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include <math.h>

#include "pzlog.h"
//#include "meshgen.h"
#include "TPZAnalyticSolution.h"
#ifdef PZ_LOG
static TPZLogger logdata("pz.mixedelasticity");
#endif

#include <fstream>
using namespace std;

TPZMixedElasticityMaterial::TPZMixedElasticityMaterial() : TPZMaterial(0) {
    fE_const = -1.; // Young modulus
    fnu_const = -1.; // Poisson coefficient
    fForce[0] = 0.; // X component of the body force
    fForce[1] = 0.; // Y component of the body force
    fForce[2] = 0.; // Z component of the body force - not used for this class
    flambda_const = 0.;
    fmu_const = 0.;

    fPlaneStress = 0;
    fPostProcIndex = 0;
    fMatrixA = 0.;
}

TPZMixedElasticityMaterial::TPZMixedElasticityMaterial(int id) : TPZMaterial(id) {
    fE_const = -1.; // Young modulus
    fnu_const = -1.; // Poisson coefficient
    fForce[0] = 0.; // X component of the body force
    fForce[1] = 0.; // Y component of the body force
    fForce[2] = 0.; // Z component of the body force - not used for this class
    flambda_const = 0.;
    fmu_const = 0.;

    fPlaneStress = 0;
    fPostProcIndex = 0;
    fMatrixA = 0.;
}

TPZMixedElasticityMaterial::TPZMixedElasticityMaterial(int num, REAL E, REAL nu, REAL fx, REAL fy, int planestress, int fdimension) : TPZMaterial(num), fDimension(fdimension) {
    this->SetElasticity(E, nu);
    fForce[0] = fx; // X component of the body force
    fForce[1] = fy; // Y component of the body force
    fPlaneStress = planestress;
    fPostProcIndex = 0;
}

TPZMixedElasticityMaterial::~TPZMixedElasticityMaterial() {
}

void TPZMixedElasticityMaterial::FillDataRequirements(TPZVec<TPZMaterialData > &datavec) {
    int nref = datavec.size();
    for (int i = 0; i < nref; i++) {
        datavec[i].SetAllRequirements(false);
        datavec[i].fNeedsNeighborSol = false;
        datavec[i].fNeedsNeighborCenter = false;
        datavec[i].fNeedsNormal = false;
    }
}

int TPZMixedElasticityMaterial::NStateVariables() const {
    return 2;
}

////////////////////////////////////////////////////////////////////

// Divergence on deformed element

void TPZMixedElasticityMaterial::ComputeDivergenceOnDeformed(TPZVec<TPZMaterialData> &datavec, TPZFMatrix<STATE> &DivergenceofPhi) {
    //itapopo conferir esse método. Foi copiado do TPZDarcyFlow3D
    int sigmaBlock = 0;

    // Getting test and basis functions
    //TPZFMatrix<REAL> phiSigmaH1 = datavec[sigmaBlock].phi; // For H1  test functions Q
    //TPZFMatrix<REAL> dphiSigmaH1 = datavec[sigmaBlock].dphi; // Derivative For H1  test functions
    TPZFMatrix<REAL> dphiSigmaH1axes = datavec[sigmaBlock].dphix; // Derivative For H1  test functions

    TPZFNMatrix<660> GradphiuH1;
    TPZAxesTools<REAL>::Axes2XYZ(dphiSigmaH1axes, GradphiuH1, datavec[sigmaBlock].axes);

    int nphiuHdiv = datavec[sigmaBlock].fVecShapeIndex.NElements();

    DivergenceofPhi.Resize(nphiuHdiv, 1);

    REAL JacobianDet = datavec[sigmaBlock].detjac;

    TPZFMatrix<REAL> Qaxes = datavec[sigmaBlock].axes;
    TPZFMatrix<REAL> QaxesT;
    TPZFMatrix<REAL> Jacobian = datavec[sigmaBlock].jacobian;
    TPZFMatrix<REAL> JacobianInverse = datavec[sigmaBlock].jacinv;

    TPZFMatrix<REAL> GradOfX;
    TPZFMatrix<REAL> GradOfXInverse;
    TPZFMatrix<REAL> VectorOnMaster;
    TPZFMatrix<REAL> VectorOnXYZ(3, 1, 0.0);
    Qaxes.Transpose(&QaxesT);
    QaxesT.Multiply(Jacobian, GradOfX);
    JacobianInverse.Multiply(Qaxes, GradOfXInverse);

    int ivectorindex = 0;
    int ishapeindex = 0;

    {
        for (int iq = 0; iq < nphiuHdiv; iq++) {
            ivectorindex = datavec[sigmaBlock].fVecShapeIndex[iq].first;
            ishapeindex = datavec[sigmaBlock].fVecShapeIndex[iq].second;

            VectorOnXYZ(0, 0) = datavec[sigmaBlock].fDeformedDirections(0, ivectorindex);
            VectorOnXYZ(1, 0) = datavec[sigmaBlock].fDeformedDirections(1, ivectorindex);
            VectorOnXYZ(2, 0) = datavec[sigmaBlock].fDeformedDirections(2, ivectorindex);

            GradOfXInverse.Multiply(VectorOnXYZ, VectorOnMaster);
            VectorOnMaster *= JacobianDet;

            /* Contravariant Piola mapping preserves the divergence */
//            DivergenceofPhi(iq, 0) = (1.0 / JacobianDet) * (dphiuH1(0, ishapeindex) * VectorOnMaster(0, 0) +
//                    dphiuH1(1, ishapeindex) * VectorOnMaster(1, 0) +
//                    dphiuH1(2, ishapeindex) * VectorOnMaster(2, 0));
        }
    }
    return;

}

void TPZMixedElasticityMaterial::ElasticityModulusTensor(TPZFMatrix<STATE> &MatrixElast, TElasticityAtPoint &elast) {
    //Matrix modulus Voigt notation:
    MatrixElast.Redim(4, 4);
    REAL mu = elast.fmu;
    REAL lambda = elast.flambda;
    if (!fPlaneStress) {
        // plane strain
        MatrixElast(Exx, Exx) = 1. / (2. * mu) - lambda / (2. * mu * (2. * lambda + 2. * mu));
        MatrixElast(Exx, Eyy) = -lambda / (2. * mu * (2. * lambda + 2. * mu));
        MatrixElast(Eyy, Exx) = MatrixElast(Exx, Eyy);
        MatrixElast(Eyy, Eyy) = MatrixElast(Exx, Exx);
        MatrixElast(Exy, Exy) = 1. / (2. * mu);
        MatrixElast(Eyx, Eyx) = MatrixElast(Exy, Exy);
    } else {
        // plane stress
        MatrixElast(Exx, Exx) = (1 - elast.fnu*elast.fnu) / elast.fE;
        MatrixElast(Exx, Eyy) = -elast.fnu * (1 + elast.fnu) / elast.fE;
        MatrixElast(Eyy, Exx) = MatrixElast(Exx, Eyy);
        MatrixElast(Eyy, Eyy) = MatrixElast(Exx, Exx);
        MatrixElast(Exy, Exy) = 1. / (2. * mu);
        MatrixElast(Eyx, Eyx) = MatrixElast(Exy, Exy);
    }
}

void TPZMixedElasticityMaterial::ComputeDeformationVector(TPZVec<STATE> &PhiStress, TPZVec<STATE> &APhiStress, TElasticityAtPoint &elast) {
    TPZFNMatrix<16, STATE> MatrixElast(4, 4, 0.);
    ElasticityModulusTensor(MatrixElast, elast);

    for (int iq = 0; iq < 4; iq++) {
        APhiStress[iq] = 0.;
        for (int jq = 0; jq < 4; jq++) {
            APhiStress[iq] += MatrixElast(iq, jq) * PhiStress[jq];
        }
    }
}

void TPZMixedElasticityMaterial::ComputeStressVector(TPZVec<STATE> &Deformation, TPZVec<STATE> &Stress, TElasticityAtPoint &elast) {
    if (fPlaneStress) {
        Stress[Exx] = elast.fE / (1 - elast.fnu * elast.fnu)*(Deformation[Exx] + elast.fnu * Deformation[Eyy]);
        Stress[Eyy] = elast.fE / (1. - elast.fnu * elast.fnu)*(elast.fnu * Deformation[Exx] + Deformation[Eyy]);
        Stress[Exy] = elast.fE / (1. + elast.fnu) * Deformation[Exy];
        Stress[Eyx] = elast.fE / (1. + elast.fnu) * Deformation[Eyx];
    } else {
        Stress[Exx] = elast.fE / ((1 + elast.fnu)*(1 - 2. * elast.fnu))*((1 - elast.fnu) * Deformation[Exx] + elast.fnu * Deformation[Eyy]);
        Stress[Eyy] = elast.fE / ((1 + elast.fnu)*(1 - 2. * elast.fnu))*((1 - elast.fnu) * Deformation[Eyy] + elast.fnu * Deformation[Exx]);
        Stress[Exy] = elast.fE / (1. + elast.fnu) * Deformation[Exy];
        Stress[Eyx] = elast.fE / (1. + elast.fnu) * Deformation[Eyx];
    }
}

void TPZMixedElasticityMaterial::Print(std::ostream &out) {
    out << "name of material : " << Name() << "\n";
    out << "properties : \n";
    out << "\tE_const   = " << fE_const << endl;
    out << "\tnu_const   = " << fnu_const << endl;
    out << "\tF   = " << fForce[0] << ' ' << fForce[1] << endl;
}

/// Transform a tensor to a Voigt notation
void TPZMixedElasticityMaterial::ToVoigt(TPZFMatrix<STATE> &S, TPZVec<STATE> &Svoigt) {
    Svoigt[Exx] = S(0, 0);
    Svoigt[Exy] = S(0, 1);
    Svoigt[Eyx] = S(1, 0);
    Svoigt[Eyy] = S(1, 1);
}

/// Transform a Voigt notation to a tensor
void TPZMixedElasticityMaterial::FromVoigt(TPZVec<STATE> &Svoigt, TPZFMatrix<STATE> &S) {
    S(0, 0) = Svoigt[Exx];
    S(0, 1) = Svoigt[Exy];
    S(1, 0) = Svoigt[Eyx];
    S(1, 1) = Svoigt[Eyy];
}

void TPZMixedElasticityMaterial::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
    if (datavec[0].fVecShapeIndex.size() == 0) {
        FillVecShapeIndex(datavec[0]);
    }

    REAL R = datavec[0].x[0];
    // Setting the phi's
    // E
    TPZFMatrix<REAL> &phiS = datavec[0].phi;
    TPZFMatrix<REAL> &dphiS = datavec[0].dphix;
    // U
    TPZFMatrix<REAL> &phiU = datavec[1].phi;
    TPZFMatrix<REAL> &dphiU = datavec[1].dphix;
    // P
    TPZFMatrix<REAL> &phiP = datavec[2].phi;

    
    TElasticityAtPoint elast(fE_const,fnu_const);
    if(fElasticity)
    {
        //TPZManVector<REAL,3> result(2);
		TPZManVector<STATE, 3> result(2);
		TPZFNMatrix<4,STATE> Dres(0,0);
        fElasticity->Execute(datavec[0].x, result, Dres);
        REAL E = result[0];
        REAL nu = result[1];
        TElasticityAtPoint modify(E,nu);
        elast = modify;
    }
    //    
//    TPZFNMatrix<220, REAL> dphiSx(fDimension, dphiS.Cols());
//    TPZAxesTools<REAL>::Axes2XYZ(dphiS, dphiSx, datavec[0].axes);

//    TPZFNMatrix<220, REAL> dphiUx(fDimension, phiU.Cols());
//    TPZAxesTools<REAL>::Axes2XYZ(dphiU, dphiUx, datavec[1].axes);

//    TPZFNMatrix<220, REAL> dphiPx(fDimension, phiP.Cols());
//    TPZAxesTools<REAL>::Axes2XYZ(dphiU, dphiPx, datavec[2].axes);

    int nshapeS, nshapeU, nshapeP;
    nshapeS = datavec[0].fVecShapeIndex.NElements();
    nshapeU = datavec[1].phi.Rows();
    nshapeP = datavec[2].phi.Rows();

    TPZManVector<STATE, 2> force(2, 0.);
    TPZFMatrix<STATE> phiSi(3, 1, 0.), phiSj(fDimension, 1, 0.);
    TPZManVector<STATE, 2> divSi1x(2, 0.), divSi1y(2, 0.);
    TPZManVector<STATE, 4> phiSi1x(4, 0.0), phiSj1x(4, 0.0), phiSi1y(4, 0.0), phiSj1y(4, 0.0), phiPj1x(4, 0.0), phiPj1y(4, 0.0);
    TPZManVector<STATE, 2> phiUi(fDimension, 0.0), phiUj(fDimension, 0.0), phiUj1x(2, 0.0), phiUj1y(2, 0.0);
    TPZFMatrix<STATE> phiPi(fDimension, 1, 0.0), phiPj(fDimension, 1, 0.0);
    TPZFNMatrix<3, REAL> ivecS(3, 1, 0.);

    force[0] = this->fForce[0];
    force[1] = this->fForce[1];
    if (this->HasForcingFunction()) {
        this->ForcingFunction()->Execute(datavec[0].x, force);
#ifdef PZ_LOG
        if (logdata.isDebugEnabled()) {
            std::stringstream sout;
            sout << " x = " << datavec[0].x << " force = " << force << std::endl;
            LOGPZ_DEBUG(logdata, sout.str())
        }
#endif
    }
    for (int i = 0; i < nshapeS; i++) {
        int iphi = datavec[0].fVecShapeIndex[i].second;
        int ivec = datavec[0].fVecShapeIndex[i].first;
        TPZFNMatrix<4> GradSi(fDimension, fDimension);

        REAL divSi = 0.;

        for (int e = 0; e < 3; e++) {

            ivecS(e, 0) = datavec[0].fDeformedDirections(e, ivec);
            phiSi(e, 0) = phiS(iphi, 0) * ivecS(e, 0);

        }
        TPZFNMatrix<3, REAL> axesvec(fDimension, 1, 0.);
        datavec[0].axes.Multiply(ivecS, axesvec);

        //calculando div(Si)
        for (int f = 0; f < fDimension; f++) {
            divSi += axesvec(f, 0) * dphiS(f, iphi);
        }

        TPZFNMatrix<4, STATE> phiTensx(2, 2, 0.), phiTensy(2, 2, 0.);

        phiTensx(0, 0) = phiSi(0, 0);
        phiTensx(0, 1) = phiSi(1, 0);

        phiTensy(1, 0) = phiSi(0, 0);
        phiTensy(1, 1) = phiSi(1, 0);
        ToVoigt(phiTensx, phiSi1x);
        ToVoigt(phiTensy, phiSi1y);

        divSi1x[0] = divSi;
        divSi1y[1] = divSi;

        // matrix K11 - (Matrix A * stress tensor) x test-function stress tensor
        for (int j = 0; j < nshapeS; j++) {
            int jphi = datavec[0].fVecShapeIndex[j].second;
            int jvec = datavec[0].fVecShapeIndex[j].first;

            for (int e = 0; e < fDimension; e++) {
                phiSj(e, 0) = phiS(jphi, 0) * datavec[0].fDeformedDirections(e, jvec);
            }

            TPZFNMatrix<4, STATE> phjTensx(2, 2, 0.), phjTensy(2, 2, 0.);

            phjTensx(0, 0) = phiSj(0, 0);
            phjTensx(0, 1) = phiSj(1, 0);

            phjTensy(1, 0) = phiSj(0, 0);
            phjTensy(1, 1) = phiSj(1, 0);
            ToVoigt(phjTensx, phiSj1x);
            ToVoigt(phjTensy, phiSj1y);

            //Multiply by Lamé parameters
            TPZManVector<STATE, 4> AphiSi1x(4, 0.0), AphiSi1y(4, 0.0);

            ComputeDeformationVector(phiSi1x, AphiSi1x,elast);
            ComputeDeformationVector(phiSi1y, AphiSi1y,elast);

            STATE valxx = InnerVec(AphiSi1x, phiSj1x);
            STATE valxy = InnerVec(AphiSi1x, phiSj1y);
            STATE valyx = InnerVec(AphiSi1y, phiSj1x);
            STATE valyy = InnerVec(AphiSi1y, phiSj1y);

            //Matrix K11
            if (fAxisSymmetric) {
                valxx /= R;
                valxy /= R;
                valyx /= R;
                valyy /= R;
            }
            ek(2 * i, 2 * j) += weight * valxx;
            ek(2 * i, 2 * j + 1) += weight * valxy;
            ek(2 * i + 1, 2 * j) += weight * valyx;
            ek(2 * i + 1, 2 * j + 1) += weight * valyy;
        }

        // matrix K21 and K12 - divergent of test-function stress tensor * displacement vector
        for (int j = 0; j < nshapeU; j++) {
            phiUj1x[0] = phiU(j, 0);
            phiUj1y[1] = phiU(j, 0);

            STATE valx = weight * InnerVec(divSi1x, phiUj1x);
            STATE valy = weight * InnerVec(divSi1y, phiUj1y);

            //position K21
            ek(2 * j + nshapeS * 2, 2 * i) += valx;
            ek(2 * j + 1 + nshapeS * 2, 2 * i + 1) += valy;

            //position K12
            ek(2 * i, 2 * j + nshapeS * 2) += valx;
            ek(2 * i + 1, 2 * j + 1 + nshapeS * 2) += valy;
        }

        // matrix K31 and K13 - test-function stress tensor x rotation tensor p
        for (int j = 0; j < nshapeP; j++) {
            TPZFNMatrix<4, STATE> phiPTensor(2, 2, 0.);
            phiPTensor(0, 1) = phiP(j, 0);
            phiPTensor(1, 0) = -phiP(j, 0);
            ToVoigt(phiPTensor, phiPj1x);

            STATE valxx = InnerVec(phiSi1x, phiPj1x);
            STATE valxy = InnerVec(phiSi1y, phiPj1x);
            //Matrix K31
            ek(j + nshapeS * 2 + nshapeU * 2, 2 * i) += weight * valxx;
            ek(j + nshapeS * 2 + nshapeU * 2, 2 * i + 1) += weight * valxy;

            //Matrix K13
            ek(2 * i, j + nshapeS * 2 + nshapeU * 2) += weight * valxx;
            ek(2 * i + 1, j + nshapeS * 2 + nshapeU * 2) += weight * valxy;
        }
    }

    for (int i = 0; i < nshapeU; i++) {
        phiUj1x[0] = phiU(i, 0);
        phiUj1y[1] = phiU(i, 0);

        // Load vector f:
        STATE factfx = -weight * phiUj1x[0] * force[0];
        STATE factfy = -weight * phiUj1y[1] * force[1];
        if (fAxisSymmetric) {
            factfx *= R;
            factfy *= R;
        }
        //if(factfx != 0) DebugStop();
        ef(nshapeS * 2 + 2 * i, 0) += factfx;
        ef(nshapeS * 2 + 2 * i + 1, 0) += factfy;

        //if(ef(nshapeS*2+2*j) != 0.) DebugStop();
        if (fAxisSymmetric)
        {
            for (int j = 0; j < nshapeU; j++) {
                ek(nshapeS * 2 + 2 * i, nshapeS * 2 + 2 * j) -= weight * phiU(i, 0) * phiU(j, 0) / (elast.fE * R);
            }
        }
    }
    ///    std::ofstream filestiff("ek.txt");
    //    ek.Print("K1 = ",filestiff,EMathematicaInput);
}

void TPZMixedElasticityMaterial::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
    TPZMaterialData::MShapeFunctionType shapetype = data.fShapeType;
    if (shapetype == data.EVecShape) {
        DebugStop();
        return;
    }

    TPZFMatrix<REAL> &dphiU = data.dphix;
    TPZFMatrix<REAL> &phi = data.phi;
    TPZFMatrix<REAL> &axes = data.axes;

    int phc, phrU, dphc, dphr, efr, efc, ekr, ekc;
    phc = phi.Cols();
    phrU = phi.Rows();
    int FirstU = 0;

    dphc = dphiU.Cols();
    dphr = dphiU.Rows();
    efr = ef.Rows();
    efc = ef.Cols();
    ekr = ek.Rows();
    ekc = ek.Cols();
    if (phc != 1 || dphr != 2 || phrU != dphc) {
        PZError << "\nTPZMixedElasticityMaterial.contr, inconsistent input data : \n" <<
                "phi.Cols() = " << phi.Cols() << " dphi.Cols() = " << dphiU.Cols() <<
                " phi.Rows = " << phi.Rows() << " dphi.Rows = " <<
                dphiU.Rows() << "\nek.Rows() = " << ek.Rows() << " ek.Cols() = "
                << ek.Cols() <<
                "\nef.Rows() = " << ef.Rows() << " ef.Cols() = "
                << ef.Cols() << "\n";
        return;
        //        PZError.show();
    }
    //TPZManVector<REAL, 2> force(2, 0.);
	TPZManVector<STATE, 2> force(2, 0.);
    force[0] = fForce[0];
    force[1] = fForce[1];
    if (fForcingFunction) { // phi(in, 0) :  node in associated forcing function
        TPZManVector<STATE, 3> res(3);
        fForcingFunction->Execute(data.x, force);
    }
    
    TElasticityAtPoint elast(fE_const,fnu_const);
    if(fElasticity)
    {
        //TPZManVector<REAL,3> result(2);
		TPZManVector<STATE, 3> result(2);
        TPZFNMatrix<4,STATE> Dres(0,0);
        fElasticity->Execute(data.x, result, Dres);
        REAL E = result[0];
        REAL nu = result[1];
        TElasticityAtPoint modify(E,nu);
        elast = modify;
    }

    REAL MuL = elast.fmu;
    REAL LambdaL = elast.flambda;
    //  ////////////////////////// Jacobian Matrix ///////////////////////////////////
    //  Contribution of domain integrals for Jacobian matrix
    //  Elasticity Block (Equation for elasticity )
    //    Elastic equation
    //    Linear strain operator
    //    Ke Matrix
    for (int iu = 0; iu < phrU; iu++) {
        for (int col = 0; col < efc; col++) {
            ef(2 * iu, col) += weight * (force[0] * phi(iu, 0)); // direcao x
            ef(2 * iu + 1, col) += weight * (force[1] * phi(iu, 0)); // direcao y <<<----
        }

        TPZManVector<REAL, 2> dv(2);
        //    Derivative for Vx
        dv[0] = dphiU(0, iu) * axes(0, 0) + dphiU(1, iu) * axes(1, 0);
        //    Derivative for Vy
        dv[1] = dphiU(0, iu) * axes(0, 1) + dphiU(1, iu) * axes(1, 1);

        for (int ju = 0; ju < phrU; ju++) {
            TPZManVector<REAL, 2> du(2);
            //    Derivative for Ux
            du[0] = dphiU(0, ju) * axes(0, 0) + dphiU(1, ju) * axes(1, 0);
            //    Derivative for Uy
            du[1] = dphiU(0, ju) * axes(0, 1) + dphiU(1, ju) * axes(1, 1);

            if (this->fPlaneStress == 1) {
                /* Plain stress state */
                ek(2 * iu + FirstU, 2 * ju + FirstU) += weight * ((4 * (MuL)*(LambdaL + MuL) / (LambdaL + 2 * MuL)) * dv[0] * du[0] + (MuL) * dv[1] * du[1]);
                ek(2 * iu + FirstU, 2 * ju + 1 + FirstU) += weight * ((2 * (MuL)*(LambdaL) / (LambdaL + 2 * MuL)) * dv[0] * du[1] + (MuL) * dv[1] * du[0]);
                ek(2 * iu + 1 + FirstU, 2 * ju + FirstU) += weight * ((2 * (MuL)*(LambdaL) / (LambdaL + 2 * MuL)) * dv[1] * du[0] + (MuL) * dv[0] * du[1]);
                ek(2 * iu + 1 + FirstU, 2 * ju + 1 + FirstU) += weight * ((4 * (MuL)*(LambdaL + MuL) / (LambdaL + 2 * MuL)) * dv[1] * du[1] + (MuL) * dv[0] * du[0]);
            } else {
                /* Plain Strain State */
                ek(2 * iu + FirstU, 2 * ju + FirstU) += weight * ((LambdaL + 2 * MuL) * dv[0] * du[0] + (MuL) * dv[1] * du[1]);
                ek(2 * iu + FirstU, 2 * ju + 1 + FirstU) += weight * (LambdaL * dv[0] * du[1] + (MuL) * dv[1] * du[0]);
                ek(2 * iu + 1 + FirstU, 2 * ju + FirstU) += weight * (LambdaL * dv[1] * du[0] + (MuL) * dv[0] * du[1]);
                ek(2 * iu + 1 + FirstU, 2 * ju + 1 + FirstU) += weight * ((LambdaL + 2 * MuL) * dv[1] * du[1] + (MuL) * dv[0] * du[0]);
            }
        }
    }
}

void TPZMixedElasticityMaterial::FillDataRequirements(TPZMaterialData &data) {
    data.fNeedsSol = true;
    data.fNeedsNormal = false;
}

void TPZMixedElasticityMaterial::FillBoundaryConditionDataRequirement(int type, TPZMaterialData &data) {
    data.fNeedsSol = false;
    data.fNeedsNormal = false;
    if (type == 4 || type == 5 || type == 6) {
        data.fNeedsNormal = true;
    }
}

void TPZMixedElasticityMaterial::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) {
    if (datavec[0].phi.Rows() != 0 && datavec[0].fShapeType != TPZMaterialData::EScalarShape) {
        DebugStop();
    }
    int ndisp = datavec[1].phi.Rows();

    TPZFNMatrix<2, STATE> v_2 = bc.Val2();
    TPZFNMatrix<4, STATE> v_1 = bc.Val1();

    // Setting forcing function
    if (bc.HasForcingFunction()) {
        TPZManVector<STATE, 3> res(2);
        TPZFNMatrix<4, STATE> tens(2, 2);
        bc.ForcingFunction()->Execute(datavec[0].x, res, tens);
        v_2(0, 0) = res[0];
        v_2(1, 0) = res[1];
    }

    // Setting the phis
    // E
    TPZFMatrix<REAL> &phiS = datavec[0].phi;

    int nshapeS;
    nshapeS = datavec[0].phi.Rows();


    //    TPZMaterialData::MShapeFunctionType shapetype = data.fShapeType;
    //    if(shapetype==data.EVecShape){
    //        ContributeVecShapeBC(data,weight,ek, ef,bc);
    //        return;
    //    }

    REAL R = datavec[0].x[0];
    if (R < 1.e-6) R = 1.e-6;

    switch (bc.Type()) {
        case 0: // Dirichlet condition
        {
            for (int iq = 0; iq < nshapeS; iq++) {
                ef(2 * iq, 0) += v_2(0, 0) * phiS(iq, 0) * weight; // forced v2 displacement
                ef(2 * iq + 1, 0) += v_2(1, 0) * phiS(iq, 0) * weight; // forced v2 displacement
            }
        }
            break;

        case 1: // Neumann condition
        {
            for (int iq = 0; iq < nshapeS; iq++) {
                for (int jq = 0; jq < nshapeS; jq++) {
                    if (fAxisSymmetric) {
                        ek(2 * iq, 2 * jq) += gBigNumber * phiS(iq, 0) * phiS(jq, 0) * weight / R;
                        ek(2 * iq + 1, 2 * jq + 1) += gBigNumber * phiS(iq, 0) * phiS(jq, 0) * weight / R;
                    } else {
                        ek(2 * iq, 2 * jq) += gBigNumber * phiS(iq, 0) * phiS(jq, 0) * weight;
                        ek(2 * iq + 1, 2 * jq + 1) += gBigNumber * phiS(iq, 0) * phiS(jq, 0) * weight;
                    }
                }
                ef(2 * iq, 0) += gBigNumber * v_2(0, 0) * phiS(iq, 0) * weight; // normal stress in x direction
                ef(2 * iq + 1, 0) += gBigNumber * v_2(1, 0) * phiS(iq, 0) * weight; // normal stress in y direction
            }
        }
            break;

        case 2: // Mixed condition
        {
            for (int iq = 0; iq < nshapeS; iq++) {
                for (int jq = 0; jq < nshapeS; jq++) {
                    if (fAxisSymmetric) {
                        ek(2 * iq, 2 * jq) += v_1(0, 0) * phiS(iq, 0) * phiS(jq, 0) * weight / R;
                        ek(2 * iq + 1, 2 * jq + 1) += v_1(1, 1) * phiS(iq, 0) * phiS(jq, 0) * weight / R;
                        ek(2 * iq + 1, 2 * jq) += v_1(1, 0) * phiS(iq, 0) * phiS(jq, 0) * weight / R;
                        ek(2 * iq, 2 * jq + 1) += v_1(0, 1) * phiS(iq, 0) * phiS(jq, 0) * weight / R;
                    } else {
                        ek(2 * iq, 2 * jq) += v_1(0, 0) * phiS(iq, 0) * phiS(jq, 0) * weight;
                        ek(2 * iq + 1, 2 * jq + 1) += v_1(1, 1) * phiS(iq, 0) * phiS(jq, 0) * weight;
                        ek(2 * iq + 1, 2 * jq) += v_1(1, 0) * phiS(iq, 0) * phiS(jq, 0) * weight;
                        ek(2 * iq, 2 * jq + 1) += v_1(0, 1) * phiS(iq, 0) * phiS(jq, 0) * weight;
                    }
                }
                ef(2 * iq, 0) += v_2(0, 0) * phiS(iq, 0) * weight; // normal stress in x direction
                ef(2 * iq + 1, 0) += v_2(1, 0) * phiS(iq, 0) * weight; // normal stress in y direction
            }
        }
            break;
        case 5:
            if (ndisp != 1) {
                DebugStop();
            }
            if (v_1(0, 0) > 1.e-1) {
                ek(2 * nshapeS, 2 * nshapeS) = 1.;
            }
            if (v_1(1, 1) > 1.e-1) {
                ek(2 * nshapeS + 1, 2 * nshapeS + 1) = 1.;
            }
            break;
        default:
            DebugStop();
            // nulo introduzindo o BIGNUMBER pelos valores da condição
    } // 1 Val1 : a leitura é 00 01 10 11
}

void TPZMixedElasticityMaterial::ContributeBC(TPZMaterialData &data, REAL weight,
        TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) {


    TPZMaterialData::MShapeFunctionType shapetype = data.fShapeType;
    if (shapetype == data.EVecShape) {
        DebugStop();
        return;
    }

    TPZFMatrix<REAL> &phi = data.phi;
    int dim = Dimension();

    const REAL BIGNUMBER = TPZMaterial::gBigNumber;

    int phr = phi.Rows();
    short in, jn;

    if (ef.Cols() != bc.NumLoadCases()) {
        DebugStop();
    }

    //        In general when the problem is needed to stablish any convention for ContributeBC implementations

    //     REAL v2[2];
    //     v2[0] = bc.Val2()(0,0);
    //     v2[1] = bc.Val2()(1,0);
    int nstate = NStateVariables();

    TPZFMatrix<STATE> &v1 = bc.Val1();

    switch (bc.Type()) {
        case 0: // Dirichlet condition
        {
            for (in = 0; in < phr; in++) {
                for (int il = 0; il < NumLoadCases(); il++) {
                    REAL v2[2];
                    v2[0] = bc.Val2(il)(0, 0);
                    v2[1] = bc.Val2(il)(1, 0);
                    ef(2 * in, il) += BIGNUMBER * v2[0] * phi(in, 0) * weight; // forced v2 displacement
                    ef(2 * in + 1, il) += BIGNUMBER * v2[1] * phi(in, 0) * weight; // forced v2 displacement
                }
                for (jn = 0; jn < phi.Rows(); jn++) {
                    ek(2 * in, 2 * jn) += BIGNUMBER * phi(in, 0) * phi(jn, 0) * weight;
                    ek(2 * in + 1, 2 * jn + 1) += BIGNUMBER * phi(in, 0) * phi(jn, 0) * weight;
                }
            }
        }
            break;

        case 1: // Neumann condition
        {
            for (in = 0; in < phr; in++) {
                for (int il = 0; il < fNumLoadCases; il++) {
                    TPZFNMatrix<2, STATE> v2 = bc.Val2(il);
                    ef(2 * in, il) += v2(0, 0) * phi(in, 0) * weight; // force in x direction
                    ef(2 * in + 1, il) += v2(1, 0) * phi(in, 0) * weight; // force in y direction
                }
            }
        }
            break;

        case 2: // Mixed Condition
        {
            for (in = 0; in < phi.Rows(); in++) {
                for (int il = 0; il < fNumLoadCases; il++) {
                    TPZFNMatrix<2, STATE> v2 = bc.Val2(il);
                    ef(2 * in, il) += v2(0, 0) * phi(in, 0) * weight; // force in x direction
                    ef(2 * in + 1, il) += v2(1, 0) * phi(in, 0) * weight; // forced in y direction
                }

                for (jn = 0; jn < phi.Rows(); jn++) {
                    ek(2 * in, 2 * jn) += bc.Val1()(0, 0) * phi(in, 0) *
                            phi(jn, 0) * weight; // peso de contorno => integral de contorno
                    ek(2 * in + 1, 2 * jn) += bc.Val1()(1, 0) * phi(in, 0) *
                            phi(jn, 0) * weight;
                    ek(2 * in + 1, 2 * jn + 1) += bc.Val1()(1, 1) * phi(in, 0) *
                            phi(jn, 0) * weight;
                    ek(2 * in, 2 * jn + 1) += bc.Val1()(0, 1) * phi(in, 0) *
                            phi(jn, 0) * weight;
                }
            } // este caso pode reproduzir o caso 0 quando o deslocamento


            case 3: // Directional Null Dirichlet - displacement is set to null in the non-null vector component direction
            for (in = 0; in < phr; in++) {
                //                ef(nstate*in+0,0) += BIGNUMBER * (0. - data.sol[0][0]) * v2[0] * phi(in,0) * weight;
                //                ef(nstate*in+1,0) += BIGNUMBER * (0. - data.sol[0][1]) * v2[1] * phi(in,0) * weight;
                for (jn = 0; jn < phr; jn++) {
                    ek(nstate * in + 0, nstate * jn + 0) += BIGNUMBER * phi(in, 0) * phi(jn, 0) * weight * bc.Val2()(0, 0);
                    ek(nstate * in + 1, nstate * jn + 1) += BIGNUMBER * phi(in, 0) * phi(jn, 0) * weight * bc.Val2()(1, 0);
                }//jn
            }//in
            break;


            case 4: // stressField Neumann condition
            {
                REAL v2[2];
                for (in = 0; in < dim; in++) {
                    v2[in] = (v1(in, 0) * data.normal[0] +
                            v1(in, 1) * data.normal[1]);
                }
                // The normal vector points towards the neighbour. The negative sign is there to
                // reflect the outward normal vector.
                for (in = 0; in < phi.Rows(); in++) {
                    ef(nstate * in + 0, 0) += v2[0] * phi(in, 0) * weight;
                    ef(nstate * in + 1, 0) += v2[1] * phi(in, 0) * weight;
                    //    cout << "normal:" << data.normal[0] << ' ' << data.normal[1] << ' ' << data.normal[2] << endl;
                    //    cout << "val2:  " << v2[0]  << endl;
                }
            }
            break;

            case 5://PRESSAO DEVE SER POSTA NA POSICAO 0 DO VETOR v2
            {
                TPZFNMatrix<2, STATE> res(2, 1, 0.);
                for (in = 0; in < phi.Rows(); in++) {
                    for (int il = 0; il < NumLoadCases(); il++) {
                        ef(nstate * in + 0, 0) += (bc.Val2(il)(0, 0) * data.normal[0]) * phi(in, 0) * weight;
                        ef(nstate * in + 1, 0) += (bc.Val2(il)(0, 0) * data.normal[1]) * phi(in, 0) * weight;
                    }
                    for (jn = 0; jn < phi.Rows(); jn++) {
                        for (int idf = 0; idf < 2; idf++) for (int jdf = 0; jdf < 2; jdf++) {
                                ek(nstate * in + idf, nstate * jn + jdf) += bc.Val1()(idf, jdf) * data.normal[idf] * data.normal[jdf] * phi(in, 0) * phi(jn, 0) * weight;
                                //BUG FALTA COLOCAR VAL2
                                //                        DebugStop();
                            }
                    }

                }
            }
            break;

            case 6://PRESSAO DEVE SER POSTA NA POSICAO 0 DO VETOR v2
            {
                TPZFNMatrix<2, STATE> res(2, 1, 0.);
                for (in = 0; in < phi.Rows(); in++) {
                    for (int il = 0; il < NumLoadCases(); il++) {
                        ef(nstate * in + 0, 0) += (bc.Val2(il)(0, 0) * data.normal[0]) * phi(in, 0) * weight;
                        ef(nstate * in + 1, 0) += (bc.Val2(il)(0, 0) * data.normal[1]) * phi(in, 0) * weight;
                    }
                    for (jn = 0; jn < phi.Rows(); jn++) {
                        for (int idf = 0; idf < 2; idf++) for (int jdf = 0; jdf < 2; jdf++) {
                                ek(nstate * in + idf, nstate * jn + jdf) += bc.Val1()(idf, jdf) * phi(in, 0) * phi(jn, 0) * weight;
                                //BUG FALTA COLOCAR VAL2
                                //                        DebugStop();
                            }
                    }

                }

            }
            break;

        } // nulo introduzindo o BIGNUMBER pelos valores da condição
    } // 1 Val1 : a leitura 00 01 10 11
}

/** Returns the variable index associated with the name. */
int TPZMixedElasticityMaterial::VariableIndex(const std::string &name) {
    if (!strcmp("displacement", name.c_str())) return 9;
    if (!strcmp("Displacement", name.c_str())) return 9; //function U ***
    if (!strcmp("DisplacementMem", name.c_str())) return 9;
    if (!strcmp("Pressure", name.c_str())) return 1;
    if (!strcmp("MaxStress", name.c_str())) return 2;
    if (!strcmp("PrincipalStress1", name.c_str())) return 3;
    if (!strcmp("PrincipalStress2", name.c_str())) return 4;
    if (!strcmp("SigmaX", name.c_str())) return 5;
    if (!strcmp("SigmaY", name.c_str())) return 6;
    if (!strcmp("TauXY", name.c_str())) return 8; //Cedric
    if (!strcmp("Strain", name.c_str())) return 11; //Philippe
    if (!strcmp("SigmaZ", name.c_str())) return 12; //Philippe
    if (!strcmp("sig_x", name.c_str())) return 5;
    if (!strcmp("sig_y", name.c_str())) return 6;
    if (!strcmp("tau_xy", name.c_str())) return 8; //Cedric
    if (!strcmp("Displacement6", name.c_str())) return 7;
    if (!strcmp("Stress", name.c_str())) return 10; //function S ***
    if (!strcmp("Flux", name.c_str())) return 10;
    if (!strcmp("J2", name.c_str())) return 20;
    if (!strcmp("I1", name.c_str())) return 21;
    if (!strcmp("J2Stress", name.c_str())) return 20;
    if (!strcmp("I1Stress", name.c_str())) return 21;
    if (!strcmp("Alpha", name.c_str())) return 22;
    if (!strcmp("PlasticSqJ2", name.c_str())) return 22;
    if (!strcmp("PlasticSqJ2El", name.c_str())) return 22;
    if (!strcmp("YieldSurface", name.c_str())) return 27;
    if (!strcmp("NormalStress", name.c_str())) return 23;
    if (!strcmp("ShearStress", name.c_str())) return 24;
    if (!strcmp("NormalStrain", name.c_str())) return 25;
    if (!strcmp("ShearStrain", name.c_str())) return 26;
    if (!strcmp("Rotation", name.c_str())) return 27; //function P ***
    if(!strcmp("Young_Modulus",name.c_str()))        return 28;
    if(!strcmp("Poisson",name.c_str()))        return 29;

    return TPZMaterial::VariableIndex(name);
}

/** Returns the number of variables associated with the variable indexed by var. */
int TPZMixedElasticityMaterial::NSolutionVariables(int var) {
    switch (var) {
        case 0:
            return 2;
        case 1:
        case 2:
            return 1;
        case 3:
        case 4:
            return 2;
        case 5:
        case 6:
        case 8:
            return 1;
        case 7:
            return 6;
        case 9:
            return 3;
        case 10: //Stress Tensor
            return 4;
        case 11: //Strain Tensor
            return 3;
            // SigZ
        case 12:
            return 1;
        case 20:
            return 1;
        case 21:
            return 1;
        case 22:
            return 1;
        case 23:
        case 24:
        case 25:
        case 26:
        case 27:
            return 3;
        case 28:
        case 29:
            return 1;
        default:
            return TPZMaterial::NSolutionVariables(var);
    }
}

void TPZMixedElasticityMaterial::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout) {
    int numbersol = data.dsol.size();
    int ipos = 0;
    if (fPostProcIndex < numbersol) {
        ipos = fPostProcIndex;
    }

    TPZVec<STATE> &Sol = data.sol[ipos];
    TPZFMatrix<STATE> &DSol = data.dsol[ipos];
    TPZFMatrix<REAL> &axes = data.axes;
    TPZFNMatrix<4, STATE> DSolxy(2, 2);

    TElasticityAtPoint elast(fE_const,fnu_const);

    REAL E(fE_const), nu(fnu_const);

    if(fElasticity)
    {
        TPZManVector<REAL,3> x = data.x;
		TPZManVector<STATE, 3> result(2);
		TPZFNMatrix<4,STATE> Dres(0,0);
        fElasticity->Execute(x, result, Dres);
        E = result[0];
        nu = result[1];
        TElasticityAtPoint modify(E,nu);
        elast = modify;
    }

    if(var == 28)
    {
        Solout[0] = E;
        return ;
    }
    if(var == 29)
    {
        Solout[0] = nu;
        return;
    }

    REAL epsx;
    REAL epsy;
    REAL epsxy;
    REAL epsz = 0.;
    REAL SigX;
    REAL SigY;
    REAL SigZ;
    REAL TauXY, aux, Sig1, Sig2, angle;

    // dudx - dudy
    DSolxy(0, 0) = DSol(0, 0) * axes(0, 0) + DSol(1, 0) * axes(1, 0);
    DSolxy(1, 0) = DSol(0, 0) * axes(0, 1) + DSol(1, 0) * axes(1, 1);
    // dvdx - dvdy
    DSolxy(0, 1) = DSol(0, 1) * axes(0, 0) + DSol(1, 1) * axes(1, 0);
    DSolxy(1, 1) = DSol(0, 1) * axes(0, 1) + DSol(1, 1) * axes(1, 1);

    epsx = DSolxy(0, 0); // du/dx
    epsy = DSolxy(1, 1); // dv/dy
    epsxy = 0.5 * (DSolxy(1, 0) + DSolxy(0, 1));

    REAL lambda = elast.flambda;
    REAL mu = elast.fmu;
    if (this->fPlaneStress == 1) {
        epsz = -lambda * (epsx + epsy) / (lambda + 2. * mu);
    } else {
        epsz = 0.;
    }
    TauXY = 2 * mu*epsxy;
#ifdef PZDEBUG
    REAL TauXY2 = elast.fE * epsxy / (1. + elast.fnu);
#ifdef REALfloat
    if (fabs(TauXY - TauXY2) > 1.e-10) {
        DebugStop();
    }
#else
    if (fabs(TauXY - TauXY2) > 1.e-6) {
        DebugStop();
    }
#endif
#endif
    if (this->fPlaneStress == 1) {
        SigX = elast.fE / (1 - elast.fnu * elast.fnu)*(epsx + elast.fnu * epsy);
        SigY = elast.fE / (1 - elast.fnu * elast.fnu)*(elast.fnu * epsx + epsy);
        SigZ = 0.;
    } else {
        SigX = elast.fE / ((1. - 2. * elast.fnu)*(1. + elast.fnu))*((1. - elast.fnu) * epsx + elast.fnu * epsy);
        SigY = elast.fE / ((1. - 2. * elast.fnu)*(1. + elast.fnu))*(elast.fnu * epsx + (1. - elast.fnu) * epsy);
        SigZ = lambda * (epsx + epsy);
    }

    switch (var) {
        case 0:
            //numvar = 2;
            Solout[0] = Sol[0];
            Solout[1] = Sol[1];
            break;
        case 7:
            //numvar = 6;
            Solout[0] = Sol[0];
            Solout[1] = Sol[1];
            Solout[2] = 0.;
            Solout[3] = 0.;
            Solout[4] = 0.;
            Solout[5] = 0.;
            break;
        case 1:
        case 2:
        case 3:
        case 4:
        case 5:
        case 6:
        case 8:
        case 10:
            // Pressure variable
            if (var == 1) {
                //numvar = 1;
                Solout[0] = SigX + SigY + SigZ;
                return;
            }
            // TauXY variable
            if (var == 8) {
                Solout[0] = TauXY;
                return;
            }
            if (var == 5) {
                Solout[0] = SigX;
                return;
            }
            if (var == 6) {
                Solout[0] = SigY;
                return;
            }
            aux = sqrt(0.25 * (SigX - SigY)*(SigX - SigY)
                    +(TauXY)*(TauXY));
            // Philippe 13/5/99
            //         if(abs(Tau) < 1.e-10 && abs(SigY-SigX) < 1.e-10) angle = 0.;
            if (fabs(TauXY) < 1.e-10 && fabs(SigY - SigX) < 1.e-10) angle = 0.;
            else angle = atan2(2 * TauXY, SigY - SigX) / 2.;
            Sig1 = 0.5 * (SigX + SigY) + aux;
            Sig2 = 0.5 * (SigX + SigY) - aux;
            if (var == 3) {
                //numvar = 2;
                Solout[0] = Sig1 * cos(angle);
                Solout[1] = Sig1 * sin(angle);
                return;
            }
            if (var == 4) {
                //numvar = 2;
                Solout[0] = -Sig2 * sin(angle);
                Solout[1] = Sig2 * cos(angle);
                return;
            }
            if (var == 2) {
                REAL sigmax;
                sigmax = (fabs(Sig1) < fabs(Sig2)) ? fabs(Sig2) : fabs(Sig1);
                Solout[0] = sigmax;
                return;
            }
            if (var == 10) {
                Solout[Exx] = SigX;
                Solout[Eyy] = SigY;
                Solout[Exy] = TauXY;
                Solout[Eyx] = TauXY;
                return;
            }
            cout << "Very critical error TPZMixedElasticityMaterial::Solution\n";
            exit(-1);
            //         Solout[0] /= 0.;
            break;
        case 9:
            Solout[0] = Sol[0];
            Solout[1] = Sol[1];
            Solout[2] = 0.;
            break;
        case 11:
            Solout[0] = epsx;
            Solout[1] = epsy;
            Solout[2] = epsxy;
            break;
        case 12:
            Solout[0] = SigZ;
            break;

        case 20:
        {
            REAL J2 = (pow(SigX + SigY, 2) - (3 * (-pow(SigX, 2) - pow(SigY, 2) + pow(SigX + SigY, 2) - 2 * pow(TauXY, 2))) / 2.) / 2.;
            Solout[0] = J2;
            break;
        }
        case 21:
        {
            REAL I1 = SigX + SigY;
            Solout[0] = I1;
            break;
        }
        case 22:
            Solout[0] = 0.;
            break;
        case 23:
            // normal stress
            Solout[0] = SigX;
            Solout[1] = SigY;
            Solout[2] = SigZ;
            break;
        case 24:
            // shear stress
            Solout[0] = TauXY;
            Solout[1] = 0.;
            Solout[2] = 0.;
            break;
        case 25:
            Solout[0] = epsx;
            Solout[1] = epsy;
            Solout[2] = epsz;
            break;
        case 26:
            Solout[0] = epsxy;
            Solout[1] = 0.;
            Solout[2] = 0.;
            break;
        case 27:
            Solout[0] = 0.;
            Solout[1] = 0.;
            Solout[2] = 0.;
            break;
        default:
            TPZMaterial::Solution(Sol, DSol, axes, var, Solout);
            break;
    }
}

/** @brief Returns the solution associated with the var index based on the finite element approximation */
void TPZMixedElasticityMaterial::Solution(TPZVec<TPZMaterialData> &data, int var, TPZVec<STATE> &Solout) {
#ifdef PZDEBUG
    if (data.size() != 3) {
        DebugStop();
    }
#endif
    TPZManVector<REAL, 3> x = data[0].x;
    REAL R = x[0];
    if (R < 1.e-6) {
        R = 1.e-6;
    }
    TPZFNMatrix<9, STATE> sigma(3, 3, 0.), sigmah(3, 3, 0.), eps(3, 3, 0.);
    int dim = Dimension();
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < 3; j++) {
            if (fAxisSymmetric) {
                sigma(i, j) = data[0].sol[0][j + i * 3] / R;
            } else {
                sigma(i, j) = data[0].sol[0][j + i * 3];
            }
        }
    }
    TPZManVector<STATE, 2> disp(2);
    for (int i = 0; i < dim; i++) {
        disp[i] = data[1].sol[0][i];
    }
    TPZFNMatrix<4, STATE> antisym(2, 2, 0.);
    antisym(0, 1) = data[2].sol[0][0];
    antisym(1, 0) = -antisym(0, 1);

    TElasticityAtPoint elast(fE_const,fnu_const);
    if(fElasticity)
    {
        //TPZManVector<REAL,3> result(2);
		TPZManVector<STATE, 3> result(2);
        TPZFNMatrix<4,STATE> Dres(0,0);
        fElasticity->Execute(x, result, Dres);
        REAL E = result[0];
        REAL nu = result[1];
        TElasticityAtPoint modify(E,nu);
        elast = modify;
    }
    
    REAL mu = elast.fmu;
    REAL E = elast.fE;
    REAL nu = elast.fnu;
    
    if(var == 28)
    {
        Solout[0] = E;
        return ;
    }
    if(var == 29)
    {
        Solout[0] = elast.fnu;
        return;
    }
    

    REAL Pressure;

    //TPZManVector<REAL, 4> SIGMA(4, 0.), EPSZ(4, 0.);
	TPZManVector<STATE, 4> SIGMA(4, 0.), EPSZ(4, 0.);

    ToVoigt(sigma, SIGMA);

    ComputeDeformationVector(SIGMA, EPSZ, elast);

    FromVoigt(EPSZ, eps);

    if (this->fPlaneStress == 1) {
        eps(2, 2) = -1. / E * nu * (sigma(0, 0) + sigma(1, 1));
    } else {
        sigma(2, 2) = nu * (sigma(0, 0) + sigma(1, 1));
        eps(2, 2) = 0;
    }

    Pressure = -1/3. * (sigma(0, 0) + sigma(1, 1) + sigma(2,2));
    sigmah(0, 1) = sigma(0, 1);
    sigmah(1, 0) = sigma(1, 0);
    sigmah(0, 0) = sigma(0, 0) + Pressure;
    sigmah(1, 1) = sigma(1, 1) + Pressure;
    sigmah(2, 2) = sigma(2, 2) + Pressure;
    // Displacement
    if (var == 9) {
        Solout[0] = disp[0];
        Solout[1] = disp[1]; //displacement is discontinuous. Can we write as it?
        return;
    }
    // Pressure
    if (var == 1) {
        Solout[0] = -1 / 3 * (sigma(0, 0) + sigma(1, 1) + sigma(2, 2));
        return;
    }
    // Sigmaz
    if (var == 12) {
        Solout[0] = sigma(2, 2);

        return;
    }
    // SigmaY
    if (var == 6) {
        Solout[0] = sigma(1, 1);
        return;
    }

    // SigmaX                
    if (var == 5) {
        Solout[0] = sigma(0, 0);
        return;
    }
    //  TauXY
    if (var == 8) {
        Solout[0] = 0.5 * (sigma(0, 1) + sigma(1, 0));
        return;
    }
    //Strain
    if (var == 11) {
        Solout[0] = eps(0, 0);
        Solout[1] = eps(1, 0);
        Solout[2] = eps(0, 1);
        Solout[3] = eps(2, 2);
        return;
    }

    //Stress
    if (var == 10) {
        Solout[Exx] = sigma(0, 0);
        Solout[Eyx] = sigma(1, 0);
        Solout[Exy] = sigma(0, 1);
        Solout[Eyy] = sigma(1, 1);
        return;
    }

    //I1
    if (var == 21) {
        Solout[0] = sigma(0, 0) + sigma(1, 1) + sigma(2, 2);
        return;
    }
    //J2
    if (var == 20) {
        Solout[0] = sigmah(1, 1) * sigmah(0, 0) + sigmah(1, 1) * sigma(2, 2)
         + sigmah(0,0) * sigmah(2,2) - sigmah(1, 0) * sigmah(1, 0) - sigmah(0, 1) * sigmah(0, 1);
        return;
    }
    //NormalStrain?

    //PrincipalStrain1
    if (var == 3) {
        Solout[0] = Pressure + sqrt(0.25 * (sigma(0, 0) - sigma(1, 1))*(sigma(0, 0) - sigma(1, 1)) + sigma(1, 2));
        return;
    }
    //Rotation
    if (var == 27) {
        Solout[0] = antisym(0, 1);
        return;
    }
}


////////////////////////////////////////////////////////////////////

STATE TPZMixedElasticityMaterial::Inner(TPZFMatrix<STATE> &S, TPZFMatrix<STATE> &T) {
    //inner product of two tensors
#ifdef DEBUG
    if (S.Rows() != S.Cols() || T.Cols() != T.Rows() || S.Rows() != T.Rows()) {
        DebugStop();
    }
#endif

    STATE Val = 0;
    for (int i = 0; i < S.Cols(); i++) {
        for (int j = 0; j < S.Cols(); j++) {
            Val += S(i, j) * T(i, j);
        }
    }
    return Val;
}

////////////////////////////////////////////////////////////////////

template <typename TVar>
TVar TPZMixedElasticityMaterial::InnerVec(const TPZVec<TVar> &S, const TPZVec<TVar> &T) {
    //inner product of two vectors
#ifdef DEBUG
    if (S.size() != T.size()) {
        DebugStop();
    }
#endif
    TVar Val = 0;
    for (int i = 0; i < S.size(); i++) {
        Val += S[i] * T[i];
    }
    return Val;
}

////////////////////////////////////////////////////////////////////
STATE TPZMixedElasticityMaterial::Tr(TPZFMatrix<REAL> &GradU) {

#ifdef DEBUG
    if (GradU.Rows() != GradU.Cols()) {
        DebugStop();
    }
#endif

    STATE Val = 0.;

    for (int i = 0; i < GradU.Rows(); i++) {
        Val += GradU(i, i);
    }

    return Val;
}


/// transform a H1 data structure to a vector data structure

void TPZMixedElasticityMaterial::FillVecShapeIndex(TPZMaterialData &data) {
    data.fDeformedDirections.Resize(fDimension, fDimension);
    data.fDeformedDirections.Identity();
    data.fVecShapeIndex.Resize(fDimension * data.phi.Rows());
    for (int d = 0; d < fDimension; d++) {
        for (int i = 0; i < data.phi.Rows(); i++) {
            data.fVecShapeIndex[i * fDimension + d].first = d;
            data.fVecShapeIndex[i * fDimension + d].second = i;
        }
    }
}

int TPZMixedElasticityMaterial::NEvalErrors() {
    return 8;
}

void TPZMixedElasticityMaterial::Errors(TPZVec<TPZMaterialData> &data, TPZVec<REAL> &errors) {
    //values[0] = 0.;
    //TPZManVector<REAL, 4> SigmaV(4, 0.), sigma_exactV(4, 0.), eps_exactV(4, 0.), EPSZV(4, 0.);
    TPZVec<STATE> u_exact;
    TPZFMatrix<STATE> du_exact;
    if(fExactSol)
    {
        fExactSol->Execute(data[0].x, u_exact, du_exact);
    }
    else DebugStop();
    
    
	TPZManVector<STATE, 4> SigmaV(4, 0.), sigma_exactV(4, 0.), eps_exactV(4, 0.), EPSZV(4, 0.);
	TPZFNMatrix<9, STATE> sigma(2, 2, 0.), eps(2, 2, 0.), grad(2, 2, 0.);
    TPZFNMatrix<4, STATE> eps_exact(2, 2, 0.);
    TPZManVector<REAL, 3> x = data[0].x;
    int dim = Dimension();
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            sigma(i, j) = data[0].sol[0][j + i * 3];
            if (fAxisSymmetric) {
                sigma(i, j) /= x[0];
            }
        }
    }
    ToVoigt(sigma, SigmaV);
    TPZManVector<STATE, 2> divSigma(2,0.);
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            divSigma[i] += data[0].dsol[0](i * 3 + j, j);
        }
    }

    TPZManVector<STATE, 2> disp(2);
    for (int i = 0; i < dim; i++) {
        disp[i] = data[1].sol[0][i];
    }

#ifdef PZ_LOG
    if (logdata.isDebugEnabled())
 {
        std::stringstream sout;
        sout << "DISP*************************** = " << disp << std::endl;
        sigma.Print("sigma************************************ = ", sout, EMathematicaInput);
        LOGPZ_DEBUG(logdata, sout.str())
    }
#endif


    eps_exact(0, 0) = du_exact(0, 0);
    eps_exact(1, 0) = eps_exact(0, 1) = 0.5 * (du_exact(0, 1) + du_exact(1, 0));
    eps_exact(1, 1) = du_exact(1, 1);
    ToVoigt(eps_exact, eps_exactV);
    
    TElasticityAtPoint elast(fE_const,fnu_const);
    if(fElasticity)
    {
        //TPZManVector<REAL,3> result(2);
		TPZManVector<STATE, 3> result(2);
        TPZFNMatrix<4,STATE> Dres(0,0);
        fElasticity->Execute(x, result, Dres);
        REAL E = result[0];
        REAL nu = result[1];
        TElasticityAtPoint modify(E,nu);
        elast = modify;
    }
    
    ComputeStressVector(eps_exactV, sigma_exactV, elast);

    STATE rotation = data[2].sol[0][0];
    STATE rotationExact = 0.5*(du_exact(1, 0)-du_exact(0, 1));
#ifdef PZDEBUG
    {
        TPZManVector<STATE, 4> eps_again(4);
        ComputeDeformationVector(sigma_exactV, eps_again, elast);
        for (int i = 0; i < 4; i++) {
            if (abs(eps_again[i] - eps_exactV[i]) > 1.e-8) {
                DebugStop();
            }
        }
    }
#endif
    //TPZFMatrix<STATE> du(dudx.Rows(),dudx.Cols());
    //du(0,0) = dudx(0,0)*axes(0,0)+dudx(1,0)*axes(1,0);
    //du(1,0) = dudx(0,0)*axes(0,1)+dudx(1,0)*axes(1,1);
    //du(0,1) = dudx(0,1)*axes(0,0)+dudx(1,1)*axes(1,0);
    //du(1,1) = dudx(0,1)*axes(0,1)+dudx(1,1)*axes(1,1);

    //tensões aproximadas : uma forma
    //gamma = du(1,0)+du(0,1);
    //sigma[0] = fPreStressXX+fEover1MinNu2*(du(0,0)+fnu*du(1,1));
    //sigma[1] = fPreStressYY+fEover1MinNu2*(fnu*du(0,0)+du(1,1));
    //sigma[2] = fPreStressXY+fE*0.5/(1.+fnu)*gamma;

    //exata
    //gamma = du_exact(1,0)+du_exact(0,1);
    //sigma_exactV[0] = fEover1MinNu2*(du_exact(0,0)+fnu*du_exact(1,1));
    //sigma_exactV[1] = fEover1MinNu2*(fnu*du_exact(0,0)+du_exact(1,1));
    //sigma_exactV[2] = fE*0.5/(1.+fnu)*gamma;
    //sigx  = (SigmaV[0] - sigma_exactV[0]);
    //sigy  = (SigmaV[1] - sigma_exactV[1]);
    //sigxy = (SigmaV[2] - sigma_exactV[2]);
    
    //L2_Error: for the displacement (u)
    errors[3] = (disp[0] - u_exact[0])*(disp[0] - u_exact[0])+(disp[1] - u_exact[1])*(disp[1] - u_exact[1]);
    //Energe norm
    //TPZManVector<REAL,4> SIGMA(4,0.) , EPSZ(4,0.);

    
    ToVoigt(sigma, SigmaV);

    ComputeDeformationVector(SigmaV, EPSZV, elast);

    //FromVoigt(EPSZV, eps);

    //SIGMA[0] = sigx;
    //SIGMA[1] = sigxy;
    //SIGMA[2] = sigxy;
    //SIGMA[3] = sigy;
    errors[0] = 0.;
    errors[1] = 0.;
    for (int i = 0; i < 4; i++) {
        //L2_Error: for the stress tensor (sigma)
        errors[0] += (SigmaV[i] - sigma_exactV[i])*(SigmaV[i] - sigma_exactV[i]);
        
        //Energy_Error: for the stress tensor (sigma)
        errors[1] += (SigmaV[i] - sigma_exactV[i])*(EPSZV[i] - eps_exactV[i]);
    }
    if (errors[1] < 0.) {
        std::cout << "I should stop \n";
    }
    //L2_Error: for the sigma_xx
    errors[7] = (SigmaV[0] - sigma_exactV[0])*(SigmaV[0] - sigma_exactV[0]);

    TPZManVector<STATE,3> divSigmaExact(fDimension,0.);
    if(HasForcingFunction())
    {
        fForcingFunction->Execute(data[0].x, divSigmaExact);
        for(int i=0; i<fDimension; i++) divSigmaExact[i] *= -1.;
    }
    //L2_Error: divergent of the stress tensor (div(sigma))
    errors[2] = pow(divSigma[0]-divSigmaExact[0],2) + pow(divSigma[1]-divSigmaExact[1],2);
    //L2_Error: for the rotation (q)
    errors[4] = pow(rotation-rotationExact,2);
    //L2_Error: for the symmetry measure (asym)
    errors[5] = pow(SigmaV[Exy]-SigmaV[Eyx],2);
    
    //Energy_Norm: for the exact displacement solution
    errors[6] = 0.;
    for(int i=0; i<4; i++)
    {
        errors[6] += eps_exactV[i]*sigma_exactV[i];
    }
    //    std::cout << "x " << data[0].x << std::endl;
    //    std::cout << "disp " << u_exact << std::endl;
    //    du_exact.Print("du ",std::cout);
    //    std::cout << "sigma_exact " << sigma_exactV << std::endl;
    //    std::cout << errors << std::endl;
    //or we can compute as this methods
    //TPZFMatrix<STATE> MatrixElast(4,4,0.);
    //ElasticityModulusTensor(MatrixElast);
    //errors[1] = 0;
    //for(int i==0; i<4； i++)
    //{
    //for(int j=0; j<4; j++)
    //{
    //errors[1]+=SIGMA[j]*MatrixElast(j,i);
    //}
    //}

    // SemiH1 norm
    //TPZFNMatrix<4,REAL> antisym(2,2,0.);
    //antisym(0,1) = data[2].sol[0][0];
    //antisym(1,0) = -antisym(0,1);
    //TPZManVector<REAL,4> P(4,0.),DU(4,0.);
    //ToVoigt(antisym,P);
    //ToVoigt(du_exact,DU);
    //errors[2]=0;
    //for(int i=0;i<4;i++)
    //{
    //   errors[2]+=(P[i]+EPSZ[i]-DU[i])*(P[i]+EPSZ[i]-DU[i]);
    //}
    //values[0] = calculo do erro estimado em norma Energia
    //values[0] = fE*(sigx*sigx + sigy*sigy + 2*fnu*sigx*sigy)/(1-fnu*fnu);
    //values[0] = (values[0] + .5*fE*sigxy*sigxy/(1+fnu));

    //values[1] : erro em norma L2 em tens�s
    //values[1] = sigx*sigx + sigy*sigy + sigxy*sigxy;

    //values[1] : erro em norma L2 em deslocamentos
    //values[1] = pow((REAL)fabs(u[0] - u_exact[0]),(REAL)2.0)+pow((REAL)fabs(u[1] - u_exact[1]),(REAL)2.0);

    //values[2] : erro estimado na norma H1
    //REAL SemiH1 =0.;
    //for(int i = 0; i < 2; i++) for(int j = 0; j < 2; j++) SemiH1 += (du(i,j) - du_exact(i,j)) * (du(i,j) - du_exact(i,j));
    //values[2] = values[1] + SemiH1;
}

void TPZMixedElasticityMaterial::Errors(TPZVec<REAL> &x, TPZVec<STATE> &u,
        TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes, 
        TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &values) {
    values[0] = 0.;
    TPZVec<REAL> sigma(3, 0.), sigma_exact(3, 0.);
    REAL sigx, sigy, sigxy, gamma;
    TPZFMatrix<STATE> du(dudx.Rows(), dudx.Cols());
    du(0, 0) = dudx(0, 0) * axes(0, 0) + dudx(1, 0) * axes(1, 0);
    du(1, 0) = dudx(0, 0) * axes(0, 1) + dudx(1, 0) * axes(1, 1);
    du(0, 1) = dudx(0, 1) * axes(0, 0) + dudx(1, 1) * axes(1, 0);
    du(1, 1) = dudx(0, 1) * axes(0, 1) + dudx(1, 1) * axes(1, 1);
    TPZManVector<STATE, 4> deform(4), deformexact(4), stress(4), stressexact(4), deformerror(4), stresserror(4);
    
    TElasticityAtPoint elast(fE_const,fnu_const);
    if(fElasticity)
    {
        //TPZManVector<REAL,3> result(2);
		TPZManVector<STATE, 3> result(2);
        TPZFNMatrix<4,STATE> Dres(0,0);
        fElasticity->Execute(x, result, Dres);
        REAL E = result[0];
        REAL nu = result[1];
        TElasticityAtPoint modify(E,nu);
        elast = modify;
    }
    
    ToVoigt(du, deform);
    ComputeStressVector(deform, stress, elast);

    ToVoigt(du_exact, deformexact);
    ComputeStressVector(deformexact, stressexact, elast);
    //exata
    for (int i = 0; i < 4; i++) {
        deformerror[i] = deform[i] - deformexact[i];
        stresserror[i] = stress[i] - stressexact[i];
    }
    //values[0] = calculo do erro estimado em norma Energia
    values[0] = 0.;
    for (int i = 0; i < 4; i++) {
        values[0] += stresserror[i] * deformerror[i];
    }

    //values[1] : erro em norma L2 em tens�s
    //values[1] = sigx*sigx + sigy*sigy + sigxy*sigxy;

    //values[1] : erro em norma L2 em deslocamentos
    values[1] = 0;
    for (int i = 0; i < 2; i++) {
        values[i] += (u[i] - u_exact[i])*(u[i] - u_exact[i]);
    }

    //values[2] : erro estimado na norma H1
    values[2] = 0.;
    for (int i = 0; i < 4; i++) {
        values[2] += deformerror[i] * deformerror[i];
    }
}

TPZMixedElasticityMaterial::TPZMixedElasticityMaterial(const TPZMixedElasticityMaterial &copy) :
TPZMaterial(copy),
fE_const(copy.fE_const),
fnu_const(copy.fnu_const),
flambda_const(copy.flambda_const),
fmu_const(copy.fmu_const),
fForce(copy.fForce),
fPlaneStress(copy.fPlaneStress),
fDimension(copy.fDimension),
fMatrixA(copy.fMatrixA),
fElasticity(copy.fElasticity),
fAxisSymmetric(copy.fAxisSymmetric){
}

int TPZMixedElasticityMaterial::ClassId() const {
    return Hash("TPZMixedElasticityMaterial") ^ TPZMaterial::ClassId() << 1;
}

void TPZMixedElasticityMaterial::Read(TPZStream &buf, void *context) {
    TPZMaterial::Read(buf, context);
    buf.Read(&fE_const, 1);
    buf.Read(&fnu_const, 1);
    fForce.Resize(2, 0.);
    buf.Read(&fForce[0], 2);
    buf.Read(&fPlaneStress, 1);
    buf.Read(&fPostProcIndex);
}

void TPZMixedElasticityMaterial::Write(TPZStream &buf, int withclassid) const {
    TPZMaterial::Write(buf, withclassid);
    buf.Write(&fE_const, 1);
    buf.Write(&fnu_const, 1);

    buf.Write(&fForce[0], 2);
    buf.Write(&fPlaneStress, 1);
    buf.Write(&fPostProcIndex);
}

