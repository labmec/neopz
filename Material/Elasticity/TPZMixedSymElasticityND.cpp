/**
 * @file
 * @brief Contains implementations of the TPZMixedSymElasticityND methods.
 */

#include "TPZMixedSymElasticityND.h"
#include "pzelmat.h"
#include "TPZBndCondT.h"
#include "pzaxestools.h"
#include "TPZMatWithMem.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include "TPZMaterialDataT.h"
#include <math.h>

#include "pzlog.h"
//#include "meshgen.h"
#include "TPZAnalyticSolution.h"
#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.mixedelasticity"));
#endif

#include <fstream>
using namespace std;

TPZMixedSymElasticityND::TPZMixedSymElasticityND() : TBase() {
    fE_const = -1.; // Young modulus
    fnu_const = -1.; // Poisson coefficient
    fForce[0] = 0.; // X component of the body force
    fForce[1] = 0.; // Y component of the body force
    fForce[2] = 0.; // Z component of the body force - not used for this class
    flambda_const = 0.;
    fmu_const = 0.;

    fPlaneStress = 0;
    fMatrixA = 0.;
}

TPZMixedSymElasticityND::TPZMixedSymElasticityND(int id, int dimension) : TBase(id), fDimension(dimension) {
    fE_const = -1.; // Young modulus
    fnu_const = -1.; // Poisson coefficient
    fForce[0] = 0.; // X component of the body force
    fForce[1] = 0.; // Y component of the body force
    fForce[2] = 0.; // Z component of the body force - not used for this class
    flambda_const = 0.;
    fmu_const = 0.;

    fPlaneStress = 0;
    fMatrixA = 0.;
}

TPZMixedSymElasticityND::TPZMixedSymElasticityND(int id, REAL E, REAL nu, REAL fx, REAL fy, int planestress, int dimension) : TBase(id), fDimension(dimension) {
    this->SetElasticity(E, nu);
    fForce[0] = fx; // X component of the body force
    fForce[1] = fy; // Y component of the body force
    fForce[2] = 0.;
    fPlaneStress = planestress;
}

TPZMixedSymElasticityND::~TPZMixedSymElasticityND() {
}

void TPZMixedSymElasticityND::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE> > &datavec) const {
    auto nref = datavec.size();
    for (int i = 0; i < nref; i++) {
        datavec[i].SetAllRequirements(false);
        datavec[i].fNeedsNeighborSol = false;
        datavec[i].fNeedsNeighborCenter = false;
        datavec[i].fNeedsNormal = false;
    }
}

int TPZMixedSymElasticityND::NStateVariables() const {
    return fDimension;
}

////////////////////////////////////////////////////////////////////

// Divergence on deformed element


void TPZMixedSymElasticityND::ElasticityModulusTensor(TPZFMatrix<STATE> &MatrixElast, TElasticityAtPoint &elast) {
    //Matrix modulus Voigt notation:
    int matdim = (fDimension == 3) ? 9 : 4;
    MatrixElast.Redim(matdim, matdim);
    REAL mu = elast.fmu;
    REAL lambda = elast.flambda;
    if(fDimension == 2)
    {
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
//            MatrixElast(Exx, Exx) = (1 - elast.fnu*elast.fnu) / elast.fE;
//            MatrixElast(Exx, Eyy) = -elast.fnu * (1 + elast.fnu) / elast.fE;
//            MatrixElast(Eyy, Exx) = MatrixElast(Exx, Eyy);
//            MatrixElast(Eyy, Eyy) = MatrixElast(Exx, Exx);
//            MatrixElast(Exy, Exy) = 1. / (2. * mu);
//            MatrixElast(Eyx, Eyx) = MatrixElast(Exy, Exy);
            MatrixElast(Exx, Exx) = 1. / elast.fE;
            MatrixElast(Exx, Eyy) = -elast.fnu / elast.fE;
            MatrixElast(Eyy, Exx) = MatrixElast(Exx, Eyy);
            MatrixElast(Eyy, Eyy) = MatrixElast(Exx, Exx);
            MatrixElast(Exy, Exy) = (1.+elast.fnu) / elast.fE;
            MatrixElast(Eyx, Eyx) = MatrixElast(Exy, Exy);
        }
    } else if (fDimension == 3)
    {
        // sig = lambda tr(E) I + 2 mu E
        // E = 1/(2 mu) (sig - lambda /(3 lambda + 2 mu) tr(sig) I)
		// NOTE: NS double-checked and this is correct
        MatrixElast(Exx, Exx) = 1./(2.*mu)*(1.-lambda/(3.*lambda+2.*mu));
        MatrixElast(Exx, Eyy) = 1./(2.*mu)*(-lambda/(3.*lambda+2.*mu));
        MatrixElast(Exx, Ezz) = MatrixElast(Exx, Eyy);
        MatrixElast(Eyy, Eyy) = MatrixElast(Exx, Exx);
        MatrixElast(Eyy, Exx) = MatrixElast(Exx, Eyy);
        MatrixElast(Eyy, Ezz) = MatrixElast(Exx, Eyy);
        MatrixElast(Ezz, Ezz) = MatrixElast(Exx, Exx);
        MatrixElast(Ezz, Exx) = MatrixElast(Exx, Eyy);
        MatrixElast(Ezz, Eyy) = MatrixElast(Exx, Eyy);
        MatrixElast(Exy, Exy) = 1. / (2. * mu);
        MatrixElast(Exz, Exz) = MatrixElast(Exy, Exy);
        MatrixElast(Eyx, Eyx) = MatrixElast(Exy, Exy);
        MatrixElast(Eyz, Eyz) = MatrixElast(Exy, Exy);
        MatrixElast(Ezx, Ezx) = MatrixElast(Exy, Exy);
        MatrixElast(Ezy, Ezy) = MatrixElast(Exy, Exy);
    }
    else
    {
        DebugStop();
    }
#ifdef LOG4CXX
    if(logdata->isDebugEnabled()){
        std::stringstream sout;
        MatrixElast.Print("A",sout,EMathematicaInput);
        LOGPZ_DEBUG(logdata,sout.str());
    }
#endif //LOG4CXX
}

void TPZMixedSymElasticityND::ComputeDeformationVector(TPZVec<STATE> &PhiStress, TPZVec<STATE> &APhiStress, TElasticityAtPoint &elast) {
    int matdim = (fDimension == 3) ? 9 : 4;
    TPZFNMatrix<81, STATE> MatrixElast(matdim, matdim, 0.);
    ElasticityModulusTensor(MatrixElast, elast);
    for (int iq = 0; iq < matdim; iq++) {
        APhiStress[iq] = 0.;
        for (int jq = 0; jq < matdim; jq++) {
            APhiStress[iq] += MatrixElast(iq, jq) * PhiStress[jq];
        }
    }
}

void TPZMixedSymElasticityND::ComputeStressVector(TPZVec<STATE> &Deformation, TPZVec<STATE> &Stress, TElasticityAtPoint &elast) {
    if(fDimension == 2)
    {
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
    } else if(fDimension == 3)
    {
        REAL mu = elast.fmu;
        REAL lambda = elast.flambda;
        STATE trE = Deformation[Exx]+Deformation[Eyy]+Deformation[Ezz];
        Stress[Exx] = lambda*trE+2.*mu*Deformation[Exx];
        Stress[Eyy] = lambda*trE+2.*mu*Deformation[Eyy];
        Stress[Ezz] = lambda*trE+2.*mu*Deformation[Ezz];
        Stress[Exy] = 2.*mu*Deformation[Exy];
        Stress[Exz] = 2.*mu*Deformation[Exz];
        Stress[Eyx] = 2.*mu*Deformation[Eyx];
        Stress[Eyz] = 2.*mu*Deformation[Eyz];
        Stress[Ezx] = 2.*mu*Deformation[Ezx];
        Stress[Ezy] = 2.*mu*Deformation[Ezy];
    }
}

void TPZMixedSymElasticityND::Print(std::ostream &out) const {
    out << "name of material : " << Name() << "\n";
    out << "properties : \n";
    out << "\tE_const   = " << fE_const << endl;
    out << "\tnu_const   = " << fnu_const << endl;
    out << "\tF   = " << fForce[0] << ' ' << fForce[1] << endl;
}

/// Transform a tensor to a Voigt notation

void TPZMixedSymElasticityND::ToVoigt(const TPZFMatrix<STATE> &S, TPZVec<STATE> &Svoigt) const {
    Svoigt[Exx] = S(0, 0);
    Svoigt[Exy] = S(0, 1);
    Svoigt[Eyx] = S(1, 0);
    Svoigt[Eyy] = S(1, 1);
    if (fDimension > 2) {
        Svoigt[Exz] = S(0, 2);
        Svoigt[Eyz] = S(1, 2);
        Svoigt[Ezx] = S(2, 0);
        Svoigt[Ezy] = S(2, 1);
        Svoigt[Ezz] = S(2, 2);
    }
}

/// Transform a Voigt notation to a tensor

void TPZMixedSymElasticityND::FromVoigt(const TPZVec<STATE> &Svoigt, TPZFMatrix<STATE> &S) const {
    S(0, 0) = Svoigt[Exx];
    S(0, 1) = Svoigt[Exy];
    S(1, 0) = Svoigt[Eyx];
    S(1, 1) = Svoigt[Eyy];
    if (fDimension > 2) {
        S(0, 2) = Svoigt[Exz];
        S(1, 2) = Svoigt[Eyz];
        S(2, 0) = Svoigt[Ezx];
        S(2, 1) = Svoigt[Ezy];
        S(2, 2) = Svoigt[Ezz];
    }
}

/** @brief Calculates the element stiffness matrix using 3 spaces - Stress tensor, displacement, and skew-symmetric tensor (for weak symmetry) */
void TPZMixedSymElasticityND::Contribute_2spaces(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    REAL R = datavec[0].x[0];
    // Setting the phi's
    // E
    TPZFMatrix<REAL> &phiS = datavec[0].phi;
    TPZFMatrix<REAL> &dphiS = datavec[0].dphix;
    // U
    TPZFMatrix<REAL> &phiU = datavec[1].phi;
//    TPZFMatrix<REAL> &dphiU = datavec[1].dphix;

    TElasticityAtPoint elast(fE_const, fnu_const);
    if (fElasticity) {
        TPZManVector<STATE, 3> result(2);
        TPZFNMatrix<4, STATE> Dres(0, 0);
        fElasticity(datavec[0].x, result, Dres);
        REAL E = result[0];
        REAL nu = result[1];
        elast = TElasticityAtPoint(E, nu);
    }

    int64_t nshapeS, nshapeU;
    int nstresses = (fDimension*(fDimension+1))/2;
    nshapeS = datavec[0].phi.Rows();
    nshapeU = datavec[1].phi.Rows();
    const int firstequation_S = 0;
    const int firstequation_U = firstequation_S + nshapeS*(fDimension*(fDimension+1))/2;
    //const int firstequation_P = firstequation_U + nshapeU*fDimension;

    int voigtdim = fDimension*fDimension;
    // matrix representing the stress associated with each shape function
    TPZFMatrix<REAL> phistress(voigtdim,nshapeS*nstresses,0.);

    for (int is = 0; is<nshapeS; is++) {
        phistress(Exx,is*nstresses) = phiS(is,0);
        phistress(Exy,is*nstresses) = 0.;
        phistress(Eyy,is*nstresses) = 0.;
        phistress(Eyx,is*nstresses) = 0.;
        phistress(Eyy,is*nstresses+1) = phiS(is,0);
        phistress(Exy,is*nstresses+1) = 0.;
        phistress(Exx,is*nstresses+1) = 0.;
        phistress(Eyx,is*nstresses+1) = 0.;
        phistress(Eyy,is*nstresses+2) = 0.;
        phistress(Exy,is*nstresses+2) = phiS(is,0);
        phistress(Exx,is*nstresses+2) = 0.;
        phistress(Eyx,is*nstresses+2) = phiS(is,0);
        if(Dimension() == 3) {
            std::cout << "Please implement me\n";
            DebugStop();
        }
    }

    TPZFMatrix<REAL> dphix(3,nshapeS);
    TPZAxesTools<REAL>::Axes2XYZ(dphiS, dphix, datavec[0].axes);
    TPZFMatrix<REAL> divstress(Dimension(),nshapeS*nstresses,0.);
    for (int is = 0; is<nshapeS; is++) {
        divstress(0,is*nstresses) = dphix(0,is);
        divstress(1,is*nstresses+1) = dphix(1,is);
        divstress(0,is*nstresses+2) = dphix(1,is);
        divstress(1,is*nstresses+2) = dphix(0,is);
        if(Dimension() == 3) {
            std::cout << "Please implement me\n";
            DebugStop();
        }

    }

    TPZManVector<STATE, 3>  force(fDimension, 0.);


    force = fForce;
    if (this->HasForcingFunction()) {
        fForcingFunction(datavec[0].x, force);
#ifdef LOG4CXX
        if (logdata->isDebugEnabled()) {
            std::stringstream sout;
            sout << " x = " << datavec[0].x << " force = " << force << std::endl;
            LOGPZ_DEBUG(logdata, sout.str())
        }
#endif
    }
    TPZFMatrix<STATE> forceMat(1,fDimension,0.);
    for(int id = 0; id<fDimension; id++) forceMat(0,id) = force[id];
    TPZFNMatrix<36,REAL> matrixelast(voigtdim, voigtdim,0.);
    ElasticityModulusTensor(matrixelast, elast);

    TPZFMatrix<REAL> phideform(voigtdim,nshapeS*nstresses,0.);
    matrixelast.Multiply(phistress, phideform);

    ek.AddContribution(0, 0, phistress, 1, phideform, 0, -weight);

    TPZFMatrix<REAL> disp(fDimension,nshapeU*fDimension);
    for (int is = 0; is<nshapeU; is++) {
        for (int d = 0; d<fDimension; d++) {
            disp(d,is*fDimension+d) = phiU(is,0);
        }
    }
    ef.AddContribution(nstresses*nshapeS, 0, disp, 1, forceMat, 1, weight);

    ek.AddContribution(0, nstresses*nshapeS, divstress, 1, disp, 0, -weight);
    ek.AddContribution(nstresses*nshapeS, 0, disp, 1, divstress, 0, -weight);
}




void TPZMixedSymElasticityND::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
    if (datavec[0].fVecShapeIndex.size() == 0) {
        FillVecShapeIndex(datavec[0]);
    }
    int64_t nspaces = datavec.size();
    if(nspaces != 2) DebugStop();
    Contribute_2spaces(datavec, weight, ek, ef);
}


void TPZMixedSymElasticityND::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) {
    if (datavec[0].phi.Rows() != 0 && datavec[0].fShapeType != TPZMaterialData::EScalarShape) {
        DebugStop();
    }
    int ndisp = datavec[1].phi.Rows();

    auto &v_2 = bc.Val2();
    auto &v_1 = bc.Val1();

    // Setting forcing function
    if (bc.HasForcingFunctionBC()) {
        TPZManVector<STATE, 3> res(fDimension);
        TPZFNMatrix<9, STATE> tens(fDimension, fDimension);
        bc.ForcingFunctionBC()(datavec[1].x, res, tens);
        STATE shear = (tens(0,1)+tens(1,0))/2.;
        tens(0,1) = shear;
        tens(1,0) = shear;
        v_2[0] = res[0];
        v_2[1] = res[1];
        if(fDimension == 3) v_2[2] = res[2];
        if(bc.Type() == 1) { // Neuman condition
            TElasticityAtPoint elast(fE_const, fnu_const);
            if (fElasticity) {
                TPZManVector<STATE, 3> result(2);
                TPZFNMatrix<4, STATE> Dres(0, 0);
                fElasticity(datavec[0].x, result, Dres);
                REAL E = result[0];
                REAL nu = result[1];
                elast = TElasticityAtPoint(E, nu);
            }

            TPZManVector<STATE,9> deform(4,0.),stressvec(4,0.);
            ToVoigt(tens, deform);
            ComputeStressVector(deform, stressvec, elast);
            TPZManVector<REAL,3> normal(2,0.);
            normal[0] = datavec[1].axes(0,1);
            normal[1] = -datavec[1].axes(0,0);
            v_2[0] = stressvec[Exx]*normal[0] + stressvec[Exy]*normal[1];
            v_2[1] = stressvec[Exy]*normal[0] + stressvec[Eyy]*normal[1];
        }
    }

    // Setting the phis
    // E
    TPZFMatrix<REAL> &phiS = datavec[0].phi;
    TPZFMatrix<REAL> &phiD = datavec[1].phi;

    int nshapeS;
    nshapeS = datavec[0].phi.Rows();


    //    TPZMaterialData::MShapeFunctionType shapetype = data.fShapeType;
    //    if(shapetype==data.EVecShape){
    //        ContributeVecShapeBC(data,weight,ek, ef,bc);
    //        return;
    //    }

    REAL R = datavec[0].x[0];
    if (R < 1.e-6) R = 1.e-6;

    int nstate = 2;
    if(fDimension == 3) nstate = 3;

    switch (bc.Type()) {
        case 0: // Dirichlet condition
        {
            for (int iq = 0; iq < nshapeS; iq++) {
                for (int idf = 0; idf < nstate; idf++) {
                    ef(nstate * iq + idf, 0) += v_2[idf] * phiS(iq, 0) * weight; // forced v2 displacement
                }
            }
            for (int iq = 0; iq < ndisp; iq++) {
                for (int jq = 0; jq < ndisp; jq++) {
                    for(int idf = 0; idf < nstate; idf++)
                    {
                        ek(nstate * iq + idf, nstate * jq + idf) += TPZMaterial::fBigNumber * phiD(iq, 0) * phiD(jq, 0) * weight;
                    }
                }
                for(int idf = 0; idf < nstate; idf++)
                {
                    ef(nstate * iq + idf, 0) += TPZMaterial::fBigNumber * v_2[idf] * phiD(iq, 0) * weight; // normal stress in x direction
                }
            }
        }
            break;

        case 1: // Neumann condition
        {
            for (int iq = 0; iq < ndisp; iq++) {
                for (int idf = 0; idf < nstate; idf++) {
                    ef(nstate * iq + idf, 0) += v_2[idf] * phiD(iq, 0) * weight; // forced v2 displacement
                }
            }
        }
            break;

       default:
            DebugStop();
            // nulo introduzindo o BIGNUMBER pelos valores da condição
    } // 1 Val1 : a leitura é 00 01 10 11

}


/** Returns the variable index associated with the name. */
int TPZMixedSymElasticityND::VariableIndex(const std::string &name) const {
    if (!strcmp("displacement", name.c_str())) return 9;
    if (!strcmp("Displacement", name.c_str())) return 9; //function U ***
    if (!strcmp("DisplacementMem", name.c_str())) return 9;
    if (!strcmp("Pressure", name.c_str())) return 1;
    if (!strcmp("MaxStress", name.c_str())) return 2;
    if (!strcmp("PrincipalStress1", name.c_str())) return 3;
    if (!strcmp("PrincipalStress2", name.c_str())) return 4;
    if (!strcmp("SigmaX", name.c_str())) return 5;
    if (!strcmp("SigmaY", name.c_str())) return 6;
    if (!strcmp("SigmaZ", name.c_str())) return 12;
    if (!strcmp("TauXY", name.c_str())) return 8;
    if (!strcmp("TauXZ", name.c_str())) return 30;
    if (!strcmp("TauYZ", name.c_str())) return 31;
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
    if (!strcmp("ExactDisplacement", name.c_str())) return 33;
    if (!strcmp("ExactStress", name.c_str())) return 34;
    if (!strcmp("ElementSigmaError", name.c_str())) return 100;

    return TPZMaterial::VariableIndex(name);
}

/** Returns the number of variables associated with the variable indexed by var. */
int TPZMixedSymElasticityND::NSolutionVariables(int var) const {
    int nstate = fDimension;
    switch (var) {
        case 0:
            return nstate;
        case 1:
        case 2:
            return 1;
        case 3:
        case 4:
            return nstate;
        case 5:
        case 6:
        case 8:
            return 1;
        case 7:
            return 6;
        case 9:
            return 3;
        case 10: //Stress Tensor
            return nstate*nstate;
        case 11: //Strain Tensor
            return nstate*nstate;
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
        case 27: // Rotation
            return Dimension() == 3 ? 3 : 1;
        case 28:
        case 29:
        case 30:
        case 31:
            return 1;
        case 33:
        case 34:
            return nstate;
        case 100:
            return 1;
        default:
            return TPZMaterial::NSolutionVariables(var);
    }
}


/** @brief Returns the solution associated with the var index based on the finite element approximation */
void TPZMixedSymElasticityND::Solution(const TPZVec<TPZMaterialDataT<STATE>> &data, int var, TPZVec<STATE> &Solout) {
#ifdef PZDEBUG
    // if (data.size() != 3) {
    //     DebugStop();
    // }
#endif
    TPZManVector<REAL, 3> x = data[0].x;
    TPZFNMatrix<9, STATE> sigma(3, 3, 0.), sigmah(3, 3, 0.), eps(3, 3, 0.);
    int dim = Dimension();
    sigma(0,0) = data[0].sol[0][ExxS];
    sigma(1,0) = data[0].sol[0][ExyS];
    sigma(0,1) = data[0].sol[0][ExyS];
    sigma(1,1) = data[0].sol[0][Eyy];
    TPZManVector<STATE, 3> disp(dim);
    for (int i = 0; i < dim; i++) {
        disp[i] = data[1].sol[0][i];
    }

    TPZFNMatrix<9, STATE> antisym(dim, dim, 0.);

    TElasticityAtPoint elast(fE_const,fnu_const);
    if(fElasticity)
    {
        //TPZManVector<REAL,3> result(2);
		TPZManVector<STATE, 3> result(dim);
        TPZFNMatrix<9,STATE> Dres(0,0);
        fElasticity(x, result, Dres);
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
	TPZManVector<STATE, 9> SIGMA(dim*dim, 0.), EPSZ(dim*dim, 0.);

    ToVoigt(sigma, SIGMA);

    ComputeDeformationVector(SIGMA, EPSZ, elast);

    FromVoigt(EPSZ, eps);

    if(dim ==2)
    {
        if (this->fPlaneStress == 1) {
            eps(2, 2) = -1. / E * nu * (sigma(0, 0) + sigma(1, 1));
        } else {
            sigma(2, 2) = nu * (sigma(0, 0) + sigma(1, 1));
            eps(2, 2) = 0;
        }
    }
    Pressure = -1/3. * (sigma(0, 0) + sigma(1, 1) + sigma(2,2));
    sigmah = sigma;
    sigmah(0, 0) += Pressure;
    sigmah(1, 1) += Pressure;
    sigmah(2, 2) += Pressure;
    // Displacement
    if (var == 9) {
        Solout[2] = 0.;
        for(int idf = 0; idf < dim; idf++)
        {
            Solout[idf] = disp[idf];
        }
        return;
    }
    // Pressure
    if (var == 1) {
        Solout[0] = Pressure;
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
    //  TauXZ
    if (var == 30) {
        Solout[0] = 0.5 * (sigma(0, 2) + sigma(2, 0));
        return;
    }
    //  TauYZ
    if (var == 31) {
        Solout[0] = 0.5 * (sigma(1, 2) + sigma(2, 1));
        return;
    }

    // Exact displacement
    if (var == 33) {
        TPZVec<STATE> u_exact(fDimension,0.);
        TPZFMatrix<STATE> du_exact(fDimension,fDimension,0.);
        if (this->fExactSol) {
            this->fExactSol(data[0].x, u_exact, du_exact);
        }
        for (int idf = 0; idf < dim; idf++) {
            Solout[idf] = u_exact[idf];
        }
        return;
    }
    // Exact stress
    if (var == 34) {

        TPZVec<STATE> u_exact(fDimension,0.);
        TPZFMatrix<STATE> du_exact(fDimension,fDimension,0.);
        if (this->fExactSol) {
            this->fExactSol(data[0].x, u_exact, du_exact);
        }
        // std::cout << "duexact = " << du_exact << std::endl;
        // For 2D only
        Solout[0] = du_exact(0,0);//Sigma x
        Solout[1] = du_exact(1,1);//Sigma y
        Solout[2] = (du_exact(1,0)+du_exact(0,1))/2;//Tau xy


        return;
    }




    //Strain
    if (var == 11) {
        Solout[Exx] = eps(0, 0);
        Solout[Eyx] = eps(1, 0);
        Solout[Exy] = eps(0, 1);
        Solout[Eyy] = eps(1, 1);
        if(dim == 3)
        {
            Solout[Exz] = eps(0, 2);
            Solout[Eyz] = eps(1, 2);
            Solout[Ezx] = eps(2, 0);
            Solout[Ezy] = eps(2, 1);
            Solout[Ezz] = eps(2, 2);
        }
        return;
    }

    //Stress
    if (var == 10) {
        Solout[Exx] = sigma(0, 0);
        Solout[Eyx] = sigma(1, 0);
        Solout[Exy] = sigma(0, 1);
        Solout[Eyy] = sigma(1, 1);
        if(dim == 3)
        {
            Solout[Exz] = sigma(0, 2);
            Solout[Eyz] = sigma(1, 2);
            Solout[Ezx] = sigma(2, 0);
            Solout[Ezy] = sigma(2, 1);
            Solout[Ezz] = sigma(2, 2);
        }
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
        if(dim == 3)
        {
            Solout[1] = antisym(0, 2);
            Solout[2] = antisym(1, 2);
        }
        return;
    }
}


////////////////////////////////////////////////////////////////////

STATE TPZMixedSymElasticityND::Inner(TPZFMatrix<STATE> &S, TPZFMatrix<STATE> &T) {
    //inner product of two tensors
#ifdef DEBUG
    if (S.Rows() != S.Cols() || T.Cols() != T.Rows() || S.Rows() != T.Rows()) {
        DebugStop();
    }
#endif

    STATE Val = 0;
    int sr = S.Rows();
    for (int i = 0; i < sr; i++) {
        for (int j = 0; j < sr; j++) {
            Val += S(i, j) * T(i, j);
        }
    }
    return Val;
}

////////////////////////////////////////////////////////////////////

template <typename TVar>
TVar TPZMixedSymElasticityND::InnerVec(const TPZVec<TVar> &S, const TPZVec<TVar> &T) {
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
STATE TPZMixedSymElasticityND::Tr(TPZFMatrix<REAL> &GradU) {
#ifdef DEBUG
    if (GradU.Rows() != GradU.Cols()) {
        DebugStop();
    }
#endif

    STATE Val = 0.;

    int grr = GradU.Rows();
    for (int i = 0; i < grr; i++) {
        Val += GradU(i, i);
    }

    return Val;
}


/// transform a H1 data structure to a vector data structure

void TPZMixedSymElasticityND::FillVecShapeIndex(TPZMaterialData &data) {
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


int TPZMixedSymElasticityND::NEvalErrors() const {
    return 5;
}

void TPZMixedSymElasticityND::Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors) {
    //values[0] = 0.;
    //TPZManVector<REAL, 4> SigmaV(4, 0.), sigma_exactV(4, 0.), eps_exactV(4, 0.), EPSZV(4, 0.);
    int dim = Dimension();
    if(dim != 2) DebugStop();
    TPZVec<STATE> u_exact(fDimension,0.);
    TPZFMatrix<STATE> du_exact(fDimension,fDimension,0.);
    TPZManVector<STATE> divsigma(fDimension,0.);
    if (this->fExactSol) {
        this->fExactSol(data[0].x, u_exact, du_exact);
    }
    if (this->fForcingFunction) {
        this->fForcingFunction(data[0].x, divsigma);
    }

    int nstate = fDimension;
    int matdim = nstate*nstate;
	TPZManVector<STATE, 9> SigmaV(matdim, 0.), sigma_exactV(matdim, 0.), eps_exactV(matdim, 0.), EPSZV(matdim, 0.);
	TPZFNMatrix<9, STATE> sigma(nstate, nstate, 0.), eps(nstate, nstate, 0.), grad(nstate, nstate, 0.);
    TPZFNMatrix<9, STATE> eps_exact(nstate, nstate, 0.);
    TPZManVector<REAL, 3> x = data[0].x;
    sigma(0, 0) = data[0].sol[0][ExxS];
    sigma(1, 0) = data[0].sol[0][ExyS];
    sigma(0, 1) = data[0].sol[0][ExyS];
    sigma(1, 1) = data[0].sol[0][EyyS];
    ToVoigt(sigma, SigmaV);
    TPZManVector<STATE, 3> divSigma(dim,0.);
    TPZFNMatrix<9,STATE> dsol = data[0].dsol[0];
    TPZFNMatrix<9,STATE> dsolxy(3, 3);
    TPZAxesTools<STATE>::Axes2XYZ(dsol, dsolxy, data[0].axes);
    divSigma[0] = dsolxy(0,ExxS)+dsolxy(1,ExyS);
    divSigma[1] = dsolxy(0,ExyS)+dsolxy(1,EyyS);

    TPZManVector<STATE, 3> disp(dim);
    for (int i = 0; i < dim; i++) {
        disp[i] = data[1].sol[0][i];
    }

#ifdef LOG4CXX
    if (logdata->isDebugEnabled())
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
    if(dim == 3)
    {
        eps_exact(2, 2) = du_exact(2, 2);
        eps_exact(2, 0) = eps_exact(0, 2) = 0.5 * (du_exact(0, 2) + du_exact(2, 0));
        eps_exact(1, 2) = eps_exact(2, 1) = 0.5 * (du_exact(2, 1) + du_exact(1, 2));
    }
    ToVoigt(eps_exact, eps_exactV);

    TElasticityAtPoint elast(fE_const,fnu_const);
    if(fElasticity)
    {
        //TPZManVector<REAL,3> result(2);
		TPZManVector<STATE, 3> result(2);
        TPZFNMatrix<4,STATE> Dres(0,0);
        fElasticity(x, result, Dres);
        REAL E = result[0];
        REAL nu = result[1];
        TElasticityAtPoint modify(E,nu);
        elast = modify;
    }

    ComputeStressVector(eps_exactV, sigma_exactV, elast);


#ifdef PZDEBUG
    {
        TPZManVector<STATE, 9> eps_again(matdim);
        ComputeDeformationVector(sigma_exactV, eps_again, elast);
        for (int i = 0; i < matdim; i++) {
            if (abs(eps_again[i] - eps_exactV[i]) > 1.e-8) {
                DebugStop();
            }
        }
    }
#endif

    //L2_Error: for the displacement (u)
    errors[3] = 0.;
    for (int idf = 0; idf < dim; idf++) {
        errors[3] += (disp[idf] - u_exact[idf])*(disp[idf] - u_exact[idf]);
    }
    //Energe norm
    //TPZManVector<REAL,4> SIGMA(4,0.) , EPSZ(4,0.);


    ToVoigt(sigma, SigmaV);

    ComputeDeformationVector(SigmaV, EPSZV, elast);

    errors[0] = 0.;
    errors[1] = 0.;
    for (int i = 0; i < matdim; i++) {
        //L2_Error: for the stress tensor (sigma)
        errors[0] += (SigmaV[i] - sigma_exactV[i])*(SigmaV[i] - sigma_exactV[i]);

        //Energy_Error: for the stress tensor (sigma)
        errors[1] += (SigmaV[i] - sigma_exactV[i])*(EPSZV[i] - eps_exactV[i]);
    }
    if (errors[1] < 0.) {
        std::cout << "I should stop \n";
    }

    TPZManVector<STATE,3> divSigmaExact(fDimension,0.);
    if(HasForcingFunction())
    {
        fForcingFunction(data[0].x, divSigmaExact);
        for(int i=0; i<fDimension; i++) divSigmaExact[i] *= -1.;
    }
    errors[2] = 0.;
    for (int idf = 0; idf < fDimension; idf++)
    {
        //L2_Error: divergent of the stress tensor (div(sigma))

        errors[2] += (divSigma[idf]-divSigmaExact[idf])*(divSigma[idf]-divSigmaExact[idf]);

    }

    //Energy_Norm: for the exact stress solution
    errors[4] = 0.;
    for(int i=0; i<matdim; i++)
    {
        errors[4] += eps_exactV[i]*sigma_exactV[i];
    }
}


TPZMixedSymElasticityND::TPZMixedSymElasticityND(const TPZMixedSymElasticityND &copy) :
TBase(copy),
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

int TPZMixedSymElasticityND::ClassId() const {
    return Hash("TPZMixedSymElasticityND") ^ TBase::ClassId() << 1;
}

void TPZMixedSymElasticityND::Read(TPZStream &buf, void *context) {
    TPZMaterial::Read(buf, context);
    buf.Read(&fE_const, 1);
    buf.Read(&fnu_const, 1);
    fForce.Resize(2, 0.);
    buf.Read(&fForce[0], 2);
    buf.Read(&fPlaneStress, 1);
}

void TPZMixedSymElasticityND::Write(TPZStream &buf, int withclassid) const {
    TPZMaterial::Write(buf, withclassid);
    buf.Write(&fE_const, 1);
    buf.Write(&fnu_const, 1);

    buf.Write(&fForce[0], 2);
    buf.Write(&fPlaneStress, 1);
}

