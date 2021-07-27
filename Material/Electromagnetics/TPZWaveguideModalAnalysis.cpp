#include "Electromagnetics/TPZWaveguideModalAnalysis.h"
#include "TPZBndCondT.h"
#include "TPZCompElHCurl.h"
#include "pzaxestools.h"

using namespace std::complex_literals;
TPZWaveguideModalAnalysis::TPZWaveguideModalAnalysis() : TBase()
{
    SetMatrixA();
}
TPZWaveguideModalAnalysis::TPZWaveguideModalAnalysis(int id) : TBase(id){
    SetMatrixA();
}

TPZWaveguideModalAnalysis::TPZWaveguideModalAnalysis(int id, 
                                                     const CSTATE ur,
                                                     const CSTATE er,
                                                     const STATE lambda,
                                                     const REAL &scale) :
    TBase(id), fUr(3), fEr(3), fScaleFactor(scale),
    fLambda(lambda)
{
    SetMatrixA();
    
    fUr = {ur,ur,ur};
    fEr = {er,er,er};
    if (lambda <0){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"Setting negative wavelength. Aborting..\n";
        DebugStop();
    }
    if (std::real(ur) <0){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"Setting negative permeability. Aborting..\n";
        DebugStop();
    }
    if (std::real(er) <0){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"Setting negative permeability. Aborting..\n";
        DebugStop();
    }
}

TPZWaveguideModalAnalysis::TPZWaveguideModalAnalysis(int id, 
                                                     const TPZVec<CSTATE> &ur,
                                                     const TPZVec<CSTATE> &er,
                                                     const STATE lambda,
                                                     const REAL &scale) :
    TBase(id), fUr(3), fEr(3), fScaleFactor(scale),
    fLambda(lambda)
{
    SetMatrixA();
    
    if (lambda <0){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"Setting negative wavelength. Aborting..\n";
        DebugStop();
    }

    SetPermeability(ur);
    SetPermittivity(er);
}

TPZWaveguideModalAnalysis* TPZWaveguideModalAnalysis::NewMaterial() const{
    return new TPZWaveguideModalAnalysis();
}


void TPZWaveguideModalAnalysis::SetWavelength(STATE lambda)
{
    if (lambda <0){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"Setting negative wavelength. Aborting..\n";
        DebugStop();
    }
    fLambda = lambda;
}


void TPZWaveguideModalAnalysis::SetPermeability(CSTATE ur)
{
    if (std::real(ur) <0){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"Setting negative permeability. Aborting..\n";
        DebugStop();
    }
    fUr = {ur,ur,ur};
}

void TPZWaveguideModalAnalysis::SetPermeability(const TPZVec<CSTATE>& ur)
{
    if(ur.size()!=3){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"\nSize of ur != 3. Aborting...\n";
        DebugStop();
    }
    for(const auto &iur : ur){
        if (std::real(iur) <0){
            PZError<<__PRETTY_FUNCTION__;
            PZError<<"Setting negative permeability. Aborting..\n";
            DebugStop();
        }
    }
    fUr = ur;
}
void TPZWaveguideModalAnalysis::SetPermittivity(CSTATE er)
{
    if (std::real(er) <0){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"Setting negative permeability. Aborting..\n";
        DebugStop();
    }
    fEr = {er,er,er};
}

void TPZWaveguideModalAnalysis::SetPermittivity(const TPZVec<CSTATE>&er)
{
    if(er.size()!=3){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"\nSize of er != 3. Aborting...\n";
        DebugStop();
    }
    for(const auto &ier : er){
        if (std::real(ier) <0){
            PZError<<__PRETTY_FUNCTION__;
            PZError<<"Setting negative permitivitty. Aborting..\n";
            DebugStop();
        }
    }
    fEr = er;
}

void TPZWaveguideModalAnalysis::SetMatrixA()
{
    TPZMatGeneralisedEigenVal::SetMatrixA();
    fCurrentContribute = [this](const auto &arg1, auto arg2, auto &arg3,
                                auto &arg4, const auto &arg5, const auto &arg6){
        this->ContributeA(arg1,arg2,arg3,arg4,arg5,arg6);
    };
    fCurrentContributeBC = [this](const auto &arg1, auto arg2, auto &arg3,
                                  auto &arg4, auto &arg5){
        this->ContributeBCA(arg1,arg2,arg3,arg4,arg5);
    };
}
void TPZWaveguideModalAnalysis::SetMatrixB()
{
    TPZMatGeneralisedEigenVal::SetMatrixB();
    fCurrentContribute = [this](const auto &arg1, auto arg2, auto &arg3,
                                auto &arg4, const auto &arg5, const auto &arg6){
        this->ContributeB(arg1,arg2,arg3,arg4, arg5, arg6);
    };
    fCurrentContributeBC = [this](const auto &arg1, auto arg2, auto &arg3,
                                  auto &arg4, auto &arg5){
        this->ContributeBCB(arg1,arg2,arg3,arg4,arg5);
    };
}

void TPZWaveguideModalAnalysis::GetPermittivity(
  [[maybe_unused]] const TPZVec<REAL> &x,TPZVec<CSTATE> &er) const
{
    er = fEr;
}

void TPZWaveguideModalAnalysis::GetPermeability(
  [[maybe_unused]] const TPZVec<REAL> &x,TPZVec<CSTATE> &ur) const
{
    ur = fUr;
}

void
TPZWaveguideModalAnalysis::Contribute(
    const TPZVec<TPZMaterialDataT<CSTATE>> &datavec,
    REAL weight,
    TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef)
{
    TPZManVector<CSTATE,3> er(3,0.),ur(3,0.);
    GetPermittivity(datavec[0].x, er);
    GetPermeability(datavec[0].x, ur);
    fCurrentContribute(datavec,weight,ek,ef,er,ur);
}

void TPZWaveguideModalAnalysis::ContributeBC(
    const TPZVec<TPZMaterialDataT<CSTATE>> &datavec,
    REAL weight,
    TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef,
    TPZBndCondT<CSTATE> &bc)
{
    fCurrentContributeBC(datavec,weight,ek,ef,bc);
}

void
TPZWaveguideModalAnalysis::ContributeA(
    const TPZVec<TPZMaterialDataT<CSTATE>> &datavec,
    REAL weight,
    TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef,
    const TPZVec<CSTATE> &er, const TPZVec<CSTATE> &ur)
{
    const CSTATE &exx = er[0];
    const CSTATE &eyy = er[1];
    const CSTATE &ezz = er[2];
    const CSTATE &uxx = ur[0];
    const CSTATE &uyy = ur[1];
    const CSTATE &uzz = ur[2];
    //get h1 functions for computing hcurl funcs
    const TPZFMatrix<REAL> &phiH1 = datavec[fH1MeshIndex].phi;
    
    TPZFNMatrix<30,REAL> phiHCurl;
    TPZHCurlAuxClass::ComputeShape(datavec[fHCurlMeshIndex].fVecShapeIndex,
                                   datavec[fHCurlMeshIndex].phi,
                                   datavec[fHCurlMeshIndex].fDeformedDirections,
                                   phiHCurl);
    
    auto & curlPhi = datavec[fHCurlMeshIndex].curlphi;
    const REAL k0 = fScaleFactor * 2*M_PI/fLambda;
    /*****************ACTUAL COMPUTATION OF CONTRIBUTION****************/

    const int nHCurlFunctions  = phiHCurl.Rows();
    const int nH1Functions  = phiH1.Rows();
    const int firstH1 = fH1MeshIndex * nHCurlFunctions;
    const int firstHCurl = fHCurlMeshIndex * nH1Functions;

    for (int iVec = 0; iVec < nHCurlFunctions; iVec++) {
        for (int jVec = 0; jVec < nHCurlFunctions; jVec++) {
            const STATE curlIzdotCurlJz =
                curlPhi(0 , iVec) * curlPhi(0 , jVec);
            STATE phiIdotPhiJx = phiHCurl(iVec , 0) * phiHCurl(jVec , 0);
            STATE phiIdotPhiJy = phiHCurl(iVec , 1) * phiHCurl(jVec , 1);

            CSTATE stiffAtt = 0.;
            stiffAtt = (1./uzz) * curlIzdotCurlJz;
            stiffAtt -= k0 * k0 * exx * phiIdotPhiJx;
            stiffAtt -= k0 * k0 * eyy * phiIdotPhiJy;
            ek( firstHCurl + iVec , firstHCurl + jVec ) += stiffAtt * weight ;
        }
    }
}

void
TPZWaveguideModalAnalysis::ContributeB(
    const TPZVec<TPZMaterialDataT<CSTATE>> &datavec,
    REAL weight,
    TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef,
    const TPZVec<CSTATE> &er, const TPZVec<CSTATE> &ur)
{
    const CSTATE &exx = er[0];
    const CSTATE &eyy = er[1];
    const CSTATE &ezz = er[2];
    const CSTATE &uxx = ur[0];
    const CSTATE &uyy = ur[1];
    const CSTATE &uzz = ur[2];
    //get h1 functions
    const TPZFMatrix<REAL> &phiH1 = datavec[fH1MeshIndex].phi;
    const TPZFMatrix<REAL> &gradPhiH1axes = datavec[fH1MeshIndex].dphix;
    TPZFNMatrix<3,REAL> gradPhiH1(3, phiH1.Rows(), 0.);
    TPZAxesTools<REAL>::Axes2XYZ(gradPhiH1axes, gradPhiH1, datavec[fH1MeshIndex].axes);
    
    TPZFNMatrix<30,REAL> phiHCurl;
    
    TPZHCurlAuxClass::ComputeShape(datavec[fHCurlMeshIndex].fVecShapeIndex,
                                   datavec[fHCurlMeshIndex].phi,
                                   datavec[fHCurlMeshIndex].fDeformedDirections,
                                   phiHCurl);
    
    const REAL k0 = fScaleFactor * 2*M_PI/fLambda;
    /*****************ACTUAL COMPUTATION OF CONTRIBUTION****************/

    const int nHCurlFunctions  = phiHCurl.Rows();
    const int nH1Functions  = phiH1.Rows();
    const int firstH1 = fH1MeshIndex * nHCurlFunctions;
    const int firstHCurl = fHCurlMeshIndex * nH1Functions;

    for (int iVec = 0; iVec < nHCurlFunctions; iVec++) {
        for (int jVec = 0; jVec < nHCurlFunctions; jVec++) {
            STATE phiIdotPhiJx = phiHCurl(iVec , 0) * phiHCurl(jVec , 0);
            STATE phiIdotPhiJy = phiHCurl(iVec , 1) * phiHCurl(jVec , 1);
            CSTATE stiffBtt = 0.;
            stiffBtt += (1./uyy) * phiIdotPhiJx;
            stiffBtt += (1./uxx) * phiIdotPhiJy;
            ek( firstHCurl + iVec , firstHCurl + jVec ) += stiffBtt * weight ;

        }
        for (int jSca = 0; jSca < nH1Functions; jSca++) {
            STATE phiVecDotGradPhiScax = phiHCurl(iVec , 0) * gradPhiH1(0,jSca);
            STATE phiVecDotGradPhiScay = phiHCurl(iVec , 1) * gradPhiH1(1,jSca);

            CSTATE stiffBzt = 0.;
            stiffBzt += (1./uyy) * phiVecDotGradPhiScax;
            stiffBzt += (1./uxx) * phiVecDotGradPhiScay;
            ek( firstHCurl + iVec , firstH1 + jSca ) += stiffBzt * weight ;
        }
    }
    for (int iSca = 0; iSca < nH1Functions; iSca++) {
        for (int jVec = 0; jVec < nHCurlFunctions; jVec++) {
            STATE phiVecDotGradPhiScax = phiHCurl(jVec , 0) * gradPhiH1(0,iSca);
            STATE phiVecDotGradPhiScay = phiHCurl(jVec , 1) * gradPhiH1(1,iSca);

            CSTATE stiffBtz = 0.;
            stiffBtz += (1./uyy) * phiVecDotGradPhiScax;
            stiffBtz += (1./uxx) * phiVecDotGradPhiScay;
            ek( firstH1 + iSca , firstHCurl +  jVec ) += stiffBtz * weight ;
        }
        for (int jSca = 0; jSca < nH1Functions; jSca++) {
            STATE gradPhiScaDotGradPhiScax = gradPhiH1(0,iSca) * gradPhiH1(0,jSca);
            STATE gradPhiScaDotGradPhiScay = gradPhiH1(1,iSca) * gradPhiH1(1,jSca);

            CSTATE stiffBzz = 0.;
            stiffBzz +=  (1./uyy) * gradPhiScaDotGradPhiScax;
            stiffBzz +=  (1./uxx) * gradPhiScaDotGradPhiScay;
            stiffBzz -=  k0 * k0 * ezz * phiH1( iSca , 0 ) * phiH1( jSca , 0 );

            ek( firstH1 + iSca , firstH1 + jSca) += stiffBzz * weight ;
        }
    }
}


void TPZWaveguideModalAnalysis::ContributeBCA(
    const TPZVec<TPZMaterialDataT<CSTATE>> &datavec,
    REAL weight,
    TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef,
    TPZBndCondT<CSTATE> &bc)
{
	
    const TPZFMatrix<REAL> &phiHCurl = datavec[fHCurlMeshIndex].phi;
    const TPZFMatrix<REAL> &phiH1 = datavec[fH1MeshIndex].phi;
    const int nHCurlFunctions  = phiHCurl.Rows();
    const int nH1Functions  = phiH1.Rows();
    const int firstH1 = fH1MeshIndex * nHCurlFunctions;
    const int firstHCurl = fHCurlMeshIndex * nH1Functions;
        
    
    const auto& BIG = TPZMaterial::fBigNumber;
    
    const CSTATE v1 = bc.Val1()(0,0);
    const CSTATE v2 = bc.Val2()[0];
    constexpr STATE tol = std::numeric_limits<STATE>::epsilon();
    if(std::abs(v2) > tol){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"\nThis method supports only homogeneous boundary conditions.\n";
        std::cout<<"Stopping now..."<<std::endl;
        DebugStop();
    }
    switch ( bc.Type() )
    {
    case 0:{
        int nshape=phiH1.Rows();
        for(int i = 0 ; i<nshape ; i++){
            for(int j=0;j<nshape;j++){
                const STATE stiff = phiH1(i,0) * phiH1(j,0) * BIG ;
                ek(firstH1+i,firstH1+j) += stiff*weight;
            }
        }
        nshape = nHCurlFunctions;
        for(int i = 0 ; i<nshape ; i++){
            for(int j=0;j<nshape;j++){
                const STATE stiff = phiHCurl(i,0) * phiHCurl(j,0) * BIG ;
                ek(firstHCurl+i,firstHCurl+j) += stiff*weight;
            }
        }
        break;
    }
    case 1:
        ///PMC condition just adds zero to both matrices. nothing to do here....
        break;
    case 2:
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"\nThis module supports only dirichlet and neumann boundary conditions.\n";
        PZError<<"Stopping now..."<<std::endl;
        DebugStop();
        break;
    }
}

void TPZWaveguideModalAnalysis::ContributeBCB(
    const TPZVec<TPZMaterialDataT<CSTATE>> &datavec,
    REAL weight,
    TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef,
    TPZBndCondT<CSTATE> &bc)
{}

int TPZWaveguideModalAnalysis::VariableIndex(const std::string &name) const
{
    if( strcmp(name.c_str(), "Et") == 0) return 0;
    if( strcmp(name.c_str(), "Ez") == 0) return 1;
    if( strcmp(name.c_str(), "Material") == 0) return 2;
    if( strcmp(name.c_str(), "POrderH1") == 0) return 3;
    if( strcmp(name.c_str(), "POrderHCurl") == 0) return 4;
    DebugStop();
    return 1;
}

int TPZWaveguideModalAnalysis::NSolutionVariables(int var) const
{
    switch (var) {
        case 0: //Et
            return 2;
            break;
        case 1://Ez
            return 1;
        case 2://material
            return 2;
        case 3://pOrderH1
            return 1;
        case 4://pOrderHCurl
            return 1;
        default:
            DebugStop();
            break;
    }
    return 1;
}

/** @brief Returns the solution associated with the var index based on the finite element approximation */
void TPZWaveguideModalAnalysis::Solution(
    const TPZVec<TPZMaterialDataT<CSTATE>> &datavec,
    int var,
    TPZVec<CSTATE> &solout)
{
    TPZManVector<CSTATE,3> er(3,0.), ur(3,0.);

    GetPermittivity(datavec[0].x, er);
    GetPermeability(datavec[0].x, ur);
    const CSTATE &exx = er[0];
    const CSTATE &eyy = er[1];
    const CSTATE &ezz = er[2];
    const CSTATE &uxx = ur[0];
    const CSTATE &uyy = ur[1];
    const CSTATE &uzz = ur[2];
    
    TPZManVector<CSTATE,3> et(3,0.);
    TPZManVector<CSTATE,1> ez(1,0.);

    et = datavec[ fHCurlMeshIndex ].sol[0];
    ez = datavec[ fH1MeshIndex ].sol[0];

    switch (var) {
        case 0:{//et
            for (int i = 0; i < et.size(); ++i) {
                et[i] /= fKz;
                et[i] = fPrintFieldRealPart ? std::real(et[i]) : std::abs(et[i]);
            }
            solout = et;
            break;
        }
        case 1:{//ez
            for (int i = 0; i < ez.size(); ++i) {
                ez[i] *= 1.0i;
                ez[i] = fPrintFieldRealPart ? std::real(ez[i]) : std::abs(ez[i]);
            }
            solout = ez;
            break;
        }

        case 2:{//material
            solout[0] = exx;
            solout[1] = eyy;
            break;
        }
        case 3:{//pOrder
            solout.Resize(1);
            solout[0] = datavec[fH1MeshIndex].p;
            break;
        }
        case 4:{//pOrder
            solout.Resize(1);
            solout[0] = datavec[fHCurlMeshIndex].p;
            break;
        }
        default:
            DebugStop();
            break;
    }
}