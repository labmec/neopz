#include "Electromagnetics/TPZAnisoWgma.h"
#include "TPZBndCondT.h"
#include "TPZCompElHCurl.h"
#include "pzaxestools.h"

using namespace std::complex_literals;
TPZAnisoWgma::TPZAnisoWgma() : TBase()
{
    SetMatrixK();
}
TPZAnisoWgma::TPZAnisoWgma(int id) : TBase(id){
    SetMatrixK();
}

TPZAnisoWgma::TPZAnisoWgma(int id, 
                           const CSTATE er,
                           const CSTATE ur,
                           const STATE lambda,
                           const REAL scale) :
    TBase(id), fScaleFactor(scale),
    fLambda(lambda)
{
    SetMatrixK();

    if (lambda <0){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"Setting negative wavelength. Aborting..\n";
        DebugStop();
    }
    SetPermittivity(er);
    SetPermeability(ur);
}

TPZAnisoWgma::TPZAnisoWgma(int id, 
                           const TPZFMatrix<CSTATE> &er,
                           const TPZFMatrix<CSTATE> &ur,
                           const STATE lambda,
                           const REAL scale) :
    TBase(id), fScaleFactor(scale),
    fLambda(lambda)
{
    SetMatrixK();
    
    if (lambda <0){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"Setting negative wavelength. Aborting..\n";
        DebugStop();
    }

    SetPermittivity(er);
    SetPermeability(ur);
    
}

TPZAnisoWgma* TPZAnisoWgma::NewMaterial() const{
    return new TPZAnisoWgma();
}


void TPZAnisoWgma::SetWavelength(STATE lambda)
{
    if (lambda <0){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"Setting negative wavelength. Aborting..\n";
        DebugStop();
    }
    fLambda = lambda;
}


void TPZAnisoWgma::SetPermeability(CSTATE ur)
{
    fUr.Redim(3,3);
    fUr.Put(0,0,ur);
    fUr.Put(1,1,ur);
    fUr.Put(2,2,ur);
}

void TPZAnisoWgma::SetPermeability(const TPZFMatrix<CSTATE>& ur)
{
    if(ur.Rows()!=3 || ur.Cols()!= 3){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"\nSize of ur != 3. Aborting...\n";
        DebugStop();
    }
    fUr = ur;
}
void TPZAnisoWgma::SetPermittivity(CSTATE er)
{
    fEr.Redim(3,3);
    fEr.Put(0,0,er);
    fEr.Put(1,1,er);
    fEr.Put(2,2,er);
}

void TPZAnisoWgma::SetPermittivity(const TPZFMatrix<CSTATE>& er)
{
    if(er.Rows()!=3 || er.Cols()!= 3){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"\nSize of er != 3. Aborting...\n";
        DebugStop();
    }
    fEr = er;
}

void TPZAnisoWgma::SetMatrixK()
{
    TPZMatQuadraticEigenVal::SetMatrixK();
    fCurrentContribute = [this](const auto &arg1, auto arg2, auto &arg3,
                                auto &arg4, const auto &arg5, const auto &arg6){
        this->ContributeK(arg1,arg2,arg3,arg4,arg5,arg6);
    };
    fCurrentContributeBC = [this](const auto &arg1, auto arg2, auto &arg3,
                                  auto &arg4, auto &arg5){
        this->ContributeBCK(arg1,arg2,arg3,arg4,arg5);
    };
}

void TPZAnisoWgma::SetMatrixL()
{
    TPZMatQuadraticEigenVal::SetMatrixL();
    fCurrentContribute = [this](const auto &arg1, auto arg2, auto &arg3,
                                auto &arg4, const auto &arg5, const auto &arg6){
        this->ContributeL(arg1,arg2,arg3,arg4,arg5,arg6);
    };
    fCurrentContributeBC = [this](const auto &arg1, auto arg2, auto &arg3,
                                  auto &arg4, auto &arg5){
        this->ContributeBCL(arg1,arg2,arg3,arg4,arg5);
    };
}


void TPZAnisoWgma::SetMatrixM()
{
    TPZMatQuadraticEigenVal::SetMatrixM();
    fCurrentContribute = [this](const auto &arg1, auto arg2, auto &arg3,
                                auto &arg4, const auto &arg5, const auto &arg6){
        this->ContributeM(arg1,arg2,arg3,arg4,arg5,arg6);
    };
    fCurrentContributeBC = [this](const auto &arg1, auto arg2, auto &arg3,
                                  auto &arg4, auto &arg5){
        this->ContributeBCM(arg1,arg2,arg3,arg4,arg5);
    };
}

void TPZAnisoWgma::GetPermittivity(
  [[maybe_unused]] const TPZVec<REAL> &x,TPZFMatrix<CSTATE> &er) const
{
    er = fEr;
}

void TPZAnisoWgma::GetPermeability(
  [[maybe_unused]] const TPZVec<REAL> &x,TPZFMatrix<CSTATE> &ur) const
{
    ur = fUr;
}

void
TPZAnisoWgma::Contribute(
    const TPZVec<TPZMaterialDataT<CSTATE>> &datavec,
    REAL weight,
    TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef)
{
    TPZFNMatrix<9,CSTATE> er(3,3,0.),ur(3,3,0.);
    GetPermittivity(datavec[0].x, er);
    GetPermeability(datavec[0].x, ur);
    ur.Decompose(ELU);
    fCurrentContribute(datavec,weight,ek,ef,er,ur);
}

void TPZAnisoWgma::ContributeBC(
    const TPZVec<TPZMaterialDataT<CSTATE>> &datavec,
    REAL weight,
    TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef,
    TPZBndCondT<CSTATE> &bc)
{
    fCurrentContributeBC(datavec,weight,ek,ef,bc);
}

void
TPZAnisoWgma::ContributeK(
    const TPZVec<TPZMaterialDataT<CSTATE>> &datavec,
    REAL weight,
    TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef,
    const TPZFMatrix<CSTATE> &er, const TPZFMatrix<CSTATE> &ur)
{
    const auto &phi_hcurl_real = datavec[fHCurlMeshIndex].phi;
    const auto &curl_phi_real = datavec[fHCurlMeshIndex].curlphi;
    const int nhcurl  = phi_hcurl_real.Rows();
    //making complex version of phi
    TPZFNMatrix<3000,CSTATE> phi_hcurl(3,nhcurl,0.);
    TPZFNMatrix<3000,CSTATE> rot_phi_hcurl(3,nhcurl,0.);
    TPZFNMatrix<3000,CSTATE> curl_phi(3,nhcurl,0.);
    for(int i = 0; i < nhcurl; i++){
        for(int x = 0; x < 2; x++){
            phi_hcurl.Put(x,i,phi_hcurl_real.Get(i,x));
        }
        rot_phi_hcurl.Put(0,i,phi_hcurl.Get(1,i));
        rot_phi_hcurl.Put(1,i,-1.*phi_hcurl.Get(0,i));
        curl_phi.Put(2,i,curl_phi_real.Get(0,i));
    }

    const auto &phi_h1_real = datavec[fH1MeshIndex].phi;
    const int nh1  = phi_h1_real.Rows();
    TPZFNMatrix<3000,REAL> grad_phi_real(3, nh1, 0.);
    {
        const TPZFMatrix<REAL> &gradPhiH1axes = datavec[fH1MeshIndex].dphix;
        TPZAxesTools<REAL>::Axes2XYZ(gradPhiH1axes, grad_phi_real, datavec[fH1MeshIndex].axes);
    }

    
    //making complex version of phi
    TPZFNMatrix<3000,CSTATE> phi_h1(3,nh1,0.);
    TPZFNMatrix<3000,CSTATE> rot_grad_phi(3,nh1,0.);
    for(int i = 0; i < nh1; i++){
        phi_h1.Put(2,i,phi_h1_real.Get(i,0));
        rot_grad_phi.Put(0,i,grad_phi_real.Get(1,i));
        rot_grad_phi.Put(1,i,-1.*grad_phi_real.Get(0,i));
    }
    
    
    const REAL k0 = fScaleFactor * 2*M_PI/fLambda;
    /*****************ACTUAL COMPUTATION OF CONTRIBUTION****************/
    const int firsthcurl = fHCurlMeshIndex * nh1;
    const int firsth1 = fH1MeshIndex * nhcurl;
    TPZFNMatrix<3000,CSTATE> tmp;

    //Att term
    tmp = curl_phi;
    ur.Substitution(&tmp);
    ek.AddContribution(firsthcurl, firsthcurl, tmp, true, curl_phi, false, weight);
    //Ctt term
    er.Multiply(phi_hcurl,tmp);
    ek.AddContribution(firsthcurl,firsthcurl,tmp,true,phi_hcurl,false, -k0*k0*weight);
    //Azz term
    tmp = rot_grad_phi;
    ur.Substitution(&tmp);
    ek.AddContribution(firsth1,firsth1,tmp,true,rot_grad_phi,false, weight);
    //Czz term
    er.Multiply(phi_h1,tmp);
    ek.AddContribution(firsth1,firsth1,tmp,true,phi_h1,false, -k0*k0*weight);
    //Ftz term
    tmp = curl_phi;
    ur.Substitution(&tmp);
    ek.AddContribution(firsthcurl,firsth1,tmp,true,rot_grad_phi,false,weight);
    //Fzt term
    tmp = rot_grad_phi;
    ur.Substitution(&tmp);
    ek.AddContribution(firsth1,firsthcurl,tmp,true, curl_phi, false,weight);
    //Gtz term
    er.Multiply(phi_hcurl,tmp);
    ek.AddContribution(firsthcurl, firsth1, tmp, true, phi_h1, false,-k0*k0*weight);
    //Gzt term
    er.Multiply(phi_h1,tmp);
    ek.AddContribution(firsth1, firsthcurl, tmp, true, phi_hcurl, false,-k0*k0*weight);
}

void
TPZAnisoWgma::ContributeL(
    const TPZVec<TPZMaterialDataT<CSTATE>> &datavec,
    REAL weight,
    TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef,
    const TPZFMatrix<CSTATE> &er, const TPZFMatrix<CSTATE> &ur)
{
    const auto &phi_hcurl_real = datavec[fHCurlMeshIndex].phi;
    const auto &curl_phi_real = datavec[fHCurlMeshIndex].curlphi;
    const int nhcurl  = phi_hcurl_real.Rows();
    //making complex version of phi
    TPZFNMatrix<3000,CSTATE> phi_hcurl(3,nhcurl,0.);
    TPZFNMatrix<3000,CSTATE> rot_phi_hcurl(3,nhcurl,0.);
    TPZFNMatrix<3000,CSTATE> curl_phi(3,nhcurl,0.);
    for(int i = 0; i < nhcurl; i++){
        for(int x = 0; x < 2; x++){
            phi_hcurl.Put(x,i,phi_hcurl_real.Get(i,x));
        }
        rot_phi_hcurl.Put(0,i,phi_hcurl.Get(1,i));
        rot_phi_hcurl.Put(1,i,-1.*phi_hcurl.Get(0,i));
        curl_phi.Put(2,i,curl_phi_real.Get(0,i));
    }

    const auto &phi_h1_real = datavec[fH1MeshIndex].phi;
    const int nh1  = phi_h1_real.Rows();
    TPZFNMatrix<3000,REAL> grad_phi_real(3, nh1, 0.);
    {
        const TPZFMatrix<REAL> &gradPhiH1axes = datavec[fH1MeshIndex].dphix;
        TPZAxesTools<REAL>::Axes2XYZ(gradPhiH1axes, grad_phi_real, datavec[fH1MeshIndex].axes);
    }

    
    TPZFNMatrix<3000,CSTATE> rot_grad_phi(3,nh1,0.);
    for(int i = 0; i < nh1; i++){
        rot_grad_phi.Put(0,i,grad_phi_real.Get(1,i));
        rot_grad_phi.Put(1,i,-1.*grad_phi_real.Get(0,i));
    }
    
    
    const REAL k0 = fScaleFactor * 2*M_PI/fLambda;
    /*****************ACTUAL COMPUTATION OF CONTRIBUTION****************/
    const int firsthcurl = fHCurlMeshIndex * nh1;
    const int firsth1 = fH1MeshIndex * nhcurl;
    TPZFNMatrix<3000,CSTATE> tmp;

    //Dtt term
    tmp = curl_phi;
    ur.Substitution(&tmp);
    ek.AddContribution(firsthcurl, firsthcurl, tmp, true, rot_phi_hcurl, false, -weight);
    tmp = rot_phi_hcurl;
    ur.Substitution(&tmp);
    ek.AddContribution(firsthcurl, firsthcurl, tmp, true, curl_phi, false, weight);
    //Atz term
    tmp = rot_phi_hcurl;
    ur.Substitution(&tmp);
    ek.AddContribution(firsthcurl,firsth1,tmp,true,rot_grad_phi,false, weight);
    //Azt term
    tmp = rot_grad_phi;
    ur.Substitution(&tmp);
    ek.AddContribution(firsth1,firsthcurl,tmp, true, rot_phi_hcurl, false, -weight);
}

void
TPZAnisoWgma::ContributeM(
    const TPZVec<TPZMaterialDataT<CSTATE>> &datavec,
    REAL weight,
    TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef,
    const TPZFMatrix<CSTATE> &er, const TPZFMatrix<CSTATE> &ur)
{
    const auto &phi_hcurl_real = datavec[fHCurlMeshIndex].phi;
    const auto &curl_phi_real = datavec[fHCurlMeshIndex].curlphi;
    const int nhcurl  = phi_hcurl_real.Rows();
    //making complex version of phi
    TPZFNMatrix<3000,CSTATE> rot_phi_hcurl(3,nhcurl,0.);
    for(int i = 0; i < nhcurl; i++){
        rot_phi_hcurl.Put(0,i,phi_hcurl_real.Get(i,1));
        rot_phi_hcurl.Put(1,i,-1.*phi_hcurl_real.Get(i,0));
    }

    /*****************ACTUAL COMPUTATION OF CONTRIBUTION****************/
    const int nh1  = datavec[fH1MeshIndex].phi.Rows();
    const int firsthcurl = fHCurlMeshIndex * nh1;
    const int firsth1 = fH1MeshIndex * nhcurl;
    TPZFNMatrix<3000,CSTATE> tmp;

    //Btt term
    tmp = rot_phi_hcurl;
    ur.Substitution(&tmp);
    ek.AddContribution(firsthcurl, firsthcurl, tmp, true, rot_phi_hcurl, false, -weight);
}


void TPZAnisoWgma::ContributeBCK(
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
        /// periodic conditions are treated at a mesh level
        break;
    default:
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"\nThis module supports only dirichlet and neumann boundary conditions.\n";
        PZError<<"Stopping now..."<<std::endl;
        DebugStop();
        break;
    }
}

int TPZAnisoWgma::VariableIndex(const std::string &name) const
{
    if( strcmp(name.c_str(), "Et_real") == 0) return 0;
    if( strcmp(name.c_str(), "Ez_real") == 0) return 1;
    if( strcmp(name.c_str(), "Et_abs") == 0) return 2;
    if( strcmp(name.c_str(), "Ez_abs") == 0) return 3;
    if( strcmp(name.c_str(), "Material") == 0) return 4;
    if( strcmp(name.c_str(), "P_H1") == 0) return 5;
    if( strcmp(name.c_str(), "P_HCurl") == 0) return 6;
    DebugStop();
    return 1;
}

int TPZAnisoWgma::NSolutionVariables(int var) const
{
    switch (var) {
        case 0: //Et_real
            return 2;
        case 1://Ez_real
            return 1;
        case 2: //Et_abs
            return 2;
        case 3://Ez_abs
            return 1;
        case 4://material
            return 2;
        case 5://pOrderH1
            return 1;
        case 6://pOrderHCurl
            return 1;
        default:
            DebugStop();
            break;
    }
    return 1;
}

/** @brief Returns the solution associated with the var index based on the finite element approximation */
void TPZAnisoWgma::Solution(
    const TPZVec<TPZMaterialDataT<CSTATE>> &datavec,
    int var,
    TPZVec<CSTATE> &solout)
{
    TPZFNMatrix<9,CSTATE> er(3,3,0.), ur(3,3,0.);

    GetPermittivity(datavec[0].x, er);
    GetPermeability(datavec[0].x, ur);
    
    TPZManVector<CSTATE,3> et(3,0.);
    TPZManVector<CSTATE,1> ez(1,0.);

    et = datavec[ fHCurlMeshIndex ].sol[0];
    ez = datavec[ fH1MeshIndex ].sol[0];

    switch (var) {
    case 0:{//et_real
        for (int i = 0; i < et.size(); ++i) {
            et[i] = std::real(et[i]);
        }
        solout = et;
        break;
    }
    case 1:{//ez_real
        for (int i = 0; i < ez.size(); ++i) {
            ez[i] = std::real(ez[i]);
        }
        solout = ez;
        break;
    }
    case 2:{//et_abs
        for (int i = 0; i < et.size(); ++i) {
            et[i] = std::abs(et[i]);
        }
        solout = et;
        break;
    }
    case 3:{//ez_abs
        for (int i = 0; i < ez.size(); ++i) {
            ez[i] = std::abs(ez[i]);
        }
        solout = ez;
        break;
    }
    case 4:{//material
        solout[0] = er.Get(0,0);
        solout[1] = er.Get(1,1);
        break;
    }
    case 5:{//pOrder
        solout.Resize(1);
        solout[0] = datavec[fH1MeshIndex].p;
        break;
    }
    case 6:{//pOrder
        solout.Resize(1);
        solout[0] = datavec[fHCurlMeshIndex].p;
        break;
    }
    default:
        DebugStop();
        break;
    }
}

#include "TPZCartesianPML.h"
template class TPZCombinedSpacesCartesianPML<TPZAnisoWgma>;