/**
 * @file
 * @brief Contains implementations of the TPZMaterialData methods.
 */

#include "pzmaterialdata.h"
#include "TPZMaterial.h"
#include "pzcompel.h"
#include "pzelmat.h"
#include <sstream>
#include "pzerror.h"
#include "TPZInterfaceEl.h"
#include "pzdiscgal.h"
#include "pzaxestools.h"
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.matrix.tpzfmatrix"));
static LoggerPtr loggerCheck(Logger::getLogger("pz.checkconsistency"));
#endif

TPZMaterialData::TPZMaterialData() : TPZRegisterClassId(&TPZMaterialData::ClassId), fShapeType(EEmpty),
    numberdualfunctions(0),normal(3,0.),x(3,0.),p(-1), fUserData(0){
    this->SetAllRequirements(false);
    this->fNeedsDeformedDirectionsFad = false;
    this->intLocPtIndex = -1;
    this->intGlobPtIndex = -1;
    this->NintPts = -1;
    this->sol.Resize(1);
    this->dsol.Resize(1);
    this->gelElId = -1;
    this->HSize = 0.;
    this->detjac = 0.;
    this->numberdualfunctions = 0;
    this->gelElId = -1;
    this->fMasterDirections = 0;
#ifdef _AUTODIFF
    this->fDeformedDirectionsFad = 0;
#endif
}

TPZMaterialData::TPZMaterialData( const TPZMaterialData &cp ) : 
TPZRegisterClassId(&TPZMaterialData::ClassId),
fShapeType(cp.fShapeType) {
    this->operator =(cp);
}

TPZMaterialData & TPZMaterialData::operator= (const TPZMaterialData &cp ){
    this->fShapeType = cp.fShapeType;
    this->fNeedsSol = cp.fNeedsSol;
    this->fNeedsNeighborSol = cp.fNeedsNeighborSol;
    this->fNeedsHSize = cp.fNeedsHSize;
    this->fNeedsNeighborCenter = cp.fNeedsNeighborCenter;
    this->fNeedsNormal = cp.fNeedsNormal;
    this->fNeedsDeformedDirectionsFad = cp.fNeedsDeformedDirectionsFad;
    this->phi = cp.phi;
    this->dphi = cp.dphi;
    this->dphix = cp.dphix;
    this->divphi = cp.divphi;
    this->curlphi = cp.curlphi;
    this->axes = cp.axes;
    this->jacobian = cp.jacobian;
    this->jacinv = cp.jacinv;
    this->normal = cp.normal;
    this->x = cp.x;
    this->p = cp.p;
    this->sol = cp.sol;
    this->dsol = cp.dsol;
    this->divsol = cp.divsol;
    this->curlsol = cp.curlsol;
    this->HSize = cp.HSize;
    this->detjac = cp.detjac;
    this->intLocPtIndex = cp.intLocPtIndex;
    this->intGlobPtIndex = cp.intGlobPtIndex;
    this->NintPts = cp.NintPts;
    this->XCenter = cp.XCenter;
    this->fMasterDirections = cp.fMasterDirections;
    this->fVecShapeIndex = cp.fVecShapeIndex;
    this->fDeformedDirections = cp.fDeformedDirections;
#ifdef _AUTODIFF
    this->fDeformedDirectionsFad = cp.fDeformedDirectionsFad;
#endif
    this->numberdualfunctions = cp.numberdualfunctions;
    this->gelElId = cp.gelElId;
    this->fUserData = cp.fUserData;
    return *this;
}

TPZMaterialData::~TPZMaterialData(){
    if(fUserData)
    {
        std::cout << "User data should be deleted and data set to zero before the destructor\n";
        DebugStop();
    }
    //NOTHING TO BE DONE!
}

void TPZMaterialData::SetAllRequirements(bool set){
    this->fNeedsSol = set;
    this->fNeedsNeighborSol = set;
    this->fNeedsHSize = set;
    this->fNeedsNeighborCenter = set;
    this->fNeedsNormal = set;
}

/** @brief Compare the object for identity with the object pointed to, eventually copy the object */
/** compare both objects bitwise for identity. Put an entry in the log file if different overwrite the calling object if the override flag is true */
bool TPZMaterialData::Compare(TPZSavable *copy, bool override)
{
    TPZMaterialData *comp = dynamic_cast<TPZMaterialData *>(copy);
    if(!comp) return false;
    bool result = true;
    bool locres;
    locres = phi.Compare(&comp->phi,override);
    if(!locres)
    {
        LOGPZ_DEBUG(loggerCheck,"phi different")
    }
    result = result && locres;
    locres = dphi.Compare(&comp->dphi,override);
    if(!locres)
    {
        LOGPZ_DEBUG(loggerCheck,"dphi different")
    }
    result = result && locres;
    locres = dphix.Compare(&comp->dphix,override);
    if(!locres)
    {
        LOGPZ_DEBUG(loggerCheck,"dphix different")
    }
    result = result && locres;
    locres = divphi.Compare(&comp->divphi,override);
    if(!locres)
    {
        LOGPZ_DEBUG(loggerCheck,"divphi different")
    }
    result = result && locres;
    locres = curlphi.Compare(&comp->curlphi,override);
    if(!locres)
    {
        LOGPZ_DEBUG(loggerCheck,"curlphi different")
    }
    result = result && locres;
    locres = axes.Compare(&comp->axes,override);
    if(!locres)
    {
        LOGPZ_DEBUG(loggerCheck,"axes different")
    }
    result = result && locres;
    locres = jacobian.Compare(&comp->jacobian,override);
    if(!locres)
    {
        LOGPZ_DEBUG(loggerCheck,"jacobian different")
    }
    result = result && locres;
    locres = jacinv.Compare(&comp->jacinv,override);
    if(!locres)
    {
        LOGPZ_DEBUG(loggerCheck,"jacinv different")
    }
    result = result && locres;
    return result;
}

// Compare the object for identity with the object pointed to, eventually copy the object
/*
 * compare both objects bitwise for identity. Put an entry in the log file if different
 * overwrite the calling object if the override flag is true
 */
bool TPZMaterialData::Compare(TPZSavable *copy, bool override) const
{
    DebugStop();
    return true;
}

/** Print the data */
void TPZMaterialData::Print(std::ostream &out) const
{
    out << "Shape function type " << ShapeFunctionType() << std::endl;
    out << "Active Approximation Space " << fActiveApproxSpace << std::endl;
    phi.Print("phi",out);
    dphi.Print("dphi",out);
    dphix.Print("dphix",out);
    out << "Number dual functions " << numberdualfunctions << std::endl;
    divphi.Print("div phi",out);
    curlphi.Print("curl phi",out);
    axes.Print("axes",out);
    jacobian.Print("jacobian",out);
    jacinv.Print("jacinv",out);
    out << "normal " << normal << std::endl;
    out << "x " << x << std::endl;
    out << "xParametric " << xParametric << std::endl;
    out << "p " << p << std::endl;
    out << "sol " << sol << std::endl;
    int nsol = dsol.size();
    for (int is=0; is<nsol; is++) {
        dsol[is].Print("dsol",out);
    }
    out << "divsol " << divsol << std::endl;
    out << "curlsol " << curlsol << std::endl;
    out << "HSize " << HSize << std::endl;
    out << "detjac " << detjac << std::endl;
    out << "XCenter " << XCenter << std::endl;
    out << "fMasterDirections" << fMasterDirections << std::endl;
    out << "fDeformedDirections" << fDeformedDirections << std::endl;
#ifdef _AUTODIFF
    if(fNeedsDeformedDirectionsFad){
        fDeformedDirectionsFad.Print(out);
    }
    else
    {
        out << "No need for directions FAD\n";
    }
#endif
    out << "gelElId " << gelElId << std::endl;
    if (fVecShapeIndex.size()) {
        out << "VecShapeIndex: ";
        for (int64_t i = 0; i < fVecShapeIndex.size(); i++) {
            out << fVecShapeIndex[i].first << '/' << fVecShapeIndex[i].second << ' ';
        }
        out << '\n';
    }
    out << "NintPts " << NintPts << std::endl;
    out << "intLocPtIndex " << intLocPtIndex << std::endl;
    out << "intGlobPtIndex " << intGlobPtIndex << std::endl;
    out << "NeedsSol " << fNeedsSol << std::endl;
    out << "fNeedsNeighborSol " << fNeedsNeighborSol << std::endl;
    out << "fNeedsHSize " << fNeedsHSize << std::endl;
    out << "fNeedsNeighborCenter " << fNeedsNeighborCenter << std::endl;
    out << "fNeedsDeformedDirectionsFad " << fNeedsDeformedDirectionsFad << std::endl;
    out << "fNeedsNormal " << fNeedsNormal << std::endl;
}

/** Print the data in a format suitable for Mathematica */
void TPZMaterialData::PrintMathematica(std::ostream &out) const
{
    phi.Print("phi = ",out,EMathematicaInput);
    dphi.Print("dphi = ",out,EMathematicaInput);
    dphix.Print("dphix = ",out,EMathematicaInput);
    axes.Print("axes = ",out,EMathematicaInput);
    jacobian.Print("jacobian = ",out,EMathematicaInput);
    jacinv.Print("jacinv = ",out,EMathematicaInput);
    out << "normal = {" << normal << "};" << std::endl;
    out << "x = {" << x << "};" << std::endl;
    out << "p = " << p << ";" << std::endl;
    out << "sol = { " << sol << "};" << std::endl;
    int nsol=dsol.size();
    for (int is=0; is<nsol; is++) {
        std::stringstream sout;
        sout << "dsol" << is << " = ";
        dsol[is].Print(sout.str().c_str(),out,EMathematicaInput);
    }
    out << "divsol = { " << divsol << "};" << std::endl;
    out << "curlsol = { " << curlsol << "};" << std::endl;

    out << "HSize = " << HSize << ";" << std::endl;
    out << "detjac = " << detjac << ";" << std::endl;
    out << "XCenter = {" << XCenter << "};" << std::endl;
    out << "fMasterDirections" << fMasterDirections << std::endl;
    out << "intLocPtIndex = " << intLocPtIndex << ";" <<std::endl;
    out << "intGlobPtIndex = " << intGlobPtIndex << ";" <<std::endl;
    out << "NintPts = " << NintPts << ";" <<std::endl;
    out << "gelElId = " << gelElId << ";" <<std::endl;
}

/** Save the element data to a stream */
void TPZMaterialData::Write(TPZStream &buf, int withclassid) const
{
    int shapetype = fShapeType;
    buf.Write(&shapetype);
    phi.Write(buf,0);
    dphi.Write(buf,0);
    dphix.Write(buf,0);
    divphi.Write(buf, 0);
    curlphi.Write(buf, 0);
    axes.Write(buf,0);
    jacobian.Write(buf,0);
    jacinv.Write(buf,0);
    buf.Write(normal.begin(),normal.size());
    buf.Write(x.begin(),x.size());
    buf.Write(&p,1);
    int nsol = sol.size();
    buf.Write(&nsol);
    for (int is=0; is<nsol; is++) {
        buf.Write(sol[is].begin(),sol[is].size());
    }
    
    nsol = dsol.size();
    buf.Write(&nsol);
    for (int is=0; is<nsol; is++) {
        dsol[is].Write(buf,0);
    }
    nsol = divsol.size();
    buf.Write(&nsol);
    for (int is=0; is<nsol; is++) {
        buf.Write(divsol[is].begin(),divsol[is].size());
    }

    nsol = curlsol.size();
    buf.Write(&nsol);
    for (int is=0; is<nsol; is++) {
        buf.Write(curlsol[is].begin(),curlsol[is].size());
    }

    buf.Write(&HSize,1);
    buf.Write(&detjac,1);
    buf.Write(XCenter.begin(),XCenter.size());
    buf.Write(&intLocPtIndex,1);
    buf.Write(&intGlobPtIndex,1);
    buf.Write(&NintPts,1);
    buf.Write(&gelElId,1);
}

/** Read the element data from a stream */
void TPZMaterialData::Read(TPZStream &buf, void *context)
{
    int shapetype;
    buf.Read(&shapetype);
    fShapeType = (MShapeFunctionType) shapetype;
    phi.Read(buf,0);
    dphi.Read(buf,0);
    dphix.Read(buf,0);
    divphi.Read(buf, 0);
    curlphi.Read(buf, 0);
    axes.Read(buf,0);
    jacobian.Read(buf,0);
    jacinv.Read(buf,0);
    buf.Read(normal);
    buf.Read(x);
    buf.Read(&p,1);
    int nsol;
    buf.Read(&nsol,1);
    for (int is=0; is<nsol; is++) {
        buf.Read(sol[is]);
    }
    buf.Read(&nsol,1);
    for (int is = 0; is<nsol; is++) {
        dsol[is].Read(buf,0);
    }
    buf.Read(&nsol,1);
    for (int is=0; is<nsol; is++) {
        buf.Read(divsol[is]);
    }

    buf.Read(&nsol,1);
    for (int is=0; is<nsol; is++) {
        buf.Read(curlsol[is]);
    }

    buf.Read(&HSize,1);
    buf.Read(&detjac,1);
    buf.Read(XCenter);
    buf.Read(&intLocPtIndex,1);
    buf.Read(&intGlobPtIndex,1);
    buf.Read(&NintPts,1);
    buf.Read(&gelElId,1);
}

int TPZMaterialData::ClassId() const{
    return Hash("TPZMaterialData");
}

/// Computes the flux values based on a Material of Hdiv approx space
void TPZMaterialData::ComputeFluxValues(TPZFMatrix<REAL> & fluxes){
    
    if (fShapeType != EVecandShape) {
        std::cout << __PRETTY_FUNCTION__ << "works only for Vec and Shape type." << std::endl;
        return;
    }
    
    int n_shape = fVecShapeIndex.size();
    fluxes.Redim(3, n_shape);
    
    for (int i = 0; i < n_shape; i++) {
        int i_v = fVecShapeIndex[i].first;
        int i_s = fVecShapeIndex[i].second;
        
        for (int j = 0; j < 3; j++) {
            fluxes(j,i) = phi(i_s,0) * fDeformedDirections(j,i_v);
        }
    }
    
}

/// Compute the divergence of the shape functions
void TPZMaterialData::ComputeFunctionDivergence()
{
    
    // Getting test and basis functions
    TPZFMatrix<REAL> dphi_s       = dphi; // Derivative For H1  test functions
    
    int n_phi_v = fVecShapeIndex.NElements();
#ifdef PZDEBUG
    if(divphi.Rows() < n_phi_v) DebugStop();
#endif
    REAL det_jac = detjac;

    int i_vec = 0;
    int i_phi_s = 0;
    
    for (int iq = 0; iq < n_phi_v; iq++)
    {
        i_vec = fVecShapeIndex[iq].first;
        i_phi_s = fVecShapeIndex[iq].second;
        divphi(iq,0) = 0.;

        int n_dir = dphi_s.Rows();
        divphi(iq,0) = 0.;
        for (int k = 0; k < n_dir; k++) {
            divphi(iq,0) +=  dphi(k,i_phi_s)*fMasterDirections(k,i_vec)/detjac;
        }
    }
        

}

/// Shape function type as a string
std::string TPZMaterialData::ShapeFunctionType() const
{
    switch(fShapeType){
        case EEmpty:
            return "NotInitialized";
        case EScalarShape:
            return "Scalar";
        case EVecandShape:
            return "Vector combined with Scalar";
        case EVecShape:
            return "Vector shape";
        default:
            DebugStop();
    }
    return "All Wrong!\n";
}



template class TPZRestoreClass<TPZMaterialData>;
