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

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.matrix.tpzfmatrix"));
static LoggerPtr loggerCheck(Logger::getLogger("pz.checkconsistency"));
#endif

TPZMaterialData::TPZMaterialData() : TPZRegisterClassId(&TPZMaterialData::ClassId), fShapeType(EEmpty), numberdualfunctions(0){
    this->SetAllRequirements(false);
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
    this->phi = cp.phi;
    this->dphi = cp.dphi;
    this->dphix = cp.dphix;
    this->axes = cp.axes;
    this->jacobian = cp.jacobian;
    this->jacinv = cp.jacinv;
    this->normal = cp.normal;
    this->x = cp.x;
    this->p = cp.p;
    this->sol = cp.sol;
    this->dsol = cp.dsol;
    this->HSize = cp.HSize;
    this->detjac = cp.detjac;
    this->intLocPtIndex = cp.intLocPtIndex;
    this->intGlobPtIndex = cp.intGlobPtIndex;
    this->NintPts = cp.NintPts;
    this->XCenter = cp.XCenter;
    this->fVecShapeIndex = cp.fVecShapeIndex;
    this->fNormalVec = cp.fNormalVec;
    this->numberdualfunctions = cp.numberdualfunctions;
    this->gelElId = cp.gelElId;
    
    return *this;
}

TPZMaterialData::~TPZMaterialData(){
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
    phi.Print("phi",out);
    dphi.Print("dphi",out);
    dphix.Print("dphix",out);
    axes.Print("axes",out);
    jacobian.Print("jacobian",out);
    jacinv.Print("jacinv",out);
    out << "normal " << normal << std::endl;
    out << "x " << x << std::endl;
    out << "p " << p << std::endl;
    out << "sol " << sol << std::endl;
    int nsol = dsol.size();
    for (int is=0; is<nsol; is++) {
        dsol[is].Print("dsol",out);
        
    }
    
    out << "HSize " << HSize << std::endl;
    out << "detjac " << detjac << std::endl;
    out << "XCenter " << XCenter << std::endl;
    out << "intLocPtIndex " << intLocPtIndex << std::endl;
    out << "intGlobPtIndex " << intGlobPtIndex << std::endl;
    out << "NintPts " << NintPts << std::endl;
    out << "gelElId " << gelElId << std::endl;
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
    
    out << "HSize = " << HSize << ";" << std::endl;
    out << "detjac = " << detjac << ";" << std::endl;
    out << "XCenter = {" << XCenter << "};" << std::endl;
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

/** @brief Computes the flux values based on a Material of Hdiv approx space */
// @TODO:: Implement a method that computes the divergence of the fluxes
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
            fluxes(j,i) = phi(i_s,0) * fNormalVec(j,i_v);
        }
    }
    
}

template class TPZRestoreClass<TPZMaterialData>;
