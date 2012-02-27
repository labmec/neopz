/**
 * \file
 * @brief Contains implementations of the TPZMaterialData methods.
 */
//$Id: pzmaterialdata.cpp,v 1.15 2011-05-11 02:18:08 phil Exp $ 

#include "pzmaterialdata.h"
#include "pzmaterial.h"
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

TPZMaterialData::TPZMaterialData() : numberdualfunctions(0){
	this->SetAllRequirements(false);
	this->intPtIndex = -1;
	//this->leftdetjac = 0.;
	//this->rightdetjac = 0.;
    this->sol.Resize(1);
    this->dsol.Resize(1);
	
}

TPZMaterialData::TPZMaterialData( const TPZMaterialData &cp ){
	this->operator =(cp);
}

TPZMaterialData & TPZMaterialData::operator= (const TPZMaterialData &cp ){
	this->fNeedsSol = cp.fNeedsSol;
	this->fNeedsNeighborSol = cp.fNeedsNeighborSol;
	this->fNeedsHSize = cp.fNeedsHSize;
	this->fNeedsNeighborCenter = cp.fNeedsNeighborCenter; 
	this->fNeedsNormal = cp.fNeedsNormal; 
	this->phi = cp.phi;
	//this-> phil = cp.phil;
	//this->phir = cp.phir;
	this->dphix = cp.dphix;
	//this->dphixl = cp.dphixl;
	//this->dphixr = cp.dphixr;
	this->axes = cp.axes;
	//this->axesleft = cp.axesleft;
	//this->axesright = cp.axesright;
	this->jacobian = cp.jacobian;
	//this->leftjac = cp.leftjac;
	//this->rightjac = cp.rightjac;
	this->jacinv = cp.jacinv;
	//this->leftjacinv = cp.leftjacinv;
	//this->rightjacinv = cp.rightjacinv;
	this->normal = cp.normal;
	this->x = cp.x;
	this->p = cp.p;
	//this->leftp = cp.leftp;
	//this->rightp = cp.rightp;
	this->sol = cp.sol;
	//this->soll = cp.soll;
	//this->solr = cp.solr;
	this->dsol = cp.dsol;
	//this->dsoll = cp.dsoll;
	//this->dsolr = cp.dsolr;
	this->HSize = cp.HSize;
	this->detjac = cp.detjac;
	//this->leftdetjac = cp.leftdetjac;
	//this->rightdetjac = cp.rightdetjac;
	this->intPtIndex = cp.intPtIndex;
    this->XCenter = cp.XCenter;
	//this->XLeftElCenter = cp.XLeftElCenter;
	//this->XRightElCenter = cp.XRightElCenter;
	this->fVecShapeIndex = cp.fVecShapeIndex;
	this->fNormalVec = cp.fNormalVec;
	this->numberdualfunctions = cp.numberdualfunctions;
	
	
	
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

/*
void TPZMaterialData::InvertLeftRightData(){
	TPZMaterialData cp(*this);
	this->leftdetjac = cp.rightdetjac;
	this->leftjac = cp.rightjac;
	this->leftjacinv = cp.rightjacinv;
	this->leftp = cp.rightp; 
	this->phil = cp.phir;
	this->dphixl = cp.dphixr;
	this->axesleft = cp.axesright;
	this->soll = cp.solr;
	this->dsoll = cp.dsolr;
	
	this->rightdetjac = cp.leftdetjac;
	this->rightjac = cp.leftjac;
	this->rightjacinv = cp.leftjacinv;
	this->rightp = cp.leftp;
	this->phir = cp.phil;
	this->dphixr = cp.dphixl;
	this->axesright = cp.axesleft;
	this->solr = cp.soll;
	this->dsolr = cp.dsoll;  
	
	this->XRightElCenter = cp.XLeftElCenter;
	this->XLeftElCenter = cp.XRightElCenter;
	
	const int n = this->normal.NElements();
	for(int i = 0; i < n; i++){
		this->normal[i] *= -1.;
	}
}
*/


/*
 * Save the element data to a stream
 */
void TPZMaterialData::Write(TPZStream &buf, int withclassid)
{
	phi.Write(buf,0);
	//phil.Write(buf,0);
	//phir.Write(buf,0);
	dphix.Write(buf,0);
	//dphixl.Write(buf,0);
	//dphixr.Write(buf,0);
	axes.Write(buf,0);
	//axesleft.Write(buf,0);
	//axesright.Write(buf,0);
	jacobian.Write(buf,0);
	//leftjac.Write(buf,0);
	//rightjac.Write(buf,0);
	jacinv.Write(buf,0);
	//leftjacinv.Write(buf,0);
	//rightjacinv.Write(buf,0);
	buf.Write(normal);
	buf.Write(x);
	buf.Write(&p,1);
	//buf.Write(&leftp,1);
	//buf.Write(&rightp,1);
    int nsol = sol.size();
    buf.Write(&nsol);
    for (int is=0; is<nsol; is++) {
        buf.Write(sol[is]);
    }
    /*
    nsol = soll.size();
    buf.Write(&nsol);
    for (int is=0; is<nsol; is++) {
        buf.Write(soll[is]);
    }
    nsol = solr.size();
    buf.Write(&nsol);
    for (int is=0; is<nsol; is++) {
        buf.Write(solr[is]);
    }
     */
    nsol = dsol.size();
    buf.Write(&nsol);
    for (int is=0; is<nsol; is++) {
        dsol[is].Write(buf,0);
    }
    /*
    nsol = dsoll.size();
    buf.Write(&nsol);
    for (int is=0; is<nsol; is++) {
        dsoll[is].Write(buf,0);
    }
    nsol = dsolr.size();
    buf.Write(&nsol);
    for (int is=0; is<nsol; is++) {
        dsolr[is].Write(buf,0);
    }
     */
	buf.Write(&HSize,1);
	buf.Write(&detjac,1);
	//buf.Write(&leftdetjac,1);
	//buf.Write(&rightdetjac,1);
    buf.Write(XCenter);
	//buf.Write(XLeftElCenter);
	//buf.Write(XRightElCenter);
	buf.Write(&intPtIndex,1);
}

/*
 * Read the element data from a stream
 */
void TPZMaterialData::Read(TPZStream &buf, void *context)
{
	phi.Read(buf,0);
	//phil.Read(buf,0);
	//phir.Read(buf,0);
	dphix.Read(buf,0);
	//dphixl.Read(buf,0);
	//dphixr.Read(buf,0);
	axes.Read(buf,0);
	//axesleft.Read(buf,0);
	//axesright.Read(buf,0);
	jacobian.Read(buf,0);
	//leftjac.Read(buf,0);
	//rightjac.Read(buf,0);
	jacinv.Read(buf,0);
	//leftjacinv.Read(buf,0);
	//rightjacinv.Read(buf,0);
	buf.Read(normal);
	buf.Read(x);
	buf.Read(&p,1);
	//buf.Read(&leftp,1);
	//buf.Read(&rightp,1);
    int nsol;
    buf.Read(&nsol,1);
    for (int is=0; is<nsol; is++) {
        buf.Read(sol[is]);
    }
    /*
    buf.Read(&nsol,1);
    for (int is=0; is<nsol; is++) {
        buf.Read(soll[is]);
    }
    buf.Read(&nsol,1);
    for (int is=0; is<nsol; is++) {
        buf.Read(solr[is]);
    }
     */
    buf.Read(&nsol,1);
    for (int is = 0; is<nsol; is++) {
        dsol[is].Read(buf,0);
    }
    /*
    for (int is = 0; is<nsol; is++) {
        dsoll[is].Read(buf,0);
    }
    for (int is = 0; is<nsol; is++) {
        dsolr[is].Read(buf,0);
    }
     */
	buf.Read(&HSize,1);
	buf.Read(&detjac,1);
	//buf.Read(&leftdetjac,1);
	//buf.Read(&rightdetjac,1);
    buf.Read(XCenter);
	//buf.Read(XLeftElCenter);
	//buf.Read(XRightElCenter);
	buf.Read(&intPtIndex,1);
	
}

//Compare the object for identity with the object pointed to, eventually copy the object
/*
 * compare both objects bitwise for identity. Put an entry in the log file if different
 * overwrite the calling object if the override flag is true
 */
bool TPZMaterialData::Compare(TPZSaveable *copy, bool override)
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
    /*
	locres = phil.Compare(&comp->phil,override);
	if(!locres)
	{
		LOGPZ_DEBUG(loggerCheck,"phil different")
	}
	result = result && locres;
	locres = phir.Compare(&comp->phir,override);
	if(!locres)
	{
		LOGPZ_DEBUG(loggerCheck,"phir different")
	}
	result = result && locres;
     */
	locres = dphix.Compare(&comp->dphix,override);
	if(!locres)
	{
		LOGPZ_DEBUG(loggerCheck,"dphix different")
	}
	result = result && locres;
    /*
	locres = dphixl.Compare(&comp->dphixl,override);
	if(!locres)
	{
		LOGPZ_DEBUG(loggerCheck,"dphil different")
	}
	result = result && locres;
	locres = dphixr.Compare(&comp->dphixr,override);
	if(!locres)
	{
		LOGPZ_DEBUG(loggerCheck,"dphixr different")
	}
	result = result && locres;
     */
	locres = axes.Compare(&comp->axes,override);
	if(!locres)
	{
		LOGPZ_DEBUG(loggerCheck,"axes different")
	}
	result = result && locres;
    /*
	locres = axesleft.Compare(&comp->axesleft,override);
	if(!locres)
	{
		LOGPZ_DEBUG(loggerCheck,"axesleft different")
	}
	result = result && locres;
	locres = axesright.Compare(&comp->axesright,override);
	if(!locres)
	{
		LOGPZ_DEBUG(loggerCheck,"axesright different")
	}
	result = result && locres;
     */
	locres = jacobian.Compare(&comp->jacobian,override);
	if(!locres)
	{
		LOGPZ_DEBUG(loggerCheck,"jacobian different")
	}
	result = result && locres;
    /*
	locres = leftjac.Compare(&comp->leftjac,override);
	if(!locres)
	{
		LOGPZ_DEBUG(loggerCheck,"left jacobian different")
	}
	result = result && locres;
	locres = rightjac.Compare(&comp->rightjac,override);
	if(!locres)
	{
		LOGPZ_DEBUG(loggerCheck,"right jacobian different")
	}
	result = result && locres;
     */
	locres = jacinv.Compare(&comp->jacinv,override);
	if(!locres)
	{
		LOGPZ_DEBUG(loggerCheck,"jacinv different")
	}
	result = result && locres;
    /*
	locres = leftjacinv.Compare(&comp->leftjacinv,override);
	if(!locres)
	{
		LOGPZ_DEBUG(loggerCheck,"leftjacinv different")
	}
	result = result && locres;
	locres = rightjacinv.Compare(&comp->rightjacinv,override);
	if(!locres)
	{
		LOGPZ_DEBUG(loggerCheck,"rightjacinv different")
	}
	result = result && locres;
     */
	
	/*
	 buf.Read(normal);
	 buf.Read(x);
	 buf.Read(&p,1);
	 buf.Read(&leftp,1);
	 buf.Read(&rightp,1);
	 buf.Read(sol);
	 buf.Read(soll);
	 buf.Read(solr);
	 dsol.Read(buf,0);
	 dsoll.Read(buf,0);
	 dsolr.Read(buf,0);
	 buf.Read(&HSize,1);
	 buf.Read(&detjac,1);
	 buf.Read(&leftdetjac,1);
	 buf.Read(&rightdetjac,1);
	 buf.Read(XLeftElCenter);
	 buf.Read(XRightElCenter);
	 buf.Read(&intPtIndex,1);
	 */
	return true;
}

// Compare the object for identity with the object pointed to, eventually copy the object
/*
 * compare both objects bitwise for identity. Put an entry in the log file if different
 * overwrite the calling object if the override flag is true
 */
bool TPZMaterialData::Compare(TPZSaveable *copy, bool override) const
{
	return true;
}

template class TPZRestoreClass<TPZMaterialData,TPZMATERIALDATAID>;

/*
 * Print the data
 */
void TPZMaterialData::Print(std::ostream &out) const
{
	phi.Print("phi",out);
	//phil.Print("phil",out);
	//phir.Print("phir",out);
	dphix.Print("dphix",out);
	//dphixl.Print("dphixl",out);
	//dphixr.Print("dphixr",out);
	axes.Print("axes",out);
	//axesleft.Print("axesleft",out);
	//axesright.Print("axesright",out);
	jacobian.Print("jacobian",out);
	//leftjac.Print("leftjac",out);
	//rightjac.Print("rightjac",out);
	jacinv.Print("jacinv",out);
	//leftjacinv.Print("leftjacinv",out);
	//rightjacinv.Print("rightjacinv",out);
	out << "normal " << normal << std::endl;
	out << "x " << x << std::endl;
	out << "p " << p << std::endl;
	//out << "leftp " << leftp << std::endl;
	//out << "rightp " << rightp << std::endl;
	out << "sol " << sol << std::endl;
	//out << "soll " << soll << std::endl;
	//out << "solr " << solr << std::endl;
    int nsol = dsol.size();
    for (int is=0; is<nsol; is++) {
        dsol[is].Print("dsol",out);

    }
    /*
    nsol = dsoll.size();
    for (int is=0; is<nsol; is++) {
        dsoll[is].Print("dsoll",out);
        
    }
    nsol = dsolr.size();
    for (int is=0; is<nsol; is++) {
        dsolr[is].Print("dsolr",out);
        
    }
     */
	out << "HSize " << HSize << std::endl;
	out << "detjac " << detjac << std::endl;
	//out << "leftdetjac " << leftdetjac << std::endl;
	//out << "rightdetjac " << rightdetjac << std::endl;
    out << "XCenter " << XCenter << std::endl;
	//out << "XLeftElCenter " << XLeftElCenter << std::endl;
	//out << "XRightElCenter " << XRightElCenter << std::endl;
	out << "intPtIndex " << intPtIndex << std::endl;
}

/*
 * Print the data in a format suitable for Mathematica
 */
void TPZMaterialData::PrintMathematica(std::ostream &out) const
{
	phi.Print("phi = ",out,EMathematicaInput);
	//phil.Print("phil = ",out,EMathematicaInput);
	//phir.Print("phir = ",out,EMathematicaInput);
	dphix.Print("dphix = ",out,EMathematicaInput);
	//dphixl.Print("dphixl = ",out,EMathematicaInput);
	//dphixr.Print("dphixr = ",out,EMathematicaInput);
	axes.Print("axes = ",out,EMathematicaInput);
	//axesleft.Print("axesleft = ",out,EMathematicaInput);
	//axesright.Print("axesright = ",out,EMathematicaInput);
	jacobian.Print("jacobian = ",out,EMathematicaInput);
	//leftjac.Print("leftjac = ",out,EMathematicaInput);
	//rightjac.Print("rightjac = ",out,EMathematicaInput);
	jacinv.Print("jacinv = ",out,EMathematicaInput);
	//leftjacinv.Print("leftjacinv = ",out,EMathematicaInput);
	//rightjacinv.Print("rightjacinv = ",out,EMathematicaInput);
	out << "normal = {" << normal << "};" << std::endl;
	out << "x = {" << x << "};" << std::endl;
	out << "p = " << p << ";" << std::endl;
	//out << "leftp = " << leftp << ";" << std::endl;
	//out << "rightp = " << rightp << ";" << std::endl;
	out << "sol = { " << sol << "};" << std::endl;
	//out << "soll = { " << soll << "};" << std::endl;
	//out << "solr = { " << solr << "};" << std::endl;
    int nsol=dsol.size();
    for (int is=0; is<nsol; is++) {
        std::stringstream sout;
        sout << "dsol" << is << " = ";
        dsol[is].Print(sout.str().c_str(),out,EMathematicaInput);
    }
    /*
    nsol=dsoll.size();
    for (int is=0; is<nsol; is++) {
        std::stringstream sout;
        sout << "dsoll" << is << " = ";
        dsoll[is].Print(sout.str().c_str(),out,EMathematicaInput);
    }
    nsol=dsolr.size();
    for (int is=0; is<nsol; is++) {
        std::stringstream sout;
        sout << "dsolr" << is << " = ";
        dsolr[is].Print(sout.str().c_str(),out,EMathematicaInput);
    }
     */
	out << "HSize = " << HSize << ";" << std::endl;
	out << "detjac = " << detjac << ";" << std::endl;
	//out << "leftdetjac =  " << leftdetjac << ";" << std::endl;
	//out << "rightdetjac =  " << rightdetjac << ";" << std::endl;
    out << "XCenter = {" << XCenter << "};" << std::endl;
	//out << "XLeftElCenter = {" << XLeftElCenter << "};" << std::endl;
	//out << "XRightElCenter = { " << XRightElCenter << "};"  << std::endl;
	out << "intPtIndex = " << intPtIndex << ";" <<std::endl;
}
