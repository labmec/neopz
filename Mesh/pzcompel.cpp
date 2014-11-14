/**
 * @file
 * @brief Contains the implementation of the TPZCompEl methods.
 */

#include "pzcompel.h"
#include "pzgeoel.h"
#include "pzgeoelside.h"
#include "pzmaterial.h"
#include "pzcmesh.h"
#include "pzbndcond.h"
#include "pzelmat.h"
#include "pzconnect.h"
#include "pzblockdiag.h"

#include "pzshapelinear.h"

#include "pzerror.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzblock.h"
#include "pzsolve.h"

#include "pzmanvector.h"
#include "pzstack.h"

#include "pztransfer.h"
#include "pztrnsform.h"
#include "pzquad.h"

#include "TPZInterfaceEl.h"
#include "TPZCompElDisc.h"
#include "pzintel.h"

#include <math.h>
#include <stdlib.h>

#include <sstream>
using namespace std;

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.tpzcompel"));
static LoggerPtr loggerSide(Logger::getLogger("pz.mesh.tpzcompelside"));
#endif

void TPZCompEl::CalcBlockDiagonal(TPZStack<long> &connectlist, TPZBlockDiagonal<STATE> & blockdiag) {
	TPZElementMatrix ek(this->Mesh(), TPZElementMatrix::EK),ef(this->Mesh(), TPZElementMatrix::EF);
	int b;
	CalcStiff(ek,ef);
	if(HasDependency()) {
		ek.ApplyConstraints();
		ef.ApplyConstraints();
		int numblock = ek.fConstrConnect.NElements();
		TPZVec<int> blocksize(numblock);
		
		for(b=0; b<numblock; b++) blocksize[b] = ek.fConstrBlock.Size(b);
		
		blockdiag.Initialize(blocksize);
		connectlist = ek.fConstrConnect;
		
		for(b=0; b<numblock; b++) {
			int blsize = blocksize[b];
			long conind = ek.fConstrConnect[b];
            TPZConnect &con = Mesh()->ConnectVec()[conind];
			if(con.HasDependency() || con.IsCondensed()) continue;
			//TPZFMatrix<REAL> ekbl(blsize,blsize);
			TPZFMatrix<STATE> ekbl(blsize,blsize);
			int r,c;
			//TPZBlock<REAL> &mbl = ek.fConstrBlock;
			TPZBlock<STATE> &mbl = ek.fConstrBlock;
			
			for(r=0; r<blsize; r++) {
				for(c=0; c<blsize; c++) {
					ekbl(r,c) = (mbl)(b,b,r,c);
				}
			}
			blockdiag.AddBlock(b,ekbl);
		}
	} else {
		int numblock = ek.fConnect.NElements();
		TPZVec<int> blocksize(numblock);
		
		for(b=0; b<numblock; b++) blocksize[b] = ek.fBlock.Size(b);
		
		blockdiag.Initialize(blocksize);
		connectlist = ek.fConnect;
		
		for(b=0; b<numblock; b++) {
			int blsize = blocksize[b];
			//TPZFMatrix<REAL> ekbl(blsize,blsize);
			TPZFMatrix<STATE> ekbl(blsize,blsize);
			long conind = ek.fConnect[b];
            TPZConnect &con = Mesh()->ConnectVec()[conind];
			if(con.HasDependency() || con.IsCondensed()) continue;
			int r, c;
			//TPZBlock<REAL> &mbl = ek.fBlock;
			TPZBlock<STATE> &mbl = ek.fBlock;

			for(r=0; r<blsize; r++) {
				for(c=0; c<blsize; c++) {
					//ekbl(r,c) = (mbl)(b,b,r,c);
					ekbl(r,c) = (mbl)(b,b,r,c);
				}
			}
			blockdiag.AddBlock(b,ekbl);
		}
	}
}

int TPZCompEl::gOrder = 2;

TPZCompEl::TPZCompEl() : fMesh(0), fIndex(-1), fReferenceIndex(-1) {
}

TPZCompEl::TPZCompEl(TPZCompMesh &mesh, TPZGeoEl *ref, long &index) {
	fMesh = &mesh;
	index = mesh.ElementVec().AllocateNewElement();
	mesh.ElementVec()[index] = this;
	fIndex = index;
	fReferenceIndex = (ref == 0) ? -1 : ref->Index();
}

TPZCompEl::TPZCompEl(TPZCompMesh &mesh, const TPZCompEl &copy) {
	fMesh = &mesh;
	long index = copy.fIndex;
	if(index >= 0) mesh.ElementVec()[index] = this;
	fIndex = index;
	fReferenceIndex = copy.fReferenceIndex;
}

TPZCompEl::TPZCompEl(TPZCompMesh &mesh, const TPZCompEl &copy, long &index) {
	fMesh = &mesh;
	index = mesh.ElementVec().AllocateNewElement();
	if(index >= 0) mesh.ElementVec()[index] = this;
	fIndex = index;
	fReferenceIndex = copy.fReferenceIndex;
}

TPZCompEl::TPZCompEl(TPZCompMesh &mesh, const TPZCompEl &copy, std::map<long,long> &gl2lcElMap)
{
	fMesh = &mesh;
	if (gl2lcElMap.find(copy.fIndex) == gl2lcElMap.end())
	{
		std::stringstream sout;
		sout << "ERROR in - " << __PRETTY_FUNCTION__
        << " original element index: " << copy.fIndex << " is not mapped!\n"
        << "Map content: ";
		std::map<long,long>::iterator it;
		for (it=gl2lcElMap.begin();it!=gl2lcElMap.end();it++) sout << " ( " << it->first << " | " << it->second << " ) ;";
		LOGPZ_ERROR (logger,sout.str().c_str());
		DebugStop();
	}
	long index = gl2lcElMap[copy.fIndex];
	if(index >= 0) mesh.ElementVec()[index] = this;
	fIndex = index;
	fReferenceIndex = copy.fReferenceIndex;
}

TPZCompEl::~TPZCompEl() {
	long index = Index();
    if (fMesh->ElementVec()[index] == this) {
        fMesh->ElementVec()[index] = 0;
        fMesh->ElementVec().SetFree(index);        
    }
#ifdef DEBUG
	TPZGeoEl *gel = Reference();
	if (gel && gel->Reference()) {
		DebugStop();
	}
#endif
    fIndex = -1;
    fReferenceIndex = -1;
    fMesh = 0;
}

MElementType TPZCompEl::Type() {
	LOGPZ_WARN(logger, "Type unknown");
	return ENoType;
}

void TPZCompEl::LoadSolution() {
	// an element without mesh is a free element
	if(!Mesh() || !HasDependency()) return;
	TPZStack<long> connectlist;
	int totalconnects;
	BuildConnectList(connectlist);
	totalconnects = connectlist.NElements();
	TPZManVector<int> dependenceorder(totalconnects);
	TPZConnect::BuildDependencyOrder(connectlist,dependenceorder,*this->Mesh());
	TPZMaterial * mat = Material();
	if(!mat) {
		LOGPZ_WARN(logger, "Exiting LoadSolution because a null material was reached.");
		return;
	}
	//int numstate = mat->NStateVariables(); 
	//TPZBlock<REAL> &block = Mesh()->Block();
	TPZBlock<STATE> &block = Mesh()->Block();
	//TPZFMatrix<REAL> &MeshSol = Mesh()->Solution();
	TPZFMatrix<STATE> &MeshSol = Mesh()->Solution();
	int maxdep = 0;
	int in;
	long iv,jv,idf;
	STATE coef;
	for(in=0;in<totalconnects;in++)
		maxdep = (maxdep < dependenceorder[in]) ? dependenceorder[in] : maxdep;
	int current_order = maxdep-1;
	while(current_order >= 0) {
		for(in=0; in<totalconnects; in++) {
			if(dependenceorder[in] != current_order) continue;
			TPZConnect *dfn = &Mesh()->ConnectVec()[connectlist[in]];
			if(!dfn->HasDependency()) continue;
			long bl = dfn->SequenceNumber();
			int nvar = block.Size(bl);
			int numstate = dfn->NState(); //numstate eh fornecida pelo connect 
			//         int numshape = nvar/numstate;
			TPZConnect::TPZDepend *dep = dfn->FirstDepend();
			long blpos = block.Position(bl);
			for(iv=0; iv<nvar; iv++) MeshSol(blpos+iv, 0) = 0.;
			while(dep) {
				long depconindex = dep->fDepConnectIndex;
				TPZConnect &depcon = Mesh()->ConnectVec()[depconindex];
				long depseq = depcon.SequenceNumber();
				int numdepvar = block.Size(depseq);
				long depseqpos = block.Position(depseq);
				for(iv=0; iv<nvar; iv+=numstate) {
					for(jv=0; jv<numdepvar; jv+=numstate) {
						coef = dep->fDepMatrix(iv/numstate,jv/numstate);
						for(idf=0; idf<numstate; idf++) MeshSol(blpos+iv+idf,0) += coef*MeshSol(depseqpos+jv+idf,0);
					}
				}
				dep = dep->fNext;
			}
		}
		current_order--;
	}
}

void TPZCompEl::SetMesh(TPZCompMesh *mesh) {
	// The mesh can be null, indicating that the element is now unused
	if(!mesh)
		LOGPZ_ERROR(logger, "TPZCompEl.SetMesh called with zero pointer.");
	fMesh = mesh;
}

TPZCompMesh *TPZCompEl::Mesh() const {
	if(!fMesh)
	{
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << " mesh address " << (void *) fMesh << " this address " << (void *) this;
		LOGPZ_WARN(logger, sout.str());
	}
	return fMesh;
}

TPZConnect &TPZCompEl::Connect(int i) const{
#ifndef NODEBUG
	if(fMesh) {
		long connectindex = ConnectIndex(i);
		if(connectindex >= 0) {
			return fMesh->ConnectVec()[connectindex];
		} else {
			LOGPZ_ERROR(logger, "TPZCompEl::Connect called for noninitialized connect\n");
			Reference()->Print(cout);
			const TPZCompEl *cel = (const TPZCompEl *) this;
			cel->Print(cout);
			DebugStop();
		}
	} else {
		LOGPZ_FATAL(logger, "Connect called for an element without mesh\n");
		DebugStop();
	}
	static TPZConnect dummy;
	return dummy;
#else
	return fMesh->ConnectVec()[ConnectIndex(i)];
#endif
}

void TPZCompEl::Print(std::ostream & out) const {
	out << "\nOutput for a computable element index: " << fIndex;
	if(this->Reference())
	{
		out << "\nCenter coordinate: ";
		TPZVec< REAL > centerMaster( this->Reference()->Dimension(),0. );
		TPZVec< REAL > centerEuclid( 3,0.);
		this->Reference()->CenterPoint(this->Reference()->NSides()-1,centerMaster);
		this->Reference()->X(centerMaster,centerEuclid);
		out << centerEuclid << std::endl;
	}
	if(this->Material())
	{
		out << "Material id " << this->Material()->Id() << "\n";
	}
	else {
		out << "No material\n";
	}
	
	out << "Number of connects = " << NConnects();
    out<< "\nConnect indexes : ";
	int nod;
	for(nod=0; nod< NConnects(); nod++)
	{
		out << ConnectIndex(nod) <<  ' ' ;
	}
	out << endl;
}

std::ostream& operator<<(std::ostream &s,TPZCompEl & el){
	el.Print(s);
	return s;
}

void TPZCompEl::PrintSolution(TPZVec<REAL> &point,char *varname,std::ostream &s) {
	TPZGeoEl *georef = Reference();
	TPZMaterial * mat = Material();
	if(!georef || !mat) {
		LOGPZ_WARN(logger, "Exiting PrintSolution should not be called for an element which doesnt have a geometric reference or material");
		Print();
		return;
	}
	int varindex = mat->VariableIndex(varname);
	int numvar = mat->NSolutionVariables(varindex);
	if(varindex == -1) {
		LOGPZ_WARN(logger, "Exiting PrintSolution should not be called for an element which has unknown variable index");
		return;
	}
	TPZManVector<STATE> sol(numvar);
	sol.Fill(0.);
	Solution(point,varindex,sol);
	for(long i=0; i<sol.NElements(); i++) {
		s << sol[i] << '\t';
	}
}

void TPZCompEl::PrintCoordinate(TPZVec<REAL> &point,int CoordinateIndex,std::ostream &s) {
	TPZGeoEl *georef = Reference();
	if(!georef) {
		LOGPZ_WARN(logger,"Exiting PrintCoordinate should not be called on an element which doesnt have a geometric reference");
		return;
	}
	TPZManVector<REAL> X(3);
	X.Fill(0.);
	georef->X(point,X);
	s << X[CoordinateIndex] << '\t';
}

void TPZCompEl::PrintTitle(char *varname,std::ostream &s) {
	TPZMaterial * mat = Material();
	TPZGeoEl *georef = Reference();
	if(!georef || !mat) {
		LOGPZ_WARN(logger,"Exiting PrintTitle should not be called for an element which doesnt have a material");
		return;
	}
	int varindex = mat->VariableIndex(varname);
	if(varindex == -1) return;
	int numvar = mat->NSolutionVariables(varindex);
	if(numvar == 0) return;
	if(numvar == 1) {
		s << varname << '\t';
		return;
	}
	for(int i=0; i<numvar; i++) s << varname << '_' << i << '\t';
}

inline void TPZCompEl::Divide(long index, TPZVec<long> &subindex, int interpolate) {
	subindex.Resize(0);
	LOGPZ_WARN(logger,"TPZCompEl::Divide called");
}

void TPZCompEl::EvaluateError(void (* /*fp*/)(const TPZVec<REAL> &loc,TPZVec<STATE> &val,TPZFMatrix<STATE> &deriv),
                              TPZVec<REAL> &/*errors*/,TPZBlock<REAL> * /*flux*/) {
	LOGPZ_WARN(logger, "EvaluateError is called.");
}

void TPZCompEl::Solution(TPZVec<REAL> &/*qsi*/,int var,TPZVec<STATE> &sol){
	if(var >= 100) {
		int ind = Index();
		if(fMesh->ElementSolution().Cols() > var-100) {
			sol[0] = fMesh->ElementSolution()(ind,var-100);
		} else {
			sol[0] = 0.;
		}
	} else {
		sol.Resize(0);
	}
}

void TPZCompEl::BuildConnectList(std::set<long> &indepconnectlist,
								 std::set<long> &depconnectlist) {
	const int ncon = this->NConnects();
	for(int i = 0; i < ncon; i++) {
		long conind = ConnectIndex(i);
		Connect(i).BuildConnectList(conind, indepconnectlist,depconnectlist,*Mesh());
	}
}

void TPZCompEl::BuildConnectList(TPZStack<long> &connectlist) {
	std::set<long> buf,buf2;
	long ncon = connectlist.NElements();
	for(long i = 0; i < ncon; i++) {
		buf.insert(connectlist[i]);
	}
	BuildConnectList(buf2);
	ncon = buf2.size();
	std::set<long>::iterator it = buf2.begin();
	for(it=buf2.begin() ; it != buf2.end(); it++) 
	{
        if(buf.find(*it) == buf.end())
        {
            connectlist.Push(*it);
            buf.insert(*it);
        }
	}
}

void TPZCompEl::BuildConnectList(std::set<long> &connectlist) {
	std::set<long> additional;
	const int ncon = this->NConnects();
	for(int i = 0; i < ncon; i++) {
		additional.insert(this->ConnectIndex(i));
	}
	TPZConnect::BuildConnectList(connectlist, additional, *this->Mesh());
}

int TPZCompEl::HasDependency() {
	int nconnects = NConnects();
	int in;
	for(in=0; in<nconnects; in++) if(Connect(in).HasDependency()){
		return 1;
	}
	return 0;
}

void TPZCompEl::SetIndex(long index) {
	fIndex = index;
	return;
}

int TPZCompEl::NEquations() {
	int numeq=0;
	for (int i=0;i<NConnects(); i++){
		TPZConnect &df = Connect(i);
		if(df.HasDependency() || df.IsCondensed() || !df.NElConnected() || df.SequenceNumber() == -1) continue;
		long seqnum = df.SequenceNumber();
		numeq += Mesh()->Block().Size(seqnum);
	}
	return numeq;
}

REAL TPZCompEl::CompareElement(int var, char *matname){
	cout << "TPZCompEl::CompareElement called!\n";
	LOGPZ_WARN(logger, "CompareElement called!");
	return 0.;
}

// TPZCompElSide //

void TPZCompEl::LoadElementReference()
{
	TPZGeoEl *ref = Reference();
	if(ref)
    {
        ref->SetReference(this);
    }
}

void TPZCompEl::CalcResidual(TPZElementMatrix &ef){
	TPZElementMatrix ek(this->Mesh(), TPZElementMatrix::EK);
	CalcStiff(ek,ef);
}

TPZGeoEl * TPZCompEl::GetRefElPatch(){
	TPZGeoEl *ref = Reference();
	if (!ref) {
		LOGPZ_ERROR(logger, "reached a null reference");
		return (0);
	}
#ifdef LOG4CXX
	std::stringstream sout;
    if (logger->isDebugEnabled())
    {
        sout << "Obtendo elemento geometrico de referencia para elemento " << Index() << endl;
        sout << "Impressao dos ancestrais\n";
        Print(sout);
    }
    
    if (logger->isDebugEnabled())
    {
        ref->Print(sout);
    }
#endif
    
	TPZStack <TPZGeoEl *> ancestors;
	ancestors.Push(ref);
	while(ref->Father()){
		ancestors.Push(ref->Father());
		ref = ref->Father();
#ifdef LOG4CXX
		ref->Print(sout);
#endif
	}
	int j;
	LOGPZ_DEBUG(logger, sout.str());
	while(ancestors.NElements()) {
		TPZGeoEl *larger = ancestors.Pop();
		for (j=0; j<larger->NSides(); j++){
			int sidedimension = larger->SideDimension(j);
			if(sidedimension == 0){
				continue;
			}
			TPZGeoElSide gelside(larger,j);
			TPZStack <TPZCompElSide> stack;
			gelside.EqualLevelCompElementList(stack,1,1);
            // the first geometric element that has a neighbour is a reference element
			if(stack.NElements()){
				//cout << " \n \n \n ==================================\n ================================\nElemento PatchReference\n";
				//aux->Print();
				return larger;
			}
		}
	}
	//cout << " \n \n \n ==================================\n ================================\nElemento PatchReference falho\n";
	//Reference()->Print();
	LOGPZ_DEBUG(logger, "Exit GetRefElPatch - Element is its own patch");
	return (Reference());
}

REAL TPZCompEl::MaximumRadiusOfEl(){
	//O elemento deve ser a envoltura convexa dos seus v�tices
	if(!this) {
		LOGPZ_ERROR(logger,"TPZCompMesh::MaximumRadiusOfEl() null element");
	}
	
	int i,j,k;
	REAL maxdist = 0.0,dist=0.0;
	TPZVec<REAL> point0(3),point1(3);
	TPZGeoNode *node0,*node1;
	int nvertices = Reference()->NNodes();
	for(i=0;i<nvertices;i++){
		for(j=i+1;j<nvertices;j++){
			node0 = Reference()->NodePtr(i);
			node1 = Reference()->NodePtr(j);
			for(k=0;k<3;k++){
				point0[k] = node0->Coord(k);
				point1[k] = node1->Coord(k);
			}
			dist = TPZGeoEl::Distance(point0,point1);
			if(dist > maxdist) maxdist = dist;
		}
	}
	return maxdist;
}

REAL TPZCompEl::LesserEdgeOfEl(){
	
	if(!this) LOGPZ_INFO(logger,"TPZCompMesh::LesserEdgeOfEl null element");
	
	int i,j,k;
	REAL mindist = 1000.0,dist=0.0;
	TPZVec<REAL> point0(3),point1(3);
	TPZGeoNode *node0,*node1;
	int nvertices = Reference()->NNodes();
	for(i=0;i<nvertices;i++){
		for(j=i+1;j<nvertices;j++){
			node0 = Reference()->NodePtr(i);
			node1 = Reference()->NodePtr(j);
			for(k=0;k<3;k++){
				point0[k] = node0->Coord(k);
				point1[k] = node1->Coord(k);
			}
			dist = TPZGeoEl::Distance(point0,point1);
			if(dist < mindist) mindist = dist;
		}
	}
	
	return mindist;
}


/**
 Save the element data to a stream
 */
void TPZCompEl::Write(TPZStream &buf, int withclassid)
{
	TPZSaveable::Write(buf,withclassid);
	buf.Write(&fIndex,1);
	buf.Write(&fReferenceIndex,1);
}

/**
 Read the element data from a stream
 */
void TPZCompEl::Read(TPZStream &buf, void *context)
{
	TPZSaveable::Read(buf,context);
	fMesh = (TPZCompMesh *) context;
	buf.Read(&fIndex,1);
	buf.Read(&fReferenceIndex,1);
}

void TPZCompEl::SetOrthogonalFunction(void (*orthogonal)(REAL x,int num,
														 TPZFMatrix<REAL> & phi,TPZFMatrix<REAL> & dphi)) {
	pzshape::TPZShapeLinear::fOrthogonal = orthogonal;
}

TPZCompElSide::TPZCompElSide() {
	fEl = 0;
	fSide = -1;
}

TPZCompElSide::TPZCompElSide(const TPZCompElSide &celside) {
	fEl = celside.fEl;
	fSide   = celside.fSide;
}

TPZCompElSide::TPZCompElSide(TPZCompEl *cel,int side) {
	fEl = cel;
	fSide   = side;
}

TPZGeoElSide TPZCompElSide::Reference() const {
	if(!fEl) {
		LOGPZ_WARN(loggerSide, "Exiting Reference - non initialized side element reached");
		return TPZGeoElSide();
	}
	TPZGeoElSide sideel(fEl->Reference(),fSide);
	return sideel;
}

void TPZCompElSide::SetSide(int side) {
	fSide = side;
	Reference().SetSide(side);
}

void TPZCompElSide::HigherLevelElementList(TPZStack<TPZCompElSide> &elvec,
                                           int onlyinterpolated, int removeduplicates) {
	TPZGeoElSide georef = Reference();
	if(!georef.Element()) {
		LOGPZ_WARN(loggerSide, "Exiting HigherLevelElementList - null reference reached");
		return;
	}
	georef.HigherLevelCompElementList2(elvec,onlyinterpolated,removeduplicates);
}

void TPZCompElSide::EqualLevelElementList(TPZStack<TPZCompElSide> &elsidevec,
                                          int onlyinterpolated, int removeduplicates) {
	TPZGeoElSide georef = Reference();
	if(!georef.Exists()) {
		LOGPZ_WARN(loggerSide, "Exiting EqualLevelElementList - null reference reached");
		return;
	}
	georef.EqualLevelCompElementList(elsidevec,onlyinterpolated,removeduplicates);
}

void TPZCompElSide::HigherDimensionElementList(TPZStack<TPZCompElSide> &elsidevec,
                                               int onlyinterpolated, int removeduplicates) {
	TPZGeoElSide georef = Reference();
	if(!georef.Exists()) {
		LOGPZ_INFO(loggerSide, "Entering HigherDimensionElementList - null reference reached");
		return;
	}
	georef.HigherDimensionElementList(elsidevec,onlyinterpolated);
	if(removeduplicates) RemoveDuplicates(elsidevec);
}


void TPZCompElSide::RemoveDuplicates(TPZStack<TPZCompElSide> &elvec) {
	RemoveConnectDuplicates(elvec);
	int i, nelem = elvec.NElements();
	for(i=0; i<nelem; i++) {
		TPZGeoElSide geli = elvec[i].Reference();
		if(!geli.Exists()) continue;
		int j;
		for(j=i+1; j<nelem; j++) {
			TPZGeoElSide gelj = elvec[j].Reference();
			if(!gelj.Exists()) continue;
			if(geli.NeighbourExists(gelj)) {
				int k;
				LOGPZ_ERROR(loggerSide, "case not identified by RemoveConnectDuplicates");
				for(k=j;k<nelem-1;k++) elvec[k] = elvec[k+1];
				elvec.Resize(nelem-1);
				nelem--;
				j--;
			}
		}
	}
}

void TPZCompElSide::ConnectedElementList(TPZStack<TPZCompElSide> &ellist,
                                         int onlyinterpolated, int removeduplicates) {
	TPZCompElSide father = LowerLevelElementList(onlyinterpolated);
	if(father.Exists()) ellist.Push(father);
	EqualLevelElementList(ellist,onlyinterpolated,removeduplicates);
	HigherLevelElementList(ellist,onlyinterpolated,removeduplicates);
}

//retorna o elem. de menor nivel ao qual estou restrito
TPZCompElSide TPZCompElSide::LowerLevelElementList(int onlyinterpolated) {
	TPZGeoElSide georef = Reference();
	if(!georef.Exists()) {
		LOGPZ_WARN(loggerSide, "Exiting LowerLevelElementList - null reference reached");
		return TPZCompElSide();
	}
	return georef.LowerLevelCompElementList2(onlyinterpolated);
}

void TPZCompElSide::ExpandConnected(TPZStack<TPZCompElSide> &expandvec,int onlyinterpolated) {
	//insidevec contem todos os elem. de maior nivel que dependem de mim,
	//obtidos com HigherLevelElementList(..)
	TPZGeoElSide georef = Reference();
	if(!georef.Exists()) {
		LOGPZ_WARN(loggerSide, "Exiting ExpandConnected - null reference reached");
		return;
	}
	TPZStack<TPZCompElSide> highdimsidevec;
	TPZStack<int> smallsides;
	TPZCompElSide compside,compsidenext;
	long exnel = expandvec.NElements();
	for (long i=0;i<exnel;i++) {
		TPZGeoElSide ref = expandvec[i].Reference();
		smallsides.Resize(0);
		ref.Element()->LowerDimensionSides(ref.Side(),smallsides);
		for (int k=0;k<smallsides.NElements();k++) {
			compside = TPZCompElSide(Element(),smallsides[k]);
			compsidenext = compside.LowerLevelElementList(onlyinterpolated);
			//      TPZGeoElSide smallgeoside = compside.Reference();
			while(compsidenext.Exists()) {
				//         TPZGeoElSide geoside = compsidenext.Reference();//linha para teste
				if (compsidenext.Reference().NeighbourExists(this->Reference())) {
					TPZCompElSide lowid = LowerIdElementList(compside,onlyinterpolated);
					long l=0;
					while (l<expandvec.NElements() && lowid.Reference() != expandvec[l].Reference()) l++;
					if(l == expandvec.NElements()) {
						expandvec.Push(lowid);
					}
					break;
				}
				compside = compsidenext;
				compsidenext = compside.LowerLevelElementList(onlyinterpolated);
			}
		}
	}
}

TPZCompElSide TPZCompElSide::LowerIdElementList(TPZCompElSide &expandvec,int onlyinterpolated) {
	if(!expandvec.Element()) {
		LOGPZ_WARN(loggerSide, "Exiting LowerIdElementList - empty list");
		return TPZCompElSide();
	}
	TPZGeoElSide gelside = expandvec.Reference();
	TPZStack<TPZGeoElSide> neighbourset;
	gelside.AllNeighbours(neighbourset);
	neighbourset.Push(gelside);
	TPZGeoElSide lowidneigh(gelside);
	int lowid = gelside.Id();
	int in = 0, nneigh = neighbourset.NElements()-1;
	while(in < nneigh) {
		TPZCompEl *ref = neighbourset[in].Reference().Element();
		if(neighbourset[in].Id() < lowid && ref && (!onlyinterpolated || dynamic_cast<TPZInterpolatedElement*>(ref)    )) {
			lowidneigh = neighbourset[in];
			lowid = lowidneigh.Id();
		}
		in++;
	}
	return lowidneigh.Reference();
}

void TPZCompElSide::RemoveConnectDuplicates(TPZStack<TPZCompElSide> &expandvec){
	long i,k;
	long nelems = expandvec.NElements();
	TPZStack<TPZCompElSide> locexpand;
	for(i=0;i<nelems;i++) locexpand.Push(expandvec[i]);
	expandvec.Resize(0);
	for(k=0;k<nelems;k++){
		TPZCompEl *kel = locexpand[k].Element();
		if(!kel) continue;
		int kside = locexpand[k].Side();
		i=k+1;
		while(i<nelems){
			TPZCompEl *iel = locexpand[i].Element();
			int iside = locexpand[i].Side();
			if(iel && kel->ConnectIndex(kside) == iel->ConnectIndex(iside))
				locexpand[i] = TPZCompElSide();
			i++;
		}
	}
	for(i=0;i<nelems;i++)
		if(locexpand[i].Element()) expandvec.Push(locexpand[i]);
}

/// Return the index of the middle side connect alon fSide
long TPZCompElSide::ConnectIndex() const
{
    if(fEl)
    {
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(fEl);
        if(intel)
        {
            int locid = intel->MidSideConnectLocId(fSide);
            // verify whether there is a connect associated with this side
            if(locid >= 0)
            {
                return intel->ConnectIndex(intel->MidSideConnectLocId(fSide));
            }
            else
            {
                return -1;
            }
        }
        else
        {
            return fEl->ConnectIndex(fSide);
        }
    }
    else return -1;
}


/** Identify the material object associated with the element */
TPZMaterial * TPZCompEl::Material() const
{
    TPZMaterial * result = 0;
	TPZGeoEl *gel = Reference();
    if(fMesh && gel) result = fMesh->FindMaterial(gel->MaterialId());
    return result;
}

/** Verify if the material associated with the element is contained in the set */
bool TPZCompEl::HasMaterial(const std::set<int> &materialids)
{
	TPZMaterial * mat = Material();
	if(!mat) return false;
	return materialids.find(mat->Id()) != materialids.end();
}

/**
 * @brief Computes solution and its derivatives in the local coordinate qsi.
 * @param qsi master element coordinate
 * @param sol finite element solution
 * @param dsol solution derivatives
 * @param axes axes associated with the derivative of the solution
 */
void TPZCompEl::ComputeSolution(TPZVec<REAL> &qsi,
                             TPZSolVec &sol, TPZGradSolVec &dsol,TPZFMatrix<REAL> &axes)
{
    std::cout << __PRETTY_FUNCTION__ << " is not implemented - bailing out\n";
    DebugStop();
}

/**
 * @brief Computes solution and its derivatives in the local coordinate qsi. \n
 * This method will function for both volumetric and interface elements
 * @param qsi master element coordinate of the interface element
 * @param normal vector
 * @param leftsol finite element solution
 * @param dleftsol solution derivatives
 * @param leftaxes axes associated with the left solution
 * @param rightsol finite element solution
 * @param drightsol solution derivatives
 * @param rightaxes axes associated with the right solution
 */
void TPZCompEl::ComputeSolution(TPZVec<REAL> &qsi,
                             TPZVec<REAL> &normal,
                             TPZSolVec &leftsol, TPZGradSolVec &dleftsol,TPZFMatrix<REAL> &leftaxes,
                             TPZSolVec &rightsol, TPZGradSolVec &drightsol,TPZFMatrix<REAL> &rightaxes)
{
    std::cout << __PRETTY_FUNCTION__ << " is not implemented - bailing out\n";
    DebugStop(); 
}

/**
 * @brief Computes solution and its derivatives in local coordinate qsi
 * @param qsi master element coordinate
 * @param phi matrix containing shape functions compute in qsi point
 * @param dphix matrix containing the derivatives of shape functions in the direction of the axes
 * @param axes [in] axes indicating the direction of the derivatives
 * @param sol finite element solution
 * @param dsol solution derivatives
 */
void TPZCompEl::ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphix,
                             const TPZFMatrix<REAL> &axes, TPZSolVec &sol, TPZGradSolVec &dsol)
{
    std::cout << __PRETTY_FUNCTION__ << " is not implemented - bailing out\n";
    DebugStop();
}

/**
 * @brief Returns the index of the pressure connect
 * @note Returns -1 if their is no pressure connect
 */
int TPZCompEl::PressureConnectIndex() const
{
    int ncon = NConnects();
    int index = -1;
    int count = 0;
    for (int ic=0; ic<ncon ; ic++) {
        long locconnectindex = ConnectIndex(ic);
        TPZConnect &c = fMesh->ConnectVec()[locconnectindex];
        if (c.LagrangeMultiplier() && ! c.IsCondensed()) {
            index = ic;
            count++;
        }
    }
    if (count > 1) {
        DebugStop();
    }
    return index;
}

/**
 * @brief Compute the integral of a variable
 */
TPZVec<STATE> TPZCompEl::IntegrateSolution(int var) const
{
    TPZManVector<STATE,3> result(0);
    TPZGeoEl *gel = Reference();
    if (!gel) {
        return result;
    }
    TPZCompEl *celnotconst = (TPZCompEl *) this;
    int matid = gel->MaterialId();
    TPZMaterial *material = Material();
    int nvar  = material->NSolutionVariables(var);
    TPZIntPoints *intrule = gel->CreateSideIntegrationRule(gel->NSides()-1, 5);
    int dim = gel->Dimension();
    TPZManVector<REAL,3> xi(dim);
    TPZMaterialData data;
    TPZFNMatrix<9,REAL> jac(dim,dim),jacinv(dim,dim),axes(dim,3);
    REAL detjac;
    TPZManVector<STATE,3> sol(nvar);
    result.Resize(nvar, 0.);
    int npoints = intrule->NPoints();
    for (int ip =0; ip<npoints; ip++) {
        REAL weight;
        intrule->Point(ip, xi, weight);
        gel->Jacobian(xi, jac, axes, detjac, jacinv);
        celnotconst->Solution(xi, var, sol);
        for (int i=0; i <nvar; i++) {
            result[i] += weight*fabs(detjac)*sol[i];
        }
    }
    delete intrule;
    return result;
}

