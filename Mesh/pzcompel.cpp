/**
 * @file
 * @brief Contains the implementation of the TPZCompEl methods.
 */
//$Id: pzcompel.cpp,v 1.52 2011-05-11 02:27:20 phil Exp $

//METHODS DEFINITION FOR CLASS ELBAS

#include "pzcompel.h"
#include "pzgeoel.h"
#include "pzmaterial.h"
#include "pzcmesh.h"
#include "pzbndcond.h"
#include "pzelmat.h"
#include "pzconnect.h"
#include "pzblockdiag.h"

#include "pzshapelinear.h"
//#include "pzshapequad.h"
//#include "pzshapetriang.h"
//#include "pzshapetetra.h"
//#include "pzshapepiram.h"
//#include "pzshapeprism.h"
//#include "pzshapecube.h"

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

// LOGPZ
//#include <LOGPZ/logger.h>
// #include <LOGPZ/basicconfigurator.h>
// #include <LOGPZ/propertyconfigurator.h>
// #include <LOGPZ/helpers/exception.h>
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.tpzcompel"));
static LoggerPtr loggerSide(Logger::getLogger("pz.mesh.tpzcompelside"));
#endif

void TPZCompEl::CalcBlockDiagonal(TPZStack<int> &connectlist, TPZBlockDiagonal<REAL> & blockdiag) {
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
			int conind = ek.fConstrConnect[b];
            TPZConnect &con = Mesh()->ConnectVec()[conind];
			if(con.HasDependency() || con.IsCondensed()) continue;
			TPZFMatrix<REAL> ekbl(blsize,blsize);
			int r,c;
			TPZBlock<REAL> &mbl = ek.fConstrBlock;
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
			TPZFMatrix<REAL> ekbl(blsize,blsize);
			int conind = ek.fConnect[b];
            TPZConnect &con = Mesh()->ConnectVec()[conind];
			if(con.HasDependency() || con.IsCondensed()) continue;
			int r,c;
			TPZBlock<REAL> &mbl = ek.fBlock;
			for(r=0; r<blsize; r++) {
				for(c=0; c<blsize; c++) {
					ekbl(r,c) = (mbl)(b,b,r,c);
				}
			}
			blockdiag.AddBlock(b,ekbl);
		}
	}
}

int TPZCompEl::gOrder = 2;

/*
 void (*TPZCompEl::fOrthogonal)(REAL x,int num,TPZFMatrix<REAL> & phi,TPZFMatrix<REAL> & dphi) = TPZCompEl::Chebyshev;
 void (*TPZShapeLinear::fOrthogonal)(REAL x,int num,TPZFMatrix<REAL> & phi,TPZFMatrix<REAL> & dphi) = TPZCompEl::fOrthogonal;
 void (*TPZShapeQuad::fOrthogonal)(REAL x,int num,TPZFMatrix<REAL> & phi,TPZFMatrix<REAL> & dphi) = TPZCompEl::fOrthogonal;
 void (*TPZShapeTriang::fOrthogonal)(REAL x,int num,TPZFMatrix<REAL> & phi,TPZFMatrix<REAL> & dphi) = TPZCompEl::fOrthogonal;
 void (*TPZShapeTetra::fOrthogonal)(REAL x,int num,TPZFMatrix<REAL> & phi,TPZFMatrix<REAL> & dphi) = TPZCompEl::fOrthogonal;
 void (*TPZShapePiram::fOrthogonal)(REAL x,int num,TPZFMatrix<REAL> & phi,TPZFMatrix<REAL> & dphi) = TPZCompEl::fOrthogonal;
 void (*TPZShapePrism::fOrthogonal)(REAL x,int num,TPZFMatrix<REAL> & phi,TPZFMatrix<REAL> & dphi) = TPZCompEl::fOrthogonal;
 void (*TPZShapeCube::fOrthogonal)(REAL x,int num,TPZFMatrix<REAL> & phi,TPZFMatrix<REAL> & dphi) = TPZCompEl::fOrthogonal;
 */

TPZCompEl::TPZCompEl() {
#ifdef LOG4CXX
	LOGPZ_DEBUG(logger,__PRETTY_FUNCTION__);
#endif
	fMesh = 0;
	fIndex = -1;
}

TPZCompEl::TPZCompEl(TPZCompMesh &mesh, TPZGeoEl *ref, int &index) {
	/*#ifdef LOG4CXX
	 {
	 std::stringstream sout;
	 sout << __PRETTY_FUNCTION__ << " Mesh add" << (void*) &mesh << " this address " << (void*) this;          LOGPZ_DEBUG(logger,sout.str());
	 }
	 #endif*/
	fMesh = &mesh;
	index = mesh.ElementVec().AllocateNewElement();
	mesh.ElementVec()[index] = this;
	fIndex = index;
	//  fReference = ref;
	fReferenceIndex = (ref == 0) ? -1 : ref->Index();
}

TPZCompEl::TPZCompEl(TPZCompMesh &mesh, const TPZCompEl &copy) {
	/*#ifdef LOG4CXX
	 LOGPZ_DEBUG(logger,__PRETTY_FUNCTION__);
	 #endif*/
	fMesh = &mesh;
	int index = copy.fIndex;
	if(index >= 0) mesh.ElementVec()[index] = this;
	fIndex = index;
	//  fReference = copy.Reference();
	fReferenceIndex = copy.fReferenceIndex;
}

TPZCompEl::TPZCompEl(TPZCompMesh &mesh, const TPZCompEl &copy, int &index) {
	/*#ifdef LOG4CXX
	 LOGPZ_DEBUG(logger,__PRETTY_FUNCTION__);
	 #endif*/
	fMesh = &mesh;
	index = mesh.ElementVec().AllocateNewElement();
	if(index >= 0) mesh.ElementVec()[index] = this;
	fIndex = index;
	//  fReference = copy.fReference;
	fReferenceIndex = copy.fReferenceIndex;
}

TPZCompEl::TPZCompEl(TPZCompMesh &mesh, const TPZCompEl &copy, std::map<int,int> &gl2lcElMap)
{
	/*#ifdef LOG4CXX
	 LOGPZ_DEBUG(logger,__PRETTY_FUNCTION__);
	 #endif*/
	fMesh = &mesh;
	if (gl2lcElMap.find(copy.fIndex) == gl2lcElMap.end())
	{
		std::stringstream sout;
		sout << "ERROR in - " << __PRETTY_FUNCTION__
        << " original element index: " << copy.fIndex << " is not mapped!\n"
        << "Map content: ";
		std::map<int,int>::iterator it;
		for (it=gl2lcElMap.begin();it!=gl2lcElMap.end();it++) sout << " ( " << it->first << " | " << it->second << " ) ;";
		LOGPZ_ERROR (logger,sout.str().c_str());
		DebugStop();
	}
	int index = gl2lcElMap[copy.fIndex];
	if(index >= 0) mesh.ElementVec()[index] = this;
	fIndex = index;
	//  fReference = copy.Reference();
	fReferenceIndex = copy.fReferenceIndex;
}

TPZCompEl::~TPZCompEl() {
	int index = Index();
    if (fMesh->ElementVec()[index] == this) {
        fMesh->ElementVec()[index] = 0;
        fMesh->ElementVec().SetFree(index);        
    }
}

MElementType TPZCompEl::Type() {
	LOGPZ_WARN(logger, "Type unknown");
	return ENoType;
}

void TPZCompEl::LoadSolution() {
	// an element without mesh is a free element
	if(!Mesh() || !HasDependency()) return;
	TPZStack<int> connectlist;
	int totalconnects;
	BuildConnectList(connectlist);
	totalconnects = connectlist.NElements();
	TPZManVector<int> dependenceorder(totalconnects);
	TPZConnect::BuildDependencyOrder(connectlist,dependenceorder,*this->Mesh());
	TPZAutoPointer<TPZMaterial> mat = Material();
	if(!mat) {
		LOGPZ_WARN(logger, "Exiting LoadSolution because a null material was reached.");
		return;
	}
	//int numstate = mat->NStateVariables(); 
	TPZBlock<REAL> &block = Mesh()->Block();
	TPZFMatrix<REAL> &MeshSol = Mesh()->Solution();
	int maxdep = 0;
	int in;
	int iv,jv,idf;
	REAL coef;
	for(in=0;in<totalconnects;in++)
		maxdep = (maxdep < dependenceorder[in]) ? dependenceorder[in] : maxdep;
	int current_order = maxdep-1;
	while(current_order >= 0) {
		for(in=0; in<totalconnects; in++) {
			if(dependenceorder[in] != current_order) continue;
			TPZConnect *dfn = &Mesh()->ConnectVec()[connectlist[in]];
			if(!dfn->HasDependency()) continue;
			int bl = dfn->SequenceNumber();
			int nvar = block.Size(bl);
			int numstate = dfn->NState(); //numstate eh fornecida pelo connect 
			//         int numshape = nvar/numstate;
			TPZConnect::TPZDepend *dep = dfn->FirstDepend();
			int blpos = block.Position(bl);
			for(iv=0; iv<nvar; iv++) MeshSol(blpos+iv, 0) = 0.;
			while(dep) {
				int depconindex = dep->fDepConnectIndex;
				TPZConnect &depcon = Mesh()->ConnectVec()[depconindex];
				int depseq = depcon.SequenceNumber();
				int numdepvar = block.Size(depseq);
				int depseqpos = block.Position(depseq);
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
		int connectindex = ConnectIndex(i);
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
	out << "Output for a computable element index: " << fIndex << std::endl;
	if(this->Reference())
	{
		out << "Center coordinate: ";
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
	
	out << "Number of connects = " << NConnects() << " Node indexes : ";
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
	TPZAutoPointer<TPZMaterial> mat = Material();
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
	TPZManVector<REAL> sol(numvar);
	sol.Fill(0.);
	Solution(point,varindex,sol);
	for(int i=0; i<sol.NElements(); i++) {
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
	TPZAutoPointer<TPZMaterial> mat = Material();
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

// void TPZCompEl::Chebyshev(REAL x,int num,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi){
//   // Quadratic or higher shape functions
//   if(num <= 0) {
//     LOGPZ_INFO(logger, "Exiting Chebyshev - num <= 0");
//     return;
//   }
//   phi.Put(0,0,1.0);
//   dphi.Put(0,0, 0.0);
//   if(num == 1) {
//     LOGPZ_INFO(logger, "Exiting Chebyshev num == 1");
//     return;
//   }
//   phi.Put(1,0, x);
//   dphi.Put(0,1, 1.0);
//   int ord;
//   for(ord = 2;ord<num;ord++) {
//     phi.Put(ord,0, 2.0*x*phi(ord-1,0) - phi(ord-2,0));
//     dphi.Put(0,ord, 2.0*x*dphi(0,ord-1) + 2.0*phi(ord-1,0) - dphi(0,ord-2));
//   }
//   LOGPZ_INFO(logger, "Exiting Chebyshev");
// }

inline void TPZCompEl::Divide(int index, TPZVec<int> &subindex, int interpolate) {
	subindex.Resize(0);
	LOGPZ_WARN(logger,"TPZCompEl::Divide called");
}

void TPZCompEl::EvaluateError(void (* /*fp*/)(TPZVec<REAL> &loc,TPZVec<REAL> &val,TPZFMatrix<REAL> &deriv),
                              TPZVec<REAL> &/*errors*/,TPZBlock<REAL> * /*flux*/) {
	LOGPZ_WARN(logger, "EvaluateError is called.");
}

void TPZCompEl::Solution(TPZVec<REAL> &/*qsi*/,int var,TPZVec<REAL> &sol){
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

void TPZCompEl::BuildConnectList(std::set<int> &indepconnectlist,
								 std::set<int> &depconnectlist) {
	const int ncon = this->NConnects();
	for(int i = 0; i < ncon; i++) {
		int conind = ConnectIndex(i);
		Connect(i).BuildConnectList(conind, indepconnectlist,depconnectlist,*Mesh());
	}
}

void TPZCompEl::BuildConnectList(TPZStack<int> &connectlist) {
	std::set<int> buf;
	int ncon = connectlist.NElements();
	for(int i = 0; i < ncon; i++) {
		buf.insert(connectlist[i]);
	}
	BuildConnectList(buf);
	ncon = buf.size();
	connectlist.Resize(ncon);
	std::set<int>::iterator it = buf.begin();
	for(int i = 0; i < ncon; i++,it++) 
	{
		connectlist[i] = *it;
	}
}

void TPZCompEl::BuildConnectList(std::set<int> &connectlist) {
	std::set<int> additional;
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

void TPZCompEl::SetIndex(int index) {
	/*   int i=0; */
	/*   int numel = fMesh->NElements(); */
	/*   TPZAdmChunkVector<TPZCompEl *> &vec = fMesh->ElementVec(); */
	/*   while ( i < numel ) { */
	/*     if (vec[i] == this) break; */
	/*     i++; */
	/*   } */
	/*   return i; */
	fIndex = index;
	return;
}

int TPZCompEl::NEquations(){
	int i=0;
	int numeq=0;
	for (;i<NConnects(); i++){
		TPZConnect &df = Connect(i);
		if(df.HasDependency() || df.IsCondensed() || !df.NElConnected() || df.SequenceNumber() == -1) continue;
		int seqnum = df.SequenceNumber();
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
	if(ref) ref->SetReference(this);
}

void TPZCompEl::CalcResidual(TPZElementMatrix &ef){
	TPZElementMatrix ek(this->Mesh(), TPZElementMatrix::EK);
	CalcStiff(ek,ef);
}

TPZGeoEl * TPZCompEl::GetRefElPatch(){
	std::stringstream sout;
	sout << "Obtendo elemento geometrico de referencia " << Index() << endl;
	Print(sout);
	TPZGeoEl *ref = Reference();
	if (!ref) {
		LOGPZ_WARN(logger, "reached a null reference");
		return (0);
	}
	ref->Print(sout);
	TPZStack <TPZGeoEl *> father;
	father.Push(ref);
	while(ref->Father()){
		father.Push(ref->Father());
		ref = ref->Father();
		ref->Print(sout);
	}
	int j;
	LOGPZ_DEBUG(logger, sout.str());
	while(father.NElements()) {
		TPZGeoEl *aux = father.Pop();
		for (j=0; j<aux->NSides(); j++){
			int sidedimension = aux->SideDimension(j);
			if(!sidedimension){
				continue;
			}
			TPZGeoElSide side(aux,j);
			TPZStack <TPZCompElSide> stack;
			side.EqualLevelCompElementList(stack,1,1);
			if(stack.NElements()){
				//cout << " \n \n \n ==================================\n ================================\nElemento PatchReference\n";
				//aux->Print();
				LOGPZ_INFO(logger, "Exing GetRefElPatch");
				return aux;
			}
		}
	}
	//cout << " \n \n \n ==================================\n ================================\nElemento PatchReference falho\n";
	//Reference()->Print();
	LOGPZ_WARN(logger, "Exing GetRefElPatch - Elemento PatchReference falho");
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
	//  georeftest = 0;
}

TPZCompElSide::TPZCompElSide(const TPZCompElSide &celside) {
	fEl = celside.fEl;
	fSide   = celside.fSide;
	//  georeftest = celside.georeftest;
}

TPZCompElSide::TPZCompElSide(TPZCompEl *cel,int side) {
	fEl = cel;
	fSide   = side;
	//  if(cel) georeftest = cel->Reference();
	//  else georeftest = 0;
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
	int exnel = expandvec.NElements();
	for (int i=0;i<exnel;i++) {
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
					int l=0;
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
	//int lowidindex = neighbourset.NElements()-1;
	//  TPZGeoElSide neighbour = gelside.Neighbour(),
	TPZGeoElSide lowidneigh(gelside);
	//if(!neighbour.Element()) return expandvec;
	int lowid = gelside.Id();
	int in = 0, nneigh = neighbourset.NElements()-1;
	while(in < nneigh) {
		TPZCompEl *ref = neighbourset[in].Reference().Element();
		if(neighbourset[in].Id() < lowid && ref && (!onlyinterpolated || dynamic_cast<TPZInterpolatedElement*>(ref)    )) {
			lowidneigh = neighbourset[in];
			lowid = lowidneigh.Id();
			//lowidindex = in;
		}
		in++;
	}
	//   while(neighbour!=gelside) {     //pode ser um viz. inativo
	//     TPZCompEl *ref = neighbour.Reference().Element();
	//     if(neighbour.Id() < lowid && ref && (!onlyinterpolated || ref->IsInterpolated())) {
	//       lowidneigh = neighbour;
	//       lowid = lowidneigh.Id();
	//     }
	//     neighbour = neighbour.Neighbour();
	//   }//for
	return lowidneigh.Reference();
}

void TPZCompElSide::RemoveConnectDuplicates(TPZStack<TPZCompElSide> &expandvec){
	int i,k;
	int nelems = expandvec.NElements();
	TPZStack<TPZCompElSide> locexpand;
	for(i=0;i<nelems;i++) locexpand.Push(expandvec[i]);
	//TPZCompElSide store[100];
	//TPZStack<TPZCompElSide> locexpand(store,nelems);
	//for(i=0;i<nelems;i++) locexpand[i] = expandvec[i];
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
int TPZCompElSide::ConnectIndex() const
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


/**
 * Identify the material object associated with the element
 */
TPZAutoPointer<TPZMaterial> TPZCompEl::Material() const
{
    TPZAutoPointer<TPZMaterial> result;
    if(fMesh && Reference()) result = fMesh->FindMaterial(Reference()->MaterialId());
    return result;
}

/**
 * Verify if the material associated with the element is contained in the set
 */
bool TPZCompEl::HasMaterial(const std::set<int> &materialids)
{
	TPZAutoPointer<TPZMaterial> mat = Material();
	if(!mat) return false;
	return materialids.find(mat->Id()) != materialids.end();
}

