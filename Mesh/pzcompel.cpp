/**
 * @file
 * @brief Contains the implementation of the TPZCompEl methods.
 */

#include "pzcompel.h"
#include "pzgeoel.h"
#include "pzgeoelside.h"
#include "TPZMaterial.h"
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

#include <algorithm>
#include <iterator>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.tpzcompel"));
static LoggerPtr loggerSide(Logger::getLogger("pz.mesh.tpzcompelside"));
#endif

void TPZCompEl::CalcBlockDiagonal(TPZStack<int64_t> &connectlist, TPZBlockDiagonal<STATE> & blockdiag) {
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
            int64_t conind = ek.fConstrConnect[b];
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
            int64_t conind = ek.fConnect[b];
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

TPZCompEl::TPZCompEl() : fMesh(0), fIndex(-1), fReferenceIndex(-1), fIntegrationRule(0) {
}

TPZCompEl::TPZCompEl(TPZCompMesh &mesh, TPZGeoEl *ref, int64_t &index) : fIntegrationRule(0) {
    fMesh = &mesh;
    index = mesh.ElementVec().AllocateNewElement();
    mesh.ElementVec()[index] = this;
    fIndex = index;
    fReferenceIndex = (ref == 0) ? -1 : ref->Index();
}

TPZCompEl::TPZCompEl(TPZCompMesh &mesh, const TPZCompEl &copy): fIntegrationRule(0) {
    fMesh = &mesh;
    int64_t index = copy.fIndex;
    if(index >= 0) mesh.ElementVec()[index] = this;
    fIndex = index;
    fReferenceIndex = copy.fReferenceIndex;
    if (copy.fIntegrationRule) {
        fIntegrationRule = copy.fIntegrationRule->Clone();
    }
}

TPZCompEl::TPZCompEl(TPZCompMesh &mesh, const TPZCompEl &copy, int64_t &index) : fIntegrationRule(0) {
    fMesh = &mesh;
    index = mesh.ElementVec().AllocateNewElement();
    if(index >= 0) mesh.ElementVec()[index] = this;
    fIndex = index;
    fReferenceIndex = copy.fReferenceIndex;
    if (copy.fIntegrationRule) {
        fIntegrationRule = copy.fIntegrationRule->Clone();
    }
}

TPZCompEl::TPZCompEl(TPZCompMesh &mesh, const TPZCompEl &copy, std::map<int64_t,int64_t> &gl2lcElMap) : fIntegrationRule(0)
{
    fMesh = &mesh;
    if (gl2lcElMap.find(copy.fIndex) == gl2lcElMap.end())
    {
        std::stringstream sout;
        sout << "ERROR in - " << __PRETTY_FUNCTION__
        << " original element index: " << copy.fIndex << " is not mapped!\n"
        << "Map content: ";
        std::map<int64_t,int64_t>::iterator it;
        for (it=gl2lcElMap.begin();it!=gl2lcElMap.end();it++) sout << " ( " << it->first << " | " << it->second << " ) ;";
        LOGPZ_ERROR (logger,sout.str().c_str());
        DebugStop();
    }
    int64_t index = gl2lcElMap[copy.fIndex];
    if(index >= 0) mesh.ElementVec()[index] = this;
    fIndex = index;
    fReferenceIndex = copy.fReferenceIndex;
    if (copy.fIntegrationRule) {
        fIntegrationRule = copy.fIntegrationRule->Clone();
    }
}

TPZCompEl::~TPZCompEl() {
    int64_t index = Index();
    if (index != -1){
        if (fMesh->ElementVec()[index] == this) {
            fMesh->ElementVec()[index] = 0;
            fMesh->ElementVec().SetFree(index);
        }
    }
#ifdef PZDEBUG
    TPZGeoEl *gel = Reference();
    if (gel && gel->Reference()) {
        DebugStop();
    }
#endif
    fIndex = -1;
    fReferenceIndex = -1;
    fMesh = 0;
    if (fIntegrationRule) {
        delete fIntegrationRule;
    }
}

MElementType TPZCompEl::Type() {
    LOGPZ_WARN(logger, "Type unknown");
    return ENoType;
}

void TPZCompEl::LoadSolution() {
    // an element without mesh is a free element
    if(!Mesh() || !HasDependency()) return;
    TPZStack<int64_t> connectlist;
    int totalconnects;
    BuildConnectList(connectlist);
    totalconnects = connectlist.NElements();
    TPZManVector<int> dependenceorder(totalconnects);
    TPZConnect::BuildDependencyOrder(connectlist,dependenceorder,*this->Mesh());
    //	TPZMaterial * mat = Material();
    //	if(!mat) {
    //		LOGPZ_WARN(logger, "Exiting LoadSolution because a null material was reached.");
    //		return;
    //	}
    //int numstate = mat->NStateVariables();
    //TPZBlock<REAL> &block = Mesh()->Block();
    TPZBlock<STATE> &block = Mesh()->Block();
    //TPZFMatrix<REAL> &MeshSol = Mesh()->Solution();
    TPZFMatrix<STATE> &MeshSol = Mesh()->Solution();
    int maxdep = 0;
    int in;
    int64_t iv,jv,idf;
    STATE coef;
    for(in=0;in<totalconnects;in++)
        maxdep = (maxdep < dependenceorder[in]) ? dependenceorder[in] : maxdep;
    int current_order = maxdep-1;
    while(current_order >= 0) {
        for(in=0; in<totalconnects; in++) {
            if(dependenceorder[in] != current_order) continue;
            TPZConnect *dfn = &Mesh()->ConnectVec()[connectlist[in]];
            if(!dfn->HasDependency()) continue;
            int64_t bl = dfn->SequenceNumber();
            int nvar = block.Size(bl);
            int numstate = dfn->NState(); //numstate eh fornecida pelo connect
            //         int numshape = nvar/numstate;
            TPZConnect::TPZDepend *dep = dfn->FirstDepend();
            int64_t blpos = block.Position(bl);
            for(iv=0; iv<nvar; iv++) MeshSol(blpos+iv, 0) = 0.;
            while(dep) {
                int64_t depconindex = dep->fDepConnectIndex;
                TPZConnect &depcon = Mesh()->ConnectVec()[depconindex];
                int64_t depseq = depcon.SequenceNumber();
                int numdepvar = block.Size(depseq);
                int64_t depseqpos = block.Position(depseq);
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
    std::list<TPZOneShapeRestraint> mylist = this->GetShapeRestraints();
    for (std::list<TPZOneShapeRestraint>::iterator it = mylist.begin(); it != mylist.end(); it++)
    {
        int64_t connectdest = it->fFaces[0].first;
        int64_t seqnumdest = Mesh()->ConnectVec()[connectdest].SequenceNumber();
        int destidf = it->fFaces[0].second;
        REAL mult = -1./it->fOrient[0];
#ifdef PZDEBUG
        STATE prevval = Mesh()->Block()(seqnumdest,0,destidf,0);
#endif
        Mesh()->Block()(seqnumdest,0,destidf,0) = 0.;
        for (int i=1; i<4; i++) {
            int64_t connectindex = it->fFaces[i].first;
            int64_t seqnum = Mesh()->ConnectVec()[connectindex].SequenceNumber();
            int idf = it->fFaces[i].second;
            STATE val = Mesh()->Block()(seqnum,0,idf,0);
            REAL multorig = it->fOrient[i];
            Mesh()->Block()(seqnumdest,0,destidf,0) += mult*multorig*val;
        }
#ifdef PZDEBUG
        STATE finalval = Mesh()->Block()(seqnumdest,0,destidf,0);
        finalval -= prevval;
#endif
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
        int64_t connectindex = ConnectIndex(i);
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
    out << "\nfReferenceIndex " << fReferenceIndex;
	if(this->Reference())
	{
		out << "\nCenter coordinate: ";
		TPZVec< REAL > centerMaster( this->Reference()->Dimension(),0. );
		TPZVec< REAL > centerEuclid( 3,0.);
		this->Reference()->CenterPoint(this->Reference()->NSides()-1,centerMaster);
		this->Reference()->X(centerMaster,centerEuclid);
		out << centerEuclid;
	}
	if(this->Material())
	{
		out << "\nMaterial id " << this->Material()->Id() << "\n";
	}
	else {
		out << "\nNo material\n";
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

void TPZCompEl::PrintSolution(TPZVec<REAL> &point,const char *varname,std::ostream &s) {
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
    for(int64_t i=0; i<sol.NElements(); i++) {
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

void TPZCompEl::PrintTitle(const char *varname,std::ostream &s) {
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

inline void TPZCompEl::Divide(int64_t index, TPZVec<int64_t> &subindex, int interpolate) {
    subindex.Resize(0);
    LOGPZ_WARN(logger,"TPZCompEl::Divide called");
}

void TPZCompEl::EvaluateError(std::function<void(const TPZVec<REAL> &loc,TPZVec<STATE> &val,TPZFMatrix<STATE> &deriv)> fp,
                              TPZVec<REAL> &/*errors*/, bool store_error) {
    LOGPZ_WARN(logger, "EvaluateError is called.");
    DebugStop();
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

void TPZCompEl::BuildConnectList(std::set<int64_t> &indepconnectlist,
                                 std::set<int64_t> &depconnectlist) {
    const int ncon = this->NConnects();
    for(int i = 0; i < ncon; i++) {
        int64_t conind = ConnectIndex(i);
        Connect(i).BuildConnectList(conind, indepconnectlist,depconnectlist,*Mesh());
    }
    std::list<TPZOneShapeRestraint> mylist = GetShapeRestraints();
    for (std::list<TPZOneShapeRestraint>::iterator it = mylist.begin(); it != mylist.end(); it++) {
        for (int i=0; i<4; i++) {
            int64_t conind = it->fFaces[i].first;
            TPZConnect &c = Mesh()->ConnectVec()[conind];
            c.BuildConnectList(conind, indepconnectlist, depconnectlist, *Mesh());
        }
    }
}

void TPZCompEl::BuildConnectList(TPZStack<int64_t> &connectlist) const {
    int64_t ncon = connectlist.NElements();
    if (ncon) {
        std::sort(&connectlist[0], &connectlist[0]+ncon);
    }
    TPZAdmChunkVector<TPZConnect> &connectvec = Mesh()->ConnectVec();
    std::list<TPZOneShapeRestraint> rest = GetShapeRestraints();
    int64_t nconloc = NConnects();
    bool hasdependency = false;
    if (ncon == 0) {
        connectlist.Resize(nconloc);
        for(int64_t i = 0; i < nconloc; i++) {
            connectlist[i] = this->ConnectIndex(i);
            if (connectlist[i] == -1) continue;
            if (connectvec[connectlist[i]].HasDependency()) {
                hasdependency = true;
            }
        }
        if (hasdependency == false && rest.size() == 0) {
            return;
        }
        std::sort(&connectlist[0], &connectlist[0]+nconloc);
        ncon = nconloc;
        nconloc = 0;
    }
    TPZManVector<int64_t> localcon(nconloc);
    for(int64_t i = 0; i < nconloc; i++) {
        localcon[i] = this->ConnectIndex(i);
        if (connectvec[localcon[i]].HasDependency()) {
            hasdependency = true;
        }
    }
    if (nconloc) {
        std::sort(&localcon[0], &localcon[0]+nconloc);
    }
    
    std::set<int64_t> buf;
    if (ncon > 0 && nconloc > 0)
    {
        std::set_union(&connectlist[0],&connectlist[0]+ncon,&localcon[0],&localcon[0]+nconloc,std::inserter(buf, buf.begin()));
    }
    else if (ncon > 0)
    {
        buf = std::set<int64_t>(&connectlist[0],&connectlist[0]+ncon);
    }
    else if(nconloc > 0)
    {
        buf = std::set<int64_t>(&localcon[0],&localcon[0]+nconloc);
    }
    std::set<int64_t> buf2;
    std::list<TPZOneShapeRestraint> mylist = GetShapeRestraints();
    for (std::list<TPZOneShapeRestraint>::iterator it = mylist.begin(); it != mylist.end(); it++) {
        for (int i=0; i<4; i++) {
            buf2.insert(it->fFaces[i].first);
        }
    }
    nconloc = NConnects();
    for(int64_t i = 0; i < nconloc; i++) {
        TPZConnect &c = Connect(i);
        if (c.HasDependency()) {
            TPZConnect::TPZDepend * dep= c.FirstDepend();
            while(dep)
            {
                buf2.insert(dep->fDepConnectIndex);
                dep = dep->fNext;
            }
        }
    }
    TPZConnect::BuildConnectList(buf, buf2, *Mesh());
    if (buf.size() != connectlist.size())
    {
        connectlist.Resize(buf.size());
        std::copy(buf.begin(), buf.end(), &connectlist[0]);
    }
}

void TPZCompEl::BuildConnectList(std::set<int64_t> &connectlist) {
    int64_t nconloc = NConnects();
    TPZManVector<int64_t> localcon(nconloc);
    for(int64_t i = 0; i < nconloc; i++) {
        localcon[i] = this->ConnectIndex(i);
    }
    if (nconloc) {
        std::sort(&localcon[0], &localcon[0]+nconloc);
    }
    
    std::set<int64_t> buf;
    if (nconloc > 0)
    {
        std::set_union(connectlist.begin(),connectlist.end(),&localcon[0],&localcon[0]+nconloc,std::inserter(buf, buf.begin()));
    }
    else
    {
        buf = connectlist;
    }
    std::set<int64_t> buf2;
    std::list<TPZOneShapeRestraint> mylist = GetShapeRestraints();
    for (std::list<TPZOneShapeRestraint>::iterator it = mylist.begin(); it != mylist.end(); it++) {
        for (int i=0; i<4; i++) {
            buf2.insert(it->fFaces[i].first);
        }
    }
    for(std::set<int64_t>::iterator it=buf.begin(); it != buf.end(); it++) {
        TPZConnect &c = Mesh()->ConnectVec()[*it];
        if (c.HasDependency()) {
            TPZConnect::TPZDepend * dep= c.FirstDepend();
            buf2.insert(dep->fDepConnectIndex);
        }
    }
    TPZConnect::BuildConnectList(buf, buf2, *Mesh());
    connectlist = buf;
}

int TPZCompEl::HasDependency() {
    int nconnects = NConnects();
    int in;
    for(in=0; in<nconnects; in++) if(Connect(in).HasDependency()){
        return 1;
    }
    if (GetShapeRestraints().size()) {
        return 1;
    }
    return 0;
}

void TPZCompEl::SetIndex(int64_t index) {
    fIndex = index;
    return;
}

int TPZCompEl::NEquations() {
    int numeq=0;
    for (int i=0;i<NConnects(); i++){
        TPZConnect &df = Connect(i);
        if(df.HasDependency() || df.IsCondensed() || !df.NElConnected() || df.SequenceNumber() == -1) continue;
        int64_t seqnum = df.SequenceNumber();
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
//    if(!this) {
//        LOGPZ_ERROR(logger,"TPZCompMesh::MaximumRadiusOfEl() null element");
//    }
    
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
    
//    if(!this) LOGPZ_INFO(logger,"TPZCompMesh::LesserEdgeOfEl null element");   ///Jorge 2017 If object exists this can not be NULL
    
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
void TPZCompEl::Write(TPZStream &buf, int withclassid) const
{
    TPZPersistenceManager::WritePointer(fMesh, &buf);
    buf.Write(&fIndex,1);
    buf.Write(&fReferenceIndex,1);
}

/**
 Read the element data from a stream
 */
void TPZCompEl::Read(TPZStream &buf, void *context)
{
    fMesh = dynamic_cast<TPZCompMesh *>(TPZPersistenceManager::GetInstance(&buf));
    buf.Read(&fIndex,1);
    buf.Read(&fReferenceIndex,1);
}

void TPZCompEl::SetOrthogonalFunction(void (*orthogonal)(REAL x,int num,
                                                         TPZFMatrix<REAL> & phi,TPZFMatrix<REAL> & dphi)) {
    pzshape::TPZShapeLinear::fOrthogonal = orthogonal;
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
    int64_t exnel = expandvec.NElements();
    for (int64_t i=0;i<exnel;i++) {
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
                    int64_t l=0;
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
    int64_t i,k;
    int64_t nelems = expandvec.NElements();
    TPZStack<TPZCompElSide> locexpand;
    for(i=0;i<nelems;i++) locexpand.Push(expandvec[i]);
    expandvec.Resize(0);
    for(k=0;k<nelems;k++){
        //TPZCompEl *kel = locexpand[k].Element();
        TPZInterpolatedElement *kel = dynamic_cast<TPZInterpolatedElement *> (locexpand[k].Element());
        if(!kel) continue;
        int kside = locexpand[k].Side();
        i=k+1;
        while(i<nelems){
            //TPZCompEl *iel = locexpand[i].Element();
            TPZInterpolatedElement *iel = dynamic_cast<TPZInterpolatedElement *> (locexpand[i].Element());
            int iside = locexpand[i].Side();
            //if(iel && kel->ConnectIndex(kside) == iel->ConnectIndex(iside))
            
            if(iel)
            {
                int a = kel->MidSideConnectLocId(kside);
                int b = iel->MidSideConnectLocId(iside);
                int64_t connecta = kel->ConnectIndex(a);
                int64_t connectb = iel->ConnectIndex(b);
                if(connecta == connectb)
                {
                    locexpand[i] = TPZCompElSide();
                }
            }
            i++;
        }
    }
    for(i=0;i<nelems;i++)
        if(locexpand[i].Element()) expandvec.Push(locexpand[i]);
}

/// Return the index of the middle side connect alon fSide
int64_t TPZCompElSide::ConnectIndex() const
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
        int64_t locconnectindex = ConnectIndex(ic);
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

/**
 * @brief Compute the integral of a variable defined by the string
 */
TPZVec<STATE> TPZCompEl::IntegrateSolution(const std::string &varname, const std::set<int> &matids)
{
    TPZMaterial *mat = Material();
    if (mat) {
        int id = mat->Id();
        if (matids.find(id) != matids.end())
        {
            int varindex = mat->VariableIndex(varname);
            if (varindex != -1) {
                return IntegrateSolution(varindex);
            }
            else
            {
                // the indicated matid does not support the variable name??
                DebugStop();
            }
        }
    }
    TPZVec<STATE> result;
    return result;
}

int TPZCompEl::ComputeIntegrationOrder() const {
    DebugStop();
	return 0;
}

void TPZCompEl::SetIntegrationRule(TPZIntPoints *intrule)
{
    if (fIntegrationRule) {
        delete fIntegrationRule;
    }
    fIntegrationRule = intrule;
}

int TPZCompEl::StaticClassId() {
    return Hash("TPZCompEl");
}


int TPZCompEl::ClassId() const{
    return StaticClassId();
}

void TPZCompEl::SetCreateFunctions(TPZCompMesh* mesh) {
    mesh->SetAllCreateFunctionsContinuous();
}

TPZGeoEl* TPZCompEl::Reference() const {
    if (fMesh == NULL || fMesh->Reference() == NULL) return NULL;
    return (fReferenceIndex == -1) ? NULL : fMesh->Reference()->ElementVec()[fReferenceIndex];
}
