/**
 * @file
 * @brief Contains the implementation of the TPZConnect methods.
 */
//$Id: pzconnect.cpp,v 1.25 2011-05-11 02:45:38 phil Exp $

#include "pzconnect.h"
#include "pzerror.h"
#include "pzvec.h"
#include "pzstack.h"
#include "pzcmesh.h"
#include "pzbndcond.h"
#include "pzsave.h"
#include "pzstream.h"
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.tpzconnect"));
#endif

using namespace std;

TPZConnect::TPZConnect() {
	fSequenceNumber = -1;
	fDependList = 0;
	fNElConnected = 0;
	//fOrder = -1;
    fFlags = 0;
}

TPZConnect::~TPZConnect() {
    if(fDependList) {
		cout << "The dependency list should be NULL\n";
		TPZDepend *ptr = fDependList;
		while(ptr) {
			int depindex = ptr->fDepConnectIndex;
			cout << "Connect index = " << depindex << endl;
			ptr->fDepMatrix.Print(" ",cout);
			ptr = ptr->fNext;
			//      nod.Print(mesh,out);
		}
		delete fDependList;
		fDependList = 0;
    }
}

void TPZConnect::operator=(const TPZConnect &copy) {
	if(this == &copy) return;
	if(fDependList) delete fDependList;
	fSequenceNumber = copy.fSequenceNumber;
	fNElConnected = copy.fNElConnected;
//	fOrder = copy.fOrder;
    fFlags = copy.fFlags;
	if(copy.fDependList) fDependList = new TPZDepend(*copy.fDependList);
}

void TPZConnect::Print(const TPZCompMesh &mesh, std::ostream & out) {
	int orde = fCompose.fOrder;
	int nstate  = fCompose.fNState;
	int nshape  = fCompose.fNShape;
	out << "TPZConnect : " << "Sequence number = " << fSequenceNumber <<"  Order = " << orde << "  NState = " << nstate << "  NShape " << nshape << " IsCondensed " << IsCondensed() << " IsLagrMult " << IsPressure();
	if(fSequenceNumber > -1)
	{
		out << "\tNumElCon = " << fNElConnected << " Block size " << mesh.Block().Size(fSequenceNumber);
		out << " Solution ";
		int ieq;
		for(ieq=0; ieq< mesh.Block().Size(fSequenceNumber); ieq++)
		{
			out << mesh.Block()(fSequenceNumber,0,ieq,0) << ' ';
		}
	}
	
	out << endl;
	if(fDependList) {
		out << "Dependency :\n";
		TPZDepend *ptr = fDependList;
		while(ptr) {
			int depindex = ptr->fDepConnectIndex;
			//TPZConnect &nod = mesh.ConnectVec()[depindex];
			out << "Connect index = " << depindex << endl;
			ptr->fDepMatrix.Print(" ",out);
			ptr = ptr->fNext;
			//      nod.Print(mesh,out);
		}
	}
}

void TPZConnect::Print(TPZCompMesh &mesh, TPZVec<REAL> &cp, std::ostream & out)
{
	out << "TPZConnect : " << "Sequence number = " << fSequenceNumber <<" Order = " << fCompose.fOrder << " NState = " << fCompose.fNState << " NShape " << fCompose.fNShape ;
	out << " coordinate " << cp;
	if(fSequenceNumber > -1)
	{
		out << "\tNumElCon = " << fNElConnected << " Block size " << mesh.Block().Size(fSequenceNumber);
		out << " Solution ";
		int ieq;
		for(ieq=0; ieq< mesh.Block().Size(fSequenceNumber); ieq++)
		{
			out << mesh.Block()(fSequenceNumber,0,ieq,0) << ' ';
		}
	}
	
	out << endl;
	if(fDependList) {
		out << "Dependency :\n";
		TPZDepend *ptr = fDependList;
		while(ptr) {
			int depindex = ptr->fDepConnectIndex;
			//TPZConnect &nod = mesh.ConnectVec()[depindex];
			out << "Connect index = " << depindex << endl;
			ptr->fDepMatrix.Print(" ",out);
			ptr = ptr->fNext;
			//      nod.Print(mesh,out);
		}
	}
}

void TPZConnect::AddDependency(int myindex, int dependindex,TPZFMatrix &depmat,int ipos,int jpos,int isize,int jsize){
	if(dependindex == myindex) return;
	TPZDepend *connect =0;
	if(dependindex == -1) {
		cout << "Chabuuuuu\n";
		int b;
		cin >> b;
	}
	if(fDependList) connect = fDependList->HasDepend(dependindex);
	if(!connect) {
		connect = new TPZDepend(dependindex,depmat,ipos,jpos,isize,jsize);
		connect->fNext = fDependList;
		fDependList = connect;
	} else {
		TPZFMatrix temp(isize,jsize);
		int i,j;
		for(i=0; i<isize; i++) for(j=0; j<jsize; j++) temp(i,j) = depmat(ipos+i,jpos+j);
		
		
		if (temp.Rows()!=connect->fDepMatrix.Rows() || temp.Cols()!=connect->fDepMatrix.Cols()){
			cout << "TPZConnect::Dependency inconsistent \t"
			<< "temp(r,c): (" << temp.Rows() << " , " << temp.Cols() << " )\t fDepMatrix(r,c): ( "
			<< connect->fDepMatrix.Rows() << " , " << connect->fDepMatrix.Cols() << " )\n";
			int a;
			cin >> a;
			return;
		}
		
		
		temp -= connect->fDepMatrix;
		REAL val = Norm(temp);
		if(val > 1.e-6) {
			cout << "TPZConnect::Dependency inconsistent\n";
		}
	}
}

void TPZConnect::RemoveDepend() {
	if(fDependList)
	{
		delete fDependList;
		fDependList = 0;
	}
}

void TPZConnect::RemoveDepend(int myindex, int dependindex) {
	if(dependindex == myindex || !fDependList) return;
	TPZDepend *dep = fDependList->HasDepend(dependindex);
	if(dep) fDependList = fDependList->RemoveDepend(dep);
}

int TPZConnect::DependencyDepth(TPZCompMesh &mesh){
	if(!fDependList) return 0;
	int maxdep = 0;
	TPZDepend *ptr = fDependList;
	while(ptr) {
		int depindex = ptr->fDepConnectIndex;
		TPZConnect &nod = mesh.ConnectVec()[depindex];
		int depth = nod.DependencyDepth(mesh);
		maxdep = depth > maxdep ? depth : maxdep;
		ptr = ptr->fNext;
	}
	maxdep++;
	return maxdep;
}

/**Adds itself and the connects from which it depends to the list
 this method will add a pointer to the current connect to connectlist if it is not already
 there, extending the list is necessary
 this method will also call AddToList for all connects from which this connect depends
 firstfree points to the first unused element of connectlist
 it is assumed that firstfree <= nodelist.capacity()*/
void TPZConnect::AddToList(int myindex, TPZCompMesh &mesh, TPZStack<int> &connectlist){
	int in=0, cap = connectlist.NElements();
	while(in < cap && connectlist[in] != myindex) in++;
	if(in == cap) connectlist.Push(myindex);
	// this inserts the node in the list and increments the pointer firstfree
	TPZDepend *dp = fDependList;
	while(dp) {
		TPZConnect &depc = mesh.ConnectVec()[dp->fDepConnectIndex];
		depc.AddToList(dp->fDepConnectIndex,mesh,connectlist);
		dp = dp->fNext;
	}
}

void TPZConnect::AddToList(int myindex, TPZCompMesh &mesh, std::set<int> &connectlist){
	connectlist.insert(myindex);
	// this inserts the node in the list and increments the pointer firstfree
	TPZDepend *dp = fDependList;
	while(dp) {
		TPZConnect &depc = mesh.ConnectVec()[dp->fDepConnectIndex];
		depc.AddToList(dp->fDepConnectIndex,mesh,connectlist);
		dp = dp->fNext;
	}
}

void TPZConnect::SetDependenceOrder(int myindex, TPZCompMesh &mesh, int CurrentOrder,TPZVec<int> &ConnectList,TPZVec<int> &DependenceOrder) {
	int in=0,cap = ConnectList.NElements();
	// identify where the current node is in the list
	while(in<cap && ConnectList[in] != myindex) in++;
	if(in== cap) {
		cout << "TPZConnect::SetDependenceOrder node not encountered in list\n";
		return;
	}
	DependenceOrder[in] = (DependenceOrder[in] < CurrentOrder) ?
    CurrentOrder : DependenceOrder[in];
	TPZDepend *dl = fDependList;
	while(dl) {
		int depindex = dl->fDepConnectIndex;
		TPZConnect &depc = mesh.ConnectVec()[depindex];
		depc.SetDependenceOrder(depindex,mesh,CurrentOrder+1,ConnectList,DependenceOrder);
		// call SetDependenceOrder recursively on the nodes from which the current node depends
		dl = dl->fNext;
		// move to the next item in the list
	}
}

TPZConnect::TPZDepend::TPZDepend(int dependindex,TPZFMatrix &depmat,int ipos,int jpos, int isize, int jsize) :
fDepMatrix(isize,jsize) {
	fDepConnectIndex = dependindex;
	int i,j;
	for(i=0; i<isize; i++) for(j=0; j<jsize; j++) fDepMatrix(i,j) = depmat(ipos+i,jpos+j);
	fNext = 0;
}

TPZConnect::TPZDepend::TPZDepend(const TPZDepend &copy) : fDepConnectIndex(copy.fDepConnectIndex),
fDepMatrix(copy.fDepMatrix), fNext(0) {
	if(copy.fNext) fNext = new TPZDepend(*copy.fNext);
}

TPZConnect::TPZDepend::TPZDepend(int connectindex) : fDepConnectIndex(connectindex),
fDepMatrix(),fNext(0)
{
}

TPZConnect::TPZDepend::~TPZDepend() {
	if(fNext) delete fNext;
}

TPZConnect::TPZDepend *TPZConnect::TPZDepend::HasDepend(int depindex) {
	TPZDepend *ptr = this;
	while(ptr && ptr->fDepConnectIndex != depindex) ptr = ptr->fNext;
	return ptr;
}

TPZConnect::TPZDepend *TPZConnect::TPZDepend::RemoveDepend(TPZDepend *ptr) {
	TPZDepend *res = this;
	if(this == ptr) {
		res = ptr->fNext;
		this->fNext = 0;
		delete this;
	} else {
		TPZDepend *run = this;
		while(run && run->fNext != ptr) run = run->fNext;
		if(run) {
			run->fNext = run->fNext->fNext;
			ptr->fNext = 0;
			delete ptr;
		}
	}
	return res;
}


int TPZConnect::NDof(TPZCompMesh &mesh) {
	if(fSequenceNumber<0) {
		PZError << "TPZConnect::NDof. Connect is inactive.\n";
		return -1;
	}
	return mesh.Block().Size(fSequenceNumber);
}

int TPZConnect::SequenceNumber() const {
	return fSequenceNumber;
}

int TPZConnect::CheckDependency(int nshape, TPZCompMesh *mesh, int nstate) {
	
	if(HasDependency()) {
		TPZConnect::TPZDepend *first = FirstDepend();
		while(first) {
			int nr = first->fDepMatrix.Rows();
			int nc = first->fDepMatrix.Cols();
			if(nr != nshape) {
				cout << "TPZConnect::CheckDependency inconsistent dependency nshape = " << nshape << " nrows " << nr << endl;
				return -1;
			}
			TPZConnect &c2 = mesh->ConnectVec()[first->fDepConnectIndex];
			if(nc*nstate != c2.NDof(*mesh)) {
				cout << "TPZConnect::CheckDependency inconsistent dependency ndof = " << c2.NDof(*mesh) << " ncols " << nc << endl;
				return -1;
			}
			int c2result = 0;
			if(c2.HasDependency()) {
				c2result = c2.CheckDependency(c2.NDof(*mesh)/nstate,mesh,nstate);
			}
			if(c2result == -1) return c2result;
			first = first->fNext;
		}
	}
	return 0;
}

void TPZConnect::ExpandShape(int cind, TPZVec<int> &connectlist, TPZVec<int> &blocksize, TPZFMatrix &phi, TPZFMatrix &dphi){
	
    if(!fDependList) return;
    int dim = dphi.Rows();
    int locind = 0;
    int ncon = connectlist.NElements();
    int eqloc = 0;
    while(locind < ncon && connectlist[locind] != cind) {
        eqloc += blocksize[locind];
        locind++;
    }
    if(locind == ncon) {
        cout << "TPZConnect::ExpandShape wrong data structure locind\n";
        return;
    }
    TPZDepend *dep = fDependList;
    while(dep) {
        int eqrem = 0;
        int remind = 0;
        while(remind < ncon && connectlist[remind] != dep->fDepConnectIndex) {
            eqrem += blocksize[remind];
            remind++;
        }
        if(remind == ncon) {
            cout << "TPZConnect::ExpandShape wrong data structure remind\n";
            return;
        }
        int rows = dep->fDepMatrix.Rows();
        int cols = dep->fDepMatrix.Cols();
        if(rows != blocksize[locind] || cols != blocksize[remind]) {
            cout << "TPZConnect::ExpandShape wrong data structure sizes\n";
            return;
        }
        int r,c,d;
        for(r=0; r<rows; r++) {
            for(c=0; c<cols; c++) {
                phi(eqrem+c,0) += phi(eqloc+r)*dep->fDepMatrix(r,c);
                for(d=0; d<dim; d++) {
                    dphi(d,eqrem+c) += dphi(d,eqloc+r)*dep->fDepMatrix(r,c);
                }
            }
        }
        dep = dep->fNext;
    }
    int r,d;
    for(r=0; r<blocksize[locind]; r++) {
        phi(eqloc+r,0) = 0.;
        for(d=0; d<dim; d++) {
            dphi(d,eqloc+r) = 0.;
        }
    }
}


/*!
 \fn TPZConnect::TPZDepend::Write(TPZStream &buf, int withclassid)
 */
void TPZConnect::TPZDepend::Write(TPZStream &buf)
{
	buf.Write(&fDepConnectIndex,1);
	fDepMatrix.Write(buf,0);
	if(fNext)
	{
		fNext->Write(buf);
	}
	else
	{
		int min1 = -1;
		buf.Write(&min1,1);
	}
}


/*!
 \fn TPZConnect::TPZDepend::Read(TPZStream &buf, void *context)
 */
void TPZConnect::TPZDepend::Read(TPZStream &buf)
{
	fDepMatrix.Read(buf,0);
	int nextindex;
	buf.Read(&nextindex,1);
	if(nextindex >= 0)
	{
		fNext = new TPZDepend(nextindex);
		fNext->Read(buf);
	}
	else
	{
		fNext = 0;
	}
}

/**
 Save the element data to a stream
 */
void TPZConnect::Write(TPZStream &buf, int withclassid)
{
	buf.Write(&fSequenceNumber,1);
	buf.Write(&fNElConnected,1);
//	buf.Write(&fOrder,1);
    buf.Write(&fFlags,1);
	if(fDependList)
	{
		fDependList->Write(buf);
	} else
	{
		int min1 = -1;
		buf.Write(&min1,1);
	}
}

/**
 Read the element data from a stream
 */
void TPZConnect::Read(TPZStream &buf, void *context)
{
	buf.Read(&fSequenceNumber,1);
	buf.Read(&fNElConnected,1);
//	buf.Read(&fOrder,1);
    buf.Read(&fFlags,1);
	int seq;
	buf.Read(&seq,1);
	if(seq >= 0)
	{
		fDependList = new TPZDepend(seq);
		fDependList->Read(buf);
	} else
	{
		fDependList = 0;
	}
}


/*!
 \fn TPZConnect::CopyFrom(TPZConnect &orig,std::map<int,int> & gl2lcIdx)
 */
void TPZConnect::CopyFrom(TPZConnect &orig,std::map<int,int> & gl2lcIdx)
{
//	fOrder = orig.fOrder;
	fSequenceNumber = orig.fSequenceNumber;
	fNElConnected = orig.fNElConnected;
    fFlags = orig.fFlags;
	TPZDepend * depend = orig.fDependList;
	bool copydepend = true;
	while ( depend )
	{
		if (gl2lcIdx.find(depend->fDepConnectIndex) != gl2lcIdx.end())
		{
			depend = depend->fNext;
		}
		else
		{
			depend = 0;
			copydepend = false;
		}
	}
	if(orig.fDependList && copydepend)
	{
		fDependList = new TPZDepend(-1);
		fDependList->CopyFrom(orig.fDependList,gl2lcIdx);
	}
}


/*!
 \fn TPZConnect::TPZDepend::CopyFrom(TPZDepend *orig,std::map<int,int>& gl2lcIdx)
 */
void TPZConnect::TPZDepend::CopyFrom(TPZDepend *orig,std::map<int,int>& gl2lcIdx)
{
	int loccondepIdx = -1;
	int origdepconIdx = orig->fDepConnectIndex;
	if (gl2lcIdx.find(origdepconIdx) != gl2lcIdx.end()) loccondepIdx = gl2lcIdx[origdepconIdx];
	else
	{
		std::stringstream sout;
		sout << "ERROR in : " << __PRETTY_FUNCTION__
		<< " trying to clone a dependency connect index: " << origdepconIdx
		<< " wich is not in mapped connect indexes!" ;
		LOGPZ_ERROR(logger,sout.str().c_str());
		return;
	}
	fDepConnectIndex = loccondepIdx;

	fDepMatrix = orig->fDepMatrix;
	
	TPZDepend *depend = orig->fNext;
	while (depend)
	{
		if (gl2lcIdx.find(depend->fDepConnectIndex) != gl2lcIdx.end())
		{
			fNext = new TPZDepend(-1);
			fNext->CopyFrom(depend,gl2lcIdx);
			depend = 0;
		}
		else
		{
			depend = depend->fNext;
		}
	}
}

void TPZConnect::BuildConnectList(int index, std::set<int> &indepconnectlist, std::set<int> &depconnectlist, TPZCompMesh &mesh){
	if(fDependList)
	{
		depconnectlist.insert(index);
		TPZDepend *dep = fDependList;
		while(dep)
		{
			TPZConnect &master = mesh.ConnectVec()[dep->fDepConnectIndex];
			master.BuildConnectList(dep->fDepConnectIndex,indepconnectlist,depconnectlist,mesh);
			dep = dep->fNext;
		}
	}
	else
	{
		indepconnectlist.insert(index);
	}
}//void

void TPZConnect::BuildConnectList(TPZStack<int> &connectlist, TPZVec<int> &ConnectIndex, TPZCompMesh &mesh){
	TPZConnect *dfn;
	int dfnindex;
	int nconnects = ConnectIndex.NElements();
	int in;
	for(in = 0; in < nconnects; in++){
		dfnindex = ConnectIndex[in];
		dfn = & (mesh.ConnectVec()[ dfnindex ]);
		dfn->AddToList(dfnindex,mesh,connectlist);
	}//for in
}//void

void TPZConnect::BuildConnectList(std::set<int> &connectlist, std::set<int> &additional, TPZCompMesh &mesh){
	TPZConnect *dfn;
	int dfnindex;
	TPZAdmChunkVector<TPZConnect> &connectvec = mesh.ConnectVec();
	//	int nconnects = additional.size();
	std::set<int>::iterator it;
	for(it = additional.begin() ; it != additional.end(); it++){
		dfnindex = *it;
		dfn = & connectvec[ dfnindex ];
		dfn->AddToList(dfnindex,mesh,connectlist);
	}//for in
}//void

void TPZConnect::BuildDependencyOrder(TPZVec<int> &connectlist, TPZVec<int> &DependenceOrder, TPZCompMesh &mesh) {
	// nodelist (input) : vector which contains pointers to all nodes which
	// are in the dependency chain of the nodes of the element
	int totalnodes = connectlist.NElements();
	DependenceOrder.Resize(totalnodes);
	DependenceOrder.Fill(0,0);
	// initialize the vector which contains the
	// dependency order to zero
	int CurrentOrder = 0;
	// order which is currently processed
	int numnodes_processed = totalnodes;
	// number of nodes processed during the current cycle
	
	while(numnodes_processed) {
		
		numnodes_processed = 0;
		int i;
		for(i=0; i<totalnodes; i++) {
			int dfnindex = connectlist[i];
			TPZConnect &dfn = mesh.ConnectVec()[dfnindex];
			if(dfn.HasDependency() && DependenceOrder[i] == CurrentOrder) {
				dfn.SetDependenceOrder(dfnindex,mesh,CurrentOrder,connectlist,DependenceOrder);
				// this method will fill in the DependenceOrder vector by recursively
				// calling SetDependenceOrder for the nodes upon which dfn depends
				numnodes_processed++;
			}
		}
		CurrentOrder++;
	}
}
