//$Id: pzconnect.cpp,v 1.5 2003-11-06 19:14:35 cesar Exp $

//METHODS DEFINITION FOR CLASS NODE


#include "pzconnect.h"
#include "pzerror.h"
#include "pzvec.h"
#include "pzstack.h"
#include "pzcmesh.h"
#include "pzbndcond.h"


TPZConnect::TPZConnect() {
  fSequenceNumber = 0;
  fDependList = 0;
  fNElConnected = 0;
  fOrder = -1;
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
  fOrder = copy.fOrder;
  if(copy.fDependList) fDependList = new TPZDepend(*copy.fDependList);
}

void TPZConnect::Print(TPZCompMesh &mesh, ostream & out) {
  out << "TPZConnect : " << "Sequence number = " << fSequenceNumber <<" Order = " << fOrder;
  if(fSequenceNumber > -1)
	  out << "\tNumElCon = " << fNElConnected << " Block size " << mesh.Block().Size(fSequenceNumber);
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
  if(copy.fNext) fNext = new TPZDepend(*fNext);
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

void TPZConnectBC::Print(TPZCompMesh &mesh,ostream &out){
  out << "Connect boundary condition :\n";
  if(fConnect) fConnect->Print(mesh,out);
  if(fBC) fBC->Print(out);
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
