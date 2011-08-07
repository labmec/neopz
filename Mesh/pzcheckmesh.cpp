/**
 * @file
 * @brief Contains the implementation of the TPZCheckMesh methods.
 */
//$Id: pzcheckmesh.cpp,v 1.7 2007-06-27 18:18:04 cesar Exp $

#include "pzcheckmesh.h"
#include "pzgeoelside.h"
#include "pzstack.h"
#include "pzintel.h"
#include "pzmaterial.h"
//#include "pzcmesh.h"
//#include "pzcompel.h"

using namespace std;

TPZCheckMesh::TPZCheckMesh(TPZCompMesh *mesh, ostream *out) {
	fMesh = mesh;
	fOut = out;
	fNState = 1;
}

void TPZCheckMesh::BuildDependList(int connect, TPZStack<int> &dependlist) {
	int nconnects = fMesh->ConnectVec().NElements();
	int ic;
	for(ic=0; ic<nconnects; ic++) {
		TPZConnect::TPZDepend *dep = fMesh->ConnectVec()[ic].FirstDepend();
		if(dep && dep->HasDepend(connect)) {
			dependlist.Push(ic);
			continue;
		}
	}
}

TPZCompElSide TPZCheckMesh::FindElement(int connect) {
	int nelem = fMesh->ElementVec().NElements();
	int el;
	for(el=0; el<nelem; el++) {
		TPZCompEl *cel = fMesh->ElementVec()[el];
		if(!cel) continue;
		int nc = cel->NConnects();
		int ic;
		for(ic=0; ic<nc; ic++) {
			int c = cel->ConnectIndex(ic);
			if(c== connect) {
				return TPZCompElSide(cel,ic);
			}
		}
	}
	TPZCompElSide nullside;
	return nullside;
}


int TPZCheckMesh::VerifyConnect(int connect) {
	
	int check = 1;
	TPZStack<int> deplist;
	BuildDependList(connect,deplist);
	int ndep = deplist.NElements();
	int id;
	for(id=0; id<ndep; id++) {
		TPZCompElSide small = FindElement(deplist[id]);
		TPZCompElSide large = small.LowerLevelElementList(1);
		if(!large.Exists()) {
			check = 0;
			(*fOut) << "VerifyConnect of " << connect << " inconsistent\n";
			continue;
		}
		TPZCompEl *cel = large.Element();
		int ncl = cel->NConnects();
		int icl;
		for(icl=0; icl<ncl; icl++) {
			if(cel->ConnectIndex(icl) == connect) break;
		}
		if(icl == ncl) {
			check = 0;
			(*fOut) << "VerifyConnect of " << connect << " inconsistent\n";
		}
	}
	return check;
}

void TPZCheckMesh::DependencyReport(int connect) {
	
	(*fOut) << "DependencyReport for connect " << connect << endl;
	TPZCompElSide large = FindElement(connect);
	DependencyReport(connect, large);
}

void TPZCheckMesh::DependencyReport(int connect, TPZCompElSide & large) {
	
	TPZStack<TPZCompElSide> small;
	large.HigherLevelElementList(small,1,1);
	int nsmall = small.NElements();
	int is;
	for(is=0; is<nsmall; is++) {
		int conind = small[is].Element()->ConnectIndex(small[is].Side());
		(*fOut) << "Connect " << conind << " may depend on " << connect << endl;
	}
	if(large.Reference().Dimension() == 1) {
		TPZStack<TPZCompElSide> equal;
		large.EqualLevelElementList(equal,1,0);
		equal.Push(large);
		TPZStack<TPZCompElSide> highdim;
		int neq = equal.NElements();
		int eq;
		for(eq=0; eq<neq; eq++) equal[eq].HigherDimensionElementList(highdim,1,0);
		large.HigherDimensionElementList(highdim,1,0);
		int nhigh = highdim.NElements();
		int ih;
		for(ih=0; ih<nhigh; ih++) {
			DependencyReport(connect,highdim[ih]);
		}
	}
}


int TPZCheckMesh::VerifyAllConnects(){
	int ncon = fMesh->NConnects();
	int i;
	int check = 1;
	for (i=0; i<ncon; i++){
		//(*fOut) << "Startint to check consistency for connect " << i << endl;
		if (!VerifyConnect(i)){
			check = 0;
			(*fOut) << "Check failed for connect: " << i << endl;
		}
	}
	return check;
}

int TPZCheckMesh::CheckDimensions()
{
	int check = 0;
	check = CheckElementShapeDimension();
	int check2 = CheckConstraintDimension();
	return (check || check2);
	
}

int TPZCheckMesh::CheckElementShapeDimension()
{
	int check = 0;
	int nelem = fMesh->ElementVec().NElements();
	int el;
	for(el=0; el<nelem; el++) {
		TPZCompEl *cel = fMesh->ElementVec()[el];
		if(!cel) continue;
		TPZInterpolatedElement *cint = dynamic_cast<TPZInterpolatedElement *> (cel);
		if(!cint) continue;
		int nc = cint->NConnects();
		int ic;
		for(ic=0; ic<nc; ic++) {
			int cind = cint->ConnectIndex(ic);
			int nshape = cint->NConnectShapeF(ic);
			int seqnum = fMesh->ConnectVec()[cind].SequenceNumber();
			int blsize = fMesh->Block().Size(seqnum);
			if(nshape*fNState != blsize) {
				cout << "TPZCheckMesh::CheckElement.. el = " << el << " connect = "
				<< cind << " nshape = " << nshape << " blsize " << blsize << endl;
				cint->Print(cout);
				check = 1;
			}
		}
	}
	return check;
}

int TPZCheckMesh::CheckConstraintDimension()
{
	int check = 0;
	int ncon = fMesh->ConnectVec().NElements();
	int ic;
	for(ic=0; ic<ncon; ic++) {
		TPZConnect &con = fMesh->ConnectVec()[ic];
		int iseqnum = con.SequenceNumber();
		int iblsize = fMesh->Block().Size(iseqnum);
		if(con.HasDependency()) {
			TPZConnect::TPZDepend *dep = con.FirstDepend();
			while(dep) {
				int r = dep->fDepMatrix.Rows();
				int c = dep->fDepMatrix.Cols();
				int jcind = dep->fDepConnectIndex;
				TPZConnect &jcon = fMesh->ConnectVec()[jcind];
				int jseqnum = jcon.SequenceNumber();
				int jblsize = fMesh->Block().Size(jseqnum);
				if(r*fNState != iblsize || c*fNState != jblsize) {
					cout << "TPZCheckMesh::CheckConstraintDi.. ic = " << ic
					<< " depends on " << jcind << " r " << r << " iblsize "
					<< iblsize  << " c " << c << " jblsize " << jblsize << endl;
					check = 1;
				}
				dep = dep->fNext;
			}
		}
	}
	return check;
}

int TPZCheckMesh::CheckConnectOrderConsistency() {
	
	int nel = fMesh->ElementVec().NElements();
	int nstate = (fMesh->MaterialVec().begin()->second)->NStateVariables();
	int iel;
	for(iel = 0; iel<nel; iel++) {
		TPZCompEl *cel = fMesh->ElementVec()[iel];
		if(!cel) continue;
		TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
		if(!intel->CheckElementConsistency()) {
			intel->Print();
			return iel;
		}
		int nc = intel->NConnects();
		int ic;
		for(ic = 0; ic<nc; ic++) {
			TPZConnect &c = intel->Connect(ic);
			int nshape = intel->NConnectShapeF(ic);
			if(c.CheckDependency(nshape, fMesh, nstate) == -1) {
				cout << "TPZCheckMesh inconsistent dependency" << endl;
				intel->Print();
				return iel;
			}
		}
	}
	return -1;
}

/*
 int TPZCheckMesh::CheckConnectOrderConsistency() {
 
 int nel = fMesh->ElementVec().NElements();
 int iel;
 for(iel = 0; iel<nel; iel++) {
 TPZCompEl *cel = fMesh->ElementVec()[iel];
 if(!cel) continue;
 TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
 if(!intel->CheckElementConsistency()) {
 intel->Print();
 return iel;
 }
 int nc = intel->NConnects();
 int ic;
 for(ic = 0; ic<nc; ic++) {
 TPZConnect &c = intel->Connect(ic);
 int nshape = intel->NConnectShapeF(ic);
 if(c.HasDependency()) {
 TPZConnect::TPZDepend *first = c.FirstDepend();
 int nr = first->fDepMatrix.Rows();
 if(nr != nshape) {
 cout << "TPZCheckMesh inconsistent dependency nshape = " << nshape << " nrows " << nr << endl;
 intel->Print();
 return iel;
 }
 }
 }
 }
 return -1;
 }
 */
