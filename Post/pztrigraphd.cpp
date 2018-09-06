/**
 * @file
 * @brief Contains the implementation of the TPZGraphElTd methods. 
 */

#include "pztrigraphd.h"
#include "pzgraphmesh.h"

using namespace std;

TPZGraphElTd::TPZGraphElTd(TPZCompEl *c, TPZGraphMesh *g) : TPZGraphEl(c,g,fConnect){
	
}

int TPZGraphElTd::NConnects(){
	return 1;
}

int TPZGraphElTd::NElements(){
	int res = fGraphMesh->Res();
	int imax = (1<<res);
	return imax*imax;
}

int TPZGraphElTd::NPoints(TPZGraphNode *n){
	int res = fGraphMesh->Res();
	int imax = (1<<res)+1;
	return ((imax+1)*imax)/2;
}

void TPZGraphElTd::Connectivity(TPZDrawStyle st){
	int res = fGraphMesh->Res();
	int imax = 1 << res;
	ostream &out = fGraphMesh->Out();
    int64_t ip = fId;
	if(st == EV3DStyle) ip++;
	TPZVec<int> co0(3,0), co1(3,0), co2(3,0);
	for(int j=0;j<imax;j++) {
		for(int i=0;i<imax-j;i++) {
			if(st == EV3DStyle) out << ip++ << " 3 ";
			if(st == EMVStyle) out << ip++ << " 1 1 1 ";
			co0[0] = i; co0[1] = j;
			co1[0] = i+1; co1[1] = j;
			co2[0] = i; co2[1]=j+1;
			out << EqNum(co0) << " " << EqNum(co1) << " " <<
			EqNum(co2) << endl;
			
			if(i <imax-j-1) {
				if(st == EV3DStyle) out << ip++ << " 3 ";
				if(st == EMVStyle) out << ip++ << " 1 1 1 ";
				co0[0] = i+1; co0[1] = j;
				co1[0] = i+1; co1[1] = j+1;
				co2[0] = i; co2[1]=j+1;
				out << EqNum(co0) << " " << EqNum(co1) << " " <<
				EqNum(co2) << endl;
			}
			
		}
	}
}

int64_t TPZGraphElTd::EqNum(TPZVec<int> &co){
	int res = fGraphMesh->Res();
	int imax = (1<<res)+1;
	int neq = fConnect->FirstPoint() + co[1];
	for(int run=0;run<co[0];run++) neq += imax-run;
	return neq;
}

void TPZGraphElTd::QsiEta(TPZVec<int> &co, int imax, TPZVec<REAL> &qsieta){
	int i,ni = co.NElements();
	for(i=0; i<ni; i++) qsieta[i] = (1.*co[i])/imax;
}


void TPZGraphElTd::FirstIJ(int no, TPZVec<int> &co, int &incr){
	int i;
	for(i=0;i<3;i++) co[i]=0;
	incr = 1;
}


void TPZGraphElTd::NextIJ(int no,TPZVec<int> &co, int /*incr*/){
	int res = fGraphMesh->Res();
	int imax;
	imax = (1 << res)+1;
	co[1]++;
	if(co[0]+co[1]>= imax) {
		co[1] = 0;
		co[0]++;
	}
}

int TPZGraphElTd::ExportType(TPZDrawStyle st){
	switch(st)
	{
		case(EVTKStyle):
			return 5;//vtk_triangle
		default:
			return -1;
	}
}

int TPZGraphElTd::NNodes()
{
	return 3;
}
