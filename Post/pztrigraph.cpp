/**
 * @file
 * @brief Contains the implementation of the TPZGraphElT methods. 
 */

#include "pztrigraph.h"
#include "pzgraphmesh.h"

using namespace std;

TPZGraphElT::TPZGraphElT(TPZCompEl *c, TPZGraphMesh *g) : TPZGraphEl(c,g,fConnects) {
}

int TPZGraphElT::NConnects(){
	return 7;
}

int TPZGraphElT::NElements(){
	int res = fGraphMesh->Res();
	int imax = (1<<res);
	return imax*imax;
}

int TPZGraphElT::NPoints(TPZGraphNode *n){
	int res = fGraphMesh->Res();
	int imax = (1<<res);                                                                       
	int in = ConnectNum(n);
	switch(in) {
		case 0:
		case 1:
		case 2:
			return 1;
		case 3:
		case 4:
		case 5:
			return imax-1;
		case 6:
			return ((imax-1)*(imax-2))/2;
	}
	return 1;
}

void TPZGraphElT::Connectivity(TPZDrawStyle st){
	int res = fGraphMesh->Res();
	int imax = 1 << res;
	ostream &out = fGraphMesh->Out();
    int64_t ip = fId;
	TPZVec<int> co0(3,0), co1(3,0), co2(3,0);
	if(st == EV3DStyle) ip++;
	//   int64_t ip = (fId+1)*imax*imax;//Cedric 22/03/99
	for(int j=0;j<imax;j++) {
		for(int i=0;i<imax-j;i++) {
			if(st == EV3DStyle) out << ip++ << " 3 ";
			if(st == EMVStyle) out << ip++ << " 1 1 1 ";
			co0[0]=i; co0[1]=j;
			co1[0]=i+1; co1[1]=j;
			co2[0]=i; co2[1]=j+1;
			if(i <imax-j-1) {
				if(st == EV3DStyle) out << ip++ << " 3 ";
				if(st == EMVStyle) out << ip++ << " 1 1 1 ";
				co0[0]=i+1; co0[1]=j;
				co1[0]=i+1; co1[1]=j+1;
				co2[0]=i; co2[1]=j+1;
			}
			out << EqNum(co0) << " " << EqNum(co1) << " " <<
			EqNum(co2) << endl;
		}
	}
}

int64_t TPZGraphElT::EqNum(TPZVec<int> &co){
	int orient[3];
	int loc;
	int res = fGraphMesh->Res();
	int imax = (1<<res);
	for(int is=0;is<3;is++) {
		orient[is] = (fConnects[is]->SequenceNumber() > fConnects[(is+1)%3]->SequenceNumber()) ? 1 : 0;
	}
	if(co[0]==0 && co[1]==0) return fConnects[0]->FirstPoint();
	if(co[0]==imax && co[1]==0) return fConnects[1]->FirstPoint();
	if(co[0]==0 && co[1]==imax) return fConnects[2]->FirstPoint();
	
	if(co[1]==0) loc = 0;
	else if(co[0]+co[1]==imax) loc=1;
	else if(co[0]==0) loc=2;
	else loc=3;
	int64_t first = fConnects[3+loc]->FirstPoint();
	int64_t neq;
	switch(loc) {
		case 0:
			neq = first + (co[0]-1) + orient[loc]*(imax-1-co[0]-co[0]+1);
			break;
		case 1:
			neq = first + (imax-1-co[0]) + orient[loc]*(co[0]-1-imax+1+co[0]);
			break;
		case 2:
			neq = first + (imax-1-co[1]) + orient[loc]*(co[1]-1-imax+1+co[1]);
			break;
		case 3:
		{
			neq = first + (co[1]-1);
			for(int run=1;run<co[0];run++) neq += imax-run-1;
			break;
		}
		default:
			neq = 0;
			break;
	}
	return neq;
}

void TPZGraphElT::QsiEta(TPZVec<int> &co, int imax, TPZVec<REAL> &qsieta){
	int i,ni = co.NElements();
	for(i=0; i<ni; i++) qsieta[i] = (1.*co[i])/imax;
}

void TPZGraphElT::FirstIJ(int no,TPZVec<int> &co, int &incr){
	int res = fGraphMesh->Res();
	int imax;
	imax = 1 << res;
	switch(no) {
		case 0:
		{
			co[0]=0;
			co[1]=0;
			incr = 1;
		}
			break;
		case 1:
		{
			co[0] = imax;
			co[1]=0;
			incr = 1;
		}
			break;
		case 2:
		{
			co[0] = 0;
			co[1] = imax;
			incr = 1;
		}
			break;
		case 3:
		{
			co[1] = 0;
			if(fConnects[0]->SequenceNumber() > fConnects[1]->SequenceNumber())
			{
				co[0] = imax-1;
				incr = -1;
			}
			else
			{
				incr = 1;
				co[0] = 1;//ifirst
			}
		}
			break;
		case 4:
		{
			if(fConnects[1]->SequenceNumber() > fConnects[2]->SequenceNumber())
			{
				incr = 1;
				co[1] = imax-1;//co[1]first
				co[0] = 1;
			}
			else
			{
				incr = -1;
				co[0] = imax-1;
				co[1] = 1;//co[1]first
			}
		}
			break;
		case 5:
		{
			co[0] = 0;
			//modified Philippe 25/7/97
			// wrong comparaison
			//				if(fConnects[0]->Id() > fConnects[3]->Id())
			if(fConnects[0]->SequenceNumber() > fConnects[2]->SequenceNumber())
			{
				incr = -1;
				co[1] = imax-1;//co[1]first
			}
			else
			{
				incr = 1;
				co[1] = 1;//co[1]first
			}
		}
			break;
		case 6:
		{
			incr = 1;
			co[0] = 1;
			co[1] = 1;
		}
	}
	return;
}

void TPZGraphElT::NextIJ(int no, TPZVec<int> &co, int incr){
	int res = fGraphMesh->Res();
	int imax;
	imax = 1 << res;
	switch(no) {
		case 0:
		case 1:
		case 2:
			return;
		case 3:
			co[0] += incr;
			break;
		case 4:
			co[0] += incr;
			co[1] -= incr;
			break;
		case 5:
			co[1]+= incr;
			break;
		case 6:
			co[1]++;
			if(co[0]+co[1]>= imax) {
				co[1] = 1;
				co[0]++;
			}
			break;
	}
}

int TPZGraphElT::ExportType(TPZDrawStyle st){
	switch(st)
	{
		case(EVTKStyle):
			return 5;//vtk_triangle
		default:
			return -1;
	}
}

int TPZGraphElT::NNodes()
{
	return 3;
}
