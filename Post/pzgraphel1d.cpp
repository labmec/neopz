/**
 * @file
 * @brief Contains the implementation of the TPZGraphEl1d methods. 
 */
#include "pzgraphel1d.h"
#include "pzgraphmesh.h"
#include "pzcompel.h"
#include "pzgeoel.h"

using namespace std;

TPZGraphEl1d::TPZGraphEl1d(TPZCompEl *ce, TPZGraphMesh *gg) : TPZGraphEl(ce,gg,fConnects)
{
}

int TPZGraphEl1d::NPoints(TPZGraphNode *n)
{
	int in = ConnectNum(n);
	int res = fGraphMesh->Res();
	switch(in)
	{
		case 0:
		case 1:
			return(1);
		case 2:
			return((1 << res)-1);
	};
	return(0);
}

int TPZGraphEl1d::NElements(){
	int res = fGraphMesh->Res();
	int imax = (1<<res);
	return imax;
}

long TPZGraphEl1d::EqNum(TPZVec<int> &co) {
	int orient;
	int res = fGraphMesh->Res();
	int imax = (1<<res);
	//orient = (fConnects[0]->Id() > fConnects[1]->Id()) ? 1 : 0;
	orient = (fConnects[0]->SequenceNumber() > fConnects[1]->SequenceNumber()) ? 1 : 0;
	if(co[0]==0) return fConnects[0]->FirstPoint();
	if(co[0]==imax) return fConnects[1]->FirstPoint();
	
	long first = fConnects[2]->FirstPoint();
	long neq;
	neq = first + (co[0]-1) + orient*(imax-1-co[0]-co[0]+1);
	return neq;
}

void TPZGraphEl1d::FirstIJ(int no,TPZVec<int> &co,int &incr) {
	
	int res = fGraphMesh->Res();
	int imax, i;
	imax = 1 << res;
	switch(no)
	{
		case 0:
		{
			for (i=0;i<3;i++) co[i]=0;
			incr = 1;
		}
			break;
		case 1:
		{
			co[0] = imax;
			for (i=1;i<3; i++) co[i]=0;
			incr = 1;
		}
			break;
		case 2:
		{
			for (i=1;i<3; i++) co[i]=0;
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
	}
	return;
}

void TPZGraphEl1d::NextIJ(int no,TPZVec<int> &co, int incr) {
	
	switch(no) {
		case 0:
		case 1:
			return;
		case 2:
			co[0] += incr;
			break;
	}
}

/*
 void TPZGraphEl1d::ComputeSequence(TPZGraphNode *n, int *ibound, int *incr)
 {
 int res = fGraphMesh->Res();
 int imax;
 imax = 1 << res;
 int in = ConnectNum(n);
 incr[0] = incr[1] = 1;
 switch(in)
 {
 case 0:
 {
 ibound[0] = 0;//ifirst
 ibound[1] = 0;//ilast
 ibound[2] = 0;//jfirst
 ibound[3] = 0;//jlast
 }
 break;
 case 1:
 {
 ibound[0] = imax;//ifirst
 ibound[1] = imax;//ilast
 ibound[2] = 0;//jfirst
 ibound[3] = 0;//jlast
 }
 break;
 case 4:
 {
 if(fConnects[0]->SequenceNumber() > fConnects[1]->SequenceNumber())
 {
 incr[0] = -1;
 ibound[0] = imax-1;//ifirst
 ibound[1] = 1;//ilast
 }
 else
 {
 incr[0] = 1;
 ibound[0] = 1;//ifirst
 ibound[1] = imax-1;//ilast
 }
 ibound[2] = 0;//jfirst
 ibound[3] = 0;//jlast
 }
 break;
 }
 return;
 }
 */

void TPZGraphEl1d::Connectivity(TPZDrawStyle st){
	int res = fGraphMesh->Res();
	int imax = 1 << res;
	ostream &out = fGraphMesh->Out();
	long ip = fId;
	TPZVec<int> co0(3,0), co1(3,0);
	if(st == EV3DStyle) ip++;
	for(int i=0;i<imax;i++) {
		if(st == EV3DStyle) out << ip << " 1 ";
		if(st == EMVStyle) out << ip << " 1 1 1 ";
		ip++;
		co0[0]=i;
		co1[0]=i+1;
		out << EqNum(co0) << " " << EqNum(co1) << endl;
	}
}

void TPZGraphEl1d::Print(ostream &out) {
	out << "TPZGraphEl1d element id = " << fId << endl;
	out << "Node numbers : ";
	int i;
	for(i=0; i<NConnects(); i++) {
		out << fConnects[i]->SequenceNumber() << " ";
	}
	out << endl << "First Equation : ";
	for(i=0; i<NConnects(); i++) {
		out << fConnects[i]->FirstPoint() << " ";
	}
	out << endl;
}

int TPZGraphEl1d::ExportType(TPZDrawStyle st){
	switch(st)
	{
		case(EVTKStyle):
			return 1;
			//		break;
		default:
			return -1;
	}
	//	return -1;
}

int TPZGraphEl1d::NNodes()
{
	return 2;	
}


