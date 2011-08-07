/**
 * @file
 * @brief Contains the implementation of the TPZGraphElQ2d methods. 
 */
#include "pzgraphelq2d.h"
#include "pzgraphmesh.h"
#include "pzcompel.h"
#include "pzgeoel.h"

using namespace std;

TPZGraphElQ2d::TPZGraphElQ2d(TPZCompEl *cel, TPZGraphMesh *gmesh) : TPZGraphEl(cel,gmesh,fConnects)
{
}

TPZGraphNode *TPZGraphElQ2d::Connect(int i) {
	return fConnects[i];
}

TPZGraphElQ2d::~TPZGraphElQ2d(void)
{
}


int TPZGraphElQ2d::NPoints(TPZGraphNode *n)
{
	int in = ConnectNum(n);
	int res = fGraphMesh->Res();
	switch(in)
	{
		case 0:
		case 1:
		case 2:
		case 3:
			return(1);
		case 4:
		case 5:
		case 6:
		case 7:
			return((1 << res)-1);
		case 8:
		{
			int a = ((1 << res)-1);
			return(a*a);
		}
	};
	return(0);
}

int TPZGraphElQ2d::NElements(){
	int res = fGraphMesh->Res();
	int imax = (1<<res);
	return imax*imax;
}

long TPZGraphElQ2d::EqNum(TPZVec<int> &co) {
	int orient[4];
	int loc;
	int res = fGraphMesh->Res();
	int imax = (1<<res);
	for(int is=0;is<4;is++) {
		orient[is] = (fConnects[is]->SequenceNumber() > fConnects[(is+1)%4]->SequenceNumber()) ? 1 : 0;
	}
	if(co[0]==0 && co[1]==0) return fConnects[0]->FirstPoint();
	if(co[0]==imax && co[1]==0) return fConnects[1]->FirstPoint();
	if(co[0]==imax && co[1]==imax) return fConnects[2]->FirstPoint();
	if(co[0]==0 && co[1]==imax) return fConnects[3]->FirstPoint();
	
	if(co[1]==0) loc = 0;
	else if(co[0]==imax) loc=1;
	else if(co[1]==imax) loc=2;
	else if(co[0]==0) loc=3;
	else loc=4;
	long first = fConnects[4+loc]->FirstPoint();
	long neq;
	switch(loc) {
		case 0:
			neq = first + (co[0]-1) + orient[loc]*(imax-1-co[0]-co[0]+1);
			break;
		case 1:
			neq = first + (co[1]-1) + orient[loc]*(imax-1-co[1]-co[1]+1);
			break;
		case 2:
			neq = first + (imax-1-co[0]) + orient[loc]*(co[0]-1-imax+1+co[0]);
			break;
		case 3:
			neq = first + (imax-1-co[1]) + orient[loc]*(co[1]-1-imax+1+co[1]);
			break;
		case 4:
			neq = first + (co[1]-1) + (co[0]-1)*(imax-1);
			break;
		default:
			neq = 0;
			break;
	}
	return neq;
}

void TPZGraphElQ2d::FirstIJ(int connect,TPZVec<int> &co, int &incr) {
	
	int res = fGraphMesh->Res();
	int imax, i;
	imax = 1 << res;
	switch(connect)
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
			for (i=1;i<3;i++) co[i]=0;
			incr = 1;
		}
			break;
		case 2:
		{
			for (i=0;i<2;i++) co[i]= imax;
			incr = 1;
		}
			break;
		case 3:
		{
			co[0] = 0;
			co[1] = imax;
			incr = 1;
		}
			break;
		case 4:
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
		case 5:
		{
			co[0] = imax;
			if(fConnects[1]->SequenceNumber() > fConnects[2]->SequenceNumber())
			{
				incr = -1;
				co[1] = imax-1;//jfirst
			}
			else
			{
				incr = 1;
				co[1] = 1;//jfirst
			}
		}
			break;
		case 6:
		{
			co[1] = imax;
			if(fConnects[3]->SequenceNumber() > fConnects[2]->SequenceNumber())
			{
				incr = -1;
				co[0] = imax-1;//ifirst
			}
			else
			{
				incr = 1;
				co[0] = 1;//ifirst
			}
		}
			break;
		case 7:
		{
			co[0] = 0;
			if(fConnects[0]->SequenceNumber() > fConnects[3]->SequenceNumber())
			{
				incr = -1;
				co[1] = imax-1;//jfirst
			}
			else
			{
				incr = 1;
				co[1] = 1;//jfirst
			}
		}
			break;
		case 8:
		{
			incr = 1;
			for (i=0;i<2;i++) co[i]= 1;
		}
	}
	return;
}

void TPZGraphElQ2d::NextIJ(int connect, TPZVec<int> &co, int incr) {
	
	int res = fGraphMesh->Res();
	int imax;
	imax = 1 << res;
	switch(connect) {
		case 0:
		case 1:
		case 2:
		case 3:
			return;
		case 4:
		case 6:
			co[0] += incr;
			break;
		case 5:
		case 7:
			co[1]+= incr;
			break;
		case 8:
			co[1]++;
			if(co[1]>= imax) {
				co[1] = 1;
				co[0]++;
			}
			break;
	}
}
/*
 void TPZGraphElQ2d::ComputeSequence(TPZGraphNode *n, int *ibound, int *incr)
 {
 int res = fGraphMesh->Res();
 int imax, jmax;
 imax = jmax = 1 << res;
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
 case 2:
 {
 ibound[0] = imax;//ifirst
 ibound[1] = imax;//ilast
 ibound[2] = jmax;//jfirst
 ibound[3] = jmax;//jlast
 }
 break;
 case 3:
 {
 ibound[0] = 0;//ifirst
 ibound[1] = 0;//ilast
 ibound[2] = jmax;//jfirst
 ibound[3] = jmax;//jlast
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
 case 5:
 {
 if(fConnects[1]->SequenceNumber() > fConnects[2]->SequenceNumber())
 {
 incr[1] = -1;
 ibound[2] = jmax-1;//jfirst
 ibound[3] = 1;//jlast
 }
 else
 {
 incr[1] = 1;
 ibound[2] = 1;//jfirst
 ibound[3] = jmax-1;//jlast
 }
 ibound[0] = imax;//ifirst
 ibound[1] = imax;//ilast
 }
 break;
 case 6:
 {
 if(fConnects[3]->SequenceNumber() > fConnects[2]->SequenceNumber())
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
 ibound[2] = jmax;//jfirst
 ibound[3] = jmax;//jlast
 }
 break;
 case 7:
 {
 if(fConnects[0]->SequenceNumber() > fConnects[3]->SequenceNumber())
 {
 incr[1] = -1;
 ibound[2] = jmax-1;//jfirst
 ibound[3] = 1;//jlast
 }
 else
 {
 incr[1] = 1;
 ibound[2] = 1;//jfirst
 ibound[3] = jmax-1;//jlast
 }
 ibound[0] = 0;//ifirst
 ibound[1] = 0;//ilast
 }
 break;
 case 8:
 {
 incr[0] = 1;
 ibound[0] =1;
 ibound[1] =imax-1;
 incr[1] = 1;
 ibound[2] = 1;
 ibound[3] = imax-1;
 //
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
 
 if(fConnects[1]->SequenceNumber() > fConnects[2]->SequenceNumber())
 {
 incr[1] = -1;
 ibound[2] = jmax-1;//jfirst
 ibound[3] = 1;//jlast
 }
 else
 {
 incr[1] = 1;
 ibound[2] = 1;//jfirst
 ibound[3] = jmax-1;//jlast
 }
 //
 }
 };
 return;
 }
 */

void TPZGraphElQ2d::Connectivity(TPZDrawStyle st){
	int res = fGraphMesh->Res();
	int imax = 1 << res;
	ostream &out = fGraphMesh->Out();
	long ip = fId;
	TPZVec<int> co0(3,0),co1(3,0),co2(3,0),co3(3,0);
	if(st == EV3DStyle) ip++;
	for(int i=0;i<imax;i++) {
		for(int j=0;j<imax;j++) {
			if(st == EV3DStyle) out << ip << " 4 ";
			if(st == EMVStyle) out << ip << " 1 1 1 ";
			ip++;
			if(st == EDXStyle) {
				co0[0]=i; co0[1]=j;
				co1[0]=i+1; co1[1]=j;
				co2[0]=i; co2[1]=j+1;
				co3[0]=i+1; co3[1]=j+1;
			} else {
				co0[0]=i; co0[1]=j;
				co1[0]=i+1; co1[1]=j;
				co2[0]=i+1; co2[1]=j+1;
				co3[0]=i; co3[1]=j+1;
			}
			out << EqNum(co0) << " " << EqNum(co1) << " " <<
			EqNum(co2) << " " << EqNum(co3) << endl;
		}
	}
}

void TPZGraphElQ2d::SetNode(int i,TPZGraphNode *gno) {
	fConnects[i] = gno;
}


int TPZGraphElQ2d::ExportType(TPZDrawStyle st){
	switch(st)
	{
		case(EVTKStyle):
			return 9;//vtk quad
			//		break;
		default:
			return -1;
	}
	//	return -1;
}

int TPZGraphElQ2d::NNodes()
{
	return 4;
}
