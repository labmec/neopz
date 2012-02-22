/**
 * @file
 * @brief Contains the implementation of the TPZGraphElQ3dd methods. 
 */

#include "pzgraphelq3dd.h"
using namespace std;


// Construction/Destruction

TPZGraphElQ3dd::TPZGraphElQ3dd(TPZCompEl *cel, TPZGraphMesh *gmesh) : TPZGraphEl(cel,gmesh,fConnect){
}

TPZGraphNode *TPZGraphElQ3dd::Connect(int i) {
	return fConnect;
}

TPZGraphElQ3dd::~TPZGraphElQ3dd(void)
{
	
}

int TPZGraphElQ3dd::NPoints(TPZGraphNode *n)
{
	int res = fGraphMesh->Res();
	int imax = (1<<res)+1;
	return imax*imax*imax;
}

int TPZGraphElQ3dd::NElements(){
	int res = fGraphMesh->Res();
	int imax = (1<<res);
	return imax*imax*imax;
}

long TPZGraphElQ3dd::EqNum(TPZVec<int> &co) {
	int res = fGraphMesh->Res();
	int imax = (1<<res)+1;
	
	return fConnect->FirstPoint()+imax*imax*co[2]+co[1]*imax+co[0];
}

void TPZGraphElQ3dd::FirstIJ(int connect, TPZVec<int> &co, int &incr) {
	int i;
	for (i=0; i<3; i++)	co[i]=0;
	incr = 1;
}

void TPZGraphElQ3dd::NextIJ(int connect,  TPZVec<int> &co, int /*incr*/) {
	
	int res = fGraphMesh->Res();
	int imax;
	imax = 1 << res;
	co[0]++;
	if(co[0]> imax) {
		co[0] = 0;
		co[1]++;
		if(co[1] > imax) {
			co[2]++;
			co[1] = 0;
		}
	}
}

void TPZGraphElQ3dd::Connectivity(TPZDrawStyle st){
	int res = fGraphMesh->Res();
	int imax = 1 << res;
	ostream &out = fGraphMesh->Out();
	long ip = fId;
	TPZVec<int> co0(3,0), co1(3,0), co2(3,0), co3(3,0), 
	co4(3,0),co5(3,0),co6(3,0),co7(3,0);
	if(st == EV3DStyle) ip++;
	for(int i=0;i<imax;i++) {
		for(int j=0;j<imax;j++) {
			for (int k=0;k<imax;k++){
				if(st == EV3DStyle) out << ip << " 8 ";
				if(st == EVTKStyle) out << "8 ";
				if(st == EMVStyle) out << ip << " 1 1 1 ";
				ip++;
				if(st == EDXStyle) {
					co0[0] = i;   co0[1] = j;   co0[2]= k; 
					co1[0] = i+1; co1[1] = j;   co1[2]= k;
					co2[0] = i;   co2[1] = j+1; co2[2]= k;
					co3[0] = i+1; co3[1] = j+1; co3[2]= k;
					co4[0] = i;   co4[1] = j;   co4[2]= k+1; 
					co5[0] = i+1; co5[1] = j;   co5[2]= k+1;
					co6[0] = i;   co6[1] = j+1; co6[2]= k+1;
					co7[0] = i+1; co7[1] = j+1; co7[2]= k+1;
				}
 				else {
					co0[0] = i;   co0[1] = j;   co0[2]= k; 
					co1[0] = i+1; co1[1] = j;   co1[2]= k;
					co2[0] = i+1; co2[1] = j+1; co2[2]= k;
					co3[0] = i;   co3[1] = j+1; co3[2]= k;
					co4[0] = i;   co4[1] = j;   co4[2]= k+1; 
					co5[0] = i+1; co5[1] = j;   co5[2]= k+1;
					co6[0] = i+1; co6[1] = j+1; co6[2]= k+1;
					co7[0] = i;   co7[1] = j+1; co7[2]= k+1;
				}
				out << EqNum(co0) << " " 
				<< EqNum(co1) << " " 
				<< EqNum(co2) << " " 
				<< EqNum(co3) << " " 
				<< EqNum(co4) << " " 
				<< EqNum(co5) << " " 
				<< EqNum(co6) << " " 
				<< EqNum(co7) << endl;
			}
		}
	}
}

void TPZGraphElQ3dd::SetNode(int i,TPZGraphNode *gno) {
	fConnect = gno;
}

int TPZGraphElQ3dd::ExportType(TPZDrawStyle st){
	switch(st)
	{
		case(EVTKStyle):
			return 12;//vtk_hexahedron
			//		break;
		default:
			return -1;
	}
	//	return -1;
}

int TPZGraphElQ3dd::NNodes()
{
	return 8;	
}
