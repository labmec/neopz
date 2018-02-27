/**
 * @file
 * @brief Contains the implementation of the TPZGraphElQ2dd methods. 
 */

#include "pzgraphelq2dd.h"
#include "pzgraphmesh.h"
#include "pzcompel.h"
#include "pzgeoel.h"

using namespace std;

TPZGraphElQ2dd::TPZGraphElQ2dd(TPZCompEl *cel, TPZGraphMesh *gmesh) : TPZGraphEl(cel,gmesh,fConnect)
{
}

TPZGraphNode *TPZGraphElQ2dd::Connect(int64_t i) {
	return fConnect;
}

TPZGraphElQ2dd::~TPZGraphElQ2dd(void)
{
}

int TPZGraphElQ2dd::NPoints(TPZGraphNode *n)
{
	int res = fGraphMesh->Res();
	int imax = (1<<res)+1;
	return imax*imax;
}

int TPZGraphElQ2dd::NElements(){
	int res = fGraphMesh->Res();
	int imax = (1<<res);
	return imax*imax;
}

int64_t TPZGraphElQ2dd::EqNum(TPZVec<int> &co) {
	int res = fGraphMesh->Res();
	int imax = (1<<res)+1;
	return fConnect->FirstPoint()+co[0]*imax+co[1];
}

void TPZGraphElQ2dd::FirstIJ(int connect,TPZVec<int> &co, int &incr) {	
	co[0]=0;
	co[1]=0;
	incr = 1;
}

void TPZGraphElQ2dd::NextIJ(int connect, TPZVec<int> &co, int /*incr*/) {
	
	int res = fGraphMesh->Res();
	int imax;
	imax = 1 << res;
	co[1]++;
	if(co[1]> imax) {
		co[1] = 0;
		co[0]++;
	}
}

void TPZGraphElQ2dd::Connectivity(TPZDrawStyle st){
	int res = fGraphMesh->Res();
	int imax = 1 << res;
	ostream &out = fGraphMesh->Out();
	int64_t ip = fId;
	if(st == EV3DStyle) ip++;
	TPZVec<int> co0(3,0), co1(3,0), co2(3,0), co3(3,0);
	for(int i=0;i<imax;i++) {
		for(int j=0;j<imax;j++) {
			if(st == EV3DStyle) out << ip << " 4 ";
			if(st == EVTKStyle) out << "4 ";
			if(st == EMVStyle) out << ip << " 1 1 1 ";
			ip++;
			if(st == EDXStyle) {
				co0[0] = i; co0[1] = j; 
				co1[0] = i+1; co1[1] = j; 
				co2[0] = i; co2[1] = j+1; 
				co3[0] = i+1; co3[1] = j+1;
			}
 			else {
				co0[0] = i; co0[1] = j; 
				co1[0] = i+1; co1[1] = j; 
				co2[0] = i+1; co2[1] = j+1; 
				co3[0] = i; co3[1] = j+1;
			}
			out << EqNum(co0) << " " << EqNum(co1) << " " <<
			EqNum(co2) << " " << EqNum(co3) << endl;
		}
	}
}

void TPZGraphElQ2dd::SetNode(int64_t i,TPZGraphNode *gno) {
	fConnect = gno;
}

int TPZGraphElQ2dd::ExportType(TPZDrawStyle st){
	switch(st)
	{
		case(EVTKStyle):
			return 9;//vtk_quad
		default:
			return -1;
	}
}

int TPZGraphElQ2dd::NNodes()
{
	return 4;
}
