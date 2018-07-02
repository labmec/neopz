/**
 * @file
 * @brief Contains the implementation of the TPZGraphEl1dd methods. 
 */

#include "pzgraphel1dd.h"
#include "pzgraphmesh.h"
#include "pzcompel.h"
#include "pzgeoel.h"

using namespace std;

TPZGraphEl1dd::TPZGraphEl1dd(TPZCompEl *ce, TPZGraphMesh *gg) : TPZGraphEl(ce,gg,fConnect)
{
}

int TPZGraphEl1dd::NPoints(TPZGraphNode *n)
{
	int res = fGraphMesh->Res();
	return((1 << res)+1);
}

int TPZGraphEl1dd::NElements(){
	int res = fGraphMesh->Res();
	int imax = (1<<res);
	return imax;
}

int64_t TPZGraphEl1dd::EqNum(TPZVec<int> &co) {
	return fConnect->FirstPoint() + co[0];
}

void TPZGraphEl1dd::FirstIJ(int no, TPZVec<int> &co, int &incr) {
    co[0]=0; // @omar:: I think the next line does not make sense
//    for (i=0;i<3;i++) co[i]=0;
	incr = 1;
}

void TPZGraphEl1dd::NextIJ(int no,TPZVec<int> &co, int incr) {
	
	co[0]+=incr;
}

void TPZGraphEl1dd::Connectivity(TPZDrawStyle st){
	int res = fGraphMesh->Res();
	int imax = 1 << res;
	ostream &out = fGraphMesh->Out();
	int64_t ip = fId;
	TPZVec<int> co0(3,0), co1(3,0);
	
	if(st == EV3DStyle) ip++;
	for(int i=0;i<imax;i++) {
		if(st == EV3DStyle) out << ip << " 1 ";
		if(st == EVTKStyle) out << "2 ";		
		if(st == EMVStyle) out << ip << " 1 1 1 ";
		ip++;
		co0[0]= i;
		co1[0]=i+1;
		out << EqNum(co0) << " " << EqNum(co1) << endl;
	}
}

void TPZGraphEl1dd::Print(ostream &out) {
	out << "TPZGraphEl1dd element id = " << fId << endl;
	out << "Node numbers : ";
	out << fConnect->SequenceNumber() << " ";
	out << endl << "First Equation : ";
	out << fConnect->FirstPoint() << " ";
	out << endl;
}

int TPZGraphEl1dd::ExportType(TPZDrawStyle st){
	switch(st)
	{
		case(EVTKStyle):
			return 3;
		default:
			return -1;
	}
}

int TPZGraphEl1dd::NNodes()
{
	return 2;
}
