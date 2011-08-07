/**
 * @file
 * @brief Contains the implementation of the TPZGraphEl methods. 
 */
#include "pzgraphel.h"
#include "pzgraphmesh.h"
#include "pzcompel.h"
#include "pzgeoel.h"
#include "pzmaterial.h"

using namespace std;

TPZGraphEl::TPZGraphEl(TPZCompEl *cel, TPZGraphMesh *gmesh, TPZGraphNode **connectvec)
{
	fCompEl = cel;
	fGraphMesh = gmesh;
	fId = cel->Index();
	TPZGraphNode *gno;
	for(int j = 0; j < cel->Reference()->NSides(); j++) {
		TPZConnect &cn = cel->Connect(j);
		int newsize = cn.SequenceNumber()+1;
		if( gmesh->NodeMap().NElements() < newsize ) {
            gmesh->NodeMap().Resize(newsize);
		}
		gno = &gmesh->NodeMap()[cn.SequenceNumber()];
		if(gno->SequenceNumber() ==-1) {
			gno->SetElement(this);
			gno->SetSequenceNumber(cn.SequenceNumber());
			gno->SetConnect(&cn);
			gno->SetGraphMesh(gmesh);
		}
		connectvec[j] = gno;
	}
	int index = gmesh->ElementList().AllocateNewElement();
	gmesh->ElementList()[index] = this;
}

TPZGraphEl::TPZGraphEl(TPZCompEl *cel, TPZGraphMesh *gmesh, TPZGraphNode *&connect) {
	fCompEl = cel;
	fGraphMesh = gmesh;
	fId = cel->Index();
	int index = gmesh->NodeMap().AllocateNewElement();
	TPZGraphNode *gno = &gmesh->NodeMap()[index];
	gno->SetElement(this);
	gno->SetSequenceNumber(index);
	gno->SetConnect(0);
	gno->SetGraphMesh(gmesh);
	connect = gno;
	
	index = gmesh->ElementList().AllocateNewElement();
   	gmesh->ElementList()[index] = this;
}

void TPZGraphEl::SetNode(int, TPZGraphNode *) {
}

TPZGraphEl::~TPZGraphEl(void)
{
}

void TPZGraphEl::QsiEta(TPZVec<int> &i, int imax, TPZVec<REAL> &qsieta)
{
	int ind,nel=i.NElements();
	for(ind=0; ind<nel; ind++)
	{
		qsieta[ind] = (-1.0+(i[ind]*2.0/imax));
	}
}


int TPZGraphEl::ConnectNum(TPZGraphNode *n) {
	int nn = NConnects();
	for(int i=0;i<nn;i++) if(n==Connect(i)) return i;
	return 0;
}


void TPZGraphEl::DrawCo(TPZGraphNode *n, TPZDrawStyle st)
{
	int in = ConnectNum(n);
	//int i,j,incr;
	int incr;
	TPZVec<int> co(3,0);
	FirstIJ(in,co,incr);
	//	ComputeSequence(n, ibound, incr);
	long ip = n->FirstPoint();
	int res = fGraphMesh->Res();
	int imax;
	imax = 1 << res;
	int np = NPoints(n);
	int point=0;
	TPZVec<REAL> qsi(3,0.),x(4,0.);
	while(point < np) {
		QsiEta(co,imax,qsi);
		fCompEl->Reference()->X(qsi,x);
		if(st == EMVStyle || st == EV3DStyle) fGraphMesh->Out() << ip++ << " ";
		fGraphMesh->Out() << x[0] << " " << x[1] << " " << x[2] << endl;
		NextIJ(in,co,incr);
		point++;
	}
}

void TPZGraphEl::DrawSolution(TPZGraphNode *n,int solind,TPZDrawStyle st) {
	TPZManVector<int> sol(1,solind);
	DrawSolution(n,sol,st);
}

void TPZGraphEl::DrawSolution(TPZGraphNode * /*n*/,TPZBlock &,TPZDrawStyle /*st*/) {
}

void TPZGraphEl::DrawSolution(TPZGraphNode *n,TPZVec<int> &solind,TPZDrawStyle st)
{
	int in = ConnectNum(n);
	//int i,j,incr;
	int incr;
	TPZVec<int> co(3,0);
	FirstIJ(in,co,incr);
	//	ComputeSequence(n, ibound, incr);
	int res = fGraphMesh->Res();
	int imax;
	int numsol = solind.NElements();
	int numvar;
	imax = 1 << res;
	int np = NPoints(n);
	int point=0;
	TPZManVector<REAL> qsi(3,0.),x(4,0.);
	TPZManVector<REAL> sol(6,0.);
	long ip = n->FirstPoint();
	while(point < np) 
	{
		QsiEta(co,imax,qsi);
		if(st == EMVStyle || st == EV3DStyle) fGraphMesh->Out() << ip++ << " ";
		for(int is=0; is<numsol; is++) 
		{
			fCompEl->Solution(qsi,solind[is],sol);
			numvar = fCompEl->Material()->NSolutionVariables(solind[is]);
			if(st == EVTKStyle)
			{
				if(numvar > 3) numvar = 3;
			}
			int iv;
			for(iv=0; iv<numvar;iv++)
			{
				if(fabs(sol[iv]) < 1.0e-20) sol[iv] = 0.0;
				fGraphMesh->Out() << sol[iv] << " ";
			}
			if((st == EMVStyle || st == EV3DStyle) && numvar ==2) fGraphMesh->Out() << 0. << " ";
			if(st == EVTKStyle && numvar != 1)
			{
				for(; iv<3; iv++) fGraphMesh->Out() << 0. << " ";
			}
		}
		
		fGraphMesh->Out() << endl;
		NextIJ(in,co,incr);
		point++;
	}
}

void TPZGraphEl::Print(ostream &out) {
	out << "TPZGraphEl element id = " << fId << endl;
	out << "Node numbers : ";
	int i;
	for(i=0; i<NConnects(); i++) {
		out << Connect(i)->SequenceNumber() << " ";
	}
	out << endl;
	for(i=0; i<NConnects(); i++) {
		out << Connect(i)->FirstPoint() << " ";
	}
	out << endl;
}

