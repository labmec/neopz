/**
 * @file
 * @brief Contains the implementation of the TPZGraphEl methods. 
 */

#include "pzgraphel.h"
#include "pzgraphmesh.h"
#include "pzcompel.h"
#include "pzgeoel.h"
#include "pzcmesh.h"
#include "TPZMaterial.h"

using namespace std;

TPZGraphEl::TPZGraphEl(TPZCompEl *cel, TPZGraphMesh *gmesh, TPZGraphNode **connectvec)
{
	fCompEl = cel;
	fGraphMesh = gmesh;
	fId = cel->Index();
	TPZGraphNode *gno;
	for(int j = 0; j < cel->Reference()->NSides(); j++) {
		TPZConnect &cn = cel->Connect(j);
		int64_t newsize = cn.SequenceNumber()+1;
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
	int64_t index = gmesh->ElementList().AllocateNewElement();
	gmesh->ElementList()[index] = this;
}

TPZGraphEl::TPZGraphEl(TPZCompEl *cel, TPZGraphMesh *gmesh, TPZGraphNode *&connect) {
	fCompEl = cel;
	fGraphMesh = gmesh;
	fId = cel->Index();
	int64_t index = gmesh->NodeMap().AllocateNewElement();
	TPZGraphNode *gno = &gmesh->NodeMap()[index];
	gno->SetElement(this);
	gno->SetSequenceNumber(index);
	gno->SetConnect(0);
	gno->SetGraphMesh(gmesh);
	connect = gno;
	
	index = gmesh->ElementList().AllocateNewElement();
   	gmesh->ElementList()[index] = this;
}

void TPZGraphEl::SetNode(int64_t, TPZGraphNode *) {
}

TPZGraphEl::~TPZGraphEl(void)
{
}

void TPZGraphEl::QsiEta(TPZVec<int> &i, int imax, TPZVec<REAL> &qsieta)
{
	int64_t ind,nel=i.NElements();
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
    int dim = this->Dimension();
	TPZVec<int> co(dim,0);
	FirstIJ(in,co,incr);
	//	ComputeSequence(n, ibound, incr);
	int64_t ip = n->FirstPoint();
	int res = fGraphMesh->Res();
	int imax;
	imax = 1 << res;
	int np = NPoints(n);
	int point=0;
	TPZManVector<REAL,3> qsi(dim,0.),x(3,0.);
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
	if(fCompEl->Mesh()->GetSolType() == EReal){
		return DrawSolutionT<STATE>(n,solind,st);
	}else{
		return DrawSolutionT<CSTATE>(n,solind,st);
	}
}
template<class TVar>
void TPZGraphEl::DrawSolutionT(TPZGraphNode *n,TPZVec<int> &solind,TPZDrawStyle st)
{
	int in = ConnectNum(n);
	//int i,j,incr;
	int incr;
    int dim = Dimension();
	TPZManVector<int,3> co(dim,0);
	FirstIJ(in,co,incr);
	//	ComputeSequence(n, ibound, incr);
	int res = fGraphMesh->Res();
	int imax;
	int64_t numsol = solind.NElements();
	int numvar;
	imax = 1 << res;
	int np = NPoints(n);
	int point=0;
	TPZManVector<REAL,4> qsi(dim,0.),x(4,0.);
	TPZManVector<TVar,10> sol(6,0.);
	int64_t ip = n->FirstPoint();
	while(point < np) 
	{
		QsiEta(co,imax,qsi);
		if(st == EMVStyle || st == EV3DStyle) fGraphMesh->Out() << ip++ << " ";
		for(int64_t is=0; is<numsol; is++) 
		{
			fCompEl->Solution(qsi,solind[is],sol);
			numvar = fCompEl->Material()->NSolutionVariables(solind[is]);
			if(st == EVTKStyle)
			{
				if(numvar > 9) numvar = 9; // Because it 3x3 tensor variables
			}
			int iv;
			for(iv=0; iv<numvar;iv++)
			{
				auto solprint = [&sol,iv](){
					if constexpr (std::is_same_v<TVar, RTVar>) {
						return sol[iv];
					} else {
						return std::real(sol[iv]);
					}
				}();
				
				if (fabs(solprint) < 1.0e-20)
					{solprint = 0.0;}
				fGraphMesh->Out() << solprint << " ";
			}
			if ((st == EMVStyle || st == EV3DStyle) && numvar == 2)
				{fGraphMesh->Out() << 0. << " ";}
			if(st == EVTKStyle && numvar != 1){
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

int TPZGraphEl::StaticClassId() {
    return Hash("TPZGraphEl");
}

int TPZGraphEl::ClassId() const {
    return StaticClassId();
}

void TPZGraphEl::Read(TPZStream &buf, void *context) {
    fCompEl = dynamic_cast<TPZCompEl *>(TPZPersistenceManager::GetInstance(&buf));
    fGraphMesh = dynamic_cast<TPZGraphMesh *>(TPZPersistenceManager::GetInstance(&buf));
    buf.Read(&fId);
}

void TPZGraphEl::Write(TPZStream &buf, int withclassid) const {
    TPZPersistenceManager::WritePointer(fCompEl, &buf);
    TPZPersistenceManager::WritePointer(fGraphMesh, &buf);
    buf.Write(&fId);
}
