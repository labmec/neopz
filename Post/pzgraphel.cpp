
/**
 * @file
 * @brief Contains the implementation of the TPZGraphEl methods. 
 */

#include "pzgraphel.h"
#include "pzgraphmesh.h"
#include "pzcompel.h"
#include "pzgeoel.h"
#include "pzmaterial.h"
#include "TPZMatElasticity2D.h"

using namespace std;

TPZGraphEl::TPZGraphEl(TPZCompEl *cel, TPZGraphMesh *gmesh, TPZGraphNode **connectvec)
{
	fCompEl = cel;
	fGraphMesh = gmesh;
	fId = cel->Index();
	TPZGraphNode *gno;
	for(int j = 0; j < cel->Reference()->NSides(); j++) {
		TPZConnect &cn = cel->Connect(j);
		long newsize = cn.SequenceNumber()+1;
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
	long index = gmesh->ElementList().AllocateNewElement();
	gmesh->ElementList()[index] = this;
}

TPZGraphEl::TPZGraphEl(TPZCompEl *cel, TPZGraphMesh *gmesh, TPZGraphNode *&connect) {
	fCompEl = cel;
	fGraphMesh = gmesh;
	fId = cel->Index();
	long index = gmesh->NodeMap().AllocateNewElement();
	TPZGraphNode *gno = &gmesh->NodeMap()[index];
	gno->SetElement(this);
	gno->SetSequenceNumber(index);
	gno->SetConnect(0);
	gno->SetGraphMesh(gmesh);
	connect = gno;
	
	index = gmesh->ElementList().AllocateNewElement();
   	gmesh->ElementList()[index] = this;
}

void TPZGraphEl::SetNode(long, TPZGraphNode *) {
}

TPZGraphEl::~TPZGraphEl(void)
{
}

void TPZGraphEl::QsiEta(TPZVec<int> &i, int imax, TPZVec<REAL> &qsieta)
{
	long ind,nel=i.NElements();
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
        
        TPZMatElasticity2D projectmaterial;
        REAL Pi = M_PI;
        
        //******* angulos colocados a mao para fazer teste *********
        REAL alpha = (60.*(Pi/180)); // azimuth
        REAL beta = (30.*(Pi/180)); // inclination
        
        int projection = 0;
        //obtem angulos da rotacao do poco (como??)
        projectmaterial.GetWellboreAngles(alpha, beta, projection);
        
//        //Print Rotated Coordinate
//        std::cout << "Coord: " << endl;
//        std::cout << "x: " << x[0] << " " << "y: " << x[1] << " " << "z: " << x[2] << endl;
        
        //************ O codigo nao esta pegando os dados do material ja inseridos no main, pois estou criando um novo construtor, como faco????  *****************//
        

        if (projection==1) {
            
            
            // Por este metodo, eh preciso rotacionar a malha geometrica para obter projecao.
//            // cria vetor normal rotacionada e coordenada projetada
//            TPZVec<REAL> nRot(3,0.),xP(3,0.);
//            
//            nRot[0] = sin(beta);
//            nRot[1] = 0;
//            nRot[2] = cos(beta);
//            
//            REAL gamma = 0.;
//            gamma = x[2]/cos(beta);
//            
//            xP[0] = x[0] - gamma*nRot[0];
//            xP[1] = x[1] - gamma*nRot[1];
//            xP[2] = x[2] - gamma*nRot[2];
//            
//            x[0] = xP[0];
//            x[1] = xP[1];
//            x[2] = xP[2];
            
            
            
//            //********* Essa transformacao assume a mesma matriz de rotacao do poco inclinado, ou seja, rotacoes no sentido horario em z e depois em y  **********//
            REAL xP = 0., yP = 0., zP = 0.;
            
            xP =  x[0]*cos(alpha)*cos(beta) + x[1]*cos(beta)*sin(alpha) - x[2]*sin(beta) + (x[2]*cos(beta) + x[0]*cos(alpha)*sin(beta) + x[1]*sin(alpha)*sin(beta))*tan(beta);
            //x*Cos(\[Alpha])*Cos(\[Beta]) + y*Cos(\[Beta])*Sin(\[Alpha]) - z*Sin(\[Beta]) + (z*Cos(\[Beta]) + x*Cos(\[Alpha])*Sin(\[Beta]) + y*Sin(\[Alpha])*Sin(\[Beta]))*Tan(\[Beta])
            
            yP = x[1]*cos(alpha) - x[0]*sin(alpha);
            //y*Cos(\[Alpha]) - x*Sin(\[Alpha])
            zP = 0;
            
            x[0] = xP;
            x[1] = yP;
            x[2] = zP;
            
            
            fGraphMesh->Out() << x[0] << " " << x[1] << " " << x[2] << endl;
        }
        
        else{
            
		fGraphMesh->Out() << x[0] << " " << x[1] << " " << x[2] << endl;
        }
        
		NextIJ(in,co,incr);
        
//        //Print Projected Coordinate
//        std::cout << "Coord Projetada: " << endl;
//        std::cout << "x: " << x[0] << " " << "y: " << x[1] << " " << "z: " << x[2] << endl;
        
		point++;
	}
}

void TPZGraphEl::DrawSolution(TPZGraphNode *n,int solind,TPZDrawStyle st) {
	TPZManVector<int> sol(1,solind);
	DrawSolution(n,sol,st);
}

void TPZGraphEl::DrawSolution(TPZGraphNode * /*n*/,TPZBlock<REAL> &,TPZDrawStyle /*st*/) {
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
	long numsol = solind.NElements();
	int numvar;
	imax = 1 << res;
	int np = NPoints(n);
	int point=0;
	TPZManVector<REAL> qsi(3,0.),x(4,0.);
	TPZManVector<STATE> sol(6,0.);
	long ip = n->FirstPoint();
	while(point < np) 
	{
		QsiEta(co,imax,qsi);
		if(st == EMVStyle || st == EV3DStyle) fGraphMesh->Out() << ip++ << " ";
		for(long is=0; is<numsol; is++) 
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
#ifdef STATE_COMPLEX //AQUIFRAN
        if(fabs(sol[iv]) < 1.0e-20) sol[iv] = 0.0;
        fGraphMesh->Out() << std::real(sol[iv]) << " ";
        //fGraphMesh->Out() << fabs(sol[iv]) << " ";
#else
				if(fabs(sol[iv]) < 1.0e-20) sol[iv] = 0.0;
        fGraphMesh->Out() << sol[iv] << " ";
#endif
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

