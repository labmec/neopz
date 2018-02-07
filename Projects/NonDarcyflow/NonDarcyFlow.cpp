
#include "pzgmesh.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"
#include "pzbndcond.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstepsolver.h"
#include "pzanalysis.h"
#include "pznonlinanalysis.h"
#include "TPZVTKGeoMesh.h"
#include "tpzautopointer.h"
#include "pzgeoquad.h"
#include "tpzgeoelrefpattern.h"
#include "TPZFrontStructMatrix.h"

#include "problemdata.h"
#include "pznondarcyanalysis.h"
#include "pznondarcyflow.h"

#include <time.h>

#include "fad.h"

REAL mypow(const REAL &a, const int &n)
{
	if(n == 0) return 1.;
	return a*mypow(a,n-1);
}

void FillGlobData();
TPZGeoMesh* CreateGMesh(const REAL q, const int nelr, const int nelc);
TPZCompMesh* CreateCMesh(TPZGeoMesh *gmesh);
void CalculatePressure(TPZCompMesh *cmesh, std::map<REAL,REAL> &mymap);
void CalculateLastPressureForMathematica(std::map<REAL,REAL> &pwmap);

int main()
{
	FillGlobData();
	const REAL q = 1.05;
	const int nelr = 40, nelc = 100;
	TPZGeoMesh *gmesh = CreateGMesh(q,nelr,nelc);
	TPZCompMesh *cmesh = CreateCMesh(gmesh);
	TPZNonDarcyAnalysis an(cmesh);
	TPZFrontStructMatrix <TPZFrontNonSym<STATE> > skyl(cmesh);	
//	TPZSkylineNSymStructMatrix skyl(cmesh);
	TPZStepSolver<STATE> step;
	step.SetDirect(ELU);
	an.SetSolver(step);
	an.SetStructuralMatrix(skyl);
	an.RunAll();
	
	std::map<REAL,REAL> pressuremap;
	CalculatePressure(cmesh,pressuremap);
	CalculateLastPressureForMathematica(pressuremap);
	
	return 0;

}

void FillGlobData()
{
	globData.fRw = 0.15557/2.;
	globData.fRe = 100.;
	globData.fH = 20.;
	
	globData.fKh = 9.86923266716013e-15;
	globData.fKv = 9.86923266716013e-16;
	
	globData.fPini = 27000000.;
	
	globData.fRhoRef = 1000.;
	globData.fPorRef = 0.2;
	globData.fViscRef = 0.00059;
	globData.fCRho = 0.e-10;
	globData.fCPor = 0.e-10;
	globData.fCVisc = 0.e-10;
	
	globData.fDeltaT = 100.e+5;
	globData.fNSteps = 2;	
}

TPZGeoMesh* CreateGMesh(const REAL q, const int nelr, const int nelc)
{
	const REAL rw = globData.fRw;
	const REAL re = globData.fRe;
	const REAL H = globData.fH;
	const int nr = nelr + 1;
	const int nc = nelc + 1;
	const REAL DeltaR = H/nelr;
	const int nnodes = nr*nc;
	
	TPZGeoMesh *gmesh = new TPZGeoMesh;
	if(q <= 1.00000001 || q > 2) DebugStop();
	const REAL a1 = (re-rw)*(q-1)/(mypow(q,nelc)-1);
	if(nelr <= 0 || nelc <= 0){
		DebugStop();
	}
	
	const int matid = 1, pressaobc1 = -1, pressaobc2 = -2;
	gmesh->NodeVec().Resize(nnodes);
	
	/// Criando Nos ----------------------------------------------
	TPZVec<REAL> coord(3,0.);
	
	REAL qacum;
	int nid = 0;
	for (int in = 0 ; in < nr ; in++){
		qacum = 0.;
		REAL lasts = 0, sacum = 0.;
		for (int jn = 0 ; jn < nc ; jn++)
		{
			coord[0] = rw + sacum + a1 * qacum;
			lasts = a1 * qacum;
			sacum+=lasts;
			//coord[1] = a1r * in;
			coord[1] = DeltaR*in;
			gmesh->NodeVec()[nid].SetNodeId(nid);
			gmesh->NodeVec()[nid].SetCoord(coord);
			nid++;
			if (jn == 0) qacum = 1;
			else qacum *= q;
		}
	}
	
	/// Criando Elementos ----------------------------------------
	TPZVec<int64_t> TopolQuad(4,0.);
	int64_t elid = 0;
	int icc = 0;
	for(int64_t i = 0 ; i < nnodes - nc ; i++){
		if(icc == nc-1){
			icc = 0;
			continue;
		}
		TopolQuad[0] = i;
		TopolQuad[1] = i+1;
		TopolQuad[2] = i+nc+1;
		TopolQuad[3] = i+nc;
        //tetra = new TPZGeoElRefPattern< pzgeom::TPZGeoTetrahedra> (index, TopolTetra, matElId, *gMesh);
		new TPZGeoElRefPattern <pzgeom::TPZGeoQuad> (elid,TopolQuad,matid,*gmesh);
		elid++;
		icc++;
	}
	
	/// Condicoes de Contorno -------------------------------------
	TPZVec<int64_t> line(2,0.);
	
	// CC na esquerda
	for (int i = 0; i < nelr; i++) {
		line[0] = nc*i;
		line[1] = nc*(i+1);
		new TPZGeoElRefPattern<pzgeom::TPZGeoLinear> (elid,line,pressaobc1,*gmesh);
		elid++;
	}
	
	// CC na direita
	for (int i = 0; i < nelr; i++) {
		line[0] = (nc*(i+1))-1;
		line[1] = (nc*(i+2))-1;
		new TPZGeoElRefPattern<pzgeom::TPZGeoLinear> (elid,line,pressaobc2,*gmesh);
		elid++;
	}
	
	
	gmesh->BuildConnectivity();
	
#ifdef DEBUG
	std::ofstream out("2Dmalha2d.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
	//std::ofstream textout("gmesh.txt");
	//gmesh->Print(textout);
#endif
	
	return gmesh;
}

TPZCompMesh* CreateCMesh(TPZGeoMesh *gmesh)
{
	
	const int matid = 1, WellBC = -1, FarfieldBC = -2;
	const int dirichlet = 0, neumann = 1;
	
	TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
	
	TPZNonDarcyFlow *mat = new TPZNonDarcyFlow(matid);
	cmesh->InsertMaterialObject(mat);
	
	// Condicoes de Contorno
	TPZFMatrix<STATE> val1(1,1,0.), val2(1,1,0.);
	
	// BC poco
	int wellBCType = dirichlet;
	val2(0,0) = 26000000.;
	TPZBndCond *bnd1 = mat->CreateBC(mat,WellBC,wellBCType,val1,val2);
	cmesh->InsertMaterialObject(bnd1);
	
	// BC farfield
	int farfieldBCType = dirichlet;
	val2(0,0) = 28000000.;
	TPZBndCond *bnd2 = mat->CreateBC(mat,FarfieldBC,farfieldBCType,val1,val2);
	cmesh->InsertMaterialObject(bnd2);
	
	int porder = 2;
	cmesh->SetDefaultOrder(porder);
	cmesh->AutoBuild();
	cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
	cmesh->ExpandSolution();
	
#ifdef DEBUG
	std::ofstream out("aMalhaC.txt");
	cmesh->Reference()->Print(out);
	cmesh->Print(out);
	out.close();
#endif
	
	return cmesh;
}

void CalculatePressure(TPZCompMesh *cmesh, std::map<REAL,REAL> &mymap)
{
	int nel = cmesh->NElements();
	
	TPZManVector<REAL,3> intpoint(2,1.),x(3,0.);
	REAL weight = 0.;
	
	for (int iel = 0 ; iel < nel ; iel++){
		TPZCompEl *cel = cmesh->ElementVec()[iel];
		if (!cel) continue;
		
		
		// conferindo se eh um elemento no sul
		TPZGeoEl *gel = cel->Reference();
		TPZGeoElSide gelside(gel,4);
		TPZGeoElSide neighside = gelside.Neighbour();
		if (gelside != neighside){
			continue;
		}
		 

		TPZInterpolationSpace *sp = dynamic_cast<TPZInterpolationSpace *> (cel);
		if (!sp) DebugStop();
		if (gel->Dimension() != 2) continue;
		const int var = 0;
		
		int np = 10;
		TPZVec <REAL> mypoints(np,0.);
		REAL delta = 2./(np-1);
		for (int i = 0; i < 10; i++) {
			mypoints[i] = -1 + delta * i;
		}
		
		for(int int_ind = 0; int_ind < np; ++int_ind){
			intpoint[0] = mypoints[int_ind];
			gel->X(intpoint,x);
			TPZManVector<REAL,3> Solout(3,0.);
			sp->Solution(intpoint,var,Solout);
			
			mymap[x[0]] = Solout[0];
		}
	}
	
}



void CalculateLastPressureForMathematica(std::map<REAL,REAL> &pwmap)
{
	std::ofstream out("LastPressureMap.txt");
	out << "LastPressureMap = " << "{";
	std::map<REAL,REAL>::const_iterator it = pwmap.begin();
	out << "{" << it->first << "," << it->second << "}";
	it++;
	for ( ; it != pwmap.end() ; it++){
		out << ",{" << it->first << "," << it->second << "}";
	}
	out << "};" << std::endl;
	
	out << "ListPlot[LastPressureMap,Joined->True]" << std::endl;
}

