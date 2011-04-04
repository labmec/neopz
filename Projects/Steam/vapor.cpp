#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/*
 *  steaminjection.cpp
 *  girkmann
 *
 *  Created by Agnaldo on 1/6/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "TPBConservacao1DMarx.h"
#include "tpbrthermaldisc.h"

#include "pzseqsolver.h"
#include <math.h>

#ifdef _AUTODIFF

void ScaleFactor(TPZFMatrix &tangentmatrix, TPZFMatrix &residualmatrix, TPZManVector<REAL> &scalevalues, TPZManVector<REAL> &statescalevalues );
void ScaleFactorSol(TPZFMatrix &residualmatrix, TPZManVector<REAL> &statescalevalues );

int mainvapor()
{
	
	//REAL domainsize = 100.;
//	int nelements = 50;
//	REAL cp = 1.;
//	REAL K = 1.;
//	REAL initialtemp = 0.;
//	REAL flux;
//	TPBRThermalDisc discrete(domainsize,nelements,cp,K,initialtemp);
//	TPZFMatrix sol(nelements+1,1,0.), nextsol(nelements+1,1,0.);
//	discrete.SetTimeStep(1.);
//	discrete.ComputeStiffness();
//	discrete.NextSolution(1., sol,nextsol,flux);
//	nextsol.Print("Next Solution", std::cout);
//	discrete.NextSolution(1., nextsol,nextsol,flux);
//	nextsol.Print("Next Solution", std::cout);
//	return 0;
	
	REAL PI = 4*atan(1.);	
	TPBrCellMarx first;
	
	//----------------- dados de entrada ------------------------------
	//dados numerico
	REAL TimeStep =10.;//(1260. s = tempo para atingir a Energia máxima)
		
	//dados da celula
	REAL rint = 0.15;
	REAL rext = 1.5;
	REAL CellSize = 1.;
	REAL LeftArea = PI*rint*rint;
	REAL RightArea = PI*rext*rext;
	REAL CellVolume = PI*(rext-rint)*(rext-rint)*CellSize;
	
	//dados da rocha
	REAL MaterialPermeability = 0.8e-12;
	PhysicalProperties minharocha(0,1);
	
	//dados de injecao
	REAL PressureWater(2.e6);
	REAL tempReservtorio = 98.0;
	REAL Temperature(tempReservtorio);
	TPZManVector<REAL> InitialSaturation(3,0.);
	InitialSaturation[TPBrCellMarx::EOil] = 0.0170871;
	InitialSaturation[TPBrCellMarx::EWater] = 1.-0.0170871;
	InitialSaturation[TPBrCellMarx::ESteam] = 0. ;
	
	TPZManVector<REAL> Massflux(3,0.);
	REAL titularidadeVapor, txInjecaoMassa, mv, mw;
	//int numX = 20;
	txInjecaoMassa =0.677507;
	//REAL taxa;
	//REAL temp1, temp2;
	//temp2 = numX-1;
	//TPZFMatrix VexpVcel(numX,2);
	
	//for (int i=0; i<numX; i++) {
		//temp1 = i;
		titularidadeVapor = /*temp1/temp2;//*/0.82;
		mv = titularidadeVapor*txInjecaoMassa;
		mw = txInjecaoMassa - mv;
		Massflux[TPBrCellMarx::EOil] = 0.0;
		Massflux[TPBrCellMarx::EWater] = mw;//0.121951;
		Massflux[TPBrCellMarx::ESteam] = mv;//0.555556
		//-----------------------------------------------------------------------------------------
		
		TPZManVector<REAL> initial(TPBrCellMarx::NUMVARS,0.), residual(TPBrCellMarx::NUMVARS,0.);
		TPZManVector<REAL> leftstate(TPBrCellMarx::NUMVARS,0.),rightstate(TPBrCellMarx::NUMVARS,0.);
		
		first.SetMaterialProperty(MaterialPermeability, minharocha);
		first.SetGeometry(CellVolume,LeftArea,RightArea,CellSize);
		first.SetCellState(PressureWater,InitialSaturation,Temperature, TimeStep);
		first.SetInjectionState(PressureWater, Massflux, leftstate);
		
		first.InitializeState(initial);
		first.TotalResidual(leftstate,initial,residual);
		
		TPZManVector<TFad<TPBrCellMarx::NUMVARS,REAL> > tangent(TPBrCellMarx::NUMVARS,0.), state(TPBrCellMarx::NUMVARS,0.);
		first.InitializeState(state);
		first.TotalResidual(leftstate, state, tangent);
		
		//-------------------------- FATORES DE ESCALA  ---------------------------------------------------
		TPZManVector <REAL> vecState(initial), scalevalues(19,0.), statescalevalues(19,0.);
		TPZFMatrix scalevaluesmatrix,statescalevaluesmatrix;
		vecState[TPBrCellMarx::EMassFluxOil] = leftstate[TPBrCellMarx::EMassFluxOil];
		vecState[TPBrCellMarx::EMassFluxWater] = leftstate[TPBrCellMarx::EMassFluxWater];
		vecState[TPBrCellMarx::EMassFluxSteam] = leftstate[TPBrCellMarx::EMassFluxSteam];
		first.ReferenceResidualValues(leftstate, scalevalues);
		first.ReferenceStateValues(vecState,statescalevalues);
		
		first.ExtractMatrix(scalevalues,scalevaluesmatrix);
		first.ExtractMatrix(statescalevalues,statescalevaluesmatrix);
	
		//----------------------------- corrigir residuo e tangente -----------------------------
       	TPZFNMatrix<20> tangentmatrix,residualmatrix,statematrix;
		first.ExtractMatrix(tangent,tangentmatrix);
		first.ExtractMatrix(residual,residualmatrix);
		first.ExtractMatrix(initial,statematrix);
		
		ScaleFactor(tangentmatrix, residualmatrix, scalevalues, statescalevalues);
	
		//----------------------------Metodo de Newton--------------------------------------	
		//cout << "\n ======= Metodo de Newton =======\n"; 
		REAL norma;
		int IterNewton =0;
		TPZVec<int> ind;
		norma = Norm(residualmatrix);
		//cout<< "Norma --> " << norma <<endl;
		while (IterNewton <100 && norma > 1.e-10) {
			IterNewton++;
			//cout<< "iter --> " << IterNewton <<endl;
			tangentmatrix.Decompose_LU(ind);
			tangentmatrix.Substitution(&residualmatrix, ind);
			//tangentmatrix.SolveDirect(residualmatrix, ELU);
			
			//residualmatrix.Print("Residualmatrix = ",cout, EMathematicaInput);
			
			// multiplicar o valor pelo scalestate
			ScaleFactorSol(residualmatrix, statescalevalues);
			
			statematrix -= residualmatrix;
			//statematrix.Print("statematrix = ",cout, EMathematicaInput);		
			
			first.ConvertState(statematrix,rightstate);
			first.ConvertState(statematrix,state);
			first.TotalResidual(leftstate, rightstate, residual);
			first.TotalResidual(leftstate, state, tangent);
			
			tangentmatrix.Zero();
			residualmatrix.Zero();
			first.ExtractMatrix(residual,residualmatrix);
			first.ExtractMatrix(tangent,tangentmatrix);
			
			// aplicar os fatores de escala novamente no tangentmatrix e residual
			ScaleFactor(tangentmatrix, residualmatrix, scalevalues, statescalevalues);
			
			norma =Norm(residualmatrix);
			//cout<< "\nNorma --> " << norma <<endl;
		}
				
		//cout << "\n Numero de Iteracoes = "<< IterNewton <<endl<<endl;
		//tangentmatrix.Print("tangentmatrix = ",cout, EMathematicaInput);
		//cout << endl;
		residualmatrix.Print("Residualmatrix = ",cout, EMathematicaInput);
		cout<<endl;
		statematrix.Print("statematrixFim = ",cout, EMathematicaInput);
		
		//first.RazaoVexpVcel(rightstate,taxa);
//		VexpVcel(i,0)=titularidadeVapor;
//		VexpVcel(i,1)=taxa;
	//}//for i
//	
	//VexpVcel.Print("VexpVcel = ",cout, EMathematicaInput);
	return EXIT_SUCCESS;
};


//----------------------------------------------------------------------------------------------------------
void ScaleFactor(TPZFMatrix &tangentmatrix, TPZFMatrix &residualmatrix, TPZManVector<REAL> &scalevalues, TPZManVector<REAL> &statescalevalues ){
	
	int numvar = TPBrCellMarx ::NUMVARS;
	for (int ir =0; ir<numvar; ir++) {
		residualmatrix(ir,0) = residualmatrix(ir,0)/scalevalues[ir] ;
	}
	
	TPZFMatrix  tangentmatrixtr;
	for (int it = 0; it<numvar; it++) {
		for (int jt = 0; jt<numvar; jt++) {
			tangentmatrix(it,jt) = tangentmatrix(it,jt)/scalevalues[it];  
		}
	}
	tangentmatrix.Transpose(&tangentmatrixtr);
	for (int itr = 0; itr<numvar; itr++) {
		for (int jtr = 0; jtr<numvar; jtr++) {
			tangentmatrixtr(itr,jtr) = tangentmatrixtr(itr,jtr)*statescalevalues[itr];  
		}
	}
	tangentmatrixtr.Transpose(&tangentmatrix);
}

//-----------------------------------------------------------------------------------------------------
void ScaleFactorSol(TPZFMatrix &residualmatrix, TPZManVector<REAL> &statescalevalues ){
	
	int numvar = TPBrCellMarx ::NUMVARS;
	for (int ir =0; ir<numvar; ir++) {
		residualmatrix(ir,0) = residualmatrix(ir,0)*statescalevalues[ir] ;
	}	
}

#endif
