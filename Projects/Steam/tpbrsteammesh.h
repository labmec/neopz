#ifndef STEAMMESHHPP
#define STEAMMESHHPP
/*
 *  SteamMesh.h
 *  PZ
 *
 *  Created by Philippe Devloo on 3/5/11.
 *  Copyright 2011 UNICAMP. All rights reserved.
 *
 */
#include "pzfmatrix.h"
#include "tpbrsteamflux.h"
#include "tpbrcellconservation.h"
/*
 Uma malha de coordenadas
 Um metodo para calcular o residuo
 Raio interno
 Raio externo
 Numero de nos
 Numero de equacoes
 
 Condicoes de entrada
 
 Pressao de saida
 
 Objeto que define as propriedades de material
 
 Fluxo entre celulas
 
 Contribuicao do interior da celula
 
 Calculo do residuo da malha
 
 Calculo da matriz tangente
 
 */
class TPBrSteamMesh
{
	
	/// number of cells of the mesh
	int fNumCells;
	
	/// fWellRadius
	REAL fWellRadius;
	
	/// Reservoir radius
	REAL fReservoirRadius;
	
	/// Geometric progression
	REAL fGeometricProgression;
	
	/// Size of first cell
	REAL fFirstCellSize;
	
	/// time step
	REAL fDelt;
	
	/// state at timestep n
	TPZFMatrix fPrevState;
	
	/// guess of state at timestep n+1
	TPZFMatrix fNextState;
	
	TPBrSteamFlux fInterface;
	
	TPBrCellConservation fCell;
	
	/// extract the state variables of the cell
	void ExtractCellState(int cell, TPZFMatrix &glob, TPZVec<REAL> &cellstate);
	
	/// extract the state variables of the interface element
	void ExtractInterfaceState(int interface, TPZFMatrix &glob, TPZVec<REAL> &interfacestate);
	
	/// assemble the contribution of the cell
	void AssembleCell(int cell, TPZFMatrix &ekcell, TPZFMatrix &efcell, TPZMatrix &glob, TPZFMatrix &res);
	
	/// assemble the contribution of the cell
	void AssembleInterface(int interface, TPZFMatrix &ekinterface, TPZFMatrix &efinterface, TPZMatrix &glob, TPZFMatrix &res);
	
	/// abstract assemble procedure
	void AssembleElement(TPZFMatrix &ek, TPZFMatrix &ef, TPZVec<int> &equation, TPZVec<int> &state, TPZMatrix &glob, TPZFMatrix &res);
	
	/// equation and state destination indexes for a cell
	void CellDestination(int cell, TPZVec<int> &equation, TPZVec<int> &state);
	
	/// equation and state destination indexes for an interface
	void FluxDestination(int interface, TPZVec<int> &equation, TPZVec<int> &state);
	
public:
	
	/// Flux access method
	TPBrSteamFlux &Interface()
	{
		return fInterface;
	}
	
	/// Cell access method
	TPBrCellConservation &Cell()
	{
		return fCell;
	}
	
	/// initialize the state variables
	TPBrSteamMesh(int numcells, REAL temperature, REAL pressure, REAL WellRadius, REAL ReservoirRadius, REAL oilsaturation);
	
	/// Define the progression of the mesh
	void SetGeometricProgression(REAL alfa)
	{
		fGeometricProgression = alfa;
		fFirstCellSize = (fReservoirRadius-fWellRadius)*(1.-alfa)/(1-pow(alfa, fNumCells));
	}
	
	REAL NodeCoord(int node)
	{
		return fWellRadius+fFirstCellSize*(1-pow(fGeometricProgression, node))/(1.-fGeometricProgression);
	}
	
	/// return the size of ith cell
	REAL CellSize(int ic)
	{
		return fFirstCellSize*pow(fGeometricProgression, (double)ic);
	}
	
	/// number of equations
	int NumEquations();
	
	/// bandwidth of the global system of equations
	int Bandwidth();
	
	/// perform an iteration returning the residual
	REAL Iterate();
	
	/// compute the tangent matrix and residual
	void ComputeTangent(TPZMatrix &tangent, TPZFMatrix &residual);
	
	/// perform a time step
	void TimeStep(REAL delt);

};
#endif