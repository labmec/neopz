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
#include "pzstack.h"

#ifdef _AUTODIFF

#include "tpbrsteamflux.h"
#include "tpbrcellconservation.h"
#include "tpbrsolutionlist.h"

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
    
    /// Height of the permeable layer
    REAL fHeight;
    
    
public:
    
    /// data structure which defines the properties of the confinement layer
    struct TPBrConfinementLayer
    {
        TPBrConfinementLayer() : fThermalCapacity(0.5643),fConductivity(0.001758),
            fDensity(2900),fHeight(100.),fNumElements(200)
        {
            
        }
        /// thermal capacity of the confinement layer [KJ/Kg]
        REAL fThermalCapacity;
        /// conductivity of the confinement layer [KJ/(m s C)]
        REAL fConductivity;
        /// density of the confinement layer [Kg/m3]
        REAL fDensity;
        /// height of the confinement layer [m]
        REAL fHeight;
        /// number of elements for the confinement layer discretization
        int fNumElements;
    };
    
    /// object defining the confinement layer configuration
    TPBrConfinementLayer fLayer;
    
	/// state at timestep n
	TPZFMatrix<REAL> fPrevState;
	
	/// guess of state at timestep n+1
	TPZFMatrix<REAL> fNextState;
private:
	TPBrSteamFlux fInterface;
	
	TPBrCellConservation fCell;
    
    TPBRSolutionList fThermal;
    
    /// initialize the data structure of the thermal problem
    void InitializeThermalProblem();
	
	/// extract the state variables of the cell
	void ExtractCellState(int cell, TPZFMatrix<REAL> &glob, TPZVec<REAL> &cellstate);
	
	/// extract the state variables of the interface element
	void ExtractInletfaceState(TPZFMatrix<REAL> &glob, TPZVec<REAL> &interfacestate);
	
	/// extract the state variables of the interface element
	void ExtractInterfaceState(int interface1, TPZFMatrix<REAL> &glob, TPZVec<REAL> &interfacestate);
	
	/// assemble the contribution of the cell
	void AssembleCell(int cell, TPZFMatrix<REAL> &ekcell, TPZFMatrix<REAL> &efcell, TPZMatrix<REAL> &glob, TPZFMatrix<REAL> &res,
                      TPZVec<REAL> &statescales);
	
	/// assemble the contribution of the cell
	void AssembleInterface(int interface1, TPZFMatrix<REAL> &ekinterface, TPZFMatrix<REAL> &efinterface, TPZMatrix<REAL> &glob, TPZFMatrix<REAL> &res, TPZVec<REAL> &statescales);
	
	/// abstract assemble procedure
	void AssembleElement(TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef, TPZVec<int> &equation, TPZVec<int> &state, TPZMatrix<REAL> &glob, TPZFMatrix<REAL> &res);
	
	/// equation and state destination indexes for a cell
	void CellDestination(int cell, TPZVec<int> &equation, TPZVec<int> &state);
	
	/// equation and state destination indexes for an interface
	void FluxDestination(int interface1, TPZVec<int> &equation, TPZVec<int> &state);
    
	/// equation and state destination indexes for an interface
	void FluxDestinationInlet(TPZVec<int> &equation, TPZVec<int> &state);
    
    /// equation number of the first variable associated with the cell
    int FirstCellEquation(int cell)
    {
        int eq = TPBrSteamFlux::NumInletVars+(2+cell)*TPBrSteamFlux::NumFluxEq+cell*TPBrCellConservation::NumCellEq;
        return eq;
    }
    
    /// equation number of the first variable associated with the interface
    int FirstInterfaceEquation(int interface1)
    {
        int eq;
        if(interface1 == 0)
        {
            eq = TPBrSteamFlux::NumInletVars;
        }
        else
        {
            eq = TPBrSteamFlux::NumInletVars+TPBrSteamFlux::NumFluxEq+(interface1-1)*(TPBrSteamFlux::NumFluxEq+TPBrCellConservation::NumCellEq);
        }
        return eq;
    }
    
public:
	
    /// compute the scales associated with the state variables
    void StateScales(TPZStack<REAL> &scales);
	
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
	TPBrSteamMesh(int numcells, REAL temperature, REAL pressure, REAL WellRadius, REAL ReservoirRadius, REAL reservoirheight, REAL oilsaturation);
    
    /// initialize the confining layer properties
    void SetConfinementLayers(REAL height, REAL thermalcapacity, REAL density, REAL conductivity);
    
    /// initialize the mass fluxes considering [kg/s]
    void SetMassFlux(REAL waterflux, REAL oilflux);
    
    /// initialize the saturation
    void SetWaterSaturation(REAL saturation);
    
    /// initialize the cell pressure
    void SetCellPressure(int cell, REAL pressure);
    
    /// initialize the temperature
    void SetTemperature(REAL temperature);
    
    /// initialize the injection data
    void SetWaterInjection(REAL massflux, REAL temperature);
    
    /// initialize the steam quality for the inlet pressure
    void SetSteamQuality(REAL quality, REAL pressurereference);
	
	/// Define the progression of the mesh
	void SetGeometricProgression(REAL alfa)
	{
		fGeometricProgression = alfa;
        if(alfa != 1.)
        {
            fFirstCellSize = (fReservoirRadius-fWellRadius)*(1.-alfa)/(1-pow(alfa, fNumCells));
        }
        else
        {
            fFirstCellSize = (fReservoirRadius-fWellRadius)/fNumCells;
        }
        InitializeThermalProblem();
	}
	
	REAL NodeCoord(int node)
	{
        if(fGeometricProgression != 1.)
        {
            return fWellRadius+fFirstCellSize*(1-pow(fGeometricProgression, node))/(1.-fGeometricProgression);
        }
        else
        {
            return fWellRadius+fFirstCellSize*node;
        }
	}
	
	/// return the size of ith cell
	REAL CellSize(int ic)
	{
		return fFirstCellSize*pow(fGeometricProgression, (REAL)ic);
	}
	
	/// number of equations
	int NumEquations();
	
	/// bandwidth of the global system of equations
	int Bandwidth();
	
	/// perform an iteration returning the residual
	REAL Iterate();
	
	/// compute the tangent matrix and residual
	void ComputeTangent(TPZMatrix<REAL> &tangent, TPZFMatrix<REAL> &residual, TPZVec<REAL> &statescales);
	
	/// compute the tangent matrix and residual for checking purposes
	void ComputeTangent2(TPZMatrix<REAL> &tangent, TPZFMatrix<REAL> &residual, TPZVec<REAL> &statescales);
	
	/// perform a time step
	void TimeStep(REAL delt);
    
    /// Update the solution of the confinement layer
    void UpdateConfinementLayer();
    
    /// set the timestep
    void SetTimeStep(REAL delt)
    {
        fDelt = delt;
    }
    
    /// Limit correction range
    REAL LimitRange(TPZFMatrix<REAL> &prevsol, TPZFMatrix<REAL> &correction);
    
    /// Project the solution to an allowable state
    void ProjectSolution(TPZFMatrix<REAL> &solution);
    
    /// print the state of the mesh
    void Print(std::ostream &out = std::cout);
    
    /// verify the consistency of the tangent matrix
    void VerifyTangent(TPZFMatrix<REAL> &direction);
    
    /// verify the conservation of energy between prevstate and nextstate
    void EnergyBalance(REAL &leftin, REAL &rightout, REAL &energychange);
    
    /// write the steam saturation to a file
    void PrintSaturations(std::ostream &out);
    
    /// write the temperature to a file
    void PrintTemperature(std::ostream &out);
    
    /// write the pressure to a file
    void PrintPressure(std::ostream &out);

    /// write all solutions for a given time
    void PrintAll(REAL time);
};

// Nothing is compiled if _AUTODIFF isnt defined
#endif

#endif