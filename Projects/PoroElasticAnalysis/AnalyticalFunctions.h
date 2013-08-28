/*
 * @file
 * @author Omar Duran Triana
 * @since 5/21/12.
 */
#ifndef ANALYTICALFUNCTIONS_HH
#define ANALYTICALFUNCTIONS_HH

#include "pzfmatrix.h"
#include "pzmatrix.h"
#include "pzgmesh.h"

//	Use this tools for Exponential Integral Fuction (Implementation defined at end of this main file) 
#include <math.h>		// required for fabsl(), expl() and logl()        
#include <float.h>		// required for LDBL_EPSILON, DBL_MAX

//	Internally Defined Routines
double      Exponential_Integral_Ei( double x );
long double xExponential_Integral_Ei( long double x );
static long double Continued_Fraction_Ei( long double x );
static long double Power_Series_Ei( long double x );
static long double Argument_Addition_Series_Ei( long double x);
//	End Use this tools for Exponential Integral Fuction (Implementation defined at end of this main file)

	
	// Differents functions for differents validations and uses
	
	// Exact Solution for Consolidation 1D available on: Finite element analysis of poro-elastic consolidation in porous media Standard and mixed approaches
	void ExactSolutionfiniteColumn1D(const TPZVec<REAL> &ptx, REAL timestep, TPZVec<STATE> &sol, TPZFMatrix<STATE> &flux);

	// Exact Solution for Consolidation 1D available on: Theory of linear Poroelasticty Wang 2000
	void ExactSolutionSemiInfiniteColumn1D(const TPZVec<REAL> &ptx, REAL timestep, TPZVec<STATE> &sol, TPZFMatrix<STATE> &flux);
	
	// Exact Solution for constant fluid injection/production on infinited radial linear elastic medium: Fluid Mass Sources and Point Forces in Linear Elastic Diffusive Solids	
	void ExactSolution2DLineSource(const TPZVec<REAL> &ptx, REAL timestep, TPZVec<STATE> &sol, TPZFMatrix<STATE> &flux);
	
	// Exact Solution 1D for transient constant fluid injection/production on infinited radial non-deformable medium: Green Book Rosa	
	void SolutionExactRosa1D(TPZVec<REAL> &ptx, REAL timestep, TPZVec<STATE> &sol, TPZFMatrix<STATE> &flux);

	// Analytical solution for flamant problem Elasticity Theory: Applications and Numerics 	
	void ExactSolutionFlamantProblem(TPZVec<REAL> &ptx, REAL timestep, TPZVec<STATE> &sol, TPZFMatrix<STATE> &flux);
	
	// Exact Solution 1D for PseudoPermanent constant fluid injection/production on infinited linear non-deformable medium: Green Book Rosa	
	void SolutionExactRosa1DPseudo(TPZVec<REAL> &ptx, REAL timestep, TPZVec<STATE> &sol, TPZFMatrix<STATE> &flux);

	// Exact Solution Madels problem ref	
	void ExactSolutionMandelsProblemwitheffect(const TPZVec<REAL> &ptx, REAL timestep, TPZVec<STATE> &sol, TPZFMatrix<STATE> &flux);
	
	// Here you can implemente diferent Forcing transient functions	
	void ForcingTimeDependFunction(TPZVec<REAL> &ptx, REAL TimeValue,int WhichStateVariable,double &StateVariable);

	// Here you can implement diferent ways to calculate initial pressure distribution	
	void InitialPressureDistribution(const TPZVec<REAL> &ptx, TPZVec<STATE> &sol);

#endif