/**
 * @file
 * @brief Experimenting with mesh for linear convection
 */

class TPZGeoMesh;
class TPZCompMesh;
class TPZGeoEl;
#include <set>
#include <pzvec.h>

template<class TVar>
class TPZFMatrix;

TPZCompMesh *CreateMeshLaxAndSod(int L,REAL &timeStep);
void InitialSolutionLaxAndSod(TPZFMatrix<REAL> &InitialSol, TPZCompMesh * cmesh);

TPZCompMesh *CreateMeshLax2D(int L,REAL &timeStep);
void InitialSolutionLax2D(TPZFMatrix<REAL> &InitialSol, TPZCompMesh * cmesh);

/** @brief For this to work PZ must be compiled with defining LinearConvection in header file pzeuler.h */
TPZCompMesh *CreateMeshLinearConvection(int L, REAL &timeStep);
void InitialSolutionLinearConvection(TPZFMatrix<REAL> &InitialSol, TPZCompMesh * cmesh);

