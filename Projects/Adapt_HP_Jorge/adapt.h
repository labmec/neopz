// PZ includes
#include <pzmaterial.h>
#include <pzbndcond.h>
#include <pzreadmeshhr.h>
#include <pzcmesh.h>
#include <pzgeoel.h>
#include <pzcompel.h>
#include <pzanalysis.h>
#include <pzskylstrmatrix.h>
#include <pzstepsolver.h>
#include <TPZRefPattern.h>
#include <TPZCompElDisc.h>
#include <TPZInterfaceEl.h>

//Material that implements the euler equation - Developed by Tiago to Olivier
#include "pzeuler.h"

// STL includes
#include <iostream>
#include <cstdlib>
#include <string>
#include <map>
#include <list>

using namespace std;

/**
 Public interface:

 ErrorEstimation
  input: a mesh with an evaluated solution
  output: a vector with the action to be done over each element

    - ProduceGradedMeshes:
      input: a mesh
      output: a vector containing the graded meshes
      -  while the higher level greater than the minimum level take the current fine mesh
        - CoarseOneLevel:
          input: the current fine mesh
          return: the coarse mesh
          - copy the given mesh to a coarse mesh
          - identify the finest level of element
          - coarse all elements whose level are equal to the finest level
          - initialize the solution of the coarse element with the average solution of the fine elements
          - return the coarse mesh
      - insert the returned coarse mesh into the vector of graded meshes
      - use the returned coarse mesh as the fine mesh for the next step
    - EvaluateAverageSolution:
      input: a mesh with an evaluated solution
      output: a vector wiht the average solution for each state variable { rho, u, v, w, P }
      - evaluate the volume of entire domain
      - for each element:
        - take the cell center solution
        - take the volume of the element
        - for each state variable:
          - produce the product state variable solution over the volume of the element
          - increase the sum of the state variable
        - divide the sum of each state variable by the volume of the domain
    - Evaluate UHat
      input: the graded mesh
      output: map containing for each level of refinement: the vector for each element
          ontaining the vector of uhat for each state variables
        - TO BE DESCRIBED !!
    - Ealuate Detail
      input:
        initial mesh
        map containing the u_hat
        vector with the average solution
      output: the vector of evaluated detail for each element
      - for each element:
        - take the level of element
        - take the uhat corresponding to the level
        - take the solution of the element
        - for each state variable:
          - evaluate the detail evaluation: d = | u - u_hat | / average solution
        - evaluate the norm of the detail for velocity: dvelocity = sqrt ( du^2 + dv^2 + dw^2 )
        - take the detail of the element: d = max { drho, dvelocity, dP }
    - Produce the action for each element
      - if the dl > Epsl: divide
      - if the dl < Epsl and dL < Epsl: coarse
       //To be done!!!

  AdaptMesh:
    input:
      current mesh
      vector containing the actio over each element
      the refinement pattern to be used for the divided elements
    output:
      current mesh adapted
    - Take the list of elements marked to divide and produce the refinement
    - take the list of elements marked to coarse
      - for each element marked to coarse
        - verify if all brothers are also marked to coarse.
          - in affirmative case: coarse the elements
    - for each element in the mesh:
      - verifiy the difference of level between the element itself and their neighbours
      - if the level is bigger than 2 for the neighbour divide refine the element and change
        its status to none ( no refine neither coarse )
 */

// Action over an element
enum EAdaptElementAction { ENone = 0, EDivide = 1, ECoarse = -1 };

void LoadDummySolution(TPZCompMesh *cmesh);

inline void DummyFunction2(TPZVec<REAL> &co, TPZVec<REAL> &val)
{
	val.Fill(co[0]);//(1-co[0])*co[0]);
}

void GetAdaptedMesh( TPZCompMesh * cmesh, double Epsl );


//Evaluate error indicator / estimator
void ErrorEstimation ( TPZCompMesh & CMesh,
                       TPZVec < EAdaptElementAction > & DivideOrCoarsen,
					   double Epsl );

//Proceed the mesh adaption
void AdaptMesh ( TPZCompMesh & CMesh,
         TPZVec < EAdaptElementAction > & DivideOrCoarsen,
         TPZAutoPointer < TPZRefPattern > & RefPattern );

/**
 * Generates a list of meshes coarsening the original mesh from the higher level
 * to the lower level
 */
void ProduceGradedMeshes ( TPZCompMesh & OriginalMesh,
               TPZVec < TPZCompMesh * > & gradedMeshVec );
/**
 * Identify the higer level of refinement in the mesh and produce a coarse mesh
 */
TPZCompMesh * CoarsenOneLevel ( TPZCompMesh & OrignalMesh );


//Must return || < u > || = Sum over elements ( fabs ( Center_Solution ) * Element_Volume ) / Domain_Volume;
void  EvaluateAverageOfSolution ( TPZCompMesh & CompMesh,
                  TPZVec < REAL > & AverageVec );

// Evaluate u^ for all level of graded meshes
void EvaluateUHat ( TPZVec < TPZCompMesh * > & gradedMeshVec,
          map < int, vector < vector < double > > > & levelToUhatVec );

/** Evaluate d =  || u - u^ || / || < u > ||
 * where || u - u^ || :
 * | u - u ^|, u = { rho, P }
 * || u - u^ ||_L2 , u = { u,v,w }
 */
void EvaluateDetail ( TPZCompMesh & CMesh,
            TPZVec < REAL > & AverageSolutionVec,
            map < int, vector < vector < double > > > & levelToElementUhatVec,
            TPZVec < REAL > & Detail );


/**
 * Return the solution for the given element in terms of:
 * Solutions: rho, u, v, w, p
 * Gradients: Grad_x_(rho),
              Grad_y_(rho),
        Grad_z_(rho),
              Grad_x_(u),
        Grad_y_(u),
        Grad_z_(u),
        Grad_x_(v),
        Grad_y_(v),
        Grad_z_(v),
        Grad_x_(w),
        Grad_y_(w),
        Grad_z_(w),
        Grad_x_(P),
        Grad_y_(P),
        Grad_z_(p)
 */
void GetSolution ( TPZCompMesh & CMesh,
           TPZCompEl * CEl,
           TPZVec < REAL > & Solutions,
           TPZVec < REAL > & Gradients );


// verify the difference of levels of refinement between neighbours
void CheckRefinementLevel ( TPZCompMesh & CMesh,
              TPZAutoPointer < TPZRefPattern > & RefPattern );

void SelectElementsByLevel ( TPZCompMesh & CMesh,
               list< int > & SelectedElements );

void RefineElements ( TPZCompMesh & CMesh,
            list < int > & RefineList,
            TPZAutoPointer < TPZRefPattern > & RefPattern );


void PrintMeshSolution ( TPZCompMesh * cmesh, std::ostream & sout);

bool CheckReferences ( TPZCompMesh & CMesh );
bool CheckElementReferences(TPZCompEl * CEl );
