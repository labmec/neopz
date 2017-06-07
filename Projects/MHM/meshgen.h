
#ifndef MESHGEN
#define MESHGEN

#include "pzvec.h"
#include "tpzautopointer.h"
#include "pzcmesh.h"
#include "pzfunction.h"

#include <string>

class TPZGeoMesh;

TPZGeoMesh *MalhaGeomFredQuadrada(int nelx, int nely, TPZVec<REAL> &x0, TPZVec<REAL> &x1, TPZVec<long> &coarseindices, int ndiv);

/// Solve the problem composed of a multiphysics mesh composed of compmeshes - applies to MHM and MHM-H(div)
void SolveProblem(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<TPZAutoPointer<TPZCompMesh> > compmeshes, std::string prefix, std::string configuration);

struct TElasticityExample1
{
    static void Force(const TPZVec<REAL> &x, TPZVec<STATE> &force)
    {
        TPZManVector<REAL,3> locforce(2);
        DivSigma(x, locforce);
        force[0] = -locforce[0];
        force[1] = -locforce[1];
    }
    
    static TPZAutoPointer<TPZFunction<STATE> > ForcingFunction();
    
    static TPZAutoPointer<TPZFunction<STATE> > DirichletFunction();
    
    static void Elasticity(const TPZVec<REAL> &x, REAL &elast, REAL &poisson);
    
    static void GradU(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu);
    
    template<class TVar>
    static void Sigma(const TPZVec<TVar> &x, TPZFMatrix<TVar> &sigma);
    
    template<class TVar>
    static void DivSigma(const TPZVec<TVar> &x, TPZVec<TVar> &divsigma);
    
    template<class TVar>
    static void uxy(const TPZVec<TVar> &x, TPZVec<TVar> &disp);
    
    static void Dirichlet(const TPZVec<REAL> &x, TPZVec<STATE> &disp)
    {
        TPZManVector<REAL,3> disploc(2,0.);
        uxy(x,disploc);
        for(int i=0; i<2; i++) disp[i] = disploc[i];
    }
    
    template<class TVar>
    static void graduxy(const TPZVec<TVar> &x, TPZFMatrix<TVar> &grad);

    template<class TVar>
    static void Elastic(const TPZVec<TVar> &x, TVar &Elast, TVar &nu);


};


#endif
