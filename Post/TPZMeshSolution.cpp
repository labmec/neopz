//
//  TPZMeshSolution.cpp
//  PZ
//
//  Created by Philippe Devloo on 30/10/16.
//
//

#include "TPZMeshSolution.h"
#include "TPZMaterial.h"
#include "pzcmesh.h"

//TPZCompMesh *fMesh;
//
//int64_t fGeoElIndex;
//
//TPZManVector<REAL,3> fLastLoc;
//
//int dimension;
//
//int fSolutionVarindex;
//
//int fGradVarindex;
//
//int fPolynomialOrder;
//
//public:

/** @brief Class constructor */
TPZMeshSolution::TPZMeshSolution(TPZCompMesh *cmesh, int materialid) : 
TPZRegisterClassId(&TPZMeshSolution::ClassId),
fMesh(cmesh), fMaterialIndex(materialid)
{
    fDimension = cmesh->Dimension();
    fLastLoc.Resize(fDimension,0.);
    TPZMaterial *material = cmesh->FindMaterial(materialid);
    if (!material) {
        DebugStop();
    }
    fSolutionVarindex = material->VariableIndex("Pressure");
    fGradVarindex = material->VariableIndex("Derivative");
    fNumSolutions = material->NSolutionVariables(fSolutionVarindex);
    int64_t nel = cmesh->Reference()->NElements();
    int64_t el;
    for (el=0; el<nel; el++) {
        TPZGeoEl *gel = fMesh->Reference()->Element(el);
        if (gel && gel->Dimension() == fDimension) {
            break;
        }
    }
    if (el==nel) {
        DebugStop();
    }
    fGeoElIndex = el;
} 

/** @brief Class destructor */
TPZMeshSolution::~TPZMeshSolution()
{
    
}

/**
 * @brief Performs function computation
 * @param x point coordinate which is suppose to be in real coordinate system but can be in master coordinate system in derived classes.
 * @param f function values
 * @param df function derivatives
 */
void TPZMeshSolution::Execute(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE> &df)
{
    TPZManVector<REAL,3> xcopy(x);
    fMesh->Reference()->FindElement(xcopy, fLastLoc, fGeoElIndex, fDimension);
    TPZGeoEl *gel = fMesh->Reference()->Element(fGeoElIndex);
    TPZCompEl *cel = gel->Reference();
    TPZManVector<STATE> sol(1), dsol(3);
    cel->Solution(fLastLoc,fSolutionVarindex,sol);
    cel->Solution(fLastLoc,fGradVarindex,dsol);
    int maxd = 3;
    if (df.Rows() < 3) {
        maxd = df.Rows();
    }
    for (int is=0; is<fNumSolutions; is++) {
        f[is] = sol[is];
        for (int d=0; d<maxd; d++) {
            df(d,is) = dsol[is*3+d];
        }
    }
}

/** @brief Print a brief statement */
void TPZMeshSolution::Print(std::ostream &out)
{
    out << __PRETTY_FUNCTION__ << std::endl;
    out << "NFunctions = " << NFunctions() << std::endl;
    out << "Polynomial Order = " << PolynomialOrder() << std::endl;
}


int TPZMeshSolution::ClassId() const{
    return Hash("TPZMeshSolution") ^ TPZFunction<STATE>::ClassId() << 1;
}
