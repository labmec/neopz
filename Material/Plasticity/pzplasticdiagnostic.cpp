//
//  pzplasticdiagnostic.cpp
//  PZ
//
//  Created by phil on 3/18/13.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "pzplasticdiagnostic.h"
#include "pzskylstrmatrix.h"
#include "pzskylmat.h"
#include "pzcmesh.h"
#include "pzcompel.h"

TPZPlasticDiagnostic::TPZPlasticDiagnostic(TPZCompMesh *cmesh) : fMesh(cmesh)
{
    
}

void TPZPlasticDiagnostic::CheckGlobal()
{
    TPZSkylineStructMatrix<STATE> strmat(fMesh);
    TPZFMatrix<STATE> referenceres, solution(fMesh->Solution());
    int neq = fMesh->NEquations();
    TPZFMatrix<STATE> rhs(neq,1,0.);
    solution.Resize(neq, 1);
    fMesh->Solution().Zero();
    TPZMatrix<STATE> *matrix =
        dynamic_cast<TPZMatrix<STATE>*>(strmat.CreateAssemble(referenceres, 0));
//    TPZFMatrix<STATE> compareres(solution.Rows(),1,0.);
//    strmat.Assemble(compareres, 0);
//    compareres -= referenceres;
//    std::cout << "Difference between residues " << Norm(compareres) << std::endl;
//    
//    int nel = fMesh->NElements();
//    for (int iel = 0; iel<nel; iel++) {
//        TPZCompEl *cel = fMesh->ElementVec()[iel];
//        if (!cel) {
//            continue;
//        }
//        TPZElementMatrix ek,ef, ef2;
//        cel->CalcStiff(ek,ef);
//        cel->CalcResidual(ef2);
//        ef.fMat -= ef2.fMat;
//        std::cout << "Element " << iel << " norm of difference " <<  Norm(ef.fMat) << std::endl;
//    }
    
    TPZManVector<STATE,10> residues(10,0.);
    for (int i=1; i<=10; i++) {
        TPZFMatrix<STATE> locsol(solution);
        locsol *= i*1./10.;
        fMesh->LoadSolution(locsol);
        strmat.Assemble(rhs, 0);
        rhs -= referenceres;
        TPZFMatrix<STATE> result;
        matrix->MultAdd(locsol, rhs, result,-1.,1.);
        residues[i-1] = Norm(result);
    }
    TPZManVector<STATE> rate(9,0.);
    for (int i=1; i<10; i++) {
        rate[i-1] = (log(residues[i]) - log(residues[i-1]))/(log(i+1.)-log(i*1.));
    }
    std::cout << "Residues " << residues << std::endl;
    std::cout << "Global convergence rates " << rate << std::endl;
    delete matrix;
}