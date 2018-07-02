//
//  ErrorOnFaces.cpp
//  PZ
//
//  Created by labmec on 02/03/18.
//
//
#include "pzgeoel.h"
#include "pzgeoelside.h"
#include "pzcompel.h"
#include "TPZMaterial.h"
#include "pzvec.h"
#include "pzstack.h"

#include "pzbndcond.h"

#include "ErrorOnFaces.h"

// Analysis contains the multiphysics mesh, Is it a Hdiv mesh? If not, how can we get the Hdiv mesh into the multiphysics mesh ?
bool ComputePressureJumpOnFaces_Hdiv(TPZAnalysis *analysis,int matid,TPZVec<long> &elIndex,TPZVec<int> &sideCoDim1,TPZVec<STATE> &PressureJump) {
    TPZCompMesh *cmesh = analysis->Mesh();
    if(!cmesh) return false;
//    if(!cmeshHdiv->IsHdiv())
//        std::cout << "The mesh in analysis is not a Hdiv mesh?";
    cmesh->LoadSolution(analysis->Solution());
    int ModelDimension = analysis->Mesh()->Dimension();
    
    TPZAdmChunkVector<TPZCompEl *> elvec = cmesh->ElementVec();
    long i, nel = elvec.NElements();
    
    // Identifying material and variable as pressure
    TPZMaterial *mat = cmesh->FindMaterial(matid);
    int varpress = mat->VariableIndex("Pressure");
    if(varpress < 0) return false;
    
    int dimvar = mat->NSolutionVariables(varpress);
    if(dimvar < 1) return false;
    TPZVec<STATE> sol(dimvar,0.);
    TPZVec<STATE> solneigh(dimvar,0.);
    
    //Initial dimensioning of the vectors
    elIndex.Resize(nel*6);
    sideCoDim1.Resize(nel*6);
    PressureJump.Resize(nel*6);
    long counter = 0L;
    
    /** Computing error for all elements with same dimension of the model */
    for (i = 0L; i<nel; i++) {
        TPZCompEl *el = (TPZCompEl *)elvec[i];
        // elements with model dimension
        if (!el || el->Dimension() != ModelDimension) continue;
        TPZGeoEl *gel = el->Reference();
        if(!gel) DebugStop();

        int j, nsides = gel->NSides();
        for(j=0;j<nsides;j++) {
            TPZGeoElSide gelside(gel,j);
            TPZStack<TPZCompElSide> neigh;
            if(gelside.Dimension()!=ModelDimension-1)
                continue;
            //if side has codimension 1 then compute the jump pressure
            elIndex[counter]=el->Index();
            sideCoDim1[counter]=j;
            //Computing pressure on center of the side
            TPZVec<REAL> pt(3,0.);
            STATE jump;
            gelside.CenterPoint(pt);
            gelside.ConnectedCompElementList(neigh, 0, 1);
            // Solution over computational side element
            el->Solution(pt,varpress,sol);

            if(!neigh.NElements()) {   // then it is boundary, we need identify if it is Neumann ou Dirischlet?
                TPZBndCond *matbc = dynamic_cast<TPZBndCond *> (cmesh->FindMaterial(gelside.Element()->MaterialId()));
                if(!matbc->Type()) // Dirichlet
                    PressureJump[counter]=0.;
                // Compute pressure when it is Neumann
                else {
                    jump = sol[0] - matbc->Val2()(0,0);
                    PressureJump[counter] = fabs(jump);
                }
            }
            else if(neigh.NElements()==1) { // Exist one neighboard, interior face with same level of refinement
                neigh[0].Element()->Solution(pt,varpress,solneigh);
                jump = sol[0] - solneigh[0];
                PressureJump[counter] = fabs(jump);
            }
            counter++;
        }
    }
    
    // Redimensioning the result vectors
    elIndex.Resize(counter);
    sideCoDim1.Resize(counter);
    PressureJump.Resize(counter);
    return true;
}

bool ComputeFluxJumpOnFaces_Hdiv(TPZAnalysis *analysis,int matid,TPZVec<long> &elIndex,TPZVec<int> &sideCoDim1,TPZVec<STATE> &PressureJump) {
    TPZCompMesh *cmesh = analysis->Mesh();
    if(!cmesh) return false;
    //    if(!cmeshHdiv->IsHdiv())
    //        std::cout << "The mesh in analysis is not a Hdiv mesh?";
    cmesh->LoadSolution(analysis->Solution());
    int ModelDimension = analysis->Mesh()->Dimension();
    
    TPZAdmChunkVector<TPZCompEl *> elvec = cmesh->ElementVec();
    long i, nel = elvec.NElements();
    
    // Identifying material and variable as pressure
    TPZMaterial *mat = cmesh->FindMaterial(matid);
    int varpress = mat->VariableIndex("Flux");
    if(varpress < 0) return false;
    
    int dimvar = mat->NSolutionVariables(varpress);
    if(dimvar != ModelDimension) return false;
    TPZVec<STATE> sol(dimvar,0.);
    TPZVec<STATE> solneigh(dimvar,0.);
    TPZVec<REAL> normal(dimvar,0.);
    TPZVec<REAL> normalneigh(dimvar,0.);
    
    //Initial dimensioning of the vectors
    elIndex.Resize(nel*6);
    sideCoDim1.Resize(nel*6);
    PressureJump.Resize(nel*6);
    long counter = 0L;
    
    /** Computing error for all elements with same dimension of the model */
    for (i = 0L; i<nel; i++) {
        TPZCompEl *el = (TPZCompEl *)elvec[i];
        // elements with model dimension
        if (!el || el->Dimension() != ModelDimension) continue;
        TPZGeoEl *gel = el->Reference();
        if(!gel) DebugStop();
        
        int j, nsides = gel->NSides();
        for(j=0;j<nsides;j++) {
            TPZGeoElSide gelside(gel,j);
            TPZStack<TPZCompElSide> neigh;
            if(gelside.Dimension()!=ModelDimension-1)
                continue;
            //if side has codimension 1 then compute the jump pressure
            elIndex[counter]=el->Index();
            sideCoDim1[counter]=j;
            //Computing pressure on center of the side
            TPZVec<REAL> pt(3,0.);
            STATE jump;
            gelside.CenterPoint(pt);
            gelside.ConnectedCompElementList(neigh, 0, 1);
            // Solution over computational side element
            el->Solution(pt,varpress,sol);
//            gelside.Normal(pt,gelside.Element(),0,normal);    /// ??

            if(!neigh.NElements()) {   // then it is boundary, we need identify if it is Neumann ou Dirischlet?
                TPZBndCond *matbc = dynamic_cast<TPZBndCond *> (cmesh->FindMaterial(gelside.Element()->MaterialId()));
                if(!matbc->Type()) // Dirichlet
                    PressureJump[counter]=0.;
                // Compute pressure when it is Neumann
                else {
                    jump = sol[0] - matbc->Val2()(0,0);
                    PressureJump[counter] = fabs(jump);
                }
            }
            else if(neigh.NElements()==1) { // Exist one neighboard, interior face with same level of refinement
                neigh[0].Element()->Solution(pt,varpress,solneigh);
                jump = sol[0] - solneigh[0];
                PressureJump[counter] = fabs(jump);
            }
            counter++;
        }
    }
    
    // Redimensioning the result vectors
    elIndex.Resize(counter);
    sideCoDim1.Resize(counter);
    PressureJump.Resize(counter);
    return true;
}
