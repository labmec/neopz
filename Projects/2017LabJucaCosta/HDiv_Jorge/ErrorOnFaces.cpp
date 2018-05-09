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
#include "pzintel.h"
#include "TPZMaterial.h"
#include "pzvec.h"
#include "pzstack.h"

#include "pzbndcond.h"

#include "ErrorOnFaces.h"

// Auxiliar function to know the maxime index of the comp elements of a mesh
long MaxCompElementsIndex(TPZCompMesh *cmesh) {
    long i, nel = cmesh->NElements();
    long maxelindex = 0;
    for(i=0L;i<nel;i++) {
        long index = cmesh->ElementVec()[i]->Index();
        if(index>maxelindex)
            maxelindex = index;
    }
    return maxelindex;
}

bool IdentifyingFaces(TPZCompMesh *cmesh,TPZStack<TPZCompElSide> &Faces, TPZStack<TPZCompElSide> &AnotherSideFaces) {

    if(!cmesh) return false;
    cmesh->Reference()->ResetReference();
    cmesh->LoadReferences();
    
    // To know the dimension of the computational elements over we search
    int ModelDimension = cmesh->Dimension();
    int MaxIndex = MaxCompElementsIndex(cmesh);
    
    TPZAdmChunkVector<TPZCompEl *> elvec = cmesh->ElementVec();
    long i, nel = elvec.NElements();
    // To check if the Face was identified already
    TPZFMatrix<int> FoundedFaces(MaxIndex+1,27,0);
    
    /** Finding faces at mesh: Face is boundary whether it has no neighboard, Face is inner if it has (only one) neighboard. */
    for (i = 0L; i<nel; i++) {
        TPZCompEl *el = (TPZCompEl *)elvec[i];
        // elements with model dimension
        if (!el || el->Dimension() != ModelDimension) continue;
        
        int j, nsides = el->Reference()->NSides();
        for(j=0;j<nsides;j++) {
            TPZCompElSide celside(el,j);
            TPZStack<TPZCompElSide> neigh;
            TPZCompElSide el_neigh;
            // Only over the sides of the codimension 1, if it is not founded
            if(celside.Reference().Dimension() != ModelDimension-1 || FoundedFaces(el->Index(),j))
                continue;
            celside.EqualLevelElementList(neigh, 1, 1);
            if(!neigh.NElements()) {
                el_neigh = celside.LowerLevelElementList(1);
                // existing or not el_neigh we can to register a face
                Faces.push_back(celside);
                AnotherSideFaces.push_back(el_neigh);
                FoundedFaces(el->Index(),j) = 1;
                if(el_neigh.Element()) {
                    FoundedFaces(el_neigh.Element()->Index(),el_neigh.Side()) = 1;
                }
            }
            else if(neigh.NElements() == 1) {
                Faces.push_back(celside);
                AnotherSideFaces.push_back(neigh[0]);
                FoundedFaces(el->Index(),j) = 1;
                FoundedFaces(neigh[0].Element()->Index(),neigh[0].Side()) = 1;
            }

        }
    }
    return true;
}

void TensorKFunction(TPZVec<REAL> &x,TPZFMatrix<REAL> &K) {
    REAL alpha = 1.;
    int p = K.Rows();
    K.Zero();
    for(int i=0;i<p;i++)
        K(i,i) = alpha;
}

/** To compute Cmin and Cmax of the elliptic equation based on tensor K */
bool ComputeCMinAndCMaxFromTensorK(void (*fp)(TPZVec<REAL> &loc,TPZFMatrix<REAL> &K),TPZCompMesh* cmesh,REAL &Cmin, REAL &Cmax) {
    TPZVec<REAL> pt(3,0.);
    REAL norm, prod;
    int dim = cmesh->Dimension();
    long nel = cmesh->NElements();
    TPZFMatrix<REAL> TensorK(dim,dim,0.);
    TPZFMatrix<REAL> Qsi(dim,1,0.), Psi(dim,1,0.);
    Cmin = Cmax = 1.;
    int i, j;

    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = cmesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZMaterial * material = cel->Material();

        if (!material) continue;
        int dimcel = cel->Dimension();
        TPZAutoPointer<TPZIntPoints> intrule = ((TPZInterpolationSpace*)cel)->GetIntegrationRule().Clone();
        
        TPZManVector<REAL,3> intpoint(dimcel);
        TPZManVector<REAL,3> pt(3,0.);
        REAL weight;

        int nintpoints = intrule->NPoints();
        
        for(int nint = 0; nint < nintpoints; nint++) {
            norm = 0.; Qsi.Zero(); Psi.Zero(); prod = 0.;
            intrule->Point(nint,intpoint,weight);

            TPZGeoEl * ref = cel->Reference();
            if(!ref) continue;
            ref->X(intpoint, pt);
            
            fp(pt,TensorK);
//            TensorKFunction(pt,TensorK);
            for(i=0; i<dim; i++) {
                norm += pt[i]*pt[i];
                for(j=0;j<dim;j++)
                    Qsi(i,0) += TensorK(i,j)*pt[j];
                prod += pt[i]*Qsi(i,0);
            }
            if(IsZero(norm)) continue;
            prod /= norm;
            
            Cmin = (Cmin < prod) ? Cmin : prod;
            Cmax = (Cmax < prod) ? prod : Cmax;
//            if(fp) {
            // if exist function to calculate Tensor K
        }

    }
    return true;
    for(long i=0;i < cmesh->Reference()->NodeVec().NElements(); i++) {
        cmesh->Reference()->NodeVec()[i].GetCoordinates(pt);
        norm = 0.;
        for(j=0;j<dim;j++) {
            norm += pt[j]*pt[j];
        }
        if(dim==1) prod = pt[0]*pt[0]*TensorK(0,0);
        else if(dim==2) prod = pt[0]*pt[0]*TensorK(0,0)+pt[0]*pt[1]*(TensorK(0,1)+TensorK(1,0))+pt[1]*pt[1]*TensorK(1,1);
        else return false;
        if(!IsZero(norm)) {
            prod /= norm;
            Cmin = (Cmin < prod) ? Cmin : prod;
            Cmax = (Cmax < prod) ? prod : Cmax;
        }
    }
    return true;
}

bool ComputePressureJumpOnFaces(TPZCompMesh *cmesh,int matid,TPZStack<TPZCompElSide> &Faces, TPZStack<TPZCompElSide> &AnotherSideFaces,STATE &Error) {
    
    if(!cmesh) return false;
    int ModelDimension = cmesh->Dimension();
    
    //    TPZAdmChunkVector<TPZCompEl *> elvec = cmesh->ElementVec();
    long i, nel = Faces.NElements();   //elvec.NElements();
    
    // Identifying material and variable as pressure
    TPZMaterial *mat = cmesh->FindMaterial(matid);
    int varpress = mat->VariableIndex("Pressure");
    if(varpress < 0) return false;
    
    int dimvar = mat->NSolutionVariables(varpress);
    if(dimvar < 1) return false;
    TPZVec<STATE> sol(dimvar,0.);
    TPZVec<STATE> solneigh(dimvar,0.);
    
    // Initializing
    REAL volEl = 0.;
    Error = 0.;
    
    /** Computing error for all elements with same dimension of the model */
    for (i = 0L; i<nel; i++) {
        TPZCompElSide celside = Faces.Pop();
        TPZCompElSide neighcelside = AnotherSideFaces.Pop();
        
        // element with higher level
        TPZCompEl *el = (TPZCompEl *)celside.Element();
        // elements with model dimension
        if (!el || el->Dimension() != ModelDimension) continue;
        TPZGeoEl *gel = el->Reference();
        if(!gel) DebugStop();
        
        // If faces is boundary made nothing
        if(!neighcelside.Element()) {
            AnotherSideFaces.Pop();
            continue;
        }
        
        //if side has codimension 1 and it is inner then compute the jump pressure
        //Computing pressure on center of the side
        TPZManVector<REAL,3> pt(3,0.);
        TPZManVector<REAL,3> pt_el(3,0.);
        TPZManVector<REAL,3> pt_n(3,0.);
        TPZManVector<REAL,3> pt_el_n(3,0.);
        volEl = celside.Reference().Area();
        
        // pt - point on face with codimension 1
        celside.Reference().CenterPoint(pt);
        
        TPZTransform<> tr;
        TPZGeoElSide geosideh(gel,gel->NSides()-1);
        tr = celside.Reference().SideToSideTransform(geosideh);
        // pt_el - point on faces with dimension ModelDimension
        tr.Apply(pt, pt_el);
        
        // Solution over computational side element
        el->Solution(pt_el,varpress,sol);
        
        // working on faces from neighboard element with commom face
        TPZCompElSide celside_n = AnotherSideFaces.Pop();
        TPZGeoElSide gelside_n = celside_n.Reference();
        gelside_n.CenterPoint(pt_n);
        TPZGeoElSide gelsideh_n(gelside_n.Element(),gelside_n.Element()->NSides()-1);
        tr = gelside_n.SideToSideTransform(gelsideh_n);
        tr.Apply(pt_n,pt_el_n);
        celside.Element()->Solution(pt_el_n,varpress,solneigh);
        Error += volEl*(sol[0] - solneigh[0]);
    }
    
    return true;
}

