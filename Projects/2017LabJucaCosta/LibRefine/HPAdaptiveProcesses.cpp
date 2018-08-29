 /**
 * @file Poisson 3D in hexahedra with shock problem
 */

#include "pzgeopoint.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "pzgeotetrahedra.h"
#include "pzgeopyramid.h"
#include "pzgeoprism.h"
#include "TPZGeoLinear.h"
#include "TPZGeoCube.h"
#include "pzgnode.h"
#include "pzgeoel.h"
#include "pzgeoelside.h"
#include "pzgeoelbc.h"
#include "tpzgeoelrefpattern.h"

#include "pzintel.h"
#include "pzcompel.h"
#include "pzcmesh.h"

#include "pzmatrix.h"
#include "pzvec.h"
#include "pzadmchunk.h"
#include "pzvec_extras.h"

#include "pzlog.h"
#include "pzcheckgeom.h"
#include "pzcheckmesh.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSSpStructMatrix.h"
#include "pzsbstrmatrix.h"
#include "pzbstrmatrix.h"
#include "pzfstrmatrix.h"
#include "pzstepsolver.h"
#include "TPZFrontStructMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZParSkylineStructMatrix.h"

#include "TPZMaterial.h"
#include "pzbndcond.h"
#include "pzelasmat.h"
#include "pzplaca.h"
#include "pzpoisson3d.h"
#include "pzmathyperelastic.h"
#include "pzmattest3d.h"
#include "pzmatplaca2.h"

#include "pzfunction.h"

#include "pzgengrid.h"
#include "TPZExtendGridDimension.h"
#include "TPZReadGIDGrid.h"
#include "TPZVTKGeoMesh.h"

#include "pzshapelinear.h"
#include "TPZRefPatternTools.h"

#include <time.h>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <cmath>

#include "HPAdaptiveProcesses.h"


/**
* Get Global L2 Error for solution and the L2 error for each element.
* Return the maxime L2 error by elements. Also return in MinErrorByElement argument the minime L2 error for all elements of the mesh.
*/
bool ProcessingErrorUAndDUKnowingExactSol(TPZAnalysis &analysis, TPZVec<REAL> &ErrorVecByIteration, int nref, TPZVec<STATE> &ErrorU, TPZVec<STATE> &ErrorDU)
{
	TPZCompMesh *cmesh = analysis.Mesh();
	cmesh->LoadSolution(analysis.Solution());
	int ModelDimension = analysis.Mesh()->Dimension();

	TPZAdmChunkVector<TPZCompEl *> elvec = cmesh->ElementVec();
	TPZManVector<REAL, 10> errors(10);
	errors.Fill(0.0);
	int64_t i, nel = elvec.NElements();
	int64_t nerrors = 0L;
	ErrorU.Resize(nel, 0.0);
	// The last position will be store the maxime value of the gradient errors
	ErrorDU.Resize(nel, 0.0);

	/** Computing error for all elements with same dimension of the model */
	for (i = 0L; i<nel; i++) {
		TPZCompEl *el = (TPZCompEl *)elvec[i];
		if (!el || el->Dimension() != ModelDimension) continue;
		errors.Fill(0.0);
		el->EvaluateError(analysis.fExact, errors, 0);
		nerrors = errors.NElements();
		REAL vol = el->VolumeOfEl();
		for (int ier = 0; ier < nerrors; ier++) {
			errors[ier] *= vol;
			ErrorVecByIteration[nerrors*nref + ier] += errors[ier];
		}

		// L2 error for each element
		ErrorU[i] = errors[1];
		ErrorDU[i] = errors[2];
	}
	return true;
}
bool ProcessingErrorUKnowingExactSol(TPZAnalysis &analysis, TPZVec<REAL> &ErrorVecByIteration, int nref, TPZVec<STATE> &ErrorU)
{
	TPZCompMesh *cmesh = analysis.Mesh();
	cmesh->LoadSolution(analysis.Solution());
	int ModelDimension = analysis.Mesh()->Dimension();

	TPZAdmChunkVector<TPZCompEl *> elvec = cmesh->ElementVec();
	TPZManVector<REAL, 10> errors(10);
	errors.Fill(0.0);
	int64_t i, nel = elvec.NElements();
	int64_t nerrors = 0L;
	ErrorU.Resize(nel, 0.0);

	/** Computing error for all elements with same dimension of the model */
	for (i = 0L; i<nel; i++) {
		TPZCompEl *el = (TPZCompEl *)elvec[i];
		if (!el || el->Dimension() != ModelDimension) continue;
		errors.Fill(0.0);
		el->EvaluateError(analysis.fExact, errors, 0);
		nerrors = errors.NElements();
		REAL vol = el->VolumeOfEl();
		for (int ier = 0; ier < nerrors; ier++) {
			errors[ier] *= vol;
			ErrorVecByIteration[nerrors*nref + ier] += errors[ier];
		}

		// L2 error for each element
		ErrorU[i] = errors[1];
	}
	return true;
}

void ApplyHPRefinement(TPZCompMesh *cmesh, TPZVec<int64_t> &PRef, int MaxPOrder,TPZVec<int64_t> &HRef,int MaxHLevel) {
	int64_t iel, nelhrefs = HRef.NElements(), nelprefs = PRef.NElements();
	int64_t nels = cmesh->NElements();

	TPZManVector<int64_t, 27> subels;
	TPZManVector<int64_t, 27> subsubels;

	// Doing P Refinement
	int pelement = 0;
	TPZGeoEl *gel = 0;
	TPZInterpolationSpace *intel;
	for (iel = 0; iel<nelprefs; iel++) {
		intel = 0;
		intel = dynamic_cast<TPZInterpolationSpace* > (cmesh->Element(PRef[iel]));
		if (!intel || intel->Dimension() != cmesh->Dimension()) continue;
		pelement = intel->GetPreferredOrder(); //->PreferredSideOrder(gel->NSides() - 1);
		if (pelement < MaxPOrder)
			intel->PRefine(pelement + 1);
	}
	cmesh->ExpandSolution();

	// Doing H Refinement
	for (iel = 0; iel<nelhrefs; iel++) {
		bool twice = false;
		if (HRef[iel]<0) {
			twice = true;
			HRef[iel] *= -1;
		}
		subels.Resize(0);
		intel = 0;
		intel = dynamic_cast<TPZInterpolatedElement* > (cmesh->Element(HRef[iel]));
		if (!intel || intel->Dimension() != cmesh->Dimension()) continue;
		gel = intel->Reference();
		if (!gel) DebugStop();
		if (gel->Level() < MaxHLevel) {
			intel->Divide(intel->Index(), subels, 1);
			cmesh->ElementVec().SetFree(HRef[iel]);
			//            intel = 0;
		}
		if (twice) {
			for (int64_t isub_el = 0; isub_el<subels.NElements(); isub_el++) {
				subsubels.Resize(0);
				TPZCompEl * isub_cel = cmesh->ElementVec()[subels[isub_el]];
				if (!isub_cel || isub_cel->Dimension() != cmesh->Dimension()) continue;
				isub_cel->Divide(subels[isub_el], subsubels);
				cmesh->ElementVec().SetFree(subels[isub_el]);
			}
			twice = false;
		}
	}
	cmesh->ExpandSolution();

	// Printing information stored
	PrintNRefinementsByType(nels, cmesh->NElements(), nelhrefs, nelprefs);
}


bool ApplyingHPAdaptiveStrategyBasedOnU_I(TPZCompMesh *cmesh,TPZVec<STATE> &ErrorU,TPZVec<STATE> &ErrorDU,TPZVec<REAL> &Tol, int MaxPOrder, int MaxHLevel,std::ostream &out) {
    if(!cmesh) return false;
    int64_t iel, nelhrefs = 0, nelprefs = 0;
    int64_t nels = cmesh->NElements();
    int64_t nelsR2=0L, nelsR1=0L, nelsR0=0L, nelsR=0L;
    
    TPZVec<int64_t> HRef(nels,0L), PRef(nels,0L);
    
    // Applying hp refinement only for elements with dimension as model dimension
    out << " Refinando malha com " << nels  << " elementos e " << cmesh->NEquations() << " equacoes.\n";

    // Applying tolerance limits to define whether the element will be h-, p-, hp-refined or not. Implementation of the hp-adaptive table.
    // Note: Some elements can to have p and h refinements. But to indicate wheter the element must to refine twice h-ref, we have changed the index by -index
    for(iel=0L;iel<nels;iel++) {
        TPZCompEl *cel = cmesh->Element(iel);
        if(!cel || cel->Dimension() != cmesh->Dimension()) continue;
        if(ErrorU[iel] > Tol[2]) {
            HRef[nelhrefs++] = -1*iel;
            nelsR2++;
        }
        else if(ErrorU[iel] > Tol[1]) {
            PRef[nelprefs++] = iel;
            HRef[nelhrefs++] = iel;
            nelsR1++;
        }
        else if(ErrorU[iel] > Tol[0]) {
            PRef[nelprefs++] = iel;
            nelsR0++;
        }
        else
            nelsR++;
    }
    out << "I - Regioes:\t" << nelsR << "\t" << nelsR0 << "\t" << nelsR1 << "\t" << nelsR2 << "\n";
	HRef.Resize(nelhrefs);
	PRef.Resize(nelprefs);

	// Doing h and p refinements
	ApplyHPRefinement(cmesh,PRef, MaxPOrder, HRef, MaxHLevel);
    
    // If no exists any element to refine, the tolerance was reached
    if(!nelhrefs && !nelprefs)
        return true;
    return false;
}
bool ApplyingHPAdaptiveStrategyBasedOnU_II(TPZCompMesh *cmesh,TPZVec<STATE> &ErrorU,TPZVec<STATE> &ErrorDU,TPZVec<REAL> &Tol, int MaxPOrder, int MaxHLevel,std::ostream &out) {
    if(!cmesh) return false;
    int64_t iel, nelhrefs = 0, nelprefs = 0;
    int64_t nels = cmesh->NElements();
    int64_t nelsR2=0L, nelsR1=0L, nelsR0=0L, nelsR=0L;
    
    TPZVec<int64_t> HRef(nels,0L), PRef(nels,0L);
    
    // Applying hp refinement only for elements with dimension as model dimension
    out << " Refinando malha com " << nels  << " elementos e " << cmesh->NEquations() << " equacoes.\n";
    
    // Applying tolerance limits to define whether the element will be h-, p-, hp-refined or not. Implementation of the hp-adaptive table.
    // Note: Some elements can to have p and h refinements. But to indicate wheter the element must to refine twice h-ref, we have changed the index by -index
    for(iel=0L;iel<nels;iel++) {
        TPZCompEl *cel = cmesh->Element(iel);
        if(!cel || cel->Dimension() != cmesh->Dimension()) continue;
        if(ErrorU[iel] > Tol[2]) {
            HRef[nelhrefs++] = iel;
            nelsR2++;
        }
        else if(ErrorU[iel] > Tol[1]) {
            nelsR1++;
            PRef[nelprefs++] = iel;
            HRef[nelhrefs++] = iel;
        }
        else if(ErrorU[iel] > Tol[0]) {
            nelsR0++;
            PRef[nelprefs++] = iel;
        }
        else
            nelsR++;
    }
    out << "II - Regioes:\t" << nelsR << "\t" << nelsR0 << "\t" << nelsR1 << "\t" << nelsR2 << "\n";
    HRef.Resize(nelhrefs);
    PRef.Resize(nelprefs);
    
    // Doing h and p refinements
    ApplyHPRefinement(cmesh,PRef, MaxPOrder, HRef, MaxHLevel);
    
    // If no exists any element to refine, the tolerance was reached
    if(!nelhrefs && !nelprefs)
        return true;
    return false;
}
bool ApplyingHPAdaptiveStrategyBasedOnU_III(TPZCompMesh *cmesh,TPZVec<STATE> &ErrorU,TPZVec<STATE> &ErrorDU,TPZVec<REAL> &Tol, int MaxPOrder, int MaxHLevel,std::ostream &out) {
    if(!cmesh) return false;
    int64_t iel, nelhrefs = 0, nelprefs = 0;
    int64_t nels = cmesh->NElements();
    int64_t nelsR2=0L, nelsR1=0L, nelsR0=0L, nelsR=0L;
    
    TPZVec<int64_t> HRef(nels,0L), PRef(nels,0L);
    
    // Applying hp refinement only for elements with dimension as model dimension
    out << " Refinando malha com " << nels  << " elementos e " << cmesh->NEquations() << " equacoes.\n";
    
    // Applying tolerance limits to define whether the element will be h-, p-, hp-refined or not. Implementation of the hp-adaptive table.
    // Note: Some elements can to have p and h refinements. But to indicate wheter the element must to refine twice h-ref, we have changed the index by -index
    for(iel=0L;iel<nels;iel++) {
        TPZCompEl *cel = cmesh->Element(iel);
        if(!cel || cel->Dimension() != cmesh->Dimension()) continue;
        if(ErrorU[iel] > Tol[2]) {
            HRef[nelhrefs++] = iel;
            PRef[nelprefs++] = iel;
            nelsR2++;
        }
        else if(ErrorU[iel] > Tol[1]) {
            nelsR1++;
            HRef[nelhrefs++] = iel;
        }
        else if(ErrorU[iel] > Tol[0]) {
            nelsR0++;
            PRef[nelprefs++] = iel;
        }
        else
            nelsR++;
    }
    out << "III - Regioes:\t" << nelsR << "\t" << nelsR0 << "\t" << nelsR1 << "\t" << nelsR2 << "\n";
    HRef.Resize(nelhrefs);
    PRef.Resize(nelprefs);
    
    // Doing h and p refinements
    ApplyHPRefinement(cmesh,PRef, MaxPOrder, HRef, MaxHLevel);
    
    // If no exists any element to refine, the tolerance was reached
    if(!nelhrefs && !nelprefs)
        return true;
    return false;
}


bool ApplyingHPAdaptiveStrategyBasedOnUAndDU_III(TPZCompMesh *cmesh,TPZVec<STATE> &ErrorU,TPZVec<STATE> &ErrorDU,TPZVec<REAL> &Tol, int MaxPOrder, int MaxHLevel,std::ostream &out) {
    if(!cmesh) return false;
    int64_t iel, nelhrefs = 0, nelprefs = 0;
    int64_t nels = cmesh->NElements();
    
    TPZVec<int64_t> HRef(nels,0L), PRef(nels,0L);

    // Applying hp refinement only for elements with dimension as model dimension
    out << " Refinando malha com " << nels << " elementos e " << cmesh->NEquations() << " equacoes.\n";
    
    // Applying tolerance limits to define whether the element will be h-, p-, hp-refined or not. Implementation of the hp-adaptive table.
    // Note: Some elements can to have p and h refinements. But to indicate wheter the element must to refine twice h-ref, we have changed the index by -index
    for(iel=0L;iel<nels;iel++) {
        TPZCompEl *cel = cmesh->Element(iel);
        if(!cel || cel->Dimension() != cmesh->Dimension()) continue;
        
        if(ErrorU[iel] < Tol[0]) {
            if(ErrorDU[iel] > Tol[2])
                PRef[nelprefs++] = iel;
        }
        else if(ErrorU[iel] < Tol[1]) {
            PRef[nelprefs++] = iel;
        }
        else if(ErrorU[iel] < Tol[2]) {
            if(ErrorDU[iel] < Tol[2])
                PRef[nelprefs++] = iel;
            else
                HRef[nelhrefs++] = iel;
        }
        else {
            HRef[nelhrefs++] = iel;
            if(ErrorDU[iel] > Tol[2])
                PRef[nelprefs++] = iel;
        }
    }
    
    HRef.Resize(nelhrefs);
    PRef.Resize(nelprefs);
    
    // Doing h and p refinements
    ApplyHPRefinement(cmesh, PRef, MaxPOrder, HRef, MaxHLevel);
    
    // If no exists any element to refine, the tolerance was reached
    if (!nelhrefs && !nelprefs)
        return true;
    return false;
}
bool ApplyingHPAdaptiveStrategyBasedOnUAndDU_IV(TPZCompMesh *cmesh,TPZVec<STATE> &ErrorU,TPZVec<STATE> &ErrorDU,TPZVec<REAL> &Tol, int MaxPOrder, int MaxHLevel,std::ostream &out) {
    if(!cmesh) return false;
    int64_t iel, nelhrefs = 0, nelprefs = 0;
    int64_t nels = cmesh->NElements();
    
    TPZVec<int64_t> HRef(nels,0L), PRef(nels,0L);
    
    // Applying hp refinement only for elements with dimension as model dimension
    out << " Refinando malha com " << nels << " elementos e " << cmesh->NEquations() << " equacoes.\n";
    
    // Applying tolerance limits to define whether the element will be h-, p-, hp-refined or not. Implementation of the hp-adaptive table.
    // Note: Some elements can to have p and h refinements. But to indicate wheter the element must to refine twice h-ref, we have changed the index by -index
    for(iel=0L;iel<nels;iel++) {
        TPZCompEl *cel = cmesh->Element(iel);
        if(!cel || cel->Dimension() != cmesh->Dimension()) continue;
        
        if(ErrorU[iel] < Tol[0]) {
            if(ErrorDU[iel] > Tol[2])
                PRef[nelprefs++] = iel;
        }
        else if(ErrorU[iel] < Tol[1]) {
            if(ErrorDU[iel] < Tol[2])
                PRef[nelprefs++] = iel;
            else
                HRef[nelhrefs++] = iel;
        }
        else if(ErrorU[iel] < Tol[2]) {
            HRef[nelprefs++] = iel;
            if(ErrorDU[iel] > Tol[2])
                PRef[nelhrefs++] = iel;
        }
        else {
            if(ErrorDU[iel] < Tol[2]) {
                HRef[nelhrefs++] = iel;
                PRef[nelprefs++] = iel;
            }
            else
                HRef[nelhrefs++] = -1*iel;
        }
    }
    
    HRef.Resize(nelhrefs);
    PRef.Resize(nelprefs);
    
    // Doing h and p refinements
    ApplyHPRefinement(cmesh, PRef, MaxPOrder, HRef, MaxHLevel);
    
    // If no exists any element to refine, the tolerance was reached
    if (!nelhrefs && !nelprefs)
        return true;
    return false;
}

bool ApplyingHPAdaptiveStrategyBasedOnUAndDU_V(TPZCompMesh *cmesh,TPZVec<STATE> &ErrorU,TPZVec<STATE> &ErrorDU,TPZVec<REAL> &Tol, int MaxPOrder, int MaxHLevel,std::ostream &out) {
	if(!cmesh) return false;
    int64_t iel, nelhrefs = 0, nelprefs = 0;
	int64_t nels = cmesh->NElements();
    
    TPZVec<int64_t> HRef(nels,0L), PRef(nels,0L);

	// Applying hp refinement only for elements with dimension as model dimension
	out << " Refinando malha com " << nels << " elementos e " << cmesh->NEquations() << " equacoes.\n";

    // Applying tolerance limits to define whether the element will be h-, p-, hp-refined or not. Implementation of the hp-adaptive table.
    // Note: Some elements can to have p and h refinements. But to indicate wheter the element must to refine twice h-ref, we have changed the index by -index
    for(iel=0L;iel<nels;iel++) {
        TPZCompEl *cel = cmesh->Element(iel);
        if(!cel || cel->Dimension() != cmesh->Dimension()) continue;
        if(ErrorU[iel] < Tol[0]) {
            if(ErrorDU[iel] > Tol[1])
                PRef[nelprefs++] = iel;
        }
        else if(ErrorU[iel] < Tol[1]) {
            if(ErrorDU[iel] < Tol[2])
                PRef[nelprefs++] = iel;
            else
                HRef[nelhrefs++] = iel;
        }
        else if(ErrorU[iel] < Tol[2]) {
            if(ErrorDU[iel] < Tol[1])
                PRef[nelprefs++] = iel;
            else {
                HRef[nelhrefs++] = iel;
                if(ErrorDU[iel] > Tol[2])
                    PRef[nelprefs++] = iel;
            }
        }
        else {
            if(ErrorDU[iel] > Tol[2])
                HRef[nelhrefs] = -1*iel;
            else if(ErrorDU[iel] > Tol[1]) {
                HRef[nelhrefs++] = iel;
                PRef[nelprefs++] = iel;
            }
        }
    }
    
	HRef.Resize(nelhrefs);
	PRef.Resize(nelprefs);

	// Doing h and p refinements
	ApplyHPRefinement(cmesh, PRef, MaxPOrder, HRef, MaxHLevel);

	// If no exists any element to refine, the tolerance was reached
	if (!nelhrefs && !nelprefs)
		return true;
	return false;
}
bool ApplyingHPAdaptiveStrategyBasedOnUAndDU_VI(TPZCompMesh *cmesh,TPZVec<STATE> &ErrorU,TPZVec<STATE> &ErrorDU,TPZVec<REAL> &Tol, int MaxPOrder, int MaxHLevel,std::ostream &out) {
    if(!cmesh) return false;
    int64_t iel, nelhrefs = 0, nelprefs = 0;
    int64_t nels = cmesh->NElements();
    
    TPZVec<int64_t> HRef(nels,0L), PRef(nels,0L);
    
    // Applying hp refinement only for elements with dimension as model dimension
    out << " Refinando malha com " << nels << " elementos e " << cmesh->NEquations() << " equacoes.\n";
    
    // Applying tolerance limits to define whether the element will be h-, p-, hp-refined or not. Implementation of the hp-adaptive table.
    // Note: Some elements can to have p and h refinements. But to indicate wheter the element must to refine twice h-ref, we have changed the index by -index
    for(iel=0L;iel<nels;iel++) {
        TPZCompEl *cel = cmesh->Element(iel);
        if(!cel || cel->Dimension() != cmesh->Dimension()) continue;
        if(ErrorU[iel] < Tol[0]) {
            if(ErrorDU[iel] > Tol[2])
                HRef[nelprefs++] = iel;
            else if(ErrorDU[iel] > Tol[1])
                PRef[nelhrefs++] = iel;
        }
        else if(ErrorU[iel] < Tol[1]) {
            if(ErrorDU[iel] > Tol[1]) {
                HRef[nelprefs++] = iel;
                if(ErrorDU[iel] > Tol[2])
                    PRef[nelhrefs++] = iel;
            }
            else
                PRef[nelhrefs++] = iel;
        }
        else if(ErrorU[iel] < Tol[2]) {
            HRef[nelhrefs++] = iel;
            if(ErrorDU[iel] > Tol[2])
                HRef[nelprefs++] *= -1;
            else if(ErrorDU[iel] > Tol[1])
                PRef[nelhrefs++] = iel;
        }
        else {
            if(ErrorDU[iel] > Tol[1])
                HRef[nelhrefs] = -1*iel;
            else {
                PRef[nelprefs++] = iel;
                HRef[nelhrefs++] = iel;
            }
        }
    }
    
    HRef.Resize(nelhrefs);
    PRef.Resize(nelprefs);
    
    // Doing h and p refinements
    ApplyHPRefinement(cmesh, PRef, MaxPOrder, HRef, MaxHLevel);
    
    // If no exists any element to refine, the tolerance was reached
    if (!nelhrefs && !nelprefs)
        return true;
    return false;
}

bool ApplyingHPAdaptiveStrategyBasedOnUAndDU_XI(TPZCompMesh *cmesh,TPZVec<STATE> &ErrorU,TPZVec<STATE> &ErrorDU,TPZVec<REAL> &Tol, int MaxPOrder, int MaxHLevel,std::ostream &out) {
    if(!cmesh) return false;
    int64_t iel, nelhrefs = 0, nelprefs = 0;
    int64_t nels = cmesh->NElements();
    
    TPZVec<int64_t> HRef(nels,0L), PRef(nels,0L);
    
    // Applying hp refinement only for elements with dimension as model dimension
    out << " Refinando malha com " << nels << " elementos e " << cmesh->NEquations() << " equacoes.\n";
    
    // Applying tolerance limits to define whether the element will be h-, p-, hp-refined or not. Implementation of the hp-adaptive table.
    // Note: Some elements can to have p and h refinements. But to indicate wheter the element must to refine twice h-ref, we have changed the index by -index
    for(iel=0L;iel<nels;iel++) {
        TPZCompEl *cel = cmesh->Element(iel);
        if(!cel || cel->Dimension() != cmesh->Dimension()) continue;
        
        if(ErrorU[iel] < Tol[0]) {
            if(ErrorDU[iel] > Tol[2])
                PRef[nelprefs++] = iel;
        }
        else if(ErrorU[iel] < Tol[1]) {
            if(ErrorDU[iel] > Tol[2])
                HRef[nelhrefs++] = iel;
            PRef[nelprefs++] = iel;
        }
        else if(ErrorU[iel] < Tol[2]) {
            if(ErrorDU[iel] < Tol[2])
                PRef[nelprefs++] = iel;
            HRef[nelhrefs++] = iel;
        }
        else {
            if(ErrorDU[iel] > Tol[2])
                HRef[nelprefs++] = -1*(iel);
            else
                HRef[nelhrefs++] = iel;
        }
    }
    
    HRef.Resize(nelhrefs);
    PRef.Resize(nelprefs);
    
    // Doing h and p refinements
    ApplyHPRefinement(cmesh, PRef, MaxPOrder, HRef, MaxHLevel);
    
    // If no exists any element to refine, the tolerance was reached
    if (!nelhrefs && !nelprefs)
        return true;
    return false;
}

void PrintNRefinementsByType(int64_t nels,int64_t newnels,int64_t hrefcounter,int64_t prefcounter,std::ostream &out) {
    out << "\n HP Refinement done, on  " << nels << " elements, given " << newnels << " elements. "<< std::endl;
    out << " Refinement type H " << hrefcounter << " elements." << std::endl;
    out << " Refinement type P " << prefcounter << " elements." << std::endl;
}


