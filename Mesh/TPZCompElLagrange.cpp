//
//  TPZCompElLagrange.cpp
//  PZ
//
//  Created by Philippe Devloo on 11/2/13.
//
//

#include "TPZCompElLagrange.h"
#include "pzmaterial.h"
#include "pzelmat.h"

TPZCompElLagrange::~TPZCompElLagrange()
{
    
}

/**
 * @brief Method for creating a copy of the element in a patch mesh
 * @param mesh Patch clone mesh
 * @param gl2lcConMap map the connects indexes from global element (original) to the local copy.
 * @param gl2lcElMap map the computational elements
 */
/**
 * Otherwise of the previous clone function, this method don't
 * copy entire mesh. Therefore it needs to map the connect index
 * from the both meshes - original and patch
 */
TPZCompEl *TPZCompElLagrange::ClonePatchEl(TPZCompMesh &mesh,
                                std::map<long,long> & gl2lcConMap,
                                std::map<long,long> & gl2lcElMap) const
{
    TPZCompElLagrange *newel = new TPZCompElLagrange(mesh,*this,gl2lcElMap);
    for (int i=0; i<2; i++) {
        newel->fIdf[i] = fIdf[i];
        std::map<long,long>::iterator it = gl2lcConMap.find(fConnect[i]);
        if (it != gl2lcConMap.end()) {
            newel->fConnect[i] = it->second;
        }
        else
        {
            DebugStop();
        }
    }
    return newel;
}

/**
 * @brief Computes the element stifness matrix and right hand side
 * @param ek element stiffness matrix
 * @param ef element load vector
 */
void TPZCompElLagrange::CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef)
{
    InitializeElementMatrix(ek, ef);
    TPZConnect &c = Connect(0);
    int blsize = c.NShape()*c.NState();
    ek.fMat(fIdf[0],fIdf[0]) = 1.;
    ek.fMat(fIdf[0],blsize+fIdf[1]) = -1.;
    ek.fMat(blsize+fIdf[1],fIdf[0]) = -1.;
    ek.fMat(blsize+fIdf[1],blsize+fIdf[1]) = 1.;
}

/**
 * @brief Computes the element right hand side
 * @param ef element load vector(s)
 */
//void TPZCompElLagrange::CalcResidual(TPZElementMatrix &ef)
//{
//    
//}

void TPZCompElLagrange::InitializeElementMatrix(TPZElementMatrix &ek, TPZElementMatrix &ef){
    int numloadcases = 1;
	int numdof = 1;
    TPZMaterial *mat = this->Material();
    if (mat)
    {
        mat->NumLoadCases();
        mat->NStateVariables();
    }
	const int ncon = this->NConnects();
    
    ek.fMesh = Mesh();
    ek.fType = TPZElementMatrix::EK;
    ef.fMesh = Mesh();
    ef.fType = TPZElementMatrix::EF;
    
	ek.fBlock.SetNBlocks(ncon);
	ef.fBlock.SetNBlocks(ncon);
	ek.fNumStateVars = numdof;
	ef.fNumStateVars = numdof;
	int i;
    int numeq=0;
	for(i=0; i<ncon; i++){
        TPZConnect &c = Connect(i);
        int nshape = c.NShape();
        int nstate = c.NState();
        
		ek.fBlock.Set(i,nshape*nstate);
		ef.fBlock.Set(i,nshape*nstate);
        numeq += nshape*nstate;
	}
	ek.fMat.Redim(numeq,numeq);
	ef.fMat.Redim(numeq,numloadcases);
	ek.fConnect.Resize(ncon);
	ef.fConnect.Resize(ncon);
	for(i=0; i<ncon; i++){
		(ef.fConnect)[i] = ConnectIndex(i);
		(ek.fConnect)[i] = ConnectIndex(i);
	}
}//void



