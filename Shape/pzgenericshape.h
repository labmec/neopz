#ifndef PZGENERICSHAPE_H
#define PZGENERICSHAPE_H

#include "pzreal.h"
#include "pzvec.h"
#include "pztrnsform.h"

template<class T>
class TPZFMatrix;


template <class TSHAPE>
TPZTransform<REAL> GetSideTransform(const int side, int trans_id);
template <class TSHAPE>
void ComputeTransforms(const TPZVec<int64_t> &id,  TPZVec<TPZTransform<REAL> > &transvec) ;

template <class TSHAPE>
void Shape(TPZVec<REAL> &pt, TPZVec<int> orders,TPZVec<TPZTransform<REAL> > &transvec, TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) ;


struct TParDefs
{
    TPZManVector<int,27> orders;
    TPZManVector<int,27> nshape;
    TPZVec<TPZTransform<REAL> > transvec;
};

#endif
