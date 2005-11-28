// -*- c++ -*-

//$Id: postprocess.h,v 1.1 2005-11-28 13:49:47 tiago Exp $

class TPZCompMesh;

void ComputeNormalFlux(std::ofstream & file, TPZCompMesh &cmesh, const int matid/*, const int npoints*/);
void ComputeGradient(TPZInterpolatedElement * cel, TPZVec<REAL> &intpoint, TPZVec<REAL> &gradU);

void ComputeSolution(std::ofstream &file, TPZCompMesh &cmesh, const int matid/*, const int npoints*/);
void ComputePointSolution(TPZInterpolatedElement * cel, TPZVec<REAL> &intpoint, REAL &solution);
