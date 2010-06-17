#ifndef VISUAL_MATRIX
#define VISUAL_MATRIX
#include "pzfmatrix.h"
#include "pzvec.h"
#include "fstream"


void VisualMatrix(TPZFMatrix &matrix, const std::string &outfilename);

void VisualMatrixDX(TPZFMatrix &matrix, const std::string &outfilename);

void VisualMatrixVTK(TPZFMatrix &matrix, const std::string &outfilename);

#endif
