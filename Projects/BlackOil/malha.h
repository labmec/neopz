//$Id: malha.h,v 1.1 2008-11-12 12:47:09 fortiago Exp $

#include "pzcmesh.h"

/** Gera malha tridimensional para problema com solucao unidimensional.
  * O problema eh injecao de agua em reservatorio saturado de oleo
  */
TPZCompMesh *Unidimensional(int h, double deltaT);

