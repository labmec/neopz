//$Id: malha.h,v 1.2 2008-11-25 13:28:15 fortiago Exp $

#include "pzcmesh.h"

/** Gera malha tridimensional para problema com solucao unidimensional.
  * O problema eh injecao de agua em reservatorio saturado de oleo
  */
TPZCompMesh *Unidimensional(int h, double deltaT);

/** Gera malha tridimensional para problema com solucao unidimensional.
  * O problema eh agua em cima e oleo em baixo que devem migrar.
  */
TPZCompMesh *UnidimensionalGravidade(int h, double deltaT);
