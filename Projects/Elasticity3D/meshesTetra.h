// -*- c++ -*-

//$Id: meshesTetra.h,v 1.1 2006-03-16 13:25:38 tiago Exp $

class TPZCompMesh;
#include "pzvec.h"

/** Viga engastada com momento concentrado na extremidade livre - Validado
 */
TPZCompMesh * VigaEngastadaTetra(int h, int p);
void MomentoExtremidadeTetra(TPZVec<REAL> &pto, TPZVec<REAL> &force);

/** Viga engastada com forca de corpo (volumetrica) - Validado
 */
TPZCompMesh * VigaEngastadaForcaVolumeTetra(int h, int p);

void SetInterpOrder(int p);
