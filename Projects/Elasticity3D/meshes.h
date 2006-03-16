// -*- c++ -*-

//$Id: meshes.h,v 1.4 2006-03-16 01:54:52 tiago Exp $

class TPZCompMesh;
class TPZFMatrix;
#include "pzvec.h"

/** Barra tracionada - Validado
 */
TPZCompMesh * BarraTracionada(int h, int p); 

/** Barra tracionada sem condicao Dirichlet - Validado
 */
TPZCompMesh * BarraTracionadaNeumann(int h, int p); 

/** Barra tracionada girada de 10 graus em torno do eixo x.
 * int cubo = 1 faz um cubo - Validado
 * int cubo = 2 faz um prisma - Validado
 */
TPZCompMesh * BarraTracionadaGirada(int h, int p);
 
/** Viga engastada com momento concentrado na extremidade livre - Validado
 */
TPZCompMesh * VigaEngastada(int h, int p);
void MomentoExtremidade(TPZVec<REAL> &pto, TPZVec<REAL> &force);

/** Viga engastada com forca de corpo (volumetrica) - Validado
 */
TPZCompMesh * VigaEngastadaForcaVolume(int h, int p);

void Rotate(TPZFMatrix & pt);

void SetPOrder(int p);
