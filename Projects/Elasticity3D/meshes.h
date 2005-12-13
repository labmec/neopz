// -*- c++ -*-

//$Id: meshes.h,v 1.1 2005-12-13 11:47:28 tiago Exp $

class TPZCompMesh;

/** Barra tracionada
 */
TPZCompMesh * BarraTracionada(int h, int p); 
 
/** Viga engastada
 */
TPZCompMesh * VigaEngastada(int h, int p);
 
/** Viga bi-engastada
 */
TPZCompMesh * VigaBiEngastada(int h, int p);
 
/** Vigas em balanćo
 */
TPZCompMesh * VigasBalanco(int h, int p);
 
/** Bloco confinado
 */
TPZCompMesh * BlocoConfinado(int h, int p);

void SetPOrder(int p);
