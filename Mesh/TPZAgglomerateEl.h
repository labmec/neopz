
#ifndef AGGLOMERATEELEMHPP
#define AGGLOMERATEELEMHPP

#include "pzreal.h"
#include "pzstack.h"
//#include "pzeltype.h"
#include "TPZCompElDisc.h"
#include <iostream>

class TPZGeoEl;
class TPZCompEl;
class TPZCompMesh;
class TPZGeoElSide;
class TPZInterfaceElement;
struct TPZElementMatrix;
using namespace std;

/**
 * TPZAgglomerateElement clase that it manages the generation of elements 
 * from the agglomeration of some geometric elements
 */

class TPZAgglomerateElement : public TPZCompElDisc { 

private:

  /**
   * indexes na malha fina dos sub-elementos computacionais aglomerados pelo atual
   */
  TPZStack<int> fIndexes;

  /**
   * malha da qual o elemento atual foi gerado
   */
  TPZCompMesh *fMotherMesh;

public:

  /** Constructor: caso o elemento é passível de ser agrupado retorna-se 
   * um novo index caso contrário index = -1
   */
  TPZAgglomerateElement(int &index,TPZCompMesh &cmesh,TPZCompMesh *finemesh);

  /** inicializa os dados caracteristicos do elemento aglomerado */
  void InitializeElement(int mat);

  /** adiciona index do sub-elemento*/
  static void AddSubElementIndex(TPZCompMesh *cmesh,int subel,int destind);

  /**Destructor do objeto*/
  ~TPZAgglomerateElement(){};

  /**
   * Type of the element 
   */  
  MElementType Type(){return ENoType;}

  /** retorna malha mae */
  TPZCompMesh *MotherMesh(){return fMotherMesh;}

  /** acumula a lista de regras de integra¢ão no elemento deformado*/
  virtual void AccumulateIntegrationRule(int degree, TPZStack<REAL> &point, TPZStack<REAL> &weight);

  /**calcula o ponto centro de massa dos elementos aglomerados */
  void CenterPoint();

  /** retorna o volume do elemento geométrico referenciado */
  REAL VolumeOfEl();

  /**
   * Calcula a restrição da solução ou resíduo dos elementos aglomerados para o
   * obtido por aglomeração - este último chamado de elemento pai
   */
  void CalcResidual(TPZFMatrix &Rhs,TPZCompElDisc *el);

  /**
   * Monta a equação diferencial do modelo sobre o elemento definido por 
   * aglomeração de elementos da malha fina
   */
  void CalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef);

  /** retorna o número de sub-elementos aglomerados */
  int NIndexes(){return fIndexes.NElements();}

  /**
   * retorna o elemento computacional da malha fina de índice index
   */
  TPZCompEl *FineElement(int index);

  /**os geométricos agrupados apontam para o computacional*/
  void SetReference();

   /** 
    * retorna o sub-elemento número sub 
    */
   TPZCompEl *SubElement(int sub);

   REAL NormalizeConst();

  static TPZCompMesh *CreateAgglomerateMesh(TPZCompMesh *finemesh,TPZVec<int> &accumlist,int numaggl);

};
#endif
