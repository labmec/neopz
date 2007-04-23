//$Id: TPZAgglomerateEl.h,v 1.27 2007-04-23 19:02:44 tiago Exp $
#ifndef AGGLOMERATEELEMHPP
#define AGGLOMERATEELEMHPP

#include "pzreal.h"
#include "pzstack.h"
#include "pzeltype.h"
#include "TPZCompElDisc.h"
#include <iostream>
#include <map>
#include <set>

#include <stdlib.h>

class TPZGeoEl;
class TPZCompEl;
class TPZCompMesh;
class TPZGeoElSide;
class TPZInterfaceElement;
struct TPZElementMatrix;
class TPZAgglomerateMesh;

/// Implements an agglomerated discontinuous element
/**
 * TPZAgglomerateElement clase that it manages the generation of elements 
 * from the agglomeration of some geometric elements
 @ingroup CompElement
 */
class TPZAgglomerateElement : public TPZCompElDisc { 

private:

  /**
   * indexes na malha fina dos sub-elementos computacionais aglomerados pelo atual
   */
  TPZStack<int> fIndexes;

  /**
   * malha à qual os elementos aglomerados fazem parte e 
   * do qual o atual foi obtido
   */
  TPZCompMesh *fMotherMesh;

  /**
   * Stores the element's inner radius.
   * It is the lessest distance between the element center point and its interface's center points.
   */
  REAL fInnerRadius;

  /**
   * Stores the number of interfaces of the element.
   */
  int fNFaces;
  
  /**
   * Material id of the agglomerated element
   */
  int fMaterialId;

public:

  /** Constructor: caso o elemento é passível de ser agrupado retorna-se 
   * um novo index caso contrário index = -1
   */
  TPZAgglomerateElement(int nummat,int &index,TPZCompMesh &aggcmesh,TPZCompMesh *finemesh);
  
  TPZAgglomerateElement();

  /** inicializa os dados caracteristicos do elemento aglomerado */
  void InitializeElement();

  /**
   * Set the inner radius value.
   */
  void SetInnerRadius(REAL InnerRadius) { fInnerRadius = InnerRadius;}

  /**
   * Returns the inner radius value.
   * Inner radius mut be set with SetInnerRadius.
   */
  REAL InnerRadius2() {return fInnerRadius;}

  /**
   * Returns the inner radius value.
   * Inner radius is the sub-element's radius average weighted by their volumes.
   */
  REAL InnerRadius(){
    int nsubel = this->NIndexes();
    REAL value = 0.;
    TPZCompElDisc * disc;
    for (int i = 0; i < nsubel; i++){
      disc = dynamic_cast<TPZCompElDisc *>(this->SubElement(i) );
#ifdef DEBUG
      if (!disc) {
	PZError << "TPZAgglomerateElement::InnerRadius FineElement must be a TPZCompElDisc" << std::endl;
	exit (-1);
      }
#endif
      value += disc->InnerRadius() * disc->VolumeOfEl();
    }
    value = value / this->VolumeOfEl();
    return value;  
  }

  /**
   * Set element's number of interfaces.
   */
  void SetNInterfaces(int nfaces) {fNFaces = nfaces; }

  /**
   * Retunrs the number of interfaces;
   */
  int NInterfaces() {return fNFaces;}

  /** adiciona index do sub-elemento*/
  static void AddSubElementIndex(TPZCompMesh *aggcmesh,int subel,int destind);

  /**Destructor do objeto*/
  ~TPZAgglomerateElement(){};

  /**
   * Type of the element 
   */  
  MElementType Type(){return EAgglomerate;}

  /** retorna malha mae */
  TPZCompMesh *MotherMesh(){return fMotherMesh;}

  /** acumula a lista de regras de integra¢ão no elemento deformado*/
  virtual void AccumulateIntegrationRule(int degree, TPZStack<REAL> &point, TPZStack<REAL> &weight);

  /** accumulate the vertices of the agglomerated elements */
  virtual void AccumulateVertices(TPZStack<TPZGeoNode *> &nodes);

  /** calcula o ponto centro de massa dos elementos aglomerados */
  void CenterPoint();

  /** devolve o centro de massa do elemento */
  virtual void CenterPoint(TPZVec<REAL> &center);

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
  int NIndexes() const {return fIndexes.NElements();}

  /**
   * retorna o elemento computacional da malha fina de índice index
   * Desnecessario porque identico a SubElement
  
  TPZCompEl *FineElement(int index);
  */

  /**return the geometric element to which this element references*/
  TPZGeoEl *CalculateReference();

//  /**os geométricos agrupados apontam para o computacional*/
/*  void SetReference(){
    int nindex = NIndexes(),i;
    for(i=0;i<nindex;i++){
      TPZCompEl *cel = SubElement(i);
      int type = cel->Type();
      //caso comp é aglomerado: chamada recursiva
      if(type == EAgglomerate){//aglomerado
	SetReference();
      } else if(type == EDiscontinuous){//descontínuo
	//o geométrico agrupado apontará para o atual computacional
	cel->Reference()->SetReference(this);
      }
    }
  }

  void SetReference2(int sub){

  TPZCompEl *cel = SubElement(sub);
  int type = cel->Type();
  //caso comp é aglomerado: chamada recursiva
  if(type == EAgglomerate){//aglomerado
    dynamic_cast<TPZAgglomerateElement *>(cel)->SetReference();
  } else if(type == EDiscontinuous){//descontínuo
    //o geométrico agrupado apontará para o atual computacional
    cel->Reference()->SetReference(this);
  }
} */

  /** \brief Computes a measure of the element
   * Computes the maximum distance in x,y and z and then returns the minimum of
   * these three measures
   */
virtual REAL LesserEdgeOfEl();


   /** 
    * retorna o sub-elemento número sub 
    */
   TPZCompEl *SubElement(int sub) const;

   REAL NormalizeConst();

  /**
   * it creates new conect that it associates the degrees of freedom of the
   * element and returns its index 
   */
  virtual int CreateMidSideConnect();

  /**
   * it returns dimension from the elements
   */
  int Dimension() const;

  /**
   * it prints the features of the element 
   */
  virtual void Print(std::ostream & out = std::cout);

  /**
   * create copy of the materials tree
   */
//  void CreateMaterialCopy(TPZCompMesh &aggcmesh);

  void ListOfDiscEl(TPZStack<TPZCompEl *> &elvec);

  void IndexesDiscSubEls(TPZStack<int> &elvec);

  /**
   * Returns the number of sides. If all the volumes agglomerated have the same number, it returns this number, else it returns -1.
   */
  int NSides();

  void CreateGraphicalElement(TPZGraphMesh &grmesh, int dimension);

//  void Solution(TPZVec<REAL> &qsi,int var,TPZManVector<REAL> &sol);

//  TPZGeoEl *FatherN(TPZGeoEl *sub,int n);

//  int NSubsOfLevels(TPZGeoEl *father,int nlevels);

//  int NSubCompEl(TPZGeoEl *father);

  static void ListOfGroupings(TPZCompMesh *finemesh,TPZVec<int> &accumlist,int nivel,int &numaggl,int dim);

//  void FineSolution(TPZVec<REAL> &x,TPZCompElDisc *disc,TPZVec<REAL> &uh);

  void Print(TPZStack<int> &listindex);

  void ProjectSolution(TPZFMatrix &projectsol);
//  void ProjectSolution2(TPZFMatrix &projectsol);


  static TPZAgglomerateMesh *CreateAgglomerateMesh(TPZCompMesh *finemesh,TPZVec<int> &accumlist,int numaggl);

  static void ComputeNeighbours(TPZCompMesh *mesh, std::map<TPZCompElDisc *,std::set<TPZCompElDisc *> > &neighbours);
  
  /**
  * returns the unique identifier for reading/writing objects to streams
  */
  virtual int ClassId() const;
  /**
  Save the element data to a stream
  */
  virtual void Write(TPZStream &buf, int withclassid);
  
  /**
  Read the element data from a stream
  */
  virtual void Read(TPZStream &buf, void *context);



};
#endif
