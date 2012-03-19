/**
 * @file
 * @brief Contains declaration of
 */
//$Id: TPZAgglomerateEl.h,v 1.29 2011-05-11 02:48:49 phil Exp $
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

/**
 * @brief Implements an agglomerated discontinuous element. \ref CompElement "Computational Element"
 * @ingroup CompElement
 */
/** TPZAgglomerateElement clase that it manages the generation of elements 
 * from the agglomeration of some geometric elements
 */
class TPZAgglomerateElement : public TPZCompElDisc { 
	
private:
	
	/**
	 * @brief Indexes into the fine mesh of computational subelements in the clusters by the current
	 */
	/** 
	 * Indexes na malha fina dos sub-elementos computacionais aglomerados pelo atual
	 */
	TPZStack<int> fIndexes;
	
	/**
	 * @brief Mesh for the clusters which elements are part and from that the current mesh was obtained
	 */
	/** Malha a qual os elementos aglomerados fazem parte e 
	 * do qual o atual foi obtido
	 */
	TPZCompMesh *fMotherMesh;
	
	/**
	 * @brief Stores the element's inner radius.
	 */ 
	/** It is the lessest distance between the element center point and its interface's center points.
	 */
	REAL fInnerRadius;
	
	/** @brief Stores the number of interfaces of the element. */
	int fNFaces;
	
	/** @brief Material id of the agglomerated element */
	int fMaterialId;
	
public:
	
	/** @brief Constructor: If the element is possible to grouped returns a new index, else returns -1. */
	TPZAgglomerateElement(int nummat,int &index,TPZCompMesh &aggcmesh,TPZCompMesh *finemesh);
	
	TPZAgglomerateElement();
	
	/** @brief Initialize the characteristics data of the clustered elements.  
	 */
	void InitializeElement();
	
	/** @brief Sets the inner radius value. */
	void SetInnerRadius(REAL InnerRadius) { fInnerRadius = InnerRadius;}
	
	/** @brief Returns the inner radius value. */
	/** Inner radius mut be set with SetInnerRadius. */
	REAL InnerRadius2() {return fInnerRadius;}
	
	/** @brief Returns the inner radius value. */ 
	/** Inner radius is the sub-element's radius average weighted by their volumes. */
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
	
	/** @brief Sets element's number of interfaces. */
	void SetNInterfaces(int nfaces) {fNFaces = nfaces; }
	
	/** @brief Retunrs the number of interfaces. */
	int NInterfaces() {return fNFaces;}
	
	/** @brief Insert the subelement index. */
	static void AddSubElementIndex(TPZCompMesh *aggcmesh,int subel,int destind);
	
	/** @brief Destructor*/
	~TPZAgglomerateElement(){};
	
	/** @brief Type of the element */  
	MElementType Type(){return EAgglomerate;}
	
	/** @brief Returns father mesh. */
	TPZCompMesh *MotherMesh(){return fMotherMesh;}
	
	/** @brief Accumulates integration rule to deformed element. */
	virtual void AccumulateIntegrationRule(int degree, TPZStack<REAL> &point, TPZStack<REAL> &weight);
	
	/** @brief Accumulate the vertices of the agglomerated elements */
	virtual void AccumulateVertices(TPZStack<TPZGeoNode *> &nodes);
	
	/** @brief Computes the center of the mass to clustered elements */
	void CenterPoint();
	
	/** @brief Returns the center of the mass */
	virtual void CenterPoint(TPZVec<REAL> &center);
	
	/** @brief Returns the volume of the geometric element referenced */
	REAL VolumeOfEl();
	
	/** @brief Computes the residual of the solution to father element from clustered subelements. */
	void CalcResidual(TPZFMatrix<REAL> &Rhs,TPZCompElDisc *el);
	
	void CalcResidual(TPZElementMatrix &ef)
	{
		std::cout << __PRETTY_FUNCTION__ << " is not implemented\n";
		exit(-1);
	}
	
	/** @brief Assembles the differential equation to model over the element defined by clustered subelements. */
	void CalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef);
	
	/** @brief Returns the number of clustered subelements. */
	int NIndexes() const {return fIndexes.NElements();}
	
	/** @brief Returns the geometric element to which this element references*/
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
	
	/** @brief Computes a measure of the element */ 
	/** Computes the maximum distance in x,y and z and then returns the minimum of these three measures */
	virtual REAL LesserEdgeOfEl();
	
	/** @brief Returns the "sub" subelement. */
	TPZCompEl *SubElement(int sub) const;
	
	REAL NormalizeConst();
	
	/** @brief It creates new conect that it associates the degrees of freedom of the element and returns its index */
	virtual int CreateMidSideConnect();
	
	/** @brief It returns dimension from the elements */
	int Dimension() const;
	
	/** @brief Prints the features of the element */
	virtual void Print(std::ostream & out = std::cout) const;

	/** @brief Returns a vector of all discontinuous elements in cluster. */
	void ListOfDiscEl(TPZStack<TPZCompEl *> &elvec);
	/** @brief Returns a vector of all indexes of the discontinuous elements in cluster. */	
	void IndexesDiscSubEls(TPZStack<int> &elvec);
	
	/** @brief Returns the number of sides. If all the volumes agglomerated have the same number, it returns this number, else it returns -1. */
	int NSides();
	
	/** Creates graphical element to postprocessing */
	void CreateGraphicalElement(TPZGraphMesh &grmesh, int dimension);
		
	static void ListOfGroupings(TPZCompMesh *finemesh,TPZVec<int> &accumlist,int nivel,int &numaggl,int dim);
	
	//  void FineSolution(TPZVec<REAL> &x,TPZCompElDisc *disc,TPZVec<REAL> &uh);
	
	void Print(TPZStack<int> &listindex);
	
	void ProjectSolution(TPZFMatrix<REAL> &projectsol);
	
	
	static TPZAgglomerateMesh *CreateAgglomerateMesh(TPZCompMesh *finemesh,TPZVec<int> &accumlist,int numaggl);
	
	static void ComputeNeighbours(TPZCompMesh *mesh, std::map<TPZCompElDisc *,std::set<TPZCompElDisc *> > &neighbours);
	
	/** @brief Returns the unique identifier for reading/writing objects to streams */
	virtual int ClassId() const;
	/*@brief Save the element data to a stream */
	virtual void Write(TPZStream &buf, int withclassid);
	
	/** @brief Read the element data from a stream */
	virtual void Read(TPZStream &buf, void *context);

};

#endif
