/**
 * @file
 * @brief Contains declaration of
 */

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
	TPZStack<int64_t> fIndexes;
	
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
	/** It is the smallest distance between the element center point and its interface's center points.
	 */
	REAL fInnerRadius;
	
	/** @brief Stores the number of interfaces of the element. */
	int fNFaces;
	
	/** @brief Material id of the agglomerated element */
	int fMaterialId;

	template<class TVar>
    void CalcStiffT(TPZElementMatrixT<TVar> &ek, TPZElementMatrixT<TVar> &ef);
public:
	
	/** @brief Constructor: If the element is possible to grouped returns a new index, else returns -1. */
	TPZAgglomerateElement(int nummat,TPZCompMesh &aggcmesh,TPZCompMesh *finemesh);
	
	TPZAgglomerateElement();
	
	/** @brief Initialize the characteristics data of the clustered elements.  
	 */
	void InitializeElement();
	
	/** @brief Sets the inner radius value. */
	void SetInnerRadius(REAL InnerRadius) override { fInnerRadius = InnerRadius;}
	
	/** @brief Returns the inner radius value. */
	/** Inner radius mut be set with SetInnerRadius. */
	REAL InnerRadius2() {return fInnerRadius;}
	
	/** @brief Returns the inner radius value. */ 
	/** Inner radius is the sub-element's radius average weighted by their volumes. */
	REAL InnerRadius() override {
		int64_t nsubel = this->NIndexes();
		REAL value = 0.;
		TPZCompElDisc * disc;
		for (int64_t i = 0; i < nsubel; i++){
			disc = dynamic_cast<TPZCompElDisc *>(this->SubElement(i) );
#ifdef PZDEBUG
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
	void SetNInterfaces(int nfaces)  override {fNFaces = nfaces; }
	
	/** @brief Retunrs the number of interfaces. */
	int NInterfaces() override {return fNFaces;}
	
	/** @brief Insert the subelement index. */
	static void AddSubElementIndex(TPZCompMesh *aggcmesh,int64_t subel,int64_t destind);
	
	/** @brief Destructor*/
	~TPZAgglomerateElement(){};
	
	/** @brief Type of the element */  
	MElementType Type() override {return EAgglomerate;}
	
	/** @brief Returns father mesh. */
	TPZCompMesh *MotherMesh(){return fMotherMesh;}
	
	/** @brief Accumulates integration rule to deformed element. */
	virtual void AccumulateIntegrationRule(int degree, TPZStack<REAL> &point, TPZStack<REAL> &weight) override;
	
	/** @brief Accumulate the vertices of the agglomerated elements */
	virtual void AccumulateVertices(TPZStack<TPZGeoNode *> &nodes) override;
	
	/** @brief Computes the center of the mass to clustered elements */
	void CenterPoint();
	
	/** @brief Returns the center of the mass */
	virtual void CenterPoint(TPZVec<REAL> &center) override;
	
	/** @brief Returns the volume of the geometric element referenced */
	REAL VolumeOfEl() override;
	
	/** @brief Computes the residual of the solution to father element from clustered subelements. */
	void CalcResidual(TPZFMatrix<REAL> &Rhs,TPZCompElDisc *el);
	
	void CalcResidual(TPZElementMatrixT<STATE> &ef) override
	{
		std::cout << __PRETTY_FUNCTION__ << " is not implemented\n";
		exit(-1);
	}
	
	/** @brief Assembles the differential equation to model over the element defined by clustered subelements. */
	void CalcStiff(TPZElementMatrixT<STATE> &ek,
				   TPZElementMatrixT<STATE> &ef) override{
		CalcStiffT(ek,ef);
	}
	
	/** @brief Returns the number of clustered subelements. */
	int64_t NIndexes() const { return fIndexes.NElements(); }
	
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
	TPZCompEl *SubElement(int64_t sub) const;
	
	REAL NormalizeConst() override;
	
	/** @brief It creates new conect that it associates the degrees of freedom of the element and returns its index */
	virtual int64_t CreateMidSideConnect() override;
	
	/** @brief It returns dimension from the elements */
	int Dimension() const override;
	
	/** @brief Prints the features of the element */
	virtual void Print(std::ostream & out = std::cout) const override;

	/** @brief Returns a vector of all discontinuous elements in cluster. */
	void ListOfDiscEl(TPZStack<TPZCompEl *> &elvec);
	/** @brief Returns a vector of all indexes of the discontinuous elements in cluster. */	
	void IndexesDiscSubEls(TPZStack<int64_t> &elvec);
	
	/** @brief Returns the number of sides. If all the volumes agglomerated have the same number, it returns this number, else it returns -1. */
	int NSides();
	
	/** Creates graphical element to postprocessing */
	void CreateGraphicalElement(TPZGraphMesh &grmesh, int dimension) override;
		
	static void ListOfGroupings(TPZCompMesh *finemesh,TPZVec<int64_t> &accumlist,int nivel,int64_t &numaggl,int dim);
	
	//  void FineSolution(TPZVec<REAL> &x,TPZCompElDisc *disc,TPZVec<REAL> &uh);
	
	void Print(TPZStack<int64_t> &listindex);

	template<class TVar>
	void ProjectSolution(TPZFMatrix<TVar> &projectsol);
	
	
	static TPZAgglomerateMesh *CreateAgglomerateMesh(TPZCompMesh *finemesh,TPZVec<int64_t> &accumlist,int64_t numaggl);
	
	static void ComputeNeighbours(TPZCompMesh *mesh, std::map<TPZCompElDisc *,std::set<TPZCompElDisc *> > &neighbours);
	
	/** @brief Returns the unique identifier for reading/writing objects to streams */
	public:
int ClassId() const override;

	/*@brief Save the element data to a stream */
	void Write(TPZStream &buf, int withclassid) const override;
	
	/** @brief Read the element data from a stream */
	void Read(TPZStream &buf, void *context) override;

};

extern template
void TPZAgglomerateElement::ProjectSolution<STATE>(TPZFMatrix<STATE> &projectsol);
extern template
void TPZAgglomerateElement::ProjectSolution<CSTATE>(TPZFMatrix<CSTATE> &projectsol);

#endif
