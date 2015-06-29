/**
 * @file TPZMatComplexExample2D.h
 * @brief Header file for class TPZMatComplexExample2D.\n
 * It implements the weak statement of the model problem from Oden's book, Chapter 1 within the PZ environment.
 */

#ifndef TPZMATCOMPLEXEXAMPLE2D_H
#define TPZMATCOMPLEXEXAMPLE2D_H

#include "pzmaterial.h"

const REAL M_C  (3*1e8); //velocidade da luz no vacuo
const REAL M_UZERO  (4*M_PI*1e-7);//permeabilidade do meio livre
const REAL M_EZERO  (8.854*1e-12);//permissividade do meio livre
const STATE imaginary(0.,1.);//unidade imaginaria
/**
 * @ingroup material
 * @brief This class implements the weak statement of the model problem from Oden's book, Chapter 1, within the PZ environment
 */
class  TPZMatComplexExample2D : public TPZMaterial
{
    
protected:
  
  //COM CERTEZA
  STATE (& fUr)( TPZVec<REAL>,REAL);
  STATE (& fEr)( TPZVec<REAL>,REAL);
  REAL fLambda;//comprimento de onda
  REAL fTheta;//angulo de incidencia da onda
  REAL fEZero;//magnitude do campo
  REAL fW;//velocidade angular da onda
  REAL fKZero;//numero de onda
   
	
public:
    TPZMatComplexExample2D(int id, REAL l, REAL t, REAL eO,STATE (& ur)( TPZVec<REAL>,REAL),STATE (& er)( TPZVec<REAL>,REAL));
  
    TPZMatComplexExample2D(int id);
  
    /** @brief Default constructor */
    TPZMatComplexExample2D();
    
    /** @brief Creates a material object based on the referred object and inserts it in the vector of material pointers of the mesh. */
	/** Upon return vectorindex contains the index of the material object within the vector */
    TPZMatComplexExample2D(const TPZMatComplexExample2D &mat);
    /** @brief Default destructor */
    virtual ~TPZMatComplexExample2D();
	
    /** @brief Returns the name of the material */
    virtual std::string Name() { return "TPZMatComplexExample2D"; }
    
    /** @brief Returns the integrable dimension of the material */
    virtual int Dimension() const {return 2;}
  
    /** @brief Change the wave incidence angle */
    void SetTheta(REAL t) {fTheta=t;}
  
    /** @brief Returns the number of state variables associated with the material */
    virtual int NStateVariables() { return 1;}
    
public:
    
    /** @name Contribute methods
	 * @{
	 */
	
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @since April 16, 2007
     */
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
		
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     * @since October 07, 2011
     */
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
  	
};

STATE urDefault( TPZVec<REAL>x,REAL l=1e-6 );//default material has ur=1
STATE erDefault( TPZVec<REAL>x,REAL l=1e-6);//default material has er=1
#endif

