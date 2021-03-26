/**
 * \file
 * @brief Contains the TPZPlacaOrthotropic class.
 */

#ifndef PLACAORTHOTROPIC
#define PLACAORTHOTROPIC
#include "pzvec.h"
#include "pzfmatrix.h"

class TPZInterpolatedElement;
class TPZCompEl;
class TPZGeoEl;

/**
 * @ingroup material
 * @brief O objeto desta classe representa uma placa do objeto multicamada
 */
/**
 * Este ultimo representa fisicamente a superposicao de varias placas. 
 * Nao eh um material propriamente.
 */
class TPZPlacaOrthotropic {
	
private:
	
	TPZGeoEl *fGeoEl;
	/** @brief Computational element related with the shell */
	TPZInterpolatedElement *fIntel;
	/**Espessura da placa*/
	REAL fH;// = fZmax - fZmin
	/**alturas mínimas e máximas da placa horizontal*/
	REAL fZMin, fZMax;

	int fTensorVar;
	
public:
	/** @brief Default constructor */
	TPZPlacaOrthotropic();
	/** @brief Constructor */
	TPZPlacaOrthotropic(TPZGeoEl *gel,REAL zmin, REAL zmax);
	/** @brief Destrutor */
	~TPZPlacaOrthotropic(){}

	/**
	 * @brief Returns the tensions tensor of the shell
	 * @param ksi point in the parametric space
	 * @param T Tension values [out]
	 */
	void Tensor(TPZVec<REAL> &ksi, TPZFMatrix<REAL> &T);

	REAL Moment(REAL zref, TPZVec<REAL> &normal, TPZVec<REAL> &direction);
	REAL Force(TPZVec<REAL> &normal, TPZVec<REAL> &direction);
	
	/**
	 * gradiente do momento na direcao indicada por graddir (a direcao e dada em espaco parametrico)
	 * a derivada e devolvida em espaco real
	 */
	REAL GradMoment(REAL zref, TPZVec<REAL> &graddir, TPZVec<REAL> &normal, TPZVec<REAL> &direction);
	
	/**
	 * gradiente da forca na direcao indicada por graddir (a direcao e dada em espaco parametrico)
	 * a derivada e devolvida em espaco real
	 */
	REAL GradForce(TPZVec<REAL> &graddir, TPZVec<REAL> &normal, TPZVec<REAL> &direction);
	
	/**
	 * gradiente do tensor na direcao indicada por graddir (a direcao e dada em espaco parametrico)
	 * a derivada e devolvida em espaco real
	 */
	void GradTensor(TPZVec<REAL> &graddir, TPZVec<REAL> &ksi,  TPZFMatrix<REAL> &gradtensor);
	
	void PrintTensors(std::ostream &out);
	
	void PrintTensors(std::ostream &out,TPZFMatrix<REAL> &tensorin,TPZFMatrix<REAL> &tensorout);
	
	REAL Height(){return fH;}
	
	REAL ZMin() { return fZMin;}
	
	REAL ZMax() { return fZMax;}
	
	void IdentifyCompEl();
	
	void Print();
	
	TPZInterpolatedElement *ComputEl() {return fIntel;}

};

#endif
