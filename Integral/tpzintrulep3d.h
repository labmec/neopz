/**
 * @file
 * @brief Contains the TPZIntRuleP3D class which defines the integration rule for pyramid.
 */

#ifndef TPZINTRULEP3D_H
#define TPZINTRULEP3D_H

#include "pzreal.h"
#include "pzvec.h"

/**
 * @ingroup integral
 * @brief Integration rule for pyramid. \ref integral "Numerical Integration"
 * @author Philippe R. B. Devloo <phil@fec.unicamp.br>
 */
class TPZIntRuleP3D {
	
	/** @brief The list can to access the constructor of the current class. */
	friend class TPZIntRuleList;
	
	/** @brief Number of integration points for this object */
    int	   fNumInt;
	/** @brief Location of the integration point Ksi */
    TPZVec<long double>	fLocationKsi;
	/** @brief Location of the integration point Eta */
    TPZVec<long double>	fLocationEta;
	/** @brief Location of the integration point ZEta */
    TPZVec<long double>	fLocationZeta;
	/** @brief Weight of the integration point */
    TPZVec<long double>	fWeight;
    /// polynomial order of the integration rule
    int fOrder;

	/** 
	 * @brief Constructor of cubature rule for pyramid 
	 * @param order Order of the polinomial will be integrated exactly with this cubature rule
	 */
	TPZIntRuleP3D(int &order);
	/** @brief Default destructor. */
	~TPZIntRuleP3D();
	

public:
	
    TPZIntRuleP3D(const TPZIntRuleP3D &copy) : fNumInt(copy.fNumInt), fLocationKsi(copy.fLocationKsi), fLocationEta(copy.fLocationEta),
        fLocationZeta(copy.fLocationZeta), fWeight(copy.fWeight), fOrder(copy.fOrder)
    {
        
    }
    
    TPZIntRuleP3D &operator=(const TPZIntRuleP3D &copy)
    {
        fNumInt = copy.fNumInt;
        fLocationKsi = copy.fLocationKsi;
        fLocationEta = copy.fLocationEta;
        fLocationZeta = copy.fLocationZeta;
        fWeight = copy.fWeight;
        fOrder = copy.fOrder;
        return *this;
    }
    
    enum {NRULESPYRAMID_ORDER = 37};
	
	/** @brief Returns number of integration points */
    int NInt() const{ return fNumInt; }
	
	/** @brief Returns location of the ith point */
	void Loc(int i, TPZVec<REAL> &pos) const;
	
	/** @brief Returns weight for the ith point */
	REAL W(int i) const;
    
    /// return the order of the polynomial order that can be integrated
    int Order()
    {
        return fOrder;
    }

	/** @brief Prints the number of integration points, all points and weights (as one dimension) */
	void Print(std::ostream & out = std::cout);

protected:

	/**
	 * @brief Computes the points and weights for pyramid cubature rule as first version of PZ 
	 * @param order Order of the polinomial will be integrated exactly with this cubature rule
	 */
	int ComputingCubatureRuleForPyramid(int order);

};

#endif
