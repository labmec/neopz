/**
 * @file
 * @brief Contains the TPZInt1d, TPZIntTriang, TPZIntQuad, TPZIntCube3D, TPZIntTetra3D, TPZIntPyram3D and TPZIntPrism3D classes \n
 * which handles the numerical integration.
 */

#ifndef INTQUADHPP
#define INTQUADHPP

#include <iostream>

#include "pzreal.h"
template<class T>
class TPZVec;
template<class T, int N>
class TPZManVector;

class TPZIntRuleT3D;
class TPZIntRuleP3D;

#include "tpzintrulet.h"
#include "tpzintpoints.h"
#include "tpzprinteg.h"

/** \addtogroup integral
 * @{
 */

//*******************************************************************
// Base Class TInt1D - 	which handles the integration
//								for 1D problems
//*******************************************************************
/** 
 * @brief Handles the numerical integration for one-dimensional problems. \ref integral "Numerical Integration"
 */
class TPZInt1d : public TPZIntPoints {
	int fOrdKsi;
	TPZGaussRule *fIntP;
public:
	enum {Dim = 1};
	TPZInt1d(int OrdK = 0,int type = 0);
	TPZInt1d(const TPZInt1d &copy ) : TPZIntPoints(copy), fOrdKsi(copy.fOrdKsi), fIntP(copy.fIntP) {
	}
    TPZInt1d &operator=(const TPZInt1d &copy)
    {
        fOrdKsi = copy.fOrdKsi;
        fIntP = copy.fIntP;
        return *this;
    }
	virtual ~TPZInt1d() {
	}
	virtual int NPoints() const override;
	virtual void Point(int ip, TPZVec<REAL> &pos, REAL &w) const override;

	virtual void SetType(int type,int order) override {
		fIntP->SetType(type,order);
	}
	virtual void SetOrder(TPZVec<int> &ord,int type = 0) override;
	virtual void GetOrder(TPZVec<int> &ord) const override;
	virtual int GetRealMaxOrder() const;
	
	/** A ordem máxima é muito maior que 36. Seu valor está em GetRealMaxOrder.
	 * Entretanto, é preciso dar algum valor aqui para sobrescrever o da classe
	 * mãe que retorna o max order do tetraedro, da pirâmide ou do triângulo,
	 * que tradicionalmente tem regras menores que linha, quadrilátero e hexaedro
	 */ 
	virtual int GetMaxOrder() const override{
		return 36;
	}
	
	virtual int Dimension() const override
	{
		return Dim;
	}

	virtual TPZIntPoints *PrismExtend(int order) override
	{
		return new TPZPrInteg<TPZInt1d>(order);
	}
	virtual TPZIntPoints *Clone() const override
	{
		return new TPZInt1d(*this);
	}
    
	void Print(std::ostream &out = std::cout) const override{
		if(fIntP) fIntP->Print(out);
	}

	/** @brief Returns the name of the cubature rule */
	void Name(std::string &name) const override{
		name = "TPZInt1D";
	}
};

//*******************************************************************
// Base Class TPZIntTriang - which handles the integration
//									for 2D problems, triangle elements
//*******************************************************************
/**
 * @brief Handles the numerical integration for two-dimensional problems using triangular elements. \ref integral "Numerical Integration"
 */
class TPZIntTriang : public TPZIntPoints{
    
protected:
	int fOrdKsi;
	TPZIntRuleT *fIntKsi;

public:
	enum {Dim = 2};
	TPZIntTriang(int OrdK = 2);
	TPZIntTriang(const TPZIntTriang &copy) : TPZIntPoints(copy), fOrdKsi(copy.fOrdKsi), fIntKsi(copy.fIntKsi) {
	}
    TPZIntTriang &operator=(const TPZIntTriang &copy)
    {
        TPZIntPoints::operator=(copy);
        fOrdKsi = copy.fOrdKsi;
        fIntKsi = copy.fIntKsi;
        return *this;
    }
	virtual ~TPZIntTriang() {
	}

	virtual int  NPoints() const override;
	virtual void Point(int ip, TPZVec<REAL> &pos, REAL &w) const override;

	virtual void SetOrder(TPZVec<int> &ord,int type = 0) override;
	virtual void GetOrder(TPZVec<int> &ord) const override;
	virtual int GetMaxOrder() const override;
	virtual int Dimension() const override
	{
		return Dim;
	}

	virtual TPZIntPoints *PrismExtend(int order) override
	{
		return new TPZPrInteg<TPZIntTriang>(order);
	}
	virtual TPZIntPoints *Clone() const override
	{
		return new TPZIntTriang(*this);
	}
	
	/** @brief Returns the name of the cubature rule */
	void Name(std::string &name) const override {
		name = "TPZIntTriang";
	}
};

//*******************************************************************
// Base Class TPZIntQuad - 	which handles the integration
//									for 2D problems, quadrilaterals elements
//*******************************************************************
/** 
 * @brief Handles the numerical integration for two-dimensional problems using quadrilateral elements. \ref integral "Numerical Integration"
 */
class TPZIntQuad : public TPZIntPoints{
	int fOrdKsi;
	int fOrdEta;
	TPZGaussRule *fIntKsi;
	TPZGaussRule *fIntEta;
public:
	enum {Dim = 2};
	/**
	 * @brief Constructor with two one dimensional rules.
	 * @param OrdK Order for \f$ \xi \f$ axe
	 * @param OrdE Order for \f$ \eta \f$ axe at master element.
	 */
	TPZIntQuad(int OrdK,int OrdE);
    
	TPZIntQuad(int OrdK = 1);
    
	/** @brief Copy constructor */
	TPZIntQuad(const TPZIntQuad &copy) : TPZIntPoints(copy), fOrdKsi(copy.fOrdKsi), fOrdEta(copy.fOrdEta), fIntKsi(copy.fIntKsi),fIntEta(copy.fIntEta) {
	}
    TPZIntQuad &operator=(const TPZIntQuad &copy)
    {
        TPZIntPoints::operator=(copy);
        fOrdKsi = copy.fOrdKsi;
        fOrdEta = copy.fOrdEta;
        fIntKsi = copy.fIntKsi;
        fIntEta = copy.fIntEta;
        return *this;
    }
	/** @brief Destructor */
	virtual ~TPZIntQuad() {
	}
	
	virtual int NPoints() const override;
	virtual void Point(int ip, TPZVec<REAL> &pos, REAL &w) const override;

    /// Set the order and the type of integration rule :
    // type 0 : Gauss Legendre
    // type 1 : Gauss Lobatto
    // type 2 : Jacobi
    // type 3 : Chebyshev
	virtual void SetType(int type,int order) override{
		fIntKsi->SetType(type,order);
		fIntEta->SetType(type,order);
	}
    
    /// Set the order and the type of integration rule :
    // type 0 : Gauss Legendre
    // type 1 : Gauss Lobatto
    // type 2 : Jacobi
    // type 3 : Chebyshev
	virtual void SetOrder(TPZVec<int> &ord,int type = 0) override;
    
    /// Return the order of the integration rule
	virtual void GetOrder(TPZVec<int> &ord) const override;
	virtual int GetRealMaxOrder() const;  
	
	/** A ordem máxima é muito maior que 36. Seu valor está em GetRealMaxOrder.
	 * Entretanto, é preciso dar algum valor aqui para sobrescrever o da classe
	 * mãe que retorna o max order do tetraedro, da pirâmide ou do triângulo,
	 * que tradicionalmente tem regras menores que linha, quadrilátero e hexaedro
	 */ 
	virtual int GetMaxOrder() const override{
		return 36;
	}
	
	virtual int Dimension() const override
	{
		return Dim;
	}

	virtual TPZIntPoints *PrismExtend(int order) override
	{
		return new TPZPrInteg<TPZIntQuad>(order);
	}
	virtual TPZIntPoints* Clone() const override
	{
		return new TPZIntQuad(*this);
	}

	/** @brief Returns the name of the cubature rule */
	void Name(std::string &name) const override{
		name = "TPZIntQuad";
	}
};


/** 
 * @brief Handles the numerical integration for three-dimensional problems using cube elements. \ref integral "Numerical Integration"
 */

class TPZIntCube3D : public TPZIntPoints{
	int fOrdKsi;
	int fOrdEta;
	int fOrdZeta;
	TPZGaussRule *fIntKsi;
	TPZGaussRule *fIntEta;
	TPZGaussRule *fIntZeta;
public:
	enum {Dim = 3};
  /**
	 * @brief Constructor with three one dimensional rules of same order in each axis.
	 * @param ord Order for every axis
	 */
	TPZIntCube3D(int OrdK = 2) : TPZIntCube3D(OrdK,OrdK,OrdK){}
	/**
	 * @brief Constructor with three one dimensional rules.
	 * @param OrdK Order for \f$ \xi \f$ axe
	 * @param OrdE Order for \f$ \eta \f$ axe
	 * @param OrdZ Order for \f$ \zeta \f$ axe at master element.
	 */
	TPZIntCube3D(int OrdK,int OrdE,int OrdZ);
	/** @brief Copy constructor */
	TPZIntCube3D(const TPZIntCube3D &copy) : TPZIntPoints(copy), fOrdKsi(copy.fOrdKsi), fOrdEta(copy.fOrdEta), fOrdZeta(copy.fOrdZeta),
		fIntKsi(copy.fIntKsi), fIntEta(copy.fIntEta), fIntZeta(copy.fIntZeta) {
	}
    
    TPZIntCube3D &operator=(const TPZIntCube3D &copy)
    {
        TPZIntPoints::operator=(copy);
        fOrdKsi = copy.fOrdKsi;
        fOrdEta = copy.fOrdEta;
        fOrdZeta = copy.fOrdZeta;
        fIntKsi = copy.fIntKsi;
        fIntEta = copy.fIntEta;
        fIntZeta = copy.fIntZeta;
        return *this;
    }

	/** @brief Destructor */
	virtual ~TPZIntCube3D() {
	}
	
	virtual int NPoints() const override;
	virtual void Point(int ip, TPZVec<REAL> &pos, REAL &w) const override;

	virtual void SetType(int type,int order) override {
		fIntKsi->SetType(type,order);
		fIntEta->SetType(type,order);
		fIntZeta->SetType(type,order);
	}
	virtual void SetOrder(TPZVec<int> &ord,int type = 0) override;
	virtual void GetOrder(TPZVec<int> &ord) const override;
	virtual int GetRealMaxOrder() const;
	virtual int GetMaxOrder() const override{
		return 20;    // With this order the number of integration points are 8000
	}
	
	virtual int Dimension() const override
	{
		return Dim;
	}

	virtual TPZIntPoints *PrismExtend(int order) override
	{
		return new TPZPrInteg<TPZIntCube3D>(order);
	}
	virtual TPZIntPoints *Clone() const override
	{
		return new TPZIntCube3D(*this);
	}

	/** @brief Returns the name of the cubature rule */
	void Name(std::string &name) const override{
		name = "TPZIntCube3D";
	}
};


/** 
 * @brief Handles the numerical integration for three-dimensional problems using tetraedra elements. \ref integral "Numerical Integration"
 */
class TPZIntTetra3D : public TPZIntPoints {
	int fOrdKsi;
	TPZIntRuleT3D *fIntKsi;
public:
	enum {Dim = 3};
	TPZIntTetra3D(int OrdK = 2);
	TPZIntTetra3D(const TPZIntTetra3D &copy) : TPZIntPoints(copy), fOrdKsi(copy.fOrdKsi), fIntKsi(copy.fIntKsi) {
	}

    TPZIntTetra3D &operator=(const TPZIntTetra3D &copy)
    {
        TPZIntPoints::operator=(copy);
        fOrdKsi = copy.fOrdKsi;
        fIntKsi = copy.fIntKsi;
        return *this;
    }
    
	virtual int NPoints() const override;
	virtual void Point(int ip, TPZVec<REAL> &pos, REAL &w) const override;

	virtual void SetOrder(TPZVec<int> &ord,int type = 0) override;
	virtual void GetOrder(TPZVec<int> &ord) const override;
	virtual int GetMaxOrder() const override;  
	virtual int Dimension() const override
	{
		return Dim;
	}

	virtual TPZIntPoints *PrismExtend(int order) override
	{
		return new TPZPrInteg<TPZIntTetra3D>(order);
	}
	virtual TPZIntPoints* Clone() const override
	{
		return new TPZIntTetra3D(*this);
	}

	/** @brief Returns the name of the cubature rule */
	void Name(std::string &name) const override {
		name = "TPZIntTetra3D";
	}
};


/** 
 * @brief Handles the numerical integration for three-dimensional problems using pyramid elements. \ref integral "Numerical Integration"
 */
class TPZIntPyram3D : public TPZIntPoints {
	int fOrdKsi;
	TPZIntRuleP3D *fIntKsi;
public:
	enum {Dim =3};
	TPZIntPyram3D(int OrdK = 2);
	TPZIntPyram3D(const TPZIntPyram3D &copy) : TPZIntPoints(copy), fOrdKsi(copy.fOrdKsi), 
			fIntKsi(copy.fIntKsi) {
	}
	
    TPZIntPyram3D &operator=(const TPZIntPyram3D &copy)
    {
        TPZIntPoints::operator=(copy);
        fOrdKsi = copy.fOrdKsi;
        fIntKsi = copy.fIntKsi;
        return *this;
    }
    
	virtual int NPoints() const override;
	virtual void Point(int ip, TPZVec<REAL> &pos, REAL &w) const override;
	
	virtual void SetOrder(TPZVec<int> &ord,int type = 0) override;
	virtual void GetOrder(TPZVec<int> &ord) const override;
	virtual int GetMaxOrder() const override;  
	virtual int Dimension() const override
	{
		return Dim;
	}

	virtual TPZIntPoints *PrismExtend(int order) override
	{
		return new TPZPrInteg<TPZIntPyram3D>(order);
	}
	TPZIntPoints *Clone() const override
	{
		return new TPZIntPyram3D(*this);
	}

	/** @brief Returns the name of the cubature rule */
	void Name(std::string &name) const override {
		name = "TPZIntPyram3D";
	}
};

/**
 * @brief Handles the numerical integration for three-dimensional problems using prism elements. \ref integral "Numerical Integration"
 * This cubature rule uses a cubature rule for one dimension (zeta) and a cubature rule for triangle (base).
 */
class TPZIntPrism3D  : public TPZIntPoints {
	int fOrdKsi,fOrdKti;
	TPZInt1d fIntRule1D;
	TPZIntTriang fIntTriang;
public:
	
	enum {Dim = 3};
	/**
	 * @brief Constructor with orders for a one dimensional rule and a cubature rule for triangle.
	 * @param OrdK Order for one dimensional cubature rule
	 * @param OrdL Order for cubature rule for triangle
	 */
	TPZIntPrism3D(int OrdK = 2,int OrdL = 2);
	/** @brief Copy constructor */
	TPZIntPrism3D(const TPZIntPrism3D &copy) : TPZIntPoints(copy), fOrdKsi(copy.fOrdKsi), fOrdKti(copy.fOrdKti), 
			fIntRule1D(copy.fIntRule1D), fIntTriang(copy.fIntTriang) {
	}
    TPZIntPrism3D &operator=(const TPZIntPrism3D &copy)
    {
        TPZIntPoints::operator=(copy);
        fOrdKsi = copy.fOrdKsi;
        fOrdKti = copy.fOrdKti;
        fIntRule1D = copy.fIntRule1D;
        fIntTriang = copy.fIntTriang;
        return *this;
    }
	/** @brief Destructor */
	virtual ~TPZIntPrism3D();
	
	int NPoints() const override;
	void Point(int ip, TPZVec<REAL> &pos, REAL &w) const override;

	virtual void SetType(int type,int order)  override{
		fIntRule1D.SetType(type,order);
	}
	void SetOrder(TPZVec<int> &ord,int type = 0) override;
	void GetOrder(TPZVec<int> &ord) const override;
	virtual int GetMaxOrder() const override;  
	virtual int Dimension() const override
	{
		return Dim;
	}

	virtual TPZIntPoints *PrismExtend(int order) override
	{
		return new TPZPrInteg<TPZIntPrism3D>(order);
	}
	TPZIntPoints *Clone() const override
	{
		return new TPZIntPrism3D(*this);
	}

	/** @brief Returns the name of the cubature rule */
	void Name(std::string &name) const override {
		name = "TPZIntPrism3D";
	}
};

/**
 * @brief Integration rule for one point. \ref integral "Numerical Integration"
 */
class TPZInt1Point : public TPZIntPoints {
	
public:
	
    enum {Dim = 0};
    TPZInt1Point(int order = 0);

	TPZInt1Point(const TPZInt1Point &copy ) : TPZIntPoints(copy) {
	}
    TPZInt1Point & operator=(const TPZInt1Point &copy )
    {
        TPZIntPoints::operator=(copy);
        return *this;
    }
    virtual ~TPZInt1Point();
	
    int NPoints() const override;
    void Point(int ip, TPZVec<REAL> &pos, REAL &w) const override;

    void SetOrder(TPZVec<int> &ord,int type = 0) override;
    void GetOrder(TPZVec<int> &ord) const override;
    int GetMaxOrder() const override;  
    int Dimension() const override {
		return Dim;
    }

    TPZIntPoints *PrismExtend(int order) override;
	TPZIntPoints *Clone() const override {
		return new TPZInt1Point(*this);
	}

	/** @brief Returns the name of the cubature rule */
	virtual void Name(std::string &name) const override {
		name = "TPZInt1Point";
	}
};

inline TPZInt1Point::TPZInt1Point(int order) : TPZIntPoints()
{
}

inline TPZInt1Point::~TPZInt1Point() {
}

inline int TPZInt1Point::NPoints() const{
	return 1;
}

inline void TPZInt1Point::Point(int ip, TPZVec<REAL> &pos, REAL &w) const {
#ifndef PZNODEBUG
	if(ip!=0) {
		std::cout << "TPZInt1Point:: Bad number point " << ip << std::endl;
		return;
	}
#endif
	w = 1.;
}

inline void  TPZInt1Point::SetOrder(TPZVec<int> &ord,int type) {
}

inline void TPZInt1Point::GetOrder(TPZVec<int> &/* ord */) const {
}

inline int TPZInt1Point::GetMaxOrder() const {
	return 1;
}

inline TPZIntPoints *TPZInt1Point::PrismExtend(int order) {
	return new TPZPrInteg<TPZInt1Point> (order);
}

/** @} */

#endif
