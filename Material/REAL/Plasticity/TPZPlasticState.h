// $Id: TPZPlasticState.h,v 1.8 2009-07-19 05:47:18 erick Exp $

#ifndef TPZPLASTICSTATE_H
#define TPZPLASTICSTATE_H

#include "fadType.h" 
#include "TPZTensor.h"
#include <iostream>

/**
 * This class holds the complete set of state variables to define a valid elastoplastic strain state
 */
template <class T>
class TPZPlasticState : public TPZSavable {

    // class member variables
	
public:
    
    //@TODO:: Rename the variables: fEpsT -> feps_t, fEpsP -> feps_p, and  fMType->fm_type
	// Tensors representing the total and plastic strain states
	TPZTensor<T> fEpsT, fEpsP;
	
	// Plastic damage variable
	T fAlpha;
    
    // Identifier for the regime of the material behaviour
    int fMType;

public:
	
	/**
	 * Default constructor - all values set to zero
	 */
    TPZPlasticState(): fEpsT(), fEpsP(), fAlpha(T(0.)), fMType(0) { }
		
	/**
	 * Constructor enabling predefinition of Alpha
	 */
	TPZPlasticState(const T & alpha):fEpsT(T(0.)), fEpsP(T(0.)), fAlpha(alpha), fMType(0) { }
	
	/**
	 * Copy constructor
	 */
	TPZPlasticState(const TPZPlasticState<T> & source):
		fEpsT(source.fEpsT), fEpsP(source.fEpsP), fAlpha(source.fAlpha), fMType(source.fMType){ }
	
	/**
	 * Default destructor
	 */
	~TPZPlasticState(){ }
	
	/**
	 * Operator =
	 */
	const TPZPlasticState<T> & operator=(const TPZPlasticState<T> & source);
	
	/**
	 * Operator -=
	 */
	const TPZPlasticState<T> & operator-=(const TPZPlasticState<T> & source);
	
	/**
	 * Operator +=
	 */
	const TPZPlasticState<T> & operator+=(const TPZPlasticState<T> & source);
	
	/**
	 * Operator *=
	 */
	const TPZPlasticState<T> & operator*=(const TPZPlasticState<T> & source);
	
	/**
	 * Similar to Operator=, but allows copies among different template specializations.
	 * When using FAD class Types as template classes, NO DERIVATIVES are copied using this functions,
	 * enabling copies from REAL to FAD types and also vice-versa.
	 */
	template <class T1>
	void CopyTo(TPZPlasticState<T1> & target) const;

	/**
	 * Operator<<
	 */
    friend std::ostream& operator<<( std::ostream& Out, const TPZPlasticState<T> & s )
	{
		s.Print(Out);
		return Out;
	}
	
	/**
	 * More complete then Operator<< because it allows derivatives supression.
	 */
	void Print(std::ostream& Out, int fadDerivatives = 1)const;
	
	// class Access Members (needed for const PlasticState access)
	
	const TPZTensor<T> & EpsT() const 
		{ return fEpsT; }
    
	const TPZTensor<T> & EpsP() const 
		{ return fEpsP; }
    
	const T & Alpha() const 
		{ return fAlpha; }
    
    const int & MType() const
        { return fMType; }

    int ClassId() const;
        
    void Read(TPZStream& buf, void* context){
        fEpsT.Read(buf,context);
        fEpsP.Read(buf,context);
        
        buf.Read(&fAlpha);
        buf.Read(&fMType);
    }
    
    void Write(TPZStream& buf, int withclassid) const{
        fEpsT.Write(buf,withclassid);
        fEpsP.Write(buf,withclassid);
	
        buf.Write(&fAlpha);
        buf.Write(&fMType);
    }
    
    void CleanUp() {
        fEpsT.Zero();
        fEpsP.Zero();
        fAlpha = T(0.);
        fMType = 0;
    }
};

template <class T>
int TPZPlasticState<T>::ClassId() const{
    return Hash("TPZPlasticState") ^ ClassIdOrHash<T>() << 1;
}

template <class T>
inline const TPZPlasticState<T> & TPZPlasticState<T>::operator=(const TPZPlasticState<T> & source)
{
	fEpsT = source.EpsT();
	fEpsP = source.EpsP();
	fAlpha = source.Alpha();
    fMType = source.MType();
		
	return *this;
}

template <class T>
inline const TPZPlasticState<T> & TPZPlasticState<T>::operator+=(const TPZPlasticState<T> & source)
{
	fEpsT += source.EpsT();
	fEpsP += source.EpsP();
	fAlpha+= source.Alpha();
	return *this;
}

template <class T>
inline const TPZPlasticState<T> & TPZPlasticState<T>::operator-=(const TPZPlasticState<T> & source)
{
	fEpsT -= source.EpsT();
	fEpsP -= source.EpsP();
	fAlpha-= source.Alpha();
	return *this;
}

template <class T>
inline const TPZPlasticState<T> & TPZPlasticState<T>::operator*=(const TPZPlasticState<T> & source)
{
	fEpsT *= source.EpsT();
	fEpsP *= source.EpsP();
	fAlpha*= source.Alpha();
	return *this;
}

template <class T>
inline void TPZPlasticState<T>::Print(std::ostream& Out, int fadDerivatives)const {
    if (fadDerivatives) {
        Out << "\tfEpsT = ";
        for (int i = 0; i < 6; i++)Out << fEpsT[i] << " ";
        Out << std::endl;
        Out << "\tfEpsP = ";
        for (int i = 0; i < 6; i++)Out << fEpsP[i] << " ";
        Out << std::endl;
        Out << "\tfAlpha = " << fAlpha << std::endl;
    } else {
        Out << "\tfEpsT = ";
        for (int i = 0; i < 6; i++)Out << TPZExtractVal::val(fEpsT[i]) << " ";
        Out << std::endl;
        Out << "\tfEpsP = ";
        for (int i = 0; i < 6; i++)Out << TPZExtractVal::val(fEpsP[i]) << " ";
        Out << std::endl;
        Out << "\tfAlpha = " << TPZExtractVal::val(fAlpha) << std::endl;
    }
    Out << "\tfMtype = " << fMType << std::endl;
}

template <class T>
template <class T1>
void TPZPlasticState<T>::CopyTo(TPZPlasticState<T1> & target) const
{
	EpsT().CopyTo(target.fEpsT);
	EpsP().CopyTo(target.fEpsP);
	target.fAlpha = TPZExtractVal::val( Alpha() );
    target.fMType = MType();
}


#endif
