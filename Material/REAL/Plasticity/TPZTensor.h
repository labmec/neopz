// $Id: TPZTensor.h,v 1.28 2010-12-04 20:41:28 diogo Exp $

#ifndef TPZTENSOR_H
#define TPZTENSOR_H

#include "pzmanvector.h"
#include "pzfmatrix.h"
#include <iostream>
#include "fadType.h"
#include <math.h>
#include "pzlog.h"
#ifndef WIN32
#include <fenv.h>//NAN DETECTOR
#pragma STDC FENV_ACCESS ON
#endif

#define _XX_ 0
#define _XY_ 1
#define _XZ_ 2
#define _YY_ 3
#define _YZ_ 4
#define _ZZ_ 5

#ifdef LOG4CXX
static LoggerPtr loggerr(Logger::getLogger("logtensor"));
#endif
/**
 Classe que implementa o comportamento de um tensor simetrico
 */
template <class T>
class TPZTensor {
public:
	
	struct TPZDecomposed
	{
		TPZManVector<TPZTensor<T>, 3> fEigenvectors;
		TPZManVector<T, 3> fEigenvalues;
		
		TPZDecomposed() : fEigenvectors(3), fEigenvalues(3,0.)
		{
			
		}
		TPZDecomposed(const TPZDecomposed &copy) : fEigenvectors(copy.fEigenvectors),fEigenvalues(copy.fEigenvalues)
		{
		}
		TPZDecomposed &operator=(const TPZDecomposed &copy)
		{
			fEigenvectors = copy.fEigenvectors;
			fEigenvalues = copy.fEigenvalues;
			return *this;
		}
		
		TPZDecomposed(const TPZTensor<T> &source) : fEigenvectors(3), fEigenvalues(3,0.)
		{
			source.EigenSystem(*this);
		}
		
		void Print(std::ostream &out) const
		{
			out << "TPZTensor::Decomposed Eigenvalues " << fEigenvalues << std::endl;
			out << "Eigenvectors:\n";
			for (int i=0; i<fEigenvectors.size(); i++) {
				fEigenvectors[i].Print(out);
			}
		}
		/**
		 * Methods for checking convergence
		 */
		/// Number of test cases implemented by this class
		static int NumCasesCheckConv()
		{
			return 3;
		}
		
		static STATE gEigval[3];
		/// Compute the tangent matrix for a particular case
		static void TangentCheckConv(TPZFMatrix<STATE> &state, TPZFMatrix<STATE> &tangent, int icase)
		{
			TPZTensor<STATE> obj;
			for (int i=0; i<6; i++)
			{
				obj[i] = state(i);
			}
			TPZTensor<STATE>::TPZDecomposed objdec(obj);
			switch (icase) {
				case 0:
				case 1:
				case 2:
				{
					tangent.Resize(1, 9);
					for (int i=0; i<3; i++) {
						gEigval[i] = objdec.fEigenvalues[i];
					}
					for (int j=0; j<6; j++) {
						tangent(0,j) = objdec.fEigenvectors[icase][j];
						tangent(0,6) = objdec.fEigenvectors[icase].XY();
						tangent(0,7) = objdec.fEigenvectors[icase].XZ();
						tangent(0,8) = objdec.fEigenvectors[icase].YZ();
					}
					STATE sum=0.,sum2=0.;
					for (int i=0; i<9; i++) {
						sum += tangent(0,i)*tangent(0,i);
					}
					sum2 = tangent(0,_XX_)+tangent(0,_YY_)+tangent(0,_ZZ_);
					tangent.Print("tangent");
					std::cout << "sum2 = " << sum2 << std::endl;
					std::cout << "sum = " << sum << std::endl;
				}
					break;
				default:
					break;
			}
		}
		/// Compute the residual for a particular state
		static void ResidualCheckConv(TPZFMatrix<STATE> &state, TPZFMatrix<STATE> &residual, int icase)
		{
			TPZTensor<STATE> obj;
			for (int i=0; i<6; i++)
			{
				obj[i] = state(i);
			}
			obj[_XY_] = obj[_XY_]*0.5+0.5*state(6);
			obj[_XZ_] = (obj[_XZ_]+state(7))*0.5;
			obj[_YZ_] = (obj[_YZ_]+state(8))*0.5;
			TPZTensor<STATE>::TPZDecomposed objdec(obj);
			residual.Resize(1, 1);
			switch (icase) {
				case 0:
				case 1:
				case 2:
				{
					residual(0,0) = objdec.fEigenvalues[icase];
					int numequal = 1;
					for (int i=0; i<3; i++) {
						if (i==icase) {
							continue;
						}
						if (gEigval[i] == gEigval[icase]) {
							numequal++;
							residual(0) += objdec.fEigenvalues[i];
						}
					}
					residual(0) /= numequal;
				}
					break;
					
				default:
					break;
			}
		}
		
		void ComputeJ2(REAL & J2)
		{
			REAL sig1,sig2,sig3;
			sig1 = this->fEigenvalues[0];
			sig2 = this->fEigenvalues[1];
			sig3 = this->fEigenvalues[2];
			J2=(pow(sig1 + (-sig1 - sig2 - sig3)/3.,2) + pow(sig2 + (-sig1 - sig2 - sig3)/3.,2) + pow((-sig1 - sig2 - sig3)/3. + sig3,2))/2.;
		}
		
		void ComputeJ3(REAL & J3)
		{
			REAL sig1,sig2,sig3,s1,s2,s3,temp;
			sig1 = this->fEigenvalues[0];
			sig2 = this->fEigenvalues[1];
			sig3 = this->fEigenvalues[2];
			s1=sig1-(1./3.)*(sig1+sig2+sig3);
			s2=sig2-(1./3.)*(sig1+sig2+sig3);
			s3=sig3-(1./3.)*(sig1+sig2+sig3);
			temp=(1./3.)*(s1*s1*s1+s2*s2*s2+s3*s3*s3);
			J3=temp;
		}
		
		void ComputeI1(REAL &I1)
		{
			REAL sig1,sig2,sig3;
			sig1 = this->fEigenvalues[0];
			sig2 = this->fEigenvalues[1];
			sig3 = this->fEigenvalues[2];
			I1=sig1+sig2+sig3;
		}
		
		void ApplyStrainComputeElasticStress(TPZTensor<REAL> &Stress,REAL &K , REAL & G)
		{
			REAL sig1,sig2,sig3,s1,s2,s3;
			sig1 = this->fEigenvalues[0];
			sig2 = this->fEigenvalues[1];
			sig3 = this->fEigenvalues[2];
			s1=sig1-(1./3.)*(sig1+sig2+sig3);
			s2=sig2-(1./3.)*(sig1+sig2+sig3);
			s3=sig3-(1./3.)*(sig1+sig2+sig3);
			
			TPZTensor<REAL> SigVol;
			SigVol.XX()=K*(sig1+sig2+sig3);
			SigVol.YY()=K*(sig1+sig2+sig3);
			SigVol.ZZ()=K*(sig1+sig2+sig3);
			
			Stress.XX()=s1;
			Stress.YY()=s2;
			Stress.ZZ()=s3;
			Stress*=(2*G);
			Stress+=SigVol;
			
			
		}
		
		
	};
	/**
	 Construtor vazio inicializando com zero
	 */
#ifdef _AUTODIFF
	TPZTensor() : fData(6, T(0.)) { }    // When T is the Fad type, the elements of the fData are initialized with T zero. See tfad.h
#else
	TPZTensor() : fData(6, T(0.)) { }
#endif
	
	/**
	 Construtor vazio inicializando com Init
	 */
	TPZTensor(const T & Init) : fData(6, Init){ }
	
	/**
	 Copy Constructor
	 */
	TPZTensor(const TPZTensor<T> & source) : fData(source.fData){ }
	
	/**
	 Construct a tensor based on its eigensystem decomposition
	 */
	TPZTensor(const TPZDecomposed &eigensystem) : fData(6, T(0.))
	{
		(*this).Scale(0.);
		for (int i=0; i<eigensystem.fEigenvectors.size(); i++) {
			Add(eigensystem.fEigenvectors[i], eigensystem.fEigenvalues[i]);
		}
	}
	
	/// Method to write to a pzstream
	void Write(TPZStream &out) const
	{
		//const TPZVec<T> *engane = &fData;
		TPZSaveable::WriteObjects(out,fData);
	}
	
	///Method to read the object from a pzstream
	void Read(TPZStream &input)
	{
		TPZSaveable::ReadObjects(input,fData);
	}
	
	operator TPZFMatrix<T>() const
	{
		TPZFMatrix<T> result(3,3);
		result(0,0) = XX();
		result(0,1) = XY();
		result(0,2) = XZ();
		result(1,0) = XY();
		result(1,1) = YY();
		result(1,2) = YZ();
		result(2,0) = XZ();
		result(2,1) = YZ();
		result(2,2) = ZZ();
		return result;
	}
    
//    operator TPZFNMatrix<6>() const
//	{
//		TPZFNMatrix<6> result(6,1,0.);
//		result(0,0) = XX();
//		result(1,0) = XY();
//		result(2,0) = XZ();
//		result(3,0) = YY();
//		result(4,0) = YZ();
//		result(5,0) = ZZ();
//		return result;
//	}
    
	
	TPZTensor(const TPZFMatrix<T> &input) : fData(6,0.)
	{
#ifdef DEBUG
		if (input.Rows() != 3 || input.Cols() != 3) {
			DebugStop();
		}
		for (int i=0; i<3; i++) {
			for (int j=0; j<3; j++) {
				if (fabs(input.GetVal(i,j) - input.GetVal(j,i)) > 1.e-15) {
					std::cout << "diff = " << input.GetVal(i,j)-input.GetVal(j,i) << std::endl;
					DebugStop();
				}
			}
		}
#endif
		fData[_XX_] = input.GetVal(0,0);
		fData[_XY_] = input.GetVal(0,1);
		fData[_XZ_] = input.GetVal(0,2);
		fData[_YY_] = input.GetVal(1,1);
		fData[_YZ_] = input.GetVal(1,2);
		fData[_ZZ_] = input.GetVal(2,2);
	}
	/**
	 * Method to print the tensor
	 */
	void Print(std::ostream &out) const;
	
	/**
	 Operator=
	 */
	const TPZTensor<T> & operator=( const TPZTensor<T> &source );
	
	/**
	 Operator+=
	 */
	const TPZTensor<T> & operator+=( const TPZTensor<T> &source );
	
	/**
	 Operator-=
	 */
	const TPZTensor<T> & operator-=( const TPZTensor<T> &source );
	
	/**
	 Operator*=
	 */
	const TPZTensor<T> & operator*=( const T &multipl );
	
	/** Identity Matrix
	 * TBASE is needed when 3rd derivatives are of interest. In such cases the FAD
	 * promotion fails.
	 */
	template <class TBASE>
	void Identity_T(TBASE &);
	
	//Specialization for TBASE = REAL
	void Identity();
	
	
	
	
	/**
	 Adiciona um tensor no tensor atual
	 @param [in] tensor tensor sendo adicionado
	 @param [in] constant fator multiplicativo
	 */
	template < class T1, class T2 >
	void Add(const TPZTensor< T1 > & tensor, const T2 & constant );
	
	/**
	 multiplica um scalar com o tensor atual
	 @param [in] multipl valor sendo multiplicada
	 @param [in] constant fator multiplicativo
	 */
	template < class T1 ,class T2>
	void Multiply(const T1 & multipl, const T2 & constant );
	
	/**
	 realiza a multiplicacao tradicional entre matrizes
	 @param [in] tensor a ser multiplicado
	 @param [out] tensor resposta
	 */
	void Multiply2(TPZTensor<T> tensor, TPZTensor<T> & resp)const;
	
	/**
	 multiplica um scalar com o tensor atual
	 @param [in] constant fator multiplicativo
	 */
	template < class T2>
	void Scale( const T2 & constant );

	/**
	 Zera o tensor
	 */
	void Zero();
  
	/**
	 Metodo que calcula os autovetores to tensor
	 @param [out] eigVec autovetores
	 */
	void EigenVector(TPZVec<TPZVec<T> > & eigVec) const;
	
	/**
	 Metodo que calcula os autovalores to tensor
	 @param [out] eigVal autovalores
	 */
	//  void EigenValue(TPZVec<T> & eigVal) const;
	
	/**
	 Metodo que calcula os invariantes to tensor
	 (atualmente somente calcula J2 e J3
	 @param [out] jInv invariantes
	 */
	void JInvariants(TPZVec<T> & jInv) const;
	
	/**
	 * Metodo que calcula o negativo do segundo invariante da parte deviatorica do tensor
	 * TBASE is needed when 3rd derivatives are of interest. In such cases the FAD
	 * promotion fails.
	 */
	template <class TBASE>
	T J2_T(TBASE &) const;
	
	// Specialization for TBASE = REAL
	T J2() const;
	
	/**
	 Metodo que calcula a derivada de J2
	 */
	void dJ2(TPZTensor<T> & Tangent) const;
	
	/**
	 * Metodo que calcula o determinante do tensor
	 */
	T Det() const;
	
	/**
	 * Metodo que calcula o determinante do tensor
	 */
	void dDet(TPZTensor<T> &grad) const;
	
	/**
	 Metodo que calcula o primeiro invariante
	 */
	T I1() const;
	
	/**
	 Metodo que calcula o segundo invariante
	 */
	T I2() const;
	
	/**
	 Metodo que calcula o terceiro invariante
	 */
	T I3() const;
	
	T Norm() const;
	
	/**
	 * Metodo que calcula o tensor deviatorio
	 * @param [out] s retorna o tensor deviatorio calculado
	 */
	void S(TPZTensor<T> &s) const;
	
	/**
	 * Return the 3rd invariant of deviatoric tensor
	 */
	T J3() const;
	
	/**
	 * The derivative of the  3rd invariant of deviatoric tensor
	 */
	void dJ3(TPZTensor<T> &deriv) const;
	
	/**
	 * @brief adjust the tensor to the given values of I1 and sqj2
	 */
	void Adjust(TPZVec<T> &sigIJ, TPZTensor<T> &result) const
	{
		TPZTensor<T> S;
		this->S(S);
		TPZTensor<T> Diag;
		Diag.Identity();
		Diag *= sigIJ[0]/T(3.);
		T J2 = this->J2();
		T sqj2 = sqrt(J2);
		REAL sqj2val = shapeFAD::val(sqj2);
		T ratio;
		if (fabs(sqj2val) > 1.e-6)
		{
			ratio = sigIJ[1]/sqj2;
		}
		else
		{
			ratio = T(1.);
		}
		S *= ratio;
		S += Diag;
		result = S;
		
	}
	
	/**
	 * Returns the Haigh-Westergaard stress representation
	 */
	void HaighWestergaard(T &LodeAngle, T &qsi, T &rho) const;
	
	/**
	 * Returns the Haigh-Westergaard stress representation and derivatives
	 */
	void HaighWestergaard(T &LodeAngle, T &qsi, T &rho, TPZTensor<T> & dLodeAngle, TPZTensor<T> & dQsi, TPZTensor<T> & dRho) const;
	
	/**
	 * Returns the tensor eigenvalues through an analytical approach
	 */
	void EigenValue(TPZTensor<T> &eigenval)const;
	
	
	/**
	 * Returns the tensor eigenvalues and derivatives through an analytical approach
	 */
	void Eigenvalue(TPZTensor<T> &eigenval,TPZTensor<T> &dSigma1,TPZTensor<T> &dSigma2,TPZTensor<T> &dSigma3)const;
	
	/**
	 * Computes the Lode angle and its derivatives
	 */
	void Lodeangle(TPZTensor<T> &GradLode,T &Lode)const;
	
	
	/**
	 * Computes the eigenvectors and eigenvalues based on (page 742 Computational methods for computational plasticity/Souza Neto)
	 */
	void EigenSystem(TPZDecomposed &eigensystem) const;
	
	/**
	 * Computes the eigenvectors and eigenvalues based on Jacobi method
	 */
	void EigenSystemJacobi(TPZDecomposed &eigensystem) const;
    
   

	
	/**
	 Mnemonical access
	 */
	inline T & XX() const {return fData[_XX_];}
	inline T & XY() const {return fData[_XY_];}
	inline T & XZ() const {return fData[_XZ_];}
	inline T & YY() const {return fData[_YY_];}
	inline T & YZ() const {return fData[_YZ_];}
	inline T & ZZ() const {return fData[_ZZ_];}
	inline T &operator[](int i) const
	{
		return fData[i];
	}
	
private:
	/**
	 Metodo que calcula a diagonal do tensor deviatorico
	 @param [out] vec diagonal
	 */
	template <class TBASE>
	void DeviatoricDiagonal_T(TPZVec<T> & vec) const;
	
	// Specialization for TBASE = REAL
	void DeviatoricDiagonal(TPZVec<T> & vec) const;
	
public:
	/**
	 Copia os valores do tensor no tensor indicado
	 Derivadas sao decartadas
	 @param [out] target onde os valores serao copiados
	 */
	template < class T1 >
	void CopyTo(TPZTensor<T1> & target) const;
	
	/**
	 Copia os valores do tensor no vetor indicado
	 Derivadas sao decartadas
	 @param [out] target onde os valores serao copiados
	 */
	void CopyTo(TPZFMatrix<REAL> & target) const;
	
	/**
	 Copia os valores do vetor para o tensor
	 Derivadas sao zeradas
	 @param [in] source onde os valores serao copiados
	 */
	void CopyFrom(const TPZFMatrix<REAL> & source);
	
	/** Converts the stress vector onto a symmetric stress tensor
	 * @param Tensor [out]
	 */
	void CopyToTensor(TPZFMatrix<REAL> & Tensor);
	
	/**
	 Initializa o valor do tensor (tensor de deformacao)
	 */
	void SetUp(const TPZVec<REAL> & Solution);
	
public:
	/**
	 Dados do tensor
	 */
	TPZManVector<T,6> fData;
};

template <class T>
template <class TBASE>
void TPZTensor<T>::Identity_T(TBASE &)
{
	fData[_XX_] = TBASE(1.);
	fData[_YY_] = TBASE(1.);
	fData[_ZZ_] = TBASE(1.);
	
	fData[_XY_] = TBASE(0.);
	fData[_XZ_] = TBASE(0.);
	fData[_YZ_] = TBASE(0.);
}

template <class T>
void TPZTensor<T>::Identity()
{
	REAL TBase;
	Identity_T(TBase);
}

template < class T >
const TPZTensor<T> & TPZTensor<T>::operator=( const TPZTensor<T> &source )
{
	fData = source.fData;
	return *this;
}

template < class T >
const TPZTensor<T> & TPZTensor<T>::operator+=( const TPZTensor<T> &source )
{
	int i;
	for(i = 0; i < 6; i++)fData[i] += source.fData[i];
	return *this;
}

template < class T >
const TPZTensor<T> & TPZTensor<T>::operator-=( const TPZTensor<T> &source )
{
	int i;
	for(i = 0; i < 6; i++)fData[i] -= source.fData[i];
	return *this;
}

template < class T >
const TPZTensor<T> & TPZTensor<T>::operator*=( const T &multipl )
{
	int i;
	for(i = 0; i < 6; i++)fData[i] *= multipl;
	return *this;
}

template < class T >
template < class T1, class T2 >
void TPZTensor<T>::Add(const TPZTensor< T1 > & tensor, const T2 & constant )
{
	int i, size=6;
	for(i=0;i<size;i++){
		fData[i]+= tensor.fData[i] * T1(constant);
	}
}

template < class T >
template < class T1, class T2 >
void TPZTensor<T>::Multiply(const T1 & multipl, const T2 & constant )
{
	int i, size=6;
	for(i=0;i<size;i++){
		fData[i] *= multipl;
		fData[i] *= constant;
	}
}

template < class T >
void TPZTensor<T>::Multiply2(TPZTensor<T> tensor, TPZTensor<T> & resp)const
{
	T XX,YY,ZZ,XY,YZ,XZ;
	XX=this->fData[_XX_];
	YY=this->fData[_YY_];
	ZZ=this->fData[_ZZ_];
	XY=this->fData[_XY_];
	XZ=this->fData[_XZ_];
	YZ=this->fData[_YZ_];
	resp.fData[_XX_] = XX*tensor.XX()+XY*tensor.XY()+XZ*tensor.XZ();
	resp.fData[_YY_] = XY*tensor.XY()+YY*tensor.YY()+YZ*tensor.YZ();
	resp.fData[_ZZ_] = XZ*tensor.XZ()+YZ*tensor.YZ()+ZZ*tensor.ZZ();
	resp.fData[_XY_] = XX*tensor.XY()+XY*tensor.YY()+XZ*tensor.YZ();
	resp.fData[_XZ_] = XX*tensor.XZ()+XY*tensor.YZ()+XZ*tensor.ZZ();
	resp.fData[_YZ_] = XY*tensor.XZ()+YY*tensor.YZ()+YZ*tensor.ZZ();
	
}


template < class T >
template < class T2 >
void TPZTensor<T>::Scale(const T2 & constant )
{
	int i, size=6;
	for(i=0;i<size;i++){
		fData[i]*= constant;
	}
}

template < class T >
void TPZTensor<T>::Zero()
{
	int i, size=6;
	for(i=0;i<size;i++){
		fData[i] = 0.;
	}
}


template < class T >
void TPZTensor<T>::SetUp(const TPZVec<REAL> & Solution){
	int i, size=6;
	for(i=0;i<size;i++)
		fData[i]=Solution[i];
	
}
template < class T >
template < class T1 >
void TPZTensor<T>::CopyTo(TPZTensor<T1> & target) const {
	int i, size=6;
	for(i=0;i<size;i++){
		target[i]=shapeFAD::val(fData[i]);
	}
	
}

template < class T >
void TPZTensor<T>::CopyTo(TPZFMatrix<REAL> & target) const
{
	int i;
	for(i = 0 ; i < 6; i++)target(i,0) = shapeFAD::val(fData[i] );
}

template < class T >
void TPZTensor<T>::CopyFrom(const TPZFMatrix<REAL> & source)
{
	int i;
	for(i = 0; i < 6; i++)fData[i] = source.Get(i,0);
}

template < class T >
template < class TBASE >
T TPZTensor<T>::J2_T(TBASE &) const {
	TPZVec<T> s(3);
	DeviatoricDiagonal_T<TBASE>(s);
	return ((s[0]*s[0]+s[1]*s[1]+s[2]*s[2])/ TBASE(2.))+ fData[_XY_]*fData[_XY_]+fData[_XZ_]*fData[_XZ_]+fData[_YZ_]*fData[_YZ_];
}

template < class T >
T TPZTensor<T>::J2() const {
	REAL TBase;
	return J2_T(TBase);
}

template < class T >
T TPZTensor<T>::I1() const {
	return fData[_XX_] + fData[_YY_] + fData[_ZZ_];
}

template < class T >
T TPZTensor<T>::I2() const {
	return -(fData[_XY_] * fData[_XY_] +
					 fData[_XZ_] * fData[_XZ_] +
					 fData[_YZ_] * fData[_YZ_])
	+ (fData[_XX_] * fData[_YY_] +
	   fData[_YY_] * fData[_ZZ_] +
	   fData[_XX_] * fData[_ZZ_] );
}

template < class T >
T TPZTensor<T>::I3() const {
	return  fData[_XX_] * fData[_YY_] * fData[_ZZ_]
	+(fData[_XY_] * fData[_XZ_] * fData[_YZ_]) * 2.
	-(fData[_XX_] * fData[_YZ_] * fData[_YZ_]+
	  fData[_YY_] * fData[_XZ_] * fData[_XZ_]+
	  fData[_ZZ_] * fData[_XY_] * fData[_XY_]);
}

template < class T >
void TPZTensor<T>::dJ2(TPZTensor<T> & Tangent) const{
	T p = I1()/T(3.);
	Tangent.fData[_XX_] = fData[_XX_] - p;
	Tangent.fData[_YY_] = fData[_YY_] - p;
	Tangent.fData[_ZZ_] = fData[_ZZ_] - p;
	Tangent.fData[_XY_] = fData[_XY_] * T(2.);
	Tangent.fData[_XZ_] = fData[_XZ_] * T(2.);
	Tangent.fData[_YZ_] = fData[_YZ_] * T(2.);
	
}

template < class T >
template < class TBASE >
void TPZTensor<T>::DeviatoricDiagonal_T(TPZVec<T> & vec) const{
	
	T p = I1() / T(TBASE(3.));
	vec[0]=fData[_XX_]-p;
	vec[1]=fData[_YY_]-p;
	vec[2]=fData[_ZZ_]-p;
}

template < class T >
void TPZTensor<T>::DeviatoricDiagonal(TPZVec<T> & vec) const{
	
	DeviatoricDiagonal_T<REAL>(vec);
}


template < class T>
T TPZTensor<T>::Norm() const
{
	int i;
	T norm = T(0.);
	// norm	= sqrt(fData[_XX_]*fData[_XX_]+fData[_YY_]*fData[_YY_]+fData[_ZZ_]*fData[_ZZ_]+T(2.)*(fData[_XY_]*fData[_XY_])+T(2.)*(fData[_XZ_]*fData[_XZ_])+T(2.)*(fData[_YZ_]*fData[_YZ_]));
	for(i=0; i<6; i++) norm += fData[i]*fData[i];
	return sqrt(norm);
	//return norm;
}

template <class T>
std::ostream &operator<<(std::ostream &out,const TPZTensor<T> &tens)
{
	for (int i=0; i<6; i++) {
		out << tens[i] << " ";
	}
	return out;
}

/**
 * Metodo que calcula o determinante do tensor
 */
template <class T>
T TPZTensor<T>::Det() const
{
	return fData[_XX_]*fData[_YY_]*fData[_ZZ_]+fData[_XY_]*fData[_XZ_]*fData[_YZ_]*2.-fData[_XZ_]*fData[_YY_]*fData[_XZ_]-
	fData[_XY_]*fData[_XY_]*fData[_ZZ_]-fData[_YZ_]*fData[_YZ_]*fData[_XX_];
}

/**
 * Metodo que calcula o determinante do tensor
 */
template <class T>
void TPZTensor<T>::dDet(TPZTensor<T> &grad) const
{
	grad.fData[_XX_] = fData[_YY_]*fData[_ZZ_]-fData[_YZ_]*fData[_YZ_];
	grad.fData[_XY_] = fData[_XZ_]*fData[_YZ_]-fData[_XY_]*fData[_ZZ_];
	grad.fData[_XZ_] = fData[_XY_]*fData[_YZ_]-fData[_XZ_]*fData[_YY_];
	grad.fData[_YY_] = fData[_XX_]*fData[_ZZ_]-fData[_XZ_]*fData[_XZ_];
	grad.fData[_YZ_] = fData[_XY_]*fData[_XZ_]-fData[_YZ_]*fData[_XX_];
	grad.fData[_ZZ_] = fData[_XX_]*fData[_YY_]-fData[_XY_]*fData[_XY_];
}

template <class T>
void TPZTensor<T>::S(TPZTensor<T> &s) const
{
	s = *this;
	TPZTensor<T> one;
	one.Identity();
	T mult = -I1()/3.;
	s.fData[_XX_] += mult;
	s.fData[_YY_] += mult;
	s.fData[_ZZ_] += mult;
}

template <class T>
T TPZTensor<T>::J3() const
{
	TPZTensor<T> s;
	S(s);
	return s.Det();
}
/**
 * The derivative of the  3rd invariant of deviatoric tensor
 */


template<class T>
void TPZTensor<T>::dJ3(TPZTensor<T> &deriv) const
{
	T state0(fData[_XX_]),state1(fData[_XY_]),statex1(fData[_XY_]),
	state2(fData[_XZ_]),statex2(fData[_XZ_]),state3(fData[_YY_]),
	state4(fData[_YZ_]),statex4(fData[_YZ_]),state5(fData[_ZZ_]);
	deriv.fData[_XX_] = (fData[_XX_]*fData[_XX_]*2.)/9. + fData[_XY_]*fData[_XY_]/3. +
	fData[_XZ_]*fData[_XZ_]/3. - (fData[_XX_]*fData[_YY_]*2.)/9. -
	fData[_YY_]*fData[_YY_]/9. - (fData[_YZ_]*fData[_YZ_]*2.)/3. - (fData[_XX_]*fData[_ZZ_]*2.)/
	9. + (fData[_YY_]*fData[_ZZ_]*4.)/9. - fData[_ZZ_]*fData[_ZZ_]/9.;
	deriv.fData[_XY_] = (state0*state1)/3. + (state1*state3)/3. -
	(state1*state5*2.)/3. + (state0*statex1)/3. +
	(state3*statex1)/3. - (state5*statex1*2.)/3. +
	state4*statex2 + state2*statex4;
	deriv.fData[_XZ_] = (state0*state2)/3. - (state2*state3*2.)/3. +
	state1*state4 + (state2*state5)/3. +
	(state0*statex2)/3. - (state3*statex2*2.)/3. +
	(state5*statex2)/3. + statex1*statex4;
	deriv.fData[_YY_] = -(fData[_XX_]*fData[_XX_])/9. +
	(fData[_XY_]*fData[_XY_])/3. - ((fData[_XZ_]*fData[_XZ_])*2.)/3. - (fData[_XX_]*fData[_YY_]*2.)/9. + ((fData[_YY_]*fData[_YY_])*2.)/9. +
	(fData[_YZ_]*fData[_YZ_])/3. + (fData[_XX_]*fData[_ZZ_]*4.)/9. - (fData[_YY_]*fData[_ZZ_]*2.)/9. -
	(fData[_ZZ_]*fData[_ZZ_])/9.;
	deriv.fData[_YZ_] = (state0*state4*(-2.))/3. + (state3*state4)/3. +
	(state4*state5)/3. + state2*statex1 +
	state1*statex2 - (state0*statex4*2.)/3. +
	(state3*statex4)/3. + (state5*statex4)/3.;
	deriv.fData[_ZZ_] = -(fData[_XX_]*fData[_XX_])/9. - (2.*(fData[_XY_]*fData[_XY_]))/3. +
	(fData[_XZ_]*fData[_XZ_])/3. + (fData[_XX_]*fData[_YY_]*4.)/9. - (fData[_YY_]*fData[_YY_])/9. +
	(fData[_YZ_]*fData[_YZ_])/3. - (fData[_XX_]*fData[_ZZ_]*2.)/9. - (fData[_YY_]*fData[_ZZ_]*2.)/
	9. + ((fData[_ZZ_]*fData[_ZZ_])*2.)/9.;
	
}

template <class T>
void TPZTensor<T>::CopyToTensor(TPZFMatrix<REAL> & Tensor)
{
	Tensor.Resize(3,3);
	Tensor(0,0) = shapeFAD::val( XX() ); //xx
	Tensor(0,1) = shapeFAD::val( XY() ); //xy
	Tensor(0,2) = shapeFAD::val( XZ() ); //xz
	
	Tensor(1,0) = shapeFAD::val( XY() ); //xy
	Tensor(1,1) = shapeFAD::val( YY() ); //yy
	Tensor(1,2) = shapeFAD::val( YZ() ); //yz
	
	Tensor(2,0) = shapeFAD::val( XZ() ); //xz
	Tensor(2,1) = shapeFAD::val( YZ() ); //yz
	Tensor(2,2) = shapeFAD::val( ZZ() ); //zz
}

template<class T>
static T Polynom(const T &x, const T &I1, const T &I2, const T &I3)
{
	T result = x*x*x- I1*x*x + I2*x - I3;
	return result;
}

template<class T>
static T DericPolynom(const T &x, const T &I1, const T &I2, const T &I3)
{
	T result = T(3.)*x*x- T(2.)*I1*x + I2;
	return result;
}

template<class T>
static T UpdateNewton(const T &x, const T &I1, const T &I2, const T &I3)
{
	T residue = Polynom(x, I1, I2, I3);
	T dres = DericPolynom(x, I1, I2, I3);
	T result = x - residue/dres;
	return result;
}

template <class T>
void TPZTensor<T>::EigenSystemJacobi(TPZDecomposed &eigensystem)const
{
	TPZFNMatrix<9,T> TensorMat(*this);
	
	long numiterations = 1000;
	REAL tol = 1.e-6;
	TPZManVector<T,3> EigenvaluesJac(3,0.);
	TPZFNMatrix<9,T> EigenvecJac(3,3,0.);
	TensorMat.SolveEigensystemJacobi(numiterations, tol, EigenvaluesJac, EigenvecJac);
	
	TPZManVector<T,3> &Eigenvalues1 = eigensystem.fEigenvalues;
	TPZManVector<TPZTensor<T>, 3> &Eigenvectors1 = eigensystem.fEigenvectors;
	Eigenvalues1 = EigenvaluesJac;
	for (int k = 0; k < 3; k++) {
		Eigenvectors1[k].XX() = EigenvecJac(0,k) * EigenvecJac(0,k);
		Eigenvectors1[k].XY() = EigenvecJac(0,k) * EigenvecJac(1,k);
		Eigenvectors1[k].XZ() = EigenvecJac(0,k) * EigenvecJac(2,k);
		Eigenvectors1[k].YY() = EigenvecJac(1,k) * EigenvecJac(1,k);
		Eigenvectors1[k].YZ() = EigenvecJac(1,k) * EigenvecJac(2,k);
		Eigenvectors1[k].ZZ() = EigenvecJac(2,k) * EigenvecJac(2,k);
	}
#ifdef LOG4CXX
	{
		std::stringstream str;
		str << "\n-------------AUTOVETORES JACOBI--------------" << std::endl;
		str << "Tensor:" << std::endl;
		this->Print(str);
		str << "Eigenvalues = " << EigenvaluesJac << std::endl;
		str << "Eigenvectors:" << std::endl;
		EigenvecJac.Print("Eigenvec:",str);
		str << "Eigenprojections:" << std::endl;
		for (int k = 0; k < 3; k++) {
			str << "EigenProjection " << k << ":" << std::endl;
			Eigenvectors1[k].Print(str);
		}
		
		str << "\n-------------FIM AUTOVETORES JACOBI--------------" << std::endl;
		LOGPZ_DEBUG(loggerr,str.str())
	}
#endif
}

template <class T>
void TPZTensor<T>::EigenSystem(TPZDecomposed &eigensystem)const
{
	
	T I1,I2,I3,R,theta,Q,x1,x2,x3,costheta,e1temp,e2temp,e3temp;
    TPZTensor<T> InternalTensor(*this);
    T norm =(InternalTensor.Norm());
    if(norm>0.&&norm>10.)
    {
        InternalTensor*=(T(1.)/norm);
    }
	I1 = (InternalTensor.I1());
	I2 = (InternalTensor.I2());
	I3 = (InternalTensor.I3());
	
    T toll=1.e-12;
	
	R=(T(-2.)*I1*I1*I1+T(9.)*I1*I2-T(27.)*I3)/T(54.);
	Q=(I1*I1-T(3.)*I2)/T(9.);
	if (Q < T(0.) || Q < toll) {
		Q = T(0.);
	}
    
	T denom = sqrt(Q*Q*Q);
	
	REAL Rval = shapeFAD::val(R);
	REAL denomval = shapeFAD::val(denom);
	TPZManVector<T,3> &Eigenvalues = eigensystem.fEigenvalues;
	TPZManVector<TPZTensor<T>, 3> &Eigenvectors = eigensystem.fEigenvectors;
	
	if (fabs(Rval) < toll && fabs(denomval) < toll && fabs(denomval) <= fabs(Rval)) {
		// three equal eigenvalues
		x1 = (this->I1())/3.;
		x2 = x1;
		x3 = x1;
        
		//Eigenvectors[0].Identity();
		//Eigenvectors[0].Scale(T(1/3.));
		//Eigenvectors[1] = Eigenvectors[0];
		//Eigenvectors[2] = Eigenvectors[1];
        
		Eigenvectors[0].XX() = 1.;
		Eigenvectors[1].YY() = 1.;
		Eigenvectors[2].ZZ() = 1.;
		
		Eigenvalues[0]=x1;
		Eigenvalues[1]=x2;
		Eigenvalues[2]=x3;
		return;
	}
	
    I1 = (this->I1());
	I2 = (this->I2());
	I3 = (this->I3());
	
	R=(T(-2.)*I1*I1*I1+T(9.)*I1*I2-T(27.)*I3)/T(54.);
	Q=(I1*I1-T(3.)*I2)/T(9.);
    
	if (Q < T(0.) || Q < toll) {
		Q = T(0.);
	}
    denom = sqrt(Q*Q*Q);
	
    Rval = shapeFAD::val(R);
	denomval = shapeFAD::val(denom);
	Eigenvalues = eigensystem.fEigenvalues;
	Eigenvectors = eigensystem.fEigenvectors;
    
	costheta=R/denom;
	
	if(shapeFAD::val(costheta) < -(1.-1.e-12))
	{
		costheta = -1.;
		theta = M_PI;
		x1 = T(2.)*sqrt(Q)+I1/T(3.);
		x1 = UpdateNewton(x1, I1, I2, I3);
		T B = x1-I1;
		x2 = -B/2.;
		x3 = -B/2.;
	}
	else if(shapeFAD::val(costheta)>(1.-1.e-12))
	{
		costheta = 1.;
		theta = 0.;
		x1 = T(-2.)*sqrt(Q)+I1/T(3.);
		x1 = UpdateNewton(x1, I1, I2, I3);
		T B = x1-I1;
		x2 = -B/2.;
		x3 = -B/2.;
	}
	else {
		theta = acos(costheta);
		x1=T(-2.)*sqrt(Q)*cos(theta/T(3.))+ (I1/T(3.));//eigenval 1
		x2=T(-2.)*sqrt(Q)*cos( ( theta +  T(M_PI*2.) )/3)+ (I1/T(3.));//eigenval 2
		x3=T(-2.)*sqrt(Q)*cos( ( theta -  T(M_PI*2.) )/3)+ (I1/T(3.));//eigenval 3
	}
	
	
	REAL valx1 = shapeFAD::val(x1), valx2 = shapeFAD::val(x2), valx3 = shapeFAD::val(x3);
	if(valx1 < valx2)
	{
		T temp = x1;
		x1 = x2;
		x2 = temp;
		REAL valtemp = valx1;
		valx1 = valx2;
		valx2 = valtemp;
	}
	if(valx1 < valx3)
	{
		T temp = x1;
		x1 = x3;
		x3 = temp;
		REAL valtemp = valx1;
		valx1 = valx3;
		valx3 = valtemp;
	}
	if(valx2 < valx3)
	{
		T temp = x2;
		x2 = x3;
		x3 = temp;
		REAL valtemp = valx2;
		valx2 = valx3;
		valx3 = valtemp;
	}
	REAL tolerance = shapeFAD::val(toll);
    tolerance = 5.e-5;
	if(valx1-valx2 > tolerance && valx2-valx3 > tolerance)
	{
		Eigenvectors.resize(3);
		TPZTensor<T> sqrsigma,sqrsigma2,sqrsigma3,sigmacopy(*this),sigmacopy2(*this),sigmacopy3(*this),Identity,Identity2,Identity3;
		Identity.XX()=T(1.);Identity.YY()=T(1.);Identity.ZZ()=T(1.);
		Identity2.XX()=T(1.);Identity2.YY()=T(1.);Identity2.ZZ()=T(1.);
		Identity3.XX()=T(1.);Identity3.YY()=T(1.);Identity3.ZZ()=T(1.);
		Multiply2(*this,sqrsigma);
		sqrsigma2 = sqrsigma;
		sqrsigma3 = sqrsigma;
		
		e1temp=x1/(T(2.)*x1*x1*x1-I1*x1*x1+I3);
		sigmacopy.Multiply(I1-x1,1);
		Identity.Multiply((I3/x1),1);
		sqrsigma.Add(sigmacopy,-1);
		sqrsigma.Add(Identity,1);
		sqrsigma.Multiply(e1temp,1);
		Eigenvectors[0]=sqrsigma;
		
		e2temp=x2/(T(2.)*x2*x2*x2-I1*x2*x2+I3);
		sigmacopy2.Multiply(I1-x2,1);
		Identity2.Multiply((I3/x2),1);
		sqrsigma2.Add(sigmacopy2,-1);
		sqrsigma2.Add(Identity2,1);
		sqrsigma2.Multiply(e2temp,1);
		Eigenvectors[1]=sqrsigma2;
		
		e3temp=x3/(T(2.)*x3*x3*x3-I1*x3*x3+I3);
		sigmacopy3.Multiply(I1-x3,1);
		Identity3.Multiply((I3/x3),1);
		sqrsigma3.Add(sigmacopy3,-1);
		sqrsigma3.Add(Identity3,1);
		sqrsigma3.Multiply(e3temp,1);
		Eigenvectors[2]=sqrsigma3;
	}
	else if(valx1-valx2 >  tolerance)
	{
		TPZTensor<T> sqrsigma,sqrsigma2,sqrsigma3,sigmacopy(*this),sigmacopy2(*this),sigmacopy3(*this),Identity,Identity2,Identity3;
		Identity.XX()=T(1.);Identity.YY()=T(1.);Identity.ZZ()=T(1.);
		Identity2.XX()=T(1.);Identity2.YY()=T(1.);Identity2.ZZ()=T(1.);
		Identity3.XX()=T(1.);Identity3.YY()=T(1.);Identity3.ZZ()=T(1.);
		Multiply2(*this,sqrsigma);
		sqrsigma2 = sqrsigma;
		sqrsigma3 = sqrsigma;
		
		e1temp=x1/(T(2.)*x1*x1*x1-I1*x1*x1+I3);
		sigmacopy.Multiply(I1-x1,1);
		Identity.Multiply((I3/x1),1);
		sqrsigma.Add(sigmacopy,-1);
		sqrsigma.Add(Identity,1);
		sqrsigma.Multiply(e1temp,1);
		
		x2 = (I1-x1)/T(2.);
		x3 = x2;
		
		Eigenvectors[0] = sqrsigma;
		
		// NATHAN ALTEROU 22/01/2014
		// O metodo soh calcula o autovetor referente ao autovalor diferente
		// Aqui calculo outros dois autovetores quaisquer que sao ortogonais ao primeiro e entre si
		int posorig = 0;
		TPZFNMatrix<3,T> epsegveFromProj(3,3,0.);
		TPZFNMatrix<9,T> EigenvecMat(3,3,0.);     
		EigenvecMat = Eigenvectors[posorig];
		for	(int k = 0 ; k < 3 ; k++){
			bool IsnoZero = false;
			for (int j = 0; j < 3; j++) {
				epsegveFromProj(j,0) = EigenvecMat(j,k);
				if (!IsZero(epsegveFromProj(j,0))) IsnoZero = true;
			}	
			if (IsnoZero) break;
		}
		
		T normvec = 0.;
		for (int i = 0; i < epsegveFromProj.Rows(); i++) {
			normvec += epsegveFromProj(i,0) * epsegveFromProj(i,0);
		}
		normvec = sqrt(normvec);
		for (int j = 0; j < 3; j++) {
			epsegveFromProj(j,0) /= normvec;
		}
		
		// Finding the smaler 2 values		
		T minabs = 0.;
		int pos = 0, pos2 = -1;
		minabs = epsegveFromProj(0,0);
		for (int i = 1; i < 3; i++) {
			T temp = epsegveFromProj(i,0);
			if (fabs(temp) < fabs(minabs)){
				minabs = temp;
				pos = i;
			} 
		}
		int k = 0;
		TPZManVector<REAL,2> LastingPos(2,-1.);
		for (int i = 0; i < 3; i++) {
			if (i == pos) continue;
			LastingPos[k++] = i;
		}
		pos2 = (fabs(epsegveFromProj(LastingPos[0],0)) < fabs(epsegveFromProj(LastingPos[1],0)) ) ? LastingPos[0] : LastingPos[1]; 
		epsegveFromProj.Resize(3, 3);
		epsegveFromProj(pos,1) = 1.;
		epsegveFromProj(pos2,2) = 1.;
		
		TPZFNMatrix <9,T> Orthog(3,3,0.),TransfToOrthog(3,3,0.);
		
		epsegveFromProj.GramSchmidt(Orthog, TransfToOrthog);
		
		//std::cout << "Autovetores:" << std::endl;
		//epsegveFromProj.Print("Eigenvectors antigo");
		//Orthog.Print("Eigenvectors Ortogonalizado");
		
		int rightvec[3] = {posorig,1,2};
		int i = -1;
		for (int k = 0; k < 3; k++) {
			i = rightvec[k];
			Eigenvectors[i].XX() = Orthog(0,k)*Orthog(0,k);
			Eigenvectors[i].XY() = Orthog(0,k)*Orthog(1,k);
			Eigenvectors[i].XZ() = Orthog(0,k)*Orthog(2,k);
			Eigenvectors[i].YY() = Orthog(1,k)*Orthog(1,k);
			Eigenvectors[i].YZ() = Orthog(1,k)*Orthog(2,k);
			Eigenvectors[i].ZZ() = Orthog(2,k)*Orthog(2,k);
			//std::cout << "i = " << i << std::endl;
			//Eigenvectors[i].Print(std::cout);
		}

#ifdef DEBUG
		{
			// Faz Transpose[Eigenvectors].DiagonalMatrix[Eigenvalues].Eigenvectors e verifica se recupero o tensor
			TPZFNMatrix <9,T> EgVecRightOrder(3,3,0.),EgVecRightOrderT(3,3,0.), DiagMat(3,3,0.),TempMat(3,3,0.);
			int rightvec[3] = {posorig,1,2};
			int i = -1;
			for (int k = 0; k < 3; k++) {
				i = rightvec[k];
				for (int j = 0; j < 3; j++) {
					EgVecRightOrder(j,i) = Orthog(j,k);
				}
			}
			EgVecRightOrder.Transpose(&EgVecRightOrderT);
			
			// montando Diagonal matrix
			TPZManVector<T,3> EgVaTemp(3,0.);
			EgVaTemp[0] = x1;
			EgVaTemp[1] = x2;
			EgVaTemp[2] = x3;
			
			int k = 0;
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3 ; j++) {
					if (i == j) DiagMat(i,j) = EgVaTemp[k++]; 
				}
			}
			
			DiagMat.Multiply(EgVecRightOrderT,TempMat);
			EgVecRightOrder.Multiply(TempMat,DiagMat); // DiagMat guardara o resultado
			TPZFNMatrix<9,T> OrigTensor(3,3,0.);
			OrigTensor = *this;
			
			REAL tol = tolerance*100.;
			
			bool IsEqual = true;
			T diff = 0.;
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					diff = fabs(OrigTensor(i,j) - DiagMat(i,j));
					if (diff > tol) {
						IsEqual = false;
					}
				}
			}
			
			
			if (!IsEqual){
				std::cout << "Nao sao Iguais!!!" << std::endl;
				OrigTensor.Print("Tensor Original",std::cout);
				DiagMat.Print("Tensor Reconstruido:",std::cout);
				DebugStop();				
			} 
		}
#endif

	}
	else if(valx2 - valx3 > tolerance)
	{
		TPZTensor<T> sqrsigma,sqrsigma2,sqrsigma3,sigmacopy(*this),sigmacopy2(*this),sigmacopy3(*this),Identity,Identity2,Identity3;
		Identity.XX()=T(1.);Identity.YY()=T(1.);Identity.ZZ()=T(1.);
		Identity2.XX()=T(1.);Identity2.YY()=T(1.);Identity2.ZZ()=T(1.);
		Identity3.XX()=T(1.);Identity3.YY()=T(1.);Identity3.ZZ()=T(1.);
		Multiply2(*this,sqrsigma);
		sqrsigma2 = sqrsigma;
		sqrsigma3 = sqrsigma;
		
		e1temp=x3/(T(2.)*x3*x3*x3-I1*x3*x3+I3);
		sigmacopy.Multiply(I1-x3,1);
		Identity.Multiply((I3/x3),1);
		sqrsigma.Add(sigmacopy,-1);
		sqrsigma.Add(Identity,1);
		sqrsigma.Multiply(e1temp,1);
		
		x1 = (I1-x3)/T(2.);
		x2 = x1;
		Eigenvectors[2]=sqrsigma;
				
		// NATHAN ALTEROU 22/01/2014
		// O metodo soh calcula o autovetor referente ao autovalor diferente
		// Aqui calculo outros dois autovetores quaisquer que sao ortogonais ao primeiro e entre si
		int posorig = 2;
		TPZFNMatrix<3,T> epsegveFromProj(3,3,0.);
		TPZFNMatrix<9,T> EigenvecMat(3,3,0.);     
		EigenvecMat = Eigenvectors[posorig];
		for	(int k = 0 ; k < 3 ; k++){
			bool IsnoZero = false;
			for (int j = 0; j < 3; j++) {
				epsegveFromProj(j,0) = EigenvecMat(j,k);
				if (!IsZero(epsegveFromProj(j,0))) IsnoZero = true;
			}	
			if (IsnoZero) break;
		}

		T normvec = 0.;
		for (int i = 0; i < epsegveFromProj.Rows(); i++) {
			normvec += epsegveFromProj(i,0) * epsegveFromProj(i,0);
		}
		normvec = sqrt(normvec);
		for (int j = 0; j < 3; j++) {
			epsegveFromProj(j,0) /= normvec;
		}
		
		// Finding the smaler 2 values		
		T minabs = 0.;
		int pos = 0, pos2 = -1;
		minabs = epsegveFromProj(0,0);
		for (int i = 1; i < 3; i++) {
			T temp = epsegveFromProj(i,0);
			if (fabs(temp) < fabs(minabs)){
				minabs = temp;
				pos = i;
			} 
		}
		int k = 0;
		TPZManVector<REAL,2> LastingPos(2,-1.);
		for (int i = 0; i < 3; i++) {
			if (i == pos) continue;
			LastingPos[k++] = i;
		}
		pos2 = (fabs(epsegveFromProj(LastingPos[0],0)) < fabs(epsegveFromProj(LastingPos[1],0)) ) ? LastingPos[0] : LastingPos[1]; 
		epsegveFromProj.Resize(3, 3);
		epsegveFromProj(pos,1) = 1.;
		epsegveFromProj(pos2,2) = 1.;
		
		TPZFNMatrix <9,T> Orthog(3,3,0.),TransfToOrthog(3,3,0.);
		
		epsegveFromProj.GramSchmidt(Orthog, TransfToOrthog);
		
		//std::cout << "Autovetores:" << std::endl;
		//epsegveFromProj.Print("Eigenvectors antigo");
		//Orthog.Print("Eigenvectors Ortogonalizado");
		
		int rightvec[3] = {posorig,0,1};
		int i = -1;
		for (int k = 0; k < 3; k++) {
			i = rightvec[k];
			Eigenvectors[i].XX() = Orthog(0,k)*Orthog(0,k);
			Eigenvectors[i].XY() = Orthog(0,k)*Orthog(1,k);
			Eigenvectors[i].XZ() = Orthog(0,k)*Orthog(2,k);
			Eigenvectors[i].YY() = Orthog(1,k)*Orthog(1,k);
			Eigenvectors[i].YZ() = Orthog(1,k)*Orthog(2,k);
			Eigenvectors[i].ZZ() = Orthog(2,k)*Orthog(2,k);
			//std::cout << "i = " << i << std::endl;
			//Eigenvectors[i].Print(std::cout);
		}

#ifdef DEBUG
		{
			// Faz Transpose[Eigenvectors].DiagonalMatrix[Eigenvalues].Eigenvectors e verifica se recupero o tensor
			TPZFNMatrix <9,T> EgVecRightOrder(3,3,0.),EgVecRightOrderT(3,3,0.), DiagMat(3,3,0.),TempMat(3,3,0.);
			int rightvec[3] = {posorig,0,1};
			int i = -1;
			for (int k = 0; k < 3; k++) {
				i = rightvec[k];
				for (int j = 0; j < 3; j++) {
					EgVecRightOrder(j,i) = Orthog(j,k);
				}
			}
			EgVecRightOrder.Transpose(&EgVecRightOrderT);

			// montando Diagonal matrix
			TPZManVector<T,3> EgVaTemp(3,0.);
			EgVaTemp[0] = x1;
			EgVaTemp[1] = x2;
			EgVaTemp[2] = x3;
			
			int k = 0;
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3 ; j++) {
					if (i == j) DiagMat(i,j) = EgVaTemp[k++]; 
				}
			}
			
			DiagMat.Multiply(EgVecRightOrderT,TempMat);
			EgVecRightOrder.Multiply(TempMat,DiagMat); // DiagMat guardara o resultado
			TPZFNMatrix<9,T> OrigTensor(3,3,0.);
			OrigTensor = *this;
			
			REAL tol = tolerance*100.;
			
			bool IsEqual = true;
			T diff = 0.;
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					diff = fabs(OrigTensor(i,j) - DiagMat(i,j));
					if (diff > tol) {
						IsEqual = false;
					}
				}
			}
			
			
			if (!IsEqual){
				std::cout << "Nao sao Iguais!!!" << std::endl;
				OrigTensor.Print("Tensor Original",std::cout);
				DiagMat.Print("Tensor Reconstruido:",std::cout);
				DebugStop();				
			} 
		}
		
#endif

	}
	else {
		
		x1 = I1/T(3.);
		x2 = x1;
		x3 = x1;
    //    Eigenvectors[0].Identity();
    //    Eigenvectors[0].Scale(T(1/3.));
    //    Eigenvectors[1] = Eigenvectors[1];
    //    Eigenvectors[2] = Eigenvectors[2];
		
		Eigenvectors[0].XX() = 1.;
		Eigenvectors[1].YY() = 1.;
		Eigenvectors[2].ZZ() = 1.;
		
		Eigenvalues[0]=x1;
		Eigenvalues[1]=x2;
		Eigenvalues[2]=x3;
	}
	Eigenvalues[0]=x1;
	Eigenvalues[1]=x2;
	Eigenvalues[2]=x3;
	
}

template <class T>
void TPZTensor<T>::HaighWestergaard(T &LodeAngle, T &qsi, T &rho) const
{/*
  T I1(this->I1()),
  J2(this->J2()),
  J3(this->J3());
  T sqrtJ2 = sqrt(J2);
  
  LodeAngle = acos( J3/J2/sqrtJ2*T( sqrt(27.)/2.) ) / T(3.);
  qsi = I1 / T( sqrt(3.) );
  rho = sqrtJ2 * T( sqrt(2./3.) );
  */
}

template <class T>
void TPZTensor<T>::HaighWestergaard(T &LodeAngle, T &qsi, T &rho,
																		TPZTensor<T> & dLodeAngle, TPZTensor<T> & dQsi, TPZTensor<T> & dRho) const
{
	//	T I1(this->I1()),
	T J2(this->J2()), J3(this->J3());
	T sqrtJ2 = sqrt(J2);
	
	TPZTensor<T> dJ2, dJ3;
	
	this->dJ2(dJ2);
	this->dJ3(dJ3);
	
	// Derivatives with respect to I1, J2 and J3
	T dLodeAngledJ2,
	dLodeAngledJ3;
	
	TPZTensor<T> TempTensor;
	
	LodeAngle = acos( J3/J2/sqrtJ2*T( sqrt(27.)/2.) ) / T(3.);
	//	T pi23 = T(2. * M_PI / 3.);
	
	T temp = sqrt( T(3.) / ( J2*J2*J2 * T(4.) - J3*J3 * T(27.)) );
	
	dLodeAngledJ2 =  temp * J3 / J2 * T(3./2.);
	dLodeAngledJ3 = -temp;
	
	//dLodeAngleJ2 = dLodeAngledJ2*dJ2 +
	//			 dLodeAngledJ3*dJ3;
	TempTensor = dJ2;
	TempTensor *= dLodeAngledJ2;
	dLodeAngle = dJ3;
	dLodeAngle *= dLodeAngledJ3;
	dLodeAngle += TempTensor;
	/*
	 T I1(this->I1()),
	 J2(this->J2()),
	 J3(this->J3());
	 T sqrtJ2 = sqrt(J2);
	 
	 TPZTensor<T> dJ2, dJ3;
	 
	 this->dJ2(dJ2);
	 this->dJ3(dJ3);
	 
	 // Derivatives with respect to I1, J2 and J3
	 T dLodeAngledJ2,
	 dLodeAngledJ3,
	 dQsidI1,
	 dRhodJ2;
	 
	 LodeAngle = acos( J3/J2/sqrtJ2*T( sqrt(27.)/2.) ) / T(3.);
	 qsi = I1 / T( sqrt(3.) );
	 rho = sqrtJ2 * T( sqrt(2./3.) );
	 
	 T temp = sqrt( T(3.) / ( J2*J2*J2 * T(4.) - J3*J3 * T(27.)) );
	 
	 dLodeAngledJ2 =  temp * J3 / J2 * T(3./2.);
	 dLodeAngledJ3 = -temp;
	 
	 dQsidI1 = T(1/sqrt(3.));
	 
	 dRhodJ2 = T(1/sqrt(2)) / sqrtJ2;
	 
	 dLodeAngle = dLodeAngledJ2*dJ2 +
	 dLodeAngledJ3*dJ3;
	 
	 dQsi = dQsidI1*dI1;
	 
	 dRho = dRhodJ2*dJ2;
	 */
}


template <class T>
void TPZTensor<T>::EigenValue(TPZTensor<T> &eigenval)const
{
	T I1(this->I1()),
	I2(this->I2()),
	I3(this->I3());
	
	//    if(I1<T(1.e-6))I1=T(1.e-6);
	//    if(I2<T(1.e-6))I2=T(1.e-6);
	//    if(I3<T(1.e-6))I3=T(1.e-6);
	
	T R,THETA,Q,temp,temp2;
	//    T verif;
	R = ( T(-2.)*(I1*I1*I1)+ T(9.)*I1*I2 - T(27.)* I3 )/54.;
	Q = ( (I1*I1) - (T(3.)*I2) )/T(9.);
	temp2 = (sqrt(Q*Q*Q));
	if(fabs(temp2)<1.e-6)temp2=1.e-6;
	temp = R/temp2;
	if(temp <= T(-1.))
	{
		temp*=0.9999999999;
	}
	if(temp >= T(1.))
	{
		temp*=0.9999999999;
	}
	THETA = acos(temp);
	
	if(THETA >= T(2*M_PI) || THETA<= T(2*M_PI))
	{
		THETA*=0.999999999;
	}
	eigenval.XX()= (T(-2.)*sqrt(Q)*cos(THETA/T(3.))) + I1/T(3.);
	eigenval.YY()= (T(-2.)*sqrt(Q)*cos((THETA+T(2.*M_PI))/T(3.))) + I1/T(3.);
	eigenval.ZZ()= (T(-2.)*sqrt(Q)*cos((THETA-T(2.*M_PI))/T(3.))) + I1/T(3.);
}

//template <class T>
//void TPZTensor<T>::EigenValue(TPZTensor<T> &eigenval,TPZTensor<T> &dSigma1,TPZTensor<T> &dSigma2,TPZTensor<T> &dSigma3)const
//{
//	T I1(this->I1()),
//	J2(this->J2()),
//	J3(this->J3());
//
//	if(fabs( shapeFAD::val(J2) ) < 1.e-6)J2 += 1.e-6;
//
//	T sqrtJ2 = sqrt(J2);
//
//	TPZTensor<T> dJ2, dJ3;
//
//	this->dJ2(dJ2);
//	this->dJ3(dJ3);
//
//	// Derivatives with respect to I1, J2 and J3
//	T dLodeAngledJ2,
//	dLodeAngledJ3;
//
//	TPZTensor<T> dLodeAngle, TempTensor;
//
//
//	T den = sqrtJ2*T( sqrt(27.)/2.);
//	T lodetemp = (J3/J2)/den;
//
//	if(shapeFAD::val(lodetemp) <= -0.999999)
//	{
//		lodetemp *= T(0.999);
//	}
//
//	if(shapeFAD::val(lodetemp) >= 0.999999)
//	{
//		lodetemp *= T(0.999);
//		DebugStop();
//	}
//
//	T acoslodetemp = acos(lodetemp);
//
//	T LodeAngle = acoslodetemp / T(3.);
//
//
//	T pi23 = T(2. * M_PI / 3.);
//
//	T denominador = ( J2*J2*J2 * T(4.) - J3*J3 * T(27.));
//	T temp;
//	if(shapeFAD::val(denominador) < 1.e-6)
//	{
//		dLodeAngledJ2 =  0.;
//	}
//	else
//	{
//		temp = sqrt( T(3.) / denominador );
//		dLodeAngledJ2 =  temp * J3 / J2 * T(3./2.);
//		dLodeAngledJ3 = -temp;
//	}
//
//
//
//	//dLodeAngleJ2 = dLodeAngledJ2*dJ2 +
//	//			 dLodeAngledJ3*dJ3;
//	TempTensor = dJ2;
//	TempTensor *= dLodeAngledJ2;
//	dLodeAngle = dJ3;
//	dLodeAngle *= dLodeAngledJ3;
//	dLodeAngle += TempTensor;
//
//	T TwoOverSqrThree = T(2./sqrt(3.));
//	T TwoOverSqrThreeJ2 = TwoOverSqrThree * sqrtJ2;
//	T I13 = I1 / T(3.);
//
//	T tempCosLode = 		cos(LodeAngle) 		  * TwoOverSqrThreeJ2;
//	T tempCosMinusLode =	cos(pi23 - LodeAngle) * TwoOverSqrThreeJ2;
//	T tempCosPlusLode =		cos(pi23 + LodeAngle) * TwoOverSqrThreeJ2;
//
//
//#ifdef LOG4CXX_PLASTICITY
//	{
//		std::stringstream sout;
//		sout << "\n  TPZTENSOR \n"<<endl;
//		sout << "\n  I1 = \n"<< I1 <<endl;
//		sout << "\n  J2 = "<<J2<<endl;
//		sout << "\n  J3 = "<<J3<<endl;
//		sout << "\n  dJ2 = "<<dJ2<<endl;
//		sout << "\n  dJ3 = "<<dJ3<<endl;
//		sout << "\n  LodeAngle= "<< LodeAngle<<endl;
//		sout << "\n  denominador = "<<denominador<<endl;
//		sout << "\n  den= "<< den << endl;
//		sout << "\n  lodetemp= "<< lodetemp << endl;
//		sout << "\n  tempCosLode= "<<tempCosLode<<endl;
//		sout << "\n tempCosMinusLode= "<<tempCosMinusLode<<endl;
//		sout << "\n  tempCosPlusLode= "<<tempCosPlusLode<<endl;
//		sout << "\n\n";
//		LOGPZ_INFO(loggerr,sout.str());
//	}
//#endif
//	eigenval.XX() = I13 + tempCosLode;
//	eigenval.YY() = I13 + tempCosMinusLode;
//	eigenval.ZZ() = I13 + tempCosPlusLode;
//	eigenval.XY() = T(0.);
//	eigenval.XZ() = T(0.);
//	eigenval.YZ() = T(0.);
//
//	T OneOverTwoJ2 = T(0.5) / J2;
//	TPZTensor<T> dI13;
//	dI13.Identity();
//	dI13 *= T(1./3.);
//
//	tempCosLode      *= OneOverTwoJ2;
//	tempCosMinusLode *= OneOverTwoJ2;
//	tempCosPlusLode  *= OneOverTwoJ2;
//
//	dSigma1 = dJ2;
//	dSigma1 *= tempCosLode;
//	dSigma1 += dI13;
//	TPZTensor<T> dLodeAngleTemp(dLodeAngle);
//	dLodeAngleTemp *= sin(LodeAngle) * TwoOverSqrThreeJ2;
//	dSigma1 -= dLodeAngleTemp;
//
//	dSigma2 = dJ2;
//	dSigma2 *= tempCosMinusLode;
//	dSigma2 += dI13;
//	dLodeAngleTemp = dLodeAngle;
//	dLodeAngleTemp *= sin(pi23 - LodeAngle) * TwoOverSqrThreeJ2;
//	dSigma2 += dLodeAngleTemp;
//
//	dSigma3 = dJ2;
//	dSigma3 *= tempCosPlusLode;
//	dSigma3 += dI13;
//	dLodeAngleTemp = dLodeAngle;
//	dLodeAngleTemp *= sin(pi23 + LodeAngle) * TwoOverSqrThreeJ2;
//	dSigma3 -= dLodeAngleTemp;
//////////////////////////////
//
//
//
//}

template <class T>
void TPZTensor<T>::Lodeangle(TPZTensor<T> &GradLode,T &Lode)const
{
	T J2t(this->J2()),
	J3t(this->J3());
	
	if(fabs( shapeFAD::val(J2t) ) < 1.e-6)J2t=T(1.e-6);
	T sqrtJ2t =sqrt(J2t);
	if(fabs( shapeFAD::val(sqrtJ2t) ) < 1.e-6)sqrtJ2t=T(1.e-6);
	
	TPZTensor<T> dJ2t, dJ3t;
	
	this->dJ2(dJ2t);
	this->dJ3(dJ3t);
	// Derivatives with respect to I1, J2 and J3
	
	TPZTensor<T> dLodeAngle, TempTensor;
	
	//QUAL DOS DOIS?
	//T theta =-asin( ( T( 3.) * sqrt( T( 3.) ) * J3t ) /( T( 2.) *  sqrt(J2t*J2t*J2t) ) )/T( 3.);
	//O GRADIENTE DO LODE ESTA EM FUNCAO DO LODE DE 0 a Pi/3
	
	
	T lodetemp =( T( 3.) * sqrt( T( 3.) ) * J3t ) /( T( 2.) *  sqrt(J2t*J2t*J2t) ) ;
	if(shapeFAD::val(lodetemp) <= -1.)
	{
		lodetemp *= T(0.999);	//	DebugStop();
		// shapeFAD::val(lodetemp) *= T(0.999);	//	DebugStop();
		//lodetemp = T(-1.);
		
	}
	
	
	//cout << "\n lodetemp "<<lodetemp<<endl;
	//cout << "\n shapeFAD::val(lodetemp);"<<shapeFAD::val(lodetemp)<<endl;
	
	if(shapeFAD::val(lodetemp) >= 1.)
	{
		lodetemp *= T(0.999);
		// shapeFAD::val(lodetemp)*= 0.999;
		// lodetemp = T(1.);
		//DebugStop();
	}
	
	//DLODE = (-2*Dj3*j2 + 3*Dj2*J3t)/(2.*pow(J2t,2.5)*sqrt(1.3333333333333333 - (9*pow(J3t,2))/pow(J2t,3)))
	//1
	T j33= T(3.)*J3t;
	dJ2t.Multiply(j33,1);
	//2
	T j22=T(2.)*J2t;
	dJ3t.Multiply(j22,1);
	//3
	dJ2t.Add(dJ3t,-1);
	//4
	//if(shapeFAD::val(J2t)<1.e-6)J2t=1.e-6;
	T checknegativeroot = ((T(9.)*J3t*J3t)/(J2t*J2t*J2t));
	if(shapeFAD::val(checknegativeroot)>=4/3.)checknegativeroot *=T(0.999);
	T denom = T(2.)*sqrt(J2t*J2t*J2t*J2t*J2t)*sqrt((T(4/3.))-checknegativeroot);
	//T denom2 = (T(2.)*pow(J2t,2.5)*sqrt(T(1.3333333333333333) - (T(9.)*pow(J3t,2.))/pow(J2t,3)));
	T oneoverden = T(1.)/denom;
	dJ2t.Multiply(oneoverden,1);
	dJ2t*=oneoverden;
	//    GradLode = dJ2t;
	
	T acoslodetemp = acos(lodetemp);
	Lode = acoslodetemp / T(3.);
	
	if(shapeFAD::val(Lode) >= (M_PI/3.)-0.0001)
	{
		Lode *= T(0.999);
	}
#ifdef MACOS
	//    feclearexcept(FE_ALL_EXCEPT);
	//	int res = fetestexcept(FE_ALL_EXCEPT);
	//	if(res)
	//	{
	//		std::cout << " \n " << __PRETTY_FUNCTION__ <<"\n NAN DETECTED \n";
	//		DebugStop();
	//	}
	feclearexcept(FE_ALL_EXCEPT);
	
	if(fetestexcept(/*FE_DIVBYZERO*/ FE_ALL_EXCEPT	)) {
		std::cout << "division by zero reported\n";
		DebugStop();
	}
	if(fetestexcept(FE_INVALID)) {
		std::cout << "invalid result reported\n";
		DebugStop();
	}
#endif
	//	T denominador;
	//	denominador = ( J2t*J2t*J2t * T(4.) - J3t*J3t * T(27.) );
	//	T temp;
	//	if(shapeFAD::val(denominador) < 1.e-6)
	//	{
	//		dLodeAngledJ2 =  0.;
	//	}
	//	else
	//	{
	//		temp = sqrt( T(3.) / denominador );
	//		dLodeAngledJ2 =  temp * J3t / J2t * T(3./2.);
	//		dLodeAngledJ3 = -temp;
	//	}
	//
	//
	//	TempTensor = dJ2t;
	//	TempTensor *= dLodeAngledJ2;
	//	GradLode = dJ3t;
	//	GradLode *= dLodeAngledJ3;
	//	GradLode += TempTensor;
	
	
	//dLodeAngleJ2 = dLodeAngledJ2*dJ2 +
	//			 dLodeAngledJ3*dJ3;
	/*
	 TempTensor = dJ2t;
	 TempTensor *= dLodeAngledJ2;
	 dLodeAngle = dJ3t;
	 dLodeAngle *= dLodeAngledJ3;
	 dLodeAngle += TempTensor;
	 GradLode = dLodeAngle;
	 */
	/*
	 
	 T temp1= -( T(9.*sqrt(3.)) * J3t) / ( T(4.)*sqrt(J2t*J2t*J2t*J2t*J2t) );
	 TPZTensor<T> resultdj2(dJ2t);
	 resultdj2.Multiply(temp1,T(1.));
	 
	 T temp2 = ( T(3.)*sqrt( T( 3. ) ) )/ ( T( 2. ) * sqrt( J2t * J2t * J2t ) );
	 TPZTensor<T> resultdj3(dJ3t);
	 resultdj3.Multiply(temp2,T(1.));
	 
	 T temp33 = ( T( 27. ) * ( J3t * J3t ) )/( T( 4. ) *  J2t * J2t * J2t);
	 
	 if(shapeFAD::val(temp33) <= -1.)
	 {
	 std::cout << "ERRO NO CALCULO DO GRADLOD EM "<< __PRETTY_FUNCTION__<< std::endl;
	 DebugStop();
	 }
	 
	 T temp3 = ( T( 3. )  *  sqrt( T( 1. ) - temp33 ));
	 
	 TPZTensor<T>  RESULT(resultdj2);
	 
	 RESULT.Add(resultdj3, T(1.));
	 
	 RESULT.Multiply( T(1.) / temp3 ,T(-1.));
	 
	 GradLode = RESULT;
	 */
	
	/*
	 T lodetemp2 = ( T( 3. ) * sqrt( T(3.) ) * J3t ) / ( T( 2. ) * sqrt( J2t*J2t*J2t ) );
	 
	 if(shapeFAD::val(lodetemp2) < -1.00000001)
	 {
	 std::cout << "\n TPZTensor LodeTemp < -1\n";
	 std::cout << __PRETTY_FUNCTION__ << std::endl;
	 DebugStop();
	 }
	 
	 if(shapeFAD::val(lodetemp2) > 1.000000001)
	 {
	 std::cout << "\n TPZTensor LodeTemp > 1\n";
	 std::cout << __PRETTY_FUNCTION__ << std::endl;
	 DebugStop();
	 }
	 
	 T lodeangle =acos(lodetemp2);
	 
	 Lode = lodeangle;
	 
	 //T DLODE = -(   (3.*sqrt(3.) * dJ3 ) / ( 2.* pow(J2,1.5) ) - (9.*sqrt(3.)*dJ2*J3)  /  (4.* pow(J2,2.5) )    ) /(3.*sqrt(1. - (27.*pow(J3,2.))/(4.*pow(J2,3.))));
	 
	 T temporaryLode1 = ( 3. * sqrt( 3. ) ) / ( 2.* ( 2. * sqrt(J2t*J2t*J2t) ) );
	 dJ3t.Multiply(temporaryLode1, T(1.));
	 T temporaryLode2 = T(-9.* sqrt(3.)) * J3t;
	 T powj2 = 4.* (sqrt(J2t*J2t*J2t*J2t*J2t));//pow(J2,2.5)
	 temporaryLode2 = temporaryLode2/powj2;
	 dJ2t.Multiply(temporaryLode2, T(1.));
	 dJ3t.Add(dJ2t,T(-1));
	 T temporaryLode3 = ( 3. * sqrt(1. - (27.*(J3t*J3t)) / (4.*(J3t*J3t*J3t)) ) );
	 dJ3t.Multiply(T(1.)/temporaryLode3, T(1.));
	 dJ3t.Multiply(T(1.),T(-1.));jghjghjghj
	 GradLode = dJ3t;
	 */
}

template <class T>
void TPZTensor<T>::Eigenvalue(TPZTensor<T> &eigenval,TPZTensor<T> &dSigma1,TPZTensor<T> &dSigma2,TPZTensor<T> &dSigma3)const
{
	T I1(this->I1()), J2(this->J2());
	//	T J3(this->J3());
	
	if(fabs( shapeFAD::val(J2) ) < 1.e-6)J2 = 1.e-6;
	//	if(fabs( shapeFAD::val(J2) ) < 1.e-6)return;
	
	T sqrtJ2 = sqrt(J2);
	if(fabs( shapeFAD::val(sqrtJ2) ) < 1.e-6)sqrtJ2 = 1.e-6;
	
	TPZTensor<T> dJ2, dJ3,dLode;
	T Lode;
	
	this->dJ2(dJ2);
	this->dJ3(dJ3);
	this->Lodeangle(dLode,Lode);
	
	T pi23 = T(2. * M_PI / 3.);
	T TwoOverSqrThree = T(2./sqrt(3.));
	T TwoOverSqrThreeJ2 = TwoOverSqrThree * sqrtJ2;
	T I13 = I1 / T(3.);
	
	T tempCosLode = 		cos(Lode) 		 * TwoOverSqrThreeJ2;
	T tempCosMinusLode =	cos(Lode - pi23) * TwoOverSqrThreeJ2;
	T tempCosPlusLode =		cos(Lode + pi23) * TwoOverSqrThreeJ2;
	
	if(shapeFAD::val(Lode) < 0.)
	{
		cout << "Lode angle é Menor que ZERO. Valido somente para sig1 > sig2 > sig3 -> 0 < theta < Pi/3 " <<endl;
		DebugStop();
	}
	if(shapeFAD::val(Lode) > M_PI/3.)
	{
		cout << "Lode angle é Maior que Pi/3. Valido somente para sig1 > sig2 > sig3 -> 0 < theta < Pi/3 " <<endl;
		DebugStop();
	}
	
	/*ORIGINAL*/
	//Valido somente para sig1 > sig2 > sig3 -> 0 < theta < Pi/3
	
	eigenval.XX() = I13 + tempCosLode;
	eigenval.YY() = I13 + tempCosMinusLode;
	eigenval.ZZ() = I13 + tempCosPlusLode;
	eigenval.XY() *= T(0.);
	eigenval.XZ() *= T(0.);
	eigenval.YZ() *= T(0.);
	
	/*	//SOUZA NETO!
	 eigenval.XX() = I13 + tempCosPlusLode;
	 eigenval.YY() = I13 + tempCosLode;
	 eigenval.ZZ() = I13 + tempCosMinusLode;
	 eigenval.XY() *= T(0.);
	 eigenval.XZ() *= T(0.);
	 eigenval.YZ() *= T(0.);
	 */
	
	
	
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "\n  TPZTENSOR \n"<<endl;
		sout << "\n  LodeAngle = \n"<<Lode <<endl;
		sout << "\n  dLodeAngle= "<<dLode<<endl;
		sout << "\n\n";
		LOGPZ_INFO(loggerr,sout.str());
	}
#endif
	
	
	
	T OneOverTwoJ2 = T(0.5) / J2;
	TPZTensor<T> dI13;
	dI13.Identity();
	dI13 *= T(1./3.);
	
	tempCosLode      *= OneOverTwoJ2;
	tempCosMinusLode *= OneOverTwoJ2;
	tempCosPlusLode  *= OneOverTwoJ2;
	
	dSigma1 = dJ2;
	dSigma1 *= tempCosLode;
	dSigma1 += dI13;
	TPZTensor<T> dLodeAngleTemp(dLode);
	dLodeAngleTemp *= sin(Lode) * TwoOverSqrThreeJ2;
	dSigma1 -= dLodeAngleTemp;
	
	dSigma2 = dJ2;
	dSigma2 *= tempCosMinusLode;
	dSigma2 += dI13;
	dLodeAngleTemp = dLode;
	dLodeAngleTemp *= sin(pi23 - Lode) * TwoOverSqrThreeJ2;
	dSigma2 += dLodeAngleTemp;
	
	dSigma3 = dJ2;
	dSigma3 *= tempCosPlusLode;
	dSigma3 += dI13;
	dLodeAngleTemp = dLode;
	dLodeAngleTemp *= sin(pi23 + Lode) * TwoOverSqrThreeJ2;
	dSigma3 -= dLodeAngleTemp;
	
	
	
	/*
	 //T dengen1 =i/3. + (DJ2*Cos(lode))/(Sqrt(3)*Sqrt(J2)) - (2*gradLode*Sqrt(J2)*Sin(lode))/Sqrt(3)
	 
	 TPZTensor<T> I;
	 I.Identity();
	 TPZTensor<T> Itemp(I);
	 Itemp.Multiply(T(1./3.), T(1.));
	 T temporary2 = cos(Lode)/ (sqrt(T(3.))*sqrt(J2));
	 TPZTensor<T> dJ2temp(dJ2);
	 dJ2temp.Multiply(temporary2, T(1.));
	 T temporary3 = -2.*(sqrt(J2)*sin(Lode))/ (sqrt(3.));
	 TPZTensor<T> dLodetemp(dLode);
	 dLodetemp.Multiply(temporary3, T(1.));
	 Itemp.Add(dJ2temp, T(1));
	 Itemp.Add(dLodetemp, T(1));
	 dSigma1 = Itemp;
	 
	 //T eigen2 = i/3. - (2*gradLode*Sqrt(J2)*Cos(lode - Pi/6.))/Sqrt(3) + (DJ2*Sin(lode))/(Sqrt(3)*Sqrt(J2))
	 
	 TPZTensor<T> Itemp2(I);
	 Itemp2.Multiply(T(1./3.), T(1.));
	 T temporary4 = sin(Lode)/ (sqrt(T(3.))*sqrt(J2));
	 TPZTensor<T> dJ2temp2(dJ2);
	 dJ2temp2.Multiply(temporary4, T(1.));
	 T temporary5 = (-2 * sqrt(J2) * cos( Lode-T(M_PI/6.) ) )/ (sqrt(3.));
	 TPZTensor<T> dLodetemp2(dLode);
	 dLodetemp2.Multiply(temporary5, T(1.));
	 Itemp2.Add(dJ2temp2, T(1));
	 Itemp2.Add(dLodetemp2, T(1));
	 dSigma2 = Itemp2;
	 
	 //T eigen2 = i/3. - (2*gradLode*Sqrt(J2)*Cos(lode + Pi/6.))/Sqrt(3) + (DJ2*Sin(lode))/(Sqrt(3)*Sqrt(J2))
	 
	 TPZTensor<T> Itemp3(I);
	 Itemp3.Multiply(T(1./3.), T(1.));
	 T temporary6 = sin(Lode)/ (sqrt(T(3.))*sqrt(J2));
	 TPZTensor<T> dJ2temp3(dJ2);
	 dJ2temp3.Multiply(temporary6, T(1.));
	 T temporary7 = (-2 * sqrt(J2) * cos( Lode+T(M_PI/6.) ) )/ (sqrt(3.));
	 TPZTensor<T> dLodetemp3(dLode);
	 dLodetemp3.Multiply(temporary7, T(1.));
	 Itemp3.Add(dJ2temp3, T(1));
	 Itemp3.Add(dLodetemp3, T(1));
	 dSigma3 = Itemp3;
	 
	 */
}

template <>
inline void TPZTensor<REAL>::Print(std::ostream &output) const {
	output << "XX = " << XX() << " XY = " << XY() << " XZ = " << XZ() << " YY = " << YY() << " YZ() = " << YZ() << " ZZ = " << ZZ() << std::endl;
}

template <>
inline void TPZTensor<TFad<6,REAL> >::Print(std::ostream &output) const {
	output << "XX = " << XX() << "\nXY = " << XY() << "\nXZ = " << XZ() << "\nYY = " << YY() << "\nYZ() = " << YZ() << "\nZZ = " << ZZ() << std::endl;
}

template <>
inline void TPZTensor<TFad<9,REAL> >::Print(std::ostream &output) const {
	output << "XX = " << XX() << "\nXY = " << XY() << "\nXZ = " << XZ() << "\nYY = " << YY() << "\nYZ() = " << YZ() << "\nZZ = " << ZZ() << std::endl;
}


template <class T>
inline void TPZTensor<T>::Print(std::ostream &output) const {
	output << __PRETTY_FUNCTION__ << " please implement me\n";
}
//template class TPZTensor<REAL>;

#endif //TPZTENSOR_H
