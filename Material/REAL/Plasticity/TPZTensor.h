// $Id: TPZTensor.h,v 1.28 2010-12-04 20:41:28 diogo Exp $

#ifndef TPZTENSOR_H
#define TPZTENSOR_H

#include "pzmanvector.h"
#include "pzfmatrix.h"
#include "pzextractval.h"
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
class TPZTensor : public TPZSavable {
public:

    struct TPZDecomposed {
        unsigned int fDistinctEigenvalues;
        TPZManVector<unsigned int, 3> fGeometricMultiplicity;
        TPZManVector<T, 3> fEigenvalues;
        TPZManVector<TPZManVector<T, 3>, 3> fEigenvectors;
        TPZManVector<TPZTensor<T>, 3> fEigentensors; // Tensors of the spectral decomposition. If there is a repeated eigenvalue, only the first occurrence should be filled.

        TPZDecomposed() : fDistinctEigenvalues(0), fGeometricMultiplicity(3, 0), fEigenvalues(3, 0.), fEigenvectors(3), fEigentensors(3) {
        }

        TPZDecomposed(const TPZDecomposed &copy) : fDistinctEigenvalues(copy.fDistinctEigenvalues), fGeometricMultiplicity(copy.fGeometricMultiplicity), fEigenvalues(copy.fEigenvalues), fEigenvectors(copy.fEigenvectors), fEigentensors(copy.fEigentensors) {
        }

        TPZDecomposed &operator=(const TPZDecomposed &copy) {
            fDistinctEigenvalues = copy.fDistinctEigenvalues;
            fGeometricMultiplicity = copy.fGeometricMultiplicity;
            fEigenvalues = copy.fEigenvalues;
            fEigenvectors = copy.fEigenvectors;
            fEigentensors = copy.fEigentensors;
            return *this;
        }

        TPZDecomposed(const TPZTensor<T> &source) : fDistinctEigenvalues(0), fGeometricMultiplicity(3, 0), fEigenvalues(3, 0.), fEigenvectors(3), fEigentensors(3) {
            source.EigenSystem(*this);
        }

        void Print(std::ostream &out) const {
#ifdef PZDEBUG
            if (fDistinctEigenvalues == 0 || fDistinctEigenvalues > 3) {
                std::stringstream str;
                str << "TPZTensor::Decomposed::Print Invalid number of distinct eigenvalues: " << fDistinctEigenvalues << std::endl;
                LOGPZ_DEBUG(loggerr, str.str());
            }
#endif
            out << "TPZTensor::Decomposed Eigenvalues (" << fDistinctEigenvalues << " distinct):" << std::endl;
            unsigned int vIndex = 0;
            unsigned int lambda = 0;
            do {
                const unsigned int geoMult(fGeometricMultiplicity[lambda]);
#ifdef PZDEBUG
                if (geoMult == 0 || geoMult > 3) {
                    std::stringstream str;
                    str << "TPZTensor::Decomposed::Print Invalid geometric multiplicity (" << geoMult << ") for eigenvalue " << fEigenvalues[lambda] << std::endl;
                    LOGPZ_DEBUG(loggerr, str.str());
                }
#endif
                out << "Eigenvalue: " << fEigenvalues[lambda] << " (" << geoMult << " times)" << std::endl;
                out << "Eigenvectors:" << std::endl;
                for (unsigned int v = 0; v < geoMult; ++v) {
                    out << fEigenvectors[vIndex][0] << " ";
                    out << fEigenvectors[vIndex][1] << " ";
                    out << fEigenvectors[vIndex][2] << std::endl;
                    ++vIndex;
                }
                out << "Eigentensor:" << std::endl;
                fEigentensors[lambda].Print(out);
                out << std::endl;
                lambda += geoMult; // Go to next distinct eigenvalue
#ifdef PZDEBUG
                if (lambda > 3) {
                    std::stringstream str;
                    str << "TPZTensor::Decomposed::Print Greater total geometric multiplicity than expected!" << std::endl;
                    LOGPZ_ERROR(loggerr, str.str());
                }
#endif
            } while (lambda < 3);
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
                    for (int i = 0; i < 3; i++) {
                        gEigval[i] = objdec.fEigenvalues[i];
                    }
                    for (int j = 0; j < 6; j++) {
                        tangent(0, j) = objdec.fEigentensors[icase][j];
                        tangent(0, 6) = objdec.fEigentensors[icase].XY();
                        tangent(0, 7) = objdec.fEigentensors[icase].XZ();
                        tangent(0, 8) = objdec.fEigentensors[icase].YZ();
                    }
                    STATE sum = 0., sum2 = 0.;
                    for (int i = 0; i < 9; i++) {
                        sum += tangent(0, i) * tangent(0, i);
                    }
                    sum2 = tangent(0, _XX_) + tangent(0, _YY_) + tangent(0, _ZZ_);
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
            obj[_XY_] = obj[_XY_]*0.5 + 0.5 * state(6);
            obj[_XZ_] = (obj[_XZ_] + state(7))*0.5;
            obj[_YZ_] = (obj[_YZ_] + state(8))*0.5;
            TPZTensor<STATE>::TPZDecomposed objdec(obj);
            residual.Resize(1, 1);
            switch (icase) {
                case 0:
                case 1:
                case 2:
                {
                    residual(0, 0) = objdec.fEigenvalues[icase];
                    int numequal = 1;
                    for (int i = 0; i < 3; i++) {
                        if (i == icase) {
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
            REAL sig1, sig2, sig3;
            sig1 = this->fEigenvalues[0];
            sig2 = this->fEigenvalues[1];
            sig3 = this->fEigenvalues[2];
            J2 = (pow(sig1 + (-sig1 - sig2 - sig3) / 3., 2) + pow(sig2 + (-sig1 - sig2 - sig3) / 3., 2) + pow((-sig1 - sig2 - sig3) / 3. + sig3, 2)) / 2.;
        }

		void ComputeJ3(REAL & J3)
		{
            REAL sig1, sig2, sig3, s1, s2, s3, temp;
            sig1 = this->fEigenvalues[0];
            sig2 = this->fEigenvalues[1];
            sig3 = this->fEigenvalues[2];
            s1 = sig1 - (1. / 3.)*(sig1 + sig2 + sig3);
            s2 = sig2 - (1. / 3.)*(sig1 + sig2 + sig3);
            s3 = sig3 - (1. / 3.)*(sig1 + sig2 + sig3);
            temp = (1. / 3.)*(s1 * s1 * s1 + s2 * s2 * s2 + s3 * s3 * s3);
            J3 = temp;
        }

        void ComputeI1(REAL &I1) {
            const REAL sig1 = this->fEigenvalues[0];
            const REAL sig2 = this->fEigenvalues[1];
            const REAL sig3 = this->fEigenvalues[2];
            I1 = sig1 + sig2 + sig3;
        }

		void ApplyStrainComputeElasticStress(TPZTensor<REAL> &Stress,REAL &K , REAL & G)
		{
            REAL sig1, sig2, sig3, s1, s2, s3;
            sig1 = this->fEigenvalues[0];
            sig2 = this->fEigenvalues[1];
            sig3 = this->fEigenvalues[2];
            s1 = sig1 - (1. / 3.)*(sig1 + sig2 + sig3);
            s2 = sig2 - (1. / 3.)*(sig1 + sig2 + sig3);
            s3 = sig3 - (1. / 3.)*(sig1 + sig2 + sig3);

            TPZTensor<REAL> SigVol;
            SigVol.XX() = K * (sig1 + sig2 + sig3);
            SigVol.YY() = K * (sig1 + sig2 + sig3);
            SigVol.ZZ() = K * (sig1 + sig2 + sig3);

            Stress.XX() = s1;
            Stress.YY() = s2;
            Stress.ZZ() = s3;
            Stress *= (2 * G);
            Stress += SigVol;
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
        for (unsigned int i = 0; i < eigensystem.fDistinctEigenvalues; i+=eigensystem.fGeometricMultiplicity[i]) {
#ifdef PZDEBUG
            if (i > 2){
                DebugStop();
            }
#endif
            Add(eigensystem.fEigentensors[i], eigensystem.fEigenvalues[i]);
        }
    }
        
        virtual int ClassId() const {
            return Hash("TPZTensor") ^ ClassIdOrHash<T>() << 1;
        }


    /// Method to write to a pzstream
        void Write(TPZStream& buf, int withclassid) const;

    ///Method to read the object from a pzstream
        void Read(TPZStream& buf, void* context);

	operator TPZFMatrix<T>() const
	{
        TPZFMatrix<T> result(3, 3);
        result(0, 0) = XX();
        result(0, 1) = XY();
        result(0, 2) = XZ();
        result(1, 0) = XY();
        result(1, 1) = YY();
        result(1, 2) = YZ();
        result(2, 0) = XZ();
        result(2, 1) = YZ();
        result(2, 2) = ZZ();
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
#ifdef PZDEBUG
        if (input.Rows() != 3 || input.Cols() != 3) {
            DebugStop();
        }
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                if (fabs(input.GetVal(i, j) - input.GetVal(j, i)) > 1.e-15) {
                    std::cout << "diff = " << input.GetVal(i, j) - input.GetVal(j, i) << std::endl;
                    DebugStop();
                }
            }
        }
#endif
        fData[_XX_] = input.GetVal(0, 0);
        fData[_XY_] = input.GetVal(0, 1);
        fData[_XZ_] = input.GetVal(0, 2);
        fData[_YY_] = input.GetVal(1, 1);
        fData[_YZ_] = input.GetVal(1, 2);
        fData[_ZZ_] = input.GetVal(2, 2);
    }
    /**
     * Method to print the tensor
     */
    void Print(std::ostream &out) const;

    /**
     Operator=
     */
    const TPZTensor<T> & operator=(const TPZTensor<T> &source);

    /**
     Operator+=
     */
    const TPZTensor<T> & operator+=(const TPZTensor<T> &source);

    /**
     Operator-=
     */
    const TPZTensor<T> & operator-=(const TPZTensor<T> &source);

    /**
     Operator*=
     */
    const TPZTensor<T> & operator*=(const T &multipl);

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
    void Add(const TPZTensor< T1 > & tensor, const T2 & constant);

    /**
     multiplica um scalar com o tensor atual
     @param [in] multipl valor sendo multiplicada
     @param [in] constant fator multiplicativo
     */
    template < class T1, class T2>
    void Multiply(const T1 & multipl, const T2 & constant);

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
    void Scale(const T2 & constant);

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
     * Metodo que calcula o negativo do segundo invariante da parte deviatorica do tensor
     * TBASE is needed when 3rd derivatives are of interest. In such cases the FAD
     * promotion fails.
     */
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
        Diag *= sigIJ[0] / T(3.);
        T J2 = this->J2();
        T sqj2 = sqrt(J2);
        REAL sqj2val = TPZExtractVal::val(sqj2);
        T ratio;
		if (fabs(sqj2val) > 1.e-6)
		{
            ratio = sigIJ[1] / sqj2;
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
    void Eigenvalue(TPZTensor<T> &eigenval, TPZTensor<T> &dSigma1, TPZTensor<T> &dSigma2, TPZTensor<T> &dSigma3)const;

    /**
     * Computes the Lode angle and its derivatives
     */
    void Lodeangle(TPZTensor<T> &GradLode, T &Lode)const;


    /**
     * Computes the eigenvectors and eigenvalues based on (page 742 Computational methods for plasticity/Souza Neto)
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
    void CopyFrom(const TPZFMatrix<T> & source);

    /** Converts the stress vector onto a symmetric stress tensor
     * @param Tensor [out]
     */
    void CopyToTensor(TPZFMatrix<T> & Tensor) const;

    /**
     Initializa o valor do tensor (tensor de deformacao)
     */
    void SetUp(const TPZVec<REAL> & Solution);

    /**
     * Updates the eigen decomposition object, filling the eigenvectors.
     * 
     * Computes the eigen projections if it wasn't done yet.
     * 
     * @param decomposed 
     */
    void ComputeEigenVectors(TPZDecomposed &eigensystem) const;

public:
    /**
     Dados do tensor
     */
    TPZManVector<T, 6> fData;

protected:

    static inline bool IsZeroVal(const T & val, REAL tol = 1e-9) {
        return (fabs(TPZExtractVal::val(val)) < tol);
    }

    bool AreEqual(const T &val1, const T &val2, const REAL tol = REAL(1e-9)) const {
        return (std::fabs(TPZExtractVal::val(val1) - TPZExtractVal::val(val2)) < tol);
    }

    /**
     * Computes the eigenprojections and eigenvalues based on (page 742 Computational methods for plasticity/Souza Neto)
     */
    void EigenProjection(const TPZVec<T> &EigenVals, int index, const TPZVec<int> &DistinctEigenvalues, TPZTensor<T> &Ei) const;
};

template <class T>
void TPZTensor<T>::Read(TPZStream& buf, void* context) {
    buf.Read(fData);
}

template <class T>
void TPZTensor<T>::Write(TPZStream& buf, int withclassid) const {
    buf.Write(fData);
}

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
    for (i = 0; i < 6; i++)fData[i] += source.fData[i];
    return *this;
}

template < class T >
const TPZTensor<T> & TPZTensor<T>::operator-=( const TPZTensor<T> &source )
{
    int i;
    for (i = 0; i < 6; i++)fData[i] -= source.fData[i];
    return *this;
}

template < class T >
const TPZTensor<T> & TPZTensor<T>::operator*=( const T &multipl )
{
    int i;
    for (i = 0; i < 6; i++)fData[i] *= multipl;
    return *this;
}

template < class T >
template < class T1, class T2 >
void TPZTensor<T>::Add(const TPZTensor< T1 > & tensor, const T2 & constant )
{
    int i, size = 6;
    for (i = 0; i < size; i++) {
        fData[i] += tensor.fData[i] * T1(constant);
    }
}

template < class T >
template < class T1, class T2 >
void TPZTensor<T>::Multiply(const T1 & multipl, const T2 & constant )
{
    int i, size = 6;
    for (i = 0; i < size; i++) {
        fData[i] *= multipl;
        fData[i] *= constant;
    }
}

template < class T >
void TPZTensor<T>::Multiply2(TPZTensor<T> tensor, TPZTensor<T> & resp)const
{
    T XX, YY, ZZ, XY, YZ, XZ;
    XX = this->fData[_XX_];
    YY = this->fData[_YY_];
    ZZ = this->fData[_ZZ_];
    XY = this->fData[_XY_];
    XZ = this->fData[_XZ_];
    YZ = this->fData[_YZ_];
    resp.fData[_XX_] = XX * tensor.XX() + XY * tensor.XY() + XZ * tensor.XZ();
    resp.fData[_YY_] = XY * tensor.XY() + YY * tensor.YY() + YZ * tensor.YZ();
    resp.fData[_ZZ_] = XZ * tensor.XZ() + YZ * tensor.YZ() + ZZ * tensor.ZZ();
    resp.fData[_XY_] = XX * tensor.XY() + XY * tensor.YY() + XZ * tensor.YZ();
    resp.fData[_XZ_] = XX * tensor.XZ() + XY * tensor.YZ() + XZ * tensor.ZZ();
    resp.fData[_YZ_] = XY * tensor.XZ() + YY * tensor.YZ() + YZ * tensor.ZZ();

}

template < class T >
template < class T2 >
void TPZTensor<T>::Scale(const T2 & constant )
{
    int i, size = 6;
    for (i = 0; i < size; i++) {
        fData[i] *= constant;
    }
}

template < class T >
void TPZTensor<T>::Zero()
{
    int i, size = 6;
    for (i = 0; i < size; i++) {
        fData[i] = 0.;
    }
}

template < class T >
void TPZTensor<T>::SetUp(const TPZVec<REAL> & Solution) {
    int i, size = 6;
    for (i = 0; i < size; i++)
        fData[i] = Solution[i];

}

template < class T >
template < class T1 >
void TPZTensor<T>::CopyTo(TPZTensor<T1> & target) const {
    int i, size = 6;
    for (i = 0; i < size; i++) {
        target[i] = TPZExtractVal::val(fData[i]);
    }

}

template < class T >
void TPZTensor<T>::CopyTo(TPZFMatrix<REAL> & target) const
{
    int i;
    for (i = 0; i < 6; i++)target(i, 0) = TPZExtractVal::val(fData[i]);
}

template < class T >
void TPZTensor<T>::CopyFrom(const TPZFMatrix<T> & source) {
    int i;
    for (i = 0; i < 6; i++)fData[i] = source.Get(i, 0);
}

template < class T >
T TPZTensor<T>::J2() const {
    TPZVec<T> s(3);
    DeviatoricDiagonal_T<T>(s);
    return - s[0] * s[1] - s[0] * s[2] - s[1] * s[2] +fData[_XY_] * fData[_XY_] + fData[_XZ_] * fData[_XZ_] + fData[_YZ_] * fData[_YZ_];
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
            fData[_XX_] * fData[_ZZ_]);
}

template < class T >
T TPZTensor<T>::I3() const {
    return fData[_XX_] * fData[_YY_] * fData[_ZZ_]
            +(fData[_XY_] * fData[_XZ_] * fData[_YZ_]) * 2.
            - (fData[_XX_] * fData[_YZ_] * fData[_YZ_] +
            fData[_YY_] * fData[_XZ_] * fData[_XZ_] +
            fData[_ZZ_] * fData[_XY_] * fData[_XY_]);
}

template < class T >
void TPZTensor<T>::dJ2(TPZTensor<T> & Tangent) const {
    T p = I1() / T(3.);
    Tangent.fData[_XX_] = fData[_XX_] - p;
    Tangent.fData[_YY_] = fData[_YY_] - p;
    Tangent.fData[_ZZ_] = fData[_ZZ_] - p;
    Tangent.fData[_XY_] = fData[_XY_] * T(2.);
    Tangent.fData[_XZ_] = fData[_XZ_] * T(2.);
    Tangent.fData[_YZ_] = fData[_YZ_] * T(2.);

}

template < class T >
template < class TBASE >
void TPZTensor<T>::DeviatoricDiagonal_T(TPZVec<T> & vec) const {

    T p = I1() / T(TBASE(3.));
    vec[0] = fData[_XX_] - p;
    vec[1] = fData[_YY_] - p;
    vec[2] = fData[_ZZ_] - p;
}

template < class T >
void TPZTensor<T>::DeviatoricDiagonal(TPZVec<T> & vec) const {

    DeviatoricDiagonal_T<REAL>(vec);
}

template < class T>
T TPZTensor<T>::Norm() const
{
    int i;
    T norm = T(0.);
    // norm	= sqrt(fData[_XX_]*fData[_XX_]+fData[_YY_]*fData[_YY_]+fData[_ZZ_]*fData[_ZZ_]+T(2.)*(fData[_XY_]*fData[_XY_])+T(2.)*(fData[_XZ_]*fData[_XZ_])+T(2.)*(fData[_YZ_]*fData[_YZ_]));
    for (i = 0; i < 6; i++) norm += fData[i] * fData[i];
    return sqrt(norm);
    //return norm;
}

template <class T>
std::ostream &operator<<(std::ostream &out,const TPZTensor<T> &tens)
{
    for (int i = 0; i < 6; i++) {
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
    return fData[_XX_] * fData[_YY_] * fData[_ZZ_] + fData[_XY_] * fData[_XZ_] * fData[_YZ_]*2. - fData[_XZ_] * fData[_YY_] * fData[_XZ_] -
            fData[_XY_] * fData[_XY_] * fData[_ZZ_] - fData[_YZ_] * fData[_YZ_] * fData[_XX_];
}

/**
 * Metodo que calcula o determinante do tensor
 */
template <class T>
void TPZTensor<T>::dDet(TPZTensor<T> &grad) const
{
    grad.fData[_XX_] = fData[_YY_] * fData[_ZZ_] - fData[_YZ_] * fData[_YZ_];
    grad.fData[_XY_] = fData[_XZ_] * fData[_YZ_] - fData[_XY_] * fData[_ZZ_];
    grad.fData[_XZ_] = fData[_XY_] * fData[_YZ_] - fData[_XZ_] * fData[_YY_];
    grad.fData[_YY_] = fData[_XX_] * fData[_ZZ_] - fData[_XZ_] * fData[_XZ_];
    grad.fData[_YZ_] = fData[_XY_] * fData[_XZ_] - fData[_YZ_] * fData[_XX_];
    grad.fData[_ZZ_] = fData[_XX_] * fData[_YY_] - fData[_XY_] * fData[_XY_];
}

template <class T>
void TPZTensor<T>::S(TPZTensor<T> &s) const
{
    s = *this;
    TPZTensor<T> one;
    one.Identity();
    T mult = -I1() / 3.;
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
    T state0(fData[_XX_]), state1(fData[_XY_]), statex1(fData[_XY_]),
            state2(fData[_XZ_]), statex2(fData[_XZ_]), state3(fData[_YY_]),
            state4(fData[_YZ_]), statex4(fData[_YZ_]), state5(fData[_ZZ_]);
    deriv.fData[_XX_] = (fData[_XX_] * fData[_XX_]*2.) / 9. + fData[_XY_] * fData[_XY_] / 3. +
            fData[_XZ_] * fData[_XZ_] / 3. - (fData[_XX_] * fData[_YY_]*2.) / 9. -
            fData[_YY_] * fData[_YY_] / 9. - (fData[_YZ_] * fData[_YZ_]*2.) / 3. - (fData[_XX_] * fData[_ZZ_]*2.) /
            9. + (fData[_YY_] * fData[_ZZ_]*4.) / 9. - fData[_ZZ_] * fData[_ZZ_] / 9.;
    deriv.fData[_XY_] = (state0 * state1) / 3. + (state1 * state3) / 3. -
            (state1 * state5 * 2.) / 3. + (state0 * statex1) / 3. +
            (state3 * statex1) / 3. - (state5 * statex1 * 2.) / 3. +
            state4 * statex2 + state2*statex4;
    deriv.fData[_XZ_] = (state0 * state2) / 3. - (state2 * state3 * 2.) / 3. +
            state1 * state4 + (state2 * state5) / 3. +
            (state0 * statex2) / 3. - (state3 * statex2 * 2.) / 3. +
            (state5 * statex2) / 3. + statex1*statex4;
    deriv.fData[_YY_] = -(fData[_XX_] * fData[_XX_]) / 9. +
            (fData[_XY_] * fData[_XY_]) / 3. - ((fData[_XZ_] * fData[_XZ_])*2.) / 3. - (fData[_XX_] * fData[_YY_]*2.) / 9. + ((fData[_YY_] * fData[_YY_])*2.) / 9. +
            (fData[_YZ_] * fData[_YZ_]) / 3. + (fData[_XX_] * fData[_ZZ_]*4.) / 9. - (fData[_YY_] * fData[_ZZ_]*2.) / 9. -
            (fData[_ZZ_] * fData[_ZZ_]) / 9.;
    deriv.fData[_YZ_] = (state0 * state4 * (-2.)) / 3. + (state3 * state4) / 3. +
            (state4 * state5) / 3. + state2 * statex1 +
            state1 * statex2 - (state0 * statex4 * 2.) / 3. +
            (state3 * statex4) / 3. + (state5 * statex4) / 3.;
    deriv.fData[_ZZ_] = -(fData[_XX_] * fData[_XX_]) / 9. - (2. * (fData[_XY_] * fData[_XY_])) / 3. +
            (fData[_XZ_] * fData[_XZ_]) / 3. + (fData[_XX_] * fData[_YY_]*4.) / 9. - (fData[_YY_] * fData[_YY_]) / 9. +
            (fData[_YZ_] * fData[_YZ_]) / 3. - (fData[_XX_] * fData[_ZZ_]*2.) / 9. - (fData[_YY_] * fData[_ZZ_]*2.) /
            9. + ((fData[_ZZ_] * fData[_ZZ_])*2.) / 9.;

}

template <class T>
void TPZTensor<T>::CopyToTensor(TPZFMatrix<T> & Tensor) const {
    Tensor.Resize(3, 3);
    Tensor(0, 0) = TPZExtractVal::val(XX()); //xx
    Tensor(0, 1) = TPZExtractVal::val(XY()); //xy
    Tensor(0, 2) = TPZExtractVal::val(XZ()); //xz

    Tensor(1, 0) = TPZExtractVal::val(XY()); //xy
    Tensor(1, 1) = TPZExtractVal::val(YY()); //yy
    Tensor(1, 2) = TPZExtractVal::val(YZ()); //yz

    Tensor(2, 0) = TPZExtractVal::val(XZ()); //xz
    Tensor(2, 1) = TPZExtractVal::val(YZ()); //yz
    Tensor(2, 2) = TPZExtractVal::val(ZZ()); //zz
}

template<class T>
static T Polynom(const T &x, const T &I1, const T &I2, const T &I3)
{
    T result = x * x * x - I1 * x * x + I2 * x - I3;
    return result;
}

template<class T>
static T DericPolynom(const T &x, const T &I1, const T &I2, const T &I3)
{
    T result = T(3.) * x * x - T(2.) * I1 * x + I2;
    return result;
}

template<class T>
static T UpdateNewton(const T &x, const T &I1, const T &I2, const T &I3)
{
    T residue = Polynom(x, I1, I2, I3);
    T dres = DericPolynom(x, I1, I2, I3);
    T result = x - residue / dres;
    return result;
}

template <class T>
void TPZTensor<T>::EigenSystemJacobi(TPZDecomposed &eigensystem)const
{
    TPZFNMatrix<9, T> TensorMat(*this);

    long numiterations = 1000;
    REAL tol = 1.e-6;

    eigensystem.fEigenvectors[0] = TPZManVector<T, 3>(3, 0.);
    eigensystem.fEigenvectors[1] = TPZManVector<T, 3>(3, 0.);
    eigensystem.fEigenvectors[2] = TPZManVector<T, 3>(3, 0.);

    TPZManVector<T, 3> Eigenvalues(3, 0.);
    TPZFNMatrix<9, T> Eigenvectors(3, 3, 0.);
    TensorMat.SolveEigensystemJacobi(numiterations, tol, Eigenvalues, Eigenvectors);


    TPZManVector<TPZTensor<T>, 3> &Eigentensors = eigensystem.fEigentensors;
    eigensystem.fEigenvalues = Eigenvalues;

    //tres autovalores iguais
    if (AreEqual(Eigenvalues[0], Eigenvalues[1]) && AreEqual(Eigenvalues[0], Eigenvalues[2])) {
        eigensystem.fDistinctEigenvalues = 1;
        for (unsigned int i = 0; i < 3; ++i) {
            eigensystem.fGeometricMultiplicity[i] = 3;
            eigensystem.fEigenvectors[i][i] = 1.;
            Eigentensors[i].Identity();
        }
    } else {
        for (unsigned int i = 0; i < 3; ++i) {
            for (unsigned int j = 0; j < 3; ++j) {
                eigensystem.fEigenvectors[i][j] = Eigenvectors(i, j);
            }
        }

        if (AreEqual(Eigenvalues[0], Eigenvalues[1]) || AreEqual(Eigenvalues[1], Eigenvalues[2]) || AreEqual(Eigenvalues[0], Eigenvalues[2])) { //dois autovalores iguais
            eigensystem.fDistinctEigenvalues = 2;

            int different = -1;
            int equals[2] = {-1, -1};
            if (AreEqual(Eigenvalues[0], Eigenvalues[1])) {
                different = 2;
                equals[0] = 0;
                equals[1] = 1;
            } else if (AreEqual(Eigenvalues[0], Eigenvalues[2])) {
                different = 1;
                equals[0] = 0;
                equals[1] = 2;
            } else if (AreEqual(Eigenvalues[1], Eigenvalues[2])) {
                different = 0;
                equals[0] = 1;
                equals[1] = 2;
            }
            eigensystem.fGeometricMultiplicity[different] = 1;
            eigensystem.fGeometricMultiplicity[equals[0]] = 2;
            eigensystem.fGeometricMultiplicity[equals[1]] = 2;

            Eigentensors[different].XX() = Eigenvectors(different, 0) * Eigenvectors(different, 0);
            Eigentensors[different].XY() = Eigenvectors(different, 0) * Eigenvectors(different, 1);
            Eigentensors[different].XZ() = Eigenvectors(different, 0) * Eigenvectors(different, 2);
            Eigentensors[different].YY() = Eigenvectors(different, 1) * Eigenvectors(different, 1);
            Eigentensors[different].YZ() = Eigenvectors(different, 1) * Eigenvectors(different, 2);
            Eigentensors[different].ZZ() = Eigenvectors(different, 2) * Eigenvectors(different, 2);

            Eigentensors[equals[0]].Identity();
            Eigentensors[equals[0]] -= Eigentensors[different];
            Eigentensors[equals[1]] = Eigentensors[equals[0]];
        } else {
            //nenhum autovalor igual
            eigensystem.fDistinctEigenvalues = 3;
            for (int i = 0; i < 3; i++) {
                eigensystem.fGeometricMultiplicity[i] = 1;
                Eigentensors[i].XX() = Eigenvectors(i, 0) * Eigenvectors(i, 0);
                Eigentensors[i].XY() = Eigenvectors(i, 0) * Eigenvectors(i, 1);
                Eigentensors[i].XZ() = Eigenvectors(i, 0) * Eigenvectors(i, 2);
                Eigentensors[i].YY() = Eigenvectors(i, 1) * Eigenvectors(i, 1);
                Eigentensors[i].YZ() = Eigenvectors(i, 1) * Eigenvectors(i, 2);
                Eigentensors[i].ZZ() = Eigenvectors(i, 2) * Eigenvectors(i, 2);
            }
        }//3 autovalores diferentes
    }

    //    for (unsigned int i = 0; i < 3; ++i) {
    //        Eigentensors[i].XX() = EigenvecJac(0, i) * EigenvecJac(0, i);
    //        Eigentensors[i].XY() = EigenvecJac(0, i) * EigenvecJac(1, i);
    //        Eigentensors[i].XZ() = EigenvecJac(0, i) * EigenvecJac(2, i);
    //        Eigentensors[i].YY() = EigenvecJac(1, i) * EigenvecJac(1, i);
    //        Eigentensors[i].YZ() = EigenvecJac(1, i) * EigenvecJac(2, i);
    //        Eigentensors[i].ZZ() = EigenvecJac(2, i) * EigenvecJac(2, i);
    //    }
#ifdef PZDEBUG
#ifdef LOG4CXX
    if (loggerr->isDebugEnabled()) {
        std::stringstream str;
        str << "\n-------------AUTOVETORES JACOBI--------------" << std::endl;
        str << "Tensor:" << std::endl;
        this->Print(str);
        str << "Eigenvalues = " << Eigenvalues << std::endl;
        str << "Eigenvectors:" << std::endl;
        Eigenvectors.Print("Eigenvec:", str);
        str << "Eigenprojections:" << std::endl;
        for (int k = 0; k < 3; k++) {
            str << "EigenProjection " << k << ":" << std::endl;
            Eigentensors[k].Print(str);
        }

        str << "\n-------------FIM AUTOVETORES JACOBI--------------" << std::endl;
        LOGPZ_DEBUG(loggerr, str.str())
    }
#endif
#endif
}

template <class T>
void TPZTensor<T>::EigenSystem(TPZDecomposed &eigensystem)const {

    TPZManVector<T, 3> &Eigenvalues = eigensystem.fEigenvalues;
    TPZManVector<TPZTensor<T>, 3> &Eigentensors = eigensystem.fEigentensors;

    TPZTensor<T> eigenValuesTensor;
    this->EigenValue(eigenValuesTensor);
    Eigenvalues[0] = eigenValuesTensor.XX();
    Eigenvalues[1] = eigenValuesTensor.YY();
    Eigenvalues[2] = eigenValuesTensor.ZZ();

    //tres autovalores iguais
    if (AreEqual(Eigenvalues[0], Eigenvalues[1]) && AreEqual(Eigenvalues[1], Eigenvalues[2])) {
        eigensystem.fDistinctEigenvalues = 1;
        eigensystem.fGeometricMultiplicity[0] = 3;
        eigensystem.fGeometricMultiplicity[1] = 3;
        eigensystem.fGeometricMultiplicity[2] = 3;
        Eigentensors[0].Identity();
        Eigentensors[1].Identity();
        Eigentensors[2].Identity();
    } else if (AreEqual(Eigenvalues[0], Eigenvalues[1]) || AreEqual(Eigenvalues[1], Eigenvalues[2]) || AreEqual(Eigenvalues[0], Eigenvalues[2])) {
        //dois autovalores iguais
        eigensystem.fDistinctEigenvalues = 2;

        int different = -1;
        int equals[2] = {-1, -1};
        if (AreEqual(Eigenvalues[0], Eigenvalues[1])) {
            different = 2;
            equals[0] = 0;
            equals[1] = 1;
        } else if (AreEqual(Eigenvalues[0], Eigenvalues[2])) {
            different = 1;
            equals[0] = 0;
            equals[1] = 2;
        } else if (AreEqual(Eigenvalues[1], Eigenvalues[2])) {
            different = 0;
            equals[0] = 1;
            equals[1] = 2;
        }
        eigensystem.fGeometricMultiplicity[different] = 1;
        eigensystem.fGeometricMultiplicity[equals[0]] = 2;
        eigensystem.fGeometricMultiplicity[equals[1]] = 2;

        TPZManVector<int, 2> DistinctEigenvalues(2);
        DistinctEigenvalues[0] = different;
        DistinctEigenvalues[1] = equals[0];

        this->EigenProjection(Eigenvalues, different, DistinctEigenvalues, Eigentensors[different]);
        Eigentensors[equals[0]].Identity();
        Eigentensors[equals[0]] -= Eigentensors[different];
        Eigentensors[equals[1]] = Eigentensors[equals[0]];
    } else {
        //3 autovalores diferentes
        eigensystem.fDistinctEigenvalues = 3;
        TPZManVector<int, 3> DistinctEigenvalues(3);
        DistinctEigenvalues[0] = 0;
        DistinctEigenvalues[1] = 1;
        DistinctEigenvalues[2] = 2;
        for (int i = 0; i < 3; i++) {
            eigensystem.fGeometricMultiplicity[i] = 1;
            this->EigenProjection(Eigenvalues, i, DistinctEigenvalues, Eigentensors[i]);
        }
    }
#ifdef PZDEBUG
    TPZTensor<T> total;
    unsigned int geometricCount = 0;
    for (unsigned int i = 0; i < 3 && geometricCount < 3; ++i) {
        total.Add(Eigentensors[i], Eigenvalues[i]);
        geometricCount += eigensystem.fGeometricMultiplicity[i];
    }
    for (unsigned int i = 0; i < 6; ++i) {
        if (!AreEqual(total[i], this->operator [](i))){
            std::cout << "Tensor decomposition error: " << std::endl;
            std::cout << "Original Tensor: " << std::endl;
            this->Print(std::cout);
            std::cout << "Decomposition: " << std::endl;
            eigensystem.Print(std::cout);
            std::cout << "Reconstruction from decomposition: " << std::endl;
            total.Print(std::cout);
            DebugStop();
        }
    }
#endif
}

template <class T>
void TPZTensor<T>::ComputeEigenVectors(TPZDecomposed &eigensystem) const {
    // Eigen decomposition not computed yet. Let's do it.
    if (eigensystem.fDistinctEigenvalues == 0) {
        this->EigenSystem(eigensystem);
    }

    eigensystem.fEigenvectors[0] = TPZManVector<T, 3>(3, 0.);
    eigensystem.fEigenvectors[1] = TPZManVector<T, 3>(3, 0.);
    eigensystem.fEigenvectors[2] = TPZManVector<T, 3>(3, 0.);
    switch (eigensystem.fDistinctEigenvalues) {
        case 1:
            eigensystem.fEigenvectors[0][0] = 1.;
            eigensystem.fEigenvectors[1][1] = 1.;
            eigensystem.fEigenvectors[2][2] = 1.;
            break;
        case 2:
        {
            // Either the two biggest or smaller eigenvalues is repeated.
            // In both cases the eigenvector in position 1 will be of a repeated eigenvalue
            // (Remember that the eigenvalues are sorted!)
#ifdef PZDEBUG
            if (eigensystem.fGeometricMultiplicity[1] != 2) {
                DebugStop();
            }
#endif
            unsigned int repeatedIndex, uniqueIndex;
            if (eigensystem.fGeometricMultiplicity[2] == 1) {
                // first and second eigenvalues are the same
#ifdef PZDEBUG
                if (eigensystem.fGeometricMultiplicity[0] != 2) {
                    DebugStop();
                }
                if (!AreEqual(eigensystem.fEigenvalues[0], eigensystem.fEigenvalues[1])) {
                    DebugStop();
                }
#endif  
                repeatedIndex = 0;
                uniqueIndex = 2;
            } else {
                // second and third eigenvalues are the same
#ifdef PZDEBUG
                if (eigensystem.fGeometricMultiplicity[2] != 2) {
                    DebugStop();
                }
                if (!AreEqual(eigensystem.fEigenvalues[1], eigensystem.fEigenvalues[2])) {
                    DebugStop();
                }
#endif  
                repeatedIndex = 2;
                uniqueIndex = 0;
            }
            TPZFMatrix<T> tempMatrix(3, 3, 0);
            tempMatrix(0, 0) = eigensystem.fEigentensors[uniqueIndex].XX();
            tempMatrix(0, 1) = eigensystem.fEigentensors[uniqueIndex].XY();
            tempMatrix(0, 2) = eigensystem.fEigentensors[uniqueIndex].XZ();
            tempMatrix(1, 0) = eigensystem.fEigentensors[uniqueIndex].XY();
            tempMatrix(1, 1) = eigensystem.fEigentensors[uniqueIndex].YY();
            tempMatrix(1, 2) = eigensystem.fEigentensors[uniqueIndex].YZ();
            tempMatrix(2, 0) = eigensystem.fEigentensors[uniqueIndex].XZ();
            tempMatrix(2, 1) = eigensystem.fEigentensors[uniqueIndex].YZ();
            tempMatrix(2, 2) = eigensystem.fEigentensors[uniqueIndex].ZZ();
            for (unsigned int i = 0; i < 3; ++i) {
                int nonZeroIndex = -1;
                for (unsigned int j = 0; j < 3; ++j) {
                    if (!IsZeroVal(tempMatrix(i, j))) {
                        nonZeroIndex = j;
                        break;
                    }
                }
                if (nonZeroIndex != -1) {
                    eigensystem.fEigenvectors[uniqueIndex][0] = tempMatrix(i, 0);
                    eigensystem.fEigenvectors[uniqueIndex][1] = tempMatrix(i, 1);
                    eigensystem.fEigenvectors[uniqueIndex][2] = tempMatrix(i, 2);
                }
                break;
            }

            // If we take the eigentensor associated with the repeated eigenvalue
            // T = v1_i v1_j + v2_i v2_j
            // We can compute the two eigenvectors writing their first components w.r.t.
            // the other two components.
            // Let l be the repeated eigenvalue. Then we know that
            // Av=lv
            // for an eigenvector v. Then
            // (A-lI)v = 0
            // and the matrix (A-lI) has rank 1 (note that l is the (only) repeated eigenvalue).
            // Therefore, all lines will be linearly dependent. In order to solve for v,
            // we only need to find a suitable line. 
            tempMatrix.Identity();
            tempMatrix *= (-eigensystem.fEigenvalues[repeatedIndex]);
            tempMatrix(0, 0) += this->XX();
            tempMatrix(0, 1) += this->XY();
            tempMatrix(0, 2) += this->XZ();
            tempMatrix(1, 0) += this->XY();
            tempMatrix(1, 1) += this->YY();
            tempMatrix(1, 2) += this->YZ();
            tempMatrix(2, 0) += this->XZ();
            tempMatrix(2, 1) += this->YZ();
            tempMatrix(2, 2) += this->ZZ();
            // First we check if the line has non-zero elements:
            for (unsigned int i = 0; i < 3; ++i) {
                int nonZeroIndex = -1;
                for (unsigned int j = 0; j < 3; ++j) {
                    if (!IsZeroVal(tempMatrix(i, j))) {
                        nonZeroIndex = j;
                        break;
                    }
                }
                if (nonZeroIndex != -1) {
                    unsigned int otherIndex1 = (nonZeroIndex + 1) % 3;
                    unsigned int otherIndex2 = (nonZeroIndex + 2) % 3;

                    eigensystem.fEigenvectors[1][otherIndex1] = 1.0;
                    eigensystem.fEigenvectors[1][otherIndex2] = 0.0;
                    eigensystem.fEigenvectors[1][nonZeroIndex] = -tempMatrix(i, otherIndex1) / tempMatrix(i, nonZeroIndex);

                    eigensystem.fEigenvectors[repeatedIndex][otherIndex1] = 0.0;
                    eigensystem.fEigenvectors[repeatedIndex][otherIndex2] = 1.0;
                    eigensystem.fEigenvectors[repeatedIndex][nonZeroIndex] = -tempMatrix(i, otherIndex2) / tempMatrix(i, nonZeroIndex);
                    // all done!
                    break;
                }
            }

        }
            break;
        case 3:
        {
            TPZFMatrix<T> tempMatrix(3, 3, 0);
            for (unsigned int lambda = 0; lambda < 3; ++lambda) {
                tempMatrix(0, 0) = eigensystem.fEigentensors[lambda].XX();
                tempMatrix(0, 1) = eigensystem.fEigentensors[lambda].XY();
                tempMatrix(0, 2) = eigensystem.fEigentensors[lambda].XZ();
                tempMatrix(1, 0) = eigensystem.fEigentensors[lambda].XY();
                tempMatrix(1, 1) = eigensystem.fEigentensors[lambda].YY();
                tempMatrix(1, 2) = eigensystem.fEigentensors[lambda].YZ();
                tempMatrix(2, 0) = eigensystem.fEigentensors[lambda].XZ();
                tempMatrix(2, 1) = eigensystem.fEigentensors[lambda].YZ();
                tempMatrix(2, 2) = eigensystem.fEigentensors[lambda].ZZ();
                for (unsigned int i = 0; i < 3; ++i) {
                    int nonZeroIndex = -1;
                    for (unsigned int j = 0; j < 3; ++j) {
                        if (!IsZeroVal(tempMatrix(i, j))) {
                            nonZeroIndex = j;
                            break;
                        }
                    }
                    if (nonZeroIndex != -1) {
                        eigensystem.fEigenvectors[lambda][0] = tempMatrix(i, 0);
                        eigensystem.fEigenvectors[lambda][1] = tempMatrix(i, 1);
                        eigensystem.fEigenvectors[lambda][2] = tempMatrix(i, 2);
                        break;
                    }
                }
            }
            break;
        }
        default:
            DebugStop();
    }
}

template <class T>
void TPZTensor<T>::EigenProjection(const TPZVec<T> &EigenVals, int index, const TPZVec<int> &DistinctEigenvalues, TPZTensor<T> &Ei) const {

    const int p = DistinctEigenvalues.NElements();
    TPZFNMatrix<9, T> local(3, 3), aux(3, 3), resultingTensor(3, 3);
    Ei.Identity();
    for (int count = 0; count < p; ++count) {
        const int j = DistinctEigenvalues[count];
        if (j == index) continue;
        local.Identity();
        local *= -1. * EigenVals[j];
        this->CopyToTensor(aux);
        local += aux;
#ifdef PZDEBUG
        if (AreEqual(EigenVals[index], EigenVals[j])) DebugStop();
#endif
        local *= 1. / (EigenVals[index] - EigenVals[j]);
        Ei.CopyToTensor(aux);
        aux.Multiply(local, resultingTensor);
        Ei[_XX_] = resultingTensor(0, 0);
        Ei[_XY_] = resultingTensor(0, 1);
        Ei[_XZ_] = resultingTensor(0, 2);
        Ei[_YY_] = resultingTensor(1, 1);
        Ei[_YZ_] = resultingTensor(1, 2);
        Ei[_ZZ_] = resultingTensor(2, 2);
    }//for j
}

template <class T>
void TPZTensor<T>::HaighWestergaard(T &LodeAngle, T &qsi, T &rho) const {/*
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
        TPZTensor<T> & dLodeAngle, TPZTensor<T> & dQsi, TPZTensor<T> & dRho) const {
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

    LodeAngle = acos(J3 / J2 / sqrtJ2 * T(sqrt(27.) / 2.)) / T(3.);
    //	T pi23 = T(2. * M_PI / 3.);

    T temp = sqrt(T(3.) / (J2 * J2 * J2 * T(4.) - J3 * J3 * T(27.)));

    dLodeAngledJ2 = temp * J3 / J2 * T(3. / 2.);
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
void TPZTensor<T>::EigenValue(TPZTensor<T> &eigenval)const {
    const T I1(this->I1()),
            I2(this->I2()),
            I3(this->I3());

    const T R((T(-2.)*(I1 * I1 * I1) + T(9.) * I1 * I2 - T(27.) * I3) / T(54.));
    const T Q(((I1 * I1) - (T(3.) * I2)) / T(9.));

    T theta(0.);
    if (IsZeroVal(Q, 1e-6)) {
        if (IsZeroVal(R, 1e-6)) {
            theta = M_PI_2;
        } else {
            DebugStop();
        }
    } else {
        if (Q < 0.) DebugStop();
        T val(R * pow(Q, T(-1.5))); // i.e. R/sqrt(Q*Q*Q);
        if (val > +1.) val = +1.;
        if (val < -1.) val = -1.;
        theta = acos(val);
    }
    T sqrtQ(0.);
    if (Q < 0.) {
        if (IsZeroVal(Q)) {
            sqrtQ = T(0.);
        } else {
            DebugStop();
        }
    } else {
        sqrtQ = sqrt(Q);
    }

    T eigenval1 = T(-2.) * sqrtQ * cos(theta / T(3.)) + I1 / T(3.);
    T eigenval2 = T(-2.) * sqrtQ * cos((theta + T(2. * M_PI)) / T(3.)) + I1 / T(3.);
    T eigenval3 = T(-2.) * sqrtQ * cos((theta - T(2. * M_PI)) / T(3.)) + I1 / T(3.);

    if (eigenval2 > eigenval1) {
        T temp(eigenval1);
        eigenval1 = eigenval2;
        eigenval2 = temp;
    }
    if (eigenval3 > eigenval1) {
        T temp(eigenval1);
        eigenval1 = eigenval3;
        eigenval3 = temp;
    }
    if (eigenval3 > eigenval2) {
        T temp(eigenval2);
        eigenval2 = eigenval3;
        eigenval3 = temp;
    }
    eigenval.XX() = eigenval1;
    eigenval.YY() = eigenval2;
    eigenval.ZZ() = eigenval3;
}

template <class T>
void TPZTensor<T>::Lodeangle(TPZTensor<T> &GradLode,T &Lode)const
{
    T J2t(this->J2()),
            J3t(this->J3());

    if (fabs(TPZExtractVal::val(J2t)) < 1.e-6)J2t = T(1.e-6);
    T sqrtJ2t = sqrt(J2t);
    if (fabs(TPZExtractVal::val(sqrtJ2t)) < 1.e-6)sqrtJ2t = T(1.e-6);

    TPZTensor<T> dJ2t, dJ3t;

    this->dJ2(dJ2t);
    this->dJ3(dJ3t);
    // Derivatives with respect to I1, J2 and J3

    TPZTensor<T> dLodeAngle, TempTensor;

    //QUAL DOS DOIS?
    //T theta =-asin( ( T( 3.) * sqrt( T( 3.) ) * J3t ) /( T( 2.) *  sqrt(J2t*J2t*J2t) ) )/T( 3.);
    //O GRADIENTE DO LODE ESTA EM FUNCAO DO LODE DE 0 a Pi/3


    T lodetemp = (T(3.) * sqrt(T(3.)) * J3t) / (T(2.) * sqrt(J2t * J2t * J2t));
	if(TPZExtractVal::val(lodetemp) <= -1.)
	{
        lodetemp *= T(0.999); //	DebugStop();
        // TPZExtractVal::val(lodetemp) *= T(0.999);	//	DebugStop();
        //lodetemp = T(-1.);

    }


    //cout << "\n lodetemp "<<lodetemp<<endl;
    //cout << "\n TPZExtractVal::val(lodetemp);"<<TPZExtractVal::val(lodetemp)<<endl;

	if(TPZExtractVal::val(lodetemp) >= 1.)
	{
        lodetemp *= T(0.999);
        // TPZExtractVal::val(lodetemp)*= 0.999;
        // lodetemp = T(1.);
        //DebugStop();
    }

    //DLODE = (-2*Dj3*j2 + 3*Dj2*J3t)/(2.*pow(J2t,2.5)*sqrt(1.3333333333333333 - (9*pow(J3t,2))/pow(J2t,3)))
    //1
    T j33 = T(3.) * J3t;
    dJ2t.Multiply(j33, 1);
    //2
    T j22 = T(2.) * J2t;
    dJ3t.Multiply(j22, 1);
    //3
    dJ2t.Add(dJ3t, -1);
    //4
    //if(TPZExtractVal::val(J2t)<1.e-6)J2t=1.e-6;
    T checknegativeroot = ((T(9.) * J3t * J3t) / (J2t * J2t * J2t));
    if (TPZExtractVal::val(checknegativeroot) >= 4 / 3.)checknegativeroot *= T(0.999);
    T denom = T(2.) * sqrt(J2t * J2t * J2t * J2t * J2t) * sqrt((T(4 / 3.)) - checknegativeroot);
    //T denom2 = (T(2.)*pow(J2t,2.5)*sqrt(T(1.3333333333333333) - (T(9.)*pow(J3t,2.))/pow(J2t,3)));
    T oneoverden = T(1.) / denom;
    dJ2t.Multiply(oneoverden, 1);
    dJ2t *= oneoverden;
    //    GradLode = dJ2t;

    T acoslodetemp = acos(lodetemp);
    Lode = acoslodetemp / T(3.);

	if(TPZExtractVal::val(Lode) >= (M_PI/3.)-0.0001)
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

    if (fetestexcept(/*FE_DIVBYZERO*/ FE_ALL_EXCEPT)) {
        std::cout << "division by zero reported\n";
        DebugStop();
    }
    if (fetestexcept(FE_INVALID)) {
        std::cout << "invalid result reported\n";
        DebugStop();
    }
#endif
    //	T denominador;
    //	denominador = ( J2t*J2t*J2t * T(4.) - J3t*J3t * T(27.) );
    //	T temp;
    //	if(TPZExtractVal::val(denominador) < 1.e-6)
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
	 
     if(TPZExtractVal::val(temp33) <= -1.)
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
	 
     if(TPZExtractVal::val(lodetemp2) < -1.00000001)
     {
     std::cout << "\n TPZTensor LodeTemp < -1\n";
     std::cout << __PRETTY_FUNCTION__ << std::endl;
     DebugStop();
     }
	 
     if(TPZExtractVal::val(lodetemp2) > 1.000000001)
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

    if (fabs(TPZExtractVal::val(J2)) < 1.e-6)J2 = 1.e-6;
    //	if(fabs( TPZExtractVal::val(J2) ) < 1.e-6)return;

    T sqrtJ2 = sqrt(J2);
    if (fabs(TPZExtractVal::val(sqrtJ2)) < 1.e-6)sqrtJ2 = 1.e-6;

    TPZTensor<T> dJ2, dJ3, dLode;
    T Lode;

    this->dJ2(dJ2);
    this->dJ3(dJ3);
    this->Lodeangle(dLode, Lode);

    T pi23 = T(2. * M_PI / 3.);
    T TwoOverSqrThree = T(2. / sqrt(3.));
    T TwoOverSqrThreeJ2 = TwoOverSqrThree * sqrtJ2;
    T I13 = I1 / T(3.);

    T tempCosLode = cos(Lode) * TwoOverSqrThreeJ2;
    T tempCosMinusLode = cos(Lode - pi23) * TwoOverSqrThreeJ2;
    T tempCosPlusLode = cos(Lode + pi23) * TwoOverSqrThreeJ2;

	if(TPZExtractVal::val(Lode) < 0.)
	{
        cout << "Lode angle é Menor que ZERO. Valido somente para sig1 > sig2 > sig3 -> 0 < theta < Pi/3 " << endl;
        DebugStop();
    }
	if(TPZExtractVal::val(Lode) > M_PI/3.)
	{
        cout << "Lode angle é Maior que Pi/3. Valido somente para sig1 > sig2 > sig3 -> 0 < theta < Pi/3 " << endl;
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
        sout << "\n  TPZTENSOR \n" << endl;
        sout << "\n  LodeAngle = \n" << Lode << endl;
        sout << "\n  dLodeAngle= " << dLode << endl;
        sout << "\n\n";
        LOGPZ_INFO(loggerr, sout.str());
    }
#endif



    T OneOverTwoJ2 = T(0.5) / J2;
    TPZTensor<T> dI13;
    dI13.Identity();
    dI13 *= T(1. / 3.);

    tempCosLode *= OneOverTwoJ2;
    tempCosMinusLode *= OneOverTwoJ2;
    tempCosPlusLode *= OneOverTwoJ2;

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
    output << "XX = " << XX() << " XY = " << XY() << " XZ = " << XZ() << " YY = " << YY() << " YZ = " << YZ() << " ZZ = " << ZZ() << std::endl;
}

template <>
inline void TPZTensor<TFad<6, REAL> >::Print(std::ostream &output) const {
    output << "XX = " << XX() << "\nXY = " << XY() << "\nXZ = " << XZ() << "\nYY = " << YY() << "\nYZ = " << YZ() << "\nZZ = " << ZZ() << std::endl;
}

template <>
inline void TPZTensor<TFad<9, REAL> >::Print(std::ostream &output) const {
    output << "XX = " << XX() << "\nXY = " << XY() << "\nXZ = " << XZ() << "\nYY = " << YY() << "\nYZ = " << YZ() << "\nZZ = " << ZZ() << std::endl;
}

template <class T>
inline void TPZTensor<T>::Print(std::ostream &output) const {
    output << __PRETTY_FUNCTION__ << " please implement me\n";
}
//template class TPZTensor<REAL>;

#endif //TPZTENSOR_H
