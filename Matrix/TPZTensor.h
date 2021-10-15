// $Id: TPZTensor.h,v 1.28 2010-12-04 20:41:28 diogo Exp $

#ifndef TPZTENSOR_H
#define TPZTENSOR_H

#include "pzmanvector.h"
#include "pzfmatrix.h"
#include "pzextractval.h"
#include <iostream>
#include <algorithm>
#include "fadType.h"
#include <math.h>
#include <array>
#include <algorithm>
#include <functional>
#include <pzvec_extras.h>

#include "pzlog.h"
#include "TPZAssert.h"
#ifndef WIN32
#include <fenv.h>//NAN DETECTOR
#endif

#define _XX_ 0
#define _XY_ 1
#define _XZ_ 2
#define _YY_ 3
#define _YZ_ 4
#define _ZZ_ 5

#ifdef PZ_LOG
static TPZLogger loggerr("logtensor");
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
        TPZManVector<TPZTensor<T>, 3> fEigentensors; // Tensors of the spectral decomposition. If there is a repeated eigenvalue, all its occurrences are filled with the same tensor.

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
#ifdef PZ_LOG
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
#ifdef PZ_LOG
                if (geoMult == 0 || geoMult > 3) {
                    std::stringstream str;
                    str << "TPZTensor::Decomposed::Print Invalid geometric multiplicity (" << geoMult << ") for eigenvalue " << fEigenvalues[lambda] << std::endl;
                    LOGPZ_DEBUG(loggerr, str.str());
                }
#endif
                out << "Eigenvalue: " << fEigenvalues[lambda] << " (" << geoMult << " time" << (geoMult > 1 ? "s" : "") << ")" << std::endl;
                out << "Eigenvector" << (geoMult > 1 ? "s" : "") << ":" << std::endl;
                if (fEigenvectors.size() == 0 || fEigenvectors[0].size() == 0) {
                    out << " Not computed." << std::endl;
                } else {
                    for (unsigned int v = 0; v < geoMult; ++v) {
                        out << fEigenvectors[vIndex][0] << " ";
                        out << fEigenvectors[vIndex][1] << " ";
                        out << fEigenvectors[vIndex][2] << std::endl;
                        ++vIndex;
                    }
                }
                out << "Eigentensor:" << std::endl;
                fEigentensors[lambda].Print(out);
                out << std::endl;
                lambda += geoMult; // Go to next distinct eigenvalue
#ifdef PZDEBUG
                if (lambda > 3) {
                    std::stringstream str;
                    str << "TPZTensor::Decomposed::Print Total geometric multiplicity greater than expected!" << std::endl;
                    LOGPZ_ERROR(loggerr, str.str());
                }
#endif
            } while (lambda < 3);
        }
        /**
         * Methods for checking convergence
         */
        /// Number of test cases implemented by this class

        static int NumCasesCheckConv() {
            return 3;
        }

        static STATE gEigval[3];
        /// Compute the tangent matrix for a particular case

        static void TangentCheckConv(TPZFMatrix<STATE> &state, TPZFMatrix<STATE> &tangent, int icase) {
            TPZTensor<STATE> obj;
            for (int i = 0; i < 6; i++) {
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

        static void ResidualCheckConv(TPZFMatrix<STATE> &state, TPZFMatrix<STATE> &residual, int icase) {
            TPZTensor<STATE> obj;
            for (int i = 0; i < 6; i++) {
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

        void ComputeJ2(T & J2) {
            T sig1, sig2, sig3;
            sig1 = this->fEigenvalues[0];
            sig2 = this->fEigenvalues[1];
            sig3 = this->fEigenvalues[2];
            J2 = (pow(sig1 + (-sig1 - sig2 - sig3) / 3., 2) + pow(sig2 + (-sig1 - sig2 - sig3) / 3., 2) + pow((-sig1 - sig2 - sig3) / 3. + sig3, 2)) / 2.;
        }

        void ComputeJ3(T & J3) {
            T sig1, sig2, sig3, s1, s2, s3, temp;
            sig1 = this->fEigenvalues[0];
            sig2 = this->fEigenvalues[1];
            sig3 = this->fEigenvalues[2];
            s1 = sig1 - (1. / 3.)*(sig1 + sig2 + sig3);
            s2 = sig2 - (1. / 3.)*(sig1 + sig2 + sig3);
            s3 = sig3 - (1. / 3.)*(sig1 + sig2 + sig3);
            temp = (1. / 3.)*(s1 * s1 * s1 + s2 * s2 * s2 + s3 * s3 * s3);
            J3 = temp;
        }

        void ComputeI1(T &I1) {
            I1 = this->fEigenvalues[0] + this->fEigenvalues[1] + this->fEigenvalues[2];
        }

        void ApplyStrainComputeElasticStress(TPZTensor<T> &Stress, REAL &K, REAL & G) {
            const T sig1 = this->fEigenvalues[0];
            const T sig2 = this->fEigenvalues[1];
            const T sig3 = this->fEigenvalues[2];
            const T trace = sig1 + sig2 + sig3;
            const T trace_3 = trace / 3.;

            TPZTensor<T> SigVol;
            SigVol.XX() = SigVol.YY() = SigVol.ZZ() = K * trace;

            Stress.XX() = sig1 - trace_3;
            Stress.YY() = sig2 - trace_3;
            Stress.ZZ() = sig3 - trace_3;
            Stress *= (2 * G);
            Stress += SigVol;
        }
    };

    /**
     Construtor vazio inicializando com zero
     */
    TPZTensor() : fData(6, T(0.)) {
    }

    /**
     Construtor inicializando com Init
     */
    TPZTensor(const T & Init) : fData(6, Init) {
    }

    /**
     Copy Constructor
     */
    TPZTensor(const TPZTensor<T> & source) : fData(source.fData) {
    }

    /**
     Construct a tensor based on its eigensystem decomposition
     */
    TPZTensor(const TPZDecomposed &eigensystem) : fData(6, T(0.)) {
        unsigned int i;
        for (i = 0; i < 3; i += eigensystem.fGeometricMultiplicity[i]) {
            Add(eigensystem.fEigentensors[i], eigensystem.fEigenvalues[i]);
        }
#ifdef PZDEBUG
        if (i != 3) {
            DebugStop();
        }
#endif
    }

    int ClassId() const override  {
        return Hash("TPZTensor") ^ ClassIdOrHash<T>() << 1;
    }

    /// Method to write to a pzstream
    void Write(TPZStream &buf, int withclassid) const override;

    ///Method to read the object from a pzstream
    void Read(TPZStream &buf, void *context) override;

    operator TPZFMatrix<T>() const {
        TPZFMatrix<T> result(3, 3);
        result(0, 0) = XX();
        result(0, 1) = result(1, 0) = XY();
        result(0, 2) = result(2, 0) = XZ();
        result(1, 1) = YY();
        result(1, 2) = result(2, 1) = YZ();
        result(2, 2) = ZZ();
        return result;
    }

    TPZTensor(const TPZFMatrix<T> &input) : fData(6, 0.) {
#ifdef PZDEBUG
        if (input.Rows() != 3 || input.Cols() != 3) {
            DebugStop();
        }
        T tol;
        ZeroTolerance(tol);
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                if (fabs(TPZExtractVal::val(input.GetVal(i, j) - input.GetVal(j, i))) > tol) {
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

    /**
     Operator+
     */
    TPZTensor<T> operator+(const TPZTensor<T> &source) const;

    /**
     Operator-
     */
    TPZTensor<T> operator-(const TPZTensor<T> &source) const;

    /**
     Operator*
     */
    TPZTensor<T> operator*(const T &multipl) const;

    T &operator()(const int64_t row, const int64_t col) {
        if (row <= col) {
            return fData[row * (5 - row) / 2 + col];
        } else {
            return fData[col * (5 - col) / 2 + row];
        }
    }

    T &operator()(const int64_t row, const int64_t col) const {
        if (row <= col) {
            return fData[row * (5 - row) / 2 + col];
        } else {
            return fData[col * (5 - col) / 2 + row];
        }
    }

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
    void Multiply(const TPZTensor<T> tensor, TPZTensor<T> & resp)const;

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
    void Adjust(TPZVec<T> &sigIJ, TPZTensor<T> &result) const {
        TPZTensor<T> S;
        this->S(S);
        TPZTensor<T> Diag;
        Diag.Identity();
        Diag *= sigIJ[0] / T(3.);
        T J2 = this->J2();
        T sqj2 = sqrt(J2);
        REAL sqj2val = TPZExtractVal::val(sqj2);
        T ratio;
        if (fabs(sqj2val) > 1.e-6) {
            ratio = sigIJ[1] / sqj2;
        } else {
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

    void ComputeEigenvalues(TPZDecomposed &eigensystem, const bool compute_eigenvectors = false) const;


    /**
     * Returns the tensor eigenvalues and derivatives through an analytical approach
     */
    void Eigenvalue(TPZTensor<T> &eigenval, TPZTensor<T> &dSigma1, TPZTensor<T> &dSigma2, TPZTensor<T> &dSigma3)const;

    /**
     * Computes the Lode angle and its derivatives
     */
    void Lodeangle(TPZTensor<T> &GradLode, T &Lode)const;

    bool IsZeroTensor(T tol = 1.e-9) const {
        REAL realTol = TPZExtractVal::val(tol);
        for (unsigned int i = 0; i < 6; ++i) {
            if (fabs(this->fData[i]) > realTol) {
                return false;
            }
        }
        return true;
    }

    bool IsDiagonal(T tol = 1.e-9) const {
        REAL realTol = TPZExtractVal::val(tol);
        if ((fabs(this->XY()) > realTol) || (fabs(this->XZ()) > realTol) || (fabs(this->YZ()) > realTol)) {
            return false;
        }
        return true;
    }

    /**
     * Computes the eigenvectors and eigenvalues based on the Geometric Tools code.
     * https://www.geometrictools.com/
     */
    void EigenSystem(TPZDecomposed &eigensystem) const;

    /**
     * Computes the eigenvectors and eigenvalues based on Jacobi method
     */
    void EigenSystemJacobi(TPZDecomposed &eigensystem) const;

    /**
     Mnemonical access
     */
    inline T & XX() const {
        return fData[_XX_];
    }

    inline T & XY() const {
        return fData[_XY_];
    }

    inline T & XZ() const {
        return fData[_XZ_];
    }

    inline T & YY() const {
        return fData[_YY_];
    }

    inline T & YZ() const {
        return fData[_YZ_];
    }

    inline T & ZZ() const {
        return fData[_ZZ_];
    }

    inline T &operator[](int i) const {
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
    void ComputeEigenvectors(TPZDecomposed &eigensystem) const;


public:
    /**
     Dados do tensor
     */
    TPZManVector<T, 6> fData;

protected:

    static inline bool IsZeroVal(const T & val, T tol = 1.e-9) {
        return (fabs(TPZExtractVal::val(val)) <= tol);
    }

    bool AreEqual(const T &val1, const T &val2, const T tol = T(1.e-9)) const {
        return (std::fabs(TPZExtractVal::val(val1) - TPZExtractVal::val(val2)) <= tol);
    }

    /**
     * Returns the tensor eigenvalues through an analytical approach
     */
    void DirectEigenValues(TPZDecomposed &eigensystem, bool compute_eigenvectors) const;

    void Precondition(REAL &conditionFactor, TPZTensor<T>& A, TPZDecomposed &decomposition) const;

    void ComputeEigenvector0(const T &eigenvalue, TPZManVector<T, 3> &eigenvector) const;
    void ComputeEigenvector1(const TPZManVector<T, 3> &eigenvector0, const T &eigenvalue1, TPZManVector<T, 3> &eigenvector1) const;
    void ComputeEigenvectorsInternal(TPZDecomposed &eigensystem) const;

    /**
     * Computes the eigenprojections and eigenvalues based on the Geometric Tools code.
     * https://www.geometrictools.com/
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
inline void TPZTensor<T>::Identity_T(TBASE &) {
    fData[_XX_] = TBASE(1.);
    fData[_YY_] = TBASE(1.);
    fData[_ZZ_] = TBASE(1.);

    fData[_XY_] = TBASE(0.);
    fData[_XZ_] = TBASE(0.);
    fData[_YZ_] = TBASE(0.);
}

template <class T>
inline void TPZTensor<T>::Identity() {
    REAL TBase;
    Identity_T(TBase);
}

template < class T >
const TPZTensor<T> & TPZTensor<T>::operator=(const TPZTensor<T> &source) {
    fData = source.fData;
    return *this;
}

template < class T >
const TPZTensor<T> & TPZTensor<T>::operator+=(const TPZTensor<T> &source) {
    int i;
    for (i = 0; i < 6; i++)fData[i] += source.fData[i];
    return *this;
}

template < class T >
const TPZTensor<T> & TPZTensor<T>::operator-=(const TPZTensor<T> &source) {
    int i;
    for (i = 0; i < 6; i++)fData[i] -= source.fData[i];
    return *this;
}

template < class T >
const TPZTensor<T> & TPZTensor<T>::operator*=(const T &multipl) {
    int i;
    for (i = 0; i < 6; i++)fData[i] *= multipl;
    return *this;
}

template < class T >
TPZTensor<T> TPZTensor<T>::operator+(const TPZTensor<T> &source) const {
    TPZTensor<T> temp(*this);
    return temp += source;
}

template < class T >
TPZTensor<T> TPZTensor<T>::operator-(const TPZTensor<T> &source) const {
    TPZTensor<T> temp(*this);
    return temp -= source;
}

template < class T >
TPZTensor<T> TPZTensor<T>::operator*(const T &multipl) const {
    TPZTensor<T> temp(*this);
    return temp *= multipl;
}

template < class T >
template < class T1, class T2 >
void TPZTensor<T>::Add(const TPZTensor< T1 > & tensor, const T2 & constant) {
    int i, size = 6;
    for (i = 0; i < size; i++) {
        fData[i] += tensor.fData[i] * T1(constant);
    }
}

template < class T >
template < class T1, class T2 >
void TPZTensor<T>::Multiply(const T1 & multipl, const T2 & constant) {
    int i, size = 6;
    for (i = 0; i < size; i++) {
        fData[i] *= multipl;
        fData[i] *= constant;
    }
}

template < class T >
void TPZTensor<T>::Multiply(const TPZTensor<T> tensor, TPZTensor<T> & resp) const {
    const T XX = this->fData[_XX_];
    const T YY = this->fData[_YY_];
    const T ZZ = this->fData[_ZZ_];
    const T XY = this->fData[_XY_];
    const T XZ = this->fData[_XZ_];
    const T YZ = this->fData[_YZ_];
    resp.fData[_XX_] = XX * tensor.XX() + XY * tensor.XY() + XZ * tensor.XZ();
    resp.fData[_YY_] = XY * tensor.XY() + YY * tensor.YY() + YZ * tensor.YZ();
    resp.fData[_ZZ_] = XZ * tensor.XZ() + YZ * tensor.YZ() + ZZ * tensor.ZZ();
    resp.fData[_XY_] = XX * tensor.XY() + XY * tensor.YY() + XZ * tensor.YZ();
    resp.fData[_XZ_] = XX * tensor.XZ() + XY * tensor.YZ() + XZ * tensor.ZZ();
    resp.fData[_YZ_] = XY * tensor.XZ() + YY * tensor.YZ() + YZ * tensor.ZZ();
}

template < class T >
template < class T2 >
void TPZTensor<T>::Scale(const T2 & constant) {
    int i, size = 6;
    const T Tconst = T(constant);
    for (i = 0; i < size; i++) {
        fData[i] *= Tconst;
    }
}

template < class T >
void TPZTensor<T>::Zero() {
    int i, size = 6;
    const T Tzero = T(0.);
    for (i = 0; i < size; i++) {
        fData[i] = Tzero;
    }
}

template < class T >
void TPZTensor<T>::SetUp(const TPZVec<REAL> & Solution) {
    int i, size = 6;
    for (i = 0; i < size; i++) {
        fData[i] = Solution[i];
    }
}

template < class T >
template < class T1 >
void TPZTensor<T>::CopyTo(TPZTensor<T1> & target) const {
    for (unsigned int i = 0; i < 6; i++) {
        target[i] = TPZExtractVal::val(fData[i]);
    }
}

template < class T >
void TPZTensor<T>::CopyTo(TPZFMatrix<REAL> & target) const {
    for (unsigned int i = 0; i < 6; i++) {
        target(i, 0) = TPZExtractVal::val(fData[i]);
    }
}

template < class T >
void TPZTensor<T>::CopyFrom(const TPZFMatrix<T> & source) {
    for (unsigned int i = 0; i < 6; i++) {
        fData[i] = source.Get(i, 0);
    }
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
T TPZTensor<T>::Norm() const {
    T norm = T(0.);
    for (unsigned int i = 0; i < 6; i++) {
        norm += fData[i] * fData[i];
    }
    norm += fData[_XY_] * fData[_XY_];
    norm += fData[_XZ_] * fData[_XZ_];
    norm += fData[_YZ_] * fData[_YZ_];
    return sqrt(norm);
}

template <class T>
std::ostream &operator<<(std::ostream &out, const TPZTensor<T> &tens) {
    tens.Print(out);
    return out;
}

/**
 * Metodo que calcula o determinante do tensor
 */
template <class T>
T TPZTensor<T>::Det() const {
    return fData[_XX_] * fData[_YY_] * fData[_ZZ_] + fData[_XY_] * fData[_XZ_] * fData[_YZ_]*2. - fData[_XZ_] * fData[_YY_] * fData[_XZ_] -
            fData[_XY_] * fData[_XY_] * fData[_ZZ_] - fData[_YZ_] * fData[_YZ_] * fData[_XX_];
}

/**
 * Metodo que calcula a derivada do determinante do tensor
 */
template <class T>
void TPZTensor<T>::dDet(TPZTensor<T> &grad) const {
    grad.fData[_XX_] = fData[_YY_] * fData[_ZZ_] - fData[_YZ_] * fData[_YZ_];
    grad.fData[_XY_] = fData[_XZ_] * fData[_YZ_] - fData[_XY_] * fData[_ZZ_];
    grad.fData[_XZ_] = fData[_XY_] * fData[_YZ_] - fData[_XZ_] * fData[_YY_];
    grad.fData[_YY_] = fData[_XX_] * fData[_ZZ_] - fData[_XZ_] * fData[_XZ_];
    grad.fData[_YZ_] = fData[_XY_] * fData[_XZ_] - fData[_YZ_] * fData[_XX_];
    grad.fData[_ZZ_] = fData[_XX_] * fData[_YY_] - fData[_XY_] * fData[_XY_];
}

template <class T>
void TPZTensor<T>::S(TPZTensor<T> &s) const {
    s = *this;
    T mult = -I1() / T(3);
    s.fData[_XX_] += mult;
    s.fData[_YY_] += mult;
    s.fData[_ZZ_] += mult;
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
T TPZTensor<T>::J2() const {
    TPZVec<T> s(3);
    DeviatoricDiagonal_T<T>(s);
    T value = -s[0] * s[1] - s[0] * s[2] - s[1] * s[2] + fData[_XY_] * fData[_XY_] + fData[_XZ_] * fData[_XZ_] + fData[_YZ_] * fData[_YZ_];
    return TPZAssert::NonNegative(value);
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

template <class T>
T TPZTensor<T>::J3() const {
    TPZTensor<T> s;
    S(s);
    return s.Det();
}

/**
 * The derivative of the  3rd invariant of deviatoric tensor
 */
template<class T>
void TPZTensor<T>::dJ3(TPZTensor<T> &deriv) const {
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
    Tensor(0, 0) = XX();
    Tensor(0, 1) = Tensor(1, 0) = XY();
    Tensor(0, 2) = Tensor(2, 0) = XZ();
    Tensor(1, 1) = YY();
    Tensor(1, 2) = Tensor(2, 1) = YZ();
    Tensor(2, 2) = ZZ();
}

template<class T>
static T Polynom(const T &x, const T &I1, const T &I2, const T &I3) {
    T result = x * x * x - I1 * x * x + I2 * x - I3;
    return result;
}

template<class T>
static T DerivPolynom(const T &x, const T &I1, const T &I2, const T &I3) {
    T result = T(3.) * x * x - T(2.) * I1 * x + I2;
    return result;
}

template<class T>
static T UpdateNewton(const T &x, const T &I1, const T &I2, const T &I3) {
    T residue = Polynom(x, I1, I2, I3);
    T dres = DerivPolynom(x, I1, I2, I3);
    return x - residue / dres;
}

template <class T>
void TPZTensor<T>::EigenSystemJacobi(TPZDecomposed &eigensystem)const {
    TPZFNMatrix<9, T> TensorMat(*this);

    int64_t numiterations = 1000;
    T tol;
    ZeroTolerance(tol);

    eigensystem.fEigenvectors[0] = TPZManVector<T, 3>(3, 0.);
    eigensystem.fEigenvectors[1] = TPZManVector<T, 3>(3, 0.);
    eigensystem.fEigenvectors[2] = TPZManVector<T, 3>(3, 0.);

    TPZManVector<T, 3> Eigenvalues(3, 0.);
    TPZFNMatrix<9, T> Eigenvectors(3, 3, 0.);
    TensorMat.SolveEigensystemJacobi(numiterations, TPZExtractVal::ref(tol), Eigenvalues, Eigenvectors);

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

#ifdef PZDEBUG
#ifdef PZ_LOG
    if (loggerr.isDebugEnabled()) {
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
void TPZTensor<T>::ComputeEigenvalues(TPZDecomposed &eigensystem, const bool compute_eigenvectors) const {
    TPZManVector<T, 3> &eigenvalues = eigensystem.fEigenvalues;
    TPZManVector<TPZManVector<T, 3>, 3> &eigenvectors = eigensystem.fEigenvectors;

    T tol = Norm()*1.e-12;

    if (this->IsZeroTensor(tol)) {
        eigensystem.fDistinctEigenvalues = 1;
        eigenvalues[0] = eigenvalues[1] = eigenvalues[2] = 0.;
        eigensystem.fGeometricMultiplicity[0] = 3;
        eigensystem.fGeometricMultiplicity[1] = 3;
        eigensystem.fGeometricMultiplicity[2] = 3;
        if (compute_eigenvectors) {
            eigenvectors[0] = TPZManVector<T, 3>(3, 0.);
            eigenvectors[0][0] = 1.0;
            eigenvectors[1] = TPZManVector<T, 3>(3, 0.);
            eigenvectors[1][1] = 1.0;
            eigenvectors[2] = TPZManVector<T, 3>(3, 0.);
            eigenvectors[2][2] = 1.0;
        }
    } else if (this->IsDiagonal(tol)) {
        eigenvalues[0] = this->XX();
        eigenvalues[1] = this->YY();
        eigenvalues[2] = this->ZZ();

        if (compute_eigenvectors) {
            eigenvectors.resize(3);
            eigenvectors[0] = TPZManVector<T, 3>(3, 0.);
            eigenvectors[0][0] = 1.0;
            eigenvectors[1] = TPZManVector<T, 3>(3, 0.);
            eigenvectors[1][1] = 1.0;
            eigenvectors[2] = TPZManVector<T, 3>(3, 0.);
            eigenvectors[2][2] = 1.0;
        }

        bool ev0eqev1 = AreEqual(eigenvalues[0], eigenvalues[1], tol);
        bool ev1eqev2 = AreEqual(eigenvalues[1], eigenvalues[2], tol);

        if (ev0eqev1 && ev1eqev2) {
            //tres autovalores iguais
            eigensystem.fDistinctEigenvalues = 1;
            eigensystem.fGeometricMultiplicity[0] = 3;
            eigensystem.fGeometricMultiplicity[1] = 3;
            eigensystem.fGeometricMultiplicity[2] = 3;
            return;
        } else if (ev0eqev1 || ev1eqev2 || AreEqual(eigenvalues[0], eigenvalues[2], tol)) {
            eigensystem.fDistinctEigenvalues = 2;
            int different = -1;
            int equals[2] = {-1, -1};
            if (ev0eqev1) {
                different = 2;
                equals[0] = 0;
                equals[1] = 1;
            } else if (ev1eqev2) {
                different = 0;
                equals[0] = 1;
                equals[1] = 2;
            } else {
                different = 1;
                equals[0] = 0;
                equals[1] = 2;
            }
            eigensystem.fGeometricMultiplicity[different] = 1;
            eigensystem.fGeometricMultiplicity[equals[0]] = 2;
            eigensystem.fGeometricMultiplicity[equals[1]] = 2;
        } else {
            eigensystem.fDistinctEigenvalues = 3;
            eigensystem.fGeometricMultiplicity[0] = 1;
            eigensystem.fGeometricMultiplicity[1] = 1;
            eigensystem.fGeometricMultiplicity[2] = 1;
        }

        if (eigenvalues[0] < eigenvalues[1]) {
            T eigenvalueTemp = eigenvalues[0];
            eigenvalues[0] = eigenvalues[1];
            eigenvalues[1] = eigenvalueTemp;
            std::swap(eigensystem.fGeometricMultiplicity[0], eigensystem.fGeometricMultiplicity[1]);
            if (compute_eigenvectors) {
                auto temp = eigenvectors[0];
                eigenvectors[0] = eigenvectors[1];
                eigenvectors[1] = temp;
            }
        }
        if (eigenvalues[1] < eigenvalues[2]) {
            T eigenvalueTemp = eigenvalues[1];
            eigenvalues[1] = eigenvalues[2];
            eigenvalues[2] = eigenvalueTemp;
            std::swap(eigensystem.fGeometricMultiplicity[1], eigensystem.fGeometricMultiplicity[2]);
            if (compute_eigenvectors) {
                auto temp = eigenvectors[1];
                eigenvectors[1] = eigenvectors[2];
                eigenvectors[2] = temp;
            }
        }
        if (eigenvalues[0] < eigenvalues[1]) {
            T eigenvalueTemp = eigenvalues[0];
            eigenvalues[0] = eigenvalues[1];
            eigenvalues[1] = eigenvalueTemp;
            std::swap(eigensystem.fGeometricMultiplicity[0], eigensystem.fGeometricMultiplicity[1]);
            if (compute_eigenvectors) {
                auto temp = eigenvectors[0];
                eigenvectors[0] = eigenvectors[1];
                eigenvectors[1] = temp;
            }
        }
    } else {
        this->DirectEigenValues(eigensystem, compute_eigenvectors);

        bool ev0eqev1 = AreEqual(eigenvalues[0], eigenvalues[1], tol);
        bool ev1eqev2 = AreEqual(eigenvalues[1], eigenvalues[2], tol);
        //tres autovalores iguais
        if (ev0eqev1 && ev1eqev2) {
            eigensystem.fDistinctEigenvalues = 1;
            eigensystem.fGeometricMultiplicity[0] = 3;
            eigensystem.fGeometricMultiplicity[1] = 3;
            eigensystem.fGeometricMultiplicity[2] = 3;
        } else if (ev0eqev1 || ev1eqev2) {
            //dois autovalores iguais
            eigensystem.fDistinctEigenvalues = 2;

            int different = -1;
            int equals = -1;
            if (ev0eqev1) {
                different = 2;
                equals = 0;
            } else {
                different = 0;
                equals = 2;
            }
            eigensystem.fGeometricMultiplicity[1] = 2;
            eigensystem.fGeometricMultiplicity[equals] = 2;
            eigensystem.fGeometricMultiplicity[different] = 1;
        } else {
            //3 autovalores diferentes
            eigensystem.fDistinctEigenvalues = 3;
            eigensystem.fGeometricMultiplicity[0] = 1;
            eigensystem.fGeometricMultiplicity[1] = 1;
            eigensystem.fGeometricMultiplicity[2] = 1;
        }
    }
}

template <class T>
void TPZTensor<T>::EigenSystem(TPZDecomposed &eigensystem)const {
    TPZManVector<T, 3> &eigenvalues = eigensystem.fEigenvalues;
    TPZManVector<TPZTensor<T>, 3> &eigentensors = eigensystem.fEigentensors;

    T precision;
    ZeroTolerance(precision);
    T tol = Norm()*precision;
    ComputeEigenvectors(eigensystem);

    TPZStack<int, 3> indices_to_add;
    for (unsigned int i = 0; i < 3; ++i) {
        if (eigensystem.fGeometricMultiplicity[i] != 1) {
            indices_to_add.push_back(i);
        }
        for (unsigned int j = 0; j < 3; ++j) {
            for (unsigned int k = j; k < 3; ++k) {
                eigentensors[i](j, k) = eigensystem.fEigenvectors[i][j] * eigensystem.fEigenvectors[i][k];
            }
        }
    }
    if (indices_to_add.size() > 0) {
        TPZTensor<T> sum;
        for (auto i : indices_to_add) {
            sum += eigentensors[i];
        }
        for (auto i : indices_to_add) {
            eigentensors[i] = sum;
        }
    }
#ifdef PZDEBUG
    TPZTensor<T> total;
    unsigned int i;
    for (i = 0; i < 3; i += eigensystem.fGeometricMultiplicity[i]) {
        total.Add(eigentensors[i], eigenvalues[i]);
    }
    if (i != 3) {
        std::cout << "Tensor decomposition error: " << std::endl;
        std::cout << "Incorrect total geometric multiplicity: " << i << std::endl;
        DebugStop();
    }
    for (unsigned int i = 0; i < 6; ++i) {
        if (!AreEqual(total[i], this->operator[](i), tol)) {
            std::cout << std::setprecision(15);
            std::cout << "Tensor decomposition error: " << std::endl;
            std::cout << "Original Tensor: ";
            this->Print(std::cout);
            std::cout << "Reconstruction from decomposition: ";
            total.Print(std::cout);
            std::cout << "Decomposition: " << std::endl;
            eigensystem.Print(std::cout);
            DebugStop();
        }
    }
#endif
}

template <class T>
void TPZTensor<T>::ComputeEigenvectors(TPZDecomposed &eigensystem) const {
    // Eigenvalues not computed yet. Let's do it.
    if (eigensystem.fDistinctEigenvalues == 0) {
        ComputeEigenvalues(eigensystem, true);
    } else {
        eigensystem.fEigenvectors[0] = TPZManVector<T, 3>(3, 0.);
        eigensystem.fEigenvectors[1] = TPZManVector<T, 3>(3, 0.);
        eigensystem.fEigenvectors[2] = TPZManVector<T, 3>(3, 0.);
        switch (eigensystem.fDistinctEigenvalues) {
            case 1:
                eigensystem.fEigenvectors[0][0] = 1.;
                eigensystem.fEigenvectors[1][1] = 1.;
                eigensystem.fEigenvectors[2][2] = 1.;
                return;
            case 2:
            case 3:
            {
                T tol = Norm()*1.e-12;
                if (IsDiagonal(tol)) {
                    for (unsigned int i = 0; i < 3; ++i) {
                        unsigned int j;
                        for (j = 0; j < 3; ++j) {
                            if (AreEqual(eigensystem.fEigenvalues[i], this->operator()(j, j), tol)) {
                                eigensystem.fEigenvectors[i][j] = 1;
                                break;
                            }
                        }
                        if (j == 3) {
                            DebugStop();
                        }
                    }
                } else {
                    REAL conditionFactor;
                    TPZDecomposed eigensystemTemp = eigensystem;
                    TPZTensor<T> A;
                    Precondition(conditionFactor, A, eigensystemTemp);
                    A.ComputeEigenvectorsInternal(eigensystemTemp);
                    eigensystem.fEigenvectors = eigensystemTemp.fEigenvectors;
                }
                break;
            }
            default:
                DebugStop();
        }
    }


    for (unsigned int lambda = 0; lambda < 3; ++lambda) {
        REAL Norm = 0.;
        for (unsigned int i = 0; i < 3; ++i) {
            Norm += pow(TPZExtractVal::val(eigensystem.fEigenvectors[lambda][i]), 2);
        }
        Norm = sqrt(Norm);
        for (unsigned int i = 0; i < 3; ++i) {
            eigensystem.fEigenvectors[lambda][i] /= Norm;
        }
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
void TPZTensor<T>::ComputeEigenvector0(const T &eigenvalue, TPZManVector<T, 3> &eigenvector) const {
    // Compute a unit-length eigenvector for eigenvalue[i0].  The matrix is
    // rank 2, so two of the rows are linearly independent.  For a robust
    // computation of the eigenvector, select the two rows whose cross product
    // has largest length of all pairs of rows.
    TPZManVector<T, 3> row0(3);
    row0[0] = XX() - eigenvalue;
    row0[1] = XY();
    row0[2] = XZ();
    TPZManVector<T, 3> row1(3);
    row1[0] = XY();
    row1[1] = YY() - eigenvalue;
    row1[2] = YZ();
    TPZManVector<T, 3> row2(3);
    row2[0] = XZ();
    row2[1] = YZ();
    row2[2] = ZZ() - eigenvalue;
    TPZManVector<T, 3> r0xr1(3);
    Cross(row0, row1, r0xr1);
    TPZManVector<T, 3> r0xr2(3);
    Cross(row0, row2, r0xr2);
    TPZManVector<T, 3> r1xr2(3);
    Cross(row1, row2, r1xr2);
    T d0 = Dot(r0xr1, r0xr1);
    T d1 = Dot(r0xr2, r0xr2);
    T d2 = Dot(r1xr2, r1xr2);

    REAL dmax = TPZExtractVal::val(d0);
    int imax = 0;
    if (TPZExtractVal::val(d1) > dmax) {
        dmax = TPZExtractVal::val(d1);
        imax = 1;
    }
    if (TPZExtractVal::val(d2) > dmax) {
        imax = 2;
    }

    if (imax == 0) {
        T inv_sqrt_val = 1. / sqrt(d0);
        eigenvector[0] = r0xr1[0] * inv_sqrt_val;
        eigenvector[1] = r0xr1[1] * inv_sqrt_val;
        eigenvector[2] = r0xr1[2] * inv_sqrt_val;
    } else if (imax == 1) {
        T inv_sqrt_val = 1. / sqrt(d1);
        eigenvector[0] = r0xr2[0] * inv_sqrt_val;
        eigenvector[1] = r0xr2[1] * inv_sqrt_val;
        eigenvector[2] = r0xr2[2] * inv_sqrt_val;
    } else {
        T inv_sqrt_val = 1. / sqrt(d2);
        eigenvector[0] = r1xr2[0] * inv_sqrt_val;
        eigenvector[1] = r1xr2[1] * inv_sqrt_val;
        eigenvector[2] = r1xr2[2] * inv_sqrt_val;
    }
}

template <class T>
void TPZTensor<T>::ComputeEigenvector1(const TPZManVector<T, 3> &eigenvector0, const T &eigenvalue1, TPZManVector<T, 3> &eigenvector1) const {
    // Robustly compute a right-handed orthonormal set { U, V, evec0 }.
    TPZManVector<T, 3> U(3), V(3);

    // The vector eigenvector0 is guaranteed to be unit-length, in which case there is no
    // need to worry about a division by zero when computing invLength.
    T invLength;
    if (fabs(eigenvector0[0]) > fabs(eigenvector0[1])) {
        // The component of maximum absolute value is either eigenvector0[0] or eigenvector0[2].
        invLength = (T) 1 / sqrt(eigenvector0[0] * eigenvector0[0] + eigenvector0[2] * eigenvector0[2]);
        U[0] = -eigenvector0[2] * invLength;
        U[1] = T(0.);
        U[2] = +eigenvector0[0] * invLength;
    } else {
        // The component of maximum absolute value is either eigenvector0[1] or eigenvector0[2].
        invLength = (T) 1 / sqrt(eigenvector0[1] * eigenvector0[1] + eigenvector0[2] * eigenvector0[2]);
        U[0] = T(0.);
        U[1] = +eigenvector0[2] * invLength;
        U[2] = -eigenvector0[1] * invLength;
    }
    Cross(eigenvector0, U, V);

    // Let e be eigenvalue1 and let v1 be a corresponding eigenvector which is a
    // solution to the linear system (A - e*I)*v1 = 0.  The matrix (A - e*I)
    // is 3x3, not invertible (so infinitely many solutions), and has rank 2
    // when eigenvalue1 and eigenvalue2 are different.  It has rank 1 when eigenvalue1 and eigenvalue2
    // are equal.  Numerically, it is difficult to compute robustly the rank
    // of a matrix.  Instead, the 3x3 linear system is reduced to a 2x2 system
    // as follows.  Define the 3x2 matrix J = [U V] whose columns are the U
    // and V computed previously.  Define the 2x1 vector X = J*v1.  The 2x2
    // system is 0 = M * X = (J^T * (A - e*I) * J) * X where J^T is the
    // transpose of J and M = J^T * (A - e*I) * J is a 2x2 matrix.  The system
    // may be written as
    //     +-                        -++-  -+       +-  -+
    //     | U^T*A*U - e  U^T*A*V     || x0 | = e * | x0 |
    //     | V^T*A*U      V^T*A*V - e || x1 |       | x1 |
    //     +-                        -++   -+       +-  -+
    // where X has row entries x0 and x1.
    TPZManVector<T, 3> AU(3);
    AU[0] = XX() * U[0] + XY() * U[1] + XZ() * U[2];
    AU[1] = XY() * U[0] + YY() * U[1] + YZ() * U[2];
    AU[2] = XZ() * U[0] + YZ() * U[1] + ZZ() * U[2];

    TPZManVector<T, 3> AV(3);
    AV[0] = XX() * V[0] + XY() * V[1] + XZ() * V[2];
    AV[1] = XY() * V[0] + YY() * V[1] + YZ() * V[2];
    AV[2] = XZ() * V[0] + YZ() * V[1] + ZZ() * V[2];

    T m00 = U[0] * AU[0] + U[1] * AU[1] + U[2] * AU[2] - eigenvalue1;
    T m01 = U[0] * AV[0] + U[1] * AV[1] + U[2] * AV[2];
    T m11 = V[0] * AV[0] + V[1] * AV[1] + V[2] * AV[2] - eigenvalue1;

    // For robustness, choose the largest-length row of M to compute the
    // eigenvector.  The 2-tuple of coefficients of U and V in the
    // assignments to eigenvector[1] lies on a circle, and U and V are
    // unit length and perpendicular, so eigenvector[1] is unit length
    // (within numerical tolerance).
    REAL absM00 = fabs(TPZExtractVal::val(m00));
    REAL absM01 = fabs(TPZExtractVal::val(m01));
    REAL absM11 = fabs(TPZExtractVal::val(m11));
    REAL maxAbsComp;
    if (absM00 >= absM11) {
        maxAbsComp = std::max(absM00, absM01);
        if (IsZeroVal(maxAbsComp)) {
            eigenvector1 = U;
        } else {
            if (absM00 >= absM01) {
                m01 /= m00;
                m00 = T(1.) / sqrt(T(1.) + m01 * m01);
                m01 *= m00;
            } else {
                m00 /= m01;
                m01 = T(1.) / sqrt(T(1.) + m00 * m00);
                m00 *= m01;
            }
            //eigenvector1 = m01*U - m00*V
            sscal(U, m01);
            eigenvector1 = U;
            saxpy(eigenvector1, V, -m00);
        }
    } else {
        maxAbsComp = std::max(absM11, absM01);
        if (IsZeroVal(maxAbsComp)) {
            eigenvector1 = U;
        } else {
            if (absM11 >= absM01) {
                m01 /= m11;
                m11 = T(1.) / sqrt(T(1.) + m01 * m01);
                m01 *= m11;
            } else {
                m11 /= m01;
                m01 = T(1.) / sqrt(T(1.) + m11 * m11);
                m11 *= m01;
            }
            //eigenvector1 = m11*U - m01*V
            sscal(U, m11);
            eigenvector1 = U;
            saxpy(eigenvector1, V, -m01);
        }
    }
}

template <class T>
void TPZTensor<T>::ComputeEigenvectorsInternal(TPZDecomposed &eigensystem) const {
    TPZManVector<T, 3> &eigenval = eigensystem.fEigenvalues;
    TPZManVector<TPZManVector<T, 3>, 3> &eigenvec = eigensystem.fEigenvectors;

    eigenvec.resize(3);
    eigenvec[0].resize(3);
    eigenvec[1].resize(3);
    eigenvec[2].resize(3);
    // Compute the eigenvectors so that the set {evec[0], evec[1], evec[2]}
    // is right handed and orthonormal.
    if (eigensystem.fGeometricMultiplicity[0] == 1) {
        ComputeEigenvector0(eigenval[0], eigenvec[0]);
        ComputeEigenvector1(eigenvec[0], eigenval[1], eigenvec[1]);
        Cross(eigenvec[0], eigenvec[1], eigenvec[2]);
    } else {
        ComputeEigenvector0(eigenval[2], eigenvec[2]);
        ComputeEigenvector1(eigenvec[2], eigenval[1], eigenvec[1]);
        Cross(eigenvec[1], eigenvec[2], eigenvec[0]);
    }
}

template <class T>
void TPZTensor<T>::Precondition(REAL &conditionFactor, TPZTensor<T>& A, TPZDecomposed &decomposition) const {
    // Precondition the matrix by factoring out the maximum absolute value
    // of the components.  This guards against floating-point overflow when
    // computing the eigenvalues.
    REAL max0 = std::max(fabs(TPZExtractVal::val(XX())), fabs(TPZExtractVal::val(XY())));
    REAL max1 = std::max(fabs(TPZExtractVal::val(XZ())), fabs(TPZExtractVal::val(YY())));
    REAL max2 = std::max(fabs(TPZExtractVal::val(YZ())), fabs(TPZExtractVal::val(ZZ())));
    conditionFactor = std::max(std::max(max0, max1), max2);

    REAL invMaxAbsElement = 1. / conditionFactor;
    for (unsigned int i = 0; i < 6; ++i) {
        A.fData[i] = this->fData[i] * invMaxAbsElement;
    }
    if (decomposition.fDistinctEigenvalues != 0) {
        decomposition.fEigenvalues[0] *= invMaxAbsElement;
        decomposition.fEigenvalues[1] *= invMaxAbsElement;
        decomposition.fEigenvalues[2] *= invMaxAbsElement;
    }
}

template <class T>
void TPZTensor<T>::DirectEigenValues(TPZDecomposed &eigensystem, bool compute_eigenvectors) const {
    TPZManVector<T, 3> &eigenval = eigensystem.fEigenvalues;

    REAL conditionFactor;
    TPZTensor<T> A;
    Precondition(conditionFactor, A, eigensystem);

    T norm = pow(A.XY(), 2.) + pow(A.XZ(), 2.) + pow(A.YZ(), 2.);
    // Compute the eigenvalues.  The acos(z) function requires |z| <= 1,
    // but will fail silently and return NaN if the input is larger than 1 in
    // magnitude.  To avoid this condition due to rounding errors, the halfDet
    // value is clamped to [-1,1].
    T traceDiv3 = A.I1() / 3.;
    T b00 = A.XX() - traceDiv3;
    T b11 = A.YY() - traceDiv3;
    T b22 = A.ZZ() - traceDiv3;
    T denom = sqrt((pow(b00, 2.) + pow(b11, 2.) + pow(b22, 2.) + norm * T(2.)) / T(6.));
    T c00 = b11 * b22 - A.YZ() * A.YZ();
    T c01 = A.XY() * b22 - A.YZ() * A.XZ();
    T c02 = A.XY() * A.YZ() - b11 * A.XZ();
    T det = (b00 * c00 - A.XY() * c01 + A.XZ() * c02) / (denom * denom * denom);
    T halfDet = det * T(0.5);
    halfDet = std::min(std::max(TPZExtractVal::val(halfDet), -1.), 1.);

    // The eigenvalues of B are ordered as beta0 <= beta1 <= beta2.  The
    // number of digits in twoThirdsPi is chosen so that, whether float or
    // double, the floating-point number is the closest to theoretical 2*pi/3.
    T angle = acos(halfDet) / T(3.);
    const T twoThirdsPi = T(2.09439510239319549);
    T beta2 = cos(angle) * T(2.);
    T beta0 = cos(angle + twoThirdsPi) * T(2.);
    T beta1 = -(beta0 + beta2);

    // The eigenvalues are ordered as alpha0 >= alpha1 >= alpha2.
    eigenval[0] = traceDiv3 + denom * beta2;
    eigenval[1] = traceDiv3 + denom * beta1;
    eigenval[2] = traceDiv3 + denom * beta0;

    if (eigenval[0] < eigenval[1]) {
        DebugStop();
    }

    if (eigenval[1] < eigenval[2]) {
        DebugStop();
    }

    if (halfDet > 0. || IsZeroVal(halfDet)) { // greatest eigenvalue has multiplicity 1
        eigensystem.fGeometricMultiplicity[0] = 1;
    } else { // lowest eigenvalue has multiplicity 1
        eigensystem.fGeometricMultiplicity[2] = 1;
    }

    if (compute_eigenvectors) {
        A.ComputeEigenvectorsInternal(eigensystem);
    }
    // The preconditioning scaled the tensor, which scales the eigenvalues.
    // Revert the scaling.
    eigenval[0] *= conditionFactor;
    eigenval[1] *= conditionFactor;
    eigenval[2] *= conditionFactor;
}

template <class T>
void TPZTensor<T>::Lodeangle(TPZTensor<T> &GradLode, T &Lode)const {
    T J2t(this->J2());
    T J3t(this->J3());

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
    if (TPZExtractVal::val(lodetemp) <= -1.) {
        lodetemp *= T(0.999); //	DebugStop();
        // TPZExtractVal::val(lodetemp) *= T(0.999);	//	DebugStop();
        //lodetemp = T(-1.);

    }


    //cout << "\n lodetemp "<<lodetemp<<endl;
    //cout << "\n TPZExtractVal::val(lodetemp);"<<TPZExtractVal::val(lodetemp)<<endl;

    if (TPZExtractVal::val(lodetemp) >= 1.) {
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

    if (TPZExtractVal::val(Lode) >= (M_PI / 3.) - 0.0001) {
        Lode *= T(0.999);
    }
}

template <class T>
void TPZTensor<T>::Eigenvalue(TPZTensor<T> &eigenval, TPZTensor<T> &dSigma1, TPZTensor<T> &dSigma2, TPZTensor<T> &dSigma3)const {
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

    if (TPZExtractVal::val(Lode) < 0.) {
        std::cout << "Lode angle é Menor que ZERO. Valido somente para sig1 > sig2 > sig3 -> 0 < theta < Pi/3 " << std::endl;
        DebugStop();
    }
    if (TPZExtractVal::val(Lode) > M_PI / 3.) {
        std::cout << "Lode angle é Maior que Pi/3. Valido somente para sig1 > sig2 > sig3 -> 0 < theta < Pi/3 " << std::endl;
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




#ifdef PZ_LOG
    {
        std::stringstream sout;
        sout << "\n  TPZTENSOR";
        sout << "\n  LodeAngle = \n" << Lode;
        sout << "\n  dLodeAngle= " << dLode;
        sout << "\n";
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
    output << "XX = " << XX() << " XY = " << XY() << " XZ = " << XZ() << " YY = " << YY() << " YZ = " << YZ() << " ZZ = " << ZZ() << std::endl;
}

#endif //TPZTENSOR_H
