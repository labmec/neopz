#ifndef TPZBFILESTREAM_H
#define TPZBFILESTREAM_H
#include <complex>              // for complex
#include <fstream>              // for string, ifstream, ofstream
#include "TPZGeneralFStream.h"  // for TPZGeneralFStream

#ifndef WIN32
#include <inttypes.h>
typedef int64_t __int64;
typedef uint64_t __uint64;
#else
typedef unsigned __int64 __uint64;
#endif

#ifdef _AUTODIFF
#include "fad.h"
#endif

class TPZBFileStream : public TPZGeneralFStream {
  private:
    std::ifstream fIn;
    std::ofstream fOut;

  public:
    TPZBFileStream();
    ~TPZBFileStream();

    void OpenRead(const std::string &fileName);
    void OpenWrite(const std::string &fileName);
    
    virtual bool AmIOpenForRead();
    
    virtual bool AmIOpenForWrite();

    void CloseRead();
    void CloseWrite();

    /** @brief Reads howMany integers from pointer location p */
    virtual void Read(int *p, int howMany) { ReadData<int>(p, howMany); }
    /** @brief Reads howMany integers from pointer location p */
    virtual void Read(unsigned int *p, int howMany) {
        ReadData<unsigned int>(p, howMany);
    }
    /** @brief Reads howMany long unsigned integers from pointer location p */
    virtual void Read(long unsigned int *p, int howMany) {//weird but necessary for working between different OSs
        __uint64 copy[howMany];
        ReadData<__uint64>(copy, howMany);
        for (int i = 0; i < howMany; i++) {
            p[i] = (long unsigned int)copy[i];
        }
    }
    /** @brief Reads howMany longs from pointer location p */
    virtual void Read(long *p, int howMany) {//weird but necessary for working between different OSs
        __int64 copy[howMany];
        ReadData<__int64>(copy, howMany);
        for (int i = 0; i < howMany; i++) {
            p[i] = (long)copy[i];
        }
    }
    /** @brief Reads howMany floating points from pointer location p */
    virtual void Read(float *p, int howMany) { ReadData<float>(p, howMany); }
    /** @brief Reads howMany floating points from pointer location p */
    virtual void Read(double *p, int howMany) { ReadData<double>(p, howMany); }
    /** @brief Reads howMany floating points from pointer location p */
    virtual void Read(long double *p, int howMany) {//weird but necessary for working between different OSs
        double copy[howMany];
        ReadData<double>(copy, howMany);
        for (int i = 0; i < howMany; i++) {
            p[i] = (long double)copy[i];
        }     
    }
    /** @brief Reads howMany unsigned chars from pointer location p */
    virtual void Read(unsigned char *p, int howMany) { ReadData<unsigned char>(p, howMany); }
    /** @brief Reads howMany chars from pointer location p */
    virtual void Read(char *p, int howMany) { ReadData<char>(p, howMany); }
    
    /** @brief Reads howMany complex-float from pointer location p */
    virtual void Read(std::complex<float> *p, int howMany) {
        ReadData<std::complex<float>>(p, howMany);
    }
    /** @brief Reads howMany complex-double from pointer location p */
    virtual void Read(std::complex<double> *p, int howMany) {
        ReadData<std::complex<double>>(p, howMany);
    }
    /** @brief Reads howMany complex-long double from pointer location p */
    virtual void Read(std::complex<long double> *p, int howMany) {//weird but necessary for working between different OSs
        std::complex<double> *copy = new std::complex<double>(sizeof(std::complex<double>)*howMany);
        ReadData<std::complex<double>>(copy, howMany);
        for (int i = 0; i < howMany; i++) {
            p[i] = (std::complex<long double>)copy[i];
        }  
        delete copy;
    }
#ifdef _AUTODIFF
    /** @brief Reads howMany TFad-REAL from pointer location p */
    virtual void Read(TFad<6,REAL> *p, int howMany) {
        ReadData<TFad<6,REAL>>(p, howMany);
    }
    /** @brief Reads howMany fad-float from pointer location p */
    virtual void Read(Fad<float> *p, int howMany) {
        ReadData<Fad<float>>(p, howMany);
    }
    /** @brief Reads howMany fad-double from pointer location p */
    virtual void Read(Fad<double> *p, int howMany) {
        ReadData<Fad<double>>(p, howMany);
    }
    /** @brief Reads howMany fad-long double from pointer location p */
    virtual void Read(Fad<long double> *p, int howMany) {//weird but necessary for working between different OSs
        Fad<double> *copy = new Fad<double>(sizeof(Fad<double>)*howMany);
        ReadData<Fad<double>>(copy, howMany);
        for (int i = 0; i < howMany; i++) {
            p[i] = (Fad<long double>)copy[i];
        }
        delete copy;
    }
#endif

    /** @brief Writes howMany integers at pointer location p */
    virtual void Write(const int *p, int howMany) {
        WriteData<int>(p,howMany);
    }
    /** @brief Writes howMany integers at pointer location p */
    virtual void Write(const unsigned int *p, int howMany) {
        WriteData<const unsigned int>(p,howMany);
    }
    /** @brief Writes howMany long unsigned integers at pointer location p */
    virtual void Write(const long unsigned int *p, int howMany) {//weird but necessary for working between different OSs
        __uint64 copy[howMany];
        for (int i = 0; i < howMany; i++) {
            copy[i] = (__uint64)p[i];
        }
        WriteData<__uint64>(copy,howMany);
    }
    /** @brief Writes howMany longs at pointer location p */
    virtual void Write(const long *p, int howMany) {//weird but necessary for working between different OSs
        __int64 copy[howMany];
        for (int i = 0; i < howMany; i++) {
            copy[i] = (__int64)p[i];
        }
        WriteData<__int64>(copy,howMany);
    }
    /** @brief Writes howMany floating points at pointer location p */
    virtual void Write(const float *p, int howMany) {
        WriteData<float>(p,howMany);
    }
    /** @brief Writes howMany floating points at pointer location p */
    virtual void Write(const double *p, int howMany) {
        WriteData<double>(p,howMany);
    }
    /** @brief Writes howMany floating points at pointer location p */
    virtual void Write(const long double *p, int howMany) {//weird but necessary for working between different OSs
        double copy[howMany];
        for (int i = 0; i < howMany; i++) {
            copy[i] = (double)p[i];
        }
        WriteData<double>(copy,howMany);
    }
    /** @brief Writes howMany unsigned chars at pointer location p */
    virtual void Write(const unsigned char *p, int howMany) {
        WriteData<unsigned char>(p,howMany);
    }
    /** @brief Writes howMany chars at pointer location p */
    virtual void Write(const char *p, int howMany) {
        WriteData<char>(p,howMany);
    }
    /** @brief Writes howMany complex-float at pointer location p */
    virtual void Write(const std::complex <float> *p, int howMany) {
        WriteData< std::complex <float> >(p,howMany);
    }
    /** @brief Writes howMany complex-double at pointer location p */
    virtual void Write(const std::complex <double> *p, int howMany) {
        WriteData< std::complex <double> >(p,howMany);
    }
    /** @brief Writes howMany complex-long double at pointer location p */
    virtual void Write(const std::complex <long double> *p, int howMany) {//weird but necessary for working between different OSs
        std::complex<double> *copy = new std::complex<double>(sizeof(std::complex<double>)*howMany);
        for (int i = 0; i < howMany; i++) {
            copy[i] = (std::complex<double>)p[i];
        }
        WriteData<std::complex<double>>(copy,howMany);
        delete copy;
    }
    
#ifdef _AUTODIFF
    /** @brief Writes howMany TFad-REAL at pointer location p */
    virtual void Write(const TFad <6,REAL> *p, int howMany) {
        WriteData< TFad <6,REAL> >(p,howMany);
    }
    /** @brief Writes howMany fad-float at pointer location p */
    virtual void Write(const Fad <float> *p, int howMany) {
        WriteData< Fad <float> >(p,howMany);
    }
    /** @brief Writes howMany fad-double at pointer location p */
    virtual void Write(const Fad <double> *p, int howMany) {
        WriteData< Fad <double> >(p,howMany);
    }
    /** @brief Writes howMany fad-long double at pointer location p */
    virtual void Write(const Fad <long double> *p, int howMany) {//weird but necessary for working between different OSs
        Fad<double> *copy = new Fad<double>(sizeof(Fad<double>)*howMany);
        for (int i = 0; i < howMany; i++) {
            copy[i] = (Fad<double>)p[i];
        }
        WriteData<Fad<double>>(copy,howMany);
        delete copy;
    }
#endif
  private:
    template <class T> void ReadData(T *p, int howMany);
    template <class T> void WriteData(const T *p, int howMany);
};
#endif // TPZBFILESTREAM_H
