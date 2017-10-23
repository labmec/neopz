#ifndef TPZBFILESTREAM_H
#define TPZBFILESTREAM_H
#include <complex>              // for complex
#include <fstream>              // for string, ifstream, ofstream
#include "TPZGeneralFStream.h"  // for TPZGeneralFStream
#ifdef _AUTODIFF
template <class T> class Fad;
template <int Num, class T> class TFad;
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
    virtual void Read(long unsigned int *p, int howMany) {
        ReadData<long unsigned int>(p, howMany);
    }
    /** @brief Reads howMany longs from pointer location p */
    virtual void Read(long *p, int howMany) { ReadData<long>(p, howMany); }
    /** @brief Reads howMany floating points from pointer location p */
    virtual void Read(float *p, int howMany) { ReadData<float>(p, howMany); }
    /** @brief Reads howMany floating points from pointer location p */
    virtual void Read(double *p, int howMany) { ReadData<double>(p, howMany); }
    /** @brief Reads howMany floating points from pointer location p */
    virtual void Read(long double *p, int howMany) {
        ReadData<long double>(p, howMany);
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
    virtual void Read(std::complex<long double> *p, int howMany) {
        ReadData<std::complex<long double>>(p, howMany);
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
    virtual void Read(Fad<long double> *p, int howMany) {
        ReadData<Fad<long double>>(p, howMany);
    }
#endif

    /** @brief Writes howMany integers at pointer location p */
    virtual void Write(const int *p, int howMany) {
        WriteData<int>(p,howMany);
    }
    /** @brief Writes howMany integers at pointer location p */
    virtual void Write(const unsigned int *p, int howMany) {
        WriteData<unsigned int>(p,howMany);
    }
    /** @brief Writes howMany long unsigned integers at pointer location p */
    virtual void Write(const long unsigned int *p, int howMany) {
        WriteData<long unsigned int>(p,howMany);
    }
    /** @brief Writes howMany longs at pointer location p */
    virtual void Write(const long *p, int howMany) {
        WriteData<long>(p,howMany);
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
    virtual void Write(const long double *p, int howMany) {
        WriteData<long double>(p,howMany);
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
    virtual void Write(const std::complex <long double> *p, int howMany) {
        WriteData< std::complex <long double> >(p,howMany);
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
    virtual void Write(const Fad <long double> *p, int howMany) {
        WriteData< Fad <long double> >(p,howMany);
    }
#endif
  private:
    template <class T> void ReadData(T *p, int howMany);
    template <class T> void WriteData(const T *p, int howMany);
};
#endif // TPZBFILESTREAM_H
