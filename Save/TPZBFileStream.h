#ifndef TPZBFILESTREAM_H
#define TPZBFILESTREAM_H
#include <complex>              // for complex
#include <fstream>              // for string, ifstream, ofstream
#include "TPZGeneralFStream.h"  // for TPZGeneralFStream
#include <stdint.h>              // for uint64_t and int64_t

template <class T> class Fad;
template <int Num, class T> class TFad;

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

	using TPZStream::Read;
    /** @brief Reads howMany integers from pointer location p */
    virtual void Read(int *p, int howMany) { ReadData<int>(p, howMany); }
    /** @brief Reads howMany integers from pointer location p */
    virtual void Read(unsigned int *p, int howMany) { ReadData<unsigned int>(p, howMany); }
    /** @brief Reads howMany uint64_t from pointer location p */
    virtual void Read(uint64_t *p, int howMany) { ReadData<uint64_t>(p, howMany); }
    /** @brief Reads howMany int64_t from pointer location p */
    virtual void Read(int64_t *p, int howMany) { ReadData<int64_t>(p, howMany); }
    /** @brief Reads howMany floating points from pointer location p */
    virtual void Read(float *p, int howMany) { ReadData<float>(p, howMany); }
    /** @brief Reads howMany floating points from pointer location p */
    virtual void Read(double *p, int howMany) { ReadData<double>(p, howMany); }
    /** @brief Reads howMany unsigned chars from pointer location p */
    virtual void Read(unsigned char *p, int howMany) { ReadData<unsigned char>(p, howMany); }
    /** @brief Reads howMany chars from pointer location p */
    virtual void Read(char *p, int howMany) { ReadData<char>(p, howMany); }
    /** @brief Reads howMany complex-float from pointer location p */
    virtual void Read(std::complex<float> *p, int howMany) { ReadData<std::complex<float>>(p, howMany); }
    /** @brief Reads howMany complex-double from pointer location p */
    virtual void Read(std::complex<double> *p, int howMany) { ReadData<std::complex<double>>(p, howMany); }
    virtual void Read(TFad<1,REAL> *p, int howMany) { ReadData<TFad<1,REAL>>(p, howMany); }
    virtual void Read(TFad<6,REAL> *p, int howMany) { ReadData<TFad<6,REAL>>(p, howMany); }
    virtual void Read(TFad<8,REAL> *p, int howMany) { ReadData<TFad<8,REAL>>(p, howMany); }
    virtual void Read(TFad<9,REAL> *p, int howMany) { ReadData<TFad<9,REAL>>(p, howMany); }
    /** @brief Reads howMany TFad-REAL from pointer location p */
    virtual void Read(TFad<10,REAL> *p, int howMany) { ReadData<TFad<10,REAL>>(p, howMany); }
    /** @brief Reads howMany TFad-REAL from pointer location p */
    virtual void Read(TFad<14,REAL> *p, int howMany) { ReadData<TFad<14,REAL>>(p, howMany); }
    /** @brief Reads howMany fad-float from pointer location p */
    virtual void Read(Fad<float> *p, int howMany) { ReadData<Fad<float>>(p, howMany); }
    /** @brief Reads howMany fad-double from pointer location p */
    virtual void Read(Fad<double> *p, int howMany) { ReadData<Fad<double>>(p, howMany); }

	using TPZStream::Write;
    /** @brief Writes howMany integers at pointer location p */
    virtual void Write(const int *p, int howMany) { WriteData<int>(p,howMany); }
    /** @brief Writes howMany integers at pointer location p */
    virtual void Write(const unsigned int *p, int howMany) { WriteData<const unsigned int>(p,howMany); }
    /** @brief Writes howMany uint64_t at pointer location p */
    virtual void Write(const uint64_t *p, int howMany) { WriteData<const uint64_t>(p,howMany); }
    /** @brief Writes howMany int64_t at pointer location p */
    virtual void Write(const int64_t *p, int howMany) { WriteData<const int64_t>(p,howMany); }
    /** @brief Writes howMany floating points at pointer location p */
    virtual void Write(const float *p, int howMany) { WriteData<float>(p,howMany); }
    /** @brief Writes howMany floating points at pointer location p */
    virtual void Write(const double *p, int howMany) { WriteData<double>(p,howMany); }
    /** @brief Writes howMany unsigned chars at pointer location p */
    virtual void Write(const unsigned char *p, int howMany) { WriteData<unsigned char>(p,howMany); }
    /** @brief Writes howMany chars at pointer location p */
    virtual void Write(const char *p, int howMany) { WriteData<char>(p,howMany); }
    /** @brief Writes howMany complex-float at pointer location p */
    virtual void Write(const std::complex <float> *p, int howMany) { WriteData< std::complex <float> >(p,howMany); }
    /** @brief Writes howMany complex-double at pointer location p */
    virtual void Write(const std::complex <double> *p, int howMany) { WriteData< std::complex <double> >(p,howMany); }
    
    /** @brief Writes howMany TFad-REAL at pointer location p */
    virtual void Write(const TFad <1,REAL> *p, int howMany) { WriteData< TFad <1,REAL> >(p,howMany); }
    /** @brief Writes howMany TFad-REAL at pointer location p */
    virtual void Write(const TFad <6,REAL> *p, int howMany) { WriteData< TFad <6,REAL> >(p,howMany); }
    /** @brief Writes howMany TFad-REAL at pointer location p */
    virtual void Write(const TFad <8,REAL> *p, int howMany) { WriteData< TFad <8,REAL> >(p,howMany); }
    /** @brief Writes howMany TFad-REAL at pointer location p */
    virtual void Write(const TFad <9,REAL> *p, int howMany) { WriteData< TFad <9,REAL> >(p,howMany); }
    /** @brief Writes howMany TFad-REAL at pointer location p */
    virtual void Write(const TFad <10,REAL> *p, int howMany) { WriteData< TFad <10,REAL> >(p,howMany); }
    /** @brief Writes howMany TFad-REAL at pointer location p */
    virtual void Write(const TFad <14,REAL> *p, int howMany) { WriteData< TFad <14,REAL> >(p,howMany); }
    /** @brief Writes howMany fad-float at pointer location p */
    virtual void Write(const Fad <float> *p, int howMany) { WriteData< Fad <float> >(p,howMany); }
    /** @brief Writes howMany fad-double at pointer location p */
    virtual void Write(const Fad <double> *p, int howMany) { WriteData< Fad <double> >(p,howMany); }
  private:
    template <class T> void ReadData(T *p, int howMany);
    template <class T> void WriteData(const T *p, int howMany);
};
#endif // TPZBFILESTREAM_H
