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

    void OpenRead(const std::string &fileName) override;
    void OpenWrite(const std::string &fileName) override;
    
    virtual bool AmIOpenForRead() override;
    
    virtual bool AmIOpenForWrite() override;

    void CloseRead() override;
    void CloseWrite() override;

	using TPZStream::Read;
    /** @brief Reads howMany integers from pointer location p */
    virtual void Read(int *p, int howMany) override{ ReadData<int>(p, howMany); }
    /** @brief Reads howMany integers from pointer location p */
    virtual void Read(unsigned int *p, int howMany) override{ ReadData<unsigned int>(p, howMany); }
    /** @brief Reads howMany uint64_t from pointer location p */
    virtual void Read(uint64_t *p, int howMany) override{ ReadData<uint64_t>(p, howMany); }
    /** @brief Reads howMany int64_t from pointer location p */
    virtual void Read(int64_t *p, int howMany)override { ReadData<int64_t>(p, howMany); }
    /** @brief Reads howMany floating points from pointer location p */
    virtual void Read(float *p, int howMany)override { ReadData<float>(p, howMany); }
    /** @brief Reads howMany floating points from pointer location p */
    virtual void Read(double *p, int howMany)override { ReadData<double>(p, howMany); }
    /** @brief Reads howMany unsigned chars from pointer location p */
    virtual void Read(unsigned char *p, int howMany)override { ReadData<unsigned char>(p, howMany); }
    /** @brief Reads howMany chars from pointer location p */
    virtual void Read(char *p, int howMany) override { ReadData<char>(p, howMany); }
    /** @brief Reads howMany complex-float from pointer location p */
    virtual void Read(std::complex<float> *p, int howMany) override { ReadData<std::complex<float>>(p, howMany); }
    /** @brief Reads howMany complex-double from pointer location p */
    virtual void Read(std::complex<double> *p, int howMany) override{ ReadData<std::complex<double>>(p, howMany); }
    virtual void Read(TFad<1,REAL> *p, int howMany) override { ReadData<TFad<1,REAL>>(p, howMany); }
    virtual void Read(TFad<6,REAL> *p, int howMany) override { ReadData<TFad<6,REAL>>(p, howMany); }
    virtual void Read(TFad<8,REAL> *p, int howMany) override { ReadData<TFad<8,REAL>>(p, howMany); }
    virtual void Read(TFad<9,REAL> *p, int howMany) override { ReadData<TFad<9,REAL>>(p, howMany); }
    /** @brief Reads howMany TFad-REAL from pointer location p */
    virtual void Read(TFad<10,REAL> *p, int howMany) override { ReadData<TFad<10,REAL>>(p, howMany); }
    /** @brief Reads howMany TFad-REAL from pointer location p */
    virtual void Read(TFad<14,REAL> *p, int howMany) override { ReadData<TFad<14,REAL>>(p, howMany); }
    /** @brief Reads howMany fad-float from pointer location p */
    virtual void Read(Fad<float> *p, int howMany) override { ReadData<Fad<float>>(p, howMany); }
    /** @brief Reads howMany fad-double from pointer location p */
    virtual void Read(Fad<double> *p, int howMany) override { ReadData<Fad<double>>(p, howMany); }

	using TPZStream::Write;
    /** @brief Writes howMany integers at pointer location p */
    virtual void Write(const int *p, int howMany) override { WriteData<int>(p,howMany); }
    /** @brief Writes howMany integers at pointer location p */
    virtual void Write(const unsigned int *p, int howMany) override { WriteData<const unsigned int>(p,howMany); }
    /** @brief Writes howMany uint64_t at pointer location p */
    virtual void Write(const uint64_t *p, int howMany) override { WriteData<const uint64_t>(p,howMany); }
    /** @brief Writes howMany int64_t at pointer location p */
    virtual void Write(const int64_t *p, int howMany) override { WriteData<const int64_t>(p,howMany); }
    /** @brief Writes howMany floating points at pointer location p */
    virtual void Write(const float *p, int howMany) override { WriteData<float>(p,howMany); }
    /** @brief Writes howMany floating points at pointer location p */
    virtual void Write(const double *p, int howMany) override { WriteData<double>(p,howMany); }
    /** @brief Writes howMany unsigned chars at pointer location p */
    virtual void Write(const unsigned char *p, int howMany) override { WriteData<unsigned char>(p,howMany); }
    /** @brief Writes howMany chars at pointer location p */
    virtual void Write(const char *p, int howMany) override { WriteData<char>(p,howMany); }
    /** @brief Writes howMany complex-float at pointer location p */
    virtual void Write(const std::complex <float> *p, int howMany) override { WriteData< std::complex <float> >(p,howMany); }
    /** @brief Writes howMany complex-double at pointer location p */
    virtual void Write(const std::complex <double> *p, int howMany) override { WriteData< std::complex <double> >(p,howMany); }
    
    /** @brief Writes howMany TFad-REAL at pointer location p */
    virtual void Write(const TFad <1,REAL> *p, int howMany)override { WriteData< TFad <1,REAL> >(p,howMany); }
    /** @brief Writes howMany TFad-REAL at pointer location p */
    virtual void Write(const TFad <6,REAL> *p, int howMany)override { WriteData< TFad <6,REAL> >(p,howMany); }
    /** @brief Writes howMany TFad-REAL at pointer location p */
    virtual void Write(const TFad <8,REAL> *p, int howMany)override { WriteData< TFad <8,REAL> >(p,howMany); }
    /** @brief Writes howMany TFad-REAL at pointer location p */
    virtual void Write(const TFad <9,REAL> *p, int howMany)override { WriteData< TFad <9,REAL> >(p,howMany); }
    /** @brief Writes howMany TFad-REAL at pointer location p */
    virtual void Write(const TFad <10,REAL> *p, int howMany)override { WriteData< TFad <10,REAL> >(p,howMany); }
    /** @brief Writes howMany TFad-REAL at pointer location p */
    virtual void Write(const TFad <14,REAL> *p, int howMany)override { WriteData< TFad <14,REAL> >(p,howMany); }
    /** @brief Writes howMany fad-float at pointer location p */
    virtual void Write(const Fad <float> *p, int howMany)override { WriteData< Fad <float> >(p,howMany); }
    /** @brief Writes howMany fad-double at pointer location p */
    virtual void Write(const Fad <double> *p, int howMany) override { WriteData< Fad <double> >(p,howMany); }
  private:
    template <class T> void ReadData(T *p, int howMany);
    template <class T> void WriteData(const T *p, int howMany);
};
#endif // TPZBFILESTREAM_H
