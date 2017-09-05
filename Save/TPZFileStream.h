#ifndef TPZFILESTREAM_H
#define TPZFILESTREAM_H

#include <complex>              // for complex, operator<<, operator>>
#include <sstream>              // for basic_stringbuf<>::int_type, basic_st...
#include "TPZGeneralFStream.h"  // for TPZGeneralFStream
#ifdef _AUTODIFF
template <class T> class Fad;
#endif//_AUTODIFF
/**
 * @ingroup save
 * @brief Implements reading from and writing to an ascii file. \ref save "Persistency"
 */
class TPZFileStream : public TPZGeneralFStream {
  private:
    std::ifstream fIn;
    std::ofstream fOut;

	template<class T>
	void ReadData(T *p, int howMany);
	
    template<class T>
	void WriteData(const T *p, int howMany);

  public:
    TPZFileStream();
    virtual ~TPZFileStream();

    void OpenRead(const std::string &fileName);
    void OpenWrite(const std::string &fileName);

    void CloseRead();
    void CloseWrite();

	bool eof() {
		return fIn.eof();
	}
	
	virtual void Write(const int *p, int howMany) {
		WriteData<int>(p,howMany);
	}
	
	virtual void Write(const unsigned int *p, int howMany) {
		WriteData<unsigned int>(p,howMany);
	}
	
	virtual void Write(const long *p, int howMany) {
		WriteData<long>(p,howMany);
	}
	
	virtual void Write(const float *p, int howMany) {
		WriteData<float>(p,howMany);
	}
	
	virtual void Write(const double *p, int howMany) {
		WriteData<double>(p,howMany);
	}
	
	virtual void Write(const long double *p, int howMany) {
		WriteData<long double>(p,howMany);
	}
	
	virtual void Write(const char *p, int howMany) {
		WriteData<char>(p,howMany);
	}
		
	virtual void Write(const std::complex <float> *p, int howMany) {
		WriteData< std::complex <float> >(p,howMany);
	}
	
	virtual void Write(const std::complex <double> *p, int howMany) {
		WriteData< std::complex <double> >(p,howMany);
	}
	
	virtual void Write(const std::complex <long double> *p, int howMany) {
		WriteData< std::complex <long double> >(p,howMany);
	}
	
#ifdef _AUTODIFF
	
	virtual void Write(const Fad <float> *p, int howMany) {
		WriteData< Fad <float> >(p,howMany);
	}
	
	virtual void Write(const Fad <double> *p, int howMany) {
		WriteData< Fad <double> >(p,howMany);
	}
	
	virtual void Write(const Fad <long double> *p, int howMany) {
		WriteData< Fad <long double> >(p,howMany);
	}
	
#endif
	
	virtual void Read(int *p, int howMany) {
		ReadData<int>(p,howMany);
	}
	
	virtual void Read(unsigned int *p, int howMany) {
		ReadData<unsigned int>(p,howMany);
	}
	
	virtual void Read(long *p, int howMany) {
		ReadData<long>(p,howMany);
	}
	
	virtual void Read(float *p, int howMany) {
		ReadData<float>(p,howMany);
	}
	
	virtual void Read(double *p, int howMany) {
		ReadData<double>(p,howMany);
	}
	
	virtual void Read(long double *p, int howMany) {
		ReadData<long double>(p,howMany);
	}
	
	virtual void Read(char *p, int howMany) {
		ReadData<char>(p,howMany);
	}
	
    virtual void Read(std::complex <float> *p, int howMany) {
		ReadData< std::complex <float> >(p,howMany);
	}
	
	virtual void Read(std::complex <double> *p, int howMany) {
		ReadData< std::complex <double> >(p,howMany);
	}
	
	virtual void Read(std::complex <long double> *p, int howMany) {
		ReadData< std::complex <long double> >(p,howMany);
	}
	
#ifdef _AUTODIFF
	
	virtual void Read(Fad <float> *p, int howMany) {
		ReadData< Fad <float> >(p,howMany);
	}
	
	virtual void Read(Fad <double> *p, int howMany) {
		ReadData< Fad <double> >(p,howMany);
	}
	
	virtual void Read(Fad <long double> *p, int howMany) {
		ReadData< Fad <long double> >(p,howMany);
	}
	
#endif
};
#endif// TPZFILESTREAM_H
