#ifndef TPZFILESTREAM_H
#define TPZFILESTREAM_H

#include "TPZStream.h"
#include <fstream>
#include <string>
/**
 * @ingroup save
 * @brief Implements reading from and writing to an ascii file. \ref save "Persistency"
 */
class TPZFileStream : public TPZStream {
  private:
    std::ifstream fIn;
    std::ofstream fOut;

	template<class T>
	void ReadData(T *p, int size);
	
    template<class T>
	void WriteData(const T *p, int size);

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
	
	virtual void Write(const std::string *p, int howMany) {
		WriteData<std::string>(p,howMany);
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
	
	virtual void Write(const Fad <float> *p, int size) {
		WriteData< Fad <float> >(p,size);
	}
	
	virtual void Write(const Fad <double> *p, int size) {
		WriteData< Fad <double> >(p,size);
	}
	
	virtual void Write(const Fad <long double> *p, int size) {
		WriteData< Fad <long double> >(p,size);
	}
	
#endif
	
	virtual void Read(int *p, int size) {
		ReadData<int>(p,size);
	}
	
	virtual void Read(unsigned int *p, int size) {
		ReadData<unsigned int>(p,size);
	}
	
	virtual void Read(long *p, int size) {
		ReadData<long>(p,size);
	}
	
	virtual void Read(float *p, int size) {
		ReadData<float>(p,size);
	}
	
	virtual void Read(double *p, int size) {
		ReadData<double>(p,size);
	}
	
	virtual void Read(long double *p, int size) {
		ReadData<long double>(p,size);
	}
	
	virtual void Read(char *p, int size) {
		ReadData<char>(p,size);
	}
	
	virtual void Read(std::string *p, int size) {
		int c;
		char buf[2560];
		for(c=0; c<size; c++)
		{
			fIn.getline(buf,2560);
			p[c] = buf;
		}
	}
	
	virtual void Read(std::complex <float> *p, int size) {
		ReadData< std::complex <float> >(p,size);
	}
	
	virtual void Read(std::complex <double> *p, int size) {
		ReadData< std::complex <double> >(p,size);
	}
	
	virtual void Read(std::complex <long double> *p, int size) {
		ReadData< std::complex <long double> >(p,size);
	}
	
#ifdef _AUTODIFF
	
	virtual void Read(Fad <float> *p, int size) {
		ReadData< Fad <float> >(p,size);
	}
	
	virtual void Read(Fad <double> *p, int size) {
		ReadData< Fad <double> >(p,size);
	}
	
	virtual void Read(Fad <long double> *p, int size) {
		ReadData< Fad <long double> >(p,size);
	}
	
#endif
};
#endif// TPZFILESTREAM_H
