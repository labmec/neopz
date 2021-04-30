#ifndef TPZFILESTREAM_H
#define TPZFILESTREAM_H

#include "TPZGeneralFStream.h" // for TPZGeneralFStream
#include <complex>             // for complex, operator<<, operator>>
#include <fstream>             //for file streams
#include <sstream>             // for basic_stringbuf<>::int_type, basic_st...
template <class T> class Fad;
template <int Num, class T> class TFad;
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

    void OpenRead(const std::string &fileName) override;
    void OpenWrite(const std::string &fileName) override;
    
    virtual bool AmIOpenForRead() override;
    virtual bool AmIOpenForWrite() override;

    void CloseRead() override;
    void CloseWrite() override;

	bool eof() {
		return fIn.eof();
	}

	using TPZStream::Write;

	virtual void Write(const int *p, int howMany) override {
		WriteData<int>(p,howMany);
	}
	
	virtual void Write(const unsigned int *p, int howMany) override {
		WriteData<unsigned int>(p,howMany);
	}
	
	virtual void Write(const uint64_t *p, int howMany) override {
		WriteData<uint64_t>(p,howMany);
	}
	
	virtual void Write(const int64_t *p, int howMany) override {
		WriteData<int64_t>(p,howMany);
	}
	
	virtual void Write(const float *p, int howMany) override {
		WriteData<float>(p,howMany);
	}
	
	virtual void Write(const double *p, int howMany) override {
		WriteData<double>(p,howMany);
	}
	
	virtual void Write(const unsigned char *p, int howMany) override {
		WriteData<unsigned char>(p,howMany);
	}
	
	virtual void Write(const char *p, int howMany) override {
		WriteData<char>(p,howMany);
	}
		
	virtual void Write(const std::complex <float> *p, int howMany) override {
		WriteData< std::complex <float> >(p,howMany);
	}
	
	virtual void Write(const std::complex <double> *p, int howMany) override {
		WriteData< std::complex <double> >(p,howMany);
	}
	

	virtual void Write(const TFad <1,REAL> *p, int howMany) override {
		WriteData< TFad <1,REAL> >(p,howMany);
	}
	
	virtual void Write(const TFad <6,REAL> *p, int howMany) override {
		WriteData< TFad <6,REAL> >(p,howMany);
	}
	
	virtual void Write(const TFad <8,REAL> *p, int howMany) override {
		WriteData< TFad <8,REAL> >(p,howMany);
	}
	
	virtual void Write(const TFad <9,REAL> *p, int howMany) override {
		WriteData< TFad <9,REAL> >(p,howMany);
	}
	
	virtual void Write(const TFad <10,REAL> *p, int howMany) override {
		WriteData< TFad <10,REAL> >(p,howMany);
	}
	
	virtual void Write(const TFad <14,REAL> *p, int howMany) override {
		WriteData< TFad <14,REAL> >(p,howMany);
	}
	
	virtual void Write(const Fad <float> *p, int howMany) override {
		WriteData< Fad <float> >(p,howMany);
	}
	
	virtual void Write(const Fad <double> *p, int howMany) override {
		WriteData< Fad <double> >(p,howMany);
	}
		

	using TPZStream::Read;

	virtual void Read(int *p, int howMany) override {
		ReadData<int>(p,howMany);
	}
	
	virtual void Read(unsigned int *p, int howMany) override {
		ReadData<unsigned int>(p,howMany);
	}
	
	virtual void Read(uint64_t *p, int howMany) override {
		ReadData<uint64_t>(p,howMany);
	}
	
	virtual void Read(int64_t *p, int howMany) override {
		ReadData<int64_t>(p,howMany);
	}
	
	virtual void Read(float *p, int howMany) override {
		ReadData<float>(p,howMany);
	}
	
	virtual void Read(double *p, int howMany) override {
		ReadData<double>(p,howMany);
	}
	
	virtual void Read(unsigned char *p, int howMany) override {
		ReadData<unsigned char>(p,howMany);
	}
	
	virtual void Read(char *p, int howMany) override {
		ReadData<char>(p,howMany);
	}
	
    virtual void Read(std::complex <float> *p, int howMany) override {
		ReadData< std::complex <float> >(p,howMany);
	}
	
	virtual void Read(std::complex <double> *p, int howMany) override {
		ReadData< std::complex <double> >(p,howMany);
	}
	

	virtual void Read(TFad <1,REAL> *p, int howMany) override {
		ReadData< TFad <1,REAL> >(p,howMany);
	}
	
	virtual void Read(TFad <6,REAL> *p, int howMany) override {
		ReadData< TFad <6,REAL> >(p,howMany);
	}
	
	virtual void Read(TFad <8,REAL> *p, int howMany) override {
		ReadData< TFad <8,REAL> >(p,howMany);
	}
	
	virtual void Read(TFad <9,REAL> *p, int howMany) override {
		ReadData< TFad <9,REAL> >(p,howMany);
	}
	
	virtual void Read(TFad <10,REAL> *p, int howMany) override {
		ReadData< TFad <10,REAL> >(p,howMany);
	}
	
	virtual void Read(TFad <14,REAL> *p, int howMany) override {
		ReadData< TFad <14,REAL> >(p,howMany);
	}
	
	virtual void Read(Fad <float> *p, int howMany) override {
		ReadData< Fad <float> >(p,howMany);
	}
	
	virtual void Read(Fad <double> *p, int howMany) override {
		ReadData< Fad <double> >(p,howMany);
	}
	
	
};
#endif// TPZFILESTREAM_H
