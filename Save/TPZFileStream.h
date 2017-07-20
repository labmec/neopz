#ifndef TPZFILESTREAM_H
#define TPZFILESTREAM_H

#include <fstream>
#include <string>
#include "TPZStream.h"
/**
 * @ingroup save
 * @brief Implements reading from and writing to an ascii file. \ref save
 * "Persistency"
 */
class TPZFileStream : public TPZStream {
 private:
  std::ifstream fIn;
  std::ofstream fOut;

//  template<class T>
//	void ReadData(T *p, int howMany);
//template<class T>
//	void  WriteData(const T *p, int howMany);

	template<class T>
	void ReadData(T *p, int howMany) {
		int c;
		char buf[100];
		if(!fIn)
		{
			DebugStop();
		}
		if(howMany)
		{
			for(c=0; c<howMany; c++) fIn >> p[c];
			fIn.getline(buf,100);
		}
	}
	
	template<class T>
	void  WriteData(const T *p, int howMany)
	{
		for(int c=0; c<howMany; c++) fOut << p[c] << std::endl;
	}
	
 public:
  TPZFileStream();
  virtual ~TPZFileStream();

  void OpenRead(const std::string &fileName);
  void OpenWrite(const std::string &fileName);

  void CloseRead();
  void CloseWrite();

  bool eof() { return fIn.eof(); }
    using TPZStream::Read;
    using TPZStream::Write;
  virtual void Read(std::string *p, int howMany) {
    int c;
    char buf[2560];
    for (c = 0; c < howMany; c++) {
      fIn.getline(buf, 2560);
      p[c] = buf;
    }
  }

  virtual void Write(const std::string *p, int howMany) {
    WriteData<std::string>(p, howMany);
  }
};


//template void TPZFileStream::ReadData<int>(int *p, int howMany);
//template void TPZFileStream::ReadData<unsigned int>(unsigned int *p, int howMany);
//template void TPZFileStream::ReadData<long>(long *p, int howMany);
//template void TPZFileStream::ReadData<float>(float *p, int howMany);
//template void TPZFileStream::ReadData<double>(double *p, int howMany);
//template void TPZFileStream::ReadData<long double>(long double *p, int howMany);
//template void TPZFileStream::ReadData<char>(char *p, int howMany);
//template void TPZFileStream::ReadData<std::complex<float> >(std::complex<float> *p, int howMany);
//template void TPZFileStream::ReadData<std::complex<double> >(std::complex<double> *p, int howMany);
//template void TPZFileStream::ReadData<std::complex<long double> >(std::complex<long double>  *p, int howMany);
//#ifdef _AUTODIFF
//template void TPZFileStream::ReadData<Fad <float> >(Fad <float>  *p, int howMany);
//template void TPZFileStream::ReadData<Fad <double> >(Fad <double>  *p, int howMany);
//template void TPZFileStream::ReadData<Fad<long double> >(Fad<long double>  *p, int howMany);
//#endif//_AUTODIFF
//
//template void TPZFileStream::WriteData<int>(const int *p, int howMany);
//template void TPZFileStream::WriteData<unsigned int>(const unsigned int *p, int howMany);
//template void TPZFileStream::WriteData<long>(const long *p, int howMany);
//template void TPZFileStream::WriteData<float>(const float *p, int howMany);
//template void TPZFileStream::WriteData<double>(const double *p, int howMany);
//template void TPZFileStream::WriteData<long double>(const long double *p, int howMany);
//template void TPZFileStream::WriteData<char>(const char *p, int howMany);
//template void TPZFileStream::WriteData<std::complex<float> >(const std::complex<float> *p, int howMany);
//template void TPZFileStream::WriteData<std::complex<double> >(const std::complex<double> *p, int howMany);
//template void TPZFileStream::WriteData<std::complex<long double> >(const std::complex<long double>  *p, int howMany);
//#ifdef _AUTODIFF
//template void TPZFileStream::WriteData<Fad <float> >(const Fad <float>  *p, int howMany);
//template void TPZFileStream::WriteData<Fad <double> >(const Fad <double>  *p, int howMany);
//template void TPZFileStream::WriteData<Fad<long double> >(const Fad<long double>  *p, int howMany);
//#endif//_AUTODIFF

#endif  // TPZFILESTREAM_H
