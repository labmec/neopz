#ifndef TPZBFILESTREAM_H
#define TPZBFILESTREAM_H

#include <fstream>
#include <string>
#include "TPZStream.h"
/**
 * @ingroup save
 * @brief Implements reading from and writing to a binary file. \ref save
 * "Persistency"
 */
class TPZBFileStream : public TPZStream {
private:
	std::ifstream fIn;
	std::ofstream fOut;
	
	/** @brief Reads howMany objects of the class T from pointer location p */
	template <class T> void ReadData(T *p, int howMany);
/** @brief Reads howMany objects of the class T from pointer location p */
	template <class T> void WriteData(const T *p, int howMany);
//	template <class T> void ReadData(T *p, int howMany) {
//		fIn.read(reinterpret_cast<char *>(p), howMany * sizeof(T));
//#ifdef PZDEBUG
//		if (fIn.bad()) {
//			PZError << "TBFileStream:Could not read from stream" << std::endl;
//			DebugStop();
//		}
//#endif
//	}
//	
//	/** @brief Reads howMany objects of the class T from pointer location p */
//	template <class T> void WriteData(const T *p, int howMany) {
//		fOut.write(reinterpret_cast<const char *>(p), howMany * sizeof(T));
//#ifdef PZDEBUG
//		if (fOut.bad()) {
//			PZError << "TBFileStream:Could not write to stream" << std::endl;
//			DebugStop();
//		}
//#endif
//	}
	
public:
	TPZBFileStream();
	~TPZBFileStream();
	
	void OpenRead(const std::string &fileName);
	void OpenWrite(const std::string &fileName);
	
	void CloseRead();
	void CloseWrite();
	
	using TPZStream::Read;
	using TPZStream::Write;
	virtual void Read(std::string *p, int howMany) {
		char buf[1000];
		int c;
		for (c = 0; c < howMany; c++) {
			int sz;
			Read(&sz, 1);
			Read(buf, sz);
			buf[sz] = 0;
			p[c] = buf;
		}
	}
	
	virtual void Write(const std::string *p, int howMany) {
		int c;
		for (c = 0; c < howMany; c++) {
			int sz = p[c].size();
			Write(&sz, 1);
			Write(p[c].c_str(), p[c].size());
		}
	}
};

//template void TPZBFileStream::ReadData<int>(int *p, int howMany);
//template void TPZBFileStream::ReadData<unsigned int>(unsigned int *p, int howMany);
//template void TPZBFileStream::ReadData<long>(long *p, int howMany);
//template void TPZBFileStream::ReadData<float>(float *p, int howMany);
//template void TPZBFileStream::ReadData<double>(double *p, int howMany);
//template void TPZBFileStream::ReadData<long double>(long double *p, int howMany);
//template void TPZBFileStream::ReadData<char>(char *p, int howMany);
//template void TPZBFileStream::ReadData<std::complex<float> >(std::complex<float> *p, int howMany);
//template void TPZBFileStream::ReadData<std::complex<double> >(std::complex<double> *p, int howMany);
//template void TPZBFileStream::ReadData<std::complex<long double> >(std::complex<long double>  *p, int howMany);
//#ifdef _AUTODIFF
//template void TPZBFileStream::ReadData<Fad <float> >(Fad <float>  *p, int howMany);
//template void TPZBFileStream::ReadData<Fad <double> >(Fad <double>  *p, int howMany);
//template void TPZBFileStream::ReadData<Fad<long double> >(Fad<long double>  *p, int howMany);
//#endif//_AUTODIFF
//
//template void TPZBFileStream::WriteData<int>(const int *p, int howMany);
//template void TPZBFileStream::WriteData<unsigned int>(const unsigned int *p, int howMany);
//template void TPZBFileStream::WriteData<long>(const long *p, int howMany);
//template void TPZBFileStream::WriteData<float>(const float *p, int howMany);
//template void TPZBFileStream::WriteData<double>(const double *p, int howMany);
//template void TPZBFileStream::WriteData<long double>(const long double *p, int howMany);
//template void TPZBFileStream::WriteData<char>(const char *p, int howMany);
//template void TPZBFileStream::WriteData<std::complex<float> >(const std::complex<float> *p, int howMany);
//template void TPZBFileStream::WriteData<std::complex<double> >(const std::complex<double> *p, int howMany);
//template void TPZBFileStream::WriteData<std::complex<long double> >(const std::complex<long double>  *p, int howMany);
//#ifdef _AUTODIFF
//template void TPZBFileStream::WriteData<Fad <float> >(const Fad <float>  *p, int howMany);
//template void TPZBFileStream::WriteData<Fad <double> >(const Fad <double>  *p, int howMany);
//template void TPZBFileStream::WriteData<Fad<long double> >(const Fad<long double>  *p, int howMany);
//#endif//_AUTODIFF

#endif  // TPZBFILESTREAM_H
