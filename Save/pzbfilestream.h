/**
 * @file
 * @brief Contains declaration of the TPZBFileStream class which implements the interface to a binary file.
 */

#ifndef STDPZBFILESTREAM_H
#define STDPZBFILESTREAM_H

#include "pzfilebuffer.h"

#include <stdio.h>

/**
 * @brief Implements the interface to a binary file. \ref save "Persistency"
 * @ingroup save 
 * @author Thiago M. N. Oliveira
 */
class TPZBFileStream : public TPZStream
{
	/** @brief Output file */
	FILE *ofd;
	/** @brief Input file */
	FILE *ifd;
	
public:
	/** @brief Default constructor */
	TPZBFileStream(){
		ofd=0;
		ifd=0;
	}
	/** @brief Default destructor */
	virtual ~TPZBFileStream() {
		if(ofd) fclose(ofd);
		if(ifd) fclose(ifd);
	}
	/** @brief Open file to write */
	void OpenWrite(const std::string &filename) {
		ofd = fopen(filename.c_str(),"wb" );
	}
	/** @brief Open file to read */
	void OpenRead(const std::string &filename) {
		ifd = fopen(filename.c_str(), "rb");
		if(!ifd)
		{
			std::cout << "could not open file " << filename << std::endl;
		}
	}
	/** @brief Writes size integers at pointer location p */
	virtual void Write(int *p, int size) {
		Writes<int>(p,size);
	}
	/** @brief Writes size floating points at pointer location p */
	virtual void Write(float *p, int size) {
		Writes<float>(p,size);
	}
	/** @brief Writes size floating points at pointer location p */
	virtual void Write(double *p, int size) {
		Writes<double>(p,size);
	}
	/** @brief Writes size floating points at pointer location p */
	virtual void Write(long double *p, int size) {
		Writes<long double>(p,size);
	}
	/** @brief Writes size chars at pointer location p */
	virtual void Write(const char *p, int size) {
		Writes<char>(p,size);
	}
	/** @brief Writes size strings at pointer location p */
	virtual void Write(std::string *p, int size) {
		int c;
		for(c=0; c<size; c++) 
		{
			int sz = p[c].size();
			Write(&sz,1);
			Write(p[c].c_str(),p[c].size());
		}
	}
	/** @brief Writes size complex-float at pointer location p */
	virtual void Write(std::complex <float> *p, int size) {
		Writes< std::complex <float> >(p,size);
	}
	/** @brief Writes size complex-double at pointer location p */
	virtual void Write(std::complex <double> *p, int size) {
		Writes< std::complex <double> >(p,size);
	}
	/** @brief Writes size complex-long double at pointer location p */
	virtual void Write(std::complex <long double> *p, int size) {
		Writes< std::complex <long double> >(p,size);
	}
	/** @brief Writes size objects of the class T at pointer location p */
	template<class T>
	void  Writes(const T *p, int size) 
	{
		fwrite(p,sizeof(T),size,ofd);
	}
	/** @brief Reads size integers from pointer location p */
	virtual void Read(int *p, int size) {
		Reads<int>(p,size);
	}
	/** @brief Reads size floating points from pointer location p */
	virtual void Read(float *p, int size) {
		Reads<float>(p,size);
	}
	/** @brief Reads size floating points from pointer location p */
	virtual void Read(double *p, int size) {
		Reads<double>(p,size);
	}
	/** @brief Reads size floating points from pointer location p */
	virtual void Read(long double *p, int size) {
		Reads<long double>(p,size);
	}
	/** @brief Reads size chars from pointer location p */
	virtual void Read(char *p, int size) {
		Reads<char>(p,size);
	}
	/** @brief Reads size strings from pointer location p */
	virtual void Read(std::string *p, int size) 
	{
		char buf[1000];
		int c;
		for(c=0; c<size; c++) 
		{
			int sz;
			Read(&sz,1);
			Read(buf,sz);
			buf[sz] = 0;
			p[c] = buf;
		}
	}
	/** @brief Reads size complex-float from pointer location p */
	virtual void Read(std::complex <float> *p, int size) {
		Reads< std::complex <float> >(p,size);
	}
	/** @brief Reads size complex-double from pointer location p */
	virtual void Read(std::complex <double> *p, int size) {
		Reads< std::complex <double> >(p,size);
	}
	/** @brief Reads size complex-long double from pointer location p */
	virtual void Read(std::complex <long double> *p, int size) {
		Reads< std::complex <long double> >(p,size);
	}
	/** @brief Reads size objects of the class T from pointer location p */
	template<class T>
	void Reads(T *p, int size)
	{
		if(ifd)
		{
			fread(p,sizeof(T),size,ifd);
		}
	}
	
};

#endif
