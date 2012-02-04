/**
 * @file
 * @brief Contains declaration of the TPZStream and TPZFileStream classes. TPZStream defines the interface for saving and reading data,\n
 * TPZFileStream implements reading from and writing to an ascii file.
 */
#ifndef PZFILEBUFFERH
#define PZFILEBUFFERH
#include <iostream>
#include <fstream>
#include <string>
#include "pzmanvector.h"
#include "pzreal.h"

/**
 * @ingroup save
 * @brief Defines the interface for saving and reading data. \ref save "Persistency"
 */
/**
 In fact, this class could use the facilities of the stream class of the std library
 This class is a subset of the functionality of the stream classes
 */
class TPZStream {
	
public:
	
	virtual ~TPZStream() {}
	
	virtual void Write(int *p, int size=1)=0;

#ifdef contar
	virtual void Write(double  *p, int size=1)=0;
#else
	virtual void Write(REAL  *p, int size=1)=0;
#endif
	
	virtual void Write(const char *p, int size=1)=0;
	
	virtual void Write(std::string *p, int size=1) = 0;
	
#ifndef ELLIPS
	void Write(TPZFlopCounter *p, int size=1) 
	{
		int i;
		for(i=0; i<size; i++) Write(&(p[i].fVal),1);
	}
#endif
	
	virtual void Read(int *p, int size=1)=0;

#ifdef contar
	virtual void Read(double *p, int size=1)=0;
#else
	virtual void Read(REAL *p, int size=1)=0;
#endif
	
#ifndef ELLIPS
	void Read(TPZFlopCounter *p, int size=1)
	{
		int i;
		for(i=0; i<size; i++)
		{
			Read(&(p[i].fVal),1);
		}
	}
#endif
	
	virtual void Read(char *p, int size=1)=0;
	
	virtual void Read(std::string *p, int size=1) = 0;
	
};

/**
 * @ingroup save
 * @brief Implements reading from and writing to an ascii file. \ref save "Persistency"
 */
class TPZFileStream : public TPZStream {
	
	std::ofstream fo;
	std::ifstream fi;
	
public:
	
	TPZFileStream(){}
	
	virtual ~TPZFileStream();
	
	void OpenWrite(const std::string &filename) {
		fo.open(filename.c_str());
	}
	
	void OpenRead(const std::string &filename) {
		fi.open(filename.c_str());
	}
	
	virtual void Write(int *p, int size) {
		Writes<int>(p,size);
	}
	
	virtual void Write(REAL *p, int size) {
		Writes<REAL>(p,size);
	}
	
	virtual void Write(const char *p, int size) {
		Writes<char>(p,size);
	}
	
	virtual void Write(std::string *p, int size) {
		Writes<std::string>(p,size);
	}
	
	template<class T>
    void  Writes(const T *p, int size) 
	{
		int c;
		for(c=0; c<size; c++) fo << p[c] << std::endl;
	}
	
	
	virtual void Read(int *p, int size) {
		Reads<int>(p,size);
	}
	
	virtual void Read(REAL *p, int size) {
		Reads<REAL>(p,size);
	}
	
	virtual void Read(char *p, int size) {
		Reads<char>(p,size);
	}
	
	virtual void Read(std::string *p, int size) {
		int c;
		char buf[2560];
		for(c=0; c<size; c++) 
		{
			fi.getline(buf,2560);
			p[c] = buf;
		}
	}
	
	template<class T>
    void Reads(T *p, int size) {
		int c;
		char buf[100];
		if(size)
		{
			for(c=0; c<size; c++) fi >> p[c];
			fi.getline(buf,100);
		}
	}
	
};

#endif
