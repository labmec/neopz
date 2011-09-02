/**
 * @file
 * @brief Contains declaration of the TPZBFileStream class which implements the interface to a binary file.
 */
//
// C++ Interface: TPZBFileStream
//
// Description: 
//
//
// Author: Thiago M. N. Oliveira <thiago@labmec.fec.unicamp.br>, (C) 2004
//
// Copyright: See COPYING file that comes with this distribution
//

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
	virtual void Write(REAL *p, int size) {
		Writes<REAL>(p,size);
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
	virtual void Read(REAL *p, int size) {
		Reads<REAL>(p,size);
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
