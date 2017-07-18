/**
 * @file
 * @brief Contains declaration of the abstract TPZStream class. TPZStream defines the interface for saving and reading data,\n
 * TPZFileStream implements reading from and writing to an ascii file.\n
 * TPZBFileStream implements reading from and writing to a binary file.\n
 * Finally, TPZBufferedStream implements reading from and writing to a buffer.
 */
#ifndef TPZSTREAM_H
#define TPZSTREAM_H
#include <string>
#include <complex>
// #include "pzmanvector.h"
#include "pzreal.h"
static unsigned long fCurrentVersion = 1;//TODO:AQUIFRANTake this away
#ifdef _AUTODIFF
#include "fad.h"
#endif


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
	
	virtual void Write(bool val)
	{
		int ival = (val == true) ? 1 : 0;
		Write(&ival);
	}
	
	virtual void Write(const int *p, int howMany=1)=0;
	
	virtual void Write(const unsigned int *p, int howMany=1)=0;
	
	virtual void Write(const long *p, int howMany=1)=0;
	
	virtual void Write(const float  *p, int howMany=1)=0;
	
	virtual void Write(const double  *p, int howMany=1)=0;
	
	virtual void Write(const long double  *p, int howMany=1)=0;
	
	virtual void Write(const char *p, int howMany=1)=0;
	
	virtual void Write(const std::string *p, int howMany=1) = 0;
	
	virtual void Write(const std::complex< float > *p, int howMany=1)=0;
	
	virtual void Write(const std::complex< double > *p, int howMany=1)=0;
	
	virtual void Write(const std::complex< long double > *p, int howMany=1)=0;
	
#ifdef _AUTODIFF
	
	virtual void Write(const Fad< float > *p, int howMany=1)=0;
	
	virtual void Write(const Fad< double > *p, int howMany=1)=0;
	
	virtual void Write(const Fad< long double > *p, int howMany=1)=0;
	
#endif
	
#ifndef ELLIPS
	void Write(TPZFlopCounter *p, int howMany=1)
	{
		int i;
		for(i=0; i<howMany; i++) Write(&(p[i].fVal),1);
	}
#endif
	
	virtual void Read(bool &val)
	{
		int ival;
		Read(&ival);
		val = (ival == 0) ? false : true;
	}
	
	virtual void Read(int *p, int howMany=1)=0;
	
	virtual void Read(unsigned int *p, int howMany=1)=0;
	
	virtual void Read(long *p, int howMany=1)=0;
	
	virtual void Read(float *p, int howMany=1)=0;
	
	virtual void Read(double *p, int howMany=1)=0;
	
	virtual void Read(long double *p, int howMany=1)=0;
	
	virtual void Read(std::complex< float > *p, int howMany=1)=0;
	
	virtual void Read(std::complex< double > *p, int howMany=1)=0;
	
	virtual void Read(std::complex< long double > *p, int howMany=1)=0;
	
#ifdef _AUTODIFF
	
	virtual void Read(Fad< float > *p, int howMany=1)=0;
	
	virtual void Read(Fad< double > *p, int howMany=1)=0;
	
	virtual void Read(Fad< long double > *p, int howMany=1)=0;
	
#endif
	
#ifndef ELLIPS
	void Read(TPZFlopCounter *p, int howMany=1)
	{
		int i;
		for(i=0; i<howMany; i++)
		{
			Read(&(p[i].fVal),1);
		}
	}
#endif
	
	virtual void Read(char *p, int howMany=1)=0;
	
	virtual void Read(std::string *p, int howMany=1) = 0;

protected:
	unsigned long fFromVersion;
	
};
#endif//TPZSTREAM_H
