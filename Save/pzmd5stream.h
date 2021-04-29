/**
 * @file
 * @brief Contains declaration of the TPZMD5Stream class which implements the interface to write and check md5 files.
 */

#ifndef PZMD5STREAM_H
#define PZMD5STREAM_H

#include "TPZStream.h"
#include "tfad.h"
#include "fad.h"

typedef struct MD5state_st MD5_CTX ;

/**
 * @brief Implements the interface to write and check MD5 files. \ref save "Persistency"
 * @ingroup save 
 * @author Edson Borin
 */
class TPZMD5Stream : public TPZStream
{

  /** @brief The MD5 signature. */
  TPZAutoPointer<MD5_CTX> md5_sig;

  int last_status; // 1 == SUCCESS

  /** @brief Return 0 if digests are equal. -2 otherwise. */ 
  int compare_digests (unsigned char* d1, unsigned char* d2, unsigned int dsz)
  {
    for (unsigned i=0; i< dsz; i++)
    {
      if (d1[i] != d2[i]) 
        return -2; // d1 != d2
    }
    return 0; // no diferences.
  }
	
public:
  /** @brief Default constructor */
  TPZMD5Stream();

  /** @brief Default destructor */
  virtual ~TPZMD5Stream();
  
  /**
   * @brief Check Stream MD5 signature against MD5 signature store on file. 
   * @return Returns 0 if ok;
   *          1 if an error ocurred when opening the file
   *          2 if an error ocurred when reading the MD5 from the file
   *         -1 if an error ocurred when computing this MD5 digest
   *         -2 if the digests are different
   *         -3 error when computing this MD5 signature (last_status != 1)
   */
  int CheckMD5(const std::string &filename) {
    FILE* fd = fopen(filename.c_str(),"rb" );
    if (!fd) return 1;
    int ret = CheckMD5(fd);
    fclose(fd);
    return ret;
  }

  /**
   * @brief Check Stream MD5 signature against MD5 signature store on file.  
   * @return Returns 0 if ok;
   *          2 if an error ocurred when reading the MD5 from the file
   *         -1 if an error ocurred when computing this MD5 digest
   *         -2 if the digests are different
   *         -3 error when computing this MD5 signature (last_status != 1)
   *          1 if is disable OPENSSL
   */
  int CheckMD5(FILE* fh);

  /**
   * @brief Write computed MD5 signature to file.
   * @return Returns 0 if ok;
   *          1 if an error ocurred when opening the file
   *          2 if an error ocurred when writing the MD5 to the file
   *         -1 if an error ocurred when computing this MD5 digest
   *         -3 error when computing this MD5 signature (last_status != 1)
   */

  int WriteMD5(const std::string &filename) {
    FILE* fd = fopen(filename.c_str(), "wb");
    if (fd == NULL) return 1;
    int ret = WriteMD5(fd);
    fclose(fd);
    return ret;
  }
  
  /**
   * @brief Write computed MD5 signature to file. 
   * @return Returns 0 if ok;
   *          2 if an error ocurred when writing the MD5 to the file
   *         -1 if an error ocurred when computing this MD5 digest
   *         -3 error when computing this MD5 signature (last_status != 1)
   *          1 if it is not enable OPENSSL
   */
  int WriteMD5(FILE* fh);

  /**
   * @brief Reset the MD5 signature. 
   * @return Returns 1 if ok, 0 otherwise. 
   */
  int ResetMD5();

  /** @brief Writes size integers at pointer location p */
  virtual void Write(const int *p, int size) override {
    Writes<int>(p,size);
  }
  /** @brief Writes size longs at pointer location p */
  virtual void Write(const int64_t *p, int size) override {
    Writes<int64_t>(p,size);
  }
  /** @brief Writes size ulongs at pointer location p */
  virtual void Write(const uint64_t *p, int size) override {
    Writes<uint64_t>(p,size);
  }
  /** @brief Writes size integers at pointer location p */
  virtual void Write(const unsigned int *p, int size) override {
    Writes<unsigned int>(p,size);
  }
  /** @brief Writes size floating points at pointer location p */
  virtual void Write(const float *p, int size) override {
    Writes<float>(p,size);
  }
  /** @brief Writes size floating points at pointer location p */
  virtual void Write(const double *p, int size) override {
    Writes<double>(p,size);
  }
  /** @brief Writes size floating points at pointer location p */
  virtual void Write(const long double *p, int size) override {
    Writes<long double>(p,size);
  }
  /** @brief Writes size chars at pointer location p */
  virtual void Write(const char *p, int size) override {
    Writes<char>(p,size);
  }
  /** @brief Writes size uchars at pointer location p */
  virtual void Write(const unsigned char *p, int size) override {
    Writes<unsigned char>(p,size);
  }
  /** @brief Writes size strings at pointer location p */
  virtual void Write(const std::string *p, int size) override {
    int c;
    for(c=0; c<size; c++) 
      {
        int sz = p[c].size();
        Write(&sz,1);
        Write(p[c].c_str(),p[c].size());
      }
  }
  /** @brief Writes size complex-float at pointer location p */
  virtual void Write(const std::complex <float> *p, int size) override {
    Writes< std::complex <float> >(p,size);
  }
  /** @brief Writes size complex-double at pointer location p */
  virtual void Write(const std::complex <double> *p, int size) override {
    Writes< std::complex <double> >(p,size);
  }
  /** @brief Writes size complex-long double at pointer location p */
  virtual void Write(const std::complex <long double> *p, int size) override {
    Writes< std::complex <long double> >(p,size);
  }
  

	
	virtual void Write(const TFad <1,REAL> *p, int howMany) override {
		Writes< TFad <1,REAL> >(p,howMany);
	}
	
	virtual void Write(const TFad <6,REAL> *p, int howMany) override {
		Writes< TFad <6,REAL> >(p,howMany);
	}
	
	virtual void Write(const TFad <8,REAL> *p, int howMany) override {
		Writes< TFad <8,REAL> >(p,howMany);
	}
	
	virtual void Write(const TFad <9,REAL> *p, int howMany) override {
		Writes< TFad <9,REAL> >(p,howMany);
	}
	
	virtual void Write(const TFad <10,REAL> *p, int howMany) override {
		Writes< TFad <10,REAL> >(p,howMany);
	}
	
	virtual void Write(const TFad <14,REAL> *p, int howMany) override {
		Writes< TFad <14,REAL> >(p,howMany);
	}
	
	virtual void Write(const Fad <float> *p, int howMany) override {
		Writes< Fad <float> >(p,howMany);
	}
	
	virtual void Write(const Fad <double> *p, int howMany) override {
		Writes< Fad <double> >(p,howMany);
	}
    
  /** @brief Writes size objects of the class T at pointer location p */
  template<class T>
  void  Writes(const T *p, int size);

  virtual void Read(TFad <1,REAL> *p, int howMany) override {
		ReadError();
	}
	
	virtual void Read(TFad <6,REAL> *p, int howMany) override {
		ReadError();
	}
	
	virtual void Read(TFad <8,REAL> *p, int howMany) override {
		ReadError();
	}
	
	virtual void Read(TFad <9,REAL> *p, int howMany) override {
		ReadError();
	}
	
	virtual void Read(TFad <10,REAL> *p, int howMany) override {
		ReadError();
	}
	
	virtual void Read(TFad <14,REAL> *p, int howMany) override {
		ReadError();
	}
	
	virtual void Read(Fad <float> *p, int howMany) override {
		ReadError();
	}
	
	virtual void Read(Fad <double> *p, int howMany) override {
		ReadError();
	}
    
  /** @brief Reads size integers from pointer location p */
  virtual void Read(int *p, int size) override {
    ReadError();
  }
  /** @brief Reads size longs from pointer location p */
  virtual void Read(int64_t *p, int size) override {
    ReadError();
  }
  /** @brief Reads size longs from pointer location p */
  virtual void Read(uint64_t *p, int size) override {
    ReadError();
  }
  /** @brief Reads size integers from pointer location p */
  virtual void Read(unsigned int *p, int size) override {
    ReadError();
  }
  /** @brief Reads size floating points from pointer location p */
  virtual void Read(float *p, int size) override {
    ReadError();
  }
  /** @brief Reads size floating points from pointer location p */
  virtual void Read(double *p, int size) override {
    ReadError();
  }
  /** @brief Reads size floating points from pointer location p */
  virtual void Read(long double *p, int size) override {
    ReadError();
  }
  /** @brief Reads size chars from pointer location p */
  virtual void Read(char *p, int size) override {
    ReadError();
  }
  /** @brief Reads size chars from pointer location p */
  virtual void Read(unsigned char *p, int size) override {
    ReadError();
  }
  /** @brief Reads size strings from pointer location p */
  virtual void Read(std::string *p, int size) override
  {
    ReadError();
  }
  /** @brief Reads size complex-float from pointer location p */
  virtual void Read(std::complex <float> *p, int size) override {
    ReadError();
  }
  /** @brief Reads size complex-double from pointer location p */
  virtual void Read(std::complex <double> *p, int size) override {
    ReadError();
  }
  /** @brief Reads size complex-long double from pointer location p */
  virtual void Read(std::complex <long double> *p, int size) override {
    ReadError();
  }

  void ReadError() 
  {
    std::cerr << "ERROR: Read methods for TPZMD5Stream object invoked." << std::endl;
    return;
  }

};

#endif // PZMD5STREAM_H
