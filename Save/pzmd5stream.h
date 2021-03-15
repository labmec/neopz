/**
 * @file
 * @brief Contains declaration of the TPZMD5Stream class which implements the interface to write and check md5 files.
 */

#ifndef PZMD5STREAM_H
#define PZMD5STREAM_H

#include "TPZStream.h"

#include <stdio.h>

#ifdef USING_OPENSSL
#include <openssl/md5.h>
#endif

/**
 * @brief Implements the interface to write and check MD5 files. \ref save "Persistency"
 * @ingroup save 
 * @author Edson Borin
 */
class TPZMD5Stream : public TPZStream
{

#ifdef USING_OPENSSL
  /** @brief The MD5 signature. */
  MD5_CTX md5_sig;

  unsigned char digest[MD5_DIGEST_LENGTH];
#endif

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
  TPZMD5Stream() 
  {
    last_status = ResetMD5();
  }

  /** @brief Default destructor */
  virtual ~TPZMD5Stream() 
    {}
  
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
  int CheckMD5(FILE* fh) 
  {
#ifdef USING_OPENSSL
    unsigned char this_digest[MD5_DIGEST_LENGTH];
    unsigned char file_digest[MD5_DIGEST_LENGTH];

    if (last_status != 1)
      return -3;
    
    /* Read file digest. */
    size_t ret = fread(file_digest, sizeof(unsigned char), MD5_DIGEST_LENGTH, fh);
    if (ret != MD5_DIGEST_LENGTH) {
      std::cerr << "fread could not read " << MD5_DIGEST_LENGTH << " items. Read only " << ret << " bytes." << std::endl;
      // Error.
      return 2;
    }

    /* Compute this digest. */
    if (MD5_Final(this_digest, &md5_sig) != 1) {
      // Error.
        return -1;
    }
    
    return compare_digests (this_digest, file_digest, MD5_DIGEST_LENGTH);
#else
    std::cerr << "Enable -DUSING_OPENSSL to use the TPZMD5Stream class." << std::endl; 
#endif
      return 1;
  }

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
  int WriteMD5(FILE* fh) 
  {
#ifdef USING_OPENSSL
    unsigned char digest[MD5_DIGEST_LENGTH];

    if (last_status != 1)
      return -3;
        
    if (MD5_Final(digest, &md5_sig) != 1) {
      return -1;
    }
    
    if (fwrite( (const void*) digest, sizeof(unsigned char), MD5_DIGEST_LENGTH, fh) < MD5_DIGEST_LENGTH)
      return 2;

    return 0; // Return OK
#else
    std::cerr << "Enable -DUSING_OPENSSL to use the TPZMD5Stream class." << std::endl; 
#endif
      return 1;
  }
  
  /**
   * @brief Reset the MD5 signature. 
   * @return Returns 1 if ok, 0 otherwise. 
   */
  int ResetMD5() 
  {
#ifdef USING_OPENSSL
    return MD5_Init(&md5_sig);
#else
    std::cerr << "Enable -DUSING_OPENSSL to use the TPZMD5Stream class." << std::endl; 
#endif
      return 0;
  }
  
  /** @brief Writes size integers at pointer location p */
  virtual void Write(const int *p, int size) {
    Writes<int>(p,size);
  }
  /** @brief Writes size longs at pointer location p */
  virtual void Write(const int64_t *p, int size) {
    Writes<int64_t>(p,size);
  }
  /** @brief Writes size ulongs at pointer location p */
  virtual void Write(const uint64_t *p, int size) {
    Writes<uint64_t>(p,size);
  }
  /** @brief Writes size integers at pointer location p */
  virtual void Write(const unsigned int *p, int size) {
    Writes<unsigned int>(p,size);
  }
  /** @brief Writes size floating points at pointer location p */
  virtual void Write(const float *p, int size) {
    Writes<float>(p,size);
  }
  /** @brief Writes size floating points at pointer location p */
  virtual void Write(const double *p, int size) {
    Writes<double>(p,size);
  }
  /** @brief Writes size floating points at pointer location p */
  virtual void Write(const long double *p, int size) {
    Writes<long double>(p,size);
  }
  /** @brief Writes size chars at pointer location p */
  virtual void Write(const char *p, int size) {
    Writes<char>(p,size);
  }
  /** @brief Writes size uchars at pointer location p */
  virtual void Write(const unsigned char *p, int size) {
    Writes<unsigned char>(p,size);
  }
  /** @brief Writes size strings at pointer location p */
  virtual void Write(const std::string *p, int size) {
    int c;
    for(c=0; c<size; c++) 
      {
        int sz = p[c].size();
        Write(&sz,1);
        Write(p[c].c_str(),p[c].size());
      }
  }
  /** @brief Writes size complex-float at pointer location p */
  virtual void Write(const std::complex <float> *p, int size) {
    Writes< std::complex <float> >(p,size);
  }
  /** @brief Writes size complex-double at pointer location p */
  virtual void Write(const std::complex <double> *p, int size) {
    Writes< std::complex <double> >(p,size);
  }
  /** @brief Writes size complex-long double at pointer location p */
  virtual void Write(const std::complex <long double> *p, int size) {
    Writes< std::complex <long double> >(p,size);
  }
  

	
	virtual void Write(const TFad <1,REAL> *p, int howMany) {
		Writes< TFad <1,REAL> >(p,howMany);
	}
	
	virtual void Write(const TFad <6,REAL> *p, int howMany) {
		Writes< TFad <6,REAL> >(p,howMany);
	}
	
	virtual void Write(const TFad <8,REAL> *p, int howMany) {
		Writes< TFad <8,REAL> >(p,howMany);
	}
	
	virtual void Write(const TFad <9,REAL> *p, int howMany) {
		Writes< TFad <9,REAL> >(p,howMany);
	}
	
	virtual void Write(const TFad <10,REAL> *p, int howMany) {
		Writes< TFad <10,REAL> >(p,howMany);
	}
	
	virtual void Write(const TFad <14,REAL> *p, int howMany) {
		Writes< TFad <14,REAL> >(p,howMany);
	}
	
	virtual void Write(const Fad <float> *p, int howMany) {
		Writes< Fad <float> >(p,howMany);
	}
	
	virtual void Write(const Fad <double> *p, int howMany) {
		Writes< Fad <double> >(p,howMany);
	}
    

  /** @brief Writes size objects of the class T at pointer location p */
  template<class T>
  void  Writes(const T *p, int size) {
#ifdef USING_OPENSSL
    if (last_status == 1)
      last_status = MD5_Update(&md5_sig, (const void*) p, sizeof(T)*size);
#else
    std::cerr << "Enable -DUSING_OPENSSL to use the TPZMD5Stream class." << std::endl; 
#endif
  }

  virtual void Read(TFad <1,REAL> *p, int howMany) {
		ReadError();
	}
	
	virtual void Read(TFad <6,REAL> *p, int howMany) {
		ReadError();
	}
	
	virtual void Read(TFad <8,REAL> *p, int howMany) {
		ReadError();
	}
	
	virtual void Read(TFad <9,REAL> *p, int howMany) {
		ReadError();
	}
	
	virtual void Read(TFad <10,REAL> *p, int howMany) {
		ReadError();
	}
	
	virtual void Read(TFad <14,REAL> *p, int howMany) {
		ReadError();
	}
	
	virtual void Read(Fad <float> *p, int howMany) {
		ReadError();
	}
	
	virtual void Read(Fad <double> *p, int howMany) {
		ReadError();
	}
    
  /** @brief Reads size integers from pointer location p */
  virtual void Read(int *p, int size) {
    ReadError();
  }
  /** @brief Reads size longs from pointer location p */
  virtual void Read(int64_t *p, int size) {
    ReadError();
  }
  /** @brief Reads size longs from pointer location p */
  virtual void Read(uint64_t *p, int size) {
    ReadError();
  }
  /** @brief Reads size integers from pointer location p */
  virtual void Read(unsigned int *p, int size) {
    ReadError();
  }
  /** @brief Reads size floating points from pointer location p */
  virtual void Read(float *p, int size) {
    ReadError();
  }
  /** @brief Reads size floating points from pointer location p */
  virtual void Read(double *p, int size) {
    ReadError();
  }
  /** @brief Reads size floating points from pointer location p */
  virtual void Read(long double *p, int size) {
    ReadError();
  }
  /** @brief Reads size chars from pointer location p */
  virtual void Read(char *p, int size) {
    ReadError();
  }
  /** @brief Reads size chars from pointer location p */
  virtual void Read(unsigned char *p, int size) {
    ReadError();
  }
  /** @brief Reads size strings from pointer location p */
  virtual void Read(std::string *p, int size) 
  {
    ReadError();
  }
  /** @brief Reads size complex-float from pointer location p */
  virtual void Read(std::complex <float> *p, int size) {
    ReadError();
  }
  /** @brief Reads size complex-double from pointer location p */
  virtual void Read(std::complex <double> *p, int size) {
    ReadError();
  }
  /** @brief Reads size complex-long double from pointer location p */
  virtual void Read(std::complex <long double> *p, int size) {
    ReadError();
  }

  void ReadError() 
  {
    std::cerr << "ERROR: Read methods for TPZMD5Stream object invoked." << std::endl;
    return;
  }

};

#endif // PZMD5STREAM_H
