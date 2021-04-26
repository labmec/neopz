/*
 * @file
 * @brief Contains implementation of TPZMD5Stream class methods.
 */

#include "pzmd5stream.h"

#ifdef USING_OPENSSL
#include <openssl/md5.h>
#else
class MD5state_st{};
#endif

TPZMD5Stream::TPZMD5Stream() 
  {
    last_status = ResetMD5();
  }
TPZMD5Stream::~TPZMD5Stream() 
  {
  }

int TPZMD5Stream::CheckMD5(FILE *fh) {
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
  if (MD5_Final(this_digest, md5_sig.operator->()) != 1) {
    // Error.
    return -1;
  }

  return compare_digests(this_digest, file_digest, MD5_DIGEST_LENGTH);
#else
  std::cerr << "Enable -DUSING_OPENSSL to use the TPZMD5Stream class." << std::endl;
  return 1;
#endif
}

int TPZMD5Stream::WriteMD5(FILE *fh) {
#ifdef USING_OPENSSL
  unsigned char digest[MD5_DIGEST_LENGTH];

  if (last_status != 1)
    return -3;

  if (MD5_Final(digest, md5_sig.operator->()) != 1) {
    return -1;
  }

  if (fwrite((const void *) digest, sizeof(unsigned char), MD5_DIGEST_LENGTH, fh) < MD5_DIGEST_LENGTH)
    return 2;

  return 0; // Return OK
#else
  std::cerr << "Enable -DUSING_OPENSSL to use the TPZMD5Stream class." << std::endl;
  return 1;
#endif
}

int TPZMD5Stream::ResetMD5() {
#ifdef USING_OPENSSL
  return MD5_Init(md5_sig.operator->());
#else
  std::cerr << "Enable -DUSING_OPENSSL to use the TPZMD5Stream class." << std::endl;
  return 0;
#endif
}

template<class T>
void  TPZMD5Stream::Writes(const T *p, int size) {
#ifdef USING_OPENSSL
  if (last_status == 1)
    last_status = MD5_Update(md5_sig.operator->(), (const void*) p, sizeof(T)*size);
#else
  std::cerr << "Enable -DUSING_OPENSSL to use the TPZMD5Stream class." << std::endl;
#endif
}

template void TPZMD5Stream::Writes<int>(const int *p, int size);
template void TPZMD5Stream::Writes<int64_t>(const int64_t* p, int size);
template void TPZMD5Stream::Writes<uint64_t>(const uint64_t* p, int size);
template void TPZMD5Stream::Writes<unsigned int>(const unsigned int* p, int size);
template void TPZMD5Stream::Writes<float>(const float* p, int size);
template void TPZMD5Stream::Writes<double>(const double* p, int size);
template void TPZMD5Stream::Writes<long double>(const long double* p, int size);
template void TPZMD5Stream::Writes<char>(const char* p, int size);
template void TPZMD5Stream::Writes<unsigned char>(const unsigned char* p, int size);

template void TPZMD5Stream::Writes<std::complex<float> >(const std::complex<float> *p, int size);
template void TPZMD5Stream::Writes<std::complex<double> >(const std::complex<double> *p, int size);
template void TPZMD5Stream::Writes<std::complex<long double> >(const std::complex<long double> *p, int size);
template void TPZMD5Stream::Writes<TFad<1, REAL> >(const TFad<1, REAL> *p, int howMany);
template void TPZMD5Stream::Writes<TFad<6, REAL> >(const TFad<6, REAL> *p, int howMany);
template void TPZMD5Stream::Writes<TFad<8, REAL> >(const TFad<8, REAL> *p, int howMany);
template void TPZMD5Stream::Writes<TFad<9, REAL> >(const TFad<9, REAL> *p, int howMany);
template void TPZMD5Stream::Writes<TFad<10, REAL> >(const TFad<10, REAL> *p, int howMany);
template void TPZMD5Stream::Writes<TFad<14, REAL> >(const TFad<14, REAL> *p, int howMany);
template void TPZMD5Stream::Writes<Fad<float> >(const Fad<float> *p, int howMany);
template void TPZMD5Stream::Writes<Fad<double> >(const Fad<double> *p, int howMany);
