#ifndef TPZCIRCBUFFEREDSTREAM_H
#define TPZCIRCBUFFEREDSTREAM_H
#include "TPZStream.h"
#include <sstream>

#include "fad.h"
#include "tfad.h"

/**
 * @brief      Class for creating a bidirectional circular buffer
 */
class TPZCircBufferedStream : public TPZStream {
  public:
    /**
     * @brief      Creates a bidirectional buffer
     */
    TPZCircBufferedStream();

    /**
     * @brief      Copy constructor. Both buffers will have the same underlying
     * stream, so this must be used with care.
     *
     * @param[in]  readBuffer  The buffer to be copied.
     */
    TPZCircBufferedStream(const TPZStream &readBuffer);

    /**
     * @brief      Destroys the object.
     */
    ~TPZCircBufferedStream();

    /**
     * @brief      Assingment operator. Both buffers will have the same
     * underlying stream, so this must be used with care.
     *
     * @param[in]  other    The other buffer
     * @return     Reference to this buffer
     */
    TPZCircBufferedStream &operator=(const TPZCircBufferedStream &other);

    /**
     * @brief      It reads all data in a buffer, consuming it
     *
     * @param      other  The other buffer
     *
     * @return     Reference to this buffer
     */
    virtual TPZCircBufferedStream &operator<<(TPZCircBufferedStream &other);

    /**
     * @brief      It reads all data in a buffer WITHOUT consuming it
     *
     * @param[in]  other  The other buffer
     *
     * @return     Reference to this buffer
     */
    virtual TPZCircBufferedStream &operator<<(const TPZCircBufferedStream &other);

    /**
     * @brief      Prints buffer info and data
     */
    void Print();    
    
    /**
     * @brief      Get all buffer data to a char* in a contiguous manner.
     * May be interesting if one needs to this data to be serialized(i.e., for MPI)
     *
     * @param[in]  dest Array where data will be stored.
     */
    void GetDataFromBuffer(char *dest) const{
        ConstReadFromBuffer(dest, fSize);
    }

	using TPZStream::Write;
    virtual void Write(const int *p, int howMany) override;
    virtual void Write(const unsigned int *p, int howMany) override;
    virtual void Write(const uint64_t *p, int howMany) override;
    virtual void Write(const int64_t *p, int howMany) override;
    virtual void Write(const float *p, int howMany) override;
    virtual void Write(const double *p, int howMany) override;
    virtual void Write(const unsigned char *p, int howMany) override;
    virtual void Write(const char *p, int howMany) override;
    virtual void Write(const std::complex<float> *p, int howMany) override;
    virtual void Write(const std::complex<double> *p, int howMany) override;

    virtual void Write(const TFad<1,REAL> *p, int howMany) override;
    virtual void Write(const TFad<6,REAL> *p, int howMany) override;
    virtual void Write(const TFad<8,REAL> *p, int howMany) override;
    virtual void Write(const TFad<9,REAL> *p, int howMany) override;
    virtual void Write(const TFad<10,REAL> *p, int howMany) override;
    virtual void Write(const TFad<14,REAL> *p, int howMany) override;
    virtual void Write(const Fad<float> *p, int howMany) override;
    virtual void Write(const Fad<double> *p, int howMany) override;

	using TPZStream::Read;
    virtual void Read(int *p, int howMany) override;
    virtual void Read(unsigned int *p, int howMany) override;
    virtual void Read(uint64_t *p, int howMany) override;
    virtual void Read(int64_t *p, int howMany) override;
    virtual void Read(float *p, int howMany) override;
    virtual void Read(double *p, int howMany) override;
    virtual void Read(long double *p, int howMany) override;
    virtual void Read(unsigned char *p, int howMany) override;
    virtual void Read(char *p, int howMany) override;
    virtual void Read(std::complex<float> *p, int howMany) override;
    virtual void Read(std::complex<double> *p, int howMany) override;

    virtual void Read(TFad<1,REAL> *p, int howMany) override;
    virtual void Read(TFad<6,REAL> *p, int howMany) override;
    virtual void Read(TFad<8,REAL> *p, int howMany) override;
    virtual void Read(TFad<9,REAL> *p, int howMany) override;
    virtual void Read(TFad<10,REAL> *p, int howMany) override;
    virtual void Read(TFad<14,REAL> *p, int howMany) override;
    virtual void Read(Fad<float> *p, int howMany) override;
    virtual void Read(Fad<double> *p, int howMany) override;

  protected:
    /**
     * @brief      Reads from buffer.
     *
     * @param      dest    The destination
     * @param[in]  nBytes  How many bytes will be read
     */
    virtual void ReadFromBuffer(char *dest, const size_t &nBytes);

    /**
     * @brief      Reads from buffer WITHOUT consuming it. Unless it is still
     * reading from its underlying stream. In this situation there is no
     * guarantee.
     *
     * @param      dest    The destination
     * @param[in]  nBytes  How many bytes will be read
     */
    virtual void ConstRead(char *dest, const size_t &nBytes) const;

    /**
     * @brief      Writes to buffer.
     *
     * @param[in]  source  The source
     * @param[in]  nBytes  How many bytes will be written
     */
    virtual void WriteToBuffer(const char *source, const size_t &nBytes);

  private:
    char *fBuffer;

    char *fFirst;

    char *fLast;

    size_t fNAllocatedBytes;

    size_t fSize;

    /**
     * @brief      Reads from buffer WITHOUT consuming it.
     *
     * @param      dest    The destination
     * @param[in]  nBytes  How many bytes will be read
     */
    virtual void ConstReadFromBuffer(char *dest, const size_t &nBytes) const;

    template <class T> void ReadData(T *p, int howMany);

    template <class T> void WriteData(const T *p, int howMany);

    static const size_t MIN_SIZE_INCREMENT = size_t(1);
};

template <typename T>
void TPZCircBufferedStream::ReadData(T *p, int howMany) {
    ReadFromBuffer(reinterpret_cast<char *> (p), howMany * sizeof (T));
}

template <typename T>
void TPZCircBufferedStream::WriteData(const T *p, int howMany) {
    WriteToBuffer(reinterpret_cast<const char *> (p), howMany * sizeof (T));
}

template void TPZCircBufferedStream::WriteData(const int *, int howMany);
template void TPZCircBufferedStream::WriteData(const unsigned int* p, int howMany);
template void TPZCircBufferedStream::WriteData(const uint64_t *p, int howMany);
template void TPZCircBufferedStream::WriteData(const int64_t *p, int howMany);
template void TPZCircBufferedStream::WriteData(const float *p, int howMany);
template void TPZCircBufferedStream::WriteData(const double *p, int howMany);
template void TPZCircBufferedStream::WriteData(const char *p, int howMany);
template void TPZCircBufferedStream::WriteData(const std::complex<float> *p, int howMany);
template void TPZCircBufferedStream::WriteData(const std::complex<double> *p, int howMany);

template void TPZCircBufferedStream::WriteData(const TFad<1,REAL> *p, int howMany);
template void TPZCircBufferedStream::WriteData(const TFad<6,REAL> *p, int howMany);
template void TPZCircBufferedStream::WriteData(const TFad<8,REAL> *p, int howMany);
template void TPZCircBufferedStream::WriteData(const TFad<9,REAL> *p, int howMany);
template void TPZCircBufferedStream::WriteData(const TFad<10,REAL> *p, int howMany);
template void TPZCircBufferedStream::WriteData(const TFad<14,REAL> *p, int howMany);
template void TPZCircBufferedStream::WriteData(const Fad<float> *p, int howMany);
template void TPZCircBufferedStream::WriteData(const Fad<double> *p, int howMany);

template void TPZCircBufferedStream::ReadData(int *p, int howMany);
template void TPZCircBufferedStream::ReadData(unsigned int *p, int howMany);
template void TPZCircBufferedStream::ReadData(uint64_t *p, int howMany);
template void TPZCircBufferedStream::ReadData(int64_t *p, int howMany);
template void TPZCircBufferedStream::ReadData(float *p, int howMany);
template void TPZCircBufferedStream::ReadData(double *p, int howMany);
template void TPZCircBufferedStream::ReadData(long double *p, int howMany);
template void TPZCircBufferedStream::ReadData(char *p, int howMany);
template void TPZCircBufferedStream::ReadData(std::complex<float> *p, int howMany);
template void TPZCircBufferedStream::ReadData(std::complex<double> *p, int howMany);

template void TPZCircBufferedStream::ReadData(TFad<1,REAL> *p, int howMany);
template void TPZCircBufferedStream::ReadData(TFad<6,REAL> *p, int howMany);
template void TPZCircBufferedStream::ReadData(TFad<8,REAL> *p, int howMany);
template void TPZCircBufferedStream::ReadData(TFad<9,REAL> *p, int howMany);
template void TPZCircBufferedStream::ReadData(TFad<10,REAL> *p, int howMany);
template void TPZCircBufferedStream::ReadData(TFad<14,REAL> *p, int howMany);
template void TPZCircBufferedStream::ReadData(Fad<float> *p, int howMany);
template void TPZCircBufferedStream::ReadData(Fad<double> *p, int howMany);

#endif // TPZCIRCBUFFEREDSTREAM_H
