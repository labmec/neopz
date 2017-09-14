#ifndef TPZCONTBUFFEREDSTREAM_H
#define TPZCONTBUFFEREDSTREAM_H
#include "TPZStream.h"
#include <sstream>

/**
 * @brief      Class for creating a bidirectional circular buffer
 */
class TPZContBufferedStream : public TPZStream {
  public:
    /**
     * @brief      Creates a buffer
     *
     */
    TPZContBufferedStream();

    /**
     * @brief      Copy constructor.
     *
     * @param[in]  readBuffer  The buffer to be copied.
     */
    TPZContBufferedStream(const TPZContBufferedStream &readBuffer);

    /**
     * @brief      Destroys the object.
     */
    ~TPZContBufferedStream();

    /**
     * @brief      Assingment operator. Both buffers will have the same
     * underlying stream, so this must be used with care.
     *
     * @param[in]  other    The other buffer
     * @return     Reference to this buffer
     */
    TPZContBufferedStream &operator=(const TPZContBufferedStream &other);

    /**
     * @brief      It reads all data in a buffer, consuming it
     *
     * @param      other  The other buffer
     *
     * @return     Reference to this buffer
     */
    virtual TPZContBufferedStream &operator<<(TPZContBufferedStream &other);

    /**
     * @brief      It reads all data in a buffer WITHOUT consuming it
     *
     * @param[in]  other  The other buffer
     *
     * @return     Reference to this buffer
     */
    virtual TPZContBufferedStream &operator<<(const TPZContBufferedStream &other);

    /**
     * @brief      Prints buffer info and data
     */
    void Print();    
    
    /**
     * @brief  Get all buffer data to a char* in a contiguous manner.
     * May be interesting if one needs to this data to be serialized(i.e., for MPI)
     *
     * @param[in]  dest Array where data will be stored.
     */
    void GetDataFromBuffer(char *dest) const;
    
    void clear();

    virtual void Write(const int *p, int howMany) {
        WriteData<int>(p, howMany);
    }

    virtual void Write(const unsigned int *p, int howMany) {
        WriteData<unsigned int>(p, howMany);
    }

    virtual void Write(const long unsigned int *p, int howMany) {
        WriteData<long unsigned int>(p, howMany);
    }

    virtual void Write(const long *p, int howMany) {
        WriteData<long>(p, howMany);
    }

    virtual void Write(const float *p, int howMany) {
        WriteData<float>(p, howMany);
    }

    virtual void Write(const double *p, int howMany) {
        WriteData<double>(p, howMany);
    }

    virtual void Write(const long double *p, int howMany) {
        WriteData<long double>(p, howMany);
    }

    virtual void Write(const char *p, int howMany) {
        WriteData<char>(p, howMany);
    }

    virtual void Write(const std::complex<float> *p, int howMany) {
        WriteData<std::complex<float>>(p, howMany);
    }

    virtual void Write(const std::complex<double> *p, int howMany) {
        WriteData<std::complex<double>>(p, howMany);
    }

    virtual void Write(const std::complex<long double> *p, int howMany) {
        WriteData<std::complex<long double>>(p, howMany);
    }

#ifdef _AUTODIFF

    virtual void Write(const Fad<float> *p, int howMany) {
        WriteData<Fad<float>>(p, howMany);
    }

    virtual void Write(const Fad<double> *p, int howMany) {
        WriteData<Fad<double>>(p, howMany);
    }

    virtual void Write(const Fad<long double> *p, int howMany) {
        WriteData<Fad<long double>>(p, howMany);
    }

#endif

    virtual void Read(int *p, int howMany) { ReadData<int>(p, howMany); }

    virtual void Read(unsigned int *p, int howMany) {
        ReadData<unsigned int>(p, howMany);
    }

    virtual void Read(long unsigned int *p, int howMany) {
        ReadData<long unsigned int>(p, howMany);
    }

    virtual void Read(long *p, int howMany) { ReadData<long>(p, howMany); }

    virtual void Read(float *p, int howMany) { ReadData<float>(p, howMany); }

    virtual void Read(double *p, int howMany) { ReadData<double>(p, howMany); }

    virtual void Read(long double *p, int howMany) {
        ReadData<long double>(p, howMany);
    }

    virtual void Read(char *p, int howMany) { ReadData<char>(p, howMany); }

    virtual void Read(std::complex<float> *p, int howMany) {
        ReadData<std::complex<float>>(p, howMany);
    }

    virtual void Read(std::complex<double> *p, int howMany) {
        ReadData<std::complex<double>>(p, howMany);
    }

    virtual void Read(std::complex<long double> *p, int howMany) {
        ReadData<std::complex<long double>>(p, howMany);
    }

#ifdef _AUTODIFF

    virtual void Read(Fad<float> *p, int howMany) {
        ReadData<Fad<float>>(p, howMany);
    }

    virtual void Read(Fad<double> *p, int howMany) {
        ReadData<Fad<double>>(p, howMany);
    }

    virtual void Read(Fad<long double> *p, int howMany) {
        ReadData<Fad<long double>>(p, howMany);
    }
    
    size_t Size() const;

#endif

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
#endif // TPZCONTBUFFEREDSTREAM_H
