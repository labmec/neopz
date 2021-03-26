#include "TPZStream.h"
#include "TPZPersistenceManager.h"

#include "fad.h"

void TPZStream::Write(const bool val) {
    int ival = (val == true) ? 1 : 0;
    Write(&ival);
}

#if defined WIN32
void TPZStream::Write(const long *p, int howMany){ //weird but necessary for working between different OSs
    int64_t *copy = new int64_t[howMany];
    for (unsigned int i = 0; i < howMany; ++i) {
        copy[i] = (int64_t) p[i];
    }
    Write(copy, howMany);
    delete[] copy;
}

void TPZStream::Write(const long unsigned int *p, int howMany){ //weird but necessary for working between different OSs
    uint64_t *copy = new uint64_t[howMany];
    for (unsigned int i = 0; i < howMany; ++i) {
        copy[i] = (uint64_t) p[i];
    }
    Write(copy, howMany);
    delete[] copy;
}
#elif defined __APPLE__

void TPZStream::Write(const long *p, int howMany){ //weird but necessary for working between different OSs
    const int64_t *alias = reinterpret_cast<const int64_t *> (p);
    Write(alias, howMany);
}

void TPZStream::Write(const long unsigned int *p, int howMany){ //weird but necessary for working between different OSs
    const uint64_t *alias = reinterpret_cast<const uint64_t *>(p);
    Write(alias, howMany);
}
#endif

/** @brief Writes howMany floating points at pointer location p */
void TPZStream::Write(const long double *p, int howMany) {//weird but necessary for working between different OSs
    double *copy = new double[howMany];
    for (unsigned int i = 0; i < howMany; ++i) {
        copy[i] = (double) p[i];
    }
    Write(copy, howMany);
    delete[] copy;
}

/** @brief Writes howMany complex-long double at pointer location p */
void TPZStream::Write(const std::complex <long double> *p, int howMany) {//weird but necessary for working between different OSs
    std::complex<double> *copy = new std::complex<double>[howMany];
    for (int i = 0; i < howMany; i++) {
        copy[i] = (std::complex<double>)p[i];
    }
    Write(copy, howMany);
    delete[] copy;
}

void TPZStream::Write(const std::string *p, int howMany) {
    int c;
    for (c = 0; c < howMany; c++) {
        int sz = p[c].size();
        Write(&sz, 1);
        Write(p[c].c_str(), sz);
    }
}

void TPZStream::Write(const TPZFlopCounter *p, int howMany) {
    int i;
    for (i = 0; i < howMany; i++)
        Write(&(p[i].fVal), 1);
}

/** @brief Writes howMany fad-long double at pointer location p */
void TPZStream::Write(const Fad <long double> *p, int howMany) {//weird but necessary for working between different OSs
    Fad<double> *copy = new Fad<double>[howMany];
    for (int i = 0; i < howMany; i++) {
        copy[i] = (Fad<double>)p[i];
    }
    Write(copy, howMany);
    delete[] copy;
}

void TPZStream::Read(bool &val) {
    int ival;
    Read(&ival);
    val = (ival == 0) ? false : true;
}

#if defined WIN32
void TPZStream::Read(long *p, int howMany) { //weird but necessary for working between different OSs
    int64_t *copy = new int64_t[howMany];
    Read(copy, howMany);
    for (unsigned int i = 0; i < howMany; ++i) {
        p[i] = (long) copy[i];
    }
    delete[] copy;
}

void TPZStream::Read(long unsigned int *p, int howMany) { //weird but necessary for working between different OSs
    uint64_t *copy = new uint64_t[howMany];
    Read(copy, howMany);
    for (unsigned int i = 0; i < howMany; ++i) {
        p[i] = (long unsigned int) copy[i];
    }
    delete[] copy;
}

#elif defined __APPLE__

void TPZStream::Read(long *p, int howMany) { //weird but necessary for working between different OSs
    int64_t *alias = reinterpret_cast<int64_t *>(p);
    Read(alias, howMany);
}

void TPZStream::Read(long unsigned int *p, int howMany) { //weird but necessary for working between different OSs
    uint64_t *alias = reinterpret_cast<uint64_t *>(p);
    Read(alias, howMany);
}
#endif

void TPZStream::Read(long double *p, int howMany) {//weird but necessary for working between different OSs
    double *copy = new double[howMany];
    Read(copy, howMany);
    for (unsigned int i = 0; i < howMany; ++i) {
        p[i] = (long double)copy[i];
    }     
    delete[] copy;
}

/** @brief Reads howMany complex-long double from pointer location p */
void TPZStream::Read(std::complex<long double> *p, int howMany) {//weird but necessary for working between different OSs
    std::complex<double> *copy = new std::complex<double>[howMany];
    Read(copy, howMany);
    for (unsigned int i = 0; i < howMany; ++i) {
        p[i] = (std::complex<long double>)copy[i];
    }
    delete[] copy;
}
    
void TPZStream::Read(std::string *p, int howMany) {
    char *temp;
    for (int c = 0; c < howMany; c++) {
        p[c].clear();
        int stringSize = -1;
        Read(&stringSize, 1);
        temp = new char[stringSize+1];
        Read(temp, stringSize);
        temp[stringSize] = 0;
        p[c] = temp;
        delete temp;
    }
}

void TPZStream::Read(TPZFlopCounter *p, int howMany) {
    int i;
    for (i = 0; i < howMany; i++) {
        Read(&(p[i].fVal), 1);
    }
}

/** @brief Reads howMany fad-long double from pointer location p */
void TPZStream::Read(Fad<long double> *p, int howMany) {//weird but necessary for working between different OSs
    Fad<double> *copy = new Fad<double>[howMany];
    Read(copy, howMany);
    for (unsigned int i = 0; i < howMany; ++i) {
        p[i] = (Fad<long double>)copy[i];
    }
    delete[] copy;
}

void TPZStream::Read(Fad<std::complex< float >> *p, int howMany)
{
    std::cout << __PRETTY_FUNCTION__ << " PLEASE IMPLEMENT ME\n";
    DebugStop();
}

void TPZStream::Read(Fad<std::complex< double >> *p, int howMany)
{
    std::cout << __PRETTY_FUNCTION__ << " PLEASE IMPLEMENT ME\n";
    DebugStop();
}

void TPZStream::Read(Fad<std::complex< long double >> *p, int howMany)
{
    std::cout << __PRETTY_FUNCTION__ << " PLEASE IMPLEMENT ME\n";
    DebugStop();
}

void TPZStream::Write(const Fad<std::complex< float > > *p, int howMany)
{
    std::cout << __PRETTY_FUNCTION__ << " PLEASE IMPLEMENT ME\n";
    DebugStop();
}

void TPZStream::Write(const Fad<std::complex< double > >*p, int howMany)
{
    std::cout << __PRETTY_FUNCTION__ << " PLEASE IMPLEMENT ME\n";
    DebugStop();
}

void TPZStream::Write(const Fad<std::complex< long double >> *p, int howMany)
{
    std::cout << __PRETTY_FUNCTION__ << " PLEASE IMPLEMENT ME\n";
    DebugStop();
}


