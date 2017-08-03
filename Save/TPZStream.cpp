#include "TPZStream.h"
#include "TPZSaveable.h"

void TPZStream::Read(std::string &vec) {
    int nel;
    this->Read(&nel, 1);
    TPZManVector<char, 1000> bufstr(nel + 1);
    if (nel)
        this->Read(&bufstr[0], nel);
    bufstr[nel] = '\0';
    vec = &bufstr[0];
}

void TPZStream::Read(TPZVec<std::string> &vec) {
    int nel;
    this->Read(&nel, 1);
    vec.resize(nel);
    for (int i = 0; i < nel; i++) {
        Read(vec[i]);
    }
}
