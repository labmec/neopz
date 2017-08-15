#include "TPZStream.h"

void TPZStream::Write(const bool val) {
    int ival = (val == true) ? 1 : 0;
    Write(&ival);
}

void TPZStream::Write(const std::string *p, int howMany) {
    int c;
    for (c = 0; c < howMany; c++) {
        int sz = p[c].size();
        Write(&sz, 1);
        Write(p[c].c_str(), p[c].size());
    }
}

#ifndef ELLIPS
void TPZStream::Write(const TPZFlopCounter *p, int howMany) {
    int i;
    for (i = 0; i < howMany; i++)
        Write(&(p[i].fVal), 1);
}
#endif

void TPZStream::Read(bool &val) {
    int ival;
    Read(&ival);
    val = (ival == 0) ? false : true;
}

void TPZStream::Read(std::string *p, int howMany) {
    char *temp;
    for (int c = 0; c < howMany; c++) {
        p[c].clear();
        int stringSize = -1;
        Read(&stringSize, 1);
        temp = new char[stringSize];
        Read(temp, stringSize);
        p[c] = temp;
        delete temp;
    }
}

#ifndef ELLIPS
void TPZStream::Read(TPZFlopCounter *p, int howMany) {
    int i;
    for (i = 0; i < howMany; i++) {
        Read(&(p[i].fVal), 1);
    }
}
#endif

