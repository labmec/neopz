#include "pzfilebuffer.h"

template <class T> void TPZStream::Write(const TPZVec<T> &vec) {
        long c, nc = vec.NElements();
        this->Write(&nc);
        for (c = 0; c < nc; c++)
            vec[c].Write(*this, 0);
    }

    template <class T> void TPZStream::Write(const std::vector<T> &vec) {
        int c, nc = vec.size();
        this->Write(&nc);
        for (c = 0; c < nc; c++)
            vec[c].Write(*this, 0);
    }

    template <class T, int EXP>
    void TPZStream::Write(const TPZChunkVector<T, EXP> &vec) {
        long c, nc = vec.NElements();
        this->Write(&nc);
        for (c = 0; c < nc; c++)
            vec[c].Write(*this, 0);
    }

    template <class T, int EXP> void TPZStream::Write(TPZAdmChunkVector<T, EXP> &vec) {
        long c, nc = vec.NElements();
        this->Write(&nc);
        for (c = 0; c < nc; c++)
            vec[c].Write(*this, 0);
        this->Write(&vec.fCompactScheme);
        Write(vec.fFree);
        Write(vec.fNFree);
    }

    template <class T> void TPZStream::Write(const std::map<T, T> &vec) {
        long sz = vec.size();
        TPZManVector<T> cp(sz * 2);
        long count = 0;
		typename std::map<T, T>::const_iterator it;
        for (it = vec.begin(); it != vec.end(); it++) {
            cp[count++] = it->first;
            cp[count++] = it->second;
        }
        Write(cp);
    }

    /**
     * Write for chunk vectors with basic elements as float, double, long
     * double, std::complex<...> .
     */
    template <class T, int EXP>
    void TPZStream::Write(const TPZAdmChunkVector<T, EXP> &vec, bool basic) {
        if (!basic) {
            DebugStop();
            return;
        }
        long c, nc = vec.NElements();
        this->Write(&nc);
        for (c = 0; c < nc; c++)
            this->Write(&vec[c]);
        this->Write(&vec.fCompactScheme);
        Write(vec.fFree, true);
        Write(vec.fNFree, true);
    }

    template <class T> void TPZStream::Write(const TPZVec<T> &vec, bool basic) {
        if (!basic) {
            DebugStop();
            return;
        }
        long c, nc = vec.NElements();
        this->Write(&nc);
        for (c = 0; c < nc; c++)
            this->Write(&vec[c]);
    }

    template <class T> void TPZStream::WritePointers(TPZVec<T *> &vec) {
        long c, nc = vec.NElements();
        int emptyPos = -1;
        this->Write(&nc);
        for (c = 0; c < nc; c++) {
            if (vec[c]) {
                vec[c]->Write(*this);
            } else {
                this->Write(&emptyPos);
            }
        }
    }

    template <class T>
    void TPZStream::WritePointers(std::map<int, TPZAutoPointer<T>> &vec) {
        int nc = vec.size(), emptyPos = -1;
        this->Write(&nc);
        typedef typename std::map<int, TPZAutoPointer<T>>::iterator vec_it;
        vec_it it;
        for (it = vec.begin(); it != vec.end(); it++) {
            int id = it->first;
            this->Write(&id);
            if (it->second) {
                it->second->Write(*this, 1);
            } else {
                this->Write(&emptyPos);
            }
        }
    }

    template <class T> void TPZStream::WritePointers(std::map<int, T *> &vec) {
        int nc = vec.size(), emptyPos = -1;
        this->Write(&nc);
        typedef typename std::map<int, T *>::iterator vec_it;
        vec_it it;
        for (it = vec.begin(); it != vec.end(); it++) {
            int id = it->first;
            this->Write(&id);
            if (it->second) {
                it->second->Write(*this, 1);
            } else {
                this->Write(&emptyPos);
            }
        }
    }
    template <class T> void TPZStream::WritePointers(std::set<T *> &vec) {
        int nel = vec.size();
        this->Write(&nel);
        typedef typename std::set<T *>::iterator vec_it;
        vec_it it;
        while (it != vec.end()) {
            it->Write(*this, 1);
            it++;
        }
    }

    template <class T, int EXP>
    void TPZStream::WritePointers(TPZChunkVector<T *, EXP> &vec) {
        long c, nc = vec.NElements();
        int emptyPos = -1;
        this->Write(&nc);
        for (c = 0; c < nc; c++) {
            T *ptr = vec[c];
            if (ptr)
                ptr->Write(*this);
            else
                this->Write(&emptyPos);
        }
    }

    template <class T, int EXP>
    void TPZStream::WritePointers(TPZAdmChunkVector<T *, EXP> &vec) {
        long c, nc = vec.NElements();
        int emptyPos = -1;
        this->Write(&nc);
        for (c = 0; c < nc; c++) {
            T *ptr = vec[c];
            if (ptr)
                ptr->Write(*this, 1);
            else
                this->Write(&emptyPos);
        }
        this->Write(&vec.fCompactScheme);
        Write(vec.fFree);
        Write(vec.fNFree);
    }

    /**
     * @brief Methods to read objects or pointer for objects.
     */

    template <class T> void TPZStream::Read(std::vector<T> &vec, void *context) {
        int c, nc;
        this->Read(&nc, 1);
        vec.resize(nc);
        for (c = 0; c < nc; c++) {
            vec[c].Read(*this, context);
        }
    }

    template <class T> void TPZStream::Read(TPZVec<T> &vec, void *context) {
        long c, nc;
        this->Read(&nc, 1);
        vec.Resize(nc);
        for (c = 0; c < nc; c++) {
            vec[c].Read(*this, context);
        }
    }

    template <int N> void TPZStream::Read(TPZManVector<REAL, N> &vec) {
        long nc;
        this->Read(&nc, 1);
        vec.Resize(nc);
        if (nc)
            this->Read(&vec[0], nc);
    }

    template <class T, int EXP>
    void TPZStream::Read(TPZChunkVector<T, EXP> &vec, void *context) {
        long c, nc;
        this->Read(&nc, 1);
        vec.Resize(nc);
        for (c = 0; c < nc; c++) {
            vec[c].Read(*this, context);
        }
    }

    template <class T, int EXP>
    void TPZStream::Read(TPZAdmChunkVector<T, EXP> &vec, void *context) {
        long c, nc;
        this->Read(&nc, 1);
        vec.Resize(nc);
        for (c = 0; c < nc; c++)
            vec[c].Read(*this, context);
        this->Read(&vec.fCompactScheme, 1);
        Read(vec.fFree);
        Read(vec.fNFree);
    }

    template <class T> void TPZStream::Read(std::map<T, T> &vec) {
        TPZManVector<T> cp;
        Read(*this, cp);
        int sz = cp.NElements();
        int i;
        for (i = 0; i < sz; i += 2) {
            vec[cp[i]] = cp[i + 1];
        }
    }

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

    template <class T> void TPZStream::ReadPointers(TPZVec<T *> &vec, void *context) {
        long c, nc;
        this->Read(&nc, 1);
        vec.Resize(nc);
        for (c = 0; c < nc; c++) {
            vec[c] = dynamic_cast<T *>(TPZSaveable::Restore(*this, context));
        }
    }

    template <class T>
    void TPZStream::ReadPointers(std::map<int, TPZAutoPointer<T>> &vec,
                      void *context) {
        int c, nc;
        this->Read(&nc, 1);
        for (c = 0; c < nc; c++) {
            int id;
            this->Read(&id, 1);
            vec[id] =
                TPZAutoPointer<T>(dynamic_cast<T *>(TPZSaveable::Restore(*this, context)));
        }
    }

    template <class T>
    void TPZStream::ReadPointers(std::map<int, T *> &vec, void *context) {
        int c, nc;
        this->Read(&nc, 1);
        for (c = 0; c < nc; c++) {
            int id;
            this->Read(&id, 1);
            vec[id] = (dynamic_cast<T *>(TPZSaveable::Restore(*this, context)));
        }
    }

    template <class T, int EXP>
    void TPZStream::ReadPointers(TPZChunkVector<T *, EXP> &vec, void *context) {
        long c, nc;
        this->Read(&nc, 1);
        vec.Resize(nc);
        for (c = 0; c < nc; c++) {
            vec[c] = dynamic_cast<T *>(TPZSaveable::Restore(*this, context));
        }
    }

    template <class T, int EXP>
    void TPZStream::ReadPointers(TPZAdmChunkVector<T *, EXP> &vec, void *context) {
        long c, nc;
        this->Read(&nc, 1);
        vec.Resize(nc);
        for (c = 0; c < nc; c++) {
            vec[c] = dynamic_cast<T *>(TPZSaveable::Restore(*this, context));
        }
        this->Read(&vec.fCompactScheme, 1);
        Read(vec.fFree);
        Read(vec.fNFree);
    }