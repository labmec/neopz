//
//  PZMatrixMarket.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/3/16.
//
//

#include "PZMatrixMarket.h"

#include "pzfmatrix.h"
#include "pzsbndmat.h"
#include "pzbndmat.h"
#include "pzskylmat.h"
#include "pzsysmp.h"


#include <fstream>

template<class T>
void TPZMatrixMarket::Read(std::string filename, TPZFMatrix<T>&fmat)
{
    int64_t nrow,ncol,nonzero=0,nsup=0,nlower=0;
    {
        std::ifstream input(filename.c_str());
        std::string buf;
        std::getline(input,buf);
        input >> nrow >> ncol >> nonzero;
        fmat.Redim(nrow,ncol);
        for (int64_t el=0; el<nonzero; el++) {
            int64_t row, col;
            T val;
            input >> row >> col >> val;
            row--;
            col--;
            if(row > col) nlower++;
            if(row< col ) nsup++;
            fmat.PutVal(row,col,val);
        }
    }
    if(nsup == 0 || nlower == 0)
    {
        std::ifstream input(filename);
        std::string buf;
        std::getline(input,buf);
        input >> nrow >> ncol >> nonzero;
        fmat.Redim(nrow,ncol);
        for (int64_t el=0; el<nonzero; el++) {
            int64_t row, col;
            T val;
            input >> row >> col >> val;
            row--;
            col--;
            int64_t swap = row;
            row = col;
            col = swap;
            fmat.PutVal(row,col,val);
        }
    }
}

template<class T>
void TPZMatrixMarket::Read(std::string filename, TPZFBMatrix<T>&fmat)
{
    int64_t nrow,ncol,nonzero=0,nsup=0,nlower=0, band=0;
    {
        std::ifstream input(filename);
        std::string buf;
        std::getline(input,buf);
        input >> nrow >> ncol >> nonzero;
        fmat.Redim(nrow,ncol);
        for (int64_t el=0; el<nonzero; el++) {
            int64_t row, col;
            T val;
            input >> row >> col >> val;
            if(row-col > band) band = row-col;
            if(col-row > band) band = col-row;
            row--;
            col--;
            if(row > col) nlower++;
            if(row< col ) nsup++;
        }
    }
    fmat.Redim(nrow,ncol);
    fmat.SetBand(band);
    {
        std::ifstream input(filename);
        std::string buf;
        std::getline(input,buf);
        input >> nrow >> ncol >> nonzero;
        for (int64_t el=0; el<nonzero; el++) {
            int64_t row, col;
            T val;
            input >> row >> col >> val;
            row--;
            col--;
            int64_t swap = row;
            row = col;
            col = swap;
            fmat.PutVal(row,col,val);
            if (nsup == 0 || nlower == 0) {
                fmat.PutVal(col,row,val);
            }
        }
    }
}

template<class T>
void TPZMatrixMarket::Read(std::string filename, TPZSBMatrix<T>&fmat)
{
    int64_t nrow,ncol,nonzero=0,nsup=0,nlower=0, band=0;
    {
        std::ifstream input(filename);
        std::string buf;
        std::getline(input,buf);
        input >> nrow >> ncol >> nonzero;
        for (int64_t el=0; el<nonzero; el++) {
            int64_t row, col;
            T val;
            input >> row >> col >> val;
            if(row-col > band) band = row-col;
            if(col-row > band) band = col-row;
            row--;
            col--;
            if(row > col) nlower++;
            if(row< col ) nsup++;
        }
    }
    fmat.Redim(nrow,ncol);
    fmat.SetBand(band);
    {
        std::ifstream input(filename);
        std::string buf;
        std::getline(input,buf);
        input >> nrow >> ncol >> nonzero;
        for (int64_t el=0; el<nonzero; el++) {
            int64_t row, col;
            T val;
            input >> row >> col >> val;
            row--;
            col--;
            fmat.PutVal(row,col,val);
        }
    }
}

template<class T>
void TPZMatrixMarket::Read(std::string filename, TPZSkylMatrix<T>&fmat)
{
    int64_t nrow,ncol,nonzero=0,nsup=0,nlower=0;
    TPZVec<int64_t> skyline(nrow,0);
    {
        std::ifstream input(filename.c_str());
        std::string buf;
        std::getline(input,buf);
        input >> nrow >> ncol >> nonzero;
        if (nrow != ncol) {
            DebugStop();
        }
        skyline.resize(nrow);
        for (int64_t el=0; el<nrow; el++) {
            skyline[el] = el;
        }
        fmat.Redim(nrow,ncol);
        for (int64_t el=0; el<nonzero; el++) {
            int64_t row, col;
            T val;
            input >> row >> col >> val;
            row--;
            col--;
            if(row < skyline[col]) skyline[col] = row;
            if(col < skyline[row]) skyline[row] = col;
            if(row > col) nlower++;
            if(row< col ) nsup++;
        }
    }
    if (nlower && nsup) {
        std::cout << "The matrix is probably not symmetric, expect trouble\n";
    }
    fmat.Redim(nrow,ncol);
    fmat.SetSkyline(skyline);
    {
        std::ifstream input(filename.c_str());
        std::string buf;
        std::getline(input,buf);
        input >> nrow >> ncol >> nonzero;
        for (int64_t el=0; el<nonzero; el++) {
            int64_t row, col;
            T val;
            input >> row >> col >> val;
            row--;
            col--;
            fmat.PutVal(row,col,val);
        }
    }
}

template<class T>
void TPZMatrixMarket::Read(std::string filename, TPZSYsmpMatrix<T> &fmat)
{
    int64_t nrow,ncol,nonzero=0,nsup=0,nlower=0;
    TPZVec<int64_t> numrowel(nrow,0);
    {
        std::ifstream input(filename);
        std::string buf;
        std::getline(input,buf);
        input >> nrow >> ncol >> nonzero;
        if (nrow != ncol) {
            DebugStop();
        }
        numrowel.resize(nrow);
        numrowel.Fill(0);
        fmat.Redim(nrow,ncol);
        for (int64_t el=0; el<nonzero; el++) {
            int64_t row, col;
            T val;
            input >> row >> col >> val;
            row--;
            col--;
            numrowel[row]++;
            if(row > col) nlower++;
            if(row< col ) nsup++;
        }
    }
    if (nlower && nsup) {
        std::cout << "The matrix is probably not symmetric, expect trouble\n";
    }
    TPZVec<int64_t> IA(nrow+1);
    TPZVec<int64_t> JA(nonzero,0);
    TPZVec<T> A(nonzero,0.);
    TPZVec<int64_t> IAcounter(nonzero);
    IA[0] = 0;
    for (int64_t i=0; i<nrow; i++) {
        IA[i+1] = IA[i]+numrowel[i];
    }
    IAcounter = IA;
    fmat.Redim(nrow,ncol);
    {
        std::ifstream input(filename);
        std::string buf;
        std::getline(input,buf);
        input >> nrow >> ncol >> nonzero;
        for (int64_t el=0; el<nonzero; el++) {
            int64_t row, col;
            T val;
            input >> row >> col >> val;
            row--;
            col--;
            JA[IAcounter[row]] = col;
            A[IAcounter[row]] = val;
            IAcounter[row]++;
        }
    }
    fmat.SetData(IA,JA,A);
}

#define IMPLEMENTMMARKET(T) \
template void TPZMatrixMarket::Read<T>(std::string filename,\
                                       TPZFMatrix<T> &fmat);\
template void TPZMatrixMarket::Read<T>(std::string filename,\
                                       TPZFBMatrix<T> &fmat);\
template void TPZMatrixMarket::Read<T>(std::string filename,\
                                       TPZSBMatrix<T> &fmat);\
template void TPZMatrixMarket::Read<T>(std::string filename,\
                                       TPZSkylMatrix<T> &fmat);\
template void TPZMatrixMarket::Read<T>(std::string filename,\
                                       TPZSYsmpMatrix<T> &fmat);

IMPLEMENTMMARKET(float)
IMPLEMENTMMARKET(double)
IMPLEMENTMMARKET(long double)
IMPLEMENTMMARKET(std::complex<float>)
IMPLEMENTMMARKET(std::complex<double>)