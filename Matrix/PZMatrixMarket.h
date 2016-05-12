//
//  PZMatrixMarket.hpp
//  PZ
//
//  Created by Philippe Devloo on 5/3/16.
//
//

#ifndef PZMatrixMarket_hpp
#define PZMatrixMarket_hpp

#include <stdio.h>
#include "pzfmatrix.h"
#include "pzsbndmat.h"
#include "pzbndmat.h"
#include "pzskylmat.h"
#include "pzsysmp.h"
#include <fstream>

class TPZMatrixMarket
{
public:
    
    template<class T>
    static void Read(std::string filename, TPZFMatrix<T> &fmat);
    
    template<class T>
    static void Read(std::string filename, TPZFBMatrix<T> &fmat);
    
    template<class T>
    static void Read(std::string filename, TPZSBMatrix<T> &fmat);
    
    template<class T>
    static void Read(std::string filename, TPZSkylMatrix<T> &fmat);
    
    template<class T>
    static void Read(std::string filename, TPZSYsmpMatrix<T> &fmat);
};

template<class T>
void TPZMatrixMarket::Read(std::string filename, TPZFMatrix<T>&fmat)
{
    long nrow,ncol,nonzero=0,nsup=0,nlower=0;
    {
        std::ifstream input(filename);
        std::string buf;
        std::getline(input,buf);
        input >> nrow >> ncol >> nonzero;
        fmat.Redim(nrow,ncol);
        for (long el=0; el<nonzero; el++) {
            long row, col;
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
        for (long el=0; el<nonzero; el++) {
            long row, col;
            T val;
            input >> row >> col >> val;
            row--;
            col--;
            long swap = row;
            row = col;
            col = swap;
            fmat.PutVal(row,col,val);
        }
    }
}

template<class T>
void TPZMatrixMarket::Read(std::string filename, TPZFBMatrix<T>&fmat)
{
    long nrow,ncol,nonzero=0,nsup=0,nlower=0, band=0;
    {
        std::ifstream input(filename);
        std::string buf;
        std::getline(input,buf);
        input >> nrow >> ncol >> nonzero;
        fmat.Redim(nrow,ncol);
        for (long el=0; el<nonzero; el++) {
            long row, col;
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
        for (long el=0; el<nonzero; el++) {
            long row, col;
            T val;
            input >> row >> col >> val;
            row--;
            col--;
            long swap = row;
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
    long nrow,ncol,nonzero=0,nsup=0,nlower=0, band=0;
    {
        std::ifstream input(filename);
        std::string buf;
        std::getline(input,buf);
        input >> nrow >> ncol >> nonzero;
        for (long el=0; el<nonzero; el++) {
            long row, col;
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
        for (long el=0; el<nonzero; el++) {
            long row, col;
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
    long nrow,ncol,nonzero=0,nsup=0,nlower=0;
    TPZVec<long> skyline(nrow,0);
    {
        std::ifstream input(filename);
        std::string buf;
        std::getline(input,buf);
        input >> nrow >> ncol >> nonzero;
        if (nrow != ncol) {
            DebugStop();
        }
        skyline.resize(nrow);
        for (long el=0; el<nrow; el++) {
            skyline[el] = el;
        }
        fmat.Redim(nrow,ncol);
        for (long el=0; el<nonzero; el++) {
            long row, col;
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
        std::ifstream input(filename);
        std::string buf;
        std::getline(input,buf);
        input >> nrow >> ncol >> nonzero;
        for (long el=0; el<nonzero; el++) {
            long row, col;
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
    long nrow,ncol,nonzero=0,nsup=0,nlower=0;
    TPZVec<long> numrowel(nrow,0);
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
        for (long el=0; el<nonzero; el++) {
            long row, col;
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
    TPZVec<long> IA(nrow+1);
    TPZVec<long> JA(nonzero,0);
    TPZVec<T> A(nonzero,0.);
    TPZVec<long> IAcounter(nonzero);
    IA[0] = 0;
    for (long i=0; i<nrow; i++) {
        IA[i+1] = IA[i]+numrowel[i];
    }
    IAcounter = IA;
    fmat.Redim(nrow,ncol);
    {
        std::ifstream input(filename);
        std::string buf;
        std::getline(input,buf);
        input >> nrow >> ncol >> nonzero;
        for (long el=0; el<nonzero; el++) {
            long row, col;
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


#endif /* PZMatrixMarket_hpp */
