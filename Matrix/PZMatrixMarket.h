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

#endif /* PZMatrixMarket_hpp */
