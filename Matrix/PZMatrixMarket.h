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

#include <string>
template<class TVar>
class TPZFMatrix;

template<class TVar>
class TPZFBMatrix;

template<class TVar>
class TPZSBMatrix;

template<class TVar>
class TPZSkylMatrix;

template<class TVar>
class TPZSYsmpMatrix;


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

#endif /* PZMatrixMarket_hpp */
