/* 
 * File:   TPZAssert.h
 * Author: thiago
 *
 * Created on 26 de Fevereiro de 2018, 13:49
 */

#ifndef TPZASSERT_H
#define TPZASSERT_H

#include "pzerror.h"
#include "pzextractval.h"

class TPZAssert {
public:
    
    template<typename T>
    static T& NonNegative(T &value) {
        //Note that IsZero uses a tolerance.
        if (value < T(0) && TPZExtractVal::IsZero(value)){
            SetValue(value, 0);
        }
        return value;
    }
    
    template<typename T>
    static T NonNegative(const T &value) {
        T copy(value);
        return NonNegative(copy);
    }
    
private:

    template<typename T, typename TValue>
    static void SetValue(T &var, const TValue value) {
        TPZExtractVal::ref(var) = value;
    }
    
};

#endif /* TPZASSERT_H */

