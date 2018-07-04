//
//  TNRFElasticMemory.h
//  pz
//
//  Created by Omar Durán on 7/2/18.
//

#ifndef TNRFElasticMemory_h
#define TNRFElasticMemory_h

#include <stdio.h>

#include "TPZTensor.h"
#include "TPZPlasticState.h"

class TNRFElasticMemory : public TPZSavable
{
    
public:
    
    /// Constructor
    TNRFElasticMemory();
    
    /// Copy constructor
    TNRFElasticMemory(const TNRFElasticMemory & other);
    
    /// assignment operator
    const TNRFElasticMemory & operator=(const TNRFElasticMemory & other);
    
    // Destructor
    virtual ~TNRFElasticMemory();
    
    // print memory
    void Print(std::ostream &out) const;
    
    /// Class Identifier
    virtual int ClassId() const;
    
};

// print memory inline
inline std::ostream &operator<<(std::ostream &out, const TNRFElasticMemory &memory)
{
    memory.Print(out);
    return out;
}


#endif /* TNRFElasticMemory_h */
