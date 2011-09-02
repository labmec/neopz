/**
 * @file
 * @brief Contains declaration of the TPZPermutation class which generates all permutations of n values.
 */
//
// C++ Interface: tpzpermutation
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef TPZPERMUTATION_H
#define TPZPERMUTATION_H

#include "pzmanvector.h"
class TPZStream;

/**
 @ingroup util
 @brief This class generates all permutations of n values. \ref util "Utility"
 @author Philippe R. B. Devloo
 */
class TPZPermutation{
public:
    TPZPermutation(int n);
    
    TPZPermutation(const TPZPermutation &copy);
	
    ~TPZPermutation();
    
    TPZPermutation &operator=(const TPZPermutation &copy);
    
    /** @brief Applies the current permutation on the vector in and produces the vector out */
    void Permute(const TPZVec<int> &in, TPZVec<int> &out) const;
    
    void operator++();
    
    void operator++(int) { operator++();}
    
    bool IsFirst();
    
    void Read(TPZStream &buf);
    
    void Write(TPZStream &buf);
	
	TPZManVector<int> Counter()
	{
		return fCounter;
	}
	
	TPZManVector<int>Order()
	{
		return fOrder;
	}
	
protected:
	
    /** @brief Variable which represents a counter for the permutations */
    TPZManVector<int> fCounter;
	
    /** @brief Variable which contains the current permutations */
    TPZManVector<int> fOrder;
};

/// Overloading operator << to print permutation object at ostream out
inline std::ostream &operator<<(std::ostream &out, TPZPermutation &obj)
{
	out << "imprimindo TPZPermutation::fCounter\n";
	for(int i = 0; i < obj.Counter().NElements(); i++)
	{
		out << obj.Counter()[i] << std::endl;
	}
	out << "\n\n";
	
	out << "imprimindo TPZPermutation::fOrder\n";
	for(int i = 0; i < obj.Order().NElements(); i++)
	{
		out << obj.Order()[i] << std::endl;
	}
	out << "\n\n";
	
	return out;
}

#endif
