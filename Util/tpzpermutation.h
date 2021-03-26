/**
 * @file
 * @brief Contains declaration of the TPZPermutation class which generates all permutations of n values.
 */
#ifndef TPZPERMUTATION_H
#define TPZPERMUTATION_H

#include "pzmanvector.h"
#include "TPZSavable.h"
class TPZStream;

/**
 * @ingroup util
 * @brief This class generates all permutations of n values. \ref util "Utility"
 * @author Philippe R. B. Devloo
 */
class TPZPermutation : public TPZSavable {
public:
	/** @brief Constructor with number of permutations */
    TPZPermutation(int n);
    /** @brief Copy constructor */
    TPZPermutation(const TPZPermutation &copy);
	/** @brief Default destructor */
    ~TPZPermutation();

    /** @brief Operator attribution */
    TPZPermutation &operator=(const TPZPermutation &copy);
    
    /** @brief Applies the current permutation on the vector in and produces the vector out */
    void Permute(const TPZVec<int> &in, TPZVec<int> &out) const;
    void Permute(const TPZVec<int64_t> &in, TPZVec<int64_t> &out) const;

    /** @brief Operator increment */
    void operator++();
    
    void operator++(int) { operator++();}
    
    bool IsFirst();
        int ClassId() const override;
    void Read(TPZStream &buf, void *context) override;
    void Write(TPZStream &buf, int withclassid) const override;
    
	/** @brief Returns the counter of the permutations */
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

/** @brief Overloading operator << to print permutation object at ostream out */
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
