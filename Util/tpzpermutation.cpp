/** 
 * @file
 * @brief Contains the implementation of the methods to TPZPermutation class.
 */
//
// C++ Implementation: tpzpermutation
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "tpzpermutation.h"
#include "pzsave.h"

TPZPermutation::TPZPermutation(int n) : fCounter(n,0), fOrder(n,-1)
{
	int i;
	for(i=0; i<n; i++) fOrder[i] = i;
}

TPZPermutation::TPZPermutation(const TPZPermutation &copy) : fCounter(copy.fCounter), fOrder(copy.fOrder)
{
}

TPZPermutation &TPZPermutation::operator=(const TPZPermutation &copy)
{
	fCounter = copy.fCounter;
	fOrder = copy.fOrder;
	return *this;
}

TPZPermutation::~TPZPermutation()
{
}

void TPZPermutation::Read(TPZStream &buf){
	TPZSaveable::ReadObjects(buf, this->fCounter);
	TPZSaveable::ReadObjects(buf, this->fOrder);
}

void TPZPermutation::Write(TPZStream &buf){
	TPZSaveable::WriteObjects(buf, this->fCounter);
	TPZSaveable::WriteObjects(buf, this->fOrder);
}

/// Applies the current permutation on the vector in and produces the vector out
void TPZPermutation::Permute(const TPZVec<int> &in, TPZVec<int> &out) const
{
	int i,n=fCounter.NElements();
	for(i=0; i<n; i++) out[fOrder[i]] = in[i];
}

void TPZPermutation::operator++()
{
	int n = fCounter.NElements();
	int i = n-2;
	while(i >= 0 && ((++fCounter[i]) %= n-i) == 0) i--;
	fOrder.Fill(-1);
	int count;
	for(i=0; i<n; i++)
	{
		count = 0;
		int index = 0;
		while(fCounter[i] > count || fOrder[index] != -1)
		{
			if(fOrder[index] == -1) count++;
			index++;
		}
		fOrder[index] = i;
	}
}

bool TPZPermutation::IsFirst()
{
	int nel = fCounter.NElements();
	int in;
	for(in=0; in<nel; in++) if(fCounter[in] != 0) return false;
	return true;
}
