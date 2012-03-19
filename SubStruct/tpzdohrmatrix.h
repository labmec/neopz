/**
 * @file
 * @brief Contains the TPZDohrMatrix class which implements a matrix divided into substructures. \n
 * Also contains the TPZDohrThreadMultData and TPZDohrThreadMultList structs.
 */
/***************************************************************************
 *   Copyright (C) 2006 by Philippe Devloo   *
 *   phil@fec.unicamp.br   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef TPZDOHRMATRIX_H
#define TPZDOHRMATRIX_H

#include "pzmatrix.h"
#include <list>
#include <sstream>
#include "tpzautopointer.h"
#include "tpzdohrsubstruct.h"
#include "tpzdohrsubstructCondense.h"
#include "tpzdohrassembly.h"

#include "tpzdohrassemblelist.h"

/**
 * @brief Implements a matrix divided into substructures. \ref matrix "Matrix" \ref substructure "Sub structure"
 * @ingroup substructure matrix
 * @author Philippe Devloo
 */
template <class TSubStruct> 
class TPZDohrMatrix : public TPZMatrix<REAL>
{
public:
	/** @brief The matrix class is a placeholder for a list of substructures */
	typedef typename std::list<TPZAutoPointer<TSubStruct> > SubsList;
private:
	SubsList fGlobal;
	
	int fNumCoarse; //n(c)
	
	/** @brief Number of threads that will be used during the matrix vector multiplication */
	int fNumThreads;
	
public:
	
	TPZAutoPointer<TPZDohrAssembly> fAssembly;
	
	TPZDohrMatrix(TPZAutoPointer<TPZDohrAssembly> dohrassembly);
	
	TPZDohrMatrix() : TPZMatrix<REAL>()
    {
        
    }
	
	TPZDohrMatrix(const TPZDohrMatrix &cp) : fGlobal(cp.fGlobal), fNumCoarse(cp.fNumCoarse), fNumThreads(cp.fNumThreads), 
	fAssembly(cp.fAssembly)
	{
	}
	
//	CLONEDEF(TPZDohrMatrix)
		virtual TPZMatrix<REAL>*Clone() const { return new TPZDohrMatrix(*this); }
	
	~TPZDohrMatrix();
	
	const SubsList &SubStructures() const
	{
		return fGlobal;
	}
	
	int NumCoarse() const
	{
		return fNumCoarse;
	}
	
	int NumThreads() const
	{
		return fNumThreads;
	}
	
	void SetNumThreads(int numthreads)
	{
		fNumThreads = numthreads;
	}
	
	/** @brief Just a method for tests */
	TPZAutoPointer<TSubStruct> GetFirstSub() {
		return (*fGlobal.begin());
	}
	/** @brief Just a method for tests */
	void Print(const char *name, std::ostream& out,const MatrixOutputFormat form = EFormatted) const
	{
		out << __PRETTY_FUNCTION__ << std::endl;
		out << name << std::endl;
		out << "Number of coarse equations " << fNumCoarse << std::endl;
		typename SubsList::const_iterator iter;
		for (iter=fGlobal.begin();iter!=fGlobal.end();iter++) {
			(*iter)->Print(out);
		}
	}
	
	void SetNumCornerEqs(int nc)
	{
		fNumCoarse = nc;
	}
	/** @brief Initialize the necessary datastructures */
	void Initialize();   
	/** @brief It adds a substruct */
	void AddSubstruct(TPZAutoPointer<TSubStruct> substruct)
	{
		fGlobal.push_back(substruct);
	}
	
	/**
	 * @brief It computes z = beta * y + alpha * opt(this)*x but z and x can not overlap in memory.
	 * @param x Is x on the above operation
	 * @param y Is y on the above operation
	 * @param z Is z on the above operation
	 * @param alpha Is alpha on the above operation
	 * @param beta Is beta on the above operation
	 * @param opt Indicates if is Transpose or not
	 * @param stride Indicates n/N where n is dimension of the right hand side vector and N is matrix dimension
	 */
	/** The only method any matrix class needs to implement */
	virtual void MultAdd(const TPZFMatrix<REAL> &x,const TPZFMatrix<REAL> &y, TPZFMatrix<REAL> &z,const REAL alpha,const REAL beta,const int opt,const int stride) const;
	
	/** @brief Adjust the residual to zero the residual of the internal connects */
	void AdjustResidual(TPZFMatrix<REAL> &res);
	
	/** @brief Add the solution corresponding to the internal residual */
	void AddInternalSolution(TPZFMatrix<REAL> &solution);
	
};

/**
 * @ingroup substructure
 * @brief .. . \ref substructure "Sub structure"
 */
template <class TSubStruct> 
struct TPZDohrThreadMultData
{
	/** @brief Default constructor */
	TPZDohrThreadMultData() : fisub(-1), fSub(0)
	{
	}
	TPZDohrThreadMultData(int isub, TPZAutoPointer<TSubStruct> submesh) : fisub(isub), fSub(submesh)
	{
	}
	/** @brief Copy constructor */
	TPZDohrThreadMultData(const TPZDohrThreadMultData<TSubStruct> &cp) : fisub(cp.fisub), fSub(cp.fSub)
	{
	}
	/** @brief Implement the attribution operator. */
	TPZDohrThreadMultData<TSubStruct> &operator=(const TPZDohrThreadMultData<TSubStruct> &cp)
	{
		fisub = cp.fisub;
		fSub = cp.fSub;
		return *this;
	}
	int fisub;
	TPZAutoPointer<TSubStruct> fSub;	
	
	bool IsValid()
	{
		return (fisub >= 0);
	}
};

/**
 * @ingroup substructure
 * @brief .. . \ref substructure "Sub structure"
 */
template <class TSubStruct> 
struct TPZDohrThreadMultList
{
	/** @brief The vector with which we will multiply */
	const TPZFMatrix<REAL> *fInput;
	/** @brief Scalar multiplication factor */
	REAL fAlpha;
	/** @brief Mutex which will enable the access protection of the list */
	pthread_mutex_t fAccessLock;
	/** @brief The data structure which defines the assemble destinations */
	TPZAutoPointer<TPZDohrAssembly> fAssembly;
	/** @brief The list of data objects which need to treated by the threads */
	std::list<TPZDohrThreadMultData<TSubStruct> > fWork;
	/** @brief The local contribution to the v2 vector */
	TPZAutoPointer<TPZDohrAssembleList> fAssemblyStructure;
	
	TPZDohrThreadMultList(const TPZFMatrix<REAL> &x, REAL alpha, TPZAutoPointer<TPZDohrAssembly> assembly, TPZAutoPointer<TPZDohrAssembleList> &assemblestruct) : fInput(&x), fAlpha(alpha), 
	fAssembly(assembly), fAssemblyStructure(assemblestruct)
	{
		pthread_mutex_init(&fAccessLock, 0);
	}
	TPZDohrThreadMultList()
	{
		pthread_mutex_destroy(&fAccessLock);
	}
	
	/** @brief The procedure which executes the lengthy process */
	static void *ThreadWork(void *voidptr);
	/** @brief Interface to add items in a thread safe way */
	void AddItem(TPZDohrThreadMultData<TSubStruct> &data)
	{
		pthread_mutex_lock(&fAccessLock);
		fWork.push_back(data);
		pthread_mutex_unlock(&fAccessLock);
	}
	/** @brief Interface to pop an item in a thread safe way */
	TPZDohrThreadMultData<TSubStruct> PopItem()
	{
		TPZDohrThreadMultData<TSubStruct> result;
		pthread_mutex_lock(&fAccessLock);
		if (fWork.size()) {
			result = *fWork.begin();
			fWork.pop_front();
		}
		pthread_mutex_unlock(&fAccessLock);
		return result;
	}
};

#endif
