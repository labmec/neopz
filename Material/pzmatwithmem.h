/**
 * \file
 * @brief Contains the TPZMatWithMem class which implements the memory features.
 */
//$Id: pzmatwithmem.h,v 1.6 2009-04-19 01:47:51 erick Exp $

#ifndef PZMATWITHMEM_H
#define PZMATWITHMEM_H

#include "pzmaterial.h"
#include "pzadmchunk.h"
#include "pzviscoelastic.h"
//#include "poroelastoplasticid.h"

/**
 * @ingroup material
 * @brief Implements an abstract class implementing the memory features.
 */
/**
 * Materials that aim to use memory shoud be derived from this class.
 */

template <class TMEM, class TFather = TPZMaterial>
class  TPZMatWithMem : public TFather
{
public:
	
	/** @brief Default constructor */
	TPZMatWithMem();
	
	/** @brief Creates a material object and inserts it in the vector of material pointers of the mesh */
	/** Upon return vectorindex contains the index of the material object within the vector */
	TPZMatWithMem(int id);
	
	/** 
	 * @brief Creates a material object based on the referred object and
	 * inserts it in the vector of material pointers of the mesh
	 */
	/** Upon return vectorindex contains the index of the material object within the vector */
	TPZMatWithMem(const TPZMatWithMem<TMEM, TFather> &mat);
	
	virtual ~TPZMatWithMem();
	
	/** @brief Returns the name of the material */
	virtual std::string Name() { return "TPZMatWithMem< >"; }
	
	/** @brief Prints out the data associated with the material */
	virtual void PrintMem(std::ostream &out = std::cout, const int memory = 0);
	
	virtual TMEM & MemItem(const int i);
	
public:
	
	/** @brief Unique identifier for serialization purposes */
	virtual int ClassId() const;
	
	virtual void Write(TPZStream &buf, int withclassid);
	
	virtual void Read(TPZStream &buf, void *context);
	
	/**
	 * @brief Pushes a new entry in the context of materials with memory
	 * @return Returning its index at the internal storage stack
	 */
	/** To be implemented only in the proper materials. */
	virtual int PushMemItem(int sourceIndex = -1);
	
	/** @brief Frees an entry in the material with memory internal history storage */
	virtual void FreeMemItem(int index);
	
	/** @brief Sets the default memory settings for initialization */
	virtual void SetDefaultMem(TMEM & defaultMem);
	
	/** @brief Sets/Unsets the internal memory data to be updated in the next assemble/contribute call */
	virtual void SetUpdateMem(int update = 1);
	
protected:
	
	
	/** @brief Material Memory */
	TPZAdmChunkVector<TMEM/*, 1024 Using Default*/> fMemory;
	
	/** @brief Default memory settings */
	TMEM fDefaultMem;
	
	/** @brief Flag to indicate wether the memory data are to be updated in an assemble loop */
	int fUpdateMem;
};

template <class TMEM, class TFather>
TPZMatWithMem<TMEM,TFather>::TPZMatWithMem() : TFather(), 
fMemory(), fDefaultMem(), fUpdateMem(0){ }	

template <class TMEM, class TFather>
TPZMatWithMem<TMEM,TFather>::TPZMatWithMem(int id) : TFather(id),
fMemory(), fDefaultMem(), fUpdateMem(0) { }

template <class TMEM, class TFather>
TPZMatWithMem<TMEM,TFather>::TPZMatWithMem(const TPZMatWithMem<TMEM, TFather> &mat) : TFather(mat),
fMemory(mat.fMemory), fDefaultMem(mat.fDefaultMem), fUpdateMem(mat.fUpdateMem){ }

template <class TMEM, class TFather>
TPZMatWithMem<TMEM,TFather>::~TPZMatWithMem(){ }

template <>
inline void TPZMatWithMem<TPZFMatrix,TPZElasticity3D>::PrintMem(std::ostream &out, const int memory)
{
	
    out << "\nTPZMatWithMem<TPZFMatrix,TPZElasticity3D> Material\n";
    out << "\n fDefaultMem = \n" << fDefaultMem;
    out << "\n fUpdateMem = " << fUpdateMem;
    int i, size = fMemory.NElements();
    out << "\n fMemory with " << size << " elements";
    if(memory)
    {
        out << "\n fMemory elements:";
        for(i = 0; i < size; i++)
        {
            out << "\n " << i << ": "; 
            fMemory[i].Print("visc",out);
        }
    }
    out << "\nEnd of TPZMatWithMem<TPZFMatrix,TPZElasticity3D >::Print\n";
    TPZElasticity3D::Print(out);
}

template <class TMEM, class TFather>
void TPZMatWithMem<TMEM,TFather>::PrintMem(std::ostream &out, const int memory)
{
	TMEM fooMEM;
	//out << "\nTPZMatWithMem< " << fooMEM.Name() << " > Material\n";
	out << "\n fDefaultMem = \n" << fDefaultMem;
	out << "\n fUpdateMem = " << fUpdateMem;
	int i, size = fMemory.NElements();
	out << "\n fMemory with " << size << " elements";
	if(memory)
	{
		out << "\n fMemory elements:";
		for(i = 0; i < size; i++)
		{
			out << "\n " << i << ": "; 
			fMemory[i].Print(out);
		}
	}
	out << "\nEnd of TPZMatWithMem< " << fooMEM.Name() << " >::Print\n";
    TFather::Print(out);
}

template <class TMEM, class TFather>
TMEM & TPZMatWithMem<TMEM,TFather>::MemItem(const int i)
{
	return fMemory[i];
}

template <class TMEM, class TFather>
int TPZMatWithMem<TMEM,TFather>::ClassId() const
{
	return -1;
}

template <class TMEM, class TFather>
void TPZMatWithMem<TMEM,TFather>::Write(TPZStream &buf, int withclassid)
{	
	TPZMaterial::Write(buf, withclassid);
	
	int size = fMemory.NElements();
    buf.Write(&size,1);
	int i;
	for(i = 0; i < size; i++)
		fMemory[i].Write(buf, withclassid);
}

template <class TMEM, class TFather>
void TPZMatWithMem<TMEM,TFather>::Read(TPZStream &buf, void *context)
{
	TPZMaterial::Read(buf, context);
	
	int i,size;
	buf.Read(&size,1);
	for(i = 0; i < size; i++)
		fMemory[i].Read(buf, context);
	
}

template <class TMEM, class TFather>
int TPZMatWithMem<TMEM,TFather>::PushMemItem(int sourceIndex)
{
	int index = fMemory.AllocateNewElement();
	if(sourceIndex < 0)
	{
		fMemory[index] = fDefaultMem;
	}else{
		fMemory[index] = fMemory[sourceIndex];
	}
	return index;
}

template <class TMEM, class TFather>
void TPZMatWithMem<TMEM,TFather>::FreeMemItem(int index)
{
	fMemory.SetFree(index);
}

template <class TMEM, class TFather>
void TPZMatWithMem<TMEM,TFather>::SetDefaultMem(TMEM & defaultMem)
{
	fDefaultMem = defaultMem;
}

template <class TMEM, class TFather>
void TPZMatWithMem<TMEM,TFather>::SetUpdateMem(int update)
{
	fUpdateMem = update;
}

#endif
