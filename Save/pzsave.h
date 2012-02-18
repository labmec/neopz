/**
 * @file
 * @brief Contains declaration of the TPZSaveable class which defines the interface to save and restore objects from TPZStream objects.
 */
#ifndef PZSAVEH
#define PZSAVEH

#include <map>
#include <vector>


#include "pzvec.h"
#include "pzmanvector.h"
#include "pzadmchunk.h"
#include "pzfilebuffer.h"
#include "pzreal.h"
#include "tpzautopointer.h"

class TPZStream;
class TPZSaveable;

/**
 * \addtogroup save
 * @{
 */
/// Identifier as saveable object
const int TPZSAVEABLEID = -1;

/// Typedef of TPZRestore_t
typedef TPZSaveable *(*TPZRestore_t)(TPZStream &,void *);


/**
 * @brief This class defines the interface to save and restore objects from TPZStream objects. \ref save "Persistency"
 */
/**
 * This class defines the interface a class needs to implement (override) in order to become persistent
 * Several static utility methods have been defined to make saving and restoring of vectors of
 * objects easier
 */
class TPZSaveable {
	
#ifndef ELLIPS
	
	/** @brief This static function garantes that the gMap object is available when needed */
	static std::map<int,TPZRestore_t> &Map() {
		static std::map<int,TPZRestore_t> gMap;
		return gMap;
	}
#endif
	
public:
	
	virtual ~TPZSaveable()
	{
	}
	
	
	/** @brief Define the class id associated with the class */
	/**
	 * This id has to be unique for all classes
	 * A non unique id is flagged at the startup of the program
	 */
	virtual int ClassId() const ;
	
	/** @brief Writes this object to the TPZStream buffer. Include the classid if withclassid = true */
	virtual void Write(TPZStream &buf, int withclassid);
	
	/** @brief read objects from the stream */
	virtual void Read(TPZStream &buf, void *context);
	
	/** @brief Compares the object for identity with the object pointed to, eventually copy the object */
	/**
	 * Compares both objects bitwise for identity. Put an entry in the log file if different
	 * overwrite \n the calling object if the override flag is true
	 */
	virtual bool Compare(TPZSaveable *copy, bool override = false);
	
	/** @brief Compares the object for identity with the object pointed to, eventually copy the object */
	/**
	 * Compares both objects bitwise for identity. Put an entry in the log file if different
	 * generate \n an interrupt if the override flag is true
	 */
	virtual bool Compare(TPZSaveable *copy, bool override = false) const;
	
	template<class T>
	static void WriteObjects(TPZStream &buf, TPZVec<T> &vec)
	{
		int c,nc = vec.NElements();
		buf.Write(&nc,1);
		for(c=0; c<nc; c++) vec[c].Write(buf,0);
	}
	
	
	template<class T>
	static void WriteObjects(TPZStream &buf, std::vector<T> &vec)
	{
		int c,nc = vec.size();
		buf.Write(&nc,1);
		for(c=0; c<nc; c++) vec[c].Write(buf,0);
	}
	
	template<class T, int EXP>
	static void WriteObjects(TPZStream &buf, TPZChunkVector<T,EXP> &vec)
	{
		int c,nc = vec.NElements();
		buf.Write(&nc,1);
		for(c=0; c<nc; c++) vec[c].Write(buf,0);
	}
	
	template<class T, int EXP>
	static void WriteObjects(TPZStream &buf, TPZAdmChunkVector<T,EXP> &vec)
	{
		int c,nc = vec.NElements();
		buf.Write(&nc,1);
		for(c=0; c<nc; c++) vec[c].Write(buf,0);
		buf.Write(&vec.fCompactScheme,1);
		WriteObjects(buf,vec.fFree);
		WriteObjects(buf,vec.fNFree);
	}
	
	static void WriteObjects(TPZStream &buf, const std::map<int,int> &vec)
	{
		int sz = vec.size();
		TPZManVector<int> cp(sz*2);
		int count = 0;
        std::map<int,int>::const_iterator it;
		for (it=vec.begin(); it!=vec.end(); it++) {
			cp[count++] = it->first;
			cp[count++] = it->second;
		}
		WriteObjects(buf, cp);
	}

    static void WriteObjects(TPZStream &buf, const std::map<REAL,REAL> &vec)
	{
		int sz = vec.size();
		TPZManVector<REAL> cp(sz*2);
		int count = 0;
		std::map<REAL,REAL>::const_iterator it;
		for (it=vec.begin(); it!=vec.end(); it++) {
			cp[count++] = it->first;
			cp[count++] = it->second;
		}
		WriteObjects(buf, cp);
	}

	
	template<class T>
	static void WriteObjectPointers(TPZStream &buf, TPZVec<T *> &vec)
	{
		int c,nc = vec.NElements(),one = -1;
		buf.Write(&nc,1);
		for(c=0; c<nc; c++)
		{
			if(vec[c])
			{
				vec[c]->Write(buf);
			} else {
				buf.Write(&one,1);
			}
		}
	}
	template<class T>
	static void WriteObjectPointers(TPZStream &buf, std::map<int, TPZAutoPointer<T> > &vec)
	{
		int nc = vec.size(),one = -1;
		buf.Write(&nc,1);
		typedef typename std::map<int, TPZAutoPointer<T> >::iterator vecit_type;
		vecit_type vecit;
		for(vecit=vec.begin(); vecit!= vec.end(); vecit++)
		{
			int id = vecit->first;
			buf.Write(&id,1);
			if(vecit->second)
			{
				vecit->second->Write(buf,1);
			} else {
				buf.Write(&one,1);
			}
		}
	}
	
	template<class T, int EXP>
	static void WriteObjectPointers(TPZStream &buf, TPZChunkVector<T *,EXP> &vec)
	{
		int c,m1=-1,nc = vec.NElements();
		buf.Write(&nc,1);
		for(c=0; c<nc; c++)
		{
			T *ptr = vec[c];
			if(ptr) ptr->Write(buf);
			else buf.Write(&m1,1);
		}
	}
	
	template<class T, int EXP>
	static void WriteObjectPointers(TPZStream &buf, TPZAdmChunkVector<T *,EXP> &vec)
	{
		int c,m1=-1,nc = vec.NElements();
		buf.Write(&nc,1);
		for(c=0; c<nc; c++)
		{
			T *ptr = vec[c];
			if(ptr) ptr->Write(buf,1);
			else buf.Write(&m1,1);
		}
		buf.Write(&vec.fCompactScheme,1);
		WriteObjects(buf,vec.fFree);
		WriteObjects(buf,vec.fNFree);
	}
	
	template<class T>
	static void ReadObjects(TPZStream &buf, TPZVec<T> &vec, void *context)
	{
		int c,nc;
		buf.Read(&nc,1);
		vec.Resize(nc);
		for(c=0; c<nc; c++)
		{
			vec[c].Read(buf,context);
		}
	}
	
	static void ReadObjects(TPZStream &buf, TPZVec<int> &vec)
	{
		int nc;
		buf.Read(&nc,1);
		vec.Resize(nc);
		if(nc) buf.Read(&vec[0],nc);
	}
	
	template<class T>
	static void ReadObjects(TPZStream &buf, std::vector<T> &vec, void *context)
	{
		int c,nc;
		buf.Read(&nc,1);
		vec.resize(nc);
		for(c=0; c<nc; c++)
		{
			vec[c].Read(buf,context);
		}
	}
	
	static void ReadObjects(TPZStream &buf, std::vector<int> &vec)
	{
		int nc;
		buf.Read(&nc,1);
		vec.resize(nc);
		if(nc) buf.Read(&vec[0],nc);
	}
	
	static void ReadObjects(TPZStream &buf, std::vector<REAL> &vec)
	{
		int nc;
		buf.Read(&nc,1);
		vec.resize(nc);
		if(nc) buf.Read(&vec[0],nc);
	}
	
	static void ReadObjects(TPZStream &buf, TPZVec<REAL> &vec)
	{
		int nc;
		buf.Read(&nc,1);
		vec.Resize(nc);
		if(nc) buf.Read(&vec[0],nc);
	}
	
	template<int N>
	static void ReadObjects(TPZStream &buf, TPZManVector<REAL,N> &vec)
	{
		int nc;
		buf.Read(&nc,1);
		vec.Resize(nc);
		if(nc) buf.Read(&vec[0],nc);
	}
	
	template<class T, int EXP>
	static void ReadObjects(TPZStream &buf, TPZChunkVector<T,EXP> &vec, void *context)
	{
		int c,nc;
		buf.Read(&nc,1);
		vec.Resize(nc);
		for(c=0; c<nc; c++)
		{
			vec[c].Read(buf,context);
		}
	}
	
	template<class T, int EXP>
	static void ReadObjects(TPZStream &buf, TPZAdmChunkVector<T,EXP> &vec, void *context)
	{
		int c,nc;
		buf.Read(&nc,1);
		vec.Resize(nc);
		for(c=0; c<nc; c++) vec[c].Read(buf,context);
		buf.Read(&vec.fCompactScheme,1);
		ReadObjects(buf,vec.fFree);
		ReadObjects(buf,vec.fNFree);
	}
	
    template<class TVar>
	static void ReadObjects(TPZStream &buf, std::map<TVar,TVar> &vec)
	{
		
		TPZManVector<TVar> cp;
		ReadObjects(buf, cp);
		int sz = cp.NElements();
		int i;
		for (i=0; i<sz; i+=2) {
			vec[cp[i]] = cp[i+1];
		}
	}
	
	template<class T>
	static void ReadObjectPointers(TPZStream &buf, TPZVec<T *> &vec, void *context)
	{
		int c,nc;
		buf.Read(&nc,1);
		vec.Resize(nc);
		for(c=0; c<nc; c++)
		{
			vec[c] = dynamic_cast<T *>(Restore(buf,context));
		}
	}
	
	template<class T>
    static void ReadObjectPointers(TPZStream &buf, std::map<int, TPZAutoPointer<T> > &vec, void *context)
	{
		int c,nc;
		buf.Read(&nc,1);
		for(c=0; c<nc; c++)
		{
			int id;
			buf.Read(&id,1);
			vec[id] = TPZAutoPointer<T>(dynamic_cast<T *>(Restore(buf,context)));
		}
	}
	
	template<class T, int EXP>
	static void ReadObjectPointers(TPZStream &buf, TPZChunkVector<T *,EXP> &vec, void *context)
	{
		int c,nc;
		buf.Read(&nc,1);
		vec.Resize(nc);
		for(c=0; c<nc; c++)
		{
			vec[c] = dynamic_cast<T *>(Restore(buf,context));
		}
	}
	
	template<class T, int EXP>
	static void ReadObjectPointers(TPZStream &buf, TPZAdmChunkVector<T *,EXP> &vec, void *context)
	{
		int c,nc;
		buf.Read(&nc,1);
		vec.Resize(nc);
		for(c=0; c<nc; c++) vec[c] = dynamic_cast<T *>(Restore(buf,context));
		buf.Read(&vec.fCompactScheme,1);
		ReadObjects(buf,vec.fFree);
		ReadObjects(buf,vec.fNFree);
	}
	
	static void WriteObjects(TPZStream &buf, TPZVec<REAL> &vec)
	{
		int nel = vec.NElements();
		buf.Write(&nel,1);
		if(nel) buf.Write(&vec[0],vec.NElements());
	}
	
	static void WriteObjects(TPZStream &buf, std::vector<REAL> &vec)
	{
		int nel = vec.size();
		buf.Write(&nel,1);
		if(nel) buf.Write(&vec[0],vec.size());
	}
	
#ifndef ELLIPS
	static void WriteObjects(TPZStream &buf, TPZVec<TPZFlopCounter> &vec)
	{
		int nel = vec.NElements();
		buf.Write(&nel,1);
		if(nel) buf.Write(&vec[0],vec.NElements());
	}
	
	static void WriteObjects(TPZStream &buf, std::vector<TPZFlopCounter> &vec)
	{
		int nel = vec.size();
		buf.Write(&nel,1);
		if(nel) buf.Write(&vec[0],vec.size());
	}
#endif
	
	static void WriteObjects(TPZStream &buf, TPZVec<int> &vec)
	{
		int nel = vec.NElements();
		buf.Write(&nel,1);
		if(nel) buf.Write(&vec[0],vec.NElements());
	}
	
	static void WriteObjects(TPZStream &buf, std::vector<int> &vec)
	{
		int nel = vec.size();
		buf.Write(&nel,1);
		if(nel) buf.Write(&vec[0],vec.size());
	}
	
	static void WriteObjects(TPZStream &buf, TPZVec<char> &vec)
	{
		int nel = vec.NElements();
		buf.Write(&nel,1);
		if(nel) buf.Write(&vec[0],vec.NElements());
	}
	
	static void WriteObjects(TPZStream &buf, std::vector<char> &vec)
	{
		int nel = vec.size();
		buf.Write(&nel,1);
		if(nel) buf.Write(&vec[0],vec.size());
	}
	
	static void Register(int classid, TPZRestore_t fun);
	
	static TPZSaveable *Restore(TPZStream &buf, void *context);
	
};

/** @brief Restores object from Map, classid is in buf */
template<class T>
TPZSaveable *Restore(TPZStream &buf, void *context) {
	T *ptr = new T;
	ptr->Read(buf,context);
	return ptr;
}

#ifndef ELLIPS
/**
 * @brief Implements an interface to register a class id and a restore function. \ref save "Persistency"
 */
/**
 * A declaration of the type "template class<classname, classid> put in .cpp file does the trick \n
 * The static object which is "automatically" created calls the proper interface of the TPZSaveable class
 */
template<class T, int N>
class TPZRestoreClass {
public:
	/** @brief Constructor */
	TPZRestoreClass()
	{
#ifdef DEBUG 
		std::string func_name = __PRETTY_FUNCTION__;
#ifndef WIN32
		std::cout << func_name << std::endl;
#endif
#endif
		TPZSaveable::Register(N,Restore);
	}
public:
	/** @brief Restores object from Map based in classid into the buf */
	static TPZSaveable *Restore(TPZStream &buf, void *context) {
		T *ptr = new T;
		ptr->Read(buf,context);
		return ptr;
	}
private:
	
	static TPZRestoreClass gRestoreObject;
};

template<>
inline TPZSaveable *TPZRestoreClass<TPZSaveable,-1>::Restore(TPZStream &buf, void *context)
{
	return 0;
}
template<class T, int N>
TPZRestoreClass<T,N> TPZRestoreClass<T,N>::gRestoreObject;

/// To restore object
template<>
inline TPZSaveable *Restore<TPZSaveable>(TPZStream &buf, void *context) {
	return 0;
}

/** @} */

#endif //ellips

#endif //PZSAVEH


