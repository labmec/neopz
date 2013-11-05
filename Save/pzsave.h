/**
 * @file
 * @brief Contains declaration of the TPZSaveable class which defines the interface to save and restore objects from TPZStream objects.
 */

#ifndef PZSAVEH
#define PZSAVEH

#include <map>
#include <vector>
#include <set>


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
/** @brief Identifier as saveable object */
const int TPZSAVEABLEID = -1;

/** @brief Typedef of TPZRestore_t */
typedef TPZSaveable *(*TPZRestore_t)(TPZStream &,void *);


/** 
 * The SAVEABLE NOTES are used to debug read and write operations on
 * TPZSaveable objects. Do not use it for other purposes.
 * To enable it, define the DEBUG_SAVEABLE macro.
 */
//#define DEBUG_SAVEABLE

#ifdef DEBUG_SAVEABLE
#define SAVEABLE_STR_NOTE(buf,str) { std::string msg(str); buf.Write(&msg,1); }
#define SAVEABLE_SKIP_NOTE(buf) { std::string str; buf.Read(&str,1); }
#else
#define SAVEABLE_STR_NOTE(buf,str)
#define SAVEABLE_SKIP_NOTE(buf)
#endif

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
	
	/** @brief Writes this object to the TPZStream buffer. Include the classid if withclassid = true */
	virtual void Write(TPZStream &buf, int withclassid) const;
	
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
	static void WriteObjects(TPZStream &buf, const TPZVec<T> &vec)
	{
		long c,nc = vec.NElements();
		buf.Write(&nc,1);
		for(c=0; c<nc; c++)
            vec[c].Write(buf,0);
	}
	
	template<class T>
	static void WriteObjects(TPZStream &buf, const std::vector<T> &vec)
	{
		int c,nc = vec.size();
		buf.Write(&nc,1);
		for(c=0; c<nc; c++) vec[c].Write(buf,0);
	}
	
	template<class T, int EXP>
	static void WriteObjects(TPZStream &buf, const TPZChunkVector<T,EXP> &vec)
	{
		long c,nc = vec.NElements();
		buf.Write(&nc,1);
		for(c=0; c<nc; c++) vec[c].Write(buf,0);
	}
	
	template<class T, int EXP>
	static void WriteObjects(TPZStream &buf, TPZAdmChunkVector<T,EXP> &vec)
	{
		long c,nc = vec.NElements();
		buf.Write(&nc,1);
		for(c=0; c<nc; c++) vec[c].Write(buf,0);
		buf.Write(&vec.fCompactScheme,1);
		WriteObjects(buf,vec.fFree);
		WriteObjects(buf,vec.fNFree);
	}
	
	static void WriteObjects(TPZStream &buf, const std::map<int,int> &vec)
	{
		long sz = vec.size();
		TPZManVector<int> cp(sz*2);
		long count = 0;
        std::map<int,int>::const_iterator it;
		for (it=vec.begin(); it!=vec.end(); it++) {
			cp[count++] = it->first;
			cp[count++] = it->second;
		}
		WriteObjects(buf, cp);
	}
	static void WriteObjects(TPZStream &buf, const std::map<long,long> &vec)
	{
		long sz = vec.size();
		TPZManVector<long> cp(sz*2);
		long count = 0;
        std::map<long,long>::const_iterator it;
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
		long c,nc = vec.NElements();
        int one = -1;
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
	
	template<class T>
	static void WriteObjectPointers(TPZStream &buf, std::map<int, T* > &vec)
	{
		int nc = vec.size(),one = -1;
		buf.Write(&nc,1);
		typedef typename std::map<int, T* >::iterator vecit_type;
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
    template<class T>
    static void WriteObjectPointers(TPZStream &buf, std::set<T* > &vec)
	{
		int nel = vec.size();
		buf.Write(&nel,1);
        typedef typename  std::set<T* >::iterator vec_it;
        vec_it it;
        while (it != vec.end()) {
            it-> Write(buf,1);
            it++;
        }
	}
	
	template<class T, int EXP>
	static void WriteObjectPointers(TPZStream &buf, TPZChunkVector<T *,EXP> &vec)
	{
		long c,nc = vec.NElements();
        int m1 = -1;
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
		long c,nc = vec.NElements();
        int m1 = -1;
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

	/**
	 * Write for chunk vectors with basic elements as float, double, long double, std::complex<...> .
	 */
	template<class T, int EXP>
	static void WriteObjects(TPZStream &buf,const TPZAdmChunkVector<T,EXP> &vec,bool basic)
	{
		if(!basic) {
			DebugStop();
			return;
		}
		long c,nc = vec.NElements();
		buf.Write(&nc,1);
		for(c=0; c<nc; c++) buf.Write(&vec[c],1);
		buf.Write(&vec.fCompactScheme,1);
		WriteObjects(buf,vec.fFree,true);
		WriteObjects(buf,vec.fNFree,true);
	}
	template<class T>
	static void WriteObjects(TPZStream &buf, const TPZVec<T> &vec,bool basic)
	{
		if(!basic) {
			DebugStop();
			return;
		}
		long c,nc = vec.NElements();
		buf.Write(&nc,1);
		for(c=0; c<nc; c++)
            buf.Write(&vec[c],1);
	}

	/**
	 * @brief Methods to read objects or pointer for objects.
	 */
	template<class T>
	static void ReadObjects(TPZStream &buf, TPZVec<T> &vec, void *context)
	{
		long c,nc;
		buf.Read(&nc,1);
		vec.Resize(nc);
		for(c=0; c<nc; c++)
		{
			vec[c].Read(buf,context);
		}
	}
    
    static void ReadObjects(TPZStream &buf, std::set<int> &vec)
	{
	  int nel;
	  buf.Read(&nel,1);
        for (int i=0; i<nel; i++)
        {
            int val;
            buf.Read(&val);
            vec.insert(val);
        }
	}
	

	static void ReadObjects(TPZStream &buf, TPZVec<int> &vec)
	{
		long nc;
		buf.Read(&nc,1);
		vec.Resize(nc);
		if(nc) buf.Read(&vec[0],nc);
	}
	static void ReadObjects(TPZStream &buf, TPZVec<long> &vec)
	{
		long nc;
		buf.Read(&nc,1);
		vec.Resize(nc);
		if(nc) buf.Read(&vec[0],nc);
	}
	
	static void ReadObjects(TPZStream &buf, TPZVec<TPZFlopCounter> &vec)
	{
		long nc;
		buf.Read(&nc,1);
		vec.Resize(nc);
        TPZVec<REAL> temp(nc);
		if(nc) buf.Read(&temp[0],nc);
        for (long ic=0; ic<nc; ic++) {
            vec[ic] = temp[ic];
        }
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
	static void ReadObjects(TPZStream &buf, std::vector<long> &vec)
	{
		long nc;
		buf.Read(&nc,1);
		vec.resize(nc);
		if(nc) buf.Read(&vec[0],nc);
	}
	
	static void ReadObjects(TPZStream &buf, std::vector<float> &vec)
	{
		int nc;
		buf.Read(&nc,1);
		vec.resize(nc);
		if(nc) buf.Read(&vec[0],nc);
	}
	
	static void ReadObjects(TPZStream &buf, TPZVec<float> &vec)
	{
		long nc;
		buf.Read(&nc,1);
		vec.Resize(nc);
		if(nc) buf.Read(&vec[0],nc);
	}
	
	static void ReadObjects(TPZStream &buf, std::vector<double> &vec)
	{
		int nc;
		buf.Read(&nc,1);
		vec.resize(nc);
		if(nc) buf.Read(&vec[0],nc);
	}
	
	static void ReadObjects(TPZStream &buf, TPZVec<double> &vec)
	{
		long nc;
		buf.Read(&nc,1);
		vec.Resize(nc);
		if(nc) buf.Read(&vec[0],nc);
	}
	
	static void ReadObjects(TPZStream &buf, std::vector<long double> &vec)
	{
		int nc;
		buf.Read(&nc,1);
		vec.resize(nc);
		if(nc) buf.Read(&vec[0],nc);
	}
	
	static void ReadObjects(TPZStream &buf, TPZVec<long double> &vec)
	{
		long nc;
		buf.Read(&nc,1);
		vec.Resize(nc);
		if(nc) buf.Read(&vec[0],nc);
	}

	template<class T, int EXP>
	static void ReadObjects(TPZStream &buf, TPZAdmChunkVector<T,EXP> &vec)
	{
		long c,nc;
		buf.Read(&nc,1);
		vec.Resize(nc);
		for(c=0; c<nc; c++) buf.Read(&vec[c],1);
		buf.Read(&vec.fCompactScheme,1);
		ReadObjects(buf,vec.fFree);
		ReadObjects(buf,vec.fNFree);
	}

	static void ReadObjects(TPZStream &buf, std::vector<std::complex<float> > &vec)
	{
		int nc;
		buf.Read(&nc,1);
		vec.resize(nc);
		if(nc) buf.Read(&vec[0],nc);
	}
	
	static void ReadObjects(TPZStream &buf, TPZVec<std::complex<float> > &vec)
	{
		long nc;
		buf.Read(&nc,1);
		vec.Resize(nc);
		if(nc) buf.Read(&vec[0],nc);
	}

	static void ReadObjects(TPZStream &buf, std::vector<std::complex<double> > &vec)
	{
		int nc;
		buf.Read(&nc,1);
		vec.resize(nc);
		if(nc) buf.Read(&vec[0],nc);
	}
	
	static void ReadObjects(TPZStream &buf, TPZVec<std::complex<double> > &vec)
	{
		long nc;
		buf.Read(&nc,1);
		vec.Resize(nc);
		if(nc) buf.Read(&vec[0],nc);
	}
	
	static void ReadObjects(TPZStream &buf, std::vector<std::complex<long double> > &vec)
	{
		int nc;
		buf.Read(&nc,1);
		vec.resize(nc);
		if(nc) buf.Read(&vec[0],nc);
	}
	
	static void ReadObjects(TPZStream &buf, TPZVec<std::complex<long double> > &vec)
	{
		long nc;
		buf.Read(&nc,1);
		vec.Resize(nc);
		if(nc) buf.Read(&vec[0],nc);
	}

	template<int N>
	static void ReadObjects(TPZStream &buf, TPZManVector<REAL,N> &vec)
	{
		long nc;
		buf.Read(&nc,1);
		vec.Resize(nc);
		if(nc) buf.Read(&vec[0],nc);
	}
	
	template<class T, int EXP>
	static void ReadObjects(TPZStream &buf, TPZChunkVector<T,EXP> &vec, void *context)
	{
		long c,nc;
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
		long c,nc;
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
		long c,nc;
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
	
	template<class T>
    static void ReadObjectPointers(TPZStream &buf, std::map<int, T* > &vec, void *context)
	{
		int c,nc;
		buf.Read(&nc,1);
		for(c=0; c<nc; c++)
		{
			int id;
			buf.Read(&id,1);
			vec[id] = (dynamic_cast<T *>(Restore(buf,context)));
		}
	}
	

	template<class T, int EXP>
	static void ReadObjectPointers(TPZStream &buf, TPZChunkVector<T *,EXP> &vec, void *context)
	{
		long c,nc;
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
		long c,nc;
		buf.Read(&nc,1);
		vec.Resize(nc);
		for(c=0; c<nc; c++) 
        {
            vec[c] = dynamic_cast<T *>(Restore(buf,context));
        }
		buf.Read(&vec.fCompactScheme,1);
		ReadObjects(buf,vec.fFree);
		ReadObjects(buf,vec.fNFree);
	}
    
	
	
	static void WriteObjects(TPZStream &buf, const TPZVec<float> &vec)
	{
		long nel = vec.NElements();
		buf.Write(&nel,1);
		if(nel) buf.Write(&vec[0],vec.NElements());
	}
	
	static void WriteObjects(TPZStream &buf, std::set<int> &vec)
	{
		int nel = vec.size();
		buf.Write(&nel,1);
        std::set<int>::iterator it = vec.begin();
        while (it != vec.end()) {
            int val = *it;
            buf.Write(&val);
            it++;
        }
	}
	
	static void WriteObjects(TPZStream &buf, const std::vector<float> &vec)
	{
		int nel = vec.size();
		buf.Write(&nel,1);
		if(nel) buf.Write(&vec[0],vec.size());
	}
	
	static void WriteObjects(TPZStream &buf, const TPZVec<double> &vec)
	{
		long nel = vec.NElements();
		buf.Write(&nel,1);
		if(nel) buf.Write(&vec[0],vec.NElements());
	}
	
	static void WriteObjects(TPZStream &buf, const TPZVec<TPZFlopCounter> &vec)
	{
		long nel = vec.NElements();
		buf.Write(&nel,1);
        TPZVec<REAL> temp(nel);
        for (int iel = 0; iel<nel; iel++) {
            temp[iel] = vec[iel];
        }
		if(nel) buf.Write(&temp[0],vec.NElements());
	}
	
	static void WriteObjects(TPZStream &buf, const std::vector<double> &vec)
	{
		int nel = vec.size();
		buf.Write(&nel,1);
		if(nel) buf.Write(&vec[0],vec.size());
	}
	
	static void WriteObjects(TPZStream &buf, const TPZVec<std::complex<double> > &vec)
	{
		long nel = vec.NElements();
		buf.Write(&nel,1);
		if(nel) buf.Write(&vec[0],vec.NElements());
	}
	
	static void WriteObjects(TPZStream &buf, const std::vector<std::complex<double> > &vec)
	{
		int nel = vec.size();
		buf.Write(&nel,1);
		if(nel) buf.Write(&vec[0],vec.size());
	}
	
	static void WriteObjects(TPZStream &buf, const TPZVec<long double> &vec)
	{
		long nel = vec.NElements();
		buf.Write(&nel,1);
		if(nel) buf.Write(&vec[0],vec.NElements());
	}
	
	static void WriteObjects(TPZStream &buf, const std::vector<long double> &vec)
	{
		int nel = vec.size();
		buf.Write(&nel,1);
		if(nel) buf.Write(&vec[0],vec.size());
	}
	
	static void WriteObjects(TPZStream &buf, const TPZVec<std::complex<long double> > &vec)
	{
		long nel = vec.NElements();
		buf.Write(&nel,1);
		if(nel) buf.Write(&vec[0],vec.NElements());
	}
	
	static void WriteObjects(TPZStream &buf, const std::vector<std::complex<long double> > &vec)
	{
		int nel = vec.size();
		buf.Write(&nel,1);
		if(nel) buf.Write(&vec[0],vec.size());
	}

	static void WriteObjects(TPZStream &buf, const TPZVec<std::complex<float> > &vec)
	{
		long nel = vec.NElements();
		buf.Write(&nel,1);
		if(nel) buf.Write(&vec[0],vec.NElements());
	}
	
	static void WriteObjects(TPZStream &buf, const std::vector<std::complex<float> > &vec)
	{
		int nel = vec.size();
		buf.Write(&nel,1);
		if(nel) buf.Write(&vec[0],vec.size());
	}

#ifndef ELLIPS
	static void WriteObjects(TPZStream &buf, TPZVec<TPZFlopCounter> &vec)
	{
		long nel = vec.NElements();
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
		long nel = vec.NElements();
		buf.Write(&nel,1);
		if(nel) buf.Write(&vec[0],vec.NElements());
	}
	static void WriteObjects(TPZStream &buf, TPZVec<long> &vec)
	{
		long nel = vec.NElements();
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
		long nel = vec.NElements();
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


