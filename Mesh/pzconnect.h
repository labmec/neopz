/**
 * @file
 * @brief Contains declaration of TPZConnect class which represents a set of shape functions associated with a computational element
 */
//$Id: pzconnect.h,v 1.19 2010-08-25 03:05:06 phil Exp $
//HEADER FILE FOR CLASS NODE

#ifndef  PZCONNECTH
#define  PZCONNECTH

#include "pzfmatrix.h"
#include "pzstack.h"
#include <iostream>
#include <set>


class TPZBndCond;
class TPZCompMesh;
template<class TVar>
class TPZBlock;
class TPZStream;



/** 
 * @brief Represents a set of shape functions associated with a computational element/side. \ref CompElement "Computational Element"
 * @ingroup CompElement
 */
/**
 * This class keeps track of information associated with an element/side such as order of interpolation
 * sequence number in the vector of blocks of equations \n
 * Objects of this class also contain the information necessary for constraints between shapefunctions
 */
class TPZConnect {
public:
    enum EConnectType {ENone = 0, EPressure = 1, ECondensed = 2};
	/** @brief Node block number */
	int		fSequenceNumber;
	
	/** @brief Number of element connected */
	int		fNElConnected;
	    
    /** @brief Flag to signal the type of connect 
     * ECondensed : the connect will be condensed before assembly 
     * EPressure : a pressure connect which should be assembled after the flux connects
     * 
     */
    union
    {
        int fFlags;
        struct
        {
            /** @brief Interpolation order of the associated shape functions */
            unsigned char fOrder;
            /** @brief Number of state variables associated with each shape function */
            unsigned char fNState;
            /** @brief Number of shape functions associated with the connect */
            unsigned char fNShape;
            /** @brief Whether the equations associated with the connect should be condensed */
            bool fIsCondensed;
            /** @brief Whether the connnect is associated with a lagrange multiplier */
            bool fIsPressure;
        } fCompose;
    };
	
public:
	/** @brief Structure to reference dependency */
	struct TPZDepend
	{
		int			fDepConnectIndex;
		TPZFNMatrix<50>	fDepMatrix;
		TPZDepend		*fNext;
		
		TPZDepend(int DepConnectIndex,TPZFMatrix<REAL> &depmat,int ipos,int jpos, int isize, int jsize);
		
		TPZDepend(const TPZDepend &copy);
		TPZDepend(int connectindex);
		
		~TPZDepend();
		TPZDepend *HasDepend(int DepConnectIndex);
		TPZDepend *RemoveDepend(TPZDepend *Ptr);
		void Write(TPZStream &buf);
		void Read(TPZStream &buf);
		
		/**
		 * @brief Copy a depend data structure to a clone depend in a clone mesh
		 * @param orig original depend to be copied
		 * @param gl2lcIdx global to local indexes map
		 */
		void CopyFrom(TPZDepend *orig , std::map<int,int>& gl2lcIdx);
	};
	
private:
	TPZDepend *fDependList;
	
public:
	/** @brief Default constructor */
	TPZConnect();
	/** @brief Default destructor */
	~TPZConnect();
	
	void operator=(const TPZConnect &con);
	
    /** @brief Reset the data of the connect */
    void Reset()
    {
        SetSequenceNumber(-1);
        SetNState(0);
        SetOrder(0);
        SetNShape(0);
        ResetElConnected();
        SetCondensed(false);
        SetPressure(false);
        if (fDependList) {
            DebugStop();
            delete fDependList;
            fDependList = 0;
        }
    }
	/**
	 * @brief Number of degrees of freedom associated with the object
	 * @note In order to compute NDof, the Connect object needs to know the mesh
	 */
	int NDof(TPZCompMesh &mesh);
    
	/**
	 * @brief Number of degrees of freedom associated with the object
	 * @note In order to compute NDof, the Connect object needs to know the mesh
	 */
	int NDof() const
    {
        return fCompose.fNShape*fCompose.fNState;
    }
    
    /**
     * @brief Number of state variables associated with the connect
     */
    unsigned char NState() const
    {
        return fCompose.fNState;
    }
    
    unsigned char NShape() const
    {
        return fCompose.fNShape;
    }

	/** @brief Returns the Sequence number of the connect object */
	/** If the \f$ sequence number == -1 \f$ this means that the node is unused */
	int SequenceNumber() const;
	
	/** @brief Set the sequence number for the global system of equations of the connect object*/
	/** If the argument \f$i==-1\f$ this means that the node is out of use*/
	void SetSequenceNumber(int i) {fSequenceNumber = i;}
	
	/** @brief Set the order of the shapefunction associated with the connect */
	void SetOrder(int order) {
#ifdef DEBUG
        if(order < 0 || order > 255)
        {
            DebugStop();
        }
#endif
		fCompose.fOrder = order;
	}
    
    /** @brief Set the number of state variables */
    void SetNState(int nstate)
    {
        if(nstate<0 || nstate > 255)
        {
            DebugStop();
        }
        fCompose.fNState = nstate;
    }
    /** @brief Set the number of shape functions associated with the connect */
    void SetNShape(int nshape)
    {
        if(nshape < 0 || nshape > 255)
        {
            DebugStop();
        }
        fCompose.fNShape = nshape;
    }
	
	/** @brief Access function to return the order associated with the connect */
	unsigned char Order() const {
		return fCompose.fOrder;
	}
    
    /** @brief Access method to return the indication whether the connect is condensed or not
     */
    bool IsCondensed() const
    {
        return fCompose.fIsCondensed;
    }
    
    /** @brief Access method to return the indication whether the connect is associated with a pressure lagrange multiplier
     */
    bool IsPressure() const
    {
        return fCompose.fIsPressure;
    }
    
    /** @brief Set the connect as a pressure connect or not */
    void SetPressure(bool flag)
    {
        fCompose.fIsPressure = flag;
    }

    /** @brief Set the connect as a condensed connect or not */
    void SetCondensed(bool flag)
    {
        fCompose.fIsCondensed = flag;
    }

	/** @brief Print the information for the connect element */
	/**
	 * The mesh argument allows the object to identify the number of variables
	 * associated with it and the solution
	 */
	void Print(const TPZCompMesh &mesh, std::ostream & out = std::cout);
	
	/** @brief Also print the center point of the side associated to the connect */
	void Print(TPZCompMesh &mesh, TPZVec<REAL> &cp, std::ostream & out = std::cout);
	
	/** @brief Initialize with zero fNElConnected */
	void ResetElConnected() { fNElConnected = 0; }
	/** @brief Increment fNElConnected */
	void IncrementElConnected() { fNElConnected++; }
	/** @brief Decrement fNElConnected */
	void DecrementElConnected() { fNElConnected--; }
	/** @brief Returns fNElConnected */
	int NElConnected() const { return fNElConnected; }
	
	/**
	 * @brief Add dependency between connects
	 * @param myindex [in] index of this connect
	 * @param dependindex [in] index of the connect this will depend upon
	 * @param depmat [in] dependency matrix which defines the relation between the connects
	 * @param ipos, jpos, isize, jsize are parameters which define the submatrix within depmat which is to be used
	 */
	void AddDependency(int myindex, int dependindex,TPZFMatrix<REAL> &depmat,int ipos,int jpos, int isize, int jsize);
	
	/**
	 * @brief Remove dependency between connects if exist
	 * @param myindex [in] index of this connect
	 * @param dependindex [in] index of the connect was exist dependency
	 */
	void RemoveDepend(int myindex, int dependindex);
	
	/** @brief Deletes all dependency information */
	void RemoveDepend();

	int DependencyDepth(TPZCompMesh &mesh);
	/** @brief Returns whether exist dependecy information */
	int HasDependency()const { return fDependList != 0; }
	
	int CheckDependency(int nshape, TPZCompMesh *mesh, int nstate);
	
	TPZDepend *FirstDepend() { return fDependList; }
	
	/** @brief Adds itself and the connects from which it depends to the list*/
	void AddToList(int myindex, TPZCompMesh &mesh, TPZStack<int> &connectlist);
	
	/** @brief Adds itself and the connects from which it depends to the list*/
	void AddToList(int myindex, TPZCompMesh &mesh, std::set<int> &connectlist);
	
	void SetDependenceOrder(int myindex, TPZCompMesh &mesh, int CurrentOrder,TPZVec<int> &connectlist,TPZVec<int> &DependenceOrder);
	
	void ExpandShape(int cind, TPZVec<int> &connectlist, TPZVec<int> &blocksize, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
	
	/** @brief Saves the element data to a stream */
	void Write(TPZStream &buf, int withclassid);
	
	/** @brief Reads the element data from a stream */
	void Read(TPZStream &buf, void *context);
	
	/**
	 * @brief Copy a connect data structure from an original connect to a new connect mapping their indexes
	 * @param orig original connect to be copied
	 * @param gl2lcIdx global to local indexes map
	 */
	void CopyFrom(TPZConnect &orig,std::map<int,int> & gl2lcIdx);
	
	/**
	 * @brief Builds the list of all connectivities related to ConnectIndex including the
	 * connects pointed to by dependent connects
	 * @note Note : this method does not reset the stack to zero. The calling
	 * method should do this
	 * @param index [in] index of the current connect
	 * @param indepconnectlist [out] set which contains the indices of independent connects
	 * @param depconnectlist [out] set which contains the indices of dependent connects
	 * @param mesh [in]
	 */
	void BuildConnectList(int index, std::set<int> &indepconnectlist, std::set<int> &depconnectlist, TPZCompMesh &mesh);
	
	/**
	 * @brief Builds the list of all connectivities related to ConnectIndex including the
	 * connects pointed to by dependent connects
	 * @note Note : this method does not reset the stack to zero. The calling
	 * method should do this
	 * @param connectlist [out] stack to receive the list
	 * @param ConnectIndex [in]
	 * @param mesh [in]
	 */
	static void BuildConnectList(TPZStack<int> &connectlist, TPZVec<int> &ConnectIndex, TPZCompMesh &mesh);
	
	/**
	 * @brief Builds the list of all connectivities related to ConnectIndex including the
	 * connects pointed to by dependent connects
	 * @note Note : this method does not reset the set to zero. The calling
	 * method should do this
	 * @param connectlist [out] set to receive the list
	 * @param additional [in] connects which should be added along with their dependencies
	 * @param mesh [in]
	 */
	static void BuildConnectList(std::set<int> &connectlist, std::set<int> &additional, TPZCompMesh &mesh);
	
	/**
	 * @brief This method builds the vector DependenceOrder which indicates in which
	 * order constrained nodes need to be processed
	 */
	/** connectlist need to be computed by BuildConnectList */
	static void BuildDependencyOrder(TPZVec<int> &connectlist, TPZVec<int> &DependenceOrder, TPZCompMesh &mesh);
	
};

/** @ingroup CompElement */
/** @brief Overload operator << to write node connect data */
inline std::ostream & operator<<(std::ostream &out,TPZConnect &con)
{
	out << "seq num: " << con.SequenceNumber()
	<< " nel con: " << con.NElConnected()
	<< " order: " << con.Order()
    << " is condensed " << con.IsCondensed()
	<< " hasdepend: " << con.HasDependency();
	return out;
}

#endif

