/**
 * @file
 * @brief Contains declaration of TPZConnect class which represents a set of shape functions associated with a computational element.
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
class TPZBlock;
class TPZStream;



/** 
 * @brief Represents a set of shape functions associated with a computational element/side. \ref interpolation "Aproximation space"
 * @ingroup interpolation
 */
/**
 * This class keeps track of information associated with an element/side such as order of interpolation
 * sequence number in the vector of blocks of equations \n
 * Objects of this class also contain the information necessary for constraints between shapefunctions
 */
class TPZConnect {
	/** @brief Node block number */
	int		fSequenceNumber;
	
	/** @brief Number of element connected */
	int		fNElConnected;
	
	/** @brief Interpolation order of the associated shape functions */
	int fOrder;
	
public:
	/** @brief Structure to reference dependency */
	struct TPZDepend
	{
		int			fDepConnectIndex;
		TPZFNMatrix<50>	fDepMatrix;
		TPZDepend		*fNext;
		
		TPZDepend(int DepConnectIndex,TPZFMatrix &depmat,int ipos,int jpos, int isize, int jsize);
		
		TPZDepend(const TPZDepend &copy);
		TPZDepend(int connectindex);
		
		~TPZDepend();
		TPZDepend *HasDepend(int DepConnectIndex);
		TPZDepend *RemoveDepend(TPZDepend *Ptr);
		void Write(TPZStream &buf);
		void Read(TPZStream &buf);
		
		/**
		 * @brief Copy a depend data structure to a clone depend in
		 * a clone mesh
		 * @param orig original depend to be copied
		 * @param gl2lcIdx global to local indexes map
		 */
		void CopyFrom(TPZDepend *orig , std::map<int,int>& gl2lcIdx);
	};
	
private:
	TPZDepend *fDependList;
	
public:
	/** @brief Constructor*/
	TPZConnect();
	/** @brief Destructor*/
	~TPZConnect();
	
	void operator=(const TPZConnect &con);
	
	/** @brief Number of degrees of freedom associated with the object
	 
	 It needs the mesh object to access this information*/
	// @note Philippe In order to compute NDof, the Connect object needs to know the mesh
	int NDof(TPZCompMesh &mesh); /*Cedric 30/09/98 -  14:05*/
	
	
	/** @brief Returns the Sequence number of the connect object.
     If the sequence number == -1 this means that the node is unused*/
	int SequenceNumber() const;
	
	/** @brief Set the sequence number for the global system of equations of the connect
     object. If the argument i==-1 this means that the node is out of use*/
	void SetSequenceNumber(int i) {fSequenceNumber = i;}
	
	/**
	 * @brief Set the order of the shapefunction associated with the connect
	 */
	void SetOrder(int order) {
		fOrder = order;
	}
	
	/**
	 * @brief Access function to return the order associated with the connect
	 */
	int Order() const {
		return fOrder;
	}
	
	/** @brief Print the information for the connect element.
	 *
	 * The mesh argument allows the object to identify the number of variables
	 * associated with it and the solution
	 */
	void Print(const TPZCompMesh &mesh, std::ostream & out = std::cout);
	
	/**
	 * @brief Also print the center point of the side associated to the connect
	 */
	void Print(TPZCompMesh &mesh, TPZVec<REAL> &cp, std::ostream & out = std::cout);
	
	/** @brief Initialize with zero fNElConnected*/
	void ResetElConnected() { fNElConnected = 0; }
	/** @brief Increment fNElConnected*/
	void IncrementElConnected() { fNElConnected++; }
	/** @brief Decrement fNElConnected*/
	void DecrementElConnected() { fNElConnected--; }
	/** @brief Returns fNElConnected*/
	int NElConnected() const { return fNElConnected; }
	
	/** 
	 * @param myindex [in] index of this connect
	 * @param dependindex [in] index of the connect this will depend upon
	 * @param depmat [in] dependency matrix which defines the relation between the connects
	 * @param ipos, jpos, isize, jsize are parameters which define the submatrix within depmat which is to be used
	 */
	void AddDependency(int myindex, int dependindex,TPZFMatrix &depmat,int ipos,int jpos, int isize, int jsize);
	
	void RemoveDepend(int myindex, int dependindex);
	
	/// delete all dependency information
	void RemoveDepend();
	
	
	int DependencyDepth(TPZCompMesh &mesh);
	
	int HasDependency()const { return fDependList != 0; }
	
	int CheckDependency(int nshape, TPZCompMesh *mesh, int nstate);
	
	
	TPZDepend *FirstDepend() { return fDependList; }
	
	/** @brief Adds itself and the connects from which it depends to the list*/
	void AddToList(int myindex, TPZCompMesh &mesh, TPZStack<int> &connectlist);
	
	/** @brief Adds itself and the connects from which it depends to the list*/
	void AddToList(int myindex, TPZCompMesh &mesh, std::set<int> &connectlist);
	
	void SetDependenceOrder(int myindex, TPZCompMesh &mesh, int CurrentOrder,TPZVec<int> &connectlist,TPZVec<int> &DependenceOrder);
	
	void ExpandShape(int cind, TPZVec<int> &connectlist, TPZVec<int> &blocksize, TPZFMatrix &phi, TPZFMatrix &dphi);
	
	/**
	 @brief Save the element data to a stream
	 */
	void Write(TPZStream &buf, int withclassid);
	
	/**
	 @brief Read the element data from a stream
	 */
	void Read(TPZStream &buf, void *context);
	
	/**
	 * @brief Copy a connect data structure from an original connect to a
	 * new connect mapping their indexes
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
	 * 
	 * connectlist need to be computed by BuildConnectList
	 */
	static void BuildDependencyOrder(TPZVec<int> &connectlist, TPZVec<int> &DependenceOrder, TPZCompMesh &mesh);
	
};


/**
 @brief Associate a degree of freedom node with a boundary condition
 such boundary condition can be dirichlet, point load or mixed boundary condition.
 @deprecated THIS CLASS IS DEPRECATED
 */
struct TPZConnectBC {
	
	TPZConnect *fConnect;
	TPZBndCond *fBC;
	
	TPZConnectBC() {
		fConnect = 0;
		fBC = 0;
	}
	TPZConnectBC(TPZConnect *nd,TPZBndCond *bc) {
		fConnect = nd;
		fBC = bc;
	}
	
	void Print(TPZCompMesh &mesh,std::ostream &out = std::cout);
	
};

/** @brief Overload operator << to write node connect data */
inline std::ostream & operator<<(std::ostream &out,TPZConnect &con)
{
	out << "seq num: " << con.SequenceNumber()
	<< " nel con: " << con.NElConnected()
	<< " order: " << con.Order()
	<< " hasdepend: " << con.HasDependency();
	return out;
}

#endif

