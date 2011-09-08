/**
 * @file
 * @brief DEPRECATED FILE. This file contains declaration of TPZConnect class which represents \n
 * a set of shape functions associated with a computational element
 */
//$Id: pzconnect.h,v 1.19 2010-08-25 03:05:06 phil Exp $
//HEADER FILE FOR CLASS NODE

#ifndef  PZCONNECTBCH
#define  PZCONNECTBCH

#include "pzfmatrix.h"
#include "pzstack.h"
#include <iostream>
#include <set>


class TPZBndCond;
class TPZCompMesh;
class TPZBlock;
class TPZStream;


/**
 * @brief Associate a degree of freedom node with a boundary condition such boundary condition can be dirichlet, \n
 * point load or mixed boundary condition.
 * @deprecated DEPRECATED BC connect CLASS.
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

/// Prints data related with connect and boundary condition
void TPZConnectBC::Print(TPZCompMesh &mesh,std::ostream &out){
	out << "Connect boundary condition :\n";
	if(fConnect) fConnect->Print(mesh,out);
	if(fBC) fBC->Print(out);
}

#endif

