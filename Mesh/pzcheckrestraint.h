/**
 * @file
 * @brief Contains declaration of TPZCheckRestraint class which verify the consistency of the restraints of shape functions along a side.
 */
//$Id: pzcheckrestraint.h,v 1.4 2005-04-25 02:31:46 phil Exp $
#ifndef PZCHECKRESTRAINTH
#define PZCHECKRESTRAINTH

#include "pzcompel.h"
#include "pzfmatrix.h"
#include "pzvec.h"

class TPZCompMesh;

/**
 * @ingroup CompMesh
 * @brief Will verify the consistency of the restraints of shape functions along a side. \ref CompMesh "Computational Mesh"
 */
class TPZCheckRestraint {
	
	TPZCompElSide fSmall;
	TPZCompElSide fLarge;
	TPZFMatrix fRestraint;
	TPZVec<int> fSmallSize,fSmallPos;
	TPZVec<int> fLargeSize,fLargePos;
	TPZVec<int> fSmallConnect;
	TPZVec<int> fLargeConnect;
	TPZCompMesh *fMesh;
	
	
public:
	
	TPZCheckRestraint(TPZCompElSide small, TPZCompElSide large);
	
	TPZFMatrix &RestraintMatrix();
	
	int CheckRestraint();
	
	void Diagnose();
	
	void Print(std::ostream &out);
	
private:
	
	int SmallConnect(int connectid);
	
	int LargeConnect(int connectid);
	
	void AddConnect(int connectindex);
	
	void AddDependency(int smallconnectid, int largeconnectid, TPZFMatrix &depend);
	
};

#endif
