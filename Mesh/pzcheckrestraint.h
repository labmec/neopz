//$Id: pzcheckrestraint.h,v 1.3 2005-02-28 22:08:50 phil Exp $
#ifndef PZCHECKRESTRAINTH
#define PZCHECKRESTRAINTH

#include "pzcompel.h"
#include "pzfmatrix.h"
#include "pzvec.h"

class TPZCompMesh;

/// this class will verify the consistency of the restraints of shape functions along a side
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

	void Print(ostream &out);

private:

	int SmallConnect(int connectid);

	int LargeConnect(int connectid);

	void AddConnect(int connectindex);

	void AddDependency(int smallconnectid, int largeconnectid, TPZFMatrix &depend);

};

#endif
