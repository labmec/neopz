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
 * @brief Will verify the consistency of the restraints of shape functions along a side. \ref CompElement "Computational Element"
 */
class TPZCheckRestraint {
	/** @brief Computational element with side included into the fLarge */
	TPZCompElSide fSmall;
	/** @brief Computational element with side including fSmall */
	TPZCompElSide fLarge;
	/** @brief Restraints matrix */
	TPZFMatrix<REAL> fRestraint;
	/** @brief Number of shape function corresponding to the connect associated with side of the small element */
	TPZVec<int> fSmallSize;
	TPZVec<int> fSmallPos;
	/** @brief Number of shape function corresponding to the connect associated with side of the large element */
	TPZVec<int> fLargeSize;
	TPZVec<int> fLargePos;
	/** @brief Stores the indexes of the connect over the side for small element */
	TPZVec<int> fSmallConnect;
	/** @brief Stores the indexes of the connect over the side for large element */
	TPZVec<int> fLargeConnect;
	 /** @brief Pointer for computational mesh containing the computational elements given */
	TPZCompMesh *fMesh;
	
	
public:
	/** @brief Constructor with small and large element with commom side */
	TPZCheckRestraint(TPZCompElSide small, TPZCompElSide large);
	/** @brief Returns the restraint matrix */
	TPZFMatrix<REAL> &RestraintMatrix();
	
	/** @brief Gets the shape functions over the sides of the small and large elements and check the matrix restraint making multiplication of matrizes */
	int CheckRestraint();
	/** @brief Get the small element and check restraints with all elements with lower dimension on side corresponding */
	void Diagnose();
	/** @brief Prints the information into the computational elements and side and geometric information also */
	void Print(std::ostream &out);
	
private:
	
	int SmallConnect(int connectid);
	
	int LargeConnect(int connectid);
	
	void AddConnect(int connectindex);
	
	void AddDependency(int smallconnectid, int largeconnectid, TPZFMatrix<REAL> &depend);
	
};

#endif
