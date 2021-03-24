/**
 * \file
 * @brief DEPRECATED FILE. This file contains implementation of the methods to Coordinate System Class.
 */
//METHODS DEFINITION FOR CLASS COSYS

#include "pzcosys.h"
#include "pzcartsys.h"
#include "pzvec.h"
using namespace std;

//***************************************
//***************************************
TPZCosys::TPZCosys() {
    fNumber = 0;
	fReference = NULL;
}

TPZCosys::TPZCosys(const TPZCosys &cp)
{
	fNumber = cp.fNumber;
	fReference = cp.fReference;
}

//***************************************
//***************************************
TPZCosys::TPZCosys(int num,  TPZCartsys* ref) {
    fNumber = num;
    fReference = ref;
}

//***************************************
//***************************************
TPZCosys::~TPZCosys(){
	//	if (fReference) delete fReference;
}

//***************************************
//***************************************
void TPZCosys::Reset() {
	fNumber = 0;
	fReference = NULL;
}

//***************************************
//***************************************
void TPZCosys::Print(std::ostream& out) {
	out << "\nCosys number " << fNumber << " and Type ";
	switch(this->Type()) {
		case cartesian: out << "cartesian\n"; break;
		case cylindric: out << "cylindric\n"; break;
		case esferic:   out << "esferic\n"; break;
	}
	if (fReference) {
		out << "Reference system:";
		fReference->Print(out);
	} else {
		out << "NULL fReference system\n";
	}
	//    out << "Origin-> " << fOrigin[0] << " " << fOrigin[1] << " " << fOrigin[2];
	out << "\nTransfer matrix:\n";
	/*        << fTr[0][0] << " " << fTr[0][1] << " " << fTr[0][2] << "\n"
	 << fTr[1][0] << " " << fTr[1][1] << " " << fTr[1][2] << "\n"
	 << fTr[2][0] << " " << fTr[2][1] << " " << fTr[2][2] << "\n";
	 */}

//***************************************
//***************************************
void TPZCosys::ToSpecific(TPZVec<REAL> &x, TPZCartsys *ref){
	ToReference(x);
	if (fReference != ref) fReference->ToSpecific(x,ref);
}

//***************************************
//***************************************

void TPZCosys::ToGlobal(TPZVec<REAL> &point) {
	ToSpecific(point,NULL);
}

//***************************************
//***************************************
void TPZCosys::FromGlobal(TPZVec<REAL> &point) {
	FromReference(point);
	if (fReference) fReference->FromGlobal(point);
}
