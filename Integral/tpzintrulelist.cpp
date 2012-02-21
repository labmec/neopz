/**
 * @file
 * @brief Contains the implementation of the TPZIntRuleList methods. It is created a list of all integration rules avaliables
 * the rule was created at first time that it is required.
 */

#include "tpzintrulelist.h"
#include "pzerror.h"
#include "tpzgaussrule.h"
#include "tpzintrulet.h"
#include "tpzintrulet3d.h"
#include "tpzintrulep3d.h"

TPZIntRuleList  gIntRuleList;

TPZIntRuleList::TPZIntRuleList() {
	
	static int	first = 1;
	
	if(first != 1) {
		PZError << "second initialization of the integration rule list\n"
        << " something fishy is going on!\n";
		DebugStop();
	}
	
	// Cleaning integration rules vectors
	fintlist.Resize(0);
	fintlistT.Resize(0);
	fintlistT3D.Resize(0);
	fintlistP3D.Resize(0);

	first++;
}

TPZIntRuleList::~TPZIntRuleList() {
	int i;
	/** Deleting dinamic allocation of the Gauss integration rules at the vector */
	for(i=0 ; i<fintlist.NElements() ; i++)   if (fintlist[i])   delete fintlist[i];
	fintlist.Resize(0);
	/** Deleting dinamic allocation of the cubature rules for triangle at the vector */
	for(i=0 ; i<fintlistT.NElements() ; i++)   if (fintlistT[i])   delete fintlistT[i];
	fintlistT.Resize(0);
	/** Deleting dinamic allocation of the cubature rules for tetrahedra at the vector */
	for(i=0 ; i<fintlistT3D.NElements() ; i++)   if (fintlistT3D[i])   delete fintlistT3D[i];
	fintlistT3D.Resize(0);
	
	/** Deleting dinamic allocation of the cubature rules for pyramid at the vector */
	for(i=0 ; i<fintlistP3D.NElements() ; i++)   if (fintlistP3D[i])   delete fintlistP3D[i];
	fintlistP3D.Resize(0);
}

TPZGaussRule* TPZIntRuleList::GetRule(int order,int type) {
	if(order < 0) {
		order = 1;
	}
	if(type == 1) {
		if(order > TPZGaussRule::NRULESLOBATTO_ORDER)
			order = 1;
	}
	else {
		if(order > TPZGaussRule::NRULESLEGENDRE_ORDER)
			order = TPZGaussRule::NRULESLEGENDRE_ORDER;
	}
	// The vectors of integration rules will be created based on order value
	if(fintlist.NElements()<order+1)
		fintlist.Resize(order+1,NULL);
	if(!fintlist[order])
		fintlist[order] = new TPZGaussRule(order,type);
	else if(fintlist[order]->Type() != type) {
		delete fintlist[order];
		fintlist[order] = new TPZGaussRule(order,type);
	}
	return fintlist[order];
}

//**************************************
TPZIntRuleT* TPZIntRuleList::GetRuleT(int order) {
	if(order < 0) order = 1;
	if(order > TPZIntRuleT::NRULESTRIANGLE_ORDER) {
		order = TPZIntRuleT::NRULESTRIANGLE_ORDER;
	}

	if(fintlistT.NElements()<order+1)
		fintlistT.Resize(order+1,NULL);
	if(!fintlistT[order])
		fintlistT[order] = new TPZIntRuleT(order);
	return fintlistT[order];
}

//**************************************
TPZIntRuleT3D* TPZIntRuleList::GetRuleT3D(int order) {
	if(order < 0) order = 1;
	if(order > TPZIntRuleT3D::NRULESTETRAHEDRA_ORDER) {
		order = TPZIntRuleT3D::NRULESTETRAHEDRA_ORDER;
	}
	
	if(fintlistT3D.NElements()<order+1)
		fintlistT3D.Resize(order+1,NULL);
	if(!fintlistT3D[order])
		fintlistT3D[order] = new TPZIntRuleT3D(order);
	return fintlistT3D[order];
}

//**************************************
TPZIntRuleP3D* TPZIntRuleList::GetRuleP3D(int order) {
	if(order < 0) order = 1;
	if(order > TPZIntRuleP3D::NRULESPYRAMID_ORDER) {
		order = TPZIntRuleP3D::NRULESPYRAMID_ORDER;
	}

	if(fintlistP3D.NElements()<order+1)
		fintlistP3D.Resize(order+1,NULL);
	if(!fintlistP3D[order])
		fintlistP3D[order] = new TPZIntRuleP3D(order);
	return fintlistP3D[order];
}
