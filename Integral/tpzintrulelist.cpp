//
// C++ Implementation: tpzintrulelist
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//

#include "tpzintrulelist.h"
#include "pzerror.h"
#include "tpzintrule.h"
#include "tpzintrulet.h"
#include "tpzintrulet3d.h"
#include "tpzintrulep3d.h"

TPZIntRuleList::TPZIntRuleList(){
	
	static int	first = 1;
	
	if(first != 1) {
		PZError << "second initialization of the integration rule list\n"
        " something fishy is going on!\n";
		//		PZError.show();
	}
	
	first++;
	
	intavail    = TPZIntRule::NUMINT_RULES;//reta, quadrilatero, cubo
	intavailT   = TPZIntRuleT::NUMINT_RULEST;//triangulo
	intavailT3D = TPZIntRuleT3D::NUMINT_RULEST3D;//tetraedro
	intavailP3D = TPZIntRuleP3D::NUMINT_RULESP3D;//piramide
	
	intlist    = new TPZIntRule*[intavail];
	intlistT   = new TPZIntRuleT*[intavailT];
	intlistT3D = new TPZIntRuleT3D*[intavailT3D];
	intlistP3D = new TPZIntRuleP3D*[intavailP3D];
	
	
	if(intlist == NULL || intlistT == NULL || intlistT3D == NULL || intlistP3D == NULL){
		PZError << "TPZIntRuleList unable to initialize a list of integration rules\n";
		//		PZError.show();
		intavail = 0;
		
	} else {
		
		int i;
		for(i = 1; i<=intavail; ++i) {
			intlist[i-1] = new TPZIntRule(i-1);
			if(intlist[i-1] == NULL) {
				PZError << "TPZIntRuleList error: some integration rules"
				"could not be initialized\n";
				//				PZError.show();
			}
		}
		
		for(i = 0; i<intavailT; ++i) {
			intlistT[i] = new TPZIntRuleT(i);
			if(intlistT[i] == NULL) {
				PZError << "TPZIntRuleList error: some integration rules"
				"for triangles could not be initialized\n";
				//				PZError.show();
			}
		}
		for(i = 1; i<=intavailT3D; ++i) {
			intlistT3D[i-1] = new TPZIntRuleT3D(i);
			if(intlistT3D[i-1] == NULL) {
				PZError << "TPZIntRuleList error: some integration rules"
				"for tetrahedros could not be initialized\n";
				//				PZError.show();
			}
		}
		for(i = 1; i<=intavailP3D; ++i) {
			intlistP3D[i-1] = new TPZIntRuleP3D(i);
			if(intlistP3D[i-1] == NULL) {
				PZError << "TPZIntRuleList error: some integration rules"
				"for pyramides could not be initialized\n";
				//				PZError.show();
			}
		}
		
		
	}
}

//***************************************
//***************************************
TPZIntRuleList::~TPZIntRuleList(){
	
	if (!intlist) return;
	int i;
	for(i=0 ; i<intavail ; ++i)    if (intlist[i])    delete intlist[i];
	for(i=0 ; i<intavailT ; ++i)   if (intlistT[i])   delete intlistT[i];
	for(i=0 ; i<intavailT3D ; ++i) if (intlistT3D[i]) delete intlistT3D[i];
	for(i=0 ; i<intavailP3D ; ++i) if (intlistP3D[i]) delete intlistP3D[i];
	delete []intlistT3D;
	delete []intlistP3D;
	delete []intlist;
	delete []intlistT;
}

//***************************************
//***************************************
TPZIntRule* TPZIntRuleList::GetRule(int fNumInt) {
	
	if (fNumInt < 0 || fNumInt >= intavail) {
		
		static bool verbose = true;
		if(verbose){  
			PZError << "\nERROR(TPZIntRuleList::getrule)-> Numint = " << fNumInt;
		}
		//		PZError.show();
		fNumInt = intavail-1;
		if(verbose){
			PZError << "\n                     precision obtained = " << fNumInt << "\n";
			PZError.flush();
			verbose = false;
		}
		//return NULL;
	}
	if(fNumInt == 0) fNumInt = 1;
	//return intlist[fNumInt-1];
	return intlist[fNumInt];
}

//**************************************
//**************************************
TPZIntRuleT* TPZIntRuleList::GetRuleT(int precision) {
	
	if (precision < 0 || precision >= intavailT) {
		static bool verbose = true;
		if(verbose){
			PZError << "\nERROR(TPZIntRuleList::getrule)-> precision required = " << precision;
		}
		precision = intavailT-1;
		if(verbose){
			PZError << "\n                                 precision obtained = " << precision;
			PZError.flush();
			verbose = false;
		}
	}
	
	return intlistT[precision];
}

//**************************************
TPZIntRuleT3D* TPZIntRuleList::GetRuleT3D(int precision) {
	
	if (precision >= intavailT3D) {
		PZError << "\nERROR(TPZIntRuleList::getrule)-> precision required = " << precision << std::endl;
		precision = intavailT3D-1;
		PZError << "\nERROR(TPZIntRuleList::getrule)-> precision gotten = " << precision << std::endl;
	}
	if (precision < 0) {
		PZError << "\nERROR(TPZIntRuleList::getrule)-> precision required = " << precision << std::endl;
		precision = 0;
		PZError << "\nERROR(TPZIntRuleList::getrule)-> precision gotten = " << precision << std::endl;
	}
	
	return intlistT3D[precision];
}


TPZIntRuleList  gIntRuleList;
