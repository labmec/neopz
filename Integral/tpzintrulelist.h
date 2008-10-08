//
// C++ Interface: tpzintrulelist
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef TPZINTRULELIST_H
#define TPZINTRULELIST_H

class TPZIntRule;
class TPZIntRuleT;
class TPZIntRuleT3D;
class TPZIntRuleP3D;
/**
This class creates instances of all integration rules for rapid selection

	@author Philippe R. B. Devloo <phil@fec.unicamp.br>
*/
class TPZIntRuleList{

  
    int		intavail;	// number of integration rules available
    int        	intavailT;  // number of integration rules available for triangles
    int        	intavailT3D;
    int        	intavailP3D;
    TPZIntRule	**intlist;	 	// pointer to an array of integration rules
    TPZIntRuleT   **intlistT; 	// pointer to an array of integration rules
    TPZIntRuleT3D **intlistT3D;
    TPZIntRuleP3D **intlistP3D;

    public :

      TPZIntRuleList();	// method which initializes all integration rules
  // should be called only once!

      ~TPZIntRuleList();

      TPZIntRule *GetRule(int numint);	// returns a pointer to an integration
  // rule with numint integration points

      TPZIntRuleT *GetRuleT(int numint); // returns a pointer to an integration
  // rule for a triangle
      TPZIntRuleT3D *GetRuleT3D(int numint);
      TPZIntRuleP3D *GetRuleP3D(int numint);
  };

  extern  TPZIntRuleList  gIntRuleList;

#endif
