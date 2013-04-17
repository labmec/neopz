//
//  pzplasticdiagnostic.h
//  PZ
//
//  Created by phil on 3/18/13.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#ifndef PZ_pzplasticdiagnostic_h
#define PZ_pzplasticdiagnostic_h

class TPZCompMesh;

class TPZPlasticDiagnostic
{

    TPZCompMesh *fMesh;
public:
    
    TPZPlasticDiagnostic(TPZCompMesh *cmesh);
    
    void CheckGlobal();
    
};


#endif
