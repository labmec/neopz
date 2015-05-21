//
//  TPZOneShapeRestraint.h
//  PZ
//
//  Created by Philippe Devloo on 5/16/15.
//
//

#ifndef __PZ__TPZOneShapeRestraint__
#define __PZ__TPZOneShapeRestraint__

#include <stdio.h>

#include <list>



/// Definition of the retraint associated with the top of the pyramid
struct TPZOneShapeRestraint
{
    /// Faces which are related. First item is the connect index, second item is the degree of freedom
    std::pair<long,int> fFaces[4];
    
    /// Orientation of each face
    int fOrient[4];
    
};


#endif /* defined(__PZ__TPZOneShapeRestraint__) */
