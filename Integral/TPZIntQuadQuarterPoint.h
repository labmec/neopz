//
//  TPZIntQuadQuarterPoint.h
//  PZ
//
//  Created by Philippe Devloo on 6/1/15.
//
//

#ifndef __PZ__TPZIntQuadQuarterPoint__
#define __PZ__TPZIntQuadQuarterPoint__

#include "pzquad.h"
#include <stdio.h>

class TPZIntQuadQuarterPoint : public TPZIntQuad
{
    /// corner associated with the integration rule
    int fCorner;
    /// number of integration points
    int fNPoints;
    
public:
    
    TPZIntQuadQuarterPoint(int ordksi) : TPZIntQuad(ordksi), fCorner(0)
    {
        fNPoints = TPZIntQuad::NPoints()*2;
    }
    
    TPZIntQuadQuarterPoint(const TPZIntQuadQuarterPoint &copy) : TPZIntQuad(copy), fCorner(copy.fCorner), fNPoints(copy.fNPoints)
    {
        
    }
    
    virtual int NPoints() const override
    {
        return fNPoints;
    }
    
    /// return the weight and position of the integration point
    virtual void Point(int ip, TPZVec<REAL> &pos, REAL &w) const override;
    
    /// set the integration order
    virtual void SetOrder(TPZVec<int> &ord,int type = 0) override;
    
    /// set the corner where the singularity resides
    void SetCorner(int corner)
    {
#ifdef PZDEBUG
        if (corner < 0 || corner > 3) {
            DebugStop();
        }
#endif
        fCorner = corner;
    }

    virtual TPZIntPoints* Clone() const override
    {
        return new TPZIntQuadQuarterPoint(*this);
    }

};

#endif /* defined(__PZ__TPZIntQuadQuarterPoint__) */
