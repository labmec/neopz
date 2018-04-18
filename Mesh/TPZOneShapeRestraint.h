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
#include "pzmanvector.h"
#include "TPZStream.h"

#include <list>



/// Definition of the retraint associated with the top of the pyramid
struct TPZOneShapeRestraint
{
    /// Faces which are related. First item is the connect index, second item is the degree of freedom
    TPZManVector<std::pair<int64_t,int>,4> fFaces;
    
    /// Orientation of each face
    TPZManVector<int,4> fOrient;
    
    TPZOneShapeRestraint() : fFaces(4), fOrient(4,-1)
    {
        for(int i=0; i<4; i++)
        {
            fFaces[i] = std::make_pair<int64_t,int>(-1, -1);
        }
    }
    TPZOneShapeRestraint(const TPZOneShapeRestraint &copy) : fFaces(copy.fFaces), fOrient(copy.fOrient)
    {
        
    }
    
    TPZOneShapeRestraint &operator=(const TPZOneShapeRestraint &copy)
    {
        fFaces = copy.fFaces;
        fOrient = copy.fOrient;
        return *this;
    }
    
    bool IsInitialized() const
    {
        return fFaces[0].first != -1;
    }
    
    void Print(std::ostream &out) const
    {
        out << "TPZOneShapeRestraint ConnectIndex/Degree of freedom : ";
        for(int i=0; i<4; i++) out << fFaces[i].first << "/" << fFaces[i].second << " ";
        out << std::endl;
        out << "Orientation of each face ";
        for(int i=0; i<4; i++) out << fOrient[i] << " ";
        out << std::endl;
    }
    void Write(TPZStream &buf) const
    {
        TPZManVector<int64_t,4> seqnums(4);
        TPZManVector<int,4> faces(4);
        for (int i=0; i<4; i++)
        {
            seqnums[i] = fFaces[i].first;
            faces[i] = fFaces[i].second;
        }
        buf.Write(&seqnums[0],4);
        buf.Write(&faces[0],4);
        buf.Write(&fOrient[0],4);
    }
    void Read(TPZStream &buf)
    {
        TPZManVector<int64_t,4> seqnums(4);
        TPZManVector<int,4> faces(4);
        buf.Read(&seqnums[0],4);
        buf.Read(&faces[0],4);
        buf.Read(&fOrient[0],4);
        for (int i=0; i<4; i++)
        {
//            fFaces[i] = std::make_pair<int64_t,int>(seqnums[i], faces[i]); // original
            fFaces[i] = std::make_pair(seqnums[i], faces[i]); // To compile in linux : Douglas 
        } 
    }
    
    
};


#endif /* defined(__PZ__TPZOneShapeRestraint__) */
