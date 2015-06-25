//
//  TPZCompMeshTools.cpp
//  PZ
//
//  Created by Philippe Devloo on 6/4/15.
//
//

#include "TPZCompMeshTools.h"
#include "pzelchdiv.h"
#include "pzshapepiram.h"
#include "TPZOneShapeRestraint.h"
#include "pzshapetriang.h"

static TPZOneShapeRestraint SetupPyramidRestraint(TPZCompEl *cel, int side);

using namespace pzshape;

void TPZCompMeshTools::AddHDivPyramidRestraints(TPZCompMesh *cmesh)
{
    long nel = cmesh->NElements();
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZCompElHDiv< pzshape::TPZShapePiram > *pyram = dynamic_cast<TPZCompElHDiv< pzshape::TPZShapePiram > *>(cel);
        if (!pyram) {
            continue;
        }
        
        if (cel->Type() != EPiramide)
        {
            DebugStop();
        }
        TPZGeoEl *gel = cel->Reference();
        /// we need to identify the connect to restraint
        /// start with face 14
        int is;
        for (is=14; is<18; is++) {
            TPZGeoElSide gelside(gel,is);
            TPZCompElSide large = gelside.LowerLevelCompElementList2(1);
            // if the side is constrained, it wont be used as a reference
            if (large) {
                continue;
            }
            // if the side has smaller elements connected to it, discard the side
            TPZStack<TPZCompElSide> celstack;
            gelside.HigherLevelCompElementList2(celstack, 1, 0);
            if (celstack.size()) {
                continue;
            }
            // see if there is a pyramid element  neighbour of the pyramid
            gelside.EqualLevelCompElementList(celstack, 1, 0);
            int localel;
            for (localel=0; localel<celstack.size(); localel++) {
                if (celstack[localel].Reference().Element()->Type() == EPiramide) break;
            }
            TPZCompElHDiv<pzshape::TPZShapePiram> *pir;
            /// the pyramid has a neighbour of pyramid type
            if (localel < celstack.size()) {
                // find out if the face is the restrained face
                pir = dynamic_cast<TPZCompElHDiv<pzshape::TPZShapePiram> *>(celstack[localel].Element());
                if (pir == NULL) {
                    DebugStop();
                }
                const int restrainedface = pir->RestrainedFace();
                const int neighbourside = celstack[localel].Side();
                if (restrainedface == neighbourside) {
                    continue;
                }
            }
            break;
        }
        if (is == 18) {
            DebugStop();
        }
        TPZOneShapeRestraint rest = SetupPyramidRestraint(pyram, is);
        pyram->AddShapeRestraint(rest);
        /// add the restraint datastructure to the neighbour
        TPZGeoElSide gelside(cel->Reference(),is);
        TPZGeoElSide neighbour = gelside.Neighbour();
        TPZCompElSide celneigh = neighbour.Reference();
        TPZInterpolatedElement *interp = dynamic_cast<TPZInterpolatedElement *>(celneigh.Element());
        interp->AddShapeRestraint(rest);
    }
    
        
        
}

static int idf[6] = {2,1,1,2,0,0};

static int nextface(int side, int inc)
{
    return ((side-14)+inc)%4+14;
}

static void GetGeoNodeIds(TPZGeoEl *gel, int side, TPZVec<long> &ids)
{
    int nsidenodes = gel->NSideNodes(side);
    if (nsidenodes != 3) {
        DebugStop();
    }
    for (int i=0; i<3; i++) {
        ids[i] = gel->SideNodePtr(side, i)->Id();
    }
}

static int SideIdf(TPZCompElHDiv<TPZShapePiram> *cel, int side)
{
    TPZManVector<long,3> ids(3);
    GetGeoNodeIds(cel->Reference(),side,ids);
    int transformid = TPZShapeTriang::GetTransformId2dT(ids);
    return idf[transformid];
}
static TPZOneShapeRestraint SetupPyramidRestraint(TPZCompEl *cel, int side)
{
    TPZCompElHDiv<TPZShapePiram> *pir = dynamic_cast<TPZCompElHDiv<TPZShapePiram> *>(cel);
    if (!pir) {
        DebugStop();
    }
    TPZOneShapeRestraint result;
    int mult = 1;
    for (int fc=0; fc<4; fc++)
    {
        int face = nextface(side, fc);
        long nodeindex = pir->SideConnectIndex(0, face);
        result.fFaces[fc] = std::make_pair(nodeindex, SideIdf(pir, face));
        result.fOrient[fc] = mult*pir->SideOrient(face-13);
        mult *= -1;
    }
    return result;
}

