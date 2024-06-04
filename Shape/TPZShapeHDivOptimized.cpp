#include "TPZShapeHDivOptimized.h"
#include "TPZShapeH1.h"
#include "TPZShapeHCurl.h"
#include "pzshapelinear.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"
#include "pzshapetetra.h"
#include "pzshapecube.h"
#include "pzshapeprism.h"
#include "pzshapepiram.h"
#include "pzgeoel.h"
#include "pzgnode.h"
#include "TPZShapeData.h"
#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.shapehdiv");
#endif

template <class TSHAPE>
TPZShapeHDivOptimized<TSHAPE>::TPZShapeHDivOptimized() {}

template struct TPZShapeHDivOptimized<pzshape::TPZShapeLinear>;

template struct TPZShapeHDivOptimized<pzshape::TPZShapeTriang>;

template struct TPZShapeHDivOptimized<pzshape::TPZShapeQuad>;

template struct TPZShapeHDivOptimized<pzshape::TPZShapeTetra>;

template struct TPZShapeHDivOptimized<pzshape::TPZShapeCube>;

template struct TPZShapeHDivOptimized<pzshape::TPZShapePrism>;