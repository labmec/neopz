#include "TPZShapeNewHDiv.h"
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
TPZShapeNewHDiv<TSHAPE>::TPZShapeNewHDiv() {}

template struct TPZShapeNewHDiv<pzshape::TPZShapeLinear>;

template struct TPZShapeNewHDiv<pzshape::TPZShapeTriang>;

template struct TPZShapeNewHDiv<pzshape::TPZShapeQuad>;

template struct TPZShapeNewHDiv<pzshape::TPZShapeTetra>;

template struct TPZShapeNewHDiv<pzshape::TPZShapeCube>;

template struct TPZShapeNewHDiv<pzshape::TPZShapePrism>;