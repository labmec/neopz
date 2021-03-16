#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzgeoelside.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"
#include "tpzgeoblend.h"

#include <TPZGenGrid2D.h>
#include "TPZMatElasticity2D.h"
#include "TPZInterfaceEl.h"


#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "tpzcompmeshreferred.h"
#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "pzanalysis.h"
#include <tpzarc3d.h>

#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "pzstrmatrix.h"
#include "TPZFrontNonSym.h"
#include "TPZFrontSym.h"
#include "TPBSpStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzbstrmatrix.h"
#include "pzl2projection.h"

#include "pzpoisson3d.h"
#include "pzpoisson3dreferred.h"

#include "pzelasmat.h"
#include "pzmultiphysicselement.h"
#include "pzmultiphysicscompel.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZSpStructMatrix.h"
#include "pzlog.h"
#include <iostream>
#include <string>
#include "TPZVTKGeoMesh.h"
#include "pzfunction.h"
#include "TPZReadGIDGrid.h"
#include "pzmultiphysicselement.h"
#include "TPZMultiphysicsInterfaceEl.h"

#ifdef USING_MATLIB // NS: Shouldnt there be a define in here? To choose using the lib matlib?
#include "MatLib.h"
USING_MATLIB_NAMESPACE

#include <cmath>
#include <set>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.elasticity"));
#endif

// Dummy Boundary Conditions
const int dirichlet = 0;
const int neumann = 1;

static bool oldmat = false;

// Defintions of Implemented Methods
TPZCompMesh *ComputationalElasticityMesh(TPZGeoMesh * gmesh,int pOrder);
void IterativeProcess(TPZAnalysis *an, std::ostream &out, int numiter = 20);

//	This Solve Different analysis
void SolveSist(TPZAnalysis *an, TPZCompMesh *fCmesh);

//	These are tools for spatial and polynomial refinement and Postprocess of solutions
void PostProcessElasticity(TPZAnalysis &an, std::string plotfile);
/**
 *  Uniformly refine a geometric mesh
 *
 *  @param gMesh Object containing the elements that will be refined
 *  @param nh    number of uniform refinements
 */
void UniformRefinement(TPZGeoMesh  *gMesh, int nh);
void UniformRefinement(TPZGeoMesh *gMesh, int nh, int MatId);
void RefinElemComp(TPZCompMesh  *cMesh, int indexEl);
void RefinUniformElemComp(TPZCompMesh  *cMesh, int ndiv);
void ReservoirPressure(const TPZVec<STATE> &x, TPZVec<STATE> &p,  TPZFMatrix<STATE> &gradp);
#endif


int main(int argc, char *argv[])
{

    std::string dirname = PZSOURCEDIR;

#ifdef USING_MATLIB
    std::cout << "Printing something\n";
    ModelDictionary::list(std:: cout);
    
    /**
     *  Create a constitutive model for MatLib
     *
     *  @param "ISOTROPIC_ELASTICITY" Example model, in this case isotropic elastic
     *
     *  @return A dynamically created object
     */
    ConstitutiveModel* model = ModelDictionary::build("ISOTROPIC_ELASTICITY");
    MaterialProperties prop("steel");
    prop.setProperty("YOUNG_MODULUS",2.1e11);
    prop.setProperty("POISSON_COEFFICIENT",0.3);
    MaterialModel mater(*model,prop);
    mater.initialize();
    

    MaterialState locState0,locState1;
    mater.initState(locState0);
    mater.initState(locState1);
    locState1.grad[0] += 1.e-3;
    ParameterSet extPar;
    MatLibMatrix K(model->nExtVar());
    std::cout << "eps=" << locState1.grad << std::endl;
    mater.updateState(extPar,locState0,locState1,1.0,K,true);
    std::cout << "sig=" << locState1.flux << std::endl;
    std::cout << "ïnternal=" << locState1.internal << std::endl;
    std::cout << "K="<< K << std::endl;
#endif
  
  return 0;
}
