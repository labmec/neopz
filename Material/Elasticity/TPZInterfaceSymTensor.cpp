#include "TPZInterfaceSymTensor.h"
#include "TPZMixedSymElasticityND.h"
// Constructor
TPZInterfaceSymTensor::TPZInterfaceSymTensor(int id, int dimension)
    : TPZLagrangeMultiplierCS<STATE>(id,dimension) {
    // Initialization code
}

// Destructor
TPZInterfaceSymTensor::~TPZInterfaceSymTensor() {
    // Cleanup code
}

// Override ContributeInterface
void TPZInterfaceSymTensor::ContributeInterface(const TPZMaterialDataT<TVar> &data,
                                 const std::map<int, TPZMaterialDataT<TVar>> &dataleft,
                                 const std::map<int, TPZMaterialDataT<TVar>> &dataright, REAL weight,
                                 TPZFMatrix<TVar> &ek,
                                                TPZFMatrix<TVar> &ef){
    if(dataleft.size() != 1 || dataright.size() != 1) DebugStop();
    const TPZMaterialDataT<STATE> &leftdata = (dataleft.begin()->second);
    const TPZMaterialDataT<STATE> &rightdata = (dataright.begin()->second);
    // we assume dataright is the displacement lagrange multiplier
    TPZManVector<REAL,3> normal(2,0.);
    {
        auto &axes = leftdata.axes;
        if(axes.Rows() != 1) DebugStop();
        normal[0] = axes(0,1);
        normal[1] = -axes(0,0);
    }

    int64_t nphitensor = leftdata.phi.Rows();

    TPZFMatrix<REAL> BMatrixL (2,nphitensor*3,0.);

    int exx = TPZMixedSymElasticityND::ExxS;
    int exy = TPZMixedSymElasticityND::ExyS;
    int eyy = TPZMixedSymElasticityND::EyyS;
    for(int p = 0; p<nphitensor; p++) {
        BMatrixL(0,3*p+exx) = leftdata.phi(p,0)*normal[0];
        BMatrixL(1,3*p+eyy) = leftdata.phi(p,0)*normal[1];
        BMatrixL(0,3*p+exy) = leftdata.phi(p,0)*normal[1];
        BMatrixL(1,3*p+exy) = leftdata.phi(p,0)*normal[0];
    }

    int64_t nphidisp = rightdata.phi.Rows();
    TPZFMatrix<REAL> BMatrixR(2,nphidisp*2,0.);
    for(int p = 0; p<nphidisp; p++) {
        BMatrixR(0,2*p) = rightdata.phi(p,0);
        BMatrixR(1,2*p+1) = rightdata.phi(p,0);
    }

    ek.AddContribution(0,3*nphitensor, BMatrixL, 1, BMatrixR, 0, weight*fMultiplier);
    ek.AddContribution(3*nphitensor, 0, BMatrixR, 1, BMatrixL, 0, weight*fMultiplier);

}


