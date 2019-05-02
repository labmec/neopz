//
// Created by natalia on 28/03/19.
//

#ifndef TRKSolution_h
#define TRKSolution_h

#include "TElastoPlasticData.h"
#include "TPZYCMohrCoulombPV.h"
#include "TPZPlasticStepPV.h"
#include "TPZElastoPlasticMem.h"
#include "TPZMatElastoPlastic2D.h"
#include "TPZTensor.h"


class TRKSolution {
    
    // @TODO:: NVB for your kind understanding please document this class.
protected:
    
    TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> m_elastoplastic_model;

    REAL m_re = -1;

    REAL m_rw = -1;

    REAL m_ure = -1;
    
    REAL m_sigma_re = -1.0;

    int m_n_points = 0;
    
    TPZManVector<TPZElastoPlasticMem,10> m_memory_vector;
    
    TPZElastoPlasticMem m_initial_state_memory;

public:
    
    /// Default constructor
    TRKSolution();

    /// Copy constructor
    TRKSolution(const TRKSolution &  other);

    /// Assignmet constructor
    TRKSolution & operator=(const TRKSolution &  other);

    /// Default destructor
    ~TRKSolution();
    
   void SetElastoPlasticModel(TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> & model);

    void SetExternalRadius(REAL re);

    REAL ExternalRadius();

    void SetWellboreRadius(REAL rw);

    REAL WellboreRadius();

    void SetRadialDisplacement(REAL m_sigma);
    
    REAL RadialRadialDisplacement();
    
    void SetRadialStress(REAL m_sigma);

    REAL RadialStress();
    
    void SetInitialStateMemory(TPZElastoPlasticMem memory);
    
    void SetNumberOfPoints(int n_points);
    
    int GetNumberOfPoints();

    void CreateMaterial();
    
    void FillPointsMemory();
    
    void F(REAL r, REAL ur, REAL sigma_r, REAL &d_ur, REAL &d_sigmar, REAL & lambda, REAL & G);
    
    void RKProcess(std::ostream &out, bool euler);

};


#endif //TRKSolution_h
