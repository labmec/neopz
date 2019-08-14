//
// Created by natalia on 17/05/19.
//

#include "TPZConstitutiveLawProcessor.h"
#include "TPZMatWithMem.h"
#include <functional>

#ifndef USING_CUDA

#include "FunctionsStressStrain.h"
#include "FunctionsTangentOperator.h"

#endif

#ifdef USING_TBB

#include "tbb/parallel_for.h"
#include "tbb/tick_count.h"

#endif

#include "Timer.h"


TPZConstitutiveLawProcessor::TPZConstitutiveLawProcessor() : fNpts(-1), fWeight(0), fMaterial(), fPlasticStrain(0, 0),
                                                             fMType(0, 0), fAlpha(0, 0) {
#ifdef USING_CUDA
    fCudaCalls = new TPZCudaCalls();

    dWeight.resize(0);
    dSigma.resize(0);
    dStrain.resize(0);
    dPlasticStrain.resize(0);
    dMType.resize(0);
    dAlpha.resize(0);
#endif
}

TPZConstitutiveLawProcessor::TPZConstitutiveLawProcessor(int npts, TPZVec<REAL> weight, TPZMaterial *material) : fNpts(
        -1), fWeight(0), fMaterial(), fPlasticStrain(0, 0), fMType(0, 0), fAlpha(0, 0) {
    SetUpDataByIntPoints(npts);
    SetWeightVector(weight);
    SetMaterial(material);
}

TPZConstitutiveLawProcessor::~TPZConstitutiveLawProcessor() {

}

TPZConstitutiveLawProcessor::TPZConstitutiveLawProcessor(const TPZConstitutiveLawProcessor &copy) {
    fNpts = copy.fNpts;
    fWeight = copy.fWeight;
    fMaterial = copy.fMaterial;
    fPlasticStrain = copy.fPlasticStrain;
    fMType = copy.fMType;
    fAlpha = copy.fAlpha;

#ifdef USING_CUDA
    dWeight = copy.dWeight;
    dSigma = copy.dSigma;
    dStrain = copy.dStrain;
    dPlasticStrain = copy.dPlasticStrain;
    dMType = copy.dMType;
    dAlpha = copy.dAlpha;
#endif
}

TPZConstitutiveLawProcessor &TPZConstitutiveLawProcessor::operator=(const TPZConstitutiveLawProcessor &copy) {
    if (&copy == this) {
        return *this;
    }

    fNpts = copy.fNpts;
    fWeight = copy.fWeight;
    fMaterial = copy.fMaterial;
    fPlasticStrain = copy.fPlasticStrain;
    fMType = copy.fMType;
    fAlpha = copy.fAlpha;

#ifdef USING_CUDA
    dWeight = copy.dWeight;
    dSigma = copy.dSigma;
    dStrain = copy.dStrain;
    dPlasticStrain = copy.dPlasticStrain;
    dMType = copy.dMType;
    dAlpha = copy.dAlpha;
#endif

    return *this;
}

void TPZConstitutiveLawProcessor::SetUpDataByIntPoints(int64_t npts) {
    fNpts = npts;

    fSigma.Resize(6 * fNpts, 1);
    fSigma.Zero();

    fStrain.Resize(6 * fNpts, 1);
    fStrain.Zero();

    fPlasticStrain.Resize(6 * fNpts, 1);
    fPlasticStrain.Zero();

    fMType.Resize(1 * fNpts, 1);
    fMType.Zero();

    fAlpha.Resize(1 * fNpts, 1);
    fAlpha.Zero();

#ifdef USING_CUDA
    dSigma.resize(6 * fNpts);
    dSigma.Zero();
    
    dStrain.resize(6 * fNpts);
    dStrain.Zero();

    dPlasticStrain.resize(6 * fNpts);
    dPlasticStrain.Zero();

    dMType.resize(1 * fNpts);
    dMType.Zero();

    dAlpha.resize(1 * fNpts);
    dAlpha.Zero();

#endif
}

void TPZConstitutiveLawProcessor::SetWeightVector(TPZVec<REAL> &weight) {
    fWeight = weight;
}

void TPZConstitutiveLawProcessor::SetMaterial(TPZMaterial *material) {
    fMaterial = dynamic_cast<TPZMatElastoPlastic2D<TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem> *>(material);
}

void TPZConstitutiveLawProcessor::De(TPZFMatrix<REAL> &De) {

    REAL lambda = fMaterial->GetPlasticModel().fER.Lambda();
    REAL mu = fMaterial->GetPlasticModel().fER.G();

    De.Zero();

    De(_XX_, _XX_) += lambda;
    De(_XX_, _YY_) += lambda;
    De(_XX_, _ZZ_) += lambda;
    De(_YY_, _XX_) += lambda;
    De(_YY_, _YY_) += lambda;
    De(_YY_, _ZZ_) += lambda;
    De(_ZZ_, _XX_) += lambda;
    De(_ZZ_, _YY_) += lambda;
    De(_ZZ_, _ZZ_) += lambda;

    int i;
    for (i = 0; i < 6; i++) De(i, i) += mu;

    De(_XX_, _XX_) += mu;
    De(_YY_, _YY_) += mu;
    De(_ZZ_, _ZZ_) += mu;
}

void TPZConstitutiveLawProcessor::ComputeSigma(TPZFMatrix<REAL> &glob_delta_strain, TPZFMatrix<REAL> &glob_sigma) {
#ifndef USING_CUDA

    int64_t rows = glob_delta_strain.Rows();
    int64_t cols = glob_delta_strain.Cols();
    glob_sigma.Resize(rows, cols);


    // Variables required for lambda's capture block
    TPZMatWithMem<TPZElastoPlasticMem> *mat = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *>(fMaterial);


    // The constitutive law is computing assuming full tensors
#ifdef USING_TBB

    tbb::parallel_for(size_t(0), size_t(fNpts), size_t(1),
                      [this, mat, & glob_delta_strain, & glob_sigma](size_t &ipts) {

                          TPZFNMatrix<6, REAL> strain(6, 1, 0.);
                          TPZFNMatrix<6, REAL> elastic_strain(6, 1, 0.);
                          TPZFNMatrix<6, REAL> sigma(6, 1, 0.);
                          TPZFNMatrix<6, REAL> plastic_strain(6, 1, 0.);
                          TPZFNMatrix<3, REAL> aux_tensor(3, 1, 0.);

                          REAL lambda = fMaterial->GetPlasticModel().fER.Lambda();
                          REAL mu = fMaterial->GetPlasticModel().fER.Mu();
                          REAL mc_phi = fMaterial->GetPlasticModel().fYC.Phi();
                          REAL mc_psi = fMaterial->GetPlasticModel().fYC.Psi();
                          REAL mc_cohesion = fMaterial->GetPlasticModel().fYC.Cohesion();
                          REAL G = mu;
                          REAL K = lambda + 2 * mu / 3;


                          REAL alpha;
                          int mtype;

                          //Get from delta strain vector
                          glob_delta_strain.GetSub(3 * ipts, 0, 3, 1, aux_tensor);

                          //Get last strain vector
                          fStrain.GetSub(6 * ipts, 0, 6, 1, strain);

                          // Translate and
                          ComposeStrain(&aux_tensor(0, 0), &strain(0, 0));

                          //Get from plastic strain vector
                          fPlasticStrain.GetSub(6 * ipts, 0, 6, 1, plastic_strain);

                          TPZFMatrix<REAL> stress_eigenvectors(9, 1, 0.);
                          TPZFMatrix<REAL> sigma_projected(3, 1, 0.);
                          TPZFMatrix<REAL> gradient(9, 1, 0.);
                          TPZFMatrix<REAL> tangent(36, 1, 0.);

                          TPZFMatrix<REAL> strain_eigenvalues(3, 1, 0.);
                          TPZFMatrix<REAL> strain_eigenvectors(36, 1, 0.);

                          // Return Mapping components
//                          ElasticStrain(&plastic_strain(0, 0), &strain(0, 0), &elastic_strain(0, 0));
//                          ComputeTrialStress(&elastic_strain(0, 0), &sigma(0, 0), mu, lambda);
//                          SpectralDecomposition(&sigma(0, 0), &sigma_projected(0, 0), &stress_eigenvectors(0, 0));
//                          SpectralDecomposition(&strain(0, 0), &strain_eigenvalues(0, 0), &strain_eigenvectors(0, 0));
//                          ProjectSigma(&sigma_projected(0, 0), &sigma_projected(0, 0), mtype, alpha, mc_phi, mc_psi,
//                                       mc_cohesion, K, G);
//                          TangentOperator(&gradient(0, 0), &sigma_projected(0, 0), &strain_eigenvalues(0, 0),
//                                          &strain_eigenvectors(0, 0), &tangent(0, 0), G, lambda);
//                          ReconstructStressTensor(&sigma_projected(0, 0), &stress_eigenvectors(0, 0), &sigma(0, 0));

                          // Update plastic strain
                          ComputeStrain(&sigma(0, 0), &elastic_strain(0, 0), mu, lambda);
                          PlasticStrain(&strain(0, 0), &elastic_strain(0, 0), &plastic_strain(0, 0));

                          //Copy to stress vector
                          TranslateStress(&sigma(0, 0), &aux_tensor(0, 0));

                          aux_tensor *= fWeight[ipts];
                          glob_sigma.PutSub(3 * ipts, 0, aux_tensor);

                          if (mat->GetUpdateMem()) {

                              fSigma.AddSub(6 * ipts, 0, sigma);

                              //Accumulate to strain vector
                              fStrain.AddSub(6 * ipts, 0, strain);

                              //Copy to plastic strain vector
                              fPlasticStrain.PutSub(6 * ipts, 0, plastic_strain);

                              //Copy to MType and Alpha vectors
                              fAlpha(ipts, 0) = alpha;
                              fMType(ipts, 0) = mtype;
                          }

                      }
    );
#else

#endif
    TPZFNMatrix<6,REAL> strain(6, 1, 0.);
    TPZFNMatrix<6,REAL> elastic_strain(6, 1, 0.);
    TPZFNMatrix<6,REAL> sigma(6, 1, 0.);
    TPZFNMatrix<6,REAL> plastic_strain(6, 1, 0.);
    TPZFNMatrix<3,REAL> aux_tensor(3, 1, 0.);


    // Lambda for evaluate flux, this is supposed to be implemented in GPU
    auto EvaluateFlux = [this, mat, & glob_delta_strain, & glob_sigma, & strain, & elastic_strain, & plastic_strain, & sigma, & aux_tensor] (int & ipts)
    {
        REAL lambda = fMaterial->GetPlasticModel().fER.Lambda();
        REAL mu =  fMaterial->GetPlasticModel().fER.Mu();
        REAL mc_phi = fMaterial->GetPlasticModel().fYC.Phi();
        REAL mc_psi = fMaterial->GetPlasticModel().fYC.Psi();
        REAL mc_cohesion = fMaterial->GetPlasticModel().fYC.Cohesion();
        REAL G = mu;
        REAL K = lambda + 2 * mu/3;

        REAL alpha;
        int mtype;

        //Get from delta strain vector
        glob_delta_strain.GetSub(3 * ipts, 0, 3, 1, aux_tensor);

        //Get last strain vector
        fStrain.GetSub(6 * ipts, 0, 6, 1, strain);

        // Translate and
        ComposeStrain(&aux_tensor(0,0), &strain(0,0));

        //Get from plastic strain vector
        fPlasticStrain.GetSub(6 * ipts, 0, 6, 1, plastic_strain);

        TPZFMatrix<REAL> stress_eigenvalues(3, 1, 0.);
        TPZFMatrix<REAL> stress_eigenvectors(9, 1, 0.);
        TPZFMatrix<REAL> sigma_projected(3, 1, 0.);
        TPZFMatrix<REAL> gradient(9, 1, 0.);

        // Return Mapping components
        ElasticStrain(&plastic_strain(0,0), &strain(0,0), &elastic_strain(0,0));
        ComputeTrialStress(&elastic_strain(0,0), &sigma(0,0), mu, lambda);
        SpectralDecomposition(&sigma(0,0), &stress_eigenvalues(0,0), &stress_eigenvectors(0,0));
        ProjectSigma(&stress_eigenvalues(0,0), &sigma_projected(0,0), mtype, alpha, mc_phi, mc_psi, mc_cohesion, K, G, &gradient(0,0), false);
        ReconstructStressTensor(&sigma_projected(0,0), &stress_eigenvectors(0,0), &sigma(0,0));

        // Update plastic strain
        ComputeStrain(&sigma(0,0), &elastic_strain(0,0), mu, lambda);
        PlasticStrain(&strain(0,0), &elastic_strain(0,0), &plastic_strain(0,0));

        //Copy to stress vector
        TranslateStress(&sigma(0,0), &aux_tensor(0,0));

        aux_tensor *= fWeight[ipts];
        glob_sigma.PutSub(3 * ipts, 0, aux_tensor);

        if (mat->GetUpdateMem()) {

            fSigma.AddSub(6 * ipts, 0, sigma);

            //Accumulate to strain vector
            fStrain.AddSub(6 * ipts, 0, strain);

            //Copy to plastic strain vector
            fPlasticStrain.PutSub(6 * ipts, 0, plastic_strain);

            //Copy to MType and Alpha vectors
            fAlpha(ipts,0) = alpha;
            fMType(ipts,0) = mtype;
        }
    };

    for (int ipts = 0; ipts < fNpts; ipts++)
    {
        EvaluateFlux(ipts);
    }

    if (mat->GetUpdateMem()) {

#ifdef USING_TBB
        tbb::parallel_for(size_t(0), size_t(fNpts), size_t(1), [this, mat](size_t &ipts) {

                              TPZFNMatrix<6, REAL> strain(6, 1, 0.);
                              TPZFNMatrix<6, REAL> plastic_strain(6, 1, 0.);
                              TPZFNMatrix<6, REAL> sigma(6, 1, 0.);

                              TPZElastoPlasticMem &memory = mat->MemItem(ipts);
                              fSigma.GetSub(6 * ipts, 0, 6, 1, sigma);
                              fStrain.GetSub(6 * ipts, 0, 6, 1, strain);
                              fPlasticStrain.GetSub(6 * ipts, 0, 6, 1, plastic_strain);
                              memory.m_sigma.CopyFrom(sigma);
                              memory.m_elastoplastic_state.m_eps_t.CopyFrom(strain);
                              memory.m_elastoplastic_state.m_eps_p.CopyFrom(plastic_strain);
                              memory.m_elastoplastic_state.m_hardening = fAlpha(ipts, 0);
                              memory.m_elastoplastic_state.m_m_type = fMType(ipts, 0);

                          }
        );
#else
        // Lambda expression to update memory
        auto UpdateElastoPlasticState = [this, mat, & strain, & plastic_strain, & sigma] (int ipts)
        {
            
            TPZElastoPlasticMem & memory = mat->MemItem(ipts);
            fSigma.GetSub(6 * ipts, 0, 6, 1, sigma);
            fStrain.GetSub(6 * ipts, 0, 6, 1, strain);
            fPlasticStrain.GetSub(6 * ipts, 0, 6, 1, plastic_strain);
            memory.m_sigma.CopyFrom(sigma);
            memory.m_elastoplastic_state.m_eps_t.CopyFrom(strain);
            memory.m_elastoplastic_state.m_eps_p.CopyFrom(plastic_strain);
            memory.m_elastoplastic_state.m_hardening = fAlpha(ipts,0);
            memory.m_elastoplastic_state.m_m_type = fMType(ipts,0);
            
        };
        
        for (int ipts = 0; ipts < fNpts; ipts++)
        {
            UpdateElastoPlasticState(ipts);
            
        }

#endif
    }

#endif
}

void TPZConstitutiveLawProcessor::ComputeSigmaDep(TPZFMatrix<REAL> &glob_delta_strain, TPZFMatrix<REAL> &glob_sigma, TPZFMatrix<REAL> &glob_dep) {

#ifndef USING_CUDA

    int64_t rows = glob_delta_strain.Rows();
    int64_t cols = glob_delta_strain.Cols();
    glob_sigma.Resize(rows, cols);
    glob_dep.Resize( 3 * 3 * fNpts, 1);

    // Variables required for lambda's capture block
    TPZMatWithMem<TPZElastoPlasticMem> *mat = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *>(fMaterial);


    // The constitutive law is computing assuming full tensors
#ifdef USING_TBB

    tbb::parallel_for(size_t(0), size_t(fNpts), size_t(1),
                      [this, mat, & glob_delta_strain, & glob_sigma](size_t &ipts) {

                          TPZFNMatrix<6, REAL> strain(6, 1, 0.);
                          TPZFNMatrix<6, REAL> elastic_strain(6, 1, 0.);
                          TPZFNMatrix<6, REAL> sigma(6, 1, 0.);
                          TPZFNMatrix<6, REAL> plastic_strain(6, 1, 0.);
                          TPZFNMatrix<3, REAL> aux_tensor(3, 1, 0.);

                          REAL lambda = fMaterial->GetPlasticModel().fER.Lambda();
                          REAL mu = fMaterial->GetPlasticModel().fER.Mu();
                          REAL mc_phi = fMaterial->GetPlasticModel().fYC.Phi();
                          REAL mc_psi = fMaterial->GetPlasticModel().fYC.Psi();
                          REAL mc_cohesion = fMaterial->GetPlasticModel().fYC.Cohesion();
                          REAL G = mu;
                          REAL K = lambda + 2 * mu / 3;


                          REAL alpha;
                          int mtype;

                          //Get from delta strain vector
                          glob_delta_strain.GetSub(3 * ipts, 0, 3, 1, aux_tensor);

                          //Get last strain vector
                          fStrain.GetSub(6 * ipts, 0, 6, 1, strain);

                          // Translate and
                          ComposeStrain(&aux_tensor(0, 0), &strain(0, 0));

                          //Get from plastic strain vector
                          fPlasticStrain.GetSub(6 * ipts, 0, 6, 1, plastic_strain);

                          TPZFMatrix<REAL> stress_eigenvectors(9, 1, 0.);
                          TPZFMatrix<REAL> sigma_projected(3, 1, 0.);
                          TPZFMatrix<REAL> gradient(9, 1, 0.);
                          TPZFMatrix<REAL> tangent(36, 1, 0.);

                          TPZFMatrix<REAL> strain_eigenvalues(3, 1, 0.);
                          TPZFMatrix<REAL> strain_eigenvectors(36, 1, 0.);

                          // Return Mapping components
//                          ElasticStrain(&plastic_strain(0, 0), &strain(0, 0), &elastic_strain(0, 0));
//                          ComputeTrialStress(&elastic_strain(0, 0), &sigma(0, 0), mu, lambda);
//                          SpectralDecomposition(&sigma(0, 0), &sigma_projected(0, 0), &stress_eigenvectors(0, 0));
//                          SpectralDecomposition(&strain(0, 0), &strain_eigenvalues(0, 0), &strain_eigenvectors(0, 0));
//                          ProjectSigma(&sigma_projected(0, 0), &sigma_projected(0, 0), mtype, alpha, mc_phi, mc_psi,
//                                       mc_cohesion, K, G);
//                          TangentOperator(&gradient(0, 0), &sigma_projected(0, 0), &strain_eigenvalues(0, 0),
//                                          &strain_eigenvectors(0, 0), &tangent(0, 0), G, lambda);
//                          ReconstructStressTensor(&sigma_projected(0, 0), &stress_eigenvectors(0, 0), &sigma(0, 0));

                          // Update plastic strain
                          ComputeStrain(&sigma(0, 0), &elastic_strain(0, 0), mu, lambda);
                          PlasticStrain(&strain(0, 0), &elastic_strain(0, 0), &plastic_strain(0, 0));

                          //Copy to stress vector
                          TranslateStress(&sigma(0, 0), &aux_tensor(0, 0));

                          aux_tensor *= fWeight[ipts];
                          glob_sigma.PutSub(3 * ipts, 0, aux_tensor);

                          if (mat->GetUpdateMem()) {

                              fSigma.AddSub(6 * ipts, 0, sigma);

                              //Accumulate to strain vector
                              fStrain.AddSub(6 * ipts, 0, strain);

                              //Copy to plastic strain vector
                              fPlasticStrain.PutSub(6 * ipts, 0, plastic_strain);

                              //Copy to MType and Alpha vectors
                              fAlpha(ipts, 0) = alpha;
                              fMType(ipts, 0) = mtype;
                          }

                      }
    );
#else

#endif
    TPZFNMatrix<6,REAL> strain(6, 1, 0.);
    TPZFNMatrix<6,REAL> elastic_strain(6, 1, 0.);
    TPZFNMatrix<6,REAL> sigma(6, 1, 0.);
    TPZFNMatrix<6,REAL> plastic_strain(6, 1, 0.);
    TPZFNMatrix<3,REAL> aux_tensor(3, 1, 0.);

    // Lambda for evaluate flux, this is supposed to be implemented in GPU
    auto EvaluateFlux = [this, mat, & glob_delta_strain, & glob_sigma, & glob_dep, & strain, & elastic_strain, & plastic_strain, & sigma, & aux_tensor] (int & ipts)
    {
        REAL lambda = fMaterial->GetPlasticModel().fER.Lambda();
        REAL mu =  fMaterial->GetPlasticModel().fER.Mu();
        REAL mc_phi = fMaterial->GetPlasticModel().fYC.Phi();
        REAL mc_psi = fMaterial->GetPlasticModel().fYC.Psi();
        REAL mc_cohesion = fMaterial->GetPlasticModel().fYC.Cohesion();
        REAL G = mu;
        REAL K = lambda + 2 * mu/3;

        REAL alpha;
        int mtype;

        //Get from delta strain vector
        glob_delta_strain.GetSub(3 * ipts, 0, 3, 1, aux_tensor);

        //Get last strain vector
        fStrain.GetSub(6 * ipts, 0, 6, 1, strain);

        // Translate and
        ComposeStrain(&aux_tensor(0,0), &strain(0,0));

        //Get from plastic strain vector
        fPlasticStrain.GetSub(6 * ipts, 0, 6, 1, plastic_strain);

        TPZFNMatrix<3,REAL> stress_eigenvalues(3, 1, 0.);
        TPZFNMatrix<9,REAL> stress_eigenvectors(9, 1, 0.);
        TPZFNMatrix<3,REAL> sigma_projected(3, 1, 0.);
        TPZFNMatrix<9,REAL> gradient(9, 1, 0.);
        TPZFNMatrix<36,REAL> dep(36, 1, 0.);

        // Return Mapping components
        ElasticStrain(&plastic_strain(0,0), &strain(0,0), &elastic_strain(0,0));
        ComputeTrialStress(&elastic_strain(0,0), &sigma(0,0), mu, lambda);
        SpectralDecomposition(&sigma(0,0), &stress_eigenvalues(0,0), &stress_eigenvectors(0,0));
        TPZFNMatrix<3,REAL> strain_eigenvalues(3, 1, 0.);
        TPZFNMatrix<9,REAL> strain_eigenvectors(9, 1, 0.);

        elastic_strain(1,0) = elastic_strain(1,0) / 2;
        SpectralDecomposition(&elastic_strain(0,0), &strain_eigenvalues(0,0), &strain_eigenvectors(0,0));
        ProjectSigma(&stress_eigenvalues(0,0), &sigma_projected(0,0), mtype, alpha, mc_phi, mc_psi, mc_cohesion, K, G, &gradient(0,0), true);
        gradient.Zero();
        gradient(0,0) = 1.;
        gradient(1,1) = 1.;
        gradient(2,2) = 1.;
        TangentOperator(&gradient(0,0), &stress_eigenvalues(0,0), &strain_eigenvalues(0,0), &strain_eigenvectors(0,0), &dep(0,0), G, lambda);

        ReconstructStressTensor(&sigma_projected(0,0), &stress_eigenvectors(0,0), &sigma(0,0));

        // Update plastic strain
        ComputeStrain(&sigma(0,0), &elastic_strain(0,0), mu, lambda);
        PlasticStrain(&strain(0,0), &elastic_strain(0,0), &plastic_strain(0,0));

        //Copy to stress vector
        TranslateStress(&sigma(0,0), &aux_tensor(0,0));

        //dep = full_dep, gradient = dep
        TranslateDep(&dep(0,0),  &gradient(0,0));

        aux_tensor *= fWeight[ipts];
        glob_sigma.PutSub(3 * ipts, 0, aux_tensor);

//        gradient *= fWeight[ipts];
        glob_dep.PutSub(3 * 3 * ipts, 0, gradient);

        if (mat->GetUpdateMem()) {

            fSigma.AddSub(6 * ipts, 0, sigma);

            //Accumulate to strain vector
            fStrain.AddSub(6 * ipts, 0, strain);

            //Copy to plastic strain vector
            fPlasticStrain.PutSub(6 * ipts, 0, plastic_strain);

            //Copy to MType and Alpha vectors
            fAlpha(ipts,0) = alpha;
            fMType(ipts,0) = mtype;
        }
    };


    for (int ipts = 0; ipts < fNpts; ipts++)
    {
        EvaluateFlux(ipts);
    }

    if (mat->GetUpdateMem()) {

#ifdef USING_TBB
        tbb::parallel_for(size_t(0), size_t(fNpts), size_t(1), [this, mat](size_t &ipts) {

                              TPZFNMatrix<6, REAL> strain(6, 1, 0.);
                              TPZFNMatrix<6, REAL> plastic_strain(6, 1, 0.);
                              TPZFNMatrix<6, REAL> sigma(6, 1, 0.);

                              TPZElastoPlasticMem &memory = mat->MemItem(ipts);
                              fSigma.GetSub(6 * ipts, 0, 6, 1, sigma);
                              fStrain.GetSub(6 * ipts, 0, 6, 1, strain);
                              fPlasticStrain.GetSub(6 * ipts, 0, 6, 1, plastic_strain);
                              memory.m_sigma.CopyFrom(sigma);
                              memory.m_elastoplastic_state.m_eps_t.CopyFrom(strain);
                              memory.m_elastoplastic_state.m_eps_p.CopyFrom(plastic_strain);
                              memory.m_elastoplastic_state.m_hardening = fAlpha(ipts, 0);
                              memory.m_elastoplastic_state.m_m_type = fMType(ipts, 0);

                          }
        );
#else
        // Lambda expression to update memory
        auto UpdateElastoPlasticState = [this, mat, & strain, & plastic_strain, & sigma] (int ipts)
        {

            TPZElastoPlasticMem & memory = mat->MemItem(ipts);
            fSigma.GetSub(6 * ipts, 0, 6, 1, sigma);
            fStrain.GetSub(6 * ipts, 0, 6, 1, strain);
            fPlasticStrain.GetSub(6 * ipts, 0, 6, 1, plastic_strain);
            memory.m_sigma.CopyFrom(sigma);
            memory.m_elastoplastic_state.m_eps_t.CopyFrom(strain);
            memory.m_elastoplastic_state.m_eps_p.CopyFrom(plastic_strain);
            memory.m_elastoplastic_state.m_hardening = fAlpha(ipts,0);
            memory.m_elastoplastic_state.m_m_type = fMType(ipts,0);

        };

        for (int ipts = 0; ipts < fNpts; ipts++)
        {
            UpdateElastoPlasticState(ipts);

        }
#endif
    }

#endif
}

#ifdef USING_CUDA
void TPZConstitutiveLawProcessor::ComputeSigma(TPZVecGPU<REAL> &glob_delta_strain, TPZVecGPU<REAL> &glob_sigma) {
    REAL lambda = fMaterial->GetPlasticModel().fER.Lambda();
    REAL mu =  fMaterial->GetPlasticModel().fER.Mu();
    REAL mc_phi = fMaterial->GetPlasticModel().fYC.Phi();
    REAL mc_psi = fMaterial->GetPlasticModel().fYC.Psi();
    REAL mc_cohesion = fMaterial->GetPlasticModel().fYC.Cohesion();

    int64_t rows = glob_delta_strain.getSize();
    glob_sigma.resize(rows);

    TPZMatWithMem<TPZElastoPlasticMem> * mat = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *>(fMaterial);

    fCudaCalls->ComputeSigma(mat->GetUpdateMem(), fNpts, glob_delta_strain.getData(), glob_sigma.getData(), lambda, mu, mc_phi, mc_psi, mc_cohesion, dPlasticStrain.getData(), 
        dMType.getData(), dAlpha.getData(), dSigma.getData(), dStrain.getData(), dWeight.getData());

        if (mat->GetUpdateMem()) {
            dPlasticStrain.get(&fPlasticStrain(0,0), 6 * fNpts);
            dSigma.get(&fSigma(0,0), 6 * fNpts);
            dStrain.get(&fStrain(0,0), 6 * fNpts);
            dMType.get(&fMType(0,0), 6 * fNpts);
            dAlpha.get(&fAlpha(0,0), 6 * fNpts);

        
#ifdef USING_TBB
        tbb::parallel_for(size_t(0),size_t(fNpts),size_t(1), [this, mat] (size_t & ipts){
            
            TPZFNMatrix<6,REAL> strain(6, 1, 0.);
            TPZFNMatrix<6,REAL> plastic_strain(6, 1, 0.);
            TPZFNMatrix<6,REAL> sigma(6, 1, 0.);
            
            TPZElastoPlasticMem & memory = mat->MemItem(ipts);
            fSigma.GetSub(6 * ipts, 0, 6, 1, sigma);
            fStrain.GetSub(6 * ipts, 0, 6, 1, strain);
            fPlasticStrain.GetSub(6 * ipts, 0, 6, 1, plastic_strain);
            memory.m_sigma.CopyFrom(sigma);
            memory.m_elastoplastic_state.m_eps_t.CopyFrom(strain);
            memory.m_elastoplastic_state.m_eps_p.CopyFrom(plastic_strain);
            memory.m_elastoplastic_state.m_hardening = fAlpha(ipts,0);
            memory.m_elastoplastic_state.m_m_type = fMType(ipts,0);
            
        }
                          );
#else
        TPZFNMatrix<6,REAL> strain(6, 1, 0.);
        TPZFNMatrix<6,REAL> plastic_strain(6, 1, 0.);
        TPZFNMatrix<6,REAL> sigma(6, 1, 0.);
        // Lambda expression to update memory
        auto UpdateElastoPlasticState = [this, mat, & strain, & plastic_strain, & sigma] (int ipts)
        {
            
            TPZElastoPlasticMem & memory = mat->MemItem(ipts);
            fSigma.GetSub(6 * ipts, 0, 6, 1, sigma);
            fStrain.GetSub(6 * ipts, 0, 6, 1, strain);
            fPlasticStrain.GetSub(6 * ipts, 0, 6, 1, plastic_strain);
            memory.m_sigma.CopyFrom(sigma);
            memory.m_elastoplastic_state.m_eps_t.CopyFrom(strain);
            memory.m_elastoplastic_state.m_eps_p.CopyFrom(plastic_strain);
            memory.m_elastoplastic_state.m_hardening = fAlpha(ipts,0);
            memory.m_elastoplastic_state.m_m_type = fMType(ipts,0);
            
        };
        
        for (int ipts = 0; ipts < fNpts; ipts++)
        {
            UpdateElastoPlasticState(ipts);
            
        }
        
#endif
    }

}
#endif

#ifdef USING_CUDA
void TPZConstitutiveLawProcessor::TransferDataToGPU() {
    dWeight.resize(fWeight.size());
    dWeight.set(&fWeight[0], fWeight.size());
}
#endif

#ifdef USING_CUDA
void TPZConstitutiveLawProcessor::De() {
    REAL lambda = fMaterial->GetPlasticModel().fER.Lambda();
    REAL mu = fMaterial->GetPlasticModel().fER.G();
    fCudaCalls->DeToDevice(lambda, mu); 
}
#endif
