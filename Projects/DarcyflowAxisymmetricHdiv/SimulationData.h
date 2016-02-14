#ifndef TPZSimulationDATAH
#define TPZSimulationDATAH
/*
 *  SimulationData.h
 *  PZ
 *
 *  Created by Omar Duran Triana on 5/21/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */


#include "tpzautopointer.h"
#include "pzfmatrix.h"
#include <math.h>


class SimulationData {
    
public:
    
    /** @brief Default Constructor */
    SimulationData();
    
    /** @brief Destructor */
    ~SimulationData();
    
    /** @brief copy constructor */
    SimulationData(const SimulationData &copy);
    
    /** @brief operator equal */
    SimulationData &operator=(const SimulationData &copy);
    
private:
    
    /**
     * @ingroup Simulation Parameters
     * @brief Define simulation parameters for Transient.
     * @since December 08, 2014
     */
    
    /** @brief Delta time - s */
    REAL fDeltaT;
    
    /** @brief Max time - s */
    REAL fMaxTime;
    
    /** @brief Time - s */
    REAL fTime;
    
    /** @brief Maximum number of newton iterations */
    int fMaxiterations;
    
    /** @brief Number of threads for assembly */
    int fnthreads;
    
    /** @brief Number of fixed jacobian */
    int fFixedJacobianIterations;
    
    /** @brief DeltaX tolerance for newton iterations */
    REAL ftoleranceDeltaX;
    
    /** @brief Residual tolerance for newton iterations */
    REAL ftoleranceResiual;
    
    /** @brief Number of uniform mesh refinement */
    int fHref;
    
    /** @brief Number of uniform mesh refinement in postprocessing */
    int fHrefpost;
    
    /** @brief Approximation order for velocity */
    int fqorder;
    
    /** @brief Approximation order for pressure */
    int fporder;
    
    /** @brief Approximation order for saturations */
    int fsorder;
    
    /** @brief Use of direct solver */
    bool fIsDirect;
    
    /** @brief Optimize band width */
    bool fOptband;
    
    /** @brief Use of Conjugated Gradient method */
    bool fIsCG;
    
    /** @brief Picard iterations */
    bool fIsPicard;
    
    /** @brief Define the use of a linear approxiamtion of S using a gradient reconstruction procedure */
    bool fUseGR;
    
    /** @brief Define the use of static condensation */
    bool fCondenseElements;
    
    /** @brief Define the use dimensionless formulation */
    bool fIsDimensionless;
    
    /** @brief State: n or n+1 temporal state */
    bool fnStep;
    
    /** @brief Number of elements, x direction */
    int fNelmx;
    
    /** @brief Number of elements, y direction  */
    int fNelmy;
    
    /** @brief Length of elements, x direction */
    REAL fLengthElementx;
    
    /** @brief Length of elements, y direction */
    REAL fLengthElementy;
    
    /** @brief Definition of the flow system one - two and  ... three phase */
    TPZStack<std::string> fSystemType;
    
    /** @brief vector that contains materials id to integrate the production */
    TPZStack<int> fMaterialsToIntegrate;
    
    /** @brief Is one-phase flow? */
    bool fIsOnePhaseQ;
    
    /** @brief Is two-phase flow? */
    bool fIsTwoPhaseQ;
    
    /** @brief Is three-phase flow? */
    bool fIsThreePhaseQ;
    
    /** @brief Using the Linear grativational segregation model */
    bool fIsLinearSegregationsQ;
    
    /** @brief Is axisymmetric analysis */
    bool fIsAxisymmetricQ;
    
    /** @brief Is Impes analysis */
    bool fIsImpesQ;
    
    /** @brief Is hydrostatic boundary condition */
    bool fIsHydrostaticBCQ;
    
    /** @brief Gravity value */
    TPZFMatrix<REAL> fGravity;
    
    /** @brief Counterclockwise rotation angle [degrees] */
    REAL fAngle;
    
    /** @brief Store time values to be reported */
    TPZManVector<REAL> fTimesToPrint;

    /** @brief Number of sub times steps per dt */
    int fn_sub_dt;
    
    /** @brief Is triangular mesh */
    bool fIsTriangularMeshQ;
    
    /** @brief Is a mesh with geometric progression */
    bool fIsMeshwithPGQ;
    
    /** @brief geometric progression ratio <= 1 */
    REAL fpg_ratio;
    
    /** @brief Compute a variable K distribution */
    bool fIsHeterogeneousQ;
    
    /** @brief Ouput directory name */
    std::string foutdirectory;
    
    /** @brief GID mesh file */
    std::string fGIDfile;
    
    /** @brief Time scale for dimensionless calculations */
    REAL ftime_scale;
    
    /** @brief Length scale for dimensionless calculations */
    REAL flength_scale;
    
    /** @brief Velocity scale for dimensionless calculations */
    REAL fvelocity_scale;
    
    
public:
    
    /** @brief void material being used for GR */
    int fMatL2;
    
    /** @brief Set Time step - s */
    void SetDeltaT(REAL DeltaT){this->fDeltaT = DeltaT;}
    
    /** @brief Set number of elements, x direction -  */
    void SetnElementsx (int Nelem){this-> fNelmx = Nelem;}
    /** @brief Get number of elements, x direction-  */
    int GetnElementsx (){return this-> fNelmx;}
    
    
    /** @brief Set number of elements, y direction -  */
    void SetnElementsy (int Nelem){this-> fNelmy = Nelem;}
    /** @brief Get number of elements, y direction -  */
    int GetnElementsy (){return this-> fNelmy;}
    
    /** @brief Set Length of elements  x direction-  */
    void SetLengthElementx (REAL LengthElementx){this-> fLengthElementx = LengthElementx;}
    /** @brief Get Length of elements  x direction- -  */
    REAL GetLengthElementx (){return this-> fLengthElementx;}
    
    /** @brief Set Length of elements  y direction -  */
    void SetLengthElementy (REAL LengthElementy){this-> fLengthElementy = LengthElementy;}
    /** @brief Get Length of elements  y direction -  */
    REAL GetLengthElementy (){return this-> fLengthElementy;}
    
    /** @brief Get Time step - s */
    REAL GetDeltaT() {return this->fDeltaT;}
    
    /** @brief Set Time - s */
    void SetTime(REAL Time) {this->fTime = Time;}
    
    /** @brief Get Time - s */
    REAL GetTime() {return this->fTime;}
    
    /** @brief Set maximum time - s */
    void SetMaxTime(REAL MaxTime) {this->fMaxTime = MaxTime;}
    
    /** @brief Get maximum time - s */
    REAL GetMaxTime() {return this->fMaxTime;}
    
    /** @brief Set tolerance for  delta X */
    void SetToleranceDX(REAL dxtol) {this->ftoleranceDeltaX = dxtol;}
    
    /** @brief Get tolerance for delta X */
    REAL GetToleranceDX() {return this->ftoleranceDeltaX;}
    
    /** @brief Set Tolerance for  Residual vector */
    void SetToleranceRes(REAL restol) {this->ftoleranceResiual = restol;}
    
    /** @brief Get Tolerance for  Residual vector */
    REAL GetToleranceRes() {return this->ftoleranceResiual;}
    
    /** @brief Set Maximum newton iterations - - */
    void SetMaxiterations(int Maxiterations){this->fMaxiterations = Maxiterations;}
    
    /** @brief Get Maximum newton iterations - - */
    int GetMaxiterations() {return this->fMaxiterations;}
    
    /** @brief Set Number of threads for assembly */
    void SetNthreads(int ntrheads){this->fnthreads = ntrheads;}
    
    /** @brief Get Number of threads for assembly */
    int GetNthreads() {return this->fnthreads;}
    
    /** @brief Set number of fixed jacobian for newton iterations - - */
    void SetFixediterations(int Fixediterations){this->fFixedJacobianIterations = Fixediterations;}
    
    /** @brief Get Maximum newton fixed jacobian iterations - - */
    int GetFixediterations() {return this->fFixedJacobianIterations;}
    
    /** @brief Set the Number of uniform mesh refinement */
    void SetHrefinement(int h){this->fHref = h;}
    
    /** @brief Get the Number of uniform mesh refinement */
    int GetHrefinement() {return this->fHref;}
    
    /** @brief Set the Number of uniform mesh refinement in postprocessing */
    void SetHPostrefinement(int h){this->fHrefpost = h;}
    
    /** @brief Get the Number of uniform mesh refinement in postprocessing */
    int GetHPostrefinement() {return this->fHrefpost;}
    
    /** @brief Set the approximation order for flux */
    void Setqorder(int qp){this->fqorder = qp;}
    
    /** @brief Get the approximation order for flux */
    int Getqorder() {return this->fqorder;}
    
    /** @brief Set the approximation order for pressure */
    void Setporder(int pp){this->fporder = pp;}
    
    /** @brief Get the approximation order for pressure */
    int Getporder() {return this->fporder;}
    
    /** @brief Set the approximation order for saturations */
    void Setsorder(int pp){this->fsorder = pp;}
    
    /** @brief Get the approximation order for saturations */
    int Getsorder() {return this->fsorder;}
        
    /** @brief Get using Broyden iterations */
    bool GetIsPicard() {return fIsPicard;}
    
    /** @brief Set the use of direct Solver */
    void SetIsDirect(bool isdirect) {fIsDirect = isdirect;}
    
    /** @brief Get the use of direct Solver */
    bool GetIsDirect() {return fIsDirect;}
    
    /** @brief Set the use of band optimization */
    void SetOptband(bool Optband) {fOptband = Optband;}
    
    /** @brief Get the use of band optimization */
    bool GetOptband() {return fOptband;}
    
    /** @brief Set the use of CG method (false is GMRES) */
    void SetIsCG(bool IsCG) {fIsCG = IsCG;}
    
    /** @brief Get the use of CG method (false is GMRES) */
    bool GetIsCG() {return fIsCG;}
    
    /** @brief Set the use of CG method (false is GMRES) */
    void SetnStep(bool nstep) {fnStep = nstep;}
    
    /** @brief Get the use of CG method (false is GMRES) */
    bool IsnStep() {return fnStep;}
    
    /** @brief Set the use of GR iterations */
    void SetGR(bool GR) {fUseGR = GR;}
    
    /** @brief Set the use of GR iterations */
    bool GetGR() {return fUseGR;}
    
    /** @brief Set the use of Static condensation */
    void SetSC(bool SC) {fCondenseElements = SC;}
    
    /** @brief Get the use of Static condensation */
    bool GetSC() {return fCondenseElements;}
    
    /** @brief Set the use of dimensionless formulation */
    void SetIsDimensionless(bool isdimensionless) {fIsDimensionless = isdimensionless;}
    
    /** @brief Get the use of dimensionless formulation */
    bool GetIsDimensionless() {return fIsDimensionless;}
    
    
    /** @brief Set gravity value */
    void SetGravity(TPZFMatrix< REAL > &gravity){
        
        fGravity = gravity;
    }
    
    /** @brief Get gravty value*/
    TPZFMatrix<REAL> GetGravity(){return fGravity;}
    
    
    /** @brief Definition of the flow system one - two and  ... three phase */
    void SetsystemType(TPZStack<std::string> SystemType){
        
        switch (SystemType.size()) {
            case 1:
            {
                fIsOnePhaseQ = true;
            }
                break;
            case 2:
            {
                fIsTwoPhaseQ = true;

            }
                break;
            case 3:
            {
                fIsThreePhaseQ = true;
            }
                break;
                
            default:
            {
                std::cout << "This code run just three-phasic systems" << std::endl;
                DebugStop();
            }
                break;
        }
        fSystemType = SystemType;
    }
    
    /** @brief Set the vector that contains materials id to integrate the production */
    void SetMaterialsToIntegrate(TPZStack<int> MaterialsToIntegrate){
        this->fMaterialsToIntegrate = MaterialsToIntegrate;
    }
    
    /** @brief Set the vector that contains materials id to integrate the production */
    TPZStack<int> MaterialsToIntegrate(){
        return this->fMaterialsToIntegrate;
    }
    
    
    /** @brief Mono-phasic system */
    bool IsOnePhaseQ() {return fIsOnePhaseQ;}
    
    /** @brief Two-phasic system */
    bool IsTwoPhaseQ() {return fIsTwoPhaseQ;}
    
    /** @brief Three-phasic system */
    bool IsThreePhaseQ() {return fIsThreePhaseQ;}
    
    /** @brief Set the use of the Linear grativational segregation model */
    void SetLinearSegregationQ(bool linearsegregation) {fIsLinearSegregationsQ = linearsegregation;}
    
    /** @brief Is it using the Linear grativational segregation model */
    bool IsLinearSegregationQ() {return fIsLinearSegregationsQ;}
    
    /** @brief Set the use of axisymmetric analysis */
    void SetAxisymmetricQ(bool IsAxisymmetric) {fIsAxisymmetricQ = IsAxisymmetric;}
    
    /** @brief Is it using axisymmetric analysis */
    bool IsAxisymmetricQ() {return fIsAxisymmetricQ;}

    /** @brief Set the use of hydrostatic boundary condition */
    void SetHydrostaticBCQ(bool HydrostaticBC) {fIsHydrostaticBCQ = HydrostaticBC;}
    
    /** @brief Is it using hydrostatic boundary condition */
    bool IsHydrostaticBCQ() {return fIsHydrostaticBCQ;}
    
    /** @brief Set the use of axisymmetric analysis */
    void SetImpesQ(bool Impes) {fIsImpesQ = Impes;
        if(Impes) {
            fIsPicard = true;
        }
    }
    
    /** @brief Is it using axisymmetric analysis */
    bool IsImpesQ() {return fIsImpesQ;}
    
    /** @brief Definition of the flow system one - two and  ... three phase */
    TPZStack<std::string> GetsystemType() {return fSystemType;}
    
    /** @brief Set counterclockwise rotation angle [degrees] */
    void SetRotationAngle(REAL angle) {this->fAngle = angle;}
    
    /** @brief Get counterclockwise rotation angle [degrees] */
    REAL GetRotationAngle() {return this->fAngle;}
    
    /** @brief Set the use of triangular mesh */
    void SetTriangularMesh(bool IsTriangularMesh) {fIsTriangularMeshQ = IsTriangularMesh;}
    
    /** @brief Is it using a triangular mesh */
    bool IsTriangularMeshQ() {return fIsTriangularMeshQ;}
    
    /** @brief Set the use of mesh with geometric progression */
    void SetMeshwithPGQ(bool pg) {this->fIsMeshwithPGQ = pg;}
    
    /** @brief Is a mesh with geometric progression */
    bool IsMeshwithPGQ() {return fIsMeshwithPGQ;}

    /** @brief Set the use of heterogeneous k */
    void SetHeterogeneousQ(bool pg) {this->fIsHeterogeneousQ = pg;}
    
    /** @brief Is a heterogeneous k */
    bool IsHeterogeneousQ() {return fIsHeterogeneousQ;}
    
    /** @brief Set geometric progression ratio */
    void SetPGRatio(REAL ratio) {this->fpg_ratio = ratio;}
    
    /** @brief Set geometric progression ratio */
    REAL GetPGRatio() {return this->fpg_ratio;}
    
    /**
     * Set time values to be reported
     */
    void SetTimes(TPZManVector<REAL> TimesToPrint){
        fTimesToPrint = TimesToPrint;
    }
    
    /**
     * Get time values to be reported
     */
    TPZManVector<REAL>GetTimes(){
        return fTimesToPrint;
    }
    
    /** @brief Set the number of sub times steps per dt */
    void SetNSubSteps(int n_sub_dt){
        fn_sub_dt = n_sub_dt;
    }
    
    /** @brief Get the number of sub times steps per dt */
    int GetNSubSteps(){
        return fn_sub_dt;
    }
    
    /** @brief Set the ouput directory name */
    void SetDirectory( std::string outdirectory ) {foutdirectory = outdirectory; }
    
    /** @brief Get the ouput directory name */
    std::string GetDirectory() {
        return this->foutdirectory;
    }
    
    /** @brief Set GID file */
    void SetGIDFile( std::string gidfile ) {
        this->fGIDfile = gidfile;
    }
    
    /** @brief Get GID file */
    std::string GetGIDFile() {
        return this->fGIDfile;
    }


    /** @brief Set the time scale for dimensionless calculations */
    void SetTime_Scale(REAL time_scale) {ftime_scale = time_scale;}
    
    /** @brief Get the time scale for dimensionless calculations */
    REAL Time_Scale() {return ftime_scale;}
    
    /** @brief Set the length scale for dimensionless calculations */
    void SetLength_Scale(REAL length_scale) {flength_scale = length_scale;}
    
    /** @brief Get the length scale for dimensionless calculations */
    REAL Length_Scale() {return flength_scale;}
    
    /** @brief Set the velocity scale for dimensionless calculations */
    void SetVelocity_Scale(REAL velocity_scale) {fvelocity_scale = velocity_scale;}
    
    /** @brief Get the velocity scale for dimensionless calculations */
    REAL Velocity_Scale() {return fvelocity_scale;}
    
};

#endif