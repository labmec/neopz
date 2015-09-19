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
    
    /** @brief Delta t - s */
    REAL fDeltaT;
    
    /** @brief Max Time - s */
    REAL fMaxTime;    
    
    /** @brief Time - s */
    REAL fTime;
    
    /** @brief Maximum number of newton iterations */
    int fMaxiterations;
    
    /** @brief Maximum number of newton iterations */
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
    
    
    /** @brief Broyden iterations */
    bool fIsBroyden;
    
    /** @brief Define the use of a linear approxiamtion of S using a gradient reconstruction procedure */
    bool fUseGR;
    
    /** @brief Define the use of static condensation */
    bool fCondenseElements;
    
    /** @brief Define the use dimensionless formulation */
    bool fIsDimensionless;
    
    /** @brief State: n or n+1 temporal state */
    bool fnStep;
    
    /** @brief State: n or n+1 temporal state */
    int fNelmx;

    /** @brief State: n or n+1 temporal state */
    int fNelmy;
    
    /** @brief State: n or n+1 temporal state */
    REAL fLengthElementx;
    
    /** @brief State: n or n+1 temporal state */
    REAL fLengthElementy;
    
    /** @brief Definition of the Top bc */
    TPZVec<REAL> fTopBCini;
    
    /** @brief Definition of the Bottom bc */
    TPZVec<REAL> fBottomBCini;
    
    /** @brief Definition of the Right bc */
    TPZVec<REAL> fRightBCini;
    
    /** @brief Definition of the Left bc */
    TPZVec<REAL> fLeftBCini;
  
    /** @brief Definition of the Top bc */
    TPZVec<REAL> fTopBC;
    
    /** @brief Definition of the Bottom bc */
    TPZVec<REAL> fBottomBC;
    
    /** @brief Definition of the Right bc */
    TPZVec<REAL> fRightBC;
    
    /** @brief Definition of the Left bc */
    TPZVec<REAL> fLeftBC;
  
     /** @brief Gravity  */
    TPZFMatrix<REAL> fGravity;
    
public:
    
    /** @brief void material being used for GR */
    int fMatL2;
    
    /** @brief Set Time step - s */
    void SetDeltaT(REAL DeltaT){this->fDeltaT = DeltaT;}
    
    /** @brief Set nElements x -  */
    void SetnElementsx (int Nelem){this-> fNelmx = Nelem;}
      /** @brief Get nElements x -  */
    int GetnElementsx (){return this-> fNelmx;}
    
    
    /** @brief Set nElements x -  */
    void SetnElementsy (int Nelem){this-> fNelmy = Nelem;}
      /** @brief Get nElements x -  */
    int GetnElementsy (){return this-> fNelmy;}
    
    /** @brief Set LengthElementx  x -  */
    void SetLengthElementx (REAL LengthElementx){this-> fLengthElementx = LengthElementx;}
      /** @brief Get LengthElementx  x -  */
    REAL GetLengthElementx (){return this-> fLengthElementx;}

    /** @brief Set LengthElementy  y -  */
    void SetLengthElementy (REAL LengthElementy){this-> fLengthElementy = LengthElementy;}
      /** @brief Get LengthElementy  y -  */
    REAL GetLengthElementy (){return this-> fLengthElementy;}
    
    /** @brief Get Time step - s */
    REAL GetDeltaT() {return this->fDeltaT;}
    
    /** @brief Set Time - s */
    void SetTime(REAL Time) {this->fTime = Time;}
    
    /** @brief Get Time - s */
    REAL GetTime() {return this->fTime;}
    
    /** @brief Set Time - s */
    void SetMaxTime(REAL MaxTime) {this->fMaxTime = MaxTime;}
    
    /** @brief Get Time - s */
    REAL GetMaxTime() {return this->fMaxTime;}    
    
    /** @brief Set Tolerance for  Delta X */
    void SetToleranceDX(REAL dxtol) {this->ftoleranceDeltaX = dxtol;}
    
    /** @brief Get Tolerance for  Delta X */
    REAL GetToleranceDX() {return this->ftoleranceDeltaX;}
    
    /** @brief Set Tolerance for  Residual vector */
    void SetToleranceRes(REAL restol) {this->ftoleranceResiual = restol;}
    
    /** @brief Get Tolerance for  Residual vector */
    REAL GetToleranceRes() {return this->ftoleranceResiual;}
    
    /** @brief Set Maximum newton iterations - - */
    void SetMaxiterations(int Maxiterations){this->fMaxiterations = Maxiterations;}
    
    /** @brief Get Maximum newton iterations - - */
    int GetMaxiterations() {return this->fMaxiterations;}
    
    /** @brief Set Maximum newton iterations - - */
    void SetFixediterations(int Fixediterations){this->fFixedJacobianIterations = Fixediterations;}
    
    /** @brief Get Maximum newton iterations - - */
    int GetFixediterations() {return this->fFixedJacobianIterations;}
    
    /** @brief Set the Number of uniform mesh refinement */
    void SetHrefinement(int h){this->fHref = h;}
    
    /** @brief Set the Number of uniform mesh refinement */
    int GetHrefinement() {return this->fHref;}
    
    /** @brief Set the Number of uniform mesh refinement in postprocessing */
    void SetHPostrefinement(int h){this->fHrefpost = h;}
    
    /** @brief Set the Number of uniform mesh refinement in postprocessing */
    int GetHPostrefinement() {return this->fHrefpost;}
    
    /** @brief Set the approximation order for velocity */
    void Setqorder(int qp){this->fqorder = qp;}
    
    /** @brief Get the approximation order for velocity */
    int Getqorder() {return this->fqorder;}
    
    /** @brief Set the approximation order for pressure */
    void Setporder(int pp){this->fporder = pp;}
    
    /** @brief Get the approximation order for pressure */
    int Getporder() {return this->fporder;}
    
    /** @brief Set the approximation order for saturations */
    void Setsorder(int pp){this->fsorder = pp;}
    
    /** @brief Get the approximation order for saturations */
    int Getsorder() {return this->fsorder;}
    
    /** @brief Using Broyden iterations */
    void SetIsBroyden(bool Broyden) {fIsBroyden = Broyden;}
    
    /** @brief Using Broyden iterations */
    bool GetIsBroyden() {return fIsBroyden;}
    
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
    
    /** @brief Using GR iterations */
    void SetGR(bool GR) {fUseGR = GR;}
    
    /** @brief Using GR iterations */
    bool GetGR() {return fUseGR;}
    
    /** @brief Using Static condensation */
    void SetSC(bool SC) {fCondenseElements = SC;}
    
    /** @brief Using Static condensation */
    bool GetSC() {return fCondenseElements;}
    
    /** @brief Using dimensionless formulation */
    void SetIsDimensionless(bool isdimensionless) {fIsDimensionless = isdimensionless;}
    
    /** @brief Using dimensionless formulation */
    bool GetIsDimensionless() {return fIsDimensionless;}
    
    /** @brief Set Top bc */
    void SetTopBC(TPZVec<REAL> topbcini,TPZVec<REAL> topbc){
        
        if (topbc.size() != 4 && topbcini.size() != 4) {
            std::cout << "The number of parameter must to be equal 4, you give me = " << topbc.size() << std::endl;
            DebugStop();
        }
        fTopBCini = topbcini;
        fTopBC = topbc;
    }
    
    /** @brief Get Top bc for initial conditions */
    TPZVec<REAL> GetTopBCini(){return fTopBCini;}
    
    /** @brief Get Top bc */
    TPZVec<REAL> GetTopBC(){return fTopBC;}
    
    /** @brief Set Top bc */
    void SetBottomBC(TPZVec<REAL> bottombcini, TPZVec<REAL> bottombc){

        if (bottombc.size() != 4 && bottombcini.size() != 4) {
            std::cout << "The number of parameter must to be equal 4, you give me = " << bottombc.size() << std::endl;
            DebugStop();
        }
        fBottomBCini = bottombcini;
        fBottomBC = bottombc;
    }
    
    /** @brief Get Top bc */
    TPZVec<REAL> GetBottomBCini(){return fBottomBCini;}

    /** @brief Get Top bc */
    TPZVec<REAL> GetBottomBC(){return fBottomBC;}
    
    /** @brief Set Top bc */
    void SetRightBC(TPZVec<REAL> rightbcini, TPZVec<REAL> rightbc){
        
        if (rightbc.size() != 4 && rightbcini.size() != 4 ) {
            std::cout << "The number of parameter must to be equal 4, you give me = " << rightbc.size() << std::endl;
            DebugStop();
        }
        fRightBCini = rightbcini;
        fRightBC = rightbc;
    }
    
    /** @brief Get Top bc */
    TPZVec<REAL> GetRightBCini(){return fRightBCini;}
    
    /** @brief Get Top bc */
    TPZVec<REAL> GetRightBC(){return fRightBC;}
    
    /** @brief Set Top bc */
    void SetLeftBC(TPZVec<REAL> leftbcini, TPZVec<REAL> leftbc){
        
        if (leftbc.size() != 4 && leftbcini.size() != 4) {
            std::cout << "The number of parameter must to be equal 4, you give me = " << leftbc.size() << std::endl;
            DebugStop();
        }
        fLeftBCini = leftbcini;
        fLeftBC = leftbc;
    }
    
    /** @brief Get Top bc */
    TPZVec<REAL> GetLeftBCini(){return fLeftBCini;}
    
    /** @brief Get Top bc */
    TPZVec<REAL> GetLeftBC(){return fLeftBC;}
    
        /** @brief Set Gravity */
    void SetGravity(TPZFMatrix< REAL > &gravity){
        
        fGravity = gravity;
    }
    
    /** @brief Get Gravty */
    TPZFMatrix<REAL> GetGravity(){return fGravity;}
    
    
};

#endif