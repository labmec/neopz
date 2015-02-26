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
    
    /** @brief State: n or n+1 temporal state */
    enum nState { n = 0, nplusone = 1 };
    nState State;
    
private:
    
    /**
     * @ingroup Simulation Parameters
     * @brief Define simulation parameters for Transient.
     * @since December 08, 2014
    */
    
	/** @brief Delta t - s */
	REAL fDeltaT;
	
	/** @brief Time - s */
	REAL fTime;

    /** @brief Maximum number of newton iterations */
    int fMaxiterations;
    
	/** @brief DeltaX tolerance for newton iterations */
	REAL ftoleranceDeltaX;
    
    /** @brief Residual tolerance for newton iterations */
    REAL ftoleranceResiual;
    
    /** @brief Broyden iterations */
    bool fIsBroyden;

	
public:

    /** @brief Set Time step - s */
    void SetDeltaT(REAL DeltaT){this->fDeltaT = DeltaT;}
    
    /** @brief Get Time step - s */
    REAL GetDeltaT() const {return this->fDeltaT;}
    
    /** @brief Set Time - s */
    void SetTime(REAL Time) {this->fTime = Time;}
    
    /** @brief Get Time - s */
    REAL GetTime() const {return this->fTime;}
    
    /** @brief Set Tolerance for  Delta X */
    void SetToleranceDX(REAL dxtol) {this->ftoleranceDeltaX = dxtol;}
    
    /** @brief Get Tolerance for  Delta X */
    REAL GetToleranceDX() const {return this->ftoleranceDeltaX;}
    
    /** @brief Set Tolerance for  Residual vector */
    void SetToleranceRes(REAL restol) {this->ftoleranceResiual = restol;}
    
    /** @brief Get Tolerance for  Residual vector */
    REAL GetToleranceRes() const {return this->ftoleranceResiual;}
    
    /** @brief Set Maximum newton iterations - - */
    void SetMaxiterations(int Maxiterations){this->fMaxiterations = Maxiterations;}
    
    /** @brief Get Maximum newton iterations - - */
    int GetMaxiterations() const {return this->fMaxiterations;}
    
    /** @brief Using Broyden iterations */
    void SetIsBroyden(bool Broyden) {fIsBroyden = Broyden;}
    
    /** @brief Using Broyden iterations */
    bool GetIsBroyden() {return fIsBroyden;}


    
};

#endif