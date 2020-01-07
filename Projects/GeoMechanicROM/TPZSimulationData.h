//
//  TPZSimulationData.h
//  PZ
//
//  Created by Omar on 8/28/16.
//
//

#ifndef TPZSimulationData_h
#define TPZSimulationData_h

#include <stdio.h>
#include "pzvec.h"
#include "pzstack.h"

class TPZSimulationData {
    
protected:
    
    /** @brief check the full coupled residue */
    bool fCheckCouplingQ;
    
    /** @brief initial REAL */
    bool fIsInitialStateQ;
    
    /** @brief current time REAL */
    bool fIsCurrentStateQ;
    
    /** @brief Mixed approximation question */
    bool fIsMixedApproxQ;
    
    /** @brief Reduced basis approximation question */
    bool fIsRBApproxQ;
    
    /** @brief Definition gravity field */
    TPZVec<REAL> fg;
    
    /** @brief Store time values to be reported */
    TPZStack< REAL , 500 > fReportingTimes;
    
    /** @brief ntime steps */
    int fn_steps;
    
    /** @brief Time step */
    REAL fdt;
    
    /** @brief Time step */
    REAL ftime;
    
    /** @brief number of corrections steps */
    int fn_corrections;
    
    /** @brief residue overal tolerance */
    REAL fepsilon_res;
    
    /** @brief correction overal tolerance */
    REAL fepsilon_cor;
    
    /** @brief number of blocks */
    TPZStack< int > fblocks;
    
public:
    
    
    /** @brief default constructor */
    TPZSimulationData();
    
    /** @brief default constructor */
    TPZSimulationData(const TPZSimulationData &copy)
    {

        this->fIsInitialStateQ = copy.fIsCurrentStateQ;
        this->fIsCurrentStateQ = copy.fIsCurrentStateQ;
        this->fg = copy.fg;
        this->fReportingTimes = copy.fReportingTimes;
        this->fn_steps = copy.fn_steps;
        this->fdt = copy.fdt;
        this->ftime = copy.ftime;
        this->fn_corrections = copy.fn_corrections;
        this->fepsilon_res  = copy.fepsilon_res;
        this->fepsilon_cor  = copy.fepsilon_cor;
    }
    
    /** @brief default constructor */
    TPZSimulationData &operator=(const TPZSimulationData &copy)
    {
        this->fIsInitialStateQ = copy.fIsCurrentStateQ;
        this->fIsCurrentStateQ = copy.fIsCurrentStateQ;
        this->fg = copy.fg;
        this->fReportingTimes = copy.fReportingTimes;
        this->fn_steps = copy.fn_steps;
        this->fdt = copy.fdt;
        this->ftime = copy.ftime;
        this->fn_corrections = copy.fn_corrections;
        this->fepsilon_res  = copy.fepsilon_res;
        this->fepsilon_cor  = copy.fepsilon_cor;
        return *this;
    }
    
    /** @brief destructor */
    ~TPZSimulationData();
    


    /** @brief Set coupling check */
    void SetCheckCouplingQQ(bool state) { fCheckCouplingQ = state; }
    
    /** @brief Get coupling check */
    bool CheckCouplingQ() {return fCheckCouplingQ;}
    
    /** @brief Set initial REAL */
    void SetInitialStateQ(bool state) { fIsInitialStateQ = state; }
    
    /** @brief Get initial REAL */
    bool IsInitialStateQ() {return fIsInitialStateQ;}
    
    /** @brief current time REAL */
    void SetCurrentStateQ(bool state) { fIsCurrentStateQ = state; }
    
    /** @brief current time REAL */
    bool IsCurrentStateQ() {return fIsCurrentStateQ;}

    /** @brief Set Mixed basis approximation question */
    void SetMixedApproxQ(bool state) { fIsMixedApproxQ = state; }
    
    /** @brief Get Mixed basis approximation question */
    bool IsMixedApproxQ() {return fIsMixedApproxQ;}
    
    /** @brief Set Reduced basis approximation question */
    void SetRBApproxQ(bool state) { fIsRBApproxQ = state; }
    
    /** @brief Get Reduced basis approximation question */
    bool IsRBApproxQ() {return fIsRBApproxQ;}

    /** @brief Setup reporting times and time step size */
    void SetTimeControls(int n_times, REAL dt);
    
    /** @brief Setup reporting times and time step size */
    void SetNumericControls(int n_corrections, REAL epsilon_res, REAL epsilon_cor);
    
    /** @brief Store time values to be reported */
    TPZStack< REAL , 500 > ReportingTimes(){
        return fReportingTimes;
    }
    
    /** @brief Set Time step */
    void Setdt(REAL dt) { fdt = dt; }
    
    /** @brief Time step */
    REAL dt() { return fdt; }
    
    /** @brief Time */
    void SetTime(REAL time) { ftime = time; }
    
    /** @brief Time */
    REAL t() { return ftime; }
    
    /** @brief number of corrections steps */
    int n_steps() { return fn_steps; }
    
    /** @brief number of corrections steps */
    int n_corrections() { return fn_corrections; }
    
    /** @brief residue overal tolerance */
    REAL epsilon_res() { return fepsilon_res; }
    
    /** @brief correction overal tolerance */
    REAL epsilon_cor() { return fepsilon_cor; }
    
    
    void SetGravity(TPZVec<REAL> &g)
    {
        fg = g;
    }
    
    TPZVec<REAL> & Gravity()
    {
        return fg;
    }
    
    void SetBlocks(TPZStack< int > & blocks)
    {
        fblocks = blocks;
    }
    
    TPZStack< int > & Blocks()
    {
        return fblocks;
    }
    
};


#endif /* TPZSimulationData_h */