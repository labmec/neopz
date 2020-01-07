//
//  TRMSimulationData.h
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#ifndef __PZ__TRMSimulationData__
#define __PZ__TRMSimulationData__

#include <stdio.h>
#include "pzreal.h"
#include "TRMPhaseProperties.h"
#include "TRMRockProperties.h"
#include "TRMPetrophysicsProperties.h"
#include "TRMSpatialPropertiesMap.h"
#include "TRMRawData.h"


class TRMSimulationData {
    
protected:

    /** @brief initial state */
    bool fIsInitialStateQ;
    
    /** @brief current time state */
    bool fIsCurrentStateQ;
    
    /** @brief one-phase flow */
    bool fIsOnePhaseQ;
    
    /** @brief two-phase flow */
    bool fIsTwoPhaseQ;
    
    /** @brief three-phase flow */
    bool fIsThreePhaseQ;
    
    /** @brief Definition of the flow system one - two or three phase */
    TPZStack<std::string> fSystemType;
    
    /** @brief Definition gravity field */
    TPZVec<STATE> fg;
    
    /** @brief Material identifier for interfaces */
    int fInterface_mat_Id;
    
    /** @brief Store time values to be reported */
    TPZStack< STATE , 500 > fReportingTimes;
    
    /** @brief Store flags values of mixed problem to be reported */
    TPZStack< bool , 500 > fReportingTimesMixedQ;
    
    /** @brief ntime steps */
    int fn_steps;
    
    /** @brief Initial time */
    STATE ftime_0;
    
    /** @brief Final time */
    STATE ftime_n;
    
    /** @brief Time step */
    STATE fdt;
    
    /** @brief Time step */
    STATE ftime;
    
    /** @brief Min time step */
    STATE fdt_min;
    
    /** @brief Max time step */
    STATE fdt_max;
    
    /** @brief Increment dt factor */
    STATE fdt_up;
    
    /** @brief Decrement dt factor */
    STATE fdt_down;
    
    /** @brief number of corrections steps */
    int fn_corrections;
    
    /** @brief number of Reduced basis */
    int fm_rb_functions;
    
    /** @brief residue overal tolerance */
    STATE fepsilon_res;
    
    /** @brief correction overal tolerance */
    STATE fepsilon_cor;

    /** @brief Use of geomechanic */
    bool fIsGeomechanicQ;
    
    /** @brief Use of quasi newton method */
    bool fIsQuasiNewtonQ;
    
    /** @brief set the use of pardiso for elliptic and parabolic operators Ae and Ap */
    bool fUsePardisoQ;
    
    /** @brief set the use p adaptation on wellbores */
    bool fIsAdataptedQ;
    
    /** @brief set the use enhanced pressure accuracy */
    bool fEnhancedPressureQ;

    
    /** @brief Use, level and resolution of MHM process */
    std::pair<bool, std::pair<int, int> > fMHMResolutionQ;
    
    /** @brief Use of increased transpor resolution transfers operators */
    std::pair<bool, int> fIncreaseTransporResolutionQ;
    
    /** @brief Use of RB method that surrogates */
    std::pair<bool, std::pair<bool, TPZStack<int> > > fReduceBasisQ;
    
    /** @brief Gmsh grid file */
    std::string fGridName;
    
    /** @brief Set SPE10 fields file */
    std::pair< std::string , std::string > fPermPorFields;
    
    /** @brief number of blocks i, j and k  */
    TPZStack<int> fNBlocks;
    
    /** @brief size of blocks dx, dy and dz  */
    TPZStack<REAL> fBlocks_sizes;
    
    /** @brief Autopointer of the RawData */
    TPZAutoPointer<TRMRawData> fRawData;
    
    /** @brief Autopointer of all the petrophysics data */
    TPZAutoPointer<TRMPetrophysicsProperties> fPetroPhysics;
    
    /** @brief phase alpha */
    TPZAutoPointer<TRMPhaseProperties> fPhase_alpha;
    
    /** @brief phase beta */
    TPZAutoPointer<TRMPhaseProperties> fPhase_beta;
    
    /** @brief phase gamma */
    TPZAutoPointer<TRMPhaseProperties> fPhase_gamma;
    
    /** @brief Stores all the rock and geomechanic properties */
    TRMRockProperties fRockProp;
    
    /** @brief Stores the spatial information given in maps */
    TPZAutoPointer<TRMSpatialPropertiesMap> fMap;
    
    /** @brief L2 projection material id for gradient reconstruction */
    int fl2_projection_material_id;
    
    /** @brief Skeleton dfault material id for MHM substructuring */
    int fSkeleton_material_id;
    
    /** @brief Define the use of linear gradient reconstruction */
    bool fUseGradientRQ;
    
    
public:
    
    
    /** @brief default constructor */
    TRMSimulationData();
    
    /** @brief default constructor */
    TRMSimulationData(const TRMSimulationData &copy)
    {
        DebugStop();
    }
    
    /** @brief default constructor */
    TRMSimulationData &operator=(const TRMSimulationData &copy)
    {
        DebugStop();
        return *this;
    }
    
    /** @brief destructor */
    ~TRMSimulationData();
    

    /** @brief Set initial state */
    void SetInitialStateQ(bool state) { fIsInitialStateQ = state; }
    
    /** @brief Get initial state */
    bool IsInitialStateQ() {return fIsInitialStateQ;}
    
    /** @brief current time state */
    void SetCurrentStateQ(bool state) { fIsCurrentStateQ = state; }
    
    /** @brief current time state */
    bool IsCurrentStateQ() {return fIsCurrentStateQ;}
    
    /** @brief Mono-phasic system */
    bool IsOnePhaseQ() {return fIsOnePhaseQ;}
    
    /** @brief Two-phasic system */
    bool IsTwoPhaseQ() {return fIsTwoPhaseQ;}
    
    /** @brief Three-phasic system */
    bool IsThreePhaseQ() {return fIsThreePhaseQ;}
    
    /** @brief Set the type of flow system one, two or three phase */
    void SetSystemType(TPZStack<std::string> &SystemType, TPZStack< TPZAutoPointer<TRMPhaseProperties> > &Phases){

        if(Phases.size() != SystemType.size()){
            DebugStop();
        }
        
        TPZAutoPointer<TRMPetrophysicsProperties> petro_physics = new TRMPetrophysicsProperties;
        petro_physics->SetSystemType(SystemType);
        switch (SystemType.size()) {
            case 1:
            {
                fIsOnePhaseQ    = true;
                SetPhaseAlpha(Phases[0]);
                petro_physics->SetPhaseAlpha(Phases[0]);
            }
                break;
            case 2:
            {
                fIsTwoPhaseQ = true;
                SetPhaseAlpha(Phases[0]);
                SetPhaseBeta(Phases[1]);
        
                petro_physics->SetPhaseAlpha(Phases[0]);
                petro_physics->SetPhaseBeta(Phases[1]);
            }
                break;
            case 3:
            {
                fIsThreePhaseQ = true;
                SetPhaseAlpha(Phases[0]);
                SetPhaseBeta(Phases[1]);
                SetPhaseGamma(Phases[2]);
                
                petro_physics->SetPhaseAlpha(Phases[0]);
                petro_physics->SetPhaseBeta(Phases[1]);
                petro_physics->SetPhaseGamma(Phases[2]);
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
        fPetroPhysics  = petro_physics;
        
    }
    
    /** @brief Get the type of flow system one, two or three phase */
    TPZStack<std::string>  SystemType(){
        return fSystemType;
    }
    
    /** @brief Set phase alpha */
    void SetPhaseAlpha(TPZAutoPointer<TRMPhaseProperties> &alpha);
    
    /** @brief Get phase alpha */
    TPZAutoPointer<TRMPhaseProperties> & AlphaProp();
    
    /** @brief Set phase beta */
    void SetPhaseBeta(TPZAutoPointer<TRMPhaseProperties> &beta);
    
    /** @brief Get phase beta */
    TPZAutoPointer<TRMPhaseProperties> & BetaProp();
    
    /** @brief Set phase gamma */
    void SetPhaseGamma(TPZAutoPointer<TRMPhaseProperties> &gamma);
    
    /** @brief Get phase gamma */
    TPZAutoPointer<TRMPhaseProperties> & GammaProp();
    
    /** @brief Setup reporting times and time step size */
    void SetTimeControls(int n_times, STATE dt, STATE dt_in, STATE dt_de, STATE dt_max, STATE dt_min, TPZStack< std::pair< STATE , bool> , 500 > ReportingTimes);
    
    /** @brief Setup reporting times and time step size */
    void SetNumericControls(int n_corrections, STATE epsilon_res, STATE epsilon_cor, bool IsQuasiNewtonQ, bool IsAdataptedQ, bool EnhancedPressureQ,  bool UsePardisoQ);

    
    /** @brief Setup Use, level and resolution of MHM process */
    void SetMHMResolution(std::pair<bool, std::pair<int, int> > MHMResolutionQ){
        fMHMResolutionQ = MHMResolutionQ;
    }
    
    /** @brief Get Use, level and resolution of MHM process */
    std::pair<bool, std::pair<int, int> > & MHMResolution(){
        return fMHMResolutionQ;
    }
    
    /** @brief Setup transpor resolution options */
    void SetTransporResolution(std::pair<bool, int> IncreaseTransporResolutionQ){
        fIncreaseTransporResolutionQ = IncreaseTransporResolutionQ;
    }
    
    /** @brief Get transpor resolution options */
    std::pair<bool, int> & TransporResolution(){
        return fIncreaseTransporResolutionQ;
    }
    
    /** @brief Setup RB method options */
    void SetReducedBasisResolution(std::pair<bool, std::pair<bool, TPZStack<int> > > ReduceBasisQ){
        fReduceBasisQ = ReduceBasisQ;
    }
    
    /** @brief Get RB method options */
    std::pair<bool, std::pair<bool, TPZStack<int> > > & ReducedBasisResolution(){
        return fReduceBasisQ;
    }
    
    /** @brief Set Gmsh grid file */
    void SetGridName(std::string GridName){
        fGridName = GridName;
    }
    
    /** @brief Get Gmsh grid file */
    std::string & GridName(){
        return fGridName;
    }
    
    /** @brief Set SPE10 fields file */
    void SetSpatialFields(TPZStack<int> NBlocks, TPZStack<REAL> Blocks_sizes, std::pair< std::string , std::string > PermPorFields){
        fNBlocks = NBlocks;
        fBlocks_sizes = Blocks_sizes;
        fPermPorFields = PermPorFields;
    }
    
    /** @brief Set SPE10 fields file */
    std::pair< std::string , std::string > & SpatialFields(){
        return fPermPorFields;
    }
    
    /** @brief Get SPE10 fields file */
    TPZStack<int> & NBlocks(){
        return fNBlocks;
    }
    
    /** @brief Get SPE10 fields file */
    TPZStack<REAL> & Blocks_sizes(){
        return fBlocks_sizes;
    }
    
    /** @brief Store time values to be reported */
    TPZStack< STATE , 500 > & ReportingTimes(){
        return fReportingTimes;
    }
    
    /** @brief Store time values to be reported */
    TPZStack< bool, 500 > & ReportingTimesMixedQ(){
        return fReportingTimesMixedQ;
    }
    
    /** @brief Initial time */
    STATE time_0() { return ftime_0; }
    
    /** @brief Final time */
    STATE time_n() { return ftime_n; }
    
    /** @brief Set Time step */
    void Setdt(STATE dt) { fdt = dt; }
    
    /** @brief Time step */
    STATE dt() { return fdt; }

    /** @brief Time */
    void SetTime(STATE time) { ftime = time; }
    
    /** @brief Time */
    STATE t() { return ftime; }
    
    /** @brief Min time step */
    STATE dt_min() { return fdt_min; }
    
    /** @brief Max time step */
    STATE dt_max() { return fdt_max; }
    
    /** @brief Increment dt factor */
    STATE dt_up() { return fdt_up; }
    
    /** @brief Decrement dt factor */
    STATE dt_down() { return fdt_down; }
    
    /** @brief number of corrections steps */
    int n_steps() { return fn_steps; }
    
    /** @brief number of corrections steps */
    int n_corrections() { return fn_corrections; }
    
    /** @brief set number of corrections steps */
    void Set_m_RB_functions(int m) { fm_rb_functions = m; }
    
    /** @brief get number of corrections steps */
    int m_RB_functions() { return fm_rb_functions; }
    
    /** @brief residue overal tolerance */
    STATE epsilon_res() { return fepsilon_res; }
    
    /** @brief correction overal tolerance */
    STATE epsilon_cor() { return fepsilon_cor; }
    
    /** @brief Set directive for the use of quasi newton method */
    void SetGeomechanicQ(bool IsGeomechanicQ) { fIsGeomechanicQ = IsGeomechanicQ;}
    
    /** @brief Get directive for the use geomechanic coupling */
    bool IsGeomechanicQ() {return fIsGeomechanicQ;}
    
    /** @brief Set the use of pardiso for elliptic and parabolic operators Ae and Ap */
    void SetUsePardisoQ(bool UsePardisoQ) { fUsePardisoQ = UsePardisoQ;}
    
    /** @brief Get the use of pardiso for elliptic and parabolic operators Ae and Ap */
    bool UsePardisoQ() {return fUsePardisoQ;}
    
    /** @brief Get directive for the use of quasi newton method */
    bool IsQuasiNewtonQ() {return fIsQuasiNewtonQ;}
    
    /** @brief Get the use p adaptation on wellbores */
    bool IsAdataptedQ() {return fIsAdataptedQ;}
    
    /** @brief Get the use enhanced pressure accuracy */
    bool IsEnhancedPressureQ() {return fEnhancedPressureQ;}
    
    /** @brief Material identifier for interfaces */
    int InterfacesMatId() { return fInterface_mat_Id; }
    
    /** @brief Set the directive for use of gradient reconstruction */
    void SetUseGradientR(bool UseGR) { fUseGradientRQ = UseGR; }
    
    /** @brief Get the directive for use of gradient reconstruction */
    int UseGradientR() { return fUseGradientRQ; }
    
    /** @brief Set autopointer of the RawData */
    void SetRawData(TPZAutoPointer<TRMRawData> &RawData);
    
    /** @brief Get autopointer of the RawData */
    TPZAutoPointer<TRMRawData>  RawData(){
        return fRawData;
    }
    
    /** @brief Stores all the petrophysics data */
    void SetPetroPhysics(TPZAutoPointer<TRMPetrophysicsProperties> &PetroPhysics)
    {
        fPetroPhysics = PetroPhysics;
    }
    
    TPZAutoPointer<TRMPetrophysicsProperties> & PetroPhysics()
    {
        return fPetroPhysics;
    }    
  
    /** @brief Stores all the rock and geomechanic properties */
    void SetRockProp(TRMRockProperties RockProp)
    {
        fRockProp = RockProp;
    }
    
    TRMRockProperties & RockProp()
    {
        return fRockProp;
    }
    
    /** @brief Stores the spatial information given in maps */
    void SetMap(TPZAutoPointer<TRMSpatialPropertiesMap> &SpatialProp)
    {
        fMap = SpatialProp;
    }
    
    TPZAutoPointer<TRMSpatialPropertiesMap> & Map()
    {
        return fMap;
    }
    
    void SetGravity(TPZVec<STATE> &g)
    {
        fg = g;
    }
    
    TPZVec<STATE> & Gravity()
    {
        return fg;
    }
    
    /** @brief L2 projection material id for gradient reconstruction */
    void SetL2_Projection_material_Id(int l2_projection){
        fl2_projection_material_id = l2_projection;
    }
    
    /** @brief L2 projection material id for gradient reconstruction */
    int L2_Projection_material_Id(){
        return fl2_projection_material_id;
    }
    
    /** @brief Set the Skeleton material id for MHM substructuring */
    void SetSkeleton_material_Id(int skeleton_id){
        fSkeleton_material_id = skeleton_id;
    }
    
    /** @brief Skeleton dfault material id for MHM substructuring */
    int Skeleton_material_Id(){
        return fSkeleton_material_id;
    }
    

};

#endif /* defined(__PZ__TRMSimulationData__) */