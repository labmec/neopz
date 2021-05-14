
#include "TPZTracerFlow.h"


TPZTracerFlow::TPZTracerFlow(){
    
}

/** @brief Constructor based on a material id */
TPZTracerFlow::TPZTracerFlow(int matid, int dimension) : TPZMaterial(matid) {
    
    m_mat_id = matid;
    m_dimension = dimension;
    m_mass_matrix_Q = false;
    m_dt = 0.0;
    m_phi = 0.0;
    m_fracture_epsilon = 0.0;
    
}

/** @brief Constructor based on a TRMMultiphase object */
TPZTracerFlow::TPZTracerFlow(const TPZTracerFlow &other) : TPZMaterial(other) {
    m_mat_id = other.m_mat_id;
    m_dimension = other.m_dimension;
    m_mass_matrix_Q = other.m_mass_matrix_Q;
    m_dt = other.m_dt;
    m_phi = other.m_phi;
    m_fracture_epsilon = other.m_fracture_epsilon;
}

TPZTracerFlow & TPZTracerFlow::operator=(const TPZTracerFlow &other){
    if (this != & other) // prevent self-assignment
    {
        TPZMaterial::operator=(other);
        m_mat_id = other.m_mat_id;
        m_dimension = other.m_dimension;
        m_mass_matrix_Q = other.m_mass_matrix_Q;
        m_dt = other.m_dt;
        m_phi = other.m_phi;
    }
    return *this;
}

/** @brief Default destructor */
TPZTracerFlow::~TPZTracerFlow(){
    
}

/** @brief Set the required data at each integration point */
void TPZTracerFlow::FillDataRequirements(TPZVec<TPZMaterialData> &datavec){
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
    }
}

/** @brief Set the required data at each integration point */
void TPZTracerFlow::FillBoundaryConditionDataRequirement(int type, TPZVec<TPZMaterialData> &datavec){
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        datavec[idata].fNeedsNormal = true;
    }
}

void TPZTracerFlow::FillDataRequirementsInterface(TPZMaterialData &data){
    data.SetAllRequirements(false);
    data.fNeedsSol = true;
    data.fNeedsNormal = true;
    if(fLinearContext == false){
        data.fNeedsNeighborSol = true;
    }
}

void TPZTracerFlow::FillDataRequirementsInterface(TPZMaterialData &data, std::map<int, TPZMaterialData> &datavec_left, std::map<int, TPZMaterialData> &datavec_right)
{
    
    data.SetAllRequirements(false);
    data.fNeedsSol = true;
    data.fNeedsNormal = true;
    for(auto &it : datavec_left)
    {
        it.second.SetAllRequirements(false);
        it.second.fNeedsSol = true;
        it.second.fNeedsNormal = true;
    }
    for(auto &it : datavec_left)
    {
        it.second.SetAllRequirements(false);
        it.second.fNeedsSol = true;
    }
    
}


/** @brief Print out the data associated with the material */
void TPZTracerFlow::Print(std::ostream &out){
    out << "\t Base class print:\n";
    out << " name of material : " << this->Name() << "\n";
    TPZMaterial::Print(out);
}

/** @brief Returns the variable index associated with the name */
int TPZTracerFlow::VariableIndex(const std::string &name){
    
    if (!strcmp("Sw", name.c_str())) return 0;
    if (!strcmp("So", name.c_str())) return 1;
    
    return TPZMaterial::VariableIndex(name);
}

/** @brief Returns the number of variables associated with varindex */
int TPZTracerFlow::NSolutionVariables(int var){
    switch(var) {
        case 0:
            return 1; // Scalar
        case 1:
            return 1; // Scalar
            
    }
    return TPZMaterial::NSolutionVariables(var);
}



// Contribute Methods being used
/** @brief Returns the solution associated with the var index */
void TPZTracerFlow::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout){
    
    int s_b    = 2;
    REAL sw = datavec[s_b].sol[0][0];
    
    Solout.Resize(this->NSolutionVariables(var));
    
    switch(var) {
        case 0:
        {
            Solout[0] = sw;
        }
            break;
        case 1:
        {
            Solout[0] = 1.0-sw;
        }
            break;
    }
    
    
}

REAL TPZTracerFlow::FractureFactor(TPZMaterialData & data){
    REAL alpha = 1.0;
    int dim = data.axes.Rows();
    switch (dim) {
        case 0:
            alpha = m_fracture_epsilon*m_fracture_epsilon*m_fracture_epsilon;
            break;
        case 1:
            alpha = m_fracture_epsilon*m_fracture_epsilon;
            break;
        case 2:
            alpha = m_fracture_epsilon;
            break;
        default:
            break;
    }
    return alpha;
}

void TPZTracerFlow::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
#ifdef PZDEBUG
    int nref =  datavec.size();
    DebugStop();
    if (nref != 3 ) {
        std::cout << " Erro. The size of the datavec is different from 3 \n";
        DebugStop();
    }
#endif
    
    
    int s_b = 2;
    
    REAL alpha = FractureFactor(datavec[s_b]);
    alpha =1.0;
    // Setting the phis
    TPZFMatrix<REAL>  &phiS =  datavec[s_b].phi;
    REAL s = datavec[s_b].sol[0][0];
    int n_phi_s = phiS.Rows();
    
    int firsts_s    = 0;
    
    for (int is = 0; is < n_phi_s; is++)
    {
        ef(is + firsts_s) += 1.0 * alpha * weight * m_phi * s * phiS(is,0); //es necesario¿?
        
        for (int js = 0; js < n_phi_s; js++)
        {
            ek(is + firsts_s, js + firsts_s) += alpha * weight * m_phi * (phiS(js,0) )* phiS(is,0);
        }
    }
    
}

void TPZTracerFlow::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    TPZFMatrix<STATE> ek_fake(ef.Rows(),ef.Rows());
    this->Contribute(datavec, weight, ek_fake, ef);
    
}


/**
 * Unique identifier for serialization purposes
 */
int TPZTracerFlow::ClassId() const{
    DebugStop();
    return  -1942;
}

/**
 * Save the element data to a stream
 */
void TPZTracerFlow::Write(TPZStream &buf, int withclassid){
    DebugStop();
}

/**
 * Read the element data from a stream
 */
void TPZTracerFlow::Read(TPZStream &buf, void *context){
    DebugStop();
}

void TPZTracerFlow::ContributeInterface(TPZMaterialData &data, std::map<int, TPZMaterialData> &datavecleft, std::map<int, TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){
    
    std::cout<<"el valor de m_mass_matrix_Q es: "<<m_mass_matrix_Q<<std::endl;
    if (m_mass_matrix_Q) {
        return;
    }
    
    int q_b = 0;
    int s_b = 2;
    
#ifdef PZDEBUG
    if(datavecleft.find(s_b) == datavecleft.end()) DebugStop();
    if(datavecright.find(q_b) == datavecright.end()) DebugStop();
#endif
    // Getting phis and solution for left material data
    TPZFMatrix<REAL>  &phiS_l =  datavecleft[s_b].phi;
    int n_phi_s_l = phiS_l.Rows();
    REAL s_l = datavecleft[s_b].sol[0][0];
    int firsts_s_l    = 0;
    
    // Getting phis and solution for right material data
    TPZFMatrix<REAL>  &phiS_r =  datavecright[s_b].phi;
    int n_phi_s_r = phiS_r.Rows();
    REAL s_r = datavecright[s_b].sol[0][0];
    int firsts_s_r    = n_phi_s_l;
    
    TPZManVector<REAL,3> n = data.normal;
    TPZManVector<STATE,3> q_l =  datavecleft[q_b].sol[0];
    REAL qn = 0.0;
    for (int i = 0; i < 3; i++) {
        qn += q_l[i]*n[i];
    }
    
    REAL beta = 0.0;
    // upwinding
    if (qn > 0.0) {
        beta = 1.0;
    }
    
    for (int is = 0; is < n_phi_s_l; is++) {
        
        ef(is + firsts_s_l) += +1.0 * m_dt * weight * (beta*s_l + (1.0-beta)*s_r)*phiS_l(is,0)*qn;
        
        for (int js = 0; js < n_phi_s_l; js++) {
            ek(is + firsts_s_l, js + firsts_s_l) += +1.0* m_dt * weight * beta * phiS_l(js,0) * phiS_l(is,0)*qn;
        }
        
        for (int js = 0; js < n_phi_s_r; js++) {
            ek(is + firsts_s_l, js + firsts_s_r) += +1.0* m_dt * weight * (1.0-beta) * phiS_r(js,0) * phiS_l(is,0)*qn;
        }
        
    }
    
    for (int is = 0; is < n_phi_s_r; is++) {
        
        ef(is + firsts_s_r) += -1.0* m_dt * weight * (beta*s_l + (1.0-beta)*s_r)*phiS_r(is,0)*qn;
        
        for (int js = 0; js < n_phi_s_l; js++) {
            ek(is + firsts_s_r, js + firsts_s_l) += -1.0* m_dt * weight * beta * phiS_l(js,0) * phiS_r(is,0)*qn;
        }
        
        for (int js = 0; js < n_phi_s_r; js++) {
            ek(is + firsts_s_r, js + firsts_s_r) += -1.0* m_dt * weight * (1.0-beta) * phiS_r(js,0) * phiS_r(is,0)*qn;
        }
        
    }
    
    
    
}

void TPZTracerFlow::ContributeInterface(TPZMaterialData &data, std::map<int, TPZMaterialData> &datavecleft, std::map<int, TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ef){
    TPZFMatrix<STATE> ek_fake(ef.Rows(),ef.Rows());
    this->ContributeInterface(data,datavecleft,datavecright, weight, ek_fake, ef);
    
}

void TPZTracerFlow::ContributeBCInterface(TPZMaterialData &data, std::map<int, TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    if (m_mass_matrix_Q) {
        return;
    }
    
    int q_b = 0;
    int s_b = 2;
    
#ifdef PZDEBUG
    if(datavecleft.find(s_b) == datavecleft.end()) DebugStop();
    if(datavecleft.find(q_b) == datavecleft.end()) DebugStop();
#endif
    // Getting phis and solution for left material data
    TPZFMatrix<REAL>  &phiS_l =  datavecleft[s_b].phi;
    int n_phi_s_l = phiS_l.Rows();
    REAL s_l = datavecleft[s_b].sol[0][0];
    int firsts_s_l    = 0;
    
    
    TPZManVector<REAL,3> n = data.normal;
    TPZManVector<STATE,3> q_l =  datavecleft[q_b].sol[0];
    REAL qn = 0.0;
    for (int i = 0; i < 3; i++) {
        qn += q_l[i]*n[i];
    }
    
    switch (bc.Type()) {
            
        case 0 :    // BC inlet
        {
            REAL s_inlet = bc.Val2()(0,0);
            if (qn > 0.0 || IsZero(qn)) {
                std::cout << "TPZTracerFlow:: Outlet flux in inlet boundary condition qn = " << qn << std::endl;
            }
            for (int is = 0; is < n_phi_s_l; is++) {
                ef(is + firsts_s_l) += +1.0* m_dt * weight * s_inlet * phiS_l(is,0)*qn;
            }
            
        }
            break;
            
        case 1 :    // BC outlet
        {
            
            if (qn > 0.0 || IsZero(qn)) {
                for (int is = 0; is < n_phi_s_l; is++) {
                    
                    ef(is + firsts_s_l) += +1.0* m_dt * weight * s_l*phiS_l(is,0)*qn;
                    
                    for (int js = 0; js < n_phi_s_l; js++) {
                        ek(is + firsts_s_l, js + firsts_s_l) += +1.0* m_dt * weight * phiS_l(js,0) * phiS_l(is,0)*qn;
                    }
                }
            }else{
                //                std::cout << "TPZTracerFlow:: Inlet flux on outlet boundary condition qn = " << qn << std::endl;
            }
            
            
        }
            break;
            
        default: std::cout << "This BC doesn't exist." << std::endl;
        {
            
            DebugStop();
        }
            break;
    }
    
    return;
    
}

void TPZTracerFlow::ContributeBCInterface(TPZMaterialData &data, std::map<int, TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    TPZFMatrix<STATE> ek_fake(ef.Rows(),ef.Rows());
    this->ContributeBCInterface(data,datavecleft, weight, ek_fake, ef, bc);
}
