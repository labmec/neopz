#include "TPZElasticity3DOtiTopo.h"
//#include <TPZMaterialDataT.h>


//void TPZElasticity3DOtiTopo::FillDataRequirements(
//                                               TPZMaterialData &data) const
//{
//    data.fNeedsNormal = true;
//}
//
////! Contribution to the integration point
//void TPZElasticity3DOtiTopo::Contribute(const TPZMaterialDataT<CSTATE> &data,
//                                     REAL weight, TPZFMatrix<CSTATE> &ef)
//{
//    //index of integration point
//    const int gp_index = data.intGlobPtIndex;
//    const auto &mem_item = this->MemItem(gp_index);
//    
//    const auto &phi = data.phi;
//    const auto beta = this->Beta();
//    const int nshape=phi.Rows();
//    
//    const auto sol = mem_item.sol;
//    const auto &dsol = mem_item.dsol;
//    const auto &n = data.normal;
//    
//    TPZManVector<CSTATE,3> er,ur;
//    GetPermittivity(data.x,er);
//    GetPermeability(data.x,ur);
//    CSTATE cGx{0};
//    switch(fMode){
//        case ModeType::TE:
//            cGx = 1./ur[1];
//            break;
//        case ModeType::TM:
//            cGx = 1./er[1];
//            break;
//    }
//    
//    //dsol \cdot normal vec //FOR NOW LET US RESTRICT TO SOURCES ALIGNED WITH Y AXIS
//    // CSTATE dsol_n = 0;
//    // for(int ix = 0; ix < 3; ix++) { dsol_n += dsol[ix] * n[ix];}
//    // std::cout<<"pt"<<std::endl;
//    // std::cout<<"\tnormal"<<n<<std::endl;
//    // std::cout<<"\tdata.x"<<data.x<<std::endl;
//    // std::cout<<"\t sol.x"<<mem_item.x<<std::endl;
//    // std::cout<<"\t   sol"<<mem_item.sol<<std::endl;
//    // std::cout<<"\tdsol.n"<<dsol_n<<std::endl;
//    
//    
//    ///Source term
//    const auto src = 2. *( -1i * beta * sol + dsol[0]);
//    
//    //Contribution
//    for(int i = 0 ; i<nshape ; i++){
//        const CSTATE load = src * cGx * phi(i,0);
//        ef(i,0) += weight * load;
//    }
//}
//

TPZElasticity3DOtiTopo::TPZElasticity3DOtiTopo(int id, STATE E, STATE nu, TPZVec<STATE> &force) : TPZElasticity3D(id, E, nu, force) {
    
}


TPZElasticity3DOtiTopo * TPZElasticity3DOtiTopo::NewMaterial() const
{
    return new TPZElasticity3DOtiTopo(*this);
}
//! Unique identifier for serialization purposes
int TPZElasticity3DOtiTopo::ClassId() const
{
    return
    Hash("TPZElasticity3DOtiTopo")
    ^
    TPZElasticity3DOtiTopo::ClassId() << 1
    ^
    TPZMatWithMem<TPZOtiTopoDensity>::ClassId() << 2;
}
//! Write to stream(serialization method)
void TPZElasticity3DOtiTopo::Write(TPZStream &buf, int withclassid) const
{
    TPZElasticity3DOtiTopo::Write(buf,withclassid);
    TPZMatWithMem<TPZOtiTopoDensity>::Write(buf, withclassid);
}
//! Read from stream(serialization method)
void TPZElasticity3DOtiTopo::Read(TPZStream &buf, void *context)
{
    TPZElasticity3DOtiTopo::Read(buf,context);
    TPZMatWithMem<TPZOtiTopoDensity>::Read(buf, context);
}

int TPZElasticity3DOtiTopo::VariableIndex(const std::string &name) const {
    if(!strcmp("topodensity",name.c_str()))     return 100;
    return TPZElasticity3D::VariableIndex(name);
}

int TPZElasticity3DOtiTopo::NSolutionVariables(int var) const {
    switch(var) {
        case 100:
            return 1;
        default:
            return TPZElasticity3D::NSolutionVariables(var);
    }
}

void TPZElasticity3DOtiTopo::Contribute(const TPZMaterialDataT<STATE> &data,
                                        REAL weight,
                                        TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) {

    TPZFMatrix<STATE> ektemp(ek.Rows(),ek.Cols(),0.);
    const int gp_index = data.intGlobPtIndex;
    TPZOtiTopoDensity &densstruct = this->MemItem(gp_index);
    const STATE den = densstruct.fDen;

    TPZElasticity3D::Contribute(data, weight, ektemp, ef);
    // RECALL TO PUT ALSO IN EF
    
    ek += ektemp*den*den*den;
}



void TPZElasticity3DOtiTopo::Solution(const TPZMaterialDataT<STATE> &data,
                                      int var, TPZVec<STATE> &Solout) {

    
    if(var < 100){
        TPZElasticity3D::Solution(data, var, Solout);
        return;
    }
    
    if(var == 100){
        const int gp_index = data.intGlobPtIndex;
        TPZOtiTopoDensity &densstruct = this->MemItem(gp_index);
        Solout[0] = densstruct.fDen;
    }
    
    
    
}
