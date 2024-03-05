#include "TPZElasticity2DOtiTopo.h"
//#include <TPZMaterialDataT.h>


//! Unique identifier for serialization purposes
int TPZOtiTopoDensity::ClassId() const
{
    return Hash("TPZOtiTopoDensity");
}

void TPZOtiTopoDensity::Write(TPZStream &buf, int withclassid) const
{
    buf.Write(&fDen);
}
//! Read from stream(serialization method)
void TPZOtiTopoDensity::Read(TPZStream &buf, void *context)
{
    buf.Read(&fDen);
}

void TPZOtiTopoDensity::Print(std::ostream &out) const
{
    out << "fDen "<<fDen<<'\n';
}

//void TPZElasticity2DOtiTopo::FillDataRequirements(
//                                               TPZMaterialData &data) const
//{
//    data.fNeedsNormal = true;
//}
//
////! Contribution to the integration point
//void TPZElasticity2DOtiTopo::Contribute(const TPZMaterialDataT<CSTATE> &data,
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

TPZElasticity2DOtiTopo::TPZElasticity2DOtiTopo(int id, STATE E, STATE nu, STATE fx, STATE fy, int planestress) : TPZElasticity2D(id, E, nu, fx, fy,planestress) {
    
}


TPZElasticity2DOtiTopo * TPZElasticity2DOtiTopo::NewMaterial() const
{
    return new TPZElasticity2DOtiTopo(*this);
}
//! Unique identifier for serialization purposes
int TPZElasticity2DOtiTopo::ClassId() const
{
    return
    Hash("TPZElasticity2DOtiTopo")
    ^
    TPZElasticity2DOtiTopo::ClassId() << 1
    ^
    TPZMatWithMem<TPZOtiTopoDensity>::ClassId() << 2;
}
//! Write to stream(serialization method)
void TPZElasticity2DOtiTopo::Write(TPZStream &buf, int withclassid) const
{
    TPZElasticity2DOtiTopo::Write(buf,withclassid);
    TPZMatWithMem<TPZOtiTopoDensity>::Write(buf, withclassid);
}
//! Read from stream(serialization method)
void TPZElasticity2DOtiTopo::Read(TPZStream &buf, void *context)
{
    TPZElasticity2DOtiTopo::Read(buf,context);
    TPZMatWithMem<TPZOtiTopoDensity>::Read(buf, context);
}

int TPZElasticity2DOtiTopo::VariableIndex(const std::string &name) const {
    if(!strcmp("topodensity",name.c_str()))     return 100;
    return TPZElasticity2D::VariableIndex(name);
}

int TPZElasticity2DOtiTopo::NSolutionVariables(int var) const {
    switch(var) {
        case 100:
            return 1;
        default:
            return TPZElasticity2D::NSolutionVariables(var);
    }
}

void TPZElasticity2DOtiTopo::Contribute(const TPZMaterialDataT<STATE> &data,
                                        REAL weight,
                                        TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) {

    TPZFMatrix<STATE> ektemp(ek.Rows(),ek.Cols(),0.);
    const int gp_index = data.intGlobPtIndex;
    TPZOtiTopoDensity &densstruct = this->MemItem(gp_index);
    const STATE den = densstruct.fDen;

    TPZElasticity2D::Contribute(data, weight, ektemp, ef);
    // RECALL TO PUT ALSO IN EF
    
    ek += ektemp*den*den*den;
}



void TPZElasticity2DOtiTopo::Solution(const TPZMaterialDataT<STATE> &data,
                                      int var, TPZVec<STATE> &Solout) {

    
    if(var < 100){
        TPZElasticity2D::Solution(data, var, Solout);
        return;
    }
    
    if(var == 100){
        const int gp_index = data.intGlobPtIndex;
        TPZOtiTopoDensity &densstruct = this->MemItem(gp_index);
        Solout[0] = densstruct.fDen;
    }
    
    
    
}
