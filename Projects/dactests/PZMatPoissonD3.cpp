//
//  PZMatPoissonD3.cpp
//  PZ
//
//  Created by Douglas Castro on 5/23/14.
//
//

#include "PZMatPoissonD3.h"
#include "pzlog.h"
#include "pzfmatrix.h"
#include "pzmaterialdata.h"
#include "pzbndcond.h"
#include "pzaxestools.h"
#include "pzlog.h"

#include <cmath>


#include <iostream>

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.tpzmatpoissonD3"));
#endif

using namespace std;

TPZMatPoissonD3::TPZMatPoissonD3():TPZDiscontinuousGalerkin(){
	
    /** Valor da funcao de carga */
    fF = 0.; //fF
    
    /** Dimensao do dominio */
    fDim = 2;
    
    /** Material id not initialized */
    fMatId = -1;
    
    /** Coeficiente que multiplica o gradiente */
    fK = 1.;
    
    fvisc = 1.;
    fInvK.Resize(1, 1);
    fTensorK.Resize(1, 1);
    fTensorK.Identity();
    fInvK.Identity();
    fPermeabilityFunction = NULL;

}

TPZMatPoissonD3::TPZMatPoissonD3(int matid, int dim):TPZDiscontinuousGalerkin(matid){
	
    if(dim<0 || dim >3){
        DebugStop();
    }
    
    /** Valor da funcao de carga */
    fF = 0.; //fF
    
    /** Dimensao do dominio */
    fDim = dim;
    
    /** Material id no initialized */
    fMatId = matid;
    
    /** Coeficiente que multiplica o gradiente */
    fK = 1.;
    
    fvisc = 1.;
    
    fInvK.Redim(dim, dim);
    fTensorK.Resize(dim, dim);
    fInvK.Identity();
    fTensorK.Identity();
    fPermeabilityFunction = NULL;
}

TPZMatPoissonD3::~TPZMatPoissonD3(){
}

TPZMatPoissonD3::TPZMatPoissonD3(const TPZMatPoissonD3 &copy):TPZDiscontinuousGalerkin(copy){
    
    this->operator=(copy);
}

TPZMatPoissonD3 & TPZMatPoissonD3::operator=(const TPZMatPoissonD3 &copy){
    
    TPZDiscontinuousGalerkin::operator = (copy);
    this->fF = copy.fF; //fF
    this->fDim = copy.fDim;
    this->fMatId = copy.fMatId;
    this->fK = copy.fK;
    
	return *this;
}


void TPZMatPoissonD3::Print(std::ostream &out) {
	out << "name of material : " << Name() << "\n";
    out << "Dimesion of problem " << fDim << endl;
	out << "Material ID  "<< fMatId << endl;
    out << "Forcing function  "<< fF << endl;
    out << "Grad Coeficient  "<< fK << endl;
	out << "Base Class properties :";
	TPZMaterial::Print(out);
	out << "\n";
}


// Contribute methods
// esse metodo esta ok
void TPZMatPoissonD3::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    
#ifdef PZDEBUG
	int nref =  datavec.size();
	if (nref != 2 ) {
        std::cout << " Erro. The size of the datavec is different from 2 \n";
		DebugStop();
	}
#endif
    REAL force = fF;
    if(fForcingFunction) {
		TPZManVector<STATE> res(1);
		fForcingFunction->Execute(datavec[1].x,res);
		force = res[0];
	}
    
    // Setting the phis
    TPZFMatrix<REAL>  &phiQ =  datavec[0].phi;
    TPZFMatrix<REAL>  &phip =  datavec[1].phi;
	//TPZFMatrix<REAL> &dphiQ = datavec[0].dphix;
    TPZFMatrix<REAL> &dphipLoc = datavec[1].dphix;
    
    TPZFNMatrix<200,REAL> dphip(3,datavec[1].dphix.Cols(),0.0);
    
    for (int ip = 0; ip<dphip.Cols(); ip++) {
        for (int d = 0; d<dphipLoc.Rows(); d++) {
            for (int j=0; j< 3; j++) {
                dphip(j,ip)+=datavec[1].axes(d,j)*dphipLoc(d,ip);
            }
        }
    }
    
    int phrq, phrp;
    phrp = phip.Rows();
    phrq = datavec[0].fVecShapeIndex.NElements();
	
	//Calculate the matrix contribution for flux. Matrix A
    // A matriz de rigidez é tal que A{ij}=\int_\Omega K^{-1} \varphi_j\cdot\varphi_i d\Omega
    // K, futuramente sera uma matriz ou funcao, deve-se ter cuidado com essa parte da inversao de K
    
    TPZFNMatrix<3,REAL> PermTensor = fTensorK;
    TPZFNMatrix<3,REAL> InvPermTensor = fInvK;
    
    if(fPermeabilityFunction){
        PermTensor.Redim(fDim,fDim);
        InvPermTensor.Redim(fDim,fDim);
        TPZFNMatrix<9,STATE> resultMat;
        TPZManVector<STATE,3> res;
        fPermeabilityFunction->Execute(datavec[1].x,res,resultMat);
        
        for(int id=0; id<fDim; id++){
            for(int jd=0; jd<fDim; jd++){
                
                PermTensor(id,jd) = resultMat(id,jd);
                InvPermTensor(id,jd) = resultMat(id+fDim,jd);
            }
        }
    }

    
    //REAL InvK = 1./fK;
    for (int iq = 0; iq<phrq; iq++)
    {
        int ivecind = datavec[0].fVecShapeIndex[iq].first;
        int ishapeind = datavec[0].fVecShapeIndex[iq].second;
        TPZFNMatrix<3> ivec(3,1);
        ivec(0,0) = datavec[0].fDeformedDirections(0,ivecind);
        ivec(1,0) = datavec[0].fDeformedDirections(1,ivecind);
        ivec(2,0) = datavec[0].fDeformedDirections(2,ivecind);
        
        TPZFNMatrix<3,REAL> jvecZ(fDim,1,0.);
        
        for (int jq=0; jq<phrq; jq++)
        {
            TPZFNMatrix<3> jvec(3,1);
            int jvecind = datavec[0].fVecShapeIndex[jq].first;
            int jshapeind = datavec[0].fVecShapeIndex[jq].second;
            jvec(0,0) = datavec[0].fDeformedDirections(0,jvecind);
            jvec(1,0) = datavec[0].fDeformedDirections(1,jvecind);
            jvec(2,0) = datavec[0].fDeformedDirections(2,jvecind);
            
            
            
            //dot product between Kinv[u]v
            jvecZ.Zero();
            for(int id=0; id<fDim; id++){
                for(int jd=0; jd<fDim; jd++){
                    jvecZ(id,0) += InvPermTensor(id,jd)*jvec(jd,0);
                }
            }
            REAL prod = 0.;
            for(int id=0; id < fDim;id++) prod += ivec(id,0)*jvecZ(id,0);
            ek(iq,jq) += fvisc*weight*phiQ(ishapeind,0)*phiQ(jshapeind,0)*prod;

            
            //dot product between u and v
//            REAL prod = 0.;
//            for(int iloc=0; iloc<3; iloc++)
//            {
//                prod += ivec(iloc,0)*jvec(iloc,0);
//                //ivec(0,0)*jvec(0,0) + ivec(1,0)*jvec(1,0) + ivec(2,0)*jvec(2,0);
//            }
//            ek(iq,jq) += InvK*weight*phiQ(ishapeind,0)*phiQ(jshapeind,0)*prod;
            
        }
    }
    
    // Coupling terms between flux and pressure. Matrix B
    // A matriz de rigidez é tal que B{ij}=\int_\Omega \nabla\phi_j\cdot\varphi_i d\Omega
    
    for(int iq=0; iq<phrq; iq++)
    {
        int ivecind = datavec[0].fVecShapeIndex[iq].first;
        int ishapeind = datavec[0].fVecShapeIndex[iq].second;
        TPZFNMatrix<3> ivec(3,1);
        ivec(0,0) = datavec[0].fDeformedDirections(0,ivecind);
        ivec(1,0) = datavec[0].fDeformedDirections(1,ivecind);
        ivec(2,0) = datavec[0].fDeformedDirections(2,ivecind);
        
        for (int jp=0; jp<phrp; jp++)
        {
            //dot product between  varphi and grad phi
            
            REAL prod = 0.;
            for(int iloc=0; iloc<3; iloc++)
            {
                prod += (  ivec(iloc,0)*phiQ(ishapeind,0)  )*dphip(iloc,jp);
            }
            
            REAL fact = weight*prod;
            // Matrix B
            ek(iq, phrq+jp) += fact;
            
            // Matrix B^T
            ek(phrq+jp,iq) += fact;
            
        }
    }
    
    //termo fonte referente a equacao da pressao
    for(int ip=0; ip<phrp; ip++){
        ef(phrq+ip,0) += (-1.)*weight*force*phip(ip,0);
    }
    
}

void TPZMatPoissonD3::ContributeBC(TPZVec<TPZMaterialData> &datavec,REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
    return;
}


//Contribute interface methods
void TPZMatPoissonD3::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, TPZVec<TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    TPZFMatrix<REAL> &phiQL = dataleft[0].phi;
	//TPZFMatrix<REAL> &phiQR = dataright[0].phi;
	TPZFMatrix<REAL> &phiPL = dataleft[1].phi;
	TPZFMatrix<REAL> &phiPR = dataright[1].phi;
    
    int pRowsL = phiPL.Rows();
    int pRowsR = phiPR.Rows();
    int QRowsL = dataleft[0].fVecShapeIndex.NElements();
    //int QRowsR = dataright[0].fVecShapeIndex.NElements();
    
    TPZManVector<REAL,3> &normal = data.normal;
    
    // Para o bloco left-left
    for (int iq = 0; iq < QRowsL; iq++)
    {
        int ivecind = dataleft[0].fVecShapeIndex[iq].first;
        int ishapeind = dataleft[0].fVecShapeIndex[iq].second;
        TPZFNMatrix<3> ivec(3,1);
        ivec(0,0) = dataleft[0].fDeformedDirections(0,ivecind);
        ivec(1,0) = dataleft[0].fDeformedDirections(1,ivecind);
        ivec(2,0) = dataleft[0].fDeformedDirections(2,ivecind);
        ivec *= phiQL(ishapeind,0);
        
        
        REAL NormalProjectioni = 0.;
        for(int iloc=0; iloc<3; iloc++)
        {
            NormalProjectioni += ivec(iloc,0)*normal[iloc];
        }
        
        for (int jp = 0; jp < pRowsL; jp++)
        {
            
            REAL fact = weight*NormalProjectioni*phiPL(jp,0);
            
            ek(iq, QRowsL+jp) += -1.0*fact; // ???????
            
            ek(QRowsL+jp, iq) += -1.0*fact; // ???????
            
        }
    }
    // Para o bloco left-right
    for (int iq = 0; iq < QRowsL; iq++)
    {
        int ivecind = dataleft[0].fVecShapeIndex[iq].first;
        int ishapeind = dataleft[0].fVecShapeIndex[iq].second;
        TPZFNMatrix<3> ivec(3,1);
        ivec(0,0) = dataleft[0].fDeformedDirections(0,ivecind);
        ivec(1,0) = dataleft[0].fDeformedDirections(1,ivecind);
        ivec(2,0) = dataleft[0].fDeformedDirections(2,ivecind);
        ivec *= phiQL(ishapeind,0);
        
        REAL NormalProjectionj = 0.;
        for(int iloc=0; iloc<fDim; iloc++)
        {
            NormalProjectionj += ivec(iloc,0)*normal[iloc];
        }
     
        
        for (int jp = 0; jp < pRowsR; jp++)
        {
        
            REAL fact = weight*NormalProjectionj*phiPR(jp,0);
            
            int pulo = QRowsL + pRowsL + QRowsL;
            
            ek(iq, pulo+jp) += 1.0*fact;
            
            ek(pulo+jp, iq) += 1.0*fact;
            
        }
    }
    
    
    //#ifdef LOG4CXX
    //	if(logdata->isDebugEnabled())
    //	{
    //		std::stringstream sout;
    //		ek.Print("ek = ",sout,EMathematicaInput);
    //		ef.Print("ef = ",sout,EMathematicaInput);
    //		LOGPZ_DEBUG(logdata,sout.str())
    //	}
    //#endif
    
}

void TPZMatPoissonD3::ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
    
    
#ifdef PZDEBUG
	int nref =  dataleft.size();
	if (nref != 2 ) {
        std::cout << " Error. This implementation needs only two computational meshes. \n";
		DebugStop();
	}
#endif
    
#ifdef PZDEBUG
	int bref =  bc.Val2().Rows();
	if (bref != 2 ) {
        std::cout << " Erro. The size of the datavec is different from 2 \n";
		DebugStop();
	}
#endif
    
    REAL Qn = bc.Val2()(0,0); // cuidado para, na hora de passar os valores de cond contorno, seguir essa ordem
    REAL Pd = 0.0; // = bc.Val2()(1,0); // fluxo normal na primeira casa e pressao na segunda
    
	TPZManVector<REAL,3> &normal = data.normal;
	//REAL n1 = normal[0];
	//REAL n2 = normal[1];
    
    //REAL v2;
    if(bc.HasForcingFunction())
    {
		TPZManVector<STATE> res(3);
		bc.ForcingFunction()->Execute(dataleft[0].x,res);
		Pd = res[0];
        Qn = res[0];
	}else
    {
        Pd = bc.Val2()(1,0);
    }
    

    // Setting the phis
    TPZFMatrix<REAL>  &phiQ =  dataleft[0].phi;
    TPZFMatrix<REAL>  &phip =  dataleft[1].phi;
	//TPZFMatrix<REAL> &dphiQ = datavec[0].dphix;
    //TPZFMatrix<REAL> &dphip = datavec[1].dphix;

    int phrq, phrp;
    phrp = phip.Rows();
    phrq = dataleft[0].fVecShapeIndex.NElements();

	//Calculate the matrix contribution for boundary conditions
    for (int iq = 0; iq<phrq; iq++)
    {
        int ivecind = dataleft[0].fVecShapeIndex[iq].first;
        int ishapeind = dataleft[0].fVecShapeIndex[iq].second;
        TPZFNMatrix<3> ivec(3,1);
        ivec(0,0) = dataleft[0].fDeformedDirections(0,ivecind);
        ivec(1,0) = dataleft[0].fDeformedDirections(1,ivecind);
        ivec(2,0) = dataleft[0].fDeformedDirections(2,ivecind);
        ivec *= phiQ(ishapeind,0);
        
        
        REAL NormalProjectioni = 0.;
        for(int iloc=0; iloc<fDim; iloc++)
        {
            NormalProjectioni += ivec(iloc,0)*normal[iloc];
        }
        
        for (int jp=0; jp<phrp; jp++)
        {
            
            REAL integration = weight*NormalProjectioni*phip(jp,0);
            
            //para a equacao do fluxo - 1o conjunto da formulacao
            ek(iq, phrq+jp) += (-1.0)*integration;
            
            // para a equacao da pressao - 2o conjunto da formulacao
            ek(phrq+jp, iq) += (-1.0)*integration;
            
        }
    }
    

    //if (bc.Type()==0){std::cout << "...." << std::endl;}

    switch (bc.Type())
    {  
        case 0:  // Dirichlet
        {
            //REAL InvK = 1./fK;
            
                        //termo fonte referente a equacao do fluxo
            for (int iq = 0; iq<phrq; iq++)
            {
                int ivecind = dataleft[0].fVecShapeIndex[iq].first;
                int ishapeind = dataleft[0].fVecShapeIndex[iq].second;
                TPZFNMatrix<3> ivec(3,1);
                ivec(0,0) = dataleft[0].fDeformedDirections(0,ivecind);
                ivec(1,0) = dataleft[0].fDeformedDirections(1,ivecind);
                ivec(2,0) = dataleft[0].fDeformedDirections(2,ivecind);
                ivec *= phiQ(ishapeind,0);
                
                
                REAL NormalProjectioni = 0.;
                for(int iloc=0; iloc<fDim; iloc++)
                {
                    NormalProjectioni += ivec(iloc,0)*normal[iloc];
                }
                
                //para a equacao do fluxo
                ef(iq,0) += (-1.0)*weight*Pd*NormalProjectioni;
            }

            // fim dirichlet

        }
            break;
        case 1:  // Neumann
        {
//            REAL InvK = 1./fK;
            for (int iq = 0; iq<phrq; iq++)
            {
                int ivecind = dataleft[0].fVecShapeIndex[iq].first;
                int ishapeind = dataleft[0].fVecShapeIndex[iq].second;
                TPZFNMatrix<3> ivec(3,1);
                ivec(0,0) = dataleft[0].fDeformedDirections(0,ivecind);
                ivec(1,0) = dataleft[0].fDeformedDirections(1,ivecind);
                ivec(2,0) = dataleft[0].fDeformedDirections(2,ivecind);
                ivec *= phiQ(ishapeind,0);
                
                
                REAL NormalProjectioni = 0.;
                for(int iloc=0; iloc<fDim; iloc++)
                {
                    NormalProjectioni += ivec(iloc,0)*normal[iloc];
                }
                ef(iq,0) += gBigNumber*weight*(Qn)*NormalProjectioni;
                //ef(iq,0) += gBigNumber*weight*(ValorPhin - Qn)*NormalProjectioni;

                for (int jq=0; jq<phrq; jq++)
                {
                    TPZFNMatrix<3> jvec(3,1);
                    int jvecind = dataleft[0].fVecShapeIndex[jq].first;
                    int jshapeind = dataleft[0].fVecShapeIndex[jq].second;
                    jvec(0,0) = dataleft[0].fDeformedDirections(0,jvecind);
                    jvec(1,0) = dataleft[0].fDeformedDirections(1,jvecind);
                    jvec(2,0) = dataleft[0].fDeformedDirections(2,jvecind);
                    
                    jvec *= phiQ(jshapeind,0);

                    REAL NormalProjectionj = 0.;
                    for(int iloc=0; iloc<fDim; iloc++)
                    {
                        NormalProjectionj += jvec(iloc,0)*normal[iloc];
                    }
                    
                    ek(iq,jq) += gBigNumber*weight*NormalProjectioni*NormalProjectionj;

                }
            }
//            //termo fonte referente a equacao da pressao no entra!!!!
//            for (int jp = 0; jp < phrp ; jp++)
//            {
//                TPZFNMatrix<3> jvec(3,1);
//                int jvecind = dataleft[0].fVecShapeIndex[jp].first;
//                int jshapeind = dataleft[0].fVecShapeIndex[jp].second;
//                jvec(0,0) = dataleft[0].fDeformedDirections(0,jvecind);
//                jvec(1,0) = dataleft[0].fDeformedDirections(1,jvecind);
//                jvec(2,0) = dataleft[0].fDeformedDirections(2,jvecind);
//                
//                jvec *= phiQ(jshapeind,0);
//                
//                REAL NormalProjectionj = 0.;
//                for(int iloc=0; iloc<fDim; iloc++)
//                {
//                    NormalProjectionj += jvec(iloc,0)*normal[iloc];
//                }
//
//                ef(jp,0) += gBigNumber*weight*Qn*NormalProjectionj;
//            }

            // fim neumann
        }
            break;
        case 3:  // Robin
        {
            std::cout << " Robin Nao implementada " << std::endl;
            DebugStop();
        }
            break;
        default:
        {
            std::cout << " Nao implementada " << std::endl;
            DebugStop();
        }
            break;
    }
    
}


void TPZMatPoissonD3::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
	DebugStop();
}


void TPZMatPoissonD3::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
	DebugStop();
}

void TPZMatPoissonD3::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec)
{
    int nref = datavec.size();
	for(int i = 0; i<nref; i++)
	{
        datavec[i].fNeedsSol = false;//true;
		datavec[i].fNeedsNormal = false;//true;
	}
}

/** Returns the variable index associated with the name */
int TPZMatPoissonD3::VariableIndex(const std::string &name){
    if(!strcmp("Flux",name.c_str()))           return 1;
    if(!strcmp("Pressure",name.c_str()))       return 2;
    if(!strcmp("GradFluxX",name.c_str()))      return 3;
    if(!strcmp("GradFluxY",name.c_str()))      return 4;
    if(!strcmp("GradFluxZ",name.c_str()))      return 5;
    if(!strcmp("ExactPressure",name.c_str()))  return 6;
    if(!strcmp("ExactFlux",name.c_str()))      return 7;
    if(!strcmp("Rhs",name.c_str()))            return 8;
    if(!strcmp("GradP",name.c_str()))          return 9;
	
	return TPZMaterial::VariableIndex(name);
}

int TPZMatPoissonD3::NSolutionVariables(int var){
    if(var == 1) return 3;
    if(var == 2) return 1;
    if(var == 3) return 3;
    if(var == 4) return 3;
    if(var == 5) return 3;
    if(var == 6) return 1;
    if(var == 7) return 3;
    if(var == 8) return 1;
    if(var == 9) return 3;
	return TPZMaterial::NSolutionVariables(var);
}

// metodo para gerar vtk
void TPZMatPoissonD3::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout){
	
	Solout.Resize( this->NSolutionVariables(var));
	
	TPZVec<STATE> SolP, SolQ;
    
    // SolQ = datavec[0].sol[0];
    SolP = datavec[1].sol[0];
    
    if(var == 1){ //function (state variable Q)
        for (int ip = 0; ip<fDim; ip++)
        {
            Solout[ip] = datavec[0].sol[0][ip];
        }
        
//		Solout[0] = datavec[0].sol[0][0];
//        Solout[1] = datavec[0].sol[0][1];
//        Solout[2] = datavec[0].sol[0][2];
		return;
	}
    
    if(var == 2){
		Solout[0] = SolP[0];//function (state variable p)
		return;
	}
    
    if(var==3){
        for (int ip = 0; ip<fDim; ip++)
        {
            Solout[ip] = datavec[0].dsol[0](ip,0);
        }
        return;
    }
    
    if(var==4){
        for (int ip = 0; ip<fDim; ip++)
        {
            Solout[ip] = datavec[0].dsol[0](ip,1);
        }
        return;
    }
    
    if(var==5){
        for (int ip = 0; ip<fDim; ip++)
        {
            Solout[ip] = datavec[0].dsol[0](ip,2);
        }
        return;
    }
    
    TPZVec<REAL> ptx(3);
	TPZVec<STATE> solExata(1);
	TPZFMatrix<STATE> flux(fDim,1);
    
    //Exact soluion
	if(var == 6){
		fForcingFunctionExact->Execute(datavec[1].x, solExata,flux);
		Solout[0] = solExata[0];
		return;
	}//var6
    
    if(var == 7){
		fForcingFunctionExact->Execute(datavec[0].x, solExata,flux);
        for (int ip = 0; ip<fDim; ip++)
        {
            Solout[ip] = flux(ip,0);
        }
		return;
	}//var7
    
    if(var == 8){
        REAL force = fF;
        if(fForcingFunction) {
            TPZManVector<STATE> res(1);
            fForcingFunction->Execute(datavec[1].x,res);
            force = res[0];
        }
        Solout[0] = force;
        
        return;
    }//var8
    if(var==9){
        
        TPZFMatrix<STATE> GradofP;
        const TPZFMatrix<STATE> &dsol = datavec[1].dsol[0];
        TPZAxesTools<STATE>::Axes2XYZ(dsol, GradofP, datavec[1].axes);
        //        int nc = GradofP.Cols();
        //        int nl = GradofP.Rows();
        
        for (int ip = 0; ip<fDim; ip++)
        {
            Solout[ip] = -1.0*GradofP(ip,0);
        }
        return;
    }
}

// metodo para computar erros
void TPZMatPoissonD3::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout){
    
    Solout.Resize( this->NSolutionVariables(var));

    
    if(var == 1){ //function (state variable Q)
        for (int ip = 0; ip<fDim; ip++)
        {
            Solout[ip] = data.sol[0][ip];
        }
        
        return;
    }
    
    if(var == 2){ //function (state variable p)
        
        TPZVec<STATE> SolP;
        SolP = data.sol[0];
        
        Solout[0] = SolP[0];
        return;
    }
    
    
}

#include "pzaxestools.h"
void TPZMatPoissonD3::Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes, int var,TPZVec<STATE> &Solout){
    
#ifndef STATE_COMPLEX
    Solout.Resize( this->NSolutionVariables( var ) );
    
    if(var == 1){
        int id;
        for(id=0 ; id<fDim; id++) {
            TPZFNMatrix<9,STATE> dsoldx;
            TPZAxesTools<STATE>::Axes2XYZ(DSol, dsoldx, axes);
            Solout[id] = dsoldx(id,0);//derivate
        }
        return;
    }
    if(var == 2) {
        Solout[0] = Sol[0];//function
        return;
    }//var == 2
    
#endif
    TPZMaterial::Solution(Sol, DSol, axes, var, Solout);
    
}//method



void TPZMatPoissonD3::FillDataRequirements(TPZVec<TPZMaterialData > &datavec)
{
	int nref = datavec.size();
	for(int i = 0; i<nref; i++)
	{
		datavec[i].SetAllRequirements(false);
        datavec[i].fNeedsSol = true;
		datavec[i].fNeedsNeighborSol = false;
		datavec[i].fNeedsNeighborCenter = false;
		datavec[i].fNeedsNormal = true;
	}
}




void TPZMatPoissonD3::ErrorsHdiv(TPZMaterialData &data,TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values){
    
    values.Fill(0.0);
    TPZVec<STATE> sol(1),dsol(3),div(1);
//    if(data.numberdualfunctions) Solution(data,2,sol);//pressao
    Solution(data,1,dsol);//fluxo
    //Solution(data,14,div);//divergente
    
#ifdef LOG4CXX
//    if(logger->isDebugEnabled()){
//        std::stringstream sout;
//        sout<< "\n";
//        sout << " Pto  " << data.x << std::endl;
//        sout<< " pressao exata " <<u_exact <<std::endl;
//        sout<< " pressao aprox " <<sol <<std::endl;
//        sout<< " ---- "<<std::endl;
//        sout<< " fluxo exato " <<du_exact(0,0)<<", " << du_exact(1,0)<<std::endl;
//        sout<< " fluxo aprox " <<dsol<<std::endl;
//        sout<< " ---- "<<std::endl;
//        if(du_exact.Rows()>fDim) sout<< " div exato " <<du_exact(2,0)<<std::endl;
//        sout<< " div aprox " <<div<<std::endl;
//        LOGPZ_DEBUG(logger,sout.str())
//    }
#endif
    
    
//    //values[0] : pressure error using L2 norm
//    if(data.numberdualfunctions){
//        REAL diffP = abs(u_exact[0]-sol[0]);
//        values[0]  = diffP*diffP;
//    }
    //values[1] : flux error using L2 norm
    for(int id=0; id<fDim; id++) {
        REAL diffFlux = abs(dsol[id] - du_exact(id,0));
        values[1]  += diffFlux*diffFlux;
    }
//    if(du_exact.Rows()>3){
//        //values[2] : divergence using L2 norm
//        REAL diffDiv = abs(div[0] - du_exact(2,0));
//        values[2]=diffDiv*diffDiv;
//        //values[3] : Hdiv norm => values[1]+values[2];
//        values[3]= values[1]+values[2];
//    }
}



void TPZMatPoissonD3::Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,
							 TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes,
							 TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values) {
	
    values.Resize(NEvalErrors());
    values.Fill(0.0);
    
    TPZManVector<STATE> sol(1),dsol(3,0.);
    Solution(u,dudx,axes,2,sol);
    Solution(u,dudx,axes,1,dsol);
    int id;
    //values[1] : eror em norma L2
    REAL  diff = fabs(sol[0] - u_exact[0]);
    values[1]  = diff*diff;
    //values[2] : erro em semi norma H1
    values[2] = 0.;
    for(id=0; id<fDim; id++) {
        diff = fabs(dsol[id] - du_exact(id,0));
        values[2]  += abs(fK)*diff*diff;
    }
    //values[0] : erro em norma H1 <=> norma Energia
    values[0]  = values[1]+values[2];
}






