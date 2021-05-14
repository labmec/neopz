/*
 * @file
 * @author Denise de Siqueira
 * @brief Contains the implementation of the TPZMaterialCoupling methods.
 */

#include "pzmaterialcoupling.h"
#include "pzmaterialdata.h"
#include "pzmatrix.h"
#include "pzlog.h"
#include "pzpoisson3d.h"
#include "pzgeoel.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.material.poisson3d.materialcoupling");
#endif

TPZMaterialCoupling::TPZMaterialCoupling(int nummat, int dim):
TPZRegisterClassId(&TPZMaterialCoupling::ClassId),
TPZMatPoisson3d( nummat, dim){
	this->NStateVariables();
	
	
}

TPZMaterialCoupling::TPZMaterialCoupling() :
TPZRegisterClassId(&TPZMaterialCoupling::ClassId){
}

TPZMaterialCoupling::~TPZMaterialCoupling() {
}

void TPZMaterialCoupling::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft,TPZMaterialData &dataright, REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){
	
	TPZFMatrix<REAL>  &phiH1 = dataright.phi;
	TPZFMatrix<REAL>  &phiHdiv = dataleft.phi;		
	int numvec=dataleft.fVecShapeIndex.NElements();
	int nrowH1=phiH1.Rows();//Funcao a direita H1
	int numdual = dataleft.numberdualfunctions;
	
	TPZFMatrix<REAL> ekCouple(numvec+numdual,nrowH1,0.);
	//vou precisar da  orientacao das normais na interface
	
	REAL leftX0=data.normal[0];
	REAL leftX1=data.normal[1];
	REAL leftX2=data.normal[2];
	
	
	for(int ilinha=0; ilinha<numvec; ilinha++) {
		int ivecind = dataleft.fVecShapeIndex[ilinha].first;
		int ishapeind = dataleft.fVecShapeIndex[ilinha].second;
		REAL prod=dataleft.fDeformedDirections(0,ivecind)*leftX0+dataleft.fDeformedDirections(1,ivecind)*leftX1+dataleft.fDeformedDirections(2,ivecind)*leftX2;
		
		
				for(int jcol=0; jcol<nrowH1; jcol++) {
						
						REAL prod1 =	phiHdiv(ishapeind,0)*phiH1(jcol,0)*prod;
						
						
#ifdef PZ_LOG
               if (logger.isDebugEnabled())
					{
							std::stringstream sout;
								sout << "prod phiHdiv[ " <<ishapeind << "]= " << phiHdiv(ishapeind,0)<< " phiH1[ "<< jcol << "] = " << phiH1(jcol,0)<< std::endl;
								LOGPZ_DEBUG(logger, sout.str().c_str());
						}
#endif
						ekCouple(ilinha,jcol)+= weight  * prod1;
						ek(ilinha,numdual+ numvec+jcol) += weight  * (prod1);
						
#ifdef PZ_LOG
                    if (logger.isDebugEnabled())
					{
				std::stringstream sout;
				sout << "-- PosJ " << numdual+ numvec+jcol<< std::endl;
				LOGPZ_DEBUG(logger, sout.str().c_str());
			}
#endif
			ek(jcol+numvec+numdual,ilinha) += weight  *(-prod1);
		}
    }
	
	ekCouple.Print("Matriz teste Acoplamento",std::cout);
#ifdef PZ_LOG
	{
		std::stringstream sout;
		ekCouple.Print("Matriz teste Acoplamento",sout);
		LOGPZ_DEBUG(logger, sout.str().c_str());
	}
#endif
	
}

void TPZMaterialCoupling::ContributeInterface2(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, 
											   REAL weight,TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef){
	
	TPZFMatrix<REAL>  &phixL = dataleft.phi;
	
	TPZFMatrix<REAL>  &phixR = dataright.phi;
	
	
	
	int numvec=dataright.fVecShapeIndex.NElements();
	int nrowR=phixR.Rows();//funcao a direita Hdiv
	int nrowL=phixL.Rows();//Funcao a esquerda H1
	int numdual = dataright.numberdualfunctions;
	
	std::cout << "numero de funcoes de Hdiv( direita ) " << nrowR<<std::endl;
	std::cout << "numero de funcoes de de pressao(direita) " << numdual<<std::endl;
	std::cout << "numero de funcoes de H1 (esquerda ) " << nrowL<<std::endl;
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
		{
				std::stringstream sout;
				sout << "numero de funcoes de Hdiv( direita ) " << nrowR<<std::endl;
				sout << "numero de funcoes de de pressao(direita) " << numdual<<std::endl;
				sout << "numero de funcoes de H1 (esquerda ) " << nrowL<<std::endl;
				LOGPZ_DEBUG(logger, sout.str().c_str());
		}
#endif
		
	for(int ir=0; ir<nrowR-1; ir++) {
		int ishapeind = dataright.fVecShapeIndex[ir].second;
		for(int jl=0; jl<nrowL; jl++) {
			REAL prod1 =	phixR(ishapeind,0)* phixL(jl);
#ifdef PZ_LOG
                    if (logger.isDebugEnabled())
			{
				std::stringstream sout;
				sout << "produto das phis " << prod1<<std::endl;
				LOGPZ_DEBUG(logger, sout.str().c_str());
			}
#endif
			ek(ir,numvec+jl) += weight  * prod1;
			ek(numvec+jl,ir) += weight  *(-prod1);
		}
    }
		
		
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
		{
				std::stringstream sout;
				ek.Print("Matriz de Acoplamento",sout);
				LOGPZ_DEBUG(logger, sout.str().c_str());
		}
#endif
}

void TPZMaterialCoupling::InitMaterialData(TPZMaterialData &data){
	return ;		
}

int TPZMaterialCoupling::ClassId() const{
    return Hash("TPZMaterialCoupling") ^ TPZMatPoisson3d::ClassId() << 1;
}
