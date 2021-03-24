/**
 * \file
 * @brief Contains implementations of the TPZGradientFlux methods.
 */

#include "pzgradientflux.h"
#include "pzmaterialdata.h"

TPZGradientFlux::TPZGradientFlux(){
	
}

TPZGradientFlux::TPZGradientFlux(const TPZGradientFlux &cp){
	
}

TPZGradientFlux::~TPZGradientFlux(){
	
}

void TPZGradientFlux::ComputeFlux(TPZVec<STATE> &solL, TPZVec<STATE> &solR, const TPZVec<REAL> &normal, TPZVec<STATE> & F){
	const int nstate = 5;
	const int dim = 3;
	F.Resize(nstate*dim);
	F.Fill(0.);
	for(int i = 0; i < nstate; i++){
		for(int id = 0; id < dim; id++){
			F[i*dim+id] = 0.5*(solL[i]+solR[i])*normal[id];
		}//for id
	}//for i
}//void

void TPZGradientFlux::ApplyLimiter(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright){
    int numbersol = dataleft.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }

	const int nstate = 5;
	const int dim = 3;
	TPZManVector<STATE,dim> dL(dim), dR(dim);
	for(int id = 0; id < dim; id++){
		dL[id] = data.x[id] - dataleft.XCenter[id];
		dR[id] = data.x[id] - dataright.XCenter[id];
	}//for id
	
	TPZManVector<STATE,dim> gradL(3), gradR(3);
	for(int is = 0; is < nstate; is++){
		for(int id = 0; id < dim; id++){
			gradL[id] = dataleft.sol[0][ nstate + is*dim +id ];
			gradR[id] = dataright.sol[0][ nstate + is*dim +id ];
		}//for id
        TPZManVector<STATE,3> stnormal(3);
        for (int i=0; i<3; i++) {
            stnormal[i] = data.normal[i];
        }
		this->ApplyVanAlbadaLimiter(dataleft.sol[0][is], dataright.sol[0][is], gradL, gradR, stnormal, dL, dR);
		//     this->ApplyMinModLimiter(data.soll[is], data.solr[is], gradL, gradR, data.normal, dL, dR);
	}//for is
}//void

void TPZGradientFlux::ApplyMinModLimiter(STATE &soll, STATE &solr,
                                         const TPZVec<STATE>& gradL, const TPZVec<STATE> &gradR,
                                         const TPZVec<STATE> &normal, 
                                         const TPZVec<STATE> &dL, const TPZVec<STATE> & dR){
	const STATE sL = this->Dot(gradL,normal);
	const STATE sR = this->Dot(gradR,normal);
	if((sL*sR) < 0.){
		return;
	}
	const STATE k = 0.5;
	if(fabs(sL) < fabs(sR)){
		soll = soll + k*this->Dot(gradL,dL);
		solr = solr + k*this->Dot(gradL,dR);
		return;
	}
	if(fabs(sL) > fabs(sR)){
		soll = soll + k*this->Dot(gradR,dL);
		solr = solr + k*this->Dot(gradR,dR);
		return;
	}
	//if (sL == sR)
	soll = soll + k*this->Dot(gradL,dL);
	solr = solr + k*this->Dot(gradR,dR);
}//void

void TPZGradientFlux::ApplyVanAlbadaLimiter(STATE &soll, STATE &solr,
                                            const TPZVec<STATE>& gradL, const TPZVec<STATE> &gradR,
                                            const TPZVec<STATE> &normal, 
                                            const TPZVec<STATE> &dL, const TPZVec<STATE> & dR){
	const double sL = this->Dot(gradL,normal);
	const double sR = this->Dot(gradR,normal);
	if(fabs(sL) < 1e-12){
		return;
	}
	const double r = sR/sL;
	double phi = (r <= 0.) ? 0. : ((r*r+r)/(1.+r*r));
	phi *= 0.5;
	soll = soll + phi*this->Dot(gradL,dL);
	solr = solr + phi*this->Dot(gradR,dR);
	
}//void
