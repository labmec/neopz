/**
 * \file
 * @brief Contains implementations of the TPZMultCamada methods.
 */

#include "tpzmultcamada.h"
#include "pzmatplaca2.h"

void TPZMultCamada::Solution(TPZVec < STATE > & Sol, TPZFMatrix<STATE> & DSol, TPZFMatrix<REAL> & axes, int var, TPZVec < STATE > & Solout) {

	if(var == 2) {
		Solout.Resize(3);
		Solout[0] = Sol[0];
		Solout[1] = Sol[1];
		Solout[2] = Sol[2];
		return;
	}
	if(var == 4) {
		Solout.Resize(1);
		Solout[0] = Sol[2];
		return;
	}
	int i, j, idf;
	if (var >= 50 && var <= 57) {
		TPZVec < STATE > SoloutAcc(3, 0.);
		REAL Ma[3];
		REAL aut1,aut2;
		Solout.Resize(3);
		for (j = 0; j < 3; j++) Ma[j] = 0.;
		TPZMaterialData data;
		data.sol[0] = Sol;
		data.dsol[0] = DSol;
		data.axes = axes;
		for (i = 0; i < fCamadas.NElements(); i++) {
			if(var < 54) fCamadas[i]->Solution(data, 50, SoloutAcc);
			else fCamadas[i]->Solution(data, 54, SoloutAcc);
			for (j = 0; j < 3; j++) {
				Ma[j] += SoloutAcc[j];
				SoloutAcc[j] = 0.;
			}
		}
		REAL sqrtval = sqrt((Ma[0]-Ma[1])*(Ma[0]-Ma[1])+Ma[2]*Ma[2]*4);
		REAL vec1[2],vec2[2],vecnorm1,vecnorm2;
		aut1 = (Ma[0]+Ma[1]-sqrtval)*0.5;
		aut2 = (Ma[0]+Ma[1]+sqrtval)*0.5;
		vec1[0] = Ma[0]-Ma[1]-sqrtval;
		vec1[1] = 2*Ma[2];
		vecnorm1 = sqrt(vec1[0]*vec1[0]+vec1[1]*vec1[1]);
		if(vecnorm1 >= 1.e-8) {
			vec1[0] /= vecnorm1;
			vec1[1] /= vecnorm1;
		} else {
			vec1[0] = 0.;
			vec1[1] = 1.;
		}
		vec2[0] = Ma[0]-Ma[1]+sqrtval;
		vec2[1] = 2.*Ma[2];
		vecnorm2 = sqrt(vec2[0]*vec2[0]+vec2[1]*vec2[1]);
		if(vecnorm2 >= 1.e-8) {
			vec2[0] /= vecnorm2;
			vec2[1] /= vecnorm2;
		} else {
			vec2[0] = 0.;
			vec2[1] = 1.;
		}
		
		// E preciso calcular os valores principais que encontram=se em SoloutLocal e orientar pelos vetores axes
		if(var == 50 || var == 54) {
			Solout.Resize(3);
			for(idf=0; idf<3; idf++) Solout[idf] = vec1[0] * axes(0,idf) + vec1[1]*axes(1,idf);
			return;
		}
		if(var == 51 || var ==55) {
			Solout.Resize(3);
			for(idf=0; idf<3; idf++) Solout[idf] = vec2[0] * axes(0,idf) + vec2[1]*axes(1,idf);
			return;
		}
		
		if(var == 52 || var == 53 || var == 56 || var == 57) {
			Solout.Resize(1);
			Solout[0] = (var == 52 || var == 56) ? aut1 : aut2;
			return;
		}
		
	}
	TPZMaterial::Solution(Sol,DSol,axes,var,Solout);
}

void TPZMultCamada::Contribute(TPZMaterialData &data,
                               REAL weight,
                               TPZFMatrix<STATE> &ek,
                               TPZFMatrix<STATE> &ef) {
	
	int i;
	for (i = 0; i < fCamadas.NElements(); i++) 
		fCamadas[i]->Contribute(data,weight,ek,ef);
}

void TPZMultCamada::ContributeBC(TPZMaterialData &data,
                                 REAL weight, 
                                 TPZFMatrix<STATE> &ek, 
                                 TPZFMatrix<STATE> &ef, 
                                 TPZBndCond &bc){
	if(fCamadas.NElements()) fCamadas[0]->ContributeBC(data,weight,ek,ef,bc);
}

int TPZMultCamada::VariableIndex(const std::string &name) {
	if(!strcmp(name.c_str(),"Displacement")) return 2;
	if(!strcmp(name.c_str(),"Deslocz")) return 4;
	if(!strcmp(name.c_str(),"M1Vec"))     return 50;// Momento na direcao do primeiro  eixo principal
	if(!strcmp(name.c_str(),"M2Vec"))     return 51;// Momento na direcao do segundo  eixo principal
	if(!strcmp(name.c_str(),"M1Scal"))     return 52;// Valor do primeiro momento principal
	if(!strcmp(name.c_str(),"M2Scal"))     return 53;// Valor do segundo momento principal
	if(!strcmp(name.c_str(),"N1Vec"))     return 54;// Momento na direcao do primeiro  eixo principal
	if(!strcmp(name.c_str(),"N2Vec"))     return 55;// Momento na direcao do segundo  eixo principal
	if(!strcmp(name.c_str(),"N1Scal"))     return 56;// Valor do primeiro momento principal
	if(!strcmp(name.c_str(),"N2Scal"))     return 57;// Valor do segundo momento principal
	
	return TPZMaterial::VariableIndex(name);
}

int TPZMultCamada::NSolutionVariables(int var) {
    if(var == 2) return 3;
    if(var == 4) return 1;
	if(var == 50 || var == 51 || var == 54 || var == 55) return 3;
	if(var == 52 || var == 53 || var == 56 || var == 57) return 1;
	return TPZMaterial::NSolutionVariables(var);
}

int TPZMultCamada::NStateVariables() const {
	if(fCamadas.NElements()) return fCamadas[0]->NStateVariables();
	return 0;
}


int TPZMultCamada::ClassId() const{
    return Hash("TPZMultCamada") ^ TPZMaterial::ClassId() << 1;
}
