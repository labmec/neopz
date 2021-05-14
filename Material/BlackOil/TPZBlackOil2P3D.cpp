/**
 * @file
 * @brief Contains implementations of the TPZBlackOil2P3D methods.
 */

#include "BlackOil/TPZBlackOil2P3D.h"
#include "TPZBndCondT.h"

using namespace std;

TPZBlackOil2P3D::EState TPZBlackOil2P3D::gState = ELastState;

/** Linear interpolation */
void TPZBlackOil2P3D::Interpolate(std::map<REAL,REAL> &dados, STATE x, STATE &y, STATE &dy){
	STATE x0, xL, y0, yL;
	std::map< REAL, REAL >::iterator w;
	w = dados.lower_bound(x);
	
	if(dados.size() < 2){
		std::cout << "Error at " << __PRETTY_FUNCTION__ << "\n";
		std::cout.flush();
	}
	
	if(w == dados.end()){
		w--;
		y = w->second;
		
		//derivada
		xL = w->first;
		yL = w->second;
		w--;
		x0 = w->first;
		y0 = w->second;
		dy = (yL-y0)/(xL-x0);
		
		return;
	}
	if (w == dados.begin()){
		y = w->second;
		
		//derivada
		x0 = w->first;
		y0 = w->second;
		w++;
		xL = w->first;
		yL = w->second;
		dy = (yL-y0)/(xL-x0);
		return;
	}
	
	xL = w->first;
	yL = w->second;
	w--;
	
	x0 = w->first;
	y0 = w->second;
	
	y = (yL-y0)*(x-x0)/(xL-x0)+y0;
	dy = (yL-y0)/(xL-x0);
	
}

void TPZBlackOil2P3D::Interpolate(std::map<REAL,REAL> &dados, BFadREAL x, BFadREAL &y){
	REAL x0, xL, y0, yL;
	std::map< REAL, REAL >::iterator w;
	w = dados.lower_bound(x.val());
	
	if(dados.size() < 2){
		std::cout << "Error at " << __PRETTY_FUNCTION__ << "\n";
		std::cout.flush();
	}
	
	if(w == dados.end()){
		w--;
		y = w->second;
		
		return;
	}
	if (w == dados.begin()){
		y = w->second;
		
		return;
	}
	
	xL = w->first;
	yL = w->second;
	w--;
	
	x0 = w->first;
	y0 = w->second;
	y = (yL-y0)*(x-x0)/(xL-x0)+y0;
}

/** Relative oil permeability \f$ Kro = Kro( Sw ) \f$ */
void TPZBlackOil2P3D::Kro(STATE So, STATE &Kro, STATE &dKroSo) {
	const int n = 8;
	STATE tabela[n][2] = {{0.12,1.},{0.2,0.8},{0.3,0.6},{0.4,0.45},{0.5,0.3},{0.6,0.2},{0.7,0.1},{0.82,0.}};
	std::map<REAL,REAL> dados;
	for(int i = 0; i < n; i++){
		dados[ 1.-tabela[i][0] ] = tabela[i][1];
	}
	this->Interpolate(dados,So,Kro,dKroSo);
}

void TPZBlackOil2P3D::Kro(BFadREAL So, BFadREAL &Kro) {
	const int n = 8;
	STATE tabela[n][2] = {{0.12,1.},{0.2,0.8},{0.3,0.6},{0.4,0.45},{0.5,0.3},{0.6,0.2},{0.7,0.1},{0.82,0.}};
	std::map<REAL,REAL> dados;
	for(int i = 0; i < n; i++){
		dados[ 1.-tabela[i][0] ] = tabela[i][1];
	}
	this->Interpolate(dados,So,Kro);
}

/** Relative water permeability \f$ Krw = Krw( Sw ) \f$ */
void TPZBlackOil2P3D::Krw(STATE So, STATE &Krw, STATE &dKrwSo) {
	const int n = 8;
	STATE tabela[n][2] = {{0.12,0.},{0.2,0.1},{0.3,0.2},{0.4,0.3},{0.5,0.4},{0.6,0.55},{0.7,0.7},{0.82,1.}};
	std::map<REAL,REAL> dados;
	for(int i = 0; i < n; i++){
		dados[ 1.-tabela[i][0] ] = tabela[i][1];
	}
	this->Interpolate(dados,So,Krw,dKrwSo);
}

void TPZBlackOil2P3D::Krw(BFadREAL So, BFadREAL &Krw){
	const int n = 8;
	STATE tabela[n][2] = {{0.12,0.},{0.2,0.1},{0.3,0.2},{0.4,0.3},{0.5,0.4},{0.6,0.55},{0.7,0.7},{0.82,1.}};
	std::map<REAL,REAL> dados;
	for(int i = 0; i < n; i++){
		dados[ 1.-tabela[i][0] ] = tabela[i][1];
	}
	this->Interpolate(dados,So,Krw);
}

/** Bo = Bo( po )
 */
void TPZBlackOil2P3D::Bo(STATE po, STATE &Bo, STATE &dBoDpo) {
	const int n = 9;
	STATE tabela[n][2] = {{14.7, 1.062}, {264.7, 1.15}, {514.7, 1.207}, {1014.7, 1.295}, {2014.7, 1.435}, {2514., 1.5}, 
		{3014.7, 1.565}, {4014.7, 1.695}, {9014.7, 1.695}};
	const STATE conv = 6894.757;
	std::map<REAL,REAL> dados;
	for(int i = 0; i < n; i++){
		dados[ tabela[i][0]*conv ] = tabela[i][1];
	}
	this->Interpolate(dados,po,Bo,dBoDpo);
}

void TPZBlackOil2P3D::Bo(BFadREAL po, BFadREAL &Bo) {
	const int n = 9;
	STATE tabela[n][2] = {{14.7, 1.062}, {264.7, 1.15}, {514.7, 1.207}, {1014.7, 1.295}, {2014.7, 1.435}, {2514., 1.5}, 
		{3014.7, 1.565}, {4014.7, 1.695}, {9014.7, 1.695}};
	const STATE conv = 6894.757;
	std::map<REAL,REAL> dados;
	for(int i = 0; i < n; i++){
		dados[ tabela[i][0]*conv ] = tabela[i][1];
	}
	this->Interpolate(dados,po,Bo);
}

/** Oil viscosity depending on pressure */
void TPZBlackOil2P3D::OilVisc(STATE po, STATE &OilVisc, STATE &dOilViscDpo){
	const int n = 9;
	STATE tabela[n][2] = {{14.7, 1.062}, {264.7, 1.15}, {514.7, 1.207}, {1014.7, 1.295}, {2014.7, 1.435}, {2514., 1.5}, 
		{3014.7, 1.565}, {4014.7, 1.695}, {9014.7, 1.6}};
	const STATE conv = 6894.757;
	const STATE convVisc = 1e-3;
	std::map<REAL,REAL> dados;
	for(int i = 0; i < n; i++){
		dados[ tabela[i][0]*conv ] = tabela[i][1]*convVisc;
	}
	this->Interpolate(dados,po,OilVisc,dOilViscDpo);
}

void TPZBlackOil2P3D::OilVisc(BFadREAL po, BFadREAL &OilVisc){
	
	const int n = 10;
	STATE tabela[n][2] = {{14.7,1.04},{264.7,0.975},{514.7,0.91},{1014.7,0.83},{2014.7,0.695},
		{2514.,0.641},{3014.7,0.594},{4014.7,0.51},{5014.7,0.449},{9014.7,0.203}};
	const STATE conv = 6894.757;
	const STATE convVisc = 1e-3;
	std::map<REAL,REAL> dados;
	for(int i = 0; i < n; i++){
		dados[ tabela[i][0]*conv ] = tabela[i][1]*convVisc;
	}
	this->Interpolate(dados,po,OilVisc);
}

/** Capilar pressure \f$ pc = pc( Sw ) \f$ */
void TPZBlackOil2P3D::CapilarPressure(STATE So, STATE &pc, STATE &DpcDSo){
	pc = 0.;
	DpcDSo = 0.;
}

void TPZBlackOil2P3D::CapilarPressure(BFadREAL So, BFadREAL &pc){
	pc = 0.;
}

/** Porosity \f$ Phi = Phi( pw ) \f$ - fizemos como Phi ( po ) */
void TPZBlackOil2P3D::Porosity(STATE po, STATE &poros, STATE &dPorosDpo){
	const STATE comp = 3.625943e-10;
	const STATE pref = 101352.93;
	const STATE porosRef = 0.22;
	poros = porosRef*exp(comp*(po-pref));
	dPorosDpo = comp*porosRef*exp(comp*(po-pref));
}

void TPZBlackOil2P3D::Porosity(BFadREAL po, BFadREAL &poros){
	const STATE comp = 3.625943e-10;
	const STATE pref = 101352.93;
	const STATE porosRef = 0.22;
	poros = porosRef*exp(comp*((po.val())-pref));
}

// Constant data

/** Densidade do oleo em condicoes padroes - kg/m3 */
STATE TPZBlackOil2P3D::OilRhoSC() const{
	return 740.75782;
}

/** Densidade da agua em condicoes padroes - kg/m3 */
STATE TPZBlackOil2P3D::WaterRhoSC() const{
	return 996.95712;
}

/** Aceleracao da gravidade */
STATE TPZBlackOil2P3D::g() const{
	return 9.81;
}

/** Bw = constante */
STATE TPZBlackOil2P3D::Bw() const{
	return 1.041;
}

/** Viscosidade da agua (constante) */
STATE TPZBlackOil2P3D::WaterVisc() const{
	return 0.31e-3;
}

/** Permeabilidade absoluta */
void TPZBlackOil2P3D::K(TPZFMatrix<REAL> &K){
	K.Resize(3,3);
	K.Zero();
	K(0,0) = 2.96077e-13;//10;
	K(1,1) = 2.96077e-13;//10;
	K(2,2) = 4.93462e-14;//11;
}


//Programa

TPZBlackOil2P3D::TPZBlackOil2P3D(int id, STATE deltaT):TBase(id){
	this->fDeltaT = deltaT;
}


TPZMaterial * TPZBlackOil2P3D::NewMaterial() const{
	return new TPZBlackOil2P3D(*this);
}

void TPZBlackOil2P3D::Contribute(const TPZMaterialDataT<STATE> &data, REAL weight,
                                 TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
	
	//un ou un+1
	STATE stateVal = 0.;
	if(gState == ELastState) stateVal = -1.;
	if(gState == ECurrentState) stateVal = +1.;
	
	//pressao e saturacao
	const BFadREAL po(data.sol[0][0],0);
	const BFadREAL So(data.sol[0][1],1);
	BFadREAL pc;
	this->CapilarPressure(So,pc);
	const BFadREAL pw = po - pc;
	const BFadREAL Sw = ((BFadREAL)1.)-So;
	
	//porosidade
	BFadREAL porosidade;
	this->Porosity(po,porosidade);
	
	//fator volume formacao
	BFadREAL Bo;
	this->Bo(po,Bo);
	const STATE Bw = this->Bw();
	
	//Equacao 1
	BFadREAL VolOp1 = (porosidade*So/Bo)/this->fDeltaT;
	
	//Equacao 2
	BFadREAL VolOp2 = (porosidade*Sw/Bw)/this->fDeltaT;
	
	ef(0,0) += -1.*weight*stateVal*VolOp1.val();
	ef(1,0) += -1.*weight*stateVal*VolOp2.val();
	
	if(gState == ECurrentState){//Last state has no tangent
		ek(0,0) += +1.*weight*stateVal*VolOp1.dx(0);
		ek(0,1) += +1.*weight*stateVal*VolOp1.dx(1);
		
		ek(1,0) += +1.*weight*stateVal*VolOp2.dx(0);
		ek(1,1) += +1.*weight*stateVal*VolOp2.dx(1);
	}
	
}//method

void TPZBlackOil2P3D::ContributeBC(const TPZMaterialDataT<STATE> &data, REAL weight,
                      TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef,
                      TPZBndCondT<STATE> &bc) {
    PZError<<__PRETTY_FUNCTION__<<'\n';
    PZError<<"Should not be called. Aborting...\n";
    DebugStop();
}//method

void TPZBlackOil2P3D::ContributeInterface(const TPZMaterialDataT<STATE> &data,
                                          const TPZMaterialDataT<STATE> &dataleft,
                                          const TPZMaterialDataT<STATE> &dataright,
                                          REAL weight,
                                          TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
	
	if(gState == ELastState) return;
	
	//calculando distancia entre centro dos elementos
	STATE dist = 0;
	for(int i = 0; i < 3; i++){
		STATE val = dataright.XCenter[i] - dataleft.XCenter[i];
		dist += val*val;
	}
	dist = sqrt(dist);
	
	//pressao e saturacao
	const BFadREAL poL(dataleft.sol[0][0],0);
	const BFadREAL SoL(dataleft.sol[0][1],1);
	const BFadREAL poR(dataright.sol[0][0],2);
	const BFadREAL SoR(dataright.sol[0][1],3);
	
	BFadREAL pcL,pcR;
	this->CapilarPressure(SoL,pcL);
	this->CapilarPressure(SoR,pcR);
	const BFadREAL pwL = poL - pcL;
	const BFadREAL pwR = poR - pcR;
	const BFadREAL SwL = ((BFadREAL)1.)-SoL;
	const BFadREAL SwR = ((BFadREAL)1.)-SoR;
	
	//Permeabilidade
	TPZFNMatrix<9> K(3,3,0.);
	this->K(K);
	
	const REAL knormal = K(0,0)*data.normal[0]*data.normal[0] + 2.*K(0,1)*data.normal[0]*data.normal[1] +
	K(1,1)*data.normal[1]*data.normal[1] + data.normal[2]*(2.*K(0,2)*data.normal[0] + 
														   2.*K(1,2)*data.normal[1] + K(2,2)*data.normal[2]);
	
	const REAL kgradZn = -K(0,2)*data.normal[0] - K(1,2)*data.normal[1] - K(2,2)*data.normal[2];
	
	// ************* Equacao 1 ******************* /
	
	//Oleo
	BFadREAL BoL,BoR;
	this->Bo(poL, BoL);
	this->Bo(poR, BoR);
	const BFadREAL GammaOleoLeft  = this->g() * this->OilRhoSC()/BoL.val();
	const BFadREAL GammaOleoRight = this->g() * this->OilRhoSC()/BoR.val();
	
	//velocidade de Darcy
	BFadREAL velocOleo = -1.*((knormal*poR-knormal*poL)/dist - (GammaOleoRight*kgradZn+GammaOleoLeft*kgradZn)/2.);
	
	//Mobilidades
	BFadREAL KroL,KroR;
	this->Kro(SoL,KroL);
	this->Kro(SoR,KroR);
	BFadREAL OilViscLeft, OilViscRight;
	this->OilVisc(poL, OilViscLeft);
	this->OilVisc(poR, OilViscRight);
	BFadREAL LambdaOleoLeft = KroL/(OilViscLeft*BoL);
	BFadREAL LambdaOleoRight = KroR/(OilViscRight*BoR);
	
	//Fluxo numerico da primeira equacao do residuo
	BFadREAL Fn1 = 0.;
	if(velocOleo.val() > 0.){
		Fn1 = -LambdaOleoLeft*velocOleo;
	}
	else{
		Fn1 = -LambdaOleoRight*velocOleo;
	}
	
	ef(0,0) += -1.*weight*( -Fn1.val() );
	ef(2,0) += -1.*weight*( +Fn1.val() );
	
#ifndef EXPLICITO
	//ek = -T (R)
	ek(0,0) += -1.*weight*( +Fn1.dx(0) );
	ek(0,1) += -1.*weight*( +Fn1.dx(1) );
	ek(0,2) += -1.*weight*( +Fn1.dx(2) );
	ek(0,3) += -1.*weight*( +Fn1.dx(3) );
	
	ek(2,0) += -1.*weight*( -Fn1.dx(0) );
	ek(2,1) += -1.*weight*( -Fn1.dx(1) );
	ek(2,2) += -1.*weight*( -Fn1.dx(2) );
	ek(2,3) += -1.*weight*( -Fn1.dx(3) );
#endif
	
	// ************* Equacao 2 ******************* /
	
	//Agua
	const REAL Bw = this->Bw();
	const REAL GammaAgua = this->g() * this->WaterRhoSC() / Bw;
	
	//velocidade de Darcy
	BFadREAL velocAgua = -1.*( (knormal*pwR-knormal*pwL)/dist - (GammaAgua*kgradZn+GammaAgua*kgradZn)/2. );
	
	//Mobilidades
	BFadREAL KrwL,KrwR;
	this->Krw(SoL,KrwL);
	this->Krw(SoR,KrwR);
	REAL WaterVisc = this->WaterVisc();
	BFadREAL LambdaAguaLeft = KrwL.val()/(WaterVisc*Bw);
	BFadREAL LambdaAguaRight = KrwR.val()/(WaterVisc*Bw);
	
	//Fluxo numerico da segunda equacao do residuo
	BFadREAL Fn2 = 0.;
	if(velocAgua.val() > 0.){
		Fn2 = -LambdaAguaLeft*velocAgua;
	}
	else{
		Fn2 = -LambdaAguaRight*velocAgua;
	}
	
	ef(1,0) += -1.*weight*( -Fn2.val() );
	ef(3,0) += -1.*weight*( +Fn2.val() );
	
#ifndef EXPLICITO
	//ek = -T (R)
	ek(1,0) += -1.*weight*( +Fn2.dx(0) );
	ek(1,1) += -1.*weight*( +Fn2.dx(1) );
	ek(1,2) += -1.*weight*( +Fn2.dx(2) );
	ek(1,3) += -1.*weight*( +Fn2.dx(3) ); 
	
	ek(3,0) += -1.*weight*( -Fn2.dx(0) );
	ek(3,1) += -1.*weight*( -Fn2.dx(1) );
	ek(3,2) += -1.*weight*( -Fn2.dx(2) );
	ek(3,3) += -1.*weight*( -Fn2.dx(3) ); 
#endif
	
}//method

void TPZBlackOil2P3D::ContributeBCInterface(
    const TPZMaterialDataT<STATE> &data,
    const TPZMaterialDataT<STATE> &dataleft,
    REAL weight,
    TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,
    TPZBndCondT<STATE> &bc)
{
	
	if(gState == ELastState) return;
	
	if(bc.Type() == 3){
		STATE vazao = -1. * (- fabs(bc.Val2()[0]) ); //em m3/m2
		ef(1,0) += weight * vazao;
	}//impondo vazao da injecao de agua
	
	if(bc.Type() == 2) return;//Parede ou simetria
    TPZMaterialDataT<STATE> dataright;
	
	if((bc.Type() == 0) || (bc.Type() == 1)){
		
		//pressao e saturacao
		dataright.sol[0][0] = bc.Val2()[0];
		dataright.sol[0][1] = bc.Val2()[1];

		TPZFNMatrix<16,STATE> auxek(4,4,0.), auxef(4,1,0.);
		this->ContributeInterface(data,dataleft,dataright,weight,auxek,auxef);
		
		for(int i = 0; i < 2; i++){
			ef(i,0) += auxef(i,0);
			for(int j = 0; j < 2; j++){
				ek(i,j) += auxek(i,j);
			}
		}
		
	}//Dirichlet na pressao e Dirichlet ou outflow na saturacao

}//method


enum ESolutionVars { ENone = 0, EWaterPressure = 1, EOilPressure, EWaterSaturation, EOilSaturation, EDarcyVelocity };

int TPZBlackOil2P3D::VariableIndex(const std::string &name) const{
	if(!strcmp("WaterPressure",name.c_str()))   return  EWaterPressure;
	if(!strcmp("OilPressure",name.c_str()))     return  EOilPressure;
	if(!strcmp("WaterSaturation",name.c_str())) return  EWaterSaturation;
	if(!strcmp("OilSaturation",name.c_str()))   return  EOilSaturation;
	if(!strcmp("DarcyVelocity",name.c_str()))   return  EDarcyVelocity;
	return TPZMaterial::VariableIndex(name);
}

int TPZBlackOil2P3D::NSolutionVariables(int var) const{
	if(var == EWaterPressure) return 1;
	if(var == EOilPressure) return 1;
	if(var == EWaterSaturation) return 1;
	if(var == EOilSaturation) return 1;
	if(var == EDarcyVelocity) return 3;
	return TPZMaterial::NSolutionVariables(var);
}

void TPZBlackOil2P3D::Solution(const TPZMaterialDataT<STATE>&data,
                               int var, TPZVec<STATE> &solout){
    auto &Sol = data.sol[0];
	const STATE po = Sol[0];
	const STATE So = Sol[1];
	
	if(var == EWaterPressure){
		BFadREAL pc;
		this->CapilarPressure(So, pc);
		solout[0] = (po-pc.val());
		return;
	}//pw
	
	if(var == EOilPressure){
		solout[0] = po;
		return;
	}//po
	
	if(var == EWaterSaturation){
		solout[0] = 1.-So;
		return;
	}//Sw
	
	if(var == EOilSaturation){
		solout[0] = So;
		return;
	}//So
	
	if(var == EDarcyVelocity){
		cout << "\nERROR AT " << __PRETTY_FUNCTION__ << " \n";
		solout.Fill(0.);
		return;
	}//Velocity
}//method

void TPZBlackOil2P3D::SolutionInterface(const TPZMaterialDataT<STATE> &data,
                                        const TPZMaterialDataT<STATE> &dataleft,
                                        const TPZMaterialDataT<STATE> &dataright,
                                        const int var,
                                        TPZVec<STATE> &Solout){
    //nothing to do here
}

void TPZBlackOil2P3D::GetSolDimensions(uint64_t &u_len,
                                        uint64_t &du_row,
                                        uint64_t &du_col) const
{
    u_len = 1;
    du_row = 1;
    du_col = 3;
}
