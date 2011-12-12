//---------------------------------------------------------------------------


#pragma hdrstop

#include "TSwxSteamInjectionData.h"

//---------------------------------------------------------------------------

#pragma package(smart_init)

/**
 * Metodos relacionados ao objeto reservatorio
 */
TSwxReservoirData::TSwxReservoirData(double h,double porosity,double rho_rocha,double c_rocha) {
	fPoros = porosity;
	fHr = h;
	fRockDensity = rho_rocha;
	fRockSpecificHeat = c_rocha;
}

TSwxReservoirData::TSwxReservoirData(double h,PhysicalProperties *rock, double temperature) {
	fPoros = rock->porosity;
	fHr = h;
	fRockDensity = rock->density;
	fRockSpecificHeat = rock->specificheat;
	fTo = temperature;
	// Preenchendo a tabela de valores para a permeabilidade relativa
	// Valores relacionados da temperatura e da pressão para o estado de saturação do fluido
	// Temperatura -> indice 0; Pressao -> indice 1; Calor Latente -> indice 7
	fTable.resize(17);   // Temos 17 valores tabulados
	// Observar: A saturação de água não pode ser ZERO, pois existe a Saturação irredutível de água, e também não pode ser UM pois existe a saturação mínima de óleo
	double WaterSaturation[17] = {0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9 };
	double RelativePermeability_Vapour[17] = {0.85, 0.71, 0.6, 0.51, 0.43, 0.36, 0.3, 0.25, 0.215, 0.187, 0.166, 0.137, 0.11, 0.1, 0.065, 0.03, 0.001};
	double RelativePermeability_Water[17] = {0.008, 0.011, 0.015, 0.023, 0.035, 0.048, 0.065, 0.083, 0.105, 0.127, 0.152, 0.185, 0.23, 0.28, 0.39, 0.56, 0.91};
    std::vector<std::vector<double> >::iterator it;
	int i;
	for(i=0,it=fTable.begin();it<fTable.end();i++,it++) {
        std::vector<double> vec;
		vec.push_back(WaterSaturation[i]);
		vec.push_back(RelativePermeability_Vapour[i]);
		vec.push_back(RelativePermeability_Water[i]);
		// Copiando o vetor preenchido na linha correspondente
		*it = vec;
	}
}

double TSwxReservoirData::Rho_C_Estrela(WaterDataInStateOfSaturation &watertable,OilData &oildata,double temp,double Sw,double So) {
	double Rho_C_Fluidos = So*oildata.getDensityToOil(temp)*oildata.getSpecificHeatToOil(temp);
	Rho_C_Fluidos += (Sw*watertable.getSaturationStateDensityToLiquidWater(temp)*watertable.getSaturationStateSpecificHeatToLiquidWater(temp));
	Rho_C_Fluidos += ((1-(Sw+So))*watertable.getSaturationStateDensityToSteam(temp)*watertable.getSaturationStateSpecificHeatToSteam(temp));
	return (((1-fPoros)*fRockDensity*fRockSpecificHeat)+fPoros*Rho_C_Fluidos);
}

double TSwxReservoirData::getRelativePermeability_Vapour(double Sw) {
	if(Sw < 0.1 || Sw > 0.9)
		return 0.0;
	return getInterpolatedValue(fTable,Sw,1);
}
double TSwxReservoirData::getRelativePermeability_Water(double Sw) {
	if(Sw < 0.1 || Sw > 0.9)
		return 0.0;
	return getInterpolatedValue(fTable,Sw,2);
}

void TSwxInjectionData::setConstantRateAndTime(double rate,double time) {
	std::pair<double,double> par(rate,time);
	fTable.push_back(par);
}

TSwxOverburdenData::TSwxOverburdenData(TSwxOverburdenData &cpy) {
	fThermalConductivity = cpy.fThermalConductivity;
	fDensity = cpy.fDensity;
	fSpecificHeat = cpy.fSpecificHeat;
}
TSwxOverburdenData::TSwxOverburdenData(double condtermica, double densidade,double calorespecif) {
	fThermalConductivity = condtermica;
	fDensity = densidade;
	fSpecificHeat = calorespecif;
}
TSwxOverburdenData::TSwxOverburdenData(PhysicalProperties *rock) {
	fThermalConductivity = ConvertWTokJperH(rock->thermalconductivity);
	fDensity = rock->density;
	fSpecificHeat = rock->specificheat;
}

double TSwxOverburdenData::getProductOfTheProperties() {
	return (fThermalConductivity*fDensity*fSpecificHeat);
}
void TSwxSteamInjectionInputData::setInjectionData(double rate,double time) {
	if(!fInjectionData.fTable.size())
		fInjectionData.setConstantRateAndTime(0.0,0.0);
	fInjectionData.setConstantRateAndTime(rate,time);
}

double TSwxSteamInjectionInputData::getSaturationWaterIntoReservoir() {
	// A saturação irredutível de água (líquida) em um reservatório é de Swirr = 0.23
	double Sw = 0.;
	// Computo a partir da qualidade do vapor
	
	if(Sw < 0.23)
		Sw = 0.23;
	if(Sw > 1.)
		return 1.;
	return Sw;
}

double TSwxSteamInjectionInputData::getSaturationWater() {
	// A saturação irredutível de água (líquida) em um reservatório é de Swirr = 0.23 
	// No presente caso não interessa se a saturação deu um valor menor a saturação irredutível
	double Sw = 0.;
	// Computo a partir da qualidade do vapor
	if(Sw < 0.)
		return 0.;
	if(Sw > 1.)
		return 1.;
	return Sw;
}
double TSwxSteamInjectionInputData::getRho_C_Estrela(int i,double So) {
	if(So < 0. || So > 1.)
		return -1.;
	// Determinando o valor da taxa de injeção da massa de água a partir da taxa de injeção de calor que é assumido constante por partes
//	double massrate = fInjectionData.fTable[i].first/(fWaterInSaturationState.getSaturationStateSpecificEnthalpyToLiquidWater(fSteamTemperature)+fQuality*fWaterInSaturationState.getSaturationStateLatentHeat(fSteamTemperature)-fWaterInSaturationState.getSpecificEnthalpyToLiquidWater(fReservoirData.getReservoirTemperature()));
	double Sw = fWaterInSaturationState.getSaturationWater(i,fQuality,So);
	return fReservoirData.Rho_C_Estrela(fWaterInSaturationState,fOilData,fSteamTemperature,Sw,So);
}

// Assume-se que a temperatura esta sendo introduzida em graus Celsius (C), entao sera convertida em graus Kelvin (K)
TSwxSteamInjectionInputData::TSwxSteamInjectionInputData(TSwxReservoirData& reservoir,TSwxOverburdenData& overrock,double temp,double quality) : fReservoirData(reservoir), fOilData(), fWaterInSaturationState(), fOverburdenData(overrock) {
	fSteamTemperature = temp;
	fQuality = quality;
}

TSwxSteamInjectionData::TSwxSteamInjectionData(TSwxReservoirData &reservoir,TSwxOverburdenData &rock_over,double temp,double quality) : fInput(reservoir,rock_over,temp,quality) {
}

double TSwxSteamInjectionData::getInjectionTemperature() {
	return fInput.getInjectionTemperature();
}

TSwxReservoirData &TSwxSteamInjectionData::getReservoirData() {
	return fInput.getReservoirData();
}

void TSwxSteamInjectionData::setInjectionData(double rate,double tempo) {
	fInput.setInjectionData(rate,tempo);
}

double TSwxSteamInjectionData::getRegionAuxiliar(int i,double tempo) {
	if(i < 0 || tempo < 0.0 || IsZero(tempo) || i+2 > fInput.getInjectionData().fTable.size())
		return 0.0;

	double tau = fInput.getInjectionData().fTable[i].second;   // isso significa que tau é o i-esimo tempo -> ti, onde i=0, 1, ..., tamanho_vetor_tempos-2
	if(tempo < tau || IsZero(tempo-tau))
		return 0.0;

	// se tempo > tau, entao calculamos o tempo - tau
	tau = tempo - tau;

	// Calculamos o A_cursiva(tau)
	double difTemp = ConvertCToK(getInjectionTemperature()) - ConvertCToK(fInput.getReservoirData().getReservoirTemperature());   // temperatura em Celsius ou Kelvin, pois a diferença é a mesma
//	double difTemp = getInjectionTemperature() - fInput.getReservoirData().getReservoirTemperature();   // temperatura em Celsius ou Kelvin, pois a diferença é a mesma
	double ResidualOilSaturation = 0.09;   // Valor copiado do livro Burger[...] pag 113
//	double ResidualOilSaturation = 0.0;
	
	double delta = 4*fInput.getPropertiesProductFromOverburden();
	double H_RhoCEstrela = fInput.getRho_C_Estrela(i,ResidualOilSaturation);
	H_RhoCEstrela *= fInput.getHReservoir();
	
	if(IsZero(H_RhoCEstrela) || IsZero(delta) || IsZero(difTemp)) 
		return 0.;

	double var_adim = sqrt(delta*tau)/H_RhoCEstrela;
	double func_adim = exp(var_adim*var_adim) * erfc(var_adim) + ((2./sqrt(M_PI))*var_adim) - 1.;

	return ((H_RhoCEstrela*func_adim)/(delta*difTemp));
}

double TSwxSteamInjectionData::getRegionOfSteamArea(double tempo) {
	if(tempo < 0.0 || IsZero(tempo) || fInput.getInjectionData().fTable.size() < 2)   // O tempo nao pode ser negativo pois ele informa o tempo de injeção a uma determinada taxa, até x horas.
		return 0.0;
	
	double QiMinus1, Qi, ti, area = 0.0;
	int counter;
	// iterators
	std::vector< std::pair<double,double> >::iterator it;
	
	for(it=fInput.getInjectionData().fTable.begin()+1,counter=0;it<fInput.getInjectionData().fTable.end();it++,counter++) {
		ti = (*it).second;
		QiMinus1 = (*(it-1)).first;
		Qi = (*it).first;
		area += ((Qi - QiMinus1)*getRegionAuxiliar(counter,tempo));
	}
	return area;
}

double TSwxSteamInjectionData::getSteamFrontPosStored(int i) {
	return fOutput.fSteamFrontPos[i].first;
}
void TSwxSteamInjectionData::setFrontPosition(double radio,double time) {
	std::pair<double,double> par(radio,time);
	fOutput.fSteamFrontPos.push_back(par);
}
void TSwxSteamInjectionData::PrintRadiosToMathematicaFile(std::ostream &out) {
	// Tabela de raios da área esquentada
	int i, NTempos = fOutput.fSteamFrontPos.size();
	out << "graphraios = {";
	for(i=0;i < NTempos;i++) {
		out << "{" << fOutput.fSteamFrontPos[i].second << "," << getSteamFrontPosStored(i) << "}";
		if(i != NTempos-1) out << ",";
	}
	out << "};" << std::endl;
	out << "gr=ListPlot[graphareas,Joined->True,Frame->True]" << std::endl;
	out << "gr=ListPlot[graphraios,Joined->True,Frame->True]" << std::endl;
}
void TSwxSteamInjectionData::PrintAreasToMathematicaFile(std::ostream &out) {
	// Tabela de áreas das regiões esquentadas para cada passo de tempo
	int i, NTempos = fOutput.fSteamFrontPos.size();
	double raio;
	out << "graphareas = {";
	for(i=0;i < NTempos;i++) {
		raio = fOutput.fSteamFrontPos[i].first;
		out << "{" << fOutput.fSteamFrontPos[i].second << "," << M_PI*raio*raio << "}";
		if(i != NTempos-1) out << ",";
	}
	out << "};" << std::endl;
}