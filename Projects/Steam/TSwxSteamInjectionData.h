//---------------------------------------------------------------------------

#ifndef TSwxSteamInjectionDataH
#define TSwxSteamInjectionDataH

#include "ThermalMethodsTables.h"
#include "ConvertionMethods.h"
#include "PropertiesTable.h"

#include <vector>
#include <fstream>
#include <math.h>

class TSwxSteamData;

class TSwxReservoirData {
	/**
	 * frw: raio do poco [m]   fre: raio de drenagem [m]   fHr: espessura do reservatorio [m]   fPe: pressao estatica [Pa]   fKh: permeabilidade horizontal [m2]
	 * fPoros: porosidade (0..1)   fSo: saturacao de Ûleo (0..1)   fTo: temperatura inicial (grau Celsius)
	 */
	double frw, fre, fHr, fPe, fKh, fPoros, fSo, fTo;

	/** Rocha da formacao:  fRockDensity: Massa especifica [kg/m3]   fRockSpecificHeat: Calor especifico [kJ/(kg C)]
	 */
	double fRockDensity, fRockSpecificHeat;

	/** Oleo:   fOilDensity: Massa especifica [kg/m3]   fOilSpecificHeat: Calor especifico [kJ/(kg C)]
	 */
	double fOilDensity, fOilSpecificHeat;

protected:
	std::vector<std::vector<double> > fTable;

public:
	// Construtor a partir das propriedades físicas de uma rocha da formação: porosidade, densidade e calor específico nas unidades acima. O h é a espessura do reservatório
	TSwxReservoirData(double h,double porosity,double rho_rocha,double c_rocha);
	// Construtor a partir da espessura do reservatório, um objeto rocha de formação e a temperatura inicial do reservatório
	TSwxReservoirData(double h,PhysicalProperties *rock, double temperature);

	// Métodos para determinar as propriedades do reservatorio, considerando a rocha e as saturacoes dos fluidos em alguma posicao radial
	
	// Método que retorna a densidade e calor específico do reservatório (rocha e fluídos) na região aquecida
	// As saturações da água, vapor e óleo é dado como média na região aquecida (de vapor)
	// A saturação do vapor é dependente da saturação da água e óleo: Sv = 1 - (Sw + So)
	double Rho_C_Estrela(WaterDataInStateOfSaturation &table,OilData &oildata,double temp,double Sw,double So);   // Unidades: [kg/m3][kJ/(kg K)] --> [kJ/(m3 K)]

	void setReservoirTemperature(double temp) { fTo = temp; }   // Unidades: C (graus Celsius),  K (graus Kelvin)

	// Métodos para recuperar dados do reservatório
	double getReservoirTemperature() { return fTo; }
	double getHReservoir() { return fHr; }

	// Método que devolve a permeabilidade relativa da rocha de formação como função da Saturação de água. Tabulados.
	// Pela existência de uma parcela irredutivel de óleo e de água, a saturação do vapor e da água não pode ser 0.0 ou 1.0
	double getRelativePermeability_Vapour(double Sw);
	double getRelativePermeability_Water(double Sw);

	// saida no Mathematica
	const void Print(std::ostream &out) {
		out << "\n/** TSwxReservoirData **/\n";
		out << "frw = " << frw << "\n" << "fre = " << fre << "\n" << "fHr = " << fHr << "\n" << "fPe = " << fPe << "\n" << "fKh = " << fKh << "\n"
				<< "fPoros = " << fPoros << "\n" << "fSo = " << fSo << "\n" << "fTo = " << fTo << "\n" << "Rock: Density = " << fRockDensity << "\t"
				<< "SpecificHeat = " << fRockSpecificHeat << "\n" << "Oil: Density = " << fOilDensity << "\t" << "SpecificHeat = " << fOilSpecificHeat << "\n";
	}
};

struct TSwxOverburdenData {
	/**
	 * fThermalConductivity: Condutividade tÈrmica [J/(s m C)]
	 * fDensity: Massa especÌfica [kg/m3]
	 * fSpecificHeat: Calor especÌfico [J/(kg C)]
	 */
	double fThermalConductivity, fDensity, fSpecificHeat;
public:
	TSwxOverburdenData(double condtermica, double densidade,double calorespecif);
	TSwxOverburdenData(PhysicalProperties *rock);
	TSwxOverburdenData(TSwxOverburdenData &cpy);
	
	// Método para recuperar o produto das propriedades físicas Condutividade_Termica * Densidade * Calor_Especifico
	double getProductOfTheProperties();
	// double getThermalCondutivity() { return fThermalConductivity; }  // Os dados são públicos, pois estamos numa estrutura

	const void Print(std::ostream &out) {
		out << "\n/** TSwxOverburdenData **/\n" << "fThermalConductivity = " << fThermalConductivity << "\n" << "fDensity = " << fDensity << "\n" << "fSpecificHeat = " << fSpecificHeat << "\n";
	}
};

struct TSwxInjectionData {
	public:
	/**
	 * fTable: tabela de injecao. Tempo [s] x Vacao [m3/2]
	 */
	std::vector< std::pair<double,double> > fTable;

	void setConstantRateAndTime(double rate,double time);
	
	const void Print(std::ostream &out) {
		out << "\n/** TSwxInjectionData **/\n" << "fTable.size() = " << fTable.size() << "\n";
		for(unsigned int i = 0; i < fTable.size(); i++) {
			out << fTable[i].first << "\t" << fTable[i].second << "\n";
		}
	}
};

/** Dados de entrada */
class TSwxSteamInjectionInputData {
	
	TSwxReservoirData fReservoirData;
	TSwxOverburdenData fOverburdenData;

	/**
	 * Objeto que armazena os métodos para o calculo das propriedades do vapor e a tabela dos valores dessas propriedades em estado de saturacao
	 */
	WaterDataInStateOfSaturation fWaterInSaturationState;
	/** Objeto com as propriedades físicas para o óleo crudo */
	OilData fOilData;
	
	double fQuality;
	double fSteamTemperature;

	TSwxInjectionData fInjectionData;

public:
	
	TSwxSteamInjectionInputData(TSwxReservoirData& reservoir,TSwxOverburdenData& overrock,double temp,double quality);

	void setInjectionData(double rate,double time);
	
	// Baseado na qualidade do vapor, retorna a saturação média esperada da água
	double getSaturationWaterIntoReservoir();   // a função considera que a saturação da água não pode ser menor da saturação irredutível de água em um reservatório
	double getSaturationWater();   // Dada a qualidade do vapor calcula a saturação da água sem se importar se o valor é menor a saturação irredutível
	
	void setInjectionTemperature(double tempv) { fSteamTemperature = tempv; }
	void setSteamQuality(double quality) { fQuality = quality; }
	
	double getRho_C_Estrela(int i,double So);  // Sw, Sv, So sao as saturacoes de agua, vapor e oleo em uma celula. De 0 ate 1.
	
	const void Print(std::ostream &out) {
		fReservoirData.Print(out);
		fOverburdenData.Print(out);
		fInjectionData.Print(out);
	}
	TSwxInjectionData &getInjectionData() { return fInjectionData; }
	TSwxReservoirData &getReservoirData() { return fReservoirData; }
	
	double getPropertiesProductFromOverburden() { return fOverburdenData.getProductOfTheProperties(); }

	double getInjectionTemperature() { return fSteamTemperature; } 
	
	double getHReservoir() { return fReservoirData.getHReservoir(); }
};

/** Resultados */
struct TSwxSteamInjectionOutputData {
	public:
	/**
	 * fSteamFrontPos: PosiÁ„o da frente vapor [m] x Tempo [s]
	 */
	std::vector< std::pair<double,double> > fSteamFrontPos;
	/**
	 * fInjectionPressure: Pres„o de injeÁ„o [Pa] x Tempo [s]
	 */
	std::vector< std::pair<double,double> > fInjectionPressure;
	/**
	 * fFracturePressure: Pres„o de quebra [Pa] x Tempo [s]
	 */
	std::vector< std::pair<double,double> > fFracturePressure;

	const void Print(std::ostream &out) {
		out << "\n/** TSwxSteamInjectionOutputData **/\n" << "fSteamFrontPos.size() = " << fSteamFrontPos.size() << "\n";
		for(unsigned int i = 0; i < fSteamFrontPos.size(); i++) {
			out << fSteamFrontPos[i].first << "\t" << fSteamFrontPos[i].second << "\n";
		}
		out << "fInjectionPressure.size() = " << fInjectionPressure.size() << "\n";
		for(unsigned int i = 0; i < fInjectionPressure.size(); i++) {
			out << fInjectionPressure[i].first << "\t" << fInjectionPressure[i].second << "\n";
		}
		out << "fFracturePressure.size() = " << fFracturePressure.size() << "\n";
		for(unsigned int i = 0; i < fFracturePressure.size(); i++) {
			out << fFracturePressure[i].first << "\t" << fFracturePressure[i].second << "\n";
		}
	}

	/** 
	 * MÈtodo para desenvolvimento da interface gr·fica enquanto n„o È integrada ao kernel
	 */
	void SampleValues(double L, int np) {
		fSteamFrontPos.resize(np);
		fInjectionPressure.resize(np);
		fFracturePressure.resize(np);
		double d = L/((double)(np)-1.);
		double x, front, ipress, fpress;
		for(int i = 0; i < np; i++) {
			x = d*(double)(i);
			front = sin(x);
			ipress = cos(x);
			fpress = sin(x)*x;
			fSteamFrontPos[i].first = x; fSteamFrontPos[i].second = front;
			fInjectionPressure[i].first = x; fInjectionPressure[i].second = ipress;
			fFracturePressure[i].first = x; fFracturePressure[i].second = fpress;
		}///for
	}///void
};

/** Dados do problema: entrada e resultados
*/
class TSwxSteamInjectionData {
	/** Dados de entrada */
	TSwxSteamInjectionInputData fInput;
	/** Resultados */
	TSwxSteamInjectionOutputData fOutput;
	
public:
	TSwxSteamInjectionData(TSwxReservoirData& reservoir,TSwxOverburdenData& rock_over,double temp,double quality);
	
	double getInjectionTemperature();
	double getRegionOfSteamArea(double tempo);
	double getRegionAuxiliar(int i,double tempo);   // Fornece a Area "cursiva" auxiliar para achar a área de vapor para taxas de injeção constantes

	double getSteamFrontPosStored(int i);
	
	TSwxReservoirData &getReservoirData();

	void setInjectionData(double rate,double tempo);

	// Preenche um valor do raio para a posição da frente do vapor em um determinado tempo
	void setFrontPosition(double radio,double time);

	void PrintRadiosToMathematicaFile(std::ostream &out);
	void PrintAreasToMathematicaFile(std::ostream &out);

	const void Print(std::ostream &out) {
		fInput.Print(out);
		fOutput.Print(out);
	}
};

//---------------------------------------------------------------------------
#endif