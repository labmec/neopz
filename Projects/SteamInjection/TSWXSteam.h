//---------------------------------------------------------------------------

#ifndef TSWXSteamH
#define TSWXSteamH

#include "TSWXInputData.h"
#include "TPZGuiInterface.h"

//---------------------------------------------------------------------------

class TSwxSteam
{
	public:

	TSwxSteam();
	void SetInputData(
			double rw,//raio do poco
			double Hr,//espessura da camada permeavel
			double Poros,//porosidade da camada permeavel
			double Tr,//temperatura da camada permeavel
			double RockDensity,//densidade da camada permeavel
			double RockSpecificHeat,//calor especifico da camada permeavel
			double re,//raio de drenagem
			double Pe,//pressao no farfield
			double Kh,//permeabilidade horizontal da camada permeavel
			double So,//saturacao inicial de oleo
			double OilDensity,//densidade do oleo
			double OilSpecificHeat,//calor especifico do oleo
			double OilVisc,//viscosidade do oleo
			double confinCondut,//condutividade termica da camada condinante
			double confinSpecificMass,//densidade da camada condinante
			double confinSpecificHeat,//calor especifico da camada condinante
			double steamTemperature,//temperatura do vapor
			double steamQuality,//titularidade do vapor
			std::map< double , std::pair< double, double > > & table,//tabela: col[0] = tempo, col[1] = taxa injecao, col[2] = fluxo de massa
			double E,//modulo de young da camada permeavel
			double nu,//poisson da camada permeavel
			double timeStep//passo de tempo (discretizacao da tabela "table")
	);

	void WriteMe(std::ofstream &outfileSI);
	void ReadMe(std::ifstream &infileSI);
	void PrintToMathematicaFile(std::map< double , std::pair<double, double> > &TempoRaioSigmaTheta, std::ostream &out);

	void getRadiusAndMaxSigmaThetaForTheseTimes(const std::vector<double> &time,
																							std::map< double , std::pair<double, double> > &Time_Radius_MaxSigmaTheta,
																							TPZGuiInterface * progressInfo);

	void getRadiusAndMaxSigmaThetaForTableTimes(std::map< double , std::pair<double, double> > &Time_Radius_MaxSigmaTheta,
																							TPZGuiInterface * progressInfo);

	TSWXInputData &getInputData()
	{
		return fInput;
	}

	//Phils equation
	double ComputeSteamPressure(double time);

	/////////////////////////////////////////////////////////////////////////////
	private:

	TSWXInputData fInput;

	double getRadiusOfSteamFront(double tempo);
	double getMassRateOfSteam(const double tempo);

	double getRegionOfSteamArea(double tempo);
	double getRegionAuxiliar(int i, double tempo);

	double m_timeStep;
};
#endif
