//---------------------------------------------------------------------------


#pragma hdrstop

#include "TSWXSteam.h"

#ifdef USING_BOOST
	#include <boost/math/special_functions/erf.hpp>
#else
	#include <math.h>
#endif
//---------------------------------------------------------------------------

#pragma package(smart_init)

TSwxSteam::TSwxSteam() : fInput()
{
//	std::ifstream readF("refpatternsVapor.txt");
//	gRefDBase.ReadRefPatternDBase(readF); // Nao fez muita diferenca refinar!
}
//---------------------------------------------------------------------------

void TSwxSteam::SetInputData(
			double rw,
			double Hr,
			double Poros,
			double Tr,
			double RockDensity,
			double RockSpecificHeat,
			double re,
			double Pe,
			double Kh,
			double So,
			double OilDensity,
			double OilSpecificHeat,
			double OilVisc,
			double confinCondut,
			double confinSpecificMass,
			double confinSpecificHeat,
			double steamTemperature,
			double steamQuality,
			std::map< double , std::pair< double, double > > & table,
			double E,
			double nu,
			double timeStep
) {
		fInput.getReservoirData().SetData(rw, re, Hr, Pe, Kh, Poros, So, Tr, RockDensity, RockSpecificHeat, OilDensity, OilSpecificHeat, OilVisc);
		fInput.getConfinementData().SetData(confinCondut, confinSpecificMass, confinSpecificHeat);
		fInput.getInjectionData().setInjectionData(table, steamQuality, steamTemperature);
		fInput.getMaxSigmaThetaData()->SetData(rw, Hr, E, nu);
		this->m_timeStep = timeStep;
		if(fabs(m_timeStep) < 1.E-12)
		{
			std::cout << "passo de tempo negativo na classe TSWXSteam!";
			DebugStop();
		}
}
//---------------------------------------------------------------------------

void TSwxSteam::WriteMe(std::ofstream &outfileSI)
{
	int size;
	double rw = fInput.getReservoirData().getWellRadius();
	double Hr = fInput.getReservoirData().getHReservoir();
	double Poros = fInput.getReservoirData().getPorosity();
	double Tr = fInput.getReservoirData().Temperature();
	double RockDensity = fInput.getReservoirData().getRockDensity();
	double RockSpecificHeat = fInput.getReservoirData().getRockSpecificHeat();
	double re = fInput.getReservoirData().getDrainageRadius();
	double Pe = fInput.getReservoirData().getStaticPressure();
	double Kh = fInput.getReservoirData().getKh();
	double So = fInput.getReservoirData().OilSaturation();
	double OilDensity = fInput.getReservoirData().OilDensity();
	double OilSpecificHeat = fInput.getReservoirData().OilSpecificHeat();
	double OilVisc = fInput.getReservoirData().OilVisc();

	double confinCondut = 1000.*fInput.getConfinementData().ThermalConductivity();
	double confinSpecificMass = fInput.getConfinementData().Density();
	double confinSpecificHeat = fInput.getConfinementData().SpecificHeat();

	double E = fInput.getMaxSigmaThetaData()->getE();
	double nu = fInput.getMaxSigmaThetaData()->getNu();

	double steamTemperature = fInput.getInjectionData().Temperature();
	double steamQuality = fInput.getInjectionData().getInjectionQuality();
	std::vector< std::pair<double, double> > table = fInput.getInjectionData().Table();
	std::map<double,double> MassRate = fInput.getInjectionData().MassRate();

	double timeStep = this->m_timeStep;

	outfileSI << "ReservoirData:\n";
	outfileSI << "rw " << rw << "\n";
	outfileSI << "Hr "  << Hr << "\n";
	outfileSI << "Poros " << Poros << "\n";
	outfileSI << "Tr " << Tr << "\n";
	outfileSI << "RockDensity " << RockDensity << "\n";
	outfileSI << "RockSpecificHeat " << RockSpecificHeat << "\n";
	outfileSI << "re " << re << "\n";
	outfileSI << "Pe " << Pe << "\n";
	outfileSI << "Kh " << Kh << "\n";
	outfileSI << "So " << So << "\n";
	outfileSI << "OilDensity " << OilDensity << "\n";
	outfileSI << "OilSpecificHeat " << OilSpecificHeat << "\n";
	outfileSI << "OilVisc " << OilVisc << "\n";
	outfileSI << "confinCondut " << confinCondut << "\n";
	outfileSI << "confinSpecificMass " << confinSpecificMass << "\n";
	outfileSI << "confinSpecificHeat " << confinSpecificHeat << "\n";
	outfileSI << "E " << E << "\n";
	outfileSI << "nu " << nu << "\n";
	outfileSI << "steamTemperature " << steamTemperature << "\n";
	outfileSI << "steamQuality " << steamQuality << "\n\n";

	outfileSI << "ConfinementData:\n";
	outfileSI << "table:\n";
	size = table.size();
	outfileSI << size << "\n";
	for(int p = 0; p < size; p++)
	{
			outfileSI << table[p].first << " " << table[p].second << "\n";
	}
	outfileSI << "\n";

	outfileSI << "massRate:\n";
	size = MassRate.size();
	outfileSI << size << "\n";
	std::map<double,double>::iterator it;
	for(it = MassRate.begin(); it != MassRate.end(); it++)
	{
			outfileSI << it->first << " " << it->second << "\n";
	}
	outfileSI << "\n";

	outfileSI << "ProcedureTimeStep:\n";
	outfileSI << "timeStep " << timeStep << "\n";
	outfileSI.close();
}
//---------------------------------------------------------------------------

void TSwxSteam::ReadMe(std::ifstream &infileSI)
{
	int size;
	string ClassData, actualData;

	double rw;
	double Hr;
	double Poros;
	double Tr;
	double RockDensity;
	double RockSpecificHeat;
	double re;
	double Pe;
	double Kh;
	double So;
	double OilDensity;
	double OilSpecificHeat;
	double OilVisc;

	double confinCondut;
	double confinSpecificMass;
	double confinSpecificHeat;

	double E;
	double nu;

	double steamTemperature;
	double steamQuality;
	std::vector< std::pair<double, double> > table;
	std::map<double,double> MassRate;
	std::map< double , std::pair< double, double > > inputTable;

	double timeStep;

	infileSI >> ClassData;
	infileSI >> actualData >> rw;
	infileSI >> actualData  >> Hr;
	infileSI >> actualData >> Poros;
	infileSI >> actualData >> Tr;
	infileSI >> actualData >> RockDensity;
	infileSI >> actualData >> RockSpecificHeat;
	infileSI >> actualData >> re;
	infileSI >> actualData >> Pe;
	infileSI >> actualData >> Kh;
	infileSI >> actualData >> So;
	infileSI >> actualData >> OilDensity;
	infileSI >> actualData >> OilSpecificHeat;
	infileSI >> actualData >> OilVisc;
	infileSI >> actualData >> confinCondut;
	infileSI >> actualData >> confinSpecificMass;
	infileSI >> actualData >> confinSpecificHeat;
	infileSI >> actualData >> E;
	infileSI >> actualData >> nu;
	infileSI >> actualData >> steamTemperature;
	infileSI >> actualData >> steamQuality;

	infileSI >> ClassData;
	infileSI >> actualData;
	infileSI >> size; table.resize(size);
	for(int p = 0; p < size; p++)
	{
			double first, second;
			infileSI >> first >> second;
			table[p] = std::make_pair(first,second);
	}

	infileSI >> actualData;
	infileSI >> size;
	for(int p = 0; p < size; p++)
	{
			double first, second;
			infileSI >> first >> second;
			MassRate[first] = second;
	}

	infileSI >> ClassData;
	infileSI >> actualData >> timeStep;

	{
		fInput.getReservoirData().SetData(rw, re, Hr, Pe, Kh, Poros, So, Tr, RockDensity, RockSpecificHeat, OilDensity, OilSpecificHeat, OilVisc);
		fInput.getConfinementData().SetData(confinCondut, confinSpecificMass, confinSpecificHeat);
		fInput.getInjectionData().setInjectionData(table, MassRate, steamQuality, steamTemperature);
		fInput.getMaxSigmaThetaData()->SetData(rw, Hr, E, nu);
		this->m_timeStep = timeStep;
	}
}
//---------------------------------------------------------------------------

void TSwxSteam::PrintToMathematicaFile(std::map< double , std::pair<double, double> > &TempoRaioSigmaTheta, std::ostream &out) {

	std::map< double , std::pair<double, double> >::iterator it;

	int i = 0, NTempos = TempoRaioSigmaTheta.size();

	//---------------------- Posicao da frente de vapor em cada instante

	out << "(*{tempo,raio}*)" << endl;
	out << "graphRaios = {";
	for(it = TempoRaioSigmaTheta.begin(); it != TempoRaioSigmaTheta.end(); it++, i++)
	{
		double instante = it->first/3600.; //em horas
		double posicao = it->second.first; //em metros
		out << "{" << instante << "," << posicao << "}";
		if(i != NTempos-1) out << ",";
	}
	out << "};" << endl << endl;
	out << "radiusGR=ListLinePlot[graphRaios, Filling -> Axis, AxesLabel -> {h, m}, AxesOrigin -> {0, 0}, PlotLabel -> \"Posicao da Frente de Vapor x Tempo\"]" << endl << endl;

	//---------------------- Sigma theta maximo em cada instante

	out << "(*{tempo,sigmathetaMax}*)" << endl;
	out << "graphSigmaMax = {";
	i = 0;
	for(it = TempoRaioSigmaTheta.begin(); it != TempoRaioSigmaTheta.end(); it++, i++)
	{
		double instante = it->first/3600.; //em horas
		double tensaoMax = it->second.second/(1.E6); //em MPa
		out << "{" << instante << "," << tensaoMax << "}";
		if(i != NTempos-1) out << ",";
	}
	out << "};" << endl << endl;
	out << "stressGR=ListLinePlot[graphSigmaMax, Filling -> Axis, AxesLabel -> {h, MPa}, AxesOrigin -> {0, 0}, PlotStyle -> Red, FillingStyle -> Opacity[0.2, Red], PlotLabel -> \"SigmaThetaMax x Tempo\"]" << endl << endl;

	out.flush();
}
//---------------------------------------------------------------------------

void TSwxSteam::getRadiusAndMaxSigmaThetaForTheseTimes(const std::vector<double> &SItime,
																									std::map< double , std::pair<double, double> > &Time_Radius_MaxSigmaTheta, TPZGuiInterface * progressInfo)
{
	double rw = fInput.getReservoirData().getWellRadius();

	Time_Radius_MaxSigmaTheta.clear();
	int nTimes = SItime.size();
	int validTimes = nTimes;
	std::map<double, double> time_r, time_maxSigma;
//	int pos = 0;
	int ini = 0;

	for(int i = 0; i < nTimes; i++)
	{

			double T = SItime[i];

			double r = getRadiusOfSteamFront(T);

			if(r < 5.*rw)

			{

				validTimes--;

				continue;

			}


			//ProgressBar da GUI

			if(progressInfo)

			{

					if(progressInfo->AmIKilled()) return;

					ini++;
					int barPos = int( 100. * ini/double(validTimes+1) + 0.5);
					progressInfo->ProgressBarPos() = barPos;

					std::stringstream infoTXT; infoTXT << "Progresso: ";
					infoTXT << barPos << "%";
					progressInfo->Message() = infoTXT.str().c_str();
					progressInfo->UpdateCaption();
			}


			time_r[T] = r;

//			double mass = getMassRateOfSteam(T);

			double DistrRightDown = ComputeSteamPressure(T);
			fInput.getMaxSigmaThetaData()->ComputeMaxSigmaTheta(T, r, DistrRightDown);

	}


	if(progressInfo)

	{

		progressInfo->ProgressBarPos() = 0;
		progressInfo->Message() = "";
		progressInfo->UpdateCaption();

	}
	time_maxSigma = fInput.getMaxSigmaThetaData()->GetSolution();


	std::map<double, double>::iterator itR, itS;

	for(itR = time_r.begin(), itS = time_maxSigma.begin(); itR != time_r.end(); itR++, itS++)

	{

			double timeR = itR->first;

			double timeS = itS->first;

			if(fabs(timeR - timeS) > 1.E-5)

			{

					std::cout << "Time dont match on TSwxSteam::getRadiusAndMaxSigmaThetaForTheseTimes(...)!\n";

					DebugStop();

			}

			double radius = itR->second;

			double sigmaMax = itS->second;


			Time_Radius_MaxSigmaTheta[timeR] = std::make_pair(radius,sigmaMax);

	}

}
//---------------------------------------------------------------------------

void TSwxSteam::getRadiusAndMaxSigmaThetaForTableTimes(std::map< double , std::pair<double, double> > &Time_Radius_MaxSigmaTheta, TPZGuiInterface * progressInfo)
{
	std::vector<double> Tabletimes = fInput.getInjectionData().getSITimes();
	std::vector<double> times(0);
	for(int tt = 0; tt < Tabletimes.size(); tt++)
	{
		double tabtime = Tabletimes[tt];
		int nTimes = int(tabtime/m_timeStep);

		for(int t = 1; t < nTimes; t++)
		{
			double stepTime = t*m_timeStep;
			times.push_back(stepTime);
		}
		times.push_back(tabtime);
	}
	getRadiusAndMaxSigmaThetaForTheseTimes(times, Time_Radius_MaxSigmaTheta, progressInfo);
}
//---------------------------------------------------------------------------

double TSwxSteam::getRegionAuxiliar(int i, double tempo)
{
	if (i < 0 || tempo < 0.0 || IsZero(tempo) || i + 2 > fInput.getInjectionData().Table().size())
	{
		return 0.0;
	}

	double tau = fInput.getInjectionData().Table()[i].first; // isso significa que tau é o i-esimo tempo -> ti, onde i=0, 1, ..., tamanho_vetor_tempos-2
	if (tempo < tau || IsZero(tempo - tau))
	{
		return 0.0;
	}

	// se tempo > tau, entao calculamos o tempo - tau
	tau = tempo - tau;

	// Calculamos o A_cursiva(tau)
	double difTemp = ConvertCToK(fInput.getInjectionData().Temperature()) - ConvertCToK(fInput.getReservoirData().Temperature());

	double ResidualOilSaturation = 0.09;//valor fixo (nao me pergunte porque)

	double delta = 4 * fInput.getConfinementData().getProductOfTheProperties();
	double H_RhoCEstrela = fInput.getRho_C_Estrela(i, ResidualOilSaturation);
	H_RhoCEstrela *= fInput.getReservoirData().getHReservoir();

	if (IsZero(H_RhoCEstrela) || IsZero(delta) || IsZero(difTemp))
	{
		return 0.;
	}

	double var_adim = sqrt(delta * tau) / H_RhoCEstrela;
	double func_adim;
#ifdef USING_BOOST
	func_adim = exp(var_adim * var_adim) * boost::math::erfc(var_adim) + ((2. / sqrt(M_PI)) * var_adim) - 1.;
#else
	func_adim = exp(var_adim * var_adim) * erfc(var_adim) + ((2. / sqrt(M_PI)) * var_adim) - 1.;
#endif
	
	return((H_RhoCEstrela * func_adim) / (delta * difTemp));
}
//---------------------------------------------------------------------------

//Phils equation
double TSwxSteam::ComputeSteamPressure(double time)
{
	double re = fInput.getReservoirData().getDrainageRadius();
	double pe = fInput.getReservoirData().getStaticPressure();
	double K = fInput.getReservoirData().getKh();
	double H = fInput.getReservoirData().getHReservoir();

	double specificMassOil = fInput.getReservoirData().OilDensity();
	double satOil = fInput.getReservoirData().OilSaturation();
	double viscOil = fInput.getReservoirData().OilVisc();
	double specificMassWater = 1000.;
	double satWater = 1. - satOil;
	double viscWater = fInput.getLiquidWaterViscosity();


	double r = getRadiusOfSteamFront(time);

	double massRate = getMassRateOfSteam(time);


	/**

	 * Valor da Pressao na Frente de Vapor ocasionada pela condensacao da agua!!!

	 */
	double Lg = log(r/re);
	double pi = M_PI;
	double DistrRightDown = pe - massRate/(2.*pi*r*H) * 1./(specificMassOil*K*satOil/viscOil + specificMassWater*K*satWater/viscWater) * Lg;

	return DistrRightDown;
}
//---------------------------------------------------------------------------

double TSwxSteam::getRadiusOfSteamFront(double tempo)
{
	double area = getRegionOfSteamArea(tempo);
	double radius = sqrt(area/M_PI);

	return radius;
}
//---------------------------------------------------------------------------

double TSwxSteam::getRegionOfSteamArea(double tempo) {
	if (tempo < 0.0 || IsZero(tempo) || fInput.getInjectionData().Table().size()	< 2) // O tempo nao pode ser negativo pois ele informa o tempo de injeção a uma determinada taxa, até x horas.
	{
		return 0.0;
	}

	double QiMinus1, Qi, ti, area = 0.0;
	int counter;
	// iterators
	std::vector<std::pair<double, double> >::iterator it = fInput.getInjectionData().Table().begin();
	it++;//comeca da segunda posicao!!!

	for (counter = 0; it != fInput.getInjectionData().Table().end();	it++, counter++)
	{
		ti = (*it).first;
		QiMinus1 = (*(it - 1)).second;
		Qi = (*it).second;
		area += ((Qi - QiMinus1) * getRegionAuxiliar(counter, tempo));
	}
	return area;
}
//---------------------------------------------------------------------------

double TSwxSteam::getMassRateOfSteam(const double tempo)
{
	std::map<double, double>::iterator it;
	double t, m;

	if(fabs(tempo) < 1.E-5)//IsZero(tempoFinal)
	{
		return 0.;
	}
	if(tempo < 0.)//Negative
	{
		return -1.;
	}
	it = fInput.getInjectionData().MassRate().end();
	it--;
	t = it->first;
	if(tempo > t)
	{
		//"Invalid Time in getMassRateOfSteam method!"
		return -999.9;
	}


	for(it = fInput.getInjectionData().MassRate().begin();
			it != fInput.getInjectionData().MassRate().end();
			it++)
	{
		t = it->first;
		m = it->second;

		if(tempo <= t)
		{
			break;
		}
	}

	return m;
}
//---------------------------------------------------------------------------

