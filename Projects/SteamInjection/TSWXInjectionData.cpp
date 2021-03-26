//---------------------------------------------------------------------------


#pragma hdrstop

#include "TSWXInjectionData.h"

//---------------------------------------------------------------------------

void TSwxInjectionData::setConstantTimeAndInjection(double time, double injection)
{
	std::pair<double, double> par(time, injection);
	fTable.push_back(par);
}

void TSwxInjectionData::setInjectionTemperature(double temperature)
{
	fSteamTemperature = temperature;
}

void TSwxInjectionData::setInjectionQuality(double quality)
{
	fQuality = quality;
}

double TSwxInjectionData::Temperature()
{
	return fSteamTemperature;
}

double TSwxInjectionData::getInjectionQuality()
{
	return fQuality;
}

void TSwxInjectionData::setInjectionData(std::map< double , std::pair< double, double > > &table,
																					double quality,
																					double temperature)
{
	if(!fTable.size())
	{
		setConstantTimeAndInjection(0., 0.);
	}
	std::map< double , std::pair< double, double > >::iterator itt;
	double taxaInj, massrate, tempoFinal;
	for(itt = table.begin(); itt != table.end(); itt++)
	{
		tempoFinal = itt->first;
		taxaInj = itt->second.first/1000.;//a formulacao contempla kJ/s
		massrate = itt->second.second;

		setConstantTimeAndInjection(tempoFinal, taxaInj);
		fMassRate[tempoFinal] = massrate;
	}

	setInjectionTemperature(temperature);
	setInjectionQuality(quality);
}

void TSwxInjectionData::setInjectionData(std::vector< std::pair<double, double> > &inj,
																				std::map<double,double> &MassRate,
																				double quality,
																				double temperature)
{
	fTable = inj;
	fMassRate = MassRate;
	setInjectionTemperature(temperature);
	setInjectionQuality(quality);
}

std::vector< std::pair<double, double> > & TSwxInjectionData::Table()
{
	return fTable;
}

std::map<double,double> & TSwxInjectionData::MassRate()
{
	return fMassRate;
}

std::vector<double> TSwxInjectionData::getSITimes()
{
	std::vector<double> times;
	int tableSize = fTable.size();
	for(int i = 0; i < tableSize; i++)
	{
		times.push_back(fTable[i].first);
	}
	return times;
}

#pragma package(smart_init)
