//---------------------------------------------------------------------------

#ifndef TSWXInjectionDataH
#define TSWXInjectionDataH

#include <map>
#include <vector>

//---------------------------------------------------------------------------

class TSwxInjectionData
{
	public:

	TSwxInjectionData()
	{

		 fTable.resize(0);
		 fMassRate.clear();
	}

	void setInjectionData(std::map< double , std::pair< double, double > > &table,
												double quality,
												double temperature);

	void setInjectionData(std::vector< std::pair<double, double> > &inj,
												std::map<double,double> &MassRate,
												double quality,
												double temperature);

	double Temperature();
	double getInjectionQuality();
	std::vector< std::pair<double, double> > & Table();
	std::map<double,double> & MassRate();
	std::vector<double> getSITimes();

	private:

	void setConstantTimeAndInjection(double time, double injection);
	void setInjectionTemperature(double temperature);
	void setInjectionQuality(double quality);

	/**
	 * fTable: tabela de injecao. Tempo [s] x Vazao [m3/s]
	 */
	std::vector< std::pair<double, double> > fTable;

	// fMassrate <instantaneous time , massrate>
	std::map<double,double> fMassRate;

	double fQuality;
	double fSteamTemperature;
};

#endif
