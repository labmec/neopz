//---------------------------------------------------------------------------


#include "TSWXSteam.h"

#include "TPZGuiInterface.h"

int main()
{
	string FileName, dirname = PZSOURCEDIR;
	FileName = dirname + "/Projects/SteamInjection/";
	FileName += "steamInjectionFILE.txt";
	
	std::ifstream input(FileName.c_str());
	
	TSwxSteam steam;
	std::map< double , std::pair<double, double> > Time_Radius_MaxSigmaTheta;
	TPZGuiInterface * progressInfo = new TPZGuiInterface;
	
	steam.ReadMe(input);
	steam.getRadiusAndMaxSigmaThetaForTableTimes(Time_Radius_MaxSigmaTheta, progressInfo);
	
	std::ofstream outMath("steamInjectionMATHEMATICA.nb");
	steam.PrintToMathematicaFile(Time_Radius_MaxSigmaTheta, outMath);

	return 0;
};


