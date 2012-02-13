/*
 *  ThermalMethodsTables.h
 *  FrenteVapor
 *
 *  Created by Jorge Calle on 4/2/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef ThermalMethodsTablesToSteamHHH
#define ThermalMethodsTablesToSteamHHH

//#define IsZero( a )    ( (a) < 1.e-8 && (a) > -1.e-8 )

#include <vector>

using namespace std;

// Devolve o valor da variável interpolada na coluna de variáveis independentes a partir do valor
// de uma variável dependente dada na coluna indexproperty
//double getIndexerFromValuesInterpolating(std::vector<std::vector<double> > &fTable,int indexproperty,double value);

// Classe que armazena a tabela de valores das propriedades físicas para algum material
// Pode ser indicado o indice da variavel independente, assumindo então todas as demais dependentes dela
// sugere-se utilizar: O índice 1 (coluna 1) para armazenar a variável da qual as propriedades físicas das outras colunas são dependentes
// O número de linhas é o total de entradas tabuladas
// O número de colunas-1 é o total de propriedades físicas tabuladas e dependentes do valor armazenado na coluna 1
class DataInTable {

protected:

	std::vector<std::vector<double> > fTable;

	// Indice da variavel independente
	int fIndependentVarIndex;

	// Método para interpolar o valor Temp, entre os valores tabulados no indice (i-1), Value(i-1) e no indice (i), Value(i)
	double getInterpolatedValue(double temp,int indexproperty);

	// Função que a partir dos pares (initialvar, initialvalue) e (finalvar, finalvalue) interpola o valor da variável dependente sobre o valor da independente - var.
	double getInterpolatedValue(double initialvar,double intialvalue,double finalvar,double finalvalue,double var);

	public:

	// Construtor com: ndata -> numero de dados associados (linhas) na tabela, nitens -> numero de valores relacionados para cada linha (colunas)
	// Primeira e segunda coluna correspondem a Temperatura e Pressão, os demais dependem do tipo de dado armazenado. Assim  nitens > 2.
	DataInTable(int ndata,int nitens);

	// Método que para uma dada temperatura obtemos o indice da linha na tabela (a temperatura é armazenada em intervalos iguais)
	int getIndex(double temp);
};

// Classe que construi uma tabela com as propriedades físicas da agua liquida e gasosa (vapor) em estado de saturação:
// Na tabela, o índice 1 (coluna 1) contem as temperaturas do estado de saturação [ C - graus Celsius ]
// e o índice 2 (coluna 2) contem a pressão associada para o estado de saturação  [ kPa - Kilo Pascal ]
// O índice 3 (coluna 3) armazena o calor latente para o correspondente estado de saturação. [ kJ/kg ]
// O calor latente para temperaturas superiores a 374 C é nulo, assim a água passa do estado líquido para o vapor diretamente.
// Os índices do 4 ao 8 armazenam propriedades físicas da água no estado aquoso (líquido)
// Os índices do 9 ao 13 armazenam propriedades físicas da água no estado gasoso (vapor)
// O índice 4(líquido) e 9(vapor) armazenam as densidades correspondentes [ kg/(m3) ]
// O índice 5(líquido) e 10(vapor) armazenam a viscosidade dinâmica [ Pa * sec ]
// O índice 6(líquido) e 11(vapor) armazenam o calor específico [ kJ/(kg * K) ]
// O índice 7(líquido) e 12(vapor) armazenam a condutividade térmica [ kJ/(sec*m*K) ]
// O índice 8(líquido) e 13(vapor) armazenam a entalphia específica [ kJ/kg ]
class WaterDataInStateOfSaturation : public DataInTable {

public:

	WaterDataInStateOfSaturation();

	// Métodos que devolvem as propriedades físicas: densidade, viscosidade, calor específico, condutividade térmica e entalphia específica 
	// a partir dos dados tabulados - para a água líquida em estado de saturação
	double getSaturationStateDensityToLiquidWater(double temp);
	double getSaturationStateSpecificHeatToLiquidWater(double temp);
	double getViscosityOfLiquidWater(double temp)//a equacao abaixo eh a funcao polinomial que interpola a viscosidade da agua no estado liquido.
	{
		//obs.: temp eh temperatura em graus Celsius.
		double visc = ( -1.7629190745815696e-20*(-147.99307632373458 + temp)*(12280.787335763045 + (-219.7423861580202 + temp)*temp)*
   (8144.360886155207 + (-166.69801220193315 + temp)*temp)*(4044.87108666885 + (-93.32992795418424 + temp)*temp)*
	 (1590.203901160502 + (-17.476286133991827 + temp)*temp)*(1064.6578436448397 + temp*(48.344763835214124 + temp)) )*1.E-3;

		return visc;
	}

	// Métodos que devolvem as propriedades físicas: densidade, viscosidade, calor específico, condutividade térmica e entalphia específica 
	// a partir dos dados tabulados - para a água gasosa (vapor) em estado de saturação
	double getSaturationStateDensityToSteam(double temp);
	double getSaturationStateViscosityToSteam(double temp);
	double getSaturationStateSpecificHeatToSteam(double temp);

	// Obtendo valores das saturações da água ou vapor dependendo da qualidade (titularidade) do vapor
	double getSaturationWater(int i, double quality, double So);
};
#endif 