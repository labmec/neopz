//---------------------------------------------------------------------------

#ifndef TSWXConfinementDataH
#define TSWXConfinementDataH
//---------------------------------------------------------------------------

class TSwxConfinementData
{

	public:

	TSwxConfinementData()
	{
	}

	void SetData(double condtermica, double densidade, double calorespecif);
	TSwxConfinementData(TSwxConfinementData &cpy);

	// Método para recuperar o produto das propriedades físicas Condutividade_Termica * Densidade * Calor_Especifico
	double getProductOfTheProperties();

	double ThermalConductivity() { return fThermalConductivity; }
	double Density() { return fDensity; }
	double SpecificHeat() { return fSpecificHeat; }

	private :
		/**
	 * fThermalConductivity: Condutividade tÈrmica [J/(s m C)]
	 * fDensity: Massa especÌfica [kg/m3]
	 * fSpecificHeat: Calor especÌfico [J/(kg C)]
	 */
	double fThermalConductivity, fDensity, fSpecificHeat;
};

#endif
