//
//  TPBrAcidFunc.h
//  PZ
//
//  Created by Philippe Devloo on 5/13/14.
//
//

#ifndef __PZ__TPBrAcidFunc__
#define __PZ__TPBrAcidFunc__

#include <iostream>
#include "pzfunction.h"


class TPBrAcidFunc : public TPZFunction<STATE>
{
    
public:
    /// porosidade inicial
    STATE fPorosidadeInicial;
    /// porosidade maxima admissivel
    STATE fPorosidadeMaximaAdmissivel;
    /// concentracao molar inicial do reagente
    STATE fConcentracaoMolar;
    /// Vazao de tratamento
    STATE fVazaoTratamento;
    /// volume de fluido de tratamento
    STATE fVolumeFluidoTratamento;
    /// tempo de injecao de tratamento
    STATE fTempodeTratamento;
    /// comprimento de trecho tratado
    STATE fComprimentoTratamento;
    /// raio do poco
    STATE fRaioPoco;
    /// massa especifica do mineral soluvel
    STATE fRhoM;
    /// massa molar do mineral soluvel
    STATE fMassaMolar;
    /// area superficial especifica dos poros
    STATE fAreaSuperficialPoros;
    /// electro valencia do fluido em dissociacao
    STATE fElectroValenciaFluido;
    /// electro valencai do mineral em dissociacao
    STATE fElectroValenciaMineral;
    /// coeficiente estequiometrico, representando quantos mols do mineral sao dissolvidos por um mol de reagente
    STATE fBetaEstequimetrico;
    /// constante de reacao
    STATE fks;
    /// coeficiente de transferencia de massa local
    STATE fkc;
    /// constante efetiva de reacao
    STATE fkeff;
    /// modulo inicial de rigidez
    STATE fEInicial;
    /// constante do modelo de reducao de rigidez
    STATE fa;
    /// constante do modelo de reducao de rigidez
    STATE fb;
    
    
public:
    
    
    
    /** @brief Class constructor */
    TPBrAcidFunc()
    {
        VoidParameters();
    }


    /** @brief Class copy constructor */
    TPBrAcidFunc(const TPBrAcidFunc &cp) : TPZFunction<STATE>(cp)
    {
        fPorosidadeInicial=cp.fPorosidadeInicial;
        fPorosidadeMaximaAdmissivel=cp.fPorosidadeMaximaAdmissivel;
        fConcentracaoMolar=cp.fConcentracaoMolar;
        fVazaoTratamento=cp.fVazaoTratamento;
        fVolumeFluidoTratamento=cp.fVolumeFluidoTratamento;
        fTempodeTratamento=cp.fTempodeTratamento;
        fComprimentoTratamento=cp.fComprimentoTratamento;
        fRaioPoco=cp.fRaioPoco;
        fRhoM=cp.fRhoM;
        fMassaMolar=cp.fMassaMolar;
        fAreaSuperficialPoros=cp.fAreaSuperficialPoros;
        fElectroValenciaFluido=cp.fElectroValenciaFluido;
        fElectroValenciaMineral=cp.fElectroValenciaMineral;
        fBetaEstequimetrico=cp.fBetaEstequimetrico;
        fks=cp.fks;
        fkc=cp.fkc;
        fkeff=cp.fkeff;
        fEInicial=cp.fEInicial;
        fa=cp.fa;
        fb=cp.fb;
    }
    

    /** @brief Class copy constructor */
    TPBrAcidFunc & operator=(const TPBrAcidFunc &cp)
    {
        TPZFunction<STATE>::operator=(cp);
        fPorosidadeInicial=cp.fPorosidadeInicial;
        fPorosidadeMaximaAdmissivel=cp.fPorosidadeMaximaAdmissivel;
        fConcentracaoMolar=cp.fConcentracaoMolar;
        fVazaoTratamento=cp.fVazaoTratamento;
        fVolumeFluidoTratamento=cp.fVolumeFluidoTratamento;
        fTempodeTratamento=cp.fTempodeTratamento;
        fComprimentoTratamento=cp.fComprimentoTratamento;
        fRaioPoco=cp.fRaioPoco;
        fRhoM=cp.fRhoM;
        fMassaMolar=cp.fMassaMolar;
        fAreaSuperficialPoros=cp.fAreaSuperficialPoros;
        fElectroValenciaFluido=cp.fElectroValenciaFluido;
        fElectroValenciaMineral=cp.fElectroValenciaMineral;
        fBetaEstequimetrico=cp.fBetaEstequimetrico;
        fks=cp.fks;
        fkc=cp.fkc;
        fkeff=cp.fkeff;
        fEInicial=cp.fEInicial;
        fa=cp.fa;
        fb=cp.fb;
        return *this;
    }
  
	/** @brief Class destructor */
	virtual ~TPBrAcidFunc()
    {
        VoidParameters();
    }
    
    /// preencha valores default para os parametros
    void StandardParameters();
    
    /// preencha valores invalidos para os parametros
    void VoidParameters();
    
    /// calcula os valores derivados das constantes fornecidas
    void CalculaDerivados();
    
	/**
	 * @brief Performs function computation
	 * @param x point coordinate which is suppose to be in real coordinate system but can be in master coordinate system in derived classes.
	 * @param f function values returns the gradient of the pressure multiplied by biot
	 * @param df pressure function multiplied by biot
	 */
	virtual void Execute(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE> &df)
    {
        DebugStop();
    }
	
	/**
	 * @brief Performs time dependent function computation
	 * @param x point coordinate which is suppose to be in real coordinate system but can be in master coordinate system in derived classes.
	 * @param ftime  time to evaluate
	 * @param f function values
	 * @param gradf function derivatives
	 */
	virtual void Execute(const TPZVec<REAL> &x, REAL ftime, TPZVec<STATE> &f, TPZFMatrix<STATE> &gradf){
        DebugStop();
    }
    
    /** @brief Execute method receiving axes. It is used in shape functions */
    virtual void Execute(const TPZVec<REAL> &x, const TPZFMatrix<REAL> &axes, TPZVec<STATE> &f, TPZFMatrix<STATE> &df){
        DebugStop();
    }
    
    /** @brief Simpler version of Execute method which does not compute function derivatives */
    virtual void Execute(const TPZVec<REAL> &x, TPZVec<STATE> &f){
        REAL r2 = x[0]*x[0]+x[1]*x[1];
        STATE alfa = fBetaEstequimetrico*fMassaMolar;
        STATE a = fPorosidadeInicial*alfa/fRhoM * fkeff * fAreaSuperficialPoros*fConcentracaoMolar*fTempodeTratamento;
        STATE c = -1185.60134*fkeff*fAreaSuperficialPoros/fVazaoTratamento *fPorosidadeInicial;
        STATE d = fComprimentoTratamento*(r2-fRaioPoco*fRaioPoco);
        STATE b = exp(c*d);
        STATE VariacaoPorosidade = a*b;
        
        STATE phiajust = MIN(fPorosidadeInicial+VariacaoPorosidade,fPorosidadeMaximaAdmissivel);
        STATE youngnovo = fEInicial*exp(fa*fPorosidadeInicial+fb*(phiajust-fPorosidadeInicial));
//        youngnovo  = 12000.;
        f[0] = youngnovo;
    }
    
	/** @brief Returns number of functions. */
	virtual int NFunctions(){
        return 1;
    }
	
	/** @brief Polynomial order of this function. */
	/** In case of non-polynomial function it can be a reasonable approximation order. */
	virtual int PolynomialOrder()
    {
        return 4;
    }
    
    /** @brief Print a brief statement */
    virtual void Print(std::ostream &out)
    {
        out << __PRETTY_FUNCTION__ << (void *) this << std::endl;
        out << "NFunctions = " << NFunctions() << std::endl;
        out << "Polynomial Order = " << PolynomialOrder() << std::endl;
        /// porosidade inicial
        out << "Porosidade inicial " << fPorosidadeInicial << std::endl;
        /// porosidade maxima admissivel
        out << "Porosidade maxima " <<  fPorosidadeMaximaAdmissivel << std::endl;
        /// area superficial especifica dos poros
        out << "Area superficial dos poros " <<  fAreaSuperficialPoros << std::endl;
        /// concentracao molar inicial do reagente
        out << "Concentracao molar " <<  fConcentracaoMolar << std::endl;
        /// massa molar do mineral soluvel
        out << "Massa molar do mineral soluvel " <<  fMassaMolar << std::endl;
        /// massa especifica do mineral soluvel
        out << "Massa especifica do mineral soluvel " << fRhoM << std::endl;
        /// electro valencia do fluido em dissociacao
        out << "Eletro valencia do fluido " << fElectroValenciaFluido << std::endl;
        /// electro valencai do mineral em dissociacao
        out << "Eletro valencia do mineral " <<  fElectroValenciaMineral << std::endl;
        /// constante de reacao
        out << "constante de reacao " << fks << std::endl;
        /// coeficiente de transferencia de massa local
        out << "coeficiente de transferencia de massa local " << fkc << std::endl;
        /// coeficiente estequiometrico, representando quantos mols do mineral sao dissolvidos por um mol de reagente
        out << "coeficiente estequiometrico " <<  fBetaEstequimetrico << std::endl;
        /// constante efetiva de reacao
        out << "constante efetiva de reacao " <<  fkeff << std::endl;
        /// Vazao de tratamento
        out << "Vazao de tratamento " <<  fVazaoTratamento << std::endl;
        /// volume de fluido de tratamento
        out << "Volume fluido tratamento " <<  fVolumeFluidoTratamento << std::endl;
        /// tempo de injecao de tratamento
        out << "Tempo de tratamento " << fTempodeTratamento << std::endl;
        /// comprimento de trecho tratado
        out << "Comprimento de tratamento " <<  fComprimentoTratamento << std::endl;
        /// raio do poco
        out << "Raio do poco " << fRaioPoco << std::endl;
        /// modulo inicial de rigidez
        out << "modulo inicial de rigidez " <<  fEInicial << std::endl;
        /// constante do modelo de reducao de rigidez
        out << "constante a do modelo de reducao de rigidez " <<  fa << std::endl;
        /// constante do modelo de reducao de rigidez
        out << "constante b do modelo de reducao de rigidez " <<  fb << std::endl;
    }
    
    /** @brief Define the class id associated with the class */
	/**
	 * This id has to be unique for all classes
	 * A non unique id is flagged at the startup of the program
	 */
	virtual int ClassId() const ;
	
	/** @brief Writes this object to the TPZStream buffer. Include the classid if withclassid = true */
	virtual void Write(TPZStream &buf, int withclassid) ;
	
	/** @brief read objects from the stream */
	virtual void Read(TPZStream &buf, void *context);
	


};

#endif /* defined(__PZ__TPBrAcidFunc__) */
