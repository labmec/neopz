//
//  TPBrAcidFunc.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/13/14.
//
//

#include "TPBrAcidFunc.h"
#include "pzmaterialid.h"


/// preencha valores default para os parametros
void TPBrAcidFunc::StandarParameters()
{
    /// porosidade inicial
    fPorosidadeInicial = 0.2;
    /// porosidade maxima admissivel
    fPorosidadeMaximaAdmissivel = 0.8;
    
    /// concentracao molar inicial do reagente
    fConcentracaoMolar = 4.4; // mol/litro
    /// Vazao de tratamento
    fVazaoTratamento = 20.; // bpm
    /// volume de fluido de tratamento
    fVolumeFluidoTratamento = 1500; // bbl
    /// tempo de injecao de tratamento
    fTempodeTratamento = -1; // calculado
    /// comprimento de trecho tratado
    fComprimentoTratamento = 100.; // m
    /// raio do poco
    fRaioPoco = 0.108; // m
    /// massa especifica do mineral soluvel
    fRhoM = 2.73; // g/cm3
    /// massa molar do mineral soluvel
    fMassaMolar = 100.087; // g/mol
    /// area superficial especifica dos poros
    fAreaSuperficialPoros = 50.; // m2/m3
    /// electro valencia do fluido em dissociacao
    fElectroValenciaFluido = 1.;
    /// electro valencai do mineral em dissociacao
    fElectroValenciaMineral = 2.;
    /// coeficiente estequiometrico, representando quantos mols do mineral sao dissolvidos por um mol de reagente
    fBetaEstequimetrico = -1.; // calculado
    /// constante de reacao
    fks = 0.2; // m/s
    /// coeficiente de transferencia de massa local
    fkc = 5.e-5; // m/s
    /// constante efetiva de reacao
    fkeff = -1.; // calculado
    /// modulo inicial de rigidez
    fEInicial = 29.459; // Mpa
    /// constante do modelo de reducao de rigidez
    fa = -8.909;
    /// constante do modelo de reducao de rigidez
    fb = 17.415;

}

void TPBrAcidFunc::CalculaDerivados()
{
    fTempodeTratamento = 60.*fVolumeFluidoTratamento/fVazaoTratamento;
    fkeff = fks*fkc/(fks+fkc);
    fBetaEstequimetrico = 0.001*fElectroValenciaFluido/fElectroValenciaMineral;
    STATE alfa = fBetaEstequimetrico*fMassaMolar;
    STATE r2 = 9.*fRaioPoco*fRaioPoco;
    STATE a = fPorosidadeInicial*alfa/fRhoM * fkeff * fAreaSuperficialPoros*fConcentracaoMolar*fTempodeTratamento;
    STATE c = -1185.60134*fkeff*fAreaSuperficialPoros/fVazaoTratamento *fPorosidadeInicial;
    STATE d = fComprimentoTratamento*(r2-fRaioPoco*fRaioPoco);
    STATE b = exp(c*d);
    STATE VariacaoPorosidade = a*b;
    
    STATE phiajust = std::min(fPorosidadeInicial+VariacaoPorosidade,fPorosidadeMaximaAdmissivel);
    STATE elast = fEInicial*exp(fa*fPorosidadeInicial+fb*phiajust);
}

/// preencha valores invalidos para os parametros
void TPBrAcidFunc::VoidParameters()
{
    /// porosidade inicial
    fPorosidadeInicial = -1;
    /// porosidade maxima admissivel
    fPorosidadeMaximaAdmissivel = -1;
    
    /// concentracao molar inicial do reagente
    fConcentracaoMolar = -1; // mol/litro
    /// Vazao de tratamento
    fVazaoTratamento = -1; // bpm
    /// volume de fluido de tratamento
    fVolumeFluidoTratamento = -1; // bbl
    /// tempo de injecao de tratamento
    fTempodeTratamento = -1; // calculado
    /// comprimento de trecho tratado
    fComprimentoTratamento = -1; // m
    /// raio do poco
    fRaioPoco = -1; // m
    /// massa especifica do mineral soluvel
    fRhoM = -1; // g/cm3
    /// massa molar do mineral soluvel
    fMassaMolar = -1; // g/mol
    /// area superficial especifica dos poros
    fAreaSuperficialPoros = -1; // m2/m3
    /// electro valencia do fluido em dissociacao
    fElectroValenciaFluido = -1;
    /// electro valencai do mineral em dissociacao
    fElectroValenciaMineral = -1;
    /// coeficiente estequiometrico, representando quantos mols do mineral sao dissolvidos por um mol de reagente
    fBetaEstequimetrico = -1.; // calculado
    /// constante de reacao
    fks = -1; // m/s
    /// coeficiente de transferencia de massa local
    fkc = -1; // m/s
    /// constante efetiva de reacao
    fkeff = -1.; // calculado
    /// modulo inicial de rigidez
    fEInicial = -1; // Mpa
    /// constante do modelo de reducao de rigidez
    fa = -1;
    /// constante do modelo de reducao de rigidez
    fb = -1;
    
}



/** @brief Define the class id associated with the class */
/**
 * This id has to be unique for all classes
 * A non unique id is flagged at the startup of the program
 */
int TPBrAcidFunc::ClassId() const
{
    return TPBrAcidFuncID;
}

/** @brief Writes this object to the TPZStream buffer. Include the classid if withclassid = true */
void TPBrAcidFunc::Write(TPZStream &buf, int withclassid) 
{
    TPZFunction<STATE>::Write(buf,withclassid);
//    buf.Write(&fRwell);
}

/** @brief read objects from the stream */
void TPBrAcidFunc::Read(TPZStream &buf, void *context)
{
    TPZFunction<STATE>::Read(buf,context);
//    buf.Read(&fRwell);
}

template class TPZRestoreClass<TPBrAcidFunc,TPBrAcidFuncID>;
