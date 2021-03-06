//
//  pznlfluidstructureData.cpp
//  PZ
//
//  Created by Cesar Lucci on 07/08/13.
//
//

#include "pznlfluidstructureData.h"
#include "pzreal.h"
#include "pzcompel.h"
#include "pzintel.h"
#include "pzcmesh.h"

InputDataStruct::InputDataStruct()
{
  fPressureMatIds_StripeId_ElastId.clear();
  fEModel = EElastic;
  fIsPropag = false;
	fNthreadsForAssemble = 0;
  fLastFracMatId = -6378;
	fDeltaC = 0.;
	fDeltaT = 0.;
	fSigmaT = 0.;
  fnElPropag = 0;
  fusingLeakOff = false;
}

InputDataStruct::~InputDataStruct()
{
  
}
void InputDataStruct::SetMohrCoulombData(REAL cohesion, REAL phiMC)
{
  fEModel = EMohrCoulomb;
  fCohesion = cohesion;
  fPhiMC = phiMC;
}

void InputDataStruct::SetSandlerData()
{
  fEModel = ESandler;
}

void InputDataStruct::SetData(REAL Lx, REAL Ly, REAL Lf, REAL Hf, REAL Lmax_edge, REAL E1, REAL Poisson1, REAL E2, REAL Poisson2, REAL XinterfaceBetween1and2,
                              REAL Fx, REAL Fy, REAL preStressXX, REAL preStressXY, REAL preStressYY,
                              int NStripes, REAL Visc, REAL SigN, REAL QinjTot, REAL Ttot, REAL maxDeltaT, int nTimes,
                              REAL Cl, REAL Pe, REAL SigmaConf, REAL Pref, REAL vsp, REAL KIc, REAL Jradius, int ndivV, int ndivH, REAL q, 
															REAL DeltaC, REAL DeltaT, REAL SigmaT,
															int NThreadsForAssemble)
{
	fNthreadsForAssemble = NThreadsForAssemble;
	fDeltaC = DeltaC;
	fDeltaT = DeltaT;
	fSigmaT = SigmaT;
	
  fLx = Lx;
  fLy = Ly;
  fLf = Lf;
  fop = 0.;
  fHf = Hf;
  fndivV = ndivV;
  fndivH = ndivH;
  fq = q;
  fLmax_edge = Lmax_edge;
  
  fE1 = E1;
  fPoisson1 = Poisson1;
  fE2 = E2;
  fPoisson2 = Poisson2;
  fXinterface = XinterfaceBetween1and2;
  
  REAL leftLimit = (Lmax_edge-1.E-10);
  //REAL rightLimit = (Lx - Lmax_edge)  -1.E-10;
  if(XinterfaceBetween1and2 < leftLimit)// || XinterfaceBetween1and2 > rightLimit) NOT GOING TO USE INTERFACE
  {
    std::cout << "A interface deve estar entre Lmax_edge e (Lx - Lmax_edge)!!!\n\n";
    DebugStop();
  }
  
  fFx = Fx;
  fFy = Fy;
  fPreStressXX = preStressXX;
  fPreStressXY = preStressXY;
  fPreStressYY = preStressYY;
  fNStripes = NStripes;
  
  fPressureMatIds_StripeId_ElastId.clear();
  
  fVisc = Visc;
  
  fLeakoffmap.clear();
  
  fSigN = SigN;
  
  REAL Qinj1asa = QinjTot / 2.;
  REAL QinjSecao = Qinj1asa / Hf;
  fQinj = QinjSecao;
  
  fTtot = Ttot;
  factTime = fminDeltaT;
  fmaxDeltaT = maxDeltaT;
  fNDeltaTsteps = nTimes;
  fminDeltaT = fmaxDeltaT/fNDeltaTsteps;
  factDeltaT = fmaxDeltaT;
  ftimeStep = 0;
  
  fCl = Cl;
  fPe = Pe;
  fSigmaConf = SigmaConf;
  fPref = Pref;
  fvsp = vsp;
  
  fKIc = KIc;
  fJradius = Jradius;
  globFractOutputData.fStepDeltaT[ftimeStep] = factDeltaT;
  globFractOutputData.fStepTAcum[ftimeStep] = 0.;
  globFractOutputData.fTxLfrac[factTime] = globFractInputData.Lf();
  globFractOutputData.fTxOp[factTime] = 0.;
}

int & InputDataStruct::GetnElPropag(){
	return fnElPropag;
}

void InputDataStruct::SetUsingLeakOff(bool usingLeakOff){
  fusingLeakOff = usingLeakOff;
}

bool InputDataStruct::GetIfUsingLeakOff(){
  return fusingLeakOff;
}

REAL InputDataStruct::DeltaC(){
	return fDeltaC;
}

REAL InputDataStruct::DeltaT(){
	return fDeltaT;
}

REAL InputDataStruct::SigmaT(){
	return fSigmaT;
}

void InputDataStruct::SetLastFracMatId(int matid)
{
  fLastFracMatId = matid;
}

int InputDataStruct::GetLastFracMatId()
{
  return fLastFracMatId;
}

REAL InputDataStruct::NthreadsForAssemble()
{
	return fNthreadsForAssemble;
}

void InputDataStruct::SetLf(REAL Lf)
{
  fLf = Lf;
}

void InputDataStruct::SetOpening(REAL op)
{
  fop = op;
}

REAL InputDataStruct::Lx()
{
  return fLx;
}

REAL InputDataStruct::Ly()
{
  return fLy;
}

REAL InputDataStruct::Lf()
{
  return fLf;
}

REAL InputDataStruct::Opening()
{
  return fop;
}

REAL InputDataStruct::Hf()
{
  return fHf;
}

int InputDataStruct::NdivV()
{
  return fndivV;
}

int InputDataStruct::NdivH()
{
  return fndivH;
}
REAL InputDataStruct::GMeshq()
{
  return fq;
}

REAL InputDataStruct::Lmax_edge()
{
  return fLmax_edge;
}

REAL InputDataStruct::E1()
{
  return fE1;
}

REAL InputDataStruct::Poisson1()
{
  return fPoisson1;
}

REAL InputDataStruct::E2()
{
  return fE2;
}

REAL InputDataStruct::Poisson2()
{
  return fPoisson2;
}

REAL InputDataStruct::Xinterface()
{
  return fXinterface;
}

REAL InputDataStruct::Fx()
{
  return fFx;
}

REAL InputDataStruct::Fy()
{
  return fFy;
}

REAL InputDataStruct::PreStressXX()
{
  return fPreStressXX;
}

REAL InputDataStruct::PreStressXY()
{
  return fPreStressXY;
}

REAL InputDataStruct::PreStressYY()
{
  return fPreStressYY;
}

bool InputDataStruct::IsElastic()
{
  return (fEModel == EElastic);
}
bool InputDataStruct::IsMohrCoulomb()
{
  return (fEModel == EMohrCoulomb);
}
bool InputDataStruct::IsSandler()
{
  return (fEModel == ESandler);
}

REAL InputDataStruct::Cohesion()
{
  return fCohesion;
}

REAL InputDataStruct::PhiMC()
{
  return fPhiMC;
}

int InputDataStruct::NStripes()
{
  return fNStripes;
}


std::map< int,std::pair<int,int> > & InputDataStruct::GetPressureMatIds_StripeId_ElastId()
{
  return fPressureMatIds_StripeId_ElastId;
}

std::map< int,std::pair<int,int> > & InputDataStruct::GetfMatID_Rec_GeoEl()
{
  return fMatID_Rec_GeoEl;
}

int InputDataStruct::StripeId(int bcId)
{
  std::map< int,std::pair<int,int> >::iterator it = fPressureMatIds_StripeId_ElastId.find(bcId);
  if(it != fPressureMatIds_StripeId_ElastId.end())
  {
    return it->second.first;
  }
  else
  {
    DebugStop();
  }
  
  return -7456;
}

int InputDataStruct::ElastId(int bcId)
{
  std::map< int,std::pair<int,int> >::iterator it = fPressureMatIds_StripeId_ElastId.find(bcId);
  if(it != fPressureMatIds_StripeId_ElastId.end())
  {
    return it->second.second;
  }
  else
  {
    DebugStop();
  }
  
  return -7456;
}

void InputDataStruct::InsertBCId_StripeId_ElastId(int BCId, int StripeId, int ElastId)
{
  fPressureMatIds_StripeId_ElastId[BCId] = std::make_pair(StripeId,ElastId);
}

bool InputDataStruct::IsBC(int matId)
{
  std::map< int,std::pair<int,int> >::iterator it = fPressureMatIds_StripeId_ElastId.find(matId);
  if(it != fPressureMatIds_StripeId_ElastId.end())
  {
    return true;
  }
  else
  {
    return false;
  }
}

REAL InputDataStruct::Visc()
{
  return fVisc;
}

std::map<int,REAL> & InputDataStruct::GetLeakoffmap()
{
  return fLeakoffmap;
}

REAL InputDataStruct::SigN()
{
  return fSigN;
}

REAL InputDataStruct::Qinj()
{
  return fQinj;
}

REAL InputDataStruct::Ttot()
{
  return fTtot;
}

REAL InputDataStruct::actTime()
{
  return factTime;
}

int InputDataStruct::actTimeStep()
{
  return ftimeStep;
}

REAL InputDataStruct::actDeltaT()
{
  return factDeltaT;
}

REAL InputDataStruct::Cl()
{
  return fCl;
}

REAL InputDataStruct::Pe()
{
  return fPe;
}

REAL InputDataStruct::SigmaConf()
{
  return fSigmaConf;
}

REAL InputDataStruct::Pref()
{
  return fPref;
}

REAL InputDataStruct::vsp()
{
  return fvsp;
}

REAL InputDataStruct::KIc()
{
  return fKIc;
}

REAL InputDataStruct::Jradius()
{
  return fJradius;
}

void InputDataStruct::SetMinDeltaT()
{
  factDeltaT = fminDeltaT;
  if(factTime + factDeltaT > fTtot)
  {
    factDeltaT = fTtot - factTime;
  }
}

void InputDataStruct::SetNextDeltaT()
{
  factDeltaT = MIN(fmaxDeltaT,(fmaxDeltaT/fNDeltaTsteps));
  if(factTime + factDeltaT > fTtot)
  {
    factDeltaT = fTtot - factTime;
  }
  globFractOutputData.fStepDeltaT[ftimeStep] = factDeltaT;
  
}

void InputDataStruct::UpdateActTime()
{
  factTime += factDeltaT;
  ftimeStep++;
  std::cout << "\n\n=============== ActTime = " << factTime << " ===============\n\n";
  globFractOutputData.fStepTAcum[ftimeStep] = factTime;
  globFractOutputData.fTxLfrac[factTime] = globFractInputData.Lf();
  globFractOutputData.fTxOp[factTime] = globFractInputData.Opening();
}

REAL InputDataStruct::MaxDeltaT()
{
  return fmaxDeltaT;
}

////////////////////////////////////////////////////////////////// Leakoff


void InputDataStruct::UpdateLeakoff(TPZCompMesh * cmesh)
{
  
  
#ifdef PZDEBUG
  if(fLeakoffmap.size() == 0)
  {//Se a fratura nao alcancou ainda a regiao elastica 2, este mapa estah vazio!!!
   //DebugStop();
  }
#endif
  
  std::map<int,REAL>::iterator it;
  
  int outVlCount = 0;
  for(int i = 0;  i < cmesh->ElementVec().NElements(); i++)
  {
    ///////////////////////
    TPZCompEl * cel = cmesh->ElementVec()[i];
    
#ifdef PZDEBUG
    if(!cel)
    {
      DebugStop();
    }
#endif
    
    TPZGeoEl * gel = cel->Reference();
    
    if(gel->Dimension() != 1)
    {
      continue;
    }
    
    TPZInterpolatedElement * sp = dynamic_cast <TPZInterpolatedElement*> (cel);
    if(!sp)
    {
      continue;
    }
    
    it = globFractInputData.GetLeakoffmap().find(gel->Id());
    
    if(it == globFractInputData.GetLeakoffmap().end())
    {
      continue;
    }
    
    TPZVec<REAL> qsi(1,0.);
    cel->Reference()->CenterPoint(cel->Reference()->NSides()-1, qsi);
    TPZMaterialData data;
    sp->InitMaterialData(data);
    
    sp->ComputeShape(qsi, data);
    sp->ComputeSolution(qsi, data);
    
    REAL pfrac = data.sol[0][0];
    ///////////////////////
    
    REAL deltaT = globFractInputData.actDeltaT();
    
    REAL VlAcum = it->second;
    REAL tStar = FictitiousTime(VlAcum, pfrac);
    REAL Vlnext = VlFtau(pfrac, tStar + deltaT);
    
    if (fusingLeakOff) {
      it->second = Vlnext;
    }
    else{
      it->second = 0.;
    }

    outVlCount++;
  }
  
#ifdef PZDEBUG
  if(outVlCount < globFractInputData.GetLeakoffmap().size())
  {
    DebugStop();
  }
#endif
}

REAL InputDataStruct::VlFtau(REAL pfrac, REAL tau)
{
  REAL Cl = globFractInputData.Cl();
  REAL sigmaConf = globFractInputData.SigmaConf();
  REAL Pe = globFractInputData.Pe();
  REAL Pref = globFractInputData.Pref();
  REAL vsp = globFractInputData.vsp();
  
  REAL Pef = pfrac + sigmaConf;
  REAL Pcalc = (Pef - Pe)/Pref;
  if(Pcalc < 0.)
  {
    Pcalc = 0.;
  }
  
  REAL Clcorr = Cl * sqrt(Pcalc);
  REAL Vl = 2. * Clcorr * sqrt(tau) + vsp;
  
  return Vl;
}

REAL InputDataStruct::FictitiousTime(REAL VlAcum, REAL pfrac)
{
  REAL Cl = globFractInputData.Cl();
  REAL sigmaConf = globFractInputData.SigmaConf();
  REAL Pe = globFractInputData.Pe();
  REAL Pref = globFractInputData.Pref();
  REAL vsp = globFractInputData.vsp();
  
  REAL tStar = 0.;
  if(VlAcum > vsp)
  {
    REAL Pef = pfrac + sigmaConf;
    REAL Pcalc = (Pef - Pe)/Pref;
    if(Pcalc < 0.)
    {
      Pcalc = 0.;
    }
    REAL Clcorr = Cl * sqrt(Pcalc);
    tStar = (VlAcum - vsp)*(VlAcum - vsp)/( (2. * Clcorr) * (2. * Clcorr) );
  }
  
  return tStar;
}

REAL InputDataStruct::QlFVl(int gelId, REAL pfrac)
{
  std::map<int,REAL>::iterator it = globFractInputData.GetLeakoffmap().find(gelId);
  if(it == globFractInputData.GetLeakoffmap().end())
  {
    globFractInputData.GetLeakoffmap()[gelId] = 0.;//Nao coloque vsp! Eh ZERO mesmo!
    it = globFractInputData.GetLeakoffmap().find(gelId);
  }
  REAL VlAcum = it->second;
  
  REAL deltaT = globFractInputData.actDeltaT();
  
  REAL tStar = FictitiousTime(VlAcum, pfrac);
  REAL Vlnext = VlFtau(pfrac, tStar + deltaT);
  REAL Ql = (Vlnext - VlAcum)/deltaT;
  
  if(fusingLeakOff){
    return Ql;
  }
  else{
    return 0.;
  }
}

REAL InputDataStruct::dQlFVl(int gelId, REAL pfrac)
{
  std::map<int,REAL>::iterator it = fLeakoffmap.find(gelId);
  if(it == fLeakoffmap.end())
  {
    fLeakoffmap[gelId] = 0.;
    it = fLeakoffmap.find(gelId);
  }
  REAL VlAcum = it->second;
  
  REAL deltaPfrac = fabs(pfrac/10000.);
  if(deltaPfrac < 1.E-10)
  {
    deltaPfrac = 1.E-10;
  }
  else if(deltaPfrac > 1.E-3)
  {
    deltaPfrac = 1.E-3;
  }
  
  REAL deltaT = globFractInputData.actDeltaT();
  /////////////////////////////////////////////////Ql maior
  REAL pfracUP = pfrac + deltaPfrac;
  REAL tStar1 = FictitiousTime(VlAcum, pfracUP);
  REAL Vlnext1 = VlFtau(pfracUP, tStar1 + deltaT);
  REAL Ql1 = (Vlnext1 - VlAcum )/deltaT;
  //...
  
  /////////////////////////////////////////////////Ql menor
  REAL pfracDOWN = pfrac - deltaPfrac;
  REAL tStar0 = FictitiousTime(VlAcum, pfracDOWN);
  REAL Vlnext0 = VlFtau(pfracDOWN, tStar0 + deltaT);
  REAL Ql0 = (Vlnext0 - VlAcum)/deltaT;
  //...
  
  REAL dQldpfrac = (Ql1-Ql0)/(2.*deltaPfrac);
  
  if(fusingLeakOff){
    return dQldpfrac;
  }
  else{
    return 0.;
  }
}

// Propagation
void InputDataStruct::SetPropagated(){
  fIsPropag = true;
}
void InputDataStruct::SetNotPropagated(){
  fIsPropag = false;
}
bool InputDataStruct::IsPropagated(){
  return fIsPropag;
}

//------------------------------------------------------------

OutputDataStruct::OutputDataStruct()
{
  fTposP.clear();
  fTposVolLeakoff.clear();
  fTAcumVolW.clear();
  fTAcumVolLeakoff.clear();
  fTKI.clear();
  fTxLfrac.clear();
  fTxOp.clear();
  fQinj1wing = 0.;
  fLfracMax = 0.;
}

OutputDataStruct::~OutputDataStruct()
{
  fTposP.clear();
  fTposVolLeakoff.clear();
  fTAcumVolW.clear();
  fTAcumVolLeakoff.clear();
  fTKI.clear();
  fTxLfrac.clear();
  fTxOp.clear();
  fQinj1wing = 0.;
  fLfracMax = 0.;
}

int OutputDataStruct::NTimes()
{
  int ntimes0 = fTposP.size();
  
#ifdef PZDEBUG
  int ntimes1 = fTposVolLeakoff.size();
  int ntimes2 = fTAcumVolW.size();
  int ntimes3 = fTAcumVolLeakoff.size();
  int ntimes4 = fTKI.size();
  if(ntimes0 != ntimes1 || ntimes0 != ntimes2 || ntimes0 != ntimes3 || ntimes0 != ntimes4)
  {
    //Todos tem que ter o mesmo tamanho!!!
    DebugStop();
  }
#endif
  
  return ntimes0;
}

void OutputDataStruct::InsertTposP(int time, std::map<REAL,REAL> & posPmap)
{
  std::map<int,posP>::iterator it = fTposP.find(time);
  if(it == fTposP.end())
  {
    posP newposP;
    newposP.fposP = posPmap;
    fTposP[time] = newposP;
  }
  else{
    DebugStop();
  }
  
}

void OutputDataStruct::InsertTposVolLeakoff(int time, REAL pos, REAL Ql)
{
  std::map<int,posVolLeakoff>::iterator it = fTposVolLeakoff.find(time);
  if(it == fTposVolLeakoff.end())
  {
    posVolLeakoff newposVolLeakoff;
    newposVolLeakoff.InsertPoint(pos, Ql);
    fTposVolLeakoff[time] = newposVolLeakoff;
  }
  else
  {
    it->second.InsertPoint(pos, Ql);
  }
}

void OutputDataStruct::InsertTAcumVolW(int time, REAL vol)
{
  fTAcumVolW[time] = vol;
}

void OutputDataStruct::InsertTAcumVolLeakoff(int time, REAL vol)
{
  fTAcumVolLeakoff[time] = vol;
}

void OutputDataStruct::InsertTKI(int time, REAL KI)
{
  fTKI[time] = KI;
}

void OutputDataStruct::SetQinj1WingAndLfracmax(REAL Qinj1wing, REAL Lfracmax)
{
  fQinj1wing = fabs(Qinj1wing);
  fLfracMax = Lfracmax;
}

void OutputDataStruct::PlotElasticVTK(TPZAnalysis * an, int anCount)
{
  std::string plotfile = "PapanastasiouEla4.vtk";
  if(anCount >= 0)
  {
    an->SetStep(anCount);
  }
  
  TPZManVector<std::string,10> scalnames(1), vecnames(1);
  scalnames[0] = "SigmaY";
  vecnames[0]  = "Displacement";
	
	const int dim = 2;
	int div = 0;
	an->DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an->PostProcess(div);
}

void OutputDataStruct::PrintMathematica(std::ofstream & outf)
{
  int sztonto = fTposP.size();
  for (int i = 0; i < sztonto; i++) {
    fTKI[i] = 1.;
  }
#ifdef PZDEBUG
  if(fTposP.size() == 0 || fTposVolLeakoff.size() == 0 || fTAcumVolW.size() == 0 || fTAcumVolLeakoff.size() == 0 || fTKI.size() == 0)
  {
    DebugStop();
  }
#endif
  
  std::cout << "fTposVolLeakoff.size() = " << fTposVolLeakoff.size() << std::endl;
  std::cout << "fTAcumVolLeakoff.size() = " << fTAcumVolLeakoff.size() << std::endl;
  
  std::map<int,posP>::iterator itTposP, itTposPLast = fTposP.end();
  itTposPLast--;
  
  std::map<int,posVolLeakoff>::iterator itTposVolLeakoff, itTposVolLeakoffLast = fTposVolLeakoff.end();
  itTposVolLeakoffLast--;
  
  std::map<int,REAL>::iterator itTAcumVolW, itTAcumVolWLast = fTAcumVolW.end();
  itTAcumVolWLast--;
  
  std::map<int,REAL>::iterator itTAcumVolLeakoff, itTAcumVolLeakoffLast = fTAcumVolLeakoff.end();
  itTAcumVolLeakoffLast--;
  
  std::map<int,REAL>::iterator itTKI, itTKILast = fTKI.end();
  itTKILast--;
  
  outf << "(* Output Fracture Propagation 1D *)\n";
  
  outf << "Nathan2014;\n\n";
  outf << "ntimes=" << NTimes() << ";\n";
  outf << "times={";
  for(itTposP = fTposP.begin(); itTposP != fTposP.end(); itTposP++)
  {
    outf << itTposP->first;
    if(itTposP != itTposPLast)
    {
      outf << ",";
    }
  }
  outf << "};\n\n";
  
  outf << "(* Position x Pressure along time *)\n";
  outf << "posvsPvsT={";
  for(itTposP = fTposP.begin(); itTposP != fTposP.end(); itTposP++)
  {
    itTposP->second.PrintMathematica(outf);
    if(itTposP != itTposPLast)
    {
      outf << ",";
    }
  }
  outf << "};\n\n";
  
  outf << "(* Position x Leakoff Volume along time *)\n";
  outf << "posvsVolleakoffvsT={";
  for(itTposVolLeakoff = fTposVolLeakoff.begin(); itTposVolLeakoff != fTposVolLeakoff.end(); itTposVolLeakoff++)
  {
    itTposVolLeakoff->second.PrintMathematica(outf);
    if(itTposVolLeakoff != itTposVolLeakoffLast)
    {
      outf << ",";
    }
  }
  outf << "};\n\n";
  
  outf << "TimeStepXDeltaT={";
  std::map<int,REAL>::iterator itstep, itstepLast = fStepTAcum.end();
  itstepLast--;
  for(itstep = fStepTAcum.begin(); itstep != fStepTAcum.end(); itstep++)
  {
    outf << "{" << itstep->first << "," << itstep->second << "}";
    if(itstep != itstepLast)
    {
      outf << ",";
    }
  }
  outf << "};\n\n";
  outf << "lastTime=TimeStepXDeltaT[[ntimes,2]];\n\n" << std::endl;
  
  outf << "(* time x W Volume *)\n";
  outf << "TvsVolW={";
  for(itTAcumVolW = fTAcumVolW.begin(); itTAcumVolW != fTAcumVolW.end(); itTAcumVolW++)
  {
    outf << "{" << "TimeStepXDeltaT[[" << itTAcumVolW->first+1 << ",2]]," << itTAcumVolW->second << "}";
    if(itTAcumVolW != itTAcumVolWLast)
    {
      outf << ",";
    }
  }
  outf << "};\n\n";
  
  outf << "(* time x Accumulated Leakoff Volume *)\n";
  outf << "TvsVolLeakoff={";
  for(itTAcumVolLeakoff = fTAcumVolLeakoff.begin(); itTAcumVolLeakoff != fTAcumVolLeakoff.end(); itTAcumVolLeakoff++)
  {
    outf << "{" << "TimeStepXDeltaT[[" << itTAcumVolLeakoff->first+1 << ",2]]," << itTAcumVolLeakoff->second << "}";
    if(itTAcumVolLeakoff != itTAcumVolLeakoffLast)
    {
      outf << ",";
    }
  }
  outf << "};\n\n";
  
  outf << "(* time x Lfrac *)\n";
  outf << "TvsLfrac={";
  std::map<REAL,REAL>::iterator itlfrac, itlfraclast = fTxLfrac.end();
  itlfraclast--;
  for(itlfrac = fTxLfrac.begin(); itlfrac != fTxLfrac.end(); itlfrac++)
  {
    outf << "{" << itlfrac->first << "," << itlfrac->second << "}";
    if(itlfrac != itlfraclast)
    {
      outf << ",";
    }
  }
  outf << "};\n\n";

  outf << "(* time x Opening *)\n";
  outf << "TvsOp={";
  std::map<REAL,REAL>::iterator itop, itoplast = fTxOp.end();
  itoplast--;
  for(itop = fTxOp.begin(); itop != fTxOp.end(); itop++)
  {
    outf << "{" << itop->first << "," << itop->second << "}";
    if(itop != itoplast)
    {
      outf << ",";
    }
  }
  outf << "};\n\n";
  
  outf << "(* Qinj 1 wing and Lfrac max *)\n";
  outf << "Qinj1wing=" << fQinj1wing << ";\n";
  outf << "LfracMax=" << fLfracMax << ";\n\n";
  
  
  REAL maxDeltaT = globFractInputData.MaxDeltaT();
  outf << "DeltaT=" << maxDeltaT << ";\n\n";
  
  outf << "maxp = Max[Transpose[Flatten[posvsPvsT, 1]][[2]]];\n";
  outf << "Manipulate[ListPlot[posvsPvsT[[t]], Joined -> True,PlotLabel ->\"Graphic A: Position x Pressure @ \" <> ToString[times[[t]]] <> \"s\",AxesLabel -> {\"position (m)\", \"pressure (Pa)\"}, Filling -> Axis,PlotRange -> {{0, LfracMax}, {0, maxp}}], {t, 1, ntimes, 1}]\n\n";
  
  outf << "posvsMeanP = Table[{Max[posvsPvsT[[o]]], Mean[posvsPvsT[[o]]][[2]]}, {o, 2,Length[posvsPvsT]}];\n";
  outf << "ListPlot[posvsMeanP, Joined -> True, PlotLabel -> \"Graphic B: Lfrac x Mean Pressure\",AxesLabel -> {\"Lfrac (m)\", \"Mean pressure (Pa)\"}, Filling ->Axis, AxesOrigin -> {Min[Transpose[posvsMeanP][[1]]],Min[Transpose[posvsMeanP][[2]]]}]\n\n";
  
  outf << "maxleakoff = Max[Transpose[Flatten[posvsVolleakoffvsT, 1]][[2]]];\n";
  outf << "Manipulate[ListPlot[posvsVolleakoffvsT[[t]], Joined -> True,PlotLabel ->\"Graphic C: Position x Leakoff penetration @ \" <> ToString[times[[t]]] <>\"s\", AxesLabel -> {\"position (m)\", \"leakoff penetration (m)\"},Filling -> Axis, PlotRange -> {{0, LfracMax}, {0, maxleakoff}}], {t, 1,ntimes, 1}]\n\n";
  
  outf << "maxinj = Qinj1wing*TimeStepXDeltaT[[ntimes,2]];\n";
  outf << "GrD = Plot[Qinj1wing*t, {t, 0, TimeStepXDeltaT[[ntimes,2]]},PlotLabel -> \"Graphic D: Time x Volume Injected\",AxesLabel -> {\"time (s)\", \"Volume injected (m3)\"},Filling -> Axis, FillingStyle -> Red,PlotRange -> {{0, lastTime}, {0, maxinj}}]\n\n";
  
  outf << "GrE = ListPlot[TvsVolW, Joined -> True,PlotLabel -> \"Graphic E: Time x Fracture Volume\",AxesLabel -> {\"time (s)\", \"Fracture volume (m3)\"},Filling -> Axis, FillingStyle -> Green,PlotRange -> {{0, lastTime}, {0, maxinj}}]\n\n";
  
  outf << "GrF = ListPlot[TvsVolLeakoff, Joined -> True,PlotLabel -> \"Graphic F: Time x Leakoff volume\",AxesLabel -> {\"time (s)\", \"Leakoff volume (m3)\"},Filling -> Axis, FillingStyle -> Blue,PlotRange -> {{0, lastTime}, {0, maxinj}}]\n\n";
  
  outf << "WplusLeakoff = {};\n";
  outf << "For[tt = 1, tt <= ntimes,\n";
  outf << "AppendTo[WplusLeakoff, {TimeStepXDeltaT[[tt,2]], TvsVolW[[tt, 2]] + TvsVolLeakoff[[tt, 2]]}];\n";
  outf << "tt++;\n";
  outf << "];\n";
  outf << "GrG = ListPlot[WplusLeakoff, Joined -> False,PlotStyle -> {Black, PointSize[0.03]},PlotLabel -> \"Graphic G: Grahics (E+F)\",AxesLabel -> {\"time (s)\", \"Vol graphics(D+E)\"},PlotRange -> All];\n";
  outf << "Show[GrD, GrE, GrF, GrG, PlotRange -> All,PlotLabel -> \"Graphic G: Grahics D, E, F and (E+F)\"]\n\n";
  
  outf << "maxlfrac = TvsLfrac[[-1,2]];\n";
  outf << "GrH = ListPlot[TvsLfrac, Joined -> True,PlotLabel -> \"Graphic H: Time x Fracture Length\",AxesLabel -> {\"time (s)\", \"Fracture length (m)\"},Filling -> Axis, FillingStyle -> Green,PlotRange -> {{0, lastTime}, {0, maxlfrac}}, PlotMarkers -> Automatic]\n\n";
  
  outf << "maxOp = TvsOp[[-1,2]];\n";
  outf << "GrI = ListPlot[TvsOp, Joined -> True,PlotLabel -> \"Graphic I: Time x Opening\",AxesLabel -> {\"time (s)\", \"Opening (m)\"},Filling -> Axis, FillingStyle -> Green,PlotRange -> {{0, lastTime}, {0, maxOp}}, PlotMarkers -> Automatic]\n\n";
  outf << std::endl;
}

InputDataStruct globFractInputData;

OutputDataStruct globFractOutputData;

