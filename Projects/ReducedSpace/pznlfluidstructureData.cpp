//
//  pznlfluidstructureData.cpp
//  PZ
//
//  Created by Cesar Lucci on 07/08/13.
//
//

#include "pznlfluidstructureData.h"
#include "pzreal.h"

InputDataStruct::InputDataStruct()
{
    
}

InputDataStruct::~InputDataStruct()
{
    
}

void InputDataStruct::SetData(REAL Lx, REAL Ly, REAL Lf, REAL Hf, REAL E, REAL Poisson, REAL Fx, REAL Fy,
                              REAL preStressXX, REAL preStressXY, REAL preStressYY,
                              int NStripes, REAL Visc, REAL SigN, REAL QinjTot, REAL Ttot, REAL maxDeltaT, int nTimes,
                              REAL Cl, REAL Pe, REAL SigmaConf, REAL Pref, REAL vsp, REAL KIc, REAL Jradius)
{
    fLx = Lx;
    fLy = Ly;
    fLf = Lf;
    fHf = Hf;
    
    fE = E;
    fPoisson = Poisson;
    fFx = Fx;
    fFy = Fy;
    fPreStressXX = preStressXX;
    fPreStressXY = preStressXY;
    fPreStressYY = preStressYY;
    fNStripes = NStripes;
    
    fVisc = Visc;
    
    fSigN = SigN;
    
    REAL Qinj1asa = QinjTot / 2.;
    REAL QinjSecao = Qinj1asa / Hf;
    fQinj = QinjSecao;
    
    fTtot = Ttot;
    factTime = fminDeltaT;
    fmaxDeltaT = maxDeltaT;
    fNDeltaTsteps = nTimes;
    fminDeltaT = fmaxDeltaT/fNDeltaTsteps;
    factDeltaT = fminDeltaT;
    
    fCl = Cl;
    fPe = Pe;
    fSigmaConf = SigmaConf;
    fPref = Pref;
    fvsp = vsp;
    
    fKIc = KIc;
    fJradius = Jradius;
}

void InputDataStruct::SetLf(REAL Lf)
{
    fLf = Lf;
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

REAL InputDataStruct::Hf()
{
    return fHf;
}

REAL InputDataStruct::E()
{
    return fE;
}

REAL InputDataStruct::Poisson()
{
    return fPoisson;
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

int InputDataStruct::NStripes()
{
    return fNStripes;
}

REAL InputDataStruct::Visc()
{
    return fVisc;
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
}

void InputDataStruct::SetNextDeltaT()
{
    factDeltaT = std::min(fmaxDeltaT,factDeltaT+fmaxDeltaT/fNDeltaTsteps);
}

void InputDataStruct::UpdateActTime()
{
    factTime += factDeltaT;
    std::cout << "\n\n=============== ActTime = " << factTime << " ===============\n\n";
}



//------------------------------------------------------------

OutputDataStruct::OutputDataStruct()
{
    fTposP.clear();
    fTposVolLeakoff.clear();
    fTAcumVolW.clear();
    fTAcumVolLeakoff.clear();
    fTKI.clear();
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
    fQinj1wing = 0.;
    fLfracMax = 0.;
}

int OutputDataStruct::NTimes()
{
    int ntimes0 = fTposP.size();
    
#ifdef DEBUG
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
    std::string plotfile = "TransientSolution.vtk";
    if(anCount >= 0)
    {
        an->SetStep(anCount);
    }
    
    TPZManVector<std::string,10> scalnames(5), vecnames(1);
	
	scalnames[0] = "DisplacementX";
	scalnames[1] = "DisplacementY";
    scalnames[2] = "SigmaX";
    scalnames[3] = "SigmaY";
    scalnames[4] = "Pressure";
    vecnames[0] = "Displacement";
	
	const int dim = 2;
	int div =0;
	an->DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an->PostProcess(div);
}

void OutputDataStruct::PrintMathematica(std::ofstream & outf)
{
#ifdef DEBUG
    if(fTposP.size() == 0 || fTposVolLeakoff.size() == 0 || fTAcumVolW.size() == 0 || fTAcumVolLeakoff.size() == 0 || fTKI.size() == 0)
    {
        DebugStop();
    }
#endif
    
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
    
    outf << "Caju2013;\n\n";
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
    
    outf << "(* time x W Volume *)\n";
    outf << "TvsVolW={";
    for(itTAcumVolW = fTAcumVolW.begin(); itTAcumVolW != fTAcumVolW.end(); itTAcumVolW++)
    {
        outf << "{" << itTAcumVolW->first << "," << itTAcumVolW->second << "}";
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
        outf << "{" << itTAcumVolLeakoff->first << "," << itTAcumVolLeakoff->second << "}";
        if(itTAcumVolLeakoff != itTAcumVolLeakoffLast)
        {
            outf << ",";
        }
    }
    outf << "};\n\n";
    
    outf << "(* time x KI *)\n";
    outf << "TvsKI={";
    for(itTKI = fTKI.begin(); itTKI != fTKI.end(); itTKI++)
    {
        outf << "{" << itTKI->first << "," << itTKI->second << "}";
        if(itTKI != itTKILast)
        {
            outf << ",";
        }
    }
    outf << "};\n\n";
    
    outf << "(* Qinj 1 wing and Lfrac max *)\n";
    outf << "Qinj1wing=" << fQinj1wing << ";\n";
    outf << "LfracMax=" << fLfracMax << ";\n\n";
    
    outf << "maxp = Max[Transpose[Flatten[posvsPvsT, 1]][[2]]];\n";
    outf << "Manipulate[ListPlot[posvsPvsT[[t]], Joined -> True,PlotLabel ->\"Graphic A: Position x Pressure @ \" <> ToString[times[[t]]] <> \"s\",AxesLabel -> {\"position (m)\", \"pressure (Pa)\"}, Filling -> Axis,PlotRange -> {{0, LfracMax}, {0, maxp}}], {t, 1, ntimes, 1}]\n\n";
    
    outf << "maxleakoff = Max[Transpose[Flatten[posvsVolleakoffvsT, 1]][[2]]];\n";
    outf << "Manipulate[ListPlot[posvsVolleakoffvsT[[t]], Joined -> True,PlotLabel ->\"Graphic B: Position x Leakoff penetration @ \" <> ToString[times[[t]]] <>\"s\", AxesLabel -> {\"position (m)\", \"leakoff penetration (m)\"},Filling -> Axis, PlotRange -> {{0, LfracMax}, {0, maxleakoff}}], {t, 1,ntimes, 1}]\n\n";
    
    outf << "maxinj = Qinj1wing*times[[ntimes]];\n";
    outf << "GrC = Plot[Qinj1wing*t, {t, 0, times[[ntimes]]},PlotLabel -> \"Graphic C: Time x Volume Injected\",AxesLabel -> {\"time (s)\", \"Volume injected (m3)\"},Filling -> Axis, FillingStyle -> Red,PlotRange -> {{0, times[[ntimes]]}, {0, maxinj}}]\n\n";
    
    outf << "GrD = ListPlot[TvsVolW, Joined -> True,PlotLabel -> \"Graphic D: Time x Fracture Volume\",AxesLabel -> {\"time (s)\", \"Fracture volume (m3)\"},Filling -> Axis, FillingStyle -> Green,PlotRange -> {{0, times[[ntimes]]}, {0, maxinj}}]\n\n";
    
    outf << "GrE = ListPlot[TvsVolLeakoff, Joined -> True,PlotLabel -> \"Graphic E: Time x Leakoff volume\",AxesLabel -> {\"time (s)\", \"Leakoff volume (m3)\"},Filling -> Axis, FillingStyle -> Blue,PlotRange -> {{0, times[[ntimes]]}, {0, maxinj}}]\n\n";
    
    outf << "WplusLeakoff = {};\n";
    outf << "For[tt = 1, tt <= ntimes,\n";
    outf << "AppendTo[WplusLeakoff, {times[[tt]], TvsVolW[[tt, 2]] + TvsVolLeakoff[[tt, 2]]}];\n";
    outf << "tt++;\n";
    outf << "];\n";
    outf << "GrF = ListPlot[WplusLeakoff, Joined -> False,PlotStyle -> {Black, PointSize[0.03]},PlotLabel -> \"Graphic F: Grahics (D+E)\",AxesLabel -> {\"time (s)\", \"Vol graphics(D+E)\"},PlotRange -> {{0, times[[ntimes]] + 1}, {0, maxinj + 1}}];\n";
    outf << "Show[GrC, GrD, GrE, GrF]\n\n";
    
    outf << "maxki = Max[Transpose[TvsKI][[2]]];\n";
    outf << "ListPlot[TvsKI, Joined -> True, PlotLabel -> \"Graphic G: Time x KI\", AxesLabel -> {\"time (s)\", \"KI (Pa.m2)\"},Filling -> Axis,PlotRange -> {{0, times[[ntimes]]}, {0, maxki}}]\n";
}

InputDataStruct globFractInputData;

OutputDataStruct globFractOutputData;

