//
//  pznlfluidstructureData.cpp
//  PZ
//
//  Created by Cesar Lucci on 07/08/13.
//
//

#include "TPZPlaneFractureData.h"
#include "TPZPlaneFractCouplingMat.h"
#include "pzreal.h"
#include "pzcompel.h"
#include "pzmultiphysicscompel.h"
#include "pzintel.h"
#include "pzbndcond.h"
#include "pzgeotriangle.h"
#include "pzgeoquad.h"


void LeakoffStorage::UpdateLeakoff(TPZCompMesh * cmesh, REAL deltaT)
{
    if(this->fDefaultLeakoffEnabled == false)
    {
        return;
    }
    
    if(this->fPressureIndependent == false)
    {
        if(fGelId_Penetration.size() == 0)
        {
            std::cout << "\n\n\nLeakoff map is empty!!!\n\n\n";
            //Estah definido Noleakoff ? Se sim, aqui deve entao ser mudado para return (e nao DebugStop)
            DebugStop();
        }
    }
    int outVlCount = 0;
    for(int i = 0;  i < cmesh->ElementVec().NElements(); i++)
    {
        TPZCompEl * cel = cmesh->ElementVec()[i];
        
        if(cel->Reference()->Dimension() != 2 || globMaterialIdGen.IsInsideFractMat(cel->Reference()->MaterialId()) == false)
        {
            continue;
        }
        std::map<int,REAL>::iterator it = fGelId_Penetration.find(cel->Reference()->Id());
        if(it == fGelId_Penetration.end())
        {
            if(this->fPressureIndependent == false)
            {
                std::cout << "\n\n\nElemento de Id = " << cel->Reference()->Id() << " e Index = " << cel->Reference()->Index() <<
                             " nao encontrado no UpdateLeakoff\n";
                std::cout << "Seria o TransferLeakoff anterior que nao o incluiu???\n";
                std::cout << "Ver método " << __PRETTY_FUNCTION__ << "\n\n\n";
                DebugStop();
            }
            else
            {
                fGelId_Penetration[cel->Reference()->Id()] = 0.;
                it = fGelId_Penetration.find(cel->Reference()->Id());
            }
        }

        TPZInterpolatedElement * sp = NULL;
        if(cel->Reference()->Type() == ETriangle)
        {
            TPZMultiphysicsCompEl<pzgeom::TPZGeoTriangle> * celmp = dynamic_cast<TPZMultiphysicsCompEl< pzgeom::TPZGeoTriangle> * >(cel);
            sp = dynamic_cast <TPZInterpolatedElement*> (celmp->ElementVec()[1]);
        }
        else
        {
            TPZMultiphysicsCompEl<pzgeom::TPZGeoQuad> * celmp = dynamic_cast<TPZMultiphysicsCompEl< pzgeom::TPZGeoQuad> * >(cel);
            sp = dynamic_cast <TPZInterpolatedElement*> (celmp->ElementVec()[1]);
        }
        
#ifdef DEBUG
        if(!sp)
        {
            DebugStop();
        }
#endif
        
        TPZVec<REAL> qsi(2,0.);
        cel->Reference()->CenterPoint(cel->Reference()->NSides()-1, qsi);
        TPZMaterialData data;
        sp->InitMaterialData(data);
        sp->ComputeShape(qsi, data);
        sp->ComputeSolution(qsi, data);

        REAL pfrac = data.sol[0][0];
        
        TPZBndCond * matbnd = dynamic_cast<TPZBndCond*> (cel->Material());

#ifdef DEBUG
        if(!matbnd)
        {
            DebugStop();
        }
#endif

        TPZPlaneFractCouplingMat * mat = dynamic_cast<TPZPlaneFractCouplingMat *>(matbnd->Material());
        
#ifdef DEBUG
        if(!mat)
        {
            DebugStop();
        }
#endif
        
        REAL VlAcum = it->second;
        REAL Cl = mat->Cl();
        REAL Pe = mat->Pe();
        REAL gradPref = mat->gradPref();
        REAL vsp = mat->vsp();
        
        REAL tStar = FictitiousTime(VlAcum, pfrac, Cl, Pe, gradPref, vsp);
        REAL Vlnext = VlFtau(pfrac, tStar + deltaT, Cl, Pe, gradPref, vsp);
        
        if(fLeakoffEnabled)
        {
            it->second = Vlnext;
        }
        else
        {
            it->second = 0.;
        }
        
        outVlCount++;
    }
    
#ifdef DEBUG
    if(outVlCount < fGelId_Penetration.size())
    {
        DebugStop();
    }
#endif
}

REAL LeakoffStorage::VlFtau(REAL pfrac, REAL tau, REAL Cl, REAL Pe, REAL gradPref, REAL vsp)
{
    REAL gradPcalc = 1.;
    if(fPressureIndependent == false)
    {
        REAL gradP = pfrac - Pe;
        gradPcalc = gradP/gradPref;
    }
    
    REAL Clcorr = Cl * sqrt(gradPcalc);
    REAL Vl = 2. * Clcorr * sqrt(tau) + vsp;
    
    return Vl;
}

REAL LeakoffStorage::FictitiousTime(REAL VlAcum, REAL pfrac, REAL Cl, REAL Pe, REAL gradPref, REAL vsp)
{
    REAL tStar = 0.;
    if(VlAcum > vsp)
    {
        REAL gradPcalc = 1.;

        if(fPressureIndependent == false)
        {
            REAL gradP = pfrac - Pe;
            gradPcalc = gradP/gradPref;
        }
        
        REAL Clcorr = Cl * sqrt(gradPcalc);
        tStar = (VlAcum - vsp)*(VlAcum - vsp)/( (2. * Clcorr) * (2. * Clcorr) );
    }
    
    return tStar;
}

REAL LeakoffStorage::QlFVl(int gelId, REAL pfrac, REAL deltaT, REAL Cl, REAL Pe, REAL gradPref, REAL vsp)
{
    if(pfrac <= Pe || fLeakoffEnabled == false)
    {
        return 0.;
    }
    
    std::map<int,REAL>::iterator it = fGelId_Penetration.find(gelId);
    if(it == fGelId_Penetration.end())
    {
        fGelId_Penetration[gelId] = 0.;//Nao coloque vsp! Eh ZERO mesmo!
        it = fGelId_Penetration.find(gelId);
    }
    
    REAL VlAcum = it->second;
    
    REAL tStar = FictitiousTime(VlAcum, pfrac, Cl, Pe, gradPref, vsp);
    REAL Vlnext = VlFtau(pfrac, tStar + deltaT, Cl, Pe, gradPref, vsp);
    REAL Ql = (Vlnext - VlAcum)/deltaT;
    
    return Ql;
}

REAL LeakoffStorage::dQlFVl(int gelId, REAL pfrac, REAL deltaT, REAL Cl, REAL Pe, REAL gradPref, REAL vsp)
{
    if(pfrac <= Pe || fLeakoffEnabled == false)
    {
        return 0.;
    }
    
    std::map<int,REAL>::iterator it = fGelId_Penetration.find(gelId);
    if(it == fGelId_Penetration.end())
    {
        fGelId_Penetration[gelId] = 0.;
        it = fGelId_Penetration.find(gelId);
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
    
    /////////////////////////////////////////////////Ql maior
    REAL pfracUP = pfrac + deltaPfrac;
    REAL tStar1 = FictitiousTime(VlAcum, pfracUP, Cl, Pe, gradPref, vsp);
    REAL Vlnext1 = VlFtau(pfracUP, tStar1 + deltaT, Cl, Pe, gradPref, vsp);
    REAL Ql1 = (Vlnext1 - VlAcum )/deltaT;
    /////////////////////////////////////////////////Ql menor
    REAL pfracDOWN = pfrac - deltaPfrac;
    REAL tStar0 = FictitiousTime(VlAcum, pfracDOWN, Cl, Pe, gradPref, vsp);
    REAL Vlnext0 = VlFtau(pfracDOWN, tStar0 + deltaT, Cl, Pe, gradPref, vsp);
    REAL Ql0 = (Vlnext0 - VlAcum)/deltaT;
    /////////////////////////////////////////////////
    
    REAL dQldpfrac = (Ql1-Ql0)/(2.*deltaPfrac);
    
    return dQldpfrac;
}

void LeakoffStorage::Printleakoff(std::ofstream & outf)
{
    std::map<int,REAL>::iterator it;
    
    for(it = fGelId_Penetration.begin(); it != fGelId_Penetration.end(); it++)
    {
        outf << "Id = " << it->first << " : Penetration = " << it->second << "\n";
    }
}


//------------------------------------------------------------


ElastReducedSolution::ElastReducedSolution()
{
    this->fElastReducedSolution.Resize(0,0);
}

ElastReducedSolution::~ElastReducedSolution()
{
    this->fElastReducedSolution.Resize(0,0);
}

void ElastReducedSolution::SetElastReducedSolution(TPZFMatrix<REAL> & ElastReducedSolution)
{
    this->fElastReducedSolution = ElastReducedSolution;
}

REAL ElastReducedSolution::GetTotalPressure(int stripe)
{
    //Como agora eh aplicado newman unitario, o alpha corresponde aa pressao aplicada!
    REAL pressAppliedEntireFract = this->fElastReducedSolution(1,0);
    
    int stripePressAppliedRow = this->GetStressAppliedSolutionRow(stripe);
    REAL pressAppliedStripe = this->fElastReducedSolution(stripePressAppliedRow,0);
    
    REAL pressApplied = pressAppliedEntireFract + pressAppliedStripe;
    
    return pressApplied;
}

REAL ElastReducedSolution::GetNetPressure(int layer, int stripe)
{
    REAL totalPressApplied = this->GetTotalPressure(stripe);
    REAL preStress = -globLayerStruct.GetLayer(layer).fSigYY;

    /** NAO COLOQUE MAX(0.,NETPRESSURE)!!! JA VI QUE NAO DA CERTO!!! */
    REAL netPressure = totalPressApplied - preStress;

    return netPressure;
}

int ElastReducedSolution::GetStressAppliedSolutionRow(int stripe)
{
    return stripe+2;
}

//------------------------------------------------------------

Output3DDataStruct::Output3DDataStruct()
{
    fQinj1wing = 0.;
    fTAcumVolW.clear();
    fTAcumVolLeakoff.clear();
    fTL.clear();
    fTHsup.clear();
    fTHinf.clear();
    fFractContour.clear();
    
    InsertTAcumVolW(0.,0.);
    InsertTAcumVolLeakoff(0,0.);
    InsertTL(0.,0.);
    InsertTHsup(0.,0.);
    InsertTHinf(0.,0.);
}

void Output3DDataStruct::SetQinj1wing(REAL Qinj1wing)
{
    fQinj1wing = Qinj1wing;
}

Output3DDataStruct::~Output3DDataStruct()
{
    fTAcumVolW.clear();
    fTAcumVolLeakoff.clear();
    fTL.clear();
    fTHsup.clear();
    fTHinf.clear();
}

int Output3DDataStruct::NTimes()
{
    int ntimes0 = fTAcumVolW.size();
    return ntimes0;
}

void Output3DDataStruct::InsertTAcumVolW(REAL time, REAL vol)
{
    fTAcumVolW[time] = vol;
}

void Output3DDataStruct::InsertTAcumVolLeakoff(REAL time, REAL vol)
{
    fTAcumVolLeakoff[time] = vol;
}

void Output3DDataStruct::InsertTL(REAL time, REAL L)
{
    fTL[time] = L;
}

void Output3DDataStruct::InsertTHsup(REAL time, REAL Hsup)
{
    fTHsup[time] = Hsup;
}

void Output3DDataStruct::InsertTHinf(REAL time, REAL Hinf)
{
    fTHinf[time] = Hinf;
}

void Output3DDataStruct::PrintMathematica(std::ofstream & outf)
{
#ifdef DEBUG
    if(fTAcumVolW.size() == 0 || fTAcumVolLeakoff.size() == 0)
    {
        DebugStop();
    }
#endif
    
    std::map<REAL,REAL>::iterator itTAcumVolW, itTAcumVolWLast = fTAcumVolW.end();
    itTAcumVolWLast--;
    
    outf << "(* Output Fracture Propagation 1D *)\n";
    outf << "Caju2013;\n\n";
    outf << "ntimes=" << NTimes() << ";\n";
    outf << "times={";
    for(itTAcumVolW = fTAcumVolW.begin(); itTAcumVolW != fTAcumVolW.end(); itTAcumVolW++)
    {
        outf << itTAcumVolW->first;
        if(itTAcumVolW != itTAcumVolWLast)
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
    
    
    std::map<REAL,REAL>::iterator itTAcumVolLeakoff, itTAcumVolLeakoffLast = fTAcumVolLeakoff.end();
    itTAcumVolLeakoffLast--;
    
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
    
    outf << "(* Qinj 1 wing *)\n";
    outf << "Qinj1wing=" << fQinj1wing << ";\n";
    
    outf << "maxinj = Qinj1wing*times[[ntimes]];\n";
    outf << "GrD = Plot[Qinj1wing*t, {t, 0, times[[ntimes]]},PlotLabel -> \"Graphic D: Time x Volume Injected\",AxesLabel -> {\"time (s)\", \"Volume injected (m3)\"},Filling -> Axis, FillingStyle -> Red,PlotRange -> {{0, times[[ntimes]]}, {0, maxinj}}]\n\n";
    
    outf << "GrE = ListPlot[TvsVolW, Joined -> True,PlotLabel -> \"Graphic E: Time x Fracture Volume\",AxesLabel -> {\"time (s)\", \"Fracture volume (m3)\"},Filling -> Axis, FillingStyle -> Green,PlotRange -> {{0, times[[ntimes]]}, {0, maxinj}}]\n\n";
    
    outf << "GrF = ListPlot[TvsVolLeakoff, Joined -> True,PlotLabel -> \"Graphic F: Time x Leakoff volume\",AxesLabel -> {\"time (s)\", \"Leakoff volume (m3)\"},Filling -> Axis, FillingStyle -> Blue,PlotRange -> {{0, times[[ntimes]]}, {0, maxinj}}]\n\n";
    
    outf << "WplusLeakoff = {};\n";
    outf << "For[tt = 1, tt <= ntimes,\n";
    outf << "AppendTo[WplusLeakoff, {times[[tt]], TvsVolW[[tt, 2]] + TvsVolLeakoff[[tt, 2]]}];\n";
    outf << "tt++;\n";
    outf << "];\n";
    outf << "GrG = ListPlot[WplusLeakoff, Joined -> False,PlotStyle -> {Black, PointSize[0.03]},PlotLabel -> \"Graphic G: Grahics (E+F)\",AxesLabel -> {\"time (s)\", \"Vol graphics(D+E)\"},PlotRange -> {{0, times[[ntimes]] + 1}, {0, maxinj + 1}}];\n";
    outf << "Show[GrD, GrE, GrF, GrG, PlotLabel -> \"Graphic G: Grahics D, E, F and (E+F)\"]\n\n\n\n";
    
    //--------------------------
    
    std::map<REAL,REAL>::iterator itT, itTlast = fTL.end();
    itTlast--;
    
    outf << "(* time x Max Lfrac *)\n";
    outf << "TvsLfracmax={";
    for(itT = fTL.begin(); itT != fTL.end(); itT++)
    {
        outf << "{" << itT->first << "," << itT->second << "}";
        if(itT != itTlast)
        {
            outf << ",";
        }
    }
    outf << "};\n\n";
    
    //--------------------------
    
    itTlast = fTHsup.end();
    itTlast--;
    
    outf << "(* time x Hsup *)\n";
    outf << "TvsHsup={";
    for(itT = fTHsup.begin(); itT != fTHsup.end(); itT++)
    {
        outf << "{" << itT->first << "," << itT->second << "}";
        if(itT != itTlast)
        {
            outf << ",";
        }
    }
    outf << "};\n\n";
    
    //--------------------------
    
    itTlast = fTHinf.end();
    itTlast--;
    
    outf << "(* time x Hinf *)\n";
    outf << "TvsHinf={";
    for(itT = fTHinf.begin(); itT != fTHinf.end(); itT++)
    {
        outf << "{" << itT->first << "," << itT->second << "}";
        if(itT != itTlast)
        {
            outf << ",";
        }
    }
    outf << "};\n\n";
    
    //--------------------------
    
    
}

TimeControl globTimeControl;

LeakoffStorage globLeakoffStorage;

MaterialIdGen globMaterialIdGen;

Output3DDataStruct globFractOutput3DData;

ElastReducedSolution globElastReducedSolution;

LayerStruct globLayerStruct;
