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
            sp = dynamic_cast <TPZInterpolatedElement*> (celmp->ElementVec()[1].Element());
        }
        else
        {
            TPZMultiphysicsCompEl<pzgeom::TPZGeoQuad> * celmp = dynamic_cast<TPZMultiphysicsCompEl< pzgeom::TPZGeoQuad> * >(cel);
            sp = dynamic_cast <TPZInterpolatedElement*> (celmp->ElementVec()[1].Element());
        }
        
#ifdef PZDEBUG
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
        
#ifdef PZDEBUG
        if(!matbnd)
        {
            DebugStop();
        }
#endif
        
        TPZPlaneFractCouplingMat * mat = dynamic_cast<TPZPlaneFractCouplingMat *>(matbnd->Material());
        
#ifdef PZDEBUG
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
    
#ifdef PZDEBUG
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
        if(gradPcalc < 0.)
        {
            return 0.;
        }
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
            if(gradPcalc <= 0.)
            {
                return 0.;
            }
        }
        
        REAL Clcorr = Cl * sqrt(gradPcalc);
        tStar = (VlAcum - vsp)*(VlAcum - vsp)/( (2. * Clcorr) * (2. * Clcorr) );
    }
    
    return tStar;
}

REAL LeakoffStorage::QlFVl(int gelId, REAL pfrac, REAL deltaT, REAL Cl, REAL Pe, REAL gradPref, REAL vsp)
{
    if(fLeakoffEnabled == false)
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
    if(fLeakoffEnabled == false)
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


//------------------------------------------------------------

//Inicializando vetor de cores
const std::string Output3DDataStruct::color[12] = {"Red","Green","Blue","Black","Gray","Cyan","Magenta","Yellow","Brown","Orange","Pink","Purple"};

Output3DDataStruct::Output3DDataStruct() : actColor(0)
{
    fQinj1wing = 0.;
    fTAcumVolW.clear();
    fTmeanW.clear();
    fTAcumVolLeakoff.clear();
    fFractContour.clear();
    fTNetPressure.clear();
    
    InsertTAcumVolW(0.,0.);
    InsertTAcumVolLeakoff(0,0.);
}

void Output3DDataStruct::SetQinj1wing(REAL Qinj1wing)
{
    fQinj1wing = Qinj1wing;
}

Output3DDataStruct::~Output3DDataStruct()
{
    fTAcumVolW.clear();
    fTmeanW.clear();
    fTAcumVolLeakoff.clear();
    fTNetPressure.clear();
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

void Output3DDataStruct::InsertTmeanW(REAL time, REAL meanW)
{
    fTmeanW[time] = meanW;
}

void Output3DDataStruct::InsertTAcumVolLeakoff(REAL time, REAL vol)
{
    fTAcumVolLeakoff[time] = vol;
}

void Output3DDataStruct::InsertTNetPressure(REAL time, REAL netpressure)
{
    fTNetPressure[time] = netpressure;
}

void Output3DDataStruct::PrintConservationMass()
{
    std::ofstream outf("000ConservationMass.txt");
    
    std::map<REAL,REAL>::iterator itTAcumVolW, itTAcumVolWLast = fTAcumVolW.end();
    itTAcumVolWLast--;
    
    outf << "(* Output Fracture Propagation 1D *)\n\n";
    
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
    outf << "Qinj1wing=" << fQinj1wing << "*60;\n";
    
    outf << "maxinj = Qinj1wing*times[[ntimes]];\n";
    outf << "GrD = Plot[Qinj1wing*t, {t, 0, times[[ntimes]]},PlotLabel -> \"Graphic D: Time x Volume Injected\",AxesLabel -> {\"time (min)\", \"Volume injected (m3)\"},Filling -> Axis, FillingStyle -> Red,PlotRange -> {{0, times[[ntimes]]}, {0, maxinj}}]\n\n";
    
    outf << "GrE = ListPlot[TvsVolW, Joined -> True,PlotLabel -> \"Graphic E: Time x Fracture Volume\",AxesLabel -> {\"time (min)\", \"Fracture volume (m3)\"},Filling -> Axis, FillingStyle -> Green,PlotRange -> {{0, times[[ntimes]]}, {0, maxinj}}]\n\n";
    
    outf << "GrF = ListPlot[TvsVolLeakoff, Joined -> True,PlotLabel -> \"Graphic F: Time x Leakoff volume\",AxesLabel -> {\"time (min)\", \"Leakoff volume (m3)\"},Filling -> Axis, FillingStyle -> Blue,PlotRange -> {{0, times[[ntimes]]}, {0, maxinj}}]\n\n";
    
    outf << "WplusLeakoff = {};\n";
    outf << "For[tt = 1, tt <= ntimes,\n";
    outf << "AppendTo[WplusLeakoff, {times[[tt]], TvsVolW[[tt, 2]] + TvsVolLeakoff[[tt, 2]]}];\n";
    outf << "tt++;\n";
    outf << "];\n";
    outf << "GrG = ListPlot[WplusLeakoff, Joined -> False,PlotStyle -> {Black, PointSize[0.03]},PlotLabel -> \"Graphic G: Grahics (E+F)\",AxesLabel -> {\"time (min)\", \"Vol graphics(D+E)\"},PlotRange -> {{0, times[[ntimes]] + 1}, {0, maxinj + 1}}];\n";
    outf << "Show[GrD, GrF, GrE, GrG, PlotLabel -> \"Graphic G: Grahics D, E, F and (E+F)\"]\n\n\n\n";
    
    
    std::map<REAL,REAL>::iterator itTNetPress, itTNetPressLast = fTNetPressure.end();
    itTNetPressLast--;
    
    outf << "(* time x NetPressure *)\n";
    outf << "TvsNetPressure={";
    for(itTNetPress = fTNetPressure.begin(); itTNetPress != fTNetPressure.end(); itTNetPress++)
    {
        outf << "{" << itTNetPress->first << "," << itTNetPress->second << "}";
        if(itTNetPress != itTNetPressLast)
        {
            outf << ",";
        }
    }
    outf << "};\n\n";
    outf << "ListPlot[TvsNetPressure, Joined -> True,PlotLabel -> \"Graphic H: Time x Net Pressure\",AxesLabel -> {\"time (min)\", \"Net pressure (MPa)\"},Filling -> Axis,FillingStyle -> Opacity[0.3,Cyan], AxesOrigin -> {0,0},PlotRange -> All]\n\n";
}

void Output3DDataStruct::PrintFractureGeometry(int num,
                                               TPZVec< std::pair<REAL,REAL> > & poligonalChain,
                                               REAL CenterTVD)
{
    if(num < 0)
    {
        std::cout << "\n\n\nnum<0. See " << __PRETTY_FUNCTION__ << ".\n\n\n";
        DebugStop();
    }
    
    std::ofstream outF("000FractContours.txt");
    
    {   //Preamble for PostProcessFractGeometry method
        outF << "(* colors = {Red,Green,Blue,Black,Gray,Cyan,Magenta,Yellow,Brown,Orange,Pink,Purple} *)\n";
        outF << "AllPolChains={};\n";
        outF << "Lgr={{0,0}};\n";
        outF << "Hsupgr={{0,0}};\n";
        outF << "Hinfgr={{0,0}};\n\n";
    }
    
    std::stringstream nmMath, nmAux;
    nmAux << "pcm={";
    
    for(int p = 0; p < poligonalChain.NElements(); p++)
    {
        globFractOutput3DData.fFractContour << "fractureDots" << p << " = {" << poligonalChain[p].first << ","
        << poligonalChain[p].second + CenterTVD << "};\n";
        nmAux << "fractureDots" << p;
        if(p < poligonalChain.NElements()-1)
        {
            nmAux << ",";
        }
    }
    nmAux << "};\n";
    nmAux << "gr" << num << "=ListPlot[pcm,Joined->True,AxesOrigin->{0,0},AspectRatio->1,PlotStyle->" << color[actColor%12] << ",AxesLabel->{\"L (m)\", \"H (m)\"}];\n";
    nmAux << "L" << num << "=Max[Transpose[pcm][[1]]];\n";
    nmAux << "Hsup" << num << "=Max[Transpose[pcm][[2]]];\n";
    nmAux << "Hinf" << num << "=-Min[Transpose[pcm][[2]]];\n";
    nmAux << "AppendTo[AllPolChains,gr" << num << "];\n";
    nmAux << "AppendTo[Lgr,{" << globTimeControl.actTime()/60. << ",L" << num << "}];\n";
    nmAux << "AppendTo[Hsupgr,{" << globTimeControl.actTime()/60. << ",Hsup" << num << "}];\n";
    nmAux << "AppendTo[Hinfgr,{" << globTimeControl.actTime()/60. << ",Hinf" << num << "}];\n";
    nmAux << "Print[\"time" << num << " = " << globTimeControl.actTime()/60. << " min\"]\n";
    nmAux << "Print[\"L" << num << " = \", L" << num << "]\n";
    nmAux << "Print[\"Hsup" << num << " = \", Hsup" << num << "]\n";
    nmAux << "Print[\"Hinf" << num << " = \", Hinf" << num << "]\n";
    nmAux << "Print[\"\"]\n\n";
    globFractOutput3DData.fFractContour << nmAux.str();
    actColor++;
    
    outF << globFractOutput3DData.fFractContour.str();
    
    outF << "grRange=Max[Max[Transpose[Lgr][[2]]],Max[Transpose[Hsupgr][[2]]],Max[Transpose[Hinfgr][[2]]]]+1;\n";
    outF << "Show[AllPolChains,PlotRange->{{0,2*grRange},{-grRange,grRange}}]\n";
    outF << "l={ListPlot[Lgr,AxesOrigin->{0,0},Filling->Axis,AxesLabel->{\"Time (min)\", \"L (green), Hsup (blue), Hinf (red) (m)\"}],";
    outF << "ListPlot[Lgr,Joined->True,AxesOrigin->{0,0},PlotStyle->Green]};\n";
    outF << "hs={ListPlot[Hsupgr,Filling->Axis,AxesOrigin->{0,0}],ListPlot[Hsupgr,Joined->True]};\n";
    outF << "hi={ListPlot[Hinfgr,PlotStyle->Red,Filling->Axis,AxesOrigin->{0,0}],ListPlot[Hinfgr,Joined->True,PlotStyle->Red]};\n";
    outF << "Show[l,hs,hi,PlotRange->All]\n\n";
    
    outF << "meanW = {{0,0},";
    std::map<REAL,REAL>::iterator itw, itwlast = fTmeanW.end();
    itwlast--;
    for(itw = fTmeanW.begin(); itw != fTmeanW.end(); itw++)
    {
        outF << "{" << itw->first << "," << itw->second << "}";
        if(itw != itwlast)
        {
            outF << ",";
        }
    }
    outF << "};\n";
    outF << "meanWgr = ListPlot[meanW, Joined -> True, PlotStyle -> Brown, Filling -> Axis,AxesLabel -> {\"Time (min)\", \"Wmed (mm)\"}];\n";
    outF << "Show[meanWgr]\n\n";
    
    outF.close();
}

TimeControl globTimeControl;

LeakoffStorage globLeakoffStorage;

MaterialIdGen globMaterialIdGen;

Output3DDataStruct globFractOutput3DData;

LayerStruct globLayerStruct;
