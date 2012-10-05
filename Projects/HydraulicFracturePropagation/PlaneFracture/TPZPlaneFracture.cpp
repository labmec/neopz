/**
 * \file
 * @brief Contains implementations of the TPZPlaneFracture methods.
 */
/*
 *  TPZPlaneFracture.cpp
 *  Crack
 *
 *  Created by Cesar Lucci on 09/08/10.
 *  Copyright 2010 LabMeC. All rights reserved.
 *
 */

#include "TPZPlaneFracture.h"

#include "pzgeopoint.h"/////MUSTDELETE
#include "TPZGeoLinear.h"
#include "TPZVTKGeoMesh.h"
#include "tpzgeoelrefpattern.h"
#include "tpzchangeel.h"
#include "pzanalysis.h"
#include "pzelast3d.h"
#include "pzelasmat.h"
#include "pzbndcond.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"

#include "adapt.h"
#include "TPZJIntegral.h"

using namespace pztopology;



/** PUBLIC METHODS */
TPZPlaneFracture::TPZPlaneFracture()
{
    std::cout << "Default constructor would not be used in this class!\n";
    DebugStop();
}
//------------------------------------------------------------------------------------------------------------

TPZPlaneFracture::TPZPlaneFracture(REAL lw, REAL bulletDepthIni, REAL bulletDepthFin, 
                                   TPZVec< std::map<REAL,REAL> > & pos_stress)
{
    fpos_stress = pos_stress;
    
    fPlaneMesh = new TPZGeoMesh;
    fFullMesh = new TPZGeoMesh;
        
    std::set<REAL> espacamentoVerticalTVD;
    std::list<REAL> espacamentoVerticalDEPTH;
    
    std::map<REAL,REAL>::iterator itM;
    
    espacamentoVerticalTVD.insert(0.);
    espacamentoVerticalTVD.insert(lw);
    espacamentoVerticalTVD.insert(bulletDepthIni);
    espacamentoVerticalTVD.insert(bulletDepthFin);
    
    int nstretches = pos_stress.NElements();
    for(int s = 0; s < nstretches; s++)
    {
        for(itM = pos_stress[s].begin(); itM != pos_stress[s].end(); itM++)
        {
            REAL pos = itM->first;
            espacamentoVerticalTVD.insert(pos);
        }
    }
    
    REAL pos0 = 0.;
    std::set<REAL>::iterator itS = espacamentoVerticalTVD.begin(); itS++;
    for(; itS != espacamentoVerticalTVD.end(); itS++)
    {
        REAL pos1 = *itS;
        REAL deltaZ = fabs(pos1 - pos0);
        
        int nrows = 1;
        REAL deltaZused = deltaZ/nrows;
        while(deltaZused > __maxLength)
        {
            nrows++;
            deltaZused = deltaZ/nrows;
        }
        
        for(int r = 1; r <= nrows; r++)
        {
            REAL z = pos0 + r*deltaZused;
            espacamentoVerticalTVD.insert(z);
        }
        pos0 = pos1;
    }
    
    for(itS = espacamentoVerticalTVD.begin(); itS != espacamentoVerticalTVD.end(); itS++)
    {
        REAL posDepth = *itS;
        espacamentoVerticalDEPTH.push_back(-posDepth);//Converting positions (TVD) in depth.
    }
    
    GeneratePlaneMesh(espacamentoVerticalDEPTH);
    GenerateFullMesh(espacamentoVerticalDEPTH);
}
//------------------------------------------------------------------------------------------------------------

TPZPlaneFracture::~TPZPlaneFracture()
{
    delete fPlaneMesh;
    delete fFullMesh;
    fpos_stress.Resize(0);
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFracture::RunThisFractureGeometry(const TPZVec<REAL> &poligonalChain, std::string vtkFile)
{
    ////CompMesh
    int porder = 2;
    TPZCompMesh * fractureCMesh = this->GetFractureCompMesh(poligonalChain, porder);
    
    int neq = fractureCMesh->NEquations();
    std::cout << "Numero de equacoes = " << neq << std::endl;
    
    ////Analysis
	TPZAnalysis an(fractureCMesh);
    
	TPZSkylineStructMatrix skylin(fractureCMesh); //caso simetrico
	TPZStepSolver<STATE> step;
	step.SetDirect(ECholesky);
    
    an.SetStructuralMatrix(skylin);
	an.SetSolver(step);
	an.Run();
    
    ////Post Processing
    TPZManVector<std::string,10> scalnames(6), vecnames(0);
    scalnames[0] = "DisplacementX";
    scalnames[1] = "DisplacementY";
    scalnames[2] = "DisplacementZ";
    scalnames[3] = "StressX";
    scalnames[4] = "StressY";
    scalnames[5] = "StressZ";
    
    int div = 0;
    const int dim = 3;
    an.DefineGraphMesh(dim,scalnames,vecnames,vtkFile);
    an.PostProcess(div,dim);
    
    /** 
     
     Post process possibilities:
     
     "Displacement"
     "state"
     "DisplacementX"
     "DisplacementY"
     "DisplacementZ"
     "PrincipalStress"
     "PrincipalStrain"
     "VonMises"
     "Stress"
     "Strain"
     "Stress1"
     "Strain1"
     "NormalStress"
     "NormalStrain"
     "StressX"
     "StressY"
     "StressZ"
     
     */
}
//------------------------------------------------------------------------------------------------------------

int TPZPlaneFracture::PointElementOnPlaneMesh(TPZGeoMesh * PlaneMesh, int & initial2DElId, TPZVec<REAL> & x, TPZVec<REAL> & qsi, int planeAxe0, int planeAxe1, int planeNormal, bool justFathers)
{
    TPZManVector<REAL,3> xproj(x);
    xproj[planeNormal] = 0.;
    
    TPZGeoEl * initialGel = PlaneMesh->ElementVec()[initial2DElId];
    qsi.Resize(2,0.);
    TPZVec<REAL> xqsi(3);
    
    #ifdef DEBUG
    if(initialGel->Dimension() != 2)
    {
        std::cout << "Given Id does not correspond to an 2D element on " << __PRETTY_FUNCTION__ << " !\n";
        DebugStop();
    }
    #endif
    
    TPZVec<REAL> dx(3);
    int count = 0;
    while(initialGel->ComputeXInverse(xproj, qsi) == false && count < PlaneMesh->NElements())
    {
        initialGel->CenterPoint(initialGel->NSides()-1, qsi);
        initialGel->X(qsi, xqsi);
        
        REAL norm = sqrt( (xproj[planeAxe0]-xqsi[planeAxe0])*(xproj[planeAxe0]-xqsi[planeAxe0]) +
                           (xproj[planeAxe1]-xqsi[planeAxe1])*(xproj[planeAxe1]-xqsi[planeAxe1]) );
        dx[planeAxe0] = (xproj[planeAxe0] - xqsi[planeAxe0])/norm;
        dx[planeNormal] = 0.;
        dx[planeAxe1] = (xproj[planeAxe1] - xqsi[planeAxe1])/norm;
        
        initialGel = CrossToNextNeighbour(initialGel, xqsi, dx, 1.E-5, planeAxe0, planeAxe1, planeNormal);
        
        count++;
    }
    
    #ifdef DEBUG
    if(count >= PlaneMesh->NElements())
    {
        std::cout.precision(15);
        std::cout << "Point ";
        std::cout << "{ " << xproj[0] << " , " << xproj[1] << " , " << xproj[2] << " } ";
        std::cout << " DO NOT belong to plane mesh domain on " << __PRETTY_FUNCTION__ << std::endl;
        DebugStop();
    }
    #endif
    
    if(justFathers == false && initialGel && initialGel->HasSubElement())
    {
        bool sonFound = false;
        TPZVec<TPZGeoEl*> sons(0);
        initialGel->GetLowerSubElements(sons);
        int nsons = sons.NElements();
        for(int s = 0; s < nsons; s++)
        {
            TPZGeoEl * actSon = sons[s];
            if(actSon->ComputeXInverse(xproj,qsi))
            {
                initialGel = actSon;
                sonFound = true;
                break;
            }
        }
        #ifdef DEBUG
        if(sonFound == false)
        {
            std::cout << "Given coordinates is in father domain, but was not on one of it sons!\nSee " << __PRETTY_FUNCTION__ << std::endl;
            std::cout << "Father Id = " << initialGel->Id() << " , Father nodes:\n";
            initialGel->PrintNodesCoordinates();
            std::cout << "Given x: ";
            std::cout.precision(15);
            std::cout << "{ " << xproj[0] << " , " << xproj[1] << " , " << xproj[2] << " }\n";
            std::cout << "Sons coordinates\n";
            for(int s = 0; s < nsons; s++)
            {
                TPZGeoEl * actSon = sons[s];
                actSon->PrintNodesCoordinates();
            }
            DebugStop();//Nao deveria chegar aqui pois um dos filhos deveria conter o ponto!
        }
        #endif
    }
	
    int id = initialGel->Id();
    
    initial2DElId = id;
    
	return id;
}


TPZGeoEl * TPZPlaneFracture::PointElementOnFullMesh(TPZVec<REAL> & x, TPZVec<REAL> & qsi, int & initial2DElId, TPZGeoMesh * fullMesh)
{
    int axe0 = 0;//axe X
    int axe1 = 2;//axe Z
    int axeNormal = 1;//axe Y
    int elFoundId = PointElementOnPlaneMesh(fullMesh, initial2DElId, x, qsi, axe0, axe1, axeNormal, true);
    
    // Da maneira com que esta classe foi construida, o elemento 2D encontrado
    // na malha fractMesh apresenta como seu dual (na malha fullMesh) o
    // elemento de mesmo id. Este dual serah o elemento de partida para a busca na direcao Y.
    TPZGeoEl * gelfullmesh = fullMesh->ElementVec()[elFoundId];
    
    TPZManVector<REAL,3> coord(x);
    while(gelfullmesh->Type() != ECube)
    {
        int side = gelfullmesh->NSides()-1;
        gelfullmesh = gelfullmesh->Neighbour(side).Element();
    }
    while(gelfullmesh->ComputeXInverse(coord, qsi) == false)
    {
        int cubeFace_in_Y_direction = 25;
        TPZGeoElSide cubeSide = gelfullmesh->Neighbour(cubeFace_in_Y_direction);
        while(cubeSide.Element()->Type() != ECube || cubeSide.Element()->Father())
        {
            cubeSide = cubeSide.Neighbour();
        }
        TPZGeoEl * nextGel = cubeSide.Element();
        if(nextGel == gelfullmesh)
        {
            std::cout << "End of domain reached without find searched element on " << __PRETTY_FUNCTION__ << std::endl;
            DebugStop();
        }
        gelfullmesh = nextGel;
        qsi.Resize(3,0.);
    }
    
    if(gelfullmesh->HasSubElement() == false)
    {
        return gelfullmesh;
    }
    else
    {
        TPZVec<TPZGeoEl*> subElements(0);
        gelfullmesh->GetLowerSubElements(subElements);
        int nCandidates = subElements.NElements();
        TPZVec< TPZVec<REAL> > subElQsi(nCandidates);
        
        for(int cd = 0; cd < nCandidates; cd++)
        {
            TPZGeoEl * cand = subElements[cd];
            if(cand->ComputeXInverse(coord, qsi) == true)
            {
                return cand;
            }
            subElQsi[cd] = qsi;
        }
        
        //if(!cand)
        {
            TPZGeoEl * cand = NULL;
            std::map<REAL,int> dist;
            for(int cd = 0; cd < nCandidates; cd++)
            {
                cand = subElements[cd];
                TPZVec<REAL> subElQsiProj(cand->Dimension());
                cand->ProjectInParametricDomain(subElQsi[cd], subElQsiProj);
                REAL distToProj = 0.;
                for(int c = 0; c < cand->Dimension(); c++)
                {
                    distToProj += (subElQsiProj[c] - subElQsi[cd][c]) * (subElQsiProj[c] - subElQsi[cd][c]);
                }
                dist[sqrt(distToProj)] = cd;
            }
            int subElPosition = dist.begin()->second;
            qsi = subElQsi[subElPosition];
            cand = subElements[subElPosition];
            
            #ifdef DEBUG
            REAL distComputed = dist.begin()->first;
            if(distComputed > 1.e-5)
            {
                DebugStop();
            }
            #endif
            
            return cand;
        }
    }
    
    DebugStop();//Nao encontrou nenhum elemento 3D.
    
    return NULL;
}


/** PRIVATE METHODS */
//------------------------------------------------------------------------------------------------------------

TPZGeoMesh * TPZPlaneFracture::GetFractureGeoMesh(const TPZVec<REAL> &poligonalChain)
{
#ifdef DEBUG
	int ncoord = poligonalChain.NElements();
	if(ncoord%2 != 0)
	{
		std::cout << "poligonalChain boundary dont have groups of 2 coordinates (x,z)!" << std::endl;
		std::cout << "See " << __PRETTY_FUNCTION__ << std::endl;
		DebugStop();
	}
#endif
	
    TPZGeoMesh * fullMesh = new TPZGeoMesh(*fFullMesh);
	TPZGeoMesh * planeMesh = new TPZGeoMesh(*fPlaneMesh);
    
    if(poligonalChain.NElements() < 4)//tem que ter no minimo 2 pontos (x1,z1,x2,z2)
    {
        delete planeMesh;
        
        return fullMesh;
    }
    
	int nelem = planeMesh->NElements();
    
	std::map< int, std::set<REAL> > elId_TrimCoords;
	std::list< std::pair<int,REAL> > elIdSequence;
	
	DetectEdgesCrossed(poligonalChain, planeMesh, elId_TrimCoords, elIdSequence);
	
	//Refining auxiliar 1D elements
	TPZVec<TPZGeoEl*> sons;
	std::map< int, std::set<REAL> >::iterator it;
	for(it = elId_TrimCoords.begin(); it != elId_TrimCoords.end(); it++)
	{
		int el1DId = it->first;
        TPZGeoEl * el1D = planeMesh->ElementVec()[el1DId];
        
		TPZAutoPointer<TPZRefPattern> linRefp = Generate1DRefPatt(it->second);
        
		el1D->SetRefPattern(linRefp);
		el1D->Divide(sons);
	}
    
	//Refining 2D and 3D elements with the intention to match the geometry of the crack boundary
	for(int el = 0; el < nelem; el++)
	{
		TPZGeoEl * gel = planeMesh->ElementVec()[el];//2D element in 2D mesh
        
		TPZAutoPointer<TPZRefPattern> elRefp = TPZRefPatternTools::PerfectMatchRefPattern(gel);
		if(elRefp)
		{
            gel = fullMesh->ElementVec()[el];//2D element in 3D mesh
            if(gel)
            {
                gel->SetRefPattern(elRefp);
                gel->Divide(sons);
                
                int innerSide = gel->NSides() - 1;
                gel = gel->Neighbour(innerSide).Element();//3D element in 3D mesh
                
                elRefp = TPZRefPatternTools::PerfectMatchRefPattern(gel);
                if(elRefp)
                {
                    gel->SetRefPattern(elRefp);
                    gel->Divide(sons);
                }
                int firstFace = 20;
                int lastFace = 25;
                for(int f = firstFace; f < lastFace; f++)//2D BC element in 3D mesh
                {
                    TPZGeoElSide hexaFace(gel,f);
                    TPZGeoElSide quadriFace = hexaFace.Neighbour();
                    if( quadriFace != hexaFace && IsBoundaryMaterial(quadriFace.Element()) )
                    {
                        elRefp = TPZRefPatternTools::PerfectMatchRefPattern(quadriFace.Element());
                        if(elRefp)
                        {
                            quadriFace.Element()->SetRefPattern(elRefp);
                            quadriFace.Element()->Divide(sons);
                        }
                    }
                }
            }
            else
            {
                DebugStop();
            }
		}
	}
	
	GenerateCrackBoundary(planeMesh, fullMesh, elIdSequence);
    SeparateElementsInMaterialSets(fullMesh);
    
    delete planeMesh;
    
#define REFINEDIRECTIONALOnce
#ifdef REFINEDIRECTIONALOnce
std::set<int> matids;
matids.insert(__1DcrackTipMat);
std::map< int,std::set<int> >::iterator itm;
for(itm = fcrackQpointsElementsIds.begin(); itm != fcrackQpointsElementsIds.end(); itm++)
{
    TPZGeoEl * gel = fullMesh->ElementVec()[itm->first];
    if(!gel->HasSubElement())
    {            
        TPZRefPatternTools::RefineDirectional(gel, matids);
    }
}
#endif

    ////4debug
//    TPZVec<REAL> qsiHere(3), qsiPZ(3);
//    TPZVec<REAL> x(3);
//    x[0] = 5.;
//    x[1] = 2.;
//    x[2] = -10.;
//    
//    int initialElId = 0;
//    
//    TPZGeoEl * gelHere = PointElementOnFullMesh(x, qsiHere, initialElId, fullMesh);
//    TPZGeoEl * gelPZ = fullMesh->FindElement(x, qsiPZ);
//
//    if(!gelPZ)
//    {
//        DebugStop();
//    }
//    TPZVec<REAL> xHere(3);
//    gelHere->X(qsiHere,xHere);
//    
//    TPZVec<REAL> xPZ(3);
//    gelPZ->X(qsiPZ, xPZ);
    //////////
    
	return fullMesh;
}
//------------------------------------------------------------------------------------------------------------

TPZCompMesh * TPZPlaneFracture::GetFractureCompMesh(const TPZVec<REAL> &poligonalChain, int porder)
{
    ////GeoMesh
    TPZGeoMesh * gmesh = this->GetFractureGeoMesh(poligonalChain);
    
    ////CompMesh
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    
    cmesh->SetDimModel(3);
    cmesh->SetAllCreateFunctionsContinuous();

    STATE young = 0.29E5;
    STATE poisson = 0.25;
    TPZVec<STATE> force(3,0.);

    TPZMaterial * materialLin = new TPZElasticity3D(__3DrockMat_linear, young, poisson, force);
    cmesh->InsertMaterialObject(materialLin); 
    
    TPZMaterial * materialQpoint = new TPZElasticity3D(__3DrockMat_quarterPoint, young, poisson, force);
    cmesh->InsertMaterialObject(materialQpoint);
    
    ////BCs    
    TPZFMatrix<STATE> k(3,3,0.), f(3,1,0.);
    int newmann = 1, mixed = 2;
    
    STATE pressureY = 5.;

    {
        k(0,0) = 1.E13;
        TPZMaterial * materialMixedLeft = new TPZElasticity3D(-300, young, poisson, force);
        TPZBndCond * mixedLeft = new TPZBndCond(materialMixedLeft,__2DleftMat, mixed, k, f);
        cmesh->InsertMaterialObject(mixedLeft);
        
        k.Zero();
        k(1,1) = 1.E13;
        TPZMaterial * materialMixedOutFracture = new TPZElasticity3D(-301, young, poisson, force);
        TPZBndCond * mixedOutFracture = new TPZBndCond(materialMixedOutFracture,__2DfractureMat_outside, mixed, k, f);
        cmesh->InsertMaterialObject(mixedOutFracture);

        k.Zero();
        k(2,2) = 1.E13;
        TPZMaterial * materialMixedTop = new TPZElasticity3D(-302, young, poisson, force);
        TPZBndCond * mixedTop = new TPZBndCond(materialMixedTop,__2DtopMat, mixed, k, f);
        cmesh->InsertMaterialObject(mixedTop);
        //
        TPZMaterial * materialMixedBottom = new TPZElasticity3D(-303, young, poisson, force);
        TPZBndCond * mixedBottom = new TPZBndCond(materialMixedBottom,__2DbottomMat, mixed, k, f);
        cmesh->InsertMaterialObject(mixedBottom);

        ///////////farField
        k.Zero();
        f(1,0) = pressureY;
        TPZMaterial * materialNewmannFarField = new TPZElasticity3D(-304, young, poisson, force);
        TPZBndCond * newmannFarfield = new TPZBndCond(materialNewmannFarField,__2DfarfieldMat, newmann, k, f);
        cmesh->InsertMaterialObject(newmannFarfield);
        
        ///////////insideFract
        f(1,0) = 0.1;
        TPZMaterial * materialNewmannInsideFract = new TPZElasticity3D(-305, young, poisson, force);
        TPZBndCond * newmannInsideFract = new TPZBndCond(materialNewmannInsideFract,__2DfractureMat_inside, newmann, k, f);
        cmesh->InsertMaterialObject(newmannInsideFract);
    }
    
    cmesh->AutoBuild();
    
    return cmesh;
}
//------------------------------------------------------------------------------------------------------------


void TPZPlaneFracture::GeneratePlaneMesh(std::list<REAL> & espacamento, REAL lengthFactor)
{
    TPZVec< TPZVec<REAL> > NodeCoord(0);
    int nrows, ncols;
    REAL Y = 0.;
    GenerateNodesAtPlaneY(espacamento, lengthFactor, NodeCoord, nrows, ncols, Y);

    int nNodesByLayer = nrows*ncols;
	
	//initializing gmesh->NodeVec()
	fPlaneMesh->NodeVec().Resize(nNodesByLayer);
	TPZVec <TPZGeoNode> Node(nNodesByLayer);
	for(int n = 0; n < nNodesByLayer; n++)
	{
		Node[n].SetNodeId(n);
		Node[n].SetCoord(NodeCoord[n]);
		fPlaneMesh->NodeVec()[n] = Node[n]; 
	}
	
	//inserting quadrilaterals
	int elId = 0;
	TPZVec <int> Topol(4);
	for(int r = 0; r < (nrows-1); r++)
	{
		for(int c = 0; c < (ncols-1); c++)
		{
			Topol[0] = ncols*r+c; Topol[1] = ncols*r+c+1; Topol[2] = ncols*(r+1)+c+1; Topol[3] = ncols*(r+1)+c;
			new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elId,Topol,__2DfractureMat_outside,*fPlaneMesh);
			elId++;
		}
	}
	
	fPlaneMesh->BuildConnectivity();
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFracture::GenerateFullMesh(std::list<REAL> & espacamento, REAL lengthFactor)
{
    TPZVec< TPZVec<REAL> > NodeCoord(0);
    int nrows, ncols;
    
    std::list<REAL>::iterator it = espacamento.end(); it--;
    REAL tickness = 7.;//4.*__maxLength;//tickness is the distance between plane of fracture and plane of farfield
    
    int nDirRef = int(log(fabs(tickness)/__maxLength)/log(2.));
    int nLayers = nDirRef + 2;
    
    REAL Y = 0.;
    GenerateNodesAtPlaneY(espacamento, lengthFactor, NodeCoord, nrows, ncols, Y);
    
    Y = tickness/pow(2.,nDirRef);
    for(int lay = 1; lay < nLayers; lay++)
    {
        GenerateNodesAtPlaneY(espacamento, lengthFactor, NodeCoord, nrows, ncols, Y);
        Y *= 2.;
    }
    int nNodesByLayer = nrows*ncols;
    int Qnodes = nNodesByLayer * nLayers;
	
	//initializing gmesh->NodeVec()
	fFullMesh->NodeVec().Resize(Qnodes);
    
    int pos = 0;
	TPZGeoNode Node;
    for(int l = 0; l < nLayers; l++)
    {
        for(int n = 0; n < nNodesByLayer; n++)
        {
            Node.SetNodeId(pos);
            Node.SetCoord(NodeCoord[pos]);
            fFullMesh->NodeVec()[pos] = Node;
            pos++;
        }
    }
	
	//inserting quadrilaterals
	int elId = 0;
	TPZVec <int> Topol(4);
	for(int r = 0; r < (nrows-1); r++)
	{
		for(int c = 0; c < (ncols-1); c++)
		{
			Topol[0] = ncols*r+c; Topol[1] = ncols*r+c+1; Topol[2] = ncols*(r+1)+c+1; Topol[3] = ncols*(r+1)+c;
			new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elId,Topol,__2DfractureMat_outside,*fFullMesh);
			elId++;
		}
	}
    
    //inserting hexaedrons
    for(int l = 0; l < (nLayers-1); l++)
    {
        for(int r = 0; r < (nrows-1); r++)
        {
            for(int c = 0; c < (ncols-1); c++)
            {
               	Topol.Resize(8);
                Topol[0] = ncols*r+c + l*nNodesByLayer;
                Topol[1] = ncols*r+c+1 + l*nNodesByLayer;
                Topol[2] = ncols*(r+1)+c+1 + l*nNodesByLayer;
                Topol[3] = ncols*(r+1)+c + l*nNodesByLayer;
                //
                Topol[4] = ncols*r+c + (l+1)*nNodesByLayer;
                Topol[5] = ncols*r+c+1 + (l+1)*nNodesByLayer;
                Topol[6] = ncols*(r+1)+c+1 + (l+1)*nNodesByLayer;
                Topol[7] = ncols*(r+1)+c + (l+1)*nNodesByLayer;
                
                new TPZGeoElRefPattern< pzgeom::TPZGeoCube > (elId,Topol,__3DrockMat_linear,*fFullMesh);
                elId++;
                
                Topol.Resize(4);
                if(l == (nLayers-2))//farfield cc
                {
                    Topol[0] = ncols*r+c + (l+1)*nNodesByLayer;
                    Topol[1] = ncols*r+c+1 + (l+1)*nNodesByLayer;
                    Topol[2] = ncols*(r+1)+c+1 + (l+1)*nNodesByLayer;
                    Topol[3] = ncols*(r+1)+c + (l+1)*nNodesByLayer;
                    
                    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elId,Topol,__2DfarfieldMat,*fFullMesh);
                    elId++;
                }
                if(c == 0)//left cc
                {
                    Topol[0] = ncols*r+c + l*nNodesByLayer;
                    Topol[1] = ncols*r+c + (l+1)*nNodesByLayer;
                    Topol[2] = ncols*(r+1)+c + (l+1)*nNodesByLayer;
                    Topol[3] = ncols*(r+1)+c + l*nNodesByLayer;
                    
                    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elId,Topol,__2DleftMat,*fFullMesh);
                    elId++;
                }
                else if(c == (ncols - 2))//right cc
                {
                    Topol[0] = ncols*r+c+1 + l*nNodesByLayer;
                    Topol[1] = ncols*(r+1)+c+1 + l*nNodesByLayer;
                    Topol[2] = ncols*(r+1)+c+1 + (l+1)*nNodesByLayer;
                    Topol[3] = ncols*r+c+1 + (l+1)*nNodesByLayer;
                    
                    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elId,Topol,__2DrightMat,*fFullMesh);
                    elId++;
                }
                if(r == 0)//top cc
                {
                    Topol[0] = ncols*r+c + l*nNodesByLayer;
                    Topol[1] = ncols*r+c+1 + l*nNodesByLayer;
                    Topol[2] = ncols*r+c+1 + (l+1)*nNodesByLayer;
                    Topol[3] = ncols*r+c + (l+1)*nNodesByLayer;
                    
                    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elId,Topol,__2DtopMat,*fFullMesh);
                    elId++;
                }
                else if(r == (nrows - 2))//bottom cc
                {
                    Topol[0] = ncols*(r+1)+c + l*nNodesByLayer;
                    Topol[1] = ncols*(r+1)+c + (l+1)*nNodesByLayer;
                    Topol[2] = ncols*(r+1)+c+1 + (l+1)*nNodesByLayer;
                    Topol[3] = ncols*(r+1)+c+1 + l*nNodesByLayer;
                    
                    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elId,Topol,__2DbottomMat,*fFullMesh);
                    elId++;
                }
            }
        }
    }
    
    // ////////////////para validar primeira simulacao completa
    Topol.Resize(1);
    
    Topol[0] = 0;
    new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (elId,Topol,__aux0DEl_Mat1,*fFullMesh);
    elId++;
    
    Topol[0] = ncols*(nrows-1);
    new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (elId,Topol,__aux0DEl_Mat2,*fFullMesh);
    elId++;
    
    Topol[0] = ncols-1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (elId,Topol,__aux0DEl_Mat3,*fFullMesh);
    elId++;
    
    Topol[0] = ncols*nrows-1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (elId,Topol,__aux0DEl_Mat4,*fFullMesh);
    elId++;
    // ////////////////////////////////////////////////////////MUSTDELETE
    
	fFullMesh->BuildConnectivity();
    fFullMesh->SetMaxElementId(fFullMesh->NElements()-1);
    fFullMesh->SetMaxNodeId(fFullMesh->NNodes()-1);
    
//    std::ofstream out("FullMesh.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(fFullMesh, out, true);
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFracture::GenerateNodesAtPlaneY(std::list<REAL> & espacamento, REAL lengthFactor, 
                                             TPZVec< TPZVec<REAL> > & NodeCoord, int & nrows, int & ncols,
                                             REAL Y)
{
    nrows = espacamento.size();
    
    std::list<REAL>::iterator it = espacamento.end(); it--;
    //REAL lastPos = fabs(*it);
    
    REAL Lx = 5.;//lengthFactor * lastPos;

    int nStretches = 1;
    REAL deltaX = Lx/nStretches;
    while(deltaX > __maxLength)
    {
        nStretches++;
        deltaX = Lx/nStretches;
    }
    ncols = nStretches+1;
    
	const int nNodesByLayer = nrows*ncols;
    int oldNNodes = NodeCoord.NElements();
    int newNNodes = oldNNodes + nNodesByLayer;
    
	NodeCoord.Resize(newNNodes);
	
	int nodeId = oldNNodes;
	
	//setting nodes coords
	for(it = espacamento.begin(); it != espacamento.end(); it++)
	{
		for(int c = 0; c < ncols; c++)
		{
            NodeCoord[nodeId].Resize(3,0.);
			NodeCoord[nodeId][0] = c*deltaX;
            NodeCoord[nodeId][1] = Y;
			NodeCoord[nodeId][2] = *it;
			nodeId++;
		}
	}
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFracture::DetectEdgesCrossed(const TPZVec<REAL> &poligonalChain, TPZGeoMesh * planeMesh,
                                          std::map< int, std::set<REAL> > &elId_TrimCoords, 
                                          std::list< std::pair<int,REAL> > &elIdSequence)
{
	int npoints = (poligonalChain.NElements())/2;
	int nelem = planeMesh->NElements();
	
    TPZVec<REAL> coord(3,0.);
    int firstPt = 0;
    coord[0] = poligonalChain[2*firstPt];
    coord[2] = poligonalChain[2*firstPt+1];
    
    int initialElId = 0;

    int axe0 = 0;//axe X
    int axe1 = 2;//axe Z
    int axeNormal = 1;//axe Y
    TPZVec<REAL> qsi(2,0.);
    int elFoundId = PointElementOnPlaneMesh(fPlaneMesh, initialElId, coord, qsi, axe0, axe1, axeNormal, true);
    
    TPZGeoEl * firstGel = planeMesh->ElementVec()[elFoundId];

	TPZGeoEl * gel = firstGel;
	TPZGeoEl * nextGel = NULL;
	
    #ifdef DEBUG
	if(!gel)
	{
		std::cout << "first point of crack tip boundary does NOT belong to any 2D element" << std::endl;
		std::cout << "(or was given a mesh with zero quantity of 2D elements)!" << std::endl;
		std::cout << "See " << __PRETTY_FUNCTION__ << std::endl;
		DebugStop();
	}
    #endif
	
	REAL alphaMin;
	bool reachNextPoint;
	int nElsCrossed, thispoint, nextpoint;
	std::map< int, std::set<REAL> >::iterator it;
	std::set<REAL> trim;
	TPZVec< TPZVec<REAL> > intersectionPoint;
	TPZManVector<REAL,3> x(3,0.), dx(3,0.);
    TPZManVector<REAL,3> xNext(3,0.), qsi2D(2,0.);
	TPZVec <int> Topol(2), edgeVec;
	
	int p;
	for(p = 0; p < (npoints-1); p++)
	{
		nElsCrossed = 0;
		alphaMin = 0.;
		thispoint = 2*p;
		nextpoint = 2*(p+1);
        
        //-----------------------coordX
        REAL xcoord = poligonalChain[thispoint+0];
        x[0] = xcoord;
        
        REAL dxcoord = poligonalChain[nextpoint+0] - poligonalChain[thispoint+0];
        dx[0] = dxcoord;
        
        REAL xnextcoord = poligonalChain[nextpoint+0];
        xNext[0] = xnextcoord;
        //-----------------------coordY
        x[1] = 0.;
        dx[1] = 0.;
        xNext[1] = 0.;
        //-----------------------coordZ
        xcoord = poligonalChain[thispoint+1];
        x[2] = xcoord;
        
        dxcoord = poligonalChain[nextpoint+1] - poligonalChain[thispoint+1];
        dx[2] = dxcoord;
        
        xnextcoord = poligonalChain[nextpoint+1];
        xNext[2] = xnextcoord;
        //-----------------------
    
		REAL norm = 0.;
		for(int c = 0; c < 3; c++)
		{
			norm += dx[c]*dx[c];
		}
		norm = sqrt(norm);		
		for(int c = 0; c < 3; c++)
		{
			dx[c] = dx[c]/norm;
		}
		
		reachNextPoint = gel->ComputeXInverse(xNext, qsi2D);
        int axe0 = 0;//axe X
        int axe1 = 2;//axe Z
        int axeNormal = 1;//axe Y
		while(reachNextPoint == false && nElsCrossed < nelem)
		{
			nextGel = CrossToNextNeighbour(gel, x, dx, alphaMin, elId_TrimCoords, elIdSequence, true, axe0, axe1, axeNormal, false);
			
			alphaMin = __smallNum;//qdo vai de vizinho a vizinho ateh chegar no proximo ponto, nao deve-se incluir os (alphaX=0)
			if(nextGel->ComputeXInverse(xNext, qsi2D))
			{
				reachNextPoint = true;
			}
			gel = nextGel;
			nElsCrossed++;
		}
		if(nElsCrossed == nelem)
		{
			//Deve ter alternado entre vizinhos!!!
			DebugStop();
		}
	}
    
    TPZGeoEl * lastGel = gel;
	
	dx[0] = -1.;//direcao oposta ao eixo x do sistema de coordenadas
    dx[1] =  0.;//direcao oposta ao eixo x do sistema de coordenadas
    dx[2] =  0.;//direcao oposta ao eixo x do sistema de coordenadas
	
	p = 0;//Fechando o inicio da fratura
	thispoint = 2*p;
	nextpoint = 2*(p+1);

    x[0] = poligonalChain[thispoint+0];
    x[1] = 0.;
    x[2] = poligonalChain[thispoint+1];

	nextGel = firstGel;
	gel = NULL;
	alphaMin = 0.;
	while(gel != nextGel)
	{
		gel = nextGel;
		nextGel = CrossToNextNeighbour(gel, x, dx, alphaMin, elId_TrimCoords, elIdSequence, false, axe0, axe1, axeNormal, true);
		alphaMin = __smallNum;
	}
	
	p = npoints;//Fechando o final da fratura
	thispoint = 2*(p-1);

    x[0] = poligonalChain[thispoint+0];
    x[1] = 0.;
    x[2] = poligonalChain[thispoint+1];
    
	nextGel = lastGel;
	gel = NULL;
	alphaMin = 0.;
	while(gel != nextGel)
	{
		gel = nextGel;
		nextGel = CrossToNextNeighbour(gel, x, dx, alphaMin, elId_TrimCoords, elIdSequence, true, axe0, axe1, axeNormal, true);
		alphaMin = __smallNum;
	}
}
//------------------------------------------------------------------------------------------------------------

TPZGeoEl * TPZPlaneFracture::CrossToNextNeighbour(TPZGeoEl * gel, TPZVec<REAL> &x, TPZVec<REAL> dx, REAL alphaMin, int planeAxe0, int planeAxe1, int planeNormal)
{
	std::map< int, std::set<REAL> >::iterator it;
	std::set<REAL> trim;
	TPZVec< TPZVec<REAL> > ExactIntersectionPoint;
	TPZVec<REAL> qsi1Dvec(1), xCrackBoundary(3);
	TPZVec<int> Topol(2), edgeVec;
	
	bool haveIntersection = EdgeIntersection(gel, x, dx, edgeVec, ExactIntersectionPoint, alphaMin, planeAxe0, planeAxe1, planeNormal);
	
	if(haveIntersection == false)
	{
		haveIntersection = EdgeIntersection(gel, x, dx, edgeVec, ExactIntersectionPoint, 0., planeAxe0, planeAxe1, planeNormal);
        
#ifdef DEBUG
		if(haveIntersection == false)
		{
            TPZManVector<REAL,3> qsi2D(2,0.);
            if(gel->ComputeXInverse(x, qsi2D))
            {
                std::cout << "The point inside element does NOT intersect its edges! EdgeIntersection method face an exeption!" << std::endl;
                std::cout << "See " << __PRETTY_FUNCTION__ << std::endl;
                
                std::cout << "Element " << gel->Id() << std::endl;
                for(int n = 0; n < gel->NNodes(); n++)
                {
                    std::cout << "Node " << n << std::endl;
                    gel->NodePtr(n)->Print(std::cout);
                    std::cout << std::endl;
                }
                std::cout << "x: " << x[0] << " , " << x[1] << " , " << x[2] << std::endl;
                std::cout << "dx: " << dx[0] << " , " << dx[1] << " , " << dx[2] << std::endl;
                std::cout << "alphaMin = " << alphaMin << std::endl << std::endl;
                
                DebugStop();
            }
            else
            {
                std::cout << "Trying to find edge intersection from an point that doesnt belong to given element!" << std::endl;
                std::cout << "See " << __PRETTY_FUNCTION__ << std::endl;
                
                DebugStop();
            }
		}
#endif
	}
    int edge = edgeVec[edgeVec.NElements()-1];
	x = ExactIntersectionPoint[ExactIntersectionPoint.NElements() - 1];
	TPZGeoElSide gelEdge(gel, edge);
	TPZGeoElSide neighEdge = gelEdge.Neighbour();
    bool neighFound = false;
	while(neighEdge != gelEdge)
	{
		if(neighEdge.Element()->Dimension() == 2 && neighEdge.Element()->Father() == NULL)
		{
            neighFound = true;
			break;
		}
		neighEdge = neighEdge.Neighbour();
	}
    if(neighFound)
    {
        gel = neighEdge.Element();
    }
    else
    {
        DebugStop();
    }
	
	return gel;
}
//------------------------------------------------------------------------------------------------------------

TPZGeoEl * TPZPlaneFracture::CrossToNextNeighbour(TPZGeoEl * gel, TPZVec<REAL> &x, TPZVec<REAL> dx, REAL alphaMin,
												  std::map< int, std::set<REAL> > &elId_TrimCoords,
                                                  std::list< std::pair<int,REAL> > &elIdSequence, 
                                                  bool pushback, int planeAxe0, int planeAxe1, int planeNormal,
                                                  bool closingFracture)
{
	bool thereIsAn1DElemAlready;
	int edge;
	std::map< int, std::set<REAL> >::iterator it;
	std::set<REAL> trim;
	TPZVec< TPZVec<REAL> > ExactIntersectionPoint, ModulatedIntersectionPoint;
	TPZVec<REAL> qsi1Dvec(1), xCrackBoundary(3);
	TPZVec<int> Topol(2), edgeVec;
	
	TPZGeoMesh * planeMesh = gel->Mesh();
	
	bool haveIntersection = EdgeIntersection(gel, x, dx, edgeVec, ExactIntersectionPoint, ModulatedIntersectionPoint, alphaMin, planeAxe0, planeAxe1, planeNormal);
	
	if(haveIntersection == false)
	{
		haveIntersection = EdgeIntersection(gel, x, dx, edgeVec, ExactIntersectionPoint, ModulatedIntersectionPoint, 0., planeAxe0, planeAxe1, planeNormal);
        
        #ifdef DEBUG
		if(haveIntersection == false)
		{
            TPZManVector<REAL,3> qsi2D(2,0.);
            if(gel->ComputeXInverse(x, qsi2D))
            {
                std::cout << "The point inside element does NOT intersect its edges! EdgeIntersection method face an exeption!" << std::endl;
                std::cout << "See " << __PRETTY_FUNCTION__ << std::endl;
                
                std::cout << "Element " << gel->Id() << std::endl;
                for(int n = 0; n < gel->NNodes(); n++)
                {
                    std::cout << "Node " << n << std::endl;
                    gel->NodePtr(n)->Print(std::cout);
                    std::cout << std::endl;
                }
                std::cout << "x: " << x[0] << " , " << x[1] << " , " << x[2] << std::endl;
                std::cout << "dx: " << dx[0] << " , " << dx[1] << " , " << dx[2] << std::endl;
                std::cout << "alphaMin = " << alphaMin << std::endl << std::endl;
                
                DebugStop();
            }
            else
            {
                std::cout << "Trying to find edge intersection from an point that doesnt belong to given element!" << std::endl;
                std::cout << "See " << __PRETTY_FUNCTION__ << std::endl;

                DebugStop();
            }
		}
        #endif
	}
	TPZVec<REAL> xLin(3);
	for(int edg = 0; edg < edgeVec.NElements(); edg++)
	{
		edge = edgeVec[edg];
		TPZVec<REAL> n0(3), n1(3);
		planeMesh->NodeVec()[gel->SideNodeIndex(edge, 0)].GetCoordinates(n0);
		planeMesh->NodeVec()[gel->SideNodeIndex(edge, 1)].GetCoordinates(n1);
		for(int c = 0; c < 3; c++)
		{
            REAL coordM = ModulatedIntersectionPoint[edg][c];
			xLin[c] = coordM;
		}
		REAL qsi1D = LinearComputeXInverse(xLin, n0, n1);
		
		TPZGeoElSide gelEdge(gel, edge);
		TPZGeoElSide neighEdge = gelEdge.Neighbour();
		thereIsAn1DElemAlready = false;
		while(neighEdge != gelEdge)
		{
			if(neighEdge.Element()->Dimension() == 1)//jah existe um elemento 1D inserido nesta aresta!
			{
				thereIsAn1DElemAlready = true;
				int neighEdgeId = neighEdge.Element()->Id();
				it = elId_TrimCoords.find(neighEdgeId);
				TPZTransform transBetweenNeigh = neighEdge.NeighbourSideTransform(gelEdge);
				qsi1D *= transBetweenNeigh.Mult()(0,0);
				it->second.insert(qsi1D);
				if(pushback)
				{
					elIdSequence.push_back(std::make_pair(neighEdgeId, qsi1D));
				}
				else // push_FRONT
				{
					elIdSequence.push_front(std::make_pair(neighEdgeId, qsi1D));
				}

				break;
			}
			neighEdge = neighEdge.Neighbour();
		}
		if(thereIsAn1DElemAlready == false)//nao existe um elemento 1D nesta aresta!
		{
			trim.clear();
			trim.insert(qsi1D);
			Topol[0] = gel->SideNodeIndex(edge, 0);
			Topol[1] = gel->SideNodeIndex(edge, 1);
			TPZGeoElRefPattern< pzgeom::TPZGeoLinear > * linGeo = 
                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear >(Topol, __aux1DEl_Mat, *planeMesh);
            
			planeMesh->BuildConnectivity();
			
			int linGeoId = linGeo->Id();
			elId_TrimCoords[linGeoId] = trim;
			
			if(pushback) // push_BACK
			{
				elIdSequence.push_back(std::make_pair(linGeoId, qsi1D));
			}
			else // push_FRONT
			{
				elIdSequence.push_front(std::make_pair(linGeoId, qsi1D));
			}
		}
	}
	x = ExactIntersectionPoint[ExactIntersectionPoint.NElements() - 1];
	TPZGeoElSide gelEdge(gel, edge);
	TPZGeoElSide neighEdge = gelEdge.Neighbour();
	bool neighFound = false;
	while(neighEdge != gelEdge)
	{
		if(neighEdge.Element()->Dimension() == 2 && neighEdge.Element()->Father() == NULL)
		{
            neighFound = true;
			break;
		}
		neighEdge = neighEdge.Neighbour();
	}
    if(neighFound)
    {
        gel = neighEdge.Element();
    }
    else if(closingFracture == false)
    {
        DebugStop();//fratura saiu para fora do dominio
    }
	
	return gel;
}
//------------------------------------------------------------------------------------------------------------

bool TPZPlaneFracture::EdgeIntersection(TPZGeoEl * gel, TPZVec<REAL> &x, TPZVec<REAL> &dx, TPZVec<int> &edge, 
                       TPZVec< TPZVec<REAL> > &ExactIntersect, REAL alphaMin,
                       int planeAxe0, int planeAxe1, int planeNormal)
{
    int nearNode;
	bool IsNearNode = TPZChangeEl::NearestNode(gel, x, nearNode, __smallNum);
	
	edge.Resize(0);
	ExactIntersect.Resize(0);
	
	int ncnodes = gel->NCornerNodes();
	
	TPZVec< TPZVec<REAL> > node(ncnodes), dnode(ncnodes);
	TPZVec<REAL> nodeCoord(3);
	for(int n = 0; n < ncnodes; n++)
	{
		node[n].Resize(3,1);
		dnode[n].Resize(3,1);
		gel->NodePtr(n)->GetCoordinates(nodeCoord);
		for(int c = 0; c < 3; c++)
		{
			(node[n])[c] = nodeCoord[c];
		}
	}
	std::list<REAL> edgeNorm;
	/** <REAL: alpha , <REAL: alphaNodemod, REAL: alphaNodesmooth, REAL: norm, int: edge that intersect> > */
    std::map<REAL, TPZVec<REAL> > alpha;
    
	REAL alphaX, alphaNodemod, alphaNodesmooth, norm;
	for(int n = 0; n < ncnodes; n++)
	{	
		norm = 0.;
		for(int c = 0; c < 3; c++)//computing the direction from node N to node N+1
		{
			if(n != (ncnodes-1) )
			{
				(dnode[n])[c] = node[n+1][c] - node[n][c];
			}
			else
			{
				dnode[n][c] = node[0][c] - node[n][c];
			}
			norm += dnode[n][c]*dnode[n][c];
		}
		norm = sqrt(norm);
		for(int c = 0; c < 3; c++)//normalizing computed direction
		{
			dnode[n][c] /= norm;
		}
		
		alphaX = ComputeAlphaX(x, dx, node[n], dnode[n], planeAxe0, planeAxe1, planeNormal);
		alphaNodemod = ComputeAlphaNode(x, dx, node[n], dnode[n], norm, true, false, planeAxe0, planeAxe1, planeNormal);
		alphaNodesmooth = ComputeAlphaNode(x, dx, node[n], dnode[n], norm, IsNearNode, true, planeAxe0, planeAxe1, planeNormal);
		if(alphaX >= alphaMin && alphaNodesmooth >= 0. && alphaNodesmooth <= norm)
		{
			int thisEdge = n+ncnodes;
			TPZVec<REAL> someData(4);
			someData[0] = alphaNodemod; someData[1] = alphaNodesmooth; someData[2] = norm; someData[3] = thisEdge;
			alpha[alphaX] = someData;
		}
	}
	
	// a aresta que serah interseccionada serah a que tiver menor alphaX não negativo, i.e.: o primeiro par do mapa alpha!
	std::map<REAL, TPZVec<REAL> >::iterator it = alpha.begin();
	if(it != alpha.end())
	{		
		alphaX = it->first;
		alphaNodemod = (it->second)[0];
		alphaNodesmooth = (it->second)[1];
		norm = (it->second)[2];
		edge.Resize(1); edge[0] = int((it->second)[3]);	
		
		ExactIntersect.Resize(1); ExactIntersect[0].Resize(3,1);
		for(int c = 0; c < 3; c++)
		{
			ExactIntersect[0][c] = node[edge[0]-ncnodes][c] + alphaNodesmooth*dnode[edge[0]-ncnodes][c];
		}
		
		if(alphaX <= alphaMin && alpha.size() > 1)
		{
			it++;
			alphaX = it->first;
			alphaNodemod = (it->second)[0];
			alphaNodesmooth = (it->second)[1];
			norm = (it->second)[2];
			edge.Resize(2); edge[1] = int((it->second)[3]);	
			
			ExactIntersect.Resize(2); ExactIntersect[1].Resize(3,1);
			for(int c = 0; c < 3; c++)
			{
				ExactIntersect[1][c] = node[edge[1]-ncnodes][c] + alphaNodesmooth*dnode[edge[1]-ncnodes][c];
			}
		}
		return true;
	}
	
	//se o ponto p estah sobre um noh e a direcao dp aponta para fora do elemento...
	if(IsNearNode)
	{
		std::cout << "Estah no noh " << nearNode << " e nao foi encontrada interseccao!" << std::endl;
		std::cout << "Tratar este caso!" << std::endl << std::endl;
	}
	
	return false;
}
//------------------------------------------------------------------------------------------------------------

bool TPZPlaneFracture::EdgeIntersection(TPZGeoEl * gel, TPZVec<REAL> &x, TPZVec<REAL> &dx, TPZVec<int> &edge,
										TPZVec< TPZVec<REAL> > &ExactIntersect,
                                        TPZVec< TPZVec<REAL> > &ModulatedIntersect, REAL alphaMin,
                                        int planeAxe0, int planeAxe1, int planeNormal)
{
	int nearNode;
	bool IsNearNode = TPZChangeEl::NearestNode(gel, x, nearNode, __smallNum);
	
	edge.Resize(0);
	ExactIntersect.Resize(0);
	ModulatedIntersect.Resize(0);
	
	int ncnodes = gel->NCornerNodes();
	
	TPZVec< TPZVec<REAL> > node(ncnodes), dnode(ncnodes);
	TPZVec<REAL> nodeCoord(3);
	for(int n = 0; n < ncnodes; n++)
	{
		node[n].Resize(3,1);
		dnode[n].Resize(3,1);
		gel->NodePtr(n)->GetCoordinates(nodeCoord);
		for(int c = 0; c < 3; c++)
		{
			(node[n])[c] = nodeCoord[c];
		}
	}
	std::list<REAL> edgeNorm;
	/** <REAL: alpha , <REAL: alphaNodemod, REAL: alphaNodesmooth, REAL: norm, int: edge that intersect> > */
    std::map<REAL, TPZVec<REAL> > alpha;
    
	REAL alphaX, alphaNodemod, alphaNodesmooth, norm;
	for(int n = 0; n < ncnodes; n++)
	{	
		norm = 0.;
		for(int c = 0; c < 3; c++)//computing the direction from node N to node N+1
		{
			if(n != (ncnodes-1) )
			{
				(dnode[n])[c] = node[n+1][c] - node[n][c];
			}
			else
			{
				dnode[n][c] = node[0][c] - node[n][c];
			}
			norm += dnode[n][c]*dnode[n][c];
		}
		norm = sqrt(norm);
		for(int c = 0; c < 3; c++)//normalizing computed direction
		{
			dnode[n][c] /= norm;
		}
		
		alphaX = ComputeAlphaX(x, dx, node[n], dnode[n], planeAxe0, planeAxe1, planeNormal);
		alphaNodemod = ComputeAlphaNode(x, dx, node[n], dnode[n], norm, true, false, planeAxe0, planeAxe1, planeNormal);
		alphaNodesmooth = ComputeAlphaNode(x, dx, node[n], dnode[n], norm, IsNearNode, true, planeAxe0, planeAxe1, planeNormal);
		if(alphaX >= alphaMin && alphaNodesmooth >= 0. && alphaNodesmooth <= norm)
		{
			int thisEdge = n+ncnodes;
			TPZVec<REAL> someData(4);
			someData[0] = alphaNodemod; someData[1] = alphaNodesmooth; someData[2] = norm; someData[3] = thisEdge;
			alpha[alphaX] = someData;
		}
	}
	
	// a aresta que serah interseccionada serah a que tiver menor alphaX não negativo, i.e.: o primeiro par do mapa alpha!
	std::map<REAL, TPZVec<REAL> >::iterator it = alpha.begin();
	if(it != alpha.end())
	{		
		alphaX = it->first;
		alphaNodemod = (it->second)[0];
		alphaNodesmooth = (it->second)[1];
		norm = (it->second)[2];
		edge.Resize(1); edge[0] = int((it->second)[3]);	
		
		ModulatedIntersect.Resize(1); ModulatedIntersect[0].Resize(3,1);
		for(int c = 0; c < 3; c++)
		{
			ModulatedIntersect[0][c] = node[edge[0]-ncnodes][c] + alphaNodemod*dnode[edge[0]-ncnodes][c];
		}
		ExactIntersect.Resize(1); ExactIntersect[0].Resize(3,1);
		for(int c = 0; c < 3; c++)
		{
			ExactIntersect[0][c] = node[edge[0]-ncnodes][c] + alphaNodesmooth*dnode[edge[0]-ncnodes][c];
		}
		
		if(alphaX <= alphaMin && alpha.size() > 1)
		{
			it++;
			alphaX = it->first;
			alphaNodemod = (it->second)[0];
			alphaNodesmooth = (it->second)[1];
			norm = (it->second)[2];
			edge.Resize(2); edge[1] = int((it->second)[3]);	
			
			ModulatedIntersect.Resize(2); ModulatedIntersect[1].Resize(3,1);
			for(int c = 0; c < 3; c++)
			{
				ModulatedIntersect[1][c] = node[edge[1]-ncnodes][c] + alphaNodemod*dnode[edge[1]-ncnodes][c];
			}
			ExactIntersect.Resize(2); ExactIntersect[1].Resize(3,1);
			for(int c = 0; c < 3; c++)
			{
				ExactIntersect[1][c] = node[edge[1]-ncnodes][c] + alphaNodesmooth*dnode[edge[1]-ncnodes][c];
			}
		}
		return true;
	}
	
	//se o ponto p estah sobre um noh e a direcao dp aponta para fora do elemento...
	if(IsNearNode)
	{
		std::cout << "Estah no noh " << nearNode << " e nao foi encontrada interseccao!" << std::endl;
		std::cout << "Tratar este caso!" << std::endl << std::endl;
	}
	
	return false;
}
//------------------------------------------------------------------------------------------------------------

REAL TPZPlaneFracture::ComputeAlphaNode(TPZVec<REAL> &x, TPZVec<REAL> &dx,
                                          TPZVec<REAL> &node, TPZVec<REAL> &dnode,
                                          REAL norm, bool modulate, bool smooth,
                                          int planeAxe0, int planeAxe1, int planeNormal)
{
	REAL fractionNumQ =	dx[planeAxe1]*node[planeAxe0] - dx[planeAxe0]*node[planeAxe1] - dx[planeAxe1]*x[planeAxe0] + dx[planeAxe0]*x[planeAxe1];
	
	REAL fractionDenomQ =	dnode[planeAxe1]*dx[planeAxe0] - dnode[planeAxe0]*dx[planeAxe1];
	
	REAL alphaNode = -1.;
	if(fabs(fractionDenomQ) > __smallNum)
	{
		alphaNode = fractionNumQ/fractionDenomQ;
		if(fabs(alphaNode) < __smallNum) alphaNode = 0.;
		if(fabs(alphaNode - norm) < __smallNum) alphaNode = norm;
	}
	
	if(modulate)
	{
		int stretchesQTD = __EdgeStretchesQTD;
		if(smooth)
		{
			stretchesQTD *= __TrimQTDmultiplier;
		}
		int nsegm = int(alphaNode/(norm/stretchesQTD) + 0.5);
		if(nsegm == 0)//Tirando a interseccao do noh inicial da aresta
		{
			nsegm = 1;
		}
		else if(nsegm == stretchesQTD)//Tirando a interseccao do noh final da aresta
		{
			nsegm = stretchesQTD - 1;
		}
		alphaNode = nsegm*(norm/stretchesQTD);//modulando o ponto para multiplo de (norm/__EdgeStretchesQTD)
	}
	
	return alphaNode;
}
//------------------------------------------------------------------------------------------------------------

REAL TPZPlaneFracture::ComputeAlphaX(TPZVec<REAL> &x, TPZVec<REAL> &dx,
                                       TPZVec<REAL> &node, TPZVec<REAL> &dnode,
                                       int planeAxe0, int planeAxe1, int planeNormal)
{	
	//computing alpha (dx vector multiplier to intersect element edges)
	REAL fractionNumP   =	dnode[planeAxe1]*node[planeAxe0] - dnode[planeAxe0]*node[planeAxe1] -
	dnode[planeAxe1]*x[planeAxe0] + dnode[planeAxe0]*x[planeAxe1];
	
	REAL fractionDenomP =	dnode[planeAxe1]*dx[planeAxe0] - dnode[planeAxe0]*dx[planeAxe1];
	
	REAL alphaX = -1.;
	if(fabs(fractionDenomP) > __smallNum)
	{
		alphaX = fractionNumP/fractionDenomP;
		if(fabs(alphaX) < __smallNum)
		{
			alphaX = 0.;
		}
	}
	
	return alphaX;
}
//------------------------------------------------------------------------------------------------------------

REAL TPZPlaneFracture::LinearComputeXInverse(TPZVec<REAL> x, TPZVec<REAL> n0, TPZVec<REAL> n1)
{
	REAL dL = 0., L = 0.;
	for(int c = 0; c < 3; c++)
	{
		L += (n1[c]-n0[c])*(n1[c]-n0[c]);
		dL += (x[c]-n0[c])*(x[c]-n0[c]);
	}
	L = sqrt(L);
	dL = sqrt(dL);
	
	#ifdef DEBUG
	if(fabs(L) < __smallNum || fabs(dL) < __smallNum || fabs(L - dL) < __smallNum)
	{
		std::cout << "n0 and n1 are coincident nodes!" << std::endl;
		std::cout << "See " << __PRETTY_FUNCTION__ << std::endl;
		DebugStop();
	}
	#endif
	
	REAL qsi = -1. + dL/L*2.;
	
	return qsi;
}
//------------------------------------------------------------------------------------------------------------

TPZAutoPointer<TPZRefPattern> TPZPlaneFracture::Generate1DRefPatt(std::set<REAL> &TrimCoord)
{
	if(TrimCoord.size() == 0)
	{
		return NULL;
	}
	
	int Qnodes = TrimCoord.size() + 2;
	int Qelements = TrimCoord.size() + 2;
	int QsubElements = Qelements - 1;
	TPZVec < TPZVec <REAL> > NodeCoord(Qnodes);
	for(int i = 0; i < Qnodes; i++) NodeCoord[i].Resize(3,0.);
	
	//1. setting initial and final nodes coordinates of 1D element mesh
	NodeCoord[0][0] = -1.;
	NodeCoord[1][0] =  1.;
	//.
	
	//2. setting intermediate nodes coordinates of 1D element mesh
	int c = 2;
	std::set<REAL>::iterator it;
	for(it = TrimCoord.begin(); it != TrimCoord.end(); it++)
	{
		REAL coord = *it;
		
		(NodeCoord[c])[0] = coord;
		c++;
	}
	//.
	
	//3. initializing internal mesh of refPattern
	TPZGeoMesh internalMesh;
	internalMesh.SetMaxNodeId(Qnodes-1);
	internalMesh.SetMaxElementId(Qelements-1);
	internalMesh.NodeVec().Resize(Qnodes);
	TPZVec <TPZGeoNode> Node(Qnodes);
	for(int n = 0; n < Qnodes; n++)
	{
		Node[n].SetNodeId(n);
		Node[n].SetCoord(NodeCoord[n]);
		internalMesh.NodeVec()[n] = Node[n]; 
	}
	//.
	
	//4. inserting 1D elements on internal mesh of refPattern
	int elId = 0;
	TPZVec <int> Topol(2);
	
	//4.1 inserting father element
	Topol[0] = 0; Topol[1] = 1;
	TPZGeoElRefPattern< pzgeom::TPZGeoLinear > * father = new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elId,Topol,__aux1DEl_Mat,internalMesh);
	elId++;
	//.
	
	//4.2 inserting subelements
	//first subelement
	Topol[0] = 0; Topol[1] = 2;
	TPZGeoElRefPattern< pzgeom::TPZGeoLinear > * son1 = new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elId,Topol,__aux1DEl_Mat,internalMesh);
	son1->SetFather(father);
	son1->SetFather(father->Index());
	elId++;
	//
	
	//last subelement
	Topol[0] = TrimCoord.size() + 1; Topol[1] = 1;
	TPZGeoElRefPattern< pzgeom::TPZGeoLinear > * son2 = new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elId,Topol,__aux1DEl_Mat,internalMesh);
	son2->SetFather(father);
	son2->SetFather(father->Index());
	elId++;
	//
	
	for(int el = 2; el < QsubElements; el++)
	{
		Topol[0] = el; Topol[1] = el+1;
		TPZGeoElRefPattern< pzgeom::TPZGeoLinear > * son = new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elId,Topol,__aux1DEl_Mat,internalMesh);
		son->SetFather(father);
		son->SetFather(father->Index());
		elId++;
	}
	//.
	//.
	
	internalMesh.BuildConnectivity();
	
	TPZAutoPointer<TPZRefPattern> refPattern = new TPZRefPattern(internalMesh);
	TPZAutoPointer<TPZRefPattern> Found = gRefDBase.FindRefPattern(refPattern);
	if(!Found)
	{
		gRefDBase.InsertRefPattern(refPattern);
		refPattern->InsertPermuted();
		
		return refPattern;
	}
	else
	{
		return Found;
	}
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFracture::UpdatePoligonalChain(TPZGeoMesh * gmesh,
                                            std::list< std::pair<int,REAL> > &elIdSequence,
                                            TPZVec<REAL> &poligonalChainUpdated)
{
	int nptos = elIdSequence.size();
	poligonalChainUpdated.Resize(2*nptos);

	TPZVec<REAL> qsi1Dvec(1), ptoCoord(3);
	int el1Did, posX, posZ, p = 0;
	REAL qsi1D;
	
	std::list< std::pair<int,REAL> >::iterator it;
	for(it = elIdSequence.begin(); it != elIdSequence.end(); it++)
	{
		el1Did = it->first;
		TPZGeoEl * el1D = gmesh->ElementVec()[el1Did];
		
		#ifdef DEBUG
		int elDim = el1D->Dimension();
		if(elDim != 1)
		{
			std::cout << "The elIdSequence supposedly would contains ids of elements exclusively 1D!" << std::endl;
			std::cout << "See " << __PRETTY_FUNCTION__ << std::endl;
			DebugStop();
		}
		#endif
		
		qsi1D = it->second;
		qsi1Dvec[0] = qsi1D;
		
		el1D->X(qsi1Dvec, ptoCoord);
		
		posX = 2*p;
		posZ = 2*p+1;
		
		poligonalChainUpdated[posX] = ptoCoord[0];
		poligonalChainUpdated[posZ] = ptoCoord[2];
		
		p++;
	}
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFracture::GenerateCrackBoundary(TPZGeoMesh * gmesh2D,
                                             TPZGeoMesh * gmesh3D,
                                             std::list< std::pair<int,REAL> > &elIdSequence)
{
    fcrackBoundaryElementsIds.Resize(0);
	TPZVec<REAL> qsi0vec(1), qsi1vec(1), node0coord(3), node1coord(3);
	TPZVec<int> Topol(2);
	int el0id, el1id, n0, n1;
	REAL qsi0, qsi1;
	std::list< std::pair<int,REAL> >::iterator crackit0, crackit1, crackitEnd;
	crackitEnd = elIdSequence.end(); crackitEnd--;
	for(crackit0 = elIdSequence.begin(); crackit0 != crackitEnd; crackit0++)
	{
		crackit1 = crackit0; crackit1++;
		
		el0id = crackit0->first;
		TPZGeoEl * el0 = gmesh2D->ElementVec()[el0id];
		qsi0 = crackit0->second;
		qsi0vec[0] = qsi0;
		el0->X(qsi0vec, node0coord);
		n0 = TPZChangeEl::NearestNode(gmesh3D, node0coord, __smallNum);
		Topol[0] = n0;
		
		el1id = crackit1->first;
		TPZGeoEl * el1 = gmesh2D->ElementVec()[el1id];
		qsi1 = crackit1->second;
		qsi1vec[0] = qsi1;
		el1->X(qsi1vec, node1coord);
		n1 = TPZChangeEl::NearestNode(gmesh3D, node1coord, __smallNum);
		Topol[1] = n1;
		
		TPZGeoEl * crack1D = new TPZGeoElRefPattern< pzgeom::TPZGeoLinear >(Topol, __1DcrackTipMat, *gmesh3D);
        int oldSize = fcrackBoundaryElementsIds.NElements();
        fcrackBoundaryElementsIds.Resize(oldSize+1);
        fcrackBoundaryElementsIds[oldSize] = crack1D->Id();
	}
    
    gmesh3D->BuildConnectivity();
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFracture::SeparateElementsInMaterialSets(TPZGeoMesh * fullMesh)
{
    int n1Dels = fcrackBoundaryElementsIds.NElements();
    std::map<int,TPZFracture2DEl> fracturedElems;
    
    fcrackQpointsElementsIds.clear();

    //Capturando subelementos que encostam no contorno da fratura
    for(int el = 0; el < n1Dels; el++)
    {
        int id = fcrackBoundaryElementsIds[el];
        TPZGeoEl * gel = fullMesh->ElementVec()[id];//1D element of crach boundary
        
        #ifdef DEBUG
        if(gel->Dimension() != 1)
        {
            DebugStop();
        }
        #endif
        
        for(int sd = 0; sd < gel->NSides(); sd++)
        {            
            TPZGeoElSide side1D(gel,sd);
            TPZGeoElSide sideNeigh = side1D.Neighbour();
            while(sideNeigh != side1D)
            {
                int sideNeighbyside = sideNeigh.Side();
                TPZGeoEl * neigh = sideNeigh.Element();
                
                if(neigh->HasSubElement() || neigh->Dimension() == 1)
                {
                    sideNeigh = sideNeigh.Neighbour();
                    continue;
                }

                //Primeiro estagio da identificacao dos elementos 2D que estao no interior da fratura.
                //brief: primeiramente sao identificados os elementos no interior da fratura que encostam no crack tip.
                //Obs.: A condicao (IsBoundaryMaterial(neigh) == false) eh para excluir os quadrilateros de condicao de contorno.
                if(sd == 2 && neigh->Dimension() == 2 && IsBoundaryMaterial(neigh) == false)
                {
                    TPZVec<REAL> neighCenterQSI(neigh->Dimension()), neighCenterX(3);
                    neigh->CenterPoint(neigh->NSides()-1, neighCenterQSI);
                    neigh->X(neighCenterQSI, neighCenterX);
                    
                    TPZVec<REAL> n0(3), n1(3);
                    fullMesh->NodeVec()[gel->SideNodeIndex(sd, 0)].GetCoordinates(n0);
                    fullMesh->NodeVec()[gel->SideNodeIndex(sd, 1)].GetCoordinates(n1);
                
                    //Como o contorno da fratura foi construido no sentido antihorario no plano x,z (normal Y > 0),
                    //interessam os elementos aa direita do elemento 1D. Portanto eh feito produto vetorial entre os vetores
                    //frac=(n1-n0) e cg_neigh=(cg-n0). O vizinho aa direita apresentarah componente em Y positiva.
                    REAL crossYcomp = n0[2]*n1[0] - n0[0]*n1[2] -
                                        n0[2]*neighCenterX[0] + n1[2]*neighCenterX[0] +
                                        n0[0]*neighCenterX[2] - n1[0]*neighCenterX[2];
            
                    if(crossYcomp > 0.)
                    {
                        neigh->SetMaterialId(__2DfractureMat_inside);
                        
                        std::map<int,TPZFracture2DEl>::iterator edgIt_temp = fracturedElems.find(neigh->Id());
                        if(edgIt_temp != fracturedElems.end())
                        {
                            TPZFracture2DEl neighFractEl = edgIt_temp->second;
                            neighFractEl.RemoveThisEdge(sideNeighbyside);
                            edgIt_temp->second = neighFractEl;
                        }
                        else
                        {
                            TPZFracture2DEl fractEl(neigh);
                            fractEl.RemoveThisEdge(sideNeighbyside);
                            fracturedElems[fractEl.Id()] = fractEl;
                        }
                    }
                }
                
                //Elementos (2D e 3D) que encostam na fratura sao inseridos no mapa fcrackQpointsElementsIds
                std::map< int , std::set<int> >::iterator it = fcrackQpointsElementsIds.find(neigh->Id());
                if(it != fcrackQpointsElementsIds.end())
                {
                    it->second.insert(sideNeighbyside);
                }
                else
                {
                    std::set<int> targetSideId;
                    targetSideId.insert(sideNeighbyside);
                    fcrackQpointsElementsIds[neigh->Id()] = targetSideId;
                }
                
                sideNeigh = sideNeigh.Neighbour();
            }
        }
    }
    
    //capturanto demais elementos no interior da fratura
    std::set<int> finishedFracturedElems;
    while(fracturedElems.size() > 0)
    {
        std::map<int,TPZFracture2DEl>::iterator edgIt = fracturedElems.begin();

        TPZFracture2DEl actEl = edgIt->second;
        
        std::set<int>::iterator sideIt;
        for(sideIt = actEl.fEdge.begin(); sideIt != actEl.fEdge.end(); sideIt++)
        {
            int side = *sideIt;
            TPZGeoElSide actElEdge(actEl.fElem2D,side);
            TPZGeoElSide neighElSide = actElEdge.Neighbour();
            
            bool wellDone = false;
            while(actElEdge != neighElSide && wellDone == false)
            {
                int sideNeighbyside = neighElSide.Side();
                TPZGeoEl * neighEl = neighElSide.Element();
                if(neighEl->Dimension() == 2 && neighEl->MaterialId() == __2DfractureMat_outside && !neighEl->HasSubElement())
                {
                    int neighElId = neighEl->Id();
                    if(finishedFracturedElems.find(neighElId) != finishedFracturedElems.end())
                    {
                        wellDone = true;
                        continue;
                    }
                    std::map<int,TPZFracture2DEl>::iterator edgIt_temp = fracturedElems.find(neighElId);
                    if(edgIt_temp != fracturedElems.end())
                    {
                        TPZFracture2DEl neighFractEl = edgIt_temp->second;
                        neighFractEl.RemoveThisEdge(sideNeighbyside);
                        
                        if(neighFractEl.IsOver())
                        {
                            finishedFracturedElems.insert(neighFractEl.Id());
                            fracturedElems.erase(edgIt_temp);
                        }
                        else
                        {
                            edgIt_temp->second = neighFractEl;
                        }
                    }
                    else
                    {
                        neighEl->SetMaterialId(__2DfractureMat_inside);

                        TPZFracture2DEl fractEl(neighEl);
                        fractEl.RemoveThisEdge(sideNeighbyside);
                        fracturedElems[fractEl.Id()] = fractEl;   
                    }
                    
                    wellDone = true;
                }
                neighElSide = neighElSide.Neighbour();
            }
        }

        finishedFracturedElems.insert(actEl.Id());
        fracturedElems.erase(edgIt);
    }
    
    TurnIntoQuarterPoint(fullMesh);
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFracture::TurnIntoQuarterPoint(TPZGeoMesh * fullMesh)
{
    //List of elements ids that will be refined, so will NOT be turned into quarterpoints (but their sons).
    //These elements must be removed from fcrackQpointsElementsIds map.
    std::set<int> elementsToRemove;
    
    std::set<int>::iterator its;
    std::map< int , std::set<int> >::iterator it;
    for(it = fcrackQpointsElementsIds.begin(); it != fcrackQpointsElementsIds.end(); it++)
    {
        int qpointId = it->first;
        TPZGeoEl * qpointEl = fullMesh->ElementVec()[qpointId];
        
        std::set<int> targetSides = it->second;
        std::set<int> targetSides0D, targetSides1D;
        
        for(its = targetSides.begin(); its != targetSides.end(); its++)
        {
            //separando os lados que encostam no cracktip por suas dimensoes
            int targetSide = *its;
            TPZGeoElSide side(qpointEl,targetSide);
            if(side.Dimension() == 0)
            {
                targetSides0D.insert(targetSide);
            }
            else if(side.Dimension() == 1)
            {
                targetSides1D.insert(targetSide);
            }
        }
        int Ntargets0D = targetSides0D.size();
        int Ntargets1D = targetSides1D.size();
        if(Ntargets0D == 1 && Ntargets1D == 0)//Soh encosta no cracktip por 1 noh apenas
        {
            if(qpointEl->Dimension() == 3)
            {
                qpointEl->SetMaterialId(__3DrockMat_quarterPoint);
            }
//PORQUE ISSO ESTAVA AQUI??? EM UM TESTE, ISSO ATRAPALHOU!!! SERAH QUE EM OUTRO ISSO FAZ FALTA???
//SE NUM FUTURO ISSO NUNCA FEZ FALTA, APAGAR DE VEZ!!!
//Caju, 06 de setembro de 2012.
//            else
//            {
//                elementsToRemove.insert(qpointEl->Id());
//            }
            int targetSide = *(targetSides0D.begin());
            TPZChangeEl::ChangeToQuarterPoint(fullMesh, qpointId, targetSide);
        }
        else if(Ntargets1D == 1 && Ntargets0D == 2)//Soh encosta no cracktip por 1 aresta apenas
        {
            if(qpointEl->Dimension() == 3)
            {
                qpointEl->SetMaterialId(__3DrockMat_quarterPoint);
            }
            int targetSide = *(targetSides1D.begin());
            TPZChangeEl::ChangeToQuarterPoint(fullMesh, qpointId, targetSide);
        }
        else//Encosta no cracktip por aresta(s) e/ou noh(s), portanto precisa refinar para isolar estes lados.
            //Isso porque o quarterpoint contempla o deslocamento dos midnodes para um lado apenas.
            //Obs.: Sao rarissimas excecoes que entram neste escopo.
        {
            bool edgesMustRefine = false;
            TPZVec<int> sidestorefine(qpointEl->NSides(),0.);
            for(int edg = qpointEl->NNodes(); edg < qpointEl->NSides(); edg++)
            {   //se 02 nohs de uma aresta aparecem na lista dos targetSides0D, e a propria aresta nao aparece na lista dos targetSides1D,
                //esta aresta deve ser particionada. Para isso, insere-se esta aresta no vetor sidestorefine
                TPZGeoElSide elSide(qpointEl,edg);
                if(elSide.Dimension() != 1)
                {
                    break;
                }
                int n0 = elSide.SideNodeLocIndex(0);
                int n1 = elSide.SideNodeLocIndex(1);
                if(targetSides0D.find(n0) != targetSides0D.end() &&
                   targetSides0D.find(n1) != targetSides0D.end() &&
                   targetSides1D.find(edg) == targetSides1D.end())
                {
                    sidestorefine[edg] = 1;
                    edgesMustRefine = true;
                }
            }
            if(edgesMustRefine == false)
            {
                //Eh o caso em que temos 2, 3 ou 4 arestas encostando no cracktip, devendo refinar "pelos nohs", e nao pelas arestas.
                //Para isso eh realizado o refinamento baricentrico pela face do plano da fratura
                TPZVec<int> sideNodeIds(targetSides0D.size());
                int pos = 0;
                for(its = targetSides0D.begin(); its != targetSides0D.end(); its++)
                {
                    sideNodeIds[pos] = *its;
                    pos++;
                }
                int sidePlaneFracture = qpointEl->WhichSide(sideNodeIds);
                sidestorefine[sidePlaneFracture] = 1;
            }
            
            TPZAutoPointer<TPZRefPattern> refp = TPZRefPatternTools::PerfectMatchRefPattern(qpointEl, sidestorefine);
            if(refp)
            {
                qpointEl->SetRefPattern(refp);
                TPZVec<TPZGeoEl*> sons;
                qpointEl->Divide(sons);
                
                //realimentando o mapa fcrackQpointsElementsIds com os subelementos que encostam no cracktip
                std::set<int> bySides;
                for(int sn = 0; sn < sons.NElements(); sn++)
                {
                    if(TouchCrackTip(sons[sn],bySides))
                    {
                        fcrackQpointsElementsIds[sons[sn]->Id()] = bySides;
                    }
                }
                elementsToRemove.insert(qpointId);
            }
            else
            {
                std::cout << "\nRefPattern NOT FOUND in " << __PRETTY_FUNCTION__ << std::endl;
                std::cout << "You should create it and add in Refinement Patterns Folder!" << std::endl;
                std::cout << "Open file QpointRefPatternNOTFOUND.vtk in Paraview to see the neighbourhood\n";
                
                std::ofstream outNotFound("QpointRefPatternNOTFOUND.vtk");
                TPZVTKGeoMesh::PrintGMeshVTKneighbourhood(qpointEl->Mesh(), qpointEl->Id(), outNotFound);
                
                DebugStop();
            }   
        }
    }
    
    for(its = elementsToRemove.begin(); its != elementsToRemove.end(); its++)
    {
        //deleting elements that was refined (those sons was turned into quarterpoints, and NOT themselves)
        it = fcrackQpointsElementsIds.find(*its);
        fcrackQpointsElementsIds.erase(it);
    }
}
//------------------------------------------------------------------------------------------------------------

bool TPZPlaneFracture::TouchCrackTip(TPZGeoEl * gel, std::set<int> &bySides)
{
    bySides.clear();
    int nsides = gel->NSides();
    for(int s = 0; s < nsides; s++)
    {
        TPZGeoElSide gelSide(gel,s);
        TPZGeoElSide sideNeigh(gelSide.Neighbour());
        while(gelSide != sideNeigh)
        {
            if(sideNeigh.Element()->MaterialId() == __1DcrackTipMat)
            {
                bySides.insert(s);
            }
            sideNeigh = sideNeigh.Neighbour();
        }
    }
    if(bySides.size() > 0)
    {
        return true;
    }
    else
    {
        return false;
    }
}
//------------------------------------------------------------------------------------------------------------

bool TPZPlaneFracture::IsBoundaryMaterial(TPZGeoEl * gel)
{
    int materialId = gel->MaterialId();
    
    if(materialId == __2DfarfieldMat)
    {
        return true;
    }
    else if(materialId == __2DleftMat)
    {
        return true;
    }
    else if(materialId == __2DrightMat)
    {
        return true;
    }
    else if(materialId == __2DtopMat)
    {
        return true;
    }
    else if(materialId == __2DbottomMat)
    {
        return true;
    }
    
    return false;
}

//** just for visualize given dots in vtk */
void TPZPlaneFracture::InsertDots4VTK(TPZGeoMesh * gmesh, const TPZVec<REAL> &fractureDots)
{
    int nDots = fractureDots.size() / 2;
    int nnodesOriginal = gmesh->NNodes();
	int Qnodes = nnodesOriginal + nDots;
    
	//initializing gmesh->NodeVec()
	gmesh->NodeVec().Resize(Qnodes);
	TPZGeoNode Node;
    TPZVec<REAL> NodeCoord(3);
    TPZVec<int> Topol(1);
    
    int elId = gmesh->NElements();
	for(int n = nnodesOriginal; n < Qnodes; n++)
	{
        Topol[0] = n;
        
        int actDot = n - nnodesOriginal;
        
        NodeCoord[0] = fractureDots[2*actDot];
        NodeCoord[1] = 0.;
        NodeCoord[2] = fractureDots[2*actDot + 1];
        
		Node.SetNodeId(n);
		Node.SetCoord(NodeCoord);
		gmesh->NodeVec()[n] = Node; 
        
        int matPoint = -100;
        new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (elId,Topol, matPoint,*gmesh);
        elId++;
	}
}
//------------------------------------------------------------------------------------------------------------


#define completeCompute

//Just 4 validation of SIF
void TPZPlaneFracture::RunModelProblemForSIFValidation(const TPZVec<REAL> &poligonalChain, std::string vtkFile, int meshDim)
{
    ////CompMesh
    REAL W = 10.;
    REAL H = 14;
    REAL a = 1.;
    REAL sigmaTraction = 5.;
    int porder = 2;
    TPZCompMesh * fractureCMesh = this->GetModelProblemForSIFValidationCompMesh(poligonalChain, porder, meshDim, W, H, a, sigmaTraction);
    
    int neq = fractureCMesh->NEquations();
    std::cout << "Numero de equacoes = " << neq << std::endl;
    
    ////Analysis
    #ifdef completeCompute
        TPZAnalysis an(fractureCMesh);
        
        TPZSkylineStructMatrix skylin(fractureCMesh); //caso simetrico
        TPZStepSolver<STATE> step;
        step.SetDirect(ECholesky);
        
        an.SetStructuralMatrix(skylin);
        an.SetSolver(step);
        an.Run();
    
        ////Post Processing
        TPZManVector<std::string,10> scalnames(0), vecnames(0);
        if(meshDim == 2)
        {
            vecnames.Resize(1);
            vecnames[0] = "displacement";
            
            scalnames.Resize(3);
            scalnames[0] = "SigmaX";
            scalnames[1] = "SigmaY";
            scalnames[2] = "TauXY";
        }
        else if(meshDim == 3)
        {
            scalnames.Resize(6);
            scalnames[0] = "DisplacementX";
            scalnames[1] = "DisplacementY";
            scalnames[2] = "DisplacementZ";
            scalnames[3] = "StressX";
            scalnames[4] = "StressY";
            scalnames[5] = "StressZ";
        }
        
        int div = 0;
        an.DefineGraphMesh(meshDim,scalnames,vecnames,vtkFile);
        an.PostProcess(div,meshDim);
    #endif
    
    TPZVec<REAL> originXYZ(3,0.), direction(3,0.);
    if(meshDim == 2)
    {
        originXYZ[0] = a;
        direction[2] = 1.;
    }
    else if(meshDim==3)
    {
        int POSmiddle1D = int(REAL(fcrackBoundaryElementsIds.NElements())/2. + 0.5);
        int middle1DId = fcrackBoundaryElementsIds[POSmiddle1D];
        
        TPZVec<REAL> originQSI(1,0.);
        TPZGeoEl * gel1D = fractureCMesh->Reference()->ElementVec()[middle1DId];
        gel1D->X(originQSI, originXYZ);
        
        direction[0] = 0.;
        direction[1] = 0.;
        direction[2] = 1.;

        originXYZ[2] = -5.;
    }

    REAL radius = 0.6;
    Path * pathMiddle = new Path(fractureCMesh, originXYZ, direction, radius, meshDim);
    
    JIntegral jInt;
    jInt.PushBackPath(pathMiddle);    
    TPZVec<REAL> Jvector(3);
    
    Jvector = jInt.IntegratePath(0);
    
    if(meshDim==2)
    {
        std::cout << "J = { " << Jvector[0] << " , " << Jvector[1] << " }\n";//Tinha que dar (+/-) 3.e-3
    }
    else
    {
        std::cout << "J = { " << Jvector[0] << " , " << Jvector[2] << " }\n";//Tinha que dar (+/-) 3.e-3
    }
}
//------------------------------------------------------------------------------------------------------------

TPZCompMesh * TPZPlaneFracture::GetModelProblemForSIFValidationCompMesh(const TPZVec<REAL> &poligonalChain, int porder, int meshDim,
                                                                        REAL W, REAL H, REAL a, REAL sigmaTraction)
{
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    TPZCompMesh * cmesh = NULL;
    if(meshDim == 2)
    {
        REAL h_2 = H/2.;
        REAL delta = 0.2;
        
        int ncolsContinuum = int((W/2. - a)/delta + 0.5);
        REAL deltaxContinuum = (W - 2.*a)/2./ncolsContinuum;
        
        int ncolsFracture = int((2.*a)/delta + 0.5);
        REAL deltaxFracture = 2.*a/ncolsFracture;
        
        int nrows = int(h_2/delta + 0.5) + 1;
        REAL deltaY = h_2/(nrows-1);
        
        int ncols = 1 + ncolsContinuum + ncolsFracture + ncolsContinuum;
        
        int NNodes = nrows*ncols;
        
        //initializing gmesh->NodeVec()
        gmesh->NodeVec().Resize(NNodes);
        TPZVec <TPZGeoNode> Node(NNodes);
        for(int r = 0; r < nrows; r++)
        {
            for(int c = 0; c < ncols; c++)
            {
                int n = r*ncols + c;
                gmesh->NodeVec()[n].SetNodeId(n);
                if(c <= ncolsContinuum)
                {
                    REAL x = -W/2. + c*deltaxContinuum;
                    gmesh->NodeVec()[n].SetCoord(0, x);
                }
                else if(c <= ncolsContinuum+ncolsFracture)
                {
                    REAL x = -W/2. + ncolsContinuum * deltaxContinuum + (c-ncolsContinuum)*deltaxFracture;
                    gmesh->NodeVec()[n].SetCoord(0, x);
                }
                else
                {
                    REAL x = -W/2. + ncolsContinuum * deltaxContinuum + ncolsFracture * deltaxFracture + (c-ncolsContinuum-ncolsFracture)*deltaxContinuum;
                    gmesh->NodeVec()[n].SetCoord(0, x);
                }
                
                REAL y = r*deltaY;
                gmesh->NodeVec()[n].SetCoord(1, y);
                
                gmesh->NodeVec()[n].SetCoord(2, 0.);
            }
        }
        
        int domainMat = 1;
        int bottomContinuumMat = -1;
        int bottomFractmat = __1DcrackTipMat;
        int topMat = -3;
        int blockXYZ = -4;
        
        //inserting quadrilaterals
        int elId = 0;
        TPZVec <int> Topol(4);
        for(int r = 0; r < (nrows-1); r++)
        {
            for(int c = 0; c < (ncols-1); c++)
            {
                Topol[0] = ncols*(r+0)+(c+0); Topol[1] = ncols*(r+0)+(c+1); Topol[2] = ncols*(r+1)+(c+1); Topol[3] = ncols*(r+1)+(c+0);
                new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elId,Topol,domainMat,*gmesh);
                elId++;
            }
        }
        
        Topol.Resize(2);
        for(int r = 0; r < (nrows-1); r++)
        {
            for(int c = 0; c < (ncols-1); c++)
            {
                if(r==0)
                {
                    Topol[0] = ncols*(r+0)+(c+0); Topol[1] = ncols*(r+0)+(c+1);
                    if(c < ncolsContinuum)
                    {
                        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elId,Topol,bottomContinuumMat,*gmesh);
                        elId++;
                    }
                    else if(c >= ncolsContinuum+ncolsFracture)
                    {
                        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elId,Topol,bottomContinuumMat,*gmesh);
                        elId++;
                    }
                }
                else if(r == (nrows-2))
                {
                    Topol[0] = ncols*(r+1)+(c+1); Topol[1] = ncols*(r+1)+(c+0);
                    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elId,Topol,topMat,*gmesh);
                    elId++;
                }
            }
        }
        Topol.Resize(1);
        Topol[0] = ncolsContinuum;
        new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (elId,Topol,bottomFractmat,*gmesh);
        elId++;

        Topol[0] =  ncolsContinuum+ncolsFracture;
        new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (elId,Topol,bottomFractmat,*gmesh);
        elId++;
        
        Topol[0] = int(ncols/2.+0.5);
        new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (elId,Topol,blockXYZ,*gmesh);
        elId++;
        
        gmesh->BuildConnectivity();
        
        int firstQuadr = ncolsContinuum-1;
        int lastQuadr = ncolsContinuum + ncolsFracture;

        TPZChangeEl::ChangeToQuarterPoint(gmesh, firstQuadr, 1);
        //
        TPZGeoEl * quad0 = gmesh->ElementVec()[firstQuadr];
        TPZGeoElSide quadSide0(quad0,4);
        TPZGeoElSide edgeElSide0(quadSide0.Neighbour());
        TPZGeoEl * edgeEl0 = edgeElSide0.Element();
        TPZChangeEl::ChangeToQuarterPoint(gmesh, edgeEl0->Id(), 1);

        TPZChangeEl::ChangeToQuarterPoint(gmesh, firstQuadr+1, 0);
        
        TPZChangeEl::ChangeToQuarterPoint(gmesh, lastQuadr-1, 1);
        
        TPZChangeEl::ChangeToQuarterPoint(gmesh, lastQuadr, 0);
        //
        TPZGeoEl * quad3 = gmesh->ElementVec()[lastQuadr];
        TPZGeoElSide quadSide3(quad3,4);
        TPZGeoElSide edgeElSide3(quadSide3.Neighbour());
        TPZGeoEl * edgeEl3 = edgeElSide3.Element();
        TPZChangeEl::ChangeToQuarterPoint(gmesh, edgeEl3->Id(), 0);
 
        
        std::set<int> matIds;
        matIds.insert(bottomFractmat);
        int nRefDir = 3;
        for(int r = 0; r < nRefDir; r++)
        {
            int nEls = gmesh->NElements();
            for(int el = 0; el < nEls; el++)
            {
                TPZGeoEl * qgel = gmesh->ElementVec()[el];
                TPZRefPatternTools::RefineDirectional(qgel, matIds);
            }
        }
        
//        std::ofstream cuco("cuco.vtk");
//        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, cuco);
        ///////////////////////////////////////////////////////////////////////////////////////////////////
        cmesh = new TPZCompMesh(gmesh);
        
        STATE young = 0.29e5;
        STATE poisson = 0.25;
        
        int planeStrain = 0;
//        int planeStress = 1;
        int planeWhat = planeStrain;
        
        TPZMaterial * materialLin = new TPZElasticityMaterial(domainMat, young, poisson, 0., 0., planeWhat);
        cmesh->InsertMaterialObject(materialLin);
        
        ////BCs
        TPZFMatrix<STATE> k(3,3,0.), f(3,1,0.);
        int newmann = 1, mixed = 2;
        
        // farfield traction and simetry
        {
            k(0,0) = 1.E13;
            TPZBndCond * dotBlocked = materialLin->CreateBC(materialLin, blockXYZ, mixed, k, f);
            cmesh->InsertMaterialObject(dotBlocked);
            
            k.Zero();
            k(1,1) = 1.E13;
            TPZBndCond * mixedContinuum = materialLin->CreateBC(materialLin, bottomContinuumMat, mixed, k, f);
            cmesh->InsertMaterialObject(mixedContinuum);

            k.Zero();
            f(1,0) = sigmaTraction;
            TPZBndCond * newmanFarfield = materialLin->CreateBC(materialLin, topMat, newmann, k, f);
            cmesh->InsertMaterialObject(newmanFarfield);
        }
    }
    else if(meshDim == 3)
    {
        cmesh = this->GetFractureCompMesh(poligonalChain,2);
    }
    cmesh->AutoBuild();
    
    return cmesh;
}
//------------------------------------------------------------------------------------------------------------

