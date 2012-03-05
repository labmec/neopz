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
#include "pzbndcond.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"

using namespace pztopology;

/** PUBLIC METHODS */
TPZPlaneFracture::TPZPlaneFracture()
{
    std::cout << "Default constructor would not be used in this class!\n";
    DebugStop();
}
//------------------------------------------------------------------------------------------------------------

TPZPlaneFracture::TPZPlaneFracture(double lw, double bulletDepthIni, double bulletDepthFin, 
                                   TPZVec< std::map<double,double> > & pos_stress)
{
    fpos_stress = pos_stress;
    
    fPlaneMesh = new TPZGeoMesh;
    fFullMesh = new TPZGeoMesh;
        
    std::set<double> espacamentoVerticalTVD;
    std::list<double> espacamentoVerticalDEPTH;
    
    std::map<double,double>::iterator itM;
    
    espacamentoVerticalTVD.insert(0.);
    espacamentoVerticalTVD.insert(lw);
    espacamentoVerticalTVD.insert(bulletDepthIni);
    espacamentoVerticalTVD.insert(bulletDepthFin);
    
    int nstretches = pos_stress.NElements();
    for(int s = 0; s < nstretches; s++)
    {
        for(itM = pos_stress[s].begin(); itM != pos_stress[s].end(); itM++)
        {
            double pos = itM->first;
            espacamentoVerticalTVD.insert(pos);
        }
    }
    
    double pos0 = 0., minEdge = 1000.;
    std::set<double>::iterator itS = espacamentoVerticalTVD.begin(); itS++;
    for(; itS != espacamentoVerticalTVD.end(); itS++)
    {
        double pos1 = *itS;
        double deltaZ = fabs(pos1 - pos0);
        
        int nrows = 1;
        double deltaZused = deltaZ/nrows;
        while(deltaZused > __maxLength)
        {
            nrows++;
            deltaZused = deltaZ/nrows;
        }
        
        minEdge = (deltaZused < minEdge) ? deltaZused : minEdge;
        
        for(int r = 1; r <= nrows; r++)
        {
            double z = pos0 + r*deltaZused;
            espacamentoVerticalTVD.insert(z);
        }
        pos0 = pos1;
    }
    minEdge = (__maxLength < minEdge) ? __maxLength : minEdge;
    fMinimumRadius = 0.9 * (minEdge/__EdgeStretchesQTD)/sqrt(2.0);
    
    #ifdef DEBUG
    if(fMinimumRadius < 0.05)
    {
        std::cout << "Cylinder radius too small on " << __PRETTY_FUNCTION__ << std::endl;
        DebugStop();
    }
    #endif
    
    for(itS = espacamentoVerticalTVD.begin(); itS != espacamentoVerticalTVD.end(); itS++)
    {
        double posDepth = *itS;
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
    
	std::map< int, std::set<double> > elId_TrimCoords;
	std::list< std::pair<int,double> > elIdSequence;
	
	DetectEdgesCrossed(poligonalChain, planeMesh, elId_TrimCoords, elIdSequence);
	
	//Refining auxiliar 1D elements
	TPZVec<TPZGeoEl*> sons;
	std::map< int, std::set<double> >::iterator it;
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

    REAL young = 1000.;
    REAL poisson = 0.2;
    TPZVec<REAL> force(3,0.);

    TPZAutoPointer<TPZMaterial> materialLin = new TPZElasticity3D(__3DrockMat_linear, young, poisson, force);
    cmesh->InsertMaterialObject(materialLin); 
    
    TPZAutoPointer<TPZMaterial> materialQpoint = new TPZElasticity3D(__3DrockMat_quarterPoint, young, poisson, force);
    cmesh->InsertMaterialObject(materialQpoint);
    
    ////BCs    
    TPZFMatrix k(3,3,0.), f(3,1,0.);
    int dirichlet = 0, newmann = 1, mista = 2;
    
    REAL pressureY = 1.;
    
//    //teste 1: compressao uniforme
//    {
//        TPZAutoPointer<TPZMaterial> materialMixedBottom = new TPZElasticity3D(-300, young, poisson, force);
//        k(0,0) = 0.;
//        k(1,1) = 0.;
//        k(2,2) = 1.E12;
//        TPZBndCond * baseDeApoio = new TPZBndCond(materialMixedBottom,__2DbottomMat, mista, k, f);
//        cmesh->InsertMaterialObject(baseDeApoio);
//        
//        TPZAutoPointer<TPZMaterial> materialMixedPoint1 = new TPZElasticity3D(-301, young, poisson, force);
//        k(0,0) = 1.E12;
//        k(1,1) = 1.E12;
//        k(2,2) = 0.;
//        TPZBndCond * pontoDeApoio1 = new TPZBndCond(materialMixedPoint1,__aux0DEl_Mat1, mista, k, f);
//        cmesh->InsertMaterialObject(pontoDeApoio1);
//        
//        k(0,0) = 0.;
//        k(1,1) = 1.E12;
//        k(2,2) = 0.;
//        TPZAutoPointer<TPZMaterial> materialMixedPoint2 = new TPZElasticity3D(-302, young, poisson, force);
//        TPZBndCond * pontoDeApoio2 = new TPZBndCond(materialMixedPoint2,__aux0DEl_Mat2, mista, k, f);
//        cmesh->InsertMaterialObject(pontoDeApoio2);
//
//        TPZAutoPointer<TPZMaterial> materialMixedPoint3 = new TPZElasticity3D(-303, young, poisson, force);
//        TPZBndCond * pontoDeApoio3 = new TPZBndCond(materialMixedPoint3,__aux0DEl_Mat3, mista, k, f);
//        cmesh->InsertMaterialObject(pontoDeApoio3);
//
//        TPZAutoPointer<TPZMaterial> materialMixedPoint4 = new TPZElasticity3D(-304, young, poisson, force);
//        TPZBndCond * pontoDeApoio4 = new TPZBndCond(materialMixedPoint4,__aux0DEl_Mat4, mista, k, f);
//        cmesh->InsertMaterialObject(pontoDeApoio4);
//        
//        TPZAutoPointer<TPZMaterial> materialNewman = new TPZElasticity3D(-400, young, poisson, force);
//        f(0,0) = 0.;
//        f(1,0) = 0.;//Pressao constante e unitaria na direcao Y>0
//        f(2,0) = -pressureY;
//        TPZBndCond * newmanOntop = new TPZBndCond(materialNewman,__2DtopMat, newmann, k, f); 
//        cmesh->InsertMaterialObject(newmanOntop);
//    }
    
//    //teste 2: mov. corpo rigido (rotacao)
    {
        f(0,0) = 0.;
        f(1,0) = 0.;
        f(2,0) = 0.;
        TPZAutoPointer<TPZMaterial> materialMixedPoint1 = new TPZElasticity3D(-301, young, poisson, force);
        TPZBndCond * pontoDeApoio1 = new TPZBndCond(materialMixedPoint1,__aux0DEl_Mat1, dirichlet, k, f);
        cmesh->InsertMaterialObject(pontoDeApoio1);

        f(0,0) = -329.8672286269283;
        f(1,0) = 0.;
        f(2,0) = 0.;
        TPZAutoPointer<TPZMaterial> materialMixedPoint2 = new TPZElasticity3D(-302, young, poisson, force);
        TPZBndCond * pontoDeApoio2 = new TPZBndCond(materialMixedPoint2,__aux0DEl_Mat2, dirichlet, k, f);
        cmesh->InsertMaterialObject(pontoDeApoio2);

        f(0,0) = 0.;
        f(1,0) = 0.;
        f(2,0) = -329.8672286269283;
        TPZAutoPointer<TPZMaterial> materialMixedPoint3 = new TPZElasticity3D(-303, young, poisson, force);
        TPZBndCond * pontoDeApoio3 = new TPZBndCond(materialMixedPoint3,__aux0DEl_Mat3, dirichlet, k, f);
        cmesh->InsertMaterialObject(pontoDeApoio3);

        f(0,0) = -329.8672286269283;
        f(1,0) = 0.;
        f(2,0) = -329.86722862692837;
        TPZAutoPointer<TPZMaterial> materialMixedPoint4 = new TPZElasticity3D(-304, young, poisson, force);
        TPZBndCond * pontoDeApoio4 = new TPZBndCond(materialMixedPoint4,__aux0DEl_Mat4, dirichlet, k, f);
        cmesh->InsertMaterialObject(pontoDeApoio4);
    }

//    //Pressure inside fracture
//    {
//f(0,0) = 0;
//f(1,0) = 0.;
//f(2,0) = 0.;
//
//TPZAutoPointer<TPZMaterial> materialDirichOut = new TPZElasticity3D(-300, young, poisson, force);
//TPZBndCond * dirichetOutside = new TPZBndCond(materialDirichOut,__2DfractureMat_outside, dirichlet, k, f);
//cmesh->InsertMaterialObject(dirichetOutside);
////
//TPZAutoPointer<TPZMaterial> materialDirichTop = new TPZElasticity3D(-301, young, poisson, force);
//TPZBndCond * dirichetTop = new TPZBndCond(materialDirichTop,__2DtopMat, dirichlet, k, f);
//cmesh->InsertMaterialObject(dirichetTop);
//
//TPZAutoPointer<TPZMaterial> materialDirichBottom = new TPZElasticity3D(-302, young, poisson, force);
//TPZBndCond * dirichetBottom = new TPZBndCond(materialDirichBottom,__2DbottomMat, dirichlet, k, f);
//cmesh->InsertMaterialObject(dirichetBottom);
//
//TPZAutoPointer<TPZMaterial> materialDirichFarField = new TPZElasticity3D(-303, young, poisson, force);
//TPZBndCond * dirichetFarfield = new TPZBndCond(materialDirichFarField,__2DfarfieldMat, dirichlet, k, f);
//cmesh->InsertMaterialObject(dirichetFarfield);    
////
////
////
////
//TPZAutoPointer<TPZMaterial> materialNewman = new TPZElasticity3D(-400, young, poisson, force);
//f(0,0) = 0.;
//f(1,0) = pressureY;//Pressao constante e unitaria na direcao Y>0
//        f(2,0) = 0.;
//        TPZBndCond * newmanInside = new TPZBndCond(materialNewman,__2DfractureMat_inside, newmann, k, f); 
//        cmesh->InsertMaterialObject(newmanInside);
//    }
    
    cmesh->AutoBuild();
    
    return cmesh;
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
	TPZStepSolver step;
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

/** PRIVATE METHODS */
//------------------------------------------------------------------------------------------------------------
void TPZPlaneFracture::GeneratePlaneMesh(std::list<double> & espacamento, double lengthFactor)
{
    TPZVec< TPZVec<REAL> > NodeCoord(0);
    int nrows, ncols;
    double Y = 0.;
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

void TPZPlaneFracture::GenerateFullMesh(std::list<double> & espacamento, double lengthFactor)
{
    TPZVec< TPZVec<REAL> > NodeCoord(0);
    int nrows, ncols;
    
    std::list<double>::iterator it = espacamento.end(); it--;
    double tickness = __maxLength;//tickness is the distance between plane of fracture and plane of farfield
    
    int nDirRef = int(log(fabs(tickness)/__maxLength)/log(2.));
    int nLayers = nDirRef + 2;
    
    double Y = 0.;
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
                Topol[1] = ncols*(r+1)+c + l*nNodesByLayer;
                Topol[2] = ncols*(r+1)+c+1 + l*nNodesByLayer;
                Topol[3] = ncols*r+c+1 + l*nNodesByLayer;
                //
                Topol[4] = ncols*r+c + (l+1)*nNodesByLayer;
                Topol[5] = ncols*(r+1)+c + (l+1)*nNodesByLayer;
                Topol[6] = ncols*(r+1)+c+1 + (l+1)*nNodesByLayer;
                Topol[7] = ncols*r+c+1 + (l+1)*nNodesByLayer;
                
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
    TPZGeoEl * pt1a = new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (elId,Topol,__aux0DEl_Mat1,*fFullMesh);
    elId++;
    
    Topol[0] = ncols*(nrows-1);
    TPZGeoEl * pt1b = new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (elId,Topol,__aux0DEl_Mat2,*fFullMesh);
    elId++;
    
    Topol[0] = ncols-1;
    TPZGeoEl * pt2a = new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (elId,Topol,__aux0DEl_Mat3,*fFullMesh);
    elId++;
    
    Topol[0] = ncols*nrows-1;
    TPZGeoEl * pt2b = new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (elId,Topol,__aux0DEl_Mat4,*fFullMesh);
    elId++;
    // ////////////////////////////////////////////////////////MUSTDELETE
    
	fFullMesh->BuildConnectivity();
    fFullMesh->SetMaxElementId(fFullMesh->NElements()-1);
    fFullMesh->SetMaxNodeId(fFullMesh->NNodes()-1);
    
    //    InsertDots4VTK(fullMesh, poligonalChain);    
//    std::ofstream out("FullMesh.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(fFullMesh, out, true);
}
//------------------------------------------------------------------------------------------------------------

void TPZPlaneFracture::GenerateNodesAtPlaneY(std::list<double> & espacamento, double lengthFactor, 
                                             TPZVec< TPZVec<REAL> > & NodeCoord, int & nrows, int & ncols,
                                             double Y)
{
    nrows = espacamento.size();
    
    std::list<double>::iterator it = espacamento.end(); it--;
    double lastPos = fabs(*it);
    
    double Lx = lengthFactor * lastPos;

    int nStretches = 1;
    double deltaX = Lx/nStretches;
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
                                          std::map< int, std::set<double> > &elId_TrimCoords, 
                                          std::list< std::pair<int,double> > &elIdSequence)
{
	int npoints = (poligonalChain.NElements())/2;
	int nelem = planeMesh->NElements();
	
    TPZGeoEl * firstGel = PointElementOnPlaneMesh(0, planeMesh, poligonalChain);

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
	
	double alphaMin;
	bool reachNextPoint;
	int nElsCrossed, thispoint, nextpoint;
	std::map< int, std::set<double> >::iterator it;
	std::set<double> trim;
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
        
        //-----------------------
        double xcoord = poligonalChain[thispoint+0];
        x[0] = xcoord;
        
        double dxcoord = poligonalChain[nextpoint+0] - poligonalChain[thispoint+0];
        dx[0] = dxcoord;
        
        double xnextcoord = poligonalChain[nextpoint+0];
        xNext[0] = xnextcoord;
        //-----------------------
        x[1] = 0.;
        dx[1] = 0.;
        xNext[1] = 0.;
        //-----------------------
        xcoord = poligonalChain[thispoint+1];
        x[2] = xcoord;
        
        dxcoord = poligonalChain[nextpoint+1] - poligonalChain[thispoint+1];
        dx[2] = dxcoord;
        
        xnextcoord = poligonalChain[nextpoint+1];
        xNext[2] = xnextcoord;
        //-----------------------
    
		double norm = 0.;
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
		while(reachNextPoint == false && nElsCrossed < nelem)
		{
			nextGel = CrossToNextNeighbour(gel, x, dx, alphaMin, elId_TrimCoords, elIdSequence, true);
			
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
		nextGel = CrossToNextNeighbour(gel, x, dx, alphaMin, elId_TrimCoords, elIdSequence, false);
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
		nextGel = CrossToNextNeighbour(gel, x, dx, alphaMin, elId_TrimCoords, elIdSequence, true);
		alphaMin = __smallNum;
	}
}
//------------------------------------------------------------------------------------------------------------

TPZGeoEl * TPZPlaneFracture::CrossToNextNeighbour(TPZGeoEl * gel, TPZVec<REAL> &x, TPZVec<REAL> dx, double alphaMin,
												  std::map< int, std::set<double> > &elId_TrimCoords,
                                                  std::list< std::pair<int,double> > &elIdSequence, 
                                                  bool pushback)
{
	bool thereIsAn1DElemAlready;
	int edge;
	std::map< int, std::set<double> >::iterator it;
	std::set<double> trim;
	TPZVec< TPZVec<REAL> > ExactIntersectionPoint, ModulatedIntersectionPoint;
	TPZVec<REAL> qsi1Dvec(1), xCrackBoundary(3);
	TPZVec<int> Topol(2), edgeVec;
	
	TPZGeoMesh * planeMesh = gel->Mesh();
	
	bool haveIntersection = EdgeIntersection(gel, x, dx, edgeVec, ExactIntersectionPoint, ModulatedIntersectionPoint, alphaMin);
	
	if(haveIntersection == false)
	{
		haveIntersection = EdgeIntersection(gel, x, dx, edgeVec, ExactIntersectionPoint, ModulatedIntersectionPoint, 0.);
        
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
            double coordM = ModulatedIntersectionPoint[edg][c];
			xLin[c] = coordM;
		}
		double qsi1D = LinearComputeXInverse(xLin, n0, n1);
		
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
	while(neighEdge != gelEdge)
	{
		if(neighEdge.Element()->Dimension() == 2)
		{
			break;
		}
		neighEdge = neighEdge.Neighbour();
	}
	gel = neighEdge.Element();
	
	return gel;
}
//------------------------------------------------------------------------------------------------------------

bool TPZPlaneFracture::EdgeIntersection(TPZGeoEl * gel, TPZVec<REAL> &x, TPZVec<REAL> &dx, TPZVec<int> &edge,
										TPZVec< TPZVec<REAL> > &ExactIntersect,
                                        TPZVec< TPZVec<REAL> > &ModulatedIntersect, double alphaMin)
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
	std::list<double> edgeNorm;
	/** <double: alpha , <double: alphaNodemod, double: alphaNodesmooth, double: norm, int: edge that intersect> > */
    std::map<double, TPZVec<REAL> > alpha;
    
	double alphaX, alphaNodemod, alphaNodesmooth, norm;
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
		
		alphaX = ComputeAlphaX(x, dx, node[n], dnode[n]);
		alphaNodemod = ComputeAlphaNode(x, dx, node[n], dnode[n], norm, true, false);
		alphaNodesmooth = ComputeAlphaNode(x, dx, node[n], dnode[n], norm, IsNearNode, true);
		if(alphaX >= alphaMin && alphaNodesmooth >= 0. && alphaNodesmooth <= norm)
		{
			int thisEdge = n+ncnodes;
			TPZVec<REAL> someData(4);
			someData[0] = alphaNodemod; someData[1] = alphaNodesmooth; someData[2] = norm; someData[3] = thisEdge;
			alpha[alphaX] = someData;
		}
	}
	
	// a aresta que serah interseccionada serah a que tiver menor alphaX não negativo, i.e.: o primeiro par do mapa alpha!
	std::map<double, TPZVec<REAL> >::iterator it = alpha.begin();
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

TPZGeoEl * TPZPlaneFracture::PointElementOnPlaneMesh(int p, TPZGeoMesh * planeMesh, const TPZVec<REAL> &poligonalChain)
{
	TPZManVector<REAL,3> x(3), qsi2D(2,0.);
	
    x[0] = poligonalChain[2*p+0];
    x[1] = 0.;
    x[2] = poligonalChain[2*p+1];
    
	TPZGeoEl * gel = NULL;
	int nelem = planeMesh->NElements();
	for(int el = 0; el < nelem; el++)//hunting the first element (that contains the first point)
	{
		TPZGeoEl * firstGel = planeMesh->ElementVec()[el];
		if(firstGel->Dimension() != 2)
		{
			continue;
		}
		if(firstGel->ComputeXInverse(x, qsi2D))
		{
			gel = firstGel;
			break;
		}
	}
    
#ifdef DEBUG
    if(!gel)
    {
        std::cout << "Point DO NOT belong to plane mesh domain on " << __PRETTY_FUNCTION__ << std::endl;
        DebugStop();
    }
#endif
	
	return gel;
}
//------------------------------------------------------------------------------------------------------------

double TPZPlaneFracture::ComputeAlphaNode(TPZVec<REAL> &x, TPZVec<REAL> &dx,
                                          TPZVec<REAL> &node, TPZVec<REAL> &dnode,
                                          double norm, bool modulate, bool smooth)
{
	double fractionNumQ =	dx[2]*node[0] - dx[0]*node[2] - dx[2]*x[0] + dx[0]*x[2];
	
	double fractionDenomQ =	dnode[2]*dx[0] - dnode[0]*dx[2];
	
	double alphaNode = -1.;
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

double TPZPlaneFracture::ComputeAlphaX(TPZVec<REAL> &x, TPZVec<REAL> &dx,
                                       TPZVec<REAL> &node, TPZVec<REAL> &dnode)
{	
	//computing alpha (dx vector multiplier to intersect element edges)
	double fractionNumP   =	dnode[2]*node[0] - dnode[0]*node[2] -
	dnode[2]*x[0] + dnode[0]*x[2];
	
	double fractionDenomP =	dnode[2]*dx[0] - dnode[0]*dx[2];
	
	double alphaX = -1.;
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

double TPZPlaneFracture::LinearComputeXInverse(TPZVec<REAL> x, TPZVec<REAL> n0, TPZVec<REAL> n1)
{
	double dL = 0., L = 0.;
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
	
	double qsi = -1. + dL/L*2.;
	
	return qsi;
}
//------------------------------------------------------------------------------------------------------------

TPZAutoPointer<TPZRefPattern> TPZPlaneFracture::Generate1DRefPatt(std::set<double> &TrimCoord)
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
	std::set<double>::iterator it;
	for(it = TrimCoord.begin(); it != TrimCoord.end(); it++)
	{
		double coord = *it;
		
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
                                            std::list< std::pair<int,double> > &elIdSequence,
                                            TPZVec<REAL> &poligonalChainUpdated)
{
	int nptos = elIdSequence.size();
	poligonalChainUpdated.Resize(2*nptos);

	TPZVec<REAL> qsi1Dvec(1), ptoCoord(3);
	int el1Did, posX, posZ, p = 0;
	double qsi1D;
	
	std::list< std::pair<int,double> >::iterator it;
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
                                             std::list< std::pair<int,double> > &elIdSequence)
{
    fcrackBoundaryElementsIds.Resize(0);
	TPZVec<REAL> qsi0vec(1), qsi1vec(1), node0coord(3), node1coord(3);
	TPZVec<int> Topol(2);
	int el0id, el1id, n0, n1;
	double qsi0, qsi1;
	std::list< std::pair<int,double> >::iterator crackit0, crackit1, crackitEnd;
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
                    double crossYcomp = n0[2]*n1[0] - n0[0]*n1[2] -
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
            else
            {
                elementsToRemove.insert(qpointEl->Id());
            }
            int targetSide = *(targetSides0D.begin());
            //TPZChangeEl::ChangeToQuarterPoint(fullMesh, qpointId, targetSide);
        }
        else if(Ntargets1D == 1 && Ntargets0D == 2)//Soh encosta no cracktip por 1 aresta apenas
        {
            if(qpointEl->Dimension() == 3)
            {
                qpointEl->SetMaterialId(__3DrockMat_quarterPoint);
            }
            int targetSide = *(targetSides1D.begin());
            //TPZChangeEl::ChangeToQuarterPoint(fullMesh, qpointId, targetSide);
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
                elementsToRemove.insert(qpointEl->Id());
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

/** to Phil - TESTAR NOVAMENTE DEPOIS DE TANTAS ALTERACOES!!! */
void TPZPlaneFracture::GetSidesCrossedByPoligonalChain(const TPZVec<REAL> &poligonalChain,
                                                       std::list< std::pair<TPZGeoElSide,double> > &sidesCrossed)
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
    
    sidesCrossed.clear();
	
	TPZGeoMesh * planeMesh = new TPZGeoMesh(*fPlaneMesh);
	
	std::map< int, std::set<double> > elId_TrimCoords;
	std::list< std::pair<int,double> > elIdSequence;
	
	DetectEdgesCrossed(poligonalChain, planeMesh, elId_TrimCoords, elIdSequence);
    
    std::list< std::pair<int,double> >::iterator it;
    
    for(it = elIdSequence.begin(); it != elIdSequence.end(); it++)
    {
        int gel1D_Id = it->first;
        double qsiNeigh, qsi = it->second;
        
        TPZGeoEl * gel1D = planeMesh->ElementVec()[gel1D_Id];
        TPZGeoElSide gel1DSide(gel1D,2);
        TPZGeoElSide neighSide = gel1DSide.Neighbour();
        while(neighSide != gel1DSide)
        {
            TPZTransform tr = gel1DSide.NeighbourSideTransform(neighSide);
            qsiNeigh = tr.Mult()(0,0)*qsi;
            
            sidesCrossed.push_back(std::make_pair(neighSide, qsiNeigh));
            neighSide = neighSide.Neighbour();
        }
    }
    
    planeMesh = NULL;
    delete planeMesh;
}
//------------------------------------------------------------------------------------------------------------

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

