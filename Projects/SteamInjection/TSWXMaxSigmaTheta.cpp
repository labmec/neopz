//$Id: TSWXMaxSigmaTheta.cpp,v 1.5 2009-07-03 19:43:07 caju Exp $

/*
 *  TSWXMaxSigmaTheta.cpp
 *  TSWXMaxSigmaTheta
 *
 *  Created by caju on 6/18/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "TSWXMaxSigmaTheta.h"
#include "TPZRefPatternTools.h"
#include "tpzgeoelrefpattern.h"


//public methods:
TSWXMaxSigmaTheta::TSWXMaxSigmaTheta()
{
}
//

TSWXMaxSigmaTheta::~TSWXMaxSigmaTheta()
{

}
//

TSWXMaxSigmaTheta::TSWXMaxSigmaTheta(double rw, double h, double E, double nu)
{
	SetData(rw, h, E, nu);
}

#define NO_Refinement
void TSWXMaxSigmaTheta::SetData(double rw, double h, double E, double nu)
{
	#ifndef NO_Refinement
	fDirectNdiv = (h < 5) ? 10 : 11;
	#else
	fDirectNdiv = 0;
	#endif

	fUnifNDiv = 0;
	fDivMAX = 6;
	fOrder = 4;

	fH = 0.;
	fB = 0.;

	fb_sol.clear();

	//applied forces
	fDistrLeftDown = 0.;
	fDistrLeftUp = 0.;
	fDistrMiddleUp = 0.;
	fDistrRightUp = 0.;
	fDistrSurface = 0.;

	frw = fabs(rw);
	fh = fabs(h);
	fE = fabs(E);
	fnu = fabs(nu);
}

void TSWXMaxSigmaTheta::ComputeMaxSigmaTheta(double t, double b, double DistrRightDown)
{
	//GeoMesh & CompMesh
	CreateGeoMesh(b);	
	CreateCompMesh(DistrRightDown);

	//Processing
	TPZAnalysis an(fCmesh);
	TPZSkylineStructMatrix full(fCmesh);
	an.SetStructuralMatrix(full);
	TPZStepSolver step;
	step.SetDirect(ECholesky);
	an.SetSolver(step);
	an.Run();
	
	//Computing Maximum SigmaTheta
	TPZGeoElSide elSide(fGmesh->ElementVec()[0], 0);
	TPZStack< TPZGeoElSide > subelements;
	TPZCompElSide celSide = elSide.Reference();
	int stop = 0;
	while (!celSide.Exists() && stop < 1000)
	{
		elSide.GetSubElements2(subelements);
		for(int i = 0; i < subelements.NElements(); i++)
		{
			elSide = subelements[i];
			celSide = elSide.Reference();
			if(celSide.Exists()) break;
		}
		stop++;
	}
	if(celSide.Exists())
	{
		TPZVec<REAL> qsi(2);
		celSide.Reference().Element()->CenterPoint(celSide.Side(),qsi);
		int var = fCmesh->ElementVec()[0]->Material()->VariableIndex("Sigmatt");
		TPZVec<REAL> sol(1);
		celSide.Element()->Solution(qsi, var, sol);
		fb_sol[t] = sol[0];
		
		//Solution visualization on Paraview (VTK)
//#define VTK
#ifdef VTK
		TPZVec<std::string> scalnames(3), vecnames(0);
		scalnames[0] = "Sigmarr";
		scalnames[1] = "Sigmazz";
		scalnames[2] = "Sigmatt";
		
		std::string plotfile("saida.vtk");
		const int dim = 2;
		an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
		an.PostProcess(0);
		std::ofstream out("malha.txt");
		an.Print("nothing",out);
#endif

//		std::ofstream outM("mat.txt");
//		this->PrintMathematica(outM);
	}
	else
	{
		std::cout << "CelSide doesnt exists on " << __PRETTY_FUNCTION__ << std::endl;
		DebugStop();
	}
}
//

//private methods:
void TSWXMaxSigmaTheta::CreateGeoMesh(double b)
{
	if(frw < 1.e-8) return;
	
	b = fabs(b);
	
	double brw = b / frw;
	double alpha   =	(brw < 120000.) ? 39.9920168733301 + 0.015966419249637404*brw - 3.3182168190175215e-7*pow(brw,2) + 3.735927711340961e-12*pow(brw,3) - 
							2.1485282981404827e-17*pow(brw,4) + 4.881566116092741e-23*pow(brw,5) + 2.372828046884787e-30*pow(brw,6) : 400.;
	fH = (alpha/2.)*fh;
	fB = alpha*fh;
	
	fB = (b > fB) ? 2.*b : fB;
	
	fGmesh = new TPZGeoMesh;
	
	const double d1 = (b < fh) ? b : fh;
	const double d2 = ((fB-b) < (fH-fh)) ? (fB-b) : (fH-fh);
	const double d = (d1 < d2) ? d1 : d2;
	
	const int h_ndiv = (int(fh/d+0.5) < fDivMAX) ? int(fh/d+0.5) : fDivMAX;
	const int hH_ndiv = (int((fH-fh)/d+0.5) < fDivMAX) ? int((fH-fh)/d+0.5) :fDivMAX;
	const int b_ndiv = (int(b/d+0.5) < fDivMAX) ? int(b/d+0.5) : fDivMAX;
	const int bB_ndiv = (int((fB-b)/d+0.5) < fDivMAX) ? int((fB-b)/d+0.5) : fDivMAX;
	
	const int nrows = (h_ndiv + hH_ndiv + 1);
	const int ncols = (b_ndiv + bB_ndiv + 1);
	
	const int Qnodes =  nrows*ncols;
	TPZVec < TPZVec <REAL> > NodeCoord(Qnodes);
	for(int i = 0; i < Qnodes; i++) NodeCoord[i].Resize(3);
	
	//setting nodes coordinates
	int nodeId = 0;
	for(int r = 0; r < nrows; r++)
	{
		for(int c = 0; c < ncols; c++)
		{
			if(c <= b_ndiv)
			{
				NodeCoord[nodeId][0] = frw + c*(b/b_ndiv);
			}
			else
			{
				NodeCoord[nodeId][0] = frw + b + (c-b_ndiv)*(fB-b)/bB_ndiv;
			}
			if(r <= h_ndiv)
			{
				NodeCoord[nodeId][1] = r*(fh/h_ndiv);
			}
			else
			{
				NodeCoord[nodeId][1] = fh + (r-h_ndiv)*(fH-fh)/hH_ndiv;
			}
			NodeCoord[nodeId][2] = 0.;
			
			nodeId++;
		}
	}
	
	//setting nodes objects in gmesh
	fGmesh->NodeVec().Resize(Qnodes);
	TPZVec <TPZGeoNode> Node(Qnodes);
	for(int n = 0; n < Qnodes; n++)
	{
		Node[n].SetNodeId(n);
		Node[n].SetCoord(NodeCoord[n]);
		fGmesh->NodeVec()[n] = Node[n]; 
	}
	
	//Quadriláteros
	int elId = 0;
	TPZVec <int> Topol(4);
	for(int r = 0; r < (nrows-1); r++)
	{
		for(int c = 0; c < (ncols-1); c++)
		{
			Topol[0] = ncols*r+c; Topol[1] = ncols*r+c+1; Topol[2] = ncols*(r+1)+c+1; Topol[3] = ncols*(r+1)+c;
			new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elId,Topol,matElId,*fGmesh);
			elId++;
		}
	}
	
	//setting bc's and materials id's
	Topol.Resize(2);
	
	//BCs - TOP and BOTTOM
	for(int c = 0; c < (ncols-1); c++)
	{
		Topol[0] = c;
		Topol[1] = c+1;
		new TPZGeoElRefPattern<TPZGeoLinear> (elId,Topol,matBOTTOMid,*fGmesh);
		elId++;
		
		Topol[0] = ncols*(nrows-1)+c;
		Topol[1] = ncols*(nrows-1)+c+1;
		new TPZGeoElRefPattern<TPZGeoLinear> (elId,Topol,matTOPid,*fGmesh);
		elId++;
	}
	
	//BCs - LEFT and RIGHT
	for(int r = 0; r < (nrows-1); r++)
	{
		if(r < h_ndiv)
		{
			Topol[0] = r*ncols;
			Topol[1] = (r+1)*ncols;
			new TPZGeoElRefPattern<TPZGeoLinear> (elId,Topol,matLeftDOWNid,*fGmesh);
			elId++;
			
			Topol[0] = (ncols-1)+r*ncols;
			Topol[1] = (ncols-1)+(r+1)*ncols;
			new TPZGeoElRefPattern<TPZGeoLinear> (elId,Topol,matRightDOWNid,*fGmesh);
			elId++;
		}
		else
		{
			Topol[0] = r*ncols;
			Topol[1] = (r+1)*ncols;
			new TPZGeoElRefPattern<TPZGeoLinear> (elId,Topol,matLeftUPid,*fGmesh);
			elId++;
			
			Topol[0] = (ncols-1)+r*ncols;
			Topol[1] = (ncols-1)+(r+1)*ncols;
			new TPZGeoElRefPattern<TPZGeoLinear> (elId,Topol,matRightUPid,*fGmesh);
			elId++;
		}
	}
	
	//BCs - MIDDLE
	for(int r = 0; r < (nrows-1); r++)
	{
		if(r < h_ndiv)
		{
			Topol[0] = b_ndiv+r*ncols;
			Topol[1] = b_ndiv+(r+1)*ncols;
			new TPZGeoElRefPattern<TPZGeoLinear> (elId,Topol,matMiddleDOWNid,*fGmesh);
			elId++;
		}
		else
		{
			Topol[0] = b_ndiv+r*ncols;
			Topol[1] = b_ndiv+(r+1)*ncols;
			new TPZGeoElRefPattern<TPZGeoLinear> (elId,Topol,matMiddleUPid,*fGmesh);
			elId++;
		}
	}
	
	fGmesh->BuildConnectivity();
	
	//Uniform Refinement
	for(int h = 0; h < fUnifNDiv; h++)
	{
		int nel = fGmesh->NElements();
		for ( int iref = 0; iref < nel; iref++ )
		{
			TPZVec<TPZGeoEl*> filhos;
			TPZGeoEl * gelP1 = fGmesh->ElementVec()[iref];
			if(!gelP1) continue;
			gelP1->Divide(filhos);
		}
	}
	
	std::set<int> matidsB, matidsL;
	matidsB.insert(matBOTTOMid);
	matidsL.insert(matLeftDOWNid);
	matidsL.insert(matLeftUPid);
	matidsL.insert(matMiddleDOWNid);
	matidsL.insert(matMiddleUPid);
	matidsL.insert(matRightDOWNid);
	matidsL.insert(matRightUPid);
	
	//Directional Refinement
	if(d/fB > 0.20) fDirectNdiv *= 2;
	for(int h = 0; h < fDirectNdiv; h++)
	{
		int nel = fGmesh->NElements();
		for ( int iref = 0; iref < nel; iref++ )
		{
			TPZVec<TPZGeoEl*> filhos;
			TPZGeoEl * gelP1 = fGmesh->ElementVec()[iref];
			if(!gelP1) continue;
			TPZRefPatternTools::RefineDirectional(gelP1,matidsB);
		}
	}
	for(int h = 0; h < fDirectNdiv; h++)
	{
		int nel = fGmesh->NElements();
		for ( int iref = 0; iref < nel; iref++ )
		{
			TPZVec<TPZGeoEl*> filhos;
			TPZGeoEl * gelP1 = fGmesh->ElementVec()[iref];
			if(!gelP1) continue;
			TPZRefPatternTools::RefineDirectional(gelP1,matidsL);
		}
	}
}
//

void TSWXMaxSigmaTheta::CreateCompMesh(double DistrRightDown)
{
	REAL   fx = fDistrSurface, fy = 0.; // fy = Gravity
	double mi = fE/(2.*(1.+fnu));
	
	TPZAutoPointer<TPZMaterial> mat = new TPZElasticityAxiMaterial(matElId, fE, fnu, fx, fy);
	
	TPZManVector<REAL> Orig(3);		Orig[0] = 0.;	Orig[1] = 0.;	Orig[2] = 0.;
	TPZManVector<REAL> AxisZ(3);		AxisZ[0] = 0.;	AxisZ[1] = 1.;	AxisZ[2] = 0.;
	TPZManVector<REAL> AxisR(3);		AxisR[0] = 1.;	AxisR[1] = 0.;	AxisR[2] = 0.;
	TPZElasticityAxiMaterial * aximat = dynamic_cast<TPZElasticityAxiMaterial*>(mat.operator->());
	aximat->SetOrigin(Orig, AxisZ, AxisR);

	///Computational Mesh
	TPZCompEl::SetgOrder(fOrder);
	fCmesh = new TPZCompMesh(fGmesh);
	fCmesh->InsertMaterialObject(mat);
	{
		///Boundary Conditions
		TPZFMatrix DistrLEFTDown1(2,2,0.), DistrLEFTDown2(2,1,0.);
		DistrLEFTDown2(0,0) = fDistrLeftDown;
		TPZAutoPointer<TPZMaterial> DistrDownLEFTBC = mat->CreateBC(mat, matLeftDOWNid, newmann, DistrLEFTDown1, DistrLEFTDown2);
		fCmesh->InsertMaterialObject(DistrDownLEFTBC);
		
		TPZFMatrix DistrLEFTUp1(2,2,0.), DistrLEFTUp2(2,1,0.);
		DistrLEFTUp2(0,0) = fDistrLeftUp;
		TPZAutoPointer<TPZMaterial> DistrUpLEFTBC = mat->CreateBC(mat, matLeftUPid, newmann, DistrLEFTUp1, DistrLEFTUp2);
		fCmesh->InsertMaterialObject(DistrUpLEFTBC);
		
		TPZFMatrix DistrRIGHTDown1(2,2,0.), DistrRIGHTDown2(2,1,0.);
		DistrRIGHTDown1(0,0) = 2.*mi/(frw+fB);
		DistrRIGHTDown2(0,0) = DistrRightDown;
		TPZAutoPointer<TPZMaterial> DistrDownRIGHTBC = mat->CreateBC(mat, matRightDOWNid, mista, DistrRIGHTDown1, DistrRIGHTDown2);
		fCmesh->InsertMaterialObject(DistrDownRIGHTBC);
		
		TPZFMatrix DistrRIGHTUp1(2,2,0.), DistrRIGHTUp2(2,1,0.);
		DistrRIGHTUp1(0,0) = 2.*mi/(frw+fB);
		DistrRIGHTUp2(0,0) = fDistrRightUp;
		TPZAutoPointer<TPZMaterial> DistrUpRIGHTBC = mat->CreateBC(mat, matRightUPid, mista, DistrRIGHTUp1, DistrRIGHTUp2);
		fCmesh->InsertMaterialObject(DistrUpRIGHTBC);
		
		TPZFMatrix DistrMIDDLEDown1(2,2,0.), DistrMIDDLEDown2(2,1,0.);
		DistrMIDDLEDown2(0,0) = fDistrMiddleDown;
		TPZAutoPointer<TPZMaterial> DistrDownMIDDLEBC = mat->CreateBC(mat, matMiddleDOWNid, mista, DistrMIDDLEDown1, DistrMIDDLEDown2);
		fCmesh->InsertMaterialObject(DistrDownMIDDLEBC);
		
		TPZFMatrix DistrMIDDLEUp1(2,2,0.), DistrMIDDLEUp2(2,1,0.);
		DistrMIDDLEUp2(0,0) = fDistrMiddleUp;
		TPZAutoPointer<TPZMaterial> DistrUpMIDDLEBC = mat->CreateBC(mat, matMiddleUPid, mista, DistrMIDDLEUp1, DistrMIDDLEUp2);
		fCmesh->InsertMaterialObject(DistrUpMIDDLEBC);
		
		//CC TOP (Mista) - impedindo desloc em y 
		double BigNum = TPZMaterial::gBigNumber;
		
		//		TPZFMatrix DistrTOP1(2,2,0.), DistrTOP2(2,1,0.);
		//		DistrTOP1(1,1) = BigNum;
		//		TPZAutoPointer<TPZMaterial> DistrTOPBC = mat->CreateBC(mat, matTOPid, mista, DistrTOP1, DistrTOP2);
		//		fCmesh->InsertMaterialObject(DistrTOPBC);
		
		//CC BOTTOM (Mista) - impedindo desloc em y 
		TPZFMatrix DistrBOTTOM1(2,2,0.), DistrBOTTOM2(2,1,0.);
		DistrBOTTOM1(1,1) = BigNum;
		TPZAutoPointer<TPZMaterial> DistrBOTTOMBC = mat->CreateBC(mat, matBOTTOMid, mista, DistrBOTTOM1, DistrBOTTOM2);
		fCmesh->InsertMaterialObject(DistrBOTTOMBC);
	}
	
	fCmesh->AutoBuild();
}
//

void TSWXMaxSigmaTheta::PrintMathematica(ostream &out)
{
	if(fb_sol.size() == 0) return;
	out.precision(8);
	
	map<double,double>::iterator p, q;
	out << "bvsSigmaThetaMAX={" << endl;
	for(p = fb_sol.begin(); p != fb_sol.end(); p++)
	{
		q = p; q++;
		out << "{" <<  p->first << " , " << p->second << "}";
		if(q != fb_sol.end()) out << ",";
		out << endl;
	}
	out << "};" << endl << endl << "ListPlot[bvsSigmaThetaMAX, Joined -> True, Filling->Axis, PlotMarkers -> Automatic, PlotRange -> All, AxesOrigin -> {0, 0}]";
}
//
