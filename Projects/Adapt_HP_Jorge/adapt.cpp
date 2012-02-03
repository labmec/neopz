/***************************************************************************
 *   Copyright (C) 2009 by Olivier Roussel   *
 *   o_roussel@yahoo.fr   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "adapt.h"
#include "TPZFakeFunction.h"
#include "pzexplfinvolanal.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.olivier.adapt"));
#endif

extern int gLMax;


int mainAdapt ( int argc, char *argv[] )
{
	// --- Start -----------------------------------------------
	
	cout << "lara: starting ..." << endl;
	
	
	// --- Define the coarse mesh ------------------------------
	
	// Read coarse mesh file
	TPZReadMeshHR Dummy ( "LaraMesh.txt" );
	
	// Create an object to store the geometrical mesh
	TPZGeoMesh * gDummyMesh = NULL;
	
	// Initialize object from read data
	gDummyMesh = Dummy.readGeoMesh();
	
	
	// --- Define the refinement patterns for tetrahedra -------
	
	// Read the refinement pattern
	string tetraFile ( "LaraRefineTetrahedron.txt" );
	
	// Create and initialize an object to store the pattern
	TPZRefPattern *tetraPattern = new TPZRefPattern ( *gDummyMesh);
	tetraPattern->SetName(tetraFile );
	
	// Include the pattern info into the mesh
	tetraPattern->InsertPermuted();
	
	// Remove this object
	delete tetraPattern;
	
	
	// --- Define the refinement patterns for triangles -------
	
	// Read the refinement pattern
	string triangleFile ( "LaraRefineTriangle.txt" );
	
	// Create and initialize an object to store the pattern
	TPZRefPattern * trianglePattern = new TPZRefPattern ( *gDummyMesh);
	trianglePattern->SetName(triangleFile );
	
	// Include the pattern info into the mesh
	trianglePattern->InsertPermuted();
	
	// Remove this object
	delete trianglePattern;
	
	// --- Divide elements ------------------------------------
	
	// Take the first element of the mesh and write it into geo
	TPZGeoEl *geo = gDummyMesh->ElementVec()[ 0 ];
	
	// Create vector of geometric elements
	// This vector returns the pointers to the children
	TPZVec < TPZGeoEl * > childVector;
	
	// Take the map of every refinement pattern for a tetrahedron defined into the mesh
	// the map below provides a relation between an index and a pointer to an object of type refinement pattern
	map < int, TPZAutoPointer < TPZRefPattern > > mapOfTetraRefPattern; //=  gDummyMesh->RefPatternList ( eltype );
	map < int, TPZAutoPointer < TPZRefPattern > >::iterator it;
	
	// this is the first refinement pattern available = uniform pattern (include 2 pyramids)
	// NOTE: this is not what we want
	it = mapOfTetraRefPattern.begin();
	
	// incresing the iterator, we have a pointer to the second refinement pattern available
	// i.e. the refinement pattern defined into the file "LaraRefineTetrahedron.txt"
	it++;
	
	// *(it->first) = index of refinement pattern;
	// it->second = pointer to refinement pattern;
	// TPZAutoPointer is an auxiliary structure to provide garbage collector facility
	TPZAutoPointer < TPZRefPattern > laraRefinePattern = it->second;
	geo->SetRefPattern ( laraRefinePattern );
	//@NOTE: must define the same refinement pattern to all children of the element!
	
	// --- Create an computational mesh ----------------------
	// Computational mesh implements the interpolation space over de partition of the domain defined in geometric mesh
	// The parameter is the geometric mesh
	TPZCompMesh * cDummyMesh = new TPZCompMesh ( gDummyMesh );
	
	// Specify the use of discontinuous interpolation space
	cDummyMesh->SetAllCreateFunctionsDiscontinuous();
	
	// --- Material / Equation data:
	// first parameter is the material identifier ( identifier is defined in element sectio of file LaraMesh.txt )
	// second parameter is the parameter gamma ( I think that is some parameter of the equation - please confirm it!)
	TPZAutoPointer < TPZMaterial > eulerMaterial = new TPZEulerEquation ( 1, 1.4 );
	
	// Insert the created material into the material data base of the computational mesh
	cDummyMesh->InsertMaterialObject ( eulerMaterial );
	
	//@NOTE: must create the boundary condition material...
	// first parameter is the material wich boundary condition will be applied
	// second parameter is material identifier ( as specified into the element definition section of the file LaraMesh.txt
	// third and fourth parameters: depends on the material... Ask to Tiago
	TPZFMatrix val1( 1, 1, 0. ), val2( 1, 1, 0. );
	TPZAutoPointer < TPZMaterial > boundaryCondition = eulerMaterial->CreateBC ( eulerMaterial, -1, 1, val1, val2 );
	
	// Method to create the interpolation spaces
	cDummyMesh->AutoBuild();
	
	int ref = 0;
	int nref = 5;
	for ( ref=0; ref < nref; ref++ )
	{
		// Creates an analysis
		TPZAnalysis eulerAnalysis ( cDummyMesh );
		
		TPZSkylineStructMatrix strskyl ( cDummyMesh );
		eulerAnalysis.SetStructuralMatrix ( strskyl );
		TPZStepSolver * direct = new TPZStepSolver;
		direct->SetDirect ( ECholesky );
		eulerAnalysis.SetSolver ( * direct );
		delete direct;
		
		eulerAnalysis.Run ();
		eulerAnalysis.PostProcess ( 0, 2 );
		
		TPZVec < EAdaptElementAction > DivideOrCoarsen;
		// Call the error evaluation and fill the decision vector for each element ( divide - coarse - none )
		ErrorEstimation ( * cDummyMesh, DivideOrCoarsen,0.0001 );
		
		AdaptMesh ( * cDummyMesh,  DivideOrCoarsen, laraRefinePattern );
		ofstream outStreamAdapted ( "refined_Lara.txt" );
		//print computational mesh;
		cDummyMesh->Print ( outStreamAdapted );
		//print geometric mesh
		gDummyMesh->Print ( outStreamAdapted );
	}
	// --- Finish ---------------------------------------------
	cout << "lara: goodbye" << endl;
	int q;
	cin >> q;
	return EXIT_SUCCESS;
}

TPZAutoPointer < TPZRefPattern > GetUsedRefinementPattern ( TPZCompMesh * CMesh )
{
	TPZGeoMesh * gmesh = CMesh->Reference();
	int iel, nel = gmesh->NElements();
	TPZAutoPointer < TPZRefPattern > refpattern;
	for ( iel = 0; iel < nel; iel++ )
	{
		TPZGeoEl * gel = gmesh->ElementVec()[ iel ];
		if ( !gel ) continue;
		if ( gel->Type() == ETetraedro )
		{
			refpattern = gel->GetRefPattern();
			if ( refpattern ) return refpattern;
		}
	}
	return refpattern;
}

/**
 * Public interface to retrieve the adapted mesh
 */
void GetAdaptedMesh ( TPZCompMesh * CMesh, double Epsl )
{
	
	cout << "\ngLMax = " << gLMax << "\n";
	cout.flush();
	
	const int nelem = CMesh->NElements();
	TPZManVector < EAdaptElementAction,1000 > DivideOrCoarsen ( nelem, ENone );
	// Call the error evaluation and fill the decision vector for each element ( divide - coarse - none )
	ErrorEstimation ( * CMesh, DivideOrCoarsen, Epsl );
	cout << " GetUsedRefinementPattern " << endl;cout.flush();
	TPZAutoPointer < TPZRefPattern > laraRefinementPattern = GetUsedRefinementPattern ( CMesh );
	cout << " AdaptMesh " << endl;cout.flush(); 
	AdaptMesh ( * CMesh,  DivideOrCoarsen, laraRefinementPattern );
	cout << " GetAdaptedMesh finished " << endl;cout.flush(); 
}

/**
 * Auxiliar method to get the solution and the gradients from the PZ solution data structure
 */
void GetSolution ( TPZCompMesh & CMesh,
				  TPZCompEl * CEl,
				  TPZVec < REAL > & Solutions,
				  TPZVec < REAL > & Gradients )
{
	if ( ! CEl ) return;
	Solutions.Resize ( 5 );
	Gradients.Resize ( 15 );
	TPZConnect connect = CEl->Connect ( 0 );
	int seqnum = connect.SequenceNumber();
	TPZBlock &block = CMesh.Block();
	int i = 0;
	for ( i = 0; i < 5; i++ )
	{
		Solutions [ i ] = block.Get ( seqnum, 0, i, 0 );
		Gradients [ 3 * i + 0 ] = block.Get ( seqnum, 0, 3 * i + 5, 0 );
		Gradients [ 3 * i + 1 ] = block.Get ( seqnum, 0, 3 * i + 6, 0 );
		Gradients [ 3 * i + 2 ] = block.Get ( seqnum, 0, 3 * i + 7, 0 );
	}
}
//-----------------------------------------------------------------------------------------------

void ErrorEstimation(TPZCompMesh & CMesh,
					  TPZVec < EAdaptElementAction > & DivideOrCoarsen,
					  double Epslref)
{
	
	int nel = CMesh.NElements();
	DivideOrCoarsen.Resize ( nel );
	DivideOrCoarsen.Fill ( ENone );
	
	// the average solutions are: rho , v (as vector ) and P
	TPZVec < REAL > AverageSolutionFineVec ( 3, 0.0 );
	//TPZVec < REAL > AverageSolutionCoarseVec ( 3, 0.0 );
	
	TPZVec < TPZCompMesh * > gradedMeshVec;
	//Produce graded mesh vector, projecting the solution
	cout <<  "ProduceGradedMeshes ... \n"; cout.flush();
	ProduceGradedMeshes ( CMesh, gradedMeshVec );
	
	//evaluate uhat for each level
	map < int, vector < vector < double > > > levelToElementUhatVec;
	cout <<  "EvaluateUHat( ... \n"; cout.flush(); 
	EvaluateUHat ( gradedMeshVec, levelToElementUhatVec );
	
	//Evaluate average solution for each state variable
	//TPZCompMesh * coarseMesh = gradedMeshVec[ 1 ];
	TPZCompMesh * fineMesh = gradedMeshVec [ 2 ];
	
	cout <<  "EvaluateAverageOfSolution( ... \n"; cout.flush();
	EvaluateAverageOfSolution ( * fineMesh, AverageSolutionFineVec );
	//EvaluateAverageOfSolution ( * coarseMesh, AverageSolutionCoarseVec );
	
	cout << "Average fine " << AverageSolutionFineVec << endl;
	
	TPZManVector < REAL,1000 > fineDetail;
	cout <<  "EvaluateDetail( ... \n"; cout.flush(); 
	EvaluateDetail ( * fineMesh , AverageSolutionFineVec, levelToElementUhatVec, fineDetail );
	
	cout <<  "Adaptando ... \n"; cout.flush();
	int i = 0;
	fineMesh->Reference()->ResetReference();
	fineMesh->LoadReferences(); 
	for ( i = 0; i < fineDetail.NElements(); i++ )
	{
		TPZCompEl * cel = fineMesh->ElementVec()[ i ];
		if (!cel || cel->Type() != EDiscontinuous || cel->NConnects() == 0) continue;
		
		int l = fineMesh->ElementVec()[i]->Reference()->Level();
		double Epsl = Epslref * pow ( 2.0 ,(double)( l - gLMax ));
		
		if ( fineDetail[ i ] > Epsl )
		{
			DivideOrCoarsen[ i ] = EDivide;
			continue;
		}
		if ( fineDetail[ i ] < Epsl )
		{
			DivideOrCoarsen[ i ] = ECoarse;
		}
	}
}


void EvaluateDetail ( TPZCompMesh & CMesh,
					 TPZVec < REAL > & AverageSolutionVec,
					 map < int, vector < vector < double > > > & levelToElementUhatVec,
					 TPZVec < REAL > & Detail )
{
	int iel = 0;
	int nel = CMesh.NElements();
	Detail.Resize ( nel, 0.0 );
#ifdef HUGE_DEBUG
#ifdef LOG4CXX
	
	map < int, vector < vector < double > > >::iterator it;
	for ( it = levelToElementUhatVec.begin(); it != levelToElementUhatVec.end(); it++ )
	{
		std::stringstream sout;
		sout << it->first << endl;
		vector < vector < double > > & solElVec = it->second;
		int nel = solElVec.size();
		for ( iel = 0; iel < nel; iel++ )
		{
			sout << "\tElement [ " << iel << " ] = \t";
			vector < double > & solution = solElVec [ iel ];
			for ( int isol = 0; isol < solution.size(); isol++ )
			{
				sout << "\t" << solution [ isol ] ;
			}
			sout << endl;
		}
		LOGPZ_DEBUG ( logger, sout.str().c_str() );
	}
#endif
#endif
	
	CMesh.Reference()->RestoreReference(&CMesh);
#ifdef LOG4CXX
	std::stringstream sout;
#endif
	for ( iel = 0; iel < nel; iel++ )
	{
		TPZCompEl * cel = CMesh.ElementVec()[ iel ];
		if ( ! cel || cel->Type() != EDiscontinuous || cel->NConnects() == 0 ) continue;
		
		unsigned int celIdx = (unsigned int)cel->Index();
		int celLevel = cel->Reference()->Level();
		if ( levelToElementUhatVec.find(celLevel) == levelToElementUhatVec.end()) continue;
		if ( celIdx < 0 || celIdx >= levelToElementUhatVec[ celLevel].size() ) continue;
		vector < double > uhat = ( levelToElementUhatVec[ celLevel] ) [ celIdx ];
		TPZVec < REAL > sol ( 5, 0.0 );
		TPZVec < REAL > grad ( 15, 0.0 );
		GetSolution ( CMesh, cel, sol, grad );
		double maxDetail = 0.0;
		
		double aux =  AverageSolutionVec[ 0 ] ;
		double dRho = 0.0;
		double dU = 0.0;
		double dV = 0.0;
		double dW = 0.0;
		double dVel = 0.0;
		double dP = 0.0;
		const double tol = 1.e-12;
		if ( fabs (aux) > tol )	dRho = fabs ( sol[ 0 ] - uhat [ 0 ] ) / aux;
		maxDetail = dRho;
		aux  = AverageSolutionVec[ 1 ];
		if ( fabs (aux) > tol ) dU = fabs ( sol[ 1 ] - uhat [ 1 ] ) / aux ;
		aux = AverageSolutionVec[ 2 ];
		if ( fabs (aux) > tol ) dV = fabs ( sol[ 2 ] - uhat [ 2 ] ) / aux ;
		aux = AverageSolutionVec[ 3 ];
		if ( fabs (aux) > tol ) dW = fabs ( sol[ 3 ] - uhat [ 3 ] ) / aux ;
		
		dVel = sqrt ( dU * dU + dV * dV + dW * dW );
		maxDetail = ( dVel > maxDetail ) ? dVel : maxDetail;
		
		aux = AverageSolutionVec[ 4 ];
		if ( fabs (aux) > tol ) dP = fabs ( sol[ 4 ] - uhat [ 4 ] ) / aux ;
		maxDetail = ( dP > maxDetail ) ? dP : maxDetail;
		
		Detail [ celIdx ] = maxDetail;
#ifdef LOG4CXX
		sout << "Detail for element [ " << celIdx << " ] = " << sol[0] << " - " << uhat[0] << " = " << dRho << " | " << dVel << " | " << dP << endl;
#endif
	}
#ifdef LOG4CXX
	LOGPZ_DEBUG(logger,sout.str().c_str());
#endif
}


void ProduceGradedMeshes ( TPZCompMesh & OriginalMesh,
						  TPZVec < TPZCompMesh * > & gradedMeshVec )
{
#ifndef NODEBUG
	{
		if ( !CheckReferences( OriginalMesh ) )
		{
			std::stringstream sout;
			OriginalMesh.Print(sout);
			OriginalMesh.Reference()->Print(sout);
			LOGPZ_ERROR(logger,sout.str().c_str());
			exit(-1);
		}
	}
#endif
	int maxLevel = -1;
	int minLevel = -1;
	int el = 0;
	int nel = OriginalMesh.NElements();
	for ( el = 0; el < nel; el++ )
	{
		TPZCompEl * cel = OriginalMesh.ElementVec()[ el ];
		if ( ! cel ) continue;
		int level = cel->Reference()->Level();
		maxLevel = level > maxLevel ? level : maxLevel;
		minLevel = level < minLevel ? level : minLevel;
	}
	int nlevels = maxLevel - minLevel;
	gradedMeshVec.Resize ( nlevels );
	int im = 0;
#ifdef HUGE_DEBUG
#ifdef LOG4CXX
	{
		TPZGeoMesh * gMesh = OriginalMesh.Reference();
		int i =0;
		std::stringstream sout;
		sout << " Antes do clone: ";
		for (i=0;i<gMesh->NElements();i++)
		{
			TPZGeoEl *gel = gMesh->ElementVec()[i];
			if (!gel) continue;
			sout << gel->Index() << " = " << gel->NumInterfaces() << "\t";
		}
		sout << endl;
		LOGPZ_DEBUG(logger,sout.str().c_str());
	}
#endif
#endif
	gradedMeshVec [0] = OriginalMesh.Clone();
	(gradedMeshVec [0])->LoadReferences();
	
#ifndef NODEBUG
	{
		if ( !CheckReferences(*gradedMeshVec[0]) )
		{
			std::stringstream sout;
			gradedMeshVec[0]->Print(sout);
			gradedMeshVec[0]->Reference()->Print(sout);
			LOGPZ_ERROR(logger,sout.str().c_str());
			exit(-1);
		}
	}
#endif
	
#ifdef HUGE_DEBUG
#ifdef LOG4CXX
	{
		TPZGeoMesh * gMesh = OriginalMesh.Reference();
		int i =0;
		std::stringstream sout;
		sout << " Depois do clone: ";
		for (i=0;i<gMesh->NElements();i++)
		{
			TPZGeoEl *gel = gMesh->ElementVec()[i];
			if (!gel) continue;
			sout << gel->Index() << " = " << gel->NumInterfaces() << "\t";
		}
		sout << endl;
		LOGPZ_DEBUG(logger,sout.str().c_str());
	}
#endif
#endif
	for ( im = 1; im < nlevels; im++ )
	{
		cout << "clonando nivel " << im << "\n"; cout.flush(); 
		gradedMeshVec[ im ] = CoarsenOneLevel ( * (gradedMeshVec [ im - 1 ]) );
#ifndef NODEBUG
		{
			if ( !CheckReferences ( *gradedMeshVec [ im ] ) )
			{
				std::stringstream sout;
				gradedMeshVec[im]->Print(sout);
				gradedMeshVec[im]->Reference()->Print(sout);
				LOGPZ_ERROR(logger,sout.str().c_str());
				exit(-1);
			}
		}
#endif
#ifdef HUGE_DEBUG
#ifdef LOG4CXX
		{
			stringstream sout;
			gradedMeshVec[im]->Print(sout);
			LOGPZ_DEBUG ( logger, sout.str().c_str() );
		}
#endif
#endif
	}
	
	//Verify the projection method
#ifdef HUGE_DEBUG
#ifdef LOG4CXX
	for ( im = 0; im < nlevels; im++ )
	{
		TPZCompMesh * cmesh = gradedMeshVec[ im ];
		stringstream sout;
		sout << "Solution for level =  " << im << " number of elements = " << cmesh->NElements() << endl; 
		PrintMeshSolution(cmesh,sout);
		LOGPZ_DEBUG(logger,sout.str().c_str());
	}
#endif
#endif
}


TPZCompMesh * CoarsenOneLevel ( TPZCompMesh & OriginalMesh )
{
	int maxLevel = 0;
#ifdef HUGE_DEBUG
#ifdef LOG4CXX
	{
		TPZGeoMesh * gMesh = OriginalMesh.Reference();
		int i =0;
		std::stringstream sout;
		sout << " Antes do clone: ";
		for (i=0;i<gMesh->NElements();i++)
		{
			TPZGeoEl *gel = gMesh->ElementVec()[i];
			if (!gel) continue;
			sout << gel->Index() << " = " << gel->NumInterfaces() << "\t";
		}
		sout << endl;
		LOGPZ_DEBUG(logger,sout.str().c_str());
	}
#endif
#endif
	TPZCompMesh * CoarseMesh = OriginalMesh.Clone();
	CoarseMesh->LoadReferences();
#ifdef HUGE_DEBUG
#ifdef LOG4CXX
	{
		TPZGeoMesh * gMesh = OriginalMesh.Reference();
		int i =0;
		std::stringstream sout;
		sout << " Depois do clone: ";
		for (i=0;i<gMesh->NElements();i++)
		{
			TPZGeoEl *gel = gMesh->ElementVec()[i];
			if (!gel) continue;
			sout << gel->Index() << " = " << gel->NumInterfaces() << "\t";
		}
		sout << endl;
		LOGPZ_DEBUG(logger,sout.str().c_str());
	}
#endif
#endif
	
	CoarseMesh->Reference()->RestoreReference(CoarseMesh);
	
	TPZAutoPointer<TPZFunction> fakefunc = new TPZFakeFunction();
#ifdef HUGE_DEBUG
#ifdef LOG4CXX
	{
		stringstream sout;
		sout << "Solution after clone ! " << endl;
		PrintMeshSolution(CoarseMesh, sout);
		LOGPZ_DEBUG ( logger, sout.str().c_str() );
	}
#endif
#endif
	
	int el = 0;
	int nel = CoarseMesh->NElements();
	for ( el = 0; el < nel; el++ )
	{
		TPZCompEl * cel = CoarseMesh->ElementVec()[ el ];
		if ( ! cel || cel->Type() != EDiscontinuous || cel->NConnects() == 0 ) continue;
		int level = cel->Reference()->Level();
		maxLevel = ( level > maxLevel ) ? level : maxLevel;
	}
	
	set < int > alreadyCoarsen;
	
	for ( el = 0; el < nel; el++ )
	{
		TPZCompEl * cel = CoarseMesh->ElementVec()[ el ];
		if ( ! cel || cel->Type() != EDiscontinuous || cel->NConnects() == 0) continue;
		int celIdx = cel->Index();
		if ( alreadyCoarsen.find(celIdx) != alreadyCoarsen.end() ) continue;
		TPZGeoEl * gel = cel->Reference();
		int level = gel->Level();
		if ( level < maxLevel ) continue;
		TPZGeoEl * father = gel->Father();
		if (!father) continue;
		int nsubel = father->NSubElements();
		TPZVec<int> subCElVec ( nsubel );
		int isub = 0;
		int isol = 0;
		
		TPZVec <REAL>  solutionVec ( 5,0. );//blocksize );
		REAL fatherVolume = father->Volume();
		
		for ( isub = 0; isub < nsubel; isub++ )
		{
			TPZGeoEl * subGel = father->SubElement ( isub );
			if (!subGel)
			{
				cout << __PRETTY_FUNCTION__ << " " << __LINE__ << " ";
				cout << " Warning: geometric subelement doesn't exist \n";
				continue;
			}
			TPZCompEl * subCel = subGel->Reference();
			if (!subCel)
			{
				cout << __PRETTY_FUNCTION__ << " " << __LINE__ << " ";
				cout << " Warning: computational subelement doesn't exist \n";
				continue;
			}
			
			subCElVec[ isub ] = subCel->Index();
			//con = subCel->Connect ( 0 );
			//seqnum = con.SequenceNumber();
			int subCelIdx = subCel->Index();
			alreadyCoarsen.insert ( subCelIdx );
			REAL sonVolume = subGel->Volume();
			TPZVec <REAL> solution (5,0.);
			TPZVec <REAL> grad (15,0.);
			GetSolution ( *CoarseMesh, subCel, solution, grad );
			for ( isol = 0; isol < 5/*blocksize*/; isol++ )
			{
				solutionVec[ isol ] += solution[isol] * fabs(sonVolume);//block.Get( seqnum,0,isol,0 ) * sonVolume;
			}
		}
		int coarseIdx = -1;
		CoarseMesh->Coarsen ( subCElVec, coarseIdx, true );
#ifndef NODEBUG
		if ( !CheckReferences(*CoarseMesh) )
		{
			cout << "Teje pego!\n";
			exit (-1);
		}
#endif
		TPZCompEl * coarseEl = CoarseMesh->ElementVec() [coarseIdx];
	    TPZCompElDisc * disc = dynamic_cast<TPZCompElDisc*>(coarseEl);
	    if(disc)
	    {
	    	disc->SetExternalShapeFunction(fakefunc);
		}
	    else
	    {
	    	cout << __PRETTY_FUNCTION__ << " Created a non discontinuous element\n";
	    }
		CoarseMesh->ExpandSolution();
		// 		CoarseMesh->CleanUpUnconnectedNodes();
		// 		CoarseMesh->AdjustBoundaryElements();
		// 		CoarseMesh->ExpandSolution();
		
		TPZConnect coarseCon = coarseEl->Connect(0);
		int coarseSeqNum = coarseCon.SequenceNumber();
		TPZBlock &coarseBlock = CoarseMesh->Block();
		int coarseblocksize = coarseBlock.Size(coarseSeqNum);
		//solution
		for ( isol = 0; isol < 5; isol++)
		{
			coarseBlock.Put( coarseSeqNum, 0, isol, 0, solutionVec [isol] / fabs(fatherVolume) );
		}
		//Gradients
		for (; isol < coarseblocksize; isol++)
		{
			coarseBlock.Put( coarseSeqNum, 0, isol, 0, 0.0 );
		}
	}
	
    CoarseMesh->ExpandSolution();
    CoarseMesh->CleanUpUnconnectedNodes();
    CoarseMesh->AdjustBoundaryElements();
    CoarseMesh->ExpandSolution();
	
	TPZExplFinVolAnal gradAnalysis ( CoarseMesh );
	TPZFMatrix solAndGrad;
	gradAnalysis.ComputeGradientForDetails( CoarseMesh->Solution(),solAndGrad );
	CoarseMesh->LoadSolution(solAndGrad);
#ifdef HUGE_DEBUG	
#ifdef LOG4CXX
	{
		stringstream sout;
		PrintMeshSolution(CoarseMesh, sout);
		LOGPZ_DEBUG ( logger, sout.str().c_str() );
	}
#endif
#endif
	return CoarseMesh;
}


void EvaluateUHat ( TPZVec < TPZCompMesh * > & gradedMeshVec, map < int, vector < vector < double > > > & levelToUhatVec )
{
	int nlevels = gradedMeshVec.NElements();
	int im = 0;		// mesh iterator
	int iel = 0;	// element iterator
	int isub = 0;	// child iterator
	int ist = 0;	// state variable iterator
	int idim = 0;	// dimension iterator
	TPZVec < REAL > solution ( 5, 0.0 );
	TPZVec < REAL > gradient ( 15, 0.0 );
	
	for(im = nlevels-1; im > 0; im-- )                  // Problema com eficiencia 
	{
		TPZCompMesh * coarseMesh = gradedMeshVec[ im ];
		coarseMesh->Reference()->RestoreReference(coarseMesh);
		int nel = coarseMesh->NElements();
		vector < vector < double > > uhatVal ( nel );
		for ( iel = 0; iel < nel; iel++ )
		{			
			TPZCompEl * cel = coarseMesh->ElementVec()[ iel ];
			if ( !cel || cel->Type() != EDiscontinuous || cel->NConnects() == 0 ) continue;
			
			TPZGeoEl * gel = cel->Reference();
			if ( !gel ) continue;
			if ( ! gel->HasSubElement() )
			{
				//copy uhat do nível anterior...
				return;
			}
			coarseMesh->Reference()->RestoreReference(coarseMesh);     
			int nsubel = gel->NSubElements();
			TPZVec < REAL > qsi(3,0.),fatherCenter (3,0.);
			gel->CenterPoint ( gel->NSides() - 1, qsi);
			gel->X(qsi,fatherCenter);
			
			
			TPZCompMesh * fineMesh = gradedMeshVec[ im - 1 ];
			fineMesh->Reference()->RestoreReference(fineMesh);
			
			for ( isub = 0; isub < nsubel; isub++ )
			{
				TPZGeoEl * subGel = gel->SubElement( isub );
				if (!subGel) continue;
				TPZCompEl * subCel = subGel->Reference();
				if (!subCel) continue;
				
				GetSolution ( * fineMesh, subCel, solution, gradient );
				
				int subcelIdx = subCel->Index();
				TPZVec < REAL > sonCenter (3,0.);
				subGel->CenterPoint ( subGel->NSides() - 1, qsi );
				subGel->X(qsi,sonCenter);    
				
				// take the solution and the gradient value from coarse mesh
				vector < double > uhat ( 5, 0.0 );
				for ( ist = 0; ist < 5; ist++)
				{
					uhat [ ist ] = solution [ ist ];
					for ( idim = 0; idim < 3; idim++ )
					{
						uhat [ ist ] += gradient [ ist * 3 + idim ] * ( fatherCenter[ 0 ] - sonCenter [ 0 ] );
					}
				}
#ifdef HUGE_DEBUG
#ifdef LOG4CXX
				{
					stringstream sout;
					sout << "Uhat for [ " << subcelIdx  << " ] = " << uhat [ 0 ] << " , "
					<< uhat [ 1 ] << " , "
					<< uhat [ 2 ] << " , "
					<< uhat [ 3 ] << " , "
					<< uhat [ 4 ] << endl;
					LOGPZ_DEBUG ( logger, sout.str().c_str() );
				}
#endif
#endif
				uhatVal [ subcelIdx ] = uhat;
			}
		}
#ifdef HUGE_DEBUG
#ifdef LOG4CXX
		{
			stringstream sout;
			sout << "UhatVAL for LEVEL " << im << ":\n";
			int iel,isol;
			sout << "Elem\t solution\n";
			for ( iel = 0; iel < uhatVal.size(); iel++ )
			{
				sout << iel << " { ";
				for ( isol = 0; isol < (uhatVal[ iel ]).size() ; isol++ )
				{
					sout << uhatVal[iel][isol] << "\t";
				}
				sout << "}\n";
			}
			LOGPZ_DEBUG ( logger, sout.str().c_str() );
		}
#endif
#endif
		levelToUhatVec [ im ] = uhatVal;
	}
}


//Must return || < u > || = Sum over elements ( fabs ( Center_Solution ) * Element_Volume ) / Domain_Volume;
void  EvaluateAverageOfSolution ( TPZCompMesh & CompMesh,
								 TPZVec < REAL > & AverageVec )
{
	double volumeOfMesh = 0.0;
	double volumeOfElement = 0.0;
	int iel = 0;
	int ist = 0;
	TPZVec < REAL > solution ( 5, 0.0 );
	TPZVec < REAL > gradient ( 15, 0.0 );
	AverageVec.Resize ( 5, 0.0 );
	
	int ncel = CompMesh.NElements();
	
	// proceed for each element: take the solution for each state variable and evaluate its norm
	// multiply the norm of state solution by the volume of the element
	// increase the average vector
	for ( iel = 0; iel < ncel; iel++ )
	{
		TPZCompEl * cel = CompMesh.ElementVec()[ iel ];
		if ( !cel || cel->Type() != EDiscontinuous || cel->NConnects() == 0 ) continue;
		volumeOfElement = fabs (cel->Reference()->Volume());
		GetSolution ( CompMesh, cel, solution, gradient );
		for ( ist = 0; ist < 5; ist++ )
		{
			AverageVec [ ist ] += fabs ( solution [ ist ] ) * volumeOfElement;
		}
		volumeOfMesh += volumeOfElement;
	}
	
	// divide the values by the volume of the entire domain...
	for ( ist = 0; ist < 5; ist++ )
	{
		AverageVec [ ist ] /= volumeOfMesh;
	}
}

// ofstream reallydividecoarsenfile("reallydividecoarsen.txt");
void AdaptMesh ( TPZCompMesh & CMesh,
				TPZVec < EAdaptElementAction > & DivideOrCoarsen,
				TPZAutoPointer < TPZRefPattern > & RefPattern )
{
	//Firt verify if some element marked to coarsen have a brother marked to refine
	int el = 0;
	int nel = CMesh.NElements();
	CMesh.Reference()->ResetReference();
	CMesh.LoadReferences();
	
#ifndef NODEBUG
	if (! CheckReferences(CMesh) )
	{
		cout << "teje pego\n";
		exit (-1);
	}
#endif
	
	
	for (el=0;el<nel;el++)
	{
		TPZCompEl *cel = CMesh.ElementVec()[el];
		if (!cel) continue;
		if ( DivideOrCoarsen[el] != ECoarse ) continue;
		TPZGeoEl *gel = cel->Reference();
		if ( gel->Level() < 3)
		{
			DivideOrCoarsen[el] = ENone;
			continue;
		}
		TPZGeoEl *father = gel->Father();
		if ( !father ) continue;
		int nsubel = father->NSubElements();
		int isub = 0;
		for ( isub=0; isub<nsubel; isub++ )
		{
			TPZGeoEl *subGel = father->SubElement(isub);
			if (!subGel) continue;
			TPZCompEl *subCel = subGel->Reference();
			if (!subCel) continue;
			int subindex = subCel->Index();
			if ( DivideOrCoarsen[subindex] == EDivide )
			{
				DivideOrCoarsen[el] = ENone;
				break;
			}
		}
	}
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Elements to divide or coarsen ";
		int el;
		for(el = 0; el<nel; el++)
		{
			if(DivideOrCoarsen[el] != ENone) sout << el << ":" << DivideOrCoarsen[el] << " ";
		}
		sout << std::endl;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	
	TPZAutoPointer<TPZFunction> fakefunc = new TPZFakeFunction();	
	//Let's divide the elements
	for (el=0;el<nel;el++)
	{
		if ( DivideOrCoarsen[el] != EDivide ) continue;
		TPZCompEl *cel = CMesh.ElementVec()[el];
		if (!cel || cel->Type() != EDiscontinuous || cel->NConnects() == 0) continue;
		int index = cel->Index();
		if (cel->Reference()->Level() > gLMax ) continue;
		TPZManVector<int,8> subIndex;
#ifdef LOG4CXX_KEEP
		{
			std::stringstream sout;
			sout << "Element to divide " << el << std::endl;
			cel->Print(sout);
			TPZGeoEl *gel = cel->Reference();
			gel->Print(sout);
			int side;
			for(side=10; side<14; side++)
			{
				TPZGeoElSide gelside(gel,side);
				TPZGeoElSide neighbour = gelside.Neighbour();
				while(neighbour != gelside)
				{
					if(neighbour.Element()->Reference())
					{
						sout << "Computational element along side " << side << std::endl;
						neighbour.Element()->Reference()->Print(sout);
					}
					neighbour = neighbour.Neighbour();
				}
			}
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
		cel->Divide ( index, subIndex, true );
		int isub = 0;
		for ( isub=0; isub< subIndex.NElements(); isub++ )
		{
			int subElIndex = subIndex [ isub ];
			TPZCompEl * subCel = CMesh.ElementVec()[ subElIndex ];
			if ( ! subCel ) continue;
			TPZGeoEl * subGel = subCel->Reference();
			subGel->SetRefPattern ( RefPattern );
			
			TPZCompElDisc * disc = dynamic_cast<TPZCompElDisc*>(subCel);
			if(disc)
			{
				disc->SetExternalShapeFunction(fakefunc);
			}
			else
			{
				cout << __PRETTY_FUNCTION__ << " Created a non discontinuous element\n";
			}
			CMesh.ExpandSolution();
		}
	}
	CMesh.ExpandSolution();
	CMesh.CleanUpUnconnectedNodes();
	CMesh.AdjustBoundaryElements();
	CMesh.ExpandSolution();
#ifdef LOG4CXX_KEEP
	{
		std::stringstream sout;
		CMesh.ExpandSolution();
		CMesh.Print(sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	int isol;
	//Let's coarsen the elements
	for (el=0;el<nel;el++)
	{
		if ( DivideOrCoarsen[el] != ECoarse ) continue;
		TPZCompEl *cel = CMesh.ElementVec()[el];
		if (!cel  || cel->Type() != EDiscontinuous || cel->NConnects() == 0 ) continue;
		if(cel->Reference()->Dimension() != 3) continue;
		TPZGeoEl *gel = cel->Reference();
		TPZGeoEl *father = gel->Father();
		if ( !father ) continue; //there are not brothers to coasen
		int nsubel = father->NSubElements();
		int isub = 0;
		TPZManVector<int> subCElVec ( nsubel,-1 );
		TPZVec<REAL> solutionVec (5,0.);
		for ( isub=0; isub<nsubel; isub++ )
		{
			TPZGeoEl *subGel = father->SubElement(isub);
			if (!subGel) continue;
			TPZCompEl *subCel = subGel->Reference();
			if (!subCel) break;
#ifdef DEBUG
		    if ( subCel->Reference() != subGel )
			{
				std::stringstream sout;
				sout << "Error: Original son " << isub << endl;
				subGel->Print(sout);
				sout << "Reference of the son\n";
				subCel->Print(sout);
				sout << "Reference of the computational element\n";
				subCel->Reference()->Print(sout);
				
				sout << "Result of the consistency test = " << CheckReferences(CMesh) << endl;
				cout << sout.str().c_str() << endl;
				continue;
			}
#endif
			int subindex = subCel->Index();
			if ( DivideOrCoarsen[subindex] == EDivide )
			{
				//        cout << "\ncedric = ";   
				//        int opa;
				//        cin >> opa;
				break;
			}
			subCElVec [isub] = subindex;
#ifdef DEBUG
			if ( (CMesh.ElementVec()[subindex])->Reference()->Father() != father )
			{
				std::stringstream sout;
				sout << "Error: Original son " << isub << endl;
				subGel->Print(sout);
				sout << "Reference of the son\n";
				CMesh.ElementVec()[subindex]->Print(sout);
				sout << "Reference of the computational element\n";
				CMesh.ElementVec()[subindex]->Reference()->Print(sout);
				
				sout << "Result of the consistency test = " << CheckReferences(CMesh) << endl;
				cout << sout.str().c_str() << endl;
				continue;
			}
#endif
			REAL sonVolume = subGel->Volume();
			TPZVec <REAL> solution (5,0.);
			TPZVec <REAL> grad (15,0.);
			GetSolution ( CMesh, subCel, solution, grad );
			for ( isol = 0; isol < 5; isol++ )
			{
				solutionVec[ isol ] += solution[isol] * fabs(sonVolume);			
			}
		}
		
		if ( isub == nsubel ) //all brothers are marked to coarse
		{
			for ( isub=0; isub<nsubel; isub++ )
			{
				DivideOrCoarsen [ subCElVec [ isub ] ] = ENone;
			}
			int newindex = -1;
			CMesh.Coarsen ( subCElVec, newindex, true );
			if (newindex < 0 )
			{
				cout << "Error trying to coarse an element - Index = -1 was returned!\n";
				stringstream sout;
				sout << "Computational list of elements " << subCElVec << endl;
				sout << " cel idx | gel idx | father idx\n";
				for (int ccc= 0; ccc<subCElVec.NElements();ccc++)
				{
					sout << CMesh.ElementVec()[ccc]->Index() << "\t" 
					<< CMesh.ElementVec()[ccc]->Reference()->Index() << "\t"
					<< CMesh.ElementVec()[ccc]->Reference()->Father()-> Index() << endl;
				}
				sout << "Father\n";
				father->Print(sout);
				sout << "filhos:\n";
				for ( isub=0; isub<nsubel; isub++ )
				{
					CMesh.ElementVec()[subCElVec [ isub ] ]->Reference()->Print(sout);
				} 
				cout << sout.str().c_str() <<endl;
				continue;
			}
			
			TPZCompEl * coarseEl = CMesh.ElementVec() [newindex];
			TPZCompElDisc * disc = dynamic_cast<TPZCompElDisc*>(coarseEl);
			if(disc)
			{
				disc->SetExternalShapeFunction(fakefunc);
			}
			else
			{
				cout << __PRETTY_FUNCTION__ << " Created a non discontinuous element\n";
			}
			// 			CMesh.ExpandSolution();
			// 			CMesh.CleanUpUnconnectedNodes();
			// 			CMesh.AdjustBoundaryElements();
			CMesh.ExpandSolution();
			
			TPZConnect coarseCon = coarseEl->Connect(0);
			int coarseSeqNum = coarseCon.SequenceNumber();
			TPZBlock &coarseBlock = CMesh.Block();
			int coarseblocksize = coarseBlock.Size(coarseSeqNum);
			//solution
			REAL fatherVolume = fabs(father->Volume());
			for ( isol = 0; isol < 5; isol++)
			{
				coarseBlock.Put( coarseSeqNum, 0, isol, 0, solutionVec [isol] / fatherVolume );
			}
			//Gradients
			
			for (; isol < coarseblocksize; isol++)
			{
				coarseBlock.Put( coarseSeqNum, 0, isol, 0, 0.0 );
			}
		}
	}
	CMesh.ExpandSolution();
	CMesh.CleanUpUnconnectedNodes();
	CMesh.AdjustBoundaryElements();
	CMesh.ExpandSolution();
	//By now, the mesh is adapted. We need to check the difference of level of refinement between neighbors
	// 	CheckRefinementLevel ( CMesh, RefPattern );
}


// Select the elements with low level of refinement (when comparing with the neighbors)
void CheckRefinementLevel ( TPZCompMesh & CMesh,
						   TPZAutoPointer < TPZRefPattern > & RefPattern )
{
	list < int > elementsToDivide;
	SelectElementsByLevel ( CMesh, elementsToDivide );
	int count = 0;
	while ( elementsToDivide.size() && count < 2 )
	{
		RefineElements ( CMesh, elementsToDivide, RefPattern );
		SelectElementsByLevel ( CMesh, elementsToDivide );
		count ++;
	}
	if (count == 10) {
		cout << " Check refinement level skipped after reach the maximum number of steps\n";
	}
}


void SelectElementsByLevel ( TPZCompMesh & CMesh,
							list < int > & SelectedElementsIdx )
{
#ifdef LOG4CXX
	{
		std::stringstream sout;
		CMesh.Print(sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	// Using the set class we are sure that we don't have any duplicates...
	set < TPZCompEl *> elementsToRefine;
	
	// counter for the elements
	int el = 0;
	
	//number of elements in computational mesh
	int nel = CMesh.NElements();
	
	// loop over the elements of the computational mesh
	// here we are looking for the interface elements.
	// Interface elements are elements of the dimension equal to dimension of elements - 1 ( i.e. faces for this case)
	// the interesting point is that interfaces have pointer to the elements on left and right side
	// Here we want to compare the difference of levels between these two elements.
	// If the difference is greater than 1 we put the element with lower level in the set of elements to divide
	for ( el = 0; el < nel; el++ )
	{
		// Take an element as a discontinuous element
		TPZInterfaceElement *cel = dynamic_cast < TPZInterfaceElement* > ( CMesh.ElementVec()[el] );
		
		// verify if the element does not exist or the cast fail proceed to the next element
		if ( ! cel ) continue;
		
		// verify if the element is an interface. Elsewhere proceed to the next element
		if ( cel->Type() != EInterface ) continue;
		
		// take the left and righ volumes associated to the interface
		TPZCompEl * cLeft = cel->LeftElement();
		TPZCompEl * cRight = cel->RightElement();
		
		if ( ! cLeft || ! cRight ) continue;
		
		// take the geometric elements associated to the volumes.
		// Note that the level of refinement is associated to the geometry
		TPZGeoEl * gLeft = cLeft->Reference();
		TPZGeoEl * gRight = cRight->Reference();
		
		// take the level of refinement of each volume
		int levelLeft = gLeft->Level();
		int levelRight = gRight->Level();
		
		// check the difference of levels
		if ( abs ( levelLeft - levelRight ) > 2 )
		{
			TPZVec < int > subIndex;
			if ( levelLeft < levelRight )
			{
				elementsToRefine.insert ( cLeft );
			}
			else
			{
				elementsToRefine.insert ( cRight );
			}
		}
	}
	
	// Move the list of elements from a set to a list
	// Note that a set is a sorted unique associative container instead the list in nor sorted neither unique...
	// In the selection we use the set to guarantee that an unique instance of an element is being selected.
	// The list is easier to insert an remove data...
	set < TPZCompEl * >::iterator it;
	for ( it = elementsToRefine.begin(); it != elementsToRefine.end(); it++ )
	{
		TPZCompEl * cel = *it;
		if (!cel  || cel->Type() != EDiscontinuous || cel->NConnects() == 0) continue;
		int index = cel->Index();
		SelectedElementsIdx.push_back (index);
	}
}


// Refine the elements on the RefineList
void RefineElements ( TPZCompMesh & CMesh,
					 list < int > & RefineList,
					 TPZAutoPointer < TPZRefPattern > & RefPattern )
{
	int siz = RefineList.size();
	TPZAutoPointer<TPZFunction> fakefunc = new TPZFakeFunction();
	while ( RefineList.size() )
	{
		// take a pointer to element
		int index = RefineList.front();
		RefineList.pop_front();
		if (index < 0 || index > CMesh.NElements() )
		{
			cout << "trying to divide element index " << index << " which is out of range ( 0 - " << CMesh.NElements() << endl;
			continue;
		}
		TPZCompEl * cel = CMesh.ElementVec()[index];
		//remove this pointer from the list
		siz = RefineList.size();
		//verify if the pointer exists
		if ( ! cel  || cel->Type() != EDiscontinuous || cel->NConnects() == 0 ) continue;
		//take index of corresponding element
		//create an vector to receive the indexes of the children
		TPZManVector<int> subIndex(8,0.);
		if (cel->Reference()->Level() > gLMax) continue;
		// refinement function...
		cel->Divide ( index, subIndex, true );
		int isub = 0;
		for ( isub=0; isub< subIndex.NElements(); isub++ )
		{
			int subElIndex = subIndex [ isub ];
			TPZCompEl * subCel = CMesh.ElementVec()[ subElIndex ];
			if ( ! subCel ) continue;
			TPZGeoEl * subGel = subCel->Reference();
			if ( !subGel) continue;
			subGel->SetRefPattern ( RefPattern );
			
			TPZCompElDisc * disc = dynamic_cast<TPZCompElDisc*>(subCel);
			if(disc)
			{
				disc->SetExternalShapeFunction(fakefunc);
			}
			else
			{
				cout << __PRETTY_FUNCTION__ << " Created a non discontinuous element\n";
			}
			CMesh.ExpandSolution();
			CMesh.CleanUpUnconnectedNodes();
			CMesh.AdjustBoundaryElements();
			CMesh.ExpandSolution();
		}
	}
	CMesh.ExpandSolution();
	CMesh.CleanUpUnconnectedNodes();
	CMesh.AdjustBoundaryElements();
	CMesh.ExpandSolution();	
}

void LoadDummySolution(TPZCompMesh *cmesh)
{
	int nel = cmesh->NElements();
	int iel;
	for(iel = 0; iel<nel; iel++)
	{
		TPZCompEl *cel = cmesh->ElementVec()[iel];
		if(!cel || !cel->Type()==EDiscontinuous || cel->NConnects() == 0 ) continue;
		TPZCompElDisc *disc = dynamic_cast<TPZCompElDisc *>(cel);
		if(!disc) continue;
		if(!disc->NConnects()) continue;
		TPZGeoEl *gel = disc->Reference();
		TPZManVector<REAL> center(3,0.);
		gel->CenterPoint(gel->NSides()-1,center);
		TPZVec<REAL> xptr (3,0.);
		gel->X( center,xptr );
		TPZManVector<REAL> val(5,0.);
		DummyFunction2(xptr,val);
		TPZConnect &c = disc->Connect(0);
		int seqnum = c.SequenceNumber();
		TPZBlock &bl = cmesh->Block();
		int i;
		for(i=0; i<5; i++)
		{
			bl(seqnum,0,i,0) = val[i];
		}
	}
#ifdef LOG4CXX
	{
		stringstream sout;
		sout << "Solution after LoadDummySolution: "; 
		PrintMeshSolution(cmesh,sout);
		LOGPZ_DEBUG(logger,sout.str().c_str());
	}
#endif
}

void PrintMeshSolution ( TPZCompMesh * cmesh, ostream & sout)
{
	sout << __PRETTY_FUNCTION__ << endl;
	int iel= 0;
	int nel = cmesh->NElements();
	int count = 0;
	for (iel=0;iel<nel;iel++)
	{
		TPZCompEl * cel = cmesh->ElementVec()[iel];
		if ( !cel || cel->Type() != EDiscontinuous || cel->NConnects() == 0 ) continue;
		count ++;
		TPZGeoEl * gel = cel->Reference();
		TPZVec <REAL> solution (5,0.0);
		TPZVec <REAL> grad (15,0.);
		TPZVec <REAL> center (3,0.);
		TPZVec<REAL> xptr (3,0.);
		gel->CenterPoint ( gel->NSides()-1, center );
		gel->X ( center, xptr );
		GetSolution ( *cmesh, cel, solution, grad );
		if ( fabs (solution [0] -xptr[0]) > 1.e-12)
		{
			sout << "ERROR___";
		}
		else
		{
			sout << "OK___";
		}
		
		sout << "elc/elg [ " << cel->Index() << " , " 
		<< gel->Index() <<  " ] = solution { " 
		<< solution << " } | center { " << xptr << " } " << endl;
	}
	sout << "number of volumes in this mesh = " << count << endl;	
}


bool CheckReferences ( TPZCompMesh & CMesh )
{
	int iel;
	for (iel = 0; iel < CMesh.NElements(); iel++)
	{
		TPZCompEl * cel = CMesh.ElementVec()[iel];
		if ( !cel ) continue;
		bool res = CheckElementReferences(cel);
		if( !res) 
		{
			return false;
		}
	}
	return true;
}


bool CheckElementReferences(TPZCompEl * CEl )
{
	TPZGeoEl * gel = CEl->Reference();
	TPZCompEl *celByRef = gel->Reference();
	if (CEl != celByRef)
	{
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << " " << __LINE__
		<< " inconsistency detected for element:\n";
		CEl->Print(sout);
		sout << "Geometric reference\n";
		gel->Print(sout);
		sout << " Computational reference of the geometric element:\n";
		celByRef->Print(sout);
		
#ifdef LOG4CXX
		LOGPZ_ERROR (logger,sout.str().c_str());
#endif
		cout << sout.str().c_str();
		return false;
	}
	return true;
}