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
	TPZRefPattern *tetraPattern = new TPZRefPattern ( gDummyMesh, tetraFile );

	// Include the pattern info into the mesh
	tetraPattern->InsertPermuted();

	// Remove this object
	delete tetraPattern;


	// --- Define the refinement patterns for triangles -------

	// Read the refinement pattern
	string triangleFile ( "LaraRefineTriangle.txt" );

	// Create and initialize an object to store the pattern
	TPZRefPattern * trianglePattern = new TPZRefPattern ( gDummyMesh, triangleFile );

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

	// Define the element type (tetrahedron)
	MElementType eltype = ETetraedro;

	// Take the map of every refinement pattern for a tetrahedron defined into the mesh
	// the map below provides a relation between an index and a pointer to an object of type refinement pattern
	map < int, TPZAutoPointer < TPZRefPattern > > mapOfTetraRefPattern = gDummyMesh->RefPatternList ( eltype );
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
	TPZCompMesh::SetAllCreateFunctionsDiscontinuous();

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
		ErrorEstimation ( * cDummyMesh, DivideOrCoarsen );

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
void GetAdaptedMesh ( TPZCompMesh * CMesh )
{
	TPZVec < EAdaptElementAction > DivideOrCoarsen;
	// Call the error evaluation and fill the decision vector for each element ( divide - coarse - none )
	ErrorEstimation ( * CMesh, DivideOrCoarsen );
	TPZAutoPointer < TPZRefPattern > laraRefinementPattern = GetUsedRefinementPattern ( CMesh );
	AdaptMesh ( * CMesh,  DivideOrCoarsen, laraRefinementPattern );
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
		Solutions [ i ] = block.Get ( seqnum + i, 0, 0, 0 );
		Gradients [ 3 * i + 5 ] = block.Get ( seqnum + 3 * i + 5, 0, 0, 0 );
		Gradients [ 3 * i + 6 ] = block.Get ( seqnum + 3 * i + 6, 0, 0, 0 );
		Gradients [ 3 * i + 7 ] = block.Get ( seqnum + 3 * i + 7, 0, 0, 0 );
	}
}
//-----------------------------------------------------------------------------------------------

void ErrorEstimation ( TPZCompMesh & CMesh,
					   TPZVec < EAdaptElementAction > & DivideOrCoarsen )
{

	int nel = CMesh.NElements();
	DivideOrCoarsen.Resize ( nel );
	DivideOrCoarsen.Fill ( ENone );

	// the average solutions are: rho , v (as vector ) and P
	TPZVec < REAL > AverageSolutionFineVec ( 3, 0.0 );
	TPZVec < REAL > AverageSolutionCoarseVec ( 3, 0.0 );

	TPZVec < TPZCompMesh * > gradedMeshVec;
	//Produce graded mesh vector, projecting the solution
	ProduceGradedMeshes ( CMesh, gradedMeshVec );

	//evaluate uhat for each level
	map < int, vector < vector < double > > > levelToElementUhatVec;
	EvaluateUHat ( gradedMeshVec, levelToElementUhatVec );

	//Evaluate average solution for each state variable
	TPZCompMesh * coarseMesh = gradedMeshVec[ 0 ];
	TPZCompMesh * fineMesh = gradedMeshVec [ 1 ];

	EvaluateAverageOfSolution ( * fineMesh, AverageSolutionFineVec );
	EvaluateAverageOfSolution ( * coarseMesh, AverageSolutionCoarseVec );

	TPZVec < REAL > fineDetail;
	TPZVec < REAL > coarseDetail;
	EvaluateDetail ( * fineMesh , AverageSolutionFineVec, levelToElementUhatVec, fineDetail );
	EvaluateDetail ( * coarseMesh, AverageSolutionCoarseVec, levelToElementUhatVec, coarseDetail );

	int i = 0;

	for ( i = 0; i < fineDetail.NElements(); i++ )
	{
		fineMesh->LoadReferences();

		int l = fineMesh->ElementVec()[i]->Reference()->Level();
		int L = l + 1;
		double NP = 1.0;
		double Epsl = pow ( 2.0 , 3.0 * (double)( l - L ) / NP );
		if ( fineDetail[ i ] > Epsl )
		{
			DivideOrCoarsen[ i ] = EDivide;
			continue;
		}
		if ( fineDetail[ i ] <= Epsl )
		{
			coarseMesh->LoadReferences();

			TPZCompEl * fineCel = fineMesh->ElementVec() [ i ];
			if ( ! fineCel )
			{
				continue;
			}
			TPZGeoEl * sonGel = fineCel->Reference();
			if ( ! sonGel )
			{
				continue;
			}
			TPZGeoEl * fatherGel = sonGel->Father();
			if ( ! fatherGel )
			{
				continue;
			}
			TPZCompEl *fatherCel = fatherGel->Reference();
			if ( ! fatherCel )
			{
				continue;
			}
			int coarseIndex = fatherCel->Index();
			if ( coarseDetail [ coarseIndex ] < Epsl )  DivideOrCoarsen[ i ] = ECoarse;
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

	CMesh.LoadReferences();

	for ( iel = 0; iel < nel; iel++ )
	{
		TPZCompEl * cel = CMesh.ElementVec()[ iel ];
		if ( ! cel || cel->Type() != EDiscontinuous ) continue;

		int celIdx = cel->Index();
		int celLevel = cel->Reference()->Level();

		vector < double > uhat = ( levelToElementUhatVec[ celLevel] ) [ celIdx ];
		TPZVec < REAL > sol ( 5, 0.0 );
		TPZVec < REAL > grad ( 15, 0.0 );
		GetSolution ( CMesh, cel, sol, grad );
		double maxDetail = 0.0;

		double dRho = fabs ( sol[ 0 ] - uhat [ 0 ] ) / AverageSolutionVec[ 0 ] ;
		maxDetail = dRho;
		double dU = fabs ( sol[ 1 ] - uhat [ 1 ] ) / AverageSolutionVec[ 1 ] ;
		double dV = fabs ( sol[ 2 ] - uhat [ 2 ] ) / AverageSolutionVec[ 2 ] ;
		double dW = fabs ( sol[ 3 ] - uhat [ 3 ] ) / AverageSolutionVec[ 3 ] ;
		double dVel = sqrt ( dU * dU + dV * dV + dW * dW );
		maxDetail = ( dVel > maxDetail ) ? dVel : maxDetail;

		double dP = fabs ( sol[ 4 ] - uhat [ 4 ] ) / AverageSolutionVec[ 4 ] ;
		maxDetail = ( dP > maxDetail ) ? dP : maxDetail;

		Detail [ celIdx ] = maxDetail;
	}
}


void ProduceGradedMeshes ( TPZCompMesh & OriginalMesh,
						   TPZVec < TPZCompMesh * > & gradedMeshVec )
{
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
	gradedMeshVec [0] = OriginalMesh.Clone();
	for ( im = 1; im < nlevels; im++ )
	{
		gradedMeshVec[ im ] = CoarsenOneLevel ( * gradedMeshVec [ im - 1 ] );
	}
}


TPZCompMesh * CoarsenOneLevel ( TPZCompMesh & OriginalMesh )
{
	int maxLevel = 0;
	TPZCompMesh * CoarseMesh = OriginalMesh.Clone();
	CoarseMesh->LoadReferences();
	TPZAutoPointer<TPZFunction> fakefunc = new TPZFakeFunction();

	int el = 0;
	int nel = CoarseMesh->NElements();
	for ( el = 0; el < nel; el++ )
	{
		TPZCompEl * cel = CoarseMesh->ElementVec()[ el ];
		if ( ! cel ) continue;
		int level = cel->Reference()->Level();
		maxLevel = ( level > maxLevel ) ? level : maxLevel;
	}
	for ( el = 0; el < nel; el++ )
	{
		TPZCompEl * cel = CoarseMesh->ElementVec()[ el ];
		if ( ! cel || cel->Type() != EDiscontinuous || cel->NConnects() == 0) continue;
		TPZGeoEl * gel = cel->Reference();
		int level = gel->Level();
		if ( level < maxLevel ) continue;
		TPZGeoEl * father = gel->Father();
		int nsubel = father->NSubElements();
		TPZVec<int> subCElVec ( nsubel );
		int isub = 0;
		int isol = 0;

		TPZConnect con = cel->Connect ( 0 );
		int seqnum = con.SequenceNumber();
		TPZBlock &block = CoarseMesh->Block();
		int blocksize = block.Size( seqnum);

		TPZVec <REAL>  solutionVec ( blocksize );

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
			con = subCel->Connect ( 0 );
			seqnum = con.SequenceNumber();
			REAL sonVolume = subGel->Volume();
			for ( isol = 0; isol < blocksize; isol++ )
			{
				solutionVec[ isol ] += block.Get( seqnum,0,isol,0 ) * sonVolume;
			}
		}
		int coarseIdx = -1;
		CoarseMesh->Coarsen ( subCElVec, coarseIdx, true );
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
		TPZConnect coarseCon = coarseEl->Connect(0);
		int coarseSeqNum = coarseCon.SequenceNumber();
		TPZBlock &coarseBlock = CoarseMesh->Block();
		int coarseblocksize = coarseBlock.Size(coarseSeqNum);
		if(coarseblocksize != blocksize)
		{
			cout << __PRETTY_FUNCTION__ << " coarseblocksize " << coarseblocksize << " blocksize " << blocksize << std::endl;
		}

		for ( isol = 0; isol < coarseblocksize; isol++)
		{
			coarseBlock.Put( coarseSeqNum, 0, isol, 0, solutionVec [isol] / fatherVolume );
		}
	}
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

	for ( im = nlevels - 1; im > 0; im-- )
	{
		TPZCompMesh * coarseMesh = gradedMeshVec[ im ];
		coarseMesh->LoadReferences();
		int nel = coarseMesh->NElements();
		vector < vector < double > > uhatVal ( nel );
		for ( iel = 0; iel < nel; iel++ )
		{
			coarseMesh->LoadReferences();
			TPZCompEl * cel = coarseMesh->ElementVec()[ iel ];
			if ( !cel || cel->Type() != EDiscontinuous ) continue;
			TPZGeoEl * gel = cel->Reference();
			if ( !gel ) continue;
			if ( ! gel->HasSubElement() )
			{
				//copy uhat do nível anterior...
				return;
			}
			int nsubel = gel->NSubElements();
			TPZVec < REAL > fatherCenter (3,0.);
			gel->CenterPoint ( gel->NSides() - 1, fatherCenter );


			TPZCompMesh * fineMesh = gradedMeshVec[ im + 1 ];
			fineMesh->LoadReferences();

			for ( isub = 0; isub < nsubel; isub++ )
			{
				TPZGeoEl * subGel = gel->SubElement( isub );
				if (!subGel) continue;
				TPZCompEl * subCel = subGel->Reference();
				if (!subCel) continue;

				GetSolution ( * fineMesh, subCel, solution, gradient );

				int subcelIdx = subCel->Index();
				TPZVec < REAL > sonCenter (3,0.);
				subGel->CenterPoint ( subGel->NSides() - 1, sonCenter );

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
				uhatVal [ subcelIdx ] = uhat;
			}
		}
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
		if ( !cel || cel->Type() != EDiscontinuous ) continue;
		volumeOfElement = cel->Reference()->Volume();
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


void AdaptMesh ( TPZCompMesh & CMesh,
				 TPZVec < EAdaptElementAction > & DivideOrCoarsen,
				 TPZAutoPointer < TPZRefPattern > & RefPattern )
{
	//Firt verify if some element marked to coarsen have a brother marked to refine
	int el = 0;
	int nel = CMesh.NElements();
	for (el=0;el<nel;el++)
	{
		TPZCompEl *cel = CMesh.ElementVec()[el];
		if (!cel) continue;
		if ( DivideOrCoarsen[el] != ECoarse ) continue;
		TPZGeoEl *gel = cel->Reference();
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

	//Let's divide the elements
	for (el=0;el<nel;el++)
	{
		if ( DivideOrCoarsen[el] != EDivide ) continue;
		TPZCompEl *cel = CMesh.ElementVec()[el];
		if (!cel) continue;
		int index = cel->Index();
		TPZVec<int> subIndex;
		cel->Divide ( index, subIndex );
		int isub = 0;
		for ( isub=0; isub< subIndex.NElements(); isub++ )
		{
			int subElIndex = subIndex [ isub ];
			TPZCompEl * subCel = CMesh.ElementVec()[ subElIndex ];
			if ( ! subCel ) continue;
			TPZGeoEl * subGel = subCel->Reference();
			subGel->SetRefPattern ( RefPattern );
		}
	}

	//Let's coarsen the elements
	for (el=0;el<nel;el++)
	{
		if ( DivideOrCoarsen[el] != ECoarse ) continue;
		TPZCompEl *cel = CMesh.ElementVec()[el];
		if (!cel) continue;
		TPZGeoEl *gel = cel->Reference();
		TPZGeoEl *father = gel->Father();
		if ( !father ) continue; //there are not brothers to coasen
		int nsubel = father->NSubElements();
		int isub = 0;
		TPZVec<int> subCElVec ( nsubel );
		for ( isub=0; isub<nsubel; isub++ )
		{
			TPZGeoEl *subGel = father->SubElement(isub);
			if (!subGel) continue;
			TPZCompEl *subCel = subGel->Reference();
			if (!subCel) continue;
			int subindex = subCel->Index();
			subCElVec [isub] = subindex;
			if ( DivideOrCoarsen[subindex] != ECoarse )
			{
				break;
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
		}
	}

	//By now, the mesh is adapted. We need to check the difference of level of refinement between neighbors
	CheckRefinementLevel ( CMesh, RefPattern );
}


// Select the elements with low level of refinement (when comparing with the neighbors)
void CheckRefinementLevel ( TPZCompMesh & CMesh,
							TPZAutoPointer < TPZRefPattern > & RefPattern )
{
	list < TPZCompEl *> elementsToDivide;
	SelectElementsByLevel ( CMesh, elementsToDivide );
	int count = 0;
	while ( elementsToDivide.size() && count < 1000 )
	{
		RefineElements ( CMesh, elementsToDivide, RefPattern );
		SelectElementsByLevel ( CMesh, elementsToDivide );
		count ++;
	}
}


void SelectElementsByLevel ( TPZCompMesh & CMesh,
						     list < TPZCompEl * > & SelectedElements )
{
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
		if ( abs ( levelLeft - levelRight ) > 1 )
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
		SelectedElements.push_back ( *it );
	}
}


// Refine the elements on the RefineList
void RefineElements ( TPZCompMesh & CMesh,
					  list < TPZCompEl * > RefineList,
					  TPZAutoPointer < TPZRefPattern > & RefPattern )
{
	while ( RefineList.size() )
	{
		// take a pointer to element
		TPZCompEl * cel = RefineList.front();
		//remove this pointer from the list
		RefineList.pop_front();
		//verify if the pointer exists
		if ( ! cel ) continue;
		//take index of corresponding element
		int index = cel->Index();
		//create an vector to receive the indexes of the children
		TPZVec<int> subIndex;
		// refinement function...
		cel->Divide ( index, subIndex );
		int isub = 0;
		for ( isub=0; isub< subIndex.NElements(); isub++ )
		{
			int subElIndex = subIndex [ isub ];
			TPZCompEl * subCel = CMesh.ElementVec()[ subElIndex ];
			if ( ! subCel ) continue;
			TPZGeoEl * subGel = subCel->Reference();
			subGel->SetRefPattern ( RefPattern );
		}
	}
}

void LoadDummySolution(TPZCompMesh *cmesh)
{
	int nel = cmesh->NElements();
	int iel;
	for(iel = 0; iel<nel; iel++)
	{
		TPZCompEl *cel = cmesh->ElementVec()[iel];
		if(!cel || !cel->Type()==EDiscontinuous) continue;
		TPZCompElDisc *disc = dynamic_cast<TPZCompElDisc *>(cel);
		if(!disc) continue;
		if(!disc->NConnects()) continue;
		TPZGeoEl *gel = disc->Reference();
		TPZManVector<REAL> center(3,0.);
		gel->CenterPoint(gel->NSides()-1,center);
		TPZManVector<REAL> val(5,0.);
		DummyFunction2(center,val);
		TPZConnect &c = disc->Connect(0);
		int seqnum = c.SequenceNumber();
		TPZBlock &bl = cmesh->Block();
		int i;
		for(i=0; i<5; i++)
		{
			bl(seqnum,0,i,0) = val[i];
		}
	}
}
