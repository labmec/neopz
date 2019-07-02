/*
 *  visual.cpp
 *  SubstructEigen
 *
 *  Created by Philippe Devloo on 03/12/08.
 *  Copyright 2008 UNICAMP. All rights reserved.
 *
 */

#include "visual.h"
#include "vtkActor.h"

/*
 void TMDIChild::RubberSheetVTK()
{
	
	TPZCompMesh * cmesh = GetRefinedCMesh();
	m_vtk->GetRenderer()->RemoveAllViewProps();
	
	int resolution = 8;
	const int nel = cmesh->NElements();
	TPZVec<REAL> qsi(3,0.), X(3,0.);
	TPZVec<REAL> sol;
	TPZFMatrix dsol;
	TPZFMatrix axes(3,3,0.);
	double minP, maxP;
	minP = 1e20;
	maxP = -1e20;
	
	
	// Create a float array which represents the points.
	vtkFloatArray * pcoords = vtkFloatArray::New();
	
	// Note that by default, an array has 1 component.
	// We have to change it to 3 for points
	pcoords->SetNumberOfComponents(3);
	// Assign each tuple. There are 5 specialized versions of SetTuple:
	// SetTuple1 SetTuple2 SetTuple3 SetTuple4 SetTuple9
	// These take 1, 2, 3, 4 and 9 components respectively.
	float pts[3];
	
	int gValue = 0;
	vtkCellArray* cells = vtkCellArray::New();
	vtkCellArray* lines = vtkCellArray::New();
	vtkDoubleArray* pressure = vtkDoubleArray::New();
	
	pressure->SetName("pressure");
	vtkPoints* points = vtkPoints::New();
	
	
	std::vector< pair<int, double> > pvalues;
	std::vector< pair< int, pair<double, double> > > dvalues;
	
	for(int iel = 0; iel < nel; iel++)
	{
		TPZCompEl * cel = cmesh->ElementVec()[iel];
		if(!cel) continue;
		TPZInterpolationSpace * sp = dynamic_cast<TPZInterpolationSpace*>(cel);
		if(!sp) continue;
		
		TPZIntPoints & intrule = sp->GetIntegrationRule();
		TPZManVector<int,3> order(sp->Dimension(),resolution);
		intrule.SetOrder(order);
		
		int intrulepoints = intrule.NPoints();
		REAL weight;
		std::vector< std::vector < int > > nodeIndices;
		REAL p=0;
		int side = 0;
		if (sp->Reference()->NCornerNodes()==3) {
			
			vtkTriangle * triang = vtkTriangle::New();
			for(side = 0 ; side < sp->Reference()->NCornerNodes(); side++)
			{
				sp->Reference()->CenterPoint(side, qsi);
				sp->ComputeSolution(qsi, sol, dsol, axes);
				sp->Reference()->X(qsi,X);
				pts[0] = X[0];
				pts[1] = X[1];
				pts[2] = 0;
				REAL p = sol[0];
				gValue = cel->Reference()->NodeIndex(side);
				points->InsertPoint(gValue, pts);
				triang->GetPointIds()->SetId(side, gValue);
				//                pressure->InsertNextValue(p);
				std::pair<int, double> tItem(gValue, p);
				pvalues.push_back(tItem);
				
				if(p > maxP) maxP = p;
				if(p < minP) minP = p;
				//            sp->Solution(qsi, 7, flow);
				//            std::pair<double, double> dpair(flow[0], flow[1]);
				//            std::pair<int, std::pair<double, double> > dItem(gValue, dpair);
				//            dvalues.push_back(dItem);
				
			}
			cells->InsertNextCell(triang);
		}
		if(sp->Reference()->NCornerNodes()==2)
		{
			vtkLine * line = vtkLine::New();
			for(side = 0 ; side < sp->Reference()->NCornerNodes(); side++)
			{
				sp->Reference()->CenterPoint(side, qsi);
				sp->ComputeSolution(qsi, sol, dsol, axes);
				sp->Reference()->X(qsi,X);
				pts[0] = X[0];
				pts[1] = X[1];
				pts[2] = 0;
				REAL p = sol[0];
				gValue = cel->Reference()->NodeIndex(side);
				points->InsertPoint(gValue, pts);
				line->GetPointIds()->SetId(side, gValue);
				//                pressure->InsertNextValue(p);
				std::pair<int, double> tItem(gValue, p);
				pvalues.push_back(tItem);
				
				if(p > maxP) maxP = p;
				if(p < minP) minP = p;
				//            if(sp->Reference()->NCornerNodes()<3) continue;
				//            sp->Solution(qsi, 7, flow);
				//            std::pair<double, double> dpair(flow[0], flow[1]);
				//            std::pair<int, std::pair<double, double> > dItem(gValue, dpair);
				//            dvalues.push_back(dItem);
				
			}
			lines->InsertNextCell(line);
		}
	}
	
	int npoints;
	npoints = points->GetNumberOfPoints();
	pressure->SetNumberOfValues(npoints);
	std::vector< std::pair< int, double > >::iterator It;
	
	for(It = pvalues.begin(); It != pvalues.end(); It++)
	{
		pressure->SetValue(It->first,It->second);
	}
	//compute vertical scale;
	double bounds[6];
	points->GetBounds(bounds);
	double lw, ll;
	lw = bounds[1]-bounds[0];
	ll = bounds[3]-bounds[2];
	
	double scale = (ll+lw)/2;
	double pw = maxP - minP;
	
	for(It = pvalues.begin(); It != pvalues.end(); It++)
	{
		double coord[3];
		points->GetPoint(It->first,coord);
		coord[2] = (It->second) * scale/pw;
		points->SetPoint(It->first, coord);
	}
	
	
	// Create the dataset. In this case, we create a vtkPolyData
	vtkPolyData* polydata = vtkPolyData::New();
	// Assign points and cells
	polydata->SetPoints(points);
	polydata->SetPolys(cells);
	polydata->SetLines(lines);
	//    SetStrips(cells);
	// Assign scalars
	polydata->GetPointData()->SetScalars(pressure);
	vtkPolyDataMapper * mapper = vtkPolyDataMapper::New();
	
	
	vtkPolyDataWriter * Writer = vtkPolyDataWriter::New();
	Writer->SetInput(polydata);
	Writer->SetFileName("PolyPhil.vtk");
	Writer->Write();
	
	mapper->SetInput(polydata);
	mapper->SetScalarRange(minP, maxP);
	
	vtkExtractEdges * edges = vtkExtractEdges::New();
	edges->SetInput(polydata);
	vtkPolyDataMapper * edgesMapper = vtkPolyDataMapper::New();
	edgesMapper->SetInput(edges->GetOutput());
	m_EdgesActor = vtkActor::New();
	m_EdgesActor->SetMapper(edgesMapper);
	
	vtkPolyDataMapper * mapContour = vtkPolyDataMapper::New();
	vtkContourFilter * contour = vtkContourFilter::New();
	contour->SetInput(polydata);
	contour->GenerateValues(10, pressure->GetRange());
	mapContour->SetInput(contour->GetOutput());
	m_ContActor = vtkActor::New();
	m_ContActor->SetMapper(mapContour);
	
	
	
	//    vtkPolyDataMapper * elevMapper = vtkPolyDataMapper::New();
	//    vtkElevationFilter * elev = vtkElevationFilter::New();
	//    elev->SetInput(polydata);
	//    elev->SetScalarRange(pressure->GetRange());
	//    elevMapper->SetInputConnection(elev->GetOutputPort());
	//    vtkActor * elevActor = vtkActor::New();
	//    elevActor->SetMapper(elevMapper);
	
	
	
	
	vtkScalarBarActor * scalarBar = vtkScalarBarActor::New();
	scalarBar->SetTitle("Pressão");
	scalarBar->SetLookupTable(mapper->GetLookupTable());
	vtkCoordinate * coord = scalarBar->GetPositionCoordinate();
	coord->SetCoordinateSystemToNormalizedViewport();
	coord->SetValue(0.8, 0.1);
	
	
	
	vtkActor * actor = vtkActor::New();
	actor->SetMapper(mapper);
	
	vtkRenderer* ren = m_vtk->GetRenderer();
	ren->AddActor(actor);
	//    ren->AddActor(elevActor);
	if(m_ShowMesh->Checked) ren->AddActor(m_EdgesActor);
	if(m_ShowIsoLines->Checked) ren->AddActor(m_ContActor);
	ren->AddActor2D(scalarBar);
	ren->ResetCamera();
	m_vtk->Invalidate();
	
	delete cmesh;
	
}
*/