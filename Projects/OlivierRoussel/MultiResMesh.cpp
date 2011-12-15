//$Id: MultiResMesh.cpp,v 1.6 2009-11-04 14:13:24 fortiago Exp $

#include "TPZFakeFunction.h"
#include "tpzautopointer.h"
#include "pzfunction.h"
#include "pzeuler.h"
#include "pzmatrix.h"
#include "tpzquadratictrig.h"
#include "tpzquadraticquad.h"
#include "tpzgeoelrefpattern.h"
#include "tpzquadraticline.h"
#include "pzgeoel.h"
#include <sstream>
#include "TPZGeoElement.h"
#include "TPZGeoCube.h"
#include "pzshapecube.h"
#include "TPZRefCube.h"
#include "pzshapelinear.h"
#include "TPZGeoLinear.h"
#include "TPZRefLinear.h"
#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "pzshapetriang.h"
#include "pzreftriangle.h"
#include "pzgeotriangle.h"
#include "pzshapeprism.h"
#include "pzrefprism.h"
#include "pzgeoprism.h"
#include "pzshapetetra.h"
#include "pzreftetrahedra.h"
#include "pzgeotetrahedra.h"
#include "pzshapepiram.h"
#include "pzrefpyram.h"
#include "pzgeopyramid.h"
#include "pzrefpoint.h"
#include "pzgeopoint.h"
#include "pzshapepoint.h"
#include "pzvec.h"
#include "pzcmesh.h"
#include "pzdebug.h"
#include "pzcheckgeom.h"
#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgeoelside.h"
#include "pzcompel.h"
#include "TPZCompElDisc.h"
#include "pzmatrix.h"
#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZParFrontStructMatrix.h"
// #include "TPZParFrontMatrix.h"
#include "TPZFrontNonSym.h"
#include "pzbdstrmatrix.h"
#include "pzblockdiag.h"
#include "TPZSpStructMatrix.h"
#include "TPZCopySolve.h"
#include "TPZStackEqnStorage.h"
#include "pzbstrmatrix.h"
#include "pzstepsolver.h"
#include "pzbndcond.h"
#include "pzpoisson3d.h"
#include "pzvisualmatrix.h"
#include "TPZGeoElement.h"
#include "TPZGeoCube.h"
#include "pzshapecube.h"
#include "TPZRefCube.h"
#include "pzshapelinear.h"
#include "TPZGeoLinear.h"
#include "TPZRefLinear.h"
#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "pzshapetriang.h"
#include "pzreftriangle.h"
#include "pzgeotriangle.h"
#include "pzshapeprism.h"
#include "pzrefprism.h"
#include "pzgeoprism.h"
#include "pzshapetetra.h"
#include "pzreftetrahedra.h"
#include "pzgeotetrahedra.h"
#include "pzshapepiram.h"
#include "pzrefpyram.h"
#include "pzgeopyramid.h"
#include "pzrefpoint.h"
#include "pzgeopoint.h"
#include "pzshapepoint.h"
#include "tpzchangeel.h"
#include "TPZInterfaceEl.h"
#include <time.h>
#include <stdio.h>
#include <set>
#include "pzeuler.h"
#include "pzreadmeshhr.h"

using namespace std;

void DefineRefinementPattern ( TPZGeoMesh * GMesh );

TPZGeoMesh * CreateCoarseMesh(int nLevel){

  // Read coarse mesh file
  TPZReadMeshHR reader ( "LaraMesh.txt" );

  // Create an object to store the geometrical mesh
  TPZGeoMesh * gmesh = NULL;

  // Initialize object from read data
  gmesh = reader.readGeoMesh();

  DefineRefinementPattern ( gmesh );

  TPZManVector<TPZGeoEl*> children;
  for(int idiv = 0; idiv < nLevel; idiv++){
    const int nel = gmesh->NElements();
    for(int iel = 0; iel < nel; iel++){
      TPZGeoEl * gel = gmesh->ElementVec()[iel];
      if(!gel->HasSubElement())
      {
    	  gel->Divide(children);
    	  for(int is = 0; is < children.NElements(); is++){
    		  children[is]->SetRefPattern( gel->GetRefPattern() );
    	  }
      }
    }///for iel
  }///for idiv

  return gmesh;

}///method

TPZCompMesh *CreateMeshMultires(TPZGeoMesh * gmesh){

  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(3);

  TPZAutoPointer<TPZMaterial> mat = new TPZEulerEquation(1,1.4);
  cmesh->InsertMaterialObject(mat);
  TPZFMatrix val1,val2;
  cmesh->InsertMaterialObject(mat->CreateBC(mat,-1,TPZEulerEquation::EFreeSlip,val1,val2));

  TPZCompMesh::SetAllCreateFunctionsDiscontinuous();
  cmesh->SetDefaultOrder(0);
  TPZCompElDisc::SetgOrder(0);

  cmesh->AutoBuild();
  TPZCompElDisc::SetTotalOrderShape(cmesh);
  TPZAutoPointer<TPZFunction> fakefunc = new TPZFakeFunction();
  for(int i = 0; i < cmesh->NElements(); i++){
    TPZCompElDisc * disc = dynamic_cast<TPZCompElDisc*>(cmesh->ElementVec()[i]);
    if(disc){
      disc->SetExternalShapeFunction(fakefunc);
    }
  }
  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();
  cmesh->ExpandSolution();

  return cmesh;

}

/** initial solution */
void InitialSolutionMultires(TPZFMatrix &InitialSol, TPZCompMesh * cmesh){
  InitialSol.Redim(cmesh->NEquations(),1);
  InitialSol.Zero();
  for(int iel = 0; iel < cmesh->NElements(); iel++){
    TPZCompEl * cel = cmesh->ElementVec()[iel];
    if(!cel) continue;
    TPZCompElDisc * disc = dynamic_cast<TPZCompElDisc*>(cel);
    if(!disc) continue;
    if(disc->NConnects() == 0) continue;
    int bl = disc->Connect(0).SequenceNumber();
    int blpos = cmesh->Block().Position(bl);
    int blocksize = cmesh->Block().Size(bl);

    TPZGeoEl * gel = cel->Reference();
    TPZVec<REAL> xi(3), xVec(3);
    gel->CenterPoint(gel->NSides()-1,xi);
    gel->X(xi,xVec);
    double x = xVec[0];
    double y = xVec[1];
    double z = xVec[2];
    double u = 0.125;

    double xCircle = 0.25;
    double yCircle = 0.5;
    double zCircle = 0.5;

    double R = 0.1;

    if( (x-xCircle)*(x-xCircle)+(y-yCircle)*(y-yCircle)+(z-zCircle)*(z-zCircle) <= R*R ) u = 1.;

    InitialSol(blpos+blocksize-20+0,0) = u;
    InitialSol(blpos+blocksize-20+1,0) = 0.;
    InitialSol(blpos+blocksize-20+2,0) = 0.;
    InitialSol(blpos+blocksize-20+3,0) = 0.;
    InitialSol(blpos+blocksize-20+4,0) = 0.;

  }///for iel

  TPZVec<REAL> celerity(3,0.);
  celerity[0] = 1.;
#ifdef LinearConvection
  TPZEulerEquation::SetLinearConvection(cmesh, celerity);
#endif
}

void DefineRefinementPattern ( TPZGeoMesh * GMesh )
{
  // --- Define the refinement patterns for tetrahedra -------

  // Read the refinement pattern
  string tetraFile ( "LaraRefineTetrahedron.txt" );

  // Create and initialize an object to store the pattern
	TPZRefPattern *tetraPattern = new TPZRefPattern (*GMesh);
	tetraPattern->SetName(tetraFile );

  // Include the pattern info into the mesh
  tetraPattern->InsertPermuted();

  // Remove this object
  delete tetraPattern;

  // --- Define the refinement patterns for triangles -------

  // Read the refinement pattern
  string triangleFile ( "LaraRefineTriangle.txt" );

  // Create and initialize an object to store the pattern
	TPZRefPattern * trianglePattern = new TPZRefPattern (*GMesh);
	trianglePattern->SetName(triangleFile );

  // Include the pattern info into the mesh
  trianglePattern->InsertPermuted();

  // Remove this object
  delete trianglePattern;

  // Take the map of every refinement pattern for a tetrahedron defined into the mesh
  // the map below provides a relation between an index and a pointer to an object of type refinement pattern
	map < int, TPZAutoPointer < TPZRefPattern > > mapOfTetraRefPattern;   // = GMesh->RefPatternList ( eltype );
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
  TPZAutoPointer < TPZRefPattern > tetraRefPattern = it->second;


  // --- Divide elements ------------------------------------
  for (int i = 0; i < GMesh->NElements(); i++ )
  {
    TPZGeoEl * gel = GMesh->ElementVec()[ i ];
    if ( ! gel || gel->Dimension() != 3 ) continue;
    gel->SetRefPattern ( tetraRefPattern );
  }

}

double ComputeTimeStep(double CFL, int Level, TPZGeoMesh * gmesh){
  double MinR = 1e12;
  for(int iel = 0; iel < gmesh->NElements(); iel++){
    TPZGeoEl * gel = gmesh->ElementVec()[iel];
    if(!gel) continue;
    if(gel->Dimension() == 3){
      int currentLevel = gel->Level();
      if(Level == currentLevel){
        double r = gel->ElementRadius();
        if(r < MinR) MinR = r;
      }
    }
  }///for

  double timeStep = CFL*MinR;
  return timeStep;
}
