
/*********************************************************************
 *
 *  Gmsh control scritp
 *	Two dimensional examples with three embeded fractures
 *
 *********************************************************************/

Include "rectangle.geo";
Include "circle.geo";

Mesh.Algorithm = 1;
Mesh.CharacteristicLengthExtendFromBoundary = 0;
General.ExpertMode = 1;

// Controls 
CircleGeometryQ = 0;
QuadrilateralMeshQ = 0;

If (CircleGeometryQ == 1)
Call CircleDomain;
Else
Call RectangleDomain;
EndIf

