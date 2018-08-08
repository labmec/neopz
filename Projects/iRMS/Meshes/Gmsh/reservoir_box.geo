// ---- iMRS reservoir geometry Gmsh scritp ----
// Creates a mesh with an inner structured-quad region and 
// an outer unstructured triangle region
//
// Created 10/01/2016 by Omar Duran
// Labmec, State University of Campinas
// --------------------------------------------

// Settings
ExpertMode = 1;
Mesh.Optimize = 1;
//Mesh.CharacteristicLengthExtendFromBoundary = 1;

Include "CadReservoir.geo";
Include "BoxReservoir.geo";
Include "BoxSideBurden.geo";
Include "BuildOmegas.geo";
Include "drill_well.geo";
Include "drill_well_box.geo";
Include "PhysicalEntities.geo";

well_index = 0;
well_p_lids = {};
well_i_lids = {};

well_p_bores = {};
well_p_regions = {};
well_p_v_regions = {};

well_i_bores = {};
well_i_regions = {};
well_i_v_regions = {};


geomechanicQ = 0;
dimension = 3;
nolinearQ = 0;
CADReservoirQ = 0;

xzQ = 0;
hexahedronsWQ = 0;
hexahedronsRQ = 1;
hexahedronsSBQ = 0;

If (nolinearQ == 1)
Mesh.ElementOrder = 2;
Mesh.SecondOrderLinear = 0;
EndIf

If (hexahedronsWQ == 1 || hexahedronsRQ == 1 || hexahedronsSBQ == 1)
Mesh.Algorithm3D = 6 ;
Else
Mesh.Algorithm3D = 1 ;
EndIf

// Gmsh allows variables; these will be used to set desired
// element sizes at various Points
cl1 = 1;
cl2 = 0.1;
cl3 = 20.0;
cl4 = 200.0;
cl5 = 5000.0;

////////////////////////////////////////////////////////////////////////////
// reservoir region geometry
////////////////////////////////////////////////////////////////////////////

// reservoir box dimensions
x_length = 100.0;
y_length = 100.0;
z_length = 100.0;


////////////////////////////////////////////////////////////////////////////
// side-burden region geometry
////////////////////////////////////////////////////////////////////////////

// side-burden box dimensions
sb_x_length = 100000.0;
sb_y_length = 100000.0;
sb_z_length = 20000.0;


sb_v_shift = 0.0;

If(dimension == 2)
sb_x_length = 100000.0;
sb_y_length = 20000.0;
sb_z_length = 10000.0;
EndIf

////////////////////////////////////////////////////////////////////////////
// reservoir rock
////////////////////////////////////////////////////////////////////////////
If(CADReservoirQ == 1)
Call ReservoirCAD;
Else
Call ReservoirBox;
EndIf

////////////////////////////////////////////////////////////////////////////
// side-burden rock
////////////////////////////////////////////////////////////////////////////
If (geomechanicQ == 1)
Call SideBurdenBox;
EndIf

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
// Computational domain, boundaries:: Tagging physical entities
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

Call DefineOmegas;
Call DrawBoundaries;
