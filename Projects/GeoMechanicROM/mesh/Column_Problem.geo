// ---- Geomechanic reduced coupling geometry Gmsh scritp ----
// Creates a mesh with an inner structured-quad region and 
// an outer unstructured triangle region
//
// Created 10/01/2016 by Omar Duran
// Labmec, State University of Campinas
// --------------------------------------------

Include "rock_box.geo";

dimension = 2;
cl = 1.0;
hexahedronsRQ = 1;

bottoms = {};
tops = {};
laterals = {};

////////////////////////////////////////////////////////////////////////////
// rock region geometry
////////////////////////////////////////////////////////////////////////////

// rock box dimensions
x_length = 0.1;
y_length = 1.0;
z_length = 0.1;

wcx = 0.0;
wcy = 0.0;
wcz = 0.0;

Call RockBox;


If (dimension == 2)

// rock
ll1 = newll; Line Loop(ll1) = {1,2,3,4};
s1  = news; Plane Surface(s1) = {ll1};
rock[] = {s1};

Transfinite Line{1,3} = 1 Using Progression 1.0;
Transfinite Line{2,4} = 11 Using Progression 1.0;

If (hexahedronsRQ == 1)
Transfinite Surface {s1};
Recombine Surface {s1};
EndIf

Physical Surface("Rock") = {rock[]};
Physical Line("bottom") = {1};
Physical Line("lateral") = {2,4};
Physical Line("top") = {3};
Physical Line("top_null") = {7};

Else

// rock
sl1 = newsl; Surface Loop(sl1) = {19,20,21,22,23,24};
v1  = newv; Volume(v1) = {sl1} ;
rock[] = {v1};

Transfinite Surface "*";
Transfinite Volume {v1};
If (hexahedronsRQ == 1)
Recombine Surface "*";
Recombine Volume {v1};
EndIf

Physical Volume("Rock") = {rock[]};
Physical Surface("bottom") = {21};
Physical Surface("lateral") = {19,20,24,22};
Physical Surface("top") = {23};
Physical Surface("top_null") = {47};


EndIf
