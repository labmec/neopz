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
// soil region geometry
////////////////////////////////////////////////////////////////////////////

// box dimensions
x_length = 2.0;
y_length = 3.0;
z_length = 1.0;

wcx = 0.0;
wcy = 0.0;
wcz = 0.0;

Call RockBox;

// rock box dimensions
x_length = 1.0;
y_length = 3.0;
z_length = 1.0;

wcx = -1.5;
wcy = 0.0;
wcz = 0.0;

Call RockBox;


If (dimension == 2)

// rock
ll1 = newll; Line Loop(ll1) = {5,-4,7,8};
s1  = news; Plane Surface(s1) = {ll1};
ll2 = newll; Line Loop(ll2) = {1,2,3,4};
s2  = news; Plane Surface(s2) = {ll2};
rock[] = {s1,s2};

Transfinite Line{7} = 8;
Transfinite Line{-3} = 8 Using Progression 1.5;
Transfinite Line{4,8} = 8 Using Progression 1.5;

If (hexahedronsRQ == 1)
Transfinite Surface {s1,s2};
Recombine Surface {s1,s2};
EndIf

Physical Surface("Rock") = {rock[]};
Physical Line("bottom") = {1,5};
Physical Line("lateral") = {8,2};
Physical Line("top") = {7};
Physical Line("top_null") = {3};

Else

// rock
sl1 = newsl; Surface Loop(sl1) = {45,44,24,43,48,47};
v1  = newv; Volume(v1) = {sl1} ;
sl2 = newsl; Surface Loop(sl2) = {19,20,21,22,24,23};
v2  = newv; Volume(v2) = {sl2} ;
rock[] = {v1,v2};


If (hexahedronsRQ == 1)
Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume {v1,v2};
Recombine Volume {v1,v2};
EndIf

Physical Volume("Rock") = {rock[]};
Physical Surface("bottom") = {21,45};
Physical Surface("lateral") = {19,20,22,43,44,48};
Physical Surface("top") = {47};
Physical Surface("top_null") = {23};


EndIf

Coherence;
