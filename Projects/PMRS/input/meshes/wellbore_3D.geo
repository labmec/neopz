
////////////////////////////////////////////////////////////////
// 3D wellbore and Reservoir
// Created 18/06/2018 by Manouchehr Sanei
// Labmec, State University of Campinas, Brazil
////////////////////////////////////////////////////////////////

// Parameteres
wr = 0.1;
or = 1.0;
orc = 2.5;
fl = 4.0;
h = 0.4;

nr = 10; 
na = 5;
nh = 1;

clr = 2.0;
cli = 1.0;
clt = 2.0;

IsquadQ = 1;

// Re
Point(1) = {-fl, -fl, 0, clr};
Point(2) = { fl, -fl, 0, clt};
Point(3) = { fl,  fl, 0, clt};
Point(4) = {-fl,  fl, 0, clr};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Wb 
Point(5)  = {0,   0, 0, cli};
Point(6)  = {0,  wr, 0, cli};
Point(7)  = {0, -wr, 0, cli};
Point(8)  = {0,  or, 0, cli};
Point(9)  = {0, -or, 0, cli};
Point(10) = {wr,  0, 0, cli};
Point(11) = {-wr, 0, 0, cli};
Point(12) = {or,  0, 0, cli};
Point(13) = {-or, 0, 0, cli};

Point(14)  = {0,  orc, 0, cli};
Point(15)  = {0, -orc, 0, cli};
Point(16) = {orc,  0, 0, cli};
Point(17) = {-orc, 0, 0, cli};

Circle(5) = {7, 5, 10};
Circle(6) = {6, 5, 11};
Circle(7) = {8, 5, 13};
Circle(8) = {9, 5, 12};

Line(9)  = {6, 8};
Line(10) = {7, 9};
Line(11) = {10, 12};
Line(12) = {11, 13};

Circle(13) = {10, 5, 6};
Circle(14) = {11, 5, 7};
Circle(15) = {13, 5, 9};
Circle(16) = {12, 5, 8};

Transfinite Line {5,6,7,8,13,14,15,16} = na; 
Transfinite Line {9,10,11,12} = nr Using Progression 1.5;  

Line Loop(1) = {1, 2, 3, 4, 7, 16, 8, 15};
Line Loop(2) = {10, 8, -11, -5};
Line Loop(3) = {7, -12, -6, 9};
Line Loop(4) = {-10, -14, 12, 15};
Line Loop(5) = {16, -9, -13, 11};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};

fixed_x_points[]={14,15};
fixed_y_points[]={16,17};
Point{fixed_x_points[],fixed_y_points[]} In Surface{1};


Transfinite Surface{2,3,4,5};

If(IsquadQ)
 Recombine Surface {1};
 Recombine Surface {2,3,4,5};
EndIf

// Height
Extrude {0, 0, h} { Surface{1,2,3,4,5}; Layers{nh}; Recombine;}

// BC
Physical Volume("Omega")     = {1,2,4,3,5};
Physical Surface("FarfieldXr") = {33};
Physical Surface("FarfieldXl") = {41};
Physical Surface("FarfieldYt") = {37};
Physical Surface("FarfieldYd") = {29};
Physical Surface("Topsurf")   = {1,2,3,4,5};
Physical Surface("Downsurf") = {58,124,80,102,146};
Physical Surface("Hole")     = {115,79,97,141};

Physical Point("fixed_x") = {fixed_x_points[]};
Physical Point("fixed_y") = {fixed_y_points[]};

Coherence Mesh;
