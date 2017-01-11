// ---- 2D Circular Cylinder Gmsh Tutorial ----
// 2D_cylinder_tutorial.geo
// Creates a mesh with an inner structured-quad region and 
// an outer unstructured tri region
//
// Created 10/01/2016 by Omar Duran
// Labmec, State University of Campinas
// --------------------------------------------

// Gmsh allows variables; these will be used to set desired
// element sizes at various Points
cl1 = 1;
cl2 = 1;
cl3 = 20;
cl4 = 20;

radius = 1.0;
outer = 20;
inner_box = 2.0*outer;
outer_box = 2.0*inner_box;
numinner = 10; 
extr = 1;

// Exterior (bounding box) of mesh
Point(1) = {-inner_box, -inner_box, 0, cl3};
Point(2) = { inner_box, -inner_box, 0, cl3};
Point(3) = { inner_box,  inner_box, 0, cl3};
Point(4) = {-inner_box,  inner_box, 0, cl3};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Circle & surrounding structured-quad region
Point(5) = {0,   0, 0, cl2};
Point(6) = {0,  radius, 0, cl2};
Point(7) = {0, -radius, 0, cl2};
Point(8) = {0,  outer, 0, cl2};
Point(9) = {0, -outer, 0, cl2};
Point(10) = {radius,  0, 0, cl2};
Point(11) = {-radius, 0, 0, cl2};
Point(12) = {outer,  0, 0, cl2};
Point(13) = {-outer, 0, 0, cl2};

Circle(5) = {7, 5, 10};
Circle(6) = {6, 5, 11};
Circle(7) = {8, 5, 13};
Circle(8) = {9, 5, 12};

Line(9)  = {6, 8};
Line(10) = {7, 9};
Line(11)  = {10, 12};
Line(12) = {11, 13};

Circle(13) = {10, 5, 6};
Circle(14) = {11, 5, 7};
Circle(15) = {13, 5, 9};
Circle(16) = {12, 5, 8};


Transfinite Line {5,6,7,8,13,14,15,16} = 4; // We want 40 points along each of these lines
Transfinite Line {9,10,11,12} = numinner Using Progression 1.5;    // And 10 points along each of these lines



//Using Progression 1.1

// Each region which to be independently meshed must have a line loop
// Regions which will be meshed with Transfinite Surface must have 4 lines
// and be labeled in CCW order, with the correct orientation of each edge
Line Loop(1) = {1, 2, 3, 4, 7, 16, 8, 15}; // Exterior
Line Loop(2) = {10, 8, -11, -5}; // RH side of quad region - note ordering
Line Loop(3) = {7, -12, -6, 9}; // LH side of quad region - note ordering
Line Loop(4) = {-10, -14, 12, 15}; // RH side of quad region - note ordering
Line Loop(5) = {16, -9, -13, 11}; // LH side of quad region - note ordering

Plane Surface(1) = {1}; // Outer unstructured region
Plane Surface(2) = {2}; // RH inner structured region
Plane Surface(3) = {3}; // LH inner structured region
Plane Surface(4) = {4}; // RH inner structured region
Plane Surface(5) = {5}; // LH inner structured region

// Mesh these surfaces in a structured manner
 Transfinite Surface{2,3,4,5};



// Turn into quads (optional, but Transfinite Surface looks best with quads)
// Recombine Surface {2,3,4,5}; 
// Turn outer region into unstructured quads (optional) 
// Recombine Surface {1};


// Change layer to increase z subdivision
//Extrude {0, 0, 10} { Surface{1,2,3,4,5}; Layers{1}; Recombine; Transfinite;}
Extrude {0, 0, 10} { Surface{1,2,3,4,5}; Layers{1};}
Transfinite Line {43,44,48,52,92,65,96,74,32,36,27,28} = 4;



////////////////////////////////////////////////////////////////////////////
/// Well box
// close well
Line Loop(6) = {5, 6, 13, 14}; // clossing the well
Plane Surface(6) = {6}; // clossing the well

Line Loop(7) = {63, 84, 105, 128}; // clossing the well
Plane Surface(7) = {7}; // clossing the well

 Transfinite Surface{6,7};
// Recombine Surface {6,7};

// outer well box

 Point(101) = {-outer_box, -outer_box, -20, cl4};
 Point(102) = { outer_box, -outer_box, -20, cl4};
 Point(103) = { outer_box,  outer_box, -20, cl4};
 Point(104) = {-outer_box,  outer_box, -20, cl4};

 Point(105) = {-outer_box, -outer_box, 30, cl4};
 Point(106) = { outer_box, -outer_box, 30, cl4};
 Point(107) = { outer_box,  outer_box, 30, cl4};
 Point(108) = {-outer_box,  outer_box, 30, cl4};

Line(201) = {101,102};
Line(202) = {102,103};
Line(203) = {103,104};
Line(204) = {104,101};

Line(205) = {105,106};
Line(206) = {106,107};
Line(207) = {107,108};
Line(208) = {108,105};

Line(209) = {105,101};
Line(210) = {106,102};
Line(211) = {107,103};
Line(212) = {108,104};

 Line Loop(301) = {201, 202, 203, 204}; // Bottom
 Line Loop(302) = {205, 206, 207, 208}; // Top
 Line Loop(303) = {201, -210, -205, 209}; // South
 Line Loop(304) = {202, -211, -206, 210}; // East
 Line Loop(305) = {203, -212, -207, 211}; // North
 Line Loop(306) = {204, -209, -208, 212}; // West

 Plane Surface(301) = {301}; // Bottom unstructured region
 Plane Surface(302) = {302}; // Top unstructured region
 Plane Surface(303) = {303}; // South unstructured region
 Plane Surface(304) = {304}; // East unstructured region
 Plane Surface(305) = {305}; // North unstructured region
 Plane Surface(306) = {306}; // West unstructured region


// Recombine Surface {301,302,303,304,305,306};

// Build the well box volume
Surface Loop(400) = {1,2,3,4,5,6,7,58,80,102,124,146,29,33,37,41,301,302,303,304,305,306};
Volume(6) = {400} ;

////////////////////////////////////////////////////////////////////////////

// Apply boundary conditions
// Note: Can change names later at top of .msh file
// Each boundary in gmsh must be labeled differently
// rename the boundaries manually in the resulting .msh file
//Physical Line("Bottom") = {1};
//Physical Line("Right")  = {2};
//Physical Line("Top")    = {3};
//Physical Line("Left")   = {4};
//Physical Line("Circle") = {5,6};
// Alternate version - make all 4 outer bounds part of the same B.C.:
//Physical Line("Char") = {1,2,3,4}; 

// IMPORTANT: "FLUID" MUST contain all fluid surfaces(2D)/volumes(3D)
//Physical Surface("FLUID") = {1,2,3};

Physical Surface("well") = {115, 79, 97,141};
Physical Surface("well_closed") = {7,6};
Physical Volume("Reservoir") = {1, 2, 4, 3, 5, 6};
Physical Volume("Reservoir_well") = {1, 2, 4, 3, 5};