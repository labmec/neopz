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
cl2 = 0.1;
cl3 = 10.0;
cl4 = 50.0;
cl5 = 250.0;

// well location
wx = 200.0;
wy = 200.0;
wz = -10.0;

// well dimensions
radius = 0.1;
outer = 10;
length = 10.0;

// reservoir box dimensions
x_length = 500.0;
y_length = 500.0;
z_length = 150.0;

// side-burden box dimensions
sb_x_length = 1000.0;
sb_y_length = 1000.0;
sb_z_length = 500.0;

// mesh controls
alpha = 1.5;
n_radial = 10;
n_azimuthal = 4;
n_axial = 8; 


// Circle & Surrounding structured-quad region
Point(1) = {0+wx, 0+wy, 0+wz, cl2};

// well bore points
Point(2) = {0+wx,  radius+wy, 0+wz, cl2};
Point(3) = {radius+wx,  0+wy, 0+wz, cl2};
Point(4) = {0+wx, -radius+wy, 0+wz, cl2};
Point(5) = {-radius+wx, 0+wy, 0+wz, cl2};

// well bore region points
Point(6) = {0+wx,  outer+wy, 0+wz, cl3};
Point(7) = {outer+wx,  0+wy, 0+wz, cl3};
Point(8) = {0+wx, -outer+wy, 0+wz, cl3};
Point(9) = {-outer+wx, 0+wy, 0+wz, cl3};


Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

Circle(5) = {6, 1, 7};
Circle(6) = {7, 1, 8};
Circle(7) = {8, 1, 9};
Circle(8) = {9, 1, 6};

Line(9)  = {2, 6};
Line(10) = {3, 7};
Line(11) = {4, 8};
Line(12) = {5, 9};

Transfinite Line {1,2,3,4,5,6,7,8} = n_azimuthal;
Transfinite Line {9,10,11,12} = n_radial Using Progression alpha;



//Using Progression 1.1

// Each region which to be independently meshed must have a line loop
// Regions which will be meshed with Transfinite Surface must have 4 lines
// and be labeled in CCW order, with the correct orientation of each edge

Line Loop(1) = {1, 10, -5, -9}; // ++ side of quad region - note ordering
Line Loop(2) = {2, 11, -6, -10}; // +- side of quad region - note ordering
Line Loop(3) = {3, 12, -7, -11}; // -- side of quad region - note ordering
Line Loop(4) = {4, 9, -8, -12}; // -+ side of quad region - note ordering

Plane Surface(1) = {1}; // ++ side of quad
Plane Surface(2) = {2}; // +- side of quad
Plane Surface(3) = {3}; // -- side of quad
Plane Surface(4) = {4}; // -+ side of quad

// Mesh these surfaces in a structured manner
 Transfinite Surface{1,2,3,4};



// Turn into quads (optional, but Transfinite Surface looks best with quads)
// Recombine Surface {1,2,3,4}; 
// Turn outer region into unstructured quads (optional) 
// Recombine Surface {1};


// Change layer to increase z subdivision
Extrude {0, 0, length} { Surface{1,2,3,4}; Layers{n_axial};}
//Extrude { {0, length, 0},{500,500,100}, -Pi/40.0} { Surface{1,2,3,4}; Layers{n_axial};}


////////////////////////////////////////////////////////////////////////////
/// Well box

well_index = 1;

// well closed regions

Line Loop(5) = {1, 2, 3, 4}; // clossing the well
Plane Surface(5) = {5}; // clossing the well

Line Loop(6) = {14, 36, 80, 58}; // clossing the well
Plane Surface(6) = {6}; // clossing the well

//well_closed[well_index] = newreg;
//Surface Loop(well_closed[well_index]) = {5,6};

//Compound Surface (7) = well_closed[well_index];

 Transfinite Surface{5,6};
// Recombine Surface {5,6};

// outer well box

 Point(1001) = {-x_length/2.0, -y_length/2.0, -z_length/2.0, cl4};
 Point(1002) = { x_length/2.0, -y_length/2.0, -z_length/2.0, cl4};
 Point(1003) = { x_length/2.0,  y_length/2.0, -z_length/2.0, cl4};
 Point(1004) = {-x_length/2.0,  y_length/2.0, -z_length/2.0, cl4};

 Point(1005) = {-x_length/2.0, -y_length/2.0, z_length/2.0, cl4};
 Point(1006) = { x_length/2.0, -y_length/2.0, z_length/2.0, cl4};
 Point(1007) = { x_length/2.0,  y_length/2.0, z_length/2.0, cl4};
 Point(1008) = {-x_length/2.0,  y_length/2.0, z_length/2.0, cl4};

Line(2001) = {1001,1002};
Line(2002) = {1002,1003};
Line(2003) = {1003,1004};
Line(2004) = {1004,1001};

Line(2005) = {1005,1006};
Line(2006) = {1006,1007};
Line(2007) = {1007,1008};
Line(2008) = {1008,1005};

Line(2009) = {1005,1001};
Line(2010) = {1006,1002};
Line(2011) = {1007,1003};
Line(2012) = {1008,1004};

 Line Loop(3001) = {2001, 2002, 2003, 2004}; // Bottom
 Line Loop(3002) = {2005, 2006, 2007, 2008}; // Top
 Line Loop(3003) = {2001, -2010, -2005, 2009}; // South
 Line Loop(3004) = {2002, -2011, -2006, 2010}; // East
 Line Loop(3005) = {2003, -2012, -2007, 2011}; // North
 Line Loop(3006) = {2004, -2009, -2008, 2012}; // West

 Plane Surface(3001) = {3001}; // Bottom unstructured region
 Plane Surface(3002) = {3002}; // Top unstructured region
 Plane Surface(3003) = {3003}; // South unstructured region
 Plane Surface(3004) = {3004}; // East unstructured region
 Plane Surface(3005) = {3005}; // North unstructured region
 Plane Surface(3006) = {3006}; // West unstructured region


// Recombine Surface {3001,3002,3003,3004,3005,3006};

// Build the well box volume

Surface Loop(4000) = {1,2,3,4,5,6,29,51,73,95,34,56,78,100,3001,3002,3003,3004,3005,3006};
Volume(6) = {4000} ;

////////////////////////////////////////////////////////////////////////////

// side-burden rocks

 Point(10001) = {-sb_x_length/2.0, -sb_y_length/2.0, -sb_z_length/2.0, cl5};
 Point(10002) = { sb_x_length/2.0, -sb_y_length/2.0, -sb_z_length/2.0, cl5};
 Point(10003) = { sb_x_length/2.0,  sb_y_length/2.0, -sb_z_length/2.0, cl5};
 Point(10004) = {-sb_x_length/2.0,  sb_y_length/2.0, -sb_z_length/2.0, cl5};

 Point(10005) = {-sb_x_length/2.0, -sb_y_length/2.0, sb_z_length/2.0, cl5};
 Point(10006) = { sb_x_length/2.0, -sb_y_length/2.0, sb_z_length/2.0, cl5};
 Point(10007) = { sb_x_length/2.0,  sb_y_length/2.0, sb_z_length/2.0, cl5};
 Point(10008) = {-sb_x_length/2.0,  sb_y_length/2.0, sb_z_length/2.0, cl5};

Line(20001) = {10001,10002};
Line(20002) = {10002,10003};
Line(20003) = {10003,10004};
Line(20004) = {10004,10001};

Line(20005) = {10005,10006};
Line(20006) = {10006,10007};
Line(20007) = {10007,10008};
Line(20008) = {10008,10005};

Line(20009) = {10005,10001};
Line(20010) = {10006,10002};
Line(20011) = {10007,10003};
Line(20012) = {10008,10004};

 Line Loop(30001) = {20001, 20002, 20003, 20004}; // Top
 Line Loop(30002) = {20005, 20006, 20007, 20008}; // Bottom
 Line Loop(30003) = {20001, -20010, -20005, 20009}; // South
 Line Loop(30004) = {20002, -20011, -20006, 20010}; // East
 Line Loop(30005) = {20003, -20012, -20007, 20011}; // North
 Line Loop(30006) = {20004, -20009, -20008, 20012}; // West

 Plane Surface(30001) = {30001}; // Top unstructured region
 Plane Surface(30002) = {30002}; // Bottom unstructured region
 Plane Surface(30003) = {30003}; // South unstructured region
 Plane Surface(30004) = {30004}; // East unstructured region
 Plane Surface(30005) = {30005}; // North unstructured region
 Plane Surface(30006) = {30006}; // West unstructured region


// Recombine Surface {30001,30002,3003,30004,30005,30006};

// Build the well box volume

Surface Loop(40000) = {3001,3002,3003,3004,3005,3006,30001,30002,30003,30004,30005,30006};
Volume(7) = {40000} ;

////////////////////////////////////////////////////////////////////////////

// Tagging boundary conditions

Physical Surface("well") = {21,43,65,87};
Physical Surface("well_closed") = {5,6};
Physical Volume("Reservoir") = {1, 2, 4, 3, 5, 6};
Physical Volume("well_region") = {1, 2, 4, 3, 5};

Physical Surface("Reservoir_bottom") = {3001};
Physical Surface("Reservoir_top") = {3002};
Physical Surface("Reservoir_South") = {3003};
Physical Surface("Reservoir_East") = {3004};
Physical Surface("Reservoir_North") = {3005};
Physical Surface("Reservoir_West") = {3006};

Physical Volume("side_burden") = {7};

Physical Surface("side_burden_bottom") = {30001};
Physical Surface("side_burden_top") = {30002};
Physical Surface("side_burden_laterals") = {30003,30004,30005,30006};


