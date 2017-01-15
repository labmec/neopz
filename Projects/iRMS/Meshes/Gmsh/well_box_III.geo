// ---- Cylinder wellbore Region Gmsh scritp ----
// 2D_cylinder_tutorial.geo
// Creates a mesh with an inner structured-quad region and 
// an outer unstructured tri region
//
// Created 10/01/2016 by Omar Duran
// Labmec, State University of Campinas
// --------------------------------------------

Include "drill_producer.geo";

// Settings
ExpertMode = 1;
Mesh.Algorithm3D = 1 ;

well_bores = {};
well_regions = {};

// Gmsh allows variables; these will be used to set desired
// element sizes at various Points
cl1 = 1;
cl2 = 0.1;
cl3 = 10.0;
cl4 = 100.0;
cl5 = 1000.0;


////////////////////////////////////////////////////////////////////////////
// Drill producer 1 
////////////////////////////////////////////////////////////////////////////

// new well data
well_index = 1;

// mesh controls on wellbore region
alpha = 1.1;
n_radial = 10;
n_azimuthal = 4;
n_axial = 10; 

// well location
wx = 0.0;
wy = -50.0;
wz = 0.0;

// Geometry well and wellbore region dimensions
radius = 0.1;
length = 100.0;
outer = 40;
angle = Pi/2.0;
beta = 0.0;

Call DrillProducer;




////////////////////////////////////////////////////////////////////////////
// Reservoir rocks
////////////////////////////////////////////////////////////////////////////

// reservoir box dimensions
x_length = 1000.0;
y_length = 1000.0;
z_length = 100.0;

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


//Recombine Surface {3001,3002,3003,3004,3005,3006};

// Build the well box volume

reservoir_region[] = {well_regions[],3001,3002,3003,3004,3005,3006};
Surface Loop(4000) = reservoir_region[];
Volume(6) = {4000} ;

//Transfinite Surface "*";
//Recombine Surface "*";
//Recombine Volume "*";

////////////////////////////////////////////////////////////////////////////
// Side-burden rocks
////////////////////////////////////////////////////////////////////////////

// side-burden box dimensions
sb_x_length = 5000.0;
sb_y_length = 5000.0;
sb_z_length = 1000.0;

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

Surface Loop(40000) = {3001,3002,3003,3004,3005,3006,30001,30002,30003,30004,30005,30006};
Volume(7) = {40000} ;


////////////////////////////////////////////////////////////////////////////
// Mark physical entities
////////////////////////////////////////////////////////////////////////////

Physical Volume("Reservoir") = {1, 2, 4, 3, 5, 6};


Physical Surface("Reservoir_bottom") = {3001};
Physical Surface("Reservoir_top") = {3002};
Physical Surface("Reservoir_South") = {3003};
Physical Surface("Reservoir_East") = {3004};
Physical Surface("Reservoir_North") = {3005};
Physical Surface("Reservoir_West") = {3006};

Physical Volume("side_burden") = {7};
Physical Surface("side_burden_bottom") = {30001};
Physical Surface("side_burden_top") = {30002};
Physical Surface("side_burden_South") = {30003};
Physical Surface("side_burden_East") = {30004};
Physical Surface("side_burden_North") = {30005};
Physical Surface("side_burden_West") = {30006};


