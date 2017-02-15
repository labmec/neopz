// ---- Cylinder wellbore Region Gmsh scritp ----
// 2D_cylinder_tutorial.geo
// Creates a mesh with an inner structured-quad region and 
// an outer unstructured tri region
//
// Created 10/01/2016 by Omar Duran
// Labmec, State University of Campinas
// --------------------------------------------


Include "drill_producer.geo";
Include "drill_injector.geo";

well_lids = {};

well_p_bores = {};
well_p_regions = {};
well_p_v_regions = {};

well_i_bores = {};
well_i_regions = {};
well_i_v_regions = {};

// Settings
ExpertMode = 1;

<<<<<<< HEAD
dimension = 3;
nolinearQ = 0;
YReservoirQ = 1;

xzQ = 0;
hexahedronsQ = 0;
hexahedronsOutQ = 0;
=======
dimension = 2;
nolinearQ = 0;

xzQ = 0;
hexahedronsQ = 1;
hexahedronsOutQ = 1;
>>>>>>> iRMS_Biot

If (nolinearQ == 1)
Mesh.ElementOrder = 2;
Mesh.SecondOrderLinear = 0;
EndIf

If (hexahedronsOutQ == 1)
n_bc_res = 100;
n_bc_sb = 1000;
EndIf

If (hexahedronsQ == 1)
Mesh.Algorithm3D = 6 ;
Else
Mesh.Algorithm3D = 1 ;
EndIf


// Gmsh allows variables; these will be used to set desired
// element sizes at various Points
cl1 = 1;
cl2 = 0.1;
cl3 = 10.0;
<<<<<<< HEAD
cl4 = 500.0;
cl5 = 100.0;
=======
cl4 = 1000.0;
cl5 = 2000.0;
>>>>>>> iRMS_Biot

////////////////////////////////////////////////////////////////////////////
// reservoir region geometry
////////////////////////////////////////////////////////////////////////////

// reservoir box dimensions
x_length = 1000.0;
y_length = 1000.0;
z_length = 200.0;

////////////////////////////////////////////////////////////////////////////
// side-burden region geometry
////////////////////////////////////////////////////////////////////////////

// side-burden box dimensions
sb_x_length = 5000.0;
sb_y_length = 5000.0;
sb_z_length = 4000.0;

////////////////////////////////////////////////////////////////////////////
// well bore regions geometry
////////////////////////////////////////////////////////////////////////////

// mesh controls on wellbore region
alpha = 4.0;
n_radial = 4;
n_azimuthal = 3;
n_axial = 4; 

// Geometry well and wellbore region dimensions
<<<<<<< HEAD
radius = 5.0;
=======
radius = 0.1;
>>>>>>> iRMS_Biot
length = 100.0;
outer = 20;
angle = Pi/2.0;
beta = 0.0;

////////////////////////////////////////////////////////////////////////////
// Drill producer 1 
////////////////////////////////////////////////////////////////////////////

// new well data
well_index = 1;

// well location
wx = 0.0;
<<<<<<< HEAD
wy = 50.0;
wz = -50.0;
Call DrillProducer;
=======
wy = 0.0;
wz = 10.0;
//Call DrillProducer;
>>>>>>> iRMS_Biot


////////////////////////////////////////////////////////////////////////////
// Drill injector 1 
////////////////////////////////////////////////////////////////////////////

// new well data
well_index = 2;

// well location
wx = 400.0;
wy = 350.0;
wz = -50.0;
<<<<<<< HEAD
Call DrillInjector;
=======
//Call DrillInjector;
>>>>>>> iRMS_Biot

////////////////////////////////////////////////////////////////////////////
// Drill injector 2 
////////////////////////////////////////////////////////////////////////////

// new well data
well_index = 3;

// well location
wx = -400.0;
wy = -400.0;
wz = -50.0;
<<<<<<< HEAD
Call DrillInjector;
=======
//Call DrillInjector;
>>>>>>> iRMS_Biot


////////////////////////////////////////////////////////////////////////////
// Drill injector 3 
////////////////////////////////////////////////////////////////////////////

// new well data
well_index = 4;

// well location
wx = -400.0;
wy = 350.0;
wz = -50.0;
//Call DrillInjector;

////////////////////////////////////////////////////////////////////////////
// Drill injector 4 
////////////////////////////////////////////////////////////////////////////

// new well data
well_index = 5;

// well location
wx = 400.0;
wy = -400.0;
wz = -50.0;
//Call DrillInjector;

If (dimension == 3)

////////////////////////////////////////////////////////////////////////////
// Reservoir rocks
////////////////////////////////////////////////////////////////////////////

<<<<<<< HEAD
If(YReservoirQ == 1)

Include "YReservoir.geo";

line_id = newl-28;
face_id = news-6;

Printf ("Mesher:: line_id = %g", line_id);
Printf ("Mesher:: face_id = %g", face_id);

bottom[] = {face_id};
top[] = {face_id + 5};
south[] = {face_id + 4};
east[] = {face_id + 3};
north[] = {face_id + 2};
west[] = {face_id + 1};


line_id = 1;
res_edges_h[] = {
line_id,line_id+1,line_id+2,line_id+3,line_id+5,line_id+8,line_id+9,line_id+11};
res_edges_v[] = {
line_id+4,line_id+6,line_id+7,line_id+10};

Transfinite Line {res_edges_h[]} = 25.0;
Transfinite Line {res_edges_v[]} = 4.0;

res_boundaries[] = {top[],bottom[],east[],west[],north[],south[]};


Compound Surface(3001) = {bottom[]}; // Bottom unstructured region
Compound Surface(3002) = {top[]}; // Top unstructured region
Compound Surface(3003) = {south[]}; // South unstructured region
Compound Surface(3004) = {east[]}; // East unstructured region
Compound Surface(3005) = {north[]}; // North unstructured region
Compound Surface(3006) = {west[]}; // West unstructured region

reservoir_region[] = {well_p_regions[],well_i_regions[],res_boundaries[]};
Surface Loop(4000) = reservoir_region[];
Volume(6) = {4000} ;

Else

=======
>>>>>>> iRMS_Biot
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

reservoir_region[] = {well_p_regions[],well_i_regions[],3001,3002,3003,3004,3005,3006};
Surface Loop(4000) = reservoir_region[];
Volume(6) = {4000} ;

<<<<<<< HEAD

EndIf


=======
>>>>>>> iRMS_Biot
////////////////////////////////////////////////////////////////////////////
// Side-burden rocks
////////////////////////////////////////////////////////////////////////////

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

<<<<<<< HEAD
Plane Surface(30001) = {30001}; // Top unstructured region
Plane Surface(30002) = {30002}; // Bottom unstructured region
Plane Surface(30003) = {30003}; // South unstructured region
Plane Surface(30004) = {30004}; // East unstructured region
Plane Surface(30005) = {30005}; // North unstructured region
Plane Surface(30006) = {30006}; // West unstructured region
=======
//Plane Surface(30001) = {30001}; // Top unstructured region
//Plane Surface(30002) = {30002}; // Bottom unstructured region
//Plane Surface(30003) = {30003}; // South unstructured region
//Plane Surface(30004) = {30004}; // East unstructured region
//Plane Surface(30005) = {30005}; // North unstructured region
//Plane Surface(30006) = {30006}; // West unstructured region
>>>>>>> iRMS_Biot

//Surface Loop(40000) = {3001,3002,3003,3004,3005,3006,30001,30002,30003,30004,30005,30006};
//Volume(7) = {40000} ;

////////////////////////////////////////////////////////////////////////////
// Converting wellbore regions to hexahedron mesh ! very experimental in 3D
////////////////////////////////////////////////////////////////////////////

If (hexahedronsQ == 1)
Transfinite Surface "*";
Transfinite Volume "*";
Recombine Surface "*";
Recombine Volume "*";
EndIf

////////////////////////////////////////////////////////////////////////////
// Mark physical entities
////////////////////////////////////////////////////////////////////////////

// Tagging boundary conditions for prodution wells
Physical Surface("well_lids") = well_lids[];

Physical Surface("producers") = well_p_bores[];
//Physical Volume("producers_region") = well_p_v_regions[];

Physical Surface("injectors") = well_i_bores[];
//Physical Volume("injectors_region") = well_i_v_regions[];

Physical Volume("Reservoir") = {6,well_p_v_regions[],well_i_v_regions[]};


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

Else

////////////////////////////////////////////////////////////////////////////
// Reservoir rocks
////////////////////////////////////////////////////////////////////////////


If (xzQ == 1)

// reservoir box dimensions on y-plane
y_length = 0.0;

Else

// reservoir box dimensions on y-plane
z_length = 0.0;

EndIf



Point(1001) = {-x_length/2.0, -y_length/2.0, -z_length/2.0, cl4};
Point(1002) = { x_length/2.0, -y_length/2.0, -z_length/2.0, cl4};
Point(1003) = { x_length/2.0, +y_length/2.0, +z_length/2.0, cl4};
Point(1004) = {-x_length/2.0, +y_length/2.0, +z_length/2.0, cl4};

Line(2001) = {1001,1002};
Line(2002) = {1002,1003};
Line(2003) = {1003,1004};
Line(2004) = {1004,1001};


Line Loop(3001) = {well_p_regions[],well_i_regions[],2001, 2002, 2003, 2004}; // Bottom
Plane Surface(3001) = {3001}; // unstructured region

////////////////////////////////////////////////////////////////////////////
// Side-burden rocks
////////////////////////////////////////////////////////////////////////////

If (xzQ == 1)

// side-burden box dimensions  on y-plane
sb_y_length = 0.0;

Else

// side-burden box dimensions  on z-plane
sb_z_length = 0.0;

EndIf



Point(10001) = {-sb_x_length/2.0, -sb_y_length/2.0, -sb_z_length/2.0, cl5};
Point(10002) = { sb_x_length/2.0, -sb_y_length/2.0, -sb_z_length/2.0, cl5};
Point(10003) = { sb_x_length/2.0,  sb_y_length/2.0, +sb_z_length/2.0, cl5};
Point(10004) = {-sb_x_length/2.0,  sb_y_length/2.0, +sb_z_length/2.0, cl5};


//Line(20001) = {10001,10002};
//Line(20002) = {10002,10003};
//Line(20003) = {10003,10004};
//Line(20004) = {10004,10001};

Line Loop(30001) = {2001, 2002, 2003, 2004, 20001, 20002, 20003, 20004};

//Plane Surface(30001) = {30001}; // unstructured region

////////////////////////////////////////////////////////////////////////////
// Converting reservoir and side burden regions to hexahedron mesh ! very experimental in 3D
////////////////////////////////////////////////////////////////////////////


If (hexahedronsOutQ == 1)
//Transfinite Line {2001,2002,2003,2004} = n_bc_res;
//Transfinite Line {20001,20002,20003,20004} = n_bc_sb;
<<<<<<< HEAD
Transfinite Surface {3001} = n_bc_res;
=======
//Transfinite Surface {3001} = n_bc_res;
Transfinite Surface {3001};
>>>>>>> iRMS_Biot
Recombine Surface "*";
Recombine Volume "*";
EndIf

////////////////////////////////////////////////////////////////////////////
// Mark physical entities
////////////////////////////////////////////////////////////////////////////

// Tagging boundary conditions for prodution wells
Physical Line("well_lids") = well_lids[];

Physical Line("producers") = well_p_bores[];

Physical Line("injectors") = well_i_bores[];

Physical Surface("Reservoir") = {3001,well_p_v_regions[],well_i_v_regions[]};
Physical Line("Reservoir_bottom") = {2001};
Physical Line("Reservoir_West") = {2002};
Physical Line("Reservoir_top") = {2003};
Physical Line("Reservoir_East") = {2004};



Physical Surface("side_burden") = {30001};
Physical Line("side_burden_bottom") = {20001};
Physical Line("side_burden_West") = {20002};
Physical Line("side_burden_top") = {20003};
Physical Line("side_burden_East") = {20004};

EndIf
