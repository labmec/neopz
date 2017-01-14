// ---- Cylinder wellbore Region Gmsh scritp ----
// 2D_cylinder_tutorial.geo
// Creates a mesh with an inner structured-quad region and 
// an outer unstructured tri region
//
// Created 10/01/2016 by Omar Duran
// Labmec, State University of Campinas
// --------------------------------------------

// Settings
ExpertMode = 1;
Mesh.Algorithm3D = 1 ;


// general structures



well_bores = {};
well_regions = {};

// Gmsh allows variables; these will be used to set desired
// element sizes at various Points
cl1 = 1;
cl2 = 0.1;
cl3 = 10.0;
cl4 = 100.0;
cl5 = 250.0;


// reservoir box dimensions
x_length = 1000.0;
y_length = 1000.0;
z_length = 100.0;

// side-burden box dimensions
sb_x_length = 5000.0;
sb_y_length = 5000.0;
sb_z_length = 1000.0;

// mesh controls on wellbore region
alpha = 1.5;
n_radial = 10;
n_azimuthal = 4;
n_axial = 2; 


// new well data
well_index = 1;

// well location
wx = 50.0;
wy = 0.0;
wz = 0.0;

// well and wellbore region dimensions
radius = 0.1;
outer = 20;
length = 10.0;
angle = 0.5;
beta = 0.5;


// Scripting

xdir = Cos(beta);
ydir = Sin(beta);

// Circle & Surrounding structured-quad region

p1=newp; Point(p1) = {0, 0, 0, cl2};

// well bore points
p2=newp; Point(p2) = {0,  radius, 0, cl2};
p3=newp; Point(p3) = {radius,  0, 0, cl2};
p4=newp; Point(p4) = {0, -radius, 0, cl2};
p5=newp; Point(p5) = {-radius, 0, 0, cl2};

// well bore region points
p6=newp; Point(p6) = {0,  outer, 0, cl3};
p7=newp; Point(p7) = {outer,  0, 0, cl3};
p8=newp; Point(p8) = {0, -outer, 0, cl3};
p9=newp; Point(p9) = {-outer, 0, 0, cl3};

// Apply translation
translate_s[] = Translate {wx, wy, wz} { 
Point {p1,p2,p3,p4,p5,p6,p7,p8,p9};
};

// Apply rotation
rotate_s[] = Rotate { { xdir, ydir,0}, {0,0,0}, angle } {
Point {p1,p2,p3,p4,p5,p6,p7,p8,p9};
};

c1[] = Point{p1};
c2[] = Point{p2};
c3[] = Point{p3};

u[] = {}; v[] = {};
For i In {0:3-1}
	u[i] =  c2[i] - c1[i];
	v[i] =  c3[i] - c1[i];
EndFor

a[] = {-u[2]*v[1]+u[1]*v[2],u[2]*v[0]-u[0]*v[2],-u[1]*v[0]+u[0]*v[1]};
norm = Sqrt( a[0]*a[0] + a[1]*a[1] + a[2]*a[2] ) / length;

l1=newl; Circle(l1) = {p2, p1, p3};
l2=newl; Circle(l2) = {p3, p1, p4};
l3=newl; Circle(l3) = {p4, p1, p5};
l4=newl; Circle(l4) = {p5, p1, p2};

l5=newl; Circle(l5) = {p6, p1, p7};
l6=newl; Circle(l6) = {p7, p1, p8};
l7=newl; Circle(l7) = {p8, p1, p9};
l8=newl; Circle(l8) = {p9, p1, p6};

l9=newl;  Line(l9)  = {p2, p6};
l10=newl; Line(l10) = {p3, p7};
l11=newl; Line(l11) = {p4, p8};
l12=newl; Line(l12) = {p5, p9};

Printf ("Mehser:: Drawing well");

Transfinite Line {l1,l2,l3,l4,l5,l6,l7,l8} = n_azimuthal;
Transfinite Line {l9,l10,l11,l12} = n_radial Using Progression alpha;


// Each region which to be independently meshed must have a line loop
// Regions which will be meshed with Transfinite Surface must have 4 lines
// and be labeled in CCW order, with the correct orientation of each edge

Line Loop(1) = {l1, l10, -l5, -l9}; // ++ side of quad region - note ordering
Line Loop(2) = {l2, l11, -l6, -l10}; // +- side of quad region - note ordering
Line Loop(3) = {l3, l12, -l7, -l11}; // -- side of quad region - note ordering
Line Loop(4) = {l4, l9, -l8, -l12}; // -+ side of quad region - note ordering

s1 = news; Plane Surface(s1) = {1}; // ++ side of quad
s2 = news; Plane Surface(s2) = {2}; // +- side of quad
s3 = news; Plane Surface(s3) = {3}; // -- side of quad
s4 = news; Plane Surface(s4) = {4}; // -+ side of quad


out[] = Extrude {a[0]/norm, a[1]/norm, a[2]/norm} { Surface{s1,s2,s3,s4}; Layers{n_axial}; };
//out[] = Extrude {a[0]/norm, a[1]/norm, a[2]/norm} { Surface{s1,s2,s3,s4}; Layers{n_axial}; Recombine; QuadTriNoNewVerts Recomblaterals;};


//N = #out[];
//Printf ("size = %g", N);
//For i In {0:N-1}
//	Printf ("out[%g] = %g", i, out[i]);
//EndFor

s1_ext = out[0];
s2_ext = out[6];
s3_ext = out[12];
s4_ext = out[18];

s1_lat = out[4];
s2_lat = out[10];
s3_lat = out[16];
s4_lat = out[22];

well_1 = out[2];
well_2 = out[8];
well_3 = out[14];
well_4 = out[20];

wing_1 = out[3];
wing_2 = out[5];
wing_3 = out[9];
wing_4 = out[15];


////////////////////////////////////////////////////////////////////////////
/// Well box

// well closed regions

Line Loop(5) = {l1, l2, l3, l4}; // clossing the well
s5 = news; Plane Surface(s5) = {5}; // clossing the well

Line Loop(6) = {l1 + 17, l2 + 38, l3 + 59, l4 + 80}; // clossing the well
s6 = news; Plane Surface(s6) = {6}; // clossing the well

well_region[] = {s1,s2,s3,s4,s5,s6,s1_ext,s2_ext,s3_ext,s4_ext,s1_lat,s2_lat,s3_lat,s4_lat};
well_bore[] = {well_1,well_2,well_3,well_4};


N = #well_region[];
M = #well_regions[];
For i In {0:N-1}
	well_regions[i + M] = well_region[i];
EndFor

//M = #well_regions[];
//For i In {0:N-1}
//	Printf ("well_regions[%g] = %g", i, well_regions[i]);
//EndFor


N = #well_bore[];
M = #well_bores[];
For i In {0:N-1}
	well_bores[i + M] = well_bore[i];
EndFor

//M = #well_bores[];
//For i In {0:N-1}
//	Printf ("well_bores[%g] = %g", i, well_bores[i]);
//EndFor


Surface Loop(1) = {s1,s2,s3,s4,s1_ext,s2_ext,s3_ext,s4_ext,well_1,well_2,well_3,well_4,s1_lat,s2_lat,s3_lat,s4_lat}; 

// wellbore region volume
vregion = newv; Volume(newv) = {1};

// Tagging boundary conditions for wells
Physical Surface("well_p") = well_bores[];
Physical Surface("well_closed") = {s5,s6};
Physical Volume("well_p_region") = {vregion};

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


//Recombine Surface {3001,3002,3003,3004,3005,3006};

// Build the well box volume

reservoir_region[] = {well_regions[],3001,3002,3003,3004,3005,3006};
Surface Loop(4000) = reservoir_region[];
Volume(6) = {4000} ;

//Transfinite Surface "*";
//Recombine Surface "*";
//Recombine Volume "*";

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
Physical Surface("side_burden_laterals") = {30003,30004,30005,30006};


