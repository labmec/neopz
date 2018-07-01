<<<<<<< HEAD
// ---- Cylinder wellbore Region Gmsh scritp ----
=======
// ---- Drill cylindrical wellbore Region Gmsh scritp ----
>>>>>>> iRMS_Biot
// 2D_cylinder_tutorial.geo
// Creates a mesh with an inner structured-quad region and 
// an outer unstructured tri region
//
// Created 10/01/2016 by Omar Duran
// Labmec, State University of Campinas
// --------------------------------------------


Macro DrillProducer

xdir = Cos(beta);
ydir = Sin(beta);

// Circle & Surrounding structured-quad region

<<<<<<< HEAD
=======
If (dimension == 3)

>>>>>>> iRMS_Biot
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

// Apply rotation
rotate_s[] = Rotate { { xdir, ydir,0}, {0,0,0}, angle } {
Point {p1,p2,p3,p4,p5,p6,p7,p8,p9};
};

<<<<<<< HEAD
// Apply translation
=======
Else

p1=newp; Point(p1) = {0, 0, 0, cl2};

If (xzQ == 1)

// well bore points
p2=newp; Point(p2) = {0,  0, radius, cl2};
p3=newp; Point(p3) = {radius,  0, 0, cl2};
p4=newp; Point(p4) = {0, 0, -radius, cl2};
p5=newp; Point(p5) = {-radius, 0, 0, cl2};

// well bore region points
p6=newp; Point(p6) = {0,  0, outer, cl3};
p7=newp; Point(p7) = {outer,  0, 0, cl3};
p8=newp; Point(p8) = {0, 0, -outer, cl3};
p9=newp; Point(p9) = {-outer, 0, 0, cl3};

Else

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

EndIf


EndIf

// Apply translation
If (dimension == 3)

>>>>>>> iRMS_Biot
translate_s[] = Translate {wx, wy, wz} { 
Point {p1,p2,p3,p4,p5,p6,p7,p8,p9};
};

<<<<<<< HEAD
=======
Else

If (xzQ == 1)

translate_s[] = Translate {wx, 0, wz} { 
Point {p1,p2,p3,p4,p5,p6,p7,p8,p9};
};

Else

translate_s[] = Translate {wx, wy, 0} { 
Point {p1,p2,p3,p4,p5,p6,p7,p8,p9};
};

EndIf

EndIf


>>>>>>> iRMS_Biot
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

Transfinite Line {l1,l2,l3,l4,l5,l6,l7,l8} = n_azimuthal;
Transfinite Line {l9,l10,l11,l12} = n_radial Using Progression alpha;


// Each region which to be independently meshed must have a line loop
// Regions which will be meshed with Transfinite Surface must have 4 lines
// and be labeled in CCW order, with the correct orientation of each edge

ll1 = newll; Line Loop(ll1) = {l1, l10, -l5, -l9}; // ++ side of quad region - note ordering
ll2 = newll; Line Loop(ll2) = {l2, l11, -l6, -l10}; // +- side of quad region - note ordering
ll3 = newll; Line Loop(ll3) = {l3, l12, -l7, -l11}; // -- side of quad region - note ordering
ll4 = newll; Line Loop(ll4) = {l4, l9, -l8, -l12}; // -+ side of quad region - note ordering

s1 = news; Plane Surface(s1) = {ll1}; // ++ side of quad
s2 = news; Plane Surface(s2) = {ll2}; // +- side of quad
s3 = news; Plane Surface(s3) = {ll3}; // -- side of quad
s4 = news; Plane Surface(s4) = {ll4}; // -+ side of quad

<<<<<<< HEAD
out[] = {};
out[] = Extrude {a[0]/norm, a[1]/norm, a[2]/norm} { Surface{s1,s2,s3,s4}; Layers{n_axial}; };
//out[] = Extrude {a[0]/norm, a[1]/norm, a[2]/norm} { Surface{s1,s2,s3,s4}; Layers{n_axial}; Recombine; QuadTriNoNewVerts Recomblaterals;};
=======
If (dimension == 3)

out[] = {};
<<<<<<< HEAD
If (hexahedronsWQ == 0)
=======
If (hexahedronsQ == 0)
>>>>>>> iRMS_Biot
out[] = Extrude {a[0]/norm, a[1]/norm, a[2]/norm} { Surface{s1,s2,s3,s4}; Layers{n_axial}; };
Else
out[] = Extrude {a[0]/norm, a[1]/norm, a[2]/norm} { Surface{s1,s2,s3,s4}; Layers{n_axial}; Recombine;};
EndIf
>>>>>>> iRMS_Biot


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

ll5 = newll; Line Loop(ll5) = {l1, l2, l3, l4}; // clossing the well
s5 = news; Plane Surface(s5) = {ll5}; // clossing the well

ll6 = newll; Line Loop(ll6) = {l1 + 21, l2 + 42, l3 + 63, l4 + 84}; // clossing the well
s6 = news; Plane Surface(s6) = {ll6}; // clossing the well

<<<<<<< HEAD
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


sl1 = newsl; Surface Loop(sl1) = {s1,s2,s3,s4,s1_ext,s2_ext,s3_ext,s4_ext,well_1,well_2,well_3,well_4,s1_lat,s2_lat,s3_lat,s4_lat}; 

// wellbore region volume
vregion = newv; Volume(newv) = {sl1};

m = 10;

// Tagging boundary conditions for wells
//well = 1 + (well_index-1)*m; Physical Surface(well) = well_bore[];
//wellVol = 2 + (well_index-1)*m; Physical Volume(wellVol) = {vregion};
//wellEnds = 3 + (well_index-1)*m; Physical Surface(wellEnds) = {s5,s6};

Printf ("Mesher:: Drilled producer well = %g", well_index);
=======
well_lid[] = {s5,s6};
well_region[] = {s1,s2,s3,s4,s5,s6,s1_ext,s2_ext,s3_ext,s4_ext,s1_lat,s2_lat,s3_lat,s4_lat};
well_bore[] = {well_1,well_2,well_3,well_4};

Else

well_lid[] = {};
well_region[] = {l5,l6,l7,l8};
well_bore[] = {l1,l2,l3,l4};

EndIf



N = #well_lid[];
M = #well_lids[];
For i In {0:N-1}
	well_lids[i + M] = well_lid[i];
EndFor

N = #well_region[];
M = #well_p_regions[];
For i In {0:N-1}
	well_p_regions[i + M] = well_region[i];
EndFor

N = #well_bore[];
M = #well_p_bores[];
For i In {0:N-1}
	well_p_bores[i + M] = well_bore[i];
EndFor

If (dimension == 3)

// wellbore region volume
sl1 = newsl; Surface Loop(sl1) = {s1,s2,s3,s4,s1_ext,s2_ext,s3_ext,s4_ext,well_1,well_2,well_3,well_4,s1_lat,s2_lat,s3_lat,s4_lat}; 
vregion = newv; Volume(vregion) = {sl1};

j=(well_index-1)*114;
well_p_v_region[] = {1+j,2+j,3+j,4+j};

<<<<<<< HEAD
If (hexahedronsWQ == 1)
=======
If (hexahedronsQ == 1)
>>>>>>> iRMS_Biot
Transfinite Volume {well_p_v_region[]};
Recombine Volume {well_p_v_region[]};
EndIf


Else

well_p_v_region[] = {s1,s2,s3,s4};

<<<<<<< HEAD
If (hexahedronsWQ == 1)
=======
If (hexahedronsQ == 1)
>>>>>>> iRMS_Biot
Transfinite Surface {well_p_v_region[]};
Recombine Surface {well_p_v_region[]};
EndIf

EndIf


N = #well_p_v_region[];
M = #well_p_v_regions[];
For i In {0:N-1}
	well_p_v_regions[i + M] = well_p_v_region[i];
EndFor

Printf ("Mesher:: Drilled producer well = %g, at position {%g,%g,%g}", well_index, wx,wy,wz);
>>>>>>> iRMS_Biot

Return
