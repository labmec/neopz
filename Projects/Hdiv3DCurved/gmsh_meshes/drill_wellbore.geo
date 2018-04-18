
// ---- Gmsh Macro ----
// ---- Vertical wellbore wih hybrid 3D mesh  ----
// Created 01/02/2018 by Omar Durán
// Labmec, University of Campinas
// --------------------------------------------

Macro MakeVerticalWellbore

l = 1.0;
r1 = inner_r/Sqrt(2.0);
r2 = outer_r/Sqrt(2.0);

p1 = newp; Point(p1) = {0,0,-h/2,l};

// interior circle
p2 = newp; Point(p2) = {r1,r1,-h/2,l};
p3 = newp; Point(p3) = {-r1,r1,-h/2,l};
p4 = newp; Point(p4) = {-r1,-r1,-h/2,l};
p5 = newp; Point(p5) = {r1,-r1,-h/2,l};
l1 = newl; Circle(l1) = {p3,p1,p2};
l2 = newl; Circle(l2) = {p2,p1,p5};
l3 = newl; Circle(l3) = {p5,p1,p4};
l4 = newl; Circle(l4) = {p4,p1,p3};

// exterior circle
pe2 = newp; Point(pe2) = {r2,r2,-h/2,l};
pe3 = newp; Point(pe3) = {-r2,r2,-h/2,l};
pe4 = newp; Point(pe4) = {-r2,-r2,-h/2,l};
pe5 = newp; Point(pe5) = {r2,-r2,-h/2,l};
le1 = newl; Circle(le1) = {pe3,p1,pe2};
le2 = newl; Circle(le2) = {pe2,p1,pe5};
le3 = newl; Circle(le3) = {pe5,p1,pe4};
le4 = newl; Circle(le4) = {pe4,p1,pe3};

pm2 = newp; Point(pm2) = {s*r1,s*r1,-h/2,l};
pm3 = newp; Point(pm3) = {-s*r1,s*r1,-h/2,l};
pm4 = newp; Point(pm4) = {-s*r1,-s*r1,-h/2,l};
pm5 = newp; Point(pm5) = {s*r1,-s*r1,-h/2,l};

lr1  = newl; Line(lr1)  = {pm2,pe2};
lr2  = newl; Line(lr2)  = {pm3,pe3};
lr3  = newl; Line(lr3)  = {pm4,pe4};
lr4  = newl; Line(lr4)  = {pm5,pe5};

lbox1  = newl; Line(lbox1)  = {pm2,pm3};
lbox2  = newl; Line(lbox2)  = {pm3,pm4};
lbox3  = newl; Line(lbox3)  = {pm4,pm5};
lbox4  = newl; Line(lbox4)  = {pm5,pm2};

// Hard coded
Line Loop(24) = {13, 10, 5, -9};
Plane Surface(25) = {24};

Line Loop(25) = {14, 11, 8, -10};
Plane Surface(26) = {25};

Line Loop(26) = {15, 12, 7, -11};
Plane Surface(27) = {26};

Line Loop(27) = {16, 9, 6, -12};
Plane Surface(28) = {27};

Line Loop(28) = {15, 16, 13, 14};
Line Loop(29) = {4, 1, 2, 3}; Plane Surface(29) = {28, 29};


If(mesh_type == 1 || mesh_type == 2) 
Extrude {0, 0, h} {
  Surface{27}; 
  Surface{28}; 
  Surface{25}; 
  Surface{26}; 
  Surface{29};
  Layers{n_vertical};
  Recombine; 
}
Else
Extrude {0, 0, h} {
  Surface{27}; 
  Surface{28}; 
  Surface{25}; 
  Surface{26}; 
  Surface{29}; 
}
EndIf

// Line Groups
box_h_edges[] = {13,14,15,16,31,53,75,97};
box_v_edges[] = {36,37,59,81};
wellbore_v_edges[] = {144,145,149,153};
wellbore_h_edges[] = {1,2,3,4,123,124,125,126};
outer_h_edges[] = {5,6,7,8,33,99,55,77};
outer_v_edges[] = {41,45,63,85};
radial_edges[]= {9,10,11,12,32,-34,54,76};

// Surface Groups
mid_surf_box[]= {38, 104, 82, 60};
mid_surf_rad[]= {86, 50, 42, 64};
mid_surf[] = {mid_surf_rad[],mid_surf_box[]};
top_bottom_reservoir_bc[]={27, 51, 28, 73, 25, 95, 26, 117};
top_bottom_wellbore_region_bc[]={29,159};
wellbore_bc[]= {146, 150, 154, 158};
reservoir_bc[]= {46, 68, 90, 112};

// Volume Groups
reservoir[]={1,2,3,4};
wellbore_region[]={5};


// Mesh types

If(mesh_type == 0)

// Meshing directives for lines
Transfinite Line {radial_edges[]} = n_radial Using Progression radial_progression; // Radial control
Transfinite Line {outer_h_edges[],box_h_edges[]} = n_azimuthal; // Azimuthal control
Transfinite Line {box_v_edges[],wellbore_v_edges[]} = n_vertical; // Vertical control

EndIf


If(mesh_type == 1)

// Meshing directives for lines
Transfinite Line {radial_edges[]} = n_radial Using Progression radial_progression; // Radial control
Transfinite Line {outer_h_edges[],box_h_edges[],wellbore_h_edges[]} = n_azimuthal; // Azimuthal control
Transfinite Line {outer_v_edges[],box_v_edges[],wellbore_v_edges[]} = n_vertical; // Vertical control

// Meshing directives for surfaces
Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";
Recombine Volume "*";

// 3D mesh algorithm (1=Delaunay, 2=New Delaunay, 4=Frontal, 5=Frontal Delaunay, 6=Frontal Hex, 7=MMG3D, 9=R-tree)
Mesh.Algorithm3D = 4;

EndIf

If(mesh_type == 2)

// Meshing directives for lines
Transfinite Line {radial_edges[]} = n_radial Using Progression radial_progression; // Radial control
Transfinite Line {outer_h_edges[],box_h_edges[],wellbore_h_edges[]} = n_azimuthal; // Azimuthal control
Transfinite Line {outer_v_edges[],box_v_edges[],wellbore_v_edges[]} = n_vertical; // Vertical control

// Meshing directives for surfaces
Transfinite Surface{mid_surf[],reservoir_bc[],wellbore_bc[]};
Recombine Surface{mid_surf[],reservoir_bc[],wellbore_bc[]};


EndIf


If(mesh_type == 3)

// Meshing directives for lines
Transfinite Line {radial_edges[]} = n_radial Using Progression radial_progression; // Radial control
Transfinite Line {outer_h_edges[],box_h_edges[]} = n_azimuthal; // Azimuthal control
Transfinite Line {outer_v_edges[],box_v_edges[],wellbore_v_edges[]} = n_vertical; // Vertical control

// Meshing directives for surfaces
Transfinite Surface{mid_surf[],reservoir_bc[],top_bottom_reservoir_bc[]};
Recombine Surface{mid_surf_rad[],top_bottom_reservoir_bc[],reservoir_bc[]};


// Meshing directives for volumes
Transfinite Volume{reservoir[]}; // Regular partition for the reservoir region
TransfQuadTri {1,5}; // Directive to force the pyramids between volumes : 1 (quads) and 5 (tri)
TransfQuadTri {2,5}; // Directive to force the pyramids between volumes : 2 (quads) and 5 (tri)
TransfQuadTri {3,5}; // Directive to force the pyramids between volumes : 3 (quads) and 5 (tri)
TransfQuadTri {4,5}; // Directive to force the pyramids between volumes : 4 (quads) and 5 (tri)

// 3D mesh algorithm (1=Delaunay, 2=New Delaunay, 4=Frontal, 5=Frontal Delaunay, 6=Frontal Hex, 7=MMG3D, 9=R-tree)
Mesh.Algorithm3D = 4;

EndIf

// Tagging
Physical Volume("reservoir") = {reservoir[],wellbore_region[]};
//Physical Volume("wellbore") = {wellbore_region[]};
//Physical Surface("outer_bc") = {reservoir_bc[]};
Physical Surface("inner_bc") = {wellbore_bc[]};
Physical Surface("non_flux_bc") = {reservoir_bc[],top_bottom_wellbore_region_bc[], top_bottom_reservoir_bc[]};


Another_entities = 0;
If(Another_entities)
Physical Surface("non_flux_res") = {top_bottom_reservoir_bc[]};
Physical Surface("non_flux_well") = {top_bottom_wellbore_region_bc[]};
Physical Surface("mid_box") = {mid_surf_box[]};
Physical Surface("mid_rad") = {mid_surf_rad[]};
Physical Line("box_h_edges") = {box_h_edges[]};
Physical Line("box_v_edges") = {box_v_edges[]};
Physical Line("wellbore_v_edges") = {wellbore_v_edges[]};
Physical Line("wellbore_h_edges") = {wellbore_h_edges[]};
Physical Line("outer_h_edges") = {outer_h_edges[]};
Physical Line("outer_v_edges") = {outer_v_edges[]};
Physical Line("radial_edges") = {radial_edges[]};
EndIf

// optimize the mesh
Mesh.Optimize = 1;
Mesh  3;


If(mesh_type == 0)
Save "./vertical_wellbore_Te.msh";
EndIf

If(mesh_type == 1)
Save "./vertical_wellbore_He.msh";
EndIf

If(mesh_type == 2)
Save "./vertical_wellbore_Pe.msh";
EndIf

If(mesh_type == 3)
Save "./vertical_wellbore_hybrid.msh";
EndIf

Return
