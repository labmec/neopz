
Mesh.CharacteristicLengthExtendFromBoundary = 1;

DFN_reservoirQ = 1;

// ---------------- Parameters ----------------
n_bc           = 4;   // number of bc elements
n_azimuthal    = 4;   // number of azimuthal elements
r_i			 = 0.1;   // Radius of the circle

width_st   = 1.0; // Width of the strucured box
height_st  = 1.0; // Height of the structured box


// ---------------- Inner circle ----------------
p0 = newp; Point(p0) = {0, 0, 0, 1e22};
p1 = newp; Point(p1) = {r_i, 0, 0, 1e22};
p2 = newp; Point(p2) = {-r_i, 0, 0, 1e22};
p3 = newp; Point(p3) = {0, -r_i, 0, 1e22};
p4 = newp; Point(p4) = {0, r_i, 0, 1e22};

i_l1 = newl; Circle(i_l1) = {p4, p0, p1};
i_l2 = newl; Circle(i_l2) = {p1, p0, p3};
i_l3 = newl; Circle(i_l3) = {p3, p0, p2};
i_l4 = newl; Circle(i_l4) = {p2, p0, p4};

wellbore_bc[] = {i_l1, i_l2, i_l3, i_l4};
ll1 = newll; Line Loop(ll1) = wellbore_bc[];
Transfinite Line{i_l1, i_l2, i_l3, i_l4} = n_azimuthal Using Progression 1;


// ---------------- Box part -----------------
p5 = newp; Point(p5) = {-width_st/2,-height_st/2, 0};
p6 = newp; Point(p6) = {width_st/2,-height_st/2, 0};
p7 = newp; Point(p7) = {width_st/2,height_st/2, 0};
p8 = newp; Point(p8) = {-width_st/2, height_st/2, 0};

///////////////////////////////////////////////////////////
// Fixed displacement points

p9 = newp; Point(p9) = {-width_st/2,0, 0};
p10 = newp; Point(p10) = {width_st/2,0, 0};
p11 = newp; Point(p11) = {0,height_st/2, 0};
p12 = newp; Point(p12) = {0,-height_st/2, 0};

fixed_uy[] = {p9,p10};
fixed_ux[] = {p11,p12};

o_l11 = newl; Line(o_l11) = {p5,p12};
o_l12 = newl; Line(o_l12) = {p12,p6};
o_l21 = newl; Line(o_l21) = {p6,p10};
o_l22 = newl; Line(o_l22) = {p10,p7};
o_l31 = newl; Line(o_l31) = {p7,p11};
o_l32 = newl; Line(o_l32) = {p11,p8};
o_l41 = newl; Line(o_l41) = {p8,p9};
o_l42 = newl; Line(o_l42) = {p9,p5};

S_bc[] = {o_l11,o_l12};
E_bc[] = {o_l21,o_l22};
N_bc[] = {o_l31,o_l32};
W_bc[] = {o_l41,o_l42};
ll2 = newll; Line Loop(ll2) = {S_bc[],E_bc[],N_bc[],W_bc[]};
Transfinite Line{S_bc[],E_bc[],N_bc[],W_bc[]} = n_bc Using Progression 1;
s1 = news; Plane Surface(s1) = {ll2,ll1};




If (DFN_reservoirQ == 1)

///////////////////////////////////////////////////////////
// Inserting fractures

// Create first Fracture
f1p0 = p1;
f1p1 = newp; Point(f1p1) = {width_st/4,0,0};
f1 = newl; Line(f1) = {f1p0,f1p1};

// Create Second Fracture
f2p0 = p2;
f2p1 = newp; Point(f2p1) = {-width_st/4,0,0};
f2 = newl; Line(f2) = {f2p0,f2p1};

// Create far Fractures
f3p0 = newp; Point(f3p0) = {width_st/4,width_st/8,0};
f3p1 = newp; Point(f3p1) = {width_st/8,width_st/4,0};
f3 = newl; Line(f3) = {f3p0,f3p1};

f4p0 = newp; Point(f4p0) = {-width_st/4,width_st/8,0};
f4p1 = newp; Point(f4p1) = {-width_st/8,width_st/4,0};
f4 = newl; Line(f4) = {f4p0,f4p1};

f5p0 = newp; Point(f5p0) = {-width_st/4,-width_st/8,0};
f5p1 = newp; Point(f5p1) = {-width_st/8,-width_st/4,0};
f5 = newl; Line(f5) = {f5p0,f5p1};

f6p0 = newp; Point(f6p0) = {width_st/4,-width_st/8,0};
f6p1 = newp; Point(f6p1) = {width_st/8,-width_st/4,0};
f6 = newl; Line(f6) = {f6p0,f6p1};

//fractures[] = {f1,f2,f3,f4,f5,f6};
fractures[] = {f3,f4,f5,f6};
fracture_tips[] = {f2p1,f1p1,f3p0,f3p1,f4p0,f4p1,f5p0,f5p1,f6p0,f6p1};
fracture_tips_bc[] = {f1p0,f2p0};

///////////////////////////////////////////////////////////
// Embedded Fractures
Line{fractures[]} In Surface {s1};

///////////////////////////////////////////////////////////
// Fracture boundaries
Point{fracture_tips[],fracture_tips_bc[]} In Surface {s1};

///////////////////////////////////////////////////////////
// Refinemente towards fractures

lc = 0.5;
r  = 0.25;

Field[1] = Attractor;
Field[1].NNodesByEdge = 5; // #attractors on the edges
Field[1].NodesList = {fracture_tips_bc[],fracture_tips[]};
Field[1].EdgesList = {fractures[]};

// Threshold field defined on the attractors
Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = lc/10; // char length inside DistMin
Field[2].LcMax = lc; // char length outside DistMax
Field[2].DistMin = 0.25*r;
Field[2].DistMax = 0.5*r;

// Define minimum of threshold and function field
Field[3] = Min;
Field[3].FieldsList = {2};
Background Field = 3;
///////////////////////////////////////////////////////////

EndIf

///////////////////////////////////////////////////////////
// Tagging surface and boundary domains
Physical Surface("Omega") = {s1};
Physical Line("Wellbore") = {wellbore_bc[]};
Physical Line("S") = {S_bc[]};
Physical Line("E") = {E_bc[]};
Physical Line("N") = {N_bc[]};
Physical Line("W") = {W_bc[]};

Physical Point("Fixed_ux") = {fixed_ux[]};
Physical Point("Fixed_uy") = {fixed_uy[]};

If (DFN_reservoirQ == 1)
Physical Line("Fractures") = {fractures[]};
Physical Point("Fractures_fixed") = {fracture_tips[]};
Physical Point("Fractures_wellbore_bc") = {fracture_tips_bc[]};
EndIf



