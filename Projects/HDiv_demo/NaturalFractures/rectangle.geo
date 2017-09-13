/*********************************************************************
 *
 *  Gmsh circle
 *	Define a macro that draw a rectangle with three embeded fractures
 *
 *********************************************************************/



Macro RectangleDomain

lx= 10.0;
ly= 10.0;

r=Sqrt(lx*lx+ly*ly);
lc=0.2*r;


p1 = newp; Point(p1) = {lx,ly,0,lc};
p2 = newp; Point(p2) = {-lx,ly,0,lc};
p3 = newp; Point(p3) = {-lx,-ly,0,lc};
p4 = newp; Point(p4) = {lx,-ly,0,lc};

l0 = newl; Line(l0) = {p1,p2};
l1 = newl; Line(l1) = {p2,p3};
l2 = newl; Line(l2) = {p3,p4};
l3 = newl; Line(l3) = {p4,p1};


ll0 = newll; Line Loop(ll0) = {l0,l1,l2,l3};
s0 = news; Plane Surface(s0) = {ll0};


// Create first Fracture
f1p0 = newp; Point(f1p0) = {r*0.25,-r*0.65,0,lc};
f1p1 = newp; Point(f1p1) = {-r*0.5,-r*0.35,0,lc};
f1 = newl; Line(f1) = {f1p0,f1p1};

// Create Second Fracture
f2p0 = newp; Point(f2p0) = {r*0.1,r*0.3,0,lc};
f2p1 = newp; Point(f2p1) = {r*0.6,r*0.5,0,lc};
f2 = newl; Line(f2) = {f2p0,f2p1};

// Create Third Fracture
f3p0 = newp; Point(f3p0) = {r*0.05,-r*0.2,0,lc};
f3p1 = newp; Point(f3p1) = {-r*0.5,r*0.4,0,lc};
f3 = newl; Line(f3) = {f3p0,f3p1};

// Embedded Fractures
Line{f1,f2,f3} In Surface {s0};

///////////////////////////////////////////////////////////
// Refinemente towards fractures
Field[1] = Attractor;
Field[1].NNodesByEdge = 10; // #attractors on the edges
Field[1].NodesList = {f1p0,f1p1,f2p0,f2p1,f3p0,f3p1};
Field[1].EdgesList = {f1,f2,f3};

// Threshold field defined on the attractors
Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = lc/10; // char length inside DistMin
Field[2].LcMax = lc; // char length outside DistMax
Field[2].DistMin = 0.1*r;
Field[2].DistMax = 0.2*r;

// Define minimum of threshold and function field
Field[3] = Min;
Field[3].FieldsList = {2};
Background Field = 3;
///////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////
// Tagging surface and boundary domains
Physical Surface("Omega") = {s0};
Physical Line("Gamma_D_inlet") = {l1};
Physical Line("Gamma_D_outlet") = {l3};
Physical Line("Gamma_N") = {l0,l2};
Physical Line("Fractures") = {f1,f2,f3};

///////////////////////////////////////////////////////////
// Coloring surface and boundary domains
Color Green{ Surface{ s0 }; }
Color Blue{ Line{ l0,l1,l2,l3 }; }
Color Red{ Line{ f1,f2,f3 }; }

If (QuadrilateralMeshQ == 1)
Recombine Surface {s0};
EndIf

Return
