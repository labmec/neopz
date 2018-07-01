// ---- iMRS reservoir geometry Gmsh scritp ----
// Creates a mesh with an inner structured-quad region and 
// an outer unstructured triangle region
//
// Created 10/01/2016 by Omar Duran
// Labmec, State University of Campinas
// --------------------------------------------

Mesh.Algorithm3D = 6 ;

IsHexaQ = 1;
IsPrismQ = 1;

xSize = 20;
ySize = 20;
zSize = 10;
esize = 0.5;
deg = 2*Pi/360;
alpha = 90*deg; // angle between hor plane and well
beta = 0*deg; // angle between x and well
n_axial = 1;

r = 5.0; // wellbore radius
// parameters of the horizontal cross-section
e = Cos(alpha);
a = r/Sin(alpha);
b = a*Sqrt(1-e^2);
esize = r;

// creating ellipse inner
p1 = newp; Point(p1) = {a*Cos(beta),a*Sin(beta),0,esize}; // first point of the arc
p2 = newp; Point(p2) = {0,0,0,esize}; // center 
p3 = newp; Point(p3) = {a/2*Cos(beta),a/2*Sin(beta),0,esize}; // point in major axis
p4 = newp; Point(p4) = {b*Sin(beta),b*Cos(beta),0,esize}; // end point of the arc
p5 = newp; Point(p5) = {-a*Cos(beta),-a*Sin(beta),0,esize}; // first point of the arc
p6 = newp; Point(p6) = {-b*Sin(beta),-b*Cos(beta),0,esize}; // end point of the arc
l1 = newl; Ellipse(l1) = {p1,p2,p3,p4};
l2 = newl; Ellipse(l2) = {p4,p2,p3,p5};
l3 = newl; Ellipse(l3) = {p5,p2,p3,p6};
l4 = newl; Ellipse(l4) = {p6,p2,p3,p1};
ll1 = newll; Line Loop(ll1) = {-l4,-l3,-l2,-l1};

r = 10.0; // wellbore region radius
// parameters of the horizontal cross-section
e = Cos(alpha);
a = r/Sin(alpha);
b = a*Sqrt(1-e^2);
esize = 2.0*r;

// creating ellipse mid
mp1 = newp; Point(mp1) = {a*Cos(beta),a*Sin(beta),0,esize}; // first point of the arc
mp2 = newp; Point(mp2) = {0,0,0,esize}; // center 
mp3 = newp; Point(mp3) = {a/2*Cos(beta),a/2*Sin(beta),0,esize}; // point in major axis
mp4 = newp; Point(mp4) = {b*Sin(beta),b*Cos(beta),0,esize}; // end point of the arc
mp5 = newp; Point(mp5) = {-a*Cos(beta),-a*Sin(beta),0,esize}; // first point of the arc
mp6 = newp; Point(mp6) = {-b*Sin(beta),-b*Cos(beta),0,esize}; // end point of the arc
ml1 = newl; Ellipse(ml1) = {mp1,mp2,mp3,mp4};
ml2 = newl; Ellipse(ml2) = {mp4,mp2,mp3,mp5};
ml3 = newl; Ellipse(ml3) = {mp5,mp2,mp3,mp6};
ml4 = newl; Ellipse(ml4) = {mp6,mp2,mp3,mp1};
mll1 = newll; Line Loop(mll1) = {-ml4,-ml3,-ml2,-ml1};

r = 50.0; // reservoir radius
// parameters of the horizontal cross-section
e = Cos(alpha);
a = r/Sin(alpha);
b = a*Sqrt(1-e^2);
esize = r;

// creating ellipse outer
rp1 = newp; Point(rp1) = {a*Cos(beta),a*Sin(beta),0,esize}; // first point of the arc
rp2 = newp; Point(rp2) = {0,0,0,esize}; // center 
rp3 = newp; Point(rp3) = {a/2*Cos(beta),a/2*Sin(beta),0,esize}; // point in major axis
rp4 = newp; Point(rp4) = {b*Sin(beta),b*Cos(beta),0,esize}; // end point of the arc
rp5 = newp; Point(rp5) = {-a*Cos(beta),-a*Sin(beta),0,esize}; // first point of the arc
rp6 = newp; Point(rp6) = {-b*Sin(beta),-b*Cos(beta),0,esize}; // end point of the arc
rl1 = newl; Ellipse(rl1) = {rp1,rp2,rp3,rp4};
rl2 = newl; Ellipse(rl2) = {rp4,rp2,rp3,rp5};
rl3 = newl; Ellipse(rl3) = {rp5,rp2,rp3,rp6};
rl4 = newl; Ellipse(rl4) = {rp6,rp2,rp3,rp1};
rll1 = newll; Line Loop(rll1) = {-rl4,-rl3,-rl2,-rl1};


// wellbore connectors
wbl1 = newl; Line(wbl1) = {p1,mp1};
wbl2 = newl; Line(wbl2) = {p4,mp4};
wbl3 = newl; Line(wbl3) = {p5,mp5};
wbl4 = newl; Line(wbl4) = {p6,mp6};

radial[] = {wbl1,wbl2,wbl3,wbl4};
aximutal[] = {l1,l2,l3,l4,ml1,ml2,ml3,ml4,rl1,rl2,rl3,rl4};

Transfinite Line {radial[]} = 4 Using Progression 2.0;
Transfinite Line {aximutal[]} = 3 Using Progression 1.0;

// reservoir connectors
resl1 = newl; Line(resl1) = {mp1,rp1};
resl2 = newl; Line(resl2) = {mp4,rp4};
resl3 = newl; Line(resl3) = {mp5,rp5};
resl4 = newl; Line(resl4) = {mp6,rp6};

// wellbore surfaces
wbll1 = newll; Line Loop(wbll1) = {-l1,wbl1,ml1,-wbl2}; wbs1 = news; Plane Surface(wbs1) = {wbll1};
wbll2 = newll; Line Loop(wbll2) = {-l2,wbl2,ml2,-wbl3}; wbs2 = news; Plane Surface(wbs2) = {wbll2};
wbll3 = newll; Line Loop(wbll3) = {-l3,wbl3,ml3,-wbl4}; wbs3 = news; Plane Surface(wbs3) = {wbll3};
wbll4 = newll; Line Loop(wbll4) = {-l4,wbl4,ml4,-wbl1}; wbs4 = news; Plane Surface(wbs4) = {wbll4};

wellbore[] = {wbs1,wbs2,wbs3,wbs4};

// reservoir surfaces
resll1 = newll; Line Loop(resll1) = {-ml1,resl1,rl1,-resl2}; ress1 = news; Plane Surface(ress1) = {resll1};
resll2 = newll; Line Loop(resll2) = {-ml2,resl2,rl2,-resl3}; ress2 = news; Plane Surface(ress2) = {resll2};
resll3 = newll; Line Loop(resll3) = {-ml3,resl3,rl3,-resl4}; ress3 = news; Plane Surface(ress3) = {resll3};
resll4 = newll; Line Loop(resll4) = {-ml4,resl4,rl4,-resl1}; ress4 = news; Plane Surface(ress4) = {resll4};

reservoir[] = {ress1,ress2,ress3,ress4};

If(IsHexaQ == 1)

Transfinite Surface "*";
Recombine Surface "*";

EndIf

If(IsPrismQ == 1)
out[] = Extrude{ zSize*Cos(alpha)*Cos(beta), zSize*Cos(alpha)*Sin(beta), zSize}{
Surface{wellbore[],reservoir[]};
Layers{n_axial}; 
Recombine;
};
Else
out[] = Extrude{ zSize*Cos(alpha)*Cos(beta), zSize*Cos(alpha)*Sin(beta), zSize}{
Surface{wellbore[],reservoir[]};
Layers{n_axial};
};
EndIf

// tagging

Physical Surface("well_p_lids") = {};
Physical Surface("well_i_lids") = {};

Physical Surface("producers") = {48,92,114,70};
Physical Surface("injectors") = {};

Physical Volume("Reservoir") = {5,6,7,8};
Physical Volume("Wellbore_p") = {1,2,3,4};
Physical Volume("Wellbore_i") = {};
Physical Surface("Reservoir_bottom") = {37,39,35,33,29,31,25,27};
Physical Surface("Reservoir_top") = {105,61,83,127,193,215,171,149};
Physical Surface("Reservoir_South") = {210};
Physical Surface("Reservoir_East") = {144};
Physical Surface("Reservoir_North") = {166};
Physical Surface("Reservoir_West") = {188};

