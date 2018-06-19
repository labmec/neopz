////////////////////////////////////////////////////////////////
// Plug
// Created 10/03/2018 by Manouchehr Sanei
// Labmec, State University of Campinas, Brazil
////////////////////////////////////////////////////////////////

lc =1.0e1;
r =0.035;
h =0.14;
r =0.25;
h =2.0;
nh = 20;
nr = 3;
Is3DQ = 0;

If(Is3DQ)

////////////////////////////////////////////////////////////////
// 3D mesh
////////////////////////////////////////////////////////////////

p1 = newp; Point(p1) = {0,0,-h/2,lc};
p2 = newp; Point(p2) = {r,0,-h/2,lc};
p3 = newp; Point(p3) = {0,r,-h/2,lc};
p4 = newp; Point(p4) = {-r,0,-h/2,lc};
p5 = newp; Point(p5) = {0,-r,-h/2,lc};

l1 = newl; Circle(l1) = {p2,p1,p3};
l2 = newl; Circle(l2) = {p3,p1,p4};
l3 = newl; Circle(l3) = {p4,p1,p5};
l4 = newl; Circle(l4) = {p5,p1,p2};

ll1 = newll; Line Loop(ll1) = {l1,l2,l3,l4};
s1 = news; Plane Surface(s1) = {ll1};

Transfinite Line {l1,l2,l3,l4} = nr;

Recombine Surface"*";
out[] = Extrude {0,0,h} {
  Surface{s1}; Layers{nh}; Recombine;
};

lateral[] = {15,19,23,27};
top[] = {28};
bottom[] = {s1};
plug[] = {1};


Else

////////////////////////////////////////////////////////////////
// 2D mesh
////////////////////////////////////////////////////////////////

p1 = newp; Point(p1) = {-r,-h/2,0,lc};
p2 = newp; Point(p2) = {r,-h/2,0,lc};
p3 = newp; Point(p3) = {r,h/2,0,lc};
p4 = newp; Point(p4) = {-r,h/2,0,lc};

l1 = newl; Line(l1) = {p1,p2};
l2 = newl; Line(l2) = {p2,p3};
l3 = newl; Line(l3) = {p3,p4};
l4 = newl; Line(l4) = {p4,p1};

ll1 = newll; Line Loop(ll1) = {l1,l2,l3,l4};
s1 = news; Plane Surface(s1) = {ll1};

Transfinite Line {l2,l4} = nh;
Transfinite Line {l1,l3} = nr;
Transfinite Surface {s1};
Recombine Surface"*";

lateral[] = {l2,l4};
top[] = {l3};
bottom[] = {l1};
plug[] = {s1};

EndIf


////////////////////////////////////////////////////////////////
// Physical tagging
////////////////////////////////////////////////////////////////

If(Is3DQ)

Physical Volume("plug") = {plug[]};
Physical Surface("bottom") = {bottom[]};
Physical Surface("top") = {top[]};
Physical Surface("lateral") = {lateral[]};

Else

Physical Surface("plug") = {plug[]};
Physical Line("bottom") = {bottom[]};
Physical Line("top") = {top[]};
Physical Line("lateral") = {lateral[]};



EndIf
