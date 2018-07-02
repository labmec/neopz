// ---- Gmsh Macro ----
// ---- Circle region  ----
// Created 05/12/2016 by Omar Duran
// Labmec, State University of Campinas
// --------------------------------------------

Macro MakeCircle

If (NonLinearQ == 1)
Mesh.ElementOrder = 2;
Mesh.SecondOrderLinear = 0;
EndIf

l = 0.5;
progr = 1.5;

h=0.5;
r1 = outer_r/Sqrt(2.0);
r2 = inner_r/Sqrt(2.0);

p1 = newp; Point(p1) = {0,0,-h/2,l};

// exterior sphere
p2 = newp; Point(p2) = {r1,r1,-h/2,l};
p3 = newp; Point(p3) = {-r1,r1,-h/2,l};
p4 = newp; Point(p4) = {-r1,-r1,-h/2,l};
p5 = newp; Point(p5) = {r1,-r1,-h/2,l};
l1 = newl; Circle(l1) = {p3,p1,p2};
l2 = newl; Circle(l2) = {p2,p1,p5};
l3 = newl; Circle(l3) = {p5,p1,p4};
l4 = newl; Circle(l4) = {p4,p1,p3};


// interior sphere
pi2 = newp; Point(pi2) = {r2,r2,-h/2,l};
pi3 = newp; Point(pi3) = {-r2,r2,-h/2,l};
pi4 = newp; Point(pi4) = {-r2,-r2,-h/2,l};
pi5 = newp; Point(pi5) = {r2,-r2,-h/2,l};
li1 = newl; Circle(li1) = {pi3,p1,pi2};
li2 = newl; Circle(li2) = {pi2,p1,pi5};
li3 = newl; Circle(li3) = {pi5,p1,pi4};
li4 = newl; Circle(li4) = {pi4,p1,pi3};

lr1  = newl; Line(lr1)  = {p2,pi2};
lr2  = newl; Line(lr2)  = {p3,pi3};
lr3  = newl; Line(lr3)  = {p4,pi4};
lr4  = newl; Line(lr4)  = {p5,pi5};

bottom  = newll; Line Loop(bottom) = {-l1,  lr2,   li1, -lr1}; // Bottom
s1  = news; Plane Surface(s1) = {bottom}; // bottom region

ritgh  = newll; Line Loop(ritgh) = {-l2,  lr1,   li2, -lr4}; // ritgh
s2  = news; Plane Surface(s2) = {ritgh}; // ritgh region

top  = newll; Line Loop(top) = {-l3,  lr4,   li3, -lr3}; // top
s3  = news; Plane Surface(s3) = {top}; // top region

left  = newll; Line Loop(left) = {-l4,  lr3,  li4, -lr2}; // left
s4  = news; Plane Surface(s4) = {left}; // left region


// define transfinite mesh
Transfinite Line {-lr1,-lr2,-lr3,-lr4} = n2 Using Progression progr;
Transfinite Line {-li1,-li2,-li3,-li4} = n1;

disc[] = {s1,s2,s3,s4};
inner[] = {li1,li2,li3,li4};
outer[] = {l1,l2,l3,l4};



// Tagging physical entities
Physical Surface("volume") = {disc[]};
Physical Line("inner") = {inner[]};
Physical Line("outer") = {outer[]};


Return
