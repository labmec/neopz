

// Create a drill pipe section

l = 13.716/2; // L = 45 feet
n_axials = 10;
n_azimuthal = 4;

r_inner = 0.1088136/2.0; // ID 4.284 inches
r_outer = 0.1270000/2.0; // OD 5 inches

p0 = newp; Point(p0) = {0,0,0};

pi1 = newp; Point(pi1) = {r_inner,0,0};
pi2 = newp; Point(pi2) = {0,r_inner,0};
pi3 = newp; Point(pi3) = {-r_inner,0,0};
pi4 = newp; Point(pi4) = {0,-r_inner,0};

li1 = newl; Circle(li1) = {pi1,p0,pi2};
li2 = newl; Circle(li2) = {pi2,p0,pi3};
li3 = newl; Circle(li3) = {pi3,p0,pi4};
li4 = newl; Circle(li4) = {pi4,p0,pi1};

po1 = newp; Point(po1) = {r_outer,0,0};
po2 = newp; Point(po2) = {0,r_outer,0};
po3 = newp; Point(po3) = {-r_outer,0,0};
po4 = newp; Point(po4) = {0,-r_outer,0};

lo1 = newl; Circle(lo1) = {po1,p0,po2};
lo2 = newl; Circle(lo2) = {po2,p0,po3};
lo3 = newl; Circle(lo3) = {po3,p0,po4};
lo4 = newl; Circle(lo4) = {po4,p0,po1};

c_innner[] = {li1,li2,li3,li4};
c_outer[] = {lo1,lo2,lo3,lo4};

Transfinite Line {c_innner[],c_outer[]} = n_azimuthal Using Progression 1.0;

lli1 = newll; Line Loop(lli1) = {li1,li2,li3,li4};
llo1 = newll; Line Loop(llo1) = {lo1,lo2,lo3,lo4};


s1 = news; Plane Surface(s1) = {lli1[],llo1[]};

Extrude {0,0,l} {
  Surface{s1}; Layers{n_axials}; Recombine;
}

Recombine Surface {11,53};

Physical Volume("Omega") = {1};
Physical Surface("Bottom") = {11};
Physical Surface("Top") = {53};


