

Macro RockBox

If (dimension == 3)

Printf ("Mesher:: Building rock box 3D.");


p1 = newp; Point(p1) = {-x_length/2.0, -y_length/2.0, -z_length/2.0, cl};
p2 = newp; Point(p2) = { x_length/2.0, -y_length/2.0, -z_length/2.0, cl};
p3 = newp; Point(p3) = { x_length/2.0,  y_length/2.0, -z_length/2.0, cl};
p4 = newp; Point(p4) = {-x_length/2.0,  y_length/2.0, -z_length/2.0, cl};

p5 = newp; Point(p5) = {-x_length/2.0, -y_length/2.0, z_length/2.0, cl};
p6 = newp; Point(p6) = { x_length/2.0, -y_length/2.0, z_length/2.0, cl};
p7 = newp; Point(p7) = { x_length/2.0,  y_length/2.0, z_length/2.0, cl};
p8 = newp; Point(p8) = {-x_length/2.0,  y_length/2.0, z_length/2.0, cl};

points[] = {p1,p2,p3,p4,p5,p6,p7,p8};

l1 = newl; Line(l1) = {p1,p2};
l2 = newl; Line(l2) = {p2,p3};
l3 = newl; Line(l3) = {p3,p4};
l4 = newl; Line(l4) = {p4,p1};

l5 = newl; Line(l5) = {p5,p6};
l6 = newl; Line(l6) = {p6,p7};
l7 = newl; Line(l7) = {p7,p8};
l8 = newl; Line(l8) = {p8,p5};

l9  = newl; Line(l9)  = {p5,p1};
l10 = newl; Line(l10) = {p6,p2};
l11 = newl; Line(l11) = {p7,p3};
l12 = newl; Line(l12) = {p8,p4};

ll1  = newll; Line Loop(ll1) = {l1,  l2,   l3, l4}; // Bottom
ll2  = newll; Line Loop(ll2) = {l5,  l6,   l7, l8}; // Top
ll3  = newll; Line Loop(ll3) = {l1, -l10, -l5, l9}; // South
ll4  = newll; Line Loop(ll4) = {l2, -l11, -l6, l10}; // East
ll5  = newll; Line Loop(ll5) = {l3, -l12, -l7, l11}; // North
ll6  = newll; Line Loop(ll6) = {l4, -l9,  -l8, l12}; // West

s1  = news; Plane Surface(s1) = {ll1}; // Bottom unstructured region
s2  = news; Plane Surface(s2) = {ll2}; // Top unstructured region
s3  = news; Plane Surface(s3) = {ll3}; // South unstructured region
s4  = news; Plane Surface(s4) = {ll4}; // East unstructured region
s5  = news; Plane Surface(s5) = {ll5}; // North unstructured region
s6  = news; Plane Surface(s6) = {ll6}; // West unstructured region

res_B[] = {s1};
res_T[] = {s2};
res_S[] = {s3};
res_E[] = {s4};
res_N[] = {s5};
res_W[] = {s6};

reservoir_boundaries[] = {res_B[],res_T[],res_S[],res_E[],res_N[],res_W[]};

Else 

// reservoir box dimensions on z-plane
Printf ("Mesher:: Building rock box 2D on z-plane.");
z_length = 0.0;


p1 = newp; Point(p1) = {-x_length/2.0, -y_length/2.0, -z_length/2.0, cl};
p2 = newp; Point(p2) = { x_length/2.0, -y_length/2.0, -z_length/2.0, cl};
p3 = newp; Point(p3) = { x_length/2.0, +y_length/2.0, +z_length/2.0, cl};
p4 = newp; Point(p4) = {-x_length/2.0, +y_length/2.0, +z_length/2.0, cl};

points[] = {p1,p2,p3,p4};

l1 = newl; Line(l1) = {p1,p2};
l2 = newl; Line(l2) = {p2,p3};
l3 = newl; Line(l3) = {p3,p4};
l4 = newl; Line(l4) = {p4,p1};

res_S[] = {l1};
res_E[] = {l2};
res_N[] = {l3};
res_W[] = {l4};
reservoir_boundaries[] = {res_S[],res_E[],res_N[],res_W[]};

lateral[] = {res_E[],res_W[]};
bottom[] = {res_S[]};
top[] = {res_N[]};

N = #lateral[];
M = #laterals[];
For i In {0:N-1}
	laterals[i + M] = lateral[i];
EndFor


N = #bottom[];
M = #bottoms[];
For i In {0:N-1}
	bottoms[i + M] = bottom[i];
EndFor

N = #top[];
M = #tops[];
For i In {0:N-1}
	tops[i + M] = top[i];
EndFor

EndIf

// Apply translation
If (dimension == 2)

wcz = 0.0;

EndIf

translate_p[] = Translate {wcx, wcy, wcz} { 
Point {points[]};
};

Return

