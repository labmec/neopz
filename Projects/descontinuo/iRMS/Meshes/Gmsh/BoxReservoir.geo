

Macro ReservoirBox

If (dimension == 3)

Printf ("Mesher:: Building reservoir box 3D.");


p1 = newp; Point(p1) = {-x_length/2.0, -y_length/2.0, -z_length/2.0, cl4};
p2 = newp; Point(p2) = { x_length/2.0, -y_length/2.0, -z_length/2.0, cl4};
p3 = newp; Point(p3) = { x_length/2.0,  y_length/2.0, -z_length/2.0, cl4};
p4 = newp; Point(p4) = {-x_length/2.0,  y_length/2.0, -z_length/2.0, cl4};

p5 = newp; Point(p5) = {-x_length/2.0, -y_length/2.0, z_length/2.0, cl4};
p6 = newp; Point(p6) = { x_length/2.0, -y_length/2.0, z_length/2.0, cl4};
p7 = newp; Point(p7) = { x_length/2.0,  y_length/2.0, z_length/2.0, cl4};
p8 = newp; Point(p8) = {-x_length/2.0,  y_length/2.0, z_length/2.0, cl4};

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


If (xzQ == 1)
// reservoir box dimensions on y-plane
Printf ("Mesher:: Building reservoir box 2D on y-plane.");
y_length = 0.0;

Else
// reservoir box dimensions on z-plane
Printf ("Mesher:: Building reservoir box 2D on z-plane.");
z_length = 0.0;

EndIf


p1 = newp; Point(p1) = {-x_length/2.0, -y_length/2.0, -z_length/2.0, cl4};
p2 = newp; Point(p2) = { x_length/2.0, -y_length/2.0, -z_length/2.0, cl4};
p3 = newp; Point(p3) = { x_length/2.0, +y_length/2.0, +z_length/2.0, cl4};
p4 = newp; Point(p4) = {-x_length/2.0, +y_length/2.0, +z_length/2.0, cl4};

l1 = newl; Line(l1) = {p1,p2};
l2 = newl; Line(l2) = {p2,p3};
l3 = newl; Line(l3) = {p3,p4};
l4 = newl; Line(l4) = {p4,p1};

res_S[] = {l1};
res_E[] = {l2};
res_N[] = {l3};
res_W[] = {l4};
reservoir_boundaries[] = {res_S[],res_E[],res_N[],res_W[]};


EndIf



Return

