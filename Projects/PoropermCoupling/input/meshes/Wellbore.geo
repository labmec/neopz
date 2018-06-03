
////////////////////////////////////////////////////////////////
// 2D wellbore
// Created 28/05/2018 by Manouchehr Sanei
// Labmec, State University of Campinas
////////////////////////////////////////////////////////////////

IsquadQ = 1;
 
Mesh.ElementOrder = 1;
Mesh.SecondOrderLinear = 0;

lf = 0.01;
lc=1;
h=4.0;
r=0.1;
nh = 20;
nv = 20;
nr = 20;

Point(1) = {-h,-h,0,lc};
Point(2) = {h,-h,0,lc};
Point(3) = {h,h,0,lc};
Point(4) = {-h,h,0,lc};

Point(5) = {0,0,0,lc};
Point(6) = {r,0,0,lc};
Point(7) = {0,r,0,lc};
Point(8) = {-r,0,0,lc};
Point(9) = {0,-r,0,lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Circle(5) = {6,5,7};
Circle(6) = {7,5,8};
Circle(7) = {8,5,9};
Circle(8) = {9,5,6};


Line Loop(1) = {1,2,3,4};
Line Loop(2) = {6,7,8,5};
Plane Surface(1) = {1,2};

fixed_y_points[]={6,8};
fixed_x_points[]={7,9};

Point{fixed_y_points[],fixed_x_points[]} In Surface{1};

Transfinite Line {2,4} = nh;
Transfinite Line {1,3} = nv;
Transfinite Line {5,6,7,8} = nr;


holes[] = {5,6,7,8};

 If(IsquadQ)
  Recombine Surface {1};
 EndIf


Physical Surface("Omega") = {1};
Physical Line("bottom") = {1};
Physical Line("top") = {3};
Physical Line("left") = {4};
Physical Line("right") = {2};  
Physical Line("holes") = {holes[]};  

Physical Point("fixed_x") = {fixed_x_points[]};
Physical Point("fixed_y") = {fixed_y_points[]};

Coherence Mesh;

