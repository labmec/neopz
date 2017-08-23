lc = 1;
lf = 1;
nLayersl=6;
nLayersc=12;

IsquadQ = 0;
 
Mesh.ElementOrder = 1;
Mesh.SecondOrderLinear = 0;

Point(4) = {0, 1, 0, lc};
Point(5) = {0, 0, 0, lc};
Point(6) = {-1, 0, 0, lc};
Point(7) = {-2, 0, 0, lc};
Point(8) = {0, 2, 0, lc};
Line (7) = {8, 4}; Transfinite Line{7} = nLayersl Using Progression lf;
Line (8) = {6, 7}; Transfinite Line{8} = nLayersl Using Progression lf;
Circle (9) = {7, 5, 8}; Transfinite Line{9} = nLayersc Using Progression lc;
Circle (10) = {4, 5, 6}; Transfinite Line{10} = nLayersc Using Progression lc;
Line Loop(11) = {7, 10, 8, 9};
Plane Surface(12) = {11};

Transfinite Surface{12} = {7, 6, 4, 8};
//Recombine Surface {12};

Physical Surface("Omega") = {12};
Physical Line("bottom") = {8};
Physical Line("top") = {7};
Physical Line("left") = {9};
Physical Line("right") = {10};

Coherence Mesh;