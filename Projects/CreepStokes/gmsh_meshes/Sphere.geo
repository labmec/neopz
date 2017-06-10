// ---- Gmsh Macro ----
// ---- Ellipsoidal wellbore region  ----
// Created 05/12/2016 by Omar Duran
// Labmec, State University of Campinas
// --------------------------------------------

Macro MakeSphere

If (NonLinearQ == 1)
Mesh.ElementOrder = 2;
Mesh.SecondOrderLinear = 0;
EndIf

l = 1.0;
progr = 1.0;

Point(1) = {0,0,0,l};

// exterior sphere
Point(2) = {r1,r1,-r1,l};
Point(3) = {-r1,r1,-r1,l};
Point(4) = {-r1,-r1,-r1,l};
Point(5) = {r1,-r1,-r1,l};
Circle(1) = {3,1,2};
Circle(2) = {2,1,5};
Circle(3) = {5,1,4};
Circle(4) = {4,1,3};
Line Loop(5) = {1,2,3,4};
Ruled Surface(6) = {5};
Rotate { {1,0,0},{0,0,0}, Pi/2 } { Duplicata{ Surface{6}; } }
Rotate { {1,0,0},{0,0,0}, Pi } { Duplicata{ Surface{6}; } }
Rotate { {1,0,0},{0,0,0}, 3*Pi/2 } { Duplicata{ Surface{6}; } }
Rotate { {0,1,0},{0,0,0}, Pi/2 } { Duplicata { Surface{6}; } }
Rotate { {0,1,0},{0,0,0}, -Pi/2 } { Duplicata { Surface{6}; } }


// interior sphere
Point(102) = {r2,r2,-r2,l};
Point(103) = {-r2,r2,-r2,l};
Point(104) = {-r2,-r2,-r2,l};
Point(105) = {r2,-r2,-r2,l};
Circle(29) = {103,1,102};
Circle(30) = {102,1,105};
Circle(31) = {105,1,104};
Circle(32) = {104,1,103};
Line Loop(33) = {29,30,31,32};
Ruled Surface(34) = {33};
Rotate { {1,0,0},{0,0,0}, Pi/2 } { Duplicata{ Surface{34}; } }
Rotate { {1,0,0},{0,0,0}, Pi } { Duplicata{ Surface{34}; } }
Rotate { {1,0,0},{0,0,0}, 3*Pi/2 } { Duplicata{ Surface{34}; } }
Rotate { {0,1,0},{0,0,0}, Pi/2 } { Duplicata { Surface{34}; } }
Rotate { {0,1,0},{0,0,0}, -Pi/2 } { Duplicata { Surface{34}; } }

// connect spheres
Line(52) = {102,2};
Line(53) = {108,8};
Line(54) = {105,5};
Line(55) = {111,11};
Line(56) = {109,9};
Line(57) = {104,4};
Line(58) = {103,3};
Line(59) = {106,6};

Line Loop(60) = {58,1,-52,-29};Plane Surface(61) = {60};
Line Loop(62) = {58,11,-59,-39};Plane Surface(63) = {62};
Line Loop(64) = {59,8,-53,-36};Plane Surface(65) = {64};
Line Loop(66) = {37,52,-9,-53};Plane Surface(67) = {66};
Line Loop(68) = {56, 21,-57,-49};Plane Surface(69) = {68};
Line Loop(70) = {31,57,-3,-54};Plane Surface(71) = {70};
Line Loop(72) = {54,19,-55,-47};Plane Surface(73) = {72};
Line Loop(74) = {55,-13,-56,41};Plane Surface(75) = {74};
Line Loop(76) = {59,16,-56,-44};Plane Surface(77) = {76};
Line Loop(78) = {58,-4,-57,32};Plane Surface(79) = {78};
Line Loop(80) = {52,2,-54,-30};Plane Surface(81) = {80};
Line Loop(82) = {42,53,-14,-55};Plane Surface(83) = {82};

// connection volumes
Surface Loop(84) = {7,61,-63,-65,67,-35}; Volume(85) = {84};
Surface Loop(86) = {34,61,-79,6,81,-71}; Volume(87) = {86};
Surface Loop(88) = {22,-79,63,77,69,-50}; Volume(89) = {88};
Surface Loop(90) = {12,83,-40,75,-77,65}; Volume(91) = {90};
Surface Loop(92) = {23,81,-67,-51,-83,73}; Volume(93) = {92};
Surface Loop(94) = {17,-71,-45,-73,-75,-69}; Volume(95) = {94};

// define transfinite mesh
Transfinite Line {53, 59, 52, 58, 55, 56, 54, 57} = n2 Using Progression progr;
Transfinite Line {1, 2, 3, 4, 9, 11,8,16,14,13,19,21,41,42,44,36,39,37, 29, 30,32, 49,31,47} = n1;

radials[] = {61,63,65,67,69,71,73,75,77,79,81,83};
inner[] = {6,7,12,17,22,23};
outer[] = {34,35,40,45,50,51};

Transfinite Surface "*";
Transfinite Volume "*";



If(IsTetraQ)

Else

If(IsPrismQ)
Recombine Surface{radials[]};
Else
Recombine Surface"*";
EndIf

EndIf


// Tagging physical entities
Physical Volume("volume") = {85,87,89,91,93,95};
Physical Surface("inner") = {inner[]};
Physical Surface("outer") = {outer[]};






Return
