Mesh.CharacteristicLengthFromPoints = 0;
Mesh.CharacteristicLengthExtendFromBoundary = 0;

Mesh.CharacteristicLengthFromPoints = 0;
Mesh.CharacteristicLengthExtendFromBoundary = 0;

// Attractor field on point 1000. This field returns the
// distance to point 1000.
Field[1] = Attractor; 
Field[1].FacesList = {well_p_bores[],well_lids[]}; // works partially 

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = 2.0;
Field[2].LcMax = 200.0;
Field[2].DistMin = 0.01;
Field[2].DistMax = 100.0;

Background Field = 2;

lc = 0.1; 
Point(1) = {0.619574, 1.427862, 0.404304, lc}; 
Point(2) = {0.823514, 1.427862, 0.442334, lc}; 
Point(3) = {0.444683, 1.427862, 1.329681, lc}; 
Point(4) = {0.909917, 1.427862, 1.076947, lc}; 
Point(5) = {0.620869, 1.413485, 0.423666, lc}; 
Point(6) = {0.461935, 1.413485, 1.264608, lc}; 
Point(7) = {0.822523, 1.413485, 0.461270, lc}; 
Point(8) = {0.899532, 1.413485, 1.026888, lc}; 
Point(9) = {0.949253, 1.601767, 1.048794, lc}; 
Point(10) = {0.379330, 1.601767, 1.358400, lc}; 
Point(11) = {0.561264, 1.601767, 0.395760, lc}; 
Point(12) = {0.868133, 1.601767, 0.452984, lc}; 
Line(13) = {9, 10}; 
Line(14) = {10, 11}; 
Line(15) = {11, 12}; 
Line(16) = {12, 9}; 
Line(17) = {1, 11}; 
Line(18) = {10, 3}; 
Line(19) = {3, 6}; 
Line(20) = {6, 5}; 
Line(21) = {5, 1}; 
Line(22) = {6, 8}; 
Line(23) = {8, 7}; 
Line(24) = {7, 5}; 
Line(25) = {3, 4}; 
Line(26) = {4, 9}; 
Line(27) = {4, 8}; 
Line(28) = {2, 12}; 
Line(29) = {7, 2}; 
Line(30) = {1, 2}; 
Line Loop(31) = {13, 14, 15, 16}; 
Line Loop(32) = {17, -14, 18, 19, 20, 21}; 
Line Loop(33) = {-20, 22, 23, 24}; 
Line Loop(34) = {25, 26, 13, 18}; 
Line Loop(35) = {25, 27, -22, -19}; 
Line Loop(36) = {28, 16, -26, 27, 23, 29}; 
Line Loop(37) = {30, 28, -15, -17}; 
Line Loop(38) = {30, -29, 24, 21}; 
Plane Surface(39) = {31}; 
Plane Surface(40) = {32}; 
Plane Surface(41) = {33}; 
Plane Surface(42) = {34}; 
Plane Surface(43) = {35}; 
Plane Surface(44) = {36}; 
Plane Surface(45) = {37}; 
Plane Surface(46) = {38}; 
Surface Loop(47) = {39, 40, 41, 42, 43, 44, 45, 46}; 
Volume(48) = {47}; 

Field[1] = Attractor; 
Field[1].EdgesList = {25}; // works 
//Field[1].FacesList = {43}; // works partially 
//Field[1].FacesList = {42}; // does not work 

Field[2] = Threshold; 
Field[2].IField = 1; 
Field[2].LcMin = lc / 4; 
Field[2].LcMax = lc; 
Field[2].DistMin = 0.05; 
Field[2].DistMax = 0.2; 

Background Field = 2; 