

IsquadQ = 0;
 
Mesh.ElementOrder = 1;
Mesh.SecondOrderLinear = 0;

//Propriedadae geomtria
rf = 0.4;

lc = 0.95;
lf=0.92;
n_bc = 28;
nLayers    = 7;


  Point(1) = {0, 0, 0, 1e+22};
  Point(2) = {8, 0, 0, 1e+22};
  Point(3) = {8, 4, 0, 1e+22};
  Point(4) = {0, 4, 0, 1e+22};

  Point(5) = {4, 0, 0, 1e+22};
  Point(6) = {4, 4, 0, 1e+22};


// Pontos nos furos

// Pontos - Furos - Linha 1

  Point(15) = {4.5, 0.5, 0, 1e+22};
  Point(16) = {4.5+rf, 0.5, 0, 1e+22};
  Point(17) = {4.5, 0.5+rf, 0, 1e+22};
  Point(18) = {4.5-rf, 0.5, 0, 1e+22};
  Point(19) = {4.5, 0.5-rf, 0, 1e+22};

  Point(20) = {5.5, 0.5, 0, 1e+22};
  Point(21) = {5.5+rf, 0.5, 0, 1e+22};
  Point(22) = {5.5, 0.5+rf, 0, 1e+22};
  Point(23) = {5.5-rf, 0.5, 0, 1e+22};
  Point(24) = {5.5, 0.5-rf, 0, 1e+22};

  Point(25) = {6.5, 0.5, 0, 1e+22};
  Point(26) = {6.5+rf, 0.5, 0, 1e+22};
  Point(27) = {6.5, 0.5+rf, 0, 1e+22};
  Point(28) = {6.5-rf, 0.5, 0, 1e+22};
  Point(29) = {6.5, 0.5-rf, 0, 1e+22};

  Point(30) = {7.5, 0.5, 0, 1e+22};
  Point(31) = {7.5+rf, 0.5, 0, 1e+22};
  Point(32) = {7.5, 0.5+rf, 0, 1e+22};
  Point(33) = {7.5-rf, 0.5, 0, 1e+22};
  Point(34) = {7.5, 0.5-rf, 0, 1e+22};

// Pontos - Furos - Linha 2

  Point(35) = {4.5, 1.5, 0, 1e+22};
  Point(36) = {4.5+rf, 1.5, 0, 1e+22};
  Point(37) = {4.5, 1.5+rf, 0, 1e+22};
  Point(38) = {4.5-rf, 1.5, 0, 1e+22};
  Point(39) = {4.5, 1.5-rf, 0, 1e+22};

  Point(40) = {5.5, 1.5, 0, 1e+22};
  Point(41) = {5.5+rf, 1.5, 0, 1e+22};
  Point(42) = {5.5, 1.5+rf, 0, 1e+22};
  Point(43) = {5.5-rf, 1.5, 0, 1e+22};
  Point(44) = {5.5, 1.5-rf, 0, 1e+22};

  Point(45) = {6.5, 1.5, 0, 1e+22};
  Point(46) = {6.5+rf, 1.5, 0, 1e+22};
  Point(47) = {6.5, 1.5+rf, 0, 1e+22};
  Point(48) = {6.5-rf, 1.5, 0, 1e+22};
  Point(49) = {6.5, 1.5-rf, 0, 1e+22};

  Point(50) = {7.5, 1.5, 0, 1e+22};
  Point(51) = {7.5+rf, 1.5, 0, 1e+22};
  Point(52) = {7.5, 1.5+rf, 0, 1e+22};
  Point(53) = {7.5-rf, 1.5, 0, 1e+22};
  Point(54) = {7.5, 1.5-rf, 0, 1e+22};

// Pontos - Furos - Linha 3

  Point(55) = {4.5, 2.5, 0, 1e+22};
  Point(56) = {4.5+rf, 2.5, 0, 1e+22};
  Point(57) = {4.5, 2.5+rf, 0, 1e+22};
  Point(58) = {4.5-rf, 2.5, 0, 1e+22};
  Point(59) = {4.5, 2.5-rf, 0, 1e+22};

  Point(60) = {5.5, 2.5, 0, 1e+22};
  Point(61) = {5.5+rf, 2.5, 0, 1e+22};
  Point(62) = {5.5, 2.5+rf, 0, 1e+22};
  Point(63) = {5.5-rf, 2.5, 0, 1e+22};
  Point(64) = {5.5, 2.5-rf, 0, 1e+22};

  Point(65) = {6.5, 2.5, 0, 1e+22};
  Point(66) = {6.5+rf, 2.5, 0, 1e+22};
  Point(67) = {6.5, 2.5+rf, 0, 1e+22};
  Point(68) = {6.5-rf, 2.5, 0, 1e+22};
  Point(69) = {6.5, 2.5-rf, 0, 1e+22};

  Point(70) = {7.5, 2.5, 0, 1e+22};
  Point(71) = {7.5+rf, 2.5, 0, 1e+22};
  Point(72) = {7.5, 2.5+rf, 0, 1e+22};
  Point(73) = {7.5-rf, 2.5, 0, 1e+22};
  Point(74) = {7.5, 2.5-rf, 0, 1e+22};

// Pontos - Furos - Linha 4

  Point(75) = {4.5, 3.5, 0, 1e+22};
  Point(76) = {4.5+rf, 3.5, 0, 1e+22};
  Point(77) = {4.5, 3.5+rf, 0, 1e+22};
  Point(78) = {4.5-rf, 3.5, 0, 1e+22};
  Point(79) = {4.5, 3.5-rf, 0, 1e+22};

  Point(80) = {5.5, 3.5, 0, 1e+22};
  Point(81) = {5.5+rf, 3.5, 0, 1e+22};
  Point(82) = {5.5, 3.5+rf, 0, 1e+22};
  Point(83) = {5.5-rf, 3.5, 0, 1e+22};
  Point(84) = {5.5, 3.5-rf, 0, 1e+22};

  Point(85) = {6.5, 3.5, 0, 1e+22};
  Point(86) = {6.5+rf, 3.5, 0, 1e+22};
  Point(87) = {6.5, 3.5+rf, 0, 1e+22};
  Point(88) = {6.5-rf, 3.5, 0, 1e+22};
  Point(89) = {6.5, 3.5-rf, 0, 1e+22};

  Point(90) = {7.5, 3.5, 0, 1e+22};
  Point(91) = {7.5+rf, 3.5, 0, 1e+22};
  Point(92) = {7.5, 3.5+rf, 0, 1e+22};
  Point(93) = {7.5-rf, 3.5, 0, 1e+22};
  Point(94) = {7.5, 3.5-rf, 0, 1e+22};


// Fronteiras

  Line(1) = {1, 5};
  Line(2) = {5, 2};
  Line(3) = {2, 3};
  Line(4) = {3, 6};  

  Line(5) = {6,4};
  Line(6) = {4,1};
  Line(7) = {5,6};

  Transfinite Line{2,4} = n_bc;
  Transfinite Line{3,7} = n_bc;

// Círculos - Furos - Linha 1

  Circle(29) = {16, 15, 17};
  Circle(30) = {17, 15, 18};
  Circle(31) = {18, 15, 19};
  Circle(32) = {19, 15, 16};

  Circle(33) = {21, 20, 22};
  Circle(34) = {22, 20, 23};
  Circle(35) = {23, 20, 24};
  Circle(36) = {24, 20, 21};

  Circle(37) = {26, 25, 27};
  Circle(38) = {27, 25, 28};
  Circle(39) = {28, 25, 29};
  Circle(40) = {29, 25, 26};

  Circle(41) = {31, 30, 32};
  Circle(42) = {32, 30, 33};
  Circle(43) = {33, 30, 34};
  Circle(44) = {34, 30, 31};

 Transfinite Line{29, 30, 31, 32} = nLayers Using Progression lf;
 Transfinite Line{33, 34, 35, 36} = nLayers Using Progression lf;
 Transfinite Line{37, 38, 39, 40} = nLayers Using Progression lf;
 Transfinite Line{41, 42, 43, 44} = nLayers Using Progression lf;

// Círculos - Furos - Linha 2

  Circle(45) = {36, 35, 37};
  Circle(46) = {37, 35, 38};
  Circle(47) = {38, 35, 39};
  Circle(48) = {39, 35, 36};

  Circle(49) = {41, 40, 42};
  Circle(50) = {42, 40, 43};
  Circle(51) = {43, 40, 44};
  Circle(52) = {44, 40, 41};

  Circle(53) = {46, 45, 47};
  Circle(54) = {47, 45, 48};
  Circle(55) = {48, 45, 49};
  Circle(56) = {49, 45, 46};

  Circle(57) = {51, 50, 52};
  Circle(58) = {52, 50, 53};
  Circle(59) = {53, 50, 54};
  Circle(60) = {54, 50, 51};

  Transfinite Line{45, 46, 47, 48} = nLayers Using Progression lf;
  Transfinite Line{49, 50, 51, 52} = nLayers Using Progression lf;
  Transfinite Line{53, 54, 55, 56} = nLayers Using Progression lf;
  Transfinite Line{57, 58, 59, 60} = nLayers Using Progression lf;

// Círculos - Furos - Linha 3

  Circle(61) = {56, 55, 57};
  Circle(62) = {57, 55, 58};
  Circle(63) = {58, 55, 59};
  Circle(64) = {59, 55, 56};

  Circle(65) = {61, 60, 62};
  Circle(66) = {62, 60, 63};
  Circle(67) = {63, 60, 64};
  Circle(68) = {64, 60, 61};

  Circle(69) = {66, 65, 67};
  Circle(70) = {67, 65, 68};
  Circle(71) = {68, 65, 69};
  Circle(72) = {69, 65, 66};

  Circle(73) = {71, 70, 72};
  Circle(74) = {72, 70, 73};
  Circle(75) = {73, 70, 74};
  Circle(76) = {74, 70, 71};

  Transfinite Line{61, 62, 63, 64} = nLayers Using Progression lf;
  Transfinite Line{65, 66, 67, 68} = nLayers Using Progression lf;
  Transfinite Line{69, 70, 71, 72} = nLayers Using Progression lf;
  Transfinite Line{73, 74, 75, 76} = nLayers Using Progression lf;

// Círculos - Furos - Linha 4

  Circle(77) = {76, 75, 77};
  Circle(78) = {77, 75, 78};
  Circle(79) = {78, 75, 79};
  Circle(80) = {79, 75, 76};

  Circle(81) = {81, 80, 82};
  Circle(82) = {82, 80, 83};
  Circle(83) = {83, 80, 84};
  Circle(84) = {84, 80, 81};

  Circle(85) = {86, 85, 87};
  Circle(86) = {87, 85, 88};
  Circle(87) = {88, 85, 89};
  Circle(88) = {89, 85, 86};

  Circle(89) = {91, 90, 92};
  Circle(90) = {92, 90, 93};
  Circle(91) = {93, 90, 94};
  Circle(92) = {94, 90, 91};

  Transfinite Line{77, 78, 79, 80} = nLayers Using Progression lf;
  Transfinite Line{81, 82, 83, 84} = nLayers Using Progression lf;
  Transfinite Line{85, 86, 87, 88} = nLayers Using Progression lf;
  Transfinite Line{89, 90, 91, 92} = nLayers Using Progression lf;

// Definição da superfície 

  Line Loop(93) = {1, 2, 3, 4, 5, 6};
  Line Loop(94) = {29, 30, 31, 32};
  Line Loop(95) = {33, 34, 35, 36};
  Line Loop(96) = {37, 38, 39, 40};
  Line Loop(97) = {41, 42, 43, 44};
  Line Loop(98) = {45, 46, 47, 48};
  Line Loop(99) = {49, 50, 51, 52};
  Line Loop(100) = {53, 54, 55, 56};
  Line Loop(101) = {57, 58, 59, 60};
  Line Loop(102) = {61, 62, 63, 64};
  Line Loop(103) = {65, 66, 67, 68};
  Line Loop(104) = {69, 70, 71, 72};
  Line Loop(105) = {73, 74, 75, 76};
  Line Loop(106) = {77, 78, 79, 80};
  Line Loop(107) = {81, 82, 83, 84};
  Line Loop(108) = {85, 86, 87, 88};
  Line Loop(109) = {89, 90, 91, 92};
  Plane Surface(110) = {93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109};


  holes[] = {29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92}; 

  If(IsquadQ)

  Recombine Surface {110};

  EndIf

 // Physical Volume("internal") = {1};
 // Extrude {0, 0, 10} {
 //  //Surface{110};
 //  Layers{1};
 //  Recombine;
 // }

  Physical Surface("Omega") = {110};
  Physical Line("bottom") = {1,2};
  Physical Line("top") = {4,5};
  Physical Line("left") = {3};
  Physical Line("right") = {6};
  
  Physical Line("holes") = {holes[]};  
  //Physical Surface("interface") = {23};

  
  Coherence Mesh;




