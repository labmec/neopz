

IsquadQ = 0;
 
Mesh.ElementOrder = 1;
Mesh.SecondOrderLinear = 0;

//Propriedadae geomtria
rf = 0.25;

//Propriedade da malha
lc = 0.1;
n_bc0 = 12;
n_bc = 6; 
n_h = 8; 
nLayers = 8;

// Pontos contorno 
  Point(1) = {0, 0, 0, 1e+22};
  Point(2) = {8, 0, 0, 1e+22};
  Point(3) = {0, 4, 0, 1e+22};
  Point(4) = {8, 4, 0, 1e+22};

// Pontos domínio dos furos

  Point(5) = {4, 0, 0, 1e+22};
  Point(6) = {5, 0, 0, 1e+22};
  Point(7) = {6, 0, 0, 1e+22};
  Point(8) = {7, 0, 0, 1e+22};
  
  Point(9) = {4, 1, 0, 1e+22};
  Point(10) = {5, 1, 0, 1e+22};
  Point(11) = {6, 1, 0, 1e+22};
  Point(12) = {7, 1, 0, 1e+22};
  Point(13) = {8, 1, 0, 1e+22};
  
  Point(14) = {4, 2, 0, 1e+22};
  Point(15) = {5, 2, 0, 1e+22};
  Point(16) = {6, 2, 0, 1e+22};
  Point(17) = {7, 2, 0, 1e+22};
  Point(18) = {8, 2, 0, 1e+22};

  Point(19) = {4, 3, 0, 1e+22};
  Point(20) = {5, 3, 0, 1e+22};
  Point(21) = {6, 3, 0, 1e+22};
  Point(22) = {7, 3, 0, 1e+22};
  Point(23) = {8, 3, 0, 1e+22};

  Point(24) = {4, 4, 0, 1e+22};
  Point(25) = {5, 4, 0, 1e+22};
  Point(26) = {6, 4, 0, 1e+22};
  Point(27) = {7, 4, 0, 1e+22};

// Pontos nos furos

  // Pontos - Contorno - Furo 1

  Point(28) = {4.5, 0.5, 0, 1e+22};
  Point(29) = {4.5+rf, 0.5, 0, 1e+22};
  Point(30) = {4.5, 0.5+rf, 0, 1e+22};
  Point(31) = {4.5-rf, 0.5, 0, 1e+22};
  Point(32) = {4.5, 0.5-rf, 0, 1e+22};

  // Pontos - Contorno - Furo 2

  Point(33) = {5.5, 0.5, 0, 1e+22};
  Point(34) = {5.5+rf, 0.5, 0, 1e+22};
  Point(35) = {5.5, 0.5+rf, 0, 1e+22};
  Point(36) = {5.5-rf, 0.5, 0, 1e+22};
  Point(37) = {5.5, 0.5-rf, 0, 1e+22};

  // Pontos - Contorno - Furo 3

  Point(38) = {6.5, 0.5, 0, 1e+22};
  Point(39) = {6.5+rf, 0.5, 0, 1e+22};
  Point(40) = {6.5, 0.5+rf, 0, 1e+22};
  Point(41) = {6.5-rf, 0.5, 0, 1e+22};
  Point(42) = {6.5, 0.5-rf, 0, 1e+22};

  // Pontos - Contorno - Furo 4

  Point(43) = {7.5, 0.5, 0, 1e+22};
  Point(44) = {7.5+rf, 0.5, 0, 1e+22};
  Point(45) = {7.5, 0.5+rf, 0, 1e+22};
  Point(46) = {7.5-rf, 0.5, 0, 1e+22};
  Point(47) = {7.5, 0.5-rf, 0, 1e+22};

  // Pontos - Contorno - Furo 5

  Point(48) = {4.5, 1.5, 0, 1e+22};
  Point(49) = {4.5+rf, 1.5, 0, 1e+22};
  Point(50) = {4.5, 1.5+rf, 0, 1e+22};
  Point(51) = {4.5-rf, 1.5, 0, 1e+22};
  Point(52) = {4.5, 1.5-rf, 0, 1e+22};

  // Pontos - Contorno - Furo 6

  Point(53) = {5.5, 1.5, 0, 1e+22};
  Point(54) = {5.5+rf, 1.5, 0, 1e+22};
  Point(55) = {5.5, 1.5+rf, 0, 1e+22};
  Point(56) = {5.5-rf, 1.5, 0, 1e+22};
  Point(57) = {5.5, 1.5-rf, 0, 1e+22};

  // Pontos - Contorno - Furo 7

  Point(58) = {6.5, 1.5, 0, 1e+22};
  Point(59) = {6.5+rf, 1.5, 0, 1e+22};
  Point(60) = {6.5, 1.5+rf, 0, 1e+22};
  Point(61) = {6.5-rf, 1.5, 0, 1e+22};
  Point(62) = {6.5, 1.5-rf, 0, 1e+22};

  // Pontos - Contorno - Furo 8

  Point(63) = {7.5, 1.5, 0, 1e+22};
  Point(64) = {7.5+rf, 1.5, 0, 1e+22};
  Point(65) = {7.5, 1.5+rf, 0, 1e+22};
  Point(66) = {7.5-rf, 1.5, 0, 1e+22};
  Point(67) = {7.5, 1.5-rf, 0, 1e+22};

  // Pontos - Contorno - Furo 9

  Point(68) = {4.5, 2.5, 0, 1e+22};
  Point(69) = {4.5+rf, 2.5, 0, 1e+22};
  Point(70) = {4.5, 2.5+rf, 0, 1e+22};
  Point(71) = {4.5-rf, 2.5, 0, 1e+22};
  Point(72) = {4.5, 2.5-rf, 0, 1e+22};

  // Pontos - Contorno - Furo 10

  Point(73) = {5.5, 2.5, 0, 1e+22};
  Point(74) = {5.5+rf, 2.5, 0, 1e+22};
  Point(75) = {5.5, 2.5+rf, 0, 1e+22};
  Point(76) = {5.5-rf, 2.5, 0, 1e+22};
  Point(77) = {5.5, 2.5-rf, 0, 1e+22};

  // Pontos - Contorno - Furo 11

  Point(78) = {6.5, 2.5, 0, 1e+22};
  Point(79) = {6.5+rf, 2.5, 0, 1e+22};
  Point(80) = {6.5, 2.5+rf, 0, 1e+22};
  Point(81) = {6.5-rf, 2.5, 0, 1e+22};
  Point(82) = {6.5, 2.5-rf, 0, 1e+22};

  // Pontos - Contorno - Furo 12

  Point(83) = {7.5, 2.5, 0, 1e+22};
  Point(84) = {7.5+rf, 2.5, 0, 1e+22};
  Point(85) = {7.5, 2.5+rf, 0, 1e+22};
  Point(86) = {7.5-rf, 2.5, 0, 1e+22};
  Point(87) = {7.5, 2.5-rf, 0, 1e+22};

  // Pontos - Contorno - Furo 13

  Point(88) = {4.5, 3.5, 0, 1e+22};
  Point(89) = {4.5+rf, 3.5, 0, 1e+22};
  Point(90) = {4.5, 3.5+rf, 0, 1e+22};
  Point(91) = {4.5-rf, 3.5, 0, 1e+22};
  Point(92) = {4.5, 3.5-rf, 0, 1e+22};

  // Pontos - Contorno - Furo 14

  Point(93) = {5.5, 3.5, 0, 1e+22};
  Point(94) = {5.5+rf, 3.5, 0, 1e+22};
  Point(95) = {5.5, 3.5+rf, 0, 1e+22};
  Point(96) = {5.5-rf, 3.5, 0, 1e+22};
  Point(97) = {5.5, 3.5-rf, 0, 1e+22};

  // Pontos - Contorno - Furo 15

  Point(98) = {6.5, 3.5, 0, 1e+22};
  Point(99) = {6.5+rf, 3.5, 0, 1e+22};
  Point(100) = {6.5, 3.5+rf, 0, 1e+22};
  Point(101) = {6.5-rf, 3.5, 0, 1e+22};
  Point(102) = {6.5, 3.5-rf, 0, 1e+22};

  // Pontos - Contorno - Furo 16

  Point(103) = {7.5, 3.5, 0, 1e+22};
  Point(104) = {7.5+rf, 3.5, 0, 1e+22};
  Point(105) = {7.5, 3.5+rf, 0, 1e+22};
  Point(106) = {7.5-rf, 3.5, 0, 1e+22};
  Point(107) = {7.5, 3.5-rf, 0, 1e+22};


// Linhas

  Line(1) = {1, 5}; Transfinite Line{1} = n_bc0;
  Line(2) = {5, 6}; Transfinite Line{2} = n_bc;
  Line(3) = {6, 7}; Transfinite Line{3} = n_bc;
  Line(4) = {7, 8}; Transfinite Line{4} = n_bc;
  Line(5) = {8, 2}; Transfinite Line{5} = n_bc;
  
  Line(6) = {1, 3}; Transfinite Line{6} = n_bc0;
  Line(7) = {5, 9}; Transfinite Line{7} = n_bc;
  Line(8) = {6, 10}; Transfinite Line{8} = n_bc;
  Line(9) = {7, 11}; Transfinite Line{9} = n_bc;
  Line(10) = {8, 12}; Transfinite Line{10} = n_bc;
  Line(11) = {2, 13}; Transfinite Line{11} = n_bc;
  
  Line(12) = {9, 10}; Transfinite Line{12} = n_bc;
  Line(13) = {10, 11}; Transfinite Line{13} = n_bc;
  Line(14) = {11, 12}; Transfinite Line{14} = n_bc;
  Line(15) = {12, 13}; Transfinite Line{15} = n_bc;
  
  Line(16) = {9, 14}; Transfinite Line{16} = n_bc;
  Line(17) = {10, 15}; Transfinite Line{17} = n_bc;
  Line(18) = {11, 16}; Transfinite Line{18} = n_bc;
  Line(19) = {12, 17}; Transfinite Line{19} = n_bc;
  Line(20) = {13, 18}; Transfinite Line{20} = n_bc;

  Line(21) = {14, 15}; Transfinite Line{21} = n_bc;
  Line(22) = {15, 16}; Transfinite Line{22} = n_bc;
  Line(23) = {16, 17}; Transfinite Line{23} = n_bc;
  Line(24) = {17, 18}; Transfinite Line{24} = n_bc;
  
  Line(25) = {14, 19}; Transfinite Line{25} = n_bc;
  Line(26) = {15, 20}; Transfinite Line{26} = n_bc;
  Line(27) = {16, 21}; Transfinite Line{27} = n_bc;
  Line(28) = {17, 22}; Transfinite Line{28} = n_bc;
  Line(29) = {18, 23}; Transfinite Line{29} = n_bc;

  Line(30) = {19, 20}; Transfinite Line{30} = n_bc;
  Line(31) = {20, 21}; Transfinite Line{31} = n_bc;
  Line(32) = {21, 22}; Transfinite Line{32} = n_bc;
  Line(33) = {22, 23}; Transfinite Line{33} = n_bc;

  Line(34) = {19, 24}; Transfinite Line{34} = n_bc;
  Line(35) = {20, 25}; Transfinite Line{35} = n_bc;
  Line(36) = {21, 26}; Transfinite Line{36} = n_bc;
  Line(37) = {22, 27}; Transfinite Line{37} = n_bc;
  Line(38) = {23, 4}; Transfinite Line{38} = n_bc;

  Line(39) = {3, 24}; Transfinite Line{39} = n_bc0;
  Line(40) = {24, 25}; Transfinite Line{40} = n_bc;
  Line(41) = {25, 26}; Transfinite Line{41} = n_bc;
  Line(42) = {26, 27}; Transfinite Line{42} = n_bc;
  Line(43) = {27, 4}; Transfinite Line{43} = n_bc;

  Line Loop(201) = {1, 2, 3, 4, 5};
  Line Loop(202) = {11, 20, 29, 38};
  Line Loop(203) = {39, 40, 41, 42, 43};
  Line Loop(204) = {6};

// Círculos - Furos - 1 a 4

  Circle(105) = {29, 28, 30}; Transfinite Line{105} = n_h;
  Circle(106) = {30, 28, 31}; Transfinite Line{106} = n_h;
  Circle(107) = {31, 28, 32}; Transfinite Line{107} = n_h;
  Circle(108) = {32, 28, 29}; Transfinite Line{108} = n_h;

  Circle(109) = {34, 33, 35}; Transfinite Line{109} = n_h;
  Circle(110) = {35, 33, 36}; Transfinite Line{110} = n_h;
  Circle(111) = {36, 33, 37}; Transfinite Line{111} = n_h;
  Circle(112) = {37, 33, 34}; Transfinite Line{112} = n_h;

  Circle(113) = {39, 38, 40}; Transfinite Line{113} = n_h;
  Circle(114) = {40, 38, 41}; Transfinite Line{114} = n_h;
  Circle(115) = {41, 38, 42}; Transfinite Line{115} = n_h;
  Circle(116) = {42, 38, 39}; Transfinite Line{116} = n_h;

  Circle(117) = {44, 43, 45}; Transfinite Line{117} = n_h;
  Circle(118) = {45, 43, 46}; Transfinite Line{118} = n_h;
  Circle(119) = {46, 43, 47}; Transfinite Line{119} = n_h;
  Circle(120) = {47, 43, 44}; Transfinite Line{120} = n_h;

  Line Loop(205) = {105, 106, 107, 108};
  Line Loop(206) = {109, 110, 111, 112};
  Line Loop(207) = {113, 114, 115, 116};
  Line Loop(208) = {117, 118, 119, 120};

// Círculos - Furos - 5 a 8


  Circle(121) = {49, 48, 50}; Transfinite Line{121} = n_h;
  Circle(122) = {50, 48, 51}; Transfinite Line{122} = n_h;
  Circle(123) = {51, 48, 52}; Transfinite Line{123} = n_h;
  Circle(124) = {52, 48, 49}; Transfinite Line{124} = n_h;

  Circle(125) = {54, 53, 55}; Transfinite Line{125} = n_h;
  Circle(126) = {55, 53, 56}; Transfinite Line{126} = n_h;
  Circle(127) = {56, 53, 57}; Transfinite Line{127} = n_h;
  Circle(128) = {57, 53, 54}; Transfinite Line{128} = n_h;

  Circle(129) = {59, 58, 60}; Transfinite Line{129} = n_h;
  Circle(130) = {60, 58, 61}; Transfinite Line{130} = n_h;
  Circle(131) = {61, 58, 62}; Transfinite Line{131} = n_h;
  Circle(132) = {62, 58, 59}; Transfinite Line{132} = n_h;

  Circle(133) = {64, 63, 65}; Transfinite Line{133} = n_h;
  Circle(134) = {65, 63, 66}; Transfinite Line{134} = n_h;
  Circle(135) = {66, 63, 67}; Transfinite Line{135} = n_h;
  Circle(136) = {67, 63, 64}; Transfinite Line{136} = n_h;

  Line Loop(209) = {121, 122, 123, 124};
  Line Loop(210) = {125, 126, 127, 128};
  Line Loop(211) = {129, 130, 131, 132};
  Line Loop(212) = {133, 134, 135, 136};

// Círculos - Furos - 8 a 12

  Circle(137) = {69, 68, 70}; Transfinite Line{137} = n_h;
  Circle(138) = {70, 68, 71}; Transfinite Line{138} = n_h;
  Circle(139) = {71, 68, 72}; Transfinite Line{139} = n_h;
  Circle(140) = {72, 68, 69}; Transfinite Line{140} = n_h;

  Circle(141) = {74, 73, 75}; Transfinite Line{141} = n_h;
  Circle(142) = {75, 73, 76}; Transfinite Line{142} = n_h;
  Circle(143) = {76, 73, 77}; Transfinite Line{143} = n_h;
  Circle(144) = {77, 73, 74}; Transfinite Line{144} = n_h;

  Circle(145) = {79, 78, 80}; Transfinite Line{145} = n_h;
  Circle(146) = {80, 78, 81}; Transfinite Line{146} = n_h;
  Circle(147) = {81, 78, 82}; Transfinite Line{147} = n_h;
  Circle(148) = {82, 78, 79}; Transfinite Line{148} = n_h;

  Circle(149) = {84, 83, 85}; Transfinite Line{149} = n_h;
  Circle(150) = {85, 83, 86}; Transfinite Line{150} = n_h;
  Circle(151) = {86, 83, 87}; Transfinite Line{151} = n_h;
  Circle(152) = {87, 83, 84}; Transfinite Line{152} = n_h;

  Line Loop(213) = {137, 138, 139, 140};
  Line Loop(214) = {141, 142, 143, 144};
  Line Loop(215) = {145, 146, 147, 148};
  Line Loop(216) = {149, 150, 151, 152};

// Círculos - Furos - 12 a 16

  Circle(153) = {89, 88, 90}; Transfinite Line{153} = n_h;
  Circle(154) = {90, 88, 91}; Transfinite Line{154} = n_h;
  Circle(155) = {91, 88, 92}; Transfinite Line{155} = n_h;
  Circle(156) = {92, 88, 89}; Transfinite Line{156} = n_h;

  Circle(157) = {94, 93, 95}; Transfinite Line{157} = n_h;
  Circle(158) = {95, 93, 96}; Transfinite Line{158} = n_h;
  Circle(159) = {96, 93, 97}; Transfinite Line{159} = n_h;
  Circle(160) = {97, 93, 94}; Transfinite Line{160} = n_h;

  Circle(161) = {99, 98, 100}; Transfinite Line{161} = n_h;
  Circle(162) = {100, 98, 101}; Transfinite Line{162} = n_h;
  Circle(163) = {101, 98, 102}; Transfinite Line{163} = n_h;
  Circle(164) = {102, 98, 99}; Transfinite Line{164} = n_h;

  Circle(165) = {104, 103, 105}; Transfinite Line{165} = n_h;
  Circle(166) = {105, 103, 106}; Transfinite Line{166} = n_h;
  Circle(167) = {106, 103, 107}; Transfinite Line{167} = n_h;
  Circle(168) = {107, 103, 104}; Transfinite Line{168} = n_h;

  Line Loop(217) = {153, 154, 155, 156};
  Line Loop(218) = {157, 158, 159, 160};
  Line Loop(219) = {161, 162, 163, 164};
  Line Loop(220) = {165, 166, 167, 168};

// Definição da superfície 
 
  Plane Surface(110) = {201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220};


  holes[] = {105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168}; 
  
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
  Physical Line("bottom") = {1,2,3,4,5};
  Physical Line("top") = {39,40,41,42,43};
  Physical Line("right") = {11,20,29,38};
  Physical Line("left") = {6};
  Physical Line("holes") = {holes[]};

  //Physical Surface("interface") = {23};


  
  Coherence Mesh;




