

IsquadQ = 0;
 
Mesh.ElementOrder = 1;
Mesh.SecondOrderLinear = 0;

//Propriedadae geomtria
df = 0.0254;
rf = df/2;

lc = 0.95;
n_bcx1=30;
n_bcy1=30;
n_bcx2=5;
n_bcy2=10;


  Point(1) = {0, 0, 0, 1e+22};
  Point(2) = {0.4318, 0, 0, 1e+22};
  Point(3) = {0.4318, 1.0414, 0, 1e+22};
  Point(4) = {0, 1.0414, 0, 1e+22};

  Point(5) = {0, 0.2794, 0, 1e+22};
  Point(6) = {0.4318, 0.2794, 0, 1e+22};


// Pontos nos furos

// Pontos - Furos - Linha 1

  Point(15) = {df+rf, df+rf, 0, 1e+22};
  Point(16) = {df+rf+rf, df+rf, 0, 1e+22};
  Point(17) = {df+rf, df+rf+rf, 0, 1e+22};
  Point(18) = {df, df+rf, 0, 1e+22};
  Point(19) = {df+rf, df, 0, 1e+22};

  Point(20) = {3*df+rf, df+rf, 0, 1e+22};
  Point(21) = {3*df+rf+rf, df+rf, 0, 1e+22};
  Point(22) = {3*df+rf, df+rf+rf, 0, 1e+22};
  Point(23) = {3*df, df+rf, 0, 1e+22};
  Point(24) = {3*df+rf, df, 0, 1e+22};

  Point(25) = {5*df+rf, df+rf, 0, 1e+22};
  Point(26) = {5*df+rf+rf, df+rf, 0, 1e+22};
  Point(27) = {5*df+rf, df+rf+rf, 0, 1e+22};
  Point(28) = {5*df, df+rf, 0, 1e+22};
  Point(29) = {5*df+rf, df, 0, 1e+22};

  Point(30) = {7*df+rf, df+rf, 0, 1e+22};
  Point(31) = {7*df+rf+rf, df+rf, 0, 1e+22};
  Point(32) = {7*df+rf, df+rf+rf, 0, 1e+22};
  Point(33) = {7*df, df+rf, 0, 1e+22};
  Point(34) = {7*df+rf, df, 0, 1e+22};

  Point(35) = {9*df+rf, df+rf, 0, 1e+22};
  Point(36) = {9*df+rf+rf, df+rf, 0, 1e+22};
  Point(37) = {9*df+rf, df+rf+rf, 0, 1e+22};
  Point(38) = {9*df, df+rf, 0, 1e+22};
  Point(39) = {9*df+rf, df, 0, 1e+22};

  Point(40) = {11*df+rf, df+rf, 0, 1e+22};
  Point(41) = {11*df+rf+rf, df+rf, 0, 1e+22};
  Point(42) = {11*df+rf, df+rf+rf, 0, 1e+22};
  Point(43) = {11*df, df+rf, 0, 1e+22};
  Point(44) = {11*df+rf, df, 0, 1e+22};

  Point(45) = {13*df+rf, df+rf, 0, 1e+22};
  Point(46) = {13*df+rf+rf, df+rf, 0, 1e+22};
  Point(47) = {13*df+rf, df+rf+rf, 0, 1e+22};
  Point(48) = {13*df, df+rf, 0, 1e+22};
  Point(49) = {13*df+rf, df, 0, 1e+22};

  Point(50) = {15*df+rf, df+rf, 0, 1e+22};
  Point(51) = {15*df+rf+rf, df+rf, 0, 1e+22};
  Point(52) = {15*df+rf, df+rf+rf, 0, 1e+22};
  Point(53) = {15*df, df+rf, 0, 1e+22};
  Point(54) = {15*df+rf, df, 0, 1e+22};

// Pontos - Furos - Linha 2

  Point(55) = {df+rf, 3*df+rf, 0, 1e+22};
  Point(56) = {df+rf+rf, 3*df+rf, 0, 1e+22};
  Point(57) = {df+rf, 3*df+rf+rf, 0, 1e+22};
  Point(58) = {df, 3*df+rf, 0, 1e+22};
  Point(59) = {df+rf, 3*df, 0, 1e+22};

  Point(60) = {3*df+rf, 3*df+rf, 0, 1e+22};
  Point(61) = {3*df+rf+rf, 3*df+rf, 0, 1e+22};
  Point(62) = {3*df+rf, 3*df+rf+rf, 0, 1e+22};
  Point(63) = {3*df, 3*df+rf, 0, 1e+22};
  Point(64) = {3*df+rf, 3*df, 0, 1e+22};

  Point(65) = {5*df+rf, 3*df+rf, 0, 1e+22};
  Point(66) = {5*df+rf+rf, 3*df+rf, 0, 1e+22};
  Point(67) = {5*df+rf, 3*df+rf+rf, 0, 1e+22};
  Point(68) = {5*df, 3*df+rf, 0, 1e+22};
  Point(69) = {5*df+rf, 3*df, 0, 1e+22};

  Point(70) = {7*df+rf, 3*df+rf, 0, 1e+22};
  Point(71) = {7*df+rf+rf, 3*df+rf, 0, 1e+22};
  Point(72) = {7*df+rf, 3*df+rf+rf, 0, 1e+22};
  Point(73) = {7*df, 3*df+rf, 0, 1e+22};
  Point(74) = {7*df+rf, 3*df, 0, 1e+22};

  Point(75) = {9*df+rf, 3*df+rf, 0, 1e+22};
  Point(76) = {9*df+rf+rf, 3*df+rf, 0, 1e+22};
  Point(77) = {9*df+rf, 3*df+rf+rf, 0, 1e+22};
  Point(78) = {9*df, 3*df+rf, 0, 1e+22};
  Point(79) = {9*df+rf, 3*df, 0, 1e+22};

  Point(80) = {11*df+rf, 3*df+rf, 0, 1e+22};
  Point(81) = {11*df+rf+rf, 3*df+rf, 0, 1e+22};
  Point(82) = {11*df+rf, 3*df+rf+rf, 0, 1e+22};
  Point(83) = {11*df, 3*df+rf, 0, 1e+22};
  Point(84) = {11*df+rf, 3*df, 0, 1e+22};

  Point(85) = {13*df+rf, 3*df+rf, 0, 1e+22};
  Point(86) = {13*df+rf+rf, 3*df+rf, 0, 1e+22};
  Point(87) = {13*df+rf, 3*df+rf+rf, 0, 1e+22};
  Point(88) = {13*df, 3*df+rf, 0, 1e+22};
  Point(89) = {13*df+rf, 3*df, 0, 1e+22};

  Point(90) = {15*df+rf, 3*df+rf, 0, 1e+22};
  Point(91) = {15*df+rf+rf, 3*df+rf, 0, 1e+22};
  Point(92) = {15*df+rf, 3*df+rf+rf, 0, 1e+22};
  Point(93) = {15*df, 3*df+rf, 0, 1e+22};
  Point(94) = {15*df+rf, 3*df, 0, 1e+22};

// Pontos - Furos - Linha 3

  Point(95) = {df+rf, 5*df+rf, 0, 1e+22};
  Point(96) = {df+rf+rf, 5*df+rf, 0, 1e+22};
  Point(97) = {df+rf, 5*df+rf+rf, 0, 1e+22};
  Point(98) = {df, 5*df+rf, 0, 1e+22};
  Point(99) = {df+rf, 5*df, 0, 1e+22};

  Point(100) = {3*df+rf, 5*df+rf, 0, 1e+22};
  Point(101) = {3*df+rf+rf, 5*df+rf, 0, 1e+22};
  Point(102) = {3*df+rf, 5*df+rf+rf, 0, 1e+22};
  Point(103) = {3*df, 5*df+rf, 0, 1e+22};
  Point(104) = {3*df+rf, 5*df, 0, 1e+22};

  Point(105) = {5*df+rf, 5*df+rf, 0, 1e+22};
  Point(106) = {5*df+rf+rf, 5*df+rf, 0, 1e+22};
  Point(107) = {5*df+rf, 5*df+rf+rf, 0, 1e+22};
  Point(108) = {5*df, 5*df+rf, 0, 1e+22};
  Point(109) = {5*df+rf, 5*df, 0, 1e+22};

  Point(110) = {7*df+rf, 5*df+rf, 0, 1e+22};
  Point(111) = {7*df+rf+rf, 5*df+rf, 0, 1e+22};
  Point(112) = {7*df+rf, 5*df+rf+rf, 0, 1e+22};
  Point(113) = {7*df, 5*df+rf, 0, 1e+22};
  Point(114) = {7*df+rf, 5*df, 0, 1e+22};

  Point(115) = {9*df+rf, 5*df+rf, 0, 1e+22};
  Point(116) = {9*df+rf+rf, 5*df+rf, 0, 1e+22};
  Point(117) = {9*df+rf, 5*df+rf+rf, 0, 1e+22};
  Point(118) = {9*df, 5*df+rf, 0, 1e+22};
  Point(119) = {9*df+rf, 5*df, 0, 1e+22};

  Point(120) = {11*df+rf, 5*df+rf, 0, 1e+22};
  Point(121) = {11*df+rf+rf, 5*df+rf, 0, 1e+22};
  Point(122) = {11*df+rf, 5*df+rf+rf, 0, 1e+22};
  Point(123) = {11*df, 5*df+rf, 0, 1e+22};
  Point(124) = {11*df+rf, 5*df, 0, 1e+22};

  Point(125) = {13*df+rf, 5*df+rf, 0, 1e+22};
  Point(126) = {13*df+rf+rf, 5*df+rf, 0, 1e+22};
  Point(127) = {13*df+rf, 5*df+rf+rf, 0, 1e+22};
  Point(128) = {13*df, 5*df+rf, 0, 1e+22};
  Point(129) = {13*df+rf, 5*df, 0, 1e+22};

  Point(130) = {15*df+rf, 5*df+rf, 0, 1e+22};
  Point(131) = {15*df+rf+rf, 5*df+rf, 0, 1e+22};
  Point(132) = {15*df+rf, 5*df+rf+rf, 0, 1e+22};
  Point(133) = {15*df, 5*df+rf, 0, 1e+22};
  Point(134) = {15*df+rf, 5*df, 0, 1e+22};

// Pontos - Furos - Linha 4

  Point(135) = {df+rf, 7*df+rf, 0, 1e+22};
  Point(136) = {df+rf+rf, 7*df+rf, 0, 1e+22};
  Point(137) = {df+rf, 7*df+rf+rf, 0, 1e+22};
  Point(138) = {df, 7*df+rf, 0, 1e+22};
  Point(139) = {df+rf, 7*df, 0, 1e+22};

  Point(140) = {3*df+rf, 7*df+rf, 0, 1e+22};
  Point(141) = {3*df+rf+rf, 7*df+rf, 0, 1e+22};
  Point(142) = {3*df+rf, 7*df+rf+rf, 0, 1e+22};
  Point(143) = {3*df, 7*df+rf, 0, 1e+22};
  Point(144) = {3*df+rf, 7*df, 0, 1e+22};

  Point(145) = {5*df+rf, 7*df+rf, 0, 1e+22};
  Point(146) = {5*df+rf+rf, 7*df+rf, 0, 1e+22};
  Point(147) = {5*df+rf, 7*df+rf+rf, 0, 1e+22};
  Point(148) = {5*df, 7*df+rf, 0, 1e+22};
  Point(149) = {5*df+rf, 7*df, 0, 1e+22};

  Point(150) = {7*df+rf, 7*df+rf, 0, 1e+22};
  Point(151) = {7*df+rf+rf, 7*df+rf, 0, 1e+22};
  Point(152) = {7*df+rf, 7*df+rf+rf, 0, 1e+22};
  Point(153) = {7*df, 7*df+rf, 0, 1e+22};
  Point(154) = {7*df+rf, 7*df, 0, 1e+22};

  Point(155) = {9*df+rf, 7*df+rf, 0, 1e+22};
  Point(156) = {9*df+rf+rf, 7*df+rf, 0, 1e+22};
  Point(157) = {9*df+rf, 7*df+rf+rf, 0, 1e+22};
  Point(158) = {9*df, 7*df+rf, 0, 1e+22};
  Point(159) = {9*df+rf, 7*df, 0, 1e+22};

  Point(160) = {11*df+rf, 7*df+rf, 0, 1e+22};
  Point(161) = {11*df+rf+rf, 7*df+rf, 0, 1e+22};
  Point(162) = {11*df+rf, 7*df+rf+rf, 0, 1e+22};
  Point(163) = {11*df, 7*df+rf, 0, 1e+22};
  Point(164) = {11*df+rf, 7*df, 0, 1e+22};

  Point(165) = {13*df+rf, 7*df+rf, 0, 1e+22};
  Point(166) = {13*df+rf+rf, 7*df+rf, 0, 1e+22};
  Point(167) = {13*df+rf, 7*df+rf+rf, 0, 1e+22};
  Point(168) = {13*df, 7*df+rf, 0, 1e+22};
  Point(169) = {13*df+rf, 7*df, 0, 1e+22};

  Point(170) = {15*df+rf, 7*df+rf, 0, 1e+22};
  Point(171) = {15*df+rf+rf, 7*df+rf, 0, 1e+22};
  Point(172) = {15*df+rf, 7*df+rf+rf, 0, 1e+22};
  Point(173) = {15*df, 7*df+rf, 0, 1e+22};
  Point(174) = {15*df+rf, 7*df, 0, 1e+22};

// Pontos - Furos - Linha 5

  Point(175) = {df+rf, 9*df+rf, 0, 1e+22};
  Point(176) = {df+rf+rf, 9*df+rf, 0, 1e+22};
  Point(177) = {df+rf, 9*df+rf+rf, 0, 1e+22};
  Point(178) = {df, 9*df+rf, 0, 1e+22};
  Point(179) = {df+rf, 9*df, 0, 1e+22};

  Point(180) = {3*df+rf, 9*df+rf, 0, 1e+22};
  Point(181) = {3*df+rf+rf, 9*df+rf, 0, 1e+22};
  Point(182) = {3*df+rf, 9*df+rf+rf, 0, 1e+22};
  Point(183) = {3*df, 9*df+rf, 0, 1e+22};
  Point(184) = {3*df+rf, 9*df, 0, 1e+22};

  Point(185) = {5*df+rf, 9*df+rf, 0, 1e+22};
  Point(186) = {5*df+rf+rf, 9*df+rf, 0, 1e+22};
  Point(187) = {5*df+rf, 9*df+rf+rf, 0, 1e+22};
  Point(188) = {5*df, 9*df+rf, 0, 1e+22};
  Point(189) = {5*df+rf, 9*df, 0, 1e+22};

  Point(190) = {7*df+rf, 9*df+rf, 0, 1e+22};
  Point(191) = {7*df+rf+rf, 9*df+rf, 0, 1e+22};
  Point(192) = {7*df+rf, 9*df+rf+rf, 0, 1e+22};
  Point(193) = {7*df, 9*df+rf, 0, 1e+22};
  Point(194) = {7*df+rf, 9*df, 0, 1e+22};

  Point(195) = {9*df+rf, 9*df+rf, 0, 1e+22};
  Point(196) = {9*df+rf+rf, 9*df+rf, 0, 1e+22};
  Point(197) = {9*df+rf, 9*df+rf+rf, 0, 1e+22};
  Point(198) = {9*df, 9*df+rf, 0, 1e+22};
  Point(199) = {9*df+rf, 9*df, 0, 1e+22};

  Point(200) = {11*df+rf, 9*df+rf, 0, 1e+22};
  Point(201) = {11*df+rf+rf, 9*df+rf, 0, 1e+22};
  Point(202) = {11*df+rf, 9*df+rf+rf, 0, 1e+22};
  Point(203) = {11*df, 9*df+rf, 0, 1e+22};
  Point(204) = {11*df+rf, 9*df, 0, 1e+22};

  Point(205) = {13*df+rf, 9*df+rf, 0, 1e+22};
  Point(206) = {13*df+rf+rf, 9*df+rf, 0, 1e+22};
  Point(207) = {13*df+rf, 9*df+rf+rf, 0, 1e+22};
  Point(208) = {13*df, 9*df+rf, 0, 1e+22};
  Point(209) = {13*df+rf, 9*df, 0, 1e+22};

  Point(210) = {15*df+rf, 9*df+rf, 0, 1e+22};
  Point(211) = {15*df+rf+rf, 9*df+rf, 0, 1e+22};
  Point(212) = {15*df+rf, 9*df+rf+rf, 0, 1e+22};
  Point(213) = {15*df, 9*df+rf, 0, 1e+22};
  Point(214) = {15*df+rf, 9*df, 0, 1e+22};

// Fronteiras

  Line(1) = {1, 2};   Transfinite Line{1} = n_bcx1;
  Line(2) = {2, 6};   Transfinite Line{2} = n_bcy1;
  Line(3) = {6, 3};   Transfinite Line{3} = n_bcy2; 
  Line(4) = {3, 4};   Transfinite Line{4} = n_bcx2;

  Line(5) = {4,5};    Transfinite Line{5} = n_bcy2;
  Line(6) = {5,1};    Transfinite Line{6} = n_bcy1;
  Line(7) = {5,6};    Transfinite Line{7} = n_bcx1;



// Círculos - Furos - Linha 1

cnum=29;
cc0=15;

  Circle(cnum++) = {cc0+1, cc0, cc0+2}; 
  Circle(cnum++) = {cc0+2, cc0, cc0+3};
  Circle(cnum++) = {cc0+3, cc0, cc0+4};
  Circle(cnum++) = {cc0+4, cc0, cc0+1};
  
  cc0=cc0+5;

  Circle(cnum++) = {cc0+1, cc0, cc0+2};
  Circle(cnum++) = {cc0+2, cc0, cc0+3};
  Circle(cnum++) = {cc0+3, cc0, cc0+4};
  Circle(cnum++) = {cc0+4, cc0, cc0+1};
  
  cc0=cc0+5;


  Circle(cnum++) = {cc0+1, cc0, cc0+2};
  Circle(cnum++) = {cc0+2, cc0, cc0+3};
  Circle(cnum++) = {cc0+3, cc0, cc0+4};
  Circle(cnum++) = {cc0+4, cc0, cc0+1};
  
  cc0=cc0+5;

  Circle(cnum++) = {cc0+1, cc0, cc0+2};
  Circle(cnum++) = {cc0+2, cc0, cc0+3};
  Circle(cnum++) = {cc0+3, cc0, cc0+4};
  Circle(cnum++) = {cc0+4, cc0, cc0+1};
  
  cc0=cc0+5;

  Circle(cnum++) = {cc0+1, cc0, cc0+2};
  Circle(cnum++) = {cc0+2, cc0, cc0+3};
  Circle(cnum++) = {cc0+3, cc0, cc0+4};
  Circle(cnum++) = {cc0+4, cc0, cc0+1};
  
  cc0=cc0+5;

  Circle(cnum++) = {cc0+1, cc0, cc0+2};
  Circle(cnum++) = {cc0+2, cc0, cc0+3};
  Circle(cnum++) = {cc0+3, cc0, cc0+4};
  Circle(cnum++) = {cc0+4, cc0, cc0+1};
  
  cc0=cc0+5;

  Circle(cnum++) = {cc0+1, cc0, cc0+2};
  Circle(cnum++) = {cc0+2, cc0, cc0+3};
  Circle(cnum++) = {cc0+3, cc0, cc0+4};
  Circle(cnum++) = {cc0+4, cc0, cc0+1};
  
  cc0=cc0+5;

  Circle(cnum++) = {cc0+1, cc0, cc0+2};
  Circle(cnum++) = {cc0+2, cc0, cc0+3};
  Circle(cnum++) = {cc0+3, cc0, cc0+4};
  Circle(cnum++) = {cc0+4, cc0, cc0+1};
  
  cc0=cc0+5;

// Círculos - Furos - Linha 2

  Circle(cnum++) = {cc0+1, cc0, cc0+2};
  Circle(cnum++) = {cc0+2, cc0, cc0+3};
  Circle(cnum++) = {cc0+3, cc0, cc0+4};
  Circle(cnum++) = {cc0+4, cc0, cc0+1};
  
  cc0=cc0+5;

  Circle(cnum++) = {cc0+1, cc0, cc0+2};
  Circle(cnum++) = {cc0+2, cc0, cc0+3};
  Circle(cnum++) = {cc0+3, cc0, cc0+4};
  Circle(cnum++) = {cc0+4, cc0, cc0+1};
  
  cc0=cc0+5;


  Circle(cnum++) = {cc0+1, cc0, cc0+2};
  Circle(cnum++) = {cc0+2, cc0, cc0+3};
  Circle(cnum++) = {cc0+3, cc0, cc0+4};
  Circle(cnum++) = {cc0+4, cc0, cc0+1};
  
  cc0=cc0+5;

  Circle(cnum++) = {cc0+1, cc0, cc0+2};
  Circle(cnum++) = {cc0+2, cc0, cc0+3};
  Circle(cnum++) = {cc0+3, cc0, cc0+4};
  Circle(cnum++) = {cc0+4, cc0, cc0+1};
  
  cc0=cc0+5;

  Circle(cnum++) = {cc0+1, cc0, cc0+2};
  Circle(cnum++) = {cc0+2, cc0, cc0+3};
  Circle(cnum++) = {cc0+3, cc0, cc0+4};
  Circle(cnum++) = {cc0+4, cc0, cc0+1};
  
  cc0=cc0+5;

  Circle(cnum++) = {cc0+1, cc0, cc0+2};
  Circle(cnum++) = {cc0+2, cc0, cc0+3};
  Circle(cnum++) = {cc0+3, cc0, cc0+4};
  Circle(cnum++) = {cc0+4, cc0, cc0+1};
  
  cc0=cc0+5;

  Circle(cnum++) = {cc0+1, cc0, cc0+2};
  Circle(cnum++) = {cc0+2, cc0, cc0+3};
  Circle(cnum++) = {cc0+3, cc0, cc0+4};
  Circle(cnum++) = {cc0+4, cc0, cc0+1};
  
  cc0=cc0+5;

  Circle(cnum++) = {cc0+1, cc0, cc0+2};
  Circle(cnum++) = {cc0+2, cc0, cc0+3};
  Circle(cnum++) = {cc0+3, cc0, cc0+4};
  Circle(cnum++) = {cc0+4, cc0, cc0+1};
  
  cc0=cc0+5;


// Círculos - Furos - Linha 3

  Circle(cnum++) = {cc0+1, cc0, cc0+2};
  Circle(cnum++) = {cc0+2, cc0, cc0+3};
  Circle(cnum++) = {cc0+3, cc0, cc0+4};
  Circle(cnum++) = {cc0+4, cc0, cc0+1};
  
  cc0=cc0+5;

  Circle(cnum++) = {cc0+1, cc0, cc0+2};
  Circle(cnum++) = {cc0+2, cc0, cc0+3};
  Circle(cnum++) = {cc0+3, cc0, cc0+4};
  Circle(cnum++) = {cc0+4, cc0, cc0+1};
  
  cc0=cc0+5;


  Circle(cnum++) = {cc0+1, cc0, cc0+2};
  Circle(cnum++) = {cc0+2, cc0, cc0+3};
  Circle(cnum++) = {cc0+3, cc0, cc0+4};
  Circle(cnum++) = {cc0+4, cc0, cc0+1};
  
  cc0=cc0+5;

  Circle(cnum++) = {cc0+1, cc0, cc0+2};
  Circle(cnum++) = {cc0+2, cc0, cc0+3};
  Circle(cnum++) = {cc0+3, cc0, cc0+4};
  Circle(cnum++) = {cc0+4, cc0, cc0+1};
  
  cc0=cc0+5;

  Circle(cnum++) = {cc0+1, cc0, cc0+2};
  Circle(cnum++) = {cc0+2, cc0, cc0+3};
  Circle(cnum++) = {cc0+3, cc0, cc0+4};
  Circle(cnum++) = {cc0+4, cc0, cc0+1};
  
  cc0=cc0+5;

  Circle(cnum++) = {cc0+1, cc0, cc0+2};
  Circle(cnum++) = {cc0+2, cc0, cc0+3};
  Circle(cnum++) = {cc0+3, cc0, cc0+4};
  Circle(cnum++) = {cc0+4, cc0, cc0+1};
  
  cc0=cc0+5;

  Circle(cnum++) = {cc0+1, cc0, cc0+2};
  Circle(cnum++) = {cc0+2, cc0, cc0+3};
  Circle(cnum++) = {cc0+3, cc0, cc0+4};
  Circle(cnum++) = {cc0+4, cc0, cc0+1};
  
  cc0=cc0+5;

  Circle(cnum++) = {cc0+1, cc0, cc0+2};
  Circle(cnum++) = {cc0+2, cc0, cc0+3};
  Circle(cnum++) = {cc0+3, cc0, cc0+4};
  Circle(cnum++) = {cc0+4, cc0, cc0+1};
  
  cc0=cc0+5;

// Círculos - Furos - Linha 4

  Circle(cnum++) = {cc0+1, cc0, cc0+2};
  Circle(cnum++) = {cc0+2, cc0, cc0+3};
  Circle(cnum++) = {cc0+3, cc0, cc0+4};
  Circle(cnum++) = {cc0+4, cc0, cc0+1};
  
  cc0=cc0+5;

  Circle(cnum++) = {cc0+1, cc0, cc0+2};
  Circle(cnum++) = {cc0+2, cc0, cc0+3};
  Circle(cnum++) = {cc0+3, cc0, cc0+4};
  Circle(cnum++) = {cc0+4, cc0, cc0+1};
  
  cc0=cc0+5;


  Circle(cnum++) = {cc0+1, cc0, cc0+2};
  Circle(cnum++) = {cc0+2, cc0, cc0+3};
  Circle(cnum++) = {cc0+3, cc0, cc0+4};
  Circle(cnum++) = {cc0+4, cc0, cc0+1};
  
  cc0=cc0+5;

  Circle(cnum++) = {cc0+1, cc0, cc0+2};
  Circle(cnum++) = {cc0+2, cc0, cc0+3};
  Circle(cnum++) = {cc0+3, cc0, cc0+4};
  Circle(cnum++) = {cc0+4, cc0, cc0+1};
  
  cc0=cc0+5;

  Circle(cnum++) = {cc0+1, cc0, cc0+2};
  Circle(cnum++) = {cc0+2, cc0, cc0+3};
  Circle(cnum++) = {cc0+3, cc0, cc0+4};
  Circle(cnum++) = {cc0+4, cc0, cc0+1};
  
  cc0=cc0+5;

  Circle(cnum++) = {cc0+1, cc0, cc0+2};
  Circle(cnum++) = {cc0+2, cc0, cc0+3};
  Circle(cnum++) = {cc0+3, cc0, cc0+4};
  Circle(cnum++) = {cc0+4, cc0, cc0+1};
  
  cc0=cc0+5;

  Circle(cnum++) = {cc0+1, cc0, cc0+2};
  Circle(cnum++) = {cc0+2, cc0, cc0+3};
  Circle(cnum++) = {cc0+3, cc0, cc0+4};
  Circle(cnum++) = {cc0+4, cc0, cc0+1};
  
  cc0=cc0+5;

  Circle(cnum++) = {cc0+1, cc0, cc0+2};
  Circle(cnum++) = {cc0+2, cc0, cc0+3};
  Circle(cnum++) = {cc0+3, cc0, cc0+4};
  Circle(cnum++) = {cc0+4, cc0, cc0+1};
  
  cc0=cc0+5;

// Círculos - Furos - Linha 5

  Circle(cnum++) = {cc0+1, cc0, cc0+2};
  Circle(cnum++) = {cc0+2, cc0, cc0+3};
  Circle(cnum++) = {cc0+3, cc0, cc0+4};
  Circle(cnum++) = {cc0+4, cc0, cc0+1};
  
  cc0=cc0+5;

  Circle(cnum++) = {cc0+1, cc0, cc0+2};
  Circle(cnum++) = {cc0+2, cc0, cc0+3};
  Circle(cnum++) = {cc0+3, cc0, cc0+4};
  Circle(cnum++) = {cc0+4, cc0, cc0+1};
  
  cc0=cc0+5;


  Circle(cnum++) = {cc0+1, cc0, cc0+2};
  Circle(cnum++) = {cc0+2, cc0, cc0+3};
  Circle(cnum++) = {cc0+3, cc0, cc0+4};
  Circle(cnum++) = {cc0+4, cc0, cc0+1};
  
  cc0=cc0+5;

  Circle(cnum++) = {cc0+1, cc0, cc0+2};
  Circle(cnum++) = {cc0+2, cc0, cc0+3};
  Circle(cnum++) = {cc0+3, cc0, cc0+4};
  Circle(cnum++) = {cc0+4, cc0, cc0+1};
  
  cc0=cc0+5;

  Circle(cnum++) = {cc0+1, cc0, cc0+2};
  Circle(cnum++) = {cc0+2, cc0, cc0+3};
  Circle(cnum++) = {cc0+3, cc0, cc0+4};
  Circle(cnum++) = {cc0+4, cc0, cc0+1};
  
  cc0=cc0+5;

  Circle(cnum++) = {cc0+1, cc0, cc0+2};
  Circle(cnum++) = {cc0+2, cc0, cc0+3};
  Circle(cnum++) = {cc0+3, cc0, cc0+4};
  Circle(cnum++) = {cc0+4, cc0, cc0+1};
  
  cc0=cc0+5;

  Circle(cnum++) = {cc0+1, cc0, cc0+2};
  Circle(cnum++) = {cc0+2, cc0, cc0+3};
  Circle(cnum++) = {cc0+3, cc0, cc0+4};
  Circle(cnum++) = {cc0+4, cc0, cc0+1};
  
  cc0=cc0+5;

  Circle(cnum++) = {cc0+1, cc0, cc0+2};
  Circle(cnum++) = {cc0+2, cc0, cc0+3};
  Circle(cnum++) = {cc0+3, cc0, cc0+4};
  Circle(cnum++) = {cc0+4, cc0, cc0+1};
  
  cc0=cc0+5;

// Definição da superfície 

  Line Loop(93) = {1, 2, 3, 4, 5, 6};
  
  llnum=94;
  llcnum=29;

  Line Loop(llnum++) = {llcnum++, llcnum++, llcnum++, llcnum++};
  Line Loop(llnum++) = {llcnum++, llcnum++, llcnum++, llcnum++};
  Line Loop(llnum++) = {llcnum++, llcnum++, llcnum++, llcnum++};
  Line Loop(llnum++) = {llcnum++, llcnum++, llcnum++, llcnum++};
  Line Loop(llnum++) = {llcnum++, llcnum++, llcnum++, llcnum++};

  Line Loop(llnum++) = {llcnum++, llcnum++, llcnum++, llcnum++};
  Line Loop(llnum++) = {llcnum++, llcnum++, llcnum++, llcnum++};
  Line Loop(llnum++) = {llcnum++, llcnum++, llcnum++, llcnum++};
  Line Loop(llnum++) = {llcnum++, llcnum++, llcnum++, llcnum++};
  Line Loop(llnum++) = {llcnum++, llcnum++, llcnum++, llcnum++};

  Line Loop(llnum++) = {llcnum++, llcnum++, llcnum++, llcnum++};
  Line Loop(llnum++) = {llcnum++, llcnum++, llcnum++, llcnum++};
  Line Loop(llnum++) = {llcnum++, llcnum++, llcnum++, llcnum++};
  Line Loop(llnum++) = {llcnum++, llcnum++, llcnum++, llcnum++};
  Line Loop(llnum++) = {llcnum++, llcnum++, llcnum++, llcnum++};

  Line Loop(llnum++) = {llcnum++, llcnum++, llcnum++, llcnum++};
  Line Loop(llnum++) = {llcnum++, llcnum++, llcnum++, llcnum++};
  Line Loop(llnum++) = {llcnum++, llcnum++, llcnum++, llcnum++};
  Line Loop(llnum++) = {llcnum++, llcnum++, llcnum++, llcnum++};
  Line Loop(llnum++) = {llcnum++, llcnum++, llcnum++, llcnum++};

  Line Loop(llnum++) = {llcnum++, llcnum++, llcnum++, llcnum++};
  Line Loop(llnum++) = {llcnum++, llcnum++, llcnum++, llcnum++};
  Line Loop(llnum++) = {llcnum++, llcnum++, llcnum++, llcnum++};
  Line Loop(llnum++) = {llcnum++, llcnum++, llcnum++, llcnum++};
  Line Loop(llnum++) = {llcnum++, llcnum++, llcnum++, llcnum++};

  Line Loop(llnum++) = {llcnum++, llcnum++, llcnum++, llcnum++};
  Line Loop(llnum++) = {llcnum++, llcnum++, llcnum++, llcnum++};
  Line Loop(llnum++) = {llcnum++, llcnum++, llcnum++, llcnum++};
  Line Loop(llnum++) = {llcnum++, llcnum++, llcnum++, llcnum++};
  Line Loop(llnum++) = {llcnum++, llcnum++, llcnum++, llcnum++};

  Line Loop(llnum++) = {llcnum++, llcnum++, llcnum++, llcnum++};
  Line Loop(llnum++) = {llcnum++, llcnum++, llcnum++, llcnum++};
  Line Loop(llnum++) = {llcnum++, llcnum++, llcnum++, llcnum++};
  Line Loop(llnum++) = {llcnum++, llcnum++, llcnum++, llcnum++};
  Line Loop(llnum++) = {llcnum++, llcnum++, llcnum++, llcnum++};

  Line Loop(llnum++) = {llcnum++, llcnum++, llcnum++, llcnum++};
  Line Loop(llnum++) = {llcnum++, llcnum++, llcnum++, llcnum++};
  Line Loop(llnum++) = {llcnum++, llcnum++, llcnum++, llcnum++};
  Line Loop(llnum++) = {llcnum++, llcnum++, llcnum++, llcnum++};
  Line Loop(llnum++) = {llcnum++, llcnum++, llcnum++, llcnum++};


  Plane Surface(110) = {93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133};


  holes[] = {29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188}; 

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
  Physical Line("bottom") = {1};
  Physical Line("top") = {4};
  Physical Line("left") = {5,6};
  Physical Line("right") = {2,3};
  
  Physical Line("holes") = {holes[]};  
  //Physical Surface("interface") = {23};

  
  Coherence Mesh;




