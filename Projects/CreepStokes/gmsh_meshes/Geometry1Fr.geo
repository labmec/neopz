

IsquadQ = 0;
 
Mesh.ElementOrder = 1;
Mesh.SecondOrderLinear = 0;

lc = 1;
n_l = 22;
n_c = 22;
n_b = 22;
pl = 1;
pc = 1;
pb = 0.95;

//nLayers    = 22;
rc=0.25;
pc1=(0.7071067812-rc)*0.7071067812;

  Point(1) = {0, 0, 0, 1e+22};
  Point(2) = {1, 0, 0, 1e+22};
  Point(3) = {1, 1, 0, 1e+22};
  Point(4) = {0, 1, 0, 1e+22};

// Pontos - Furos - Linha 1

  Point(15) = {pc1, pc1, 0, lc};
  Point(16) = {pc1+rc*0.7071067812*2, pc1, 0, lc};
  Point(17) = {pc1+rc*0.7071067812*2, pc1+rc*0.7071067812*2, 0, lc};
  Point(18) = {pc1, pc1+rc*0.7071067812*2, 0, lc};
  Point(19) = {0.5, 0.5, 0, lc};

// Fronteiras

  Line(1) = {1, 2};  Transfinite Line{1} = n_l Using Progression pl;
  Line(2) = {2, 3};  Transfinite Line{2} = n_l Using Progression pl;  
  Line(3) = {3, 4};  Transfinite Line{3} = n_l Using Progression pl; 
  Line(4) = {4, 1};  Transfinite Line{4} = n_l Using Progression pl;
 
// Círculos - Furos - Linha 1

  Circle(5) = {15, 19, 16}; Transfinite Line{5} = n_c Using Progression pc;
  Circle(6) = {16, 19, 17}; Transfinite Line{6} = n_c Using Progression pc;
  Circle(7) = {17, 19, 18}; Transfinite Line{7} = n_c Using Progression pc;
  Circle(8) = {18, 19, 15}; Transfinite Line{8} = n_c Using Progression pc;

// Block Lines

  Line(9) = {1, 15};  Transfinite Line{9} = n_b Using Progression pb;
  Line(10) = {2, 16};  Transfinite Line{10} = n_b Using Progression pb;  
  Line(11) = {3, 17};  Transfinite Line{11} = n_b Using Progression pb; 
  Line(12) = {4, 18};  Transfinite Line{12} = n_b Using Progression pb;


// surfaces

  Line Loop(13) = {9, 5, -10, -1};
  Plane Surface(14) = {13};
  
  Line Loop(15) = {10, 6, -11, -2};
  Plane Surface(16) = {15};

  Line Loop(17) = {11, 7, -12, -3};
  Plane Surface(18) = {17};

  Line Loop(19) = {12, 8, -9, -4};
  Plane Surface(20) = {19};

  Transfinite Surface {14};
  Transfinite Surface {16};
  Transfinite Surface {18};
  Transfinite Surface {20};

 // Recombine Surface {14};
 // Recombine Surface {16};
 // Recombine Surface {18};
 // Recombine Surface {20};

  holes[] = {5, 6, 7, 8}; 
  
//  If(IsquadQ)
//  Recombine Surface {110};
// EndIf


 // Physical Volume("internal") = {1};
 // Extrude {0, 0, 10} {
 //  //Surface{110};
 //  Layers{1};
 //  Recombine;
 // }

  Physical Surface("Omega") = {14,16,18,20};
  Physical Line("bottom") = {1};
  Physical Line("top") = {3};
  Physical Line("right") = {2};
  Physical Line("left") = {4};
  
  Physical Line("holes") = {holes[]};  
  //Physical Surface("interface") = {23};

  
  Coherence Mesh;




