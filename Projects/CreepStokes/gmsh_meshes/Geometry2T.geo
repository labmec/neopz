

IsquadQ = 1;
 
Mesh.ElementOrder = 2;
Mesh.SecondOrderLinear = 0;

lc = 1;
n_bc = 3;
nLayers    = 8;


  Point(1) = {0, 0, 0, 1e+22};
  Point(2) = {1, 0, 0, 1e+22};
  Point(3) = {1, 1, 0, 1e+22};
  Point(4) = {0, 1, 0, 1e+22};

// Pontos - Furos - Linha 1

//  Point(15) = {0.25, 0.25, 0, lc};
//  Point(16) = {0.75, 0.25, 0, lc};
//  Point(17) = {0.75, 0.75, 0, lc};
//  Point(18) = {0.25, 0.75, 0, lc};

// Fronteiras

  Line(1) = {1, 2};
  Line(2) = {2, 3};
  Line(3) = {3, 4};
  Line(4) = {4, 1};  

  Transfinite Line{1,3} = n_bc;
  Transfinite Line{2,4} = n_bc;

// Quadrados - Furos - Linha 1

//  Line(29) = {15, 16};
//  Line(30) = {16, 17};
//  Line(31) = {17, 18};
//  Line(32) = {18, 15};


   Transfinite Line{29, 30, 31, 32} = nLayers Using Progression 1;



// Definição da superfície 

  Line Loop(93) = {1, 2, 3, 4};
  Line Loop(94) = {29, 30, 31, 32};
 
  Plane Surface(110) = {93, 94};
 
  holes[] = {29, 30, 31, 32}; 

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
  Physical Line("top") = {3};
  Physical Line("right") = {2};
  Physical Line("left") = {4};
  
  Physical Line("holes") = {holes[]};  
  //Physical Surface("interface") = {23};

  
  Coherence Mesh;




