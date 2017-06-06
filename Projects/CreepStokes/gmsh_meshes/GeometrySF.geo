

IsquadQ = 0;
 
Mesh.ElementOrder = 2;
Mesh.SecondOrderLinear = 0;

lc = 0.1;
n_bc = 10;
nLayers    = 12;


  Point(1) = {0, 0, 0, 1e+22};
  Point(2) = {8, 0, 0, 1e+22};
  Point(3) = {8, 4, 0, 1e+22};
  Point(4) = {0, 4, 0, 1e+22};

  Point(5) = {4, 0, 0, 1e+22};
  Point(6) = {4, 4, 0, 1e+22};


// Fronteiras

  Line(1) = {1, 5};
  Line(2) = {5, 2};
  Line(3) = {2, 3};
  Line(4) = {3, 6};  

  Line(5) = {6,4};
  Line(6) = {4,1};
  Line(7) = {5,6};

  Transfinite Line{1,5} = n_bc;
  Transfinite Line{6,7} = n_bc;

  Transfinite Line{2,4} = n_bc;
  Transfinite Line{3,7} = n_bc;

// Círculos - Furos - Linha 1


// Definição da superfície 

  Line Loop(93) = {1, 2, 3, 4, 5, 6};

  Plane Surface(110) = {93};



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
  
  //Physical Line("holes") = {holes[]};  
  //Physical Surface("interface") = {23};

  
  Coherence Mesh;




