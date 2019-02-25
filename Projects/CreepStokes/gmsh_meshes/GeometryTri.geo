

IsquadQ = 0;
 
Mesh.ElementOrder = 2;
Mesh.SecondOrderLinear = 0;

lc = 1;
n_bc = 3;
nLayers    = 9;


  Point(1) = {0, -1, 0, 1e+22};
  Point(2) = {2, -1, 0, 1e+22};
  Point(3) = {2, 1, 0, 1e+22};
  Point(4) = {0, 1, 0, 1e+22};

// Fronteiras

  Line(1) = {1, 2};
  Line(2) = {2, 3};
  Line(3) = {3, 4};
  Line(4) = {4, 1};  

  Transfinite Line{1,3} = n_bc;
  Transfinite Line{2,4} = n_bc;

// Quadrados - Furos - Linha 1

   Transfinite Line{1, 2, 3, 4} = nLayers Using Progression 1;

// Definição da superfície 

  Line Loop(1) = {1, 2, 3, 4};
 
  Plane Surface(1) = {1};
 
  Transfinite Surface {1};
  
  If(IsquadQ)

  Recombine Surface {1};

  EndIf

  Physical Surface("Omega") = {1};
  Physical Line("bottom") = {1};
  Physical Line("top") = {3};
  Physical Line("left") = {4};
  Physical Line("right") = {2};
  
  
  Coherence Mesh;




