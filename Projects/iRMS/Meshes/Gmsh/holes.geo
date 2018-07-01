/*********************************************************************
 *
 *  Gmsh holes by using a macro
 *
 *  Characteristic lengths, arrays of variables, macros, loops
 *
 *********************************************************************/

// We start by defining some target mesh sizes:

lcar1 = .25;
lcar2 = .01;
lcar3 = .25;

// If we wanted to change these mesh sizes globally (without changing the above
// definitions), we could give a global scaling factor for all characteristic
// lengths on the command line with the `-clscale' option (or with
// `Mesh.CharacteristicLengthFactor' in an option file). For example, with:
//
// > gmsh t5.geo -clscale 1
//
// this input file produces a mesh of approximately 1,300 nodes and 11,000
// tetrahedra. With
//
// > gmsh t5.geo -clscale 0.2
//
// the mesh counts approximately 350,000 nodes and 2.1 million tetrahedra. You
// can check mesh statistics in the graphical user interface with the
// `Tools->Statistics' menu.

// We proceed by defining some elementary entities describing a truncated cube:

Point(1) = {0,0,0,lcar1};     Point(2) = {1,0,0,lcar1};
Point(3) = {1,1,0,lcar1};     Point(4) = {0,1,0,lcar1};
Point(5) = {0,0,1,lcar1};     Point(6) = {1,0,1,lcar1};
Point(7) = {1,1,1,lcar1};     Point(8) = {0,1,1,lcar1};

Line(1) = {1,2};    Line(2) = {2,3};   Line(3) = {3,4};  Line(4) = {4,1};
Line(5) = {5,6};    Line(6) = {6,7};   Line(7) = {7,8};  Line(8) = {8,5};
Line(9) = {5,1};    Line(10) = {6,2};  Line(11) = {7,3}; Line(12) = {8,4};

Line Loop(1) = {1,2,3,4};    Plane Surface(1) = {1};
Line Loop(2) = {5,6,7,8};    Plane Surface(2) = {2};
Line Loop(3) = {1,-10,-5,9};   Plane Surface(3) = {3};
Line Loop(4) = {2,-11,-6,10};   Plane Surface(4) = {4};
Line Loop(5) = {3,-12,-7,11};   Plane Surface(5) = {5};
Line Loop(6) = {-4,-12,8,9};   Plane Surface(6) = {6};



// Instead of using included files, we now use a user-defined macro in order
// to carve some holes in the cube:

Macro CheeseHole

  // In the following commands we use the reserved variable name `newp', which
  // automatically selects a new point number. This number is chosen as the
  // highest current point number, plus one. (Note that, analogously to `newp',
  // the variables `newl', `news', `newv' and `newreg' select the highest number
  // amongst currently defined curves, surfaces, volumes and `any entities other
  // than points', respectively.)

  p1 = newp; Point(p1) = {x,  y,  z,  lcar2} ;
  p2 = newp; Point(p2) = {x+r,y,  z,  lcar2} ;
  p3 = newp; Point(p3) = {x,  y+r,z,  lcar2} ;
  p4 = newp; Point(p4) = {x,  y,  z+r,lcar2} ;
  p5 = newp; Point(p5) = {x-r,y,  z,  lcar2} ;
  p6 = newp; Point(p6) = {x,  y-r,z,  lcar2} ;
  p7 = newp; Point(p7) = {x,  y,  z-r,lcar2} ;

  c1 = newreg; Circle(c1) = {p2,p1,p7}; 
  c2 = newreg; Circle(c2) = {p7,p1,p5};
  c3 = newreg; Circle(c3) = {p5,p1,p4}; 
  c4 = newreg; Circle(c4) = {p4,p1,p2};
  c5 = newreg; Circle(c5) = {p2,p1,p3}; 
  c6 = newreg; Circle(c6) = {p3,p1,p5};
  c7 = newreg; Circle(c7) = {p5,p1,p6}; 
  c8 = newreg; Circle(c8) = {p6,p1,p2};
  c9 = newreg; Circle(c9) = {p7,p1,p3}; 
  c10 = newreg; Circle(c10) = {p3,p1,p4};
  c11 = newreg; Circle(c11) = {p4,p1,p6}; 
  c12 = newreg; Circle(c12) = {p6,p1,p7};

  // We need non-plane surfaces to define the spherical holes. Here we use ruled
  // surfaces, which can have 3 or 4 sides:

  l1 = newreg; Line Loop(l1) = {c5,c10,c4};   Ruled Surface(newreg) = {l1};
  l2 = newreg; Line Loop(l2) = {c9,-c5,c1};   Ruled Surface(newreg) = {l2};
  l3 = newreg; Line Loop(l3) = {c12,-c8,-c1}; Ruled Surface(newreg) = {l3};
  l4 = newreg; Line Loop(l4) = {c8,-c4,c11};  Ruled Surface(newreg) = {l4};
  l5 = newreg; Line Loop(l5) = {-c10,c6,c3};  Ruled Surface(newreg) = {l5};
  l6 = newreg; Line Loop(l6) = {-c11,-c3,c7}; Ruled Surface(newreg) = {l6};
  l7 = newreg; Line Loop(l7) = {-c2,-c7,-c12};Ruled Surface(newreg) = {l7};
  l8 = newreg; Line Loop(l8) = {-c6,-c9,c2};  Ruled Surface(newreg) = {l8};

  // We then store the surface loops identification numbers in a list for later
  // reference (we will need these to define the final volume):

  holes[t] = newreg;

  Surface Loop(holes[t]) = {l8+1,l5+1,l1+1,l2+1,l3+1,l7+1,l6+1,l4+1};

  // Solid hole
  //solidholes = newreg ;
  //Volume(solidholes) = holes[t] ;

  // Define a physical volume for each hole:
  //Physical Volume ("wellbore_region") = solidholes ;

Return

// We can use a `For' loop to generate one hole in the cube:

x = 0.25 ; y = 0.5 ; z = 0.25 ; r = 0.01 ;


For t In {1:1}

  x += 0.166 ;
  z += 0.166 ;

  // We call the `CheeseHole' macro:

  Call CheeseHole;

  // We also print some variables on the terminal (note that, since all
  // variables are treated internally as floating point numbers, the format
  // string should only contain valid floating point format specifiers like
  // `%g', `%f', '%e', etc.):

  //Printf("Hole %g (center = {%g,%g,%g}, radius = %g) has number %g!",t, x, y, z, r, thehole) ;

EndFor



// We can then define the surface loop for the exterior surface of the cube:

theloops[0] = newreg ;

Surface Loop(holes[0]) = {1,2,3,4,5,6};

// The volume of the cube, without the 5 holes, is now defined by 6 surface
// loops: the first surface loop defines the exterior surface; the surface loops
// other than the first one define holes.  (Again, to reference an array of
// variables, its identifier is followed by square brackets):

Volume(1) = {holes[]} ;

 // transform to hexahedra mesh

//  Transfinite Line "*" = 10 Using Bump 0.25;
//  Transfinite Surface "*";
//  Recombine Surface "*";
//  Transfinite Volume "*";

// We finally define a physical volume for the elements discretizing the cube,
// without the holes (whose elements were already tagged with numbers 1 to 5 in
// the `For' loop):

Physical Volume("reservoir") = {1};
Physical Surface("res_top") = {2};
Physical Surface("res_bottom") = {1};
Physical Surface("res_sides") = {3, 4, 5, 6};
Physical Surface("well") = {holes[1]};

  Printf("Hole %g (center = {%g,%g,%g}, radius = %g) has number %g!",t, x, y, z, r, holes) ;



// We could make only part of the model visible to only mesh this subset:
//
// Hide "*";
// Recursive Show { Volume{129}; }
// Mesh.MeshOnlyVisible = 1;
//+
//+
Coherence;
