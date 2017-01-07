/*********************************************************************
 *
 *  Gmsh holes by using a macro
 *
 *  Characteristic lengths, arrays of variables, macros, loops
 *
 *********************************************************************/

ft = 0.3048;
xSize = 10*ft;
ySize = 10*ft;
zSize = 10*ft;
esize = 0.125*ft;
deg = 2*Pi/360;
alpha = 85*deg; // angle between hor plane and well
beta = 0*deg; // angle between x and well
r = 3*ft; // well radius
// parameters of the horizontal cross-section
e = Cos(alpha);
a = r/Sin(alpha);
b = a*Sqrt(1-e^2);

// creating ellipse
Point(1) = {a*Cos(beta),a*Sin(beta),0,esize}; // first point of the arc
Point(2) = {0,0,0,esize}; // center 
Point(3) = {a/2*Cos(beta),a/2*Sin(beta),0,esize}; // point in major axis
Point(4) = {b*Sin(beta),b*Cos(beta),0,esize}; // end point of the arc
Point(5) = {-a*Cos(beta),-a*Sin(beta),0,esize}; // first point of the arc
Point(6) = {-b*Sin(beta),-b*Cos(beta),0,esize}; // end point of the arc
Ellipse(1) = {1,2,3,4};
Ellipse(2) = {4,2,3,5};
Ellipse(3) = {5,2,3,6};
Ellipse(4) = {6,2,3,1};
// Line Loop(5) = {1,2,3,4};
Line Loop(5) = {-4,-3,-2,-1};

out[] = Extrude{
zSize*Cos(alpha)*Cos(beta),
zSize*Cos(alpha)*Sin(beta),
zSize}
{
// Surface{6};
// Loop{5};
Line{-4,-3,-2,-1};
};
// // Points
// Bottom
Point(31) = {xSize/2,ySize/2,0,esize};
Point(32) = {-xSize/2,ySize/2,0,esize};
Point(33) = {-xSize/2,-ySize/2,0,esize};
Point(34) = {xSize/2,-ySize/2,0,esize};
// Top
Point(35) = {xSize/2,ySize/2,zSize,esize};
Point(36) = {-xSize/2,ySize/2,zSize,esize};
Point(37) = {-xSize/2,-ySize/2,zSize,esize};
Point(38) = {xSize/2,-ySize/2,zSize,esize};
// Lines
// Bottom
Line(31) = {31,32};
Line(32) = {32,33};
Line(33) = {33,34};
Line(34) = {34,31};

Line Loop(35) = {31,32,33,34};
surf2 = news;
Plane Surface(surf2) = {35,5}; // bottom surf
// // Top
Line(36) = {35,36};
Line(37) = {36,37};
Line(38) = {37,38};
Line(39) = {38,35};
Line Loop(40) = {36,37,38,39};
Line Loop(41) = {out[0],out[5],out[10],out[15]}; // upper ellipse
// Plane Surface(666) = {41};
surf3 = news;
Plane Surface(surf3) = {40,41}; //top surf
// Front
Line(42) = {34,38};
Line(43) = {38,37};
Line(44) = {37,33};
Line(45) = {33,34};
Line Loop(46) = {42,43,44,45};
surf4 = news;
Plane Surface(surf4) = {46};
// Left
Line(47) = {33,37};
Line(48) = {37,36};
Line(49) = {36,32};
Line(50) = {32,33};
Line Loop(51) = {47,48,49,50};
surf5 = news;
Plane Surface(surf5) = {51};
// Back
Line(52) = {32,36};
Line(53) = {36,35};
Line(54) = {35,31};
Line(55) = {31,32};
Line Loop(56) = {52,53,54,55};
surf6 = news;
Plane Surface(surf6) = {56};
// Right
Line(57) = {31,35};
Line(58) = {35,38};
Line(59) = {38,34};
Line(60) = {34,31};
Line Loop(61) = {57,58,59,60};
surf7 = news;
Plane Surface(surf7) = {61};
surf8 = news;
// Surface Loop(surf8) = {surf2,surf5,surf4,surf3,surf6,surf7,-9,-13,-17,-21} ;
Surface Loop(surf8) = {surf2,-9,-13,-17,-21,surf3,surf5,surf4,surf7,surf6} ;
Volume(1) = {surf8} ;
//