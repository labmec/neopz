

alfa = 0.2;
beta = 0.5;

wcx = 0.0;
wcy = 0.0;
wcz = 0.0;

wl = 1.0;

rw = 0.1;
rw_cell= 0.1;

LensA = 2.0;
LensB = 4.0; 
wr_cell= 2.0;

xdir = Cos(beta);
ydir = Sin(beta);
zdir = Cos(alfa);


ap1 = newp; Point(ap1) = {0, 0, 0, rw_cell};   	// Sphere center
ap2 = newp; Point(ap2) = {rw, 0, 0, rw_cell};  //Sphere X radius
ap3 = newp; Point(ap3) = {0, rw, 0, rw_cell};  //Sphere y radius
ap4 = newp; Point(ap4) = {0, 0, rw, rw_cell};  //Sphere z radius
ap[] = {ap1,ap2,ap3,ap4};

bp1 = newp; Point(bp1) = {0, 0, -wl, rw_cell};   	// Sphere center
bp2 = newp; Point(bp2) = {rw, 0, -wl, rw_cell};  //Sphere X radius
bp3 = newp; Point(bp3) = {0, rw, -wl, rw_cell};  //Sphere y radius
bp4 = newp; Point(bp4) = {0, 0, -rw-wl, rw_cell};  //Sphere z radius
bp[] = {bp1,bp2,bp3,bp4};


wr_pc = newp; Point(wr_pc) = {0, 0, -0.5*wl, wr_cell};   	// Sphere center
wr_pc1 = newp; Point(wr_pc1) = {LensA, 0, -0.5*wl, wr_cell};  //Sphere X radius
wr_pc2 = newp; Point(wr_pc2) = {0, LensA, -0.5*wl, wr_cell};  //Sphere y radius
wr_ap1 = newp; Point(wr_ap1) = {0, 0, LensB-0.5*wl, wr_cell};  //Sphere z radius
wr_bp1 = newp; Point(wr_bp1) = {0, 0, -LensB-0.5*wl, wr_cell};  //Sphere z radius
wr_p[] = {wr_pc,wr_pc1,wr_pc2,wr_ap1,wr_bp1};


// Apply rotation
rotate_p[] = Rotate { { xdir, ydir,0}, {0,0,0}, alfa } {
Point {ap[],bp[],wr_p[]};
};


c1[] = Point{ap1};
c2[] = Point{ap2};
c3[] = Point{ap3};

u[] = {}; v[] = {};
For i In {0:3-1}
	u[i] =  c2[i] - c1[i];
	v[i] =  c3[i] - c1[i];
EndFor

a[] = {-u[2]*v[1]+u[1]*v[2],u[2]*v[0]-u[0]*v[2],-u[1]*v[0]+u[0]*v[1]};
norm = Sqrt( a[0]*a[0] + a[1]*a[1] + a[2]*a[2] );


al1 = newl; Circle(al1) = {ap2, ap1, ap3};    //Sphere XY arc
al2 = newl; Ellipse(al2) = {ap2, ap1, ap4, ap4};   //Sphere XZ arc
al3 = newl; Ellipse(al3) = {ap3, ap1, ap4, ap4};   //Sphere YZ arc
all1 = newll; Line Loop(all1) = {al1, al3, -al2} ;  //Sphere slice edge 
as1 = news; Ruled Surface (as1) = {all1};    //Sphere slice surface

al[] = {al1,al2,al3};
as[] = {as1};

bl1 = newl; Circle(bl1) = {bp2, bp1, bp3};    //Sphere XY arc
bl2 = newl; Ellipse(bl2) = {bp2, bp1, bp4, bp4};   //Sphere XZ arc
bl3 = newl; Ellipse(bl3) = {bp3, bp1, bp4, bp4};   //Sphere YZ arc
bll1 = newll; Line Loop(bll1) = {bl1, bl3, -bl2} ;  //Sphere slice edge 
bs1 = news; Ruled Surface (bs1) = {bll1};    //Sphere slice surface

bl[] = {bl1,bl2,bl3};
bs[] = {bs1};


cl1 = newl; Line(cl1) = {ap2,bp2};
cl2 = newl; Line(cl2) = {ap3,bp3};
cll1 = newll; Line Loop(cll1) = {al1,cl2,-bl1,-cl1} ;  //Sphere slice edge 
cs1 = news; Ruled Surface (cs1) = {cll1};    //Sphere slice surface

cl[] = {cl1,cl2};
cs[] = {cs1};


wr_alc = newl; Circle(wr_alc) = {wr_pc1, wr_pc, wr_pc2};    //Sphere XY arc

wr_al1 = newl; Ellipse(wr_al1) = {wr_pc1, wr_pc, wr_ap1, wr_ap1};   //Sphere XZ arc
wr_al2 = newl; Ellipse(wr_al2) = {wr_pc2, wr_pc, wr_ap1, wr_ap1};   //Sphere YZ arc
wr_all1 = newll; Line Loop(wr_all1) = {wr_alc, wr_al2, -wr_al1} ;  //Sphere slice edge 
wr_as1 = news; Ruled Surface (wr_as1) = {wr_all1};    //Sphere slice surface

wr_bl1 = newl; Ellipse(wr_bl1) = {wr_pc1, wr_pc, wr_bp1, wr_bp1};   //Sphere XZ arc
wr_bl2 = newl; Ellipse(wr_bl2) = {wr_pc2, wr_pc, wr_bp1, wr_bp1};   //Sphere YZ arc
wr_bll1 = newll; Line Loop(wr_bll1) = {wr_alc, wr_bl2, -wr_bl1} ;  //Sphere slice edge 
wr_bs1 = news; Ruled Surface (wr_bs1) = {wr_bll1};    //Sphere slice surface


axis[] = {a[0]/norm, a[1]/norm, a[2]/norm};

lidst1[] = Rotate {{axis[0],axis[1],axis[2]},{0,0,0},Pi/2} {Duplicata{Surface{as1,bs1};}}; 
lidst2[] = Rotate {{axis[0],axis[1],axis[2]},{0,0,0},Pi} {Duplicata{Surface{as1,bs1};}}; 
lidst3[] = Rotate {{axis[0],axis[1],axis[2]},{0,0,0},3*Pi/2} {Duplicata{Surface{as1,bs1};}};


wellbore_t1[] = Rotate {{axis[0],axis[1],axis[2]},{0,0,0},Pi/2} {Duplicata{Surface{cs1};}}; 
wellbore_t2[] = Rotate {{axis[0],axis[1],axis[2]},{0,0,0},Pi} {Duplicata{Surface{cs1};}}; 
wellbore_t3[] = Rotate {{axis[0],axis[1],axis[2]},{0,0,0},3*Pi/2} {Duplicata{Surface{cs1};}};


wellregion_t1[] = Rotate {{axis[0],axis[1],axis[2]},{0,0,0},Pi/2} {Duplicata{Surface{wr_as1,wr_bs1};}}; 
wellregion_t2[] = Rotate {{axis[0],axis[1],axis[2]},{0,0,0},Pi} {Duplicata{Surface{wr_as1,wr_bs1};}}; 
wellregion_t3[] = Rotate {{axis[0],axis[1],axis[2]},{0,0,0},3*Pi/2} {Duplicata{Surface{wr_as1,wr_bs1};}};





