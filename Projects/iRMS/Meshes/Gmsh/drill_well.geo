
// ---- Gmsh Macro ----
// ---- drill a cylindrical wellbore and ellipsoidal wellbore region  ----
// Created 10/01/2016 by Omar Duran
// Labmec, State University of Campinas
// --------------------------------------------

Macro DrillWell

xdir = Cos(beta);
ydir = Sin(beta);

If(dimension == 3)

////////////////////////////////////////////////////////////////////////////
// 3D model
////////////////////////////////////////////////////////////////////////////

ap1 = newp; Point(ap1) = {0, 0, +0.5*wl, rw_cell};   	// Sphere center
ap2 = newp; Point(ap2) = {rw, 0, +0.5*wl, rw_cell};  //Sphere X radius
ap3 = newp; Point(ap3) = {0, rw, +0.5*wl, rw_cell};  //Sphere y radius
ap4 = newp; Point(ap4) = {0, 0, rw+0.5*wl, rw_cell};  //Sphere z radius
ap[] = {ap1,ap2,ap3,ap4};

bp1 = newp; Point(bp1) = {0, 0, -wl+0.5*wl, rw_cell};   	// Sphere center
bp2 = newp; Point(bp2) = {rw, 0, -wl+0.5*wl, rw_cell};  //Sphere X radius
bp3 = newp; Point(bp3) = {0, rw, -wl+0.5*wl, rw_cell};  //Sphere y radius
bp4 = newp; Point(bp4) = {0, 0, -rw-wl+0.5*wl, rw_cell};  //Sphere z radius
bp[] = {bp1,bp2,bp3,bp4};


wr_pc = newp; Point(wr_pc) = {0, 0, 0, wr_cell};   	// Sphere center
wr_pc1 = newp; Point(wr_pc1) = {wbr, 0, 0, wr_cell};  //Sphere X radius
wr_pc2 = newp; Point(wr_pc2) = {0, wbr, 0, wr_cell};  //Sphere y radius
wr_ap1 = newp; Point(wr_ap1) = {0, 0, ela, wr_cell};  //Sphere z radius
wr_bp1 = newp; Point(wr_bp1) = {0, 0, -ela, wr_cell};  //Sphere z radius
wr_p[] = {wr_pc,wr_pc1,wr_pc2,wr_ap1,wr_bp1};


// Apply translation
If (dimension == 3)

translate_p[] = Translate {wcx, wcy, wcz} { 
Point {ap[],bp[],wr_p[]};
};

Else

If (xzQ == 1)

translate_p[] = Translate {wcx, 0, wcz} { 
Point {ap[],bp[],wr_p[]};
};

Else

translate_p[] = Translate {wcx, wcy, 0} { 
Point {ap[],bp[],wr_p[]};
};

EndIf

EndIf

// Apply rotation
rotate_p[] = Rotate { { xdir, ydir,0}, {wcx, wcy, wcz}, alfa } {
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

lidst1[] = Rotate {{axis[0],axis[1],axis[2]},{wcx,wcy,wcz},Pi/2} {Duplicata{Surface{as1,bs1};}}; 
lidst2[] = Rotate {{axis[0],axis[1],axis[2]},{wcx,wcy,wcz},Pi} {Duplicata{Surface{as1,bs1};}}; 
lidst3[] = Rotate {{axis[0],axis[1],axis[2]},{wcx,wcy,wcz},3*Pi/2} {Duplicata{Surface{as1,bs1};}};

wellbore_t1[] = Rotate {{axis[0],axis[1],axis[2]},{wcx,wcy,wcz},Pi/2} {Duplicata{Surface{cs1};}}; 
wellbore_t2[] = Rotate {{axis[0],axis[1],axis[2]},{wcx,wcy,wcz},Pi} {Duplicata{Surface{cs1};}}; 
wellbore_t3[] = Rotate {{axis[0],axis[1],axis[2]},{wcx,wcy,wcz},3*Pi/2} {Duplicata{Surface{cs1};}};

wellregion_t1[] = Rotate {{axis[0],axis[1],axis[2]},{wcx,wcy,wcz},Pi/2} {Duplicata{Surface{wr_as1,wr_bs1};}}; 
wellregion_t2[] = Rotate {{axis[0],axis[1],axis[2]},{wcx,wcy,wcz},Pi} {Duplicata{Surface{wr_as1,wr_bs1};}}; 
wellregion_t3[] = Rotate {{axis[0],axis[1],axis[2]},{wcx,wcy,wcz},3*Pi/2} {Duplicata{Surface{wr_as1,wr_bs1};}};


well_lid[] = {as1,bs1,lidst1[],lidst2[],lidst3[]};
well_bore[] = {cs1,wellbore_t1[],wellbore_t2[],wellbore_t3[]};
well_region[] = {wr_as1,wr_bs1,wellregion_t1[],wellregion_t2[],wellregion_t3[]};

// wellbore region volume
wr_sl1 = newsl; Surface Loop(wr_sl1) = {well_lid[],well_bore[],well_region[]}; 
wr_v1 = newv; Volume(wr_v1) = {wr_sl1};
well_v_region[] = {wr_v1};


Else

////////////////////////////////////////////////////////////////////////////
// 2D model
////////////////////////////////////////////////////////////////////////////

// mesh controls on wellbore region on 2D
alpha = 2.0;
n_radial = 5;
n_azimuthal = 3;

apc = newp; Point(apc) = {0, 0, 0, rw_cell};   	// Sphere center

If (xzQ == 1)

ap1 = newp; Point(ap1) = {rw, 0, 0, rw_cell};  //Sphere X radius
ap2 = newp; Point(ap2) = {0, 0, rw, rw_cell};  //Sphere z radius

bp1 = newp; Point(bp1) = {wbr, 0, 0, wr_cell};  //Sphere X radius
bp2 = newp; Point(bp2) = {0, 0, wbr, wr_cell};  //Sphere z radius

Else

ap1 = newp; Point(ap1) = {rw, 0, 0, rw_cell};  //Sphere X radius
ap2 = newp;Point(ap2) = {0, rw, 0, rw_cell};  //Sphere z radius

bp1 = newp; Point(bp1) = {wbr, 0, 0, wr_cell};  //Sphere X radius
bp2 = newp; Point(bp2) = {0, wbr, 0, wr_cell};  //Sphere z radius

EndIf

ap[] = {ap1,ap2,apc};
bp[] = {bp1,bp2};

If (xzQ == 1)
wcy = 0;
Else
wcz = 0;
EndIf

translate_p[] = Translate {wcx, wcy, wcz} { 
Point {ap[],bp[]};
};

al1 = newl; Circle(al1) = {ap1, apc, ap2};    //Sphere XY arc
bl1 = newl; Circle(bl1) = {bp1, apc, bp2};    //Sphere XY arc
cl1 = newl; Line(cl1) = {ap2, bp2};    //Line XY arc

If (xzQ == 1)
axis[] = {0, -1, 0};

Else
axis[] = {0, 0, 1};

EndIf


wellbore_t1[] = Rotate {{axis[0],axis[1],axis[2]},{wcx,wcy,wcz},Pi/2} {Duplicata{Line{al1};}};
wellbore_t2[] = Rotate {{axis[0],axis[1],axis[2]},{wcx,wcy,wcz},Pi} {Duplicata{Line{al1};}};
wellbore_t3[] = Rotate {{axis[0],axis[1],axis[2]},{wcx,wcy,wcz},3*Pi/2} {Duplicata{Line{al1};}};

wellregion_t1[] = Rotate {{axis[0],axis[1],axis[2]},{wcx,wcy,wcz},Pi/2} {Duplicata{Line{bl1};}};
wellregion_t2[] = Rotate {{axis[0],axis[1],axis[2]},{wcx,wcy,wcz},Pi} {Duplicata{Line{bl1};}};
wellregion_t3[] = Rotate {{axis[0],axis[1],axis[2]},{wcx,wcy,wcz},3*Pi/2} {Duplicata{Line{bl1};}};

radial_t1[] = Rotate {{axis[0],axis[1],axis[2]},{wcx,wcy,wcz},Pi/2} {Duplicata{Line{cl1};}};
radial_t2[] = Rotate {{axis[0],axis[1],axis[2]},{wcx,wcy,wcz},Pi} {Duplicata{Line{cl1};}};
radial_t3[] = Rotate {{axis[0],axis[1],axis[2]},{wcx,wcy,wcz},3*Pi/2} {Duplicata{Line{cl1};}};


well_bore[] = {al1,wellbore_t1[0],wellbore_t2[0],wellbore_t3[0]};
well_region[] = {bl1,wellregion_t1[0],wellregion_t2[0],wellregion_t3[0]};
radial[] = {cl1,radial_t1[0],radial_t2[0],radial_t3[0]};
well_lid[] = {};



// wellbore region volume
wr_ll1 = newll; Line Loop(wr_ll1) = {-al1,radial_t3[0],bl1,-cl1}; 
wr_ll2 = newll; Line Loop(wr_ll2) = {-wellbore_t1[0],cl1,wellregion_t1[0],-radial_t1[0]}; 
wr_ll3 = newll; Line Loop(wr_ll3) = {-wellbore_t2[0],radial_t1[0],wellregion_t2[0],-radial_t2[0]}; 
wr_ll4 = newll; Line Loop(wr_ll4) = {-wellbore_t3[0],radial_t2[0],wellregion_t3[0],-radial_t3[0]}; 
wr_s1 = news; Plane Surface(wr_s1) = {wr_ll1};
wr_s2 = news; Plane Surface(wr_s2) = {wr_ll2};
wr_s3 = news; Plane Surface(wr_s3) = {wr_ll3};
wr_s4 = news; Plane Surface(wr_s4) = {wr_ll4};
well_v_region[] = {wr_s1,wr_s2,wr_s3,wr_s4};

If (hexahedronsWQ == 1)
Transfinite Surface {well_v_region[]};
Recombine Surface {well_v_region[]};
EndIf

Transfinite Line {well_bore[],well_region[]} = n_azimuthal;
Transfinite Line {radial[]} = n_radial Using Progression alpha;

EndIf



////////////////////////////////////////////////////////////////////////////
// Transfering indexes
////////////////////////////////////////////////////////////////////////////

N = #well_lid[];
If(IsInjectorQ == 1)
M = #well_i_lids[];
For i In {0:N-1}
	well_i_lids[i + M] = well_lid[i];
EndFor
Else
M = #well_p_lids[];
For i In {0:N-1}
	well_p_lids[i + M] = well_lid[i];
EndFor
EndIf

N = #well_bore[];
If(IsInjectorQ == 1)
M = #well_i_regions[];
For i In {0:N-1}
	well_i_bores[i + M] = well_bore[i];
EndFor
Else
M = #well_p_regions[];
For i In {0:N-1}
	well_p_bores[i + M] = well_bore[i];
EndFor
EndIf

N = #well_region[];
If(IsInjectorQ == 1)
M = #well_i_regions[];
For i In {0:N-1}
	well_i_regions[i + M] = well_region[i];
EndFor
Else
M = #well_p_regions[];
For i In {0:N-1}
	well_p_regions[i + M] = well_region[i];
EndFor
EndIf

N = #well_v_region[];
If(IsInjectorQ == 1)
M = #well_i_v_regions[];
For i In {0:N-1}
	well_i_v_regions[i + M] = well_v_region[i];
EndFor
Else
M = #well_p_v_regions[];
For i In {0:N-1}
	well_p_v_regions[i + M] = well_v_region[i];
EndFor
EndIf

well_index = well_index + 1;

If(IsInjectorQ == 1)
Printf ("Mesher:: Drilled injector  well = %g, at position {%g,%g,%g}", well_index, wcx,wcy,wcz);
Else
Printf ("Mesher:: Drilled productor well = %g, at position {%g,%g,%g}", well_index, wcx,wcy,wcz);
EndIf

Return
