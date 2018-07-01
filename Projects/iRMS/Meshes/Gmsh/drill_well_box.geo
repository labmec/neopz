
// ---- Gmsh Macro ----
// ---- drill a cylindrical wellbore and ellipsoidal wellbore region  ----
// Created 10/01/2016 by Omar Duran
// Labmec, State University of Campinas
// --------------------------------------------

Macro DrillWellBox

xdir = Cos(beta);
ydir = Sin(beta);

If(dimension == 3)

// mesh controls on wellbore region on 2D
alpha = 1.5;
n_radial = 4;
n_azimuthal = 4;
n_longitudinal = 4;
s=1.5;

////////////////////////////////////////////////////////////////////////////
// 3D model
////////////////////////////////////////////////////////////////////////////

ap0 = newp; Point(ap0) = {0, 0, +0.5*wl, rw_cell};   	// Circle center
ap1 = newp; Point(ap1) = {-rw, -rw, +0.5*wl, rw_cell};
ap2 = newp; Point(ap2) = {+rw, -rw, +0.5*wl, rw_cell};
ap3 = newp; Point(ap3) = {+rw, +rw, +0.5*wl, rw_cell};
ap4 = newp; Point(ap4) = {-rw, +rw, +0.5*wl, rw_cell};

a1 = newp; Point(a1) = {-wbr, -wbr, +s*wl, wr_cell}; 
a2 = newp; Point(a2) = {+wbr, -wbr, +s*wl, wr_cell}; 
a3 = newp; Point(a3) = {+wbr, +wbr, +s*wl, wr_cell}; 
a4 = newp; Point(a4) = {-wbr, +wbr, +s*wl, wr_cell}; 

ap[] = {ap0,ap1,ap2,ap3,ap4,a1,a2,a3,a4};

bp0 = newp; Point(bp0) = {0, 0, -wl+0.5*wl, rw_cell};   	// Circle center
bp1 = newp; Point(bp1) = {-rw, -rw, -wl+0.5*wl, rw_cell};
bp2 = newp; Point(bp2) = {+rw, -rw, -wl+0.5*wl, rw_cell};
bp3 = newp; Point(bp3) = {+rw, +rw, -wl+0.5*wl, rw_cell};
bp4 = newp; Point(bp4) = {-rw, +rw, -wl+0.5*wl, rw_cell};

b1 = newp; Point(b1) = {-wbr, -wbr, -s*wl, wr_cell}; 
b2 = newp; Point(b2) = {+wbr, -wbr, -s*wl, wr_cell}; 
b3 = newp; Point(b3) = {+wbr, +wbr, -s*wl, wr_cell}; 
b4 = newp; Point(b4) = {-wbr, +wbr, -s*wl, wr_cell}; 

bp[] = {bp0,bp1,bp2,bp3,bp4,b1,b2,b3,b4};

// Apply translation
If (dimension == 3)

translate_p[] = Translate {wcx, wcy, wcz} { 
Point {ap[],bp[]};
};

Else

If (xzQ == 1)

translate_p[] = Translate {wcx, 0, wcz} { 
Point {ap[],bp[]};
};

Else

translate_p[] = Translate {wcx, wcy, 0} { 
Point {ap[],bp[]};
};

EndIf

EndIf

// Apply rotation
rotate_p[] = Rotate { { xdir, ydir,0}, {wcx, wcy, wcz}, alfa } {
Point {ap[],bp[]};
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


////////////////////////////////////////////////////////////////////////////
// 3D model:: Lines
////////////////////////////////////////////////////////////////////////////

// south

sca1 = newl; Circle(sca1) = {ap1, ap0, ap2};
sca2 = newl; Circle(sca2) = {ap2, ap0, ap3};
sca3 = newl; Circle(sca3) = {ap3, ap0, ap4};
sca4 = newl; Circle(sca4) = {ap4, ap0, ap1};


sl1 = newl; Line(sl1) = {a1,a2};
sl2 = newl; Line(sl2) = {a2,a3};
sl3 = newl; Line(sl3) = {a3,a4};
sl4 = newl; Line(sl4) = {a4,a1};

sli1 = newl; Line(sli1) = {a1,ap1};
sli2 = newl; Line(sli2) = {a2,ap2};
sli3 = newl; Line(sli3) = {a3,ap3};
sli4 = newl; Line(sli4) = {a4,ap4};


// north

ncb1 = newl; Circle(ncb1) = {bp1, bp0, bp2};
ncb2 = newl; Circle(ncb2) = {bp2, bp0, bp3};
ncb3 = newl; Circle(ncb3) = {bp3, bp0, bp4};
ncb4 = newl; Circle(ncb4) = {bp4, bp0, bp1};

nl1 = newl; Line(nl1) = {b1,b2};
nl2 = newl; Line(nl2) = {b2,b3};
nl3 = newl; Line(nl3) = {b3,b4};
nl4 = newl; Line(nl4) = {b4,b1};

nli1 = newl; Line(nli1) = {b1,bp1};
nli2 = newl; Line(nli2) = {b2,bp2};
nli3 = newl; Line(nli3) = {b3,bp3};
nli4 = newl; Line(nli4) = {b4,bp4};

// laterals
l1 = newl; Line(l1) = {a1,b1};
l2 = newl; Line(l2) = {a2,b2};
l3 = newl; Line(l3) = {a3,b3};
l4 = newl; Line(l4) = {a4,b4};

li1 = newl; Line(li1) = {ap1,bp1};
li2 = newl; Line(li2) = {ap2,bp2};
li3 = newl; Line(li3) = {ap3,bp3};
li4 = newl; Line(li4) = {ap4,bp4};

radial_lines[] = {-nli1,-nli2,-nli3,-nli4,-sli1,-sli2,-sli3,-sli4};
frontal_lines[] = {sl1,sl2,sl3,sl4,sli1,sli2,sli3,sli4,nl1,nl2,nl3,nl4,nli1,nli2,nli3,nli4};
lateral_lines[] = {l1,l2,l3,l4,li1,li2,li3,li4}; 

Transfinite Line {lateral_lines[]} = n_longitudinal;
Transfinite Line {frontal_lines[]} = n_azimuthal;
Transfinite Line {radial_lines[]} = n_radial Using Progression alpha;

////////////////////////////////////////////////////////////////////////////
// 3D model:: Surfaces
////////////////////////////////////////////////////////////////////////////

// south

sal0 = newll; Line Loop(sal0) = {sca1, sca2, sca3, sca4}; sa0 = news; Plane Surface (sa0) = {sal0};

sal1 = newll; Line Loop(sal1) = {sl1, sli2, -sca1,-sli1}; sa1 = news; Ruled Surface (sa1) = {sal1};
sal2 = newll; Line Loop(sal2) = {sl2, sli3, -sca2,-sli2}; sa2 = news; Ruled Surface (sa2) = {sal2};
sal3 = newll; Line Loop(sal3) = {sl3, sli4, -sca3,-sli3}; sa3 = news; Ruled Surface (sa3) = {sal3};
sal4 = newll; Line Loop(sal4) = {sl4, sli1, -sca4,-sli4}; sa4 = news; Ruled Surface (sa4) = {sal4};

sal5 = newll; Line Loop(sal5) = {sl1, sl2, sl3, sl4}; sa5 = news; Plane Surface (sa5) = {sal5};


// north

nal0 = newll; Line Loop(nal0) = {ncb1, ncb2, ncb3, ncb4}; na0 = news; Plane Surface (na0) = {nal0};

nal1 = newll; Line Loop(nal1) = {nl1, nli2, -ncb1,-nli1}; na1 = news; Ruled Surface (na1) = {nal1};
nal2 = newll; Line Loop(nal2) = {nl2, nli3, -ncb2,-nli2}; na2 = news; Ruled Surface (na2) = {nal2};
nal3 = newll; Line Loop(nal3) = {nl3, nli4, -ncb3,-nli3}; na3 = news; Ruled Surface (na3) = {nal3};
nal4 = newll; Line Loop(nal4) = {nl4, nli1, -ncb4,-nli4}; na4 = news; Ruled Surface (na4) = {nal4};

nal5 = newll; Line Loop(nal5) = {nl1, nl2, nl3, nl4}; na5 = news; Plane Surface (na5) = {nal5};

// laterals

ell1 = newll; Line Loop(ell1) = {nl2, -l3, -sl2, l2}; se1 = news; Plane Surface (se1) = {ell1};
tll2 = newll; Line Loop(tll2) = {nl3, -l4, -sl3, l3}; st2 = news; Plane Surface (st2) = {tll2};
wll3 = newll; Line Loop(wll3) = {nl4, -l1, -sl4, l4}; sw3 = news; Plane Surface (sw3) = {wll3};
bll4 = newll; Line Loop(bll4) = {nl1, -l2, -sl1, l1}; sb4 = news; Plane Surface (sb4) = {bll4};

eill1 = newll; Line Loop(eill1) = {ncb2, -li3, -sca2, li2}; sei1 = news; Ruled Surface (sei1) = {eill1};
till2 = newll; Line Loop(till2) = {ncb3, -li4, -sca3, li3}; sti2 = news; Ruled Surface (sti2) = {till2};
will3 = newll; Line Loop(will3) = {ncb4, -li1, -sca4, li4}; swi3 = news; Ruled Surface (swi3) = {will3};
bill4 = newll; Line Loop(bill4) = {ncb1, -li2, -sca1, li1}; sbi4 = news; Ruled Surface (sbi4) = {bill4};

////////////////////////////////////////////////////////////////////////////
// 3D model:: grouping
////////////////////////////////////////////////////////////////////////////

internal_surfaces[] = {na1,na2,na3,na4,sa1,sa2,sa3,sa4};
well_lid[] = {sa0,na0};
well_bore[] = {sei1,sti2,swi3,sbi4};
well_region[] = {sa5,se1,st2,sw3,sb4,na5};

// wellbore region volume
wr_sl1 = newsl; Surface Loop(wr_sl1) = {well_lid[],well_bore[],well_region[]}; 
wr_v1 = newv; Volume(wr_v1) = {wr_sl1};
well_v_region[] = {wr_v1};

Transfinite Surface {well_region[],well_bore[],well_lid[],internal_surfaces[]};
If (hexahedronsWQ == 1)
Recombine Surface {well_region[],well_bore[],well_lid[],internal_surfaces[]};
Recombine Volume {well_v_region[]};
EndIf

Else

////////////////////////////////////////////////////////////////////////////
// 2D model
////////////////////////////////////////////////////////////////////////////

// mesh controls on wellbore region on 2D
alpha = 1.5;
n_radial = 10;
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
M = #well_lids[];
For i In {0:N-1}
	well_lids[i + M] = well_lid[i];
EndFor

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
