
Macro DefineOmegas

If (dimension == 3)

// reservoir
sl1 = newsl; Surface Loop(sl1) = {well_p_regions[],well_i_regions[],reservoir_boundaries[]};
v1  = newv; Volume(v1) = {sl1} ;
reservoir[] = {v1};


If (hexahedronsRQ == 1)
Transfinite Surface {reservoir_boundaries[]};
Recombine Surface {reservoir_boundaries[]};
Transfinite Volume {v1};
Recombine Volume {v1};
EndIf

// side-burden
If (geomechanicQ == 1)
sl2 = newsl; Surface Loop(sl2) = {reservoir_boundaries[],sb_boundaries[]};
v2  = newv; Volume(v2) = {sl2} ;
side_burden[] = {v2};

If (hexahedronsSBQ == 1)
Transfinite Surface {sb_boundaries[]};
Recombine Surface {sb_boundaries[]};
Transfinite Volume {v2};
Recombine Volume {v2};
EndIf

EndIf


Else

// reservoir
ll1 = newll; Line Loop(ll1) = {well_p_regions[],well_i_regions[],reservoir_boundaries[]};
s1  = news; Plane Surface(s1) = {ll1};
reservoir[] = {s1};

If (hexahedronsRQ == 1)
Transfinite Surface {s1};
Recombine Surface {s1};
EndIf

// side-burden
If (geomechanicQ == 1)
ll2 = newll; Line Loop(ll2) = {reservoir_boundaries[],sb_boundaries[]};
s2  = news; Plane Surface(s2) = {ll2};
side_burden[] = {s2};

If (hexahedronsSBQ == 1)
Transfinite Surface {s2};
Recombine Surface {s2};
EndIf

EndIf


EndIf

Return

