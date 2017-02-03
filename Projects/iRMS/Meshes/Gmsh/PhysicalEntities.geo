
Macro DrawBoundaries

If (dimension == 3)

// Tagging boundary conditions for prodution wells
Physical Surface("well_lids") = well_lids[];

Physical Surface("producers") = well_p_bores[];

Physical Surface("injectors") = well_i_bores[];


Physical Volume("Reservoir") = {v1,well_p_v_regions[],well_i_v_regions[]};
Physical Surface("Reservoir_bottom") = {res_B[]};
Physical Surface("Reservoir_top") = {res_T[]};
Physical Surface("Reservoir_South") = {res_S[]};
Physical Surface("Reservoir_East") = {res_E[]};
Physical Surface("Reservoir_North") = {res_N[]};
Physical Surface("Reservoir_West") = {res_W[]};


// side-burden
If (geomechanicQ == 1)
Physical Volume("side_burden") = {v2};
Physical Surface("side_burden_bottom") = {sb_B[]};
Physical Surface("side_burden_top") = {sb_T[]};
Physical Surface("side_burden_South") = {sb_S[]};
Physical Surface("side_burden_East") = {sb_E[]};
Physical Surface("side_burden_North") = {sb_N[]};
Physical Surface("side_burden_West") = {sb_W[]};
EndIf



Else

// Tagging boundary conditions for prodution wells
Physical Line("well_lids") = well_lids[];

Physical Line("producers") = well_p_bores[];

Physical Line("injectors") = well_i_bores[];

Physical Surface("Reservoir") = {s1,well_p_v_regions[],well_i_v_regions[]};
Physical Line("Reservoir_south") = {res_S[]};
Physical Line("Reservoir_West") = {res_E[]};
Physical Line("Reservoir_north") = {res_N[]};
Physical Line("Reservoir_East") = {res_W[]};

// side-burden
If (geomechanicQ == 1)
Physical Surface("side_burden") = {s2};
Physical Line("side_burden_south") = {sb_S[]};
Physical Line("side_burden_West") = {sb_E[]};
Physical Line("side_burden_north") = {sb_N[]};
Physical Line("side_burden_East") = {sb_W[]};
EndIf

EndIf

Return

