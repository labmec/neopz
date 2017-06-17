
YShapeReservoirQ = 0;

Macro ReservoirCAD

//Geometry.Tolerance = 0.01; // Geometrical tolerance
Geometry.OCCSewFaces = 1; // Sew faces in STEP, IGES and BRep models

If(dimension == 3)

If(YShapeReservoirQ == 1)

Merge "YReservoir.igs";

s1  = 1; // Top unstructured region
s9  = 9; // Bottom unstructured region
s8  = 8; // West unstructured region
s6  = 6; // East unstructured region
s7  = 7; // North unstructured region
s2[]  = {2,3,4,5}; // South unstructured region

res_T[] = {s1};
res_B[] = {s9};
res_W[] = {s8};
res_E[] = {s6};
res_N[] = {s7};
res_S[] = {s2[]};


line_id = 1;
res_edges_h[] = {1,2,3,4,5,6,7,9,12,14,16,18,20,21};
res_edges_v[] = {8,10,11,13,15,17,19};

Transfinite Line {res_edges_h[]} = 15.0;
Transfinite Line {res_edges_v[]} = 5.0;

reservoir_boundaries[] = {res_B[],res_T[],res_S[],res_E[],res_N[],res_W[]};

Else

Merge "Reservoir.igs";

s1  = 1; // Bottom unstructured region
s2  = 2; // Top unstructured region
s3  = 3; // South unstructured region
s4  = 4; // East unstructured region
s5  = 5; // North unstructured region
s6  = 6; // West unstructured region

res_B[] = {s1};
res_W[] = {s2};
res_N[] = {s3};
res_E[] = {s4};
res_S[] = {s5};
res_T[] = {s6};

res_edges_h[] = {1,2,3,4,6,9,10,12};
res_edges_v[] = {5,7,8,11};

Transfinite Line {res_edges_h[]} = 25.0;
Transfinite Line {res_edges_v[]} = 2.0;

reservoir_boundaries[] = {res_B[],res_T[],res_S[],res_E[],res_N[],res_W[]};

EndIf


Else

Include "Reservoir_2D.geo";

s1  = 1; // North unstructured region
s2  = 2; // South unstructured region
s3  = 3; // West unstructured region
s4  = 4; // East unstructured region

res_N[] = {s1};
res_S[] = {s2};
res_E[] = {s3};
res_W[] = {s4};

res_edges_h[] = {res_S[],res_N[]};
res_edges_v[] = {res_W[],res_E[]};

Transfinite Line {res_edges_h[]} = 35.0;
Transfinite Line {res_edges_v[]} = 2.0;

reservoir_boundaries[] = {res_N[],res_E[],-res_S[],-res_W[]};

EndIf


Return


