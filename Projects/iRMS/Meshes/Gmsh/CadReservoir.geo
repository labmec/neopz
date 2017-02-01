
Macro ReservoirCAD

Geometry.Tolerance = 0.01; // Geometrical tolerance
Geometry.OCCSewFaces = 1; // Sew faces in STEP, IGES and BRep models

Merge "Reservoir.igs";
//Merge "YReservoir.igs";

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

line_id = 1;
res_edges_h[] = {
line_id,line_id+1,line_id+2,line_id+3,line_id+5,line_id+8,line_id+9,line_id+11};
res_edges_v[] = {
line_id+4,line_id+6,line_id+7,line_id+10};

Transfinite Line {res_edges_h[]} = 10.0;
Transfinite Line {res_edges_v[]} = 2.0;

reservoir_boundaries[] = {res_B[],res_T[],res_S[],res_E[],res_N[],res_W[]};

Return

