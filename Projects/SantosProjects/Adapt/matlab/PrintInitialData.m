function PrintInitialData(md, solutionfile)
% print the solution for each node
% the sequence is the same of the mesh nodes

% open a file for writing
fid = fopen(solutionfile, 'w');

% print a title, followed by a blank line
% fprintf(fid, 'solution\n\n');

% nodes data
nnodes = md.mesh.numberofvertices;

% number of nodes
fprintf(fid, '%i\n', nnodes);

% surface
surface = md.geometry.surface;
for i = 1:nnodes
    value = surface(i);
    fprintf(fid, '%.12e\n', value);
end

% base
base = md.geometry.base;
for i = 1:nnodes
    value = base(i);
    fprintf(fid, '%.12e\n', value);
end

% bed
bed = md.geometry.bed;
for i = 1:nnodes
    value = bed(i);
    fprintf(fid, '%.12e\n', value);
end

% pressure
pressure = md.initialization.pressure;
for i = 1:nnodes
    value = pressure(i);
    fprintf(fid, '%.12e\n', value);
end

% temperature
temperature = md.initialization.temperature;
for i = 1:nnodes
    value = temperature(i);
    fprintf(fid, '%.12e\n', value);
end

% vx
vx = md.initialization.vx;
for i = 1:nnodes
    value = vx(i);
    fprintf(fid, '%.12e\n', value);
end

% vy
vy = md.initialization.vy;
for i = 1:nnodes
    value = vy(i);
    fprintf(fid, '%.12e\n', value);
end

% MaskGroundediceLevelset
maskLevelSet = md.mask.groundedice_levelset;
for i = 1:nnodes
    value = maskLevelSet(i);
    fprintf(fid, '%.12e\n', value);
end


fclose(fid);

end
