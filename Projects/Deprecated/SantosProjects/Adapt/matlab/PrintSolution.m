function PrintSolution(md, solutionfile)
% print the solution for each node
% the sequence is the same of the mesh nodes

% open a file for writing
fid = fopen(solutionfile, 'w');

% print a title, followed by a blank line
% fprintf(fid, 'solution\n\n');

% nodes data
nnodes = md.mesh.numberofvertices;
if(md.mesh.dimension==3)
    nnodes = md.mesh.numberofvertices2d;
end
% number of nodes
fprintf(fid, '%i\n', nnodes);

% last time with solution
LastTime = size(md.results.TransientSolution);    
LastTime = LastTime(2);

% surface
surface = md.results.TransientSolution(LastTime).Surface;
for i = 1:nnodes
    value = surface(i);
    fprintf(fid, '%.12e\n', value);
end

% base
base = md.results.TransientSolution(LastTime).Base;
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
pressure = md.results.TransientSolution(LastTime).Pressure;
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
vx = md.results.TransientSolution(LastTime).Vx;
for i = 1:nnodes
    value = vx(i);
    fprintf(fid, '%.12e\n', value);
end

% vy
vy = md.results.TransientSolution(LastTime).Vy;
for i = 1:nnodes
    value = vy(i);
    fprintf(fid, '%.12e\n', value);
end

% MaskGroundediceLevelset
maskLevelSet = md.results.TransientSolution(LastTime).MaskGroundediceLevelset;
for i = 1:nnodes
    value = maskLevelSet(i);
    fprintf(fid, '%.12e\n', value);
end


fclose(fid);

end
