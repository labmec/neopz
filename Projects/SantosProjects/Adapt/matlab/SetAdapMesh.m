function mdOut = SetAdapMesh(md, meshfile, datafile, parfile)
%Set the new mesh (refined) into the model md

[x,y,elements,segments,segmentmarkers] = ReadNewMesh(md, meshfile);
[surface, base, bed, pressure, temperature, vx, vy, masklevelset] = ReadInitialData(datafile); 
%[edgenodes] = FatherEdgeNodes( elements, fatherindex, x, y, md );

nnodes = length(x);
ndata = length(surface);

if nnodes ~= ndata
    return;
end

NewModel = model;

% geometry
NewModel.mesh = mesh2d();
NewModel.mesh.x = x;
NewModel.mesh.y = y;
NewModel.mesh.elements = elements;
NewModel.mesh.segments = segments;
NewModel.mesh.segmentmarkers = segmentmarkers;

% connectivity
% Fill in rest of fields:
NewModel.mesh.numberofelements = size(NewModel.mesh.elements,1);
NewModel.mesh.numberofvertices = length(NewModel.mesh.x);
NewModel.mesh.vertexonboundary = zeros(NewModel.mesh.numberofvertices,1);
NewModel.mesh.vertexonboundary(NewModel.mesh.segments(:,1:2)) = 1;

% Now, build the connectivity tables for this mesh.
NewModel.mesh.vertexconnectivity = NodeConnectivity(NewModel.mesh.elements,NewModel.mesh.numberofvertices);
NewModel.mesh.elementconnectivity = ElementConnectivity(NewModel.mesh.elements,NewModel.mesh.vertexconnectivity);

% Itapopo testando aqui
NewModel = setmask(NewModel,'','');
NewModel = parameterize(NewModel, parfile);
% NewModel = parameterize(NewModel, './Exp_Par/Mismip.par');


% interpolating data
%LastTime = size(md.results.TransientSolution);    
%LastTime = LastTime(2);

NewModel.geometry.surface = surface;%InterpData( edgenodes, md.results.TransientSolution(LastTime).Surface , md, NewModel); 
NewModel.geometry.base = base;%InterpData( edgenodes, md.results.TransientSolution(LastTime).Base, md, NewModel); 
NewModel.geometry.bed = bed;%InterpData( edgenodes, md.geometry.bed, md, NewModel); 
NewModel.geometry.thickness = NewModel.geometry.surface - NewModel.geometry.base;

% common settings
%disp('!!!!!!!!!!!!!!!     USAR BC.par       !!!!!!!!!!!!!!!' );
% itapopo testando em outro lugar esta chamada
%NewModel = setmask(NewModel,'','');
%NewModel = parameterize(NewModel, './Exp_Par/Mismip.par');% parfle 'BCdebug.par');%itapopo usar BC.par, inclusive no MAIN!!!

% initialization from last md.results time step
NewModel.initialization.pressure = zeros(NewModel.mesh.numberofvertices,1);
NewModel.initialization.vx = zeros(NewModel.mesh.numberofvertices,1);
NewModel.initialization.vy = zeros(NewModel.mesh.numberofvertices,1);
NewModel.initialization.vz = zeros(NewModel.mesh.numberofvertices,1);
NewModel.initialization.vel = zeros(NewModel.mesh.numberofvertices,1);
NewModel.initialization.temperature = zeros(NewModel.mesh.numberofvertices,1);

NewModel.initialization.pressure = pressure;%InterpData( edgenodes, md.results.TransientSolution(LastTime).Pressure, md, NewModel );
NewModel.initialization.vx = vx; %InterpData( edgenodes, md.results.TransientSolution(LastTime).Vx, md, NewModel );
NewModel.initialization.vy = vy; %InterpData( edgenodes, md.results.TransientSolution(LastTime).Vy, md, NewModel );
NewModel.initialization.vel = sqrt(power(vx,2)+power(vy,2));%InterpData( edgenodes, md.results.TransientSolution(LastTime).Vel, md, NewModel );
NewModel.initialization.temperature = temperature;%InterpData( edgenodes, md.initialization.temperature, md, NewModel );

% groundedice_levelset from last md.results time step
%NewModel.mask.groundedice_levelset = InterpData( edgenodes, md.results.TransientSolution(LastTime).MaskGroundediceLevelset, md, NewModel);

%NewModel.mask.groundedice_levelset = NewModel.geometry.thickness + (NewModel.materials.rho_water/NewModel.materials.rho_ice) * NewModel.geometry.bed;
NewModel.mask.groundedice_levelset = masklevelset;

pos = find(NewModel.mask.groundedice_levelset>0.);
pos2 = find(abs( NewModel.geometry.base(pos)-NewModel.geometry.bed(pos) ) > 10^-10 );
posSize = length(pos2);
for i = 1:posSize
    point = pos(pos2(i));
    NewModel.geometry.base(point) = NewModel.geometry.bed(point);%itapopo
    NewModel.mask.groundedice_levelset(point) = 0.;
end
%itapopo
NewModel.geometry.thickness = NewModel.geometry.surface - NewModel.geometry.base;
%NewModel.mask.groundedice_levelset = NewModel.geometry.thickness + (NewModel.materials.rho_water/NewModel.materials.rho_ice) * NewModel.geometry.bed;
%itapopo

mdOut = NewModel;

end

