function mdOut = SetAdapMeshAndInterpData(md, solutionenum, meshfile, parfile)
%Set the new mesh (refined) into the model md

[x,y,elements,segments,segmentmarkers] = ReadNewMesh(md, meshfile);

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

if ~strncmp(fliplr(EnumToString(solutionenum)),fliplr('Solution'),8),
	error(['solutionenum ' EnumToString(solutionenum) ' not supported!']);
end

if strncmp(EnumToString(solutionenum),'Stressbalance',13),
    
	vx              = md.initialization.vx;
    vy              = md.initialization.vy;
    vz              = md.initialization.vz;
    vel             = md.initialization.vel;
    pressure        = md.initialization.pressure;
    temperature     = md.initialization.temperature;
    surface         = md.geometry.surface;
    base            = md.geometry.base;
    bed             = md.geometry.bed;
    thickness       = md.geometry.thickness;
    masklevelset    = md.mask.groundedice_levelset;
    
elseif strncmp(EnumToString(solutionenum),'Transient',9),

	vx              = md.results.TransientSolution(end).Vx;
	vy              = md.results.TransientSolution(end).Vy;
    vz              = zeros(md.mesh.numberofvertices,1);%md.results.TransientSolution(end).Vz;itapopo
	vel             = md.results.TransientSolution(end).Vel;
    pressure        = md.results.TransientSolution(end).Pressure;
    temperature     = md.initialization.temperature;%md.results.TransientSolution(end).Temperature;itapopo
    surface         = md.results.TransientSolution(end).Surface;
    base            = md.results.TransientSolution(end).Base;
    bed             = md.geometry.bed;
    thickness       = md.results.TransientSolution(end).Thickness;
	masklevelset    = md.results.TransientSolution(end).MaskGroundediceLevelset;
    
else
    error(['solutionenum ' EnumToString(solutionenum) ' not supported!']);
end

% interpolating geometry data
NewModel.geometry.surface           = InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,surface,NewModel.mesh.x,NewModel.mesh.y);
NewModel.geometry.base              = InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,base,NewModel.mesh.x,NewModel.mesh.y); 
NewModel.geometry.bed               = InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,bed,NewModel.mesh.x,NewModel.mesh.y); 
NewModel.geometry.thickness         = InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,thickness,NewModel.mesh.x,NewModel.mesh.y);

% interpolating initializaton data
NewModel.initialization.pressure    = InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,pressure,NewModel.mesh.x,NewModel.mesh.y);
NewModel.initialization.vx          = InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,vx,NewModel.mesh.x,NewModel.mesh.y);
NewModel.initialization.vy          = InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,vy,NewModel.mesh.x,NewModel.mesh.y);
NewModel.initialization.vz          = InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,vz,NewModel.mesh.x,NewModel.mesh.y);
NewModel.initialization.vel         = InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,vel,NewModel.mesh.x,NewModel.mesh.y);
NewModel.initialization.temperature = InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,temperature,NewModel.mesh.x,NewModel.mesh.y);

% interpolation masklevelset data
NewModel.mask.groundedice_levelset = InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,masklevelset,NewModel.mesh.x,NewModel.mesh.y);

% verifying masklevel set
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
%itapopo

%copy other data
NewModel.miscellaneous          = md.miscellaneous;
%NewModel.flowequation          = md.flowequation;
NewModel.timestepping           = md.timestepping;
NewModel.settings               = md.settings;
NewModel.stressbalance.maxiter  = md.stressbalance.maxiter;
NewModel.stressbalance.abstol   = md.stressbalance.abstol;
NewModel.stressbalance.restol   = md.stressbalance.restol;
NewModel.verbose        = md.verbose;
NewModel.cluster        = md.cluster;

mdOut = NewModel;

end
