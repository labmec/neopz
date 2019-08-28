% main to print the first mesh and to run with new adpated meshes

%initial mesh (mesh0)
if isRun == 0
    
    % set the model
    md = model;
    md.miscellaneous.name = 'Mismip2D';
    %md=triangle(md, 'SquareDomain.exp',75000);%500000 para 4 elementos
    md=triangle(md, expfile, resolution);%500000 para 4 elementos
    
    % initial geometry
    md.geometry.thickness = 100.*ones(md.mesh.numberofvertices,1); 
    md.geometry.base      = -100. - (1./1000.) * abs(md.mesh.x);
    md.geometry.bed       = md.geometry.base;
    md.geometry.surface   = md.geometry.bed + md.geometry.thickness;

    %md=parameterize(md,'BCdebug.par'); 
    md=parameterize(md, parfile);
    
    %initial velocity
    md.initialization.vx = ones(md.mesh.numberofvertices,1);
    md.initialization.vy = ones(md.mesh.numberofvertices,1);
    md.initialization.vz = ones(md.mesh.numberofvertices,1);
    md.initialization.vel = sqrt(2)*ones(md.mesh.numberofvertices,1);
    md.initialization.pressure = md.constants.g*md.materials.rho_ice*md.geometry.thickness;

    %set the number of proccessor on my host
    %md.cluster.np = 4;

    %set the cluster LabMec-Mac-Pro.local
    md.cluster=generic();
    md.cluster.np = 8;
    md.cluster.name = 'LabMec-Mac-Pro.local';
    md.cluster.login = 'santos';

    % first solve
    md.timestepping.start_time = 0.;
    md.timestepping.final_time = 1; %1
    md.timestepping.time_step = 0.125; %1;
    md.settings.output_frequency = 1;
    md=solve(md,TransientSolutionEnum);

    % print mesh0 and its solution
    %meshname = 'mesh0.txt';
    PrintMesh4Adapt(md, meshfile);
    PrintSolution(md, solutionfile);

    % print results to verify
    LastTime = size(md.results.TransientSolution);    
    LastTime = LastTime(2);
    %plotmodel(md,'data',md.results.TransientSolution(LastTime).MaskGroundediceLevelset);

else

    %set the new mesh
    md = SetAdapMesh(md, meshfile, datafile, parfile);
    
    %plotmodel(md,'data',md.mask.groundedice_levelset);
    
    %set the number of proccessor on my host
    %md.cluster.np = 4;
    %set the cluster LabMec-Mac-Pro.local
    md.cluster=generic();
    md.cluster.np = 10;
    md.cluster.name = 'LabMec-Mac-Pro.local';
    md.cluster.login = 'santos';

    % solve
    md.timestepping.start_time = 0; %md.timestepping.final_time;
    md.timestepping.final_time = 5000; %50 p 2x // md.timestepping.start_time + 500;
    md.timestepping.time_step = 1.;%1.; %0.25 p 2x
    md.settings.output_frequency = 1000;
    md = solve(md,TransientSolutionEnum);
    
    % print results to verify
    LastTime = size(md.results.TransientSolution);    
    LastTime = LastTime(2);
    %plotmodel(md,'data','BC');
    %plotmodel(md,'data',md.results.TransientSolution(LastTime).MaskGroundediceLevelset);
    
    %print the solution
    PrintSolution(md, solutionfile);
    
end
