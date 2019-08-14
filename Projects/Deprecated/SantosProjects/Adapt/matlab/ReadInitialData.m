function [ surface, base, bed, pressure, temperature, vx, vy, masklevelset ] = ReadInitialData( datafile )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% open a file for reading
fid = fopen(datafile, 'r');

% read the number of points
NPoints = fscanf(fid, '%d\n',1);

%read surface
surface = zeros(NPoints, 1);

for i = 1:NPoints
    
    value = fscanf(fid,'%e\t',1);
    
    surface(i) = value;
  
end

%read base
base = zeros(NPoints, 1);

for i = 1:NPoints
    
    value = fscanf(fid,'%e\t',1);
    
    base(i) = value;
  
end

%read bed
bed = zeros(NPoints, 1);

for i = 1:NPoints
    
    value = fscanf(fid,'%e\t',1);
    
    bed(i) = value;
  
end

%read pressure
pressure = zeros(NPoints, 1);

for i = 1:NPoints
    
    value = fscanf(fid,'%e\t',1);
    
    pressure(i) = value;
  
end

%read temperature
temperature = zeros(NPoints, 1);

for i = 1:NPoints
    
    value = fscanf(fid,'%e\t',1);
    
    temperature(i) = value;
  
end

%read vx
vx = zeros(NPoints, 1);

for i = 1:NPoints
    
    value = fscanf(fid,'%e\t',1);
    
    vx(i) = value;
  
end

%read vy
vy = zeros(NPoints, 1);

for i = 1:NPoints
    
    value = fscanf(fid,'%e\t',1);
    
    vy(i) = value;
  
end

%read masklevelset
masklevelset = zeros(NPoints, 1);

for i = 1:NPoints
    
    value = fscanf(fid,'%e\t',1);
    
    masklevelset(i) = value;
  
end

fclose(fid);

end

