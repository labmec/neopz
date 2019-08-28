function [x, y, elements, segments, segmentmarkers ] = ReadNewMesh(meshfile)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% open a file for reading
fid = fopen(meshfile, 'r');

% read the number of points
NPoints = fscanf(fid, '%d\n',1);

%read the x and y values
fX = zeros(NPoints, 1);
fY = zeros(NPoints, 1);

for i = 1:NPoints
    xv = fscanf(fid,'%e\t',1);
    yv = fscanf(fid,'%e\n',1);
    
    fX(i)=xv;
    fY(i)=yv;    
end

%read the number of elements
NElements = fscanf(fid, '%d\n',1);

%read the elements
fElements = zeros(NElements,3);

for i = 1:NElements
    N1 = fscanf(fid, '%d\t', 1);
    N2 = fscanf(fid, '%d\t', 1);
    N3 = fscanf(fid, '%d\n', 1);
    
    fElements(i,:) = [N1 N2 N3]; 
end

%read the number of segments
NSegments = fscanf(fid, '%d/n',1);

%read the segmentmarkers and segments
fSegmentmarkers = zeros(NSegments,1);
fSegments = zeros(NSegments, 3);

for i = 1:NSegments
    N1 = fscanf(fid, '%d/t',1);
    N2 = fscanf(fid, '%d/t',1);
    ID = fscanf(fid, '%d/n',1);
    
    fSegmentmarkers(i) = i;
    fSegments(i,:) = [N1 N2 ID];

end



%Set the x and y coords
x = fX;
y = fY;

%Set the elements
elements = fElements;

%Set the segments and segmentmarkers
segments = fSegments;
segmentmarkers = fSegmentmarkers;

fclose(fid);

end

