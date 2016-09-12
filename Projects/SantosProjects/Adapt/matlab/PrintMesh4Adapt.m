function PrintMesh4Adapt(md, meshname)
%print the mesh 2D for adaptivity

% open a file for writing
fid = fopen(meshname, 'w');

% print a title, followed by a blank line
%fprintf(fid, 'mesh\n\n');

% nodes data
nnodes = md.mesh.numberofvertices;
x = md.mesh.x;
y = md.mesh.y;

%elements data
nelements = md.mesh.numberofelements;
elements = md.mesh.elements;

%segments data
SizeSegmentMatrix = size(md.mesh.segments);
nsegments = SizeSegmentMatrix(1); %number of rows
segments = md.mesh.segments;

% number of nodes
fprintf(fid, '%i\n', nnodes);
% x y coords
for i = 1:nnodes
    xValue = x(i);
    yValue = y(i);
    fprintf(fid, '%.12e\t%.12e\n', xValue, yValue);
end


%number of elements
fprintf(fid, '%i\n', nelements);
%nodes ID of the elements
for i = 1:nelements
    OneElement = elements(i,:);
    fprintf(fid, '%i\t%i\t%i\n', OneElement);
end

%number of segments (1D elements on boundary)
fprintf(fid, '%i\n', nsegments);
%node n, node n+1, element ID
for i = 1:nsegments
    OneSegment = segments(i,:);
    fprintf(fid, '%i\t%i\t%i\n', OneSegment);
end

%printing elements to refine
%first test: uniform
%fprintf(fid, '%i\n', nelements);
%IDs of elements to refine
%for i = 1:nelements
%    fprintf(fid, '%i\n', i);
%end

%flagelements = FlagElementsToRefine(md);
%flagelements = find(flagelements == 1);
%nFlag = length(flagelements);
%fprintf(fid, '%i\n', nFlag);
%IDs of elements to refine
%for i = 1:nFlag
 %   fprintf(fid, '%i\n', flagelements(i));
%end

fclose(fid);


end
