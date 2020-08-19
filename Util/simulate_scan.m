function [scan] = simulate_scan(Shape, Camera)
% Generate a synthetic scan of a given shape from a fixed camera position
% Input arguments:
%   Shape:  input shape
%   Camera: the camera configuration
% Output argument:
%   meshPoints: coordinates on the original mesh, the last row stores the
%               2D saliency score

% Rectify the camera
Shape.faceNors = compute_face_normal(Shape);
axis_z = Camera.origin  - Camera.lookAt;
axis_z = axis_z/norm(axis_z);
Camera.upVec = Camera.upVec - axis_z*(axis_z'*Camera.upVec);
axis_y = Camera.upVec;
axis_x = cross(axis_y, axis_z);

% Perform ray tracing to find the 3D mesh points
Height = Camera.nHeight_inner;
Width = Camera.nWidth_inner;
%
cols = kron(1:Width, ones(1,Height));
rows = kron(ones(1,Width), 1:Height);
coordX = ((2*cols-1) - Width)/min(Height, Width);
coordY = (Height - (2*rows-1))/min(Height, Width);
%


points = axis_x*coordX*Camera.scale + axis_y*coordY*Camera.scale +...
    Camera.lookAt*ones(1, Width*Height);

meshPoints = unproject(double(Shape.vertexPoss), double(Shape.faceVIds),...
    Camera.origin, points);

% Perform ray triangle intersection to obtain the 3d locations
ids = find(meshPoints(4,:) < 1e5);
meshPoints = meshPoints(:, ids);
scan.points = meshpoints2surfels(Shape, meshPoints);
scan.points = [scan.points; ones(3, size(scan.points,2))/2];

function [surfels] = meshpoints2surfels(Shape, meshPoints)
%
faceids = meshPoints(1,:);
surfelNors = Shape.faceNors(:, faceids);
%
v1_pos = Shape.vertexPoss(:, Shape.faceVIds(1, faceids));
v2_pos = Shape.vertexPoss(:, Shape.faceVIds(2, faceids));
v3_pos = Shape.vertexPoss(:, Shape.faceVIds(3, faceids));
tus = meshPoints(2, :);
tvs = meshPoints(3, :);
surfelPoss = v1_pos.*(ones(3,1)*tus) + v2_pos.*(ones(3,1)*tvs) + v3_pos.*(ones(3,1)*(1-tus-tvs));
surfels = [surfelPoss;surfelNors]; 

function [faceNors] = compute_face_normal(Shape)
%
v1_pos = Shape.vertexPoss(:, Shape.faceVIds(1,:));
v2_pos = Shape.vertexPoss(:, Shape.faceVIds(2,:));
v3_pos = Shape.vertexPoss(:, Shape.faceVIds(3,:));
e12 = v1_pos - v2_pos;
e13 = v1_pos - v3_pos;
faceNors = cross(e12, e13);
norms = sqrt(sum(faceNors.*faceNors));
faceNors = faceNors./(ones(3,1)*norms);
