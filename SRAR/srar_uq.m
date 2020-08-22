function [CovM] = srar_uq(scans, poses_opt, LatentSurf, Para)
% This function performs uncertainty quantification under the formulation
% of simultaneous registration and reconstruction
% Input parameters:
%       'scans': input scans, where each scan is given by a collection of
%                surfels (i.e., position + normal)
%       'poses_opt': optimized scan poses 
%       'LatentSurf': a datastructure that encodes the latent surface and
%                     optimized correspondences between the input scans 
%                     and the latent surface
%       'Para.uq_numEigens': Number of leading eigenvectors we extract from
%                         the predicted covariance matrix
% Output parameters:
%        'CovM': A convariance matrix data structure that encodes 
%                the predicted covariance matrix. Note that for large-scale
%                problems, we only store the diagonal blocks and the
%                leading eigenvectors

% LatentSurf.surfels: The optimized surfels that encode this latent surface
%                     where each surfel is given by a position and a normal
% LatentSurf.corres:  This is a cell array whose size is identical to the
%                      number of scans
% LatentSurf.corres{scan_id}: is an array of 2 x numPoints
%                                where the first element of each row is the
%                                point index, and the corresponding second
%                                element is the index of the corresponding
%                                surfel of the latent surface

% Step 1: Pre-compute all the relevant points and store them in a data
% array
numpoints_total = size(LatentSurf.corres, 2);
datapoints_opt = zeros(3, numpoints_total); % global variable
scanindices = LatentSurf.corres(1,:); % global variable
surfelindices = LatentSurf.corres(3,:); %global variable
%
ids = 1:(numpoints_total-1);
offsets = find(LatentSurf.corres(1, ids)~=LatentSurf.corres(1, ids+1));
offsets = [0, offsets, numpoints_total];
%
for scan_id = 1 : length(scans)
    t = poses_opt(1:3, scan_id);
    R = reshape(poses_opt(4:12, scan_id), [3,3]);
    left = offsets(scan_id)+1;
    right = offsets(scan_id+1);
    pointids = LatentSurf.corres(2, left:right);
    pointPoss = scans{scan_id}.points(1:3, pointids);
    % Transform the points
    datapoints_opt(:, left:right) =...
        R*pointPoss + t*ones(1, right-left+1);
end

% Step 2: predict the sample variance

footPoss = LatentSurf.surfels(1:3, surfelindices); %global variable
footNors = LatentSurf.surfels(4:6, surfelindices); %global variable
signedDist = sum((datapoints_opt - footPoss).*footNors);
sigma = sqrt(mean(signedDist.*signedDist));

% Step 3: form the Jacobi matrix for simultaneous registration and
% reconstruction
J_rows = ones(9,1)*(1:numpoints_total);
num_scans = length(scans);
num_surfels = size(LatentSurf.surfels, 2);
J_cols = zeros(9, numpoints_total);
J_cols(1:6, :) = 6*ones(6,1)*scanindices + (-5:0)'*ones(1, numpoints_total);
J_cols(7:9, :) = 6*num_scans + (3*ones(3,1)*surfelindices...
    + (-2:0)'*ones(1, numpoints_total));

% Buffer the values of the elements
J_vals = -ones(9, numpoints_total);
J_vals(1:3, :) = footNors;
J_vals(4, :) = datapoints_opt(2,:).*footNors(3,:)...
    - datapoints_opt(3,:).*footNors(2,:);
J_vals(5, :) = datapoints_opt(3,:).*footNors(1,:)...
    - datapoints_opt(1,:).*footNors(3,:);
J_vals(6, :) = datapoints_opt(1,:).*footNors(2,:)...
    - datapoints_opt(2,:).*footNors(1,:);
% tangent directions 
[footTansI, footTansII] = coordinateFrame(footNors);
% tangent I
J_vals(7, :) = sum(footTansI.*datapoints_opt);
% tangent II
J_vals(8, :) = sum(footTansII.*datapoints_opt);
%
dim = 6*num_scans + 3*num_surfels;
J = sparse(J_rows, J_cols, J_vals, numpoints_total, dim);
H = J'*J;
%
ids0 = 7:(6*num_scans);
ids1 = (6*num_scans+1):dim;
A = H(ids0, ids0);
B = H(ids0, ids1);
C = H(ids1, ids1);
% invert the matrix C (Note that the C++ implementation is faster)
for surfelid = 1 : num_surfels
    ids = (3*surfelid-2):(3*surfelid);
    C(ids,ids) = inv(C(ids,ids));
end
% Predict the covariance matrix
Cov_inv = inv(A - B*C*B')*sigma*sigma;

% Extract the diagonal blocks which are stable
CovM.Diag{1} = zeros(6,6);
for id = 2:num_scans
    offs = (6*id-11):(6*id-6);
    CovM.Diag{id} = Cov_inv(offs, offs);
end

% Extract the leading eigen-vectors 
[eigenVecs, eigenVals] = eig(full(Cov_inv));
[s,order] = sort(-diag(eigenVals));
ids = order((length(order)-Para.uq_numEigens+1):length(order));
CovM.eigenVecs = [zeros(6, Para.uq_numEigens);eigenVecs(:, ids)];
CovM.eigenVals = sqrt(diag(eigenVals(ids,ids)));

function [footTansI, footTansII] = coordinateFrame(footNors)
% Complete the coordinate frame 
num_surfels = size(footNors, 2);
footTansI = ones(3, num_surfels);
[s,ids] = max(abs(footNors));
footTansI((0:3:(3*num_surfels-3))+ids) = 0;
inner = sum(footTansI.*footNors);
footTansI = footTansI - (ones(3,1)*inner).*footNors;
norms = sqrt(sum(footTansI.*footTansI));
footTansI = footTansI./(ones(3,1)*norms);
footTansII = cross(footNors, footTansI);
