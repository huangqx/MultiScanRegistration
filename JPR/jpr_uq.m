function [CovM] = jpr_uq(scans, poses_opt, NNStruct, Para)
% Performs uncertainty quantification under the formulation of minimizing
% distances between pairs of overlaping scans
% Input parameters:
%       'scans': the input scans, where each scan is a collection of
%                surfels (positions + normals) 
%       'poses_opt': optimized scan poses 
%       'NNStruct': a data structure that encodes the nearest neighbors
%                   after scan registration
% Output parameters:
%        'CovM': A convariance matrix data structure that encodes 
%                the predicted covariance matrix. Note that for large-scale
%                problems, we only store the diagonal blocks and the
%                leading eigenvectors

%Transform the input scans
num_scans = length(scans);
for scanid = 1 : num_scans
    t = poses_opt(1:3, scanid);
    R = reshape(poses_opt(4:12, scanid), [3,3]);
    scans{scanid}.points(1:3,:) = R*scans{scanid}.points(1:3,:)...
        + t*ones(1, size(scans{scanid}.points,2));
    scans{scanid}.points(4:6,:) = R*scans{scanid}.points(4:6,:);
end
%
numcorres = size(NNStruct.corres, 2);
%
rowsJ = ones(12,1)*(1:numcorres);
colsJ = [6*ones(6,1)*NNStruct.corres(1,:) + [-5:0]'*ones(1,numcorres);
    6*ones(6,1)*NNStruct.corres(2,:) + [-5:0]'*ones(1,numcorres)];
valsJ = zeros(12, numcorres);
sumSqrDis = 0;
sumW = 0;
for cid = 1 : numcorres
    sid_surf = NNStruct.corres(1, cid);
    tid_surf = NNStruct.corres(2, cid);
    sid_pt = NNStruct.corres(3, cid);
    tid_pt = NNStruct.corres(4, cid);
    w = sqrt(NNStruct.corres(5, cid));
    spos = scans{sid_surf}.points(1:3, sid_pt);
    snor = scans{sid_surf}.points(4:6, sid_pt);
    tpos = scans{tid_surf}.points(1:3, tid_pt);
    tnor = scans{tid_surf}.points(4:6, tid_pt);
    dis = (spos-tpos)'*tnor;
    sumSqrDis = sumSqrDis + (dis*w)*(dis*w);
    sumW = sumW + w*w;
    valsJ(1:3, cid) = tnor*w;
    valsJ(4:6, cid) = cross(spos, tnor)*w;
    valsJ(7:9, cid) = -tnor*w;
    valsJ(10:12,cid) = -cross(tpos,tnor)*w;
end
sqrSigma = sumSqrDis/sumW;
dim = 6*num_scans;
J = sparse(rowsJ, colsJ, valsJ, numcorres, dim);
H = J'*J;
H = H(7:dim, 7:dim);
Cov_inv = inv(full(H))*sqrSigma;

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