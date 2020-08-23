function [poses_opt, CovM] = jpr_main(scans, Para)
% This function performs multi-scan registration among a collection of
% roughly aligned scans under the formulation of minimizng distances
% between overlaping scans. The details can be found in the function 'jpr_opt'
% Input parameters:
%       'scans':  a cell that encodes the input scans. Each scan is a 
%                 collection of points, where each point has a position
%                 and normal
%        'Para':  parameters that specify the joint pairwise registration
%                 procedure
% Output parameters:
%   'poses_opt': The optimized poses after performing scan alignment
%   'VariancePred.pointwise_diag': Pointwise variance derived from the
%      diagonal blocks of a predicted covariance matrix. This is a cell
%      whose size is the same as the number of scans, and each cell is a
%      vector whose size is the same as the corresponding scan
%   'VariancePred.pointwise_eigen': Pointwise variance derived from a
%      leading eigenvector of a predicted covariance matrix. The default
%      eigenvector is the first eigenvector, but you can configurate the
%      parameter to change that. Its size is similar to that of
%      'VariancePred.pointwise_diag' Note that for visualization purpose,
%      one can further dump the predicted pointwise variance onto a
%      reconstructed mesh

% Step 1: Performs simultaneous registration and reconstruction
[poses_opt, NNStruct] = jpr_opt(scans, Para);

% Step 2: Extract the predicted convariance matrix
CovM = jpr_uq(scans, poses_opt, NNStruct, Para);

% Step 3: Extract point-wise variance
for scanid = 1 : length(scans)
    VariancePred.pointwise_diag{scanid} = extract_pvar(scans{scanid},...
        CovM.Diag{scanid});
    VariancePred.pointwise_eigen{scanid} = extract_pvarII(scans{scanid},...
        CovM.eigenVecs((6*scanid-5):(6*scanid-3),Para.uq_focused_eigenId),...
        CovM.eigenVecs((6*scanid-2):(6*scanid),Para.uq_focused_eigenId));
end
