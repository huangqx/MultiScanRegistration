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
%        'poses_opt': The optimized poses after performing scan alignment
%        'CovM': A convariance matrix data structure that encodes 
%                the predicted covariance matrix. Note that for large-scale
%                problems, we only store the diagonal blocks and the
%                leading eigenvectors

% Step 1: Performs simultaneous registration and reconstruction
[poses_opt, NNStruct] = jpr_opt(scans, Para);

% Step 2: Extract the predicted convariance matrix
CovM = jpr_uq(scans, poses_opt, NNStruct);