function [CovM] = jpr_uq(scans, poses_opt, NNStruct)
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