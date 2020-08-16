function [CovM] = srar_uq(scans, poses_opt, LatentSurf)
% This function performs uncertainty quantification under the formulation
% of simultaneous registration and reconstruction
% Input parameters:
%       'scans': input scans, where each scan is given by a collection of
%                surfels (i.e., position + normal)
%       'poses_opt': optimized scan poses 
%       'LatentSurf': a datastructure that encodes the latent surface and
%                     optimized correspondences between the input scans 
%                     and the latent surface
% Output parameters:
%        'CovM': A convariance matrix data structure that encodes 
%                the predicted covariance matrix. Note that for large-scale
%                problems, we only store the diagonal blocks and the
%                leading eigenvectors