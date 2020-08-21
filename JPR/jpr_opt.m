function [poses_opt, NNStruct] = jpr_opt(scans, Para)
% This function performs joint registration of a collection of scans by
% minimizing a columulative sum of distances between overlapping scans
% \min_{\{T_i, 1\leq i \leq N\}} \sum\limits_{(i,j)\in \set{E}} d^2(T_i(S_i), T_j(S_j))
% subject to T_1 = (I_3, 0)
% The distance function is defined as the sum of weighted squared
% point-2-plane distances
% d^2(T_i(S_i), T_j(S_j))
%:= \sum_{k} w_{ijk} ((T_j^{-1}T_i(p_{ijk}) - q_{ijk})^T n_{ijk})^2
% where the correspondence weight is computed using reweighted scheme
% For more details, please refer to the survey paper by
%
% 'Registration of 3D Point Clouds and Meshes: A Survey From Rigid to
% Non-Rigid' Gary K.L. Tam1, Zhi-Quan Cheng, Yu-Kun Lai, Frank C. Langbein,
% Yonghuai Liu, David Marshall,Ralph R. Martin, Xian-Fang Sun,
% and Paul L. Rosin. PAMI'13
% Note that the entire registration algorithm is implemented in C++, and
% this matlab file is a wrapper
%
% Input arguments:
%       'scans': input scans, where each scan is given by a collection of
%                surfels (i.e., position + normal)
%       'Para': includes all the parameters of joint pairwise registration
%       'Para.jpr_max_distance': The maximum distance between corresponding
%                                points of each correspondence
%       'Para.jpr_overlap_ratio': The minimum percentage of shared
%                                 correspondences for each pair of scans
%                                 to be considered overlapping
%       'Para.jpr_weight_normal': The contribution of normal when computing
%                                 closest points
%       'Para.jpr_weight_color':  The contribution of color when computing
%                                 closest points
%       'Para.jpr_down_sampling_rate': We randomly select every this number
%                                      for a sample point
%       'Para.jpr_num_reweighting_iters': Number of iterations for
%                                         calculating closest points
%       'Para.jpr_num_gauss_newton_iters': Number of Gauss Newton
%                                          iterations when fixing the correspondences
%       'Para.jpr_weight_point2planeDis': the weight for the point-2-plane
%                                         disance
% Output arguments:
%       'poses_opt': optimized scan poses 
%       'NNStruct': a data structure that encodes the nearest neighbors
%                   after scan registration
paras = [Para.jpr_max_distance, Para.jpr_overlap_ratio,...
    Para.jpr_weight_normal, Para.jpr_weight_color,...
    Para.jpr_down_sampling_rate, Para.jpr_num_reweighting_iters,...
    Para.jpr_num_gauss_newton_iters, Para.jpr_weight_point2planeDis];
%
numscans = length(scans);
numpoints = 0;
point_offsets = zeros(1, numscans+1);
for id = 1 : numscans
    numpoints_new = size(scans{id}.points, 2);
    numpoints = numpoints + numpoints_new;
    point_offsets(id+1) = numpoints;
end
points_all = zeros(9, numpoints);
numpoints = 0;
for id = 1 : numscans
    numpoints_new = size(scans{id}.points, 2);
    points_all(:, (numpoints+1):(numpoints+numpoints_new)) = scans{id}.points;
    numpoints = numpoints + numpoints_new;   
end
%
[poses_opt, NNStruct.corres] = jpr_interface(points_all, point_offsets, paras);