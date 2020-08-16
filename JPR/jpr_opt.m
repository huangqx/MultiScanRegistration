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
% Output arguments:
%       'poses_opt': optimized scan poses 
%       'NNStruct': a data structure that encodes the nearest neighbors
%                   after scan registration