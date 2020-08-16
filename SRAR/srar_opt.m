function [poses_opt, LatentSurf] = srar_opt(scans, Para)
% This function performs joint registration of a collection of scans by
% Minimizing their distances to a dynamic latent surface M:
% \min_{\{T_i, 1\leq i \leq N\}, M} \sum\limits_{i=1}^N d^2(T_i(S_i), M)
% subject to T_1 = (I_3, 0)

% This paper defines the dynamic surface as a collection of moving surfels
% and use the point-2-plane distance measure to align the each input scan
% with the dynamic latent surface M
% d^2(T_i(S_i), M)
%:= \sum_{k} (T_i(p_{ik}) - q_{c_{ij}})^T n_{k_{c_{ij}})^2
% where $c_{ij}$ is the index of the corresponding point of $p_{ik}$ on the
% latent surface. For more details, please refer to the following paper:
%
% 'High Quality Pose Estimation by Aligning Multiple Scans to a Latent Map'
% Qixing Huang and Dragomir Anguelov. ICRA 2010.
% 
% Note that the entire registration algorithm is implemented in C++, and
% this matlab file is a wrapper

% Input arguments:
%       'scans': input scans, where each scan is given by a collection of
%                surfels (i.e., position + normal)
%       'Para': includes all the parameters of joint pairwise registration
% Output arguments:
%       'poses_opt': optimized scan poses 
%       'LatentSurf': a datastructure that encodes the latent surface and
%                     optimized correspondences between the input scans 
%                     and the latent surface