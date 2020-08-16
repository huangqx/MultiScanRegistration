function [poses_opt, CovM] = msg_main(scans, Para)
% The main function of this repository. It aligns a collection of partially
% overlapping scans. The output consists of optimized scan poses and the
% predicted covariance matrix among the scan poses. The predicted covariance 
% matrix utilizes the local parameterizion defined around the optimized 
% scan poses
% Input arguments:
%       'scan': input scans, where each scan consists of a collection of
%               surfels, and each surfel is given by a position + normal
%       'Para': parameters used in scan registration. Please refer to the
%               individual functions for details.
%       'Para.using_jpr': If Para.using_jpr == 1, then we perform
%               multi-scan registration by minimizing distances between
%               pairs of overlapping scans. Otherwise, we use the approach
%               of simultaneous registration and reconstruction which
%               minimizes the distances between the input scans and a
%               dynamic latent surface
if Para.using_jpr == 1
    [poses_opt, CovM] = jpr_main(scans, Para);
else
    [poses_opt, CovM] = srar_main(scans, Para);
end