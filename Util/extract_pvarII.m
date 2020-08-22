function [pointwise_var] = extract_pvarII(scan, bar_c, c)
% Compute the pointwise variance derived from a leading eigenvectors of the
% predicted covariance matrix
numpoints = size(scan.points, 2);
dev_vec = bar_c*ones(1, numpoints)...
    + cross(c*ones(1,numpoints), scan.points(1:3,:));
pointwise_var = sqrt(sum(dev_vec.*dev_vec));