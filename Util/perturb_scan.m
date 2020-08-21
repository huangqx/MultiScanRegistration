function [scan_output, R_gt, t_gt] = perturb_scan(scan_input, paras)
% Perturb the position and the normal of each input scan
eps_pos = paras(3); % noise level for perturbing the point positions
eps_nor = paras(4); % noise level for perturbing the point normals
eps_barc = paras(1); % noise level for perturbing the scan translation
eps_c = paras(2); % noise level for perturbing the scan rotation
numpoints = size(scan_input.points, 2);
% Perturbing the scan along the normal direction
dpos = (scan_input.points(4:6,:).*(ones(3,1)*(2*rand(1, numpoints)-1)))*eps_pos;
% Perturbing surface normals
dnor = (2*rand(3, numpoints)-1)*eps_nor;
scan_output = scan_input;
scan_output.points(1:3, :) = scan_output.points(1:3, :) + dpos;
scan_output.points(4:6, :) = scan_output.points(4:6, :) + dnor;
norms = sqrt(sum(scan_output.points(4:6, :).*scan_output.points(4:6, :)));
scan_output.points(4:6, :) = scan_output.points(4:6, :)./(ones(3,1)*norms);
t = (2*rand(3,1)-1)*eps_barc;
c = (2*rand(3,1)-1)*eps_c;
R = expm([0 -c(3) c(2);
    c(3) 0 -c(1);
    -c(2) c(1) 0]);
scan_output.points(1:3,:) = R*scan_output.points(1:3,:) + t*ones(1, numpoints);
scan_output.points(4:6,:) = R*scan_output.points(4:6,:);
R_gt = R';
t_gt = -R'*t;
