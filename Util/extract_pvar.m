function [pointwise_var] = extract_pvar(scan, C_diag)
% extract pointwise variance based on scan.points and the predicted
% diagonal 6x6 matrix C_diag
C_diag = full(C_diag);
pointPoss = scan.points(1:3, :);
%
p = pointPoss(:, 1);
c1 = trace(C_diag(1:3,1:3));
v1 = [C_diag(6,2)-C_diag(5,3),C_diag(4,3)-C_diag(6,1),C_diag(5,1)-C_diag(4,2)]';
C2 = C_diag(4:6,4:6);
c2 = trace(C2);
%
pointwise_var = sqrt(c1 + sum(pointPoss.*pointPoss)*c2 - sum(pointPoss.*(C2*pointPoss))...
    -2*v1'*pointPoss);

