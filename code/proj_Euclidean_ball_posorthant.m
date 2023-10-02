function [out] = proj_Euclidean_ball_posorthant(x,N,nUsers,nAPs)
%proj_Euclidean_ball_posorthant Project x onto the feasible set,
% nAPs : the number of APs
% nUsers : the number of users
% x : nAPs x  nUsers matrix of power control coefficients
% x(m,k) is associated with the mth AP and the kth user
% N : number of transmit antennas per AP

%   Detailed explanation goes here
out = zeros(nAPs,nUsers);
for iAP=1:nAPs
    out(iAP,:) = 1/sqrt(N)/max(norm(max(x(iAP,:),0)),...
        1/sqrt(N))*max(x(iAP,:),0);
end
end

