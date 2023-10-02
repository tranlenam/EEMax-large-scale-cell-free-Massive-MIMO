function [out] = computeSumRate(nAPs,nTx,nUsers,Gammaa,BETAA,Phii_cf,Pd,mytheta)
%computeSumRate Compute the achievable sum rate
% mytheta is an nAPx x nUsers matrix
% mytheta(m,k) is the power control coefficient between the m-th AP m and the k-th user 
% NOTE \mytheta(m,)=\sqrt{\eta(m,k)*\bar{gamma(m,k)}}
%   Detailed explanation goes here
out = 0;

% Gammaa means sqrt(Gammaa)in the paper
%mytheta = reshape(mytheta,nUsers,[])';
for iUser =1:nUsers
    sig = (Pd)*(nTx^2)*((Gammaa(:,iUser))'*mytheta(:,iUser))^2; 
    interference = computeInterference(nAPs,nTx,nUsers,mytheta,Gammaa,BETAA,Phii_cf,iUser,Pd);
    out = out+log2(1+sig/(interference+1));
end

end

