function [out] = computeEnergyEfficiency(xi,So,B,tau,T,nAPs,nTx,nUsers,Gammaa,BETAA,Phii_cf,Pd,mytheta,P_fix_bar,c,myalpha)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[~,SumRate] = gradPenalty(tau,T,nAPs,nTx,nUsers,Gammaa,BETAA,Phii_cf,Pd,mytheta,So);
% replace Gammaa by sqrt(Gammaa);
out = B*(1-tau/T)*computeSumRate(nAPs,nTx,nUsers,sqrt(Gammaa),BETAA,Phii_cf,Pd,mytheta)...
    /computeTotalPower(P_fix_bar,c,myalpha,mytheta,nUsers)- xi*SumRate;

end

