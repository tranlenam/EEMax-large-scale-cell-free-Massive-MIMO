function [out] = gradEnergyEfficiency(xi,So,B,tau,T,nAPs,nTx,nUsers,Gammaa,BETAA,Phii_cf,Pd,mytheta,P_fix_bar,c,myalpha)
%UNTITLED2 Summary of this function goes here

[gradRate,~] = gradPenalty(tau,T,nAPs,nTx,nUsers,Gammaa,BETAA,Phii_cf,Pd,mytheta,So);

[gradSumRate,SumRate]=computegradSumRate(nAPs,nTx,nUsers,Gammaa,...
    BETAA,Phii_cf,Pd,mytheta);

out = B*(1-tau/T)*(computeTotalPower(P_fix_bar,c,myalpha,mytheta,nUsers)...
    *gradSumRate -SumRate*gradtotalpower(c,myalpha,mytheta))...
    /(computeTotalPower(P_fix_bar,c,myalpha,mytheta,nUsers)^2)- xi*gradRate;

end

