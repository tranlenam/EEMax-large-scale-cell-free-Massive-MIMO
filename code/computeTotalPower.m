function out = computeTotalPower(P_fix_bar,c,myalpha,mytheta,K)
%computeTotalPower computes the total power
%   Detailed explanation goes here
out =  P_fix_bar+c*norm(sqrt(myalpha).*mytheta,'fro')^2;
end

