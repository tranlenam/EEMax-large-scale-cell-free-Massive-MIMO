function [v,stepsize_theta] = linesearch_theta(xi,So,B,tau,T,M,Nm,K,Gammaa,BETAA,Phii_cf,Pd,...
    P_fix_bar,noise_p,myalpha,mybeta,myrho,mytheta_now,mytheta_prev,v_prev)

maxLSiteratios = 20; % maximum number of iterations for line search

s = v_prev-mytheta_prev; s=s(:);
r = gradEnergyEfficiency(xi,So,B,tau,T,M,Nm,K,Gammaa,BETAA,Phii_cf,Pd,v_prev,...
    P_fix_bar,Pd*noise_p*Nm,myalpha) - gradEnergyEfficiency(xi,So,B,tau,T,M,Nm,K,Gammaa,BETAA,Phii_cf,Pd,mytheta_prev,...
    P_fix_bar,Pd*noise_p*Nm,myalpha);
r = r(:);

stepsize_theta = abs(max((s'*s)/(s'*r),(s'*r)/(r'*r)));

Fmytheta  = computeEnergyEfficiency(xi,So,B,tau,T,M,Nm,K,(Gammaa),BETAA,Phii_cf,Pd,mytheta_now,...
    P_fix_bar,Pd*noise_p*Nm,myalpha);
gradFmytheta = gradEnergyEfficiency(xi,So,B,tau,T,M,Nm,K,Gammaa,BETAA,Phii_cf,Pd,mytheta_now,...
    P_fix_bar,Pd*noise_p*Nm,myalpha);

for iLineSearch = 1:maxLSiteratios
    %         iLineSearch
    v = proj_Euclidean_ball_posorthant(mytheta_now+stepsize_theta*gradFmytheta,Nm,K,M);
    Fv = computeEnergyEfficiency(xi,So,B,tau,T,M,Nm,K,(Gammaa),BETAA,Phii_cf,Pd,v,...
        P_fix_bar,Pd*noise_p*Nm,myalpha);
    if(Fv>=(Fmytheta+myrho*norm(v-mytheta_now)^2))
        break
    else
        stepsize_theta = stepsize_theta*mybeta;
    end
end
end

