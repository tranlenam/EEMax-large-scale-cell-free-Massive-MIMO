function [z,stepsize_y] = linesearch_y(xi,So,B,tau,T,M,Nm,K,Gammaa,BETAA,Phii_cf,Pd,...
    P_fix_bar,noise_p,myalpha,mybeta,myrho,y,y_prev,z_prev)

maxLSiteratios = 20; % maximum number of iterations for line search
z = z_prev;
s = z-y_prev; s=s(:);
r = gradEnergyEfficiency(xi,So,B,tau,T,M,Nm,K,Gammaa,BETAA,Phii_cf,Pd,z,...
    P_fix_bar,Pd*noise_p*Nm,myalpha) - gradEnergyEfficiency(xi,So,B,tau,T,M,Nm,K,Gammaa,BETAA,Phii_cf,Pd,y_prev,...
    P_fix_bar,Pd*noise_p*Nm,myalpha);
r = r(:);
stepsize_y = abs( max((s'*s)/(s'*r),(s'*r)/(r'*r)));

Fy = computeEnergyEfficiency(xi,So,B,tau,T,M,Nm,K,(Gammaa),BETAA,Phii_cf,Pd,y,...
    P_fix_bar,Pd*noise_p*Nm,myalpha); % compute the objective at y
gradFy = gradEnergyEfficiency(xi,So,B,tau,T,M,Nm,K,Gammaa,BETAA,Phii_cf,Pd,y,...
    P_fix_bar,Pd*noise_p*Nm,myalpha); % compute the gradient at y

for iLineSearch = 1:maxLSiteratios
    z = proj_Euclidean_ball_posorthant(y+stepsize_y*gradFy,Nm,K,M); % projected gradient step from y
    Fz = computeEnergyEfficiency(xi,So,B,tau,T,M,Nm,K,(Gammaa),BETAA,Phii_cf,Pd,z,...
        P_fix_bar,Pd*noise_p*Nm,myalpha);
    
    if(Fz>=(Fy+myrho*norm(z-y)^2)) % if an improved solution is obtained
        break
    else
        stepsize_y = stepsize_y*mybeta; % if not, reduce step size
    end
    
end

end

