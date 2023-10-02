function [gradSumRate,SumRate] = computegradSumRate(nAPs,nTx,nUsers,Gammaa,BETAA,Phii_cf,Pd,mytheta)
Gammaa = sqrt(Gammaa);
% gradSumRate = zeros(nAPs*nUsers,1);
gradSumRate = zeros(nAPs,nUsers);
SumRate = 0;
for iUser =1:nUsers
    c1 = Gammaa(:,iUser)'*mytheta(:,iUser);
    sig = (Pd)*(nTx)^2*(c1)^2; 
    gradsig = 2*((Pd)*(nTx)^2)*c1*Gammaa(:,iUser);
    gradinterference = zeros(nAPs,nUsers);
    
    interference = 0;
    for jUser = 1:nUsers
        
        % the second sum in the intererence term
        
        interference = interference + (Pd)*(nTx)*((BETAA(:,iUser))'*(mytheta(:,jUser)).^2);
        
        gradinterference(:,jUser) = gradinterference(:,jUser)+2*(Pd)*(nTx)*(BETAA(:,iUser)).*...
            mytheta(:,jUser); % why /nTx?
        
        if(jUser~=iUser)
            gamma_tilde = abs((Phii_cf(:,iUser)'*Phii_cf(:,jUser)))...
                *(Gammaa(:,jUser))./BETAA(:,jUser)...
                .*BETAA(:,iUser);
            c = gamma_tilde'*mytheta(:,jUser);
            
            gradinterference(:,jUser) = gradinterference(:,jUser)+...
                2*(Pd)*(nTx)^2*(c)*gamma_tilde;
            interference = interference + (Pd)*(nTx)^2*(c)^2;
        end
    end    
    
    mubar = interference+1;
    mu = sig+mubar;
    
    % gradient of the signal part
    gradSumRate(:,iUser) = gradSumRate(:,iUser)+ gradsig/mu;
    gradSumRate = gradSumRate + ((gradinterference)/mu-gradinterference/mubar);

SumRate = SumRate+log2(1+sig/(interference+1));

end
gradSumRate = gradSumRate*log2(exp(1));
end

