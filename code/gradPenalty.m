function [gradSumRate,SumRate] = gradPenalty(tau,T,nAPs,nTx,nUsers,Gammaa,BETAA,Phii_cf,Pd,mytheta,So)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%So is Kx1 matrix of SE constraints
% gradSumRate = (MxK)
% SumRate (real value)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Gammaa = sqrt(Gammaa); % Sqrt if gamma (MxK)
gradSumRate = zeros(nAPs,nUsers); %MxK
SumRate = 0;
ak = zeros(nUsers,1); % a parameter
for iUser =1:nUsers  %loop all users
    ak(iUser,1) = sqrt((2^(So(iUser,1)*T/(T-tau)) -1)/(Pd*nTx^2));   
    c1 = Gammaa(:,iUser)'*mytheta(:,iUser); % The second part
%     sig = (Pd)*(nTx*c1)^2; % Numerator of Eq. 3

    gradRate = zeros(nAPs,nUsers); %MxK; reset gradRate
    gradinterference = zeros(nAPs,nUsers);
    interference = 0;
    for jUser = 1:nUsers
        % the first sum in the intererence term
        gradinterference(:,jUser) = gradinterference(:,jUser)+2*Pd*nTx*(BETAA(:,iUser)).*...
            mytheta(:,jUser);
        interference = interference + Pd*nTx*((BETAA(:,iUser))'*(mytheta(:,jUser)).^2);
        if(jUser~=iUser)
            gamma_tilde = abs((Phii_cf(:,iUser)'*Phii_cf(:,jUser)))...
                *(Gammaa(:,jUser))./BETAA(:,jUser)...
                .*BETAA(:,iUser);  %Eq. 20
            c = gamma_tilde'*mytheta(:,jUser);
        
            gradinterference(:,jUser) = gradinterference(:,jUser)+...
            2*Pd*(nTx^2)*(c)*gamma_tilde;
            interference = interference+ Pd*(nTx*c)^2;
        end      
    end

    % gradient of rate
    gradRate = gradRate + ak(iUser,1)*gradinterference/(2*sqrt(interference + 1));
    gradRate =  gradRate  - Gammaa(:,iUser);

    gradSumRate = gradSumRate - 2*(max(0,(ak(iUser,1)*sqrt(interference + 1) - c1)))*gradRate;
    % sum of rate's gradient  
    SumRate = SumRate + (max(0,(ak(iUser,1)*sqrt(interference + 1) - c1)))^2;
end

end