function interference = computeInterference(nAPs,nTx,nUsers,mytheta,Gammaan,BETAAn,Phii_cf,index,Pd)
interference=0;
firstterm = 0;
c=(BETAAn(:,index));
for jUser =1:nUsers
    if (jUser~=index)
        gamma_tilde=abs(Phii_cf(:,index)'*Phii_cf(:,jUser))*(Gammaan(:,jUser))./BETAAn(:,jUser)...
            .*(c);
        
        firstterm = firstterm+(gamma_tilde'*mytheta(:,jUser))^2;
    end
    interference = interference+(Pd)*nTx*(c)'*(mytheta(:,jUser).^2);
end

interference = interference + (Pd)*(nTx^2)*firstterm;

end

